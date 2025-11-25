import csv
import os

from typing import Any, Dict, List, Optional

from qgis.PyQt.QtCore import Qt, QCoreApplication, QUrl
from qgis.PyQt.QtGui import QIcon, QDesktopServices
from qgis.PyQt.QtWidgets import QAction, QDialog, QLabel, QMessageBox, QProgressBar
from qgis.core import QgsApplication, Qgis, QgsProject, QgsRasterLayer, QgsVectorLayer

from .processing_provider import LinkscapeProcessingProvider


class LinkscapePlugin:
    def __init__(self, iface):
        self.iface = iface
        self.plugin_dir = os.path.dirname(__file__)
        self.action = None
        self.provider = None

    def initGui(self):
        self.action = QAction("Run Linkscape", self.iface.mainWindow())
        icon_path = os.path.join(self.plugin_dir, "icon.png")
        icon = QIcon(icon_path)
        if not icon.isNull():
            self.action.setIcon(icon)
        self.action.triggered.connect(self.run)
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("&Linkscape (SORUS)", self.action)

        if self.provider is None:
            self.provider = LinkscapeProcessingProvider(self)
            QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        self.iface.removeToolBarIcon(self.action)
        self.iface.removePluginMenu("&Linkscape (SORUS)", self.action)
        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None

    def run(self):
        from .linkscape_dialog import LinkscapeDialog

        dialog = LinkscapeDialog(self.iface)
        if dialog.exec_() != QDialog.Accepted:
            return

        layer = dialog.get_selected_layer()
        layer_type = dialog.get_layer_type()
        output_dir = dialog.get_output_directory()
        params = dialog.get_parameters()
        strategy = dialog.get_strategy()
        use_temporary = dialog.use_temporary_output()

        if layer is None:
            expected_label = "raster" if layer_type == "raster" else "vector"
            QMessageBox.warning(
                self.iface.mainWindow(),
                "Linkscape",
                f"No {expected_label} layer selected for analysis.",
            )
            return
        if layer_type == "vector" or isinstance(layer, QgsVectorLayer):
            feat_count = 0
            try:
                feat_count = layer.featureCount()
            except Exception:
                pass
            if feat_count <= 1:
                QMessageBox.warning(
                    self.iface.mainWindow(),
                    "Linkscape",
                    "Vector mode requires multiple features (one per patch). "
                    "Your layer currently contains a single feature. Please split it into individual patches "
                    "using the Multipart to Singleparts tool or similar, then try again.",
                )
                return

        progress_widget = None
        progress_bar: Optional[QProgressBar] = None
        progress_label: Optional[QLabel] = None
        QCoreApplication.processEvents()

        def progress_callback(value: int, message: Optional[str] = None) -> None:
            if progress_bar is not None:
                progress_bar.setMaximum(100)
                progress_bar.setValue(max(0, min(100, int(value))))
            if message and progress_label is not None:
                progress_label.setText(message)
            QCoreApplication.processEvents()

        results = []
        def _format_number(val):
            if isinstance(val, float):
                return f"{val:,.2f}"
            return f"{val:,}"

        def _extract_metrics(result: dict) -> Dict[str, Any]:
            stats = result.get("stats", {}) or {}
            title = stats.get("layer_name", result["strategy"].replace("_", " ").title())
            area_unit = stats.get("area_units_label") or ("px" if layer_type != "vector" else "")
            connected_area = (
                stats.get("total_connected_area_display")
                or stats.get("total_connected_size")
                or stats.get("final_patch_area_display")
                or stats.get("final_patch_size")
            )
            largest_area = (
                stats.get("largest_group_area_display")
                or stats.get("largest_group_size")
                or stats.get("final_patch_size")
            )
            return {
                "output": title,
                "area_unit": area_unit,
                "corridors_created": stats.get("corridors_used"),
                "patches_connected": stats.get("patches_connected") or stats.get("patches_merged"),
                "total_contiguous_area": connected_area,
                "largest_contiguous_area": largest_area,
                "budget_used": stats.get("budget_used_display") or stats.get("budget_used"),
                "budget_total": stats.get("budget_total_display") or stats.get("budget_total"),
                "budget_units": area_unit or ("px" if layer_type != "vector" else ""),
                "note": stats.get("note"),
            }

        def _format_result_summary(result: dict) -> str:
            metrics = _extract_metrics(result)
            details: List[str] = []

            if metrics["corridors_created"] is not None:
                details.append(f"Corridors created: {_format_number(metrics['corridors_created'])}")

            if metrics["patches_connected"]:
                details.append(f"Patches connected: {_format_number(metrics['patches_connected'])}")

            area_unit = metrics["area_unit"]
            connected_area = metrics["total_contiguous_area"]
            if connected_area is not None:
                suffix = f" {area_unit}" if area_unit else ""
                details.append(f"Total contiguous area: {_format_number(connected_area)}{suffix}")

            largest_area = metrics["largest_contiguous_area"]
            if largest_area is not None:
                suffix = f" {area_unit}" if area_unit else " px"
                details.append(f"Largest contiguous area: {_format_number(largest_area)}{suffix}")

            budget_used = metrics["budget_used"]
            budget_total = metrics["budget_total"]
            budget_unit = metrics["budget_units"]
            if budget_used is not None:
                if budget_total:
                    details.append(
                        f"Budget used: {_format_number(budget_used)}/{_format_number(budget_total)}"
                        f"{f' {budget_unit}' if budget_unit else ''}"
                    )
                else:
                    details.append(
                        f"Budget used: {_format_number(budget_used)}{f' {budget_unit}' if budget_unit else ''}"
                    )

            if not details:
                return metrics["output"]
            return metrics["output"] + "\n  - " + "\n  - ".join(details)

        try:
            if layer_type == "vector" or isinstance(layer, QgsVectorLayer):
                from .analysis_vector import VectorAnalysisError, run_vector_analysis

                try:
                    progress_widget = self.iface.messageBar().createMessage("Linkscape", "Running vector analysis…")
                    progress_bar = QProgressBar()
                    progress_bar.setMaximum(100)
                    progress_bar.setValue(0)
                    progress_widget.layout().addWidget(progress_bar)
                    progress_label = progress_widget.findChild(QLabel)
                    if progress_label is not None:
                        progress_label.setText("Running vector analysis…")
                    self.iface.messageBar().pushWidget(progress_widget, Qgis.Info)

                    results = run_vector_analysis(
                        layer,
                        output_dir,
                        params,
                        strategy=strategy,
                        temporary=use_temporary,
                        iface=self.iface,
                        progress_cb=progress_callback,
                    )
                except VectorAnalysisError as exc:
                    if progress_widget:
                        self.iface.messageBar().popWidget(progress_widget)
                        progress_widget = None
                        progress_bar = None
                        progress_label = None
                    QMessageBox.critical(
                        self.iface.mainWindow(),
                        "Linkscape",
                        str(exc),
                    )
                    return
            else:
                from .analysis_raster import RasterAnalysisError, run_raster_analysis

                try:
                    progress_widget = self.iface.messageBar().createMessage("Linkscape", "Running raster analysis…")
                    progress_bar = QProgressBar()
                    progress_bar.setMaximum(100)
                    progress_bar.setValue(0)
                    progress_widget.layout().addWidget(progress_bar)
                    progress_label = progress_widget.findChild(QLabel)
                    if progress_label is not None:
                        progress_label.setText("Running raster analysis…")
                    self.iface.messageBar().pushWidget(progress_widget, Qgis.Info)

                    results = run_raster_analysis(
                        layer,
                        output_dir,
                        params,
                        strategy=strategy,
                        temporary=use_temporary,
                        iface=self.iface,
                        progress_cb=progress_callback,
                    )
                except RasterAnalysisError as exc:
                    if progress_widget:
                        self.iface.messageBar().popWidget(progress_widget)
                        progress_widget = None
                        progress_bar = None
                        progress_label = None
                    QMessageBox.critical(
                        self.iface.mainWindow(),
                        "Linkscape",
                        str(exc),
                    )
                    return
        except Exception as exc:  # noqa: BLE001
            if progress_widget:
                self.iface.messageBar().popWidget(progress_widget)
                progress_widget = None
                progress_bar = None
                progress_label = None
            QgsApplication.messageLog().logMessage(
                f"Unexpected error: {exc}",
                "Linkscape",
                level=Qgis.Critical,
            )
            QMessageBox.critical(
                self.iface.mainWindow(),
                "Linkscape",
                "An unexpected error occurred. See logs for details.",
            )
            return
        finally:
            if progress_widget:
                self.iface.messageBar().popWidget(progress_widget)
                progress_widget = None
                progress_bar = None
                progress_label = None
        if results:
            summary_lines = [_format_result_summary(result) for result in results]
            metrics_rows = [_extract_metrics(result) for result in results]
            message = "Linkscape analysis complete:\n" + "\n\n".join(summary_lines)
            if layer_type != "vector":
                message += (
                    "\n\n*Note* Value of raster pixels indicates the total area of the new contiguous "
                    "patch the corridor creates."
                )

            log_text = message.replace("Linkscape analysis complete:\n", "")
            log_dir = output_dir or os.path.dirname(layer.source())
            os.makedirs(log_dir, exist_ok=True)
            strategy_suffix = strategy.replace(" ", "").lower() if strategy else ""
            prefix = "ls"
            if strategy_suffix == "most_connectivity":
                summary_basename = f"{prefix}_mostconn_summary"
            elif strategy_suffix == "largest_patch":
                summary_basename = f"{prefix}_largestpatch_summary"
            else:
                summary_basename = f"{prefix}_summary"

            log_path = os.path.join(log_dir, f"{summary_basename}.txt")
            with open(log_path, "w", encoding="utf-8") as fh:
                fh.write("Linkscape Run Summary\n")
                fh.write("=" * 70 + "\n\n")
                fh.write(log_text.strip() + "\n")
                if layer_type != "vector":
                    fh.write(
                        "\n*Note* Value of raster pixels indicates the total area of the new contiguous "
                        "patch the corridor creates.\n"
                    )

            csv_path = os.path.join(log_dir, f"{summary_basename}.csv")
            with open(csv_path, "w", newline="", encoding="utf-8") as csvfile:
                fieldnames = [
                    "output",
                    "corridors_created",
                    "patches_connected",
                    "total_contiguous_area",
                    "largest_contiguous_area",
                    "budget_used",
                    "budget_total",
                    "budget_units",
                ]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for row in metrics_rows:
                    writer.writerow({key: row.get(key, "") for key in fieldnames})

            QMessageBox.information(
                self.iface.mainWindow(),
                "Linkscape",
                message,
            )
            if os.path.exists(csv_path):
                QgsApplication.messageLog().logMessage(
                    f"Summary written to {csv_path}", "Linkscape", Qgis.Info
                )
                csv_uri = (
                    f"file://{csv_path}?encoding=UTF-8&type=csv&detectTypes=yes&geomType=none"
                    "&subsetIndex=no&watchFile=no"
                )
                summary_layer_name = f"Linkscape Summary ({os.path.basename(csv_path)})"
                summary_layer = QgsVectorLayer(csv_uri, summary_layer_name, "delimitedtext")
                if summary_layer.isValid():
                    QgsProject.instance().addMapLayer(summary_layer)
        else:
            QMessageBox.information(
                self.iface.mainWindow(),
                "Linkscape",
                "Analysis completed, but no corridor outputs were generated.",
            )
