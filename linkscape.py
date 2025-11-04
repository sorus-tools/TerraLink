import os

from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction, QDialog, QMessageBox
from qgis.core import QgsApplication, Qgis, QgsRasterLayer, QgsVectorLayer

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

        results = []

        try:
            if layer_type == "vector" or isinstance(layer, QgsVectorLayer):
                from .analysis_vector import VectorAnalysisError, run_vector_analysis

                if isinstance(layer, QgsVectorLayer) and layer.crs().isGeographic():
                    QMessageBox.information(
                        self.iface.mainWindow(),
                        "Linkscape",
                        "Linkscape will temporarily project this layer to a local UTM CRS to ensure accurate area measurements.",
                    )

                try:
                    results = run_vector_analysis(
                        layer,
                        output_dir,
                        params,
                        strategy=strategy,
                        temporary=use_temporary,
                        iface=self.iface,
                    )
                except VectorAnalysisError as exc:
                    QMessageBox.critical(
                        self.iface.mainWindow(),
                        "Linkscape",
                        str(exc),
                    )
                    return
            else:
                from .analysis_raster import RasterAnalysisError, run_raster_analysis

                try:
                    results = run_raster_analysis(
                        layer,
                        output_dir,
                        params,
                        strategy=strategy,
                        temporary=use_temporary,
                        iface=self.iface,
                    )
                except RasterAnalysisError as exc:
                    QMessageBox.critical(
                        self.iface.mainWindow(),
                        "Linkscape",
                        str(exc),
                    )
                    return
        except Exception as exc:  # noqa: BLE001
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

        if results:
            summary_lines = [
                (
                    f"{result['stats'].get('layer_name', result['strategy'].replace('_', ' ').title())}: "
                    f"{result.get('output_path') or 'Temporary layer'}"
                )
                for result in results
            ]
            QMessageBox.information(
                self.iface.mainWindow(),
                "Linkscape",
                "Linkscape analysis complete:\n" + "\n".join(summary_lines),
            )
        else:
            QMessageBox.information(
                self.iface.mainWindow(),
                "Linkscape",
                "Analysis completed, but no corridor outputs were generated.",
            )
