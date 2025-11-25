import os
from html import escape

from typing import Dict, List, Optional, Union

from qgis.PyQt import uic
from qgis.PyQt.QtCore import Qt
from qgis.PyQt.QtGui import QTextOption
from qgis.PyQt.QtWidgets import (
    QFileDialog,
    QDialog,
    QMessageBox,
    QTextBrowser,
    QHBoxLayout,
    QWidget,
    QTextEdit,
    QSplitter,
    QSizePolicy,
)
from qgis.core import QgsProject, QgsRasterLayer, QgsVectorLayer, QgsWkbTypes


FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'linkscape_dialog_base.ui'))


class LinkscapeDialog(QDialog, FORM_CLASS):
    def __init__(self, iface, parent=None):
        super(LinkscapeDialog, self).__init__(parent)
        self.setupUi(self)
        self.iface = iface
        self._help_browser: Optional[QTextBrowser] = None
        self._inject_help_panel()

        self.output_dir_line.setPlaceholderText("Select output folder…")

        self._parameters = {}
        self.available_layers = {}
        self._current_layer_type = "raster"
        self._selected_layer_id: Optional[str] = None

        self._connect_signals()
        self._configure_defaults()
        self._on_layer_type_changed(self.layer_type_combo.currentText())

    # ------------------------------------------------------------------
    # UI setup helpers
    # ------------------------------------------------------------------

    def _inject_help_panel(self) -> None:
        """Add a right-hand help panel fed by a markdown file."""
        base_layout = self.layout()
        if base_layout is None:
            return

        contents_margins = base_layout.contentsMargins()
        contents_spacing = base_layout.spacing()

        # Wrap the existing layout in a container widget so we can place it beside the help panel.
        left_container = QWidget()
        left_container.setLayout(base_layout)

        help_browser = QTextBrowser()
        help_browser.setObjectName("help_browser")
        help_browser.setOpenExternalLinks(True)
        help_browser.setMinimumWidth(220)
        help_browser.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        help_browser.setLineWrapMode(QTextEdit.WidgetWidth)
        help_browser.setWordWrapMode(QTextOption.WordWrap)
        help_browser.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        splitter = QSplitter(Qt.Horizontal)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(left_container)
        splitter.addWidget(help_browser)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)

        wrapper = QHBoxLayout()
        wrapper.setContentsMargins(contents_margins)
        wrapper.setSpacing(contents_spacing)
        wrapper.addWidget(splitter)

        self.setLayout(wrapper)
        self._help_browser = help_browser
        self._load_help_content()

    def _connect_signals(self):
        self.output_browse_button.clicked.connect(self._choose_output_dir)
        self.input_layer_combo.currentIndexChanged.connect(self._on_layer_changed)
        self.layer_type_combo.currentTextChanged.connect(self._on_layer_type_changed)
        self.temporary_output_checkbox.toggled.connect(self._on_temporary_toggled)
        self.vector_units_combo.currentTextChanged.connect(self._update_vector_units_labels)
        self.strategy_combo.currentIndexChanged.connect(self._on_strategy_changed)
        self.range_checkbox.toggled.connect(self._update_range_state)
        self.obstacle_enable_checkbox.toggled.connect(self._update_obstacle_controls)
        self.obstacle_range_checkbox.toggled.connect(self._update_obstacle_range_state)
        self.vector_obstacle_enable_checkbox.toggled.connect(self._update_vector_obstacle_controls)

    def _configure_defaults(self):
        self._vector_unit_system = "metric"
        self.temporary_output_checkbox.setChecked(False)
        self.strategy_combo.setItemData(0, "most_connectivity")
        self.strategy_combo.setItemData(1, "largest_patch")
        self.vector_units_combo.setItemData(0, "metric")
        self.vector_units_combo.setItemData(1, "imperial")
        idx = self.pixel_neighborhood_combo.findText("8")
        if idx >= 0:
            self.pixel_neighborhood_combo.setCurrentIndex(idx)
        self._on_temporary_toggled(self.temporary_output_checkbox.isChecked())
        self._update_vector_units_labels(self.vector_units_combo.currentText())
        self._update_range_state(self.range_checkbox.isChecked())
        self._update_obstacle_controls(self.obstacle_enable_checkbox.isChecked())
        self._update_vector_obstacle_controls(self.vector_obstacle_enable_checkbox.isChecked())

    def _load_help_content(self) -> None:
        """Load markdown help content into the right-hand panel."""
        if self._help_browser is None:
            return
        help_path = os.path.join(os.path.dirname(__file__), "linkscape_help.md")
        if not os.path.exists(help_path):
            self._help_browser.setPlainText("Linkscape\n\nHelp file not found.")
            return
        try:
            with open(help_path, "r", encoding="utf-8") as fh:
                text = fh.read()
            # Try to render Markdown
            try:
                import markdown  # type: ignore
                html = markdown.markdown(text, extensions=["extra"])
                self._help_browser.setHtml(html)
            except ImportError:
                # Fallback if markdown lib is missing: Use plain text (wrap mode handles long lines).
                self._help_browser.setPlainText(text)
            except Exception as e:
                # Fallback for other errors
                self._help_browser.setPlainText(f"Error loading help:\n{e}\n\n{text}")

        except Exception as exc:  # pragma: no cover - UI only
            self._help_browser.setPlainText(f"Linkscape help could not be loaded:\n{exc}")

    def populate_layers(self):
        self.available_layers.clear()
        self.input_layer_combo.clear()

        project = QgsProject.instance()
        selected_type = self._current_layer_type
        if selected_type == "vector":
            target_cls = QgsVectorLayer
        else:
            target_cls = QgsRasterLayer

        layers = [
            layer
            for layer in project.mapLayers().values()
            if isinstance(layer, target_cls) and layer.isValid()
        ]

        if not layers:
            label = "No raster layers available" if selected_type == "raster" else "No vector layers available"
            self.input_layer_combo.addItem(label, None)
            self.input_layer_combo.setEnabled(False)
            self._selected_layer_id = None
            return

        self.input_layer_combo.setEnabled(True)
        for lyr in layers:
            self.available_layers[lyr.id()] = lyr
            self.input_layer_combo.addItem(lyr.name(), lyr.id())

        active_layer = self.iface.activeLayer() if self.iface else None
        if isinstance(active_layer, target_cls) and active_layer.id() in self.available_layers:
            index = self.input_layer_combo.findData(active_layer.id())
            if index >= 0:
                self.input_layer_combo.setCurrentIndex(index)
        else:
            self.input_layer_combo.setCurrentIndex(0)

        self._on_layer_changed(self.input_layer_combo.currentIndex())
        self._populate_vector_obstacle_layers()

    def _populate_vector_obstacle_layers(self):
        combo = self.vector_obstacle_layer_combo
        combo.blockSignals(True)
        combo.clear()
        combo.addItem("None", "")
        project = QgsProject.instance()
        for layer in project.mapLayers().values():
            if not isinstance(layer, QgsVectorLayer) or not layer.isValid():
                continue
            if QgsWkbTypes.geometryType(layer.wkbType()) != QgsWkbTypes.PolygonGeometry:
                continue
            combo.addItem(layer.name(), layer.id())
        combo.blockSignals(False)

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _choose_output_dir(self):
        start_dir = self.output_dir_line.text() or os.path.expanduser("~")
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory", start_dir)
        if directory:
            self.output_dir_line.setText(directory)

    def _on_layer_changed(self, index: int):
        layer = self._layer_from_index(index)
        if layer:
            self._selected_layer_id = layer.id()
            if not self.output_dir_line.text():
                self.output_dir_line.setText(os.path.dirname(layer.source()))
            if (
                self._current_layer_type == "vector"
                and self.vector_output_name_line.text().strip() in ("", "optimized_corridors.gpkg", "linkscape_corridors.gpkg")
            ):
                base_name = f"{layer.name()}_linkscape.gpkg"
                self.vector_output_name_line.setText(base_name)

    def _on_temporary_toggled(self, checked: bool):
        self.output_dir_line.setEnabled(not checked)
        self.output_browse_button.setEnabled(not checked)
        if checked:
            self.output_dir_line.clear()

    def _update_vector_units_labels(self, _: str):
        new_units = self.vector_units_combo.currentData()
        if not hasattr(self, "_vector_unit_system"):
            self._vector_unit_system = "metric"

        if new_units == self._vector_unit_system:
            # Still update labels in case they were not initialised
            pass
        else:
            if new_units == "imperial":
                # Convert existing values from metric to imperial
                self.vector_min_corridor_width_spin.setValue(
                    round(self.vector_min_corridor_width_spin.value() * 3.280839895, 4)
                )
                self.vector_max_search_spin.setValue(
                    round(self.vector_max_search_spin.value() * 3.280839895, 4)
                )
                self.vector_min_patch_size_spin.setValue(
                    round(self.vector_min_patch_size_spin.value() * 2.471053814, 4)
                )
                self.vector_budget_spin.setValue(
                    round(self.vector_budget_spin.value() * 2.471053814, 4)
                )
                if self.vector_max_corridor_area_spin.value() != 0.0:
                    self.vector_max_corridor_area_spin.setValue(
                        round(self.vector_max_corridor_area_spin.value() * 2.471053814, 4)
                    )
            else:
                # Convert imperial to metric
                self.vector_min_corridor_width_spin.setValue(
                    round(self.vector_min_corridor_width_spin.value() * 0.3048, 4)
                )
                self.vector_max_search_spin.setValue(
                    round(self.vector_max_search_spin.value() * 0.3048, 4)
                )
                self.vector_min_patch_size_spin.setValue(
                    round(self.vector_min_patch_size_spin.value() * 0.404685642, 4)
                )
                self.vector_budget_spin.setValue(
                    round(self.vector_budget_spin.value() * 0.404685642, 4)
                )
                if self.vector_max_corridor_area_spin.value() != 0.0:
                    self.vector_max_corridor_area_spin.setValue(
                        round(self.vector_max_corridor_area_spin.value() * 0.404685642, 4)
                    )
            self._vector_unit_system = new_units

        if self.vector_units_combo.currentData() == "imperial":
            self.label_vector_min_width.setText("Min Corridor Width (ft):")
            self.label_vector_max_search.setText("Max Search Distance (ft):")
            self.label_vector_min_patch.setText("Min Patch Size (ac):")
            self.label_vector_budget.setText("Budget (ac):")
            self.label_vector_max_area.setText("Max Corridor Area (ac):")
        else:
            self.label_vector_min_width.setText("Min Corridor Width (m):")
            self.label_vector_max_search.setText("Max Search Distance (m):")
            self.label_vector_min_patch.setText("Min Patch Size (ha):")
            self.label_vector_budget.setText("Budget (ha):")
            self.label_vector_max_area.setText("Max Corridor Area (ha):")

    def _on_strategy_changed(self, _index: int):
        # Placeholder for future dynamic behaviour.
        return

    def _update_range_state(self, checked: bool):
        self.patch_value_line.setVisible(not checked)
        self.label_patch_value.setVisible(not checked)
        self.label_range_lower.setVisible(checked)
        self.range_lower_line.setVisible(checked)
        self.label_range_upper.setVisible(checked)
        self.range_upper_line.setVisible(checked)

    def _update_obstacle_controls(self, enabled: bool):
        is_enabled = bool(enabled)
        self.label_obstacle_value.setEnabled(is_enabled)
        self.obstacle_value_line.setEnabled(is_enabled)
        self.obstacle_range_checkbox.setEnabled(is_enabled)
        self._update_obstacle_range_state(self.obstacle_range_checkbox.isChecked())
        if not is_enabled:
            self.obstacle_value_line.clear()

    def _update_obstacle_range_state(self, checked: bool):
        is_enabled = self.obstacle_enable_checkbox.isChecked()
        show_range = bool(checked) and is_enabled

        self.label_obstacle_lower.setVisible(show_range)
        self.obstacle_range_lower_line.setVisible(show_range)
        self.label_obstacle_upper.setVisible(show_range)
        self.obstacle_range_upper_line.setVisible(show_range)

        self.label_obstacle_lower.setEnabled(show_range)
        self.obstacle_range_lower_line.setEnabled(show_range)
        self.label_obstacle_upper.setEnabled(show_range)
        self.obstacle_range_upper_line.setEnabled(show_range)

        value_enabled = is_enabled and not show_range
        self.label_obstacle_value.setVisible(not show_range)
        self.obstacle_value_line.setVisible(not show_range)
        self.label_obstacle_value.setEnabled(value_enabled)
        self.obstacle_value_line.setEnabled(value_enabled)

        if not is_enabled:
            self.obstacle_range_lower_line.clear()
            self.obstacle_range_upper_line.clear()

    def _update_vector_obstacle_controls(self, enabled: bool):
        is_enabled = bool(enabled)
        self.label_vector_obstacle_layer.setEnabled(is_enabled)
        self.vector_obstacle_layer_combo.setEnabled(is_enabled)
        self.label_vector_obstacle_resolution.setEnabled(is_enabled)
        self.vector_obstacle_resolution_spin.setEnabled(is_enabled)
        if not is_enabled:
            self.vector_obstacle_layer_combo.setCurrentIndex(0)

    def _on_layer_type_changed(self, text: str):
        self._current_layer_type = (text or "Raster").lower()
        self._update_group_visibility()
        self.populate_layers()
        if self._current_layer_type == "raster":
            self._update_range_state(self.range_checkbox.isChecked())

    def _update_group_visibility(self):
        is_raster = self._current_layer_type == "raster"
        self.patchGroup.setVisible(is_raster)
        self.paramGroup.setVisible(is_raster)
        self.obstacleGroup.setVisible(is_raster)
        self.vectorGroup.setVisible(not is_raster)
        self.vectorObstacleGroup.setVisible(not is_raster)

    # ------------------------------------------------------------------
    # Data extraction helpers
    # ------------------------------------------------------------------

    def _layer_from_index(self, index: int) -> Optional[Union[QgsRasterLayer, QgsVectorLayer]]:
        layer_id = self.input_layer_combo.itemData(index)
        if not layer_id:
            return None
        return self.available_layers.get(layer_id)

    def get_selected_layer(self) -> Optional[QgsRasterLayer]:
        layer_id = self._parameters.get("layer_id") or self._selected_layer_id
        if not layer_id:
            return None
        layer = QgsProject.instance().mapLayer(layer_id)
        if isinstance(layer, (QgsRasterLayer, QgsVectorLayer)) and layer.isValid():
            return layer
        return self.available_layers.get(layer_id)

    def get_layer_id(self) -> Optional[str]:
        return self._parameters.get("layer_id")

    def get_output_directory(self) -> str:
        return self._parameters.get("output_dir", "")

    def get_parameters(self) -> Dict:
        return self._parameters.get("params", {})

    def get_strategy(self) -> str:
        return self._parameters.get("strategy", "most_connectivity")

    def use_temporary_output(self) -> bool:
        return self._parameters.get("use_temporary_output", False)

    def get_layer_type(self) -> str:
        return self._parameters.get("layer_type", self._current_layer_type)

    # ------------------------------------------------------------------
    # Validation and acceptance
    # ------------------------------------------------------------------

    def _parse_values(self, text: str) -> List[float]:
        if not text.strip():
            return []
        values = []
        for token in text.split(","):
            token = token.strip()
            if not token:
                continue
            try:
                values.append(float(token))
            except ValueError as exc:
                raise ValueError(f"Invalid value '{token}' in Ecosystem Values.") from exc
        return values

    def _collect_parameters(self) -> Dict:
        params = self._collect_vector_parameters() if self._current_layer_type == "vector" else self._collect_raster_parameters()
        params["strategy"] = self.strategy_combo.currentData() or "most_connectivity"
        return params

    def _collect_raster_parameters(self) -> Dict:
        connectivity = int(self.pixel_neighborhood_combo.currentText())
        use_range = self.range_checkbox.isChecked()
        params = {
            "patch_connectivity": connectivity,
            "patch_mode": "range" if use_range else "value",
            "patch_values": [],
            "range_lower": None,
            "range_upper": None,
            "value_tolerance": 1e-6,
            "nodata_fallback": -9999.0,
            "min_patch_size": self.min_patch_size_spin.value(),
            "budget_pixels": self.budget_spin.value(),
            "max_search_distance": self.max_search_spin.value(),
            "max_corridor_area": self.max_corridor_area_spin.value() or None,
            "min_corridor_width": self.min_corridor_width_spin.value(),
            "allow_bottlenecks": self.allow_bottlenecks_checkbox.isChecked(),
            "obstacle_enabled": self.obstacle_enable_checkbox.isChecked(),
            "obstacle_mode": "value",
            "obstacle_values": [],
            "obstacle_range_lower": None,
            "obstacle_range_upper": None,
        }

        if use_range:
            lower_text = self.range_lower_line.text().strip()
            upper_text = self.range_upper_line.text().strip()
            if not lower_text or not upper_text:
                raise ValueError("Provide both lower and upper values for the range.")
            try:
                lower = float(lower_text)
                upper = float(upper_text)
            except ValueError as exc:
                raise ValueError("Range bounds must be numeric.") from exc
            if lower > upper:
                lower, upper = upper, lower
            params["range_lower"] = lower
            params["range_upper"] = upper
        else:
            values = self._parse_values(self.patch_value_line.text())
            if not values:
                raise ValueError("Enter at least one patch value.")
            params["patch_values"] = values

        if params["obstacle_enabled"]:
            obstacle_is_range = self.obstacle_range_checkbox.isChecked()
            params["obstacle_mode"] = "range" if obstacle_is_range else "value"
            if obstacle_is_range:
                lower_text = self.obstacle_range_lower_line.text().strip()
                upper_text = self.obstacle_range_upper_line.text().strip()
                if not lower_text or not upper_text:
                    raise ValueError("Provide both lower and upper values for the obstacle range.")
                try:
                    lower = float(lower_text)
                    upper = float(upper_text)
                except ValueError as exc:
                    raise ValueError("Obstacle range bounds must be numeric.") from exc
                if lower > upper:
                    lower, upper = upper, lower
                params["obstacle_range_lower"] = lower
                params["obstacle_range_upper"] = upper
            else:
                obstacle_values = self._parse_values(self.obstacle_value_line.text())
                if not obstacle_values:
                    raise ValueError("Enter at least one obstacle value.")
                params["obstacle_values"] = obstacle_values

        return params

    def _collect_vector_parameters(self) -> Dict:
        output_name = self.vector_output_name_line.text().strip() or "linkscape_corridors.gpkg"
        units = self.vector_units_combo.currentData() or "metric"

        width_value = self.vector_min_corridor_width_spin.value()
        search_value = self.vector_max_search_spin.value()
        min_patch_value = self.vector_min_patch_size_spin.value()
        budget_value = self.vector_budget_spin.value()
        max_area_value = self.vector_max_corridor_area_spin.value()
        obstacle_enabled = self.vector_obstacle_enable_checkbox.isChecked()
        obstacle_layer_id = self.vector_obstacle_layer_combo.currentData()
        resolution_value = self.vector_obstacle_resolution_spin.value()

        if units == "imperial":
            width_value *= 0.3048
            search_value *= 0.3048
            min_patch_value *= 0.404685642
            budget_value *= 0.404685642
            max_area_value = max_area_value * 0.404685642 if max_area_value else None
            resolution_value *= 0.3048
        else:
            max_area_value = max_area_value or None

        if obstacle_enabled and not obstacle_layer_id:
            raise ValueError("Select an obstacle layer to enable obstacle avoidance.")

        return {
            "min_corridor_width": width_value,
            "max_corridor_area": max_area_value,
            "min_patch_size": min_patch_value,
            "budget_area": budget_value,
            "max_search_distance": search_value,
            "output_name": output_name,
            "unit_system": units,
            "obstacle_enabled": obstacle_enabled and bool(obstacle_layer_id),
            "obstacle_layer_id": obstacle_layer_id if obstacle_enabled else None,
            "grid_resolution": resolution_value,
        }

    def accept(self):
        layer = self._layer_from_index(self.input_layer_combo.currentIndex())
        expected_cls = QgsRasterLayer if self._current_layer_type == "raster" else QgsVectorLayer
        if layer is None or not isinstance(layer, expected_cls) or not layer.isValid():
            msg = "Select a valid raster layer before running." if self._current_layer_type == "raster" else \
                "Select a valid vector layer before running."
            QMessageBox.warning(self, "Linkscape", msg)
            return

        use_temporary = self.temporary_output_checkbox.isChecked()
        output_dir = self.output_dir_line.text().strip()
        if not use_temporary and not output_dir:
            QMessageBox.warning(self, "Linkscape", "Select an output directory.")
            return

        try:
            params = self._collect_parameters()
        except ValueError as exc:
            QMessageBox.warning(self, "Linkscape", str(exc))
            return

        strategy = params.pop("strategy", "most_connectivity")

        if not use_temporary:
            try:
                os.makedirs(output_dir, exist_ok=True)
            except OSError as exc:  # noqa: PERF203
                QMessageBox.warning(
                    self,
                    "Linkscape",
                    f"Unable to create or access the output directory:\n{exc}",
                )
                return
        else:
            output_dir = ""

        self._parameters = {
            "layer_id": layer.id(),
            "layer_type": self._current_layer_type,
            "output_dir": output_dir,
            "params": params,
            "strategy": strategy,
            "use_temporary_output": use_temporary,
        }

        super(LinkscapeDialog, self).accept()
