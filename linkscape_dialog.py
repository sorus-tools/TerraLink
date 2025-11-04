import os

from typing import Dict, List, Optional, Union

from qgis.PyQt import uic
from qgis.PyQt.QtWidgets import QFileDialog, QDialog, QMessageBox
from qgis.core import QgsProject, QgsRasterLayer, QgsVectorLayer


FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'linkscape_dialog_base.ui'))


class LinkscapeDialog(QDialog, FORM_CLASS):
    def __init__(self, iface, parent=None):
        super(LinkscapeDialog, self).__init__(parent)
        self.setupUi(self)
        self.iface = iface

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

    def _connect_signals(self):
        self.output_browse_button.clicked.connect(self._choose_output_dir)
        self.input_layer_combo.currentIndexChanged.connect(self._on_layer_changed)
        self.layer_type_combo.currentTextChanged.connect(self._on_layer_type_changed)
        self.temporary_output_checkbox.toggled.connect(self._on_temporary_toggled)
        self.vector_units_combo.currentTextChanged.connect(self._update_vector_units_labels)
        self.strategy_combo.currentIndexChanged.connect(self._on_strategy_changed)
        self.range_checkbox.toggled.connect(self._update_range_state)

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
        self.vectorGroup.setVisible(not is_raster)

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

        return params

    def _collect_vector_parameters(self) -> Dict:
        output_name = self.vector_output_name_line.text().strip() or "linkscape_corridors.gpkg"
        units = self.vector_units_combo.currentData() or "metric"

        width_value = self.vector_min_corridor_width_spin.value()
        search_value = self.vector_max_search_spin.value()
        min_patch_value = self.vector_min_patch_size_spin.value()
        budget_value = self.vector_budget_spin.value()
        max_area_value = self.vector_max_corridor_area_spin.value()

        if units == "imperial":
            width_value *= 0.3048
            search_value *= 0.3048
            min_patch_value *= 0.404685642
            budget_value *= 0.404685642
            max_area_value = max_area_value * 0.404685642 if max_area_value else None
        else:
            max_area_value = max_area_value or None

        return {
            "min_corridor_width": width_value,
            "max_corridor_area": max_area_value,
            "min_patch_size": min_patch_value,
            "budget_area": budget_value,
            "max_search_distance": search_value,
            "output_name": output_name,
            "unit_system": units,
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
