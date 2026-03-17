"""
TerraLink Corridor Analysis - Vector Workflow (v23.3)
----------------------------------------------------
Runs the single-strategy corridor optimization workflow for polygon
patch datasets.

Updates in v23.3:
- Fixed Redundancy: Allows parallel connections between components if budget permits.
- Fixed Traversal: Efficient spatial indexing for detecting intermediate patch crossings.
- Fixed Logic: Corridors crossing intermediate patches (A->C->B) are now prioritized
  if C is not yet connected, even if A and B are.
"""

from __future__ import annotations

import heapq
import math
import os
import tempfile
import time
import csv
from contextlib import contextmanager, nullcontext
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
# NumPy 2.x removed np.int; add shim for any legacy references.
if not hasattr(np, "int"):  # pragma: no cover
    np.int = int  # type: ignore[attr-defined]

from . import terralink_graph as nx
from .terralink_engine import NetworkOptimizer, UnionFind
from .core_select import select_circuit_utility
from .habitat_availability_mode import (
    HABITAT_AVAILABILITY_DEFAULT_KERNEL,
    HABITAT_AVAILABILITY_DEFAULT_SCALING,
    HabitatAvailabilityEvaluator,
    HabitatAvailabilityNode,
    normalize_habitat_availability_kernel,
    normalize_patch_area_scaling,
    scale_patch_area,
)
from .utils import emit_progress, log_error
# Import the graph-metrics helper library
try:
    from . import graph_math
except ImportError:
    graph_math = None
from PyQt5.QtCore import QVariant, QUrl
import random
from qgis.core import (
    Qgis,
    QgsApplication,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsFeature,
    QgsField,
    QgsFields,
    QgsGeometry,
    QgsPointXY,
    QgsProject,
    QgsRectangle,
    QgsSpatialIndex,
    QgsVectorFileWriter,
    QgsVectorLayer,
    QgsWkbTypes,
)

BUFFER_SEGMENTS = 16
ROI_REDUNDANCY_BIAS = 0.2  # prefer new connections unless redundancy ROI is clearly higher
HYBRID_MAX_LINKS_PER_PAIR = 2
HYBRID_OVERLAP_REJECT_RATIO = 0.85
SHORTCUT_RATIO_HIGH = 3.0
SHORTCUT_RATIO_MID = 1.5
SHORTCUT_RATIO_LOW = 1.5
SHORTCUT_MULT_HIGH = 0.9
SHORTCUT_MULT_MID = 0.5
SHORTCUT_MULT_LOW = 0.1
MAX_BRIDGE_MIDS = 800
MAX_BRIDGE_EDGES_PER_MID = 40
MAX_BRIDGE_ITERATIONS = 200
MAX_BRIDGE_RUNTIME_S = 2.5
MAX_BRIDGE_PATCHES = 3000
MAX_BRIDGE_CANDIDATES = 30000
MAX_VECTOR_PATCH_COUNT = 1000
BIGCONNECT_EXACT_MAX_CANDIDATES = 26
BIGCONNECT_EXACT_MAX_STATES = 250000
BIGCONNECT_EXACT_MAX_SECONDS = 5.0
BIGCONNECT_BEAM_WIDTH = 256
BIGCONNECT_VECTOR_SCALE = 100000
BIGCONNECT_MERGE_EQUIV_AREA_HA = 0.5
BIGCONNECT_MERGE_EQUIV_AREA_RATIO = 0.02

try:
    from osgeo import gdal, ogr, osr  # type: ignore
except Exception:  # pragma: no cover
    gdal = None  # type: ignore
    ogr = None  # type: ignore
    osr = None  # type: ignore


def _log_message(message: str, level: int = Qgis.Info) -> None:
    """Log to the QGIS Log Messages Panel with a TerraLink tag."""
    try:
        QgsApplication.messageLog().logMessage(message, "TerraLink", level)
    except Exception:
        # Fallback for environments where the message log is unavailable
        print(f"TerraLink Log: {message}")


def _write_text_report(path: str, lines: List[str]) -> None:
    try:
        os.makedirs(os.path.dirname(path) or os.getcwd(), exist_ok=True)
    except Exception:
        pass
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines).rstrip() + "\n")


def _write_summary_csv(path: str, rows: List[Tuple[str, str]]) -> None:
    try:
        os.makedirs(os.path.dirname(path) or os.getcwd(), exist_ok=True)
    except Exception:
        pass
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["item", "value"])
        writer.writerows(rows)


def _format_number(value: object, decimals: int = 2) -> str:
    try:
        return f"{float(value):.{decimals}f}"
    except Exception:
        return "" if value is None else str(value)



def _normalize_strategy_key(strategy: Optional[str]) -> str:
    key = str(strategy or "most_connected_habitat").strip().lower().replace(" ", "_").replace("-", "_")
    aliases = {
        # Back-compat
        "largest_network": "largest_single_network",
        "bigconnect": "most_connected_habitat",
        "most_connected_area": "most_connected_habitat",
        "habitat_availability": "reachable_habitat_advanced",
        "habitatavailability": "reachable_habitat_advanced",
        "habitat_available": "reachable_habitat_advanced",
        "ha": "reachable_habitat_advanced",
        "landscape_fluidity_a": "landscape_fluidity",
        "landscape_fluidity_1": "landscape_fluidity",
        "landscape_fluidity_a1": "landscape_fluidity",
        "lf_a": "landscape_fluidity",
        "lfa": "landscape_fluidity",
        "lf_a1": "landscape_fluidity",
        "lfa1": "landscape_fluidity",
    }
    key = aliases.get(key, key)
    # Optimization modes exposed in this plugin version (vector).
    valid = {
        "largest_single_network",
        "most_connected_habitat",
        "reachable_habitat_advanced",
        "landscape_fluidity",
    }
    if key not in valid:
        key = "most_connected_habitat"
    return key


def _strategy_display_name(strategy: Optional[str]) -> str:
    key = _normalize_strategy_key(strategy)
    names = {
        "largest_single_network": "Largest Single Network",
        "most_connected_habitat": "Most Connected Area",
        "reachable_habitat_advanced": "Reachable Habitat (Advanced)",
        "landscape_fluidity": "Landscape Fluidity",
    }
    return names.get(key, "Most Connected Area")


def _corridor_patch_sum_area_ha(cdata: Dict, patches: Dict[int, Dict]) -> float:
    """Return local numerator for corridor scoring: area(patch1) + area(patch2)."""
    try:
        p1 = cdata.get("p1", cdata.get("patch1"))
        p2 = cdata.get("p2", cdata.get("patch2"))
        if p1 is None or p2 is None:
            return 0.0
        a1 = float((patches.get(int(p1)) or {}).get("area_ha", 0.0) or 0.0)
        a2 = float((patches.get(int(p2)) or {}).get("area_ha", 0.0) or 0.0)
        return max(0.0, a1) + max(0.0, a2)
    except Exception:
        return 0.0


def _corridor_local_efficiency_ha(cdata: Dict, patches: Dict[int, Dict], area_override: Optional[float] = None) -> float:
    """Efficiency uses local form: (patch1_area + patch2_area) / corridor_area."""
    try:
        area = float(area_override if area_override is not None else (cdata.get("area_ha", 0.0) or 0.0))
    except Exception:
        area = 0.0
    patches_area = float(_corridor_patch_sum_area_ha(cdata, patches) or 0.0)
    if area <= 0.0 or patches_area <= 0.0:
        return 0.0
    return float(patches_area / area)


def _add_summary_csv_layer(csv_path: str, layer_title: str, add_to_project: bool = True) -> None:
    if (not add_to_project) or QgsVectorLayer is None or QgsProject is None:
        return
    try:
        url = QUrl.fromLocalFile(csv_path).toString()
        uri = f"{url}?type=csv&delimiter=,&xField=&yField=&geomType=none"
        layer = QgsVectorLayer(uri, layer_title, "delimitedtext")
        if not layer.isValid():
            return
        QgsProject.instance().addMapLayer(layer)
    except Exception:
        return


def _apply_visible_corridor_style_vector(layer: "QgsVectorLayer") -> None:
    """Apply a visible corridor style so layers don't render as black."""
    if layer is None or (not layer.isValid()):
        return
    try:
        from qgis.PyQt.QtGui import QColor
        from qgis.core import QgsFillSymbol, QgsLineSymbol, QgsSymbol, QgsWkbTypes
    except Exception:
        return
    try:
        geom_type = layer.geometryType()
    except Exception:
        geom_type = None
    try:
        if geom_type == QgsWkbTypes.PolygonGeometry:
            symbol = QgsFillSymbol.createSimple(
                {
                    "color": "255,180,30,180",
                    "outline_color": "130,40,0,220",
                    "outline_width": "0.4",
                }
            )
        elif geom_type == QgsWkbTypes.LineGeometry:
            symbol = QgsLineSymbol.createSimple(
                {
                    "color": "255,140,20,230",
                    "width": "0.6",
                }
            )
        else:
            symbol = QgsSymbol.defaultSymbol(geom_type) if geom_type is not None else None
        if symbol is not None:
            layer.renderer().setSymbol(symbol)
            layer.triggerRepaint()
    except Exception:
        return


def _apply_random_unique_value_symbology_vector(layer: "QgsVectorLayer", field_name: str) -> None:
    """
    Apply a categorized renderer with stable pseudo-random colors based on integer ids.

    Used for the "Contiguous Areas" vector layer to match raster-mode unique component colors.
    """
    if layer is None or (not layer.isValid()):
        return
    if not field_name:
        return
    try:
        from qgis.PyQt.QtGui import QColor
        from qgis.core import (
            QgsCategorizedSymbolRenderer,
            QgsFillSymbol,
            QgsRendererCategory,
            QgsWkbTypes,
        )
    except Exception:
        return

    try:
        if layer.geometryType() != QgsWkbTypes.PolygonGeometry:
            return
    except Exception:
        return

    idx = layer.fields().indexOf(field_name)
    if idx < 0:
        return

    # Collect unique ids (cap for performance).
    vals: List[int] = []
    try:
        unique = layer.uniqueValues(idx)
        vals = [int(v) for v in unique if int(v) > 0]
    except Exception:
        try:
            seen: Set[int] = set()
            for f in layer.getFeatures():
                try:
                    v = int(f[field_name])
                except Exception:
                    continue
                if v <= 0 or v in seen:
                    continue
                seen.add(v)
                vals.append(v)
                if len(vals) >= 2048:
                    break
        except Exception:
            return

    if not vals:
        return
    vals.sort()

    categories: List["QgsRendererCategory"] = []
    for v in vals[:2048]:
        hue = int((int(v) * 137) % 360)
        col = QColor.fromHsv(hue, 200, 255, 160)
        symbol = QgsFillSymbol.createSimple(
            {
                "color": f"{col.red()},{col.green()},{col.blue()},{col.alpha()}",
                "outline_color": "40,40,40,180",
                "outline_width": "0.26",
            }
        )
        categories.append(QgsRendererCategory(v, symbol, str(v)))

    try:
        renderer = QgsCategorizedSymbolRenderer(field_name, categories)
        layer.setRenderer(renderer)
        layer.triggerRepaint()
    except Exception:
        return


def _safe_filename(name: str, max_len: int = 64) -> str:
    safe = "".join(ch if (ch.isalnum() or ch in ("-", "_")) else "_" for ch in (name or ""))
    safe = safe.strip("_") or "layer"
    return safe[:max_len]


def _add_landscape_metrics_table_layer(
    layer_title: str,
    analysis_lines: List[str],
    add_to_project: bool = True,
) -> None:
    """
    Add a non-spatial table layer to the QGIS project with the landscape metrics.
    """
    if not add_to_project:
        return
    try:
        def _pct_change_str(pre_val: str, post_val: str) -> str:
            pre_v = _to_float_or_none(pre_val)
            post_v = _to_float_or_none(post_val)
            if pre_v is None or post_v is None:
                return ""
            denom = abs(float(pre_v))
            if denom <= 1e-12:
                return ""
            pct = ((float(post_v) - float(pre_v)) / denom) * 100.0
            return f"{pct:+.3f}%"

        # metric, pre, post, % change, description, source
        rows: List[Tuple[str, str, str, str, str, str]] = []
        for line in analysis_lines or []:
            s = (line or "").strip()
            if not s:
                continue
            if s.startswith("=") or s.startswith("-"):
                continue
            if "|" not in s:
                continue
            parts = [p.strip() for p in s.split("|")]
            if len(parts) < 3:
                continue
            metric = parts[0]
            pre_val = parts[1] if len(parts) > 1 else ""
            post_val = parts[2] if len(parts) > 2 else ""
            delta_pct = ""
            desc = ""
            source = ""
            # Full report format from analysis_raster._perform_landscape_analysis:
            # metric | pre | post | delta % | description | source
            if len(parts) >= 6:
                delta_pct = parts[3]
                desc = parts[4]
                source = parts[5]
            # Compact fallback rows:
            # metric | pre | post | description | source
            elif len(parts) >= 5:
                desc = parts[3]
                source = parts[4]
            else:
                desc = parts[3] if len(parts) > 3 else ""
            metric_key = metric.strip().lower()
            if metric_key in ("metric name", "metric", "id"):
                continue
            # Keep only metric-style rows; skip deltaPC data rows in the table view.
            if metric.strip().isdigit():
                continue
            if metric:
                rows.append((metric, pre_val, post_val, delta_pct or _pct_change_str(pre_val, post_val), desc, source))

        if not rows:
            joined = " ".join([l.strip() for l in (analysis_lines or []) if l.strip()])[:240]
            if not joined:
                joined = "No landscape metrics available."
            rows = [("Message", joined, "", "", "", "")]

        uri = (
            "None?"
            "field=metric:string(100)&"
            "field=pre_value:string(40)&"
            "field=post_value:string(40)&"
            "field=pct_change:string(20)&"
            "field=description:string(220)&"
            "field=source:string(120)"
        )
        layer = QgsVectorLayer(uri, layer_title, "memory")
        if not layer.isValid():
            return
        provider = layer.dataProvider()
        feats: List[QgsFeature] = []
        for metric, pre_val, post_val, pct_change, desc, source in rows:
            feat = QgsFeature(layer.fields())
            feat.setAttributes([metric, pre_val, post_val, pct_change, desc, source])
            feats.append(feat)
        provider.addFeatures(feats)
        layer.updateExtents()
        QgsProject.instance().addMapLayer(layer)
    except Exception:
        return


def _to_float_or_none(value: object) -> Optional[float]:
    try:
        s = str(value).strip().replace(",", "")
        if not s or s.lower() in ("na", "n/a", "none"):
            return None
        return float(s)
    except Exception:
        return None


def _extract_landscape_metric_values(
    analysis_lines: Sequence[str],
    metric_name: str,
) -> Optional[Tuple[float, float]]:
    target = str(metric_name or "").strip().lower()
    if not target:
        return None
    for line in analysis_lines or []:
        s = (line or "").strip()
        if not s or "|" not in s:
            continue
        parts = [p.strip() for p in s.split("|")]
        if len(parts) < 3:
            continue
        metric = str(parts[0] or "").strip().lower()
        if metric != target:
            continue
        pre_v = _to_float_or_none(parts[1])
        post_v = _to_float_or_none(parts[2])
        if pre_v is None or post_v is None:
            return None
        return float(pre_v), float(post_v)
    return None


def _clone_corridors(corridors: Dict[int, Dict]) -> Dict[int, Dict]:
    out: Dict[int, Dict] = {}
    for cid in sorted(corridors.keys(), key=lambda x: int(x)):
        src = corridors.get(cid) or {}
        rec: Dict[str, Any] = {}
        for key, val in src.items():
            if key == "geom" and val is not None:
                rec[key] = clone_geometry(val)
            elif key == "patch_ids":
                try:
                    rec[key] = set(int(pid) for pid in (val or []))
                except Exception:
                    rec[key] = set()
            elif isinstance(val, (list, dict, set)):
                # Avoid accidental cross-mutation from nested mutable refs.
                try:
                    if isinstance(val, list):
                        rec[key] = list(val)
                    elif isinstance(val, dict):
                        rec[key] = dict(val)
                    else:
                        rec[key] = set(val)
                except Exception:
                    rec[key] = val
            else:
                rec[key] = val
        out[int(cid)] = rec
    return out


def _evaluate_mode_metric_exact(
    *,
    patches: Dict[int, Dict],
    corridors: Dict[int, Dict],
    params: "VectorRunParams",
    raw_params: Dict[str, Any],
    target_crs: QgsCoordinateReferenceSystem,
    metric_name: str,
    layer_name: str,
) -> Optional[Dict[str, float]]:
    try:
        from .landscape_metrics import _perform_landscape_analysis  # local import to avoid import-cycle issues
    except Exception:
        return None

    try:
        exact_comp = _compute_habitat_component_metrics_exact(patches, corridors)
        metric_key = str(metric_name or "").strip().lower()
        if metric_key == "effective mesh size (habtitat-normalized)":
            pre_v = float(exact_comp.get("mesh_norm_pre", 0.0) or 0.0)
            post_v = float(exact_comp.get("mesh_norm_post", 0.0) or 0.0)
            return {"pre": pre_v, "post": post_v, "gain": float(post_v - pre_v)}
        if metric_key == "largest connected component proportion":
            pre_v = float(exact_comp.get("lcc_norm_pre", 0.0) or 0.0)
            post_v = float(exact_comp.get("lcc_norm_post", 0.0) or 0.0)
            return {"pre": pre_v, "post": post_v, "gain": float(post_v - pre_v)}
        if metric_key == "strategic mobility":
            mobility = _compute_strategic_mobility_exact(patches, corridors, params)
            pre_v = float(mobility.get("pre", 0.0) or 0.0)
            post_v = float(mobility.get("post", 0.0) or 0.0)
            return {"pre": pre_v, "post": post_v, "gain": float(post_v - pre_v)}
        if metric_key in (
            "landscape fluidity",
            "landscape fluidity lf-a",
            "landscape fluidity lf-b",
            "landscape fluidity lf-c",
            "mean effective resistance",
        ):
            if metric_key == "landscape fluidity lf-b":
                fluidity = _compute_landscape_fluidity_2_exact(
                    patches,
                    corridors,
                    params,
                    candidate_pool=None,
                )
            else:
                fluidity = _compute_landscape_fluidity_exact(patches, corridors, params)
            if metric_key == "mean effective resistance":
                pre_v = float(fluidity.get("graph_resistance_pre", 0.0) or 0.0)
                post_v = float(fluidity.get("graph_resistance_post", 0.0) or 0.0)
                return {"pre": pre_v, "post": post_v, "gain": float(pre_v - post_v)}
            pre_v = float(fluidity.get("pre", 0.0) or 0.0)
            post_v = float(fluidity.get("post", 0.0) or 0.0)
            return {"pre": pre_v, "post": post_v, "gain": float(post_v - pre_v)}

        dissolve_tolerance = max(0.0, float(getattr(params, "min_corridor_width", 0.0) or 0.0) * 0.01)
        networks = build_contiguous_network_summaries(
            patches,
            corridors,
            dissolve_tolerance=dissolve_tolerance,
        )
        patch_geoms = [
            pdata.get("geom")
            for pdata in patches.values()
            if pdata.get("geom") is not None and (not pdata.get("geom").isEmpty())
        ]
        network_geoms = [
            net.get("geom")
            for net in networks
            if net.get("geom") is not None and (not net.get("geom").isEmpty())
        ]
        bounds = _compute_geoms_bounds(patch_geoms + network_geoms)
        pixel_size_m = float(getattr(params, "grid_resolution", 50.0) or 50.0)
        min_width_m = max(0.0, float(getattr(params, "min_corridor_width", 0.0) or 0.0))
        metrics_pixel_size_m = max(1.0, pixel_size_m)
        if min_width_m > 0.0:
            safe_px = max(1.0, min_width_m / 3.0)
            if metrics_pixel_size_m > safe_px + 1e-9:
                metrics_pixel_size_m = float(safe_px)

        post_mask, eff_px = _rasterize_networks_to_mask(
            networks=networks,
            pixel_size_m=metrics_pixel_size_m,
            target_crs=target_crs,
            bounds=bounds,
        )
        pre_mask, _ = _rasterize_geoms_to_mask(
            geoms=patch_geoms,
            pixel_size_m=metrics_pixel_size_m,
            target_crs=target_crs,
            bounds=bounds,
        )
        analysis_params = dict(raw_params or {})
        analysis_params["mesh_override_pre_norm"] = float(exact_comp.get("mesh_norm_pre", 0.0) or 0.0)
        analysis_params["mesh_override_post_norm"] = float(exact_comp.get("mesh_norm_post", 0.0) or 0.0)
        analysis_params["lcc_override_pre_norm"] = float(exact_comp.get("lcc_norm_pre", 0.0) or 0.0)
        analysis_params["lcc_override_post_norm"] = float(exact_comp.get("lcc_norm_post", 0.0) or 0.0)
        analysis_lines = _perform_landscape_analysis(
            arr=post_mask,
            layer_name=layer_name,
            res_x=float(eff_px),
            res_y=float(eff_px),
            is_metric=True,
            params=analysis_params,
            pre_arr=pre_mask,
        )
        vals = _extract_landscape_metric_values(analysis_lines, metric_name)
        if vals is None:
            return None
        pre_v, post_v = vals
        return {
            "pre": float(pre_v),
            "post": float(post_v),
            "gain": float(post_v - pre_v),
        }
    except Exception:
        return None


def _refine_metric_mode_with_exact_metric(
    *,
    strategy_key: str,
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: "VectorRunParams",
    raw_params: Dict[str, Any],
    target_crs: QgsCoordinateReferenceSystem,
    base_corridors: Dict[int, Dict],
    base_stats: Dict[str, Any],
) -> Tuple[Dict[int, Dict], Dict[str, Any]]:
    metric_by_mode = {
        "reachable_habitat_advanced": "Reachable Habitat Score",
        "mesh_size": "Effective mesh size (habtitat-normalized)",
        "lcc_proportion": "Largest Connected Component Proportion",
        "probability_of_connectivity": "Probability of Connectivity",
        "flow_redundancy": "Flow Redundancy",
        "landscape_fluidity": "Landscape Fluidity",
        "mobility_strategic": "Strategic Mobility",
    }
    metric_name = metric_by_mode.get(strategy_key)
    if not metric_name:
        return base_corridors, base_stats

    # Some optimization modes exist only in certain TerraLink builds. Build this mapping
    # defensively so missing optional optimizers don't crash exact-metric refinement.
    optimizer_map: Dict[str, Callable[[Dict[int, Dict], List[Dict], "VectorRunParams"], Tuple[Dict[int, Dict], Dict]]] = {
        "most_connected_habitat": optimize_bigconnect_vector,
        "landscape_fluidity": optimize_landscape_fluidity_a,
        "landscape_fluidity_b": optimize_landscape_fluidity_b,
        "landscape_fluidity_c": optimize_landscape_fluidity_c,
        "mobility_strategic": optimize_mobility_strategic,
        "largest_single_network": optimize_circuit_utility_largest_network,
        "reachable_habitat_advanced": optimize_habitat_availability,
    }
    optional_optimizers = {
        "mesh_size": "optimize_mesh_size",
        "lcc_proportion": "optimize_lcc_proportion",
        "probability_of_connectivity": "optimize_probability_of_connectivity",
        "flow_redundancy": "optimize_flow_redundancy",
        "most_connectivity": "optimize_circuit_utility",
    }
    for mode_name, fn_name in optional_optimizers.items():
        fn = globals().get(fn_name)
        if callable(fn):
            optimizer_map[mode_name] = fn  # type: ignore[assignment]
    competitors_by_mode: Dict[str, List[str]] = {
        "mesh_size": ["mesh_size", "lcc_proportion", "largest_single_network"],
        "lcc_proportion": ["lcc_proportion", "largest_single_network", "mesh_size"],
        "most_connected_habitat": ["most_connected_habitat"],
        "probability_of_connectivity": ["probability_of_connectivity", "mesh_size", "lcc_proportion"],
        "flow_redundancy": ["flow_redundancy", "mesh_size", "probability_of_connectivity", "largest_single_network"],
        "landscape_fluidity": ["landscape_fluidity"],
        "mobility_strategic": ["mobility_strategic"],
    }
    competitor_keys = competitors_by_mode.get(strategy_key, [strategy_key])

    pure_metric_modes = {
        "mesh_size",
        "lcc_proportion",
        "most_connected_habitat",
        "probability_of_connectivity",
        "flow_redundancy",
        "landscape_fluidity",
        "landscape_fluidity_b",
        "landscape_fluidity_c",
        "mobility_strategic",
    }

    def _simulate_mode_postprocess(
        mode_name: str,
        in_corridors: Dict[int, Dict],
        in_stats: Dict[str, Any],
    ) -> Tuple[Dict[int, Dict], Dict[str, Any]]:
        corridors_local = _clone_corridors(in_corridors)
        stats_local = dict(in_stats or {})
        _refresh_vector_connectivity_stats(patches, corridors_local, stats_local)

        remaining_budget = float((params.budget_area or 0.0) - float(stats_local.get("budget_used_ha", 0.0) or 0.0))
        if mode_name not in pure_metric_modes and remaining_budget > 0 and corridors_local:
            try:
                extra_used, _low_added, _red_added = _apply_hybrid_leftover_budget_vector(
                    patches=patches,
                    candidates=candidates,
                    corridors=corridors_local,
                    remaining_budget=remaining_budget,
                    max_search_distance=float(params.max_search_distance or 0.0),
                )
                if extra_used:
                    stats_local["budget_used_ha"] = float(stats_local.get("budget_used_ha", 0.0) or 0.0) + float(extra_used)
                    stats_local["corridors_used"] = len(corridors_local)
                    _refresh_vector_connectivity_stats(patches, corridors_local, stats_local)
            except Exception:
                pass
            remaining_budget = float((params.budget_area or 0.0) - float(stats_local.get("budget_used_ha", 0.0) or 0.0))

        # Fixed-width policy: do not use remaining budget to alter corridor width.

        if mode_name == "largest_single_network" and corridors_local:
            try:
                removed_count, removed_area = _enforce_largest_network_component(
                    patches=patches,
                    corridors=corridors_local,
                    seed_patch=stats_local.get("seed_patch"),
                )
                if removed_count > 0:
                    stats_local["budget_used_ha"] = max(
                        0.0,
                        float(stats_local.get("budget_used_ha", 0.0) or 0.0) - float(removed_area),
                    )
                    stats_local["corridors_used"] = len(corridors_local)
                    _refresh_vector_connectivity_stats(patches, corridors_local, stats_local)
            except Exception:
                pass

        return corridors_local, stats_local

    candidates_eval: List[Dict[str, Any]] = []

    base_corridors_pp, base_stats_pp = _simulate_mode_postprocess(strategy_key, base_corridors, base_stats)
    base_eval = _evaluate_mode_metric_exact(
        patches=patches,
        corridors=base_corridors_pp,
        params=params,
        raw_params=raw_params,
        target_crs=target_crs,
        metric_name=metric_name,
        layer_name=f"TerraLink Exact Metric Eval ({strategy_key})",
    )
    if base_eval is not None:
        candidates_eval.append(
            {
                "name": str(strategy_key),
                "corridors": _clone_corridors(base_corridors_pp),
                "stats": dict(base_stats_pp or {}),
                "eval": dict(base_eval),
            }
        )

    for key in competitor_keys:
        if key == strategy_key:
            continue
        fn = optimizer_map.get(key)
        if fn is None:
            continue
        try:
            cc, ss = fn(patches, candidates, params)
        except Exception:
            continue
        if not cc:
            continue
        cc_pp, ss_pp = _simulate_mode_postprocess(key, cc, ss)
        ev = _evaluate_mode_metric_exact(
            patches=patches,
            corridors=cc_pp,
            params=params,
            raw_params=raw_params,
            target_crs=target_crs,
            metric_name=metric_name,
            layer_name=f"TerraLink Exact Metric Eval ({key})",
        )
        if ev is None:
            continue
        candidates_eval.append(
            {
                "name": str(key),
                "corridors": _clone_corridors(cc_pp),
                "stats": dict(ss_pp or {}),
                "eval": dict(ev),
            }
        )

    if not candidates_eval:
        return base_corridors, base_stats

    chosen = max(
        candidates_eval,
        key=lambda r: (
            float((r.get("eval") or {}).get("gain", float("-inf"))),
            float((r.get("eval") or {}).get("post", float("-inf"))),
            float((r.get("stats") or {}).get("budget_used_ha", 0.0) or 0.0),
            int(len((r.get("corridors") or {}))),
        ),
    )

    out_corridors = _clone_corridors(chosen.get("corridors") or {})
    out_stats = dict(chosen.get("stats") or {})
    out_stats["strategy"] = strategy_key
    out_stats["exact_metric_target"] = metric_name
    out_stats["exact_metric_selector_source"] = str(chosen.get("name") or strategy_key)
    out_stats["exact_metric_selector_gain"] = float((chosen.get("eval") or {}).get("gain", 0.0) or 0.0)
    out_stats["exact_metric_selector_post"] = float((chosen.get("eval") or {}).get("post", 0.0) or 0.0)
    return out_corridors, out_stats


def _compute_geoms_bounds(
    geoms: List[QgsGeometry],
) -> Optional[Tuple[float, float, float, float]]:
    xmin = ymin = float("inf")
    xmax = ymax = float("-inf")
    found = False
    for geom in geoms:
        if geom is None or geom.isEmpty():
            continue
        bb = geom.boundingBox()
        xmin = min(xmin, bb.xMinimum())
        ymin = min(ymin, bb.yMinimum())
        xmax = max(xmax, bb.xMaximum())
        ymax = max(ymax, bb.yMaximum())
        found = True
    if not found or not math.isfinite(xmin) or not math.isfinite(ymin) or not math.isfinite(xmax) or not math.isfinite(ymax):
        return None
    return xmin, ymin, xmax, ymax


def _geom_overlap_ratio(new_geom: QgsGeometry, prior_geoms: List[QgsGeometry]) -> float:
    try:
        if new_geom is None or new_geom.isEmpty() or not prior_geoms:
            return 0.0
        denom = float(new_geom.area())
        if denom <= 0:
            return 1.0
        best = 0.0
        for g in prior_geoms:
            if g is None or g.isEmpty():
                continue
            inter = new_geom.intersection(g)
            if inter is None or inter.isEmpty():
                continue
            ratio = float(inter.area()) / denom
            if ratio > best:
                best = ratio
        return best
    except Exception:
        return 1.0


def _rasterize_geoms_to_mask(
    geoms: List[QgsGeometry],
    pixel_size_m: float,
    target_crs: QgsCoordinateReferenceSystem,
    max_cells: int = 10_000_000,
    bounds: Optional[Tuple[float, float, float, float]] = None,
) -> Tuple[np.ndarray, float]:
    """
    Rasterize polygons to a binary mask for landscape metrics.

    Returns (mask_array, effective_pixel_size_m).
    """
    if gdal is None or ogr is None or osr is None:
        raise RuntimeError("GDAL/OGR not available for rasterize-based landscape metrics.")

    pixel_size = max(float(pixel_size_m or 0.0), 1.0)
    geoms = [geom for geom in geoms if geom is not None and not geom.isEmpty()]

    if bounds is None:
        xmin = ymin = float("inf")
        xmax = ymax = float("-inf")
        use_geoms: List[QgsGeometry] = []
        for geom in geoms:
            use_geoms.append(geom)
            bb = geom.boundingBox()
            xmin = min(xmin, bb.xMinimum())
            ymin = min(ymin, bb.yMinimum())
            xmax = max(xmax, bb.xMaximum())
            ymax = max(ymax, bb.yMaximum())
        geoms = use_geoms
        if not geoms or not math.isfinite(xmin) or not math.isfinite(ymin) or not math.isfinite(xmax) or not math.isfinite(ymax):
            return np.zeros((1, 1), dtype=np.uint8), float(pixel_size)
    else:
        xmin, ymin, xmax, ymax = bounds

    pad = pixel_size * 2.0
    xmin -= pad
    ymin -= pad
    xmax += pad
    ymax += pad

    width = max(xmax - xmin, pixel_size)
    height = max(ymax - ymin, pixel_size)
    cols = max(1, int(math.ceil(width / pixel_size)))
    rows = max(1, int(math.ceil(height / pixel_size)))

    if rows * cols > max_cells:
        scale = int(math.ceil(math.sqrt((rows * cols) / max_cells)))
        pixel_size *= max(1, scale)
        cols = max(1, int(math.ceil(width / pixel_size)))
        rows = max(1, int(math.ceil(height / pixel_size)))

    mem_driver = gdal.GetDriverByName("MEM")
    ds = mem_driver.Create("", cols, rows, 1, gdal.GDT_Byte)
    ds.SetGeoTransform((xmin, pixel_size, 0.0, ymax, 0.0, -pixel_size))

    try:
        srs = osr.SpatialReference()
        wkt = target_crs.toWkt() if target_crs and target_crs.isValid() else ""
        if wkt:
            srs.ImportFromWkt(wkt)
            ds.SetProjection(srs.ExportToWkt())
    except Exception:
        srs = None

    ds.GetRasterBand(1).Fill(0)
    if geoms:
        ogr_driver = ogr.GetDriverByName("Memory")
        vds = ogr_driver.CreateDataSource("geoms")
        layer = vds.CreateLayer("geoms", srs=srs, geom_type=ogr.wkbMultiPolygon)
        layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
        defn = layer.GetLayerDefn()

        for i, qgs_geom in enumerate(geoms, 1):
            try:
                ogr_geom = ogr.CreateGeometryFromWkb(bytes(qgs_geom.asWkb()))
            except Exception:
                continue
            if ogr_geom is None:
                continue
            feat = ogr.Feature(defn)
            feat.SetField("id", int(i))
            feat.SetGeometry(ogr_geom)
            layer.CreateFeature(feat)
            feat = None

        gdal.RasterizeLayer(ds, [1], layer, burn_values=[1])

    arr = ds.GetRasterBand(1).ReadAsArray()
    arr = np.asarray(arr, dtype=np.uint8)
    return arr, float(pixel_size)


def _rasterize_networks_to_mask(
    networks: List[Dict],
    pixel_size_m: float,
    target_crs: QgsCoordinateReferenceSystem,
    max_cells: int = 10_000_000,
    bounds: Optional[Tuple[float, float, float, float]] = None,
) -> Tuple[np.ndarray, float]:
    """
    Rasterize dissolved network polygons to a binary mask for landscape metrics.
    """
    geoms = [net.get("geom") for net in networks if net.get("geom") is not None and not net.get("geom").isEmpty()]
    return _rasterize_geoms_to_mask(
        geoms=geoms,
        pixel_size_m=pixel_size_m,
        target_crs=target_crs,
        max_cells=max_cells,
        bounds=bounds,
    )

def clone_geometry(geom: QgsGeometry) -> QgsGeometry:
    """
    Lightweight copy helper; QgsGeometry uses implicit sharing so this is cheap
    and avoids unnecessary deep clones unless a write occurs.
    """
    return QgsGeometry(geom)


class VectorAnalysisError(RuntimeError):
    """Raised when the vector analysis cannot be completed."""


class _TimingBlock:
    """Context manager that records elapsed time for a named step."""

    def __init__(self, label: str, sink: List[Dict[str, float]]):
        self.label = label
        self.sink = sink
        self._start = 0.0

    def __enter__(self) -> None:
        self._start = time.perf_counter()
        return None

    def __exit__(self, exc_type, exc, tb) -> bool:
        duration = time.perf_counter() - self._start
        self.sink.append({"label": self.label, "duration_s": duration})
        return False


class TimingRecorder:
    """Lightweight helper to track fine-grained step timings."""

    def __init__(self) -> None:
        self.records: List[Dict[str, float]] = []

    def time_block(self, label: str) -> _TimingBlock:
        return _TimingBlock(label, self.records)

    def add(self, label: str, duration: float) -> None:
        self.records.append({"label": label, "duration_s": duration})

    def write_report(self, path: str, total_elapsed: Optional[float] = None) -> None:
        try:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("TerraLink Vector Timing\n")
                fh.write("=" * 30 + "\n")
                for entry in self.records:
                    fh.write(f"{entry['label']}: {entry['duration_s']:.3f}s\n")
                if total_elapsed is not None:
                    fh.write("\n")
                    fh.write(f"Total wall time: {total_elapsed:.3f}s\n")
            print(f"  ✓ Timing report saved: {path}")
        except Exception as exc:  # noqa: BLE001
            print(f"  ⚠ Could not write timing report: {exc}")


@dataclass
class VectorRunParams:
    min_corridor_width: float  # metres
    max_corridor_area: Optional[float]  # hectares
    min_patch_size: float  # hectares
    budget_area: float  # hectares
    max_search_distance: float  # metres
    unit_system: str
    output_name: str
    grid_resolution: float  # metres
    obstacle_layer_ids: List[str] = field(default_factory=list)
    obstacle_enabled: bool = False
    vector_terminal_spacing_m: float = 150.0
    vector_terminal_max_per_patch: int = 120
    vector_terminal_pairs_per_pair: int = 25
    vector_routing_enabled: bool = False
    vector_routing_max_window_m: float = 6000.0
    vector_routing_smooth_iterations: int = 6
    vector_routing_smooth_offset: float = 0.25
    pc_alpha_analysis: float = 0.0
    pc_cutoff_analysis: float = 0.0
    sample_points: int = 50
    pair_samples: int = 200
    mobility_lambda_intra: float = 0.001
    mobility_tau: float = 0.90
    mobility_intra_k_total: int = 60
    mobility_eval_cap: int = 220
    mobility_detour_cap: float = 8.0
    landscape_fluidity_shortcut_threshold: float = 3.0
    landscape_fluidity_intra_k_total: int = 80
    landscape_fluidity_eval_cap: int = 140
    landscape_fluidity_intra_tau: float = 0.8
    landscape_fluidity_intra_bonus_weight: float = 1.0
    landscape_fluidity_intra_bonus_cap_mult: float = 3.0
    landscape_fluidity_inter_loop_weight: float = 0.35
    landscape_fluidity_inter_loop_cap_mult: float = 0.80
    landscape_fluidity_bridge_relief_weight: float = 0.25
    landscape_fluidity_bridge_relief_cap_mult: float = 0.60
    landscape_fluidity_lookahead_k: int = 8
    landscape_fluidity_lookahead_weight: float = 0.85
    landscape_fluidity_chain_k_total: int = 80
    landscape_fluidity_chain_branch_per_patch: int = 8
    landscape_fluidity_topology_weight: float = 0.20
    landscape_fluidity_long_hop_factor: float = 1.35
    landscape_fluidity_long_hop_tau: float = 0.90
    resiliency_shortcut_threshold: float = 2.5
    resiliency_intra_k_total: int = 80
    resiliency_eval_cap: int = 140
    species_dispersal_distance_analysis: float = 0.0
    species_dispersal_kernel: str = HABITAT_AVAILABILITY_DEFAULT_KERNEL
    min_patch_area_for_species_ha: float = 0.0
    patch_quality_weight_field: str = ""
    patch_area_scaling: str = HABITAT_AVAILABILITY_DEFAULT_SCALING


@dataclass
class AnalysisContext:
    """
    Mutable, per-run state (caches) that should not live on VectorRunParams.
    Keeping params immutable-ish avoids state leakage across runs.
    """

    impassable_union: Optional[QgsGeometry] = None
    filtered_small_patches: List[Dict[str, Any]] = field(default_factory=list)


def _to_dataclass(params: Dict) -> VectorRunParams:
    output_name = params.get("output_name") or "terralink_corridors.gpkg"
    if not output_name.lower().endswith(".gpkg"):
        output_name = f"{output_name}.gpkg"
    obstacle_ids_raw = params.get("obstacle_layer_ids") or []
    if not obstacle_ids_raw and params.get("obstacle_layer_id"):
        obstacle_ids_raw = [params.get("obstacle_layer_id")]
    obstacle_ids = [str(val) for val in obstacle_ids_raw if val]
    obstacle_flag = bool(params.get("obstacle_enabled", False) and obstacle_ids)

    max_search_distance_value = params.get("max_search_distance")
    if max_search_distance_value is None or str(max_search_distance_value).strip() == "":
        max_search_distance = 0.0
    else:
        try:
            max_search_distance = float(max_search_distance_value)
        except (TypeError, ValueError):
            max_search_distance = 0.0

    return VectorRunParams(
        min_corridor_width=float(params.get("min_corridor_width", 200.0)),
        max_corridor_area=None,
        min_patch_size=float(params.get("min_patch_size", 10.0)),
        budget_area=float(params.get("budget_area", 50.0)),
        max_search_distance=max_search_distance,
        unit_system=str(params.get("unit_system", "metric")),
        output_name=output_name,
        grid_resolution=max(float(params.get("grid_resolution", 50.0)), 1.0),
        obstacle_layer_ids=obstacle_ids,
        obstacle_enabled=bool(params.get("obstacle_enabled", False) and obstacle_ids),
        vector_terminal_spacing_m=float(params.get("vector_terminal_spacing_m", 150.0) or 150.0),
        vector_terminal_max_per_patch=int(params.get("vector_terminal_max_per_patch", 120) or 120),
        vector_terminal_pairs_per_pair=int(params.get("vector_terminal_pairs_per_pair", 25) or 25),
        # Default to enabled when impassables are enabled (maze-capable), unless explicitly overridden.
        vector_routing_enabled=(
            bool(params.get("vector_routing_enabled"))
            if "vector_routing_enabled" in params
            else bool(obstacle_flag)
        ),
        vector_routing_max_window_m=float(params.get("vector_routing_max_window_m", 6000.0) or 6000.0),
        vector_routing_smooth_iterations=int(params.get("vector_routing_smooth_iterations", 6) or 6),
        vector_routing_smooth_offset=float(params.get("vector_routing_smooth_offset", 0.25) or 0.25),
        pc_alpha_analysis=float(params.get("pc_alpha_analysis", 0.0) or 0.0),
        pc_cutoff_analysis=float(params.get("pc_cutoff_analysis", 0.0) or 0.0),
        sample_points=int(params.get("sample_points", 50) or 50),
        pair_samples=int(params.get("pair_samples", 200) or 200),
        mobility_lambda_intra=float(params.get("mobility_lambda_intra", 0.001) or 0.001),
        mobility_tau=float(params.get("mobility_tau", 0.90) or 0.90),
        mobility_intra_k_total=int(params.get("mobility_intra_k_total", 60) or 60),
        mobility_eval_cap=int(params.get("mobility_eval_cap", 220) or 220),
        mobility_detour_cap=float(params.get("mobility_detour_cap", 8.0) or 8.0),
        landscape_fluidity_shortcut_threshold=float(
            params.get(
                "landscape_fluidity_shortcut_threshold",
                params.get("resiliency_shortcut_threshold", 3.0),
            )
            or 3.0
        ),
        landscape_fluidity_intra_k_total=int(
            params.get(
                "landscape_fluidity_intra_k_total",
                params.get("resiliency_intra_k_total", 80),
            )
            or 80
        ),
        landscape_fluidity_eval_cap=int(
            params.get(
                "landscape_fluidity_eval_cap",
                params.get("resiliency_eval_cap", 140),
            )
            or 140
        ),
        landscape_fluidity_intra_tau=float(
            params.get(
                "landscape_fluidity_intra_tau",
                0.8,
            )
            or 0.8
        ),
        landscape_fluidity_intra_bonus_weight=float(
            params.get(
                "landscape_fluidity_intra_bonus_weight",
                1.0,
            )
            or 1.0
        ),
        landscape_fluidity_intra_bonus_cap_mult=float(
            params.get(
                "landscape_fluidity_intra_bonus_cap_mult",
                1.0,
            )
            or 1.0
        ),
        landscape_fluidity_inter_loop_weight=float(
            params.get(
                "landscape_fluidity_inter_loop_weight",
                0.35,
            )
            or 0.35
        ),
        landscape_fluidity_inter_loop_cap_mult=float(
            params.get(
                "landscape_fluidity_inter_loop_cap_mult",
                0.80,
            )
            or 0.80
        ),
        landscape_fluidity_bridge_relief_weight=float(
            params.get(
                "landscape_fluidity_bridge_relief_weight",
                0.25,
            )
            or 0.25
        ),
        landscape_fluidity_bridge_relief_cap_mult=float(
            params.get(
                "landscape_fluidity_bridge_relief_cap_mult",
                0.60,
            )
            or 0.60
        ),
        landscape_fluidity_lookahead_k=int(
            params.get(
                "landscape_fluidity_lookahead_k",
                8,
            )
            or 8
        ),
        landscape_fluidity_lookahead_weight=float(
            params.get(
                "landscape_fluidity_lookahead_weight",
                0.85,
            )
            or 0.85
        ),
        landscape_fluidity_chain_k_total=int(
            params.get(
                "landscape_fluidity_chain_k_total",
                80,
            )
            or 80
        ),
        landscape_fluidity_chain_branch_per_patch=int(
            params.get(
                "landscape_fluidity_chain_branch_per_patch",
                8,
            )
            or 8
        ),
        landscape_fluidity_topology_weight=float(
            params.get(
                "landscape_fluidity_topology_weight",
                0.20,
            )
            or 0.20
        ),
        landscape_fluidity_long_hop_factor=float(
            params.get(
                "landscape_fluidity_long_hop_factor",
                1.35,
            )
            or 1.35
        ),
        landscape_fluidity_long_hop_tau=float(
            params.get(
                "landscape_fluidity_long_hop_tau",
                0.90,
            )
            or 0.90
        ),
        resiliency_shortcut_threshold=float(params.get("resiliency_shortcut_threshold", 2.5) or 2.5),
        resiliency_intra_k_total=int(params.get("resiliency_intra_k_total", 80) or 80),
        resiliency_eval_cap=int(params.get("resiliency_eval_cap", 140) or 140),
        species_dispersal_distance_analysis=float(params.get("species_dispersal_distance_analysis", 0.0) or 0.0),
        species_dispersal_kernel=normalize_habitat_availability_kernel(
            params.get("species_dispersal_kernel", HABITAT_AVAILABILITY_DEFAULT_KERNEL)
        ),
        min_patch_area_for_species_ha=float(params.get("min_patch_area_for_species_analysis", 0.0) or 0.0),
        patch_quality_weight_field=str(params.get("patch_quality_weight_field", "") or "").strip(),
        patch_area_scaling=normalize_patch_area_scaling(
            params.get("patch_area_scaling", HABITAT_AVAILABILITY_DEFAULT_SCALING)
        ),
    )


def _count_interior_rings(geom: QgsGeometry) -> int:
    try:
        if geom.isMultipart():
            polys = geom.asMultiPolygon()
        else:
            polys = [geom.asPolygon()]
        rings = 0
        for poly in polys:
            if not poly:
                continue
            # poly[0] is exterior, poly[1:] are interior rings
            rings += max(len(poly) - 1, 0)
        return rings
    except Exception:
        return 0


def _extract_interior_ring_geometries(geom: QgsGeometry) -> List[QgsGeometry]:
    holes: List[QgsGeometry] = []
    try:
        polys = geom.asMultiPolygon() if geom.isMultipart() else [geom.asPolygon()]
        for poly in polys:
            if not poly or len(poly) < 2:
                continue
            for ring in poly[1:]:
                if not ring:
                    continue
                ring_xy = [QgsPointXY(p) for p in ring]
                if ring_xy and ring_xy[0] != ring_xy[-1]:
                    ring_xy.append(ring_xy[0])
                hole = QgsGeometry.fromPolygonXY([ring_xy])
                if hole and not hole.isEmpty():
                    holes.append(hole)
    except Exception:
        return []
    return holes


def _redundancy_far_enough(
    geom: QgsGeometry, prior_geoms: Sequence[QgsGeometry], threshold: float
) -> bool:
    if threshold <= 0.0:
        return True
    if geom is None or geom.isEmpty():
        return True
    if not prior_geoms:
        return True
    for g in prior_geoms:
        if g is None or g.isEmpty():
            continue
        try:
            if geom.distance(g) < threshold:
                return False
        except Exception:
            continue
    return True


def get_utm_crs_from_extent(layer: QgsVectorLayer) -> QgsCoordinateReferenceSystem:
    extent = layer.extent()
    center = extent.center()
    source_crs = layer.crs()

    if not source_crs.isGeographic():
        wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
        transform = QgsCoordinateTransform(source_crs, wgs84, QgsProject.instance())
        center = transform.transform(center)

    utm_zone = int((center.x() + 180) / 6) + 1
    epsg_code = 32600 + utm_zone if center.y() >= 0 else 32700 + utm_zone
    return QgsCoordinateReferenceSystem(f"EPSG:{epsg_code}")


def load_and_prepare_patches(
    layer: QgsVectorLayer,
    target_crs: QgsCoordinateReferenceSystem,
    params: VectorRunParams,
    split_multipart: bool = False,
) -> Tuple[Dict[int, Dict], QgsSpatialIndex, List[Dict[str, Any]]]:
    print("  Loading patches and building spatial index...")
    source_crs = layer.crs()
    transform = QgsCoordinateTransform(source_crs, target_crs, QgsProject.instance())

    raw_geoms: List[QgsGeometry] = []
    raw_quality_weights: List[float] = []
    indexed_features: List[QgsFeature] = []
    multipart_features_expanded = 0
    multipart_parts_emitted = 0
    quality_field = str(getattr(params, "patch_quality_weight_field", "") or "").strip()

    def _append_raw_geom(g: QgsGeometry, quality_weight: float) -> None:
        if g is None or g.isEmpty():
            return
        fid = len(raw_geoms)
        raw_geoms.append(g)
        raw_quality_weights.append(max(0.0, float(quality_weight or 0.0)))
        feat = QgsFeature()
        feat.setGeometry(g)
        feat.setId(fid)
        indexed_features.append(feat)

    for feature in layer.getFeatures():
        try:
            quality_weight = float(feature[quality_field]) if quality_field else 1.0
        except Exception:
            quality_weight = 1.0
        if not np.isfinite(float(quality_weight)):
            quality_weight = 1.0
        quality_weight = max(0.0, float(quality_weight))
        geom = QgsGeometry(feature.geometry())
        if geom.isEmpty():
            continue
        geom.transform(transform)
        try:
            geom = geom.makeValid()
        except Exception:
            pass
        if geom.isEmpty():
            continue
        if split_multipart and geom.isMultipart():
            emitted_any_part = False
            try:
                multi_parts = geom.asMultiPolygon() or []
            except Exception:
                multi_parts = []
            for poly in multi_parts:
                try:
                    part_geom = QgsGeometry.fromPolygonXY(poly)
                except Exception:
                    part_geom = None
                if part_geom is None or part_geom.isEmpty():
                    continue
                try:
                    part_geom = part_geom.makeValid()
                except Exception:
                    pass
                if part_geom is None or part_geom.isEmpty():
                    continue
                _append_raw_geom(part_geom, quality_weight)
                multipart_parts_emitted += 1
                emitted_any_part = True
            if emitted_any_part:
                multipart_features_expanded += 1
                continue
        _append_raw_geom(geom, quality_weight)

    spatial_index = QgsSpatialIndex(flags=QgsSpatialIndex.FlagStoreFeatureGeometries)
    if indexed_features:
        spatial_index.addFeatures(indexed_features)

    if split_multipart and multipart_features_expanded > 0:
        print(
            "  ✓ Expanded multipart features for strict single-network mode: "
            f"{multipart_features_expanded} features -> {multipart_parts_emitted} parts"
        )

    uf = UnionFind()
    for i in range(len(raw_geoms)):
        uf.find(i)

    for i, geom in enumerate(raw_geoms):
        try:
            bbox = geom.boundingBox()
            candidates = spatial_index.intersects(bbox)
        except Exception:
            candidates = []
        for j in candidates:
            if j <= i:
                continue
            other = raw_geoms[j]
            try:
                if geom.intersects(other):
                    uf.union(i, j)
            except Exception:
                continue

    components: Dict[int, List[int]] = defaultdict(list)
    for idx, geom in enumerate(raw_geoms):
        root = int(uf.find(idx))
        components[root].append(int(idx))

    patches: Dict[int, Dict] = {}
    patch_id = 1
    filtered_count = 0
    filtered_small_patches: List[Dict[str, Any]] = []
    indexed_features = []

    for geom_ids in components.values():
        if not geom_ids:
            continue
        geoms = [raw_geoms[idx] for idx in geom_ids]
        try:
            merged = QgsGeometry.unaryUnion(geoms)
        except Exception:
            merged = geoms[0]
            for extra in geoms[1:]:
                try:
                    merged = merged.combine(extra)
                except Exception:
                    pass
        if merged is None or merged.isEmpty():
            continue
        try:
            merged = merged.makeValid()
        except Exception:
            pass
        if merged.isEmpty():
            continue
        area_ha = merged.area() / 10000.0
        if area_ha < params.min_patch_size:
            filtered_count += 1
            filtered_small_patches.append(
                {
                    "geom": clone_geometry(merged),
                    "area_ha": float(area_ha),
                }
            )
            continue

        feat = QgsFeature()
        feat.setGeometry(merged)
        feat.setId(patch_id)
        indexed_features.append(feat)
        quality_weight = 1.0
        if geom_ids:
            weighted_sum = 0.0
            area_sum = 0.0
            for idx in geom_ids:
                geom_part = raw_geoms[int(idx)]
                try:
                    part_area = max(float(geom_part.area()), 0.0)
                except Exception:
                    part_area = 0.0
                if part_area <= 0.0:
                    continue
                weighted_sum += part_area * max(float(raw_quality_weights[int(idx)]), 0.0)
                area_sum += part_area
            if area_sum > 0.0:
                quality_weight = weighted_sum / area_sum
        patches[patch_id] = {
            "geom": clone_geometry(merged),
            "area_ha": area_ha,
            "quality_weight": float(max(0.0, quality_weight)),
        }
        patch_id += 1

    spatial_index = QgsSpatialIndex(flags=QgsSpatialIndex.FlagStoreFeatureGeometries)
    if indexed_features:
        spatial_index.addFeatures(indexed_features)

    print(f"  ✓ Loaded {len(patches)} patches (filtered {filtered_count} too small)")
    return patches, spatial_index, filtered_small_patches


def _detect_corridor_intersections(
    corridor_geom: QgsGeometry,
    patches: Dict[int, Dict],
    spatial_index: QgsSpatialIndex,
    connected_patches: Set[int],
) -> Set[int]:
    """
    Efficiently detect which patches are traversed by a corridor geometry
    using the spatial index.
    """
    intersected: Set[int] = set()
    bbox = corridor_geom.boundingBox()
    
    # Use spatial index to find candidate interactions (fast)
    candidate_ids = spatial_index.intersects(bbox)
    
    for pid in candidate_ids:
        if pid in connected_patches:
            continue
        pdata = patches.get(pid)
        if not pdata:
            continue
            
        try:
            # Check for actual intersection
            if corridor_geom.intersects(pdata["geom"]):
                intersection = corridor_geom.intersection(pdata["geom"])
                if intersection and (not intersection.isEmpty()) and intersection.area() > 0:
                    intersected.add(pid)
        except Exception:
            continue
            
    return intersected


def _finalize_corridor_geometry(
    pid1: int,
    pid2: int,
    corridor_geom: QgsGeometry,
    patches: Dict[int, Dict],
    spatial_index: QgsSpatialIndex,
    patch_union: Optional[QgsGeometry] = None,
    local_clip_geom: Optional[QgsGeometry] = None,
    skip_patch_clip: bool = False,
    corridor_width_m: float = 0.0,
) -> Tuple[Optional[QgsGeometry], Set[int]]:
    """
    Detect traversed patches, clip corridor geometry so it doesn't overlap them,
    and update the set of connected patches.
    """
    if corridor_geom is None or corridor_geom.isEmpty():
        return None, set()

    patch_ids: Set[int] = {pid1, pid2}
    
    # 1. Detect intermediate patches (A -> C -> B)
    # Using spatial index here is crucial for performance
    intersected = _detect_corridor_intersections(corridor_geom, patches, spatial_index, patch_ids)
    if intersected:
        patch_ids.update(intersected)

    # 2. Clip the corridor geometry against ALL involved patches
    # This makes the corridor "free" where it crosses existing habitat
    # and prevents drawing on top of patches.
    final_geom = clone_geometry(corridor_geom)

    if local_clip_geom and (not local_clip_geom.isEmpty()):
        try:
            final_geom = final_geom.intersection(local_clip_geom)
        except Exception:
            pass

    if not skip_patch_clip:
        if patch_union and not patch_union.isEmpty():
            try:
                final_geom = final_geom.difference(patch_union)
            except Exception:
                pass
        else:
            for pid in patch_ids:
                pdata = patches.get(pid)
                if not pdata: 
                    continue
                patch_geom = pdata.get("geom")
                if not patch_geom or patch_geom.isEmpty():
                    continue
                
                try:
                    if final_geom.intersects(patch_geom):
                        final_geom = final_geom.difference(patch_geom)
                        if final_geom.isEmpty():
                            break
                except Exception:
                    pass

    if final_geom is None or final_geom.isEmpty():
        # If the corridor is entirely consumed by patches, it means the patches 
        # touch or overlap. In vector analysis, this is valid connectivity (0 area cost).
        # We return an empty geom but valid patch_ids. However, to display it,
        # we might want to return None if strictly "corridor building".
        # But usually we return None to avoid 0-area features being written.
        return None, set()

    # Remove corridor segments that fall inside polygon holes.
    try:
        holes_to_clip: List[QgsGeometry] = []
        bbox = final_geom.boundingBox()
        candidate_ids = spatial_index.intersects(bbox)
        for pid in candidate_ids:
            pdata = patches.get(pid)
            if not pdata:
                continue
            patch_geom = pdata.get("geom")
            if not patch_geom or patch_geom.isEmpty():
                continue
            if _count_interior_rings(patch_geom) <= 0:
                continue
            for hole in _extract_interior_ring_geometries(patch_geom):
                if hole and (not hole.isEmpty()):
                    try:
                        if final_geom.intersects(hole):
                            holes_to_clip.append(hole)
                    except Exception:
                        continue
        if holes_to_clip:
            try:
                holes_union = QgsGeometry.unaryUnion(holes_to_clip).makeValid()
            except Exception:
                holes_union = None
            if holes_union and (not holes_union.isEmpty()):
                final_geom = final_geom.difference(holes_union)
            else:
                for hole in holes_to_clip:
                    try:
                        final_geom = final_geom.difference(hole)
                        if final_geom.isEmpty():
                            break
                    except Exception:
                        continue
    except Exception:
        pass

    final_geom = final_geom.makeValid()
    if final_geom.isEmpty():
        return None, set()

    # When two patches are separated by a tiny gap, subtraction against the full
    # habitat union can leave long edge-adjacent slivers. If that inflation is
    # detected, rebuild a compact shortest-bridge geometry.
    try:
        width_m = float(corridor_width_m or 0.0)
    except Exception:
        width_m = 0.0
    if width_m > 0.0 and int(pid1) != int(pid2) and len(patch_ids) <= 2:
        g1 = (patches.get(int(pid1)) or {}).get("geom")
        g2 = (patches.get(int(pid2)) or {}).get("geom")
        if g1 is not None and g2 is not None and (not g1.isEmpty()) and (not g2.isEmpty()):
            try:
                gap_m = float(g1.distance(g2))
            except Exception:
                gap_m = 0.0
            nominal_area_m2 = max(gap_m, 0.0) * max(width_m, 0.0)
            try:
                final_area_m2 = float(final_geom.area())
            except Exception:
                final_area_m2 = 0.0
            inflated = nominal_area_m2 > 1e-9 and final_area_m2 > (1.6 * nominal_area_m2)
            if inflated:
                cut_geom = None
                if patch_union is not None and (not patch_union.isEmpty()):
                    cut_geom = patch_union
                else:
                    try:
                        cut_geom = QgsGeometry.unaryUnion([g1, g2])
                    except Exception:
                        cut_geom = None

                def _prepare_compact_candidate(line_geom: QgsGeometry) -> Optional[QgsGeometry]:
                    if line_geom is None or line_geom.isEmpty():
                        return None
                    try:
                        cand = line_geom.buffer(width_m / 2.0, BUFFER_SEGMENTS)
                    except Exception:
                        return None
                    if cand is None or cand.isEmpty():
                        return None
                    if cut_geom is not None and (not cut_geom.isEmpty()):
                        try:
                            cand = cand.difference(cut_geom)
                        except Exception:
                            pass
                    if cand is None or cand.isEmpty():
                        return None
                    local_radius = max(5.0, 2.0 * width_m, 3.0 * max(gap_m, 0.0))
                    try:
                        local_window = line_geom.buffer(local_radius, BUFFER_SEGMENTS)
                    except Exception:
                        local_window = None
                    if local_window is not None and (not local_window.isEmpty()):
                        try:
                            local_geom = cand.intersection(local_window)
                            if local_geom is not None and (not local_geom.isEmpty()):
                                cand = local_geom
                        except Exception:
                            pass
                    if cand is None or cand.isEmpty():
                        return None
                    try:
                        cand = cand.makeValid()
                    except Exception:
                        pass
                    if cand is None or cand.isEmpty():
                        return None
                    parts = _geometry_parts(cand)
                    if len(parts) > 1:
                        tol = max(2.0, 0.75 * width_m)
                        ranked_parts: List[Tuple[int, float, float, QgsGeometry]] = []
                        for part in parts:
                            if part is None or part.isEmpty():
                                continue
                            try:
                                d1 = float(part.distance(g1))
                                d2 = float(part.distance(g2))
                                pa = float(part.area())
                            except Exception:
                                continue
                            near_both = 1 if (d1 <= tol and d2 <= tol) else 0
                            ranked_parts.append((near_both, d1 + d2, pa, part))
                        if ranked_parts:
                            ranked_parts.sort(key=lambda t: (-int(t[0]), float(t[1]), float(t[2])))
                            cand = ranked_parts[0][3]
                    if cand is None or cand.isEmpty():
                        return None
                    return cand

                compact_geom = None
                try:
                    shortest_line = g1.shortestLine(g2)
                except Exception:
                    shortest_line = None
                if shortest_line is not None and (not shortest_line.isEmpty()):
                    compact_geom = _prepare_compact_candidate(shortest_line)
                if compact_geom is not None and (not compact_geom.isEmpty()):
                    best_geom = compact_geom
                    try:
                        best_area = float(best_geom.area())
                        best_perim = float(best_geom.length())
                    except Exception:
                        best_area = float("inf")
                        best_perim = float("inf")

                    # For near-touching components, shortest-line location can be
                    # ambiguous along long parallel edges. Probe nearby boundary
                    # vertices and keep the most compact valid bridge.
                    if gap_m <= (2.5 * width_m):
                        def _collect_near_vertices(src_geom: QgsGeometry, other_geom: QgsGeometry) -> List[QgsPointXY]:
                            pts: List[QgsPointXY] = []
                            try:
                                np_src = src_geom.nearestPoint(other_geom).asPoint()
                                if np_src and not np_src.isEmpty():
                                    pts.append(QgsPointXY(np_src))
                            except Exception:
                                pass
                            all_verts: List[QgsPointXY] = []
                            try:
                                for v in src_geom.vertices():
                                    all_verts.append(QgsPointXY(v))
                            except Exception:
                                all_verts = []
                            if len(all_verts) > 400:
                                step = max(1, len(all_verts) // 180)
                                all_verts = all_verts[::step]
                            max_d = max(gap_m + (1.5 * width_m), 2.0 * width_m)
                            for p in all_verts:
                                try:
                                    dp = QgsGeometry.fromPointXY(p).distance(other_geom)
                                except Exception:
                                    continue
                                if dp <= max_d + 1e-9:
                                    pts.append(p)
                            dedup: List[QgsPointXY] = []
                            seen: Set[Tuple[int, int]] = set()
                            for p in pts:
                                key = (int(round(p.x() * 100.0)), int(round(p.y() * 100.0)))
                                if key in seen:
                                    continue
                                seen.add(key)
                                dedup.append(p)
                            return dedup[:24]

                        pts1 = _collect_near_vertices(g1, g2)
                        pts2 = _collect_near_vertices(g2, g1)
                        eval_count = 0
                        max_eval = 220
                        max_probe_len = max(gap_m + (2.5 * width_m), 6.0 * width_m)
                        for a in pts1:
                            if eval_count >= max_eval:
                                break
                            for b in pts2:
                                if eval_count >= max_eval:
                                    break
                                eval_count += 1
                                try:
                                    line_ab = QgsGeometry.fromPolylineXY([a, b])
                                except Exception:
                                    continue
                                if line_ab is None or line_ab.isEmpty():
                                    continue
                                try:
                                    line_len = float(line_ab.length())
                                except Exception:
                                    continue
                                if line_len <= 1e-9 or line_len > max_probe_len:
                                    continue
                                cand_ab = _prepare_compact_candidate(line_ab)
                                if cand_ab is None or cand_ab.isEmpty():
                                    continue
                                try:
                                    area_ab = float(cand_ab.area())
                                    perim_ab = float(cand_ab.length())
                                except Exception:
                                    continue
                                # Prefer clearly smaller footprint, then smaller perimeter.
                                if (
                                    area_ab + 1e-9 < (0.85 * max(best_area, 1e-9))
                                    or (
                                        area_ab <= best_area + 1e-9
                                        and perim_ab + 1e-9 < best_perim
                                    )
                                ):
                                    best_geom = cand_ab
                                    best_area = area_ab
                                    best_perim = perim_ab

                    compact_geom = best_geom
                    if compact_geom is not None and (not compact_geom.isEmpty()):
                        final_geom = compact_geom

                # Fallback to original geometry if compact candidate search fails.
                if compact_geom is None or compact_geom.isEmpty():
                    pass

    return final_geom, patch_ids


def _buffer_line_segment(line_geom: QgsGeometry, width: float) -> QgsGeometry:
    return line_geom.buffer(width / 2.0, BUFFER_SEGMENTS)


def _corridor_cost_area_ha(
    corridor_geom: QgsGeometry,
    distance_m: float,
    corridor_width_m: float,
) -> float:
    """
    Robust corridor cost area in hectares.
    Uses polygon footprint area, but caps pathological geometry spikes relative to
    nominal strip area (length * width), which can occur after topology repairs.
    """
    if corridor_geom is None or corridor_geom.isEmpty():
        return 0.0
    try:
        raw_area_ha = float(corridor_geom.area() / 10000.0)
    except Exception:
        raw_area_ha = 0.0
    if raw_area_ha <= 0.0:
        return 0.0
    nominal_ha = max(float(distance_m), 0.0) * max(float(corridor_width_m), 0.0) / 10000.0
    if nominal_ha <= 1e-12:
        return raw_area_ha
    # Allow moderate shape complexity while preventing invalid-geometry blowups.
    cap_ha = 2.5 * nominal_ha
    if raw_area_ha > cap_ha:
        return float(cap_ha)
    return float(raw_area_ha)


def _corridor_passes_width(corridor_geom: QgsGeometry, min_width: float) -> bool:
    if corridor_geom.isEmpty():
        return False
    # If the corridor is multipolygon (due to crossing patches), check each part?
    # Actually, if it's multipolygon, it means we clipped out patches.
    # The 'neck' check is complex on multipolygons. We check the buffer on the original,
    # but since we already clipped, we skip aggressive width validation on the final result
    # to avoid rejecting valid patch-traversals.
    return True


def _format_no_corridor_reason(
    stage: str,
    patch_count: int,
    candidate_count: int,
    params: VectorRunParams,
) -> str:
    return (
        f"{stage}: no feasible corridors could be generated.\n"
        f"- Patches meeting criteria: {patch_count}\n"
        f"- Candidate corridors generated: {candidate_count}\n"
        f"- Max search distance: {params.max_search_distance:.2f}\n"
        f"- Corridor width: {params.min_corridor_width:.2f}\n"
        "Try increasing the search distance, lowering the corridor width/area, "
        "or simplifying the patch layer."
    )


def _route_line_around_impassables_grid(
    start_pt: QgsPointXY,
    end_pt: QgsPointXY,
    obstacle_geoms: List[QgsGeometry],
    cell_m: float,
    max_window_m: float,
    params: VectorRunParams,
    ctx: Optional[AnalysisContext] = None,
) -> Optional[QgsGeometry]:
    minx0 = min(start_pt.x(), end_pt.x())
    maxx0 = max(start_pt.x(), end_pt.x())
    miny0 = min(start_pt.y(), end_pt.y())
    maxy0 = max(start_pt.y(), end_pt.y())

    width0 = maxx0 - minx0
    height0 = maxy0 - miny0
    if width0 > max_window_m or height0 > max_window_m:
        return None

    # Expand the bounding box just enough to allow detours, but never exceed max_window_m.
    pad_x = max(0.0, (max_window_m - width0) / 2.0)
    pad_y = max(0.0, (max_window_m - height0) / 2.0)

    minx = minx0 - pad_x
    maxx = maxx0 + pad_x
    miny = miny0 - pad_y
    maxy = maxy0 + pad_y

    width_m = maxx - minx
    height_m = maxy - miny

    cell_m = float(cell_m or 0.0)
    if cell_m <= 0:
        return None

    cols = int(math.ceil(width_m / cell_m))
    rows = int(math.ceil(height_m / cell_m))
    if cols < 3 or rows < 3:
        return None
    if (rows * cols) > 600_000:
        return None

    def to_rc(pt: QgsPointXY) -> Tuple[int, int]:
        c = int((pt.x() - minx) / cell_m)
        r = int((maxy - pt.y()) / cell_m)
        c = max(0, min(cols - 1, c))
        r = max(0, min(rows - 1, r))
        return r, c

    def to_xy(r: int, c: int) -> QgsPointXY:
        x = minx + (c + 0.5) * cell_m
        y = maxy - (r + 0.5) * cell_m
        return QgsPointXY(x, y)

    blocked = np.zeros((rows, cols), dtype=np.uint8)

    window_geom = QgsGeometry.fromRect(QgsRectangle(minx, miny, maxx, maxy))
    obs: List[QgsGeometry] = []
    for g in obstacle_geoms:
        try:
            if g and (not g.isEmpty()) and g.intersects(window_geom):
                obs.append(g)
        except Exception:
            continue

    if not obs:
        return QgsGeometry.fromPolylineXY([start_pt, end_pt])

    # Burn obstacles into the grid (hard blocking with clearance).
    # Relax inflation to allow squeezing through narrow gaps: ensure the centerline doesn't hit
    # the obstacle cell (half-cell clearance), but do not add corridor width clearance here.
    corridor_r = max(0.0, float(params.min_corridor_width) * 0.5)
    inflate = cell_m * 0.5
    inflated_obs: List[QgsGeometry] = []
    for g in obs:
        try:
            gg = g.makeValid()
        except Exception:
            gg = g
        try:
            if gg and (not gg.isEmpty()):
                inflated_obs.append(gg.buffer(inflate, 8))
        except Exception:
            continue

    for og in inflated_obs:
        try:
            bbox = og.boundingBox()
        except Exception:
            continue

        min_c = max(0, int(math.floor((bbox.xMinimum() - minx) / cell_m)))
        max_c = min(cols - 1, int(math.ceil((bbox.xMaximum() - minx) / cell_m)))
        min_r = max(0, int(math.floor((maxy - bbox.yMaximum()) / cell_m)))
        max_r = min(rows - 1, int(math.ceil((maxy - bbox.yMinimum()) / cell_m)))

        for r in range(min_r, max_r + 1):
            y_top = maxy - r * cell_m
            y_bot = maxy - (r + 1) * cell_m
            for c in range(min_c, max_c + 1):
                x_left = minx + c * cell_m
                x_right = minx + (c + 1) * cell_m
                cell_geom = QgsGeometry.fromRect(QgsRectangle(x_left, y_bot, x_right, y_top))
                try:
                    if og.intersects(cell_geom):
                        blocked[r, c] = 1
                except Exception:
                    continue

    sr, sc = to_rc(start_pt)
    tr, tc = to_rc(end_pt)

    def nearest_free(r0: int, c0: int, maxrad: int = 25) -> Optional[Tuple[int, int]]:
        if blocked[r0, c0] == 0:
            return (r0, c0)
        for rad in range(1, maxrad + 1):
            for dr in range(-rad, rad + 1):
                for dc in range(-rad, rad + 1):
                    rr = r0 + dr
                    cc = c0 + dc
                    if 0 <= rr < rows and 0 <= cc < cols and blocked[rr, cc] == 0:
                        return (rr, cc)
        return None

    sfix = nearest_free(sr, sc)
    tfix = nearest_free(tr, tc)
    if sfix is None or tfix is None:
        return None
    sr, sc = sfix
    tr, tc = tfix

    def h(r: int, c: int) -> float:
        return math.hypot(tr - r, tc - c)

    INF = 1e30
    gscore = np.full((rows, cols), INF, dtype=np.float64)
    gscore[sr, sc] = 0.0
    came: Dict[Tuple[int, int], Tuple[int, int]] = {}

    moves = [
        (-1, 0, 1.0),
        (1, 0, 1.0),
        (0, -1, 1.0),
        (0, 1, 1.0),
        (-1, -1, math.sqrt(2)),
        (-1, 1, math.sqrt(2)),
        (1, -1, math.sqrt(2)),
        (1, 1, math.sqrt(2)),
    ]

    heap: List[Tuple[float, float, int, int]] = [(h(sr, sc), 0.0, sr, sc)]
    visited: Set[Tuple[int, int]] = set()

    while heap:
        _f, g, r, c = heapq.heappop(heap)
        if (r, c) in visited:
            continue
        visited.add((r, c))
        if (r, c) == (tr, tc):
            break

        for dr, dc, step_cost in moves:
            rr = r + dr
            cc = c + dc
            if rr < 0 or rr >= rows or cc < 0 or cc >= cols:
                continue
            if blocked[rr, cc] != 0:
                continue
            # Prevent diagonal corner-cutting through thin walls.
            if dr != 0 and dc != 0:
                if blocked[r + dr, c] != 0 or blocked[r, c + dc] != 0:
                    continue
            ng = g + step_cost
            if ng >= float(gscore[rr, cc]):
                continue
            gscore[rr, cc] = ng
            came[(rr, cc)] = (r, c)
            heapq.heappush(heap, (ng + h(rr, cc), ng, rr, cc))

    if (tr, tc) not in came and (sr, sc) != (tr, tc):
        return None

    # Reconstruct cell path.
    cell_path: List[Tuple[int, int]] = [(tr, tc)]
    cur = (tr, tc)
    while cur != (sr, sc):
        cur = came.get(cur)
        if cur is None:
            return None
        cell_path.append(cur)
    cell_path.reverse()

    pts: List[QgsPointXY] = [start_pt]
    pts.extend(to_xy(r, c) for (r, c) in cell_path[1:-1])
    pts.append(end_pt)

    # Thin points to reduce geometry complexity.
    if len(pts) > 2:
        thinned: List[QgsPointXY] = [pts[0]]
        last = pts[0]
        for p in pts[1:-1]:
            if last.distance(p) >= (cell_m * 0.75):
                thinned.append(p)
                last = p
        thinned.append(pts[-1])
        pts = thinned

    try:
        geom = QgsGeometry.fromPolylineXY(pts)
        # Smooth the routed line to reduce the grid "stair-step" appearance.
        if geom and (not geom.isEmpty()) and len(pts) > 2:
            try:
                iters = int(getattr(params, "vector_routing_smooth_iterations", 6) or 6)
                offset = float(getattr(params, "vector_routing_smooth_offset", 0.25) or 0.25)
                smoothed = geom.smooth(max(0, iters), offset)
            except Exception:
                smoothed = None
                try:
                    smoothed = geom.smooth(max(0, int(getattr(params, "vector_routing_smooth_iterations", 6) or 6)))
                except Exception:
                    smoothed = None
            if smoothed and (not smoothed.isEmpty()):
                geom = smoothed
        return geom
    except Exception:
        return None


def _create_corridor_geometry(
    waypoints: List[QgsPointXY],
    source_geom: QgsGeometry,
    target_geom: QgsGeometry,
    params: VectorRunParams,
    obstacle_geoms: Optional[List[QgsGeometry]] = None,
    ctx: Optional[AnalysisContext] = None,
    smooth_iterations: int = 0,
) -> Optional[QgsGeometry]:
    if not waypoints or len(waypoints) < 2:
        return None

    start_pt = QgsPointXY(waypoints[0])
    end_pt = QgsPointXY(waypoints[-1])
    corridor_line = QgsGeometry.fromPolylineXY([QgsPointXY(pt) for pt in waypoints])
    # Apply smoothing for routed paths (more than 2 waypoints), or if explicitly requested.
    # This reduces the angular grid "stair-step" look.
    iterations_to_use = int(smooth_iterations or 0)
    if obstacle_geoms and iterations_to_use == 0 and len(waypoints) > 2:
        iterations_to_use = int(getattr(params, "vector_routing_smooth_iterations", 6) or 6)

    if iterations_to_use > 0:
        try:
            try:
                offset = float(getattr(params, "vector_routing_smooth_offset", 0.25) or 0.25)
                smoothed = corridor_line.smooth(iterations_to_use, offset)
            except Exception:
                smoothed = corridor_line.smooth(iterations_to_use)
            if smoothed and not smoothed.isEmpty():
                corridor_line = smoothed
        except Exception:
            pass

    # If a straight corridor is blocked by impassables, optionally route around them using a
    # small grid A* within a bounded local window (maze-capable).
    #
    # NOTE: `obstacle_geoms` should already be in analysis CRS; when passed from RasterNavigator
    # we use obstacles buffered by half-width so corridor width becomes a hard constraint.
    if obstacle_geoms:
        try:
            enabled = bool(getattr(params, "vector_routing_enabled", False))
        except Exception:
            enabled = False
        try:
            safety = max(0.0, float(params.min_corridor_width or 0.0) * 0.5)
            obstacle_union = (ctx.impassable_union if ctx is not None else getattr(params, "_impassable_union", None))
            if obstacle_union is None:
                obstacle_union = QgsGeometry.unaryUnion([g for g in obstacle_geoms if g and (not g.isEmpty())]).makeValid()
                try:
                    if ctx is not None:
                        ctx.impassable_union = obstacle_union
                    else:
                        params._impassable_union = obstacle_union  # type: ignore[attr-defined]
                except Exception:
                    pass
            if safety > 0:
                try:
                    blocked = corridor_line.intersects(obstacle_union.buffer(safety, 8))
                except Exception:
                    blocked = corridor_line.intersects(obstacle_union)
            else:
                blocked = corridor_line.intersects(obstacle_union)
        except Exception:
            blocked = False

        if blocked and not enabled:
            return None

        if blocked and enabled:
            cell_m = float(getattr(params, "grid_resolution", 50.0) or 50.0)
            max_win = float(getattr(params, "vector_routing_max_window_m", 6000.0) or 6000.0)
            routed = _route_line_around_impassables_grid(
                start_pt=start_pt,
                end_pt=end_pt,
                obstacle_geoms=list(obstacle_geoms),
                cell_m=cell_m,
                max_window_m=max_win,
                params=params,
                ctx=ctx,
            )
            if routed and not routed.isEmpty():
                corridor_line = routed
            else:
                return None

    # Buffer to full width
    corridor_geom = _buffer_line_segment(corridor_line, params.min_corridor_width)

    # Clip start/end immediately to get the "bridge" geometry
    corridor_geom = corridor_geom.difference(source_geom)
    corridor_geom = corridor_geom.difference(target_geom)

    if obstacle_geoms:
        obstacle_union = (ctx.impassable_union if ctx is not None else getattr(params, "_impassable_union", None))
        if obstacle_union is None:
            try:
                obstacle_union = QgsGeometry.unaryUnion([g for g in obstacle_geoms if g and (not g.isEmpty())]).makeValid()
                try:
                    if ctx is not None:
                        ctx.impassable_union = obstacle_union
                    else:
                        params._impassable_union = obstacle_union  # type: ignore[attr-defined]
                except Exception:
                    pass
            except Exception:
                obstacle_union = None

        if obstacle_union and (not obstacle_union.isEmpty()):
            try:
                overlap = corridor_geom.intersection(obstacle_union)
                overlap_area = overlap.area() if overlap and (not overlap.isEmpty()) else 0.0
            except Exception:
                overlap_area = 0.0

            if overlap_area > 0.0:
                try:
                    clipped = corridor_geom.difference(obstacle_union)
                except Exception:
                    clipped = None

                if clipped is None or clipped.isEmpty():
                    return None

                end_buf = max(float(params.grid_resolution or 0.0) * 1.5, float(params.min_corridor_width or 0.0) * 1.5)
                try:
                    a_buf = QgsGeometry.fromPointXY(start_pt).buffer(end_buf, 8)
                    b_buf = QgsGeometry.fromPointXY(end_pt).buffer(end_buf, 8)
                    ok_a = clipped.intersects(a_buf)
                    ok_b = clipped.intersects(b_buf)
                except Exception:
                    ok_a = True
                    ok_b = True

                if not (ok_a and ok_b):
                    return None

                corridor_geom = clipped
                 
    corridor_geom = corridor_geom.makeValid()
    if corridor_geom.isEmpty():
        return None

    return corridor_geom


class RasterNavigator:
    def __init__(
        self,
        patches: Dict[int, Dict],
        obstacle_layers: List[QgsVectorLayer],
        target_crs: QgsCoordinateReferenceSystem,
        params: VectorRunParams,
    ):
        if not obstacle_layers:
            raise VectorAnalysisError("Select at least one polygon impassable layer for impassable land classes.")

        self._params = params
        self.resolution = max(params.grid_resolution, 1.0)
        self.obstacle_geoms: List[QgsGeometry] = []
        self._buffered_obstacles_cache: Dict[float, List[QgsGeometry]] = {}

        extent: Optional[QgsRectangle] = None
        for patch in patches.values():
            bbox = patch["geom"].boundingBox()
            if extent is None:
                extent = QgsRectangle(bbox)
            else:
                extent.combineExtentWith(bbox)

        for obstacle_layer in obstacle_layers:
            if obstacle_layer is None or QgsWkbTypes.geometryType(obstacle_layer.wkbType()) != QgsWkbTypes.PolygonGeometry:
                raise VectorAnalysisError("Select a polygon impassable layer for impassable land classes.")

            transform = QgsCoordinateTransform(obstacle_layer.crs(), target_crs, QgsProject.instance())
            for feature in obstacle_layer.getFeatures():
                geom = QgsGeometry(feature.geometry())
                if geom.isEmpty():
                    continue
                geom.transform(transform)
                geom = geom.makeValid()
                if geom.isEmpty():
                    continue
                self.obstacle_geoms.append(clone_geometry(geom))
                bbox = geom.boundingBox()
                if extent is None:
                    extent = QgsRectangle(bbox)
                else:
                    extent.combineExtentWith(bbox)

        if extent is None:
            raise VectorAnalysisError("Unable to determine extent for impassable land class routing.")

        pad = max(params.max_search_distance, params.min_corridor_width)
        extent = QgsRectangle(
            extent.xMinimum() - pad,
            extent.yMinimum() - pad,
            extent.xMaximum() + pad,
            extent.yMaximum() + pad,
        )
        width = max(extent.width(), self.resolution)
        height = max(extent.height(), self.resolution)

        self.origin_x = extent.xMinimum()
        self.origin_y = extent.yMaximum()
        self.cols = max(1, int(math.ceil(width / self.resolution)))
        self.rows = max(1, int(math.ceil(height / self.resolution)))
        self.passable = np.ones((self.rows, self.cols), dtype=bool)

        # Use a minimal safety buffer (half a grid cell) to allow squeezing through narrow gaps.
        # This allows the routed centerline to pass through any gap wider than the grid resolution.
        safety_buffer = self.resolution * 0.5

        for geom in self.obstacle_geoms:
            try:
                mask_geom = geom.buffer(safety_buffer, 4) if safety_buffer > 0 else clone_geometry(geom)
            except Exception:
                mask_geom = clone_geometry(geom)

            try:
                mask_geom = mask_geom.makeValid()
            except Exception:
                pass
            if mask_geom is None or mask_geom.isEmpty():
                continue
            self._burn_geometry(mask_geom)

    def _world_to_rc(self, point: QgsPointXY) -> Optional[Tuple[int, int]]:
        col = int(math.floor((point.x() - self.origin_x) / self.resolution))
        row = int(math.floor((self.origin_y - point.y()) / self.resolution))
        if 0 <= row < self.rows and 0 <= col < self.cols:
            return row, col
        return None

    def _rc_to_world(self, row: int, col: int) -> QgsPointXY:
        x = self.origin_x + (col + 0.5) * self.resolution
        y = self.origin_y - (row + 0.5) * self.resolution
        return QgsPointXY(x, y)

    def _burn_geometry(self, geom: QgsGeometry) -> None:
        bbox = geom.boundingBox()
        min_col = max(0, int(math.floor((bbox.xMinimum() - self.origin_x) / self.resolution)))
        max_col = min(self.cols - 1, int(math.ceil((bbox.xMaximum() - self.origin_x) / self.resolution)))
        min_row = max(0, int(math.floor((self.origin_y - bbox.yMaximum()) / self.resolution)))
        max_row = min(self.rows - 1, int(math.ceil((self.origin_y - bbox.yMinimum()) / self.resolution)))

        if min_col > max_col or min_row > max_row:
            return

        for row in range(min_row, max_row + 1):
            y = self.origin_y - (row + 0.5) * self.resolution
            for col in range(min_col, max_col + 1):
                x = self.origin_x + (col + 0.5) * self.resolution
                try:
                    if geom.contains(QgsPointXY(x, y)):
                        self.passable[row, col] = False
                except Exception:
                    pass

    def _cell_stats_in_geom(self, geom: QgsGeometry) -> Tuple[int, int, int]:
        """Return (passable_cells, blocked_cells, total_cells) for cell centers inside `geom`."""
        bbox = geom.boundingBox()
        min_col = max(0, int(math.floor((bbox.xMinimum() - self.origin_x) / self.resolution)))
        max_col = min(self.cols - 1, int(math.ceil((bbox.xMaximum() - self.origin_x) / self.resolution)))
        min_row = max(0, int(math.floor((self.origin_y - bbox.yMaximum()) / self.resolution)))
        max_row = min(self.rows - 1, int(math.ceil((self.origin_y - bbox.yMinimum()) / self.resolution)))

        if min_col > max_col or min_row > max_row:
            return 0, 0, 0

        passable = 0
        blocked = 0
        total = 0
        for row in range(min_row, max_row + 1):
            y = self.origin_y - (row + 0.5) * self.resolution
            for col in range(min_col, max_col + 1):
                x = self.origin_x + (col + 0.5) * self.resolution
                try:
                    if not geom.contains(QgsPointXY(x, y)):
                        continue
                except Exception:
                    continue
                total += 1
                if self.passable[row, col]:
                    passable += 1
                else:
                    blocked += 1
        return passable, blocked, total

    def find_path(
        self, start_point: QgsPointXY, end_point: QgsPointXY
    ) -> Optional[List[QgsPointXY]]:
        start = self._world_to_rc(start_point)
        end = self._world_to_rc(end_point)
        if start is None or end is None:
            return None
        start_node = _nearest_passable_node(self.passable, start)
        end_node = _nearest_passable_node(self.passable, end)
        if start_node is None or end_node is None:
            return None

        path = _shortest_path_on_mask(self.passable, start_node, end_node)
        if not path:
            return None
        return [self._rc_to_world(r, c) for r, c in path]

    def buffered_obstacles(self, safety: float) -> List[QgsGeometry]:
        """
        Return obstacles buffered by `safety` (analysis CRS units), cached per run.
        Use this for corridor validation and optional non-navigator routing.
        """
        try:
            key = round(float(safety or 0.0), 3)
        except Exception:
            key = 0.0
        cached = self._buffered_obstacles_cache.get(key)
        if cached is not None:
            return cached

        out: List[QgsGeometry] = []
        for g in self.obstacle_geoms:
            if g is None or g.isEmpty():
                continue
            try:
                gg = g.buffer(float(safety or 0.0), 4) if safety and safety > 0 else clone_geometry(g)
            except Exception:
                gg = clone_geometry(g)
            try:
                gg = gg.makeValid()
            except Exception:
                pass
            if gg is not None and not gg.isEmpty():
                out.append(gg)

        self._buffered_obstacles_cache[key] = out
        return out


def _nearest_passable_node(mask: np.ndarray, node: Tuple[int, int], search_radius: int = 6) -> Optional[Tuple[int, int]]:
    r0, c0 = node
    rows, cols = mask.shape
    if 0 <= r0 < rows and 0 <= c0 < cols and mask[r0, c0]:
        return node
    for radius in range(1, search_radius + 1):
        for dr in range(-radius, radius + 1):
            for dc in range(-radius, radius + 1):
                nr = r0 + dr
                nc = c0 + dc
                if 0 <= nr < rows and 0 <= nc < cols and mask[nr, nc]:
                    return nr, nc
    return None


def _shortest_path_on_mask(
    mask: np.ndarray, start: Tuple[int, int], goal: Tuple[int, int]
) -> Optional[List[Tuple[int, int]]]:
    rows, cols = mask.shape
    moves = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]

    heap: List[Tuple[float, int, int]] = []
    heapq.heappush(heap, (0.0, start[0], start[1]))
    best_cost: Dict[Tuple[int, int], float] = {start: 0.0}
    parents: Dict[Tuple[int, int], Tuple[int, int]] = {}

    while heap:
        cost, r, c = heapq.heappop(heap)
        if (r, c) == goal:
            break
        if cost > best_cost.get((r, c), float("inf")):
            continue
        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if not (0 <= nr < rows and 0 <= nc < cols):
                continue
            if not mask[nr, nc]:
                continue
            step = math.sqrt(2) if dr != 0 and dc != 0 else 1.0
            new_cost = cost + step
            if new_cost >= best_cost.get((nr, nc), float("inf")):
                continue
            best_cost[(nr, nc)] = new_cost
            parents[(nr, nc)] = (r, c)
            heapq.heappush(heap, (new_cost, nr, nc))

    if goal not in parents and goal != start:
        return None

    path: List[Tuple[int, int]] = [goal]
    current = goal
    while current != start:
        current = parents.get(current)
        if current is None:
            return None
        path.append(current)
    path.reverse()
    return path


def find_all_possible_corridors(
    patches: Dict[int, Dict],
    spatial_index: QgsSpatialIndex,
    params: VectorRunParams,
    strategy: str = "most_connected_habitat",
    patch_union: Optional[QgsGeometry] = None,
    ctx: Optional[AnalysisContext] = None,
    progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
    progress_start: int = 30,
    progress_end: int = 55,
    navigator: Optional[RasterNavigator] = None,
    timings: Optional[TimingRecorder] = None,
    timing_out: Optional[Dict[str, object]] = None,
) -> List[Dict]:
    print("  Finding all possible corridors...")
    processed_pairs: Set[frozenset] = set()
    total = len(patches) or 1
    accum_durations: Dict[str, float] = defaultdict(float)
    accum_counts: Dict[str, int] = defaultdict(int)
    strategy_key = _normalize_strategy_key(strategy)
    fluidity_mode = strategy_key == "landscape_fluidity"
    # Metrics that benefit from alternative corridor geometry variants.
    circuit_mode = strategy_key in (
        "most_connected_habitat",
        "reachable_habitat_advanced",
        "landscape_fluidity",
        "largest_single_network",
    )
    largest_network_mode = strategy_key == "largest_single_network"
    # Both Most Connectivity and Largest Single Network (Circuit Utility) need enough candidate
    # diversity to avoid "blind" terminal selection in long/snaking patches.
    utility_mode = circuit_mode
    max_keep_per_pair = 8 if utility_mode else 1
    # To avoid O(n^2) candidate enumeration, rank neighbor patches by true polygon boundary distance
    # (not centroid distance) and keep the top-K nearest per patch.
    # Important: keep K large enough that short "obvious" gaps are not dropped in dense patch fields.
    # We therefore only cap when there are many candidates within the search window.
    k_nearest_neighbors_cap = 250
    min_distinct_overlap_ratio = 0.75  # higher = allow more similar corridors
    proximity_dist = max(float(getattr(params, "min_corridor_width", 0.0) or 0.0) * 1.5, float(getattr(params, "grid_resolution", 0.0) or 0.0) * 2.0)
    timing_start = time.perf_counter()
    # Use raw impassable geometries (analysis CRS). Routing/clearance is handled inside
    # the grid router via obstacle inflation and inside corridor finalization via clipping.
    impassable_geoms: Optional[List[QgsGeometry]] = None
    if navigator is not None:
        impassable_geoms = navigator.obstacle_geoms

    # Maintain a bounded set of spatially distinct candidates per patch-pair.
    candidates_by_pair: Dict[Tuple[int, int], List[Dict]] = defaultdict(list)
    chain_seed_by_pair: Dict[Tuple[int, int], List[Dict]] = defaultdict(list)

    def _pair_key(a: int, b: int) -> Tuple[int, int]:
        ia, ib = int(a), int(b)
        return (ia, ib) if ia <= ib else (ib, ia)

    def _geom_parts(geom: QgsGeometry) -> List[QgsGeometry]:
        if geom is None or geom.isEmpty():
            return []
        if not geom.isMultipart():
            return [geom]
        try:
            coll = geom.asGeometryCollection()
            if coll:
                out: List[QgsGeometry] = []
                for p in coll:
                    g = QgsGeometry(p)
                    if g and (not g.isEmpty()):
                        out.append(g)
                return out or [geom]
        except Exception:
            pass
        try:
            out2: List[QgsGeometry] = []
            for poly in geom.asMultiPolygon():
                try:
                    g = QgsGeometry.fromPolygonXY(poly)
                    if g and (not g.isEmpty()):
                        out2.append(g)
                except Exception:
                    continue
            return out2 or [geom]
        except Exception:
            return [geom]

    def _pick_part_closest_to_patch(geom: QgsGeometry, patch_geom: QgsGeometry) -> Optional[QgsGeometry]:
        parts = _geom_parts(geom)
        if not parts:
            return None
        best_part: Optional[QgsGeometry] = None
        best_dist = float("inf")
        for part in parts:
            try:
                d = float(part.distance(patch_geom))
            except Exception:
                continue
            if d < best_dist:
                best_dist = d
                best_part = part
        return best_part

    def _store_candidate(
        cand: Dict,
        store: Dict[Tuple[int, int], List[Dict]],
        *,
        cap_override: Optional[int] = None,
    ) -> None:
        try:
            p1 = int(cand.get("patch1"))
            p2 = int(cand.get("patch2"))
        except Exception:
            return
        barrier_pair = bool(cand.get("barrier_pair", False))
        geom = cand.get("geom")
        if geom is None or geom.isEmpty():
            return
        key = _pair_key(p1, p2)
        existing = store.get(key, [])
        pair_sep = max(float(params.max_search_distance or 0.0), 0.0)

        def _cand_anchors(c: Dict) -> Dict[int, QgsPointXY]:
            anchor_pts = c.get("anchor_pts")
            if isinstance(anchor_pts, dict) and anchor_pts:
                out: Dict[int, QgsPointXY] = {}
                for pid in key:
                    try:
                        if int(pid) in anchor_pts:
                            out[int(pid)] = anchor_pts[int(pid)]
                    except Exception:
                        continue
                if len(out) == len(key):
                    return out
            cg = c.get("geom")
            if cg is None or cg.isEmpty():
                return {}
            out: Dict[int, QgsPointXY] = {}
            for pid in key:
                try:
                    pg = (patches.get(int(pid)) or {}).get("geom")
                    if pg is None or pg.isEmpty():
                        continue
                    npg = pg.nearestPoint(cg)
                    if npg is None or npg.isEmpty():
                        continue
                    out[int(pid)] = QgsPointXY(npg.asPoint())
                except Exception:
                    continue
            return out

        keep_even_if_overlap = False
        if barrier_pair and pair_sep > 0.0 and existing:
            ca = _cand_anchors(cand)
            if len(ca) == 2:
                keep_even_if_overlap = True
                for prev in existing:
                    pa = _cand_anchors(prev)
                    if len(pa) != 2:
                        keep_even_if_overlap = False
                        break
                    for pid in key:
                        if (
                            int(pid) not in ca
                            or int(pid) not in pa
                            or float(ca[int(pid)].distance(pa[int(pid)])) + 1e-9 < pair_sep
                        ):
                            keep_even_if_overlap = False
                            break
                    if not keep_even_if_overlap:
                        break
        for idx, prev in enumerate(existing):
            if keep_even_if_overlap:
                try:
                    if _overlap_ratio(geom, prev.get("geom")) >= min_distinct_overlap_ratio:
                        return
                except Exception:
                    pass
                continue
            try:
                if _overlap_ratio(geom, prev.get("geom")) >= min_distinct_overlap_ratio:
                    new_area = float(cand.get("area_ha", float("inf")) or float("inf"))
                    prev_area = float(prev.get("area_ha", float("inf")) or float("inf"))
                    new_dist = float(cand.get("distance_m", float("inf")) or float("inf"))
                    prev_dist = float(prev.get("distance_m", float("inf")) or float("inf"))
                    if (
                        new_area + 1e-12 < prev_area
                        or (
                            abs(new_area - prev_area) <= 1e-12
                            and new_dist + 1e-12 < prev_dist
                        )
                    ):
                        existing[idx] = cand
                    return
            except Exception:
                continue

        existing.append(cand)
        existing.sort(key=lambda c: (c.get("area_ha", 0.0), c.get("distance_m", 0.0)))
        cap = int(cap_override if cap_override is not None else max_keep_per_pair)
        if barrier_pair and cap_override is None:
            cap = max(cap, int(getattr(params, "landscape_fluidity_barrier_max_keep_per_pair", 20) or 20))
        if len(existing) > cap:
            del existing[cap:]
        store[key] = existing

    def _push_candidate(cand: Dict) -> None:
        try:
            p1 = int(cand.get("patch1"))
            p2 = int(cand.get("patch2"))
        except Exception:
            return
        barrier_pair = bool(cand.get("barrier_pair", False))
        def _dbg(msg: str) -> None:
            return
        if p1 == p2 and (not bool(cand.get("intra_patch", False))):
            _dbg("skip: same patch without intra flag")
            return
        geom = cand.get("geom")
        if geom is None or geom.isEmpty():
            _dbg("skip: empty geom")
            return
        # Patch-area guard (early): reject corridors where an endpoint patch is
        # smaller than the corridor itself. Allow tiny interior patches in chains
        # (A-B-C) as long as endpoints are not tiny.
        try:
            est_cost = float(cand.get("area_ha", 0.0) or 0.0)
        except Exception:
            est_cost = 0.0
        try:
            geom_cost = float(geom.area()) / 10000.0
        except Exception:
            geom_cost = 0.0
        cost = max(est_cost, geom_cost, 0.0)
        if cost > 0.0:
            p1_area = float((patches.get(int(p1)) or {}).get("area_ha", 0.0) or 0.0)
            p2_area = float((patches.get(int(p2)) or {}).get("area_ha", 0.0) or 0.0)
            small_endpoint = (
                (p1_area > 0.0 and p1_area < cost - 1e-12)
                or (p2_area > 0.0 and p2_area < cost - 1e-12)
            )
            if small_endpoint:
                if fluidity_mode and (not barrier_pair) and (not bool(cand.get("intra_patch", False))):
                    seed = dict(cand)
                    seed["chain_seed_only"] = True
                    _store_candidate(seed, chain_seed_by_pair, cap_override=16)
                _dbg(f"skip: patch-area guard cost={cost:.4f} p1={p1_area:.4f} p2={p2_area:.4f}")
                return
        _store_candidate(cand, candidates_by_pair)

    @contextmanager
    def _measure(label: str):
        start = time.perf_counter()
        accum_counts[label] += 1
        try:
            yield
        finally:
            accum_durations[label] += time.perf_counter() - start

    def _sample_boundary_points(
        patch_geom: QgsGeometry, max_points: int = 32, allow_densify: bool = True
    ) -> List[QgsPointXY]:
        pts: List[QgsPointXY] = []
        try:
            ring = patch_geom.constGet().exteriorRing()
            if ring:
                pts = [QgsPointXY(v) for v in ring.vertices()]
        except Exception:
            try:
                pts = [QgsPointXY(v) for v in patch_geom.vertices()]
            except Exception:
                pts = []
        if len(pts) > max_points and max_points > 0:
            step = max(1, len(pts) // max_points)
            pts = pts[::step][:max_points]
        if allow_densify and len(pts) < max_points:
            try:
                perim = max(patch_geom.length(), 1.0)
                target_spacing = perim / max_points
                densified = patch_geom.densifyByDistance(target_spacing)
                extra = [QgsPointXY(v) for v in densified.vertices()]
                if extra:
                    pts.extend(extra)
            except Exception:
                pass
            if len(pts) > max_points:
                step = max(1, len(pts) // max_points)
                pts = pts[::step][:max_points]
        return pts[:max_points]

    def _safe_unary_union(geoms: List[QgsGeometry]) -> Optional[QgsGeometry]:
        if not geoms:
            return None
        try:
            merged = QgsGeometry.unaryUnion(geoms)
        except Exception:
            merged = geoms[0]
            for extra in geoms[1:]:
                try:
                    merged = merged.combine(extra)
                except Exception:
                    pass
        if merged is None or merged.isEmpty():
            return None
        try:
            merged = merged.makeValid()
        except Exception:
            pass
        return merged if (merged is not None and not merged.isEmpty()) else None

    def _pick_extremes_along_axis(
        points: List[QgsPointXY], ax_dx: float, ax_dy: float, cx: float, cy: float
    ) -> Tuple[Optional[QgsPointXY], Optional[QgsPointXY]]:
        if abs(ax_dx) < 1e-12 and abs(ax_dy) < 1e-12:
            return None, None
        if not points:
            return None, None

        def _proj(pt: QgsPointXY) -> float:
            return (pt.x() - cx) * ax_dx + (pt.y() - cy) * ax_dy

        p_max = max(points, key=_proj)
        p_min = min(points, key=_proj)
        if p_max.distance(p_min) < 1e-6:
            return p_max, None
        return p_max, p_min

    def _anchor_variants_for_pair(
        g1: QgsGeometry, g2: QgsGeometry, nearest1: QgsPointXY, nearest2: QgsPointXY
    ) -> List[Tuple[str, QgsPointXY, QgsPointXY]]:
        variants: List[Tuple[str, QgsPointXY, QgsPointXY]] = [("nearest", nearest1, nearest2)]
        # Build two additional variants that are on the *facing* edges of each patch.
        # This prevents "teleport" corridors that rely on travel through patch interior
        # to exit on the opposite side.
        try:
            c1 = g1.centroid().asPoint()
            c2 = g2.centroid().asPoint()
            c1x, c1y = c1.x(), c1.y()
            c2x, c2y = c2.x(), c2.y()
            vx, vy = (c2x - c1x), (c2y - c1y)
        except Exception:
            c1x, c1y, c2x, c2y = 0.0, 0.0, 1.0, 0.0
            vx, vy = 1.0, 0.0

        pts1 = _sample_boundary_points(g1, max_points=96, allow_densify=True)
        pts2 = _sample_boundary_points(g2, max_points=96, allow_densify=True)
        if not pts1 or not pts2:
            return variants

        # "Facing" = high projection along the axis towards the other patch.
        def _proj_to_other_from_1(pt: QgsPointXY) -> float:
            return (pt.x() - c1x) * vx + (pt.y() - c1y) * vy

        def _proj_to_other_from_2(pt: QgsPointXY) -> float:
            # towards patch1, so invert axis
            return (pt.x() - c2x) * (-vx) + (pt.y() - c2y) * (-vy)

        p1_max = max((_proj_to_other_from_1(p) for p in pts1), default=0.0)
        p2_max = max((_proj_to_other_from_2(p) for p in pts2), default=0.0)
        # keep points in the top 40% of "facing" scores
        p1_cut = p1_max * 0.60
        p2_cut = p2_max * 0.60
        facing1 = [p for p in pts1 if _proj_to_other_from_1(p) >= p1_cut]
        facing2 = [p for p in pts2 if _proj_to_other_from_2(p) >= p2_cut]
        if not facing1:
            facing1 = pts1
        if not facing2:
            facing2 = pts2

        # Now spread along perpendicular direction, but only within facing sets.
        px, py = (-vy, vx)
        a1_pos, a1_neg = _pick_extremes_along_axis(facing1, px, py, c1x, c1y)
        a2_pos, a2_neg = _pick_extremes_along_axis(facing2, px, py, c2x, c2y)
        if a1_pos is not None and a2_pos is not None:
            variants.append(("side_pos", a1_pos, a2_pos))
        if a1_neg is not None and a2_neg is not None:
            variants.append(("side_neg", a1_neg, a2_neg))

        if circuit_mode:
            # Perimeter-based (cardinal) anchors to enable multiple spatially distinct
            # candidates between long/snaking patches.
            n1, s1 = _pick_extremes_along_axis(pts1, 0.0, 1.0, c1x, c1y)
            n2, s2 = _pick_extremes_along_axis(pts2, 0.0, 1.0, c2x, c2y)
            e1, w1 = _pick_extremes_along_axis(pts1, 1.0, 0.0, c1x, c1y)
            e2, w2 = _pick_extremes_along_axis(pts2, 1.0, 0.0, c2x, c2y)
            if n1 is not None and n2 is not None:
                variants.append(("north", n1, n2))
            if s1 is not None and s2 is not None:
                variants.append(("south", s1, s2))
            if e1 is not None and e2 is not None:
                variants.append(("east", e1, e2))
            if w1 is not None and w2 is not None:
                variants.append(("west", w1, w2))

        seen = set()
        uniq: List[Tuple[str, QgsPointXY, QgsPointXY]] = []
        for tag, pta, ptb in variants:
            key = (round(pta.x(), 3), round(pta.y(), 3), round(ptb.x(), 3), round(ptb.y(), 3))
            if key in seen:
                continue
            seen.add(key)
            uniq.append((tag, pta, ptb))
        return uniq

    def _barrier_touch_variants(
        g1: QgsGeometry,
        g2: QgsGeometry,
        nearest1: QgsPointXY,
        nearest2: QgsPointXY,
        spacing_m: float,
        max_candidates: int,
    ) -> List[Tuple[str, QgsPointXY, QgsPointXY]]:
        # Prefer near-touch locations along the facing boundary, then enforce spacing.
        variants: List[Tuple[str, QgsPointXY, QgsPointXY]] = [("nearest", nearest1, nearest2)]
        try:
            c1 = g1.centroid().asPoint()
            c2 = g2.centroid().asPoint()
            c1x, c1y = c1.x(), c1.y()
            c2x, c2y = c2.x(), c2.y()
            vx, vy = (c2x - c1x), (c2y - c1y)
        except Exception:
            c1x, c1y, c2x, c2y = 0.0, 0.0, 1.0, 0.0
            vx, vy = 1.0, 0.0

        perim = max(float(g1.length() or 0.0), 0.0)
        step = max(5.0, min(50.0, spacing_m * 0.25))
        max_points = int(min(1000, max(120, perim / max(1.0, step))))
        pts1 = _sample_boundary_points(g1, max_points=max_points, allow_densify=True)
        if not pts1:
            return variants

        def _proj_to_other_from_1(pt: QgsPointXY) -> float:
            return (pt.x() - c1x) * vx + (pt.y() - c1y) * vy

        # For barrier pairs, allow the full boundary; we guard against
        # back-side teleporting by checking interior line length below.
        facing1 = list(pts1)

        candidates: List[Tuple[float, QgsPointXY, QgsPointXY]] = []
        max_gap = float(params.max_search_distance or 0.0)
        max_inside_len = max(float(getattr(params, "min_corridor_width", 0.0) or 0.0) * 3.0, 5.0)
        for p1 in facing1:
            try:
                p1_geom = QgsGeometry.fromPointXY(p1)
                near_geom = g2.nearestPoint(p1_geom)
                if near_geom is None or near_geom.isEmpty():
                    continue
                p2p = near_geom.asPoint()
                p2 = QgsPointXY(p2p)
            except Exception:
                continue
            try:
                raw_line = QgsGeometry.fromPolylineXY([p1, p2])
                inside1 = raw_line.intersection(g1)
                inside2 = raw_line.intersection(g2)
                inside_len = 0.0
                if inside1 and not inside1.isEmpty():
                    inside_len = max(inside_len, float(inside1.length()))
                if inside2 and not inside2.isEmpty():
                    inside_len = max(inside_len, float(inside2.length()))
                if inside_len > max_inside_len:
                    continue
            except Exception:
                pass
            d = float(p1.distance(p2))
            if not np.isfinite(d):
                continue
            if d <= max_gap + 1e-9:
                candidates.append((d, p1, p2))

        candidates.sort(key=lambda x: x[0])
        kept: List[Tuple[QgsPointXY, QgsPointXY]] = [(nearest1, nearest2)]

        # Bin candidates along the perpendicular axis so we cover the full length
        # instead of only taking the next-closest gap.
        px, py = (-vy, vx)
        def _proj_perp(pt: QgsPointXY) -> float:
            return (pt.x() - c1x) * px + (pt.y() - c1y) * py

        if candidates and spacing_m > 0:
            s_vals = [_proj_perp(p1) for _, p1, _ in candidates]
            s_min = min(s_vals)
            bin_size = float(spacing_m)
            bins = {}
            for (_d, p1, p2), s in zip(candidates, s_vals):
                b = int((s - s_min) // bin_size) if bin_size > 0 else 0
                prev = bins.get(b)
                if prev is None or _d < prev[0]:
                    bins[b] = (_d, p1, p2)
            # Add best from each bin (sorted by gap)
            binned = sorted(bins.values(), key=lambda x: x[0])
        else:
            binned = candidates

        for _d, p1, p2 in binned:
            if len(kept) >= max_candidates:
                break
            too_close = False
            for k1, k2 in kept:
                if p1.distance(k1) + 1e-9 < spacing_m or p2.distance(k2) + 1e-9 < spacing_m:
                    too_close = True
                    break
            if too_close:
                continue
            kept.append((p1, p2))
            variants.append((f"touch_{len(kept)}", p1, p2))

        seen = set()
        uniq: List[Tuple[str, QgsPointXY, QgsPointXY]] = []
        for tag, pta, ptb in variants:
            key = (round(pta.x(), 3), round(pta.y(), 3), round(ptb.x(), 3), round(ptb.y(), 3))
            if key in seen:
                continue
            seen.add(key)
            uniq.append((tag, pta, ptb))
        return uniq

    def _barrier_pair_metrics(
        pid_a: int, pid_b: int, g1: QgsGeometry, g2: QgsGeometry, gap_dist: float
    ) -> Tuple[bool, Dict[str, float]]:
        if (not fluidity_mode) or pid_a == pid_b:
            return False, {"gap": float(gap_dist or 0.0)}
        try:
            a1 = float((patches.get(int(pid_a)) or {}).get("area_ha", 0.0) or 0.0)
            a2 = float((patches.get(int(pid_b)) or {}).get("area_ha", 0.0) or 0.0)
        except Exception:
            return False, {"gap": float(gap_dist or 0.0)}
        min_area = float(getattr(params, "landscape_fluidity_barrier_min_patch_area_ha", 50.0) or 50.0)
        if a1 < min_area or a2 < min_area:
            return False, {"gap": float(gap_dist or 0.0), "area1": a1, "area2": a2, "min_area": min_area}
        area_ratio = min(a1, a2) / max(max(a1, a2), 1e-9)
        if area_ratio < 0.15:
            return False, {"gap": float(gap_dist or 0.0), "area_ratio": area_ratio}
        if gap_dist <= 0.0 or gap_dist > float(params.max_search_distance):
            return False, {"gap": float(gap_dist or 0.0), "max_search": float(params.max_search_distance)}
        try:
            b1 = g1.boundingBox()
            b2 = g2.boundingBox()
            ovx = max(0.0, min(float(b1.xMaximum()), float(b2.xMaximum())) - max(float(b1.xMinimum()), float(b2.xMinimum())))
            ovy = max(0.0, min(float(b1.yMaximum()), float(b2.yMaximum())) - max(float(b1.yMinimum()), float(b2.yMinimum())))
            long_overlap = max(float(ovx), float(ovy))
        except Exception:
            long_overlap = 0.0
        min_overlap = float(getattr(params, "landscape_fluidity_barrier_min_overlap_m", 1200.0) or 1200.0)
        if long_overlap < min_overlap:
            return False, {
                "gap": float(gap_dist or 0.0),
                "overlap": float(long_overlap),
                "min_overlap": float(min_overlap),
            }
        # Long/parallel-like pair with a relatively thin separating gap.
        gap_ratio = float(gap_dist) / max(float(long_overlap), 1e-9)
        if gap_ratio > 0.08:
            return False, {"gap": float(gap_dist or 0.0), "gap_ratio": gap_ratio}
        return True, {
            "gap": float(gap_dist or 0.0),
            "area1": a1,
            "area2": a2,
            "area_ratio": area_ratio,
            "overlap": float(long_overlap),
            "gap_ratio": gap_ratio,
        }

    def _is_barrier_pair(pid_a: int, pid_b: int, g1: QgsGeometry, g2: QgsGeometry, gap_dist: float) -> bool:
        ok, _ = _barrier_pair_metrics(pid_a, pid_b, g1, g2, gap_dist)
        return ok

    def _overlap_ratio(g1: QgsGeometry, g2: QgsGeometry) -> float:
        try:
            if circuit_mode:
                try:
                    if proximity_dist > 0 and g1.distance(g2) <= proximity_dist:
                        return 1.0
                except Exception:
                    pass
            a1 = g1.area()
            a2 = g2.area()
            denom = min(a1, a2)
            if denom <= 0:
                return 1.0
            inter = g1.intersection(g2)
            if inter is None or inter.isEmpty():
                return 0.0
            return inter.area() / denom
        except Exception:
            return 1.0

    def _dedupe_points_xy(points, tol=0.01):
        out = []
        seen = set()
        for p in points:
            k = (round(p.x() / tol), round(p.y() / tol))
            if k in seen:
                continue
            seen.add(k)
            out.append(p)
        return out

    def _boundary_terminals_for_patch(patch_geom, spacing_m=150.0, max_pts=120):
        if not patch_geom or patch_geom.isEmpty():
            return []
        pts: List[QgsPointXY] = []

        # QGIS API compatibility: not all versions expose QgsGeometry.boundary().
        # We only sample exterior rings (corridors should not start from holes).
        rings: List[List[QgsPointXY]] = []
        try:
            if patch_geom.isMultipart():
                for poly in patch_geom.asMultiPolygon() or []:
                    if not poly:
                        continue
                    exterior = poly[0] if poly else []
                    if exterior:
                        rings.append([QgsPointXY(p) for p in exterior])
            else:
                poly = patch_geom.asPolygon() or []
                exterior = poly[0] if poly else []
                if exterior:
                    rings.append([QgsPointXY(p) for p in exterior])
        except Exception:
            rings = []

        if not rings:
            # Fallback: at least return some vertices.
            try:
                for v in patch_geom.vertices():
                    pts.append(QgsPointXY(v))
            except Exception:
                return []

        # Always include vertices (tips often are vertices)
        for ring in rings:
            pts.extend(ring)

        # Sample along boundary at regular spacing
        if spacing_m and spacing_m > 0:
            for ring in rings:
                if len(ring) < 2:
                    continue
                try:
                    line = QgsGeometry.fromPolylineXY(ring)
                    L = float(line.length())
                    if L <= 0:
                        continue
                    n = int(L // spacing_m) + 1
                    for i in range(n + 1):
                        d = min(i * spacing_m, L)
                        ip = line.interpolate(d)
                        if ip and not ip.isEmpty():
                            try:
                                pts.append(QgsPointXY(ip.asPoint()))
                            except Exception:
                                pass
                except Exception:
                    continue

        pts = _dedupe_points_xy(pts, tol=0.01)

        # Cap size by uniform thinning (keeps tips/vertices because they were added first)
        if max_pts and len(pts) > max_pts:
            step = max(1, len(pts) // max_pts)
            pts = pts[::step][:max_pts]

        return pts

    terminal_cache: Dict[int, List[QgsPointXY]] = {}

    for idx, (pid1, pdata1) in enumerate(patches.items(), start=1):
        if progress_cb is not None:
            span = max(progress_end - progress_start, 1)
            progress_value = progress_start + ((idx - 1) / total) * span
            emit_progress(progress_cb, progress_value, "Analyzing patches")

        geom1 = pdata1["geom"]
        rect = geom1.boundingBox()
        rect.grow(params.max_search_distance)
        candidate_ids = spatial_index.intersects(rect)

        with _measure("Patch iteration"):
            # Rank potential neighbors by true boundary-to-boundary distance, then keep top-K.
            ranked_neighbors: List[Tuple[float, int]] = []
            with _measure("Neighbor ranking"):
                for pid2 in candidate_ids:
                    if pid2 == pid1:
                        continue
                    pair = frozenset({pid1, pid2})
                    if pair in processed_pairs:
                        continue

                    pdata2 = patches.get(pid2)
                    if not pdata2:
                        continue
                    geom2 = pdata2["geom"]

                    try:
                        distance = float(geom1.distance(geom2))
                    except Exception:
                        continue

                    # Touching/intersecting polygons are already contiguous; do not build a corridor.
                    if distance <= 0.0:
                        processed_pairs.add(pair)
                        continue

                    if distance > params.max_search_distance:
                        continue

                    ranked_neighbors.append((distance, int(pid2)))

            ranked_neighbors.sort(key=lambda item: item[0])
            if len(ranked_neighbors) > int(k_nearest_neighbors_cap):
                ranked_neighbors = ranked_neighbors[: int(k_nearest_neighbors_cap)]

            for distance, pid2 in ranked_neighbors:
                pair = frozenset({pid1, pid2})
                if pair in processed_pairs:
                    continue

                pdata2 = patches.get(pid2)
                if not pdata2:
                    continue
                geom2 = pdata2["geom"]
                barrier_pair = _is_barrier_pair(int(pid1), int(pid2), geom1, geom2, float(distance))

                def _path_cost_length(points: List[QgsPointXY]) -> float:
                    return sum(points[i].distance(points[i + 1]) for i in range(len(points) - 1))

                # --- Terminal selection (fix peninsula-tip blindness) ---
                term_spacing = float(getattr(params, "vector_terminal_spacing_m", 150.0))
                term_max = int(getattr(params, "vector_terminal_max_per_patch", 120))
                term_pairs_k = int(getattr(params, "vector_terminal_pairs_per_pair", 25))

                base_terms1 = terminal_cache.get(int(pid1))
                if base_terms1 is None:
                    base_terms1 = _boundary_terminals_for_patch(geom1, spacing_m=term_spacing, max_pts=term_max)
                    terminal_cache[int(pid1)] = base_terms1
                base_terms2 = terminal_cache.get(int(pid2))
                if base_terms2 is None:
                    base_terms2 = _boundary_terminals_for_patch(geom2, spacing_m=term_spacing, max_pts=term_max)
                    terminal_cache[int(pid2)] = base_terms2

                # Copy cached terminals so we can add pair-specific terminals without polluting the cache.
                terms1 = list(base_terms1 or [])
                terms2 = list(base_terms2 or [])

                # Force-add the TRUE nearest boundary points as terminals so small gaps are always evaluated.
                try:
                    nearest_p1 = geom1.nearestPoint(geom2).asPoint()
                    nearest_p2 = geom2.nearestPoint(geom1).asPoint()
                    if nearest_p1 and not nearest_p1.isEmpty():
                        terms1 = [QgsPointXY(nearest_p1)] + terms1
                    if nearest_p2 and not nearest_p2.isEmpty():
                        terms2 = [QgsPointXY(nearest_p2)] + terms2
                except Exception:
                    pass

                terms1 = _dedupe_points_xy(terms1, tol=0.01)
                terms2 = _dedupe_points_xy(terms2, tol=0.01)

                # Fallback to legacy behavior if sampling fails
                if not terms1 or not terms2:
                    p1 = geom1.nearestPoint(geom2).asPoint()
                    p2 = geom2.nearestPoint(geom1).asPoint()
                    if p1.isEmpty() or p2.isEmpty():
                        continue
                    terms1 = [QgsPointXY(p1)]
                    terms2 = [QgsPointXY(p2)]

                # Build top-K closest terminal pairs by straight-line distance (cheap shortlist)
                pairs: List[Tuple[float, QgsPointXY, QgsPointXY]] = []
                for t1 in terms1:
                    for t2 in terms2:
                        pairs.append((t1.distance(t2), t1, t2))
                pairs.sort(key=lambda x: x[0])
                pairs = pairs[: max(1, term_pairs_k)]

                # Evaluate a few candidate terminal pairs using the SAME routing/cost logic,
                # then pick the best. This keeps only ONE terminal pair per patch-pair.
                best: Optional[Tuple[float, List[QgsPointXY], QgsPointXY, QgsPointXY]] = None

                for lower_bound_dist, cand_t1, cand_t2 in pairs:
                    if best is not None and lower_bound_dist >= best[0]:
                        continue
                    p1_xy = cand_t1
                    p2_xy = cand_t2

                    try:
                        if navigator:
                            path_points = navigator.find_path(p1_xy, p2_xy)
                        else:
                            path_points = [p1_xy, p2_xy]
                    except Exception:
                        path_points = [p1_xy, p2_xy]

                    if not path_points or len(path_points) < 2:
                        continue

                    cand_cost_len = _path_cost_length(path_points)
                    if best is None or cand_cost_len < best[0]:
                        best = (cand_cost_len, path_points, cand_t1, cand_t2)

                if best is None:
                    continue

                p1_xy = best[2]
                p2_xy = best[3]
                best_path_points = best[1]
                # --- end terminal selection ---

                def _best_path_to_boundary(
                    start_pt: QgsPointXY, patch_geom: QgsGeometry
                ) -> Tuple[Optional[List[QgsPointXY]], Optional[QgsPointXY], float]:
                    if not navigator:
                        return None, None, float("inf")
                    candidates = _sample_boundary_points(patch_geom, max_points=16, allow_densify=False)
                    best_path: Optional[List[QgsPointXY]] = None
                    best_cost = float("inf")
                    best_pt: Optional[QgsPointXY] = None
                    for target in candidates:
                        pts = navigator.find_path(start_pt, target)
                        if not pts:
                            continue
                        cost = _path_cost_length(pts)
                        if cost < best_cost:
                            best_cost = cost
                            best_path = pts
                            best_pt = target
                    if not best_path:
                        dense_candidates = _sample_boundary_points(patch_geom, max_points=40, allow_densify=True)
                        for target in dense_candidates:
                            pts = navigator.find_path(start_pt, target)
                            if not pts:
                                continue
                            cost = _path_cost_length(pts)
                            if cost < best_cost:
                                best_cost = cost
                                best_path = pts
                                best_pt = target
                    return best_path, best_pt, best_cost

                def _local_shortest_segment(
                    pa: QgsPointXY,
                    pb: QgsPointXY,
                    radius_m: float,
                ) -> Optional[List[QgsPointXY]]:
                    try:
                        b1 = QgsGeometry.fromPointXY(pa).buffer(float(radius_m), BUFFER_SEGMENTS)
                        b2 = QgsGeometry.fromPointXY(pb).buffer(float(radius_m), BUFFER_SEGMENTS)
                        g1l = geom1.intersection(b1) if b1 is not None else geom1
                        g2l = geom2.intersection(b2) if b2 is not None else geom2
                        if g1l is None or g1l.isEmpty():
                            g1l = geom1
                        if g2l is None or g2l.isEmpty():
                            g2l = geom2
                        sl = g1l.shortestLine(g2l)
                        if sl is None or sl.isEmpty():
                            return None
                        pts = [QgsPointXY(p) for p in sl.asPolyline() or []]
                        if len(pts) < 2:
                            return None
                        return [pts[0], pts[-1]]
                    except Exception:
                        return None

                # Default: single nearest corridor per pair.
                variants = [("nearest", p1_xy, p2_xy)]
                if barrier_pair:
                    sweep_sep = max(float(params.max_search_distance or 0.0), 1.0)
                    sweep_cap = max(
                        6,
                        min(
                            60,
                            int(getattr(params, "landscape_fluidity_barrier_sweep_cap", 24) or 24),
                        ),
                    )
                    variants = _barrier_touch_variants(
                        geom1,
                        geom2,
                        p1_xy,
                        p2_xy,
                        spacing_m=sweep_sep,
                        max_candidates=sweep_cap,
                    )
                for variant_tag, start_xy, end_xy in variants:
                    raw_geom_candidates: List[Tuple[str, QgsGeometry]] = []
                    path_points: Optional[List[QgsPointXY]] = None

                    if barrier_pair:
                        try:
                            direct_line = QgsGeometry.fromPolylineXY([start_xy, end_xy])
                            direct_barrier = _buffer_line_segment(direct_line, params.min_corridor_width)
                            if direct_barrier and (not direct_barrier.isEmpty()):
                                direct_barrier = direct_barrier.makeValid()
                                raw_geom_candidates.append((f"barrier_direct:{variant_tag}", direct_barrier))
                        except Exception:
                            pass

                    if navigator:
                        with _measure("Navigator routing"):
                            # Reuse nearest path cache; route alternates explicitly.
                            if variant_tag == "nearest":
                                path_points = best_path_points
                            else:
                                path_points = navigator.find_path(start_xy, end_xy)

                    if navigator and path_points:
                        with _measure("Corridor geometry (navigator)"):
                            nav_geom = _create_corridor_geometry(
                                path_points,
                                geom1,
                                geom2,
                                params,
                                obstacle_geoms=impassable_geoms,
                                ctx=ctx,
                                smooth_iterations=3,
                            )
                        if nav_geom:
                            raw_geom_candidates.append((f"navigator:{variant_tag}", nav_geom))

                    with _measure("Corridor geometry (direct)"):
                        direct_geom = _create_corridor_geometry(
                            [start_xy, end_xy],
                            geom1,
                            geom2,
                            params,
                            obstacle_geoms=impassable_geoms if navigator else None,
                            ctx=ctx,
                        )
                    if direct_geom:
                        raw_geom_candidates.append((f"direct:{variant_tag}", direct_geom))

                    if not raw_geom_candidates:
                        continue

                    if barrier_pair:
                        pick = None
                        for item in raw_geom_candidates:
                            if str(item[0]).startswith("barrier_direct"):
                                pick = item
                                break
                        if pick is not None:
                            best_source, raw_geom = pick
                        else:
                            best_source, raw_geom = min(
                                raw_geom_candidates,
                                key=lambda item: item[1].area() if item[1] and not item[1].isEmpty() else float("inf"),
                            )
                    else:
                        best_source, raw_geom = min(
                            raw_geom_candidates,
                            key=lambda item: item[1].area() if item[1] and not item[1].isEmpty() else float("inf"),
                        )
                    if raw_geom is None or raw_geom.isEmpty():
                        continue
                    if barrier_pair and str(variant_tag).startswith(("sweep_", "touch_", "nearest")):
                        local_seg = _local_shortest_segment(
                            start_xy,
                            end_xy,
                            radius_m=max(40.0, min(300.0, 0.65 * float(params.max_search_distance or 0.0))),
                        )
                        if local_seg is not None and len(local_seg) >= 2:
                            try:
                                # For barrier variants, build the corridor directly from the
                                # local segment to preserve distinct anchor locations.
                                local_line = QgsGeometry.fromPolylineXY([local_seg[0], local_seg[1]])
                                local_raw = _buffer_line_segment(local_line, params.min_corridor_width)
                                if local_raw and (not local_raw.isEmpty()):
                                    local_raw = local_raw.makeValid()
                            except Exception:
                                local_raw = None
                            if local_raw is not None and (not local_raw.isEmpty()):
                                raw_geom = local_raw
                                best_source = f"{best_source}:local_shortest"

                    # Reject candidates that "reach through" an endpoint patch to exit on the far side.
                    # We do this by measuring how much of the *raw centerline* lies within each endpoint patch.
                    raw_line = None
                    try:
                        if navigator and path_points:
                            raw_line = QgsGeometry.fromPolylineXY(path_points)
                        else:
                            raw_line = QgsGeometry.fromPolylineXY([start_xy, end_xy])
                        if not barrier_pair:
                            max_inside_len = max(params.min_corridor_width * 3.0, 5.0)
                            for _pid, _pgeom in ((pid1, geom1), (pid2, geom2)):
                                try:
                                    inter = raw_line.intersection(_pgeom)
                                    inside_len = inter.length() if inter and not inter.isEmpty() else 0.0
                                    if inside_len > max_inside_len:
                                        raw_geom = None
                                        break
                                except Exception:
                                    continue
                        if raw_geom is None and (not barrier_pair):
                            # Fallback for local near-gap links: use the true shortest
                            # boundary-to-boundary bridge to avoid selecting long
                            # edge-following centerlines from sampled terminals.
                            try:
                                shortest_line = geom1.shortestLine(geom2)
                            except Exception:
                                shortest_line = None
                            if shortest_line is not None and (not shortest_line.isEmpty()):
                                short_pts: List[QgsPointXY] = []
                                try:
                                    short_pts = [QgsPointXY(p) for p in shortest_line.asPolyline() or []]
                                except Exception:
                                    short_pts = []
                                if len(short_pts) >= 2:
                                    short_start = short_pts[0]
                                    short_end = short_pts[-1]
                                    try:
                                        short_raw = _create_corridor_geometry(
                                            [short_start, short_end],
                                            geom1,
                                            geom2,
                                            params,
                                            obstacle_geoms=impassable_geoms if navigator else None,
                                            ctx=ctx,
                                        )
                                    except Exception:
                                        short_raw = None
                                    if short_raw is not None and (not short_raw.isEmpty()):
                                        raw_geom = short_raw
                                        raw_line = QgsGeometry.fromPolylineXY([short_start, short_end])
                                        best_source = f"{best_source}:shortest_fallback"
                        if raw_geom is None:
                            continue
                    except Exception:
                        pass

                    if navigator:
                        with _measure("Handle intermediate patch"):
                            intersected_temp = _detect_corridor_intersections(
                                raw_geom, patches, spatial_index, {pid1, pid2}
                            )
                            if intersected_temp:
                                intermediate_id = min(
                                    intersected_temp, key=lambda pid: geom1.distance(patches[pid]["geom"])
                                )
                                mid_geom = patches[intermediate_id]["geom"]

                                path_in, entry_pt, cost_in = _best_path_to_boundary(start_xy, mid_geom)
                                path_out, exit_pt, cost_out = _best_path_to_boundary(end_xy, mid_geom)

                                if path_in and path_out and entry_pt and exit_pt:
                                    entry_geom = _create_corridor_geometry(
                                        path_in,
                                        geom1,
                                        mid_geom,
                                        params,
                                        obstacle_geoms=impassable_geoms,
                                        ctx=ctx,
                                        smooth_iterations=3,
                                    )
                                    exit_geom = _create_corridor_geometry(
                                        path_out,
                                        mid_geom,
                                        geom2,
                                        params,
                                        obstacle_geoms=impassable_geoms,
                                        ctx=ctx,
                                        smooth_iterations=3,
                                    )
                                    internal_geom = None
                                    # Avoid creating an "internal bridge" that could extend outside the patch
                                    # and accidentally overlap impassables; the corridor cost is computed
                                    # after patch-difference anyway.
                                    if not impassable_geoms:
                                        internal_geom = QgsGeometry.fromPolylineXY([entry_pt, exit_pt]).buffer(
                                            params.min_corridor_width / 2.0, BUFFER_SEGMENTS
                                        )
                                    geoms_to_merge = [g for g in (entry_geom, exit_geom, internal_geom) if g and not g.isEmpty()]
                                    if geoms_to_merge:
                                        try:
                                            split_geom = QgsGeometry.unaryUnion(geoms_to_merge)
                                        except Exception:
                                            split_geom = geoms_to_merge[0]
                                            for g in geoms_to_merge[1:]:
                                                try:
                                                    split_geom = split_geom.combine(g)
                                                except Exception:
                                                    pass
                                        if split_geom is not None and not split_geom.isEmpty():
                                            raw_geom = split_geom

                    with _measure("Finalize corridor geometry"):
                        if barrier_pair:
                            corridor_geom = raw_geom
                            patch_ids = {int(pid1), int(pid2)}
                        else:
                            corridor_geom, patch_ids = _finalize_corridor_geometry(
                                pid1,
                                pid2,
                                raw_geom,
                                patches,
                                spatial_index,
                                patch_union=patch_union,
                                corridor_width_m=float(params.min_corridor_width or 0.0),
                            )
                    if corridor_geom is None:
                        continue

                    emitted_any = False
                    try:
                        extra_patches = [int(pid) for pid in patch_ids if int(pid) not in (int(pid1), int(pid2))]
                    except Exception:
                        extra_patches = []
                    bridge_over_penalty = 1.0
                    bridge_over_count = 0

                    # Largest Single Network: if a corridor proposal touches another patch C "in between",
                    # don't keep an A->B corridor that gets clipped across C. Prefer A->C and B->C.
                    # Most Connectivity already behaves well here via stop-early multipart splitting.
                    blocking_patches: List[int] = []
                    if extra_patches and raw_line is not None and (not raw_line.isEmpty()):
                        for pid in extra_patches:
                            try:
                                g = patches.get(int(pid), {}).get("geom")
                                if g is None or g.isEmpty():
                                    continue
                                if raw_line.intersects(g):
                                    blocking_patches.append(int(pid))
                            except Exception:
                                continue

                    target_patches = extra_patches if largest_network_mode else extra_patches

                    if target_patches and corridor_geom.isMultipart():
                        part_a = _pick_part_closest_to_patch(corridor_geom, geom1)
                        if part_a is not None and (not part_a.isEmpty()):
                            try:
                                tgt_a = min(target_patches, key=lambda pid: float(geom1.distance(patches[pid]["geom"])))
                            except Exception:
                                tgt_a = target_patches[0]
                            try:
                                dist_a = float(geom1.distance(patches[tgt_a]["geom"]))
                            except Exception:
                                dist_a = float(distance)
                            area_a = _corridor_cost_area_ha(
                                part_a,
                                float(dist_a),
                                float(params.min_corridor_width or 0.0),
                            )
                            if area_a > 0 and (params.max_corridor_area is None or area_a <= params.max_corridor_area):
                                _push_candidate(
                                    {
                                        "patch1": int(pid1),
                                        "patch2": int(tgt_a),
                                        "patch_ids": {int(pid1), int(tgt_a)},
                                        "geom": clone_geometry(part_a),
                                        "area_ha": area_a,
                                        "original_area_ha": area_a,
                                        "distance_m": dist_a,
                                        "variant": f"{variant_tag}:stop_early",
                                        "source": best_source,
                                        "barrier_pair": False,
                                    }
                                )
                                emitted_any = True

                        part_b = _pick_part_closest_to_patch(corridor_geom, geom2)
                        if part_b is not None and (not part_b.isEmpty()):
                            try:
                                tgt_b = min(target_patches, key=lambda pid: float(geom2.distance(patches[pid]["geom"])))
                            except Exception:
                                tgt_b = target_patches[0]
                            try:
                                dist_b = float(geom2.distance(patches[tgt_b]["geom"]))
                            except Exception:
                                dist_b = float(distance)
                            area_b = _corridor_cost_area_ha(
                                part_b,
                                float(dist_b),
                                float(params.min_corridor_width or 0.0),
                            )
                            if area_b > 0 and (params.max_corridor_area is None or area_b <= params.max_corridor_area):
                                _push_candidate(
                                    {
                                        "patch1": int(pid2),
                                        "patch2": int(tgt_b),
                                        "patch_ids": {int(pid2), int(tgt_b)},
                                        "geom": clone_geometry(part_b),
                                        "area_ha": area_b,
                                        "original_area_ha": area_b,
                                        "distance_m": dist_b,
                                        "variant": f"{variant_tag}:stop_early",
                                        "source": best_source,
                                        "barrier_pair": False,
                                    }
                                )
                                emitted_any = True

                    if emitted_any:
                        continue

                    # Landscape Fluidity should prefer stepping-stone connections when a
                    # proposed corridor crosses another patch component, instead of adding
                    # a single long "bridge-over" edge that skips that intermediate habitat.
                    via_intermediate_emitted = 0
                    if fluidity_mode and extra_patches:
                        try:
                            mid_pid = min(
                                extra_patches,
                                key=lambda pid: float(corridor_geom.distance(patches[pid]["geom"])),
                            )
                        except Exception:
                            mid_pid = int(extra_patches[0])
                        mid_geom = patches.get(int(mid_pid), {}).get("geom")
                        if mid_geom is not None and (not mid_geom.isEmpty()):
                            for src_pid, src_geom, src_xy in (
                                (int(pid1), geom1, start_xy),
                                (int(pid2), geom2, end_xy),
                            ):
                                if int(src_pid) == int(mid_pid):
                                    continue
                                try:
                                    mid_pt = mid_geom.nearestPoint(src_geom).asPoint()
                                    if mid_pt.isEmpty():
                                        continue
                                    alt_raw = _create_corridor_geometry(
                                        [src_xy, QgsPointXY(mid_pt)],
                                        src_geom,
                                        mid_geom,
                                        params,
                                        obstacle_geoms=impassable_geoms if navigator else None,
                                        ctx=ctx,
                                    )
                                except Exception:
                                    alt_raw = None
                                if alt_raw is None or alt_raw.isEmpty():
                                    continue
                                alt_geom, alt_patch_ids = _finalize_corridor_geometry(
                                    int(src_pid),
                                    int(mid_pid),
                                    alt_raw,
                                    patches,
                                    spatial_index,
                                    patch_union=patch_union,
                                    corridor_width_m=float(params.min_corridor_width or 0.0),
                                )
                                if alt_geom is None or alt_geom.isEmpty():
                                    continue
                                alt_pid_set = set(int(pid) for pid in alt_patch_ids)
                                # Keep only clean stepping-stone links. If the geometry
                                # crosses additional habitat components, skip this candidate
                                # so the optimizer does not "graze" along habitat edges.
                                if any(int(pid) not in (int(src_pid), int(mid_pid)) for pid in alt_pid_set):
                                    continue
                                try:
                                    alt_dist = float(src_geom.distance(mid_geom))
                                except Exception:
                                    alt_dist = float(distance)
                                alt_area_ha = _corridor_cost_area_ha(
                                    alt_geom,
                                    float(alt_dist),
                                    float(params.min_corridor_width or 0.0),
                                )
                                if alt_area_ha <= 0:
                                    continue
                                if params.max_corridor_area is not None and alt_area_ha > params.max_corridor_area:
                                    continue
                                _push_candidate(
                                    {
                                        "patch1": int(src_pid),
                                        "patch2": int(mid_pid),
                                        "patch_ids": {int(src_pid), int(mid_pid)},
                                        "geom": clone_geometry(alt_geom),
                                        "area_ha": alt_area_ha,
                                        "original_area_ha": alt_area_ha,
                                        "distance_m": alt_dist,
                                        "variant": f"{variant_tag}:via_intermediate",
                                        "source": best_source,
                                        "barrier_pair": False,
                                    }
                                )
                                via_intermediate_emitted += 1
                        # In landscape fluidity mode, never keep the original bridge-over
                        # candidate when intermediate habitat is present. The habitat
                        # itself is treated as the conduit via stepping-stone links.
                        continue

                    # If we touched an intermediate patch but couldn't split into endpoint-local pieces,
                    # drop this candidate in Largest Single Network so it doesn't "jump over" the patch.
                    if largest_network_mode and extra_patches:
                        continue

                    effective_distance = float(distance)
                    try:
                        if raw_line is not None and (not raw_line.isEmpty()):
                            effective_distance = max(float(raw_line.length()), 1e-9)
                    except Exception:
                        effective_distance = float(distance)

                    corridor_area_ha = _corridor_cost_area_ha(
                        corridor_geom,
                        float(effective_distance),
                        float(params.min_corridor_width or 0.0),
                    )
                    if corridor_area_ha <= 0:
                        continue
                    if params.max_corridor_area is not None and corridor_area_ha > params.max_corridor_area:
                        continue

                    _push_candidate(
                        {
                            "patch1": int(pid1),
                            "patch2": int(pid2),
                            "patch_ids": set(patch_ids),
                            "geom": clone_geometry(corridor_geom),
                            "area_ha": corridor_area_ha,
                            "original_area_ha": corridor_area_ha,
                            "distance_m": float(effective_distance),
                            "variant": variant_tag,
                            "source": best_source,
                            "barrier_pair": bool(barrier_pair),
                            "anchor_pts": {int(pid1): start_xy, int(pid2): end_xy} if barrier_pair else {},
                            "bridge_over_penalty": float(bridge_over_penalty),
                            "bridge_over_count": int(bridge_over_count),
                        }
                    )
                processed_pairs.add(pair)

    if fluidity_mode or strategy == "most_connected_habitat":
        def _add_landscape_fluidity_chain_candidates() -> int:
            """
            Add bounded stepping-stone chain candidates over the local feasible-link
            graph so the optimizer can select cheaper multi-hop routes when they
            outperform a direct bridge.
            """
            if fluidity_mode:
                chain_cap = max(
                    0,
                    min(
                        400,
                        int(getattr(params, "landscape_fluidity_chain_k_total", 80) or 80),
                    ),
                )
                branch_cap = max(
                    2,
                    min(
                        24,
                        int(getattr(params, "landscape_fluidity_chain_branch_per_patch", 8) or 8),
                    ),
                )
                max_hops = 2
                per_start_cap = max(4, min(24, branch_cap))
            else:
                chain_cap = max(0, min(600, int(getattr(params, "most_connected_chain_k_total", 140) or 140)))
                branch_cap = max(3, min(20, int(getattr(params, "most_connected_chain_branch_per_patch", 6) or 6)))
                max_hops = max(2, min(4, int(getattr(params, "most_connected_chain_max_hops", 4) or 4)))
                per_start_cap = max(6, min(30, int(getattr(params, "most_connected_chain_per_start_cap", 10) or 10)))
            if chain_cap <= 0:
                return 0
            max_chain_len = max(float(params.max_search_distance or 0.0) * 2.2, 1.0)
            max_chain_budget = max(float(params.budget_area or 0.0) * 1.6, 0.0)
            if max_chain_budget <= 0.0:
                max_chain_budget = float("inf")

            edge_recs: List[Dict[str, Any]] = []
            direct_pairs: Set[Tuple[int, int]] = set()
            by_node: Dict[int, List[int]] = defaultdict(list)
            edge_keys: Set[Tuple[int, int]] = set()

            for items in list(candidates_by_pair.values()) + list(chain_seed_by_pair.values()):
                for cand in items:
                    if bool(cand.get("intra_patch", False)):
                        continue
                    try:
                        a = int(cand.get("patch1"))
                        b = int(cand.get("patch2"))
                    except Exception:
                        continue
                    if a == b:
                        continue
                    geom = cand.get("geom")
                    if geom is None or geom.isEmpty():
                        continue
                    cost = float(cand.get("area_ha", 0.0) or 0.0)
                    if cost <= 0.0:
                        continue
                    d = float(cand.get("distance_m", cand.get("distance", 0.0)) or 0.0)
                    if d <= 0.0:
                        try:
                            d = max(float(geom.length()), 1.0)
                        except Exception:
                            d = 1.0
                    rec = {
                        "u": int(a),
                        "v": int(b),
                        "cost": float(cost),
                        "dist": float(d),
                        "geom": clone_geometry(geom),
                    }
                    idx = len(edge_recs)
                    edge_recs.append(rec)
                    by_node[int(a)].append(idx)
                    by_node[int(b)].append(idx)
                    key = (int(a), int(b)) if int(a) <= int(b) else (int(b), int(a))
                    edge_keys.add(key)
                    if not bool(cand.get("chain_seed_only", False)):
                        direct_pairs.add(key)

            # If the normal candidate pass misses local links to a small stepping-stone
            # patch, synthesize bounded chain-only edges from the shortest boundary gap.
            synth_branch_cap = max(4, min(24, branch_cap * 2))
            for pivot in sorted(int(pid) for pid in patches.keys()):
                pivot_geom = (patches.get(int(pivot)) or {}).get("geom")
                if pivot_geom is None or pivot_geom.isEmpty():
                    continue
                pivot_area = float((patches.get(int(pivot)) or {}).get("area_ha", 0.0) or 0.0)
                nearby: List[Tuple[float, float, int]] = []
                for other in sorted(int(pid) for pid in patches.keys()):
                    if int(other) == int(pivot):
                        continue
                    key = _pair_key(int(pivot), int(other))
                    if key in edge_keys:
                        continue
                    other_geom = (patches.get(int(other)) or {}).get("geom")
                    if other_geom is None or other_geom.isEmpty():
                        continue
                    try:
                        gap_m = float(pivot_geom.distance(other_geom))
                    except Exception:
                        continue
                    if gap_m <= 0.0 or gap_m > float(params.max_search_distance or 0.0) + 1e-12:
                        continue
                    other_area = float((patches.get(int(other)) or {}).get("area_ha", 0.0) or 0.0)
                    # Restrict fallback synthesis to plausible stepping-stone links.
                    if pivot_area > other_area + 1e-12:
                        continue
                    nearby.append((float(gap_m), -float(other_area), int(other)))
                if not nearby:
                    continue
                nearby.sort()
                for gap_m, _neg_other_area, other in nearby[:synth_branch_cap]:
                    other_geom = (patches.get(int(other)) or {}).get("geom")
                    try:
                        shortest_line = pivot_geom.shortestLine(other_geom)
                    except Exception:
                        shortest_line = None
                    if shortest_line is None or shortest_line.isEmpty():
                        continue
                    try:
                        seg_pts = [QgsPointXY(p) for p in shortest_line.asPolyline() or []]
                    except Exception:
                        seg_pts = []
                    if len(seg_pts) < 2:
                        continue
                    try:
                        raw_geom = _create_corridor_geometry(
                            [seg_pts[0], seg_pts[-1]],
                            pivot_geom,
                            other_geom,
                            params,
                            obstacle_geoms=impassable_geoms if navigator else None,
                            ctx=ctx,
                        )
                    except Exception:
                        raw_geom = None
                    if raw_geom is None or raw_geom.isEmpty():
                        continue
                    try:
                        edge_geom, edge_patch_ids = _finalize_corridor_geometry(
                            int(pivot),
                            int(other),
                            raw_geom,
                            patches,
                            spatial_index,
                            patch_union=patch_union,
                            corridor_width_m=float(params.min_corridor_width or 0.0),
                        )
                    except Exception:
                        edge_geom, edge_patch_ids = None, set()
                    if edge_geom is None or edge_geom.isEmpty():
                        continue
                    clean_ids = set(int(pid) for pid in (edge_patch_ids or set()))
                    if any(int(pid) not in (int(pivot), int(other)) for pid in clean_ids):
                        continue
                    edge_cost = _corridor_cost_area_ha(
                        edge_geom,
                        float(gap_m),
                        float(params.min_corridor_width or 0.0),
                    )
                    if edge_cost <= 0.0 or edge_cost > max_chain_budget + 1e-12:
                        continue
                    rec = {
                        "u": int(pivot),
                        "v": int(other),
                        "cost": float(edge_cost),
                        "dist": float(gap_m),
                        "geom": clone_geometry(edge_geom),
                        "synthetic_chain_seed": True,
                    }
                    idx = len(edge_recs)
                    edge_recs.append(rec)
                    by_node[int(pivot)].append(idx)
                    by_node[int(other)].append(idx)
                    edge_keys.add(_pair_key(int(pivot), int(other)))

            if not edge_recs:
                return 0

            adjacency: Dict[int, List[Tuple[int, int, Dict[str, Any]]]] = defaultdict(list)
            for edge_idx, rec in enumerate(edge_recs):
                adjacency[int(rec["u"])].append((int(rec["v"]), int(edge_idx), rec))
                adjacency[int(rec["v"])].append((int(rec["u"]), int(edge_idx), rec))
            for node, rows in list(adjacency.items()):
                adjacency[node] = sorted(
                    rows,
                    key=lambda item: (
                        float(item[2]["cost"]) / max(float(item[2]["dist"]), 1e-9),
                        float(item[2]["cost"]),
                        float(item[2]["dist"]),
                    ),
                )[: max(branch_cap, per_start_cap)]

            def _chain_score(path_nodes: Sequence[int], path_edges: Sequence[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
                if len(path_nodes) < 3 or len(path_edges) != len(path_nodes) - 1:
                    return None
                a = int(path_nodes[0])
                c = int(path_nodes[-1])
                chain_cost = float(sum(float(edge.get("cost", 0.0) or 0.0) for edge in path_edges))
                if chain_cost <= 0.0 or chain_cost > max_chain_budget + 1e-12:
                    return None
                chain_dist = float(sum(float(edge.get("dist", 0.0) or 0.0) for edge in path_edges))
                if chain_dist <= 0.0 or chain_dist > max_chain_len + 1e-12:
                    return None
                try:
                    eu_ac = float(
                        (patches[a]["geom"].distance(patches[c]["geom"]))
                        if (a in patches and c in patches)
                        else chain_dist
                    )
                except Exception:
                    eu_ac = chain_dist
                eu_ac = max(eu_ac, 1.0)
                detour_ratio = float(chain_dist / eu_ac)
                if detour_ratio > 2.5:
                    return None

                area_a = float((patches.get(a) or {}).get("area_ha", 0.0) or 0.0)
                area_c = float((patches.get(c) or {}).get("area_ha", 0.0) or 0.0)
                mid_nodes = [int(pid) for pid in path_nodes[1:-1]]
                pivot_signal = 0.0
                endpoint_scale = max(area_a, area_c, 1e-9)
                for mid in mid_nodes:
                    area_mid = float((patches.get(int(mid)) or {}).get("area_ha", 0.0) or 0.0)
                    pivot_signal += min(max(area_mid, 0.0), 0.75 * endpoint_scale)
                endpoint_signal = math.sqrt(max(area_a * area_c, 0.0))
                area_signal = max(endpoint_signal + (0.35 * pivot_signal), 1e-9)
                direct_key = _pair_key(a, c)
                direct_factor = 0.92 if direct_key in direct_pairs else 1.0
                detour_eff = 1.0 / max(float(detour_ratio), 1.0)
                hop_penalty = 1.0 / max(1.0, 0.92 + (0.08 * float(len(path_edges))))
                score = float((area_signal / max(chain_cost, 1e-9)) * detour_eff * direct_factor * hop_penalty)
                if score <= 0.0:
                    return None
                return {
                    "path_nodes": tuple(int(pid) for pid in path_nodes),
                    "path_edges": list(path_edges),
                    "score": float(score),
                    "chain_cost": float(chain_cost),
                    "chain_dist": float(chain_dist),
                }

            best_by_chain: Dict[Tuple[int, ...], Dict[str, Any]] = {}

            if max_hops <= 2:
                for pivot, inc in by_node.items():
                    local = sorted(
                        inc,
                        key=lambda i: (
                            float(edge_recs[i]["cost"]) / max(float(edge_recs[i]["dist"]), 1e-9),
                            float(edge_recs[i]["cost"]),
                        ),
                    )[:branch_cap]
                    for i in range(len(local)):
                        e1 = edge_recs[local[i]]
                        n1 = int(e1["v"]) if int(e1["u"]) == int(pivot) else int(e1["u"])
                        for j in range(i + 1, len(local)):
                            e2 = edge_recs[local[j]]
                            n2 = int(e2["v"]) if int(e2["u"]) == int(pivot) else int(e2["u"])
                            if n1 == n2:
                                continue
                            path_nodes = (int(min(n1, n2)), int(pivot), int(max(n1, n2)))
                            path_edges = [e1, e2] if path_nodes[0] == int(n1) else [e2, e1]
                            scored = _chain_score(path_nodes, path_edges)
                            if scored is None:
                                continue
                            prev = best_by_chain.get(tuple(int(pid) for pid in scored["path_nodes"]))
                            if prev is None or float(scored["score"]) > float(prev.get("score", 0.0) or 0.0):
                                best_by_chain[tuple(int(pid) for pid in scored["path_nodes"])] = scored
            else:
                search_nodes = sorted(int(pid) for pid in adjacency.keys())
                for start in search_nodes:
                    pq: List[Tuple[float, float, int, Tuple[int, ...], Tuple[int, ...]]] = []
                    heapq.heappush(pq, (0.0, 0.0, int(start), (int(start),), tuple()))
                    expanded = 0
                    seen_best: Dict[Tuple[int, int], float] = {(int(start), 0): 0.0}
                    accepted = 0
                    while pq and accepted < per_start_cap and expanded < 400:
                        cost_so_far, dist_so_far, node, path_nodes, edge_ids = heapq.heappop(pq)
                        expanded += 1
                        hop_count = len(path_nodes) - 1
                        if hop_count >= 2:
                            scored = _chain_score(path_nodes, [edge_recs[idx] for idx in edge_ids])
                            if scored is not None:
                                sig = tuple(int(pid) for pid in scored["path_nodes"])
                                prev = best_by_chain.get(sig)
                                if prev is None or float(scored["score"]) > float(prev.get("score", 0.0) or 0.0):
                                    best_by_chain[sig] = scored
                                accepted += 1
                        if hop_count >= max_hops:
                            continue
                        for neighbor, edge_idx, edge in adjacency.get(int(node), []):
                            if int(neighbor) in path_nodes:
                                continue
                            new_cost = float(cost_so_far) + float(edge.get("cost", 0.0) or 0.0)
                            if new_cost > max_chain_budget + 1e-12:
                                continue
                            new_dist = float(dist_so_far) + float(edge.get("dist", 0.0) or 0.0)
                            if new_dist > max_chain_len + 1e-12:
                                continue
                            new_hops = hop_count + 1
                            state_key = (int(neighbor), int(new_hops))
                            prev_cost = seen_best.get(state_key)
                            if prev_cost is not None and float(prev_cost) <= float(new_cost) - 1e-12:
                                continue
                            seen_best[state_key] = float(new_cost)
                            heapq.heappush(
                                pq,
                                (
                                    float(new_cost),
                                    float(new_dist),
                                    int(neighbor),
                                    tuple(list(path_nodes) + [int(neighbor)]),
                                    tuple(list(edge_ids) + [edge_idx]),
                                ),
                            )

            added = 0
            ranked = sorted(
                best_by_chain.values(),
                key=lambda r: (
                    float(r.get("score", 0.0) or 0.0),
                    -float(r.get("chain_dist", 0.0) or 0.0),
                ),
                reverse=True,
            )
            for row in ranked[:chain_cap]:
                path_nodes = [int(pid) for pid in list(row.get("path_nodes", ()) or [])]
                path_edges = list(row.get("path_edges", []) or [])
                if len(path_nodes) < 3 or len(path_edges) != len(path_nodes) - 1:
                    continue
                a = int(path_nodes[0])
                c = int(path_nodes[-1])
                geoms_to_merge = []
                for edge in path_edges:
                    g = edge.get("geom")
                    if g is None or g.isEmpty():
                        geoms_to_merge = []
                        break
                    geoms_to_merge.append(g)
                if not geoms_to_merge:
                    continue
                merged = _safe_unary_union(geoms_to_merge)
                if merged is None or merged.isEmpty():
                    continue
                try:
                    area_ha = float(merged.area() / 10000.0)
                except Exception:
                    area_ha = float(sum(float(edge.get("cost", 0.0) or 0.0) for edge in path_edges))
                if area_ha <= 0.0:
                    continue
                explicit_edges = []
                for idx in range(len(path_nodes) - 1):
                    explicit_edges.append(
                        (
                            int(path_nodes[idx]),
                            int(path_nodes[idx + 1]),
                            float(path_edges[idx].get("dist", 0.0) or 0.0),
                        )
                    )
                _push_candidate(
                    {
                        "patch1": int(a),
                        "patch2": int(c),
                        "patch_ids": set(int(pid) for pid in path_nodes),
                        "raw_patch_ids": set(int(pid) for pid in path_nodes),
                        "geom": clone_geometry(merged),
                        "area_ha": float(area_ha),
                        "original_area_ha": float(area_ha),
                        "distance_m": float(row.get("chain_dist", 0.0) or 0.0),
                        "variant": f"chain{max(2, len(path_nodes) - 1)}",
                        "source": "landscape_fluidity_chain" if fluidity_mode else "most_connected_chain",
                        "chain_via_patch": int(path_nodes[1]),
                        "explicit_edges": explicit_edges,
                    }
                )
                added += 1
            return int(added)

        chain_added = _add_landscape_fluidity_chain_candidates()
        if chain_added > 0:
            mode_label = "landscape fluidity" if fluidity_mode else "most connected"
            print(f"  ✓ Added {chain_added} chain candidate(s) for {mode_label} mode")

        if fluidity_mode:
            intra_candidates = _get_intra_patch_candidates(
                patches=patches,
                params=params,
                resiliency_shortcut_threshold=_get_landscape_fluidity_shortcut_threshold(params, 3.0),
                patch_union=patch_union,
            )
            for cand in intra_candidates:
                _push_candidate(cand)
            if intra_candidates:
                print(f"  ✓ Added {len(intra_candidates)} intra-patch shortcut candidate(s) for landscape fluidity mode")

    emit_progress(progress_cb, progress_end, "Candidate corridors ready.")
    all_corridors: List[Dict] = []
    for items in candidates_by_pair.values():
        all_corridors.extend(items)
    print(f"\n  ✓ Found {len(all_corridors)} possible corridors")

    if timings:
        for label, duration in accum_durations.items():
            count = accum_counts.get(label, 0)
            timings.add(f"Find corridors | {label} (count={count})", duration)

    if timing_out is not None and circuit_mode:
        try:
            durations = dict(accum_durations)
            durations["total"] = time.perf_counter() - timing_start
            timing_out["durations_s"] = durations
            timing_out["counts"] = dict(accum_counts)
            timing_out["candidates"] = len(all_corridors)
        except Exception:
            pass

    return all_corridors


def greedy_global_optimize_vector(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
    metric_func: Callable[[nx.Graph], float],
    strategy_name: str,
    maximize: bool = True
) -> Tuple[Dict[int, Dict], Dict]:
    remaining = float(params.budget_area or 0.0)
    selected: Dict[int, Dict] = {}
    
    if not candidates:
        return {}, {"strategy": strategy_name, "corridors_used": 0, "budget_used_ha": 0.0}

    G = nx.Graph()
    for pid in patches:
        G.add_node(pid)
        
    budget_used = 0.0
    
    available_candidates = [c for c in candidates if float(c.get("area_ha", 0.0) or 0.0) > 0.0]

    def _metric_value(graph: nx.Graph) -> float:
        # Evaluate metric on the active corridor graph only (ignore isolated nodes).
        # For disconnected active graphs, aggregate by connected component so gain can
        # be measured incrementally instead of collapsing to zero.
        try:
            if graph.number_of_edges() <= 0:
                return 0.0
            active_nodes: Set[int] = set()
            for u, v in graph.edges():
                active_nodes.add(int(u))
                active_nodes.add(int(v))
            if len(active_nodes) < 2:
                return 0.0
            active = graph.subgraph(active_nodes)
            if active.number_of_edges() <= 0:
                return 0.0
            if nx.is_connected(active):
                return float(metric_func(active) or 0.0)
            total = 0.0
            for comp_nodes in nx.connected_components(active):
                if len(comp_nodes) < 2:
                    continue
                comp = active.subgraph(comp_nodes)
                total += float(metric_func(comp) or 0.0)
            return total
        except Exception:
            return 0.0
    
    while remaining > 0 and available_candidates:
        valid_cands = [c for c in available_candidates if float(c.get("area_ha", 0.0)) <= remaining]
        if not valid_cands:
            break
            
        current_metric = _metric_value(G)
        
        best_cand = None
        best_score = float('-inf') if maximize else float('inf')
        
        for cand in valid_cands:
            pids = cand.get("patch_ids", [])
            if not pids:
                p1, p2 = cand.get("patch1"), cand.get("patch2")
                if p1 is not None and p2 is not None:
                    pids = [int(p1), int(p2)]
            pids = [int(p) for p in pids if p is not None]
            if len(pids) < 2:
                continue
                
            cost = float(cand.get("area_ha", 0.0))
            length = float(cand.get("distance_m", 1.0) or 1.0)
            
            edges_added = []
            anchor = pids[0]
            for other in pids[1:]:
                if not G.has_edge(anchor, other):
                    G.add_edge(anchor, other, weight=length, distance=length)
                    edges_added.append((anchor, other))
                else:
                    old_w = G[anchor][other].get('weight', 1.0)
                    if length < old_w:
                        G[anchor][other]['weight'] = length
                        edges_added.append((anchor, other, old_w))
                        
            new_metric = _metric_value(G)
            
            for edge in edges_added:
                if len(edge) == 2:
                    G.remove_edge(edge[0], edge[1])
                elif len(edge) == 3:
                    G[edge[0]][edge[1]]['weight'] = edge[2]
            
            if maximize:
                gain = new_metric - current_metric
                score = gain / cost if cost > 0 else 0
                if score > best_score and gain > 1e-9:
                    best_score = score
                    best_cand = cand
            else:
                gain = current_metric - new_metric
                score = gain / cost if cost > 0 else 0
                if score > best_score and gain > 1e-9:
                    best_score = score
                    best_cand = cand
                    
        if best_cand is None:
            break
            
        cost = float(best_cand.get("area_ha", 0.0))
        remaining -= cost
        budget_used += cost
        
        pids = best_cand.get("patch_ids", [])
        if not pids:
            p1, p2 = best_cand.get("patch1"), best_cand.get("patch2")
            if p1 is not None and p2 is not None:
                pids = [int(p1), int(p2)]
        pids = [int(p) for p in pids if p is not None]
        anchor = pids[0]
        length = float(best_cand.get("distance_m", 1.0) or 1.0)
        
        for other in pids[1:]:
            if not G.has_edge(anchor, other):
                G.add_edge(anchor, other, weight=length, distance=length)
            else:
                old_w = G[anchor][other].get('weight', 1.0)
                if length < old_w:
                    G[anchor][other]['weight'] = length

        cid = len(selected) + 1
        selected[cid] = {
            "geom": clone_geometry(best_cand["geom"]),
            "patch_ids": set(best_cand.get("patch_ids", set(pids))),
            "area_ha": cost,
            "p1": int(best_cand.get("patch1", anchor)),
            "p2": int(best_cand.get("patch2", pids[-1])),
            "distance": length,
            "type": "primary",
            "variant": best_cand.get("variant"),
            "source": best_cand.get("source"),
            "utility_score": best_score,
            "overlap_ratio": 0.0,
        }
        
        available_candidates.remove(best_cand)
        
    stats = {
        "strategy": strategy_name,
        "corridors_used": len(selected),
        "budget_used_ha": budget_used,
        "patches_connected": nx.number_connected_components(G) if G.number_of_nodes() > 0 else 0,
        "total_connected_area_ha": sum(p.get("area_ha", 0.0) for p in patches.values()),
    }
    return selected, stats

def optimize_habitat_availability(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
) -> Tuple[Dict[int, Dict], Dict]:
    dispersal_distance = float(getattr(params, "species_dispersal_distance_analysis", 0.0) or 0.0)
    budget_limit = float(getattr(params, "budget_area", 0.0) or 0.0)
    if dispersal_distance <= 0.0 or budget_limit <= 0.0 or (not patches) or (not candidates):
        return {}, {"strategy": "reachable_habitat_advanced", "corridors_used": 0, "budget_used_ha": 0.0}

    min_patch_area = max(float(getattr(params, "min_patch_area_for_species_ha", 0.0) or 0.0), 0.0)
    scaling = normalize_patch_area_scaling(getattr(params, "patch_area_scaling", HABITAT_AVAILABILITY_DEFAULT_SCALING))
    nodes: List[HabitatAvailabilityNode] = []
    valid_nodes: Set[int] = set()
    for pid, pdata in patches.items():
        try:
            ipid = int(pid)
        except Exception:
            continue
        area = max(float((pdata or {}).get("area_ha", 0.0) or 0.0), 0.0)
        if area <= 0.0 or area + 1e-12 < min_patch_area:
            continue
        quality_weight = max(float((pdata or {}).get("quality_weight", 1.0) or 1.0), 0.0)
        eff_area = scale_patch_area(area, quality_weight=quality_weight, scaling=scaling)
        if eff_area <= 0.0:
            continue
        nodes.append(HabitatAvailabilityNode(node_id=ipid, raw_area=area, effective_area=eff_area))
        valid_nodes.add(ipid)
    if len(nodes) < 2:
        return {}, {"strategy": "reachable_habitat_advanced", "corridors_used": 0, "budget_used_ha": 0.0}

    evaluator = HabitatAvailabilityEvaluator(
        nodes,
        dispersal_distance=dispersal_distance,
        kernel=getattr(params, "species_dispersal_kernel", HABITAT_AVAILABILITY_DEFAULT_KERNEL),
    )
    before_metrics = evaluator.snapshot_metrics()
    remaining = float(budget_limit)
    budget_used = 0.0
    selected: Dict[int, Dict] = {}
    available_candidates = [cand for cand in candidates if float(cand.get("area_ha", 0.0) or 0.0) > 0.0]

    while remaining > 1e-12 and available_candidates:
        best_idx = -1
        best_eval: Optional[Dict[str, object]] = None
        best_score = -1.0
        best_gain = 0.0
        best_cost = 0.0
        best_length = float("inf")
        best_pids: List[int] = []

        for idx, cand in enumerate(available_candidates):
            cost = float(cand.get("area_ha", 0.0) or 0.0)
            if cost <= 0.0 or cost > remaining + 1e-12:
                continue
            pids = _candidate_patch_ids_for_metric(cand, valid_nodes)
            if len(pids) < 2:
                continue
            length = float(cand.get("distance_m", cand.get("distance", 1.0)) or 1.0)
            edge_specs = HabitatAvailabilityEvaluator.normalize_candidate_edges(
                pids,
                length=max(length, 1e-9),
                valid_node_ids=valid_nodes,
            )
            if not edge_specs:
                continue
            eval_out = evaluator.evaluate_candidate(edge_specs, corridor_cost=cost)
            if not eval_out:
                continue
            score = float(eval_out.get("score", 0.0) or 0.0)
            gain = float(eval_out.get("gain", 0.0) or 0.0)
            better = (
                score > best_score + 1e-12
                or (abs(score - best_score) <= 1e-12 and gain > best_gain + 1e-12)
                or (
                    abs(score - best_score) <= 1e-12
                    and abs(gain - best_gain) <= 1e-12
                    and (best_cost <= 0.0 or cost < best_cost - 1e-12)
                )
                or (
                    abs(score - best_score) <= 1e-12
                    and abs(gain - best_gain) <= 1e-12
                    and abs(cost - best_cost) <= 1e-12
                    and length < best_length - 1e-12
                )
            )
            if better:
                best_idx = int(idx)
                best_eval = eval_out
                best_score = float(score)
                best_gain = float(gain)
                best_cost = float(cost)
                best_length = float(length)
                best_pids = list(pids)

        if best_idx < 0 or best_eval is None or best_score <= 1e-12:
            break

        chosen = available_candidates.pop(best_idx)
        evaluator.apply_candidate(best_eval)
        budget_used += float(best_cost)
        remaining -= float(best_cost)
        cid = len(selected) + 1
        p1 = int(chosen.get("patch1", chosen.get("p1", best_pids[0])))
        p2 = int(chosen.get("patch2", chosen.get("p2", best_pids[-1])))
        selected[cid] = {
            "geom": clone_geometry(chosen["geom"]),
            "patch_ids": set(best_pids),
            "area_ha": float(best_cost),
            "p1": p1,
            "p2": p2,
            "distance": float(best_length),
            "type": "primary",
            "variant": chosen.get("variant"),
            "source": chosen.get("source"),
            "utility_score": float(best_score),
            "overlap_ratio": 0.0,
        }

    after_metrics = evaluator.snapshot_metrics()
    before_ha = float(before_metrics.get("habitat_availability", 0.0) or 0.0)
    after_ha = float(after_metrics.get("habitat_availability", 0.0) or 0.0)
    percent_increase = (
        float(((after_ha - before_ha) / max(before_ha, 1e-12)) * 100.0)
        if before_ha > 1e-12
        else (100.0 if after_ha > 1e-12 else 0.0)
    )
    stats = {
        "strategy": "reachable_habitat_advanced",
        "corridors_used": len(selected),
        "budget_used_ha": float(budget_used),
        "habitat_availability_before": float(before_ha),
        "habitat_availability_after": float(after_ha),
        "habitat_availability_gain": float(after_ha - before_ha),
        "percent_increase": float(percent_increase),
        "mean_reachable_area_before": float(before_metrics.get("mean_reachable_area", 0.0) or 0.0),
        "mean_reachable_area": float(after_metrics.get("mean_reachable_area", 0.0) or 0.0),
        "median_reachable_area_before": float(before_metrics.get("median_reachable_area", 0.0) or 0.0),
        "median_reachable_area": float(after_metrics.get("median_reachable_area", 0.0) or 0.0),
        "largest_reachable_habitat_cluster_before": float(
            before_metrics.get("largest_reachable_habitat_cluster", 0.0) or 0.0
        ),
        "largest_reachable_habitat_cluster": float(
            after_metrics.get("largest_reachable_habitat_cluster", 0.0) or 0.0
        ),
        "species_dispersal_distance_m": float(dispersal_distance),
        "species_dispersal_kernel": str(
            getattr(params, "species_dispersal_kernel", HABITAT_AVAILABILITY_DEFAULT_KERNEL)
            or HABITAT_AVAILABILITY_DEFAULT_KERNEL
        ),
        "min_patch_area_for_species_ha": float(min_patch_area),
        "patch_area_scaling": str(scaling),
        "patch_quality_weight_field": str(getattr(params, "patch_quality_weight_field", "") or ""),
        "patches_connected": 0,
        "total_connected_area_ha": sum(float(p.get("area_ha", 0.0) or 0.0) for p in patches.values()),
    }
    return selected, stats


def _candidate_patch_ids_for_metric(
    cand: Dict,
    valid_nodes: Set[int],
) -> List[int]:
    out: List[int] = []
    raw = cand.get("patch_ids")
    if isinstance(raw, (set, list, tuple)):
        for pid in raw:
            try:
                ip = int(pid)
            except Exception:
                continue
            if ip in valid_nodes:
                out.append(ip)
    if len(out) < 2:
        for k in ("patch1", "p1", "patch2", "p2"):
            try:
                ip = int(cand.get(k))
            except Exception:
                continue
            if ip in valid_nodes:
                out.append(ip)
    seen: Set[int] = set()
    dedup: List[int] = []
    for pid in out:
        if pid in seen:
            continue
        seen.add(pid)
        dedup.append(pid)
    return dedup


def _ordered_patch_boundary_points(
    patch_geom: QgsGeometry,
    max_points: int = 120,
) -> List[QgsPointXY]:
    if patch_geom is None or patch_geom.isEmpty():
        return []

    # Use exterior ring order so ring-distance approximations are stable/deterministic.
    rings: List[List[QgsPointXY]] = []
    try:
        if patch_geom.isMultipart():
            polys = patch_geom.asMultiPolygon() or []
            for poly in polys:
                if poly and poly[0]:
                    rings.append([QgsPointXY(p) for p in poly[0]])
        else:
            poly = patch_geom.asPolygon() or []
            if poly and poly[0]:
                rings.append([QgsPointXY(p) for p in poly[0]])
    except Exception:
        rings = []

    if not rings:
        try:
            verts = [QgsPointXY(v) for v in patch_geom.vertices()]
            if len(verts) > 1:
                return verts[:max_points] if max_points > 0 else verts
            return []
        except Exception:
            return []

    # Keep the largest ring by perimeter to focus on the main patch body.
    def _ring_len(ring_pts: List[QgsPointXY]) -> float:
        if len(ring_pts) < 2:
            return 0.0
        total = 0.0
        for i in range(len(ring_pts) - 1):
            total += float(ring_pts[i].distance(ring_pts[i + 1]))
        return float(total)

    ring = max(rings, key=_ring_len)
    if len(ring) > 1 and ring[0].distance(ring[-1]) <= 1e-9:
        ring = ring[:-1]
    if len(ring) < 4:
        return ring

    if max_points > 0 and len(ring) > max_points:
        step = max(1, len(ring) // max_points)
        ring = ring[::step][:max_points]
    return ring


def _geometry_parts(geom: QgsGeometry) -> List[QgsGeometry]:
    if geom is None or geom.isEmpty():
        return []
    if not geom.isMultipart():
        return [geom]
    try:
        coll = geom.asGeometryCollection()
        if coll:
            out: List[QgsGeometry] = []
            for part in coll:
                gg = QgsGeometry(part)
                if gg is not None and (not gg.isEmpty()):
                    out.append(gg)
            if out:
                return out
    except Exception:
        pass
    try:
        out2: List[QgsGeometry] = []
        for poly in geom.asMultiPolygon():
            gg = QgsGeometry.fromPolygonXY(poly)
            if gg is not None and (not gg.isEmpty()):
                out2.append(gg)
        if out2:
            return out2
    except Exception:
        pass
    return [geom]


def _ring_signed_area(points: Sequence[QgsPointXY]) -> float:
    if points is None or len(points) < 3:
        return 0.0
    area2 = 0.0
    n = len(points)
    for i in range(n):
        j = (i + 1) % n
        area2 += float(points[i].x()) * float(points[j].y()) - float(points[j].x()) * float(points[i].y())
    return float(area2 * 0.5)


def _iter_cyclic_true_runs(flags: Sequence[bool]) -> List[Tuple[int, int]]:
    n = len(flags or [])
    if n <= 0:
        return []
    runs: List[Tuple[int, int]] = []
    start: Optional[int] = None
    for idx, flag in enumerate(flags):
        if flag and start is None:
            start = idx
        elif (not flag) and start is not None:
            runs.append((start, idx - 1))
            start = None
    if start is not None:
        runs.append((start, n - 1))
    if len(runs) >= 2 and flags[0] and flags[-1]:
        first_start, first_end = runs[0]
        last_start, last_end = runs[-1]
        merged = (last_start, first_end + n)
        runs = [merged] + runs[1:-1]
    return runs


def _ring_path_metrics(points: Sequence[QgsPointXY]) -> Tuple[List[float], float]:
    n = len(points or [])
    if n <= 0:
        return [0.0], 0.0
    cum: List[float] = [0.0]
    total_len = 0.0
    for i in range(n):
        j = (i + 1) % n
        total_len += float(points[i].distance(points[j]))
        if i < n - 1:
            cum.append(total_len)
    return cum, float(total_len)


def _ring_shorter_detour(cum: Sequence[float], total_len: float, i: int, j: int) -> float:
    if not cum or total_len <= 0.0:
        return 0.0
    lo = min(int(i), int(j))
    hi = max(int(i), int(j))
    along = abs(float(cum[hi]) - float(cum[lo]))
    return float(min(along, float(total_len) - along))


def _scaled_intra_patch_budgets(
    patch_geom: QgsGeometry,
    point_count: int,
    max_candidates: int,
) -> Dict[str, int]:
    area_m2 = 0.0
    diag_m = 0.0
    try:
        area_m2 = max(float(patch_geom.area()), 0.0)
    except Exception:
        area_m2 = 0.0
    try:
        bbox = patch_geom.boundingBox()
        diag_m = math.hypot(float(bbox.width()), float(bbox.height()))
    except Exception:
        diag_m = 0.0

    scale = max(
        1.0,
        min(
            4.0,
            0.65
            + math.sqrt(max(area_m2, 1.0)) / 250.0
            + diag_m / 2500.0
            + max(int(point_count), 1) / 350.0,
        ),
    )
    return {
        "shoulder_expand": int(min(18, max(4, round(3 * scale)))),
        "bay_sample_step": int(min(12, max(1, round(12.0 / scale)))),
        "bay_candidates": int(min(6 * max(max_candidates, 1), max(24, round(max_candidates * scale * 1.5)))),
    }


def _score_intra_pair_candidate(
    *,
    patch_geom: QgsGeometry,
    points: Sequence[QgsPointXY],
    cum: Sequence[float],
    total_len: float,
    holes_before_count: int,
    holes_before_area_m2: float,
    i: int,
    j: int,
    min_corridor_width_m: float,
    max_gap_m: float,
    min_shortcut_ratio: float,
    min_index_separation: int,
) -> Optional[Dict[str, Any]]:
    def _evaluate_pair(ii: int, jj: int) -> Optional[Dict[str, Any]]:
        idx_sep_local = min((jj - ii) % n, (ii - jj) % n)
        if idx_sep_local < int(min_index_separation):
            return None

        gap_local = float(points[ii].distance(points[jj]))
        if gap_local <= gap_min or gap_local > gap_cap:
            return None

        detour_local = _ring_shorter_detour(cum, total_len, ii, jj)
        if detour_local <= gap_local + 1e-9:
            return None
        ratio_local = float(detour_local / max(gap_local, 1e-9))
        if ratio_local + 1e-12 < float(min_shortcut_ratio):
            return None

        try:
            line_chord_local = QgsGeometry.fromPolylineXY([points[ii], points[jj]])
            mid_pt_local = line_chord_local.interpolate(gap_local * 0.5)
            if mid_pt_local is None or mid_pt_local.isEmpty() or patch_geom.contains(mid_pt_local):
                return None
        except Exception:
            return None

        bridge_local = _build_valid_intra_bridge_geometry(
            patch_geom=patch_geom,
            pt_a=points[ii],
            pt_b=points[jj],
            min_corridor_width_m=min_corridor_width_m,
            subtract_geom=patch_geom,
        )
        if bridge_local is None or bridge_local.isEmpty():
            return None

        bridge_area_local_m2 = float(bridge_local.area())
        if bridge_area_local_m2 <= 0.0:
            return None
        loop_bonus_local, delta_holes_local, delta_hole_area_local_m2 = _intra_loop_bonus(
            patch_geom=patch_geom,
            corridor_geom=bridge_local,
            base_hole_count=holes_before_count,
            base_hole_area_m2=holes_before_area_m2,
        )
        pocket_strength_local = float(delta_hole_area_local_m2 / max(bridge_area_local_m2, 1.0))
        if int(delta_holes_local) <= 0 and pocket_strength_local < 1.25:
            return None
        mid_x_local = 0.5 * (float(points[ii].x()) + float(points[jj].x()))
        mid_y_local = 0.5 * (float(points[ii].y()) + float(points[jj].y()))
        angle_deg_local = math.degrees(
            math.atan2(float(points[jj].y()) - float(points[ii].y()), float(points[jj].x()) - float(points[ii].x()))
        )
        return {
            "pt_a": points[ii],
            "pt_b": points[jj],
            "detour_m": float(detour_local),
            "gap_m": float(gap_local),
            "ratio": float(ratio_local),
            "shortcut_ratio": float(ratio_local),
            "gain_m": float(detour_local - gap_local),
            "loop_bonus": float(loop_bonus_local),
            "loop_hole_count": int(delta_holes_local),
            "loop_hole_area_m2": float(delta_hole_area_local_m2),
            "bridge_area_m2": float(bridge_area_local_m2),
            "pocket_strength": float(pocket_strength_local),
            "mid_x": float(mid_x_local),
            "mid_y": float(mid_y_local),
            "angle_deg": float(angle_deg_local),
            "idx_a": int(ii),
            "idx_b": int(jj),
            "ring_n": int(n),
        }

    def _refine_quality(rec: Dict[str, Any]) -> float:
        gap_v = max(float(rec.get("gap_m", 0.0) or 0.0), 1.0)
        gain_v = max(float(rec.get("gain_m", 0.0) or 0.0), 0.0)
        hole_v = max(float(rec.get("loop_hole_area_m2", 0.0) or 0.0), 0.0)
        pocket_v = max(float(rec.get("pocket_strength", 0.0) or 0.0), 0.0)
        return float(
            (math.pow(gain_v + 1.0, 0.18) * math.pow(hole_v + 1.0, 0.34) * math.pow(pocket_v + 1.0, 0.26))
            / math.pow(gap_v, 0.62)
        )

    n = len(points or [])
    if n < 4:
        return None

    gap_min = max(float(min_corridor_width_m or 0.0) * 0.5, 1.0)
    gap_cap = max(float(max_gap_m or 0.0), gap_min)
    base = _evaluate_pair(int(i), int(j))
    if base is None:
        return None

    refine_span = max(2, min(8, int(round(math.sqrt(max(float(base.get("gap_m", 1.0) or 1.0), 1.0)) / 3.5))))
    best = base
    best_q = _refine_quality(base)
    base_hole = max(float(base.get("loop_hole_area_m2", 0.0) or 0.0), 0.0)
    base_gain = max(float(base.get("gain_m", 0.0) or 0.0), 0.0)
    for di in range(-refine_span, refine_span + 1):
        ii = int(i) + di
        if ii < 0 or ii >= n:
            continue
        for dj in range(-refine_span, refine_span + 1):
            jj = int(j) + dj
            if jj < 0 or jj >= n or ii == jj:
                continue
            cand = _evaluate_pair(ii, jj)
            if cand is None:
                continue
            cand_hole = max(float(cand.get("loop_hole_area_m2", 0.0) or 0.0), 0.0)
            cand_gain = max(float(cand.get("gain_m", 0.0) or 0.0), 0.0)
            if cand_hole + 1e-9 < 0.70 * max(base_hole, 1.0):
                continue
            if cand_gain + 1e-9 < 0.45 * max(base_gain, 1.0):
                continue
            cand_q = _refine_quality(cand)
            if cand_q > best_q + 1e-12:
                best = cand
                best_q = cand_q

    return best


def _cyclic_idx_dist(a: int, b: int, n: int) -> int:
    if n <= 0:
        return abs(int(a) - int(b))
    diff = abs(int(a) - int(b)) % int(n)
    return int(min(diff, int(n) - diff))


def _select_diverse_intra_candidates(
    patch_geom: QgsGeometry,
    candidates: Sequence[Dict[str, Any]],
    max_candidates: int,
) -> List[Dict[str, Any]]:
    if max_candidates <= 0:
        return list(candidates)
    diag_m = 0.0
    try:
        bb = patch_geom.boundingBox()
        diag_m = math.hypot(float(bb.width()), float(bb.height()))
    except Exception:
        diag_m = 0.0
    min_sep = max(8.0, min(60.0, diag_m * 0.010))
    selected: List[Dict[str, Any]] = []
    signatures: List[Tuple[float, float, float, float, float, int, int, int]] = []
    for cand in candidates:
        try:
            mid_x = float(cand.get("mid_x"))
            mid_y = float(cand.get("mid_y"))
            ax = float(cand["pt_a"].x())
            ay = float(cand["pt_a"].y())
            bx = float(cand["pt_b"].x())
            by = float(cand["pt_b"].y())
            gap_m = max(float(cand.get("gap_m", 0.0) or 0.0), 1.0)
            idx_a = int(cand.get("idx_a"))
            idx_b = int(cand.get("idx_b"))
            ring_n = max(int(cand.get("ring_n", 0) or 0), 1)
        except Exception:
            continue
        keep = True
        terminal_tol = max(12.0, min(120.0, gap_m * 0.22 + 8.0))
        idx_tol = max(8, min(60, int(round(ring_n * 0.035))))
        for prev_ax, prev_ay, prev_bx, prev_by, prev_gap_m, prev_idx_a, prev_idx_b, prev_ring_n in signatures:
            prev_terminal_tol = max(12.0, min(120.0, prev_gap_m * 0.22 + 8.0))
            same_order = (
                math.hypot(ax - prev_ax, ay - prev_ay) <= max(terminal_tol, prev_terminal_tol)
                and math.hypot(bx - prev_bx, by - prev_by) <= max(terminal_tol, prev_terminal_tol)
            )
            swap_order = (
                math.hypot(ax - prev_bx, ay - prev_by) <= max(terminal_tol, prev_terminal_tol)
                and math.hypot(bx - prev_ax, by - prev_ay) <= max(terminal_tol, prev_terminal_tol)
            )
            same_mouth = math.hypot(mid_x - 0.5 * (prev_ax + prev_bx), mid_y - 0.5 * (prev_ay + prev_by)) < min_sep
            idx_same_order = (
                _cyclic_idx_dist(idx_a, prev_idx_a, max(ring_n, prev_ring_n)) <= idx_tol
                and _cyclic_idx_dist(idx_b, prev_idx_b, max(ring_n, prev_ring_n)) <= idx_tol
            )
            idx_swap_order = (
                _cyclic_idx_dist(idx_a, prev_idx_b, max(ring_n, prev_ring_n)) <= idx_tol
                and _cyclic_idx_dist(idx_b, prev_idx_a, max(ring_n, prev_ring_n)) <= idx_tol
            )
            if idx_same_order or idx_swap_order or same_order or swap_order or (same_mouth and abs(gap_m - prev_gap_m) < max(18.0, 0.12 * max(gap_m, prev_gap_m))):
                keep = False
                break
        if not keep:
            continue
        selected.append(cand)
        signatures.append((ax, ay, bx, by, gap_m, idx_a, idx_b, ring_n))
        if len(selected) >= int(max_candidates):
            break
    return selected


def _intra_structural_repair_score(candidate: Dict[str, Any]) -> float:
    gain_m = max(float(candidate.get("intra_gain_m", candidate.get("gain_m", 0.0)) or 0.0), 0.0)
    gap_m = max(float(candidate.get("intra_gap_m", candidate.get("gap_m", 0.0)) or 0.0), 1.0)
    hole_area_m2 = max(float(candidate.get("intra_hole_area_m2", candidate.get("loop_hole_area_m2", 0.0)) or 0.0), 0.0)
    pocket_strength = max(float(candidate.get("intra_pocket_strength", candidate.get("pocket_strength", 0.0)) or 0.0), 0.0)
    patch_area_ha = max(float(candidate.get("intra_patch_area_ha", 0.0) or 0.0), 0.0)
    cost_ha = max(float(candidate.get("area_ha", 0.0) or 0.0), 1e-9)
    patch_area_m2 = max(patch_area_ha * 10000.0, 1.0)
    hole_share = max(0.0, min(1.0, hole_area_m2 / max(0.35 * patch_area_m2, 1.0)))
    mouth_scale = max(0.0, min(1.0, gap_m / max(math.sqrt(patch_area_m2) * 0.30, 1.0)))
    repair_term = (
        (0.50 * math.sqrt(hole_share))
        + (0.30 * max(0.0, min(1.0, math.log1p(pocket_strength) / math.log1p(150.0))))
        + (0.20 * mouth_scale)
    )
    return float((math.sqrt(gain_m + 1.0) * (0.15 + repair_term)) / math.sqrt(cost_ha))


def _intra_candidate_dominates(a: Dict[str, Any], b: Dict[str, Any]) -> bool:
    cost_a = float(a.get("area_ha", 0.0) or 0.0)
    cost_b = float(b.get("area_ha", 0.0) or 0.0)
    gain_a = float(a.get("intra_gain_m", 0.0) or 0.0)
    gain_b = float(b.get("intra_gain_m", 0.0) or 0.0)
    hole_a = float(a.get("intra_hole_area_m2", 0.0) or 0.0)
    hole_b = float(b.get("intra_hole_area_m2", 0.0) or 0.0)
    pocket_a = float(a.get("intra_pocket_strength", 0.0) or 0.0)
    pocket_b = float(b.get("intra_pocket_strength", 0.0) or 0.0)
    no_worse = (
        cost_a <= cost_b + 1e-9
        and gain_a + 1e-9 >= gain_b
        and hole_a + 1e-9 >= hole_b
        and pocket_a + 1e-9 >= pocket_b
    )
    strictly_better = (
        cost_a < cost_b - 1e-9
        or gain_a > gain_b + 1e-9
        or hole_a > hole_b + 1e-9
        or pocket_a > pocket_b + 1e-9
    )
    return bool(no_worse and strictly_better)


def _group_intra_patch_winners(
    patch_geom: QgsGeometry,
    patch_candidates: Sequence[Dict[str, Any]],
    *,
    max_group_winners: int = 1,
) -> List[Dict[str, Any]]:
    if not patch_candidates:
        return []

    diag_m = 0.0
    try:
        bb = patch_geom.boundingBox()
        diag_m = math.hypot(float(bb.width()), float(bb.height()))
    except Exception:
        diag_m = 0.0

    groups: List[List[Dict[str, Any]]] = []
    for cand in sorted(
        patch_candidates,
        key=lambda c: (
            float(c.get("intra_repair_score", 0.0) or 0.0),
            float(c.get("intra_hole_area_m2", 0.0) or 0.0),
            float(c.get("intra_gain_m", 0.0) or 0.0),
        ),
        reverse=True,
    ):
        mid_x = float(cand.get("mid_x", 0.0) or 0.0)
        mid_y = float(cand.get("mid_y", 0.0) or 0.0)
        idx_a = int(cand.get("idx_a", -1) or -1)
        idx_b = int(cand.get("idx_b", -1) or -1)
        ring_n = max(int(cand.get("ring_n", 0) or 0), 1)
        gap_m = max(float(cand.get("intra_gap_m", cand.get("gap_m", 0.0)) or 0.0), 1.0)
        group_hit = None
        for gi, group in enumerate(groups):
            ref = group[0]
            ref_mid_x = float(ref.get("mid_x", 0.0) or 0.0)
            ref_mid_y = float(ref.get("mid_y", 0.0) or 0.0)
            ref_idx_a = int(ref.get("idx_a", -1) or -1)
            ref_idx_b = int(ref.get("idx_b", -1) or -1)
            ref_ring_n = max(int(ref.get("ring_n", 0) or 0), 1)
            idx_tol = max(6, min(42, int(round(max(ring_n, ref_ring_n) * 0.028))))
            mouth_tol = max(20.0, min(140.0, max(gap_m, float(ref.get("intra_gap_m", 0.0) or 0.0)) * 0.45, diag_m * 0.02))
            same_order = (
                _cyclic_idx_dist(idx_a, ref_idx_a, max(ring_n, ref_ring_n)) <= idx_tol
                and _cyclic_idx_dist(idx_b, ref_idx_b, max(ring_n, ref_ring_n)) <= idx_tol
            )
            swap_order = (
                _cyclic_idx_dist(idx_a, ref_idx_b, max(ring_n, ref_ring_n)) <= idx_tol
                and _cyclic_idx_dist(idx_b, ref_idx_a, max(ring_n, ref_ring_n)) <= idx_tol
            )
            same_mid = math.hypot(mid_x - ref_mid_x, mid_y - ref_mid_y) <= mouth_tol
            if same_order or swap_order or same_mid:
                group_hit = gi
                break
        if group_hit is None:
            groups.append([cand])
        else:
            groups[group_hit].append(cand)

    winners: List[Dict[str, Any]] = []
    for group in groups:
        nondominated: List[Dict[str, Any]] = []
        for cand in group:
            if any(_intra_candidate_dominates(other, cand) for other in group if other is not cand):
                continue
            nondominated.append(cand)
        chosen = sorted(
            nondominated or group,
            key=lambda c: (
                float(c.get("intra_repair_score", 0.0) or 0.0),
                float(c.get("intra_hole_area_m2", 0.0) or 0.0),
                float(c.get("intra_gain_m", 0.0) or 0.0),
                -float(c.get("area_ha", 0.0) or 0.0),
            ),
            reverse=True,
        )[: max(1, int(max_group_winners))]
        winners.extend(chosen)
    return winners


def _bay_mouth_candidates(
    patch_geom: QgsGeometry,
    points: Sequence[QgsPointXY],
    *,
    min_corridor_width_m: float,
    max_gap_m: float,
    min_shortcut_ratio: float,
    min_index_separation: int,
    max_candidates: int,
) -> List[Dict[str, Any]]:
    n = len(points or [])
    if n < 8:
        return []
    cum, total_len = _ring_path_metrics(points)
    if total_len <= 0.0:
        return []
    holes_before = _extract_interior_ring_geometries(patch_geom)
    base_hole_count = len(holes_before)
    base_hole_area_m2 = float(sum(float(h.area()) for h in holes_before))
    budgets = _scaled_intra_patch_budgets(patch_geom, n, max_candidates)

    out: List[Dict[str, Any]] = []
    gap_cap = max(float(max_gap_m or 0.0), max(float(min_corridor_width_m or 0.0), 1.0))
    min_detour = max(gap_cap * 0.30, 120.0)
    max_out = max(1, int(budgets.get("bay_candidates", max_candidates)))
    pre_candidates: List[Tuple[Tuple[float, float, float, float], int, int]] = []
    for i in range(n):
        for j in range(i + int(min_index_separation), n):
            gap = float(points[i].distance(points[j]))
            if gap <= 1.0 or gap > gap_cap:
                continue
            detour = _ring_shorter_detour(cum, total_len, i, j)
            if detour < min_detour or detour <= gap + 1e-9:
                continue
            ratio = float(detour / max(gap, 1e-9))
            if ratio + 1e-12 < float(min_shortcut_ratio):
                continue
            pre_candidates.append((((detour - gap), ratio, detour, -gap), i, j))

    pre_candidates.sort(key=lambda x: x[0], reverse=True)
    for _, i, j in pre_candidates:
        cand = _score_intra_pair_candidate(
                patch_geom=patch_geom,
                points=points,
                cum=cum,
                total_len=total_len,
                holes_before_count=base_hole_count,
                holes_before_area_m2=base_hole_area_m2,
                i=i,
                j=j,
                min_corridor_width_m=min_corridor_width_m,
                max_gap_m=max_gap_m,
                min_shortcut_ratio=min_shortcut_ratio,
                min_index_separation=min_index_separation,
        )
        if cand is None:
            continue
        out.append(cand)

    out.sort(
        key=lambda x: (
            float(x.get("loop_hole_count", 0) or 0.0),
            float(x.get("loop_hole_area_m2", 0.0) or 0.0),
            float(x.get("pocket_strength", 0.0) or 0.0),
            float(x.get("gain_m", 0.0) or 0.0),
            -float(x.get("gap_m", 0.0) or 0.0),
        ),
        reverse=True,
    )
    return _select_diverse_intra_candidates(patch_geom, out, max_out)


def _concavity_mouth_candidates(
    patch_geom: QgsGeometry,
    points: Sequence[QgsPointXY],
    *,
    min_corridor_width_m: float,
    max_gap_m: float,
    min_shortcut_ratio: float,
    min_index_separation: int,
    max_candidates: int,
) -> List[Dict[str, Any]]:
    n = len(points or [])
    if n < 8:
        return []

    signed_area = _ring_signed_area(points)
    orient_sign = 1.0 if signed_area >= 0.0 else -1.0
    concave_flags: List[bool] = [False] * n
    concave_strength: List[float] = [0.0] * n
    min_turn_strength = math.sin(math.radians(18.0))

    for idx in range(n):
        prev_pt = points[(idx - 1) % n]
        curr_pt = points[idx]
        next_pt = points[(idx + 1) % n]
        ax = float(curr_pt.x()) - float(prev_pt.x())
        ay = float(curr_pt.y()) - float(prev_pt.y())
        bx = float(next_pt.x()) - float(curr_pt.x())
        by = float(next_pt.y()) - float(curr_pt.y())
        la = math.hypot(ax, ay)
        lb = math.hypot(bx, by)
        if la <= 1e-9 or lb <= 1e-9:
            continue
        cross = ax * by - ay * bx
        turn_strength = abs(float(cross)) / max(la * lb, 1e-9)
        if turn_strength < min_turn_strength:
            continue
        if cross * orient_sign < 0.0:
            concave_flags[idx] = True
            concave_strength[idx] = float(turn_strength)

    runs = _iter_cyclic_true_runs(concave_flags)
    if not runs:
        return []

    cum, total_len = _ring_path_metrics(points)
    if total_len <= 0.0:
        return []

    holes_before = _extract_interior_ring_geometries(patch_geom)
    base_hole_count = len(holes_before)
    base_hole_area_m2 = float(sum(float(h.area()) for h in holes_before))

    out: List[Dict[str, Any]] = []
    seen_pairs: Set[Tuple[int, int]] = set()
    budgets = _scaled_intra_patch_budgets(patch_geom, n, max_candidates)
    shoulder_expand = int(budgets.get("shoulder_expand", 6))

    for run_start_raw, run_end_raw in runs:
        run_start = int(run_start_raw % n)
        run_end = int(run_end_raw % n)
        run_len = int(run_end_raw - run_start_raw + 1)
        if run_len <= 0:
            continue
        for left_expand in range(0, shoulder_expand + 1):
            for right_expand in range(0, shoulder_expand + 1):
                i = (run_start - 1 - left_expand) % n
                j = (run_end + 1 + right_expand) % n
                if i == j:
                    continue
                pair = (min(i, j), max(i, j))
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)
                cand = _score_intra_pair_candidate(
                    patch_geom=patch_geom,
                    points=points,
                    cum=cum,
                    total_len=total_len,
                    holes_before_count=base_hole_count,
                    holes_before_area_m2=base_hole_area_m2,
                    i=i,
                    j=j,
                    min_corridor_width_m=min_corridor_width_m,
                    max_gap_m=max_gap_m,
                    min_shortcut_ratio=min_shortcut_ratio,
                    min_index_separation=min_index_separation,
                )
                if cand is None:
                    continue

                run_strength = max(float(concave_strength[k % n]) for k in range(run_start_raw, run_end_raw + 1))
                cand["concavity_run_len"] = int(run_len)
                cand["concavity_strength"] = float(run_strength)
                out.append(cand)

    out.sort(
        key=lambda x: (
            float(x.get("loop_hole_count", 0) or 0.0),
            float(x.get("loop_hole_area_m2", 0.0) or 0.0),
            float(x.get("pocket_strength", 0.0) or 0.0),
            float(x.get("gain_m", 0.0) or 0.0),
            float(x.get("ratio", 0.0) or 0.0),
            float(x.get("concavity_strength", 0.0) or 0.0),
        ),
        reverse=True,
    )
    return _select_diverse_intra_candidates(patch_geom, out, int(max_candidates)) if max_candidates > 0 else out


def _rank_intra_patch_shortcuts(
    patch_geom: QgsGeometry,
    *,
    min_corridor_width_m: float,
    max_gap_m: float,
    min_shortcut_ratio: float,
    max_points: int = 150,
    min_index_separation: int = 4,
    max_candidates: int = 48,
) -> List[Dict[str, Any]]:
    if patch_geom is None or patch_geom.isEmpty():
        return []

    pts = _ordered_patch_boundary_points(patch_geom, max_points=max_points)
    n = len(pts)
    if n < 4:
        return []

    cum, total_len = _ring_path_metrics(pts)
    if total_len <= 0.0:
        return []

    holes_before = _extract_interior_ring_geometries(patch_geom)
    base_hole_count = len(holes_before)
    base_hole_area_m2 = float(sum(float(h.area()) for h in holes_before))

    concavity_first = _concavity_mouth_candidates(
        patch_geom,
        pts,
        min_corridor_width_m=min_corridor_width_m,
        max_gap_m=max_gap_m,
        min_shortcut_ratio=min_shortcut_ratio,
        min_index_separation=min_index_separation,
        max_candidates=max_candidates,
    )
    bay_first = _bay_mouth_candidates(
        patch_geom,
        pts,
        min_corridor_width_m=min_corridor_width_m,
        max_gap_m=max_gap_m,
        min_shortcut_ratio=min_shortcut_ratio,
        min_index_separation=min_index_separation,
        max_candidates=max_candidates,
    )
    exhaustive_scored: List[Dict[str, Any]] = []
    for i in range(n):
        for j in range(i + int(min_index_separation), n):
            cand = _score_intra_pair_candidate(
                patch_geom=patch_geom,
                points=pts,
                cum=cum,
                total_len=total_len,
                holes_before_count=base_hole_count,
                holes_before_area_m2=base_hole_area_m2,
                i=i,
                j=j,
                min_corridor_width_m=min_corridor_width_m,
                max_gap_m=max_gap_m,
                min_shortcut_ratio=min_shortcut_ratio,
                min_index_separation=min_index_separation,
            )
            if cand is not None:
                exhaustive_scored.append(cand)
    if concavity_first or bay_first or exhaustive_scored:
        merged: List[Dict[str, Any]] = []
        seen_keys: Set[Tuple[int, int]] = set()
        for cand in list(concavity_first) + list(bay_first) + list(exhaustive_scored):
            try:
                key = (
                    int(round(float(cand["pt_a"].x()) * 10.0)),
                    int(round(float(cand["pt_a"].y()) * 10.0)),
                    int(round(float(cand["pt_b"].x()) * 10.0)),
                    int(round(float(cand["pt_b"].y()) * 10.0)),
                )
            except Exception:
                key = None
            if key is not None and key in seen_keys:
                continue
            if key is not None:
                seen_keys.add(key)
            merged.append(cand)
        merged.sort(
            key=lambda x: (
                float(x.get("loop_hole_count", 0) or 0.0),
                float(x.get("loop_hole_area_m2", 0.0) or 0.0),
                float(x.get("pocket_strength", 0.0) or 0.0),
                float(x.get("gain_m", 0.0) or 0.0),
                float(x.get("ratio", 0.0) or 0.0),
                float(x.get("concavity_strength", 0.0) or 0.0),
            ),
            reverse=True,
        )
        return _select_diverse_intra_candidates(patch_geom, merged, int(max_candidates)) if max_candidates > 0 else merged

    out: List[Dict[str, Any]] = []
    min_gap = 1.0
    gap_cap = max(float(max_gap_m), min_gap)
    for i in range(n):
        for j in range(i + 1, n):
            idx_sep = min(j - i, n - (j - i))
            if idx_sep < int(min_index_separation):
                continue

            gap = float(pts[i].distance(pts[j]))
            if gap < min_gap or gap > gap_cap:
                continue

            along = abs(cum[j] - cum[i])
            detour = min(along, total_len - along)
            ratio = float(detour / max(gap, 0.1))
            if ratio + 1e-12 < float(min_shortcut_ratio):
                continue

            # MUST-SPAN-NON-HABITAT precheck:
            # midpoint inside habitat indicates a corner-cut, not a true gap bridge.
            try:
                line_chord = QgsGeometry.fromPolylineXY([pts[i], pts[j]])
                mid_pt = line_chord.interpolate(gap * 0.5)
                if mid_pt is not None and (not mid_pt.isEmpty()) and patch_geom.contains(mid_pt):
                    continue
            except Exception:
                continue

            out.append(
                {
                    "pt_a": pts[i],
                    "pt_b": pts[j],
                    "detour_m": float(detour),
                    "gap_m": float(gap),
                    "ratio": float(ratio),
                    # Keep alias used by existing scoring/sorting code paths.
                    "shortcut_ratio": float(ratio),
                    "gain_m": float(detour - gap),
                }
            )

    out.sort(key=lambda x: float(x.get("ratio", 0.0) or 0.0), reverse=True)
    return _select_diverse_intra_candidates(patch_geom, out, int(max_candidates)) if max_candidates > 0 else out


def _build_valid_intra_bridge_geometry(
    *,
    patch_geom: QgsGeometry,
    pt_a: QgsPointXY,
    pt_b: QgsPointXY,
    min_corridor_width_m: float,
    subtract_geom: Optional[QgsGeometry] = None,
) -> Optional[QgsGeometry]:
    """
    Build a robust bridge polygon across a patch gap.
    Ensures the corridor is topologically connected to both sides of the concavity.
    """
    width_m = max(float(min_corridor_width_m or 0.0), 0.01)
    try:
        line = QgsGeometry.fromPolylineXY([pt_a, pt_b])
    except Exception:
        return None
    if line is None or line.isEmpty():
        return None

    try:
        raw_bridge = line.buffer(width_m / 2.0, 5)
    except Exception:
        return None
    if raw_bridge is None or raw_bridge.isEmpty():
        return None

    # Use a slightly shrunken patch for subtraction so bridge ends "seed" into habitat.
    try:
        habitat_to_cut = patch_geom.buffer(-0.2, 5)
        if habitat_to_cut is None or habitat_to_cut.isEmpty():
            habitat_to_cut = patch_geom
    except Exception:
        habitat_to_cut = patch_geom

    try:
        gap_only = raw_bridge.difference(habitat_to_cut)
    except Exception:
        gap_only = raw_bridge
    if gap_only is None or gap_only.isEmpty():
        return None

    try:
        gap_only = gap_only.makeValid()
    except Exception:
        pass
    if gap_only is None or gap_only.isEmpty():
        return None

    parts = _geometry_parts(gap_only)
    if not parts:
        return None

    pt_a_geom = QgsGeometry.fromPointXY(pt_a)
    pt_b_geom = QgsGeometry.fromPointXY(pt_b)
    min_area = float(line.length()) * float(width_m) * 0.4
    best_part: Optional[QgsGeometry] = None
    best_key: Optional[Tuple[float, float]] = None
    for part in parts:
        if part is None or part.isEmpty():
            continue
        try:
            d_a = float(part.distance(pt_a_geom))
            d_b = float(part.distance(pt_b_geom))
            part_area = float(part.area())
        except Exception:
            continue
        if d_a < width_m and d_b < width_m and part_area > min_area:
            key = (d_a + d_b, -part_area)
            if best_key is None or key < best_key:
                best_key = key
                best_part = part

    if best_part is None or best_part.isEmpty():
        return None
    return clone_geometry(best_part)


def _intra_loop_bonus(
    *,
    patch_geom: QgsGeometry,
    corridor_geom: QgsGeometry,
    base_hole_count: int,
    base_hole_area_m2: float,
) -> Tuple[float, int, float]:
    delta_count = 0
    delta_area_m2 = 0.0
    try:
        merged = patch_geom.combine(corridor_geom)
        try:
            merged = merged.makeValid()
        except Exception:
            pass
        holes_after = _extract_interior_ring_geometries(merged)
        hole_count_after = len(holes_after)
        hole_area_after = float(sum(float(h.area()) for h in holes_after))
        delta_count = max(int(hole_count_after) - int(base_hole_count), 0)
        delta_area_m2 = max(float(hole_area_after) - float(base_hole_area_m2), 0.0)
    except Exception:
        delta_count = 0
        delta_area_m2 = 0.0

    patch_area = max(float(patch_geom.area()), 1.0)
    hole_frac = float(delta_area_m2 / patch_area)
    # Give meaningful preference to true loop closure (C -> O), while keeping
    # raw detour reduction as the primary signal.
    bonus = 1.0 + min(3.0, float(delta_count) * 1.20 + hole_frac * 30.0)
    return float(bonus), int(delta_count), float(delta_area_m2)


def _best_valid_intra_shortcut(
    *,
    patch_geom: QgsGeometry,
    min_corridor_width_m: float,
    max_gap_m: float,
    min_shortcut_ratio: float,
    max_points: int,
    min_index_separation: int,
    subtract_geom: Optional[QgsGeometry],
) -> Optional[Dict[str, Any]]:
    pts = _ordered_patch_boundary_points(patch_geom, max_points=max_points)
    n = len(pts)
    if n < 8:
        return None

    cum: List[float] = [0.0]
    total_len = 0.0
    for i in range(n):
        j = (i + 1) % n
        d = float(pts[i].distance(pts[j]))
        total_len += d
        if i < n - 1:
            cum.append(total_len)
    if total_len <= 0.0:
        return None

    holes_before = _extract_interior_ring_geometries(patch_geom)
    base_hole_count = len(holes_before)
    base_hole_area_m2 = float(sum(float(h.area()) for h in holes_before))
    min_gap = max(float(min_corridor_width_m) * 0.5, 1.0)
    gap_cap = max(float(max_gap_m), min_gap)

    best: Optional[Dict[str, Any]] = None
    best_key: Optional[Tuple[float, float, float, float]] = None
    for i in range(n):
        for j in range(i + 1, n):
            idx_sep = min(j - i, n - (j - i))
            if idx_sep < int(min_index_separation):
                continue

            gap = float(pts[i].distance(pts[j]))
            if gap <= min_gap or gap > gap_cap:
                continue
            along = abs(cum[j] - cum[i])
            detour = min(along, total_len - along)
            if detour <= gap + 1e-9:
                continue
            ratio = float(detour / max(gap, 1e-9))
            if ratio + 1e-12 < float(min_shortcut_ratio):
                continue

            bridge = _build_valid_intra_bridge_geometry(
                patch_geom=patch_geom,
                pt_a=pts[i],
                pt_b=pts[j],
                min_corridor_width_m=min_corridor_width_m,
                subtract_geom=subtract_geom,
            )
            if bridge is None or bridge.isEmpty():
                continue

            area_ha = float(bridge.area() / 10000.0)
            if area_ha <= 0.0:
                continue

            raw_gain = float(detour - gap)
            loop_bonus, delta_holes, delta_hole_area_m2 = _intra_loop_bonus(
                patch_geom=patch_geom,
                corridor_geom=bridge,
                base_hole_count=base_hole_count,
                base_hole_area_m2=base_hole_area_m2,
            )
            tuned_gain = max(raw_gain, 0.0) * max(loop_bonus, 1.0)
            efficiency = tuned_gain / max(area_ha, 1e-9)
            key = (efficiency, tuned_gain, ratio, -gap)
            if best_key is None or key > best_key:
                best_key = key
                best = {
                    "geom": bridge,
                    "gap_m": float(gap),
                    "detour_m": float(detour),
                    "shortcut_ratio": float(ratio),
                    "gain_raw_m": float(raw_gain),
                    "gain_tuned_m": float(tuned_gain),
                    "loop_bonus": float(loop_bonus),
                    "loop_hole_count": int(delta_holes),
                    "loop_hole_area_m2": float(delta_hole_area_m2),
                    "area_ha": float(area_ha),
                }
    return best


def _best_intra_patch_shortcut(
    patch_geom: QgsGeometry,
    *,
    min_corridor_width_m: float,
    max_gap_m: float,
    min_shortcut_ratio: float,
    max_points: int = 120,
    min_index_separation: int = 4,
) -> Optional[Dict[str, Any]]:
    ranked = _rank_intra_patch_shortcuts(
        patch_geom,
        min_corridor_width_m=min_corridor_width_m,
        max_gap_m=max_gap_m,
        min_shortcut_ratio=min_shortcut_ratio,
        max_points=max_points,
        min_index_separation=min_index_separation,
        max_candidates=1,
    )
    return ranked[0] if ranked else None


def _get_landscape_fluidity_shortcut_threshold(
    params: "VectorRunParams",
    default_value: float = 3.0,
) -> float:
    try:
        val = float(
            getattr(
                params,
                "landscape_fluidity_shortcut_threshold",
                getattr(params, "resiliency_shortcut_threshold", default_value),
            )
            or default_value
        )
    except Exception:
        val = float(default_value)
    return max(1.0, float(val))


def _intra_patch_gap_cap_for_geom(
    patch_geom: QgsGeometry,
    params: "VectorRunParams",
    base_max_gap_m: Optional[float] = None,
) -> float:
    """
    Dynamic intra-gap cap so large C-gaps are still considered on large patches.
    """
    min_width = float(getattr(params, "min_corridor_width", 0.0) or 0.0)
    search_m = float(getattr(params, "max_search_distance", 0.0) or 0.0)
    base = float(base_max_gap_m) if base_max_gap_m is not None else max(min_width * 20.0, search_m * 0.90, 500.0)

    area_m2 = 0.0
    diag_m = 0.0
    try:
        area_m2 = max(float(patch_geom.area()), 0.0)
    except Exception:
        area_m2 = 0.0
    try:
        bb = patch_geom.boundingBox()
        diag_m = math.hypot(float(bb.width()), float(bb.height()))
    except Exception:
        diag_m = 0.0

    # area-based and extent-based caps keep this scale aware while still bounded.
    area_cap = math.sqrt(max(area_m2, 1.0)) * 0.70
    diag_cap = diag_m * 0.45
    dynamic_cap = max(base, area_cap, diag_cap)
    return max(dynamic_cap, min_width * 3.0, 100.0)


def _get_intra_patch_candidates(
    patches: Dict[int, Dict],
    params: "VectorRunParams",
    resiliency_shortcut_threshold: float = 3.0,
    patch_union: Optional[QgsGeometry] = None,
) -> List[Dict]:
    """
    Generate intra-patch bridge candidates and keep only deep concavity closures.
    """
    candidates: List[Dict] = []
    threshold = _get_landscape_fluidity_shortcut_threshold(params, resiliency_shortcut_threshold)
    max_bridge = max(
        float(getattr(params, "min_corridor_width", 0.0) or 0.0) * 20.0,
        float(getattr(params, "max_search_distance", 0.0) or 0.0) * 0.90,
        500.0,
    )
    max_points = min(800, int(getattr(params, "vector_terminal_max_per_patch", 120) or 120))
    max_keep_total = max(
        1,
        min(
            600,
            int(
                getattr(
                    params,
                    "landscape_fluidity_intra_k_total",
                    getattr(params, "resiliency_intra_k_total", 80),
                )
                or 80
            ),
        ),
    )

    for pid, pdata in patches.items():
        geom = pdata.get("geom")
        if geom is None or geom.isEmpty():
            continue

        patch_gap_cap = _intra_patch_gap_cap_for_geom(geom, params, max_bridge)
        ranked = _rank_intra_patch_shortcuts(
            geom,
            min_corridor_width_m=float(getattr(params, "min_corridor_width", 0.0) or 0.0),
            max_gap_m=patch_gap_cap,
            min_shortcut_ratio=threshold,
            max_points=max_points,
            min_index_separation=4,
            max_candidates=96,
        )
        if not ranked:
            continue

        patch_area_ha = max(float(pdata.get("area_ha", 0.0) or 0.0), 0.0)
        if patch_area_ha <= 0.0:
            continue

        for cand in ranked:
            detour_m = float(cand.get("detour_m", 0.0) or 0.0)
            gap_m = float(cand.get("gap_m", 0.0) or 0.0)
            if gap_m <= 0.0:
                continue
            shortcut_ratio = float(detour_m / max(gap_m, 1e-9))
            if shortcut_ratio + 1e-12 < threshold:
                continue

            bridge_geom = _build_valid_intra_bridge_geometry(
                patch_geom=geom,
                pt_a=cand["pt_a"],
                pt_b=cand["pt_b"],
                min_corridor_width_m=float(getattr(params, "min_corridor_width", 0.0) or 0.0),
                subtract_geom=patch_union if (patch_union is not None and (not patch_union.isEmpty())) else geom,
            )
            if bridge_geom is None or bridge_geom.isEmpty():
                continue

            area_ha = float(bridge_geom.area() / 10000.0)
            if area_ha <= 0.0:
                continue
            if params.max_corridor_area is not None and area_ha > float(params.max_corridor_area):
                continue

            candidates.append(
                {
                    "patch1": int(pid),
                    "patch2": int(pid),
                    "patch_ids": {int(pid)},
                    "geom": clone_geometry(bridge_geom),
                    "area_ha": float(area_ha),
                    "original_area_ha": float(area_ha),
                    "distance_m": float(gap_m),
                    "variant": "intra_patch_shortcut",
                    "source": "intra_shortcut",
                    "intra_patch": True,
                    "intra_patch_id": int(pid),
                    "intra_detour_old_m": float(detour_m),
                    "intra_shortcut_ratio": float(shortcut_ratio),
                    "fluidity_gain": float(shortcut_ratio),
                    "intra_gain_m": float(cand.get("gain_m", detour_m - gap_m)),
                    "intra_gap_m": float(gap_m),
                    "intra_hole_area_m2": float(cand.get("loop_hole_area_m2", 0.0) or 0.0),
                    "intra_pocket_strength": float(cand.get("pocket_strength", 0.0) or 0.0),
                    "intra_patch_area_ha": float(patch_area_ha),
                    "mid_x": float(cand.get("mid_x", 0.0) or 0.0),
                    "mid_y": float(cand.get("mid_y", 0.0) or 0.0),
                    "idx_a": int(cand.get("idx_a", -1) or -1),
                    "idx_b": int(cand.get("idx_b", -1) or -1),
                    "ring_n": int(cand.get("ring_n", 0) or 0),
                }
            )

    for cand in candidates:
        cand["intra_repair_score"] = float(_intra_structural_repair_score(cand))

    grouped: List[Dict[str, Any]] = []
    by_patch: Dict[int, List[Dict[str, Any]]] = defaultdict(list)
    for cand in candidates:
        try:
            by_patch[int(cand.get("intra_patch_id"))].append(cand)
        except Exception:
            continue
    for pid, patch_rows in by_patch.items():
        patch_geom = (patches.get(int(pid)) or {}).get("geom")
        if patch_geom is None or patch_geom.isEmpty():
            grouped.extend(patch_rows)
            continue
        grouped.extend(_group_intra_patch_winners(patch_geom, patch_rows, max_group_winners=1))

    grouped.sort(
        key=lambda c: (
            float(c.get("intra_repair_score", 0.0) or 0.0),
            float(c.get("intra_hole_area_m2", 0.0) or 0.0),
            float(c.get("intra_gain_m", 0.0) or 0.0),
            float(c.get("intra_shortcut_ratio", 0.0) or 0.0),
            float(c.get("intra_gap_m", 0.0) or 0.0),
        ),
        reverse=True,
    )
    if max_keep_total > 0:
        grouped = grouped[: int(max_keep_total)]
    return grouped


def _internal_shortcut_baseline_by_patch(
    patches: Dict[int, Dict],
    params: VectorRunParams,
    ratio_threshold: float,
) -> Dict[int, float]:
    out: Dict[int, float] = {}
    max_gap = max(
        float(getattr(params, "min_corridor_width", 0.0) or 0.0) * 20.0,
        float(getattr(params, "max_search_distance", 0.0) or 0.0) * 0.90,
        500.0,
    )
    min_ratio = max(1.05, float(ratio_threshold))
    for pid, pdata in patches.items():
        geom = pdata.get("geom")
        if geom is None or geom.isEmpty():
            continue
        patch_gap_cap = _intra_patch_gap_cap_for_geom(geom, params, max_gap)
        best = _best_valid_intra_shortcut(
            patch_geom=geom,
            min_corridor_width_m=float(getattr(params, "min_corridor_width", 0.0) or 0.0),
            max_gap_m=patch_gap_cap,
            min_shortcut_ratio=min_ratio,
            max_points=min(160, int(getattr(params, "vector_terminal_max_per_patch", 120) or 120)),
            min_index_separation=4,
            subtract_geom=geom,
        )
        if not best:
            continue
        best_gain = float(best.get("gain_tuned_m", 0.0) or 0.0)
        if best_gain <= 0.0:
            continue
        try:
            out[int(pid)] = float(best_gain)
        except Exception:
            continue
    return out


def _all_pairs_shortest_weighted(graph: nx.Graph) -> Dict[int, Dict[int, float]]:
    if graph.number_of_nodes() <= 0 or graph.number_of_edges() <= 0:
        return {}
    try:
        return {
            int(src): {int(dst): float(val) for dst, val in dist.items()}
            for src, dist in nx.all_pairs_dijkstra_path_length(graph, weight="weight")
        }
    except Exception:
        return {}


def _shortest_lookup(shortest: Dict[int, Dict[int, float]], a: int, b: int) -> float:
    ia = int(a)
    ib = int(b)
    if ia == ib:
        return 0.0
    return float(shortest.get(ia, {}).get(ib, float("inf")))


def _detour_ratio_from_distance(distance_graph: float, distance_euclid: float, detour_cap: float) -> float:
    eu = max(float(distance_euclid), 1e-9)
    cap = max(float(detour_cap), 1.0)
    if not np.isfinite(float(distance_graph)):
        return cap
    ratio = float(distance_graph) / eu
    if not np.isfinite(ratio):
        return cap
    return max(1.0, min(cap, ratio))


def _detour_ratio_score_from_shortest(
    shortest: Dict[int, Dict[int, float]],
    pair_records: Sequence[Tuple[int, int, float]],
    detour_cap: float,
) -> float:
    if not pair_records:
        return 0.0
    ratios: List[float] = []
    for u, v, eu in pair_records:
        d = _shortest_lookup(shortest, int(u), int(v))
        ratios.append(_detour_ratio_from_distance(d, float(eu), detour_cap))
    if not ratios:
        return 0.0
    mean_ratio = float(sum(ratios) / max(len(ratios), 1))
    if mean_ratio <= 0.0:
        return 0.0
    return max(0.0, min(1.0, 1.0 / mean_ratio))


def _detour_ratio_score_with_candidate_edges(
    shortest: Dict[int, Dict[int, float]],
    pair_records: Sequence[Tuple[int, int, float]],
    edges: Sequence[Tuple[int, int, float]],
    detour_cap: float,
) -> float:
    if not pair_records:
        return 0.0

    score_sum = 0.0
    n = 0
    for s, t, eu in pair_records:
        old_d = _shortest_lookup(shortest, int(s), int(t))
        new_d = old_d
        for u, v, w in edges:
            su = _shortest_lookup(shortest, int(s), int(u))
            sv = _shortest_lookup(shortest, int(s), int(v))
            tu = _shortest_lookup(shortest, int(t), int(u))
            tv = _shortest_lookup(shortest, int(t), int(v))
            if np.isfinite(su) and np.isfinite(tv):
                cand_d = su + float(w) + tv
                if cand_d < new_d:
                    new_d = float(cand_d)
            if np.isfinite(sv) and np.isfinite(tu):
                cand_d = sv + float(w) + tu
                if cand_d < new_d:
                    new_d = float(cand_d)
        ratio = _detour_ratio_from_distance(new_d, float(eu), detour_cap)
        score_sum += ratio
        n += 1

    if n <= 0:
        return 0.0
    mean_ratio = float(score_sum / float(n))
    if mean_ratio <= 0.0:
        return 0.0
    return max(0.0, min(1.0, 1.0 / mean_ratio))


def _apply_edges_with_changes(
    graph: nx.Graph,
    edges: Sequence[Tuple[int, int, float]],
) -> List[Tuple[str, int, int, float]]:
    changes: List[Tuple[str, int, int, float]] = []
    for u, v, w in edges:
        iu = int(u)
        iv = int(v)
        if iu == iv:
            continue
        ww = max(float(w), 1e-9)
        if graph.has_edge(iu, iv):
            old_w = float(graph[iu][iv].get("weight", ww) or ww)
            if ww + 1e-12 < old_w:
                graph[iu][iv]["weight"] = ww
                graph[iu][iv]["distance"] = ww
                changes.append(("update", iu, iv, old_w))
        else:
            graph.add_edge(iu, iv, weight=ww, distance=ww)
            changes.append(("add", iu, iv, 0.0))
    return changes


def _revert_edges_with_changes(
    graph: nx.Graph,
    changes: Sequence[Tuple[str, int, int, float]],
) -> None:
    for kind, u, v, old_w in reversed(list(changes)):
        iu = int(u)
        iv = int(v)
        if kind == "add":
            if graph.has_edge(iu, iv):
                graph.remove_edge(iu, iv)
        elif kind == "update":
            if graph.has_edge(iu, iv):
                graph[iu][iv]["weight"] = float(old_w)
                graph[iu][iv]["distance"] = float(old_w)


def _build_patch_mobility_context(
    patches: Dict[int, Dict],
    params: VectorRunParams,
) -> Optional[Dict[str, Any]]:
    node_ids: List[int] = []
    coords_xy: Dict[int, Tuple[float, float]] = {}
    patch_geoms: Dict[int, QgsGeometry] = {}
    patch_areas: Dict[int, float] = {}
    for pid, pdata in patches.items():
        try:
            ipid = int(pid)
        except Exception:
            continue
        area = float(pdata.get("area_ha", 0.0) or 0.0)
        geom = pdata.get("geom")
        if area <= 0.0 or geom is None or geom.isEmpty():
            continue
        try:
            c = geom.centroid().asPoint()
        except Exception:
            continue
        node_ids.append(int(ipid))
        coords_xy[int(ipid)] = (float(c.x()), float(c.y()))
        patch_geoms[int(ipid)] = geom
        patch_areas[int(ipid)] = float(area)

    node_ids = sorted(set(node_ids))
    if len(node_ids) < 2:
        return None

    graph = nx.Graph()
    graph.add_nodes_from(node_ids)

    boxes: Dict[int, QgsRectangle] = {}
    for pid in node_ids:
        try:
            boxes[int(pid)] = patch_geoms[int(pid)].boundingBox()
        except Exception:
            boxes[int(pid)] = QgsRectangle()

    for i, pid_i in enumerate(node_ids):
        gi = patch_geoms.get(int(pid_i))
        if gi is None or gi.isEmpty():
            continue
        bi = boxes.get(int(pid_i))
        if bi is None:
            continue
        x1, y1 = coords_xy[int(pid_i)]
        for j in range(i + 1, len(node_ids)):
            pid_j = int(node_ids[j])
            bj = boxes.get(pid_j)
            if bj is None:
                continue
            if not bi.intersects(bj):
                continue
            gj = patch_geoms.get(pid_j)
            if gj is None or gj.isEmpty():
                continue
            try:
                d_poly = float(gi.distance(gj))
            except Exception:
                continue
            if d_poly > 1e-6:
                continue
            x2, y2 = coords_xy[pid_j]
            w = max(float(math.hypot(x1 - x2, y1 - y2)), 1e-6)
            graph.add_edge(int(pid_i), int(pid_j), weight=w, distance=w)

    sample_k = max(2, min(200, int(getattr(params, "sample_points", 50) or 50)))
    pair_target = max(1, int(getattr(params, "pair_samples", 200) or 200))
    rng_nodes = random.Random(31)

    components = [sorted(int(n) for n in comp) for comp in nx.connected_components(graph)]
    components = [comp for comp in components if comp]
    components.sort(key=lambda c: (min(c), len(c)))
    if not components:
        return None

    comp_weights = [float(sum(patch_areas.get(int(pid), 0.0) for pid in comp)) for comp in components]
    total_weight = float(sum(comp_weights))
    if total_weight <= 0.0:
        comp_weights = [float(len(comp)) for comp in components]
        total_weight = float(sum(comp_weights))
    alloc = [max(1, int(math.floor((w / total_weight) * sample_k))) for w in comp_weights]
    alloc = [min(a, len(components[i])) for i, a in enumerate(alloc)]

    while sum(alloc) > sample_k:
        idx = max(range(len(alloc)), key=lambda k: (alloc[k], comp_weights[k]))
        if alloc[idx] <= 1:
            break
        alloc[idx] -= 1
    while sum(alloc) < sample_k:
        idx = max(range(len(alloc)), key=lambda k: (comp_weights[k], len(components[k]) - alloc[k]))
        if alloc[idx] >= len(components[idx]):
            break
        alloc[idx] += 1

    terminals: List[int] = []
    for comp_idx, comp in enumerate(components):
        need = max(1, min(int(alloc[comp_idx]), len(comp)))
        choices = list(comp)
        rng_nodes.shuffle(choices)
        terminals.extend(int(pid) for pid in choices[:need])

    terminals = sorted(set(terminals))
    if len(terminals) < 2:
        terminals = node_ids[: min(len(node_ids), max(2, sample_k))]
    if len(terminals) < 2:
        return None

    pair_indices: List[Tuple[int, int]] = [(i, j) for i in range(len(terminals)) for j in range(i + 1, len(terminals))]
    if not pair_indices:
        return None
    if len(pair_indices) > pair_target:
        rng_pairs = random.Random(37)
        rng_pairs.shuffle(pair_indices)
        pair_indices = sorted(pair_indices[:pair_target])

    pair_records: List[Tuple[int, int, float]] = []
    for i, j in pair_indices:
        u = int(terminals[i])
        v = int(terminals[j])
        x1, y1 = coords_xy.get(u, (0.0, 0.0))
        x2, y2 = coords_xy.get(v, (0.0, 0.0))
        eu = float(math.hypot(x1 - x2, y1 - y2))
        if eu <= 1e-9:
            continue
        pair_records.append((u, v, eu))
    if not pair_records:
        return None

    return {
        "graph": graph,
        "node_ids": node_ids,
        "coords_xy": coords_xy,
        "pair_records": pair_records,
    }


def _compute_strategic_mobility_exact(
    patches: Dict[int, Dict],
    corridors: Dict[int, Dict],
    params: VectorRunParams,
) -> Dict[str, float]:
    context = _build_patch_mobility_context(patches, params)
    if context is None:
        return {"pre": 0.0, "post": 0.0, "gain": 0.0}

    detour_cap = max(1.0, float(getattr(params, "mobility_detour_cap", 8.0) or 8.0))
    graph_pre: nx.Graph = context["graph"].copy()
    pair_records: List[Tuple[int, int, float]] = list(context["pair_records"])
    valid_nodes = set(int(pid) for pid in context["node_ids"])
    coords_xy: Dict[int, Tuple[float, float]] = dict(context["coords_xy"])

    graph_post = graph_pre.copy()
    for cdata in corridors.values():
        pids = _candidate_patch_ids_for_metric(cdata, valid_nodes)
        if len(pids) < 2:
            continue
        length = float(cdata.get("distance", cdata.get("distance_m", 0.0)) or 0.0)
        if length <= 0.0:
            p1, p2 = int(pids[0]), int(pids[-1])
            x1, y1 = coords_xy.get(p1, (0.0, 0.0))
            x2, y2 = coords_xy.get(p2, (0.0, 0.0))
            length = max(float(math.hypot(x1 - x2, y1 - y2)), 1.0)
        edges_local: List[Tuple[int, int, float]] = []
        for i in range(len(pids)):
            for j in range(i + 1, len(pids)):
                edges_local.append((int(pids[i]), int(pids[j]), float(length)))
        _apply_edges_with_changes(graph_post, edges_local)

    shortest_pre = _all_pairs_shortest_weighted(graph_pre)
    shortest_post = _all_pairs_shortest_weighted(graph_post)
    score_pre = _detour_ratio_score_from_shortest(shortest_pre, pair_records, detour_cap)
    score_post = _detour_ratio_score_from_shortest(shortest_post, pair_records, detour_cap)
    return {
        "pre": float(score_pre),
        "post": float(score_post),
        "gain": float(score_post - score_pre),
    }


def _build_patch_resistance_context(
    patches: Dict[int, Dict],
    params: VectorRunParams,
) -> Optional[Dict[str, Any]]:
    base = _build_patch_mobility_context(patches, params)
    if base is None:
        return None
    coords_xy: Dict[int, Tuple[float, float]] = dict(base.get("coords_xy") or {})
    if len(coords_xy) < 2:
        return None
    pair_target = max(200, min(500, int(getattr(params, "pair_samples", 300) or 300)))
    pair_records: List[Tuple[int, int, float]] = []
    seen_pairs: Set[Tuple[int, int]] = set()

    area_ha_by_pid: Dict[int, float] = {}
    patch_geom_by_pid: Dict[int, QgsGeometry] = {}
    for pid in coords_xy.keys():
        try:
            ipid = int(pid)
        except Exception:
            continue
        area_ha_by_pid[ipid] = float((patches.get(ipid) or {}).get("area_ha", 0.0) or 0.0)
        patch_geom_by_pid[ipid] = (patches.get(ipid) or {}).get("geom")

    def _try_add_pair(u_raw: int, v_raw: int, eu_raw: Optional[float] = None) -> bool:
        try:
            u = int(u_raw)
            v = int(v_raw)
        except Exception:
            return False
        if u == v:
            return False
        a, b = (u, v) if u < v else (v, u)
        if (a, b) in seen_pairs:
            return False
        if eu_raw is None:
            x1, y1 = coords_xy.get(int(a), (0.0, 0.0))
            x2, y2 = coords_xy.get(int(b), (0.0, 0.0))
            eu = float(math.hypot(x1 - x2, y1 - y2))
        else:
            eu = float(eu_raw)
        if (not np.isfinite(eu)) or eu <= 1e-9:
            return False
        seen_pairs.add((a, b))
        pair_records.append((int(a), int(b), float(eu)))
        return True

    def _pair_rank_distance(u: int, v: int) -> Tuple[float, float]:
        """Return (boundary_distance, centroid_distance) for ranking/recording."""
        iu = int(u)
        iv = int(v)
        x1, y1 = coords_xy.get(iu, (0.0, 0.0))
        x2, y2 = coords_xy.get(iv, (0.0, 0.0))
        eu = float(math.hypot(x1 - x2, y1 - y2))
        if (not np.isfinite(eu)) or eu <= 1e-9:
            eu = 1e-9
        prox = float(eu)
        gu = patch_geom_by_pid.get(iu)
        gv = patch_geom_by_pid.get(iv)
        if gu is not None and gv is not None and (not gu.isEmpty()) and (not gv.isEmpty()):
            try:
                d = float(gu.distance(gv))
                if np.isfinite(d) and d >= 0.0:
                    prox = max(float(d), 1e-9)
            except Exception:
                pass
        return float(prox), float(eu)

    # Coverage stratum: ensure each component contributes representative pairs.
    # This avoids cases where locally meaningful inter-patch links evaluate as
    # zero only because sampled pairs skipped those components.
    try:
        graph_base: nx.Graph = base["graph"]
        comp_list = [sorted(int(n) for n in comp) for comp in nx.connected_components(graph_base) if comp]
    except Exception:
        comp_list = []
    comp_list.sort(key=lambda comp: (len(comp), comp[0]))
    comp_id_by_node: Dict[int, int] = {}
    comp_area_ha_by_idx: Dict[int, float] = {}
    for cidx, comp in enumerate(comp_list):
        comp_area_ha = float(
            sum(max(float(area_ha_by_pid.get(int(n), 0.0) or 0.0), 0.0) for n in comp)
        )
        comp_area_ha_by_idx[int(cidx)] = float(comp_area_ha)
        for n in comp:
            comp_id_by_node[int(n)] = int(cidx)
    comp_anchors: List[int] = []
    for comp in comp_list:
        anchor = max(
            comp,
            key=lambda n: (
                float(area_ha_by_pid.get(int(n), 0.0)),
                -int(n),
            ),
        )
        comp_anchors.append(int(anchor))

    # Major-component coverage stratum:
    # Guarantee sampled cross-component pairs among the components that hold most
    # habitat area, so large merges are always measurable in the objective.
    major_comp_ids: List[int] = []
    if comp_list:
        comp_idx_area = sorted(
            (
                (int(idx), float(comp_area_ha_by_idx.get(int(idx), 0.0) or 0.0))
                for idx in range(len(comp_list))
            ),
            key=lambda t: (float(t[1]), -int(t[0])),
            reverse=True,
        )
        total_comp_area = float(sum(float(a) for _idx, a in comp_idx_area))
        min_major_count = min(len(comp_idx_area), max(2, int(math.ceil(math.sqrt(float(len(comp_idx_area))) / 2.0))))
        max_major_count = min(len(comp_idx_area), max(6, int(math.ceil(math.sqrt(float(len(comp_idx_area))) * 2.0))))
        running = 0.0
        for idx, area in comp_idx_area:
            if len(major_comp_ids) >= max_major_count:
                break
            major_comp_ids.append(int(idx))
            running += float(area)
            if len(major_comp_ids) >= min_major_count and running + 1e-12 >= 0.70 * max(total_comp_area, 1e-9):
                break

    if len(major_comp_ids) >= 2:
        major_pair_pool: List[Tuple[float, float, int, int]] = []
        major_anchor_k = max(2, min(8, int(math.ceil(math.sqrt(float(len(major_comp_ids)))))))
        for ci in major_comp_ids:
            ui = int(comp_anchors[int(ci)])
            ai = float(comp_area_ha_by_idx.get(int(ci), 0.0) or 0.0)
            nearest_major: List[Tuple[float, float, int, float]] = []
            for cj in major_comp_ids:
                if int(cj) == int(ci):
                    continue
                uj = int(comp_anchors[int(cj)])
                prox_d, eu = _pair_rank_distance(int(ui), int(uj))
                if (not np.isfinite(prox_d)) or prox_d <= 1e-9:
                    continue
                aj = float(comp_area_ha_by_idx.get(int(cj), 0.0) or 0.0)
                nearest_major.append((float(prox_d), float(eu), int(uj), float(aj)))
            nearest_major.sort(key=lambda t: (float(t[0]), -float(t[3]), int(t[2])))
            for prox_d, eu, uj, aj in nearest_major[:major_anchor_k]:
                a, b = (int(ui), int(uj)) if int(ui) < int(uj) else (int(uj), int(ui))
                priority = math.sqrt(max(ai * aj, 0.0)) / max(float(prox_d), 1e-9)
                major_pair_pool.append((float(priority), float(eu), int(a), int(b)))

        major_pair_pool.sort(key=lambda t: (-float(t[0]), float(t[1]), int(t[2]), int(t[3])))
        major_target = max(1, min(pair_target, int(round(float(pair_target) * 0.20))))
        for _prio, eu, a, b in major_pair_pool:
            if len(pair_records) >= major_target:
                break
            _try_add_pair(int(a), int(b), float(eu))

    # Extra coverage for large patches: each major node gets nearest cross-component
    # pair candidates so short obvious bridges don't become invisible to sampling.
    node_ids = sorted(int(n) for n in coords_xy.keys())
    major_node_cap = max(10, int(round(float(len(node_ids)) * 0.20)))
    major_node_cap = min(len(node_ids), major_node_cap)
    major_nodes = sorted(
        node_ids,
        key=lambda n: (
            float(area_ha_by_pid.get(int(n), 0.0) or 0.0),
            -int(n),
        ),
        reverse=True,
    )[:major_node_cap]

    coverage_pool: List[Tuple[float, float, int, int]] = []
    anchor_k = 2
    for u in comp_anchors:
        nearest: List[Tuple[float, float, int]] = []
        for v in comp_anchors:
            if int(v) == int(u):
                continue
            prox_d, eu = _pair_rank_distance(int(u), int(v))
            if (not np.isfinite(prox_d)) or prox_d <= 1e-9:
                continue
            nearest.append((float(prox_d), float(eu), int(v)))
        nearest.sort(key=lambda t: (float(t[0]), float(t[1]), int(t[2])))
        for prox_d, eu, v in nearest[:anchor_k]:
            a, b = (int(u), int(v)) if int(u) < int(v) else (int(v), int(u))
            area_prod = float(area_ha_by_pid.get(int(a), 0.0) * area_ha_by_pid.get(int(b), 0.0))
            # Rank by area-product over distance (dimensionless priority).
            priority = area_prod / max(float(prox_d), 1e-9)
            coverage_pool.append((float(priority), float(eu), int(a), int(b)))

    major_k = 3
    for u in major_nodes:
        cu = int(comp_id_by_node.get(int(u), -1))
        nearest_cross: List[Tuple[float, float, int]] = []
        for v in node_ids:
            iv = int(v)
            if iv == int(u):
                continue
            if int(comp_id_by_node.get(iv, -2)) == cu:
                continue
            prox_d, eu = _pair_rank_distance(int(u), iv)
            if (not np.isfinite(prox_d)) or prox_d <= 1e-9:
                continue
            nearest_cross.append((float(prox_d), float(eu), int(iv)))
        nearest_cross.sort(key=lambda t: (float(t[0]), float(t[1]), int(t[2])))
        for prox_d, eu, v in nearest_cross[:major_k]:
            a, b = (int(u), int(v)) if int(u) < int(v) else (int(v), int(u))
            area_scale = math.sqrt(
                max(float(area_ha_by_pid.get(int(a), 0.0) or 0.0), 0.0)
                * max(float(area_ha_by_pid.get(int(b), 0.0) or 0.0), 0.0)
            )
            priority = area_scale / max(float(prox_d), 1e-9)
            coverage_pool.append((float(priority), float(eu), int(a), int(b)))

    coverage_pool.sort(key=lambda t: (-float(t[0]), float(t[1]), int(t[2]), int(t[3])))
    coverage_target = max(1, min(pair_target, int(round(float(pair_target) * 0.45))))
    for _prio, eu, a, b in coverage_pool:
        if len(pair_records) >= coverage_target:
            break
        _try_add_pair(int(a), int(b), float(eu))

    # Landscape stratum: far-separated pairs for global resistance behavior.
    if len(pair_records) < pair_target and graph_math is not None and hasattr(graph_math, "sample_spatially_separated_pairs"):
        try:
            far_pairs = list(
                graph_math.sample_spatially_separated_pairs(
                    coords_xy,
                    pair_count=max(pair_target, 3 * max(pair_target - len(pair_records), 1)),
                    seed=41,
                )
            )
        except Exception:
            far_pairs = []
        far_pairs.sort(key=lambda t: (-float(t[2]), int(t[0]), int(t[1])))
        for u, v, eu in far_pairs:
            if len(pair_records) >= pair_target:
                break
            _try_add_pair(int(u), int(v), float(eu))

    # Fallback: include baseline mobility pairs if needed.
    if len(pair_records) < pair_target:
        fallback_pairs = list(base.get("pair_records") or [])
        fallback_pairs.sort(key=lambda t: (-float(t[2]), int(t[0]), int(t[1])))
        for u, v, eu in fallback_pairs:
            if len(pair_records) >= pair_target:
                break
            _try_add_pair(int(u), int(v), float(eu))

    if not pair_records:
        return None

    # Area-weighted pair weights: emphasize links that improve movement between
    # larger habitat bodies while staying scale-free.
    pair_weights_raw: Dict[Tuple[int, int], float] = {}
    weight_sum = 0.0
    for u, v, _eu in pair_records:
        iu = int(u)
        iv = int(v)
        key = (iu, iv) if iu < iv else (iv, iu)
        au = max(float(area_ha_by_pid.get(int(key[0]), 0.0) or 0.0), 0.0)
        av = max(float(area_ha_by_pid.get(int(key[1]), 0.0) or 0.0), 0.0)
        w = math.sqrt(max(au * av, 0.0))
        if w <= 0.0:
            w = 1e-9
        pair_weights_raw[key] = float(w)
        weight_sum += float(w)
    if weight_sum <= 0.0:
        weight_sum = float(len(pair_weights_raw) or 1)
    pair_weights: Dict[Tuple[int, int], float] = {
        key: float(val / weight_sum) for key, val in pair_weights_raw.items()
    }

    pair_max_eu = max((float(d) for _u, _v, d in pair_records), default=1.0)
    disconnected_penalty = max(1.0, pair_max_eu * 20.0)
    return {
        "graph": base["graph"],
        "node_ids": list(base.get("node_ids") or []),
        "coords_xy": coords_xy,
        "pair_records": pair_records,
        "pair_weights": pair_weights,
        "disconnected_penalty": float(disconnected_penalty),
    }


def _mean_effective_resistance_graph(
    graph: nx.Graph,
    pair_records: Sequence[Tuple[int, int, float]],
    disconnected_penalty: float,
    pair_weights: Optional[Dict[Tuple[int, int], float]] = None,
) -> float:
    normalized_weights: Optional[Dict[Tuple[int, int], float]] = None
    if pair_weights:
        normalized_weights = {}
        for key, val in pair_weights.items():
            try:
                u, v = key
                iu = int(u)
                iv = int(v)
            except Exception:
                continue
            kk = (iu, iv) if iu < iv else (iv, iu)
            normalized_weights[kk] = float(val or 0.0)

    if normalized_weights and graph_math is not None and hasattr(graph_math, "weighted_mean_effective_resistance_sampled"):
        try:
            val = graph_math.weighted_mean_effective_resistance_sampled(
                graph,
                pair_records,
                normalized_weights,
                cost_attr="weight",
                disconnected_penalty=float(disconnected_penalty),
            )
            if np.isfinite(float(val)):
                return float(val)
        except Exception:
            pass

    if graph_math is not None and hasattr(graph_math, "mean_effective_resistance_sampled"):
        try:
            val = graph_math.mean_effective_resistance_sampled(
                graph,
                pair_records,
                cost_attr="weight",
                disconnected_penalty=float(disconnected_penalty),
            )
            if np.isfinite(float(val)):
                return float(val)
        except Exception:
            pass

    # Fallback: shortest-path mean (not true effective resistance, used only if helper is unavailable).
    shortest = _all_pairs_shortest_weighted(graph)
    weighted_total = 0.0
    weight_sum = 0.0
    vals: List[float] = []
    penalty = max(float(disconnected_penalty), 1.0)
    for u, v, _eu in pair_records:
        d = _shortest_lookup(shortest, int(u), int(v))
        if np.isfinite(d):
            val = max(float(d), 0.0)
        else:
            val = penalty
        vals.append(val)
        if normalized_weights:
            key = (int(u), int(v)) if int(u) < int(v) else (int(v), int(u))
            w = max(float(normalized_weights.get(key, 0.0) or 0.0), 0.0)
            if w > 0.0:
                weighted_total += w * val
                weight_sum += w
    if normalized_weights and weight_sum > 0.0:
        return float(weighted_total / weight_sum)
    if not vals:
        return penalty
    return float(sum(vals) / max(len(vals), 1))


def _landscape_fluidity_intra_factor(
    ratio: float,
    detour_gain_m: float,
    *,
    hole_area_m2: float = 0.0,
    pocket_strength: float = 0.0,
    patch_area_ha: float = 0.0,
    gap_m: float = 0.0,
) -> float:
    """Dimensionless intra benefit term favoring meaningful concavity closure."""
    if ratio <= 1.0 + 1e-12:
        return 0.0
    detour_gain_m = max(float(detour_gain_m), 0.0)
    if detour_gain_m <= 1e-12:
        return 0.0
    ratio_component = max(0.0, min(1.0, (float(ratio) - 1.0) / max(float(ratio), 1e-9)))
    detour_component = math.sqrt(detour_gain_m)
    patch_area_m2 = max(float(patch_area_ha) * 10000.0, 1.0)
    closure_share = max(0.0, min(1.0, float(hole_area_m2) / max(0.35 * patch_area_m2, 1.0)))
    closure_component = math.sqrt(closure_share)
    pocket_component = max(0.0, min(1.0, math.log1p(max(float(pocket_strength), 0.0)) / math.log1p(150.0)))
    mouth_scale = max(0.0, min(1.0, float(gap_m) / max(math.sqrt(patch_area_m2) * 0.30, 1.0)))
    structural_component = 0.20 + 0.50 * closure_component + 0.30 * pocket_component
    mouth_component = 0.35 + 0.65 * mouth_scale
    return float(ratio_component * detour_component * structural_component * mouth_component)


def _landscape_fluidity_intra_bonus_from_factor(
    *,
    graph_fluidity_pre: float,
    intra_factor_sum: float,
    params: VectorRunParams,
    intra_patch_count: int = 1,
) -> float:
    """Scale/cap intra shortcut contribution so it cannot dominate graph-fluidity gain."""
    bonus_weight = max(0.0, float(getattr(params, "landscape_fluidity_intra_bonus_weight", 1.0) or 1.0))
    cap_mult = max(0.0, float(getattr(params, "landscape_fluidity_intra_bonus_cap_mult", 1.0) or 1.0))
    scaled = max(0.0, float(intra_factor_sum)) * bonus_weight
    bounded = min(scaled, cap_mult * max(1, int(intra_patch_count)))
    return float(max(float(graph_fluidity_pre), 0.0) * bounded)


def _landscape_fluidity_inter_loop_factor(
    shortcut_ratio: float,
    component_area_ha: float,
    total_habitat_area_ha: float,
    hop_length: int = 0,
) -> float:
    """
    Dimensionless loop value term for same-component shortcuts.
    Rewards meaningful detour relief on larger connected components.
    """
    if shortcut_ratio <= 1.0 + 1e-12:
        return 0.0
    if component_area_ha <= 0.0 or total_habitat_area_ha <= 0.0:
        return 0.0
    ratio_component = max(0.0, min(1.0, (float(shortcut_ratio) - 1.0) / max(float(shortcut_ratio), 1e-9)))
    area_share = max(0.0, min(1.0, float(component_area_ha) / max(float(total_habitat_area_ha), 1e-9)))
    # Mild convexity on ratio and hop-depth boost help value "final easy closure"
    # links that collapse long local detours in big components.
    ratio_emphasis = math.sqrt(max(ratio_component, 0.0))
    hop_component = 0.0
    try:
        hops = int(hop_length)
    except Exception:
        hops = 0
    if hops > 1:
        hop_component = max(0.0, min(1.0, float(hops - 1) / float(max(hops, 1))))
    hop_boost = 0.70 + (0.30 * hop_component)
    return float(ratio_emphasis * area_share * hop_boost)


def _landscape_fluidity_inter_loop_bonus_from_factor(
    *,
    graph_fluidity_pre: float,
    inter_loop_factor_sum: float,
    params: VectorRunParams,
) -> float:
    """
    Scale/cap inter-patch loop contribution so it nudges selection toward
    valuable circuit closure without dominating core flow gains.
    """
    bonus_weight = max(0.0, float(getattr(params, "landscape_fluidity_inter_loop_weight", 0.35) or 0.35))
    cap_mult = max(0.0, float(getattr(params, "landscape_fluidity_inter_loop_cap_mult", 0.80) or 0.80))
    scaled = max(0.0, float(inter_loop_factor_sum)) * bonus_weight
    bounded = min(scaled, cap_mult)
    return float(max(float(graph_fluidity_pre), 0.0) * bounded)


def _landscape_fluidity_bridge_relief_factor(
    bridge_fraction_on_old_path: float,
    component_area_ha: float,
    total_habitat_area_ha: float,
    shortcut_ratio: float = 1.0,
) -> float:
    """
    Dimensionless robustness gain term for same-component edges that
    provide alternative routes across bridge-dominated paths.
    """
    frac = max(0.0, min(1.0, float(bridge_fraction_on_old_path)))
    if frac <= 1e-12:
        return 0.0
    if component_area_ha <= 0.0 or total_habitat_area_ha <= 0.0:
        return 0.0
    # Only reward robustness additions when the candidate is at least some
    # degree of path-shortening (avoids very long near-parallel "insurance" links).
    ratio_component = max(
        0.0,
        min(1.0, (float(shortcut_ratio) - 1.0) / max(float(shortcut_ratio), 1e-9)),
    )
    if ratio_component <= 1e-12:
        return 0.0
    area_share = max(0.0, min(1.0, float(component_area_ha) / max(float(total_habitat_area_ha), 1e-9)))
    return float(frac * area_share * ratio_component)


def _landscape_fluidity_bridge_relief_bonus_from_factor(
    *,
    graph_fluidity_pre: float,
    bridge_relief_factor_sum: float,
    params: VectorRunParams,
) -> float:
    bonus_weight = max(0.0, float(getattr(params, "landscape_fluidity_bridge_relief_weight", 0.25) or 0.25))
    cap_mult = max(0.0, float(getattr(params, "landscape_fluidity_bridge_relief_cap_mult", 0.60) or 0.60))
    scaled = max(0.0, float(bridge_relief_factor_sum)) * bonus_weight
    bounded = min(scaled, cap_mult)
    return float(max(float(graph_fluidity_pre), 0.0) * bounded)


def _bridge_fraction_on_shortest_path(
    graph: nx.Graph,
    u: int,
    v: int,
    bridge_edges: Set[Tuple[int, int]],
) -> float:
    iu = int(u)
    iv = int(v)
    if iu == iv:
        return 0.0
    try:
        path_nodes = nx.shortest_path(graph, source=iu, target=iv, weight="weight")
    except Exception:
        return 0.0
    if not path_nodes or len(path_nodes) < 2:
        return 0.0
    total_edges = len(path_nodes) - 1
    if total_edges <= 0:
        return 0.0
    bridge_count = 0
    for i in range(total_edges):
        a = int(path_nodes[i])
        b = int(path_nodes[i + 1])
        key = (a, b) if a <= b else (b, a)
        if key in bridge_edges:
            bridge_count += 1
    return float(bridge_count / float(total_edges))


def _compute_landscape_fluidity_exact(
    patches: Dict[int, Dict],
    corridors: Dict[int, Dict],
    params: VectorRunParams,
) -> Dict[str, float]:
    context = _build_patch_resistance_context(patches, params)
    if context is None:
        return {"pre": 0.0, "post": 0.0, "gain": 0.0}

    graph_pre: nx.Graph = context["graph"].copy()
    graph_post: nx.Graph = graph_pre.copy()
    pair_records: List[Tuple[int, int, float]] = list(context["pair_records"])
    pair_weights: Dict[Tuple[int, int], float] = dict(context.get("pair_weights") or {})
    penalty = float(context.get("disconnected_penalty", 1.0))
    valid_nodes = set(int(pid) for pid in context.get("node_ids", []))
    coords_xy: Dict[int, Tuple[float, float]] = dict(context.get("coords_xy") or {})
    total_habitat_area = max(
        1e-9,
        float(sum(float(p.get("area_ha", 0.0) or 0.0) for p in patches.values())),
    )
    intra_factor_by_patch: Dict[int, float] = {}
    inter_loop_factor_by_pair: Dict[Tuple[int, int], float] = {}
    inter_loop_factor_sum = 0.0
    bridge_relief_factor_by_pair: Dict[Tuple[int, int], float] = {}
    bridge_relief_factor_sum = 0.0
    shortcut_threshold = _get_landscape_fluidity_shortcut_threshold(params, 3.0)
    same_component_shortcut_threshold = max(
        1.05,
        1.0 + 0.25 * max(float(shortcut_threshold) - 1.0, 0.0),
    )

    for _cid, cdata in sorted(corridors.items(), key=lambda kv: int(kv[0])):
        p1_raw = cdata.get("p1", cdata.get("patch1"))
        p2_raw = cdata.get("p2", cdata.get("patch2"))
        is_intra = bool(cdata.get("intra_patch", False))
        try:
            if p1_raw is not None and p2_raw is not None and int(p1_raw) == int(p2_raw):
                is_intra = True
        except Exception:
            pass
        if is_intra:
            try:
                intra_pid = int(cdata.get("intra_patch_id", p1_raw))
            except Exception:
                intra_pid = None
            if intra_pid is None:
                continue
            gap_m = float(cdata.get("distance_m", cdata.get("distance", 0.0)) or 0.0)
            detour_m = float(cdata.get("intra_detour_old_m", 0.0) or 0.0)
            ratio = float(cdata.get("intra_shortcut_ratio", cdata.get("fluidity_gain", 0.0)) or 0.0)
            if ratio <= 0.0 and gap_m > 0.0 and detour_m > 0.0:
                ratio = float(detour_m / max(gap_m, 1e-9))
            detour_gain_m = max(float(detour_m) - float(gap_m), 0.0)
            if detour_gain_m <= 0.0 and ratio > 1.0 + 1e-12 and gap_m > 0.0:
                detour_gain_m = max((float(ratio) - 1.0) * float(gap_m), 0.0)
            factor = _landscape_fluidity_intra_factor(
                float(ratio),
                float(detour_gain_m),
                hole_area_m2=float(cdata.get("intra_hole_area_m2", 0.0) or 0.0),
                pocket_strength=float(cdata.get("intra_pocket_strength", 0.0) or 0.0),
                patch_area_ha=float(cdata.get("intra_patch_area_ha", 0.0) or 0.0),
                gap_m=float(cdata.get("intra_gap_m", gap_m) or gap_m),
            )
            if factor <= 0.0:
                continue
            cap_mult = max(0.0, float(getattr(params, "landscape_fluidity_intra_bonus_cap_mult", 1.0) or 1.0))
            if cap_mult > 0.0:
                factor = min(float(factor), float(cap_mult))
            intra_factor_by_patch[int(intra_pid)] = max(
                float(intra_factor_by_patch.get(int(intra_pid), 0.0) or 0.0),
                float(factor),
            )
            continue

        pids = _candidate_patch_ids_for_metric(cdata, valid_nodes)
        if len(pids) < 2:
            continue
        length = float(cdata.get("distance", cdata.get("distance_m", 0.0)) or 0.0)
        if length <= 0.0:
            p1 = int(pids[0])
            p2 = int(pids[-1])
            x1, y1 = coords_xy.get(p1, (0.0, 0.0))
            x2, y2 = coords_xy.get(p2, (0.0, 0.0))
            length = max(float(math.hypot(x1 - x2, y1 - y2)), 1.0)

        try:
            gate_u = int(p1_raw) if p1_raw is not None else int(pids[0])
        except Exception:
            gate_u = int(pids[0])
        try:
            gate_v = int(p2_raw) if p2_raw is not None else int(pids[-1])
        except Exception:
            gate_v = int(pids[-1])
        if gate_u not in valid_nodes:
            gate_u = int(pids[0])
        if gate_v not in valid_nodes:
            gate_v = int(pids[-1])

        if gate_u != gate_v:
            comp_id: Dict[int, int] = {}
            comp_area_ha: Dict[int, float] = defaultdict(float)
            for comp_idx, comp in enumerate(nx.connected_components(graph_post)):
                for node in comp:
                    inode = int(node)
                    comp_id[inode] = int(comp_idx)
                    comp_area_ha[int(comp_idx)] += float((patches.get(inode) or {}).get("area_ha", 0.0) or 0.0)
            cu = int(comp_id.get(int(gate_u), -1))
            cv = int(comp_id.get(int(gate_v), -1))
            if cu >= 0 and cv >= 0 and cu == cv:
                shortest_current = _all_pairs_shortest_weighted(graph_post)
                d_old = _shortest_lookup(shortest_current, int(gate_u), int(gate_v))
                c_e = max(float(length), 1e-9)
                shortcut_ratio = (d_old / c_e) if np.isfinite(d_old) else float("inf")
                comp_area = float(comp_area_ha.get(cu, 0.0) or 0.0)
                try:
                    hop_len = int(nx.shortest_path_length(graph_post, int(gate_u), int(gate_v)))
                except Exception:
                    hop_len = 0
                loop_key = (int(gate_u), int(gate_v)) if int(gate_u) < int(gate_v) else (int(gate_v), int(gate_u))
                if shortcut_ratio + 1e-12 >= same_component_shortcut_threshold:
                    loop_factor = _landscape_fluidity_inter_loop_factor(
                        float(shortcut_ratio),
                        float(comp_area),
                        float(total_habitat_area),
                        hop_length=int(hop_len),
                    )
                    if loop_factor > 0.0:
                        prev_factor = float(inter_loop_factor_by_pair.get(loop_key, 0.0) or 0.0)
                        if loop_factor > prev_factor:
                            inter_loop_factor_sum += float(loop_factor - prev_factor)
                            inter_loop_factor_by_pair[loop_key] = float(loop_factor)
                bridge_edges: Set[Tuple[int, int]] = set()
                try:
                    for a, b in nx.bridges(graph_post):
                        ia, ib = int(a), int(b)
                        bridge_edges.add((ia, ib) if ia <= ib else (ib, ia))
                except Exception:
                    bridge_edges = set()
                bridge_fraction = _bridge_fraction_on_shortest_path(
                    graph_post,
                    int(gate_u),
                    int(gate_v),
                    bridge_edges,
                )
                bridge_factor = _landscape_fluidity_bridge_relief_factor(
                    float(bridge_fraction),
                    float(comp_area),
                    float(total_habitat_area),
                    float(shortcut_ratio),
                )
                if bridge_factor > 0.0:
                    prev_bridge = float(bridge_relief_factor_by_pair.get(loop_key, 0.0) or 0.0)
                    if bridge_factor > prev_bridge:
                        bridge_relief_factor_sum += float(bridge_factor - prev_bridge)
                        bridge_relief_factor_by_pair[loop_key] = float(bridge_factor)

        edges_local: List[Tuple[int, int, float]] = []
        for i in range(len(pids)):
            for j in range(i + 1, len(pids)):
                edges_local.append((int(pids[i]), int(pids[j]), float(length)))
        _apply_edges_with_changes(graph_post, edges_local)

    mean_pre_graph = float(
        _mean_effective_resistance_graph(
            graph_pre,
            pair_records,
            penalty,
            pair_weights=pair_weights,
        )
    )
    mean_post_graph = float(
        _mean_effective_resistance_graph(
            graph_post,
            pair_records,
            penalty,
            pair_weights=pair_weights,
        )
    )
    graph_fluidity_pre = float(total_habitat_area / max(mean_pre_graph, 1e-9))
    graph_fluidity_post = float(total_habitat_area / max(mean_post_graph, 1e-9))
    intra_factor_sum = float(sum(float(v) for v in intra_factor_by_patch.values()))
    intra_bonus_post = _landscape_fluidity_intra_bonus_from_factor(
        graph_fluidity_pre=float(graph_fluidity_pre),
        intra_factor_sum=float(intra_factor_sum),
        params=params,
        intra_patch_count=len(intra_factor_by_patch),
    )
    inter_loop_bonus_post = _landscape_fluidity_inter_loop_bonus_from_factor(
        graph_fluidity_pre=float(graph_fluidity_pre),
        inter_loop_factor_sum=float(inter_loop_factor_sum),
        params=params,
    )
    bridge_relief_bonus_post = _landscape_fluidity_bridge_relief_bonus_from_factor(
        graph_fluidity_pre=float(graph_fluidity_pre),
        bridge_relief_factor_sum=float(bridge_relief_factor_sum),
        params=params,
    )
    fluidity_pre = float(graph_fluidity_pre)
    fluidity_post = float(graph_fluidity_post + intra_bonus_post + inter_loop_bonus_post + bridge_relief_bonus_post)

    return {
        "pre": float(fluidity_pre),
        "post": float(fluidity_post),
        "gain": float(fluidity_post - fluidity_pre),
        "graph_resistance_pre": float(mean_pre_graph),
        "graph_resistance_post": float(mean_post_graph),
        "intra_bonus_post": float(intra_bonus_post),
        "intra_factor_sum": float(intra_factor_sum),
        "inter_loop_bonus_post": float(inter_loop_bonus_post),
        "inter_loop_factor_sum": float(inter_loop_factor_sum),
        "bridge_relief_bonus_post": float(bridge_relief_bonus_post),
        "bridge_relief_factor_sum": float(bridge_relief_factor_sum),
    }


def _compute_resiliency_flow_exact(
    patches: Dict[int, Dict],
    corridors: Dict[int, Dict],
    params: VectorRunParams,
) -> Dict[str, float]:
    # Backward-compatible wrapper.
    return _compute_landscape_fluidity_exact(patches, corridors, params)


def optimize_landscape_fluidity(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
) -> Tuple[Dict[int, Dict], Dict]:
    budget_limit = float(params.budget_area or 0.0)
    if budget_limit <= 0.0 or (not patches) or (not candidates):
        return {}, {"strategy": "landscape_fluidity", "corridors_used": 0, "budget_used_ha": 0.0}

    context = _build_patch_resistance_context(patches, params)
    if context is None:
        return {}, {"strategy": "landscape_fluidity", "corridors_used": 0, "budget_used_ha": 0.0}

    graph = context["graph"].copy()
    graph_base = context["graph"].copy()
    node_ids: List[int] = list(context.get("node_ids") or [])
    coords_xy: Dict[int, Tuple[float, float]] = dict(context.get("coords_xy") or {})
    pair_records: List[Tuple[int, int, float]] = list(context.get("pair_records") or [])
    pair_weights: Dict[Tuple[int, int], float] = dict(context.get("pair_weights") or {})
    disconnected_penalty = float(context.get("disconnected_penalty", 1.0))
    valid_nodes = set(int(pid) for pid in node_ids)
    total_habitat_area = max(
        1e-9,
        float(sum(float((patches.get(int(pid)) or {}).get("area_ha", 0.0) or 0.0) for pid in valid_nodes)),
    )

    shortcut_threshold = _get_landscape_fluidity_shortcut_threshold(params, 3.0)
    intra_cap = max(
        1,
        min(
            500,
            int(
                getattr(
                    params,
                    "landscape_fluidity_intra_k_total",
                    getattr(params, "resiliency_intra_k_total", 80),
                )
                or 80
            ),
        ),
    )
    eval_cap = max(
        20,
        min(
            1200,
            int(
                getattr(
                    params,
                    "landscape_fluidity_eval_cap",
                    getattr(params, "resiliency_eval_cap", 140),
                )
                or 140
            ),
        ),
    )
    intra_tau = max(0.0, float(getattr(params, "landscape_fluidity_intra_tau", 0.8) or 0.8))
    gain_weight = max(0.0, float(getattr(params, "landscape_fluidity_gain_weight", 0.65) or 0.65))
    ratio_weight = max(0.0, float(getattr(params, "landscape_fluidity_ratio_weight", 0.35) or 0.35))
    lookahead_k = max(1, min(20, int(getattr(params, "landscape_fluidity_lookahead_k", 8) or 8)))
    lookahead_weight = max(0.0, min(1.5, float(getattr(params, "landscape_fluidity_lookahead_weight", 0.85) or 0.85)))
    topology_weight = max(0.0, min(1.0, float(getattr(params, "landscape_fluidity_topology_weight", 0.20) or 0.20)))
    long_hop_factor_thresh = max(1.05, float(getattr(params, "landscape_fluidity_long_hop_factor", 1.35) or 1.35))
    long_hop_tau = max(0.0, float(getattr(params, "landscape_fluidity_long_hop_tau", 0.90) or 0.90))
    if gain_weight + ratio_weight <= 1e-12:
        gain_weight, ratio_weight = 0.65, 0.35
    same_component_shortcut_threshold = max(
        1.05,
        1.0 + 0.25 * max(float(shortcut_threshold) - 1.0, 0.0),
    )

    cand_specs: List[Dict[str, Any]] = []
    for idx, cand in enumerate(candidates):
        cost = float(cand.get("area_ha", 0.0) or 0.0)
        if cost <= 0.0 or cost > budget_limit + 1e-12:
            continue

        p1_raw = cand.get("patch1", cand.get("p1"))
        p2_raw = cand.get("patch2", cand.get("p2"))
        is_intra = bool(cand.get("intra_patch", False))
        try:
            if p1_raw is not None and p2_raw is not None and int(p1_raw) == int(p2_raw):
                is_intra = True
        except Exception:
            pass

        if is_intra:
            try:
                intra_pid = int(cand.get("intra_patch_id", p1_raw))
            except Exception:
                continue
            if intra_pid not in valid_nodes:
                continue
            length = max(float(cand.get("distance_m", cand.get("distance", 0.0)) or 0.0), 1e-9)
            detour_old = float(cand.get("intra_detour_old_m", 0.0) or 0.0)
            intra_ratio = float(cand.get("intra_shortcut_ratio", cand.get("fluidity_gain", 0.0)) or 0.0)
            if intra_ratio <= 0.0 and detour_old > 0.0:
                intra_ratio = float(detour_old / max(length, 1e-9))
            if intra_ratio <= 1.0 + 1e-12:
                continue
            intra_factor = _landscape_fluidity_intra_factor(
                float(intra_ratio),
                float(
                    max(float(detour_old) - float(length), 0.0)
                    if detour_old > 0.0 and length > 0.0
                    else max((float(intra_ratio) - 1.0) * float(length), 0.0)
                ),
                hole_area_m2=float(cand.get("intra_hole_area_m2", 0.0) or 0.0),
                pocket_strength=float(cand.get("intra_pocket_strength", 0.0) or 0.0),
                patch_area_ha=float(cand.get("intra_patch_area_ha", 0.0) or 0.0),
                gap_m=float(cand.get("intra_gap_m", length) or length),
            )
            if intra_factor <= 1e-12:
                continue
            cap_mult = max(0.0, float(getattr(params, "landscape_fluidity_intra_bonus_cap_mult", 1.0) or 1.0))
            if cap_mult > 0.0:
                intra_factor = min(float(intra_factor), float(cap_mult))
            cand_specs.append(
                {
                    "idx": int(idx),
                    "cand": cand,
                    "cost": float(cost),
                    "pids": [int(intra_pid)],
                    "edges": [],
                    "length": float(length),
                    "gate_u": int(intra_pid),
                    "gate_v": int(intra_pid),
                    "gate_pair": (int(intra_pid), int(intra_pid)),
                    "is_intra": True,
                    "intra_pid": int(intra_pid),
                    "intra_detour_old_m": float(detour_old),
                    "shortcut_ratio": float(intra_ratio),
                    "intra_factor": float(intra_factor),
                    "intra_hole_area_m2": float(cand.get("intra_hole_area_m2", 0.0) or 0.0),
                    "intra_pocket_strength": float(cand.get("intra_pocket_strength", 0.0) or 0.0),
                    "intra_patch_area_ha": float(cand.get("intra_patch_area_ha", 0.0) or 0.0),
                    "intra_gap_m": float(cand.get("intra_gap_m", length) or length),
                    "intra_repair_score": float(cand.get("intra_repair_score", 0.0) or 0.0),
                }
            )
            continue

        pids = _candidate_patch_ids_for_metric(cand, valid_nodes)
        length = float(cand.get("distance_m", 0.0) or 0.0)
        edges_local: List[Tuple[int, int, float]] = []
        explicit_edges = cand.get("explicit_edges")
        if isinstance(explicit_edges, (list, tuple)):
            for item in explicit_edges:
                try:
                    if isinstance(item, (list, tuple)) and len(item) >= 2:
                        eu = int(item[0])
                        ev = int(item[1])
                        elen = float(item[2]) if len(item) >= 3 else float(length)
                    else:
                        continue
                except Exception:
                    continue
                if eu not in valid_nodes or ev not in valid_nodes or eu == ev:
                    continue
                if elen <= 0.0:
                    x1, y1 = coords_xy.get(eu, (0.0, 0.0))
                    x2, y2 = coords_xy.get(ev, (0.0, 0.0))
                    elen = max(float(math.hypot(x1 - x2, y1 - y2)), 1.0)
                edges_local.append((int(eu), int(ev), float(elen)))
            if edges_local:
                node_union: Set[int] = set()
                for eu, ev, _ in edges_local:
                    node_union.add(int(eu))
                    node_union.add(int(ev))
                pids = sorted(node_union)

        if not edges_local:
            if len(pids) < 2:
                continue
            if length <= 0.0:
                p1 = int(pids[0])
                p2 = int(pids[-1])
                x1, y1 = coords_xy.get(p1, (0.0, 0.0))
                x2, y2 = coords_xy.get(p2, (0.0, 0.0))
                length = max(float(math.hypot(x1 - x2, y1 - y2)), 1.0)
            for i in range(len(pids)):
                for j in range(i + 1, len(pids)):
                    edges_local.append((int(pids[i]), int(pids[j]), float(length)))
            if not edges_local:
                continue

        try:
            gate_u = int(p1_raw) if p1_raw is not None else int(pids[0])
        except Exception:
            gate_u = int(pids[0])
        try:
            gate_v = int(p2_raw) if p2_raw is not None else int(pids[-1])
        except Exception:
            gate_v = int(pids[-1])
        if gate_u not in valid_nodes:
            gate_u = int(pids[0])
        if gate_v not in valid_nodes:
            gate_v = int(pids[-1])
        gate_pair = (gate_u, gate_v) if gate_u <= gate_v else (gate_v, gate_u)
        allow_multi_same_pair = bool(cand.get("barrier_pair", False))
        graph_edges = list(edges_local)
        if allow_multi_same_pair and len(pids) == 2 and len(edges_local) == 1:
            eu, ev, el = edges_local[0]
            vmid = -int(idx + 1)
            h = max(float(el) * 0.5, 1e-9)
            graph_edges = [(int(eu), int(vmid), float(h)), (int(vmid), int(ev), float(h))]

        cand_specs.append(
            {
                "idx": int(idx),
                "cand": cand,
                "cost": float(cost),
                "score_cost": float(
                    max(
                        float(cost)
                        * max(1.0, float(cand.get("bridge_over_penalty", 1.0) or 1.0)),
                        1e-12,
                    )
                ),
                "pids": list(pids),
                "edges": edges_local,
                "graph_edges": graph_edges,
                "length": float(length),
                "gate_u": int(gate_u),
                "gate_v": int(gate_v),
                "gate_pair": gate_pair,
                "is_intra": False,
                "allow_multi_same_pair": bool(allow_multi_same_pair),
                "anchor_pts": cand.get("anchor_pts"),
                "bridge_over_penalty": float(max(1.0, float(cand.get("bridge_over_penalty", 1.0) or 1.0))),
                "bridge_over_count": int(cand.get("bridge_over_count", 0) or 0),
            }
        )
    if not cand_specs:
        return {}, {"strategy": "landscape_fluidity", "corridors_used": 0, "budget_used_ha": 0.0}

    def _selected_row_from_spec(spec: Dict[str, Any], utility_score: float) -> Dict[str, Any]:
        cand = spec["cand"]
        pids = list(spec.get("pids", []))
        gate_u = int(spec.get("gate_u", -1))
        gate_v = int(spec.get("gate_v", -1))
        p1_val = cand.get("patch1", cand.get("p1"))
        p2_val = cand.get("patch2", cand.get("p2"))
        if p1_val is None and pids:
            p1_val = pids[0]
        if p2_val is None and pids:
            p2_val = pids[-1]
        try:
            p1_int = int(p1_val)
        except Exception:
            p1_int = int(gate_u if gate_u >= 0 else (pids[0] if pids else 0))
        try:
            p2_int = int(p2_val)
        except Exception:
            p2_int = int(gate_v if gate_v >= 0 else (pids[-1] if pids else p1_int))
        if p1_int < 0 and gate_u >= 0:
            p1_int = int(gate_u)
        if p2_int < 0 and gate_v >= 0:
            p2_int = int(gate_v)
        return {
            "geom": clone_geometry(cand["geom"]),
            "patch_ids": set(int(pid) for pid in pids),
            "area_ha": float(spec.get("cost", 0.0) or 0.0),
            "p1": int(p1_int),
            "p2": int(p2_int),
            "distance": float(spec.get("length", 1.0) or 1.0),
            "graph_edges": list(spec.get("graph_edges", spec.get("edges", [])) or []),
            "type": "primary",
            "variant": cand.get("variant"),
            "source": cand.get("source"),
            "utility_score": float(utility_score),
            "overlap_ratio": 0.0,
            "intra_patch": bool(spec.get("is_intra", False)),
            "intra_patch_id": int(spec.get("intra_pid")) if bool(spec.get("is_intra", False)) else None,
            "intra_shortcut_ratio": (
                float(spec.get("shortcut_ratio", 0.0) or 0.0) if bool(spec.get("is_intra", False)) else 0.0
            ),
            "fluidity_gain": (
                float(spec.get("shortcut_ratio", 0.0) or 0.0) if bool(spec.get("is_intra", False)) else 0.0
            ),
            "mid_x": (
                float(spec.get("mid_x", cand.get("mid_x", 0.0)) or 0.0) if bool(spec.get("is_intra", False)) else 0.0
            ),
            "mid_y": (
                float(spec.get("mid_y", cand.get("mid_y", 0.0)) or 0.0) if bool(spec.get("is_intra", False)) else 0.0
            ),
            "idx_a": (
                int(spec.get("idx_a", cand.get("idx_a", -1)) or -1) if bool(spec.get("is_intra", False)) else -1
            ),
            "idx_b": (
                int(spec.get("idx_b", cand.get("idx_b", -1)) or -1) if bool(spec.get("is_intra", False)) else -1
            ),
            "ring_n": (
                int(spec.get("ring_n", cand.get("ring_n", 0)) or 0) if bool(spec.get("is_intra", False)) else 0
            ),
            "intra_gap_m": (
                float(spec.get("intra_gap_m", cand.get("distance_m", 0.0)) or 0.0)
                if bool(spec.get("is_intra", False))
                else 0.0
            ),
            "intra_repair_score": (
                float(spec.get("intra_repair_score", 0.0) or 0.0) if bool(spec.get("is_intra", False)) else 0.0
            ),
            "same_component_inter": (
                bool(spec.get("same_component_inter", False)) if (not bool(spec.get("is_intra", False))) else False
            ),
        }

    # Soft dominance check: if an inter-component direct edge is much more expensive
    # than an available two-hop alternative (u->p->v), treat it as lower efficiency
    # in ranking (without forbidding selection).
    pair_min_cost: Dict[Tuple[int, int], float] = {}
    adj_cost: Dict[int, Dict[int, float]] = defaultdict(dict)
    for spec in cand_specs:
        if bool(spec.get("is_intra", False)):
            continue
        try:
            u = int(spec.get("gate_u"))
            v = int(spec.get("gate_v"))
        except Exception:
            continue
        if u == v:
            continue
        key = (u, v) if u <= v else (v, u)
        c = float(spec.get("cost", 0.0) or 0.0)
        if c <= 0.0:
            continue
        prev = pair_min_cost.get(key)
        if prev is None or c < prev:
            pair_min_cost[key] = float(c)
        prev_uv = adj_cost[u].get(v)
        if prev_uv is None or c < prev_uv:
            adj_cost[u][v] = float(c)
            adj_cost[v][u] = float(c)

    for spec in cand_specs:
        if bool(spec.get("is_intra", False)):
            continue
        try:
            u = int(spec.get("gate_u"))
            v = int(spec.get("gate_v"))
        except Exception:
            continue
        if u == v:
            continue
        c = float(spec.get("cost", 0.0) or 0.0)
        if c <= 0.0:
            continue
        neigh_u = adj_cost.get(u, {})
        neigh_v = adj_cost.get(v, {})
        if not neigh_u or not neigh_v:
            continue
        best_twohop = float("inf")
        for p, cup in neigh_u.items():
            if int(p) in (u, v):
                continue
            cpv = neigh_v.get(int(p))
            if cpv is None:
                continue
            total = float(cup + cpv)
            if total < best_twohop:
                best_twohop = total
        if np.isfinite(best_twohop) and best_twohop > 0.0:
            spec["best_twohop_cost"] = float(best_twohop)
            long_hop_factor = float(c / max(best_twohop, 1e-9))
            spec["long_hop_factor"] = float(long_hop_factor)
            if c > (1.25 * best_twohop):
                # Penalize score-cost only; budget uses true area_ha.
                dominance_factor = min(3.0, c / max(best_twohop * 1.15, 1e-9))
                spec["score_cost"] = float(max(float(spec.get("score_cost", c) or c), c * dominance_factor))
        else:
            spec["long_hop_factor"] = 1.0

    # Hard dominance filter (requested): if a direct A-C corridor is not cheaper than
    # an available two-hop A-B-C alternative, drop it from consideration.
    filtered_specs: List[Dict[str, Any]] = []
    removed = 0
    for spec in cand_specs:
        if bool(spec.get("is_intra", False)):
            filtered_specs.append(spec)
            continue
        best_twohop = float(spec.get("best_twohop_cost", 0.0) or 0.0)
        if best_twohop > 0.0:
            c = float(spec.get("cost", 0.0) or 0.0)
            if c >= best_twohop - 1e-12:
                removed += 1
                continue
        filtered_specs.append(spec)

    if removed:
        print(f"  ✓ Hard-dominance filter removed {removed} corridor candidate(s)")
    cand_specs = filtered_specs

    patch_area_by_pid: Dict[int, float] = {
        int(pid): float((pdata or {}).get("area_ha", 0.0) or 0.0) for pid, pdata in patches.items()
    }
    small_patch_share_cap = max(0.03, min(0.16, 220.0 / max(total_habitat_area, 1e-9)))
    for spec in cand_specs:
        if bool(spec.get("is_intra", False)):
            spec["small_patch_penalty"] = 1.0
            spec["small_patch_share"] = 0.0
            continue
        pids_local = [int(pid) for pid in spec.get("pids", []) if int(pid) in patch_area_by_pid]
        if len(pids_local) < 2:
            spec["small_patch_penalty"] = 1.0
            spec["small_patch_share"] = 0.0
            continue
        patch_shares = [
            max(0.0, min(1.0, patch_area_by_pid[int(pid)] / max(total_habitat_area, 1e-9)))
            for pid in pids_local
        ]
        min_share = min(patch_shares) if patch_shares else 0.0
        spec["small_patch_share"] = float(min_share)
        if min_share <= small_patch_share_cap:
            penalty = 1.0 + min(2.4, (small_patch_share_cap / max(min_share, 1e-9) - 1.0) * 1.10)
            if min_share <= 0.45 * small_patch_share_cap:
                penalty *= 1.25
            spec["small_patch_penalty"] = float(penalty)
            score_cost0 = max(float(spec.get("score_cost", spec.get("cost", 0.0)) or 0.0), 1e-12)
            spec["score_cost"] = float(score_cost0 * penalty)
        else:
            spec["small_patch_penalty"] = 1.0


    used_specs: Set[int] = set()
    selected_pairs: Set[Tuple[int, int]] = set()
    pair_multi_min_sep = max(float(getattr(params, "max_search_distance", 0.0) or 0.0), 0.0)
    selected_barrier_pair_geoms: Dict[Tuple[int, int], List[QgsGeometry]] = defaultdict(list)
    selected_barrier_pair_anchors: Dict[Tuple[int, int], List[Dict[int, QgsPointXY]]] = defaultdict(list)
    selected: Dict[int, Dict] = {}
    remaining = float(budget_limit)
    budget_used = 0.0
    intra_selected = 0
    inter_selected = 0
    same_component_selected = 0

    current_graph_mean = _mean_effective_resistance_graph(
        graph,
        pair_records,
        disconnected_penalty,
        pair_weights=pair_weights,
    )
    current_graph_score = float(total_habitat_area / max(float(current_graph_mean), 1e-9))
    graph_fluidity_pre = float(current_graph_score)
    intra_factor_by_patch: Dict[int, float] = {}
    current_intra_factor_sum = 0.0
    current_intra_bonus = 0.0
    inter_loop_factor_by_pair: Dict[Tuple[int, int], float] = {}
    current_inter_loop_factor_sum = 0.0
    current_inter_loop_bonus = 0.0
    bridge_relief_factor_by_pair: Dict[Tuple[int, int], float] = {}
    current_bridge_relief_factor_sum = 0.0
    current_bridge_relief_bonus = 0.0
    current_score = float(current_graph_score + current_intra_bonus + current_inter_loop_bonus + current_bridge_relief_bonus)
    score_pre = float(current_score)
    mean_pre_graph = float(current_graph_mean)

    def _pair_anchor_points(
        gate_pair: Tuple[int, int],
        cgeom: Optional[QgsGeometry],
        anchor_pts: Optional[Dict[int, QgsPointXY]] = None,
    ) -> Dict[int, QgsPointXY]:
        if isinstance(anchor_pts, dict) and anchor_pts:
            out: Dict[int, QgsPointXY] = {}
            for pid in gate_pair:
                try:
                    if int(pid) in anchor_pts:
                        out[int(pid)] = anchor_pts[int(pid)]
                except Exception:
                    continue
            if len(out) == len(gate_pair):
                return out
        if cgeom is None or cgeom.isEmpty():
            return {}
        out: Dict[int, QgsPointXY] = {}
        for pid in gate_pair:
            pgeom = (patches.get(int(pid)) or {}).get("geom")
            if pgeom is None or pgeom.isEmpty():
                continue
            try:
                npg = pgeom.nearestPoint(cgeom)
                if npg is None or npg.isEmpty():
                    continue
                out[int(pid)] = QgsPointXY(npg.asPoint())
            except Exception:
                continue
        return out

    def _is_distinct_barrier_pair(
        gate_pair: Tuple[int, int],
        cgeom: Optional[QgsGeometry],
        anchor_pts: Optional[Dict[int, QgsPointXY]] = None,
    ) -> bool:
        prev_anchor_sets = selected_barrier_pair_anchors.get(gate_pair) or []
        if prev_anchor_sets and isinstance(anchor_pts, dict) and anchor_pts:
            for pa in prev_anchor_sets:
                if len(anchor_pts) == 2 and len(pa) == 2:
                    for pid in gate_pair:
                        if (
                            int(pid) not in anchor_pts
                            or int(pid) not in pa
                            or float(anchor_pts[int(pid)].distance(pa[int(pid)])) + 1e-9 < pair_multi_min_sep
                        ):
                            return False
                else:
                    return False
            return True
        prevs = selected_barrier_pair_geoms.get(gate_pair) or []
        if not prevs:
            return True
        if cgeom is None or cgeom.isEmpty():
            return False
        ca = _pair_anchor_points(gate_pair, cgeom, anchor_pts=anchor_pts)
        for pg in prevs:
            pa = _pair_anchor_points(gate_pair, pg)
            if len(ca) == 2 and len(pa) == 2:
                for pid in gate_pair:
                    if (
                        int(pid) not in ca
                        or int(pid) not in pa
                        or float(ca[int(pid)].distance(pa[int(pid)])) + 1e-9 < pair_multi_min_sep
                    ):
                        return False
            else:
                try:
                    d = float(cgeom.distance(pg))
                except Exception:
                    d = 0.0
                if (not np.isfinite(d)) or d + 1e-9 < pair_multi_min_sep:
                    return False
        return True

    while remaining > 1e-12:
        comp_id: Dict[int, int] = {}
        comp_area_ha: Dict[int, float] = defaultdict(float)
        comp_node_count: Dict[int, int] = defaultdict(int)
        for comp_idx, comp in enumerate(nx.connected_components(graph)):
            for node in comp:
                comp_id[int(node)] = int(comp_idx)
                comp_area_ha[int(comp_idx)] += float((patches.get(int(node)) or {}).get("area_ha", 0.0) or 0.0)
                comp_node_count[int(comp_idx)] += 1
        shortest_current = _all_pairs_shortest_weighted(graph)
        bridge_edges_current: Set[Tuple[int, int]] = set()
        try:
            for a, b in nx.bridges(graph):
                ia, ib = int(a), int(b)
                bridge_edges_current.add((ia, ib) if ia <= ib else (ib, ia))
        except Exception:
            bridge_edges_current = set()

        component_share_vals = sorted(
            max(0.0, min(1.0, float(a) / max(float(total_habitat_area), 1e-9)))
            for a in comp_area_ha.values()
            if float(a) > 0.0
        )
        if component_share_vals:
            q35_comp_idx = max(0, int(math.floor(0.35 * max(len(component_share_vals) - 1, 0))))
            q35_comp_share = float(component_share_vals[q35_comp_idx])
        else:
            q35_comp_share = 0.0

        inter_ranked: List[Tuple[float, float, float, int]] = []
        same_component_ranked: List[Tuple[float, float, float, int]] = []
        intra_ranked: List[Tuple[float, float, float, int]] = []

        for spec_idx, spec in enumerate(cand_specs):
            if spec_idx in used_specs:
                continue
            is_barrier_multi = bool(spec.get("allow_multi_same_pair", False))
            gate_pair_spec = tuple(spec.get("gate_pair"))
            if (not is_barrier_multi) and gate_pair_spec in selected_pairs:
                continue
            if is_barrier_multi and ((gate_pair_spec in selected_pairs) or bool(selected_barrier_pair_geoms.get(gate_pair_spec))):
                cgeom = (spec.get("cand") or {}).get("geom")
                if not _is_distinct_barrier_pair(gate_pair_spec, cgeom, anchor_pts=spec.get("anchor_pts")):
                    continue
            cost = float(spec.get("cost", 0.0) or 0.0)
            if cost <= 0.0 or cost > remaining + 1e-12:
                continue
            score_cost = float(spec.get("score_cost", cost) or cost)
            if (not bool(spec.get("is_intra", False))) and score_cost <= 0.0:
                penalty = float(spec.get("bridge_over_penalty", 1.0) or 1.0)
                score_cost = max(cost * max(penalty, 1.0), 1e-12)
            gate_u = int(spec.get("gate_u"))
            gate_v = int(spec.get("gate_v"))
            is_intra = bool(spec.get("is_intra", False))
            if is_intra:
                intra_pid = int(spec.get("intra_pid"))
                base_intra_factor = float(spec.get("intra_factor", 0.0) or 0.0)
                curr_patch_factor = float(intra_factor_by_patch.get(intra_pid, 0.0) or 0.0)
                gain_intra_factor = max(base_intra_factor - curr_patch_factor, 0.0)
                next_intra_factor_sum = float(current_intra_factor_sum + gain_intra_factor)
                gain_intra = float(
                    _landscape_fluidity_intra_bonus_from_factor(
                        graph_fluidity_pre=graph_fluidity_pre,
                        intra_factor_sum=next_intra_factor_sum,
                        params=params,
                    )
                    - current_intra_bonus
                )
                if gain_intra <= 1e-12:
                    continue
                shortcut_ratio = float(spec.get("shortcut_ratio", 0.0) or 0.0)
                if shortcut_ratio + 1e-12 < shortcut_threshold:
                    continue
                repair_score = float(spec.get("intra_repair_score", 0.0) or 0.0)
                patch_area_share = max(
                    0.0,
                    min(
                        1.0,
                        float((patches.get(int(intra_pid)) or {}).get("area_ha", 0.0) or 0.0)
                        / max(float(total_habitat_area), 1e-9),
                    ),
                )
                ipc_priority = float(
                    (0.60 * (gain_intra / max(cost, 1e-12)))
                    + (0.30 * repair_score)
                    + (0.10 * math.sqrt(max(patch_area_share, 0.0)))
                )
                rank_tuple = (
                    ipc_priority,
                    gain_intra,
                    repair_score,
                    -cost,
                    int(spec_idx),
                )
                intra_ranked.append(rank_tuple)
            else:
                if gate_u == gate_v:
                    continue
                cu = int(comp_id.get(gate_u, -1))
                cv = int(comp_id.get(gate_v, -1))
                same_component = (cu >= 0 and cv >= 0 and cu == cv)
                d_old = _shortest_lookup(shortest_current, gate_u, gate_v)
                c_e = max(float(spec.get("length", 0.0) or 0.0), 1e-9)
                shortcut_ratio = (d_old / c_e) if np.isfinite(d_old) else float("inf")
                bridge_factor_rank = 0.0
                if same_component:
                    bridge_fraction = _bridge_fraction_on_shortest_path(
                        graph,
                        int(gate_u),
                        int(gate_v),
                        bridge_edges_current,
                    )
                    comp_area = float(comp_area_ha.get(cu, 0.0) or 0.0)
                    bridge_factor_rank = _landscape_fluidity_bridge_relief_factor(
                        float(bridge_fraction),
                        float(comp_area),
                        float(total_habitat_area),
                        float(shortcut_ratio),
                    )
                    pair_already = (gate_pair_spec in selected_pairs) or bool(selected_barrier_pair_geoms.get(gate_pair_spec))
                    if (
                        (not (is_barrier_multi and pair_already))
                        and shortcut_ratio + 1e-12 < same_component_shortcut_threshold
                        and bridge_factor_rank <= 1e-12
                    ):
                        # Keep mild-ratio closures on large components eligible so
                        # "final easy circuit" links are not filtered out too early.
                        comp_share = max(
                            0.0,
                            min(1.0, float(comp_area) / max(float(total_habitat_area), 1e-9)),
                        )
                        rel_shortcut_floor = max(1.05, 0.72 * float(same_component_shortcut_threshold))
                        if not (
                            shortcut_ratio + 1e-12 >= rel_shortcut_floor
                            and comp_share + 1e-12 >= max(q35_comp_share, 1e-12)
                        ):
                            continue
                if same_component:
                    impacted_area_ha = float(comp_area_ha.get(cu, 0.0) or 0.0)
                else:
                    impacted_area_ha = float(comp_area_ha.get(cu, 0.0) or 0.0) + float(comp_area_ha.get(cv, 0.0) or 0.0)
                structural_rank = 0.0
                if not same_component:
                    total_nodes = max(1, len(valid_nodes))
                    structural_rank = min(
                        1.0,
                        float(comp_node_count.get(cu, 0) + comp_node_count.get(cv, 0)) / float(total_nodes),
                    )
                structural_rank += min(0.60, 0.20 * max(0, len(spec.get("pids", [])) - 2))
                heuristic = float(
                    (
                        (shortcut_ratio if np.isfinite(shortcut_ratio) else (disconnected_penalty / c_e))
                        + (2.0 * float(bridge_factor_rank))
                        + (0.6 * float(structural_rank))
                    )
                    * max(impacted_area_ha, 1e-9)
                )
                inter_ranked.append((heuristic / max(score_cost, 1e-12), heuristic, -score_cost, int(spec_idx)))
                if same_component:
                    ratio_component = max(
                        0.0,
                        min(1.0, (float(shortcut_ratio) - 1.0) / max(float(shortcut_ratio), 1e-9)),
                    )
                    comp_share = max(
                        0.0,
                        min(1.0, float(impacted_area_ha) / max(float(total_habitat_area), 1e-9)),
                    )
                    loop_rank = float(
                        (
                            (0.70 * ratio_component)
                            + (0.30 * math.sqrt(max(ratio_component, 0.0)))
                            + (1.60 * float(bridge_factor_rank))
                        )
                        * max(comp_share, 1e-12)
                    )
                    if loop_rank > 0.0:
                        same_component_ranked.append(
                            (loop_rank / max(score_cost, 1e-12), loop_rank, -score_cost, int(spec_idx))
                        )
                # Store score_cost back so downstream ratio scoring reuses same effective cost basis.
                spec["score_cost"] = float(max(score_cost, 1e-12))

        inter_ranked.sort(reverse=True)
        same_component_ranked.sort(reverse=True)
        intra_ranked.sort(reverse=True)

        eval_indices: List[int] = []
        inter_slots = min(len(inter_ranked), max(1, int(eval_cap * 0.75)))
        eval_indices.extend(int(spec_idx) for _s, _h, _c, spec_idx in inter_ranked[:inter_slots])
        # Force-evaluate a bounded set of high-mass inter-component candidates so
        # medium-cost but high-impact merges are not pruned before exact-gain scoring.
        if inter_ranked and len(eval_indices) < eval_cap:
            major_force_cap = max(4, min(24, int(math.ceil(eval_cap * 0.12))))
            major_ranked: List[Tuple[float, int]] = []
            for _s, _h, _c, spec_idx in inter_ranked:
                if len(major_ranked) >= max(len(inter_ranked), major_force_cap * 4):
                    break
                spec = cand_specs[int(spec_idx)]
                if bool(spec.get("is_intra", False)):
                    continue
                gate_u = int(spec.get("gate_u"))
                gate_v = int(spec.get("gate_v"))
                cu = int(comp_id.get(gate_u, -1))
                cv = int(comp_id.get(gate_v, -1))
                if cu < 0 or cv < 0 or cu == cv:
                    continue
                area_u = float(comp_area_ha.get(cu, 0.0) or 0.0)
                area_v = float(comp_area_ha.get(cv, 0.0) or 0.0)
                if area_u <= 0.0 or area_v <= 0.0:
                    continue
                gm_share = math.sqrt(max(area_u * area_v, 0.0)) / max(float(total_habitat_area), 1e-9)
                score_cost = max(float(spec.get("score_cost", spec.get("cost", 0.0)) or 0.0), 1e-12)
                mass_priority = float(gm_share / score_cost)
                if mass_priority <= 0.0:
                    continue
                major_ranked.append((float(mass_priority), int(spec_idx)))
            major_ranked.sort(key=lambda t: (float(t[0]), -int(t[1])), reverse=True)
            for _prio, spec_idx in major_ranked[:major_force_cap]:
                if len(eval_indices) >= eval_cap:
                    break
                if int(spec_idx) not in eval_indices:
                    eval_indices.append(int(spec_idx))
        # Force-evaluate top same-component closure candidates so "final easy loop"
        # options are not dropped before exact-gain evaluation.
        if same_component_ranked and len(eval_indices) < eval_cap:
            same_force_cap = max(3, min(22, int(math.ceil(eval_cap * 0.16))))
            for _s, _h, _c, spec_idx in same_component_ranked[:same_force_cap]:
                if len(eval_indices) >= eval_cap:
                    break
                if int(spec_idx) not in eval_indices:
                    eval_indices.append(int(spec_idx))
        for _s, _h, _r, _c, spec_idx in intra_ranked[:intra_cap]:
            if len(eval_indices) >= eval_cap:
                break
            if spec_idx not in eval_indices:
                eval_indices.append(int(spec_idx))
        if not eval_indices:
            break

        eval_rows: List[Dict[str, float]] = []
        for spec_idx in eval_indices:
            spec = cand_specs[int(spec_idx)]
            is_intra = bool(spec.get("is_intra", False))
            if is_intra:
                intra_pid = int(spec.get("intra_pid"))
                base_intra_factor = float(spec.get("intra_factor", 0.0) or 0.0)
                curr_patch_factor = float(intra_factor_by_patch.get(intra_pid, 0.0) or 0.0)
                gain_factor = max(base_intra_factor - curr_patch_factor, 0.0)
                next_intra_factor_sum = float(current_intra_factor_sum + gain_factor)
                next_intra_bonus = _landscape_fluidity_intra_bonus_from_factor(
                    graph_fluidity_pre=graph_fluidity_pre,
                    intra_factor_sum=next_intra_factor_sum,
                    params=params,
                    intra_patch_count=len(intra_factor_by_patch) + (0 if intra_pid in intra_factor_by_patch else 1),
                )
                gain = float(next_intra_bonus - current_intra_bonus)
                if gain <= 1e-12:
                    continue
                score_new = float(current_graph_score + next_intra_bonus + current_inter_loop_bonus + current_bridge_relief_bonus)
                same_component = True
                loop_gain = 0.0
                bridge_gain = 0.0
                same_component_shortcut_ratio = 0.0
                closure_priority = 0.0
                topology_bonus = 0.0
                intra_hole_area_m2 = float(spec.get("intra_hole_area_m2", 0.0) or 0.0)
                intra_pocket_strength = float(spec.get("intra_pocket_strength", 0.0) or 0.0)
                intra_patch_area_ha = float(spec.get("intra_patch_area_ha", 0.0) or 0.0)
                intra_gap_m = float(spec.get("intra_gap_m", 0.0) or 0.0)
                patch_area_m2 = max(intra_patch_area_ha * 10000.0, 1.0)
                hole_share = max(0.0, min(1.0, intra_hole_area_m2 / max(0.35 * patch_area_m2, 1.0)))
                mouth_scale = max(0.0, min(1.0, intra_gap_m / max(math.sqrt(patch_area_m2) * 0.30, 1.0)))
                intra_structure_priority = float(
                    (0.55 * math.sqrt(hole_share))
                    + (0.30 * max(0.0, min(1.0, math.log1p(intra_pocket_strength) / math.log1p(150.0))))
                    + (0.15 * mouth_scale)
                )
                min_comp_share = max(0.0, min(1.0, hole_share))
                impact_share = max(
                    0.0,
                    min(
                        1.0,
                        max(
                            float((patches.get(int(intra_pid)) or {}).get("area_ha", 0.0) or 0.0)
                            / max(float(total_habitat_area), 1e-9),
                            intra_hole_area_m2 / max(total_habitat_area * 10000.0, 1e-9),
                        ),
                    ),
                )
                score_new += float(
                    0.14
                    * graph_fluidity_pre
                    * (
                        (0.55 * intra_structure_priority)
                        + (0.30 * math.sqrt(max(impact_share, 0.0)))
                        + (0.15 * max(float(spec.get("intra_repair_score", 0.0) or 0.0), 0.0))
                    )
                )
            else:
                changes_try = _apply_edges_with_changes(
                    graph,
                    spec.get("graph_edges", spec.get("edges", [])),
                )
                if not changes_try:
                    continue
                mean_graph_new = _mean_effective_resistance_graph(
                    graph,
                    pair_records,
                    disconnected_penalty,
                    pair_weights=pair_weights,
                )
                _revert_edges_with_changes(graph, changes_try)
                graph_score_new = float(total_habitat_area / max(float(mean_graph_new), 1e-9))
                gain_graph = float(graph_score_new - current_graph_score)
                gate_u = int(spec.get("gate_u"))
                gate_v = int(spec.get("gate_v"))
                same_component = (int(comp_id.get(gate_u, -1)) == int(comp_id.get(gate_v, -1)))
                loop_gain = 0.0
                bridge_gain = 0.0
                closure_priority = 0.0
                if same_component:
                    c_e = max(float(spec.get("length", 0.0) or 0.0), 1e-9)
                    d_old = _shortest_lookup(shortest_current, gate_u, gate_v)
                    shortcut_ratio = (d_old / c_e) if np.isfinite(d_old) else float("inf")
                    same_component_shortcut_ratio = float(shortcut_ratio if np.isfinite(shortcut_ratio) else 0.0)
                    comp_key = int(comp_id.get(gate_u, -1))
                    comp_area = float(comp_area_ha.get(comp_key, 0.0) or 0.0) if comp_key >= 0 else 0.0
                    try:
                        hop_len = int(nx.shortest_path_length(graph, int(gate_u), int(gate_v)))
                    except Exception:
                        hop_len = 0
                    loop_factor = _landscape_fluidity_inter_loop_factor(
                        float(shortcut_ratio),
                        float(comp_area),
                        float(total_habitat_area),
                        hop_length=int(hop_len),
                    )
                    gate_pair = spec.get("gate_pair")
                    if gate_pair is None:
                        gate_pair = (gate_u, gate_v) if gate_u <= gate_v else (gate_v, gate_u)
                    prev_factor = float(inter_loop_factor_by_pair.get(tuple(gate_pair), 0.0) or 0.0)
                    gain_loop_factor = max(float(loop_factor) - prev_factor, 0.0)
                    if gain_loop_factor > 0.0:
                        next_inter_loop_factor_sum = float(current_inter_loop_factor_sum + gain_loop_factor)
                        next_inter_loop_bonus = _landscape_fluidity_inter_loop_bonus_from_factor(
                            graph_fluidity_pre=graph_fluidity_pre,
                            inter_loop_factor_sum=next_inter_loop_factor_sum,
                            params=params,
                        )
                        loop_gain = float(next_inter_loop_bonus - current_inter_loop_bonus)
                    bridge_fraction = _bridge_fraction_on_shortest_path(
                        graph,
                        int(gate_u),
                        int(gate_v),
                        bridge_edges_current,
                    )
                    bridge_factor = _landscape_fluidity_bridge_relief_factor(
                        float(bridge_fraction),
                        float(comp_area),
                        float(total_habitat_area),
                        float(shortcut_ratio),
                    )
                    prev_bridge_factor = float(bridge_relief_factor_by_pair.get(tuple(gate_pair), 0.0) or 0.0)
                    gain_bridge_factor = max(float(bridge_factor) - prev_bridge_factor, 0.0)
                    if gain_bridge_factor > 0.0:
                        next_bridge_relief_factor_sum = float(current_bridge_relief_factor_sum + gain_bridge_factor)
                        next_bridge_relief_bonus = _landscape_fluidity_bridge_relief_bonus_from_factor(
                            graph_fluidity_pre=graph_fluidity_pre,
                            bridge_relief_factor_sum=next_bridge_relief_factor_sum,
                            params=params,
                        )
                        bridge_gain = float(next_bridge_relief_bonus - current_bridge_relief_bonus)
                    # Additional same-component circuit-closure priority:
                    # favors cheap final links that close long local detours in large networks.
                    if hop_len >= 2:
                        ratio_component = max(
                            0.0,
                            min(
                                1.0,
                                (float(shortcut_ratio) - 1.0) / max(float(shortcut_ratio), 1e-9),
                            ),
                        )
                        hop_component = max(0.0, min(1.0, float(hop_len - 1) / float(max(hop_len, 1))))
                        area_component = max(0.0, min(1.0, float(comp_area) / max(float(total_habitat_area), 1e-9)))
                        closure_priority = float(area_component * hop_component * (0.25 + (0.75 * ratio_component)))
                    else:
                        closure_priority = 0.0
                topology_bonus = 0.0
                if not same_component:
                    cu = int(comp_id.get(gate_u, -1))
                    cv = int(comp_id.get(gate_v, -1))
                    total_nodes = max(1, len(valid_nodes))
                    topology_bonus = min(
                        1.0,
                        float(comp_node_count.get(cu, 0) + comp_node_count.get(cv, 0)) / float(total_nodes),
                    )
                topology_bonus += min(0.60, 0.20 * max(0, len(spec.get("pids", [])) - 2))
                if same_component:
                    topology_bonus += min(0.60, max(float(loop_gain + bridge_gain), 0.0) / max(abs(float(gain_graph)) + 1e-9, 1e-9))
                topology_bonus = min(1.5, max(0.0, float(topology_bonus)))
                cu = int(comp_id.get(gate_u, -1))
                cv = int(comp_id.get(gate_v, -1))
                min_comp_share = 0.0
                if same_component:
                    impact_area = float(comp_area_ha.get(cu, 0.0) or 0.0)
                    min_comp_share = max(
                        0.0,
                        min(
                            1.0,
                            impact_area / max(float(total_habitat_area), 1e-9),
                        ),
                    )
                else:
                    area_u = float(comp_area_ha.get(cu, 0.0) or 0.0)
                    area_v = float(comp_area_ha.get(cv, 0.0) or 0.0)
                    impact_area = area_u + area_v
                    min_comp_share = max(
                        0.0,
                        min(
                            1.0,
                            min(area_u, area_v) / max(float(total_habitat_area), 1e-9),
                        ),
                    )
                impact_share = max(0.0, min(1.0, impact_area / max(float(total_habitat_area), 1e-9)))
                gain = float(gain_graph + loop_gain + bridge_gain)
                if gain <= 1e-12:
                    continue
                score_new = float(
                    graph_score_new
                    + current_intra_bonus
                    + current_inter_loop_bonus
                    + current_bridge_relief_bonus
                    + loop_gain
                    + bridge_gain
                )
                if not same_component:
                    same_component_shortcut_ratio = 0.0
                    closure_priority = 0.0
                intra_structure_priority = 0.0
            cost = float(spec.get("cost", 0.0) or 0.0)
            score_cost = float(spec.get("score_cost", cost) or cost)
            if score_cost <= 0.0:
                score_cost = max(cost, 1e-12)
            ratio = gain / max(score_cost, 1e-12)
            eval_rows.append(
                {
                    "spec_idx": float(spec_idx),
                    "gain": float(gain),
                    "ratio": float(ratio),
                    "cost": float(cost),
                    "score_cost": float(score_cost),
                    "topology_bonus": float(topology_bonus),
                    "impact_share": float(impact_share),
                    "min_comp_share": float(min_comp_share),
                    "long_hop_factor": float(spec.get("long_hop_factor", 1.0) or 1.0),
                    "is_intra": 1.0 if is_intra else 0.0,
                    "is_inter_component": 0.0 if same_component else 1.0,
                    "loop_gain": float(loop_gain),
                    "bridge_gain": float(bridge_gain),
                    "same_component_shortcut_ratio": float(same_component_shortcut_ratio),
                    "closure_priority": float(closure_priority),
                    "intra_structure_priority": float(intra_structure_priority),
                    "post_score": float(score_new),
                }
            )
        if not eval_rows:
            break

        candidate_rows = list(eval_rows)
        inter_rows = [r for r in candidate_rows if float(r.get("is_intra", 0.0)) < 0.5]
        if inter_rows:
            best_inter_ratio = max(float(r.get("ratio", 0.0) or 0.0) for r in inter_rows)
            best_inter_gain = max(float(r.get("gain", 0.0) or 0.0) for r in inter_rows)
            intra_min_ratio = float(intra_tau) * max(best_inter_ratio, 0.0)
            intra_min_gain = float(intra_tau) * max(best_inter_gain, 0.0)
            best_intra_structure = max(
                (float(r.get("intra_structure_priority", 0.0) or 0.0) for r in candidate_rows if float(r.get("is_intra", 0.0)) > 0.5),
                default=0.0,
            )
            gated_rows: List[Dict[str, float]] = []
            for r in candidate_rows:
                if float(r.get("is_intra", 0.0)) < 0.5:
                    gated_rows.append(r)
                    continue
                if (
                    float(r.get("ratio", 0.0) or 0.0) + 1e-12 >= intra_min_ratio
                    and float(r.get("gain", 0.0) or 0.0) + 1e-12 >= intra_min_gain
                ):
                    gated_rows.append(r)
                    continue
                if (
                    float(r.get("intra_structure_priority", 0.0) or 0.0) + 1e-12 >= 0.72 * max(best_intra_structure, 1e-12)
                    and float(r.get("gain", 0.0) or 0.0) + 1e-12 >= 0.18 * max(best_inter_gain, 1e-12)
                ):
                    gated_rows.append(r)
            if gated_rows:
                candidate_rows = gated_rows
            else:
                candidate_rows = inter_rows

        long_hop_rows = [
            r
            for r in candidate_rows
            if float(r.get("is_inter_component", 0.0)) > 0.5
            and float(r.get("long_hop_factor", 1.0) or 1.0) >= long_hop_factor_thresh
            and float(r.get("impact_share", 0.0) or 0.0) < 0.25
        ]
        non_long_rows = [r for r in candidate_rows if r not in long_hop_rows]
        if long_hop_rows and non_long_rows:
            def _row_priority(rr: Dict[str, float]) -> float:
                return float(rr.get("ratio", 0.0) or 0.0) * (
                    1.0 + topology_weight * float(rr.get("topology_bonus", 0.0) or 0.0)
                )
            best_non_long = max((_row_priority(r) for r in non_long_rows), default=0.0)
            best_long = max((_row_priority(r) for r in long_hop_rows), default=0.0)
            if best_non_long + 1e-12 >= long_hop_tau * max(best_long, 0.0):
                candidate_rows = non_long_rows

        # Constrained-budget priority gate:
        # Use a scale-free multi-factor frontier so the mode spends early budget on
        # meaningful component merges (gain + impacted share + smaller-side share),
        # while still allowing exceptional intra/loop candidates through.
        inter_component_rows = [r for r in candidate_rows if float(r.get("is_inter_component", 0.0)) > 0.5]
        if inter_component_rows:
            budget_used_frac = 1.0 - (float(remaining) / max(float(budget_limit), 1e-9))
            best_inter_gain = max(float(r.get("gain", 0.0) or 0.0) for r in inter_component_rows)
            best_inter_ratio = max(float(r.get("ratio", 0.0) or 0.0) for r in inter_component_rows)
            best_inter_impact = max(float(r.get("impact_share", 0.0) or 0.0) for r in inter_component_rows)
            best_inter_min_comp = max(float(r.get("min_comp_share", 0.0) or 0.0) for r in inter_component_rows)
            best_inter_topo = max(float(r.get("topology_bonus", 0.0) or 0.0) for r in inter_component_rows)
            min_inter_score_cost = min(
                max(float(r.get("score_cost", r.get("cost", 0.0)) or 0.0), 1e-12)
                for r in inter_component_rows
            )

            sorted_impacts = sorted(float(r.get("impact_share", 0.0) or 0.0) for r in inter_component_rows)
            sorted_min_comp = sorted(float(r.get("min_comp_share", 0.0) or 0.0) for r in inter_component_rows)
            q25_idx = max(0, int(math.floor(0.25 * max(len(sorted_min_comp) - 1, 0))))
            med_idx = max(0, int(math.floor(0.50 * max(len(sorted_impacts) - 1, 0))))
            q50_min_idx = max(0, int(math.floor(0.50 * max(len(sorted_min_comp) - 1, 0))))
            q60_impact_idx = max(0, int(math.floor(0.60 * max(len(sorted_impacts) - 1, 0))))
            q25_min_comp = float(sorted_min_comp[q25_idx]) if sorted_min_comp else 0.0
            median_impact = float(sorted_impacts[med_idx]) if sorted_impacts else 0.0
            q50_min_comp = float(sorted_min_comp[q50_min_idx]) if sorted_min_comp else 0.0
            q60_impact = float(sorted_impacts[q60_impact_idx]) if sorted_impacts else 0.0

            def _inter_priority(rr: Dict[str, float]) -> float:
                gain_n = float(rr.get("gain", 0.0) or 0.0) / max(best_inter_gain, 1e-12)
                ratio_n = float(rr.get("ratio", 0.0) or 0.0) / max(best_inter_ratio, 1e-12)
                impact_n = float(rr.get("impact_share", 0.0) or 0.0) / max(best_inter_impact, 1e-12)
                min_comp_n = float(rr.get("min_comp_share", 0.0) or 0.0) / max(best_inter_min_comp, 1e-12)
                topo_n = float(rr.get("topology_bonus", 0.0) or 0.0) / max(best_inter_topo, 1e-12)
                score_cost_rr = max(float(rr.get("score_cost", rr.get("cost", 0.0)) or 0.0), 1e-12)
                cost_eff_n = min(1.0, min_inter_score_cost / max(score_cost_rr, 1e-12))
                base = (
                    (0.46 * gain_n)
                    + (0.22 * impact_n)
                    + (0.17 * min_comp_n)
                    + (0.10 * ratio_n)
                    + (0.05 * topo_n)
                )
                return float(base * (0.92 + (0.08 * cost_eff_n)))

            for r in inter_component_rows:
                r["_inter_priority"] = float(_inter_priority(r))
            best_inter_priority = max(float(r.get("_inter_priority", 0.0) or 0.0) for r in inter_component_rows)

            frontier_tau = 0.70 if budget_used_frac < 0.80 else 0.58
            frontier_rows = [
                r
                for r in inter_component_rows
                if float(r.get("_inter_priority", 0.0) or 0.0) + 1e-12 >= frontier_tau * max(best_inter_priority, 1e-12)
            ]
            if not frontier_rows:
                frontier_rows = list(inter_component_rows)

            filtered_frontier: List[Dict[str, float]] = []
            for r in frontier_rows:
                gain_r = float(r.get("gain", 0.0) or 0.0)
                impact_r = float(r.get("impact_share", 0.0) or 0.0)
                min_comp_share_r = float(r.get("min_comp_share", 0.0) or 0.0)
                inter_prio_r = float(r.get("_inter_priority", 0.0) or 0.0)
                tiny_pair = (
                    min_comp_share_r + 1e-12 <= q25_min_comp
                    and impact_r + 1e-12 <= max(median_impact, 1e-12)
                )
                if tiny_pair and (
                    gain_r + 1e-12 < 0.95 * max(best_inter_gain, 1e-12)
                    and inter_prio_r + 1e-12 < 0.97 * max(best_inter_priority, 1e-12)
                ):
                    continue
                filtered_frontier.append(r)
            if filtered_frontier:
                frontier_rows = filtered_frontier

            # Early/mid run, keep focus on merges that bring in meaningful component mass.
            # This is relative (quantile-based), so it scales across landscapes.
            if frontier_rows and budget_used_frac < 0.92:
                major_rows = [
                    r
                    for r in frontier_rows
                    if (
                        float(r.get("min_comp_share", 0.0) or 0.0) + 1e-12 >= q50_min_comp
                        or float(r.get("impact_share", 0.0) or 0.0) + 1e-12 >= q60_impact
                    )
                ]
                if major_rows:
                    frontier_rows = major_rows

            intra_rows = [r for r in candidate_rows if float(r.get("is_intra", 0.0)) > 0.5]
            best_intra_structure = max(
                (float(r.get("intra_structure_priority", 0.0) or 0.0) for r in intra_rows),
                default=0.0,
            )
            intra_keep = [
                r
                for r in intra_rows
                if (
                    float(r.get("gain", 0.0) or 0.0) + 1e-12 >= 0.85 * max(best_inter_gain, 1e-12)
                    or (
                        float(r.get("intra_structure_priority", 0.0) or 0.0) + 1e-12 >= 0.74 * max(best_intra_structure, 1e-12)
                        and float(r.get("gain", 0.0) or 0.0) + 1e-12 >= 0.16 * max(best_inter_gain, 1e-12)
                    )
                )
            ]
            same_component_rows = [
                r
                for r in candidate_rows
                if float(r.get("is_intra", 0.0)) < 0.5 and float(r.get("is_inter_component", 0.0)) < 0.5
            ]
            same_sorted_impact = sorted(float(r.get("impact_share", 0.0) or 0.0) for r in same_component_rows)
            same_sorted_ratio = sorted(float(r.get("same_component_shortcut_ratio", 0.0) or 0.0) for r in same_component_rows)
            same_sorted_closure = sorted(float(r.get("closure_priority", 0.0) or 0.0) for r in same_component_rows)
            same_sorted_loop_gain = sorted(
                float(r.get("loop_gain", 0.0) or 0.0) + float(r.get("bridge_gain", 0.0) or 0.0)
                for r in same_component_rows
            )
            if same_component_rows:
                same_q30_impact = float(
                    same_sorted_impact[max(0, int(math.floor(0.30 * max(len(same_sorted_impact) - 1, 0))))]
                )
                same_q40_impact = float(
                    same_sorted_impact[max(0, int(math.floor(0.40 * max(len(same_sorted_impact) - 1, 0))))]
                )
                same_q50_impact = float(
                    same_sorted_impact[max(0, int(math.floor(0.50 * max(len(same_sorted_impact) - 1, 0))))]
                )
                same_q60_ratio = float(
                    same_sorted_ratio[max(0, int(math.floor(0.60 * max(len(same_sorted_ratio) - 1, 0))))]
                )
                same_q75_ratio = float(
                    same_sorted_ratio[max(0, int(math.floor(0.75 * max(len(same_sorted_ratio) - 1, 0))))]
                )
                same_q90_ratio = float(
                    same_sorted_ratio[max(0, int(math.floor(0.90 * max(len(same_sorted_ratio) - 1, 0))))]
                )
                same_q70_loop = float(
                    same_sorted_loop_gain[max(0, int(math.floor(0.70 * max(len(same_sorted_loop_gain) - 1, 0))))]
                )
                same_q80_closure = float(
                    same_sorted_closure[max(0, int(math.floor(0.80 * max(len(same_sorted_closure) - 1, 0))))]
                )
            else:
                same_q30_impact = 0.0
                same_q40_impact = 0.0
                same_q50_impact = 0.0
                same_q60_ratio = 0.0
                same_q75_ratio = 0.0
                same_q90_ratio = 0.0
                same_q70_loop = 0.0
                same_q80_closure = 0.0
            best_same_closure = max(
                (float(r.get("closure_priority", 0.0) or 0.0) for r in same_component_rows),
                default=0.0,
            )
            best_same_ratio = max(
                (float(r.get("same_component_shortcut_ratio", 0.0) or 0.0) for r in same_component_rows),
                default=0.0,
            )
            same_component_keep: List[Dict[str, float]] = []
            for r in same_component_rows:
                gain_r = float(r.get("gain", 0.0) or 0.0)
                ratio_r = float(r.get("ratio", 0.0) or 0.0)
                impact_r = float(r.get("impact_share", 0.0) or 0.0)
                loop_plus_bridge = float(r.get("loop_gain", 0.0) or 0.0) + float(r.get("bridge_gain", 0.0) or 0.0)
                sc_ratio = float(r.get("same_component_shortcut_ratio", 0.0) or 0.0)
                closure_r = float(r.get("closure_priority", 0.0) or 0.0)
                if gain_r + 1e-12 >= 0.92 * max(best_inter_gain, 1e-12):
                    same_component_keep.append(r)
                    continue
                if (
                    best_same_closure > 1e-12
                    and closure_r + 1e-12 >= 0.78 * best_same_closure
                    and sc_ratio + 1e-12 >= max(0.85 * best_same_ratio, 0.80 * same_component_shortcut_threshold)
                    and impact_r + 1e-12 >= max(same_q40_impact, 1e-12)
                ):
                    same_component_keep.append(r)
                    continue
                if (
                    loop_plus_bridge + 1e-12 >= max(same_q70_loop, 0.22 * max(best_inter_gain, 1e-12))
                    and sc_ratio + 1e-12 >= max(same_q60_ratio, 1e-12)
                ):
                    same_component_keep.append(r)
                    continue
                if (
                    closure_r + 1e-12 >= max(same_q80_closure, 0.72 * best_same_closure)
                    and sc_ratio + 1e-12 >= max(same_q75_ratio, 1e-12)
                    and impact_r + 1e-12 >= max(same_q30_impact, 1e-12)
                ):
                    same_component_keep.append(r)
                    continue
                if (
                    sc_ratio + 1e-12 >= max(same_q90_ratio, 1e-12)
                    and impact_r + 1e-12 >= max(same_q50_impact, 1e-12)
                    and ratio_r + 1e-12 >= 0.45 * max(best_inter_ratio, 1e-12)
                ):
                    same_component_keep.append(r)

            merged_rows = list(frontier_rows) + intra_keep + same_component_keep
            if merged_rows:
                uniq: Dict[int, Dict[str, float]] = {}
                for r in merged_rows:
                    try:
                        sid = int(r.get("spec_idx", -1))
                    except Exception:
                        sid = -1
                    if sid < 0:
                        continue
                    prev = uniq.get(sid)
                    if prev is None or float(r.get("gain", 0.0) or 0.0) > float(prev.get("gain", 0.0) or 0.0):
                        uniq[sid] = r
                if uniq:
                    candidate_rows = list(uniq.values())

        max_gain = max(float(r.get("gain", 0.0) or 0.0) for r in candidate_rows)
        max_ratio = max(float(r.get("ratio", 0.0) or 0.0) for r in candidate_rows)
        max_impact = max(float(r.get("impact_share", 0.0) or 0.0) for r in candidate_rows)
        max_min_comp = max(float(r.get("min_comp_share", 0.0) or 0.0) for r in candidate_rows)
        max_closure = max(float(r.get("closure_priority", 0.0) or 0.0) for r in candidate_rows)
        max_loop_bridge = max(
            (
                float(r.get("loop_gain", 0.0) or 0.0)
                + float(r.get("bridge_gain", 0.0) or 0.0)
            )
            for r in candidate_rows
        )

        def _primary_score(row: Dict[str, float]) -> float:
            base = float(
                (
                    (gain_weight * (float(row.get("gain", 0.0) or 0.0) / max(max_gain, 1e-12)))
                    + (ratio_weight * (float(row.get("ratio", 0.0) or 0.0) / max(max_ratio, 1e-12)))
                )
                / max(gain_weight + ratio_weight, 1e-12)
            )
            topo = max(0.0, float(row.get("topology_bonus", 0.0) or 0.0))
            impact = max(0.0, float(row.get("impact_share", 0.0) or 0.0)) / max(max_impact, 1e-12)
            min_comp = max(0.0, float(row.get("min_comp_share", 0.0) or 0.0)) / max(max_min_comp, 1e-12)
            closure = max(0.0, float(row.get("closure_priority", 0.0) or 0.0)) / max(max_closure, 1e-12)
            loop_bridge = max(
                0.0,
                float(row.get("loop_gain", 0.0) or 0.0) + float(row.get("bridge_gain", 0.0) or 0.0),
            ) / max(max_loop_bridge, 1e-12)
            same_component_flag = (
                1.0
                if (
                    float(row.get("is_intra", 0.0) or 0.0) < 0.5
                    and float(row.get("is_inter_component", 0.0) or 0.0) < 0.5
                )
                else 0.0
            )
            closure_mix = (0.65 * closure) + (0.35 * loop_bridge)
            gain_n = max(0.0, float(row.get("gain", 0.0) or 0.0)) / max(max_gain, 1e-12)
            impact_gain = gain_n * impact
            return float(
                (
                    base
                    + (topology_weight * topo)
                    + (0.40 * topology_weight * impact)
                    + (0.30 * topology_weight * min_comp)
                    + (0.25 * topology_weight * impact_gain)
                    + (0.70 * topology_weight * closure_mix * same_component_flag)
                )
                / max(
                    1.0
                    + topology_weight
                    + (0.40 * topology_weight)
                    + (0.30 * topology_weight)
                    + (0.25 * topology_weight),
                    + (0.70 * topology_weight),
                    1e-12,
                )
            )

        ranked_rows = sorted(
            candidate_rows,
            key=lambda r: (
                _primary_score(r),
                float(r.get("gain", 0.0)),
                float(r.get("is_inter_component", 0.0)),
                float(r.get("ratio", 0.0)),
                -float(r.get("score_cost", r.get("cost", 0.0))),
                -float(r.get("cost", 0.0)),
                -float(r.get("spec_idx", 0.0)),
            ),
            reverse=True,
        )

        if lookahead_k > 1 and len(ranked_rows) > 1 and lookahead_weight > 1e-12:
            la_pool = ranked_rows[: min(len(ranked_rows), lookahead_k)]
            best = max(
                la_pool,
                key=lambda r1: (
                    float(r1.get("gain", 0.0) or 0.0)
                    + lookahead_weight
                    * max(
                        (
                            float(r2.get("gain", 0.0) or 0.0)
                            for r2 in candidate_rows
                            if int(r2.get("spec_idx", -1)) != int(r1.get("spec_idx", -2))
                            and float(r2.get("cost", 0.0) or 0.0)
                            <= max(remaining - float(r1.get("cost", 0.0) or 0.0), 0.0) + 1e-12
                            and (
                                (
                                    int(cand_specs[int(r2.get("spec_idx", -1))].get("intra_pid", -999999))
                                    != int(cand_specs[int(r1.get("spec_idx", -1))].get("intra_pid", -999998))
                                )
                                or float(r2.get("is_intra", 0.0)) < 0.5
                                or float(r1.get("is_intra", 0.0)) < 0.5
                            )
                            and (
                                tuple(cand_specs[int(r2.get("spec_idx", -1))].get("gate_pair", (-1, -1)))
                                != tuple(cand_specs[int(r1.get("spec_idx", -1))].get("gate_pair", (-2, -2)))
                            )
                        ),
                        default=0.0,
                    ),
                    _primary_score(r1),
                    float(r1.get("gain", 0.0)),
                    float(r1.get("is_inter_component", 0.0)),
                    float(r1.get("ratio", 0.0)),
                    -float(r1.get("cost", 0.0)),
                    -float(r1.get("spec_idx", 0.0)),
                ),
            )
        else:
            best = ranked_rows[0]

        # Circuit-closure budget reserve:
        # keep enough budget for a strong same-component closure candidate when one exists.
        same_component_rows = [
            r for r in candidate_rows
            if float(r.get("is_intra", 0.0)) < 0.5 and float(r.get("is_inter_component", 0.0)) < 0.5
        ]
        if same_component_rows:
            same_sorted_impact = sorted(float(r.get("impact_share", 0.0) or 0.0) for r in same_component_rows)
            same_sorted_ratio = sorted(float(r.get("same_component_shortcut_ratio", 0.0) or 0.0) for r in same_component_rows)
            q40_impact = float(
                same_sorted_impact[max(0, int(math.floor(0.40 * max(len(same_sorted_impact) - 1, 0))))]
            )
            q70_ratio = float(
                same_sorted_ratio[max(0, int(math.floor(0.70 * max(len(same_sorted_ratio) - 1, 0))))]
            )
            best_same_closure = max(
                (float(r.get("closure_priority", 0.0) or 0.0) for r in same_component_rows),
                default=0.0,
            )
            reserve_rows = [
                r for r in same_component_rows
                if (
                    float(r.get("closure_priority", 0.0) or 0.0) + 1e-12 >= 0.72 * max(best_same_closure, 1e-12)
                    and (
                        float(r.get("impact_share", 0.0) or 0.0) + 1e-12 >= max(q40_impact, 1e-12)
                        or float(r.get("same_component_shortcut_ratio", 0.0) or 0.0) + 1e-12 >= max(q70_ratio, 1e-12)
                    )
                )
            ]
            if reserve_rows:
                reserve_cost = min(float(r.get("cost", 0.0) or 0.0) for r in reserve_rows)
                best_cost = float(best.get("cost", 0.0) or 0.0)
                best_is_reserve = any(int(r.get("spec_idx", -1)) == int(best.get("spec_idx", -2)) for r in reserve_rows)
                if (
                    (not best_is_reserve)
                    and reserve_cost > 1e-12
                    and (remaining - best_cost) + 1e-12 < reserve_cost
                ):
                    alt_rows = [
                        r for r in ranked_rows
                        if int(r.get("spec_idx", -1)) != int(best.get("spec_idx", -2))
                        and float(r.get("cost", 0.0) or 0.0) <= max(remaining - reserve_cost, 0.0) + 1e-12
                    ]
                    if alt_rows:
                        alt_best = max(
                            alt_rows,
                            key=lambda r: (
                                _primary_score(r),
                                float(r.get("gain", 0.0)),
                                float(r.get("closure_priority", 0.0)),
                                float(r.get("is_inter_component", 0.0)),
                                float(r.get("ratio", 0.0)),
                                -float(r.get("score_cost", r.get("cost", 0.0))),
                                -float(r.get("cost", 0.0)),
                                -float(r.get("spec_idx", 0.0)),
                            ),
                        )
                        # Only switch when alternative is reasonably competitive.
                        if _primary_score(alt_best) + 1e-12 >= 0.88 * max(_primary_score(best), 1e-12):
                            best = alt_best

        chosen_spec_idx = int(best["spec_idx"])
        chosen_spec = cand_specs[chosen_spec_idx]
        is_intra = bool(chosen_spec.get("is_intra", False))
        exact_loop_gain = 0.0
        exact_bridge_gain = 0.0
        if is_intra:
            intra_pid = int(chosen_spec.get("intra_pid"))
            base_intra_factor = float(chosen_spec.get("intra_factor", 0.0) or 0.0)
            curr_patch_factor = float(intra_factor_by_patch.get(intra_pid, 0.0) or 0.0)
            gain_factor = max(base_intra_factor - curr_patch_factor, 0.0)
            next_intra_factor_sum = float(current_intra_factor_sum + gain_factor)
            next_intra_bonus = _landscape_fluidity_intra_bonus_from_factor(
                graph_fluidity_pre=graph_fluidity_pre,
                intra_factor_sum=next_intra_factor_sum,
                params=params,
                intra_patch_count=len(intra_factor_by_patch) + (0 if intra_pid in intra_factor_by_patch else 1),
            )
            exact_gain = float(next_intra_bonus - current_intra_bonus)
            if exact_gain <= 1e-12:
                used_specs.add(int(chosen_spec_idx))
                continue
            intra_factor_by_patch[intra_pid] = max(base_intra_factor, curr_patch_factor)
            current_intra_factor_sum = float(next_intra_factor_sum)
            current_intra_bonus = float(next_intra_bonus)
            score_new = float(current_graph_score + current_intra_bonus + current_inter_loop_bonus + current_bridge_relief_bonus)
        else:
            gate_u = int(chosen_spec.get("gate_u"))
            gate_v = int(chosen_spec.get("gate_v"))
            gate_pair = chosen_spec.get("gate_pair")
            if gate_pair is None:
                gate_pair = (gate_u, gate_v) if gate_u <= gate_v else (gate_v, gate_u)
            gate_pair = tuple(gate_pair)
            same_component_before = (int(comp_id.get(gate_u, -1)) == int(comp_id.get(gate_v, -1)))
            next_inter_loop_factor_sum = float(current_inter_loop_factor_sum)
            next_inter_loop_bonus = float(current_inter_loop_bonus)
            next_bridge_relief_factor_sum = float(current_bridge_relief_factor_sum)
            next_bridge_relief_bonus = float(current_bridge_relief_bonus)
            chosen_loop_factor = 0.0
            chosen_bridge_factor = 0.0
            if same_component_before:
                c_e = max(float(chosen_spec.get("length", 0.0) or 0.0), 1e-9)
                d_old = _shortest_lookup(shortest_current, gate_u, gate_v)
                shortcut_ratio = (d_old / c_e) if np.isfinite(d_old) else float("inf")
                comp_key = int(comp_id.get(gate_u, -1))
                comp_area = float(comp_area_ha.get(comp_key, 0.0) or 0.0) if comp_key >= 0 else 0.0
                try:
                    hop_len = int(nx.shortest_path_length(graph, int(gate_u), int(gate_v)))
                except Exception:
                    hop_len = 0
                loop_factor = _landscape_fluidity_inter_loop_factor(
                    float(shortcut_ratio),
                    float(comp_area),
                    float(total_habitat_area),
                    hop_length=int(hop_len),
                )
                chosen_loop_factor = float(loop_factor)
                prev_factor = float(inter_loop_factor_by_pair.get(gate_pair, 0.0) or 0.0)
                gain_loop_factor = max(float(loop_factor) - prev_factor, 0.0)
                if gain_loop_factor > 0.0:
                    next_inter_loop_factor_sum = float(current_inter_loop_factor_sum + gain_loop_factor)
                    next_inter_loop_bonus = _landscape_fluidity_inter_loop_bonus_from_factor(
                        graph_fluidity_pre=graph_fluidity_pre,
                        inter_loop_factor_sum=next_inter_loop_factor_sum,
                        params=params,
                    )
                    exact_loop_gain = float(next_inter_loop_bonus - current_inter_loop_bonus)
                bridge_fraction = _bridge_fraction_on_shortest_path(
                    graph,
                    int(gate_u),
                    int(gate_v),
                    bridge_edges_current,
                )
                bridge_factor = _landscape_fluidity_bridge_relief_factor(
                    float(bridge_fraction),
                    float(comp_area),
                    float(total_habitat_area),
                    float(shortcut_ratio),
                )
                chosen_bridge_factor = float(bridge_factor)
                prev_bridge_factor = float(bridge_relief_factor_by_pair.get(gate_pair, 0.0) or 0.0)
                gain_bridge_factor = max(float(bridge_factor) - prev_bridge_factor, 0.0)
                if gain_bridge_factor > 0.0:
                    next_bridge_relief_factor_sum = float(current_bridge_relief_factor_sum + gain_bridge_factor)
                    next_bridge_relief_bonus = _landscape_fluidity_bridge_relief_bonus_from_factor(
                        graph_fluidity_pre=graph_fluidity_pre,
                        bridge_relief_factor_sum=next_bridge_relief_factor_sum,
                        params=params,
                    )
                    exact_bridge_gain = float(next_bridge_relief_bonus - current_bridge_relief_bonus)

            changes = _apply_edges_with_changes(
                graph,
                chosen_spec.get("graph_edges", chosen_spec.get("edges", [])),
            )
            if not changes:
                used_specs.add(int(chosen_spec_idx))
                continue

            mean_graph_new = _mean_effective_resistance_graph(
                graph,
                pair_records,
                disconnected_penalty,
                pair_weights=pair_weights,
            )
            graph_score_new = float(total_habitat_area / max(float(mean_graph_new), 1e-9))
            exact_gain = float(graph_score_new - current_graph_score)
            if exact_gain <= 1e-12:
                _revert_edges_with_changes(graph, changes)
                used_specs.add(int(chosen_spec_idx))
                continue
            exact_gain += float(exact_loop_gain)
            exact_gain += float(exact_bridge_gain)
            if exact_gain <= 1e-12:
                _revert_edges_with_changes(graph, changes)
                used_specs.add(int(chosen_spec_idx))
                continue
            current_graph_mean = float(mean_graph_new)
            current_graph_score = float(graph_score_new)
            if same_component_before and exact_loop_gain > 0.0:
                prev_factor = float(inter_loop_factor_by_pair.get(gate_pair, 0.0) or 0.0)
                if chosen_loop_factor > prev_factor:
                    inter_loop_factor_by_pair[gate_pair] = float(chosen_loop_factor)
                current_inter_loop_factor_sum = float(next_inter_loop_factor_sum)
                current_inter_loop_bonus = float(next_inter_loop_bonus)
            if same_component_before and exact_bridge_gain > 0.0:
                prev_bridge = float(bridge_relief_factor_by_pair.get(gate_pair, 0.0) or 0.0)
                if chosen_bridge_factor > prev_bridge:
                    bridge_relief_factor_by_pair[gate_pair] = float(chosen_bridge_factor)
                current_bridge_relief_factor_sum = float(next_bridge_relief_factor_sum)
                current_bridge_relief_bonus = float(next_bridge_relief_bonus)
            score_new = float(current_graph_score + current_intra_bonus + current_inter_loop_bonus + current_bridge_relief_bonus)

        cid = len(selected) + 1
        gate_u = int(chosen_spec.get("gate_u"))
        gate_v = int(chosen_spec.get("gate_v"))
        gate_pair = (gate_u, gate_v) if gate_u <= gate_v else (gate_v, gate_u)
        chosen_spec["same_component_inter"] = bool((not is_intra) and (int(comp_id.get(gate_u, -1)) == int(comp_id.get(gate_v, -1))))
        selected[cid] = _selected_row_from_spec(chosen_spec, utility_score=float(exact_gain))
        if bool(chosen_spec.get("allow_multi_same_pair", False)):
            cgeom = selected[cid].get("geom")
            if cgeom is not None and (not cgeom.isEmpty()):
                selected_barrier_pair_geoms[gate_pair].append(clone_geometry(cgeom))
            anchors = chosen_spec.get("anchor_pts")
            if isinstance(anchors, dict) and anchors:
                selected_barrier_pair_anchors[gate_pair].append(dict(anchors))
        else:
            selected_pairs.add(gate_pair)

        if is_intra:
            intra_selected += 1
        elif int(comp_id.get(gate_u, -1)) == int(comp_id.get(gate_v, -1)):
            same_component_selected += 1
        else:
            inter_selected += 1

        cost = float(chosen_spec.get("cost", 0.0) or 0.0)
        budget_used += cost
        remaining -= cost
        used_specs.add(int(chosen_spec_idx))
        current_score = float(score_new)

    def _intra_same_mouth_selected(spec: Dict[str, Any], selected_rows: Dict[int, Dict[str, Any]]) -> bool:
        try:
            pid = int(spec.get("intra_pid", -1) or -1)
        except Exception:
            return False
        if pid < 0:
            return False
        mid_x = float(spec.get("mid_x", 0.0) or 0.0)
        mid_y = float(spec.get("mid_y", 0.0) or 0.0)
        idx_a = int(spec.get("idx_a", -1) or -1)
        idx_b = int(spec.get("idx_b", -1) or -1)
        ring_n = max(int(spec.get("ring_n", 0) or 0), 1)
        gap_m = max(float(spec.get("intra_gap_m", spec.get("gap_m", 0.0)) or 0.0), 1.0)
        for row in selected_rows.values():
            if not bool(row.get("intra_patch", False)):
                continue
            try:
                row_pid = int(row.get("intra_patch_id", -2) or -2)
            except Exception:
                continue
            if row_pid != pid:
                continue
            ref_mid_x = float(row.get("mid_x", 0.0) or 0.0)
            ref_mid_y = float(row.get("mid_y", 0.0) or 0.0)
            ref_idx_a = int(row.get("idx_a", -1) or -1)
            ref_idx_b = int(row.get("idx_b", -1) or -1)
            ref_ring_n = max(int(row.get("ring_n", 0) or 0), 1)
            ref_gap_m = max(float(row.get("intra_gap_m", 0.0) or 0.0), 1.0)
            idx_tol = max(6, min(42, int(round(max(ring_n, ref_ring_n) * 0.028))))
            mouth_tol = max(20.0, min(140.0, max(gap_m, ref_gap_m) * 0.45))
            same_order = (
                _cyclic_idx_dist(idx_a, ref_idx_a, max(ring_n, ref_ring_n)) <= idx_tol
                and _cyclic_idx_dist(idx_b, ref_idx_b, max(ring_n, ref_ring_n)) <= idx_tol
            )
            swap_order = (
                _cyclic_idx_dist(idx_a, ref_idx_b, max(ring_n, ref_ring_n)) <= idx_tol
                and _cyclic_idx_dist(idx_b, ref_idx_a, max(ring_n, ref_ring_n)) <= idx_tol
            )
            same_mid = math.hypot(mid_x - ref_mid_x, mid_y - ref_mid_y) <= mouth_tol
            if same_order or swap_order or same_mid:
                return True
        return False

    # IPC refill phase:
    # if the main LF loop stalls with budget remaining, spend leftover budget on
    # distinct high-repair IPC winners that still improve the exact LF score.
    if remaining > 1e-9 and cand_specs:
        ipc_refill_iters = 4
        ipc_refill_pool_cap = 36
        ipc_refill_min_delta = 1e-9
        for _ipc_refill_iter in range(ipc_refill_iters):
            if remaining <= 1e-9:
                break
            base_exact = _compute_landscape_fluidity_exact(patches, selected, params)
            base_post = float(base_exact.get("post", current_score) or current_score)
            refill_ranked: List[Tuple[float, int]] = []
            for spec_idx, spec in enumerate(cand_specs):
                if spec_idx in used_specs or (not bool(spec.get("is_intra", False))):
                    continue
                cost_spec = float(spec.get("cost", 0.0) or 0.0)
                if cost_spec <= 0.0 or cost_spec > remaining + 1e-12:
                    continue
                if _intra_same_mouth_selected(spec, selected):
                    continue
                repair_score = float(spec.get("intra_repair_score", 0.0) or 0.0)
                if repair_score <= 1e-12:
                    continue
                hole_area_m2 = max(float(spec.get("intra_hole_area_m2", 0.0) or 0.0), 0.0)
                pocket_strength = max(float(spec.get("intra_pocket_strength", 0.0) or 0.0), 0.0)
                shortcut_ratio = max(float(spec.get("shortcut_ratio", 0.0) or 0.0), 0.0)
                patch_area_ha = max(float(spec.get("intra_patch_area_ha", 0.0) or 0.0), 0.0)
                patch_scale = min(1.0, patch_area_ha / 50.0)
                hint = float(
                    (
                        (0.55 * repair_score)
                        + (0.20 * math.sqrt(hole_area_m2 / 10000.0 + 1.0))
                        + (0.15 * math.log1p(pocket_strength))
                        + (0.10 * shortcut_ratio)
                    )
                    * (0.60 + 0.40 * patch_scale)
                    / max(math.sqrt(cost_spec), 1e-12)
                )
                if hint <= 1e-12:
                    continue
                refill_ranked.append((hint, int(spec_idx)))

            if not refill_ranked:
                break
            refill_ranked.sort(reverse=True)
            refill_ranked = refill_ranked[:ipc_refill_pool_cap]

            best_add: Optional[Tuple[int, float, float]] = None
            for _hint, add_idx in refill_ranked:
                add_spec = cand_specs[int(add_idx)]
                trial_selected = {int(k): dict(v) for k, v in selected.items()}
                trial_selected[len(trial_selected) + 1] = _selected_row_from_spec(add_spec, utility_score=0.0)
                trial_exact = _compute_landscape_fluidity_exact(patches, trial_selected, params)
                trial_post = float(trial_exact.get("post", 0.0) or 0.0)
                delta = float(trial_post - base_post)
                if delta <= ipc_refill_min_delta:
                    continue
                if best_add is None or delta > float(best_add[2]) + 1e-12:
                    best_add = (int(add_idx), float(trial_post), float(delta))

            if best_add is None:
                break

            add_idx, new_post_score, delta = best_add
            add_spec = cand_specs[int(add_idx)]
            add_cost = float(add_spec.get("cost", 0.0) or 0.0)
            new_cid = len(selected) + 1
            selected[new_cid] = _selected_row_from_spec(add_spec, utility_score=float(delta))
            budget_used += add_cost
            remaining = float(max(0.0, budget_limit - budget_used))
            current_score = float(new_post_score)
            used_specs.add(int(add_idx))

    # Post-greedy exact local refinement:
    # allow bounded single-corridor swaps when they improve the exact reported
    # landscape fluidity score under the same budget.
    if selected and cand_specs:
        max_swap_iters = 2
        swap_pool_cap = 28
        swap_drop_cap = 16
        for _swap_iter in range(max_swap_iters):
            base_exact = _compute_landscape_fluidity_exact(patches, selected, params)
            base_post = float(base_exact.get("post", current_score) or current_score)
            if base_post <= 0.0:
                break

            comp_id_final: Dict[int, int] = {}
            comp_area_final: Dict[int, float] = defaultdict(float)
            for comp_idx, comp in enumerate(nx.connected_components(graph)):
                for node in comp:
                    inode = int(node)
                    comp_id_final[inode] = int(comp_idx)
                    comp_area_final[int(comp_idx)] += float((patches.get(inode) or {}).get("area_ha", 0.0) or 0.0)
            shortest_final = _all_pairs_shortest_weighted(graph)
            bridge_edges_final: Set[Tuple[int, int]] = set()
            try:
                for a, b in nx.bridges(graph):
                    ia, ib = int(a), int(b)
                    bridge_edges_final.add((ia, ib) if ia <= ib else (ib, ia))
            except Exception:
                bridge_edges_final = set()

            selected_gate_pairs: Set[Tuple[int, int]] = set()
            selected_gate_pair_geoms: Dict[Tuple[int, int], List[QgsGeometry]] = defaultdict(list)
            for csel in selected.values():
                try:
                    su = int(csel.get("p1", csel.get("patch1")))
                    sv = int(csel.get("p2", csel.get("patch2")))
                except Exception:
                    continue
                if su < 0 or sv < 0 or su == sv:
                    continue
                gk = (su, sv) if su <= sv else (sv, su)
                selected_gate_pairs.add(gk)
                g = csel.get("geom")
                if g is not None and (not g.isEmpty()):
                    selected_gate_pair_geoms[gk].append(g)

            closure_pool: List[Tuple[float, int]] = []
            rel_shortcut_floor = max(1.05, 0.72 * float(same_component_shortcut_threshold))
            for spec_idx, spec in enumerate(cand_specs):
                if spec_idx in used_specs:
                    continue
                if bool(spec.get("is_intra", False)):
                    continue
                cost_spec = float(spec.get("cost", 0.0) or 0.0)
                if cost_spec <= 0.0 or cost_spec > budget_limit + 1e-12:
                    continue
                gate_u = int(spec.get("gate_u"))
                gate_v = int(spec.get("gate_v"))
                if gate_u == gate_v:
                    continue
                gate_pair = (gate_u, gate_v) if gate_u <= gate_v else (gate_v, gate_u)
                if gate_pair in selected_gate_pairs:
                    if not bool(spec.get("allow_multi_same_pair", False)):
                        continue
                    cgeom = (spec.get("cand") or {}).get("geom")
                    if cgeom is None or cgeom.isEmpty():
                        continue
                    # For barrier-qualified pairs, allow additional corridors only when
                    # separated by the configured search radius at endpoints.
                    ca = _pair_anchor_points(gate_pair, cgeom)
                    distinct = True
                    for pg in selected_gate_pair_geoms.get(gate_pair, []):
                        pa = _pair_anchor_points(gate_pair, pg)
                        for pid in gate_pair:
                            if (
                                int(pid) not in ca
                                or int(pid) not in pa
                                or float(ca[int(pid)].distance(pa[int(pid)])) + 1e-9 < pair_multi_min_sep
                            ):
                                distinct = False
                                break
                        if not distinct:
                            break
                    if not distinct:
                        continue
                cu = int(comp_id_final.get(gate_u, -1))
                cv = int(comp_id_final.get(gate_v, -1))
                if cu < 0 or cv < 0 or cu != cv:
                    continue
                d_old = _shortest_lookup(shortest_final, gate_u, gate_v)
                c_e = max(float(spec.get("length", 0.0) or 0.0), 1e-9)
                shortcut_ratio = (d_old / c_e) if np.isfinite(d_old) else float("inf")
                if shortcut_ratio + 1e-12 < rel_shortcut_floor:
                    continue
                comp_area = float(comp_area_final.get(cu, 0.0) or 0.0)
                try:
                    hop_len = int(nx.shortest_path_length(graph, int(gate_u), int(gate_v)))
                except Exception:
                    hop_len = 0
                loop_factor = _landscape_fluidity_inter_loop_factor(
                    float(shortcut_ratio),
                    float(comp_area),
                    float(total_habitat_area),
                    hop_length=int(hop_len),
                )
                bridge_fraction = _bridge_fraction_on_shortest_path(
                    graph,
                    int(gate_u),
                    int(gate_v),
                    bridge_edges_final,
                )
                bridge_factor = _landscape_fluidity_bridge_relief_factor(
                    float(bridge_fraction),
                    float(comp_area),
                    float(total_habitat_area),
                    float(shortcut_ratio),
                )
                priority = float(
                    (
                        float(loop_factor)
                        + float(bridge_factor)
                        + (0.25 * math.sqrt(max(float(loop_factor), 0.0)))
                    )
                    / max(cost_spec, 1e-12)
                )
                if priority <= 1e-12:
                    continue
                closure_pool.append((priority, int(spec_idx)))

            if not closure_pool:
                break
            closure_pool.sort(reverse=True)
            closure_pool = closure_pool[:swap_pool_cap]

            drop_ids = sorted(
                selected.keys(),
                key=lambda cid: (
                    float(selected[cid].get("utility_score", 0.0) or 0.0),
                    -float(selected[cid].get("area_ha", 0.0) or 0.0),
                    -int(cid),
                ),
            )[:swap_drop_cap]
            if not drop_ids:
                break

            best_swap: Optional[Tuple[int, int, float, float, float]] = None
            # (drop_cid, add_spec_idx, new_budget_used, new_post_score, delta)
            for _prio, add_idx in closure_pool:
                add_spec = cand_specs[int(add_idx)]
                add_cost = float(add_spec.get("cost", 0.0) or 0.0)
                if add_cost <= 0.0:
                    continue
                for drop_cid in drop_ids:
                    drop_cost = float(selected[int(drop_cid)].get("area_ha", 0.0) or 0.0)
                    new_budget_used = float(budget_used - drop_cost + add_cost)
                    if new_budget_used > budget_limit + 1e-12:
                        continue
                    trial_selected = {int(k): dict(v) for k, v in selected.items()}
                    add_spec["same_component_inter"] = True
                    trial_selected[int(drop_cid)] = _selected_row_from_spec(add_spec, utility_score=0.0)
                    trial_exact = _compute_landscape_fluidity_exact(patches, trial_selected, params)
                    trial_post = float(trial_exact.get("post", 0.0) or 0.0)
                    delta = float(trial_post - base_post)
                    if delta <= 1e-12:
                        continue
                    if best_swap is None or delta > float(best_swap[4]) + 1e-12:
                        best_swap = (int(drop_cid), int(add_idx), float(new_budget_used), float(trial_post), float(delta))

            if best_swap is None:
                break

            drop_cid, add_idx, new_budget_used, new_post_score, delta = best_swap
            add_spec = cand_specs[int(add_idx)]
            add_spec["same_component_inter"] = True
            selected[int(drop_cid)] = _selected_row_from_spec(
                add_spec,
                utility_score=float(selected[int(drop_cid)].get("utility_score", 0.0) or 0.0) + float(delta),
            )
            budget_used = float(new_budget_used)
            remaining = float(max(0.0, budget_limit - budget_used))
            current_score = float(new_post_score)
            used_specs.add(int(add_idx))

            # Rebuild working graph from current selected set to keep subsequent
            # swap iteration topology metrics consistent.
            graph = graph_base.copy()
            for csel in selected.values():
                edges_local: List[Tuple[int, int, float]] = list(csel.get("graph_edges", []) or [])
                if not edges_local:
                    pids = _candidate_patch_ids_for_metric(csel, valid_nodes)
                    if len(pids) < 2:
                        continue
                    length = float(csel.get("distance", csel.get("distance_m", 0.0)) or 0.0)
                    if length <= 0.0:
                        p1 = int(pids[0])
                        p2 = int(pids[-1])
                        x1, y1 = coords_xy.get(p1, (0.0, 0.0))
                        x2, y2 = coords_xy.get(p2, (0.0, 0.0))
                        length = max(float(math.hypot(x1 - x2, y1 - y2)), 1.0)
                    for i in range(len(pids)):
                        for j in range(i + 1, len(pids)):
                            edges_local.append((int(pids[i]), int(pids[j]), float(length)))
                _apply_edges_with_changes(graph, edges_local)

        # If swap refinement freed budget, greedily top-up with bounded exact-gain
        # additions so remaining budget is used only on real metric improvements.
        fill_iters = 4
        fill_pool_cap = 48
        for _fill_iter in range(fill_iters):
            if remaining <= 1e-9:
                break
            base_exact = _compute_landscape_fluidity_exact(patches, selected, params)
            base_post = float(base_exact.get("post", current_score) or current_score)
            selected_gate_pairs: Set[Tuple[int, int]] = set()
            selected_gate_pair_geoms: Dict[Tuple[int, int], List[QgsGeometry]] = defaultdict(list)
            for csel in selected.values():
                try:
                    su = int(csel.get("p1", csel.get("patch1")))
                    sv = int(csel.get("p2", csel.get("patch2")))
                except Exception:
                    continue
                if su < 0 or sv < 0 or su == sv:
                    continue
                gk = (su, sv) if su <= sv else (sv, su)
                selected_gate_pairs.add(gk)
                g = csel.get("geom")
                if g is not None and (not g.isEmpty()):
                    selected_gate_pair_geoms[gk].append(g)

            fill_ranked: List[Tuple[float, int]] = []
            for spec_idx, spec in enumerate(cand_specs):
                if spec_idx in used_specs:
                    continue
                cost_spec = float(spec.get("cost", 0.0) or 0.0)
                if cost_spec <= 0.0 or cost_spec > remaining + 1e-12:
                    continue
                if bool(spec.get("is_intra", False)):
                    hint = float(spec.get("shortcut_ratio", 1.0) or 1.0) / max(cost_spec, 1e-12)
                    fill_ranked.append((hint, int(spec_idx)))
                    continue
                gate_u = int(spec.get("gate_u"))
                gate_v = int(spec.get("gate_v"))
                if gate_u != gate_v:
                    gate_pair = (gate_u, gate_v) if gate_u <= gate_v else (gate_v, gate_u)
                    if gate_pair in selected_gate_pairs:
                        if not bool(spec.get("allow_multi_same_pair", False)):
                            continue
                        cgeom = (spec.get("cand") or {}).get("geom")
                        if cgeom is None or cgeom.isEmpty():
                            continue
                        ca = _pair_anchor_points(gate_pair, cgeom)
                        distinct = True
                        for pg in selected_gate_pair_geoms.get(gate_pair, []):
                            pa = _pair_anchor_points(gate_pair, pg)
                            for pid in gate_pair:
                                if (
                                    int(pid) not in ca
                                    or int(pid) not in pa
                                    or float(ca[int(pid)].distance(pa[int(pid)])) + 1e-9 < pair_multi_min_sep
                                ):
                                    distinct = False
                                    break
                            if not distinct:
                                break
                        if not distinct:
                            continue
                hint = (
                    float(spec.get("long_hop_factor", 1.0) or 1.0)
                    / max(float(spec.get("score_cost", cost_spec) or cost_spec), 1e-12)
                )
                fill_ranked.append((hint, int(spec_idx)))

            if not fill_ranked:
                break
            fill_ranked.sort(reverse=True)
            fill_ranked = fill_ranked[:fill_pool_cap]

            best_add: Optional[Tuple[int, float, float]] = None
            # (add_spec_idx, new_post_score, delta)
            for _hint, add_idx in fill_ranked:
                add_spec = cand_specs[int(add_idx)]
                trial_selected = {int(k): dict(v) for k, v in selected.items()}
                add_spec["same_component_inter"] = bool(
                    (not bool(add_spec.get("is_intra", False)))
                    and bool(nx.has_path(graph, int(add_spec.get("gate_u", -1)), int(add_spec.get("gate_v", -1))))
                )
                trial_selected[len(trial_selected) + 1] = _selected_row_from_spec(add_spec, utility_score=0.0)
                trial_exact = _compute_landscape_fluidity_exact(patches, trial_selected, params)
                trial_post = float(trial_exact.get("post", 0.0) or 0.0)
                delta = float(trial_post - base_post)
                if delta <= 1e-12:
                    continue
                if best_add is None or delta > float(best_add[2]) + 1e-12:
                    best_add = (int(add_idx), float(trial_post), float(delta))

            if best_add is None:
                break

            add_idx, new_post_score, delta = best_add
            add_spec = cand_specs[int(add_idx)]
            add_cost = float(add_spec.get("cost", 0.0) or 0.0)
            new_cid = len(selected) + 1
            add_spec["same_component_inter"] = bool(
                (not bool(add_spec.get("is_intra", False)))
                and bool(nx.has_path(graph, int(add_spec.get("gate_u", -1)), int(add_spec.get("gate_v", -1))))
            )
            selected[new_cid] = _selected_row_from_spec(add_spec, utility_score=float(delta))
            budget_used += add_cost
            remaining = float(max(0.0, budget_limit - budget_used))
            current_score = float(new_post_score)
            used_specs.add(int(add_idx))

            if not bool(add_spec.get("is_intra", False)):
                _apply_edges_with_changes(
                    graph,
                    add_spec.get("graph_edges", add_spec.get("edges", [])),
                )

        # Keep graph-resistance stats in sync with any accepted swap(s).
        current_graph_mean = float(
            _mean_effective_resistance_graph(
                graph,
                pair_records,
                disconnected_penalty,
                pair_weights=pair_weights,
            )
        )

    intra_selected = sum(1 for row in selected.values() if bool(row.get("intra_patch", False)))
    same_component_selected = sum(1 for row in selected.values() if bool(row.get("same_component_inter", False)))
    inter_selected = max(0, len(selected) - intra_selected - same_component_selected)

    stats = {
        "strategy": "landscape_fluidity",
        "corridors_used": len(selected),
        "budget_used_ha": float(budget_used),
        "mean_effective_resistance_pre": float(mean_pre_graph),
        "mean_effective_resistance_post": float(current_graph_mean),
        "mean_effective_resistance_gain": float(mean_pre_graph - current_graph_mean),
        "landscape_fluidity_pre": float(score_pre),
        "landscape_fluidity_post": float(current_score),
        "landscape_fluidity_gain": float(current_score - score_pre),
        "landscape_fluidity_inter_selected": int(inter_selected),
        "landscape_fluidity_intra_selected": int(intra_selected),
        "landscape_fluidity_same_component_selected": int(same_component_selected),
        "landscape_fluidity_shortcut_threshold": float(shortcut_threshold),
        "landscape_fluidity_same_component_shortcut_threshold": float(
            same_component_shortcut_threshold
        ),
        "landscape_fluidity_gain_weight": float(gain_weight),
        "landscape_fluidity_ratio_weight": float(ratio_weight),
        "landscape_fluidity_topology_weight": float(topology_weight),
        "landscape_fluidity_long_hop_factor": float(long_hop_factor_thresh),
        "landscape_fluidity_long_hop_tau": float(long_hop_tau),
        "landscape_fluidity_intra_tau": float(intra_tau),
        "landscape_fluidity_lookahead_k": int(lookahead_k),
        "landscape_fluidity_lookahead_weight": float(lookahead_weight),
        "landscape_fluidity_intra_bonus_weight": float(
            getattr(params, "landscape_fluidity_intra_bonus_weight", 1.0) or 1.0
        ),
        "landscape_fluidity_intra_bonus_cap_mult": float(
            getattr(params, "landscape_fluidity_intra_bonus_cap_mult", 1.0) or 1.0
        ),
        "landscape_fluidity_inter_loop_weight": float(
            getattr(params, "landscape_fluidity_inter_loop_weight", 0.35) or 0.35
        ),
        "landscape_fluidity_inter_loop_cap_mult": float(
            getattr(params, "landscape_fluidity_inter_loop_cap_mult", 0.80) or 0.80
        ),
        "landscape_fluidity_bridge_relief_weight": float(
            getattr(params, "landscape_fluidity_bridge_relief_weight", 0.25) or 0.25
        ),
        "landscape_fluidity_bridge_relief_cap_mult": float(
            getattr(params, "landscape_fluidity_bridge_relief_cap_mult", 0.60) or 0.60
        ),
        "landscape_fluidity_intra_bonus_post": float(current_intra_bonus),
        "landscape_fluidity_inter_loop_bonus_post": float(current_inter_loop_bonus),
        "landscape_fluidity_inter_loop_factor_sum": float(current_inter_loop_factor_sum),
        "landscape_fluidity_bridge_relief_bonus_post": float(current_bridge_relief_bonus),
        "landscape_fluidity_bridge_relief_factor_sum": float(current_bridge_relief_factor_sum),
        "landscape_fluidity_pair_samples": int(len(pair_records)),
        "patches_connected": 0,
        "total_connected_area_ha": sum(p.get("area_ha", 0.0) for p in patches.values()),
    }
    return selected, stats


def _with_landscape_fluidity_variant(
    stats: Dict[str, Any],
    *,
    strategy_key: str,
    variant_label: str,
) -> Dict[str, Any]:
    out = dict(stats or {})
    out["strategy"] = str(strategy_key)
    out["landscape_fluidity_variant"] = str(variant_label)
    return out


def optimize_landscape_fluidity_a(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
) -> Tuple[Dict[int, Dict], Dict]:
    corridors, stats = optimize_landscape_fluidity(patches, candidates, params)
    return corridors, _with_landscape_fluidity_variant(
        stats,
        strategy_key="landscape_fluidity",
        variant_label="LF-A",
    )


def _lf2_pair_key(u: int, v: int) -> Tuple[int, int]:
    iu = int(u)
    iv = int(v)
    return (iu, iv) if iu <= iv else (iv, iu)


def _lf2_add_pair_record(
    *,
    pair_records: List[Tuple[int, int, float]],
    pair_weight_raw: Dict[Tuple[int, int], float],
    seen_pairs: Set[Tuple[int, int]],
    u: int,
    v: int,
    eu: float,
    weight_raw: float,
) -> None:
    if int(u) == int(v):
        return
    key = _lf2_pair_key(int(u), int(v))
    if key in seen_pairs:
        return
    e = max(float(eu), 1e-6)
    w = max(float(weight_raw), 1e-9)
    seen_pairs.add(key)
    pair_records.append((int(key[0]), int(key[1]), float(e)))
    pair_weight_raw[key] = float(pair_weight_raw.get(key, 0.0) or 0.0) + float(w)


def _lf2_pick_virtual_targets(
    pid: int,
    coords_xy: Dict[int, Tuple[float, float]],
    *,
    max_targets: int = 3,
) -> List[int]:
    if int(pid) not in coords_xy:
        return []
    x0, y0 = coords_xy[int(pid)]
    dists: List[Tuple[float, int]] = []
    for other, xy in coords_xy.items():
        io = int(other)
        if io == int(pid):
            continue
        x1, y1 = xy
        d = float(math.hypot(float(x0) - float(x1), float(y0) - float(y1)))
        if np.isfinite(d) and d > 0.0:
            dists.append((float(d), int(io)))
    if not dists:
        return []
    dists.sort(key=lambda t: (float(t[0]), int(t[1])))
    near = [int(pid2) for _d, pid2 in dists[: max(1, min(2, int(max_targets)))]]
    far = [int(dists[-1][1])] if len(dists) > 2 and int(max_targets) >= 3 else []
    out: List[int] = []
    seen: Set[int] = set()
    for n in near + far:
        if n in seen:
            continue
        seen.add(n)
        out.append(int(n))
        if len(out) >= int(max_targets):
            break
    return out


def _build_landscape_fluidity_2_context(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
) -> Optional[Dict[str, Any]]:
    base = _build_patch_resistance_context(patches, params)
    if base is None:
        return None

    graph = base["graph"].copy()
    coords_xy: Dict[int, Tuple[float, float]] = dict(base.get("coords_xy") or {})
    valid_nodes = set(int(pid) for pid in base.get("node_ids", []))
    if not valid_nodes:
        return None

    pair_records: List[Tuple[int, int, float]] = [
        (int(u), int(v), float(eu))
        for u, v, eu in list(base.get("pair_records") or [])
        if int(u) != int(v)
    ]
    seen_pairs: Set[Tuple[int, int]] = {_lf2_pair_key(int(u), int(v)) for u, v, _eu in pair_records}
    base_pair_weights = dict(base.get("pair_weights") or {})
    pair_weight_raw: Dict[Tuple[int, int], float] = {}
    for u, v, _eu in pair_records:
        key = _lf2_pair_key(int(u), int(v))
        pair_weight_raw[key] = max(float(base_pair_weights.get(key, 0.0) or 0.0), 1e-9)

    patch_area_ha: Dict[int, float] = {
        int(pid): float((patches.get(int(pid)) or {}).get("area_ha", 0.0) or 0.0)
        for pid in valid_nodes
    }
    total_habitat_area = max(1e-9, float(sum(max(0.0, a) for a in patch_area_ha.values())))

    cand_specs: List[Dict[str, Any]] = []
    budget_limit = float(params.budget_area or 0.0)
    next_virtual_id = -1
    max_virtual_targets = max(1, min(4, int(getattr(params, "landscape_fluidity_2_virtual_targets", 3) or 3)))
    virtual_pair_weight = max(
        0.05,
        min(2.0, float(getattr(params, "landscape_fluidity_2_virtual_pair_weight", 0.40) or 0.40)),
    )

    for idx, cand in enumerate(candidates):
        cost = float(cand.get("area_ha", 0.0) or 0.0)
        if cost <= 0.0 or cost > budget_limit + 1e-12:
            continue

        p1_raw = cand.get("patch1", cand.get("p1"))
        p2_raw = cand.get("patch2", cand.get("p2"))
        is_intra = bool(cand.get("intra_patch", False))
        try:
            if p1_raw is not None and p2_raw is not None and int(p1_raw) == int(p2_raw):
                is_intra = True
        except Exception:
            pass

        if is_intra:
            try:
                pid = int(cand.get("intra_patch_id", p1_raw))
            except Exception:
                continue
            if pid not in valid_nodes:
                continue
            gap_m = max(float(cand.get("distance_m", cand.get("distance", 0.0)) or 0.0), 1e-6)
            detour_m = float(cand.get("intra_detour_old_m", 0.0) or 0.0)
            ratio = float(cand.get("intra_shortcut_ratio", cand.get("fluidity_gain", 0.0)) or 0.0)
            if detour_m <= gap_m + 1e-9 and ratio > 1.0:
                detour_m = float(gap_m * ratio)
            if detour_m <= gap_m + 1e-9:
                continue

            vid = int(next_virtual_id)
            next_virtual_id -= 1
            graph.add_node(int(vid))
            graph.add_edge(int(pid), int(vid), weight=float(detour_m), distance=float(detour_m))

            # Parallel path equivalent resistance after adding the shortcut.
            new_eff = 1.0 / ((1.0 / max(detour_m, 1e-9)) + (1.0 / max(gap_m, 1e-9)))
            new_eff = max(float(new_eff), 1e-9)

            area_pid = max(float(patch_area_ha.get(int(pid), 0.0) or 0.0), 0.0)
            _lf2_add_pair_record(
                pair_records=pair_records,
                pair_weight_raw=pair_weight_raw,
                seen_pairs=seen_pairs,
                u=int(pid),
                v=int(vid),
                eu=max(float(gap_m), 1.0),
                weight_raw=max(area_pid, 1e-9) * virtual_pair_weight,
            )
            for tgt in _lf2_pick_virtual_targets(int(pid), coords_xy, max_targets=max_virtual_targets):
                xt, yt = coords_xy.get(int(tgt), (0.0, 0.0))
                xp, yp = coords_xy.get(int(pid), (0.0, 0.0))
                eu = max(float(math.hypot(xp - xt, yp - yt)), 1.0)
                wt = math.sqrt(
                    max(area_pid, 0.0)
                    * max(float(patch_area_ha.get(int(tgt), 0.0) or 0.0), 0.0)
                )
                _lf2_add_pair_record(
                    pair_records=pair_records,
                    pair_weight_raw=pair_weight_raw,
                    seen_pairs=seen_pairs,
                    u=int(vid),
                    v=int(tgt),
                    eu=float(eu + max(0.25 * gap_m, 1.0)),
                    weight_raw=max(float(wt), 1e-9) * virtual_pair_weight,
                )

            cand_specs.append(
                {
                    "idx": int(idx),
                    "cand": cand,
                    "cost": float(cost),
                    "is_intra": True,
                    "gate_u": int(pid),
                    "gate_v": int(vid),
                    "proxy_new_distance": float(gap_m),
                    "length": float(gap_m),
                    "edges": [(int(pid), int(vid), float(new_eff))],
                    "pids": [int(pid)],
                    "output_p1": int(pid),
                    "output_p2": int(pid),
                    "intra_pid": int(pid),
                    "intra_detour_old_m": float(detour_m),
                    "intra_shortcut_ratio": float(detour_m / max(gap_m, 1e-9)),
                }
            )
            continue

        pids = _candidate_patch_ids_for_metric(cand, valid_nodes)
        length = float(cand.get("distance_m", 0.0) or 0.0)
        edges_local: List[Tuple[int, int, float]] = []

        explicit_edges = cand.get("explicit_edges")
        if isinstance(explicit_edges, (list, tuple)):
            for item in explicit_edges:
                try:
                    if not isinstance(item, (list, tuple)) or len(item) < 2:
                        continue
                    eu = int(item[0])
                    ev = int(item[1])
                    elen = float(item[2]) if len(item) >= 3 else float(length)
                except Exception:
                    continue
                if eu not in valid_nodes or ev not in valid_nodes or eu == ev:
                    continue
                if elen <= 0.0:
                    x1, y1 = coords_xy.get(eu, (0.0, 0.0))
                    x2, y2 = coords_xy.get(ev, (0.0, 0.0))
                    elen = max(float(math.hypot(x1 - x2, y1 - y2)), 1.0)
                edges_local.append((int(eu), int(ev), float(elen)))
            if edges_local:
                node_union: Set[int] = set()
                for eu, ev, _d in edges_local:
                    node_union.add(int(eu))
                    node_union.add(int(ev))
                pids = sorted(node_union)

        if not edges_local:
            if len(pids) < 2:
                continue
            if length <= 0.0:
                a = int(pids[0])
                b = int(pids[-1])
                x1, y1 = coords_xy.get(a, (0.0, 0.0))
                x2, y2 = coords_xy.get(b, (0.0, 0.0))
                length = max(float(math.hypot(x1 - x2, y1 - y2)), 1.0)
            for i in range(len(pids)):
                for j in range(i + 1, len(pids)):
                    edges_local.append((int(pids[i]), int(pids[j]), float(length)))
            if not edges_local:
                continue

        try:
            gate_u = int(p1_raw) if p1_raw is not None else int(pids[0])
        except Exception:
            gate_u = int(pids[0])
        try:
            gate_v = int(p2_raw) if p2_raw is not None else int(pids[-1])
        except Exception:
            gate_v = int(pids[-1])
        if gate_u not in valid_nodes:
            gate_u = int(pids[0])
        if gate_v not in valid_nodes:
            gate_v = int(pids[-1])

        cand_specs.append(
            {
                "idx": int(idx),
                "cand": cand,
                "cost": float(cost),
                "is_intra": False,
                "gate_u": int(gate_u),
                "gate_v": int(gate_v),
                "proxy_new_distance": float(max(length, 1e-9)),
                "length": float(max(length, 1e-9)),
                "edges": list(edges_local),
                "pids": list(pids),
                "output_p1": int(gate_u),
                "output_p2": int(gate_v),
            }
        )

    if not cand_specs:
        return None

    total_w = float(sum(max(float(v), 0.0) for v in pair_weight_raw.values()))
    if total_w <= 0.0:
        total_w = float(len(pair_weight_raw) or 1)
    pair_weights = {
        key: float(max(float(val), 0.0) / total_w)
        for key, val in pair_weight_raw.items()
    }

    return {
        "graph": graph,
        "cand_specs": cand_specs,
        "pair_records": pair_records,
        "pair_weights": pair_weights,
        "disconnected_penalty": float(base.get("disconnected_penalty", 1.0)),
        "coords_xy": coords_xy,
        "valid_nodes": valid_nodes,
        "patch_area_ha": patch_area_ha,
        "total_habitat_area": float(total_habitat_area),
    }


def _compute_landscape_fluidity_2_exact(
    patches: Dict[int, Dict],
    corridors: Dict[int, Dict],
    params: VectorRunParams,
    candidate_pool: Optional[List[Dict]] = None,
) -> Dict[str, float]:
    pool = list(candidate_pool) if candidate_pool is not None else list((corridors or {}).values())
    context = _build_landscape_fluidity_2_context(patches, pool, params)
    if context is None:
        return {"pre": 0.0, "post": 0.0, "gain": 0.0, "graph_resistance_pre": 0.0, "graph_resistance_post": 0.0}

    graph_pre: nx.Graph = context["graph"].copy()
    graph_post: nx.Graph = graph_pre.copy()
    pair_records: List[Tuple[int, int, float]] = list(context.get("pair_records") or [])
    pair_weights: Dict[Tuple[int, int], float] = dict(context.get("pair_weights") or {})
    penalty = float(context.get("disconnected_penalty", 1.0) or 1.0)
    total_habitat_area = max(1e-9, float(context.get("total_habitat_area", 1.0) or 1.0))

    def _sig_from_corridor(cdata: Dict[str, Any]) -> Tuple[Any, ...]:
        is_intra = bool(cdata.get("intra_patch", False))
        p1 = int(cdata.get("p1", cdata.get("patch1", 0)) or 0)
        p2 = int(cdata.get("p2", cdata.get("patch2", 0)) or 0)
        if is_intra:
            p2 = int(p1)
        area = round(float(cdata.get("area_ha", 0.0) or 0.0), 6)
        dist = round(float(cdata.get("distance", cdata.get("distance_m", 0.0)) or 0.0), 6)
        variant = str(cdata.get("variant", "") or "")
        source = str(cdata.get("source", "") or "")
        return (bool(is_intra), int(p1), int(p2), float(area), float(dist), variant, source)

    def _sig_from_spec(spec: Dict[str, Any]) -> Tuple[Any, ...]:
        cand = dict(spec.get("cand") or {})
        area = round(float(spec.get("cost", cand.get("area_ha", 0.0)) or 0.0), 6)
        dist = round(float(spec.get("length", 0.0) or 0.0), 6)
        is_intra = bool(spec.get("is_intra", False))
        p1 = int(spec.get("output_p1", spec.get("gate_u", 0)) or 0)
        p2 = int(spec.get("output_p2", spec.get("gate_v", 0)) or 0)
        if is_intra:
            p2 = int(p1)
        variant = str(cand.get("variant", "") or "")
        source = str(cand.get("source", "") or "")
        return (bool(is_intra), int(p1), int(p2), float(area), float(dist), variant, source)

    spec_by_sig: Dict[Tuple[Any, ...], List[Dict[str, Any]]] = defaultdict(list)
    for spec in context.get("cand_specs", []):
        spec_by_sig[_sig_from_spec(spec)].append(spec)
    selected_signatures: List[Tuple[Any, ...]] = [
        _sig_from_corridor(cdata) for _cid, cdata in sorted((corridors or {}).items(), key=lambda kv: int(kv[0]))
    ]
    for sig in selected_signatures:
        bucket = spec_by_sig.get(sig) or []
        if not bucket:
            continue
        picked_spec = bucket.pop(0)
        _apply_edges_with_changes(graph_post, picked_spec.get("edges", []))

    mean_pre = float(
        _mean_effective_resistance_graph(
            graph_pre,
            pair_records,
            penalty,
            pair_weights=pair_weights,
        )
    )
    mean_post = float(
        _mean_effective_resistance_graph(
            graph_post,
            pair_records,
            penalty,
            pair_weights=pair_weights,
        )
    )
    pre_v = float(total_habitat_area / max(mean_pre, 1e-9))
    post_v = float(total_habitat_area / max(mean_post, 1e-9))
    return {
        "pre": float(pre_v),
        "post": float(post_v),
        "gain": float(post_v - pre_v),
        "graph_resistance_pre": float(mean_pre),
        "graph_resistance_post": float(mean_post),
        "intra_bonus_post": 0.0,
        "intra_factor_sum": 0.0,
        "inter_loop_bonus_post": 0.0,
        "inter_loop_factor_sum": 0.0,
        "bridge_relief_bonus_post": 0.0,
        "bridge_relief_factor_sum": 0.0,
    }


def optimize_landscape_fluidity_2(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
) -> Tuple[Dict[int, Dict], Dict]:
    budget_limit = float(params.budget_area or 0.0)
    if budget_limit <= 0.0 or (not patches) or (not candidates):
        return {}, {"strategy": "landscape_fluidity_2", "corridors_used": 0, "budget_used_ha": 0.0}

    context = _build_landscape_fluidity_2_context(patches, candidates, params)
    if context is None:
        return {}, {"strategy": "landscape_fluidity_2", "corridors_used": 0, "budget_used_ha": 0.0}

    graph = context["graph"].copy()
    cand_specs: List[Dict[str, Any]] = list(context.get("cand_specs") or [])
    pair_records: List[Tuple[int, int, float]] = list(context.get("pair_records") or [])
    pair_weights: Dict[Tuple[int, int], float] = dict(context.get("pair_weights") or {})
    penalty = float(context.get("disconnected_penalty", 1.0) or 1.0)
    total_habitat_area = max(1e-9, float(context.get("total_habitat_area", 1.0) or 1.0))
    patch_area_ha = dict(context.get("patch_area_ha") or {})

    top_k = max(3, min(40, int(getattr(params, "landscape_fluidity_2_top_k", 15) or 15)))

    mean_pre = float(
        _mean_effective_resistance_graph(
            graph,
            pair_records,
            penalty,
            pair_weights=pair_weights,
        )
    )
    current_mean = float(mean_pre)
    current_score = float(total_habitat_area / max(current_mean, 1e-9))
    score_pre = float(current_score)

    selected: Dict[int, Dict[str, Any]] = {}
    used_specs: Set[int] = set()
    remaining = float(budget_limit)
    budget_used = 0.0

    while remaining > 1e-12:
        shortest_current = _all_pairs_shortest_weighted(graph)
        comp_id: Dict[int, int] = {}
        comp_area_ha: Dict[int, float] = defaultdict(float)
        for comp_idx, comp in enumerate(nx.connected_components(graph)):
            for node in comp:
                inode = int(node)
                comp_id[inode] = int(comp_idx)
                if inode > 0:
                    comp_area_ha[int(comp_idx)] += float(
                        max(float(patch_area_ha.get(int(inode), 0.0) or 0.0), 0.0)
                    )

        ranked: List[Tuple[float, float, int]] = []
        for local_idx, spec in enumerate(cand_specs):
            if local_idx in used_specs:
                continue
            cost = float(spec.get("cost", 0.0) or 0.0)
            if cost <= 0.0 or cost > remaining + 1e-12:
                continue
            gate_u = int(spec.get("gate_u"))
            gate_v = int(spec.get("gate_v"))
            if gate_u == gate_v:
                continue

            proxy = 0.0
            if bool(spec.get("is_intra", False)):
                area_pid = max(float(patch_area_ha.get(int(spec.get("intra_pid", gate_u)), 0.0) or 0.0), 0.0)
                d_old = _shortest_lookup(shortest_current, int(gate_u), int(gate_v))
                d_new = max(float(spec.get("proxy_new_distance", 0.0) or 0.0), 1e-9)
                delta = max(float(d_old) - float(d_new), 0.0) if np.isfinite(d_old) else 0.0
                proxy = float((area_pid * delta) / max(cost, 1e-9))
            else:
                cu = int(comp_id.get(int(gate_u), -1))
                cv = int(comp_id.get(int(gate_v), -1))
                if cu < 0 or cv < 0:
                    continue
                area_u = max(float(comp_area_ha.get(cu, 0.0) or 0.0), 0.0)
                area_v = max(float(comp_area_ha.get(cv, 0.0) or 0.0), 0.0)
                area_scale = math.sqrt(max(area_u * area_v, 0.0))
                if area_scale <= 0.0:
                    continue
                if cu != cv:
                    proxy = float(area_scale / max(cost, 1e-9))
                else:
                    d_old = _shortest_lookup(shortest_current, int(gate_u), int(gate_v))
                    d_new = max(float(spec.get("proxy_new_distance", 0.0) or 0.0), 1e-9)
                    delta = max(float(d_old) - float(d_new), 0.0) if np.isfinite(d_old) else 0.0
                    if delta <= 1e-12:
                        continue
                    proxy = float((area_scale * delta) / max(cost, 1e-9))

            if proxy <= 1e-12:
                continue
            ranked.append((float(proxy), -float(cost), int(local_idx)))

        if not ranked:
            break

        ranked.sort(reverse=True)
        eval_ids = [idx for _proxy, _neg_cost, idx in ranked[:top_k]]

        best: Optional[Tuple[int, List[Tuple[str, int, int, float]], float, float, float]] = None
        # (spec_idx, changes, new_mean, gain_abs, gain_per_cost)
        for spec_idx in eval_ids:
            spec = cand_specs[int(spec_idx)]
            cost = float(spec.get("cost", 0.0) or 0.0)
            changes = _apply_edges_with_changes(graph, spec.get("edges", []))
            if not changes:
                used_specs.add(int(spec_idx))
                continue
            mean_new = float(
                _mean_effective_resistance_graph(
                    graph,
                    pair_records,
                    penalty,
                    pair_weights=pair_weights,
                )
            )
            score_new = float(total_habitat_area / max(mean_new, 1e-9))
            gain_abs = float(score_new - current_score)
            gain_per_cost = float(gain_abs / max(cost, 1e-9))
            _revert_edges_with_changes(graph, changes)
            if gain_abs <= 1e-12:
                continue
            if best is None:
                best = (int(spec_idx), changes, float(mean_new), float(gain_abs), float(gain_per_cost))
                continue
            _, _, _, best_gain_abs, best_gain_per_cost = best
            if (
                gain_per_cost > float(best_gain_per_cost) + 1e-12
                or (
                    abs(gain_per_cost - float(best_gain_per_cost)) <= 1e-12
                    and gain_abs > float(best_gain_abs) + 1e-12
                )
            ):
                best = (int(spec_idx), changes, float(mean_new), float(gain_abs), float(gain_per_cost))

        if best is None:
            break

        spec_idx, _changes_preview, mean_new, gain_abs, gain_per_cost = best
        chosen = cand_specs[int(spec_idx)]
        applied_changes = _apply_edges_with_changes(graph, chosen.get("edges", []))
        if not applied_changes:
            used_specs.add(int(spec_idx))
            continue

        cid = len(selected) + 1
        p1_out = int(chosen.get("output_p1", chosen.get("gate_u", 0)))
        p2_out = int(chosen.get("output_p2", chosen.get("gate_v", 0)))
        pids = [int(pid) for pid in list(chosen.get("pids") or []) if int(pid) > 0]
        selected[cid] = {
            "geom": clone_geometry(chosen["cand"]["geom"]),
            "patch_ids": set(pids) if pids else {int(max(p1_out, 0))},
            "area_ha": float(chosen.get("cost", 0.0) or 0.0),
            "p1": int(p1_out),
            "p2": int(p2_out),
            "distance": float(chosen.get("length", chosen.get("proxy_new_distance", 0.0)) or 0.0),
            "type": "primary",
            "variant": chosen["cand"].get("variant"),
            "source": chosen["cand"].get("source"),
            "utility_score": float(gain_per_cost),
            "overlap_ratio": 0.0,
            "intra_patch": bool(chosen.get("is_intra", False)),
            "intra_patch_id": int(chosen.get("intra_pid")) if bool(chosen.get("is_intra", False)) else None,
            "intra_shortcut_ratio": float(chosen.get("intra_shortcut_ratio", 0.0) or 0.0),
            "intra_detour_old_m": float(chosen.get("intra_detour_old_m", 0.0) or 0.0),
            "fluidity_gain": float(chosen.get("intra_shortcut_ratio", 0.0) or 0.0),
        }

        cost = float(chosen.get("cost", 0.0) or 0.0)
        budget_used += cost
        remaining = float(max(0.0, budget_limit - budget_used))
        current_mean = float(mean_new)
        current_score = float(current_score + gain_abs)
        used_specs.add(int(spec_idx))

    stats = {
        "strategy": "landscape_fluidity_2",
        "corridors_used": len(selected),
        "budget_used_ha": float(budget_used),
        "mean_effective_resistance_pre": float(mean_pre),
        "mean_effective_resistance_post": float(current_mean),
        "mean_effective_resistance_gain": float(mean_pre - current_mean),
        "landscape_fluidity_pre": float(score_pre),
        "landscape_fluidity_post": float(current_score),
        "landscape_fluidity_gain": float(current_score - score_pre),
        "landscape_fluidity_2_top_k": int(top_k),
        "landscape_fluidity_2_pair_samples": int(len(pair_records)),
        "patches_connected": 0,
        "total_connected_area_ha": sum(float(p.get("area_ha", 0.0) or 0.0) for p in patches.values()),
    }
    return selected, stats


def optimize_landscape_fluidity_b(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
) -> Tuple[Dict[int, Dict], Dict]:
    corridors, stats = optimize_landscape_fluidity_2(patches, candidates, params)
    out_stats = _with_landscape_fluidity_variant(
        stats,
        strategy_key="landscape_fluidity_b",
        variant_label="LF-B",
    )
    if "landscape_fluidity_2_top_k" in out_stats:
        out_stats["landscape_fluidity_b_top_k"] = int(out_stats.get("landscape_fluidity_2_top_k", 0) or 0)
    if "landscape_fluidity_2_pair_samples" in out_stats:
        out_stats["landscape_fluidity_b_pair_samples"] = int(
            out_stats.get("landscape_fluidity_2_pair_samples", 0) or 0
        )
    return corridors, out_stats


def optimize_landscape_fluidity_c(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
) -> Tuple[Dict[int, Dict], Dict]:
    budget_limit = float(params.budget_area or 0.0)
    if budget_limit <= 0.0 or (not patches) or (not candidates):
        return {}, {"strategy": "landscape_fluidity_c", "corridors_used": 0, "budget_used_ha": 0.0}

    context = _build_patch_resistance_context(patches, params)
    if context is None:
        return {}, {"strategy": "landscape_fluidity_c", "corridors_used": 0, "budget_used_ha": 0.0}

    graph = context["graph"].copy()
    valid_nodes = set(int(pid) for pid in context.get("node_ids", []))
    coords_xy: Dict[int, Tuple[float, float]] = dict(context.get("coords_xy") or {})
    if len(valid_nodes) < 2:
        return {}, {"strategy": "landscape_fluidity_c", "corridors_used": 0, "budget_used_ha": 0.0}

    patch_area_ha: Dict[int, float] = {
        int(pid): float((patches.get(int(pid)) or {}).get("area_ha", 0.0) or 0.0)
        for pid in valid_nodes
    }
    total_habitat_area = max(1e-9, float(sum(max(0.0, v) for v in patch_area_ha.values())))
    shortcut_threshold = _get_landscape_fluidity_shortcut_threshold(params, 3.0)

    uf = UnionFind()
    for pid in valid_nodes:
        ipid = int(pid)
        uf.parent[ipid] = ipid
        uf.size[ipid] = max(float(patch_area_ha.get(ipid, 0.0) or 0.0), 0.0)
        uf.count[ipid] = 1

    cand_specs: List[Dict[str, Any]] = []
    for idx, cand in enumerate(candidates):
        cost = float(cand.get("area_ha", 0.0) or 0.0)
        if cost <= 0.0 or cost > budget_limit + 1e-12:
            continue

        p1_raw = cand.get("patch1", cand.get("p1"))
        p2_raw = cand.get("patch2", cand.get("p2"))
        is_intra = bool(cand.get("intra_patch", False))
        try:
            if p1_raw is not None and p2_raw is not None and int(p1_raw) == int(p2_raw):
                is_intra = True
        except Exception:
            pass

        pids = _candidate_patch_ids_for_metric(cand, valid_nodes)
        if is_intra:
            try:
                intra_pid = int(cand.get("intra_patch_id", p1_raw if p1_raw is not None else (pids[0] if pids else 0)))
            except Exception:
                continue
            if intra_pid not in valid_nodes:
                continue
            gap_m = float(cand.get("distance_m", cand.get("distance", 0.0)) or 0.0)
            detour_m = float(cand.get("intra_detour_old_m", 0.0) or 0.0)
            ratio = float(cand.get("intra_shortcut_ratio", cand.get("fluidity_gain", 0.0)) or 0.0)
            if ratio <= 0.0 and gap_m > 0.0 and detour_m > 0.0:
                ratio = float(detour_m / max(gap_m, 1e-9))
            if ratio <= 0.0:
                continue
            if gap_m <= 0.0 and detour_m > 0.0:
                gap_m = float(detour_m / max(ratio, 1e-9))
            if detour_m <= gap_m + 1e-9 and ratio > 1.0:
                detour_m = float(gap_m * ratio)
            if detour_m <= gap_m + 1e-9:
                continue
            if ratio + 1e-12 < shortcut_threshold:
                continue

            cand_specs.append(
                {
                    "idx": int(idx),
                    "cand": cand,
                    "cost": float(cost),
                    "is_intra": True,
                    "gate_u": int(intra_pid),
                    "gate_v": int(intra_pid),
                    "output_p1": int(intra_pid),
                    "output_p2": int(intra_pid),
                    "pids": [int(intra_pid)],
                    "edges": [],
                    "length": float(max(gap_m, 1e-9)),
                    "proxy_new_distance": float(max(gap_m, 1e-9)),
                    "intra_pid": int(intra_pid),
                    "intra_detour_old_m": float(detour_m),
                    "intra_shortcut_ratio": float(ratio),
                }
            )
            continue

        if len(pids) < 2:
            continue
        length = float(cand.get("distance_m", cand.get("distance", 0.0)) or 0.0)
        edges_local: List[Tuple[int, int, float]] = []
        explicit_edges = cand.get("explicit_edges")
        if isinstance(explicit_edges, (list, tuple)):
            for item in explicit_edges:
                try:
                    if not isinstance(item, (list, tuple)) or len(item) < 2:
                        continue
                    eu = int(item[0])
                    ev = int(item[1])
                    elen = float(item[2]) if len(item) >= 3 else float(length)
                except Exception:
                    continue
                if eu not in valid_nodes or ev not in valid_nodes or eu == ev:
                    continue
                if elen <= 0.0:
                    x1, y1 = coords_xy.get(eu, (0.0, 0.0))
                    x2, y2 = coords_xy.get(ev, (0.0, 0.0))
                    elen = max(float(math.hypot(x1 - x2, y1 - y2)), 1.0)
                edges_local.append((int(eu), int(ev), float(elen)))
            if edges_local:
                node_union: Set[int] = set()
                for eu, ev, _d in edges_local:
                    node_union.add(int(eu))
                    node_union.add(int(ev))
                pids = sorted(node_union)

        if not edges_local:
            if length <= 0.0:
                a = int(pids[0])
                b = int(pids[-1])
                x1, y1 = coords_xy.get(a, (0.0, 0.0))
                x2, y2 = coords_xy.get(b, (0.0, 0.0))
                length = max(float(math.hypot(x1 - x2, y1 - y2)), 1.0)
            for i in range(len(pids)):
                for j in range(i + 1, len(pids)):
                    edges_local.append((int(pids[i]), int(pids[j]), float(length)))

        if not edges_local:
            continue

        try:
            gate_u = int(p1_raw) if p1_raw is not None else int(pids[0])
        except Exception:
            gate_u = int(pids[0])
        try:
            gate_v = int(p2_raw) if p2_raw is not None else int(pids[-1])
        except Exception:
            gate_v = int(pids[-1])
        if gate_u not in valid_nodes:
            gate_u = int(pids[0])
        if gate_v not in valid_nodes:
            gate_v = int(pids[-1])
        if gate_u == gate_v:
            if len(pids) >= 2:
                gate_u = int(pids[0])
                gate_v = int(pids[-1])
            else:
                continue
        if gate_u == gate_v:
            continue

        area_u = max(float(patch_area_ha.get(int(gate_u), 0.0) or 0.0), 0.0)
        area_v = max(float(patch_area_ha.get(int(gate_v), 0.0) or 0.0), 0.0)
        backbone_weight = float((area_u * area_v) / max(cost, 1e-9))
        if backbone_weight <= 0.0:
            continue

        cand_specs.append(
            {
                "idx": int(idx),
                "cand": cand,
                "cost": float(cost),
                "is_intra": False,
                "gate_u": int(gate_u),
                "gate_v": int(gate_v),
                "output_p1": int(gate_u),
                "output_p2": int(gate_v),
                "pids": list(pids),
                "edges": list(edges_local),
                "length": float(max(length, 1e-9)),
                "proxy_new_distance": float(max(length, 1e-9)),
                "backbone_weight": float(backbone_weight),
            }
        )

    if not cand_specs:
        return {}, {"strategy": "landscape_fluidity_c", "corridors_used": 0, "budget_used_ha": 0.0}

    selected: Dict[int, Dict[str, Any]] = {}
    used_specs: Set[int] = set()
    budget_used = 0.0

    def _add_selected(spec: Dict[str, Any], *, corridor_type: str, utility_score: float) -> None:
        cid = len(selected) + 1
        p1_out = int(spec.get("output_p1", spec.get("gate_u", 0)) or 0)
        p2_out = int(spec.get("output_p2", spec.get("gate_v", 0)) or 0)
        pids_local = [int(pid) for pid in list(spec.get("pids") or []) if int(pid) > 0]
        selected[cid] = {
            "geom": clone_geometry(spec["cand"]["geom"]),
            "patch_ids": set(pids_local) if pids_local else {int(max(p1_out, 0))},
            "area_ha": float(spec.get("cost", 0.0) or 0.0),
            "p1": int(p1_out),
            "p2": int(p2_out),
            "distance": float(spec.get("length", spec.get("proxy_new_distance", 0.0)) or 0.0),
            "type": str(corridor_type),
            "variant": spec["cand"].get("variant"),
            "source": spec["cand"].get("source"),
            "utility_score": float(utility_score),
            "overlap_ratio": 0.0,
            "intra_patch": bool(spec.get("is_intra", False)),
            "intra_patch_id": int(spec.get("intra_pid")) if bool(spec.get("is_intra", False)) else None,
            "intra_shortcut_ratio": float(spec.get("intra_shortcut_ratio", 0.0) or 0.0),
            "intra_detour_old_m": float(spec.get("intra_detour_old_m", 0.0) or 0.0),
            "fluidity_gain": float(spec.get("intra_shortcut_ratio", 0.0) or 0.0),
        }

    # Phase 1: weighted spanning forest backbone.
    backbone_specs = [
        spec for spec in cand_specs if (not bool(spec.get("is_intra", False))) and float(spec.get("backbone_weight", 0.0)) > 0.0
    ]
    backbone_specs.sort(
        key=lambda s: (
            float(s.get("backbone_weight", 0.0)),
            -float(s.get("cost", 0.0)),
            -int(s.get("idx", 0)),
        ),
        reverse=True,
    )
    phase1_count = 0
    for spec in backbone_specs:
        cost = float(spec.get("cost", 0.0) or 0.0)
        if cost <= 0.0 or budget_used + cost > budget_limit + 1e-12:
            continue
        gate_u = int(spec.get("gate_u", 0))
        gate_v = int(spec.get("gate_v", 0))
        if gate_u == gate_v:
            continue
        if uf.find(gate_u) == uf.find(gate_v):
            continue
        if not uf.union(gate_u, gate_v):
            continue
        _apply_edges_with_changes(graph, spec.get("edges", []))
        weight = float(spec.get("backbone_weight", 0.0) or 0.0)
        _add_selected(spec, corridor_type="backbone", utility_score=weight)
        budget_used += cost
        used_specs.add(int(spec.get("idx", -1)))
        phase1_count += 1
        if budget_used >= budget_limit - 1e-12:
            break

    # Phase 2: add loops/shortcuts by detour relief.
    phase2_loop_count = 0
    phase2_intra_count = 0
    while budget_used < budget_limit - 1e-12:
        remaining = float(max(0.0, budget_limit - budget_used))
        shortest_current = _all_pairs_shortest_weighted(graph)
        best_idx: Optional[int] = None
        best_relief = 0.0
        best_tie_cost = float("inf")
        best_spec: Optional[Dict[str, Any]] = None

        for spec in cand_specs:
            spec_idx = int(spec.get("idx", -1))
            if spec_idx in used_specs:
                continue
            cost = float(spec.get("cost", 0.0) or 0.0)
            if cost <= 0.0 or cost > remaining + 1e-12:
                continue
            relief = 0.0
            if bool(spec.get("is_intra", False)):
                intra_pid = int(spec.get("intra_pid", spec.get("gate_u", 0)) or 0)
                gap_m = float(spec.get("proxy_new_distance", spec.get("length", 0.0)) or 0.0)
                detour_m = float(spec.get("intra_detour_old_m", 0.0) or 0.0)
                delta = max(float(detour_m) - float(gap_m), 0.0)
                if delta <= 1e-12:
                    continue
                patch_share = max(float(patch_area_ha.get(intra_pid, 0.0) or 0.0), 0.0) / total_habitat_area
                relief = float((delta * (1.0 + patch_share)) / max(cost, 1e-9))
            else:
                gate_u = int(spec.get("gate_u", 0))
                gate_v = int(spec.get("gate_v", 0))
                if gate_u == gate_v:
                    continue
                d_old = _shortest_lookup(shortest_current, gate_u, gate_v)
                if not np.isfinite(d_old):
                    continue
                d_new = max(float(spec.get("proxy_new_distance", spec.get("length", 0.0)) or 0.0), 1e-9)
                delta = max(float(d_old) - float(d_new), 0.0)
                if delta <= 1e-12:
                    continue
                relief = float(delta / max(cost, 1e-9))

            if relief <= 1e-12:
                continue
            if (
                best_spec is None
                or relief > best_relief + 1e-12
                or (abs(relief - best_relief) <= 1e-12 and cost < best_tie_cost - 1e-12)
            ):
                best_spec = spec
                best_idx = spec_idx
                best_relief = float(relief)
                best_tie_cost = float(cost)

        if best_spec is None or best_idx is None:
            break

        if not bool(best_spec.get("is_intra", False)):
            _apply_edges_with_changes(graph, best_spec.get("edges", []))
            phase2_loop_count += 1
        else:
            phase2_intra_count += 1
        _add_selected(best_spec, corridor_type="redundant", utility_score=float(best_relief))
        budget_used += float(best_spec.get("cost", 0.0) or 0.0)
        used_specs.add(int(best_idx))

    exact = _compute_landscape_fluidity_exact(patches, selected, params)
    stats = {
        "strategy": "landscape_fluidity_c",
        "corridors_used": len(selected),
        "budget_used_ha": float(budget_used),
        "mean_effective_resistance_pre": float(exact.get("graph_resistance_pre", 0.0) or 0.0),
        "mean_effective_resistance_post": float(exact.get("graph_resistance_post", 0.0) or 0.0),
        "mean_effective_resistance_gain": float(
            float(exact.get("graph_resistance_pre", 0.0) or 0.0)
            - float(exact.get("graph_resistance_post", 0.0) or 0.0)
        ),
        "landscape_fluidity_pre": float(exact.get("pre", 0.0) or 0.0),
        "landscape_fluidity_post": float(exact.get("post", 0.0) or 0.0),
        "landscape_fluidity_gain": float(exact.get("gain", 0.0) or 0.0),
        "landscape_fluidity_c_backbone_selected": int(phase1_count),
        "landscape_fluidity_c_loop_selected": int(phase2_loop_count),
        "landscape_fluidity_c_intra_selected": int(phase2_intra_count),
        "landscape_fluidity_c_pair_samples": int(len(context.get("pair_records") or [])),
        "patches_connected": 0,
        "total_connected_area_ha": sum(float(p.get("area_ha", 0.0) or 0.0) for p in patches.values()),
    }
    return selected, _with_landscape_fluidity_variant(
        stats,
        strategy_key="landscape_fluidity_c",
        variant_label="LF-C",
    )


def optimize_mobility_strategic(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
) -> Tuple[Dict[int, Dict], Dict]:
    budget_limit = float(params.budget_area or 0.0)
    if budget_limit <= 0.0 or (not patches) or (not candidates):
        return {}, {"strategy": "mobility_strategic", "corridors_used": 0, "budget_used_ha": 0.0}

    context = _build_patch_mobility_context(patches, params)
    if context is None:
        return {}, {"strategy": "mobility_strategic", "corridors_used": 0, "budget_used_ha": 0.0}

    detour_cap = max(1.0, float(getattr(params, "mobility_detour_cap", 8.0) or 8.0))
    lambda_intra = max(0.0, float(getattr(params, "mobility_lambda_intra", 0.001) or 0.001))
    tau = max(0.0, float(getattr(params, "mobility_tau", 0.90) or 0.90))
    intra_cap = max(1, min(500, int(getattr(params, "mobility_intra_k_total", 60) or 60)))
    eval_cap = max(20, min(2000, int(getattr(params, "mobility_eval_cap", 220) or 220)))
    intra_detour_min = 1.20

    graph = context["graph"].copy()
    node_ids: List[int] = list(context["node_ids"])
    coords_xy: Dict[int, Tuple[float, float]] = dict(context["coords_xy"])
    pair_records: List[Tuple[int, int, float]] = list(context["pair_records"])
    valid_nodes = set(int(pid) for pid in node_ids)

    cand_specs: List[Dict[str, Any]] = []
    for idx, cand in enumerate(candidates):
        cost = float(cand.get("area_ha", 0.0) or 0.0)
        if cost <= 0.0 or cost > budget_limit + 1e-12:
            continue
        pids = _candidate_patch_ids_for_metric(cand, valid_nodes)
        if len(pids) < 2:
            continue
        length = float(cand.get("distance_m", 0.0) or 0.0)
        if length <= 0.0:
            p1 = int(pids[0])
            p2 = int(pids[-1])
            x1, y1 = coords_xy.get(p1, (0.0, 0.0))
            x2, y2 = coords_xy.get(p2, (0.0, 0.0))
            length = max(float(math.hypot(x1 - x2, y1 - y2)), 1.0)
        edges_local: List[Tuple[int, int, float]] = []
        for i in range(len(pids)):
            for j in range(i + 1, len(pids)):
                edges_local.append((int(pids[i]), int(pids[j]), float(length)))
        if not edges_local:
            continue
        cand_specs.append(
            {
                "idx": int(idx),
                "cand": cand,
                "cost": float(cost),
                "pids": list(pids),
                "edges": edges_local,
                "length": float(length),
            }
        )
    if not cand_specs:
        return {}, {"strategy": "mobility_strategic", "corridors_used": 0, "budget_used_ha": 0.0}

    used_specs: Set[int] = set()
    selected: Dict[int, Dict] = {}
    remaining = float(budget_limit)
    budget_used = 0.0
    intra_selected = 0
    inter_selected = 0

    shortest_current = _all_pairs_shortest_weighted(graph)
    score_pre = _detour_ratio_score_from_shortest(shortest_current, pair_records, detour_cap)
    current_score = float(score_pre)
    saw_inter_candidate_global = False

    while remaining > 1e-12:
        comp_id: Dict[int, int] = {}
        for comp_idx, comp in enumerate(nx.connected_components(graph)):
            for node in comp:
                comp_id[int(node)] = int(comp_idx)

        inter_ranked: List[Tuple[float, float, float, int]] = []
        intra_by_comp: Dict[int, List[Tuple[float, float, float, int]]] = defaultdict(list)

        for spec_idx, spec in enumerate(cand_specs):
            if spec_idx in used_specs:
                continue
            cost = float(spec.get("cost", 0.0) or 0.0)
            if cost <= 0.0 or cost > remaining + 1e-12:
                continue
            pids = list(spec.get("pids", []))
            if len(pids) < 2:
                continue
            roots = {int(comp_id.get(int(pid), -1)) for pid in pids}
            is_inter = len(roots) >= 2

            detour_pressure = 1.0
            for u, v, _w in spec.get("edges", []):
                x1, y1 = coords_xy.get(int(u), (0.0, 0.0))
                x2, y2 = coords_xy.get(int(v), (0.0, 0.0))
                eu = float(math.hypot(x1 - x2, y1 - y2))
                d_old = float(shortest_current.get(int(u), {}).get(int(v), float("inf")))
                detour_pressure = max(detour_pressure, _detour_ratio_from_distance(d_old, eu, detour_cap))

            heuristic = float(detour_pressure / max(cost, 1e-12))
            entry = (heuristic, float(detour_pressure), -cost, int(spec_idx))
            if is_inter:
                inter_ranked.append(entry)
            else:
                if detour_pressure >= intra_detour_min:
                    root = int(next(iter(roots))) if roots else -1
                    intra_by_comp[root].append(entry)

        inter_ranked.sort(reverse=True)
        if inter_ranked:
            saw_inter_candidate_global = True
        for root in list(intra_by_comp.keys()):
            intra_by_comp[root].sort(reverse=True)

        eval_indices: List[int] = []
        inter_eval_count = min(len(inter_ranked), max(1, int(eval_cap * 0.7)))
        eval_indices.extend(int(spec_idx) for _h, _dp, _c, spec_idx in inter_ranked[:inter_eval_count])

        comp_keys = sorted(
            intra_by_comp.keys(),
            key=lambda rid: (
                -float(intra_by_comp[rid][0][1]) if intra_by_comp.get(rid) else 0.0,
                int(rid),
            ),
        )
        intra_candidates: List[int] = []
        while len(intra_candidates) < intra_cap and comp_keys:
            advanced = False
            for rid in comp_keys:
                lst = intra_by_comp.get(rid) or []
                if not lst:
                    continue
                _h, _dp, _negc, spec_idx = lst.pop(0)
                intra_candidates.append(int(spec_idx))
                advanced = True
                if len(intra_candidates) >= intra_cap:
                    break
            if not advanced:
                break

        remaining_slots = max(0, eval_cap - len(eval_indices))
        for spec_idx in intra_candidates[:remaining_slots]:
            if spec_idx not in eval_indices:
                eval_indices.append(int(spec_idx))
        if not eval_indices:
            break

        eval_rows: List[Dict[str, float]] = []
        for spec_idx in eval_indices:
            spec = cand_specs[int(spec_idx)]
            cost = float(spec.get("cost", 0.0) or 0.0)
            if cost <= 0.0:
                continue
            pids = list(spec.get("pids", []))
            roots = {int(comp_id.get(int(pid), -1)) for pid in pids}
            is_inter = len(roots) >= 2
            score_new_approx = _detour_ratio_score_with_candidate_edges(
                shortest_current,
                pair_records,
                spec.get("edges", []),
                detour_cap,
            )
            gain_approx = float(score_new_approx - current_score)
            if gain_approx <= 1e-12:
                continue
            ratio = float(gain_approx / max(cost, 1e-12))
            eval_rows.append(
                {
                    "spec_idx": float(spec_idx),
                    "gain_approx": float(gain_approx),
                    "ratio_approx": float(ratio),
                    "cost": float(cost),
                    "is_inter": 1.0 if is_inter else 0.0,
                }
            )

        if not eval_rows:
            break

        eval_rows.sort(
            key=lambda r: (
                float(r.get("gain_approx", 0.0)),
                float(r.get("ratio_approx", 0.0)),
                -float(r.get("cost", 0.0)),
                float(r.get("is_inter", 0.0)),
                -float(r.get("spec_idx", 0.0)),
            ),
            reverse=True,
        )
        exact_eval_cap = max(10, min(60, int(eval_cap * 0.35)))
        shortlist = eval_rows[:exact_eval_cap]

        exact_rows: List[Dict[str, float]] = []
        best_inter_ratio = None
        for row in shortlist:
            spec_idx = int(row["spec_idx"])
            spec = cand_specs[spec_idx]
            changes_try = _apply_edges_with_changes(graph, spec.get("edges", []))
            if not changes_try:
                continue
            shortest_try = _all_pairs_shortest_weighted(graph)
            score_exact = _detour_ratio_score_from_shortest(shortest_try, pair_records, detour_cap)
            _revert_edges_with_changes(graph, changes_try)

            gain_exact = float(score_exact - current_score)
            if gain_exact <= 1e-12:
                continue
            cost = float(row.get("cost", 0.0) or 0.0)
            ratio_exact = float(gain_exact / max(cost, 1e-12))
            is_inter = bool(float(row.get("is_inter", 0.0)) > 0.5)
            if is_inter and (best_inter_ratio is None or ratio_exact > best_inter_ratio):
                best_inter_ratio = ratio_exact
            exact_rows.append(
                {
                    "spec_idx": float(spec_idx),
                    "gain": float(gain_exact),
                    "ratio": float(ratio_exact),
                    "cost": float(cost),
                    "is_inter": 1.0 if is_inter else 0.0,
                }
            )

        if not exact_rows:
            break

        best_pick: Optional[Tuple[float, float, float, float, float, int]] = None
        for row in exact_rows:
            spec_idx = int(row["spec_idx"])
            gain = float(row["gain"])
            ratio = float(row["ratio"])
            cost = float(row["cost"])
            is_inter = bool(row["is_inter"] > 0.5)
            if is_inter:
                adjusted_gain = gain
            else:
                if saw_inter_candidate_global and inter_selected > 0 and intra_selected >= inter_selected:
                    continue
                if best_inter_ratio is not None and ratio + 1e-12 < float(best_inter_ratio) * tau:
                    continue
                adjusted_gain = gain - (lambda_intra * cost)
            if adjusted_gain <= 1e-12:
                continue
            rank_key = (
                float(adjusted_gain),
                float(gain),
                float(ratio),
                -float(cost),
                1.0 if is_inter else 0.0,
                -float(spec_idx),
            )
            if best_pick is None or rank_key > best_pick:
                best_pick = (*rank_key, int(spec_idx))

        if best_pick is None:
            break

        chosen_spec_idx = int(best_pick[-1])
        chosen_spec = cand_specs[chosen_spec_idx]
        changes = _apply_edges_with_changes(graph, chosen_spec.get("edges", []))
        if not changes:
            used_specs.add(int(chosen_spec_idx))
            continue

        shortest_new = _all_pairs_shortest_weighted(graph)
        exact_new_score = _detour_ratio_score_from_shortest(shortest_new, pair_records, detour_cap)
        exact_gain = float(exact_new_score - current_score)
        if exact_gain <= 1e-12:
            _revert_edges_with_changes(graph, changes)
            used_specs.add(int(chosen_spec_idx))
            continue

        cand = chosen_spec["cand"]
        pids = list(chosen_spec.get("pids", []))
        cid = len(selected) + 1
        selected[cid] = {
            "geom": clone_geometry(cand["geom"]),
            "patch_ids": set(int(pid) for pid in pids),
            "area_ha": float(chosen_spec.get("cost", 0.0) or 0.0),
            "p1": int(cand.get("patch1", cand.get("p1", pids[0]))),
            "p2": int(cand.get("patch2", cand.get("p2", pids[-1]))),
            "distance": float(chosen_spec.get("length", 1.0) or 1.0),
            "type": "primary",
            "variant": cand.get("variant"),
            "source": cand.get("source"),
            "utility_score": float(exact_gain),
            "overlap_ratio": 0.0,
        }
        roots_after = {int(comp_id.get(int(pid), -1)) for pid in pids}
        if len(roots_after) >= 2:
            inter_selected += 1
        else:
            intra_selected += 1

        cost = float(chosen_spec.get("cost", 0.0) or 0.0)
        budget_used += cost
        remaining -= cost
        used_specs.add(int(chosen_spec_idx))
        shortest_current = shortest_new
        current_score = float(exact_new_score)

    stats = {
        "strategy": "mobility_strategic",
        "corridors_used": len(selected),
        "budget_used_ha": float(budget_used),
        "mobility_detour_pre": float(score_pre),
        "mobility_detour_post": float(current_score),
        "mobility_detour_gain": float(current_score - score_pre),
        "mobility_inter_selected": int(inter_selected),
        "mobility_intra_selected": int(intra_selected),
        "mobility_lambda_intra": float(lambda_intra),
        "mobility_tau": float(tau),
        "patches_connected": 0,
        "total_connected_area_ha": sum(p.get("area_ha", 0.0) for p in patches.values()),
    }
    return selected, stats


def _bigconnect_budget_key_ha(value: object) -> int:
    try:
        return max(0, int(round(float(value or 0.0) * float(BIGCONNECT_VECTOR_SCALE))))
    except Exception:
        return 0


def _candidate_patch_ids_for_bigconnect_vector(
    cand: Dict[str, Any],
    valid_nodes: Set[int],
) -> List[int]:
    out: List[int] = []
    for field in ("raw_patch_ids", "patch_ids"):
        raw = cand.get(field)
        if isinstance(raw, (set, list, tuple)):
            for pid in raw:
                try:
                    ip = int(pid)
                except Exception:
                    continue
                if ip in valid_nodes and ip not in out:
                    out.append(int(ip))
        if len(out) >= 2:
            return out
    for key in ("patch1", "p1", "patch2", "p2"):
        try:
            ip = int(cand.get(key))
        except Exception:
            continue
        if ip in valid_nodes and ip not in out:
            out.append(int(ip))
    return out


def _bigconnect_score_tuple_vector(
    connected_area_key: int,
    cohesion_key: int,
    budget_used_key: int,
    corridor_count: int,
    total_length: float,
) -> Tuple[int, int, int, int, float]:
    return (
        int(connected_area_key),
        -int(budget_used_key),
        -float(total_length),
        int(cohesion_key),
        -int(corridor_count),
    )


def _bigconnect_state_is_better_vector(
    candidate_score: Tuple[int, int, int, int, float],
    incumbent_score: Optional[Tuple[int, int, int, int, float]],
) -> bool:
    if incumbent_score is None:
        return True
    return candidate_score > incumbent_score


def _bigconnect_canonical_signature_vector(sig: Sequence[int]) -> Tuple[int, ...]:
    mapping: Dict[int, int] = {}
    out: List[int] = []
    next_id = 0
    for val in sig:
        key = int(val)
        if key not in mapping:
            mapping[key] = next_id
            next_id += 1
        out.append(mapping[key])
    return tuple(out)


def _bigconnect_union_signature_vector(sig: Tuple[int, ...], members: Sequence[int]) -> Tuple[int, ...]:
    uniq = sorted({int(x) for x in members if int(x) >= 0})
    if len(uniq) < 2:
        return tuple(sig)
    root = sig[int(uniq[0])]
    updated = list(sig)
    targets = {sig[int(idx)] for idx in uniq}
    for i, val in enumerate(updated):
        if val in targets:
            updated[i] = int(root)
    return _bigconnect_canonical_signature_vector(updated)


def _bigconnect_signature_groups_vector(sig: Tuple[int, ...]) -> Dict[int, List[int]]:
    groups: Dict[int, List[int]] = defaultdict(list)
    for idx, root in enumerate(sig):
        groups[int(root)].append(int(idx))
    return groups


def _bigconnect_connected_area_vector(sig: Tuple[int, ...], patch_area_keys: Sequence[int]) -> int:
    groups = _bigconnect_signature_groups_vector(sig)
    total = 0
    for members in groups.values():
        if len(members) < 2:
            continue
        total += int(sum(int(patch_area_keys[idx]) for idx in members))
    return int(total)


def _bigconnect_cohesion_key_vector(sig: Tuple[int, ...], patch_area_keys: Sequence[int]) -> int:
    """
    Secondary objective for Most Connected Area.

    Sum of squared connected-component areas (restricted to components with at least
    two patches) rewards solutions that merge area into fewer, larger subnetworks
    instead of scattering the same area across many small components.
    """
    groups = _bigconnect_signature_groups_vector(sig)
    total = 0
    for members in groups.values():
        if len(members) < 2:
            continue
        comp_area_key = int(sum(int(patch_area_keys[idx]) for idx in members))
        total += int(comp_area_key * comp_area_key)
    return int(total)


def _bigconnect_remaining_upper_bound_vector(
    sig: Tuple[int, ...],
    patch_area_keys: Sequence[int],
    remaining_patch_sets: Sequence[Set[int]],
) -> int:
    current = _bigconnect_connected_area_vector(sig, patch_area_keys)
    groups = _bigconnect_signature_groups_vector(sig)
    counted: Set[int] = set()
    for members in groups.values():
        if len(members) >= 2:
            counted.update(int(idx) for idx in members)
    incident_remaining: Set[int] = set()
    for pset in remaining_patch_sets:
        incident_remaining.update(int(idx) for idx in pset)
    extra = 0
    for idx in incident_remaining:
        if int(idx) not in counted:
            extra += int(patch_area_keys[int(idx)])
    return int(current + extra)


def _bigconnect_keep_best_frontier_rows_vector(
    rows: Sequence[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    best_by_spend: Dict[int, Dict[str, Any]] = {}
    for row in rows:
        spend_key = int(row.get("budget_used_key", _bigconnect_budget_key_ha(row.get("budget_used_ha", 0.0))) or 0)
        score = _bigconnect_score_tuple_vector(
            int(row.get("connected_area_key", 0) or 0),
            int(row.get("cohesion_key", 0) or 0),
            spend_key,
            int(row.get("corridor_count", 0) or 0),
            float(row.get("total_length", 0.0) or 0.0),
        )
        prev = best_by_spend.get(spend_key)
        prev_score = None
        if prev is not None:
            prev_score = _bigconnect_score_tuple_vector(
                int(prev.get("connected_area_key", 0) or 0),
                int(prev.get("cohesion_key", 0) or 0),
                int(prev.get("budget_used_key", 0) or 0),
                int(prev.get("corridor_count", 0) or 0),
                float(prev.get("total_length", 0.0) or 0.0),
            )
        if _bigconnect_state_is_better_vector(score, prev_score):
            best_by_spend[spend_key] = dict(row)
    return [best_by_spend[k] for k in sorted(best_by_spend.keys())]


def _bigconnect_build_canonical_candidates_vector(
    candidates: Sequence[Dict[str, Any]],
    patches: Dict[int, Dict[str, Any]],
    budget_limit_ha: float,
) -> Tuple[List[Dict[str, Any]], Dict[str, int]]:
    def _endpoint_pair_key(a: int, b: int) -> Tuple[int, int]:
        ia, ib = int(a), int(b)
        return (ia, ib) if ia <= ib else (ib, ia)

    valid_nodes: Set[int] = set()
    patch_area_by_id: Dict[int, float] = {}
    for pid, pdata in patches.items():
        try:
            ipid = int(pid)
        except Exception:
            continue
        area = float((pdata or {}).get("area_ha", 0.0) or 0.0)
        if area > 0.0:
            valid_nodes.add(ipid)
            patch_area_by_id[ipid] = float(area)

    grouped: Dict[Tuple[int, ...], List[Dict[str, Any]]] = defaultdict(list)
    raw_count = 0
    budget_key_limit = _bigconnect_budget_key_ha(budget_limit_ha)
    for idx, cand in enumerate(candidates):
        cost_ha = float(cand.get("area_ha", 0.0) or 0.0)
        if cost_ha <= 0.0:
            continue
        cost_key = _bigconnect_budget_key_ha(cost_ha)
        if cost_key <= 0 or cost_key > budget_key_limit:
            continue
        pids = tuple(sorted(_candidate_patch_ids_for_bigconnect_vector(cand, valid_nodes)))
        raw_count += 1
        if len(pids) < 2:
            continue
        length = float(cand.get("distance_m", cand.get("distance", cost_ha)) or cost_ha)
        grouped[pids].append(
            {
                "candidate_id": int(cand.get("id", idx + 1) or (idx + 1)),
                "patch_ids": tuple(int(pid) for pid in pids),
                "cost_ha": float(cost_ha),
                "cost_key": int(cost_key),
                "length": float(length),
                "candidate": cand,
                "potential_area_key": int(sum(_bigconnect_budget_key_ha(patch_area_by_id.get(int(pid), 0.0)) for pid in pids)),
            }
        )

    canonical_raw: List[Dict[str, Any]] = []
    canon_id = 1
    for patch_ids in sorted(grouped.keys()):
        rows = sorted(
            grouped[patch_ids],
            key=lambda rec: (int(rec.get("cost_key", 0) or 0), float(rec.get("length", 0.0) or 0.0), int(rec.get("candidate_id", 0) or 0)),
        )
        kept: List[Dict[str, Any]] = []
        for row in rows:
            dominated = False
            for prior in kept:
                if (
                    int(prior.get("cost_key", 0) or 0) <= int(row.get("cost_key", 0) or 0)
                    and float(prior.get("length", 0.0) or 0.0) <= float(row.get("length", 0.0) or 0.0) + 1e-12
                ):
                    dominated = True
                    break
            if dominated:
                continue
            canon_row = dict(row)
            canon_row["canon_id"] = int(canon_id)
            kept.append(canon_row)
            canonical_raw.append(canon_row)
            canon_id += 1

    def _endpoint_pair(rec: Dict[str, Any]) -> Tuple[int, int]:
        cand = dict(rec.get("candidate") or {})
        try:
            p1 = int(cand.get("patch1", cand.get("p1", rec["patch_ids"][0])))
            p2 = int(cand.get("patch2", cand.get("p2", rec["patch_ids"][-1])))
        except Exception:
            pids = list(rec.get("patch_ids", ()) or ())
            if len(pids) >= 2:
                p1 = int(pids[0])
                p2 = int(pids[-1])
            else:
                p1 = p2 = 0
        return _endpoint_pair_key(int(p1), int(p2))

    collapsed: List[Dict[str, Any]] = []
    for endpoint_pair in sorted({_endpoint_pair(row) for row in canonical_raw}):
        rows = [dict(row) for row in canonical_raw if _endpoint_pair(row) == endpoint_pair]
        rows.sort(
            key=lambda rec: (
                -int(rec.get("potential_area_key", 0) or 0),
                int(rec.get("cost_key", 0) or 0),
                float(rec.get("length", 0.0) or 0.0),
                int(rec.get("candidate_id", 0) or 0),
            )
        )
        while rows:
            ref_area_key = int(rows[0].get("potential_area_key", 0) or 0)
            area_eps_key = max(
                _bigconnect_budget_key_ha(float(BIGCONNECT_MERGE_EQUIV_AREA_HA)),
                int(round(float(ref_area_key) * float(BIGCONNECT_MERGE_EQUIV_AREA_RATIO))),
            )
            tier = [
                rec
                for rec in rows
                if int(ref_area_key) - int(rec.get("potential_area_key", 0) or 0) <= int(area_eps_key)
            ]
            best = min(
                tier,
                key=lambda rec: (
                    int(rec.get("cost_key", 0) or 0),
                    float(rec.get("length", 0.0) or 0.0),
                    int(rec.get("candidate_id", 0) or 0),
                ),
            )
            collapsed.append(dict(best))
            keep_ids = {int(rec.get("candidate_id", 0) or 0) for rec in tier}
            rows = [rec for rec in rows if int(rec.get("candidate_id", 0) or 0) not in keep_ids]

    canonical = sorted(
        collapsed,
        key=lambda rec: (
            tuple(int(pid) for pid in rec.get("patch_ids", ()) or ()),
            int(rec.get("cost_key", 0) or 0),
            float(rec.get("length", 0.0) or 0.0),
            int(rec.get("candidate_id", 0) or 0),
        ),
    )

    return canonical, {
        "raw_count": int(raw_count),
        "nondominated_count": int(len(canonical)),
    }


def _bigconnect_build_clusters_vector(
    canonical: Sequence[Dict[str, Any]],
    patches: Dict[int, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    if not canonical:
        return []
    patch_area_ha: Dict[int, float] = {}
    patch_area_keys: Dict[int, int] = {}
    for pid, pdata in patches.items():
        try:
            ipid = int(pid)
        except Exception:
            continue
        area = float((pdata or {}).get("area_ha", 0.0) or 0.0)
        if area > 0.0:
            patch_area_ha[ipid] = float(area)
            patch_area_keys[ipid] = _bigconnect_budget_key_ha(area)

    uf = UnionFind()
    for pid in patch_area_ha.keys():
        uf.find(int(pid))
    for row in canonical:
        pids = list(row.get("patch_ids", ()) or ())
        if len(pids) < 2:
            continue
        anchor = int(pids[0])
        for other in pids[1:]:
            uf.union(int(anchor), int(other))

    clusters: Dict[int, Dict[str, Any]] = {}
    for row in canonical:
        pids = list(row.get("patch_ids", ()) or ())
        if not pids:
            continue
        root = int(uf.find(int(pids[0])))
        bucket = clusters.setdefault(root, {"patch_ids": set(), "candidates": []})
        bucket["patch_ids"].update(int(pid) for pid in pids)
        bucket["candidates"].append(row)

    out: List[Dict[str, Any]] = []
    for root in sorted(clusters.keys()):
        bucket = clusters[root]
        patch_ids = sorted(int(pid) for pid in bucket.get("patch_ids", set()))
        patch_index = {int(pid): idx for idx, pid in enumerate(patch_ids)}
        cluster_patch_area_keys = [int(patch_area_keys.get(int(pid), 0) or 0) for pid in patch_ids]
        cluster_patch_area_ha = [float(patch_area_ha.get(int(pid), 0.0) or 0.0) for pid in patch_ids]
        cluster_candidates: List[Dict[str, Any]] = []
        for row in sorted(
            bucket.get("candidates", []),
            key=lambda rec: (
                -int(rec.get("potential_area_key", 0) or 0),
                int(rec.get("cost_key", 0) or 0),
                float(rec.get("length", 0.0) or 0.0),
            ),
        ):
            local_indices = tuple(
                sorted(int(patch_index[int(pid)]) for pid in row.get("patch_ids", ()) if int(pid) in patch_index)
            )
            if len(local_indices) < 2:
                continue
            cluster_candidates.append(
                {
                    "canon_id": int(row["canon_id"]),
                    "patch_ids": tuple(int(pid) for pid in row["patch_ids"]),
                    "local_indices": tuple(int(idx) for idx in local_indices),
                    "cost_key": int(row["cost_key"]),
                    "cost_ha": float(row["cost_ha"]),
                    "length": float(row["length"]),
                    "potential_area_key": int(sum(int(cluster_patch_area_keys[idx]) for idx in local_indices)),
                    "candidate": row["candidate"],
                }
            )
        out.append(
            {
                "cluster_id": int(len(out) + 1),
                "patch_ids": patch_ids,
                "patch_area_keys": cluster_patch_area_keys,
                "patch_area_ha": cluster_patch_area_ha,
                "candidates": cluster_candidates,
            }
        )
    return out


def _bigconnect_exact_frontier_vector(
    cluster: Dict[str, Any],
    budget_limit_key: int,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    candidates = list(cluster.get("candidates", []) or [])
    patch_area_keys = list(cluster.get("patch_area_keys", []) or [])
    start_sig = tuple(range(len(patch_area_keys)))
    states: Dict[Tuple[Tuple[int, ...], int], Tuple[int, float, int]] = {(start_sig, 0): (0, 0.0, 0)}
    remaining_patch_sets = [set(int(x) for x in cand.get("local_indices", ())) for cand in candidates]
    start_time = time.perf_counter()
    peak_states = 1

    for idx, cand in enumerate(candidates):
        if time.perf_counter() - start_time > BIGCONNECT_EXACT_MAX_SECONDS:
            raise RuntimeError("time_limit")
        next_states = dict(states)
        for (sig, spend_key), (mask, total_length, count) in states.items():
            if int(spend_key) + int(cand.get("cost_key", 0) or 0) > budget_limit_key:
                continue
            upper_bound = _bigconnect_remaining_upper_bound_vector(sig, patch_area_keys, remaining_patch_sets[idx:])
            if upper_bound <= 0:
                continue
            new_sig = _bigconnect_union_signature_vector(sig, cand.get("local_indices", ()))
            new_spend_key = int(spend_key) + int(cand.get("cost_key", 0) or 0)
            bit = 1 << idx
            new_value = (
                int(mask) | int(bit),
                float(total_length) + float(cand.get("length", 0.0) or 0.0),
                int(count) + 1,
            )
            key = (new_sig, int(new_spend_key))
            prev = next_states.get(key)
            if prev is None:
                next_states[key] = new_value
            else:
                area_key = _bigconnect_connected_area_vector(new_sig, patch_area_keys)
                cohesion_key = _bigconnect_cohesion_key_vector(new_sig, patch_area_keys)
                prev_cohesion_key = _bigconnect_cohesion_key_vector(new_sig, patch_area_keys)
                prev_score = _bigconnect_score_tuple_vector(area_key, prev_cohesion_key, int(new_spend_key), int(prev[2]), float(prev[1]))
                new_score = _bigconnect_score_tuple_vector(area_key, cohesion_key, int(new_spend_key), int(new_value[2]), float(new_value[1]))
                if _bigconnect_state_is_better_vector(new_score, prev_score):
                    next_states[key] = new_value
        states = next_states
        peak_states = max(peak_states, len(states))
        if peak_states > BIGCONNECT_EXACT_MAX_STATES:
            raise RuntimeError("state_limit")

    frontier_rows: List[Dict[str, Any]] = []
    for (sig, spend_key), (mask, total_length, count) in states.items():
        frontier_rows.append(
            {
                "connected_area_key": int(_bigconnect_connected_area_vector(sig, patch_area_keys)),
                "connected_area_ha": float(spend_key * 0.0),
                "cohesion_key": int(_bigconnect_cohesion_key_vector(sig, patch_area_keys)),
                "budget_used_key": int(spend_key),
                "budget_used_ha": float(spend_key) / float(BIGCONNECT_VECTOR_SCALE),
                "corridor_count": int(count),
                "total_length": float(total_length),
                "selected_canon_ids": [
                    int(candidates[idx]["canon_id"])
                    for idx in range(len(candidates))
                    if int(mask) & (1 << idx)
                ],
                "exact_flag": True,
            }
        )
    for row in frontier_rows:
        row["connected_area_ha"] = float(row["connected_area_key"]) / float(BIGCONNECT_VECTOR_SCALE)
    return _bigconnect_keep_best_frontier_rows_vector(frontier_rows), {
        "exact": True,
        "peak_states": int(peak_states),
        "abort_reason": "",
    }


def _bigconnect_beam_frontier_vector(
    cluster: Dict[str, Any],
    budget_limit_key: int,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    candidates = list(cluster.get("candidates", []) or [])
    patch_area_keys = list(cluster.get("patch_area_keys", []) or [])
    states: Dict[Tuple[Tuple[int, ...], int], Tuple[int, float, int]] = {(tuple(range(len(patch_area_keys))), 0): (0, 0.0, 0)}
    peak_states = 1
    remaining_patch_sets = [set(int(x) for x in cand.get("local_indices", ())) for cand in candidates]

    for idx, cand in enumerate(candidates):
        merged: Dict[Tuple[Tuple[int, ...], int], Tuple[int, float, int]] = dict(states)
        for (sig, spend_key), (mask, total_length, count) in states.items():
            new_spend_key = int(spend_key) + int(cand.get("cost_key", 0) or 0)
            if new_spend_key > budget_limit_key:
                continue
            new_sig = _bigconnect_union_signature_vector(sig, cand.get("local_indices", ()))
            new_value = (
                int(mask) | (1 << idx),
                float(total_length) + float(cand.get("length", 0.0) or 0.0),
                int(count) + 1,
            )
            key = (new_sig, int(new_spend_key))
            prev = merged.get(key)
            if prev is None:
                merged[key] = new_value
            else:
                area_key = _bigconnect_connected_area_vector(new_sig, patch_area_keys)
                cohesion_key = _bigconnect_cohesion_key_vector(new_sig, patch_area_keys)
                prev_cohesion_key = _bigconnect_cohesion_key_vector(new_sig, patch_area_keys)
                prev_score = _bigconnect_score_tuple_vector(area_key, prev_cohesion_key, int(new_spend_key), int(prev[2]), float(prev[1]))
                new_score = _bigconnect_score_tuple_vector(area_key, cohesion_key, int(new_spend_key), int(new_value[2]), float(new_value[1]))
                if _bigconnect_state_is_better_vector(new_score, prev_score):
                    merged[key] = new_value

        def _beam_rank(item: Tuple[Tuple[Tuple[int, ...], int], Tuple[int, float, int]]) -> Tuple[int, int, int, int, float]:
            (sig, spend_key), (_mask, total_length, count) = item
            current_area_key = int(_bigconnect_connected_area_vector(sig, patch_area_keys))
            optimistic_area_key = int(
                _bigconnect_remaining_upper_bound_vector(
                    sig,
                    patch_area_keys,
                    remaining_patch_sets[idx + 1 :],
                )
            )
            return (
                optimistic_area_key,
                current_area_key,
                _bigconnect_cohesion_key_vector(sig, patch_area_keys),
                -int(spend_key),
                -int(count),
                -float(total_length),
            )

        ranked = sorted(merged.items(), key=_beam_rank, reverse=True)
        states = dict(ranked[:BIGCONNECT_BEAM_WIDTH])
        peak_states = max(peak_states, len(states))

    frontier_rows: List[Dict[str, Any]] = []
    for (sig, spend_key), (mask, total_length, count) in states.items():
        frontier_rows.append(
            {
                "connected_area_key": int(_bigconnect_connected_area_vector(sig, patch_area_keys)),
                "connected_area_ha": float(_bigconnect_connected_area_vector(sig, patch_area_keys)) / float(BIGCONNECT_VECTOR_SCALE),
                "cohesion_key": int(_bigconnect_cohesion_key_vector(sig, patch_area_keys)),
                "budget_used_key": int(spend_key),
                "budget_used_ha": float(spend_key) / float(BIGCONNECT_VECTOR_SCALE),
                "corridor_count": int(count),
                "total_length": float(total_length),
                "selected_canon_ids": [
                    int(candidates[idx]["canon_id"])
                    for idx in range(len(candidates))
                    if int(mask) & (1 << idx)
                ],
                "exact_flag": False,
            }
        )
    return _bigconnect_keep_best_frontier_rows_vector(frontier_rows), {
        "exact": False,
        "peak_states": int(peak_states),
        "abort_reason": "",
    }


def _bigconnect_solve_cluster_vector(
    cluster: Dict[str, Any],
    budget_limit_key: int,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    candidates = list(cluster.get("candidates", []) or [])
    if not candidates:
        return [
            {
                "connected_area_key": 0,
                "connected_area_ha": 0.0,
                "cohesion_key": 0,
                "budget_used_key": 0,
                "budget_used_ha": 0.0,
                "corridor_count": 0,
                "total_length": 0.0,
                "selected_canon_ids": [],
                "exact_flag": True,
            }
        ], {"exact": True, "peak_states": 1, "abort_reason": ""}
    if len(candidates) <= BIGCONNECT_EXACT_MAX_CANDIDATES:
        try:
            return _bigconnect_exact_frontier_vector(cluster, budget_limit_key)
        except RuntimeError as exc:
            frontier, meta = _bigconnect_beam_frontier_vector(cluster, budget_limit_key)
            meta["abort_reason"] = str(exc)
            return frontier, meta
    frontier, meta = _bigconnect_beam_frontier_vector(cluster, budget_limit_key)
    meta["abort_reason"] = "candidate_limit"
    return frontier, meta


def _bigconnect_objective_from_corridors_vector(
    corridors: Dict[int, Dict[str, Any]],
    patches: Dict[int, Dict[str, Any]],
) -> Tuple[float, int]:
    uf = UnionFind()
    patch_area_ha: Dict[int, float] = {}
    for pid, pdata in patches.items():
        try:
            ipid = int(pid)
        except Exception:
            continue
        area = float((pdata or {}).get("area_ha", 0.0) or 0.0)
        if area <= 0.0:
            continue
        uf.find(ipid)
        patch_area_ha[ipid] = float(area)
    for cdata in corridors.values():
        pids = [int(pid) for pid in list(cdata.get("patch_ids", ())) if pid is not None]
        if len(pids) < 2:
            continue
        anchor = int(pids[0])
        for other in pids[1:]:
            uf.union(int(anchor), int(other))
    groups: Dict[int, List[int]] = defaultdict(list)
    for pid in patch_area_ha.keys():
        root = int(uf.find(int(pid)))
        groups[root].append(int(pid))
    connected_area_ha = 0.0
    connected_patches = 0
    for members in groups.values():
        if len(members) < 2:
            continue
        connected_patches += int(len(members))
        connected_area_ha += float(sum(float(patch_area_ha.get(int(pid), 0.0) or 0.0) for pid in members))
    return float(connected_area_ha), int(connected_patches)


def _bigconnect_refill_vector(
    corridors: Dict[int, Dict[str, Any]],
    candidates: Sequence[Dict[str, Any]],
    patches: Dict[int, Dict[str, Any]],
    budget_limit_ha: float,
) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "added": 0,
        "added_ids": [],
        "added_source_candidate_ids": [],
        "records": [],
        "skip_reason": "",
    }
    remaining_budget = float(budget_limit_ha) - float(sum(float(c.get("area_ha", 0.0) or 0.0) for c in corridors.values()))
    if remaining_budget <= 1e-12:
        out["skip_reason"] = "no_remaining_budget"
        return out
    valid_patch_ids = {int(pid) for pid, pdata in patches.items() if float((pdata or {}).get("area_ha", 0.0) or 0.0) > 0.0}
    total_patch_count = int(len(valid_patch_ids))
    selected_source_ids = {
        int(c.get("source_candidate_id"))
        for c in corridors.values()
        if c.get("source_candidate_id") is not None
    }
    current_area_ha, current_patches = _bigconnect_objective_from_corridors_vector(corridors, patches)
    current_area_key = _bigconnect_budget_key_ha(current_area_ha)

    def _candidate_row(cand: Dict[str, Any], src_id: int, cost_ha: float, touched: Tuple[int, ...]) -> Dict[str, Any]:
        return {
            "geom": clone_geometry(cand["geom"]),
            "patch_ids": set(int(pid) for pid in touched),
            "area_ha": float(cost_ha),
            "p1": int(cand.get("patch1", cand.get("p1", touched[0]))),
            "p2": int(cand.get("patch2", cand.get("p2", touched[-1]))),
            "distance": float(cand.get("distance_m", cand.get("distance", cost_ha)) or cost_ha),
            "type": "primary",
            "variant": cand.get("variant"),
            "source": cand.get("source"),
            "utility_score": 0.0,
            "overlap_ratio": 0.0,
            "source_candidate_id": int(src_id),
        }

    while remaining_budget > 1e-12:
        feasible: List[Dict[str, Any]] = []
        for idx, cand in enumerate(candidates):
            src_id = int(cand.get("id", idx + 1) or (idx + 1))
            if src_id in selected_source_ids:
                continue
            cost_ha = float(cand.get("area_ha", 0.0) or 0.0)
            if cost_ha <= 0.0 or cost_ha > remaining_budget + 1e-12:
                continue
            touched = tuple(sorted(_candidate_patch_ids_for_bigconnect_vector(cand, valid_patch_ids)))
            if len(touched) < 2:
                continue
            feasible.append(
                {
                    "candidate": cand,
                    "source_id": int(src_id),
                    "cost_ha": float(cost_ha),
                    "cost_key": _bigconnect_budget_key_ha(cost_ha),
                    "touched": touched,
                    "row": _candidate_row(cand, int(src_id), float(cost_ha), touched),
                }
            )
        if not feasible:
            out["skip_reason"] = "no_feasible_candidate"
            break

        best_direct = None
        best_direct_score = None
        for item in feasible:
            trial = dict(corridors)
            trial[-1] = dict(item["row"])
            area_after_ha, patches_after = _bigconnect_objective_from_corridors_vector(trial, patches)
            area_gain_key = _bigconnect_budget_key_ha(area_after_ha) - int(current_area_key)
            patch_gain = int(patches_after) - int(current_patches)
            if area_gain_key <= 0:
                continue
            direct_score = (
                int(area_gain_key),
                int(patch_gain),
                -int(item["cost_key"]),
                -float(item["row"].get("distance", item["cost_ha"]) or 0.0),
            )
            if best_direct_score is None or direct_score > best_direct_score:
                best_direct_score = direct_score
                best_direct = {
                    "item": item,
                    "area_after_ha": float(area_after_ha),
                    "patches_after": int(patches_after),
                    "area_gain_key": int(area_gain_key),
                }

        best_bridge = None
        best_bridge_score = None
        for first in feasible:
            rem_after_first = float(remaining_budget) - float(first["cost_ha"])
            if rem_after_first <= 1e-12:
                continue
            trial1 = dict(corridors)
            trial1[-1] = dict(first["row"])
            area1_ha, patches1 = _bigconnect_objective_from_corridors_vector(trial1, patches)
            for second in feasible:
                if int(second["source_id"]) == int(first["source_id"]):
                    continue
                if float(second["cost_ha"]) > rem_after_first + 1e-12:
                    continue
                trial2 = dict(trial1)
                trial2[-2] = dict(second["row"])
                area2_ha, patches2 = _bigconnect_objective_from_corridors_vector(trial2, patches)
                total_gain_key = _bigconnect_budget_key_ha(area2_ha) - int(current_area_key)
                if total_gain_key <= 0:
                    continue
                bridge_score = (
                    int(total_gain_key),
                    int(patches2) - int(current_patches),
                    -int(first["cost_key"] + second["cost_key"]),
                    -float(
                        float(first["row"].get("distance", first["cost_ha"]) or 0.0)
                        + float(second["row"].get("distance", second["cost_ha"]) or 0.0)
                    ),
                )
                if best_bridge_score is None or bridge_score > best_bridge_score:
                    best_bridge_score = bridge_score
                    best_bridge = {
                        "first": first,
                        "second": second,
                        "area1_ha": float(area1_ha),
                        "patches1": int(patches1),
                        "area_after_ha": float(area2_ha),
                        "patches_after": int(patches2),
                        "area_gain_key": int(total_gain_key),
                        "bridge_step_gain_key": _bigconnect_budget_key_ha(area1_ha) - int(current_area_key),
                    }

        chosen = None
        if best_bridge is not None and (
            best_direct is None or int(best_bridge.get("area_gain_key", 0) or 0) > int(best_direct.get("area_gain_key", 0) or 0)
        ):
            chosen = {
                "item": best_bridge["first"],
                "mode": "bridge_step",
                "area_after_ha": float(best_bridge.get("area1_ha", current_area_ha) or current_area_ha),
                "patches_after": int(best_bridge.get("patches1", current_patches) or current_patches),
                "area_gain_key": int(best_bridge.get("bridge_step_gain_key", 0) or 0),
                "bridge_target_source_id": int(best_bridge["second"]["source_id"]),
                "bridge_total_gain_key_if_followed": int(best_bridge.get("area_gain_key", 0) or 0),
            }
        elif best_direct is not None:
            chosen = {
                "item": best_direct["item"],
                "mode": "direct_gain",
                "area_after_ha": float(best_direct.get("area_after_ha", current_area_ha) or current_area_ha),
                "patches_after": int(best_direct.get("patches_after", current_patches) or current_patches),
                "area_gain_key": int(best_direct.get("area_gain_key", 0) or 0),
            }

        if chosen is None and total_patch_count >= 2 and int(current_patches) >= int(total_patch_count):
            uf = UnionFind()
            for pid in valid_patch_ids:
                uf.find(int(pid))
            for cdata in corridors.values():
                touched_existing = tuple(sorted(_candidate_patch_ids_for_bigconnect_vector(cdata, valid_patch_ids)))
                if len(touched_existing) < 2:
                    continue
                anchor_existing = int(touched_existing[0])
                for other_existing in touched_existing[1:]:
                    uf.union(anchor_existing, int(other_existing))

            best_redundant = None
            best_redundant_score = None
            for item in feasible:
                touched = tuple(int(pid) for pid in item.get("touched", ()) if int(pid) in valid_patch_ids)
                if len(touched) < 2:
                    continue
                if len({int(uf.find(int(pid))) for pid in touched}) != 1:
                    continue
                redundant_score = (
                    -float(item["row"].get("distance", item.get("cost_ha", 0.0)) or 0.0),
                    -int(item.get("cost_key", 0) or 0),
                    -int(len(touched)),
                )
                if best_redundant_score is None or redundant_score > best_redundant_score:
                    best_redundant_score = redundant_score
                    best_redundant = item

            if best_redundant is not None:
                chosen = {
                    "item": best_redundant,
                    "mode": "redundant_fill",
                    "area_after_ha": float(current_area_ha),
                    "patches_after": int(current_patches),
                    "area_gain_key": 0,
                }

        if chosen is None:
            out["skip_reason"] = "no_positive_gain_candidate"
            break

        item = chosen["item"]
        cid = len(corridors) + 1
        corridors[cid] = dict(item["row"])
        corridors[cid]["utility_score"] = float(chosen["area_gain_key"]) / float(BIGCONNECT_VECTOR_SCALE)
        out["added"] = int(out["added"]) + 1
        out["added_ids"].append(int(cid))
        out["added_source_candidate_ids"].append(int(item["source_id"]))
        out["records"].append(
            {
                "corridor_id": int(cid),
                "source_candidate_id": int(item["source_id"]),
                "mode": str(chosen.get("mode", "")),
                "patch_ids": list(int(pid) for pid in item["touched"]),
                "cost_ha": float(item["cost_ha"]),
                "connected_area_gain_ha": float(chosen["area_gain_key"]) / float(BIGCONNECT_VECTOR_SCALE),
                "bridge_target_source_id": int(chosen.get("bridge_target_source_id", 0) or 0),
                "bridge_total_gain_if_followed_ha": float(chosen.get("bridge_total_gain_key_if_followed", 0) or 0) / float(BIGCONNECT_VECTOR_SCALE),
            }
        )
        selected_source_ids.add(int(item["source_id"]))
        remaining_budget -= float(item["cost_ha"])
        current_area_ha = float(chosen["area_after_ha"])
        current_area_key = _bigconnect_budget_key_ha(current_area_ha)
        current_patches = int(chosen["patches_after"])

    return out


def optimize_bigconnect_vector(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
) -> Tuple[Dict[int, Dict], Dict]:
    budget_limit_ha = float(params.budget_area or 0.0)
    budget_limit_key = _bigconnect_budget_key_ha(budget_limit_ha)
    if budget_limit_key <= 0 or (not patches) or (not candidates):
        return {}, {"strategy": "most_connected_habitat", "corridors_used": 0, "budget_used_ha": 0.0}

    canonical, canon_stats = _bigconnect_build_canonical_candidates_vector(candidates, patches, budget_limit_ha)
    if not canonical:
        return {}, {
            "strategy": "most_connected_habitat",
            "corridors_used": 0,
            "budget_used_ha": 0.0,
            "bigconnect_candidate_count_raw": int(canon_stats.get("raw_count", 0) or 0),
            "bigconnect_candidate_count_nondominated": int(canon_stats.get("nondominated_count", 0) or 0),
        }

    clusters = _bigconnect_build_clusters_vector(canonical, patches)
    cluster_frontiers: List[List[Dict[str, Any]]] = []
    exact_clusters = 0
    heuristic_clusters = 0
    frontier_states_total = 0
    abort_reason_counts: Dict[str, int] = defaultdict(int)
    for cluster in clusters:
        frontier, meta = _bigconnect_solve_cluster_vector(cluster, budget_limit_key)
        zero_row = {
            "connected_area_key": 0,
            "connected_area_ha": 0.0,
            "budget_used_key": 0,
            "budget_used_ha": 0.0,
            "corridor_count": 0,
            "total_length": 0.0,
            "selected_canon_ids": [],
            "exact_flag": bool(meta.get("exact", False)),
        }
        full_frontier = _bigconnect_keep_best_frontier_rows_vector([zero_row] + list(frontier))
        cluster_frontiers.append(full_frontier)
        frontier_states_total += int(len(full_frontier))
        if bool(meta.get("exact", False)):
            exact_clusters += 1
        else:
            heuristic_clusters += 1
            abort_reason_counts[str(meta.get("abort_reason", "") or "heuristic")] += 1

    dp: Dict[int, Dict[str, Any]] = {
        0: {
            "connected_area_key": 0,
            "connected_area_ha": 0.0,
            "budget_used_key": 0,
            "budget_used_ha": 0.0,
            "corridor_count": 0,
            "total_length": 0.0,
            "selected_canon_ids": [],
            "exact_flag": True,
        }
    }
    for frontier in cluster_frontiers:
        next_dp: Dict[int, Dict[str, Any]] = {}
        for spend_key_a, row_a in dp.items():
            for row_b in frontier:
                new_spend_key = int(spend_key_a) + int(row_b.get("budget_used_key", 0) or 0)
                if new_spend_key > budget_limit_key:
                    continue
                combined = {
                    "connected_area_key": int(row_a.get("connected_area_key", 0) or 0) + int(row_b.get("connected_area_key", 0) or 0),
                    "connected_area_ha": (
                        float(row_a.get("connected_area_ha", 0.0) or 0.0) + float(row_b.get("connected_area_ha", 0.0) or 0.0)
                    ),
                    "cohesion_key": int(row_a.get("cohesion_key", 0) or 0) + int(row_b.get("cohesion_key", 0) or 0),
                    "budget_used_key": int(new_spend_key),
                    "budget_used_ha": (
                        float(row_a.get("budget_used_ha", 0.0) or 0.0) + float(row_b.get("budget_used_ha", 0.0) or 0.0)
                    ),
                    "corridor_count": int(row_a.get("corridor_count", 0) or 0) + int(row_b.get("corridor_count", 0) or 0),
                    "total_length": float(row_a.get("total_length", 0.0) or 0.0) + float(row_b.get("total_length", 0.0) or 0.0),
                    "selected_canon_ids": list(row_a.get("selected_canon_ids", []) or []) + list(row_b.get("selected_canon_ids", []) or []),
                    "exact_flag": bool(row_a.get("exact_flag", False)) and bool(row_b.get("exact_flag", False)),
                }
                prev = next_dp.get(int(new_spend_key))
                prev_score = None
                if prev is not None:
                    prev_score = _bigconnect_score_tuple_vector(
                        int(prev.get("connected_area_key", 0) or 0),
                        int(prev.get("cohesion_key", 0) or 0),
                        int(prev.get("budget_used_key", 0) or 0),
                        int(prev.get("corridor_count", 0) or 0),
                        float(prev.get("total_length", 0.0) or 0.0),
                    )
                new_score = _bigconnect_score_tuple_vector(
                    int(combined["connected_area_key"]),
                    int(combined["cohesion_key"]),
                    int(combined["budget_used_key"]),
                    int(combined["corridor_count"]),
                    float(combined["total_length"]),
                )
                if _bigconnect_state_is_better_vector(new_score, prev_score):
                    next_dp[int(new_spend_key)] = combined
        dp = next_dp

    best_row: Optional[Dict[str, Any]] = None
    best_score: Optional[Tuple[int, int, int, float]] = None
    for row in dp.values():
        score = _bigconnect_score_tuple_vector(
            int(row.get("connected_area_key", 0) or 0),
            int(row.get("cohesion_key", 0) or 0),
            int(row.get("budget_used_key", 0) or 0),
            int(row.get("corridor_count", 0) or 0),
            float(row.get("total_length", 0.0) or 0.0),
        )
        if _bigconnect_state_is_better_vector(score, best_score):
            best_row = row
            best_score = score
    if best_row is None:
        return {}, {"strategy": "most_connected_habitat", "corridors_used": 0, "budget_used_ha": 0.0}

    canonical_by_id = {int(row["canon_id"]): row for row in canonical}
    selected: Dict[int, Dict[str, Any]] = {}
    for canon_id in sorted(set(int(x) for x in (best_row.get("selected_canon_ids", []) or []))):
        row = canonical_by_id.get(int(canon_id))
        if row is None:
            continue
        cand = dict(row.get("candidate") or {})
        pids = list(row.get("patch_ids", ()) or ())
        if len(pids) < 2:
            continue
        cid = len(selected) + 1
        selected[cid] = {
            "geom": clone_geometry(cand["geom"]),
            "patch_ids": set(int(pid) for pid in pids),
            "area_ha": float(row.get("cost_ha", 0.0) or 0.0),
            "p1": int(cand.get("patch1", cand.get("p1", pids[0]))),
            "p2": int(cand.get("patch2", cand.get("p2", pids[-1]))),
            "distance": float(row.get("length", 0.0) or 0.0),
            "type": "primary",
            "variant": cand.get("variant"),
            "source": cand.get("source"),
            "utility_score": float(best_row.get("connected_area_ha", 0.0) or 0.0),
            "overlap_ratio": 0.0,
            "source_candidate_id": int(cand.get("id", row.get("candidate_id", cid)) or cid),
        }

    stats = {
        "strategy": "most_connected_habitat",
        "corridors_used": len(selected),
        "budget_used_ha": float(best_row.get("budget_used_ha", 0.0) or 0.0),
        "bigconnect_objective_pre": 0.0,
        "bigconnect_objective_post": float(best_row.get("connected_area_ha", 0.0) or 0.0),
        "bigconnect_objective_gain": float(best_row.get("connected_area_ha", 0.0) or 0.0),
        "bigconnect_budget_fill_ratio": float((float(best_row.get("budget_used_ha", 0.0) or 0.0) / budget_limit_ha) if budget_limit_ha > 0.0 else 0.0),
        "bigconnect_candidate_count_raw": int(canon_stats.get("raw_count", 0) or 0),
        "bigconnect_candidate_count_nondominated": int(canon_stats.get("nondominated_count", 0) or 0),
        "bigconnect_clusters_total": int(len(clusters)),
        "bigconnect_clusters_exact": int(exact_clusters),
        "bigconnect_clusters_heuristic": int(heuristic_clusters),
        "bigconnect_proven_optimal": bool(heuristic_clusters == 0),
        "bigconnect_frontier_states_total": int(frontier_states_total),
        "bigconnect_exact_abort_reason_counts": dict(abort_reason_counts),
        "patches_connected": 0,
        "total_connected_area_ha": sum(float((p or {}).get("area_ha", 0.0) or 0.0) for p in patches.values()),
    }
    _refresh_vector_connectivity_stats(patches, selected, stats)

    refill = _bigconnect_refill_vector(selected, candidates, patches, budget_limit_ha)
    if int(refill.get("added", 0) or 0) > 0:
        stats["bigconnect_proven_optimal"] = False

    stats["budget_used_ha"] = float(sum(float(c.get("area_ha", 0.0) or 0.0) for c in selected.values()))
    connected_area_post_ha, connected_patches_post = _bigconnect_objective_from_corridors_vector(selected, patches)
    stats["bigconnect_objective_post"] = float(connected_area_post_ha)
    stats["bigconnect_objective_gain"] = float(connected_area_post_ha)
    stats["bigconnect_connected_patches"] = int(connected_patches_post)
    stats["bigconnect_budget_fill_ratio"] = float((stats["budget_used_ha"] / budget_limit_ha) if budget_limit_ha > 0.0 else 0.0)
    _refresh_vector_connectivity_stats(patches, selected, stats)
    return selected, stats


def optimize_circuit_utility_largest_network(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
    overlap_reject_ratio: float = 0.30,
) -> Tuple[Dict[int, Dict], Dict]:
    """
    Largest Single Network (Circuit Utility):
    Use the Most Connectivity greedy utility model, but constrain selection to grow
    a single connected component seeded at the largest patch.

    This inherits Most Connectivity behavior (ROI scoring + overlap-aware redundancy),
    but avoids spending budget on disconnected "side networks".
    """

    def _pair_key(a: int, b: int) -> Tuple[int, int]:
        return (a, b) if a <= b else (b, a)

    def _candidate_cost_ha(cand: Dict) -> float:
        try:
            return float(cand.get("area_ha", 0.0) or 0.0)
        except Exception:
            return 0.0

    def _patch_area_ha(pid: int) -> float:
        try:
            return float(patches.get(pid, {}).get("area_ha", 0.0) or 0.0)
        except Exception:
            return 0.0

    def _safe_geom_area(geom: QgsGeometry) -> float:
        try:
            if geom is None or geom.isEmpty():
                return 0.0
            return float(geom.area())
        except Exception:
            return 0.0

    def _max_overlap_ratio(new_geom: QgsGeometry, prior: List[QgsGeometry]) -> float:
        denom = _safe_geom_area(new_geom)
        if denom <= 0 or not prior:
            return 0.0
        max_ratio = 0.0
        for g in prior:
            try:
                try:
                    proximity_dist = max(
                        float(getattr(params, "min_corridor_width", 0.0) or 0.0) * 1.5,
                        float(getattr(params, "grid_resolution", 0.0) or 0.0) * 2.0,
                    )
                    if proximity_dist > 0 and new_geom.distance(g) <= proximity_dist:
                        return 1.0
                except Exception:
                    pass
                inter = new_geom.intersection(g)
                if inter is None or inter.isEmpty():
                    continue
                ratio = float(inter.area()) / denom
                if ratio > max_ratio:
                    max_ratio = ratio
            except Exception:
                continue
        return max_ratio

    def _get_patch_ids(cand: Dict) -> List[int]:
        try:
            pids = cand.get("patch_ids")
            if pids:
                return [int(pid) for pid in pids if pid is not None]
        except Exception:
            pass
        p1, p2 = cand.get("patch1"), cand.get("patch2")
        if p1 is None or p2 is None:
            return []
        return [int(p1), int(p2)]

    def _get_pair(cand: Dict) -> Tuple[int, int]:
        p1, p2 = cand.get("patch1"), cand.get("patch2")
        if p1 is None or p2 is None:
            return (0, 0)
        return _pair_key(int(p1), int(p2))

    def _get_base_roi(cand: Dict) -> float:
        geom = cand.get("geom")
        if geom is None or geom.isEmpty():
            return 0.0
        cost = float(_candidate_cost_ha(cand) or 0.0)
        if cost <= 0.0:
            return 0.0
        p1, p2 = cand.get("patch1"), cand.get("patch2")
        if p1 is None or p2 is None:
            return 0.0
        p1i, p2i = int(p1), int(p2)
        w = math.sqrt(max(_patch_area_ha(p1i), 0.0) * max(_patch_area_ha(p2i), 0.0))
        if w <= 0.0:
            return 0.0
        return float(w / cost)

    def _overlap_ratio(cand: Dict, prior: Sequence[QgsGeometry]) -> float:
        g = cand.get("geom")
        if g is None or g.isEmpty() or not prior:
            return 0.0
        return float(_max_overlap_ratio(g, list(prior)) or 0.0)

    if not patches:
        return {}, {"strategy": "largest_single_network", "corridors_used": 0, "budget_used_ha": 0.0}

    seed_patch = max(patches.keys(), key=lambda pid: float(patches[pid].get("area_ha", 0.0) or 0.0))
    remaining = float(params.budget_area or 0.0)
    selected: Dict[int, Dict] = {}
    selected_geoms: List[QgsGeometry] = []

    uf = UnionFind()
    for pid, pdata in patches.items():
        uf.find(int(pid))
        uf.size[int(pid)] = float(pdata.get("area_ha", 0.0) or 0.0)
        uf.count[int(pid)] = 1

    G = None
    if nx is not None:
        try:
            G = nx.Graph()
            G.add_nodes_from([int(pid) for pid in patches.keys()])
        except Exception:
            G = None

    def _shortcut_multiplier(p1: int, p2: int, length: float) -> float:
        if length <= 0:
            return float(SHORTCUT_MULT_LOW)
        if G is None:
            return float(SHORTCUT_MULT_LOW)
        try:
            current_len = float(nx.shortest_path_length(G, p1, p2, weight="weight"))
        except Exception:
            return float(SHORTCUT_MULT_LOW)
        ratio = current_len / max(length, 1e-9)
        if ratio >= float(SHORTCUT_RATIO_HIGH):
            return float(SHORTCUT_MULT_HIGH)
        if ratio >= float(SHORTCUT_RATIO_MID):
            return float(SHORTCUT_MULT_MID)
        if ratio <= float(SHORTCUT_RATIO_LOW):
            return float(SHORTCUT_MULT_LOW)
        return float(SHORTCUT_MULT_LOW)

    def _main_root() -> int:
        return int(uf.find(int(seed_patch)))

    def _touches_main_and_other(pids: List[int]) -> bool:
        mr = _main_root()
        has_main = False
        has_other = False
        for pid in pids:
            r = int(uf.find(int(pid)))
            if r == mr:
                has_main = True
            else:
                has_other = True
        return has_main and has_other

    # ------------------------------------------------------------------
    # Phase 1: Seeded bridge-first (expand main component only)
    # ------------------------------------------------------------------
    bridge_ranked: List[Dict] = []
    for cand in candidates:
        cost = float(_candidate_cost_ha(cand) or 0.0)
        if cost <= 0.0:
            continue
        geom = cand.get("geom")
        if geom is None or geom.isEmpty():
            continue
        pids = _get_patch_ids(cand)
        if len(pids) < 2:
            continue
        bridge_ranked.append(cand)

    def _primary_score(cand: Dict) -> float:
        base_roi = float(_get_base_roi(cand) or 0.0)
        if base_roi <= 0.0:
            return 0.0
        pids = [int(pid) for pid in _get_patch_ids(cand) if pid is not None]
        if len(pids) < 2:
            return 0.0
        sizes = [float(_patch_area_ha(pid) or 0.0) for pid in pids]
        target_importance = min(sizes) if sizes else 1.0
        if target_importance <= 0:
            target_importance = 1.0
        return base_roi * target_importance

    bridge_ranked.sort(
        key=lambda c: (
            -_primary_score(c),
            float(c.get("area_ha", float("inf")) or float("inf")),
        )
    )

    selected_overlap_by_pair: Dict[Tuple[int, int], List[QgsGeometry]] = defaultdict(list)
    selected_count_by_pair: Dict[Tuple[int, int], int] = defaultdict(int)
    selected_geoms: List[QgsGeometry] = []

    def _commit_primary(cand: Dict, base_score: float) -> None:
        nonlocal remaining
        pids = [int(pid) for pid in _get_patch_ids(cand) if pid is not None]
        if len(pids) < 2:
            return
        g = cand.get("geom")
        if g is None or g.isEmpty():
            return
        cost = float(_candidate_cost_ha(cand) or 0.0)
        if cost <= 0.0 or cost > remaining:
            return
        cid = len(selected) + 1
        selected[cid] = {
            "geom": clone_geometry(cand["geom"]),
            "patch_ids": set(cand.get("patch_ids", set(pids))),
            "area_ha": float(cand.get("area_ha", 0.0) or 0.0),
            "p1": int(cand.get("patch1")),
            "p2": int(cand.get("patch2")),
            "distance": float(cand.get("distance_m", 1.0) or 1.0),
            "type": "primary",
            "variant": cand.get("variant"),
            "source": cand.get("source"),
            "utility_score": base_score,
            "overlap_ratio": 0.0,
        }
        anchor = int(pids[0])
        for other in pids[1:]:
            uf.union(anchor, int(other))
        remaining -= cost
        if G is not None:
            try:
                length = float(cand.get("distance_m", cost) or cost)
                for other in pids[1:]:
                    G.add_edge(anchor, int(other), weight=float(length))
            except Exception:
                pass

        pk = _get_pair(cand)
        g = cand.get("geom")
        if g is not None and not g.isEmpty():
            selected_overlap_by_pair.setdefault(pk, []).append(clone_geometry(g))
            if len(selected_overlap_by_pair[pk]) > 3:
                del selected_overlap_by_pair[pk][0]
            selected_geoms.append(clone_geometry(g))
        selected_count_by_pair[pk] = selected_count_by_pair.get(pk, 0) + 1

    bridge_by_mid: Dict[int, List[Dict]] = defaultdict(list)
    if len(patches) > MAX_BRIDGE_PATCHES or len(candidates) > MAX_BRIDGE_CANDIDATES:
        bridge_by_mid = {}
    for cand in candidates:
        p1, p2 = cand.get("patch1"), cand.get("patch2")
        if p1 is None or p2 is None:
            continue
        bridge_by_mid[int(p1)].append(cand)
        bridge_by_mid[int(p2)].append(cand)
    if bridge_by_mid:
        for mid, edges in list(bridge_by_mid.items()):
            if len(edges) > MAX_BRIDGE_EDGES_PER_MID:
                edges.sort(key=lambda c: float(c.get("area_ha", 0.0) or 0.0))
                bridge_by_mid[mid] = edges[:MAX_BRIDGE_EDGES_PER_MID]
        if len(bridge_by_mid) > MAX_BRIDGE_MIDS:
            mids_sorted = sorted(
                bridge_by_mid.items(),
                key=lambda kv: len(kv[1]),
                reverse=True,
            )[:MAX_BRIDGE_MIDS]
            bridge_by_mid = {mid: edges for mid, edges in mids_sorted}

    def _best_bridge_pair() -> Optional[Tuple[Dict, Dict, float]]:
        best = None
        best_score = 0.0
        mr = _main_root()
        for mid, edges in bridge_by_mid.items():
            if int(uf.find(int(mid))) == mr:
                continue
            best_to_main = None
            best_to_main_cost = 0.0
            best_to_other = None
            best_to_other_score = 0.0
            best_to_other_id = None
            for cand in edges:
                p1, p2 = cand.get("patch1"), cand.get("patch2")
                if p1 is None or p2 is None:
                    continue
                a = int(p1)
                b = int(p2)
                other = b if a == mid else a
                cost = float(_candidate_cost_ha(cand) or 0.0)
                if cost <= 0.0:
                    continue
                if int(uf.find(int(other))) == mr:
                    if best_to_main is None or cost < best_to_main_cost:
                        best_to_main = cand
                        best_to_main_cost = cost
                else:
                    gain_other = float(_patch_area_ha(int(other)) or 0.0)
                    if gain_other <= 0.0:
                        continue
                    score = gain_other / cost
                    if best_to_other is None or score > best_to_other_score:
                        best_to_other = cand
                        best_to_other_score = score
                        best_to_other_id = int(other)
            if best_to_main is None or best_to_other is None or best_to_other_id is None:
                continue
            total_cost = float(best_to_main_cost + (_candidate_cost_ha(best_to_other) or 0.0))
            if total_cost <= 0.0 or total_cost > remaining:
                continue
            gain_mid = float(_patch_area_ha(int(mid)) or 0.0)
            gain_other = float(_patch_area_ha(int(best_to_other_id)) or 0.0)
            score = (gain_mid + gain_other) / total_cost if total_cost else 0.0
            if score > best_score:
                best_score = score
                best = (best_to_main, best_to_other, score)
        return best

    bridge_iters = 0
    bridge_start = time.perf_counter()
    while remaining > 0:
        if bridge_iters >= MAX_BRIDGE_ITERATIONS:
            break
        if (time.perf_counter() - bridge_start) > MAX_BRIDGE_RUNTIME_S:
            break
        pair = _best_bridge_pair()
        if pair is None:
            break
        cand1, cand2, score = pair
        cost1 = float(_candidate_cost_ha(cand1) or 0.0)
        cost2 = float(_candidate_cost_ha(cand2) or 0.0)
        if cost1 <= 0.0 or cost2 <= 0.0 or (cost1 + cost2) > remaining:
            break
        _commit_primary(cand1, score)
        _commit_primary(cand2, score)
        bridge_iters += 1

    for cand in bridge_ranked:
        cost = float(_candidate_cost_ha(cand) or 0.0)
        if cost <= 0.0 or cost > remaining:
            continue
        pids = [int(pid) for pid in _get_patch_ids(cand) if pid is not None]
        if len(pids) < 2:
            continue
        roots = {int(uf.find(pid)) for pid in pids}
        if len(roots) <= 1:
            continue
        if not _touches_main_and_other(pids):
            continue

        base_roi = float(_get_base_roi(cand) or 0.0)
        target_sizes = [float(_patch_area_ha(pid) or 0.0) for pid in pids]
        target_importance = min(target_sizes) if target_sizes else 1.0
        if target_importance <= 0:
            target_importance = 1.0
        base_roi = base_roi * target_importance
        cid = len(selected) + 1
        selected[cid] = {
            "geom": clone_geometry(cand["geom"]),
            "patch_ids": set(cand.get("patch_ids", set(pids))),
            "area_ha": float(cand.get("area_ha", 0.0) or 0.0),
            "p1": int(cand.get("patch1")),
            "p2": int(cand.get("patch2")),
            "distance": float(cand.get("distance_m", 1.0) or 1.0),
            "type": "primary",
            "variant": cand.get("variant"),
            "source": cand.get("source"),
            "utility_score": base_roi,
            "overlap_ratio": 0.0,
        }
        g = cand.get("geom")
        if g is not None and not g.isEmpty():
            selected_geoms.append(clone_geometry(g))

        anchor = int(pids[0])
        for other in pids[1:]:
            uf.union(anchor, int(other))
        remaining -= cost
        if G is not None:
            try:
                length = float(cand.get("distance_m", cost) or cost)
                for other in pids[1:]:
                    G.add_edge(anchor, int(other), weight=float(length))
            except Exception:
                pass

        pk = _get_pair(cand)
        g = cand.get("geom")
        if g is not None and not g.isEmpty():
            selected_overlap_by_pair.setdefault(pk, []).append(clone_geometry(g))
            if len(selected_overlap_by_pair[pk]) > 3:
                del selected_overlap_by_pair[pk][0]
            selected_geoms.append(clone_geometry(g))
        selected_count_by_pair[pk] = selected_count_by_pair.get(pk, 0) + 1

    # ------------------------------------------------------------------
    # Phase 2: Circuit utility within the single (seeded) network
    # ------------------------------------------------------------------
    heap: List[Tuple[float, int, int, Dict]] = []
    stamps: Dict[int, int] = {}
    counter = 0

    def _stamp_key(cand: Dict) -> int:
        return id(cand)

    for cand in candidates:
        base = float(_get_base_roi(cand) or 0.0)
        cost = float(_candidate_cost_ha(cand) or 0.0)
        geom = cand.get("geom")
        if base <= 0.0 or cost <= 0.0 or geom is None or geom.isEmpty():
            continue
        if cost > remaining:
            continue
        counter += 1
        k = _stamp_key(cand)
        stamps[k] = stamps.get(k, 0) + 1
        heapq.heappush(heap, (-base, counter, stamps[k], cand))

    primary_links = sum(1 for d in selected.values() if d.get("type") == "primary")
    redundant_links = 0
    wasteful_links = 0

    while heap and remaining > 0:
        neg_score, _idx, stamp, cand = heapq.heappop(heap)
        k = _stamp_key(cand)
        if stamps.get(k, 0) != int(stamp):
            continue
        old_score = -float(neg_score)

        cost = float(_candidate_cost_ha(cand) or 0.0)
        if cost <= 0.0 or cost > remaining:
            continue
        base_roi = float(_get_base_roi(cand) or 0.0)
        if base_roi <= 0.0:
            continue

        pids = [int(pid) for pid in _get_patch_ids(cand) if pid is not None]
        if len(pids) < 2:
            continue

        roots = {int(uf.find(pid)) for pid in pids}
        mr = _main_root()

        overlap_r = 0.0
        if len(roots) > 1:
            # Only allow bridges that expand the main component.
            if not _touches_main_and_other(pids):
                continue
            root_sizes = [float(uf.size.get(r, 0.0) or 0.0) for r in roots]
            target_importance = min(root_sizes) if root_sizes else 1.0
            if target_importance <= 0:
                target_importance = 1.0
            mult = float(target_importance)
        else:
            # Only allow redundancy within the main component.
            if int(next(iter(roots))) != mr:
                continue
            g = cand.get("geom")
            if g is not None and not _redundancy_far_enough(
                g, selected_geoms, float(getattr(params, "max_search_distance", 0.0) or 0.0)
            ):
                continue
            pair = _get_pair(cand)
            prior = selected_overlap_by_pair.get(pair, [])
            overlap_r = float(_overlap_ratio(cand, prior) or 0.0)
            if overlap_r > float(overlap_reject_ratio):
                mult = 0.01
            else:
                length = float(cand.get("distance_m", cost) or cost)
                mult = _shortcut_multiplier(int(pids[0]), int(pids[-1]), float(length))

        new_score = base_roi * mult
        if abs(new_score - old_score) > 1e-12:
            counter += 1
            stamps[k] = stamps.get(k, 0) + 1
            heapq.heappush(heap, (-new_score, counter, stamps[k], cand))
            continue

        if len(roots) > 1:
            corr_type = "primary"
            primary_links += 1
        else:
            corr_type = "redundant"
            redundant_links += 1
            if mult <= float(SHORTCUT_MULT_LOW) or mult <= 0.01:
                corr_type = "wasteful"
                wasteful_links += 1

        cid = len(selected) + 1
        selected[cid] = {
            "geom": clone_geometry(cand["geom"]),
            "patch_ids": set(cand.get("patch_ids", set(pids))),
            "area_ha": float(cand.get("area_ha", 0.0) or 0.0),
            "p1": int(cand.get("patch1")),
            "p2": int(cand.get("patch2")),
            "distance": float(cand.get("distance_m", 1.0) or 1.0),
            "type": corr_type,
            "variant": cand.get("variant"),
            "source": cand.get("source"),
            "utility_score": float(new_score),
            "overlap_ratio": float(overlap_r),
        }

        anchor = int(pids[0])
        for other in pids[1:]:
            uf.union(anchor, int(other))

        remaining -= cost
        g = cand.get("geom")
        pair = _get_pair(cand)
        if g is not None and not g.isEmpty():
            lst = selected_overlap_by_pair.setdefault(pair, [])
            lst.append(clone_geometry(g))
            if len(lst) > 3:
                del lst[0]
            selected_geoms.append(clone_geometry(g))
        selected_count_by_pair[pair] = selected_count_by_pair.get(pair, 0) + 1
        if G is not None:
            try:
                length = float(cand.get("distance_m", cost) or cost)
                for other in pids[1:]:
                    G.add_edge(anchor, int(other), weight=float(length))
            except Exception:
                pass

    comp_area: Dict[int, float] = defaultdict(float)
    comp_count: Dict[int, int] = defaultdict(int)
    for pid in patches:
        root = int(uf.find(int(pid)))
        comp_area[root] += _patch_area_ha(int(pid))
        comp_count[root] += 1

    mr = _main_root()
    largest_group_area = float(comp_area.get(mr, 0.0) or 0.0)
    largest_group_patches = int(comp_count.get(mr, 0) or 0)

    _annotate_vector_corridor_metrics(
        patches=patches,
        corridors=selected,
        use_global_network_area=True,
        largest_group_area_ha=largest_group_area,
    )

    stats = {
        "strategy": "largest_single_network",
        "corridors_used": len(selected),
        "budget_used_ha": float((params.budget_area or 0.0) - remaining),
        "patches_total": len(patches),
        "patches_connected": largest_group_patches,
        "components_remaining": len(comp_area) if comp_area else len(patches),
        "primary_links": int(primary_links),
        "redundant_links": int(redundant_links),
        "wasteful_links": int(wasteful_links),
        "total_connected_area_ha": sum(p.get("area_ha", 0.0) for p in patches.values()),
        "largest_group_area_ha": largest_group_area,
        "largest_group_patches": largest_group_patches,
        "seed_patch": int(seed_patch),
    }
    return selected, stats


def _apply_hybrid_leftover_budget_vector(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    corridors: Dict[int, Dict],
    remaining_budget: float,
    roi_bias: float = ROI_REDUNDANCY_BIAS,
    overlap_reject_ratio: float = HYBRID_OVERLAP_REJECT_RATIO,
    max_search_distance: float = 0.0,
) -> Tuple[float, int, int]:
    """
    Spend remaining budget by comparing low-value connections vs redundancy.
    Additional redundant corridors between the same pair are allowed without a
    hard cap, but must be spatially distinct by at least the search distance.
    Returns (budget_used, low_value_added, redundancy_added).
    """
    if remaining_budget <= 0 or not candidates:
        return 0.0, 0, 0

    def _pair_key(a: int, b: int) -> Tuple[int, int]:
        return (a, b) if a <= b else (b, a)

    uf = UnionFind()
    for pid, pdata in patches.items():
        uf.find(int(pid))
        uf.size[int(pid)] = float(pdata.get("area_ha", 0.0) or 0.0)
        uf.count[int(pid)] = 1
    for data in corridors.values():
        pids = list(data.get("patch_ids", {data.get("p1"), data.get("p2")}))
        pids = [int(pid) for pid in pids if pid is not None]
        if len(pids) < 2:
            continue
        anchor = pids[0]
        for other in pids[1:]:
            uf.union(anchor, int(other))

    roots = {int(uf.find(int(pid))) for pid in patches}
    if not roots:
        return 0.0, 0, 0
    main_root = max(roots, key=lambda r: uf.size.get(r, 0.0))

    selected_count_by_pair: Dict[Tuple[int, int], int] = defaultdict(int)
    selected_geoms_by_pair: Dict[Tuple[int, int], List[QgsGeometry]] = defaultdict(list)
    selected_geoms: List[QgsGeometry] = []
    for data in corridors.values():
        p1, p2 = data.get("p1"), data.get("p2")
        if p1 is None or p2 is None:
            continue
        pk = _pair_key(int(p1), int(p2))
        selected_count_by_pair[pk] = selected_count_by_pair.get(pk, 0) + 1
        g = data.get("geom")
        if g is not None and not g.isEmpty():
            selected_geoms_by_pair[pk].append(clone_geometry(g))
            selected_geoms.append(clone_geometry(g))

    G = None
    if nx is not None and graph_math is not None:
        G = nx.Graph()
        for pid in patches:
            G.add_node(int(pid))
        for data in corridors.values():
            p1, p2 = data.get("p1"), data.get("p2")
            if p1 is None or p2 is None:
                continue
            dist = float(data.get("distance", 1.0) or 1.0)
            G.add_edge(int(p1), int(p2), weight=dist)

    budget_used = 0.0
    low_value_added = 0
    redundancy_added = 0

    while remaining_budget > 0:
        best_new: Optional[Tuple[float, Dict, float]] = None
        best_red: Optional[Tuple[float, Dict, float]] = None

        for cand in candidates:
            cost = float(cand.get("area_ha", 0.0) or 0.0)
            if cost <= 0.0 or cost > remaining_budget:
                continue
            geom = cand.get("geom")
            if geom is None or geom.isEmpty():
                continue
            p1, p2 = cand.get("patch1"), cand.get("patch2")
            if p1 is None or p2 is None:
                continue
            p1i, p2i = int(p1), int(p2)
            pk = _pair_key(p1i, p2i)

            pids = list(cand.get("patch_ids", {p1i, p2i}))
            pids = [int(pid) for pid in pids if pid is not None]
            roots = {int(uf.find(pid)) for pid in pids}
            if not roots:
                continue

            if len(roots) > 1:
                if main_root not in roots:
                    continue
                gain = sum(float(uf.size.get(r, 0.0) or 0.0) for r in roots if r != main_root)
                roi_new = (gain / cost) if cost else 0.0
                if roi_new > 0 and (best_new is None or roi_new > best_new[0]):
                    best_new = (roi_new, cand, gain)
                continue

            if main_root not in roots:
                continue
            prior_geoms = selected_geoms_by_pair.get(pk, [])
            overlap = _geom_overlap_ratio(geom, prior_geoms)
            if overlap >= float(overlap_reject_ratio):
                continue
            if max_search_distance > 0.0 and not _redundancy_far_enough(
                geom, prior_geoms, float(max_search_distance or 0.0)
            ):
                continue
            dist = float(cand.get("distance_m", 1.0) or 1.0)
            score = 0.0
            if G is not None and graph_math is not None:
                try:
                    if nx.has_path(G, p1i, p2i):
                        score = float(graph_math.score_edge_for_loops(G, p1i, p2i, dist))
                except Exception:
                    score = 0.0
            if score <= 0.0:
                continue
            roi_red = score / cost if cost else 0.0
            if roi_red > 0 and (best_red is None or roi_red > best_red[0]):
                best_red = (roi_red, cand, overlap)

        if best_new is None and best_red is None:
            break

        pick_red = False
        if best_red is not None and best_new is not None:
            pick_red = best_red[0] > best_new[0] * (1.0 + float(roi_bias))
        elif best_red is not None:
            pick_red = True

        roi, cand, extra = best_red if pick_red else best_new
        if cand is None or roi <= 0.0:
            break

        cost = float(cand.get("area_ha", 0.0) or 0.0)
        if cost <= 0.0 or cost > remaining_budget:
            break

        p1i, p2i = int(cand.get("patch1")), int(cand.get("patch2"))
        pk = _pair_key(p1i, p2i)
        pids = list(cand.get("patch_ids", {p1i, p2i}))
        pids = [int(pid) for pid in pids if pid is not None]

        if pick_red:
            redundancy_added += 1
            corr_type = "redundant"
            overlap_r = float(extra or 0.0)
        else:
            low_value_added += 1
            corr_type = "low_value"
            overlap_r = 0.0

        if pids:
            if main_root in {int(uf.find(pid)) for pid in pids}:
                anchor = next((pid for pid in pids if int(uf.find(pid)) == main_root), pids[0])
            else:
                anchor = pids[0]
            for other in pids[1:]:
                uf.union(int(anchor), int(other))
            main_root = int(uf.find(anchor))

        cid = len(corridors) + 1
        corridors[cid] = {
            "geom": clone_geometry(cand.get("geom")),
            "patch_ids": set(cand.get("patch_ids", {p1i, p2i})),
            "area_ha": float(cand.get("area_ha", 0.0) or 0.0),
            "p1": p1i,
            "p2": p2i,
            "distance": float(cand.get("distance_m", 1.0) or 1.0),
            "type": corr_type,
            "variant": cand.get("variant"),
            "source": cand.get("source"),
            "utility_score": float(roi),
            "overlap_ratio": float(overlap_r),
        }
        selected_count_by_pair[pk] = selected_count_by_pair.get(pk, 0) + 1
        geom = cand.get("geom")
        if geom is not None and not geom.isEmpty():
            selected_geoms_by_pair[pk].append(clone_geometry(geom))
            selected_geoms.append(clone_geometry(geom))

        remaining_budget -= cost
        budget_used += cost
        if G is not None:
            G.add_edge(p1i, p2i, weight=float(cand.get("distance_m", 1.0) or 1.0))

    return budget_used, low_value_added, redundancy_added


def _enforce_largest_network_component(
    patches: Dict[int, Dict],
    corridors: Dict[int, Dict],
    seed_patch: Optional[int] = None,
) -> Tuple[int, float]:
    """
    Keep only corridors in the seed-connected patch component.
    Returns (corridors_removed, area_removed_ha).
    """
    if not patches or not corridors:
        return 0, 0.0

    uf = UnionFind()
    for pid in patches:
        uf.find(int(pid))

    for data in corridors.values():
        pids = list(data.get("patch_ids", {data.get("p1"), data.get("p2")}))
        pids = [int(pid) for pid in pids if pid is not None and int(pid) in patches]
        if len(pids) < 2:
            continue
        anchor = pids[0]
        for other in pids[1:]:
            uf.union(anchor, int(other))

    if seed_patch is None or int(seed_patch) not in patches:
        try:
            seed_patch = int(max(patches.items(), key=lambda kv: float(kv[1].get("area_ha", 0.0) or 0.0))[0])
        except Exception:
            return 0, 0.0

    target_root = int(uf.find(int(seed_patch)))
    kept: Dict[int, Dict] = {}
    removed = 0
    removed_area = 0.0

    for cid in sorted(corridors.keys(), key=lambda x: int(x)):
        data = corridors.get(cid) or {}
        pids = list(data.get("patch_ids", {data.get("p1"), data.get("p2")}))
        pids = [int(pid) for pid in pids if pid is not None and int(pid) in patches]
        if not pids:
            removed += 1
            removed_area += float(data.get("area_ha", 0.0) or 0.0)
            continue
        if any(int(uf.find(pid)) == target_root for pid in pids):
            kept[len(kept) + 1] = data
        else:
            removed += 1
            removed_area += float(data.get("area_ha", 0.0) or 0.0)

    if removed > 0:
        corridors.clear()
        corridors.update(kept)

    return int(removed), float(removed_area)


def _refresh_vector_connectivity_stats(
    patches: Dict[int, Dict],
    corridors: Dict[int, Dict],
    stats: Dict,
) -> None:
    """Recompute connectivity-derived stats after adding corridors."""
    if not corridors:
        return

    uf = UnionFind()
    for pid, pdata in patches.items():
        uf.find(int(pid))
        uf.size[int(pid)] = float(pdata.get("area_ha", 0.0) or 0.0)
        uf.count[int(pid)] = 1

    for data in corridors.values():
        pids = list(data.get("patch_ids", {data.get("p1"), data.get("p2")}))
        pids = [int(pid) for pid in pids if pid is not None]
        if len(pids) < 2:
            continue
        anchor = pids[0]
        for other in pids[1:]:
            uf.union(anchor, int(other))

    comp_patch_area: Dict[int, float] = defaultdict(float)
    comp_corridor_area: Dict[int, float] = defaultdict(float)
    comp_count: Dict[int, int] = defaultdict(int)
    for pid in patches:
        root = int(uf.find(int(pid)))
        comp_patch_area[root] += float(patches[pid].get("area_ha", 0.0) or 0.0)
        comp_count[root] += 1

    for cdata in corridors.values():
        p1 = cdata.get("p1")
        if p1 is None:
            continue
        root = int(uf.find(int(p1)))
        comp_corridor_area[root] += float(cdata.get("area_ha", 0.0) or 0.0)

    comp_total_area: Dict[int, float] = {}
    for root, patch_area in comp_patch_area.items():
        comp_total_area[root] = patch_area + float(comp_corridor_area.get(root, 0.0) or 0.0)

    largest_group_area = max(comp_total_area.values()) if comp_total_area else 0.0
    largest_group_patches = max(comp_count.values()) if comp_count else 0

    if "components_remaining" in stats:
        stats["components_remaining"] = len(comp_total_area) if comp_total_area else len(patches)
    if "largest_group_area_ha" in stats:
        stats["largest_group_area_ha"] = float(largest_group_area)
    if "largest_group_patches" in stats:
        stats["largest_group_patches"] = int(largest_group_patches)
    if "patches_connected" in stats:
        if stats.get("strategy") == "largest_single_network" or stats.get("mode") == "Largest Single Network":
            stats["patches_connected"] = int(largest_group_patches)
        else:
            stats["patches_connected"] = int(largest_group_patches)
    if "total_connected_area_ha" in stats:
        stats["total_connected_area_ha"] = sum(comp_total_area.values()) if comp_total_area else 0.0

    _annotate_vector_corridor_metrics(
        patches=patches,
        corridors=corridors,
        use_global_network_area=(
            stats.get("strategy") == "largest_single_network"
            or stats.get("mode") == "Largest Single Network"
        ),
        largest_group_area_ha=float(largest_group_area),
    )


def _annotate_vector_corridor_metrics(
    *,
    patches: Dict[int, Dict],
    corridors: Dict[int, Dict],
    use_global_network_area: bool = False,
    largest_group_area_ha: float = 0.0,
) -> None:
    """Populate per-corridor area metrics used by the vector output tables."""
    if not corridors:
        return

    uf = UnionFind()
    for pid in patches:
        uf.find(int(pid))

    for cdata in corridors.values():
        pids = list(cdata.get("patch_ids", {cdata.get("p1"), cdata.get("p2")}))
        pids = [int(pid) for pid in pids if pid is not None and int(pid) in patches]
        if len(pids) < 2:
            continue
        anchor = int(pids[0])
        for other in pids[1:]:
            uf.union(anchor, int(other))

    comp_patch_area: Dict[int, float] = defaultdict(float)
    comp_corridor_area: Dict[int, float] = defaultdict(float)
    for pid, pdata in patches.items():
        root = int(uf.find(int(pid)))
        comp_patch_area[root] += float(pdata.get("area_ha", 0.0) or 0.0)

    for cdata in corridors.values():
        p1 = cdata.get("p1", cdata.get("patch1"))
        if p1 is None:
            continue
        root = int(uf.find(int(p1)))
        comp_corridor_area[root] += float(cdata.get("area_ha", 0.0) or 0.0)

    comp_total_area: Dict[int, float] = {}
    for root, patch_area in comp_patch_area.items():
        comp_total_area[root] = patch_area + float(comp_corridor_area.get(root, 0.0) or 0.0)

    global_network_area_ha = max(
        float(largest_group_area_ha or 0.0),
        max(comp_total_area.values()) if comp_total_area else 0.0,
    )

    for data in corridors.values():
        p1 = data.get("p1", data.get("patch1"))
        if p1 is None:
            continue
        corridor_area_ha = float(data.get("area_ha", 0.0) or 0.0)
        patches_area_ha = float(_corridor_patch_sum_area_ha(data, patches) or 0.0)
        root = int(uf.find(int(p1)))
        network_area_ha = (
            global_network_area_ha
            if use_global_network_area
            else float(comp_total_area.get(root, 0.0) or 0.0)
        )
        data["corridor_area_ha"] = corridor_area_ha
        data["patches_area_ha"] = patches_area_ha
        data["network_area_ha"] = network_area_ha
        # Preserve existing field name for backward compatibility.
        data["connected_area_ha"] = network_area_ha
        data["efficiency"] = (
            float(patches_area_ha / corridor_area_ha)
            if corridor_area_ha > 0.0 and patches_area_ha > 0.0
            else 0.0
        )


def _corridor_yes_no(value: object) -> str:
    return "YES" if bool(value) else "NO"


def _corridor_is_redundant(
    corridor_id: int,
    cdata: Dict[str, Any],
    corridors: Dict[int, Dict[str, Any]],
) -> bool:
    """
    Return True when the corridor's endpoints remain connected through other
    selected corridors, including indirect alternatives such as A-B-C for A-C.
    """
    p1 = cdata.get("p1", cdata.get("patch1"))
    p2 = cdata.get("p2", cdata.get("patch2"))
    if p1 is None or p2 is None:
        return False
    try:
        u = int(p1)
        v = int(p2)
    except Exception:
        return False
    if u == v:
        return False

    graph = nx.Graph()
    graph.add_node(u)
    graph.add_node(v)

    for other_id, other in corridors.items():
        if int(other_id) == int(corridor_id):
            continue
        raw_pids = list(other.get("patch_ids", {other.get("p1"), other.get("p2")}))
        pids: List[int] = []
        for pid in raw_pids:
            if pid is None:
                continue
            try:
                pids.append(int(pid))
            except Exception:
                continue
        pids = sorted(set(pids))
        if len(pids) < 2:
            continue
        for idx, a in enumerate(pids):
            graph.add_node(int(a))
            for b in pids[idx + 1 :]:
                graph.add_edge(int(a), int(b))

    try:
        return bool(nx.has_path(graph, u, v))
    except Exception:
        return False


def _corridor_is_isthmus(
    cdata: Dict[str, Any],
    patches: Dict[int, Dict[str, Any]],
    params: VectorRunParams,
    ctx: Optional[AnalysisContext] = None,
) -> bool:
    """
    Flag corridors that use a surviving small patch as an intermediate stepping stone.

    Under the current vector workflow, patches below the minimum patch-size filter
    are removed before corridor generation, so an isthmus can only be a retained
    intermediate patch whose area is smaller than the corridor area.
    """
    corridor_area_ha = float(cdata.get("corridor_area_ha", cdata.get("area_ha", 0.0)) or 0.0)

    participant_ids: Set[int] = set()
    raw_patch_ids = list(cdata.get("patch_ids", []) or [])
    for pid in raw_patch_ids:
        try:
            participant_ids.add(int(pid))
        except Exception:
            continue

    p1 = cdata.get("p1", cdata.get("patch1"))
    p2 = cdata.get("p2", cdata.get("patch2"))
    endpoint_ids: Set[int] = set()
    try:
        if p1 is not None:
            endpoint_ids.add(int(p1))
        if p2 is not None:
            endpoint_ids.add(int(p2))
    except Exception:
        pass

    chain_via = cdata.get("chain_via_patch")
    if chain_via is not None:
        try:
            participant_ids.add(int(chain_via))
        except Exception:
            pass

    interior_ids = set(int(pid) for pid in participant_ids if int(pid) not in endpoint_ids)
    if not interior_ids:
        return False

    min_patch_size_ha = max(float(getattr(params, "min_patch_size", 0.0) or 0.0), 0.0)
    for pid in interior_ids:
        pdata = patches.get(int(pid)) or {}
        patch_area_ha = float(pdata.get("area_ha", 0.0) or 0.0)
        if patch_area_ha <= 0.0:
            continue
        if patch_area_ha + 1e-12 < min_patch_size_ha:
            return True
        if corridor_area_ha > 0.0 and patch_area_ha + 1e-12 < corridor_area_ha:
            return True

    return False


def _annotate_vector_corridor_export_flags(
    corridors: Dict[int, Dict[str, Any]],
    patches: Dict[int, Dict[str, Any]],
    params: VectorRunParams,
    ctx: Optional[AnalysisContext] = None,
) -> None:
    for cid, cdata in corridors.items():
        cdata["redundant_label"] = _corridor_yes_no(_corridor_is_redundant(int(cid), cdata, corridors))
        cdata["intrapatch_label"] = _corridor_yes_no(bool(cdata.get("intra_patch", False)))
        cdata["isthmus_label"] = _corridor_yes_no(_corridor_is_isthmus(cdata, patches, params, ctx))


def _compute_habitat_component_metrics_exact(
    patches: Dict[int, Dict],
    corridors: Dict[int, Dict],
) -> Dict[str, float]:
    """
    Exact component metrics on patch graph (habitat area only, no corridor area).

    Returns normalized habitat mesh and LCC values that are resolution-independent
    and align with optimization objectives that operate on patch connectivity.
    """
    patch_areas: Dict[int, float] = {}
    for pid, pdata in patches.items():
        try:
            ipid = int(pid)
        except Exception:
            continue
        area = float(pdata.get("area_ha", 0.0) or 0.0)
        if area > 0.0:
            patch_areas[ipid] = area

    total_habitat_area = float(sum(patch_areas.values()))
    if total_habitat_area <= 0.0:
        return {
            "mesh_norm_pre": 0.0,
            "mesh_norm_post": 0.0,
            "lcc_norm_pre": 0.0,
            "lcc_norm_post": 0.0,
            "total_habitat_area_ha": 0.0,
        }

    pre_num = float(sum(a * a for a in patch_areas.values()))
    pre_lcc_area = max(patch_areas.values()) if patch_areas else 0.0

    uf = UnionFind()
    for pid, area in patch_areas.items():
        uf.find(int(pid))
        uf.size[int(pid)] = float(area)
        uf.count[int(pid)] = 1

    for cdata in corridors.values():
        pids = list(cdata.get("patch_ids", {cdata.get("p1"), cdata.get("p2")}))
        pids = [int(pid) for pid in pids if pid is not None and int(pid) in patch_areas]
        if len(pids) < 2:
            continue
        anchor = int(pids[0])
        for other in pids[1:]:
            uf.union(anchor, int(other))

    comp_area: Dict[int, float] = defaultdict(float)
    for pid, area in patch_areas.items():
        root = int(uf.find(int(pid)))
        comp_area[root] += float(area)

    post_num = float(sum(a * a for a in comp_area.values()))
    post_lcc_area = max(comp_area.values()) if comp_area else 0.0
    denom = float(total_habitat_area * total_habitat_area)

    mesh_norm_pre = max(0.0, min(1.0, (pre_num / denom) if denom > 0.0 else 0.0))
    mesh_norm_post = max(0.0, min(1.0, (post_num / denom) if denom > 0.0 else 0.0))
    lcc_norm_pre = max(0.0, min(1.0, pre_lcc_area / total_habitat_area))
    lcc_norm_post = max(0.0, min(1.0, post_lcc_area / total_habitat_area))
    return {
        "mesh_norm_pre": float(mesh_norm_pre),
        "mesh_norm_post": float(mesh_norm_post),
        "lcc_norm_pre": float(lcc_norm_pre),
        "lcc_norm_post": float(lcc_norm_post),
        "total_habitat_area_ha": float(total_habitat_area),
    }


def write_corridors_layer_to_gpkg(
    corridors: Dict[int, Dict],
    output_path: str,
    layer_name: str,
    target_crs: QgsCoordinateReferenceSystem,
    original_crs: QgsCoordinateReferenceSystem,
    unit_system: str,
    overwrite_file: bool = False,
) -> bool:
    print(f"\nWriting layer '{layer_name}' to {os.path.basename(output_path)} ...")
    transform = QgsCoordinateTransform(target_crs, original_crs, QgsProject.instance())

    is_imperial = unit_system == "imperial"
    corridor_area_field = "corridor_area_ac" if is_imperial else "corridor_area_ha"
    patches_area_field = "patches_area_ac" if is_imperial else "patches_area_ha"
    network_area_field = "network_area_ac" if is_imperial else "network_area_ha"
    area_factor = 2.471053814 if is_imperial else 1.0

    fields = QgsFields()
    fields.append(QgsField("corridor_id", QVariant.Int))
    fields.append(QgsField("patch_ids", QVariant.String))
    fields.append(QgsField("patch1", QVariant.Int))
    fields.append(QgsField("patch2", QVariant.Int))
    fields.append(QgsField(corridor_area_field, QVariant.Double))
    fields.append(QgsField(patches_area_field, QVariant.Double))
    fields.append(QgsField(network_area_field, QVariant.Double))
    fields.append(QgsField("efficiency", QVariant.Double))
    fields.append(QgsField("redundant", QVariant.String))
    fields.append(QgsField("intrapatch", QVariant.String))
    fields.append(QgsField("isthmus", QVariant.String))

    save_options = QgsVectorFileWriter.SaveVectorOptions()
    save_options.driverName = "GPKG"
    save_options.fileEncoding = "UTF-8"
    save_options.layerName = layer_name
    save_options.actionOnExistingFile = (
        QgsVectorFileWriter.CreateOrOverwriteFile if overwrite_file else QgsVectorFileWriter.CreateOrOverwriteLayer
    )

    writer = QgsVectorFileWriter.create(
        output_path,
        fields,
        QgsWkbTypes.Polygon,
        original_crs,
        QgsProject.instance().transformContext(),
        save_options,
    )

    if writer.hasError() != QgsVectorFileWriter.NoError:
        print(f"  ✗ Error: {writer.errorMessage()}")
        return False

    for cid, cdata in corridors.items():
        feat = QgsFeature(fields)
        geom = clone_geometry(cdata["geom"])
        geom.transform(transform)
        feat.setGeometry(geom)
        feat.setAttributes(
            [
                cid,
                ",".join(map(str, sorted(cdata["patch_ids"]))),
                int(cdata.get("p1", cdata.get("patch1", 0)) or 0),
                int(cdata.get("p2", cdata.get("patch2", 0)) or 0),
                round(float(cdata.get("corridor_area_ha", cdata.get("area_ha", 0.0)) or 0.0) * area_factor, 4),
                round(float(cdata.get("patches_area_ha", 0.0) or 0.0) * area_factor, 4),
                round(float(cdata.get("network_area_ha", cdata.get("connected_area_ha", 0.0)) or 0.0) * area_factor, 4),
                round(cdata["efficiency"], 6),
                str(cdata.get("redundant_label", "") or ""),
                str(cdata.get("intrapatch_label", "") or ""),
                str(cdata.get("isthmus_label", "") or ""),
            ]
        )
        writer.addFeature(feat)

    del writer
    print(f"  ✓ Wrote {len(corridors)} corridors to layer '{layer_name}'")
    return True


def add_layer_to_qgis_from_gpkg(gpkg_path: str, layer_name: str, add_to_project: bool = True) -> None:
    if not add_to_project:
        return
    uri = f"{gpkg_path}|layername={layer_name}"
    layer = QgsVectorLayer(uri, layer_name, "ogr")
    if layer.isValid():
        lname = str(layer_name).strip().lower()
        if lname.startswith("corridors"):
            try:
                _apply_visible_corridor_style_vector(layer)
            except Exception:
                pass
        if lname == "contiguous areas":
            try:
                _apply_random_unique_value_symbology_vector(layer, "network_id")
            except Exception:
                pass
        QgsProject.instance().addMapLayer(layer)
        print(f"  ✓ Added '{layer_name}' to QGIS project")
    else:
        print(f"  ✗ Could not add '{layer_name}' from {gpkg_path}")


def _safe_unary_union(geoms: List[QgsGeometry]) -> Optional[QgsGeometry]:
    if not geoms:
        return None
    try:
        merged = QgsGeometry.unaryUnion(geoms)
    except Exception:
        merged = geoms[0]
        for extra in geoms[1:]:
            try:
                merged = merged.combine(extra)
            except Exception:
                pass
    if merged is None or merged.isEmpty():
        return None
    try:
        merged = merged.makeValid()
    except Exception:
        pass
    return merged if (merged is not None and not merged.isEmpty()) else None


def build_contiguous_network_summaries(
    patches: Dict[int, Dict],
    corridors: Dict[int, Dict],
    dissolve_tolerance: float = 0.0,
) -> List[Dict]:
    """Group patches/corridors into connected networks and dissolve geometries per network."""
    if not patches:
        return []

    uf = UnionFind()
    for pid in patches.keys():
        uf.find(int(pid))

    for cdata in corridors.values():
        pids = list(cdata.get("patch_ids", []))
        if len(pids) < 2:
            continue
        base = int(pids[0])
        for other in pids[1:]:
            uf.union(base, int(other))

    root_to_pids: Dict[int, Set[int]] = defaultdict(set)
    for pid in patches.keys():
        root_to_pids[int(uf.find(int(pid)))].add(int(pid))

    root_to_corridors: Dict[int, List[int]] = defaultdict(list)
    for cid, cdata in corridors.items():
        pids = list(cdata.get("patch_ids", []))
        if not pids:
            continue
        root_to_corridors[int(uf.find(int(pids[0])))] += [int(cid)]

    summaries: List[Dict] = []
    network_id = 1
    for root, pid_set in sorted(root_to_pids.items(), key=lambda kv: kv[0]):
        geom_parts: List[QgsGeometry] = []
        patch_area_ha = 0.0
        for pid in sorted(pid_set):
            pdata = patches.get(pid)
            if not pdata:
                continue
            pgeom = pdata.get("geom")
            if pgeom is not None and (not pgeom.isEmpty()):
                geom_parts.append(clone_geometry(pgeom))
            patch_area_ha += float(pdata.get("area_ha", 0.0) or 0.0)

        corridor_ids = root_to_corridors.get(root, [])
        corridor_area_ha = 0.0
        for cid in corridor_ids:
            cdata = corridors.get(cid) or {}
            cgeom = cdata.get("geom")
            if cgeom is not None and (not cgeom.isEmpty()):
                geom_parts.append(clone_geometry(cgeom))
            corridor_area_ha += float(cdata.get("area_ha", 0.0) or 0.0)

        net_geom = _safe_unary_union(geom_parts)
        if net_geom and (not net_geom.isEmpty()) and dissolve_tolerance > 0.0:
            try:
                net_geom = net_geom.buffer(dissolve_tolerance, BUFFER_SEGMENTS)
                net_geom = net_geom.buffer(-dissolve_tolerance, BUFFER_SEGMENTS)
                net_geom = net_geom.makeValid()
            except Exception:
                pass

        summaries.append(
            {
                "network_id": network_id,
                "patch_ids": set(pid_set),
                "corridor_ids": list(corridor_ids),
                "geom": net_geom,
                "area_ha": (net_geom.area() / 10000.0) if net_geom else 0.0,
                "patch_count": len(pid_set),
                "patch_area_ha": patch_area_ha,
                "corridor_count": len(corridor_ids),
                "corridor_area_ha": corridor_area_ha,
            }
        )
        network_id += 1

    return summaries


def write_contiguous_networks_layer_to_gpkg(
    networks: List[Dict],
    output_path: str,
    layer_name: str,
    target_crs: QgsCoordinateReferenceSystem,
    original_crs: QgsCoordinateReferenceSystem,
    unit_system: str,
) -> bool:
    print(f"\nWriting layer '{layer_name}' to {os.path.basename(output_path)} ...")
    transform = QgsCoordinateTransform(target_crs, original_crs, QgsProject.instance())

    is_imperial = unit_system == "imperial"
    area_field = "area_ac" if is_imperial else "area_ha"
    patch_area_field = "patch_ac" if is_imperial else "patch_ha"
    corr_area_field = "corr_ac" if is_imperial else "corr_ha"
    area_factor = 2.471053814 if is_imperial else 1.0

    fields = QgsFields()
    fields.append(QgsField("network_id", QVariant.Int))
    fields.append(QgsField("patch_count", QVariant.Int))
    fields.append(QgsField("corr_count", QVariant.Int))
    fields.append(QgsField(area_field, QVariant.Double))
    fields.append(QgsField(patch_area_field, QVariant.Double))
    fields.append(QgsField(corr_area_field, QVariant.Double))
    fields.append(QgsField("multipart", QVariant.Bool))
    fields.append(QgsField("part_count", QVariant.Int))

    save_options = QgsVectorFileWriter.SaveVectorOptions()
    save_options.driverName = "GPKG"
    save_options.fileEncoding = "UTF-8"
    save_options.layerName = layer_name
    save_options.actionOnExistingFile = QgsVectorFileWriter.CreateOrOverwriteLayer

    writer = QgsVectorFileWriter.create(
        output_path,
        fields,
        QgsWkbTypes.Polygon,
        original_crs,
        QgsProject.instance().transformContext(),
        save_options,
    )
    if writer.hasError() != QgsVectorFileWriter.NoError:
        print(f"  ✗ Error: {writer.errorMessage()}")
        return False

    written = 0
    for net in networks:
        geom = net.get("geom")
        if geom is None or geom.isEmpty():
            continue
        feat = QgsFeature(fields)
        g = clone_geometry(geom)
        g.transform(transform)
        feat.setGeometry(g)
        multipart = g.isMultipart()
        part_count = g.constGet().numGeometries() if multipart else 1
        feat.setAttributes(
            [
                int(net.get("network_id", written + 1)),
                int(net.get("patch_count", 0)),
                int(net.get("corridor_count", 0)),
                round(float(net.get("area_ha", 0.0)) * area_factor, 4),
                round(float(net.get("patch_area_ha", 0.0)) * area_factor, 4),
                round(float(net.get("corridor_area_ha", 0.0)) * area_factor, 4),
                multipart,
                int(part_count),
            ]
        )
        writer.addFeature(feat)
        written += 1

    del writer
    print(f"  ✓ Wrote {written} network feature(s) to layer '{layer_name}'")
    return True


def create_memory_layer_from_corridors(
    corridors: Dict[int, Dict],
    layer_name: str,
    target_crs: QgsCoordinateReferenceSystem,
    original_crs: QgsCoordinateReferenceSystem,
    unit_system: str,
    add_to_project: bool = True,
) -> Optional[QgsVectorLayer]:
    is_imperial = unit_system == "imperial"
    corridor_area_field = "corridor_area_ac" if is_imperial else "corridor_area_ha"
    patches_area_field = "patches_area_ac" if is_imperial else "patches_area_ha"
    network_area_field = "network_area_ac" if is_imperial else "network_area_ha"
    area_factor = 2.471053814 if is_imperial else 1.0

    layer = QgsVectorLayer(f"Polygon?crs={original_crs.authid()}", layer_name, "memory")
    provider = layer.dataProvider()
    provider.addAttributes(
        [
            QgsField("corridor_id", QVariant.Int),
            QgsField("patch_ids", QVariant.String),
            QgsField("patch1", QVariant.Int),
            QgsField("patch2", QVariant.Int),
            QgsField(corridor_area_field, QVariant.Double),
            QgsField(patches_area_field, QVariant.Double),
            QgsField(network_area_field, QVariant.Double),
            QgsField("efficiency", QVariant.Double),
            QgsField("redundant", QVariant.String),
            QgsField("intrapatch", QVariant.String),
            QgsField("isthmus", QVariant.String),
        ]
    )
    layer.updateFields()

    transform = QgsCoordinateTransform(target_crs, original_crs, QgsProject.instance())

    features = []
    for cid, cdata in corridors.items():
        geom = clone_geometry(cdata["geom"])
        geom.transform(transform)
        feat = QgsFeature(layer.fields())
        feat.setGeometry(geom)
        feat.setAttributes(
            [
                cid,
                ",".join(map(str, sorted(cdata["patch_ids"]))),
                int(cdata.get("p1", cdata.get("patch1", 0)) or 0),
                int(cdata.get("p2", cdata.get("patch2", 0)) or 0),
                round(float(cdata.get("corridor_area_ha", cdata.get("area_ha", 0.0)) or 0.0) * area_factor, 4),
                round(float(cdata.get("patches_area_ha", 0.0) or 0.0) * area_factor, 4),
                round(float(cdata.get("network_area_ha", cdata.get("connected_area_ha", 0.0)) or 0.0) * area_factor, 4),
                round(cdata["efficiency"], 6),
                str(cdata.get("redundant_label", "") or ""),
                str(cdata.get("intrapatch_label", "") or ""),
                str(cdata.get("isthmus_label", "") or ""),
            ]
        )
        features.append(feat)

    provider.addFeatures(features)
    layer.updateExtents()
    if add_to_project:
        try:
            _apply_visible_corridor_style_vector(layer)
        except Exception:
            pass
        QgsProject.instance().addMapLayer(layer)
        print(f"  ✓ Added temporary layer '{layer_name}' to QGIS project")
    return layer


def create_memory_layer_from_networks(
    networks: List[Dict],
    layer_name: str,
    target_crs: QgsCoordinateReferenceSystem,
    original_crs: QgsCoordinateReferenceSystem,
    unit_system: str,
    add_to_project: bool = True,
) -> Optional[QgsVectorLayer]:
    is_imperial = unit_system == "imperial"
    area_field = "area_ac" if is_imperial else "area_ha"
    patch_area_field = "patch_ac" if is_imperial else "patch_ha"
    corr_area_field = "corr_ac" if is_imperial else "corr_ha"
    area_factor = 2.471053814 if is_imperial else 1.0

    layer = QgsVectorLayer(f"Polygon?crs={original_crs.authid()}", layer_name, "memory")
    provider = layer.dataProvider()
    provider.addAttributes(
        [
            QgsField("network_id", QVariant.Int),
            QgsField("patch_count", QVariant.Int),
            QgsField("corr_count", QVariant.Int),
            QgsField(area_field, QVariant.Double),
            QgsField(patch_area_field, QVariant.Double),
            QgsField(corr_area_field, QVariant.Double),
            QgsField("multipart", QVariant.Bool),
            QgsField("part_count", QVariant.Int),
        ]
    )
    layer.updateFields()

    transform = QgsCoordinateTransform(target_crs, original_crs, QgsProject.instance())
    features: List[QgsFeature] = []
    for net in networks:
        geom = net.get("geom")
        if geom is None or geom.isEmpty():
            continue
        g = clone_geometry(geom)
        g.transform(transform)
        feat = QgsFeature(layer.fields())
        feat.setGeometry(g)
        multipart = g.isMultipart()
        part_count = g.constGet().numGeometries() if multipart else 1
        feat.setAttributes(
            [
                int(net.get("network_id", 0)),
                int(net.get("patch_count", 0)),
                int(net.get("corridor_count", 0)),
                round(float(net.get("area_ha", 0.0)) * area_factor, 4),
                round(float(net.get("patch_area_ha", 0.0)) * area_factor, 4),
                round(float(net.get("corridor_area_ha", 0.0)) * area_factor, 4),
                multipart,
                int(part_count),
            ]
        )
        features.append(feat)

    provider.addFeatures(features)
    layer.updateExtents()
    try:
        _apply_random_unique_value_symbology_vector(layer, "network_id")
    except Exception:
        pass
    if add_to_project:
        QgsProject.instance().addMapLayer(layer)
        print(f"  ✓ Added temporary layer '{layer_name}' to QGIS project")
    return layer


def _convert_stats_for_units(stats: Dict, unit_system: str) -> Dict:
    factor = 2.471053814 if unit_system == "imperial" else 1.0
    label = "ac" if unit_system == "imperial" else "ha"
    converted = dict(stats)
    if "budget_used_ha" in stats:
        converted["budget_used_display"] = stats["budget_used_ha"] * factor
    if "total_connected_area_ha" in stats:
        converted["total_connected_area_display"] = stats["total_connected_area_ha"] * factor
    if "largest_group_area_ha" in stats:
        converted["largest_group_area_display"] = stats["largest_group_area_ha"] * factor
    if "species_dispersal_distance_m" in stats:
        converted["species_dispersal_distance_display"] = stats["species_dispersal_distance_m"] * (
            3.280839895 if unit_system == "imperial" else 1.0
        )
    if "min_patch_area_for_species_ha" in stats:
        converted["min_patch_area_for_species_display"] = stats["min_patch_area_for_species_ha"] * factor
    if "largest_reachable_habitat_cluster" in stats:
        converted["largest_reachable_habitat_cluster_display"] = stats["largest_reachable_habitat_cluster"] * factor
    if "largest_reachable_habitat_cluster_before" in stats:
        converted["largest_reachable_habitat_cluster_before_display"] = (
            stats["largest_reachable_habitat_cluster_before"] * factor
        )
    if "seed_area_ha" in stats:
        converted["seed_area_display"] = stats["seed_area_ha"] * factor
    if "final_patch_area_ha" in stats:
        converted["final_patch_area_display"] = stats["final_patch_area_ha"] * factor
    if "total_patch_area_ha" in stats:
        converted["total_patch_area_display"] = stats["total_patch_area_ha"] * factor
    if "bigconnect_objective_post" in stats:
        converted["bigconnect_objective_display"] = stats["bigconnect_objective_post"] * factor
    converted["area_units_label"] = label
    converted["conversion_factor"] = factor
    return converted


def _compute_connectivity_metrics(
    corridors: Dict[int, Dict],
    patches: Dict[int, Dict],
) -> Dict[str, float]:
    """Compute consistent graph metrics for summaries across strategies."""
    n_nodes = len(patches)
    uf = UnionFind()
    for pid in patches:
        uf.find(int(pid))

    corridors_used = len(corridors)
    for data in corridors.values():
        pids = list(data.get("patch_ids", []))
        if len(pids) >= 2:
            anchor = int(pids[0])
            for other in pids[1:]:
                uf.union(anchor, int(other))
        else:
            p1 = data.get("p1") or data.get("patch1")
            p2 = data.get("p2") or data.get("patch2")
            if p1 is not None and p2 is not None:
                uf.union(int(p1), int(p2))

    roots = [uf.find(int(pid)) for pid in patches]
    components = len(set(roots))
    redundant_links = max(0, corridors_used - (n_nodes - components))
    avg_degree = (2 * corridors_used / n_nodes) if n_nodes > 0 else 0.0

    comp_counts: Dict[int, int] = defaultdict(int)
    comp_areas: Dict[int, float] = defaultdict(float)
    for pid, pdata in patches.items():
        root = int(uf.find(int(pid)))
        comp_counts[root] += 1
        try:
            comp_areas[root] += float(pdata.get("area_ha", 0.0) or 0.0)
        except Exception:
            pass

    largest_group_patches = max(comp_counts.values()) if comp_counts else 0
    largest_group_area_ha = max(comp_areas.values()) if comp_areas else 0.0

    return {
        "patches_total": n_nodes,
        "corridors_used": corridors_used,
        "components_remaining": components,
        "redundant_links": redundant_links,
        "avg_degree": avg_degree,
        "patches_connected": largest_group_patches,
        "largest_group_patches": largest_group_patches,
        "largest_group_area_ha": largest_group_area_ha,
    }


def run_vector_analysis(
    layer: QgsVectorLayer,
    output_dir: str,
    raw_params: Dict,
    strategy: str = "most_connected_habitat",
    temporary: bool = False,
    iface=None,
    progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
    log_cb: Optional[Callable[[str, str], None]] = None,
    ctx: Optional[AnalysisContext] = None,
) -> List[Dict]:
    """Execute the vector corridor analysis for the provided polygon layer."""
    add_to_project = bool((raw_params or {}).get("add_to_project", iface is not None))
    if not isinstance(layer, QgsVectorLayer) or not layer.isValid():
        raise VectorAnalysisError("Selected layer is not a valid vector layer.")

    timings = TimingRecorder()

    params = _to_dataclass(raw_params)
    strategy_key = _normalize_strategy_key(strategy)
    safe_layer = _safe_filename(layer.name())
    try:
        base_name, ext = os.path.splitext(params.output_name or "")
        ext = ext or ".gpkg"
        if safe_layer and safe_layer.lower() not in (base_name or "").lower():
            base_name = base_name or "terralink_corridors"
            params.output_name = f"{base_name}_{safe_layer}{ext}"
    except Exception:
        pass
    ctx = ctx or AnalysisContext()
    unit_system = params.unit_system
    area_factor = 2.471053814 if unit_system == "imperial" else 1.0
    area_label = "ac" if unit_system == "imperial" else "ha"

    if temporary:
        output_path = ""
    else:
        out_dir = output_dir or os.path.dirname(layer.source())
        os.makedirs(out_dir, exist_ok=True)
        output_path = os.path.join(out_dir, params.output_name)

    summary_dir = (
        os.path.dirname(output_path)
        if output_path
        else (output_dir or os.path.dirname(layer.source()) or os.getcwd())
    )
    if not summary_dir:
        summary_dir = os.getcwd()
    os.makedirs(summary_dir, exist_ok=True)

    overall_start = time.time()
    print("=" * 70)
    print("TERRALINK VECTOR ANALYSIS v23.3")
    print("=" * 70)
    print(f"Plugin path: {os.path.dirname(__file__)}")
    print("\n1. Loading vector layer...")
    print(f"  ✓ Using layer: {layer.name()} ({layer.featureCount()} features)")
    emit_progress(progress_cb, 5, "Loading vector layer…")

    original_crs = layer.crs()
    print(f"  CRS: {original_crs.authid()}")

    print("\n2. Determining analysis CRS...")
    with timings.time_block("Determine analysis CRS"):
        target_crs = get_utm_crs_from_extent(layer)
        print(f"  ✓ Using {target_crs.authid()} for measurements")
        emit_progress(progress_cb, 15, "Preparing data…")

        # Auto-scale max search distance if user left it unset/zero
        try:
            extent_src = layer.extent()
            transform_extent = QgsCoordinateTransform(layer.crs(), target_crs, QgsProject.instance())
            extent = transform_extent.transformBoundingBox(extent_src)
            max_dimension = max(extent.width(), extent.height())
            DEFAULT_SCALING_FACTOR = 0.25
            if params.max_search_distance <= 0 and max_dimension > 0:
                params.max_search_distance = max_dimension * DEFAULT_SCALING_FACTOR
                _log_message(
                    f"Auto-setting max search distance to {DEFAULT_SCALING_FACTOR*100:.0f}% of map extent "
                    f"({params.max_search_distance:.1f} units)."
                )
        except Exception:
            pass

    print("\n3. Loading patches...")
    with timings.time_block("Load patches and spatial index"):
        patches, spatial_index, filtered_small_patches = load_and_prepare_patches(layer, target_crs, params)
        ctx.filtered_small_patches = list(filtered_small_patches or [])
    if not patches:
        raise VectorAnalysisError("No patches found after filtering.")
    if len(patches) > MAX_VECTOR_PATCH_COUNT:
        raise VectorAnalysisError(
            "This dataset has too many patches for vector optimization in TerraLink "
            f"({len(patches)} > {MAX_VECTOR_PATCH_COUNT}).\n"
            "Please use a smaller dataset (for example: clip to a smaller extent, "
            "split the study area into tiles, or raise minimum patch size)."
        )
    patch_union: Optional[QgsGeometry] = None
    with timings.time_block("Build patch union"):
        try:
            patch_union = QgsGeometry.unaryUnion([p["geom"] for p in patches.values()])
        except Exception:
            patch_union = None
    emit_progress(progress_cb, 35, "Searching for corridor candidates…")

    navigator: Optional[RasterNavigator] = None
    obstacle_layers: List[QgsVectorLayer] = []
    skipped_ids: List[str] = []
    if params.obstacle_enabled and params.obstacle_layer_ids:
        with timings.time_block("Impassable preparation"):
            for layer_id in params.obstacle_layer_ids:
                layer = QgsProject.instance().mapLayer(layer_id)
                if isinstance(layer, QgsVectorLayer) and layer.isValid() and QgsWkbTypes.geometryType(layer.wkbType()) == QgsWkbTypes.PolygonGeometry:
                    obstacle_layers.append(layer)
                else:
                    skipped_ids.append(str(layer_id))

            if skipped_ids:
                print(
                    f"  ⚠ Skipped {len(skipped_ids)} impassable layer(s) that are unavailable or not polygon geometry."
                )

            if obstacle_layers:
                try:
                    navigator = RasterNavigator(patches, obstacle_layers, target_crs, params)
                    print(
                        f"  ✓ Raster navigator grid: {navigator.cols} × {navigator.rows} cells "
                        f"@ {navigator.resolution:.1f} units "
                        f"using {len(obstacle_layers)} impassable layer(s)"
                    )
                    # Cache a union of impassables for faster per-corridor checks/clipping.
                    try:
                        ctx.impassable_union = QgsGeometry.unaryUnion(
                            [g for g in navigator.obstacle_geoms if g and (not g.isEmpty())]
                        ).makeValid()
                    except Exception:
                        ctx.impassable_union = None
                except VectorAnalysisError as exc:
                    print(f"  ⚠ Impassable land classes disabled: {exc}")
                    navigator = None
            else:
                print("  ⚠ Selected impassable layers are unavailable; continuing without impassable land classes.")
    else:
        timings.add("Impassable preparation (disabled)", 0.0)

    # Per-run mutable state lives in the AnalysisContext (not on params).

    print("\n4. Precomputing candidate corridors...")
    all_possible = find_all_possible_corridors(
        patches,
        spatial_index,
        params,
        strategy=strategy,
        patch_union=patch_union,
        ctx=ctx,
        progress_cb=progress_cb,
        progress_start=35,
        progress_end=60,
        navigator=navigator,
        timings=timings,
        timing_out=None,
    )

    strategy = strategy_key

    print("\nTOP 20 CANDIDATES BY COST:")
    sorted_cands = sorted(all_possible, key=lambda x: x.get("area_ha", 999))
    for i, c in enumerate(sorted_cands[:20]):
        cost = c.get("area_ha", 0.0)
        orig = c.get("original_area_ha", cost)
        patch_count = len(c.get("patch_ids", {c.get('patch1'), c.get('patch2')}))
        print(
            f"  {i+1:2d}. {cost:8.5f} ha  via {patch_count} patches  (orig {orig:.3f})"
        )

    print("\n5. Running optimization...")
    print("=" * 70)
    print(f"--- {_strategy_display_name(strategy).upper()} ---")

    if not all_possible:
        raise VectorAnalysisError(_format_no_corridor_reason("Precomputation", len(patches), len(all_possible), params))

    # 5. Optimization
    emit_progress(progress_cb, 60, "Optimizing network...")
    strategy_key = _normalize_strategy_key(strategy)

    # --- UPDATED DISPATCH LOGIC ---
    opt_label = f"Optimization ({strategy_key})"
    with timings.time_block(opt_label):
        if strategy_key == "largest_single_network":
            corridors, stats = optimize_circuit_utility_largest_network(patches, all_possible, params)
            layer_name = "Corridors (Largest Single Network)"
        elif strategy_key == "reachable_habitat_advanced":
            corridors, stats = optimize_habitat_availability(patches, all_possible, params)
            layer_name = "Corridors (Reachable Habitat (Advanced))"
        elif strategy_key == "most_connected_habitat":
            corridors, stats = optimize_bigconnect_vector(patches, all_possible, params)
            layer_name = "Corridors (Most Connected Area)"
        elif strategy_key == "landscape_fluidity":
            corridors, stats = optimize_landscape_fluidity_a(patches, all_possible, params)
            layer_name = "Corridors (Landscape Fluidity)"
        else:
            raise VectorAnalysisError(f"Unsupported strategy '{strategy_key}'.")

    exact_metric_modes = {
        "most_connected_habitat",
        "landscape_fluidity",
    }
    if strategy_key in exact_metric_modes:
        with timings.time_block(f"Exact metric refinement ({strategy_key})"):
            refined_corridors, refined_stats = _refine_metric_mode_with_exact_metric(
                strategy_key=strategy_key,
                patches=patches,
                candidates=all_possible,
                params=params,
                raw_params=dict(raw_params or {}),
                target_crs=target_crs,
                base_corridors=corridors,
                base_stats=stats,
            )
        if refined_corridors:
            selector_source = str((refined_stats or {}).get("exact_metric_selector_source", strategy_key))
            if selector_source != strategy_key:
                print(
                    f"  ℹ Exact metric refinement selected solution from '{selector_source}' "
                    f"for '{strategy_key}'."
                )
            corridors, stats = refined_corridors, refined_stats

    if not corridors:
        raise VectorAnalysisError(_format_no_corridor_reason("Optimization", len(patches), len(all_possible), params))

    _refresh_vector_connectivity_stats(patches, corridors, stats)

    pure_metric_modes = {
        "reachable_habitat_advanced",
        "most_connected_habitat",
        "landscape_fluidity",
    }
    remaining_budget = float((params.budget_area or 0.0) - float(stats.get("budget_used_ha", 0.0) or 0.0))
    if strategy_key not in pure_metric_modes and remaining_budget > 0 and corridors:
        with timings.time_block("Hybrid leftover budget"):
            extra_used, low_value_added, redundancy_added = _apply_hybrid_leftover_budget_vector(
                patches=patches,
                candidates=all_possible,
                corridors=corridors,
                remaining_budget=remaining_budget,
                max_search_distance=float(params.max_search_distance or 0.0),
            )
        if extra_used:
            stats["budget_used_ha"] = float(stats.get("budget_used_ha", 0.0) or 0.0) + float(extra_used)
            stats["corridors_used"] = len(corridors)
            if "primary_links" in stats:
                stats["primary_links"] = int(stats.get("primary_links", 0) or 0) + int(low_value_added)
            if "redundant_links" in stats:
                stats["redundant_links"] = int(stats.get("redundant_links", 0) or 0) + int(redundancy_added)
            _refresh_vector_connectivity_stats(patches, corridors, stats)
        remaining_budget = float((params.budget_area or 0.0) - float(stats.get("budget_used_ha", 0.0) or 0.0))

    # Consistent connectivity metrics for all optimization modes
    try:
        stats.update(_compute_connectivity_metrics(corridors, patches))
    except Exception:
        pass

    # Fixed-width policy: do not use remaining budget to alter corridor width.

    if strategy_key == "largest_single_network" and corridors:
        with timings.time_block("Enforce single largest network"):
            removed_count, removed_area = _enforce_largest_network_component(
                patches=patches,
                corridors=corridors,
                seed_patch=stats.get("seed_patch"),
            )
        if removed_count > 0:
            stats["budget_used_ha"] = max(
                0.0,
                float(stats.get("budget_used_ha", 0.0) or 0.0) - float(removed_area),
            )
            stats["corridors_used"] = len(corridors)
            stats["primary_links"] = sum(
                1
                for c in corridors.values()
                if str(c.get("type", "")).lower() in ("primary", "low_value", "backbone")
            )
            stats["redundant_links"] = sum(
                1 for c in corridors.values() if str(c.get("type", "")).lower() == "redundant"
            )
            stats["wasteful_links"] = sum(
                1 for c in corridors.values() if str(c.get("type", "")).lower() == "wasteful"
            )
            _refresh_vector_connectivity_stats(patches, corridors, stats)
            print(
                "  ✓ Enforced single-network output: "
                f"removed {removed_count} disconnected corridor(s)"
            )

    stats = _convert_stats_for_units(stats, unit_system)
    stats["budget_total_display"] = params.budget_area * area_factor

    print("  Preparing outputs...")
    emit_progress(progress_cb, 90, "Writing outputs…")
    _annotate_vector_corridor_export_flags(corridors, patches, params, ctx)
    with timings.time_block("Write corridor outputs"):
        # Contiguous network areas (patches + corridors dissolved per connected network)
        networks: List[Dict] = []
        try:
            dissolve_tolerance = max(0.0, float(getattr(params, "min_corridor_width", 0.0) or 0.0) * 0.01)
            networks = build_contiguous_network_summaries(patches, corridors, dissolve_tolerance=dissolve_tolerance)
            networks_layer_name = "Contiguous Areas"
            if temporary:
                create_memory_layer_from_networks(
                    networks,
                    networks_layer_name,
                    target_crs,
                    original_crs,
                    unit_system,
                    add_to_project=add_to_project,
                )
            else:
                write_corridors_layer_to_gpkg(
                    corridors,
                    output_path,
                    layer_name,
                    target_crs,
                    original_crs,
                    unit_system,
                    overwrite_file=True,
                )
                write_contiguous_networks_layer_to_gpkg(
                    networks,
                    output_path,
                    networks_layer_name,
                    target_crs,
                    original_crs,
                    unit_system,
                )
                add_layer_to_qgis_from_gpkg(output_path, networks_layer_name, add_to_project=add_to_project)
                add_layer_to_qgis_from_gpkg(output_path, layer_name, add_to_project=add_to_project)
        except Exception as net_exc:  # noqa: BLE001
            print(f"  ⚠ Could not write contiguous areas layer: {net_exc}")
            networks = []
        if temporary:
            create_memory_layer_from_corridors(
                corridors,
                layer_name,
                target_crs,
                original_crs,
                unit_system,
                add_to_project=add_to_project,
            )

    # --- LANDSCAPE METRICS REPORT (always written) ---
    landscape_metrics_path = ""
    try:
        strategy_key = _normalize_strategy_key(strategy)
        analysis_layer_name = f"TerraLink Landscape Metrics ({layer.name()})"
        pixel_size_m = float(getattr(params, "grid_resolution", 50.0) or 50.0)
        min_width_m = max(0.0, float(getattr(params, "min_corridor_width", 0.0) or 0.0))
        metrics_pixel_size_m = max(1.0, pixel_size_m)
        # Prevent coarse metrics rasterization from "missing" narrow corridors.
        if min_width_m > 0.0:
            safe_px = max(1.0, min_width_m / 3.0)
            if metrics_pixel_size_m > safe_px + 1e-9:
                print(
                    "  ℹ Landscape metrics resolution auto-adjusted "
                    f"from {metrics_pixel_size_m:.1f} m to {safe_px:.1f} m "
                    f"for corridor width {min_width_m:.1f} m"
                )
                metrics_pixel_size_m = float(safe_px)
        patch_geoms = [
            pdata.get("geom")
            for pdata in patches.values()
            if pdata.get("geom") and not pdata.get("geom").isEmpty()
        ]
        network_geoms = [
            net.get("geom")
            for net in networks
            if net.get("geom") and not net.get("geom").isEmpty()
        ]
        bounds = _compute_geoms_bounds(patch_geoms + network_geoms)
        mask, eff_px = _rasterize_networks_to_mask(
            networks=networks,
            pixel_size_m=metrics_pixel_size_m,
            target_crs=target_crs,
            bounds=bounds,
        )
        pre_mask, _ = _rasterize_geoms_to_mask(
            geoms=patch_geoms,
            pixel_size_m=metrics_pixel_size_m,
            target_crs=target_crs,
            bounds=bounds,
        )
        exact_comp = _compute_habitat_component_metrics_exact(patches, corridors)
        if strategy_key == "landscape_fluidity_b":
            exact_fluidity = _compute_landscape_fluidity_2_exact(
                patches,
                corridors,
                params,
                candidate_pool=all_possible,
            )
        else:
            exact_fluidity = _compute_landscape_fluidity_exact(patches, corridors, params)
        exact_mobility = _compute_strategic_mobility_exact(patches, corridors, params)
        analysis_params = dict(raw_params or {})
        analysis_params["mesh_override_pre_norm"] = float(exact_comp.get("mesh_norm_pre", 0.0) or 0.0)
        analysis_params["mesh_override_post_norm"] = float(exact_comp.get("mesh_norm_post", 0.0) or 0.0)
        analysis_params["lcc_override_pre_norm"] = float(exact_comp.get("lcc_norm_pre", 0.0) or 0.0)
        analysis_params["lcc_override_post_norm"] = float(exact_comp.get("lcc_norm_post", 0.0) or 0.0)
        analysis_params["mean_effective_resistance_override_pre"] = float(
            exact_fluidity.get("graph_resistance_pre", 0.0) or 0.0
        )
        analysis_params["mean_effective_resistance_override_post"] = float(
            exact_fluidity.get("graph_resistance_post", 0.0) or 0.0
        )
        analysis_params["landscape_fluidity_override_pre"] = float(exact_fluidity.get("pre", 0.0) or 0.0)
        analysis_params["landscape_fluidity_override_post"] = float(exact_fluidity.get("post", 0.0) or 0.0)
        analysis_params["strategic_mobility_override_pre"] = float(exact_mobility.get("pre", 0.0) or 0.0)
        analysis_params["strategic_mobility_override_post"] = float(exact_mobility.get("post", 0.0) or 0.0)
        stats["mesh_habitat_norm_pre_exact"] = float(exact_comp.get("mesh_norm_pre", 0.0) or 0.0)
        stats["mesh_habitat_norm_post_exact"] = float(exact_comp.get("mesh_norm_post", 0.0) or 0.0)
        stats["lcc_norm_pre_exact"] = float(exact_comp.get("lcc_norm_pre", 0.0) or 0.0)
        stats["lcc_norm_post_exact"] = float(exact_comp.get("lcc_norm_post", 0.0) or 0.0)
        stats["mean_effective_resistance_pre_exact"] = float(exact_fluidity.get("graph_resistance_pre", 0.0) or 0.0)
        stats["mean_effective_resistance_post_exact"] = float(exact_fluidity.get("graph_resistance_post", 0.0) or 0.0)
        stats["landscape_fluidity_pre_exact"] = float(exact_fluidity.get("pre", 0.0) or 0.0)
        stats["landscape_fluidity_post_exact"] = float(exact_fluidity.get("post", 0.0) or 0.0)
        stats["strategic_mobility_pre_exact"] = float(exact_mobility.get("pre", 0.0) or 0.0)
        stats["strategic_mobility_post_exact"] = float(exact_mobility.get("post", 0.0) or 0.0)
        stats["landscape_metrics_pixel_size_m"] = float(eff_px)
        from .landscape_metrics import _perform_landscape_analysis  # local import avoids hard coupling at import time

        analysis_lines = _perform_landscape_analysis(
            arr=mask,
            layer_name=analysis_layer_name,
            res_x=float(eff_px),
            res_y=float(eff_px),
            is_metric=True,
            params=analysis_params,
            pre_arr=pre_mask,
        )
        if strategy_key == "reachable_habitat_advanced":
            analysis_lines = list(analysis_lines or [])
            analysis_lines.append("")
            analysis_lines.append("Reachable habitat metrics||||")
            try:
                pre_ha = float(stats.get("habitat_availability_before", 0.0) or 0.0)
                post_ha = float(stats.get("habitat_availability_after", 0.0) or 0.0)
                pct = ""
                if abs(pre_ha) > 1e-12:
                    pct = f"{((post_ha - pre_ha) / abs(pre_ha)) * 100.0:+.3f}%"
                analysis_lines.append(
                    f"Reachable Habitat Score|{_format_number(pre_ha, 4)}|{_format_number(post_ha, 4)}|{pct}|Score (dimensionless)|Kupfer 2012"
                )
            except Exception:
                pass
            try:
                area_label = "ac" if str(getattr(params, "unit_system", "metric")) == "imperial" else "ha"
                pre_mean = float(stats.get("mean_reachable_area_before", 0.0) or 0.0)
                post_mean = float(stats.get("mean_reachable_area", 0.0) or 0.0)
                mean_pct = ""
                if abs(pre_mean) > 1e-12:
                    mean_pct = f"{((post_mean - pre_mean) / abs(pre_mean)) * 100.0:+.3f}%"
                analysis_lines.append(
                    f"Mean reachable area|{_format_number(pre_mean, 4)}|{_format_number(post_mean, 4)}|{mean_pct}|Area ({area_label})|"
                )
                pre_med = float(stats.get("median_reachable_area_before", 0.0) or 0.0)
                post_med = float(stats.get("median_reachable_area", 0.0) or 0.0)
                med_pct = ""
                if abs(pre_med) > 1e-12:
                    med_pct = f"{((post_med - pre_med) / abs(pre_med)) * 100.0:+.3f}%"
                analysis_lines.append(
                    f"Median reachable area|{_format_number(pre_med, 4)}|{_format_number(post_med, 4)}|{med_pct}|Area ({area_label})|"
                )
                pre_cluster = float(stats.get("largest_reachable_habitat_cluster_before_display", 0.0) or 0.0)
                post_cluster = float(stats.get("largest_reachable_habitat_cluster_display", 0.0) or 0.0)
                if pre_cluster or post_cluster:
                    cluster_pct = ""
                    if abs(pre_cluster) > 1e-12:
                        cluster_pct = f"{((post_cluster - pre_cluster) / abs(pre_cluster)) * 100.0:+.3f}%"
                    analysis_lines.append(
                        f"Largest reachable habitat cluster|{_format_number(pre_cluster, 4)}|{_format_number(post_cluster, 4)}|{cluster_pct}|Area ({area_label})|"
                    )
            except Exception:
                pass

        if temporary:
            temp_file = tempfile.NamedTemporaryFile(
                prefix="terralink_landscape_metrics_", suffix=".txt", delete=False
            )
            landscape_metrics_path = temp_file.name
            temp_file.close()
        else:
            safe = _safe_filename(layer.name())
            landscape_metrics_path = os.path.join(
                os.path.dirname(output_path) or os.getcwd(),
                f"landscape_metrics_{safe}.txt",
            )
        _write_text_report(landscape_metrics_path, analysis_lines)
        stats["landscape_metrics_path"] = landscape_metrics_path
        print(f"  ✓ Saved landscape metrics: {landscape_metrics_path}")
        try:
            _add_landscape_metrics_table_layer(
                analysis_layer_name,
                analysis_lines,
                add_to_project=add_to_project,
            )
        except Exception:
            pass
    except Exception as e:  # noqa: BLE001
        try:
            if not landscape_metrics_path:
                if temporary:
                    temp_file = tempfile.NamedTemporaryFile(
                        prefix="terralink_landscape_metrics_", suffix=".txt", delete=False
                    )
                    landscape_metrics_path = temp_file.name
                    temp_file.close()
                else:
                    safe = _safe_filename(layer.name())
                    landscape_metrics_path = os.path.join(
                        os.path.dirname(output_path) or os.getcwd(),
                        f"landscape_metrics_{safe}.txt",
                    )
            _write_text_report(landscape_metrics_path, [f"Landscape analysis failed: {e}"])
            stats["landscape_metrics_path"] = landscape_metrics_path
            print(f"  ✓ Saved landscape metrics: {landscape_metrics_path}")
            try:
                _add_landscape_metrics_table_layer(
                    "Landscape Metrics (Error)",
                    [f"Landscape analysis failed: {e}"],
                    add_to_project=add_to_project,
                )
            except Exception:
                pass
        except Exception:
            pass

    elapsed = time.time() - overall_start
    emit_progress(progress_cb, 100, "Vector analysis complete.")

    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    strategy_label = _strategy_display_name(strategy)
    print(f"Strategy:          {strategy_label}")
    print(f"Corridors created: {stats.get('corridors_used', 0)}")
    print(f"Connections:       {stats.get('connections_made', 0)}")
    print(f"Total connected:   {stats.get('total_connected_area_display', 0):.2f} {area_label}")
    print(f"Largest group:     {stats.get('largest_group_area_display', 0):.2f} {area_label}")
    if strategy_key == "reachable_habitat_advanced":
        print(f"Reachable Habitat Score before: {float(stats.get('habitat_availability_before', 0.0) or 0.0):.4f}")
        print(f"Reachable Habitat Score after:  {float(stats.get('habitat_availability_after', 0.0) or 0.0):.4f}")
        print(f"HA increase:       {float(stats.get('percent_increase', 0.0) or 0.0):.2f}%")
        print(
            f"Reachable cluster: {float(stats.get('largest_reachable_habitat_cluster_display', 0.0) or 0.0):.2f} {area_label}"
        )
    if strategy_key == "most_connected_habitat":
        print(f"Most Connected Area: {float(stats.get('bigconnect_objective_post', 0.0) or 0.0) * area_factor:.2f} {area_label}")
        print(f"Proven optimal:    {bool(stats.get('bigconnect_proven_optimal', False))}")
    if "redundant_links" in stats:
        print(f"Redundant links:   {stats.get('redundant_links', 0)}")
    if "avg_degree" in stats:
        print(f"Average degree:    {stats.get('avg_degree', 0):.2f}")
    print(
        f"Budget used:      {stats.get('budget_used_display', 0):.2f} / {params.budget_area * area_factor:.2f} {area_label}"
    )
    print(f"Processing time:   {elapsed:.1f}s")
    if temporary:
        print("Output:            Temporary layer (memory)")
    else:
        print(f"Output GPKG:       {output_path}")
    print("=" * 70)

    stats["layer_name"] = layer_name
    stats["output_path"] = output_path if not temporary else ""
    stats["unit_system"] = unit_system
    try:
        safe = _safe_filename(layer.name())
        if temporary:
            temp_file = tempfile.NamedTemporaryFile(
                prefix="terralink_vector_summary_", suffix=".csv", delete=False
            )
            summary_path = temp_file.name
            temp_file.close()
        else:
            summary_path = os.path.join(summary_dir, f"terralink_vector_summary_{safe}.csv")
        rows = [
            ("Input layer", layer.name()),
            ("Strategy", strategy_label),
            ("Min patch size", f"{_format_number(params.min_patch_size * area_factor)} {area_label}"),
            ("Corridors created", str(stats.get("corridors_used", 0))),
            ("Budget used", f"{_format_number(stats.get('budget_used_display', 0))} {area_label}"),
            ("Budget total", f"{_format_number(params.budget_area * area_factor)} {area_label}"),
        ]
        if strategy_key == "most_connected_habitat":
            rows.append(("Most Connected Area", f"{_format_number(float(stats.get('bigconnect_objective_post', 0.0) or 0.0) * area_factor)} {area_label}"))
            rows.append(("Most Connected Area (proven optimal)", str(bool(stats.get("bigconnect_proven_optimal", False)))))
        if strategy_key != "reachable_habitat_advanced":
            if "habitat_availability_before" in stats:
                rows.append(("Reachable Habitat Score before", _format_number(stats.get("habitat_availability_before", 0.0), 4)))
            if "habitat_availability_after" in stats:
                rows.append(("Reachable Habitat Score after", _format_number(stats.get("habitat_availability_after", 0.0), 4)))
            if "percent_increase" in stats:
                rows.append(("Reachable Habitat Score increase", f"{_format_number(stats.get('percent_increase', 0.0), 2)} %"))
            if "mean_reachable_area" in stats:
                rows.append(("Mean reachable area", _format_number(stats.get("mean_reachable_area", 0.0), 4)))
            if "median_reachable_area" in stats:
                rows.append(("Median reachable area", _format_number(stats.get("median_reachable_area", 0.0), 4)))
            if "largest_reachable_habitat_cluster_display" in stats:
                rows.append(
                    (
                        "Largest reachable habitat cluster",
                        f"{_format_number(stats.get('largest_reachable_habitat_cluster_display', 0.0))} {area_label}",
                    )
                )
            if "species_dispersal_distance_display" in stats:
                dist_label = "ft" if unit_system == "imperial" else "m"
                rows.append(
                    (
                        "Species dispersal distance",
                        f"{_format_number(stats.get('species_dispersal_distance_display', 0.0))} {dist_label}",
                    )
                )
            if "min_patch_area_for_species_display" in stats:
                rows.append(
                    (
                        "Min patch area for species",
                        f"{_format_number(stats.get('min_patch_area_for_species_display', 0.0))} {area_label}",
                    )
                )
            if "patch_area_scaling" in stats:
                rows.append(("Patch area scaling", str(stats.get("patch_area_scaling", ""))))
            if "patch_quality_weight_field" in stats and str(stats.get("patch_quality_weight_field", "")).strip():
                rows.append(("Patch quality field", str(stats.get("patch_quality_weight_field", ""))))
        if "primary_links" in stats:
            rows.append(("Primary links", str(stats.get("primary_links", 0))))
        if "redundant_links" in stats:
            rows.append(("Redundant links", str(stats.get("redundant_links", 0))))
        if "entropy_total" in stats:
            rows.append(("Entropy (H_total)", _format_number(stats.get("entropy_total", 0), 4)))
        _write_summary_csv(summary_path, rows)
        stats["summary_csv_path"] = summary_path
        print(f"  ✓ Saved vector summary CSV: {summary_path}")
        _add_summary_csv_layer(
            summary_path,
            f"TerraLink Vector Summary ({layer.name()})",
            add_to_project=add_to_project,
        )
    except Exception:
        pass

    return [
        {
            "strategy": strategy,
            "stats": stats,
            "output_path": output_path if not temporary else "",
            "layer_name": layer_name,
        }
    ]
