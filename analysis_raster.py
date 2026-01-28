"""
TerraLink Corridor Analysis - Raster Workflow (v23.0)
-----------------------------------------------------
Runs the selected raster optimization workflow for corridor analysis.

The logic is adapted from the standalone raster script and packaged so it can
be invoked from the QGIS plugin with user-supplied parameters.
"""

from __future__ import annotations

import heapq
import math
import os
import random
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
# NumPy 2.x removed np.int; add shim for legacy references.
if not hasattr(np, "int"):  # pragma: no cover
    np.int = int  # type: ignore[attr-defined]
try:
    from scipy import ndimage

    HAS_NDIMAGE = True
except Exception:  # pragma: no cover
    ndimage = None  # type: ignore
    HAS_NDIMAGE = False
try:
    from osgeo import gdal
except ImportError:  # pragma: no cover
    gdal = None  # type: ignore

try:
    from qgis.core import (
        Qgis,
        QgsApplication,
        QgsFeature,
        QgsField,
        QgsPalettedRasterRenderer,
        QgsProject,
        QgsRasterLayer,
        QgsUnitTypes,
        QgsVectorLayer,
    )
except ImportError:  # pragma: no cover
    Qgis = None  # type: ignore
    QgsApplication = None  # type: ignore
    QgsFeature = None  # type: ignore
    QgsField = None  # type: ignore
    QgsPalettedRasterRenderer = None  # type: ignore
    QgsProject = None  # type: ignore
    QgsRasterLayer = None  # type: ignore
    QgsVectorLayer = None  # type: ignore

from PyQt5.QtCore import QVariant
from PyQt5.QtGui import QColor

# Optional OpenCV import for faster connected component labeling
try:
    import cv2

    HAS_CV2 = True
except ImportError:
    HAS_CV2 = False

# Import graph-metrics helper library for entropy/robustness calculations
try:
    from . import graph_math
    import networkx as nx

    HAS_GRAPH_MATH = True
except ImportError:
    graph_math = None  # type: ignore
    nx = None  # type: ignore
    HAS_GRAPH_MATH = False

GTIFF_OPTIONS = ["COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"]

PIXEL_COUNT_WARNING_THRESHOLD = 40_000_000  # warn when raster exceeds ~40 million pixels
PIXEL_SIZE_WARNING_THRESHOLD = 10.0  # warn when pixel size < 10 map units
PIXEL_FINE_CRITICAL_COUNT = 2_000_000  # only block fine rasters when they are also large
ROI_REDUNDANCY_BIAS = 0.2  # prefer new connections unless redundancy ROI is clearly higher
# Heuristic work-unit cap to keep Circuit Utility runs under ~3 minutes on fast hardware.
MAX_CANDIDATE_SEARCH_WORK = 6000.0

def _log_message(message: str, level: int = Qgis.Info) -> None:
    """Log to the QGIS Log Messages Panel with a TerraLink tag."""
    try:
        QgsApplication.messageLog().logMessage(message, "TerraLink", level)
    except Exception:
        print(f"TerraLink Log: {message}")


def _safe_filename(name: str, max_len: int = 64) -> str:
    safe = "".join(ch if (ch.isalnum() or ch in ("-", "_")) else "_" for ch in (name or ""))
    safe = safe.strip("_") or "layer"
    return safe[:max_len]


def _apply_random_unique_value_symbology(layer: "QgsRasterLayer", values: List[int]) -> None:
    """
    Apply a paletted renderer with stable pseudo-random colors for integer component ids.
    """
    if layer is None or not getattr(layer, "isValid", lambda: False)():
        return
    if QgsPalettedRasterRenderer is None:
        return
    try:
        vals = [int(v) for v in values if int(v) > 0]
    except Exception:
        return
    if not vals:
        return
    max_classes = 2048
    if len(vals) > max_classes:
        vals = vals[:max_classes]
        try:
            _log_message(f"[Symbology] Too many components; using first {max_classes} classes.")
        except Exception:
            pass

    classes = []
    try:
        # Make background transparent.
        classes.append(QgsPalettedRasterRenderer.Class(0, QColor(0, 0, 0, 0), "NoData"))
    except Exception:
        pass
    for v in vals:
        # Stable color based on id (so reruns look consistent).
        hue = int((v * 137) % 360)
        col = QColor.fromHsv(hue, 200, 255, 255)
        classes.append(QgsPalettedRasterRenderer.Class(v, col, str(v)))

    try:
        renderer = QgsPalettedRasterRenderer(layer.dataProvider(), 1, classes)
        layer.setRenderer(renderer)
        layer.triggerRepaint()
    except Exception:
        return


def _apply_random_unique_value_symbology_from_array(layer: "QgsRasterLayer", arr: np.ndarray) -> None:
    """Apply random unique-value colors for the non-zero unique values present in `arr`."""
    try:
        if arr is None or arr.size == 0:
            return
        vals = np.unique(arr)
        vals = [int(v) for v in vals.tolist() if int(v) > 0]
        _apply_random_unique_value_symbology(layer, vals)
    except Exception:
        return


def _write_text_report(path: str, lines: List[str]) -> None:
    try:
        os.makedirs(os.path.dirname(path) or os.getcwd(), exist_ok=True)
    except Exception:
        pass
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines).rstrip() + "\n")


def _add_landscape_metrics_table_layer(layer_title: str, analysis_lines: List[str]) -> None:
    """
    Add a non-spatial table layer to the QGIS project with the landscape metrics.

    This is a UI-friendly alternative to opening text files or showing popups.
    """
    if QgsVectorLayer is None or QgsProject is None or QgsFeature is None or QgsField is None:
        return
    try:
        rows: List[Tuple[str, str, str, str]] = []
        in_table = False
        for line in analysis_lines or []:
            if (not in_table) and ("METRIC NAME" in line) and ("|" in line):
                in_table = True
                continue
            if not in_table:
                continue
            s = (line or "").strip()
            if not s:
                continue
            if s.startswith("="):
                break
            if s.startswith("-"):
                continue
            if "|" not in s:
                continue
            parts = [p.strip() for p in s.split("|")]
            if len(parts) < 3:
                continue
            if len(parts) >= 4:
                metric, pre_val, post_val, interp = parts[0], parts[1], parts[2], parts[3]
            else:
                metric, pre_val, post_val, interp = parts[0], parts[1], "", parts[2]
            if metric:
                rows.append((metric, pre_val, post_val, interp))

        if not rows:
            joined = " ".join([l.strip() for l in (analysis_lines or []) if l.strip()])[:240]
            if not joined:
                joined = "No landscape metrics available."
            rows = [("Message", joined, "", "")]

        uri = "None?field=metric:string(80)&field=pre_value:string(40)&field=post_value:string(40)&field=interpretation:string(120)"
        layer = QgsVectorLayer(uri, layer_title, "memory")
        if not layer.isValid():
            return
        provider = layer.dataProvider()
        feats: List[QgsFeature] = []
        for metric, pre_val, post_val, interp in rows:
            feat = QgsFeature(layer.fields())
            feat.setAttributes([metric, pre_val, post_val, interp])
            feats.append(feat)
        provider.addFeatures(feats)
        layer.updateExtents()
        QgsProject.instance().addMapLayer(layer)
    except Exception:
        return


def _emit_progress(
    progress_cb: Optional[Callable[[int, Optional[str]], None]],
    value: float,
    message: Optional[str] = None,
) -> None:
    if progress_cb is None:
        return
    try:
        progress_cb(int(max(0, min(100, value))), message)
    except Exception:
        pass


def _map_units_to_meters(units: int) -> Optional[float]:
    if units == QgsUnitTypes.DistanceMeters:
        return 1.0
    if units == QgsUnitTypes.DistanceFeet:
        return 0.3048
    if units == QgsUnitTypes.DistanceFeetUS:
        return 0.3048006096
    return None


def _apply_raster_display_units(
    stats: Dict,
    *,
    layer: QgsRasterLayer,
    unit_system: str,
    pixel_w: float,
    pixel_h: float,
) -> None:
    stats["raster_units"] = unit_system
    if unit_system == "pixels":
        stats["raster_area_label"] = "px"
        stats["raster_dist_label"] = "px"
        stats["min_patch_size_display"] = stats.get("min_patch_size")
        stats["min_corridor_width_display"] = stats.get("min_corridor_width")
        stats["max_search_distance_display"] = stats.get("max_search_distance")
        stats["budget_used_display"] = stats.get("budget_used")
        stats["budget_total_display"] = stats.get("budget_total")
        return

    units_to_m = _map_units_to_meters(layer.crs().mapUnits())
    if units_to_m is None:
        return

    pixel_area_m2 = abs(pixel_w) * abs(pixel_h) * units_to_m * units_to_m
    pixel_size_m = max(abs(pixel_w), abs(pixel_h)) * units_to_m
    if pixel_area_m2 <= 0 or pixel_size_m <= 0:
        return

    if unit_system == "metric":
        area_label = "ha"
        dist_label = "m"
        area_factor = 10000.0
        dist_factor = 1.0
    else:
        area_label = "ac"
        dist_label = "ft"
        area_factor = 4046.8564224
        dist_factor = 0.3048

    stats["raster_area_label"] = area_label
    stats["raster_dist_label"] = dist_label
    try:
        stats["min_patch_size_display"] = float(stats.get("min_patch_size", 0)) * pixel_area_m2 / area_factor
        stats["min_corridor_width_display"] = float(stats.get("min_corridor_width", 0)) * (
            pixel_size_m / dist_factor
        )
        stats["max_search_distance_display"] = float(stats.get("max_search_distance", 0)) * (
            pixel_size_m / dist_factor
        )
        stats["budget_used_display"] = float(stats.get("budget_used", 0)) * pixel_area_m2 / area_factor
        stats["budget_total_display"] = float(stats.get("budget_total", 0)) * pixel_area_m2 / area_factor
    except Exception:
        return


class RasterAnalysisError(RuntimeError):
    """Raised when the raster analysis cannot be completed."""


@dataclass
class RasterRunParams:
    patch_connectivity: int
    patch_mode: str
    patch_values: List[float]
    range_lower: Optional[float]
    range_upper: Optional[float]
    obstacle_enabled: bool
    obstacle_mode: str
    obstacle_values: List[float]
    obstacle_range_lower: Optional[float]
    obstacle_range_upper: Optional[float]
    value_tolerance: float
    nodata_fallback: float
    min_patch_size: int
    budget_pixels: int
    max_search_distance: int
    max_corridor_area: Optional[int]
    min_corridor_width: int
    allow_bottlenecks: bool


class UnionFind:
    """Union-Find data structure for tracking connected components."""

    def __init__(self):
        self.parent: Dict[int, int] = {}
        self.size: Dict[int, int] = {}
        self.count: Dict[int, int] = {}

    def find(self, x: int) -> int:
        if x not in self.parent:
            self.parent[x] = x
            self.size[x] = 0
            self.count[x] = 0
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: int, b: int) -> int:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return ra
        sa, sb = self.size.get(ra, 0), self.size.get(rb, 0)
        if sa < sb:
            ra, rb = rb, ra
            sa, sb = sb, sa
        self.parent[rb] = ra
        self.size[ra] = sa + sb
        self.count[ra] = self.count.get(ra, 0) + self.count.get(rb, 0)
        return ra

    def get_size(self, x: int) -> int:
        return self.size.get(self.find(x), 0)

    def get_count(self, x: int) -> int:
        return self.count.get(self.find(x), 0)


def label_components_numpy(binary_array: np.ndarray, connectivity: int = 8) -> Tuple[np.ndarray, int]:
    """Label connected components using numpy (no external dependencies)."""
    rows, cols = binary_array.shape
    uf = UnionFind()
    if connectivity == 4:
        neighbors = [(-1, 0), (0, -1)]
    else:
        neighbors = [(-1, -1), (-1, 0), (-1, 1), (0, -1)]

    for i in range(rows):
        for j in range(cols):
            if binary_array[i, j]:
                cur = (i, j)
                uf.find(cur)
                for di, dj in neighbors:
                    ni, nj = i + di, j + dj
                    if 0 <= ni < rows and 0 <= nj < cols and binary_array[ni, nj]:
                        uf.union(cur, (ni, nj))

    root_to_label: Dict[Tuple[int, int], int] = {}
    next_label = 1
    labeled = np.zeros_like(binary_array, dtype=np.int32)
    for i in range(rows):
        for j in range(cols):
            if binary_array[i, j]:
                root = uf.find((i, j))
                if root not in root_to_label:
                    root_to_label[root] = next_label
                    next_label += 1
                labeled[i, j] = root_to_label[root]

    return labeled, next_label - 1


def label_components_opencv(binary_array: np.ndarray, connectivity: int = 8) -> Tuple[np.ndarray, int]:
    """Label connected components using OpenCV when available."""
    cv_conn = 4 if connectivity == 4 else 8
    n_labels, labeled = cv2.connectedComponents(binary_array.astype(np.uint8), connectivity=cv_conn)
    return labeled.astype(np.int32), n_labels - 1


def label_patches(binary_array: np.ndarray, connectivity: int = 8) -> Tuple[np.ndarray, int]:
    """Label connected components using the fastest available approach."""
    if HAS_CV2:
        print("  Using OpenCV for labeling...")
        return label_components_opencv(binary_array, connectivity)
    print("  Using numpy for labeling...")
    return label_components_numpy(binary_array, connectivity)


def read_band(
    band: gdal.Band,
    rows: int,
    cols: int,
    progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
    progress_start: int = 0,
    progress_end: int = 10,
) -> np.ndarray:
    """Read a raster band as a numpy array with incremental progress updates."""
    data = np.empty((rows, cols), dtype=np.float32)
    chunk_rows = max(1, min(1024, rows // 50 or 1))
    span = max(progress_end - progress_start, 1)

    for start_row in range(0, rows, chunk_rows):
        this_rows = min(chunk_rows, rows - start_row)
        buf = band.ReadRaster(0, start_row, cols, this_rows, cols, this_rows, gdal.GDT_Float32)
        if not buf:
            raise RasterAnalysisError("Failed to read raster data chunk.")
        arr = np.frombuffer(buf, dtype=np.float32, count=cols * this_rows).reshape(this_rows, cols)
        data[start_row : start_row + this_rows] = arr

        if progress_cb is not None:
            ratio = (start_row + this_rows) / max(rows, 1)
            progress_value = progress_start + ratio * span
            _emit_progress(progress_cb, progress_value, "Reading raster data…")

    return data


def write_raster(path: str, arr: np.ndarray, gt: Tuple[float, ...], proj: str, nodata: float = 0) -> None:
    """Write a numpy array out to GeoTIFF."""
    rows, cols = arr.shape
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, cols, rows, 1, gdal.GDT_Int32, options=GTIFF_OPTIONS)
    if ds is None:
        raise RasterAnalysisError(f"Unable to create output dataset: {path}")
    ds.SetGeoTransform(gt)
    ds.SetProjection(proj)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(int(nodata))

    if arr.dtype != np.int32:
        arr = arr.astype(np.int32)

    buf = np.ascontiguousarray(arr).tobytes()
    band.WriteRaster(0, 0, cols, rows, buf, cols, rows)
    band.FlushCache()
    ds = None


def define_habitat(data: np.ndarray, nodata_mask: np.ndarray, params: RasterRunParams) -> np.ndarray:
    """Identify patch pixels based on the selected configuration."""
    valid = ~nodata_mask
    patches = np.zeros(data.shape, dtype=np.uint8)
    mode = params.patch_mode.lower()
    tol = params.value_tolerance

    if mode == "value" and params.patch_values:
        for val in params.patch_values:
            patches |= (np.abs(data - val) < tol) & valid
        print(f"  Patch = values {params.patch_values}")
    elif mode == "range" and params.range_lower is not None and params.range_upper is not None:
        patches = ((data >= params.range_lower) & (data <= params.range_upper)) & valid
        print(
            "  Patch = range "
            f"{params.range_lower:.4f} - {params.range_upper:.4f}"
        )
    else:
        raise RasterAnalysisError("Patch configuration did not yield any valid pixels.")

    return patches


def define_obstacles(data: np.ndarray, nodata_mask: np.ndarray, patch_mask: np.ndarray, params: RasterRunParams) -> np.ndarray:
    """Create a boolean mask for impassable pixels corridors must avoid."""
    if not params.obstacle_enabled:
        return np.zeros(data.shape, dtype=bool)

    mask = np.zeros(data.shape, dtype=bool)
    tol = params.value_tolerance
    mode = params.obstacle_mode.lower()

    if mode == "range" and params.obstacle_range_lower is not None and params.obstacle_range_upper is not None:
        lower = min(params.obstacle_range_lower, params.obstacle_range_upper)
        upper = max(params.obstacle_range_lower, params.obstacle_range_upper)
        mask = (data >= lower) & (data <= upper)
    elif mode == "value" and params.obstacle_values:
        for val in params.obstacle_values:
            mask |= np.abs(data - val) < tol
    else:
        return np.zeros(data.shape, dtype=bool)

    mask &= ~nodata_mask
    return mask


def find_shortest_corridor(
    start_patches: Set[int],
    labels: np.ndarray,
    habitat: np.ndarray,
    max_width: int,
    connectivity: int,
    obstacle_mask: Optional[np.ndarray] = None,
    passable_mask: Optional[np.ndarray] = None,
    params: Optional[RasterRunParams] = None,
    target_patch_id: Optional[int] = None,
    target_boundary_center: Optional[Tuple[int, int]] = None,
    target_boundary_radius: Optional[int] = None,
    ghost_pixels: Optional[Set[Tuple[int, int]]] = None,
    ghost_factor: float = 10.0,
    start_positions_override: Optional[List[Tuple[int, int]]] = None,
) -> List[Tuple[frozenset, frozenset, int, float]]:
    """
    Dijkstra search to find shortest corridors connecting start_patches to other patches.
    Returns a list of (path_pixels, habitat_pixels, target_patch, length).

    `habitat_pixels` records any habitat cells (labels==0 but habitat==True) traversed along
    the chosen path. These are not "corridor cost" pixels, but are used to ensure the
    contiguous output raster includes connected habitat (e.g., small patches used as bridges).
    """
    rows, cols = labels.shape
    SQRT2 = 1.4142135623730951  # Pre-computed math.sqrt(2)
    if connectivity == 4:
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    else:
        moves = [
            (-1, -1),
            (-1, 0),
            (-1, 1),
            (0, -1),
            (0, 1),
            (1, -1),
            (1, 0),
            (1, 1),
        ]

    # Convert ghost_pixels set to numpy mask for faster lookup
    ghost_mask: Optional[np.ndarray] = None
    if ghost_pixels:
        ghost_mask = np.zeros((rows, cols), dtype=bool)
        for gr, gc in ghost_pixels:
            if 0 <= gr < rows and 0 <= gc < cols:
                ghost_mask[gr, gc] = True

    start_positions: List[Tuple[int, int]] = start_positions_override or []
    if not start_positions_override:
        for i in range(rows):
            for j in range(cols):
                if habitat[i, j]:
                    continue
                if obstacle_mask is not None and obstacle_mask[i, j]:
                    continue
                if passable_mask is not None and not passable_mask[i, j]:
                    continue
                for di, dj in moves:
                    ni, nj = i + di, j + dj
                    if 0 <= ni < rows and 0 <= nj < cols and labels[ni, nj] in start_patches:
                        start_positions.append((i, j))
                        break

    if not start_positions:
        return []

    heap: List[Tuple[float, int, int]] = []
    best_cost: Dict[Tuple[int, int], float] = {}
    parents: Dict[Tuple[int, int], Tuple[Optional[Tuple[int, int]], Optional[str]]] = {}
    results: List[Tuple[frozenset, frozenset, int, float]] = []
    visited_targets: Set[int] = set()

    def _reconstruct(state: Tuple[int, int]) -> Tuple[frozenset, frozenset]:
        path: Set[Tuple[int, int]] = set()
        habitat_out: Set[Tuple[int, int]] = set()
        current: Optional[Tuple[int, int]] = state
        while current is not None:
            prev, kind = parents.get(current, (None, None))
            r0, c0 = current
            if kind == "path":
                path.add((r0, c0))
            elif kind == "habitat":
                habitat_out.add((r0, c0))
            current = prev
        return frozenset(path), frozenset(habitat_out)

    for r, c in start_positions:
        if obstacle_mask is not None and obstacle_mask[r, c]:
            continue
        if passable_mask is not None and not passable_mask[r, c]:
            continue
        state = (r, c)
        parents[state] = (None, "path")
        heapq.heappush(heap, (0.0, r, c))
        best_cost[state] = 0.0

    nodes_explored = 0
    t_start = time.perf_counter()

    while heap:
        cost, r, c = heapq.heappop(heap)
        state = (r, c)
        prev_best = best_cost.get(state)
        if prev_best is not None and cost > prev_best:
            continue
        if cost > max_width:
            continue
        nodes_explored += 1

        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                lbl = labels[nr, nc]
                if lbl > 0 and lbl not in start_patches:
                    if target_patch_id is not None:
                        if int(lbl) == int(target_patch_id):
                            if target_boundary_center is not None and target_boundary_radius is not None:
                                tr, tc = target_boundary_center
                                if abs(int(nr) - int(tr)) > int(target_boundary_radius) or abs(int(nc) - int(tc)) > int(target_boundary_radius):
                                    continue
                            path_set, habitat_set = _reconstruct(state)
                            return [(path_set, habitat_set, int(lbl), cost)]
                    else:
                        if lbl not in visited_targets:
                            visited_targets.add(lbl)
                            path_set, habitat_set = _reconstruct(state)
                            results.append((path_set, habitat_set, int(lbl), cost))
                    continue

                if obstacle_mask is not None and obstacle_mask[nr, nc]:
                    continue
                if habitat[nr, nc]:
                    if lbl == 0:
                        move_cost = 0.0
                        include_kind = "habitat"  # record connected habitat; no corridor cost
                    else:
                        continue
                else:
                    if passable_mask is not None and not passable_mask[nr, nc]:
                        continue
                    move_cost = SQRT2 if dr != 0 and dc != 0 else 1.0
                    if ghost_mask is not None and ghost_mask[nr, nc]:
                        move_cost *= float(ghost_factor)
                    include_kind = "path"
                new_cost = cost + move_cost
                if new_cost > max_width:
                    continue
                next_state = (nr, nc)
                prev_best = best_cost.get(next_state)
                if prev_best is not None and prev_best <= new_cost:
                    continue
                best_cost[next_state] = new_cost
                parents[next_state] = (state, include_kind)
                heapq.heappush(heap, (new_cost, nr, nc))

    return results



def _compute_start_positions_by_patch(
    labels: np.ndarray,
    habitat: np.ndarray,
    connectivity: int,
    obstacle_mask: Optional[np.ndarray] = None,
    passable_mask: Optional[np.ndarray] = None,
    params: Optional[RasterRunParams] = None,
) -> Dict[int, List[Tuple[int, int]]]:
    """Precompute Dijkstra start positions for each patch id."""
    rows, cols = labels.shape
    if connectivity == 4:
        moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    else:
        moves = [
            (-1, -1),
            (-1, 0),
            (-1, 1),
            (0, -1),
            (0, 1),
            (1, -1),
            (1, 0),
            (1, 1),
        ]

    by_patch: Dict[int, List[Tuple[int, int]]] = defaultdict(list)
    for r in range(rows):
        for c in range(cols):
            if habitat[r, c]:
                continue
            if obstacle_mask is not None and obstacle_mask[r, c]:
                continue
            if passable_mask is not None and not passable_mask[r, c]:
                continue

            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    pid = int(labels[nr, nc])
                    if pid > 0:
                        by_patch[pid].append((r, c))
    return by_patch


def _compute_boundary_coords_by_label(labels: np.ndarray, connectivity: int) -> Dict[int, np.ndarray]:
    """
    Compute boundary pixel coordinates for each positive label id.

    Boundary pixel = a patch pixel that touches background or a different label
    under the given connectivity (4 or 8).
    """
    rows, cols = labels.shape
    boundary = np.zeros_like(labels, dtype=bool)

    def _mark_pair(a: np.ndarray, b: np.ndarray, a_view: Tuple[slice, slice], b_view: Tuple[slice, slice]) -> None:
        mism = a != b
        if mism.any():
            boundary[a_view] |= mism & (a > 0)
            boundary[b_view] |= mism & (b > 0)

    # 4-neighborhood comparisons (mark both sides)
    _mark_pair(labels[:-1, :], labels[1:, :], (slice(0, rows - 1), slice(0, cols)), (slice(1, rows), slice(0, cols)))
    _mark_pair(labels[:, :-1], labels[:, 1:], (slice(0, rows), slice(0, cols - 1)), (slice(0, rows), slice(1, cols)))

    if connectivity != 4:
        # Diagonals
        _mark_pair(labels[:-1, :-1], labels[1:, 1:], (slice(0, rows - 1), slice(0, cols - 1)), (slice(1, rows), slice(1, cols)))
        _mark_pair(labels[:-1, 1:], labels[1:, :-1], (slice(0, rows - 1), slice(1, cols)), (slice(1, rows), slice(0, cols - 1)))

    ys, xs = np.nonzero(boundary)
    if len(ys) == 0:
        return {}
    ids = labels[ys, xs].astype(np.int32, copy=False)
    keep = ids > 0
    ys = ys[keep]
    xs = xs[keep]
    ids = ids[keep]
    if ids.size == 0:
        return {}

    order = np.argsort(ids, kind="mergesort")
    ids_s = ids[order]
    ys_s = ys[order]
    xs_s = xs[order]

    out: Dict[int, np.ndarray] = {}
    unique_ids, starts = np.unique(ids_s, return_index=True)
    for i, pid in enumerate(unique_ids):
        start = int(starts[i])
        end = int(starts[i + 1]) if i + 1 < len(starts) else int(ids_s.size)
        coords = np.stack([ys_s[start:end], xs_s[start:end]], axis=1).astype(np.int32, copy=False)
        out[int(pid)] = coords
    return out


def _subsample_coords_grid(
    coords: np.ndarray,
    block_size: int,
    max_points: Optional[int] = None,
) -> np.ndarray:
    """
    Subsample coordinates by keeping at most one point per block in a coarse grid.
    This preserves local coverage (good for "obvious gaps") without exploding point count.
    """
    if coords is None or coords.size == 0:
        return coords
    bs = max(1, int(block_size))
    rr = (coords[:, 0] // bs).astype(np.int64, copy=False)
    cc = (coords[:, 1] // bs).astype(np.int64, copy=False)
    keys = (rr << 32) ^ cc
    _, idx = np.unique(keys, return_index=True)
    sub = coords[np.sort(idx)]
    if max_points is not None and sub.shape[0] > int(max_points):
        stride = max(1, int(sub.shape[0] // int(max_points)))
        sub = sub[::stride]
    return sub


def _filter_start_positions_near(
    starts: List[Tuple[int, int]],
    center: Tuple[int, int],
    radius: int,
) -> List[Tuple[int, int]]:
    if not starts:
        return []
    cr, cc = int(center[0]), int(center[1])
    rad = max(1, int(radius))
    out: List[Tuple[int, int]] = []
    for r, c in starts:
        if abs(int(r) - cr) <= rad and abs(int(c) - cc) <= rad:
            out.append((r, c))
    return out


def _estimate_search_work_units(estimated_searches: int, total_pixels: int) -> float:
    """Estimate search work units for circuit candidate generation."""
    if estimated_searches <= 0:
        return 0.0
    # Scale by ~250k-pixel blocks (heuristic for keeping runs under a few minutes).
    pixels_factor = max(1.0, float(total_pixels) / 250_000.0)
    return float(estimated_searches) * pixels_factor


def _maybe_block_large_candidate_search(
    *,
    strategy_label: str,
    estimated_searches: int,
    candidate_pairs: int,
    total_pixels: int,
    min_patch_size: int,
    max_search_distance: int,
    min_corridor_width: int,
) -> None:
    """Raise if estimated candidate search is likely to exceed practical runtime."""
    work_units = _estimate_search_work_units(estimated_searches, total_pixels)
    if work_units <= MAX_CANDIDATE_SEARCH_WORK:
        return
    msg = (
        f"Corridor search too large for {strategy_label} (est. {estimated_searches:,} searches across "
        f"{candidate_pairs:,} pairs for {total_pixels:,} raster pixels; "
        f"min patch {min_patch_size}px, max search {max_search_distance}px, "
        f"min width {min_corridor_width}px). "
        "Try increasing min patch size, lowering max search distance, "
        "reducing corridor width, clipping the study area, or resampling to coarser resolution."
    )
    raise RasterAnalysisError(msg)


def _generate_boundary_pair_seeds(
    labels: np.ndarray,
    patch_ids: List[int],
    connectivity: int,
    max_search_distance: int,
    min_corridor_width: int,
    max_seeds_per_pair: int = 8,
) -> Dict[Tuple[int, int], List[Tuple[Tuple[int, int], Tuple[int, int], float]]]:
    """
    Build a "gap-aware" list of nearby patch pairs using boundary pixels and spatial hashing.

    Returns a dict: (p1, p2) -> [(seedA_rc, seedB_rc, euclid_dist_px), ...]
    """
    max_d = max(1, int(max_search_distance))
    max_d2 = float(max_d * max_d)

    # Boundary pixels per patch (then subsample for speed but keep local coverage).
    t0 = time.perf_counter()
    boundary_by_pid = _compute_boundary_coords_by_label(labels, connectivity=connectivity)
    _ = t0  # local timing handled by caller if desired

    # Subsample settings: smaller blocks preserve micro-gaps; cap prevents explosion.
    block_size = 2
    cap_per_patch = 2000
    pts: List[Tuple[int, int, int]] = []
    for pid in patch_ids:
        coords = boundary_by_pid.get(int(pid))
        if coords is None or coords.size == 0:
            continue
        sub = _subsample_coords_grid(coords, block_size=block_size, max_points=cap_per_patch)
        for r, c in sub:
            pts.append((int(pid), int(r), int(c)))
    if not pts:
        return {}

    # Spatial bins
    bin_size = max(8, min(64, max_d))
    bins: Dict[Tuple[int, int], List[Tuple[int, int, int]]] = defaultdict(list)
    for pid, r, c in pts:
        bins[(r // bin_size, c // bin_size)].append((pid, r, c))

    # Collect closest boundary-to-boundary seeds per pair.
    # We keep a small candidate list per pair, then enforce spatial diversity later.
    pair_cands: Dict[Tuple[int, int], List[Tuple[float, Tuple[int, int], Tuple[int, int]]]] = defaultdict(list)

    neighbor_offsets = [(dr, dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1)]
    for (br, bc), cell_pts in bins.items():
        # Compare points in this bin to points in neighboring bins.
        neighbor_pts: List[Tuple[int, int, int]] = []
        for dr, dc in neighbor_offsets:
            neighbor_pts.extend(bins.get((br + dr, bc + dc), []))
        if not neighbor_pts:
            continue

        for pid1, r1, c1 in cell_pts:
            for pid2, r2, c2 in neighbor_pts:
                if pid2 <= pid1:
                    continue
                dr = float(r1 - r2)
                dc = float(c1 - c2)
                d2 = dr * dr + dc * dc
                if d2 > max_d2:
                    continue
                pk = (pid1, pid2) if pid1 <= pid2 else (pid2, pid1)
                pair_cands[pk].append((d2, (r1, c1), (r2, c2)))
                # Soft cap per pair while collecting to keep this step bounded.
                if len(pair_cands[pk]) > 200:
                    pair_cands[pk].sort(key=lambda t: t[0])
                    pair_cands[pk] = pair_cands[pk][:50]

    if not pair_cands:
        return {}

    # Enforce "boundary states" diversity: avoid picking seeds that are basically the same spot.
    # The diversity radius scales with corridor width to avoid adjacent-parallel redundancies.
    diversity_radius = max(4, int(math.ceil(float(min_corridor_width) * 1.5)))
    out: Dict[Tuple[int, int], List[Tuple[Tuple[int, int], Tuple[int, int], float]]] = {}
    for pk, items in pair_cands.items():
        items.sort(key=lambda t: t[0])
        chosen: List[Tuple[Tuple[int, int], Tuple[int, int], float]] = []
        for d2, a, b in items:
            if len(chosen) >= int(max_seeds_per_pair):
                break
            ok = True
            for a0, b0, _ in chosen:
                if abs(a[0] - a0[0]) <= diversity_radius and abs(a[1] - a0[1]) <= diversity_radius:
                    ok = False
                    break
                if abs(b[0] - b0[0]) <= diversity_radius and abs(b[1] - b0[1]) <= diversity_radius:
                    ok = False
                    break
            if not ok:
                continue
            chosen.append((a, b, float(math.sqrt(d2))))
        if chosen:
            out[pk] = chosen
    return out


def find_all_possible_corridors(
    labels: np.ndarray,
    habitat: np.ndarray,
    patch_sizes: Dict[int, int],
    max_width: int,
    max_area: Optional[int],
    connectivity: int,
    min_corridor_width: int,
    obstacle_mask: Optional[np.ndarray] = None,
    passable_mask: Optional[np.ndarray] = None,
    strategy: str = "circuit_utility",
    progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
    progress_start: int = 45,
    progress_end: int = 75,
    params: Optional[RasterRunParams] = None,
    timing_out: Optional[Dict[str, object]] = None,
) -> List[Dict]:
    """Find all possible corridors between patch pairs."""
    print("  Finding all possible corridors...")
    rows, cols = labels.shape
    total_pixels = int(rows * cols)
    all_corridors: List[Dict] = []
    processed_pairs: Set[frozenset] = set()
    offsets = _corridor_offsets(min_corridor_width)

    unique_patches = [p for p in patch_sizes.keys() if p > 0]
    total = len(unique_patches) or 1
    span = max(progress_end - progress_start, 1)

    timing_durations: Dict[str, float] = defaultdict(float)
    timing_counts: Dict[str, int] = defaultdict(int)
    timing_start = time.perf_counter()

    def _record_timing(label: str, start_s: float) -> None:
        try:
            timing_durations[label] += time.perf_counter() - start_s
            timing_counts[label] += 1
        except Exception:
            pass

    strategy_key = (strategy or "circuit_utility").lower()
    if strategy_key not in ("largest_network", "circuit_utility"):
        strategy_key = "circuit_utility"
    circuit_mode = strategy_key == "circuit_utility"
    # Largest Single Network mode keeps the fast multi-target search, but adds a targeted micro-gap pass
    # to catch obvious missed bridges (boundary-aware + gap-aware).
    micro_gap_mode = strategy_key in ("largest_network",)

    # Circuit mode: keep multiple "boundary-state" variants per pair (no ghosting).
    max_keep_per_pair = 4 if circuit_mode else 1
    overlap_keep_ratio = 0.75

    start_positions_by_patch: Dict[int, List[Tuple[int, int, int, bool, bool]]] = {}
    if circuit_mode or micro_gap_mode:
        t0 = time.perf_counter()
        start_positions_by_patch = _compute_start_positions_by_patch(
            labels,
            habitat,
            connectivity,
            obstacle_mask=obstacle_mask,
            passable_mask=passable_mask,
            params=params,
        )
        _record_timing("precompute_start_positions", t0)

    pair_buffers: Dict[Tuple[int, int], List[frozenset]] = defaultdict(list)
    pair_buffers_expanded: Dict[Tuple[int, int], List[frozenset]] = defaultdict(list)
    pair_counts: Dict[Tuple[int, int], int] = defaultdict(int)

    # Proximity-aware overlap expansion: treat near-parallel corridors as overlap/thickening.
    prox_radius = max(1, int(math.ceil(float(min_corridor_width) / 2.0)))
    expand_offsets = [(dr, dc) for dr in range(-prox_radius, prox_radius + 1) for dc in range(-prox_radius, prox_radius + 1)]

    def _expand_pixels(pixels: Set[Tuple[int, int]]) -> frozenset:
        expanded: Set[Tuple[int, int]] = set()
        for r, c in pixels:
            for dr, dc in expand_offsets:
                expanded.add((r + dr, c + dc))
        return frozenset(expanded)

    def _pair_key(a: int, b: int) -> Tuple[int, int]:
        return (a, b) if a <= b else (b, a)

    def _max_overlap_ratio(new_buf: Set[Tuple[int, int]], prior_expanded: List[frozenset]) -> float:
        if not new_buf or not prior_expanded:
            return 0.0
        denom = float(len(new_buf))
        if denom <= 0:
            return 0.0
        best = 0.0
        new_f = frozenset(new_buf)
        for prev in prior_expanded:
            try:
                ratio = len(new_f.intersection(prev)) / denom
            except Exception:
                continue
            if ratio > best:
                best = ratio
        return best

    def _add_candidate(
        p1: int,
        p2: int,
        path_pixels: frozenset,
        habitat_pixels: frozenset,
        path_len: float,
        tag: str,
    ) -> None:
        if max_area is not None and path_len > max_area:
            return
        t0 = time.perf_counter()
        buffered = _inflate_corridor_pixels(set(path_pixels), offsets, rows, cols, obstacle_mask=obstacle_mask)
        _record_timing("inflate_pixels", t0)
        area_px = len(buffered)
        pk = _pair_key(p1, p2)
        if circuit_mode:
            if pair_counts[pk] >= max_keep_per_pair:
                return
            t1 = time.perf_counter()
            overlap = _max_overlap_ratio(buffered, pair_buffers_expanded.get(pk, []))
            _record_timing("overlap_check", t1)
            if overlap >= overlap_keep_ratio:
                return
        pair_counts[pk] += 1
        buf_f = frozenset(buffered)
        pair_buffers[pk].append(buf_f)
        pair_buffers_expanded[pk].append(_expand_pixels(buffered))
        all_corridors.append(
            {
                "patch1": p1,
                "patch2": p2,
                "pixels": path_pixels,
                "habitat_pixels": habitat_pixels,
                "length": path_len,
                "area": area_px,
                "buffered_pixels": frozenset(buffered),
                "variant": tag,
                # Used by Most Connectivity overlap logic to treat near-parallel corridors as overlap.
                "_prox_radius": int(prox_radius),
            }
        )

    if circuit_mode:
        # ------------------------------------------------------------------
        # MOST CONNECTIVITY: boundary-aware + gap-aware candidate generation
        #
        # Instead of "search until you touch any patch", we:
        #  - build boundary-to-boundary seed pairs for nearby patch pairs
        #  - run targeted A->B searches (terminate only when reaching B)
        #  - keep spatially distinct variants; treat near-parallel as overlap
        # ------------------------------------------------------------------
        patch_ids = [int(p) for p in unique_patches]
        t_seed = time.perf_counter()
        pair_seeds = _generate_boundary_pair_seeds(
            labels=labels,
            patch_ids=patch_ids,
            connectivity=connectivity,
            max_search_distance=int(max_width),
            min_corridor_width=int(min_corridor_width),
            max_seeds_per_pair=8,
        )
        _record_timing("pair_seed_build", t_seed)
        if pair_seeds:
            est_searches = sum(len(seeds) for seeds in pair_seeds.values())
            _maybe_block_large_candidate_search(
                strategy_label="Circuit Utility",
                estimated_searches=int(est_searches),
                candidate_pairs=len(pair_seeds),
                total_pixels=total_pixels,
                min_patch_size=int(params.min_patch_size) if params is not None else 0,
                max_search_distance=int(max_width),
                min_corridor_width=int(min_corridor_width),
            )
        if not pair_seeds:
            try:
                _log_message(
                    "[Circuit] No boundary pair seeds found; falling back to baseline LCP candidate generation.",
                    Qgis.Warning if Qgis is not None else 1,
                )
            except Exception:
                pass

        if not pair_seeds:
            # Baseline fallback: old behavior (no ghosting), used only when boundary seeds fail.
            for idx, patch_id in enumerate(unique_patches):
                if progress_cb is not None:
                    pre_value = progress_start + (idx / total) * span
                    _emit_progress(progress_cb, pre_value, "Analyzing patches")

                start_override = start_positions_by_patch.get(int(patch_id)) if circuit_mode else None
                t0 = time.perf_counter()
                results = find_shortest_corridor(
                    {int(patch_id)},
                    labels,
                    habitat,
                    int(max_width),
                    connectivity,
                    obstacle_mask=obstacle_mask,
                    passable_mask=passable_mask,
                    params=params,
                    start_positions_override=start_override,
                )
                _record_timing("lcp_search", t0)
                if not results:
                    continue
                for path_pixels, habitat_pixels, target_id, path_len in results:
                    if int(patch_id) > int(target_id):
                        continue
                    _add_candidate(
                        int(patch_id),
                        int(target_id),
                        path_pixels,
                        habitat_pixels,
                        path_len,
                        tag="lcp",
                    )
        else:
            total_tasks = sum(len(v) for v in pair_seeds.values()) or 1
            done = 0

            # Start filtering radius scales with patch width to force "closest edge" behavior.
            start_filter_radius = max(6, int(math.ceil(float(min_corridor_width) * 2.0)))
            target_filter_radius = max(6, int(math.ceil(float(min_corridor_width) * 2.0)))

            for (p1i, p2i), seeds in sorted(pair_seeds.items()):
                if p1i not in patch_sizes or p2i not in patch_sizes:
                    continue
                pk = _pair_key(p1i, p2i)
                for s_idx, (seed_a, seed_b, _d) in enumerate(seeds, 1):
                    if pair_counts[pk] >= max_keep_per_pair:
                        break
                    if progress_cb is not None:
                        value = progress_start + (done / total_tasks) * span
                        _emit_progress(progress_cb, value, "Analyzing patches")

                    start_override = start_positions_by_patch.get(int(p1i), [])
                    filtered = _filter_start_positions_near(start_override, seed_a, radius=start_filter_radius)
                    if filtered:
                        start_override_use = filtered
                    else:
                        start_override_use = start_override

                    t0 = time.perf_counter()
                    res = find_shortest_corridor(
                        {int(p1i)},
                        labels,
                        habitat,
                        int(max_width),
                        connectivity,
                        obstacle_mask=obstacle_mask,
                        passable_mask=passable_mask,
                        params=params,
                        target_patch_id=int(p2i),
                        target_boundary_center=seed_b,
                        target_boundary_radius=target_filter_radius,
                        start_positions_override=start_override_use,
                    )
                    _record_timing("target_search", t0)
                    if (not res) and s_idx == 1 and pair_counts[pk] == 0:
                        # Fallback for pairs whose "closest boundary" is blocked by obstacles:
                        # allow reaching the target patch anywhere for at least one baseline candidate.
                        t_fb = time.perf_counter()
                        res = find_shortest_corridor(
                            {int(p1i)},
                            labels,
                            habitat,
                            int(max_width),
                            connectivity,
                            obstacle_mask=obstacle_mask,
                            passable_mask=passable_mask,
                            params=params,
                            target_patch_id=int(p2i),
                            start_positions_override=start_override_use,
                        )
                        _record_timing("target_search_fallback", t_fb)
                    if res:
                        path_pixels, habitat_pixels, _tid, path_len = res[0]
                        _add_candidate(
                            p1i, p2i, path_pixels, habitat_pixels, path_len, tag=f"seed_{s_idx}"
                        )
                    done += 1

        _emit_progress(progress_cb, progress_end, "Corridor candidates ready.")
    else:
        for idx, patch_id in enumerate(unique_patches):
            if progress_cb is not None:
                pre_value = progress_start + (idx / total) * span
                _emit_progress(
                    progress_cb,
                    pre_value,
                    "Analyzing patches",
                )

            t0 = time.perf_counter()
            start_override = start_positions_by_patch.get(int(patch_id)) if micro_gap_mode else None
            results = find_shortest_corridor(
                {patch_id},
                labels,
                habitat,
                max_width,
                connectivity,
                obstacle_mask=obstacle_mask,
                passable_mask=passable_mask,
                params=params,
                start_positions_override=start_override,
            )
            _record_timing("lcp_search", t0)
            if not results:
                continue

            for path_pixels, habitat_pixels, target_id, path_len in results:
                pair = frozenset({patch_id, target_id})
                if pair in processed_pairs:
                    continue
                processed_pairs.add(pair)
                _add_candidate(
                    int(patch_id),
                    int(target_id),
                    path_pixels,
                    habitat_pixels,
                    path_len,
                    tag="lcp",
                )

            if progress_cb is not None:
                post_value = progress_start + ((idx + 1) / total) * span
                _emit_progress(progress_cb, post_value, "Analyzing patches")

        # ------------------------------------------------------------------
        # MICRO-GAP AUGMENTATION (boundary-aware + gap-aware)
        #
        # Goal: catch "obvious" near-boundary gaps that the multi-target search
        # can miss when another patch is touched first.
        # ------------------------------------------------------------------
        if micro_gap_mode:
            try:
                micro_gap_dist = min(int(max_width), max(12, int(min_corridor_width) * 3))
                if micro_gap_dist > 0:
                    patch_ids = [int(p) for p in unique_patches]
                    t_seed = time.perf_counter()
                    micro_seeds = _generate_boundary_pair_seeds(
                        labels=labels,
                        patch_ids=patch_ids,
                        connectivity=connectivity,
                        max_search_distance=int(micro_gap_dist),
                        min_corridor_width=int(min_corridor_width),
                        max_seeds_per_pair=2,
                    )
                    _record_timing("micro_seed_build", t_seed)
                    try:
                        max_micro_pairs_total = 4000
                        if len(micro_seeds) > max_micro_pairs_total:
                            items = list(micro_seeds.items())
                            items.sort(key=lambda kv: float(kv[1][0][2]) if kv[1] else 1e18)
                            micro_seeds = dict(items[:max_micro_pairs_total])
                            _log_message(
                                f"[Micro-gap] Pairs capped to {max_micro_pairs_total} (of {len(items)}) for performance.",
                                Qgis.Info if Qgis is not None else 0,
                            )
                    except Exception:
                        pass

                    # Nudge progress forward (don't spam patch-level logs).
                    micro_pairs = len(micro_seeds) or 1
                    done_pairs = 0
                    micro_value_start = max(progress_start, min(progress_end - 1, progress_start + int(span * 0.90)))
                    micro_span = max(1, progress_end - micro_value_start)

                    start_filter_radius = max(6, int(math.ceil(float(min_corridor_width) * 2.0)))
                    target_filter_radius = max(6, int(math.ceil(float(min_corridor_width) * 2.0)))

                    for (p1i, p2i), seeds in sorted(micro_seeds.items()):
                        pair_fs = frozenset({int(p1i), int(p2i)})
                        if pair_fs in processed_pairs:
                            continue
                        if p1i not in patch_sizes or p2i not in patch_sizes:
                            continue

                        if progress_cb is not None:
                            value = micro_value_start + int((done_pairs / micro_pairs) * micro_span)
                            _emit_progress(progress_cb, value, "Analyzing patches")

                        best: Optional[Tuple[frozenset, frozenset, float, str]] = None
                        start_override = start_positions_by_patch.get(int(p1i), [])
                        for s_idx, (seed_a, seed_b, _d) in enumerate(seeds, 1):
                            filtered = _filter_start_positions_near(start_override, seed_a, radius=start_filter_radius)
                            start_override_use = filtered if filtered else start_override

                            t0 = time.perf_counter()
                            res = find_shortest_corridor(
                                {int(p1i)},
                                labels,
                                habitat,
                                int(max_width),
                                connectivity,
                                obstacle_mask=obstacle_mask,
                                passable_mask=passable_mask,
                                params=params,
                                target_patch_id=int(p2i),
                                target_boundary_center=seed_b,
                                target_boundary_radius=target_filter_radius,
                                start_positions_override=start_override_use,
                            )
                            _record_timing("micro_target_search", t0)
                            if (not res) and s_idx == 1:
                                # Fallback: allow reaching the target patch anywhere (once).
                                t_fb = time.perf_counter()
                                res = find_shortest_corridor(
                                    {int(p1i)},
                                    labels,
                                    habitat,
                                    int(max_width),
                                    connectivity,
                                    obstacle_mask=obstacle_mask,
                                    passable_mask=passable_mask,
                                    params=params,
                                    target_patch_id=int(p2i),
                                    start_positions_override=start_override_use,
                                )
                                _record_timing("micro_target_search_fallback", t_fb)

                            if not res:
                                continue
                            path_pixels, habitat_pixels, _tid, path_len = res[0]
                            tag = f"micro_gap_seed_{s_idx}"
                            if best is None or float(path_len) < float(best[2]):
                                best = (path_pixels, habitat_pixels, float(path_len), tag)

                        if best is not None:
                            path_pixels, habitat_pixels, path_len, tag = best
                            _add_candidate(
                                int(p1i),
                                int(p2i),
                                path_pixels,
                                habitat_pixels,
                                float(path_len),
                                tag=tag,
                            )
                            processed_pairs.add(pair_fs)
                        done_pairs += 1
            except Exception:
                pass

    if not circuit_mode:
        _emit_progress(progress_cb, progress_end, "Corridor candidates ready.")
    print(f"\n  ✓ Found {len(all_corridors)} possible corridors")
    timing_durations["total"] = time.perf_counter() - timing_start
    if timing_out is not None and circuit_mode:
        timing_out["durations_s"] = dict(timing_durations)
        timing_out["counts"] = dict(timing_counts)
        timing_out["candidates"] = len(all_corridors)
    return all_corridors


def _corridor_offsets(min_corridor_width: int) -> List[Tuple[int, int]]:
    """Precompute offsets used to inflate corridors to the minimum width."""
    width = max(1, int(min_corridor_width))
    if width <= 1:
        return [(0, 0)]

    radius = max(0.0, width / 2.0)
    max_offset = int(math.ceil(radius))
    radius_sq = radius * radius

    offsets: List[Tuple[int, int]] = []
    for dr in range(-max_offset, max_offset + 1):
        for dc in range(-max_offset, max_offset + 1):
            if dr * dr + dc * dc <= radius_sq + 1e-9:
                offsets.append((dr, dc))

    if not offsets:
        offsets.append((0, 0))
    return offsets


def _inflate_corridor_pixels(
    pixels: Set[Tuple[int, int]],
    offsets: List[Tuple[int, int]],
    rows: int,
    cols: int,
    obstacle_mask: Optional[np.ndarray] = None,
) -> Set[Tuple[int, int]]:
    """Apply offsets to centerline pixels, respecting impassables when provided."""
    inflated: Set[Tuple[int, int]] = set()
    use_mask = obstacle_mask is not None

    for r, c in pixels:
        for dr, dc in offsets:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                if use_mask and obstacle_mask[nr, nc]:
                    continue
                inflated.add((nr, nc))
    return inflated


def _shift_mask(mask: np.ndarray, dr: int, dc: int) -> np.ndarray:
    """Return a shifted copy of mask aligned so (r, c) reads (r+dr, c+dc) from original."""
    rows, cols = mask.shape
    shifted = np.zeros_like(mask, dtype=bool)

    if abs(dr) >= rows or abs(dc) >= cols:
        return shifted

    if dr >= 0:
        src_r = slice(dr, rows)
        dst_r = slice(0, rows - dr)
    else:
        src_r = slice(0, rows + dr)
        dst_r = slice(-dr, rows)

    if dc >= 0:
        src_c = slice(dc, cols)
        dst_c = slice(0, cols - dc)
    else:
        src_c = slice(0, cols + dc)
        dst_c = slice(-dc, cols)

    shifted[dst_r, dst_c] = mask[src_r, src_c]
    return shifted


def _erode_mask(mask: np.ndarray, offsets: List[Tuple[int, int]]) -> np.ndarray:
    """Morphologically erode a mask using the provided offsets."""
    if not offsets:
        return mask.copy()

    eroded = mask.copy()
    for dr, dc in offsets:
        if dr == 0 and dc == 0:
            continue
        shifted = _shift_mask(mask, dr, dc)
        eroded &= shifted
        if not eroded.any():
            break
    return eroded


def _build_passable_mask(
    habitat: np.ndarray,
    obstacle_mask: np.ndarray,
    min_corridor_width: int,
    allow_bottlenecks: bool,
) -> np.ndarray:
    """Derive which pixels can host corridor centerlines under width constraints."""
    habitat_bool = habitat.astype(bool)
    obstacle_bool = obstacle_mask.astype(bool) if obstacle_mask is not None else np.zeros_like(habitat_bool)
    base_passable = (~habitat_bool) & (~obstacle_bool)

    if allow_bottlenecks or min_corridor_width <= 1:
        passable = base_passable
    else:
        offsets = _corridor_offsets(min_corridor_width)
        # Ignore habitat when evaluating clearance so corridors can still touch patches.
        true_obstacles = obstacle_bool & (~habitat_bool)
        clearance_space = ~true_obstacles
        clearance_ok = _erode_mask(clearance_space, offsets)
        passable = base_passable & clearance_ok

    return passable


def create_output_raster(
    labels: np.ndarray,
    corridors: Dict[int, Dict],
    min_corridor_width: int,
    obstacle_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Create an output raster with corridors marked by connected size."""
    output = np.zeros_like(labels, dtype=np.int32)
    rows, cols = labels.shape
    use_mask = obstacle_mask is not None

    for corridor_data in corridors.values():
        score = corridor_data["connected_size"]
        buffered = corridor_data.get("buffered_pixels")
        if buffered:
            for r, c in buffered:
                if 0 <= r < rows and 0 <= c < cols:
                    if use_mask and obstacle_mask[r, c]:
                        continue
                    if score > output[r, c]:
                        output[r, c] = score
            continue
        offsets = _corridor_offsets(min_corridor_width)
        for r, c in corridor_data["pixels"]:
            for dr, dc in offsets:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    if use_mask and obstacle_mask[nr, nc]:
                        continue
                    if score > output[nr, nc]:
                        output[nr, nc] = score
    return output


def create_combined_raster(
    labels: np.ndarray,
    corridors: Dict[int, Dict],
    min_corridor_width: int,
    obstacle_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Create a raster where habitat patches and corridors are merged and labeled by connectivity."""
    rows, cols = labels.shape
    use_mask = obstacle_mask is not None
    combined = np.zeros_like(labels, dtype=np.int32)

    # Start with patches (existing labels)
    combined = labels.copy().astype(np.int32)

    # Burn corridors as a unique label (negative), then re-label components
    next_label = int(labels.max()) + 1
    corridor_mask = np.zeros_like(labels, dtype=bool)
    for corridor_data in corridors.values():
        buffered = corridor_data.get("buffered_pixels")
        if buffered:
            for r, c in buffered:
                if 0 <= r < rows and 0 <= c < cols:
                    if use_mask and obstacle_mask is not None and obstacle_mask[r, c]:
                        continue
                    corridor_mask[r, c] = True
            continue
        offsets = _corridor_offsets(min_corridor_width)
        for r, c in corridor_data["pixels"]:
            for dr, dc in offsets:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    if use_mask and obstacle_mask is not None and obstacle_mask[nr, nc]:
                        continue
                    corridor_mask[nr, nc] = True

    combined[corridor_mask] = -1  # temporary mark corridors

    if HAS_NDIMAGE and ndimage is not None:
        habitat_corr = combined != 0
        structure = np.ones((3, 3), dtype=np.uint8)
        labeled, _ = ndimage.label(habitat_corr, structure=structure)
        labeled[combined == 0] = 0
        return labeled.astype(np.int32)

    # Fallback: return patches with corridor mask burned as a unique label
    combined[corridor_mask] = labels.max() + 1 if labels.size > 0 else 1
    return combined.astype(np.int32)


def create_contiguous_network_area_raster(
    labels: np.ndarray,
    corridors: Dict[int, Dict],
    min_corridor_width: int,
    obstacle_mask: Optional[np.ndarray] = None,
    labels_all: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Build a raster where every habitat/corridor cell value equals the total area (pixels)
    of its connected component.

    This includes all habitat patches in `labels` (even if not part of the corridor network).

    Additionally, corridor pathfinding may traverse habitat pixels that are not in `labels`
    (e.g., small filtered-out patches); if `labels_all` is provided and a traversed habitat
    pixel belongs to a filtered-out patch in `labels_all`, the entire patch footprint is
    included in the component before computing its area.
    """
    rows, cols = labels.shape
    if rows == 0 or cols == 0:
        return labels.astype(np.int32)

    patch_ids = np.unique(labels[labels > 0]).astype(int)
    if patch_ids.size == 0:
        return np.zeros_like(labels, dtype=np.int32)

    uf = UnionFind()
    for pid in patch_ids.tolist():
        uf.find(int(pid))
    for corr in corridors.values():
        pids = list(corr.get("patch_ids", []))
        if len(pids) >= 2:
            uf.union(int(pids[0]), int(pids[1]))

    root_to_idx: Dict[int, int] = {}
    next_idx = 1
    for pid in patch_ids.tolist():
        root = int(uf.find(int(pid)))
        if root not in root_to_idx:
            root_to_idx[root] = next_idx
            next_idx += 1

    max_label = int(labels.max())
    lookup = np.zeros(max_label + 1, dtype=np.int32)
    for pid in patch_ids.tolist():
        lookup[int(pid)] = root_to_idx[int(uf.find(int(pid)))]

    network_labels = lookup[labels.astype(np.int32)]

    use_mask = obstacle_mask is not None
    for corr in corridors.values():
        pids = list(corr.get("patch_ids", []))
        if len(pids) < 1:
            continue
        comp_idx = root_to_idx.get(int(uf.find(int(pids[0]))))
        if comp_idx is None:
            continue
        buffered = corr.get("buffered_pixels")
        if buffered:
            for r, c in buffered:
                if 0 <= r < rows and 0 <= c < cols:
                    if use_mask and obstacle_mask is not None and obstacle_mask[r, c]:
                        continue
                    network_labels[r, c] = comp_idx
        else:
            offsets = _corridor_offsets(min_corridor_width)
            for r, c in corr.get("pixels", frozenset()):
                for dr, dc in offsets:
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        if use_mask and obstacle_mask is not None and obstacle_mask[nr, nc]:
                            continue
                        network_labels[nr, nc] = comp_idx

        # Include connected habitat pixels used by the corridor path.
        # This ensures small "bridge" patches are present in the contiguous output.
        touched_old_ids: Set[int] = set()
        for r, c in corr.get("habitat_pixels", frozenset()) or frozenset():
            if 0 <= r < rows and 0 <= c < cols:
                if use_mask and obstacle_mask is not None and obstacle_mask[r, c]:
                    continue
                network_labels[r, c] = comp_idx
                if labels_all is not None:
                    try:
                        if int(labels[r, c]) == 0:
                            old_id = int(labels_all[r, c])
                            if old_id > 0:
                                touched_old_ids.add(old_id)
                    except Exception:
                        pass

        if labels_all is not None and touched_old_ids:
            for old_id in touched_old_ids:
                try:
                    mask = (labels_all == int(old_id)) & (labels == 0)
                    network_labels[mask] = comp_idx
                except Exception:
                    continue

    # Convert component ids -> component areas in pixels (per your requested semantics).
    counts = np.bincount(network_labels.ravel())
    out = counts[network_labels]
    out[network_labels == 0] = 0
    return out.astype(np.int32)


def _compute_connectivity_metrics(corridors: Dict[int, Dict], patch_sizes: Dict[int, int]) -> Dict[str, float]:
    """Compute consistent graph metrics for summaries across strategies."""
    n_nodes = len(patch_sizes)
    uf = UnionFind()
    for pid in patch_sizes:
        uf.find(pid)

    corridors_used = len(corridors)
    for data in corridors.values():
        p1 = data.get("p1") or data.get("patch1")
        p2 = data.get("p2") or data.get("patch2")
        if p1 is not None and p2 is not None:
            uf.union(int(p1), int(p2))
        else:
            pids = list(data.get("patch_ids", []))
            if len(pids) >= 2:
                uf.union(int(pids[0]), int(pids[1]))

    components = len({uf.find(pid) for pid in patch_sizes})
    redundant_links = max(0, corridors_used - (n_nodes - components))
    avg_degree = (2 * corridors_used / n_nodes) if n_nodes > 0 else 0.0

    return {
        "patches_total": n_nodes,
        "corridors_used": corridors_used,
        "components_remaining": components,
        "redundant_links": redundant_links,
        "avg_degree": avg_degree,
    }
def _build_corridor_stats(
    selected: Dict[int, Dict],
    patch_sizes: Dict[int, int],
) -> Tuple[Dict[int, int], Dict[int, int]]:
    """Compute component sizes and counts for selected corridors."""
    uf = UnionFind()
    for pid, size in patch_sizes.items():
        uf.find(pid)
        uf.size[pid] = size
        uf.count[pid] = 1

    for corr in selected.values():
        pids = list(corr.get("patch_ids", []))
        if len(pids) >= 2:
            uf.union(int(pids[0]), int(pids[1]))

    comp_sizes: Dict[int, int] = {}
    comp_counts: Dict[int, int] = {}
    for pid in patch_sizes.keys():
        root = uf.find(pid)
        comp_sizes[root] = comp_sizes.get(root, 0) + int(patch_sizes.get(pid, 0) or 0)
        comp_counts[root] = comp_counts.get(root, 0) + 1

    for corr in selected.values():
        pids = list(corr.get("patch_ids", []))
        if not pids:
            continue
        root = uf.find(int(pids[0]))
        comp_sizes[root] = comp_sizes.get(root, 0) + int(corr.get("area", 0) or 0)

    return comp_sizes, comp_counts


def _redundancy_far_enough_pixels(
    pixels: Sequence[Tuple[int, int]], prior_sets: Sequence[object], threshold: int
) -> bool:
    if threshold <= 0:
        return True
    if not pixels or not prior_sets:
        return True
    cell = max(1, int(threshold))
    index: Dict[Tuple[int, int], List[Tuple[int, int]]] = defaultdict(list)
    for obj in prior_sets:
        if not isinstance(obj, (set, frozenset)):
            continue
        for r, c in obj:
            index[(int(r) // cell, int(c) // cell)].append((int(r), int(c)))
    if not index:
        return True
    thr2 = int(threshold) * int(threshold)
    for r, c in pixels:
        br, bc = int(r) // cell, int(c) // cell
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                for pr, pc in index.get((br + dr, bc + dc), []):
                    if (int(r) - pr) ** 2 + (int(c) - pc) ** 2 < thr2:
                        return False
    return True


def optimize_resilient(
    candidates: List[Dict], patch_sizes: Dict[int, int], budget_pixels: int
) -> Tuple[Dict[int, Dict], Dict]:
    """
    Lightweight resilient strategy: MST-style on corridor cost within budget.
    """
    uf = UnionFind()
    for pid, size in patch_sizes.items():
        uf.find(pid)
        uf.size[pid] = size
        uf.count[pid] = 1

    selected: Dict[int, Dict] = {}
    remaining = max(0, int(budget_pixels))
    budget_used = 0
    redundant_links = 0
    backbone_graph_edges: List[Tuple[int, int, float]] = []

    for cand in sorted(candidates, key=lambda c: c.get("area", 0)):
        p1, p2 = cand.get("patch1"), cand.get("patch2")
        if p1 is None or p2 is None:
            continue
        cost = int(cand.get("area", 0))
        if cost <= 0 or cost > remaining:
            continue
        if uf.find(p1) == uf.find(p2):
            continue

        uf.union(p1, p2)
        remaining -= cost
        budget_used += cost
        comp_root = uf.find(p1)
        conn_size = uf.size.get(comp_root, cost)
        backbone_graph_edges.append((int(p1), int(p2), float(cand.get("length", cost) or cost)))

        cid = len(selected) + 1
        selected[cid] = {
            "pixels": cand.get("pixels", set()),
            "habitat_pixels": cand.get("habitat_pixels", frozenset()),
            "patch_ids": {p1, p2},
            "area": cost,
            "connected_size": conn_size,
            "buffered_pixels": cand.get("buffered_pixels", frozenset()),
            "length": float(cand.get("length", cost) or cost),
            "type": "backbone",
        }

    if remaining > 0 and selected:
        loop_candidates: List[Tuple[float, Dict]] = []
        selected_pairs = {tuple(sorted(map(int, c.get("patch_ids", [])))) for c in selected.values()}

        if nx is not None and graph_math is not None:
            G_backbone = nx.Graph()
            for pid in patch_sizes:
                G_backbone.add_node(int(pid))
            for a, b, w in backbone_graph_edges:
                G_backbone.add_edge(int(a), int(b), weight=float(w))

            for cand in candidates:
                p1, p2 = cand.get("patch1"), cand.get("patch2")
                if p1 is None or p2 is None:
                    continue
                p1i, p2i = int(p1), int(p2)
                if uf.find(p1i) != uf.find(p2i):
                    continue
                pair = tuple(sorted((p1i, p2i)))
                if pair in selected_pairs:
                    continue
                cost = int(cand.get("area", 0))
                if cost <= 0 or cost > remaining:
                    continue
                length = float(cand.get("length", 0) or 0.0)
                if length <= 0:
                    length = float(cost)
                try:
                    score = float(graph_math.score_edge_for_loops(G_backbone, p1i, p2i, length))
                except Exception:
                    score = 0.0
                efficiency = score / float(cost) if cost else 0.0
                loop_candidates.append((efficiency, cand))

            loop_candidates.sort(key=lambda t: t[0], reverse=True)

            for _eff, cand in loop_candidates:
                cost = int(cand.get("area", 0))
                if cost <= 0 or cost > remaining:
                    continue
                p1, p2 = int(cand.get("patch1")), int(cand.get("patch2"))
                pair = tuple(sorted((p1, p2)))
                if pair in selected_pairs:
                    continue
                remaining -= cost
                budget_used += cost
                redundant_links += 1
                selected_pairs.add(pair)
                length = float(cand.get("length", 0) or 0.0)
                if length <= 0:
                    length = float(cost)
                try:
                    G_backbone.add_edge(p1, p2, weight=length)
                except Exception:
                    pass

                cid = len(selected) + 1
                selected[cid] = {
                    "pixels": cand.get("pixels", set()),
                    "habitat_pixels": cand.get("habitat_pixels", frozenset()),
                    "patch_ids": {p1, p2},
                    "area": cost,
                    "connected_size": uf.size.get(uf.find(p1), cost),
                    "buffered_pixels": cand.get("buffered_pixels", frozenset()),
                    "length": float(cand.get("length", cost) or cost),
                    "type": "strategic_loop",
                }
        else:
            for cand in sorted(candidates, key=lambda c: c.get("area", 0)):
                if remaining <= 0:
                    break
                p1, p2 = cand.get("patch1"), cand.get("patch2")
                if p1 is None or p2 is None:
                    continue
                if uf.find(p1) != uf.find(p2):
                    continue
                pair = tuple(sorted((int(p1), int(p2))))
                if pair in selected_pairs:
                    continue
                cost = int(cand.get("area", 0))
                if cost <= 0 or cost > remaining:
                    continue
                remaining -= cost
                budget_used += cost
                redundant_links += 1
                selected_pairs.add(pair)
                cid = len(selected) + 1
                selected[cid] = {
                    "pixels": cand.get("pixels", set()),
                    "habitat_pixels": cand.get("habitat_pixels", frozenset()),
                    "patch_ids": {p1, p2},
                    "area": cost,
                    "connected_size": uf.size.get(uf.find(p1), cost),
                    "buffered_pixels": cand.get("buffered_pixels", frozenset()),
                    "length": float(cand.get("length", cost) or cost),
                    "type": "redundant",
                }

    comp_sizes, comp_counts = _build_corridor_stats(selected, patch_sizes)
    largest_group_size = max(comp_sizes.values()) if comp_sizes else 0
    largest_group_patches = max(comp_counts.values()) if comp_counts else 0

    stats = {
        "strategy": "largest_network",
        "corridors_used": len(selected),
        "budget_used": budget_used,
        "budget_total": budget_pixels,
        "patches_total": len(patch_sizes),
        "patches_connected": len(patch_sizes),
        "components_remaining": len(comp_sizes) if comp_sizes else len(patch_sizes),
        "largest_group_size": largest_group_size,
        "largest_group_patches": largest_group_patches,
        "total_connected_size": sum(comp_sizes.values()) if comp_sizes else 0,
        "redundant_links": redundant_links,
        "avg_degree": (2 * len(selected) / max(len(patch_sizes), 1)) if patch_sizes else 0,
    }
    _refresh_raster_component_stats(selected, patch_sizes, stats)
    return selected, stats


def optimize_largest_network(
    candidates: List[Dict],
    patch_sizes: Dict[int, int],
    budget_pixels: int,
    max_search_distance: int = 0,
) -> Tuple[Dict[int, Dict], Dict]:
    """
    Grow corridors outward from the largest patch, preferring highest gain per cost.

    This is intentionally different from MST-by-cost: it prioritizes "micro-gaps" that
    add a large patch/component to the growing network, instead of spending budget on
    very cheap links to tiny patches.
    """
    if not patch_sizes:
        return {}, {}
    seed = max(patch_sizes.items(), key=lambda kv: kv[1])[0]
    visited: Set[int] = {seed}
    visited_size = int(patch_sizes.get(int(seed), 0) or 0)
    remaining = max(0, int(budget_pixels))
    budget_used = 0
    selected: Dict[int, Dict] = {}
    selected_pixels: List[object] = []
    redundant_links = 0
    backbone_graph_edges: List[Tuple[int, int, float]] = []

    adjacency: Dict[int, List[Dict]] = defaultdict(list)
    for cand in candidates:
        adjacency[cand.get("patch1")].append(cand)
        adjacency[cand.get("patch2")].append(cand)

    bridge_by_mid: Dict[int, List[Dict]] = defaultdict(list)
    for cand in candidates:
        p1, p2 = cand.get("patch1"), cand.get("patch2")
        if p1 is None or p2 is None:
            continue
        bridge_by_mid[int(p1)].append(cand)
        bridge_by_mid[int(p2)].append(cand)

    # Heap entry: (-roi, -gain, cost, counter, cand)
    heap: List[Tuple[float, int, int, int, Dict]] = []
    counter = 0
    for cand in adjacency.get(seed, []):
        cost = int(cand.get("area", 0) or 0)
        if cost <= 0:
            continue
        p1, p2 = cand.get("patch1"), cand.get("patch2")
        if p1 is None or p2 is None:
            continue
        in1, in2 = p1 in visited, p2 in visited
        if not (in1 ^ in2):
            continue
        cand_pixels = cand.get("buffered_pixels") or cand.get("pixels", set())
        if not _redundancy_far_enough_pixels(
            list(cand_pixels),
            selected_pixels,
            int(max_search_distance or 0),
        ):
            continue
        new_node = int(p2) if in1 else int(p1)
        gain = int(patch_sizes.get(new_node, 0) or 0)
        if gain <= 0:
            continue
        roi = gain / float(cost)
        counter += 1
        heapq.heappush(heap, (-roi, -gain, cost, counter, cand))

    def _best_bridge_pair() -> Optional[Tuple[Dict, Dict, int]]:
        best = None
        best_score = 0.0
        for mid, edges in bridge_by_mid.items():
            if int(mid) in visited:
                continue
            best_to_visited = None
            best_to_visited_cost = 0
            best_to_unvisited = None
            best_to_unvisited_score = 0.0
            best_to_unvisited_other = None
            for cand in edges:
                p1, p2 = cand.get("patch1"), cand.get("patch2")
                if p1 is None or p2 is None:
                    continue
                a = int(p1)
                b = int(p2)
                other = b if a == mid else a
                cost = int(cand.get("area", 0) or 0)
                if cost <= 0:
                    continue
                if other in visited:
                    if best_to_visited is None or cost < best_to_visited_cost:
                        best_to_visited = cand
                        best_to_visited_cost = cost
                else:
                    gain_other = int(patch_sizes.get(int(other), 0) or 0)
                    if gain_other <= 0:
                        continue
                    score = gain_other / float(cost)
                    if best_to_unvisited is None or score > best_to_unvisited_score:
                        best_to_unvisited = cand
                        best_to_unvisited_score = score
                        best_to_unvisited_other = int(other)
            if best_to_visited is None or best_to_unvisited is None or best_to_unvisited_other is None:
                continue
            total_cost = int(best_to_visited_cost + int(best_to_unvisited.get("area", 0) or 0))
            if total_cost <= 0 or total_cost > remaining:
                continue
            gain_mid = int(patch_sizes.get(int(mid), 0) or 0)
            gain_other = int(patch_sizes.get(int(best_to_unvisited_other), 0) or 0)
            score = (gain_mid + gain_other) / float(total_cost) if total_cost else 0.0
            if score > best_score:
                best_score = score
                best = (best_to_visited, best_to_unvisited, int(best_to_unvisited_other))
        return best

    while remaining > 0:
        bridge_pick = _best_bridge_pair()
        if bridge_pick is None:
            break
        cand1, cand2, other_node = bridge_pick
        cost1 = int(cand1.get("area", 0) or 0)
        cost2 = int(cand2.get("area", 0) or 0)
        if cost1 <= 0 or cost2 <= 0 or (cost1 + cost2) > remaining:
            break
        p1, p2 = cand1.get("patch1"), cand1.get("patch2")
        if p1 is None or p2 is None:
            break
        mid_node = int(p1) if int(p1) not in visited else int(p2)
        gain_mid = int(patch_sizes.get(int(mid_node), 0) or 0)

        remaining -= cost1
        budget_used += cost1
        visited.add(int(mid_node))
        if gain_mid > 0:
            visited_size += gain_mid
        backbone_graph_edges.append((int(p1), int(p2), float(cand1.get("length", cost1) or cost1)))
        cid = len(selected) + 1
        selected[cid] = {
            "pixels": cand1.get("pixels", set()),
            "habitat_pixels": cand1.get("habitat_pixels", frozenset()),
            "patch_ids": {int(p1), int(p2)},
            "area": cost1,
            "connected_size": visited_size,
            "buffered_pixels": cand1.get("buffered_pixels", frozenset()),
            "type": "backbone",
        }
        selected_pixels.append(cand1.get("buffered_pixels") or cand1.get("pixels", set()))

        p1b, p2b = cand2.get("patch1"), cand2.get("patch2")
        if p1b is None or p2b is None:
            break
        remaining -= cost2
        budget_used += cost2
        visited.add(int(other_node))
        gain_other = int(patch_sizes.get(int(other_node), 0) or 0)
        if gain_other > 0:
            visited_size += gain_other
        backbone_graph_edges.append((int(p1b), int(p2b), float(cand2.get("length", cost2) or cost2)))
        cid = len(selected) + 1
        selected[cid] = {
            "pixels": cand2.get("pixels", set()),
            "habitat_pixels": cand2.get("habitat_pixels", frozenset()),
            "patch_ids": {int(p1b), int(p2b)},
            "area": cost2,
            "connected_size": visited_size,
            "buffered_pixels": cand2.get("buffered_pixels", frozenset()),
            "type": "backbone",
        }
        selected_pixels.append(cand2.get("buffered_pixels") or cand2.get("pixels", set()))

        for new_node in (int(mid_node), int(other_node)):
            for nxt in adjacency.get(new_node, []):
                n1, n2 = nxt.get("patch1"), nxt.get("patch2")
                if n1 is None or n2 is None:
                    continue
                touches = (n1 in visited) ^ (n2 in visited)
                if not touches:
                    continue
                cost2 = int(nxt.get("area", 0) or 0)
                if cost2 <= 0:
                    continue
                new2 = int(n2) if (n1 in visited) else int(n1)
                if new2 in visited:
                    continue
                gain2 = int(patch_sizes.get(new2, 0) or 0)
                if gain2 <= 0:
                    continue
                roi2 = gain2 / float(cost2)
                counter += 1
                heapq.heappush(heap, (-roi2, -gain2, cost2, counter, nxt))

    while heap and remaining > 0:
        _neg_roi, _neg_gain, cost, _, cand = heapq.heappop(heap)
        if cost <= 0 or cost > remaining:
            continue
        p1, p2 = cand.get("patch1"), cand.get("patch2")
        in1, in2 = p1 in visited, p2 in visited
        if in1 and in2:
            continue
        if not in1 and not in2:
            continue

        new_node = int(p2) if in1 else int(p1)
        if new_node in visited:
            continue
        gain = int(patch_sizes.get(int(new_node), 0) or 0)
        if gain <= 0:
            continue
        visited.add(new_node)
        visited_size += gain
        remaining -= cost
        budget_used += cost
        comp_size = visited_size
        backbone_graph_edges.append((int(p1), int(p2), float(cand.get("length", cost) or cost)))

        cid = len(selected) + 1
        selected[cid] = {
            "pixels": cand.get("pixels", set()),
            "habitat_pixels": cand.get("habitat_pixels", frozenset()),
            "patch_ids": {p1, p2},
            "area": cost,
            "connected_size": comp_size,
            "buffered_pixels": cand.get("buffered_pixels", frozenset()),
            "type": "backbone",
        }
        selected_pixels.append(cand.get("buffered_pixels") or cand.get("pixels", set()))

        for nxt in adjacency.get(new_node, []):
            n1, n2 = nxt.get("patch1"), nxt.get("patch2")
            if n1 is None or n2 is None:
                continue
            touches = (n1 in visited) ^ (n2 in visited)
            if not touches:
                continue
            cost2 = int(nxt.get("area", 0) or 0)
            if cost2 <= 0:
                continue
            new2 = int(n2) if (n1 in visited) else int(n1)
            if new2 in visited:
                continue
            gain2 = int(patch_sizes.get(new2, 0) or 0)
            if gain2 <= 0:
                continue
            roi2 = gain2 / float(cost2)
            counter += 1
            heapq.heappush(heap, (-roi2, -gain2, cost2, counter, nxt))

    if remaining > 0 and selected and visited:
        selected_pairs = {tuple(sorted(map(int, c.get("patch_ids", [])))) for c in selected.values()}

        if nx is not None and graph_math is not None:
            G_backbone = nx.Graph()
            for pid in visited:
                G_backbone.add_node(int(pid))
            for a, b, w in backbone_graph_edges:
                if int(a) in visited and int(b) in visited:
                    G_backbone.add_edge(int(a), int(b), weight=float(w))

            loop_candidates: List[Tuple[float, Dict]] = []
            for cand in candidates:
                p1, p2 = cand.get("patch1"), cand.get("patch2")
                if p1 is None or p2 is None:
                    continue
                p1i, p2i = int(p1), int(p2)
                if p1i not in visited or p2i not in visited:
                    continue
                pair = tuple(sorted((p1i, p2i)))
                if pair in selected_pairs:
                    continue
                cost = int(cand.get("area", 0))
                if cost <= 0 or cost > remaining:
                    continue
                cand_pixels = cand.get("buffered_pixels") or cand.get("pixels", set())
                if not _redundancy_far_enough_pixels(
                    list(cand_pixels),
                    selected_pixels,
                    int(max_search_distance or 0),
                ):
                    continue
                length = float(cand.get("length", 0) or 0.0)
                if length <= 0:
                    length = float(cost)
                try:
                    if not nx.has_path(G_backbone, p1i, p2i):
                        continue
                    score = float(graph_math.score_edge_for_loops(G_backbone, p1i, p2i, length))
                except Exception:
                    continue
                efficiency = score / float(cost) if cost else 0.0
                loop_candidates.append((efficiency, cand))

            loop_candidates.sort(key=lambda t: t[0], reverse=True)
            for _eff, cand in loop_candidates:
                cost = int(cand.get("area", 0))
                if cost <= 0 or cost > remaining:
                    continue
                p1, p2 = int(cand.get("patch1")), int(cand.get("patch2"))
                pair = tuple(sorted((p1, p2)))
                if pair in selected_pairs:
                    continue
                remaining -= cost
                budget_used += cost
                redundant_links += 1
                selected_pairs.add(pair)
                length = float(cand.get("length", 0) or 0.0)
                if length <= 0:
                    length = float(cost)
                try:
                    G_backbone.add_edge(p1, p2, weight=length)
                except Exception:
                    pass

                cid = len(selected) + 1
                selected[cid] = {
                    "pixels": cand.get("pixels", set()),
                    "habitat_pixels": cand.get("habitat_pixels", frozenset()),
                    "patch_ids": {p1, p2},
                    "area": cost,
                    "connected_size": visited_size,
                    "buffered_pixels": cand.get("buffered_pixels", frozenset()),
                    "type": "strategic_loop",
                }
                selected_pixels.append(cand.get("buffered_pixels") or cand.get("pixels", set()))
        else:
            for cand in sorted(candidates, key=lambda c: c.get("area", 0)):
                if remaining <= 0:
                    break
                p1, p2 = cand.get("patch1"), cand.get("patch2")
                if p1 is None or p2 is None:
                    continue
                if int(p1) not in visited or int(p2) not in visited:
                    continue
                pair = tuple(sorted((int(p1), int(p2))))
                if pair in selected_pairs:
                    continue
                cost = int(cand.get("area", 0))
                if cost <= 0 or cost > remaining:
                    continue
                cand_pixels = cand.get("buffered_pixels") or cand.get("pixels", set())
                if not _redundancy_far_enough_pixels(
                    list(cand_pixels),
                    selected_pixels,
                    int(max_search_distance or 0),
                ):
                    continue
                remaining -= cost
                budget_used += cost
                redundant_links += 1
                selected_pairs.add(pair)
                cid = len(selected) + 1
                selected[cid] = {
                    "pixels": cand.get("pixels", set()),
                    "habitat_pixels": cand.get("habitat_pixels", frozenset()),
                    "patch_ids": {p1, p2},
                    "area": cost,
                    "connected_size": visited_size,
                    "buffered_pixels": cand.get("buffered_pixels", frozenset()),
                    "type": "redundant",
                }
                selected_pixels.append(cand.get("buffered_pixels") or cand.get("pixels", set()))

    comp_sizes, comp_counts = _build_corridor_stats(selected, patch_sizes)
    largest_group_size = max(comp_sizes.values()) if comp_sizes else 0
    largest_group_patches = max(comp_counts.values()) if comp_counts else 0

    stats = {
        "strategy": "largest_network",
        "corridors_used": len(selected),
        "budget_used": budget_used,
        "budget_total": budget_pixels,
        "patches_total": len(patch_sizes),
        "patches_connected": len(visited),
        "components_remaining": len(comp_sizes) if comp_sizes else len(patch_sizes),
        "largest_group_size": largest_group_size,
        "largest_group_patches": largest_group_patches,
        "total_connected_size": sum(comp_sizes.values()) if comp_sizes else 0,
        "redundant_links": redundant_links,
        "avg_degree": (2 * len(selected) / max(len(patch_sizes), 1)) if patch_sizes else 0,
    }
    _refresh_raster_component_stats(selected, patch_sizes, stats)
    return selected, stats


def optimize_circuit_utility(
    candidates: List[Dict],
    patch_sizes: Dict[int, int],
    budget_pixels: int,
    max_search_distance: int = 0,
    overlap_reject_ratio: float = 0.30,
) -> Tuple[Dict[int, Dict], Dict]:
    """
    Most Connectivity (Utility) strategy for raster candidates.

    Uses a weighted greedy loop over corridor candidates:
      score = (sqrt(size_p1 * size_p2) / cost) * multiplier

    multiplier:
      1.0  if candidate connects different components
      0.5  if already connected but spatially distinct
      0.01 if overlaps existing corridor for the pair by > overlap_reject_ratio
    """

    from .core_select import select_circuit_utility

    def _pair_key(a: int, b: int) -> Tuple[int, int]:
        return (a, b) if a <= b else (b, a)

    # Proximity-aware expansion: treat "parallel adjacent" corridors as overlap too.
    prox_radius = 1
    try:
        for cand in candidates:
            pr = int(cand.get("_prox_radius", 0) or 0)
            if pr > prox_radius:
                prox_radius = pr
    except Exception:
        prox_radius = 1
    prox_radius = max(1, int(prox_radius))
    expand_offsets = [(dr, dc) for dr in range(-prox_radius, prox_radius + 1) for dc in range(-prox_radius, prox_radius + 1)]

    def _expand_pixels(pixels: frozenset) -> frozenset:
        expanded: Set[Tuple[int, int]] = set()
        for r, c in pixels:
            for dr, dc in expand_offsets:
                expanded.add((r + dr, c + dc))
        return frozenset(expanded)

    def _get_patch_ids(cand: Dict) -> List[int]:
        p1, p2 = cand.get("patch1"), cand.get("patch2")
        if p1 is None or p2 is None:
            return []
        return [int(p1), int(p2)]

    def _get_pair(cand: Dict) -> Tuple[int, int]:
        p1, p2 = cand.get("patch1"), cand.get("patch2")
        if p1 is None or p2 is None:
            return (0, 0)
        return _pair_key(int(p1), int(p2))

    def _get_cost(cand: Dict) -> float:
        try:
            return float(int(cand.get("area", 0) or 0))
        except Exception:
            return 0.0

    def _get_length(cand: Dict) -> float:
        try:
            length = float(cand.get("length", 0) or 0.0)
        except Exception:
            length = 0.0
        if length <= 0.0:
            try:
                length = float(_get_cost(cand) or 0.0)
            except Exception:
                length = 0.0
        return length

    def _get_patch_size(pid: int) -> float:
        try:
            return float(patch_sizes.get(int(pid), 0) or 0.0)
        except Exception:
            return 0.0

    def _get_base_roi(cand: Dict) -> float:
        try:
            pids = _get_patch_ids(cand)
            if len(pids) < 2:
                return 0.0
            s1 = int(patch_sizes.get(int(pids[0]), 0) or 0)
            s2 = int(patch_sizes.get(int(pids[1]), 0) or 0)
            cost = float(int(cand.get("area", 0) or 0))
            if s1 <= 0 or s2 <= 0 or cost <= 0:
                return 0.0
            return float(math.sqrt(float(s1) * float(s2)) / cost)
        except Exception:
            return 0.0

    def _overlap_obj(cand: Dict) -> object:
        new_buf = cand.get("buffered_pixels") or frozenset(cand.get("pixels", []))
        try:
            return _expand_pixels(frozenset(new_buf))
        except Exception:
            return frozenset()

    def _overlap_ratio(cand: Dict, prior: Sequence[object]) -> float:
        new_buf = cand.get("buffered_pixels") or frozenset(cand.get("pixels", []))
        new_buf = frozenset(new_buf)
        denom = float(len(new_buf))
        if denom <= 0 or not prior:
            return 0.0
        best = 0.0
        for prev in prior:
            try:
                ratio = len(new_buf.intersection(prev)) / denom
            except Exception:
                continue
            if ratio > best:
                best = ratio
        return best

    def _redundancy_distance_ok(cand: Dict, prior: Sequence[object]) -> bool:
        threshold = int(max_search_distance or 0)
        if threshold <= 0:
            return True
        cand_pixels = cand.get("buffered_pixels") or cand.get("pixels", set())
        if not cand_pixels:
            return True
        return _redundancy_far_enough_pixels(list(cand_pixels), prior, threshold)

    picks, base_stats = select_circuit_utility(
        candidates,
        budget=float(max(0, int(budget_pixels))),
        get_patch_ids=_get_patch_ids,
        get_pair_key=_get_pair,
        get_cost=_get_cost,
        get_base_roi=_get_base_roi,
        get_length=_get_length,
        get_patch_size=_get_patch_size,
        overlap_ratio=_overlap_ratio,
        overlap_obj=_overlap_obj,
        redundancy_distance_ok=_redundancy_distance_ok,
        overlap_reject_ratio=float(overlap_reject_ratio),
        max_prior_per_pair=3,
        diminishing_base=0.5,
    )

    selected: Dict[int, Dict] = {}

    uf = UnionFind()
    for pid, size in patch_sizes.items():
        uf.find(int(pid))
        uf.size[int(pid)] = int(size)
        uf.count[int(pid)] = 1

    for pick in picks:
        cand = pick.candidate
        pids = _get_patch_ids(cand)
        if len(pids) < 2:
            continue
        p1, p2 = int(pids[0]), int(pids[1])
        cost = int(_get_cost(cand) or 0.0)
        if cost <= 0:
            continue
        uf.union(p1, p2)
        comp_root = int(uf.find(p1))
        conn_size = int(uf.size.get(comp_root, 0) or 0)
        cid = len(selected) + 1
        selected[cid] = {
            "pixels": cand.get("pixels", set()),
            "habitat_pixels": cand.get("habitat_pixels", frozenset()),
            "patch_ids": {p1, p2},
            "area": cost,
            "connected_size": conn_size,
            "buffered_pixels": cand.get("buffered_pixels", frozenset()),
            "type": pick.corr_type,
            "utility_score": float(pick.score),
            "overlap_ratio": float(pick.overlap_ratio),
        }

    comp_sizes, comp_counts = _build_corridor_stats(selected, patch_sizes)
    for corr in selected.values():
        pids = list(corr.get("patch_ids", []))
        if not pids:
            continue
        root = uf.find(int(pids[0]))
        corr["connected_size"] = int(comp_sizes.get(root, corr.get("connected_size", 0)) or 0)

    largest_group_size = max(comp_sizes.values()) if comp_sizes else 0
    largest_group_patches = max(comp_counts.values()) if comp_counts else 0

    stats = {
        "strategy": "circuit_utility",
        "corridors_used": len(selected),
        "budget_used": int(float(base_stats.get("budget_used", 0.0) or 0.0)),
        "budget_total": budget_pixels,
        "patches_total": len(patch_sizes),
        "patches_connected": largest_group_patches,
        "components_remaining": len(comp_sizes) if comp_sizes else len(patch_sizes),
        "largest_group_size": largest_group_size,
        "largest_group_patches": largest_group_patches,
        "primary_links": int(base_stats.get("primary_links", 0) or 0),
        "redundant_links": int(base_stats.get("redundant_links", 0) or 0),
        "wasteful_links": int(base_stats.get("wasteful_links", 0) or 0),
    }
    return selected, stats


def _refresh_raster_component_stats(
    corridors: Dict[int, Dict],
    patch_sizes: Dict[int, int],
    stats: Dict,
) -> None:
    """Update connectivity fields in stats after adding corridors."""
    uf = UnionFind()
    for pid, size in patch_sizes.items():
        uf.find(int(pid))
        uf.size[int(pid)] = int(size)
        uf.count[int(pid)] = 1
    for corr in corridors.values():
        pids = list(corr.get("patch_ids", []))
        if len(pids) >= 2:
            uf.union(int(pids[0]), int(pids[1]))

    comp_sizes: Dict[int, int] = {}
    comp_counts: Dict[int, int] = {}
    for pid in patch_sizes.keys():
        root = uf.find(int(pid))
        comp_sizes[root] = comp_sizes.get(root, 0) + int(patch_sizes.get(pid, 0) or 0)
        comp_counts[root] = comp_counts.get(root, 0) + 1

    for corr in corridors.values():
        pids = list(corr.get("patch_ids", []))
        if not pids:
            continue
        root = uf.find(int(pids[0]))
        comp_sizes[root] = comp_sizes.get(root, 0) + int(corr.get("area", 0) or 0)

    for corr in corridors.values():
        pids = list(corr.get("patch_ids", []))
        if not pids:
            continue
        root = uf.find(int(pids[0]))
        corr["connected_size"] = int(comp_sizes.get(root, corr.get("connected_size", 0)) or 0)
    if comp_sizes:
        stats["largest_group_size"] = max(comp_sizes.values())
        stats["total_connected_size"] = sum(comp_sizes.values())
    if comp_counts:
        stats["largest_group_patches"] = max(comp_counts.values())
        stats["patches_connected"] = stats.get("largest_group_patches", 0)
        stats["components_remaining"] = len(comp_counts)
    elif patch_sizes:
        stats["components_remaining"] = len(patch_sizes)


def _min_unselected_candidate_cost(candidates: List[Dict], corridors: Dict[int, Dict]) -> Optional[int]:
    """Return the minimum corridor cost among candidates not already selected."""
    if not candidates:
        return None
    selected_pairs: Set[Tuple[int, int]] = set()
    for data in corridors.values():
        pids = list(data.get("patch_ids", []))
        if len(pids) >= 2:
            selected_pairs.add(tuple(sorted((int(pids[0]), int(pids[1])))))
    min_cost: Optional[int] = None
    for cand in candidates:
        try:
            p1 = int(cand.get("patch1", 0) or 0)
            p2 = int(cand.get("patch2", 0) or 0)
        except Exception:
            continue
        if p1 <= 0 or p2 <= 0:
            continue
        pair = (p1, p2) if p1 <= p2 else (p2, p1)
        if pair in selected_pairs:
            continue
        cost = int(cand.get("area", 0) or 0)
        if cost <= 0:
            continue
        if min_cost is None or cost < min_cost:
            min_cost = cost
    return min_cost


def _apply_hybrid_leftover_budget_raster(
    candidates: List[Dict],
    corridors: Dict[int, Dict],
    patch_sizes: Dict[int, int],
    remaining_budget: int,
    roi_bias: float = ROI_REDUNDANCY_BIAS,
    max_search_distance: int = 0,
) -> Tuple[int, int, int]:
    """
    Spend remaining budget by comparing low-value connections vs redundancy.
    Returns (budget_used, low_value_added, redundancy_added).
    """
    if remaining_budget <= 0 or not candidates:
        return 0, 0, 0

    uf = UnionFind()
    for pid, size in patch_sizes.items():
        uf.find(int(pid))
        uf.size[int(pid)] = int(size)
        uf.count[int(pid)] = 1
    for data in corridors.values():
        pids = list(data.get("patch_ids", []))
        if len(pids) >= 2:
            uf.union(int(pids[0]), int(pids[1]))

    roots = {int(uf.find(int(pid))) for pid in patch_sizes}
    if not roots:
        return 0, 0, 0
    main_root = max(roots, key=lambda r: uf.size.get(r, 0))

    selected_pairs: Set[Tuple[int, int]] = set()
    selected_pixels: List[object] = []
    for data in corridors.values():
        pids = list(data.get("patch_ids", []))
        if len(pids) >= 2:
            selected_pairs.add(tuple(sorted((int(pids[0]), int(pids[1])))))
        pix = data.get("buffered_pixels") or data.get("pixels", set())
        if pix:
            selected_pixels.append(pix)

    G = None
    if nx is not None and graph_math is not None:
        G = nx.Graph()
        for pid in patch_sizes:
            G.add_node(int(pid))
        for data in corridors.values():
            pids = list(data.get("patch_ids", []))
            if len(pids) < 2:
                continue
            length = float(data.get("length", data.get("area", 1) or 1) or 1)
            G.add_edge(int(pids[0]), int(pids[1]), weight=length)

    budget_used = 0
    low_value_added = 0
    redundancy_added = 0

    while remaining_budget > 0:
        best_new = None
        best_red = None

        for cand in candidates:
            cost = int(cand.get("area", 0) or 0)
            if cost <= 0 or cost > remaining_budget:
                continue
            p1 = cand.get("patch1")
            p2 = cand.get("patch2")
            if p1 is None or p2 is None:
                continue
            p1i, p2i = int(p1), int(p2)
            pair = (p1i, p2i) if p1i <= p2i else (p2i, p1i)
            if pair in selected_pairs:
                continue
            r1, r2 = int(uf.find(p1i)), int(uf.find(p2i))

            if r1 != r2:
                if (r1 == main_root) ^ (r2 == main_root):
                    gain = int(uf.size.get(r2 if r1 == main_root else r1, 0) or 0)
                    roi_new = (gain / cost) if cost else 0.0
                    if roi_new > 0 and (best_new is None or roi_new > best_new[0]):
                        best_new = (roi_new, cand, gain)
                continue

            if r1 != main_root:
                continue
            length = float(cand.get("length", cost) or cost)
            score = 0.0
            if G is not None and graph_math is not None:
                try:
                    if nx.has_path(G, p1i, p2i):
                        score = float(graph_math.score_edge_for_loops(G, p1i, p2i, length))
                except Exception:
                    score = 0.0
            if score <= 0.0:
                continue
            cand_pixels = cand.get("buffered_pixels") or cand.get("pixels", set())
            if not _redundancy_far_enough_pixels(
                list(cand_pixels),
                selected_pixels,
                int(max_search_distance or 0),
            ):
                continue
            roi_red = score / cost if cost else 0.0
            if roi_red > 0 and (best_red is None or roi_red > best_red[0]):
                best_red = (roi_red, cand, score)

        if best_new is None and best_red is None:
            break

        pick_red = False
        if best_red is not None and best_new is not None:
            pick_red = best_red[0] > best_new[0] * (1.0 + float(roi_bias))
        elif best_red is not None:
            pick_red = True

        roi, cand, _extra = best_red if pick_red else best_new
        if cand is None or roi <= 0:
            break

        cost = int(cand.get("area", 0) or 0)
        if cost <= 0 or cost > remaining_budget:
            break
        p1i, p2i = int(cand.get("patch1")), int(cand.get("patch2"))
        pair = (p1i, p2i) if p1i <= p2i else (p2i, p1i)

        if pick_red:
            redundancy_added += 1
            corr_type = "redundant"
        else:
            low_value_added += 1
            corr_type = "low_value"

        r1 = int(uf.find(p1i))
        r2 = int(uf.find(p2i))
        if r1 != r2:
            if r1 == main_root:
                uf.union(p1i, p2i)
                main_root = int(uf.find(p1i))
            elif r2 == main_root:
                uf.union(p2i, p1i)
                main_root = int(uf.find(p2i))
            else:
                uf.union(p1i, p2i)
                main_root = int(uf.find(p1i))

        comp_root = int(uf.find(p1i))
        conn_size = int(uf.size.get(comp_root, 0) or 0)
        cid = len(corridors) + 1
        corridors[cid] = {
            "pixels": cand.get("pixels", set()),
            "habitat_pixels": cand.get("habitat_pixels", frozenset()),
            "patch_ids": {p1i, p2i},
            "area": cost,
            "connected_size": conn_size,
            "buffered_pixels": cand.get("buffered_pixels", frozenset()),
            "length": float(cand.get("length", cost) or cost),
            "type": corr_type,
        }
        selected_pairs.add(pair)
        pix = cand.get("buffered_pixels") or cand.get("pixels", set())
        if pix:
            selected_pixels.append(pix)
        remaining_budget -= cost
        budget_used += cost
        if G is not None:
            G.add_edge(p1i, p2i, weight=float(cand.get("length", cost) or cost))

    return budget_used, low_value_added, redundancy_added


def _thicken_corridors_raster(
    corridors: Dict[int, Dict],
    remaining_budget: int,
    params: "RasterRunParams",
    patch_sizes: Dict[int, int],
    obstacle_mask: Optional[np.ndarray],
    labels_shape: Tuple[int, int],
    max_width_factor: int = 3,
) -> int:
    """
    Use remaining budget to thicken corridors up to max_width_factor (in steps of 1x).
    Prioritize corridors adjacent to the largest patch.
    """
    if remaining_budget <= 0 or not corridors or max_width_factor <= 1:
        return 0

    try:
        largest_patch = max(patch_sizes.items(), key=lambda kv: kv[1])[0]
    except Exception:
        largest_patch = None

    ordered = list(corridors.keys())
    if largest_patch is not None:
        ordered = [cid for cid in ordered if largest_patch in corridors[cid].get("patch_ids", set())] + [
            cid for cid in ordered if largest_patch not in corridors[cid].get("patch_ids", set())
        ]

    rows, cols = labels_shape
    budget_used = 0

    for cid in ordered:
        if remaining_budget <= 0:
            break
        cdata = corridors.get(cid) or {}
        base_pixels = cdata.get("_thicken_base_pixels")
        if base_pixels is None:
            base_pixels = frozenset(cdata.get("pixels", set()))
            cdata["_thicken_base_pixels"] = base_pixels
        if not base_pixels:
            continue

        current_buf = cdata.get("buffered_pixels")
        if current_buf is None:
            offsets = _corridor_offsets(params.min_corridor_width)
            current_buf = frozenset(
                _inflate_corridor_pixels(set(base_pixels), offsets, rows, cols, obstacle_mask=obstacle_mask)
            )
            cdata["buffered_pixels"] = current_buf

        current_factor = int(cdata.get("_thicken_factor", 1))
        for factor in range(current_factor + 1, max_width_factor + 1):
            width = max(1, int(params.min_corridor_width * factor))
            offsets = _corridor_offsets(width)
            new_buf = frozenset(
                _inflate_corridor_pixels(set(base_pixels), offsets, rows, cols, obstacle_mask=obstacle_mask)
            )
            added = len(new_buf) - len(current_buf)
            if added <= 0:
                cdata["_thicken_factor"] = factor
                current_buf = new_buf
                continue
            if added > remaining_budget:
                return budget_used
            remaining_budget -= added
            budget_used += added
            current_buf = new_buf
            cdata["buffered_pixels"] = current_buf
            cdata["area"] = len(current_buf)
            cdata["_thicken_factor"] = factor

    return budget_used


def _to_dataclass(params: Dict) -> RasterRunParams:
    """Convert raw parameter dict into the expected dataclass."""
    return RasterRunParams(
        patch_connectivity=int(params.get("patch_connectivity", 4)),
        patch_mode=str(params.get("patch_mode", "value")).lower(),
        patch_values=list(params.get("patch_values", [])),
        range_lower=params.get("range_lower"),
        range_upper=params.get("range_upper"),
        obstacle_enabled=bool(params.get("obstacle_enabled", False)),
        obstacle_mode=str(params.get("obstacle_mode", "value")).lower(),
        obstacle_values=list(params.get("obstacle_values", [])),
        obstacle_range_lower=params.get("obstacle_range_lower"),
        obstacle_range_upper=params.get("obstacle_range_upper"),
        value_tolerance=float(params.get("value_tolerance", 1e-6)),
        nodata_fallback=float(params.get("nodata_fallback", -9999)),
        min_patch_size=int(params.get("min_patch_size", 30)),
        budget_pixels=int(params.get("budget_pixels", 10)),
        max_search_distance=int(params.get("max_search_distance", 50)),
        max_corridor_area=None,
        min_corridor_width=int(params.get("min_corridor_width", 1)),
        allow_bottlenecks=bool(params.get("allow_bottlenecks", True)),
    )


def run_raster_analysis(
    layer: QgsRasterLayer,
    output_dir: str,
    raw_params: Dict,
    strategy: str = "circuit_utility",
    temporary: bool = False,
    iface=None,
    progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
) -> List[Dict]:
    """Execute the raster corridor analysis for the provided layer."""
    if not isinstance(layer, QgsRasterLayer) or not layer.isValid():
        raise RasterAnalysisError("Selected layer is not a valid raster layer.")

    params = _to_dataclass(raw_params)
    overall_start = time.time()
    src_path = layer.source()
    summary_dir = output_dir or os.path.dirname(src_path) or os.getcwd()
    if not summary_dir:
        summary_dir = os.getcwd()
    os.makedirs(summary_dir, exist_ok=True)
    return _run_raster_analysis_core(
        layer,
        output_dir,
        raw_params,
        strategy,
        temporary,
        iface,
        progress_cb,
        params,
        overall_start,
    )


def _add_raster_run_summary_layer(
    *,
    input_layer_name: str,
    strategy: str,
    stats: Dict,
    out_path: str,
    corridor_path: Optional[str] = None,
    combined_path: Optional[str] = None,
) -> None:
    """Add a small in-project summary table for the raster run."""
    if QgsVectorLayer is None or QgsProject is None or QgsFeature is None or QgsField is None:
        return
    try:
        title = f"TerraLink Raster Summary ({input_layer_name})"
        summary_layer = QgsVectorLayer("NoGeometry", title, "memory")
        if not summary_layer.isValid():
            return
        provider = summary_layer.dataProvider()
        provider.addAttributes([QgsField("item", QVariant.String), QgsField("value", QVariant.String)])
        summary_layer.updateFields()

        def _add(item: str, value: object) -> None:
            feat = QgsFeature(summary_layer.fields())
            feat.setAttributes([str(item), "" if value is None else str(value)])
            provider.addFeature(feat)

        _add("Input layer", input_layer_name)
        _add("Strategy", str(strategy))
        rows = stats.get("raster_rows")
        cols = stats.get("raster_cols")
        total_px = stats.get("raster_pixels_total")
        if rows is not None and cols is not None:
            if total_px is None:
                size_val = f"{int(rows):,} x {int(cols):,}"
            else:
                size_val = f"{int(rows):,} x {int(cols):,} ({int(total_px):,})"
            _add("Raster size (px)", size_val)
        raw_patches = stats.get("patches_total_raw")
        filtered_patches = stats.get("patches_total_filtered")
        if raw_patches is not None or filtered_patches is not None:
            raw_display = "" if raw_patches is None else f"{int(raw_patches):,}"
            filt_display = "" if filtered_patches is None else f"{int(filtered_patches):,}"
            patches_val = raw_display or filt_display
            if raw_display and filt_display:
                patches_val = f"{raw_display} / {filt_display}"
            _add("Patches (raw / filtered)", patches_val)
        area_label = stats.get("raster_area_label", "px")
        dist_label = stats.get("raster_dist_label", "px")

        def _fmt(value: object) -> str:
            try:
                return f"{float(value):.2f}"
            except Exception:
                return "" if value is None else str(value)

        min_patch_disp = stats.get("min_patch_size_display", stats.get("min_patch_size"))
        if min_patch_disp is not None:
            _add(f"Min patch size ({area_label})", _fmt(min_patch_disp) if area_label != "px" else min_patch_disp)
        max_search_disp = stats.get("max_search_distance_display", stats.get("max_search_distance"))
        if max_search_disp is not None:
            _add(
                f"Max search distance ({dist_label})",
                _fmt(max_search_disp) if dist_label != "px" else max_search_disp,
            )
        min_width_disp = stats.get("min_corridor_width_display", stats.get("min_corridor_width"))
        if min_width_disp is not None:
            _add(
                f"Min corridor width ({dist_label})",
                _fmt(min_width_disp) if dist_label != "px" else min_width_disp,
            )
        candidate_pairs = stats.get("candidate_pairs")
        possible_pairs = stats.get("possible_pairs")
        if candidate_pairs is not None:
            if possible_pairs:
                _add("Candidate pairs", f"{int(candidate_pairs):,} / {int(possible_pairs):,}")
            else:
                _add("Candidate pairs", f"{int(candidate_pairs):,}")
        if stats.get("candidate_corridors") is not None:
            _add("Candidate corridors", stats.get("candidate_corridors"))
        candidate_time = stats.get("candidate_search_s")
        if candidate_time is not None:
            _add("Candidate search time (s)", f"{float(candidate_time):.2f}")
        _add("Corridors created", stats.get("corridors_used", 0))
        budget_used_disp = stats.get("budget_used_display", stats.get("budget_used", 0))
        budget_total_disp = stats.get("budget_total_display", stats.get("budget_total", 0))
        _add(
            f"Budget used ({area_label})",
            _fmt(budget_used_disp) if area_label != "px" else budget_used_disp,
        )
        _add(
            f"Budget total ({area_label})",
            _fmt(budget_total_disp) if area_label != "px" else budget_total_disp,
        )
        _add("Contiguous areas raster", out_path or "")
        if corridor_path:
            _add("Corridor raster", corridor_path)
        if combined_path:
            _add("Combined raster", combined_path)

        QgsProject.instance().addMapLayer(summary_layer)
    except Exception:
        return


def _perform_landscape_analysis(
    arr: np.ndarray,
    layer_name: str,
    res_x: float,
    res_y: float,
    is_metric: bool,
    params: Optional[RasterRunParams] = None,
    pre_arr: Optional[np.ndarray] = None,
) -> List[str]:
    """Calculate landscape metrics for the output raster (contiguous networks)."""
    if not HAS_NDIMAGE or ndimage is None:
        return ["Landscape analysis skipped: scipy.ndimage not found."]
    results: List[str] = []

    if is_metric:
        area_unit, dist_unit = "ha", "km"
        pixel_area_ha = (res_x * res_y) / 10000.0
        dist_factor = 1000.0
        gyrate_unit = "m"
    else:
        area_unit, dist_unit = "pixels", "pixels"
        pixel_area_ha, dist_factor = 1.0, 1.0
        gyrate_unit = "pixels"

    def _metrics_from_mask(mask: np.ndarray) -> Tuple[Dict[str, float], bool, np.ndarray, int]:
        hab_pixels = int(np.sum(mask))
        if hab_pixels <= 0:
            empty_labels = np.zeros_like(mask, dtype=int)
            return {
                "total_area": 0.0,
                "num_patches": 0.0,
                "mesh": 0.0,
                "lps": 0.0,
                "split": 0.0,
                "pd": 0.0,
                "lpi": 0.0,
                "total_edge": 0.0,
                "msi": 0.0,
                "ed": 0.0,
                "mean_para": 0.0,
                "frac": 0.0,
                "ai": 0.0,
                "cai": 0.0,
                "mean_gyrate": 0.0,
            }, True, empty_labels, 0

        s = ndimage.generate_binary_structure(2, 2)
        labeled_array, num_patches = ndimage.label(mask, structure=s)

        patch_ids = np.arange(1, num_patches + 1)
        patch_pixel_counts = np.bincount(labeled_array.ravel())[1:] if num_patches > 0 else np.array([], dtype=float)

        kernel = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        neighbor_count = ndimage.convolve(mask.astype(float), kernel, mode="constant", cval=0)
        edge_map = (4.0 - neighbor_count) * mask

        patch_areas_ha = patch_pixel_counts * pixel_area_ha
        patch_perims_m = ndimage.sum(edge_map, labeled_array, index=patch_ids) * (res_x if is_metric else 1.0) if num_patches > 0 else np.array([], dtype=float)
        patch_areas_m2 = patch_pixel_counts * (res_x * res_y if is_metric else 1.0)
        total_area_ha = hab_pixels * pixel_area_ha
        total_edge_m = float(np.sum(edge_map)) * (res_x if is_metric else 1.0)

        pd = (num_patches / total_area_ha) * 100 if total_area_ha > 0 else 0.0
        mesh = float(np.sum(patch_areas_ha**2) / total_area_ha) if total_area_ha > 0 else 0.0
        split = float((total_area_ha**2) / np.sum(patch_areas_ha**2)) if np.sum(patch_areas_ha**2) > 0 else 0.0
        lps_ha = float(np.max(patch_areas_ha)) if num_patches > 0 else 0.0
        lpi = (lps_ha / total_area_ha) * 100 if total_area_ha > 0 else 0.0

        denom_ai = (2 * hab_pixels - 2 * int(np.sqrt(hab_pixels)) - 2)
        ai = ((hab_pixels * 4 - np.sum(edge_map)) / 2) / denom_ai * 100 if denom_ai > 0 else 0.0

        ed = total_edge_m / total_area_ha if total_area_ha > 0 else 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            msi = float(np.mean(0.25 * patch_perims_m / np.sqrt(patch_areas_m2))) if num_patches > 0 else 0.0
            mean_para = float(np.mean(patch_perims_m / patch_areas_m2)) if num_patches > 0 else 0.0
            frac_vals = (2 * np.log(0.25 * patch_perims_m)) / np.log(patch_areas_m2) if num_patches > 0 else np.array([], dtype=float)
            frac = float(np.nanmean(frac_vals)) if num_patches > 0 and frac_vals.size else 0.0

        core_mask = ndimage.binary_erosion(mask)
        cai = (float(np.sum(core_mask)) / hab_pixels) * 100 if hab_pixels > 0 else 0.0

        gyrations = []
        for p_id in patch_ids:
            if num_patches > 1000 and p_id > 200:
                break
            coords = np.argwhere(labeled_array == p_id)
            if len(coords) < 2:
                continue
            gyrations.append(
                np.sqrt(np.mean(np.sum((coords - coords.mean(axis=0)) ** 2, axis=1))) * (res_x if is_metric else 1.0)
            )
        mean_gyrate = float(np.mean(gyrations)) if gyrations else 0.0

        return {
            "total_area": float(total_area_ha),
            "num_patches": float(num_patches),
            "mesh": float(mesh),
            "lps": float(lps_ha),
            "split": float(split),
            "pd": float(pd),
            "lpi": float(lpi),
            "total_edge": float(total_edge_m / dist_factor),
            "msi": float(msi),
            "ed": float(ed),
            "mean_para": float(mean_para),
            "frac": float(frac),
            "ai": float(ai),
            "cai": float(cai),
            "mean_gyrate": float(mean_gyrate),
        }, False, labeled_array, int(num_patches)

    pre_mask = (pre_arr if pre_arr is not None else arr) > 0
    post_mask = arr > 0

    pre_metrics, pre_empty, pre_labels, _ = _metrics_from_mask(pre_mask)
    post_metrics, post_empty, post_labels, post_count = _metrics_from_mask(post_mask)

    connected_groups = 0
    isolated_groups = 0
    if post_count > 0 and np.any(pre_labels):
        for group_id in range(1, post_count + 1):
            overlap = np.unique(pre_labels[post_labels == group_id])
            overlap = overlap[overlap > 0]
            if overlap.size == 1:
                isolated_groups += 1
            elif overlap.size >= 2:
                connected_groups += 1

    header = "=" * 100
    results.append(header)
    results.append(f" LANDSCAPE ANALYSIS: {layer_name}")
    results.append(header)
    if pre_empty:
        results.append("Note: pre-corridor habitat mask contains no habitat pixels.")
    if post_empty:
        results.append("Note: post-corridor habitat mask contains no habitat pixels.")
    results.append(f"{'METRIC NAME':<30} | {'PRE':<15} | {'POST':<15} | {'INTERPRETATION'}")
    results.append("-" * 100)

    pre_np = f"{pre_metrics.get('num_patches', 0.0):,.0f}"
    post_np = f"{post_metrics.get('num_patches', 0.0):,.0f}"
    results.append(
        f"{'Num. Patches (NP)':<30} | {pre_np:<15} | {post_np:<15} | count"
    )
    connected_post = f"{connected_groups}"
    results.append(
        f"{'Num. Connected Groups':<30} | {'0':<15} | {connected_post:<15} | "
        "Post = networks"
    )

    specs = [
        ("Total Area", "total_area", "{:,.2f}", area_unit),
        ("Eff. Mesh Size", "mesh", "{:,.2f}", f"{area_unit} (effective connected habitat at landscape scale)"),
        ("Largest Patch Size", "lps", "{:,.2f}", area_unit),
        ("Splitting Index", "split", "{:,.2f}", "Higher = more fragmented"),
        ("Largest Patch Index (LPI)", "lpi", "{:,.2f}", "% of landscape"),
        ("Mean Shape Index (MSI)", "msi", "{:,.2f}", "1 = compact, higher = irregular"),
        ("Edge Density (ED)", "ed", "{:,.2f}", f"{dist_unit} edge per {area_unit}"),
        ("Aggregation Index (AI)", "ai", "{:,.2f}", "0 = Dispersed, 100 = Clumped"),
        ("Core Area Index (CAI)", "cai", "{:,.2f}", "% core habitat (excluding edge effects)"),
    ]

    for label, key, fmt, interp in specs:
        pre_val = fmt.format(pre_metrics.get(key, 0.0))
        post_val = fmt.format(post_metrics.get(key, 0.0))
        results.append(f"{label:<30} | {pre_val:<15} | {post_val:<15} | {interp}")

    results.append(header)
    return results


def _run_raster_analysis_core(
    layer: QgsRasterLayer,
    output_dir: str,
    raw_params: Dict,
    strategy: str = "circuit_utility",
    temporary: bool = False,
    iface=None,
    progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
    params: Optional[RasterRunParams] = None,
    overall_start: Optional[float] = None,
) -> List[Dict]:
    """Core raster corridor logic (called by `run_raster_analysis`)."""
    if params is None:
        raise RasterAnalysisError("Raster parameters not supplied.")
    if overall_start is None:
        overall_start = time.time()

    ds = gdal.Open(layer.source())
    if ds is None:
        raise RasterAnalysisError(f"Cannot open raster source: {layer.source()}")

    rows, cols = ds.RasterYSize, ds.RasterXSize
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    print("=" * 70)
    print("LINKSCAPE RASTER ANALYSIS v23.0")
    print("=" * 70)
    print("\n1. Loading raster...")
    print(f"  ✓ Using layer: {layer.name()}")
    total_pixels = rows * cols
    print(f"  Size: {rows:,} x {cols:,} = {total_pixels:,} pixels")

    pixel_w = abs(gt[1])
    pixel_h = abs(gt[5]) if gt[5] != 0 else pixel_w
    pixel_size = max(pixel_w, pixel_h)
    unit_system = (raw_params or {}).get("raster_units", "pixels")

    warnings: List[str] = []
    blockers: List[str] = []

    if total_pixels >= PIXEL_COUNT_WARNING_THRESHOLD:
        blockers.append(
            "Large raster detected (>{:,} pixels).".format(
                PIXEL_COUNT_WARNING_THRESHOLD
            )
        )
    if 0 < pixel_size < PIXEL_SIZE_WARNING_THRESHOLD:
        msg = (
            f"High-resolution data detected (≈{pixel_size:.2f} map units per pixel). "
            "Consider resampling to a coarser resolution for faster corridor modelling."
        )
        if total_pixels >= PIXEL_FINE_CRITICAL_COUNT:
            blockers.append(msg)
        else:
            warnings.append(msg + " Proceeding because raster size is small.")

    if warnings:
        warning_text = " ".join(warnings)
        if iface and hasattr(iface, "messageBar"):
            try:
                iface.messageBar().pushMessage("TerraLink", warning_text, level=Qgis.Warning)
            except Exception:
                print(f"WARNING: {warning_text}")
        else:
            print(f"WARNING: {warning_text}")

    if blockers:
        warning_text = " ".join(blockers)
        if iface and hasattr(iface, "messageBar"):
            try:
                iface.messageBar().pushWarning("TerraLink", warning_text)
            except Exception:
                print(f"WARNING: {warning_text}")
        else:
            print(f"WARNING: {warning_text}")
        raise RasterAnalysisError(
            f"Raster is too large/fine for TerraLink to process efficiently. {warning_text} "
            "Please resample to a coarser resolution or process in smaller chunks."
        )
    _emit_progress(progress_cb, 5, "Loading raster data…")

    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    if nodata is None:
        nodata = params.nodata_fallback

    print("  Reading data...")
    data = read_band(
        band,
        rows,
        cols,
        progress_cb=progress_cb,
        progress_start=5,
        progress_end=18,
    )
    nodata_mask = np.abs(data - nodata) < params.value_tolerance if nodata is not None else np.zeros_like(
        data, dtype=bool
    )
    _emit_progress(progress_cb, 20, "Defining habitat patches…")

    print("\n2. Identifying patch pixels...")
    patch_mask = define_habitat(data, nodata_mask, params)
    habitat_mask = patch_mask.astype(np.uint8)
    patch_pixels = int(np.sum(habitat_mask))
    if patch_pixels == 0:
        raise RasterAnalysisError("No patch pixels found with the current configuration.")
    print(f"  ✓ Patch pixels: {patch_pixels:,}")
    _emit_progress(progress_cb, 25, "Processing habitat patches…")

    obstacle_mask = define_obstacles(data, nodata_mask, habitat_mask, params)
    if params.obstacle_enabled:
        obstacle_pixels = int(np.sum(obstacle_mask))
        if obstacle_pixels:
            print(f"  ✓ Impassable pixels: {obstacle_pixels:,}")
        else:
            print("  ⚠ Impassable land class configuration matched no pixels; proceeding without impassables.")
    else:
        print("  ✓ Impassable land classes disabled.")
    obstacle_pixels = int(np.sum(obstacle_mask))

    # NOTE: Background ("matrix") is already treated as passable by default: passable_mask
    # is derived from (~habitat) & (~impassable). Do not clear impassable pixels here.
    _emit_progress(progress_cb, 35, "Labeling patches…")

    print("\n3. Labeling patches...")
    t0 = time.time()
    labels, n_patches = label_patches(habitat_mask, params.patch_connectivity)
    print(f"  ✓ Patches: {n_patches:,} in {time.time() - t0:.2f}s")
    raw_patch_count = int(n_patches)

    # Preserve the full label map (including small patches) so the contiguous output
    # can include "bridge patches" that were filtered out by min_patch_size but were
    # later connected via a corridor.
    labels_all = labels.copy()

    if params.min_patch_size > 0:
        unique_labels, counts = np.unique(labels[labels > 0], return_counts=True)
        valid_labels = unique_labels[counts >= params.min_patch_size]
        new_labels = np.zeros_like(labels)
        for new_id, old_id in enumerate(valid_labels, 1):
            new_labels[labels == old_id] = new_id
        labels = new_labels
        print(f"  ✓ After filter: {len(valid_labels):,} patches")

    patch_mask = (labels > 0).astype(np.uint8)
    final_patch_ids = np.unique(labels[labels > 0])

    obstacle_mask = obstacle_mask.astype(bool) | patch_mask.astype(bool)

    passable_mask = _build_passable_mask(
        patch_mask,
        obstacle_mask,
        params.min_corridor_width,
        params.allow_bottlenecks,
    )
    passable_cells = int(np.sum(passable_mask))

    unique_labels, counts = np.unique(labels[labels > 0], return_counts=True)
    patch_sizes = dict(zip(unique_labels.tolist(), counts.tolist()))
    if not patch_sizes:
        raise RasterAnalysisError("No valid patches remain after filtering.")
    filtered_patch_count = len(patch_sizes)
    _emit_progress(progress_cb, 45, "Searching for corridors…")

    print("\n4. Finding possible corridors...")
    candidate_start = time.perf_counter()
    candidates = find_all_possible_corridors(
        labels,
        habitat_mask,
        patch_sizes,
        params.max_search_distance,
        params.max_corridor_area,
        params.patch_connectivity,
        min_corridor_width=params.min_corridor_width,
        obstacle_mask=obstacle_mask,
        passable_mask=passable_mask,
        strategy=strategy,
        progress_cb=progress_cb,
        progress_start=45,
        progress_end=75,
        params=params,
        timing_out=None,
    )
    candidate_elapsed = time.perf_counter() - candidate_start
    candidate_count = len(candidates)
    candidate_pairs: Set[Tuple[int, int]] = set()
    for cand in candidates:
        try:
            p1 = int(cand.get("patch1", 0))
            p2 = int(cand.get("patch2", 0))
        except Exception:
            continue
        if p1 <= 0 or p2 <= 0:
            continue
        candidate_pairs.add((p1, p2) if p1 <= p2 else (p2, p1))
    candidate_pairs_count = len(candidate_pairs)
    possible_pairs = (filtered_patch_count * (filtered_patch_count - 1)) // 2 if filtered_patch_count > 1 else 0

    if not candidates:
        raise RasterAnalysisError("No feasible corridors found with the current configuration.")
    _emit_progress(progress_cb, 78, "Optimizing corridor selection…")

    strategy = (strategy or "circuit_utility").lower()
    if strategy not in ("largest_network", "circuit_utility"):
        strategy = "circuit_utility"
    safe_layer = _safe_filename(layer.name())
    default_filename = f"terralink_contiguous_{safe_layer}.tif"
    strategy_map = {
        "largest_network": (
            optimize_largest_network,
            default_filename,
            "Contiguous Areas (Largest Single Network)",
        ),
        "circuit_utility": (
            optimize_circuit_utility,
            default_filename,
            "Contiguous Areas (Most Connectivity)",
        ),
    }

    if strategy not in strategy_map:
        raise RasterAnalysisError(f"Unsupported strategy '{strategy}'.")

    optimize_func, default_filename, layer_name = strategy_map[strategy]

    print("\n5. Running optimization...")
    print("=" * 70)
    print(f"--- {strategy.replace('_', ' ').upper()} ---")

    corridors, stats = optimize_func(
        candidates,
        patch_sizes,
        params.budget_pixels,
        params.max_search_distance,
    )
    if not corridors:
        raise RasterAnalysisError("Selected optimization did not produce any corridors.")
    remaining_budget = int(max(0, int(params.budget_pixels) - int(stats.get("budget_used", 0) or 0)))
    if remaining_budget > 0:
        min_extra_cost = _min_unselected_candidate_cost(candidates, corridors)
        if min_extra_cost is not None and remaining_budget < int(min_extra_cost):
            stats["remaining_budget_unusable"] = remaining_budget
            stats["min_unselected_corridor_cost"] = int(min_extra_cost)
            print(
                "  Remaining budget too small to add new corridors; skipping extra selection."
            )
        else:
            extra_used, low_value_added, redundancy_added = _apply_hybrid_leftover_budget_raster(
                candidates=candidates,
                corridors=corridors,
                patch_sizes=patch_sizes,
                remaining_budget=remaining_budget,
                max_search_distance=int(params.max_search_distance or 0),
            )
            if extra_used:
                stats["budget_used"] = int(stats.get("budget_used", 0) or 0) + int(extra_used)
                stats["corridors_used"] = len(corridors)
                if "primary_links" in stats:
                    stats["primary_links"] = int(stats.get("primary_links", 0) or 0) + int(low_value_added)
                if "redundant_links" in stats:
                    stats["redundant_links"] = int(stats.get("redundant_links", 0) or 0) + int(redundancy_added)
                _refresh_raster_component_stats(corridors, patch_sizes, stats)
            remaining_budget = int(max(0, int(params.budget_pixels) - int(stats.get("budget_used", 0) or 0)))
    if remaining_budget > 0:
        thickened = _thicken_corridors_raster(
            corridors=corridors,
            remaining_budget=remaining_budget,
            params=params,
            patch_sizes=patch_sizes,
            obstacle_mask=obstacle_mask,
            labels_shape=labels.shape,
            max_width_factor=3,
        )
        if thickened:
            stats["budget_used"] = int(stats.get("budget_used", 0) or 0) + int(thickened)
    _emit_progress(progress_cb, 85, "Rendering output raster…")

    print("  Creating corridor and contiguous network area rasters...")
    corridor_output = create_output_raster(
        labels, corridors, params.min_corridor_width, obstacle_mask=obstacle_mask
    )
    network_output = create_contiguous_network_area_raster(
        labels,
        corridors,
        params.min_corridor_width,
        obstacle_mask=obstacle_mask,
        labels_all=labels_all,
    )

    if temporary:
        temp_file = tempfile.NamedTemporaryFile(prefix="terralink_contiguous_", suffix=".tif", delete=False)
        out_path = temp_file.name
        temp_file.close()
        temp_corr = tempfile.NamedTemporaryFile(prefix="terralink_corridors_", suffix=".tif", delete=False)
        corridor_path = temp_corr.name
        temp_corr.close()
    else:
        out_dir = output_dir or os.path.dirname(src_path)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, default_filename)
        corridor_path = os.path.join(out_dir, f"corridors_{default_filename}")

    try:
        # Corridor raster (connected_size per corridor cell)
        write_raster(corridor_path, corridor_output, gt, proj, nodata=0)
        print(f"  ✓ Saved corridor raster: {corridor_path}")
        corr_layer = None
        try:
            corr_layer_name = layer_name.replace("Contiguous Areas", "Corridors")
            corr_layer = QgsRasterLayer(corridor_path, corr_layer_name)
            if corr_layer.isValid():
                pass
        except Exception:
            pass

        write_raster(out_path, network_output, gt, proj, nodata=0)
        print(f"  ✓ Saved: {out_path}")

        # --- LANDSCAPE METRICS REPORT (always written) ---
        landscape_metrics_path = ""
        analysis_lines: List[str] = []
        analysis_title = f"TerraLink Landscape Metrics ({layer.name()})"
        try:
            is_metric = (layer.crs().mapUnits() == QgsUnitTypes.DistanceMeters)
            analysis_lines = _perform_landscape_analysis(
                network_output,
                analysis_title,
                abs(gt[1]),
                abs(gt[5]),
                is_metric,
                params=params,
                pre_arr=patch_mask,
            )
            safe = _safe_filename(layer.name())
            landscape_metrics_path = os.path.join(
                os.path.dirname(out_path) or os.getcwd(),
                f"landscape_metrics_{safe}.txt",
            )
            _write_text_report(landscape_metrics_path, analysis_lines)
            print(f"  ✓ Saved landscape metrics: {landscape_metrics_path}")
            try:
                _add_landscape_metrics_table_layer(analysis_title, analysis_lines)
            except Exception:
                pass
        except Exception as e:
            # Always write a report, even if computation fails.
            try:
                safe = _safe_filename(layer.name())
                landscape_metrics_path = os.path.join(
                    os.path.dirname(out_path) or os.getcwd(),
                    f"landscape_metrics_{safe}.txt",
                )
                analysis_lines = [f"Landscape analysis failed: {e}"]
                _write_text_report(landscape_metrics_path, analysis_lines)
                print(f"  ✓ Saved landscape metrics: {landscape_metrics_path}")
                try:
                    _add_landscape_metrics_table_layer(analysis_title, analysis_lines)
                except Exception:
                    pass
            except Exception:
                print(f"  ⚠ Landscape analysis failed: {e}")
        net_layer = QgsRasterLayer(out_path, layer_name)
        if net_layer.isValid():
            try:
                _apply_random_unique_value_symbology_from_array(net_layer, network_output)
            except Exception:
                pass
            QgsProject.instance().addMapLayer(net_layer)
            print("  ✓ Added to project")
        if corr_layer is not None and corr_layer.isValid():
            QgsProject.instance().addMapLayer(corr_layer)
            print("  ✓ Added corridor raster to project")
    except Exception as net_exc:  # noqa: BLE001
        print(f"  ⚠ Could not save/add output raster: {net_exc}")

    _emit_progress(progress_cb, 95, "Finishing up…")

    elapsed = time.time() - overall_start
    stats = dict(stats)
    try:
        if landscape_metrics_path:
            stats["landscape_metrics_path"] = landscape_metrics_path
    except Exception:
        pass
    stats["output_path"] = out_path if not temporary else ""
    try:
        stats["corridor_output_path"] = corridor_path if not temporary else ""
    except Exception:
        stats["corridor_output_path"] = ""
    stats["layer_name"] = layer_name
    stats["budget_total"] = params.budget_pixels
    stats["patches_total"] = stats.get("patches_total", len(patch_sizes))
    stats["habitat_pixels_total"] = sum(patch_sizes.values())
    stats["raster_rows"] = rows
    stats["raster_cols"] = cols
    stats["raster_pixels_total"] = total_pixels
    stats["patches_total_raw"] = raw_patch_count
    stats["patches_total_filtered"] = filtered_patch_count
    stats["min_patch_size"] = params.min_patch_size
    stats["max_search_distance"] = params.max_search_distance
    stats["min_corridor_width"] = params.min_corridor_width
    stats["candidate_corridors"] = candidate_count
    stats["candidate_pairs"] = candidate_pairs_count
    stats["possible_pairs"] = possible_pairs
    stats["candidate_search_s"] = candidate_elapsed
    stats["passable_cells"] = passable_cells
    stats["obstacle_pixels"] = obstacle_pixels
    stats.update(_compute_connectivity_metrics(corridors, patch_sizes))
    _apply_raster_display_units(
        stats,
        layer=layer,
        unit_system=unit_system,
        pixel_w=pixel_w,
        pixel_h=pixel_h,
    )

    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    stats_strategy = stats.get("strategy", strategy)
    print(f"Strategy:          {stats_strategy.replace('_', ' ').title()}")
    print(f"Corridors created: {stats.get('corridors_used', 0)}")
    print(f"Connected groups:  {stats.get('components_remaining', 0)}")
    print(f"Corridors used:    {stats.get('corridors_used', 0)}")
    print(f"Avg degree:        {stats.get('avg_degree', 0):.2f}")
    print(f"Redundant links:   {stats.get('redundant_links', 0)}")
    if "connections_made" in stats:
        print(f"Connections:       {stats.get('connections_made', 0)}")
    if "seed_id" in stats:
        print(f"Seed patch:        {stats.get('seed_id')}")
    print(f"Final size:        {stats.get('final_patch_size', 0):,} px")
    area_label = stats.get("raster_area_label", "px")
    if area_label == "px":
        print(f"Corridor budget:   {params.budget_pixels} px")
        print(f"Budget used:       {stats.get('budget_used', 0)} px")
    else:
        print(f"Corridor budget:   {stats.get('budget_total_display', 0):.2f} {area_label}")
        print(f"Budget used:       {stats.get('budget_used_display', 0):.2f} {area_label}")
    if "entropy_total" in stats:
        print("\nNetwork Metrics:")
        print(f"  Entropy (H_total): {stats.get('entropy_total', 0):.4f}")
        if "entropy_movement" in stats:
            print(f"  Movement entropy:  {stats.get('entropy_movement', 0):.4f}")
        if "robustness_rho2" in stats:
            print(f"  Robustness (ρ₂):   {stats.get('robustness_rho2', 0):.4f}")
        if "topology_penalty" in stats:
            print(f"  Topology penalty:  {stats.get('topology_penalty', 0):.1f}")
        print(f"  Mode:              {stats.get('mode', 'N/A')}")
    print(f"Processing time:   {elapsed:.1f}s")
    if temporary:
        print("Output:            Temporary raster layer")
    else:
        print(f"Output GeoTIFF:    {out_path}")
    print("=" * 70)

    # Add a lightweight in-project summary table instead of writing/opening text files.
    try:
        _add_raster_run_summary_layer(
            input_layer_name=layer.name(),
            strategy=strategy.replace("_", " ").title(),
            stats=stats,
            out_path=out_path if not temporary else "",
            corridor_path=corridor_path if "corridor_path" in locals() and not temporary else None,
            combined_path=None,
        )
    except Exception:
        pass

    ds = None
    _emit_progress(progress_cb, 100, "Raster analysis complete.")
    return [{"strategy": strategy, "stats": stats, "output_path": out_path if not temporary else ""}]
