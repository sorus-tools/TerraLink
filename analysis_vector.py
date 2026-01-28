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
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
# NumPy 2.x removed np.int; add shim for any legacy references.
if not hasattr(np, "int"):  # pragma: no cover
    np.int = int  # type: ignore[attr-defined]

import networkx as nx  # Required for graph-based optimization/metrics
from .terralink_engine import NetworkOptimizer, UnionFind
from .core_select import select_circuit_utility
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



def _add_summary_csv_layer(csv_path: str, layer_title: str) -> None:
    if QgsVectorLayer is None or QgsProject is None:
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


def _safe_filename(name: str, max_len: int = 64) -> str:
    safe = "".join(ch if (ch.isalnum() or ch in ("-", "_")) else "_" for ch in (name or ""))
    safe = safe.strip("_") or "layer"
    return safe[:max_len]


def _add_landscape_metrics_table_layer(layer_title: str, analysis_lines: List[str]) -> None:
    """
    Add a non-spatial table layer to the QGIS project with the landscape metrics.
    """
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


@dataclass
class AnalysisContext:
    """
    Mutable, per-run state (caches) that should not live on VectorRunParams.
    Keeping params immutable-ish avoids state leakage across runs.
    """

    impassable_union: Optional[QgsGeometry] = None


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
) -> Tuple[Dict[int, Dict], QgsSpatialIndex]:
    print("  Loading patches and building spatial index...")
    source_crs = layer.crs()
    transform = QgsCoordinateTransform(source_crs, target_crs, QgsProject.instance())

    raw_geoms: List[QgsGeometry] = []
    indexed_features: List[QgsFeature] = []

    for feature in layer.getFeatures():
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
        fid = len(raw_geoms)
        raw_geoms.append(geom)
        feat = QgsFeature()
        feat.setGeometry(geom)
        feat.setId(fid)
        indexed_features.append(feat)

    spatial_index = QgsSpatialIndex(flags=QgsSpatialIndex.FlagStoreFeatureGeometries)
    if indexed_features:
        spatial_index.addFeatures(indexed_features)

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

    components: Dict[int, List[QgsGeometry]] = defaultdict(list)
    for idx, geom in enumerate(raw_geoms):
        root = int(uf.find(idx))
        components[root].append(geom)

    patches: Dict[int, Dict] = {}
    patch_id = 1
    filtered_count = 0
    indexed_features = []

    for geoms in components.values():
        if not geoms:
            continue
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
            continue

        feat = QgsFeature()
        feat.setGeometry(merged)
        feat.setId(patch_id)
        indexed_features.append(feat)
        patches[patch_id] = {
            "geom": clone_geometry(merged),
            "area_ha": area_ha,
        }
        patch_id += 1

    spatial_index = QgsSpatialIndex(flags=QgsSpatialIndex.FlagStoreFeatureGeometries)
    if indexed_features:
        spatial_index.addFeatures(indexed_features)

    print(f"  ✓ Loaded {len(patches)} patches (filtered {filtered_count} too small)")
    return patches, spatial_index


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

    return final_geom, patch_ids


def _buffer_line_segment(line_geom: QgsGeometry, width: float) -> QgsGeometry:
    return line_geom.buffer(width / 2.0, BUFFER_SEGMENTS)


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
        f"- Min corridor width: {params.min_corridor_width:.2f}\n"
        "Try increasing the search distance, lowering the minimum corridor width/area, "
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
    strategy: str = "circuit_utility",
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
    strategy_key = (strategy or "circuit_utility").lower()
    if strategy_key not in ("largest_network", "circuit_utility"):
        strategy_key = "circuit_utility"
    circuit_mode = strategy_key == "circuit_utility"
    largest_network_mode = strategy_key == "largest_network"
    # Both Most Connectivity and Largest Single Network (Circuit Utility) need enough candidate
    # diversity to avoid "blind" terminal selection in long/snaking patches.
    utility_mode = strategy_key in ("circuit_utility", "largest_network")
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

    def _push_candidate(cand: Dict) -> None:
        try:
            p1 = int(cand.get("patch1"))
            p2 = int(cand.get("patch2"))
        except Exception:
            return
        if p1 == p2:
            return
        geom = cand.get("geom")
        if geom is None or geom.isEmpty():
            return
        key = _pair_key(p1, p2)
        existing = candidates_by_pair.get(key, [])

        for prev in existing:
            try:
                if _overlap_ratio(geom, prev.get("geom")) >= min_distinct_overlap_ratio:
                    return
            except Exception:
                continue

        existing.append(cand)
        existing.sort(key=lambda c: (c.get("area_ha", 0.0), c.get("distance_m", 0.0)))
        if len(existing) > max_keep_per_pair:
            del existing[max_keep_per_pair:]
        candidates_by_pair[key] = existing

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

                # We evaluate multiple terminals above, but still emit only ONE corridor per patch-pair
                # to avoid parallel duplicate corridors between the same patch IDs.
                variants = [("nearest", p1_xy, p2_xy)]
                for variant_tag, start_xy, end_xy in variants:
                    raw_geom_candidates: List[Tuple[str, QgsGeometry]] = []
                    path_points: Optional[List[QgsPointXY]] = None

                    if navigator:
                        with _measure("Navigator routing"):
                            # Always reuse the best evaluated path for the selected terminals.
                            path_points = best_path_points

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

                    best_source, raw_geom = min(
                        raw_geom_candidates,
                        key=lambda item: item[1].area() if item[1] and not item[1].isEmpty() else float("inf"),
                    )
                    if raw_geom is None or raw_geom.isEmpty():
                        continue

                    # Reject candidates that "reach through" an endpoint patch to exit on the far side.
                    # We do this by measuring how much of the *raw centerline* lies within each endpoint patch.
                    raw_line = None
                    try:
                        if navigator and path_points:
                            raw_line = QgsGeometry.fromPolylineXY(path_points)
                        else:
                            raw_line = QgsGeometry.fromPolylineXY([start_xy, end_xy])
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
                        corridor_geom, patch_ids = _finalize_corridor_geometry(
                            pid1, pid2, raw_geom, patches, spatial_index, patch_union=patch_union
                        )
                    if corridor_geom is None:
                        continue

                    emitted_any = False
                    try:
                        extra_patches = [int(pid) for pid in patch_ids if int(pid) not in (int(pid1), int(pid2))]
                    except Exception:
                        extra_patches = []

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
                            area_a = float(part_a.area() / 10000.0)
                            if area_a > 0 and (params.max_corridor_area is None or area_a <= params.max_corridor_area):
                                try:
                                    dist_a = float(geom1.distance(patches[tgt_a]["geom"]))
                                except Exception:
                                    dist_a = float(distance)
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
                                    }
                                )
                                emitted_any = True

                        part_b = _pick_part_closest_to_patch(corridor_geom, geom2)
                        if part_b is not None and (not part_b.isEmpty()):
                            try:
                                tgt_b = min(target_patches, key=lambda pid: float(geom2.distance(patches[pid]["geom"])))
                            except Exception:
                                tgt_b = target_patches[0]
                            area_b = float(part_b.area() / 10000.0)
                            if area_b > 0 and (params.max_corridor_area is None or area_b <= params.max_corridor_area):
                                try:
                                    dist_b = float(geom2.distance(patches[tgt_b]["geom"]))
                                except Exception:
                                    dist_b = float(distance)
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
                                    }
                                )
                                emitted_any = True

                    if emitted_any:
                        continue

                    # If we touched an intermediate patch but couldn't split into endpoint-local pieces,
                    # drop this candidate in Largest Single Network so it doesn't "jump over" the patch.
                    if largest_network_mode and extra_patches:
                        continue

                    corridor_area_ha = float(corridor_geom.area() / 10000.0)
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
                            "distance_m": float(distance),
                            "variant": variant_tag,
                            "source": best_source,
                        }
                    )
                processed_pairs.add(pair)

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


def optimize_largest_network_strategy(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
    allow_loops: bool = True,
    mode: str = "resilient",
) -> Tuple[Dict[int, Dict], Dict]:
    """
    Graph-based Strategy: Entropy Minimization

    Modes:
      - resilient: MST + optional strategic loops (multi-network allowed).
      - largest_network: Grow the largest single network (seed = largest patch), then optional loops.
    """
    if graph_math is None:
        raise VectorAnalysisError("NetworkX is required for Largest Single Network optimization but could not be imported.")

    mode = (mode or "resilient").lower()
    selected: Dict[int, Dict] = {}
    selected_geoms: List[QgsGeometry] = []
    selected_pairs: Set[Tuple[int, int]] = set()
    budget_used = 0.0
    remaining = params.budget_area

    # ------------------------------------------------------------------
    # Phase 1: Backbone construction
    # ------------------------------------------------------------------
    def _add_corridor(cand: Dict, corr_type: str, apply_distance_guard: bool = True) -> None:
        nonlocal budget_used, remaining
        try:
            a = int(cand.get("patch1"))
            b = int(cand.get("patch2"))
        except Exception:
            return
        geom = cand.get("geom")
        if geom is None or geom.isEmpty():
            return
        if apply_distance_guard:
            if not _redundancy_far_enough(
                geom, selected_geoms, float(getattr(params, "max_search_distance", 0.0) or 0.0)
            ):
                return
        pk = (a, b) if a <= b else (b, a)
        # Largest Single Network should not contain parallel corridors between the same patch pair.
        # (Rewires explicitly delete the old graph edge before adding the replacement.)
        if pk in selected_pairs:
            return
        base_roi = float(_get_base_roi(cand) or 0.0)
        target_sizes = [float(_patch_area_ha(pid) or 0.0) for pid in pids]
        target_importance = min(target_sizes) if target_sizes else 1.0
        if target_importance <= 0:
            target_importance = 1.0
        base_roi = base_roi * target_importance
        cid = len(selected) + 1
        selected[cid] = {
            "geom": clone_geometry(cand["geom"]),
            "patch_ids": cand.get("patch_ids", {cand.get("patch1"), cand.get("patch2")}),
            "area_ha": cand.get("area_ha", 0.0),
            "p1": cand.get("patch1"),
            "p2": cand.get("patch2"),
            "distance": cand.get("distance_m", 1.0),
            "type": corr_type,
            "variant": cand.get("variant"),
            "source": cand.get("source"),
        }
        selected_pairs.add(pk)
        budget_used += cand.get("area_ha", 0.0)
        remaining -= cand.get("area_ha", 0.0)
        selected_geoms.append(clone_geometry(geom))

    backbone_graph_edges: List[Tuple[int, int, float, float]] = []

    if mode == "largest_network":
        print("  Strategy: Largest Single Network")
        # Seed with largest patch by area and grow outward (Prim-style) to avoid hub-and-spoke
        seed_patch = max(patches.keys(), key=lambda pid: patches[pid]["area_ha"])
        visited: Set[int] = {seed_patch}

        adjacency: Dict[int, List[Dict]] = defaultdict(list)
        for cand in candidates:
            p1, p2 = cand.get("patch1"), cand.get("patch2")
            adjacency[p1].append(cand)
            adjacency[p2].append(cand)

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

        heap: List[Tuple[float, float, int, Dict]] = []
        counter = 0
        for cand in adjacency.get(seed_patch, []):
            counter += 1
            heapq.heappush(
                heap,
                (cand.get("distance_m", float("inf")), cand.get("area_ha", float("inf")), counter, cand),
            )

        def _best_bridge_pair() -> Optional[Tuple[Dict, Dict, int]]:
            best = None
            best_score = 0.0
            for mid, edges in bridge_by_mid.items():
                if mid in visited:
                    continue
                best_to_visited = None
                best_to_visited_cost = 0.0
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
                    cost = float(cand.get("area_ha", 0.0) or 0.0)
                    if cost <= 0.0:
                        continue
                    if other in visited:
                        if best_to_visited is None or cost < best_to_visited_cost:
                            best_to_visited = cand
                            best_to_visited_cost = cost
                    else:
                        gain_other = float(patches.get(int(other), {}).get("area_ha", 0.0) or 0.0)
                        if gain_other <= 0.0:
                            continue
                        score = gain_other / cost
                        if best_to_unvisited is None or score > best_to_unvisited_score:
                            best_to_unvisited = cand
                            best_to_unvisited_score = score
                            best_to_unvisited_other = int(other)
                if best_to_visited is None or best_to_unvisited is None or best_to_unvisited_other is None:
                    continue
                total_cost = float(best_to_visited_cost + (best_to_unvisited.get("area_ha", 0.0) or 0.0))
                if total_cost <= 0.0 or total_cost > remaining:
                    continue
                gain_mid = float(patches.get(int(mid), {}).get("area_ha", 0.0) or 0.0)
                gain_other = float(patches.get(int(best_to_unvisited_other), {}).get("area_ha", 0.0) or 0.0)
                score = (gain_mid + gain_other) / total_cost if total_cost else 0.0
                if score > best_score:
                    best_score = score
                    best = (best_to_visited, best_to_unvisited, int(best_to_unvisited_other))
            return best

        bridge_iters = 0
        bridge_start = time.perf_counter()
        while remaining > 0:
            if bridge_iters >= MAX_BRIDGE_ITERATIONS:
                break
            if (time.perf_counter() - bridge_start) > MAX_BRIDGE_RUNTIME_S:
                break
            bridge_pick = _best_bridge_pair()
            if bridge_pick is None:
                break
            cand1, cand2, other_node = bridge_pick
            cost1 = float(cand1.get("area_ha", 0.0) or 0.0)
            cost2 = float(cand2.get("area_ha", 0.0) or 0.0)
            if cost1 <= 0.0 or cost2 <= 0.0 or (cost1 + cost2) > remaining:
                break
            _add_corridor(cand1, "backbone", apply_distance_guard=False)
            _add_corridor(cand2, "backbone", apply_distance_guard=False)
            p1, p2 = int(cand1.get("patch1")), int(cand1.get("patch2"))
            mid_node = p1 if p1 not in visited else p2
            visited.add(int(mid_node))
            visited.add(int(other_node))
            backbone_graph_edges.append((p1, p2, cost1, cand1.get("distance_m", 1.0)))
            backbone_graph_edges.append((int(cand2.get("patch1")), int(cand2.get("patch2")), cost2, cand2.get("distance_m", 1.0)))

            for new_node in (int(mid_node), int(other_node)):
                for nxt in adjacency.get(new_node, []):
                    n_p1, n_p2 = nxt.get("patch1"), nxt.get("patch2")
                    if (n_p1 in visited and n_p2 in visited) or (n_p1 not in visited and n_p2 not in visited):
                        continue
                    counter += 1
                    heapq.heappush(
                        heap,
                        (nxt.get("distance_m", float("inf")), nxt.get("area_ha", float("inf")), counter, nxt),
                    )
            bridge_iters += 1

        while heap and remaining > 0:
            _dist, cost, _idx, cand = heapq.heappop(heap)
            if cost > remaining:
                continue
            p1, p2 = cand.get("patch1"), cand.get("patch2")
            in1, in2 = p1 in visited, p2 in visited
            if in1 and in2:
                continue  # already connected; skip to avoid redundant spokes
            if not in1 and not in2:
                continue  # does not touch current network

            new_node = p2 if in1 else p1
            _add_corridor(cand, "backbone", apply_distance_guard=False)
            backbone_graph_edges.append((p1, p2, cost, cand.get("distance_m", 1.0)))
            visited.add(new_node)

            for nxt in adjacency.get(new_node, []):
                n_p1, n_p2 = nxt.get("patch1"), nxt.get("patch2")
                if (n_p1 in visited and n_p2 in visited) or (n_p1 not in visited and n_p2 not in visited):
                    continue
                counter += 1
                heapq.heappush(
                    heap,
                    (nxt.get("distance_m", float("inf")), nxt.get("area_ha", float("inf")), counter, nxt),
                )

    else:
        print(f"  Strategy: Resilient (loops allowed={allow_loops})")
        _log_message("--- Phase 1: Building Backbone ---")
        mst_candidates = sorted(
            candidates,
            key=lambda x: (
                x.get("distance_m", float("inf")),
                x.get("area_ha", float("inf")),
            ),
        )
        uf = UnionFind()
        for pid in patches:
            uf.find(pid)

        for cand in mst_candidates:
            cost = cand.get("area_ha", 0.0)
            p1 = cand.get("patch1")
            p2 = cand.get("patch2")

            if cost > remaining:
                continue

            if uf.union(p1, p2):
                _add_corridor(cand, "backbone")
                backbone_graph_edges.append((p1, p2, cost, cand.get("distance_m", 1.0)))

    _log_message(f"  Backbone complete. Budget used: {budget_used:.2f}/{params.budget_area:.2f}")

    # ------------------------------------------------------------------
    # Phase 2a: Parallel Reinforcement (Optional)
    # ------------------------------------------------------------------
    def _pair_key(p1: int, p2: int) -> Tuple[int, int]:
        return (p1, p2) if p1 <= p2 else (p2, p1)

    def _overlap_ratio(g1: QgsGeometry, g2: QgsGeometry) -> float:
        try:
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

    # Largest Single Network should not add a second corridor between the same backbone patch pair.
    # (Users expect redundancy via alternate patch pairs, not parallel duplicates.)
    if mode != "largest_network" and remaining > 0 and selected and backbone_graph_edges:
        _log_message("--- Phase 2a: Reinforcing Existing Links ---")
        # Prefer spending remaining budget on redundant (spatially distinct) links
        # between already-backbone-connected patch pairs, rather than creating many
        # long-distance loops across the landscape.
        max_total_links_per_pair = 2  # 1 backbone + 1 redundant
        overlap_reject_ratio = 0.85

        backbone_pairs: Set[Tuple[int, int]] = {_pair_key(p1, p2) for p1, p2, _c, _d in backbone_graph_edges}

        selected_geoms_by_pair: Dict[Tuple[int, int], List[QgsGeometry]] = defaultdict(list)
        for data in selected.values():
            p1 = data.get("p1")
            p2 = data.get("p2")
            if p1 is None or p2 is None:
                continue
            g = data.get("geom")
            if g is None or g.isEmpty():
                continue
            selected_geoms_by_pair[_pair_key(int(p1), int(p2))].append(g)

        reinforce_ranked: List[Tuple[float, float, Dict]] = []
        for cand in candidates:
            p1, p2 = cand.get("patch1"), cand.get("patch2")
            if p1 is None or p2 is None:
                continue
            pk = _pair_key(int(p1), int(p2))
            if pk not in backbone_pairs:
                continue
            if len(selected_geoms_by_pair.get(pk, [])) >= max_total_links_per_pair:
                continue
            cost = float(cand.get("area_ha", 0.0))
            if cost <= 0 or cost > remaining:
                continue
            g = cand.get("geom")
            if g is None or g.isEmpty():
                continue
            prior_geoms = selected_geoms_by_pair.get(pk, [])
            overlap = max((_overlap_ratio(g, pg) for pg in prior_geoms), default=0.0)
            if overlap >= overlap_reject_ratio:
                continue
            dist = float(cand.get("distance_m", 1.0))
            # Reward: low cost, low overlap, short distance.
            score = (1.0 - overlap) / (cost + 1e-9) * (1.0 / (dist + 1e-6))
            reinforce_ranked.append((score, overlap, cand))

        reinforce_ranked.sort(key=lambda t: t[0], reverse=True)

        reinforced = 0
        for _score, overlap, cand in reinforce_ranked:
            cost = float(cand.get("area_ha", 0.0))
            if cost > remaining:
                continue
            p1 = int(cand.get("patch1"))
            p2 = int(cand.get("patch2"))
            pk = _pair_key(p1, p2)
            if len(selected_geoms_by_pair.get(pk, [])) >= max_total_links_per_pair:
                continue
            _add_corridor(cand, "reinforcement")
            geom = cand.get("geom")
            if geom is not None:
                selected_geoms_by_pair.setdefault(pk, []).append(clone_geometry(geom))
            reinforced += 1

        _log_message(f"  Added {reinforced} redundant link(s) across backbone pairs.")

    # ------------------------------------------------------------------
    # Phase 2: Strategic Loops (Optional)
    # ------------------------------------------------------------------
    if allow_loops and remaining > 0 and selected:
        _log_message("--- Phase 2: Adding Strategic Loops ---")

        G_backbone = nx.Graph()
        for pid in patches:
            G_backbone.add_node(pid)
        for p1, p2, _, dist in backbone_graph_edges:
            G_backbone.add_edge(p1, p2, weight=dist)

        selected_geoms: List[QgsGeometry] = [
            clone_geometry(d["geom"])
            for d in selected.values()
            if d.get("geom") is not None and not d.get("geom").isEmpty()
        ]
        redundant_distance = float(getattr(params, "max_search_distance", 0.0) or 0.0)

        loop_candidates = []
        for cand in candidates:
            p1, p2 = cand.get("patch1"), cand.get("patch2")
            if not nx.has_path(G_backbone, p1, p2):
                continue

            cost = float(cand.get("area_ha", 0.0) or 0.0)
            if cost <= 0:
                continue

            cand_dist = float(cand.get("distance_m", 1.0) or 1.0)
            cand_geom = cand.get("geom")
            if cand_geom is None or cand_geom.isEmpty():
                continue

            # Rewire: if there is already a direct backbone graph edge, and this candidate is significantly shorter,
            # swap it in (acyclic replacement), paying only the cost delta.
            if G_backbone.has_edge(p1, p2):
                existing_dist = float(G_backbone[p1][p2].get("weight", 1.0) or 1.0)
                if cand_dist >= existing_dist * 0.85:
                    continue

                old_cid = None
                old_cost = 0.0
                for cid, data in selected.items():
                    if data.get("type") != "backbone":
                        continue
                    a = data.get("p1")
                    b = data.get("p2")
                    if (a == p1 and b == p2) or (a == p2 and b == p1):
                        old_cid = cid
                        old_cost = float(data.get("area_ha", 0.0) or 0.0)
                        break
                if old_cid is None:
                    continue

                extra_needed = max(0.0, cost - old_cost)
                if extra_needed > remaining:
                    continue

                # Refund old graph-edge cost then add new graph-edge cost (net = extra_needed).
                del selected[old_cid]
                budget_used -= old_cost
                remaining += old_cost
                try:
                    pk_old = _pair_key(int(p1), int(p2))
                    selected_pairs.discard(pk_old)
                except Exception:
                    pass
                _add_corridor(cand, "backbone")
                G_backbone[p1][p2]["weight"] = cand_dist

                for i, (a, b, c, d) in enumerate(backbone_graph_edges):
                    if (a == p1 and b == p2) or (a == p2 and b == p1):
                        backbone_graph_edges[i] = (p1, p2, cost, cand_dist)
                        break

                continue

            if redundant_distance > 0.0 and not _redundancy_far_enough(
                cand_geom, selected_geoms, redundant_distance
            ):
                continue

            # Otherwise, treat it as a loop candidate (shortcuts after connectivity exists).
            if cost <= remaining:
                score = graph_math.score_edge_for_loops(G_backbone, p1, p2, cand_dist)
                efficiency = score / cost if cost else 0.0
                loop_candidates.append((efficiency, cand))

        loop_candidates.sort(key=lambda x: x[0], reverse=True)

        loops_added = 0
        for _eff, cand in loop_candidates:
            cost = cand.get("area_ha", 0.0)
            if cost <= remaining:
                _add_corridor(cand, "strategic_loop")
                p1, p2 = cand.get("patch1"), cand.get("patch2")
                backbone_graph_edges.append((p1, p2, cost, cand.get("distance_m", 1.0)))
                g = cand.get("geom")
                if g is not None and not g.isEmpty():
                    selected_geoms.append(clone_geometry(g))
                loops_added += 1

        _log_message(f"  Added {loops_added} strategic loops.")

    # ------------------------------------------------------------------
    # Phase 3: Metrics
    # ------------------------------------------------------------------
    final_corridors_list = []
    for cid, data in selected.items():
        final_corridors_list.append(
            {"id": cid, "patch1": data.get("p1"), "patch2": data.get("p2"), "distance_m": data.get("distance", 1.0)}
        )

    G_final = graph_math.build_graph_from_corridors(patches, final_corridors_list)
    entropy_stats = graph_math.calculate_total_entropy(G_final)
    rho2 = graph_math.calculate_two_edge_connectivity(G_final)

    stats = {
        "strategy": "largest_network",
        "corridors_used": len(selected),
        "budget_used_ha": budget_used,
        "patches_connected": G_final.number_of_nodes(),
        "entropy_total": entropy_stats["H_total"],
        "robustness_rho2": rho2,
        "mode": "Largest Single Network",
        "total_connected_area_ha": sum(p["area_ha"] for p in patches.values()),
        "largest_group_area_ha": 0.0,
    }

    comps = list(nx.connected_components(G_final))
    if comps:
        max_area = max(sum(patches[pid]["area_ha"] for pid in comp) for comp in comps)
        stats["largest_group_area_ha"] = max_area

    for data in selected.values():
        connected_area = stats.get("largest_group_area_ha", 0.0)
        data["connected_area_ha"] = connected_area
        area = data.get("area_ha", 0.0)
        data["efficiency"] = (connected_area / area) if area else 0.0

    return selected, stats


def optimize_circuit_utility(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
    overlap_reject_ratio: float = 0.30,
) -> Tuple[Dict[int, Dict], Dict]:
    """
    Most Connectivity (Utility) strategy:
    Weighted greedy selection by marginal utility per cost.

    Utility Score:
      score = (sqrt(area_p1 * area_p2) / cost) * multiplier

    Multiplier:
      1.0  if candidate connects different components
      0.5  if already connected but spatially distinct
      0.01 if overlaps existing corridor for the pair by > overlap_reject_ratio
    """
    def _pair_key(a: int, b: int) -> Tuple[int, int]:
        return (a, b) if a <= b else (b, a)

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
                    proximity_dist = max(float(getattr(params, "min_corridor_width", 0.0) or 0.0) * 1.5, float(getattr(params, "grid_resolution", 0.0) or 0.0) * 2.0)
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

    def _get_cost(cand: Dict) -> float:
        return float(_candidate_cost_ha(cand) or 0.0)

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

    def _overlap_obj(cand: Dict) -> object:
        g = cand.get("geom")
        if g is None or g.isEmpty():
            return None
        return clone_geometry(g)

    def _overlap_ratio(cand: Dict, prior: Sequence[object]) -> float:
        g = cand.get("geom")
        if g is None or g.isEmpty() or not prior:
            return 0.0
        prior_geoms: List[QgsGeometry] = [pg for pg in prior if isinstance(pg, QgsGeometry)]
        return float(_max_overlap_ratio(g, prior_geoms) or 0.0)

    def _redundancy_distance_ok(cand: Dict, prior: Sequence[object]) -> bool:
        threshold = float(getattr(params, "max_search_distance", 0.0) or 0.0)
        if threshold <= 0.0:
            return True
        g = cand.get("geom")
        if g is None or g.isEmpty():
            return True
        prior_geoms = [pg for pg in prior if isinstance(pg, QgsGeometry)]
        return _redundancy_far_enough(g, prior_geoms, threshold)

    remaining = float(params.budget_area or 0.0)
    selected: Dict[int, Dict] = {}
    selected_geoms: List[QgsGeometry] = []

    uf = UnionFind()
    for pid, pdata in patches.items():
        uf.find(int(pid))
        uf.size[int(pid)] = float(pdata.get("area_ha", 0.0) or 0.0)
        uf.count[int(pid)] = 1

    def _get_length(cand: Dict) -> float:
        try:
            length = float(cand.get("distance_m", 0.0) or 0.0)
        except Exception:
            length = 0.0
        if length <= 0.0:
            length = float(_get_cost(cand) or 0.0)
        return length

    def _get_patch_size(pid: int) -> float:
        return float(_patch_area_ha(pid) or 0.0)

    picks, base_stats = select_circuit_utility(
        candidates,
        budget=remaining,
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

    for pick in picks:
        cand = pick.candidate
        pids = [int(pid) for pid in _get_patch_ids(cand) if pid is not None]
        if len(pids) < 2:
            continue
        cost = float(_candidate_cost_ha(cand) or 0.0)
        if cost <= 0.0:
            continue
        anchor = int(pids[0])
        for other in pids[1:]:
            uf.union(anchor, int(other))

        cid = len(selected) + 1
        selected[cid] = {
            "geom": clone_geometry(cand["geom"]),
            "patch_ids": set(cand.get("patch_ids", set(pids))),
            "area_ha": float(cand.get("area_ha", 0.0) or 0.0),
            "p1": int(cand.get("patch1")),
            "p2": int(cand.get("patch2")),
            "distance": float(cand.get("distance_m", 1.0) or 1.0),
            "type": pick.corr_type,
            "variant": cand.get("variant"),
            "source": cand.get("source"),
            "utility_score": float(pick.score),
            "overlap_ratio": float(pick.overlap_ratio),
        }

    # Component areas and per-corridor connected area
    comp_area: Dict[int, float] = defaultdict(float)
    comp_count: Dict[int, int] = defaultdict(int)
    for pid in patches:
        root = int(uf.find(int(pid)))
        comp_area[root] += _patch_area_ha(int(pid))
        comp_count[root] += 1
    largest_group_area = max(comp_area.values()) if comp_area else 0.0
    largest_group_patches = max(comp_count.values()) if comp_count else 0

    for data in selected.values():
        root = int(uf.find(int(data.get("p1"))))
        connected_area = float(comp_area.get(root, 0.0))
        data["connected_area_ha"] = connected_area
        area = float(data.get("area_ha", 0.0) or 0.0)
        data["efficiency"] = (connected_area / area) if area else 0.0

    budget_used = float(base_stats.get("budget_used", 0.0) or 0.0)
    stats = {
        "strategy": "circuit_utility",
        "corridors_used": len(selected),
        "budget_used_ha": budget_used,
        "patches_total": len(patches),
        "patches_connected": largest_group_patches,
        "components_remaining": len(comp_area) if comp_area else len(patches),
        "primary_links": int(base_stats.get("primary_links", 0) or 0),
        "redundant_links": int(base_stats.get("redundant_links", 0) or 0),
        "wasteful_links": int(base_stats.get("wasteful_links", 0) or 0),
        "total_connected_area_ha": sum(p.get("area_ha", 0.0) for p in patches.values()),
        "largest_group_area_ha": largest_group_area,
        "largest_group_patches": largest_group_patches,
    }
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
        return {}, {"strategy": "largest_network", "corridors_used": 0, "budget_used_ha": 0.0}

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

    for data in selected.values():
        data["connected_area_ha"] = largest_group_area
        area = float(data.get("area_ha", 0.0) or 0.0)
        data["efficiency"] = (largest_group_area / area) if area else 0.0

    stats = {
        "strategy": "largest_network",
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


def _thicken_corridors(
    corridors: Dict[int, Dict],
    remaining_budget: float,
    params: VectorRunParams,
    patches: Dict[int, Dict],
    max_width_factor: float = 3.0,
) -> float:
    """
    Use remaining budget to widen existing corridors up to a max factor of the base width.
    Prioritize corridors adjacent to the largest patch.
    Returns additional budget used.
    """
    if remaining_budget <= 0 or not corridors:
        return 0.0

    print(f"\n  Thickening corridors with remaining budget: {remaining_budget:.4f} ha")

    try:
        largest_patch = max(patches.items(), key=lambda kv: kv[1].get("area_ha", 0.0))[0]
    except Exception:
        largest_patch = None

    ordered = list(corridors.keys())
    if largest_patch is not None:
        ordered = [cid for cid in ordered if largest_patch in corridors[cid].get("patch_ids", set())] + [
            cid for cid in ordered if largest_patch not in corridors[cid].get("patch_ids", set())
        ]

    budget_used = 0.0
    max_width_factor = max(1.0, float(max_width_factor))

    for cid in ordered:
        if remaining_budget <= 0:
            break
        cdata = corridors.get(cid) or {}
        base_geom = cdata.get("_thicken_base_geom")
        if base_geom is None:
            base_geom = clone_geometry(cdata.get("geom"))
            cdata["_thicken_base_geom"] = base_geom
        if base_geom is None or base_geom.isEmpty():
            continue

        current_factor = float(cdata.get("_thicken_factor", 1.0))
        start_factor = int(math.floor(current_factor))
        for factor in range(start_factor + 1, int(math.floor(max_width_factor)) + 1):
            buffer_dist = (params.min_corridor_width * (factor - 1)) / 2.0
            if buffer_dist <= 0:
                continue
            try:
                thickened_geom = base_geom.buffer(buffer_dist, BUFFER_SEGMENTS)
            except Exception:
                continue
            if thickened_geom is None or thickened_geom.isEmpty():
                continue
            new_area = thickened_geom.area() / 10000.0
            old_area = float(cdata.get("area_ha", new_area))
            added_area = new_area - old_area
            if added_area <= 0:
                cdata["_thicken_factor"] = float(factor)
                continue
            if added_area > remaining_budget:
                print("  Budget exhausted; stopping.")
                return budget_used

            remaining_budget -= added_area
            budget_used += added_area
            cdata["geom"] = thickened_geom
            cdata["area_ha"] = new_area
            cdata["_thicken_factor"] = float(factor)
            if cdata.get("connected_area_ha") is not None:
                conn = float(cdata.get("connected_area_ha") or 0.0)
                cdata["efficiency"] = (conn / new_area) if new_area else 0.0
            print(f"    - Corridor {cid} thickened to {factor:.0f}x width (Cost: {added_area:.4f} ha)")

    if budget_used:
        print(f"  ✓ Finalized thickening. Total extra budget used: {budget_used:.4f} ha")
    return budget_used


def _apply_hybrid_leftover_budget_vector(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    corridors: Dict[int, Dict],
    remaining_budget: float,
    roi_bias: float = ROI_REDUNDANCY_BIAS,
    max_links_per_pair: int = HYBRID_MAX_LINKS_PER_PAIR,
    overlap_reject_ratio: float = HYBRID_OVERLAP_REJECT_RATIO,
    max_search_distance: float = 0.0,
) -> Tuple[float, int, int]:
    """
    Spend remaining budget by comparing low-value connections vs redundancy.
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
            if selected_count_by_pair.get(pk, 0) >= max_links_per_pair:
                continue
            prior_geoms = selected_geoms_by_pair.get(pk, [])
            overlap = _geom_overlap_ratio(geom, prior_geoms)
            if overlap >= float(overlap_reject_ratio):
                continue
            if max_search_distance > 0.0 and not _redundancy_far_enough(
                geom, selected_geoms, float(max_search_distance or 0.0)
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
        if stats.get("strategy") == "largest_network" or stats.get("mode") == "Largest Single Network":
            stats["patches_connected"] = len(patches)
        else:
            stats["patches_connected"] = int(largest_group_patches)
    if "total_connected_area_ha" in stats:
        stats["total_connected_area_ha"] = sum(comp_total_area.values()) if comp_total_area else 0.0

    use_global = stats.get("strategy") == "largest_network" or stats.get("mode") == "Largest Single Network"
    for data in corridors.values():
        if use_global:
            connected_area = largest_group_area
        else:
            root = int(uf.find(int(data.get("p1"))))
            connected_area = float(comp_total_area.get(root, 0.0) or 0.0)
        data["connected_area_ha"] = float(connected_area)
        area = float(data.get("area_ha", 0.0) or 0.0)
        data["efficiency"] = (connected_area / area) if area else 0.0

    if graph_math is None:
        return

    final_corridors = []
    for cid, data in corridors.items():
        p1, p2 = data.get("p1"), data.get("p2")
        if p1 is None or p2 is None:
            continue
        final_corridors.append(
            {"id": cid, "patch1": int(p1), "patch2": int(p2), "distance_m": float(data.get("distance", 1.0) or 1.0)}
        )

    try:
        G_final = graph_math.build_graph_from_corridors(patches, final_corridors)
        if "entropy_total" in stats:
            stats["entropy_total"] = graph_math.calculate_total_entropy(G_final).get("H_total", 0.0)
        if "robustness_rho2" in stats:
            stats["robustness_rho2"] = graph_math.calculate_two_edge_connectivity(G_final)
    except Exception:
        return


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
    area_field = "area_ac" if is_imperial else "area_ha"
    conn_field = "conn_area_ac" if is_imperial else "conn_area_ha"
    area_factor = 2.471053814 if is_imperial else 1.0

    fields = QgsFields()
    fields.append(QgsField("corridor_id", QVariant.Int))
    fields.append(QgsField("patch_ids", QVariant.String))
    fields.append(QgsField(area_field, QVariant.Double))
    fields.append(QgsField(conn_field, QVariant.Double))
    fields.append(QgsField("efficiency", QVariant.Double))
    fields.append(QgsField("multipart", QVariant.Bool))
    fields.append(QgsField("segment_count", QVariant.Int))

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
        multipart = geom.isMultipart()
        segment_count = geom.constGet().numGeometries() if multipart else 1
        feat.setAttributes(
            [
                cid,
                ",".join(map(str, sorted(cdata["patch_ids"]))),
                round(cdata["area_ha"] * area_factor, 4),
                round(cdata["connected_area_ha"] * area_factor, 2),
                round(cdata["efficiency"], 6),
                multipart,
                segment_count,
            ]
        )
        writer.addFeature(feat)

    del writer
    print(f"  ✓ Wrote {len(corridors)} corridors to layer '{layer_name}'")
    return True


def add_layer_to_qgis_from_gpkg(gpkg_path: str, layer_name: str) -> None:
    uri = f"{gpkg_path}|layername={layer_name}"
    layer = QgsVectorLayer(uri, layer_name, "ogr")
    if layer.isValid():
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
) -> Optional[QgsVectorLayer]:
    is_imperial = unit_system == "imperial"
    area_field = "area_ac" if is_imperial else "area_ha"
    conn_field = "conn_area_ac" if is_imperial else "conn_area_ha"
    area_factor = 2.471053814 if is_imperial else 1.0

    layer = QgsVectorLayer(f"Polygon?crs={original_crs.authid()}", layer_name, "memory")
    provider = layer.dataProvider()
    provider.addAttributes(
        [
            QgsField("corridor_id", QVariant.Int),
            QgsField("patch_ids", QVariant.String),
            QgsField(area_field, QVariant.Double),
            QgsField(conn_field, QVariant.Double),
            QgsField("efficiency", QVariant.Double),
            QgsField("multipart", QVariant.Bool),
            QgsField("segment_count", QVariant.Int),
        ]
    )
    layer.updateFields()

    transform = QgsCoordinateTransform(target_crs, original_crs, QgsProject.instance())

    features = []
    for cid, cdata in corridors.items():
        geom = clone_geometry(cdata["geom"])
        geom.transform(transform)
        multipart = geom.isMultipart()
        segment_count = geom.constGet().numGeometries() if multipart else 1
        feat = QgsFeature(layer.fields())
        feat.setGeometry(geom)
        feat.setAttributes(
            [
                cid,
                ",".join(map(str, sorted(cdata["patch_ids"]))),
                round(cdata["area_ha"] * area_factor, 4),
                round(cdata["connected_area_ha"] * area_factor, 2),
                round(cdata["efficiency"], 6),
                multipart,
                segment_count,
            ]
        )
        features.append(feat)

    provider.addFeatures(features)
    layer.updateExtents()
    QgsProject.instance().addMapLayer(layer)
    print(f"  ✓ Added temporary layer '{layer_name}' to QGIS project")
    return layer


def create_memory_layer_from_networks(
    networks: List[Dict],
    layer_name: str,
    target_crs: QgsCoordinateReferenceSystem,
    original_crs: QgsCoordinateReferenceSystem,
    unit_system: str,
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
    if "seed_area_ha" in stats:
        converted["seed_area_display"] = stats["seed_area_ha"] * factor
    if "final_patch_area_ha" in stats:
        converted["final_patch_area_display"] = stats["final_patch_area_ha"] * factor
    if "total_patch_area_ha" in stats:
        converted["total_patch_area_display"] = stats["total_patch_area_ha"] * factor
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
    strategy: str = "circuit_utility",
    temporary: bool = False,
    iface=None,
    progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
    log_cb: Optional[Callable[[str, str], None]] = None,
    ctx: Optional[AnalysisContext] = None,
) -> List[Dict]:
    """Execute the vector corridor analysis for the provided polygon layer."""
    if not isinstance(layer, QgsVectorLayer) or not layer.isValid():
        raise VectorAnalysisError("Selected layer is not a valid vector layer.")

    timings = TimingRecorder()

    params = _to_dataclass(raw_params)
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
    print("LINKSCAPE VECTOR ANALYSIS v23.3")
    print("=" * 70)
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
        patches, spatial_index = load_and_prepare_patches(layer, target_crs, params)
    if not patches:
        raise VectorAnalysisError("No patches found after filtering.")
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

    strategy = (strategy or "circuit_utility").lower()
    if strategy not in ("largest_network", "circuit_utility"):
        strategy = "circuit_utility"

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
    print(f"--- {strategy.replace('_', ' ').upper()} ---")

    if not all_possible:
        raise VectorAnalysisError(_format_no_corridor_reason("Precomputation", len(patches), len(all_possible), params))

    # 5. Optimization
    emit_progress(progress_cb, 60, "Optimizing network...")
    strategy_key = strategy.lower().replace(" ", "_")

    # --- UPDATED DISPATCH LOGIC ---
    opt_label = f"Optimization ({strategy_key})"
    with timings.time_block(opt_label):
        if strategy_key == "largest_network":
            corridors, stats = optimize_circuit_utility_largest_network(patches, all_possible, params)
            layer_name = "Corridors (Largest Single Network)"
        elif strategy_key == "circuit_utility":
            corridors, stats = optimize_circuit_utility(patches, all_possible, params)
            layer_name = "Corridors (Most Connectivity)"
        else:
            raise VectorAnalysisError(f"Unsupported strategy '{strategy_key}'.")

    if not corridors:
        raise VectorAnalysisError(_format_no_corridor_reason("Optimization", len(patches), len(all_possible), params))

    _refresh_vector_connectivity_stats(patches, corridors, stats)

    remaining_budget = float((params.budget_area or 0.0) - float(stats.get("budget_used_ha", 0.0) or 0.0))
    if remaining_budget > 0 and corridors:
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

    if remaining_budget > 0 and corridors:
        with timings.time_block("Thicken corridors"):
            extra_used = _thicken_corridors(
                corridors=corridors,
                remaining_budget=remaining_budget,
                params=params,
                patches=patches,
                max_width_factor=3.0,
            )
        if extra_used:
            stats["budget_used_ha"] = float(stats.get("budget_used_ha", 0.0) or 0.0) + float(extra_used)

    stats = _convert_stats_for_units(stats, unit_system)
    stats["budget_total_display"] = params.budget_area * area_factor

    print("  Preparing outputs...")
    emit_progress(progress_cb, 90, "Writing outputs…")
    with timings.time_block("Write corridor outputs"):
        # Contiguous network areas (patches + corridors dissolved per connected network)
        networks: List[Dict] = []
        try:
            dissolve_tolerance = max(0.0, float(getattr(params, "min_corridor_width", 0.0) or 0.0) * 0.01)
            networks = build_contiguous_network_summaries(patches, corridors, dissolve_tolerance=dissolve_tolerance)
            networks_layer_name = "Contiguous Areas"
            if temporary:
                create_memory_layer_from_networks(networks, networks_layer_name, target_crs, original_crs, unit_system)
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
                add_layer_to_qgis_from_gpkg(output_path, networks_layer_name)
                add_layer_to_qgis_from_gpkg(output_path, layer_name)
        except Exception as net_exc:  # noqa: BLE001
            print(f"  ⚠ Could not write contiguous areas layer: {net_exc}")
            networks = []
        if temporary:
            create_memory_layer_from_corridors(corridors, layer_name, target_crs, original_crs, unit_system)

    # --- LANDSCAPE METRICS REPORT (always written) ---
    landscape_metrics_path = ""
    try:
        strategy_key = (strategy or "circuit_utility").lower()
        if strategy_key not in ("largest_network", "circuit_utility"):
            strategy_key = "circuit_utility"
        analysis_layer_name = f"TerraLink Landscape Metrics ({layer.name()})"
        pixel_size_m = float(getattr(params, "grid_resolution", 50.0) or 50.0)
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
            pixel_size_m=pixel_size_m,
            target_crs=target_crs,
            bounds=bounds,
        )
        pre_mask, _ = _rasterize_geoms_to_mask(
            geoms=patch_geoms,
            pixel_size_m=pixel_size_m,
            target_crs=target_crs,
            bounds=bounds,
        )
        from .analysis_raster import _perform_landscape_analysis  # local import avoids hard coupling at import time

        analysis_lines = _perform_landscape_analysis(
            arr=mask,
            layer_name=analysis_layer_name,
            res_x=float(eff_px),
            res_y=float(eff_px),
            is_metric=True,
            params=None,
            pre_arr=pre_mask,
        )

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
            _add_landscape_metrics_table_layer(analysis_layer_name, analysis_lines)
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
                    "Landscape Metrics (Error)", [f"Landscape analysis failed: {e}"]
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
    strategy_key = (strategy or "circuit_utility").lower()
    if strategy_key not in ("largest_network", "circuit_utility"):
        strategy_key = "circuit_utility"
    if strategy_key == "largest_network":
        strategy_label = "Largest Single Network"
    else:
        strategy_label = "Most Connectivity"
    print(f"Strategy:          {strategy_label}")
    print(f"Corridors created: {stats.get('corridors_used', 0)}")
    print(f"Connections:       {stats.get('connections_made', 0)}")
    print(f"Total connected:   {stats.get('total_connected_area_display', 0):.2f} {area_label}")
    print(f"Largest group:     {stats.get('largest_group_area_display', 0):.2f} {area_label}")
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
        if "primary_links" in stats:
            rows.append(("Primary links", str(stats.get("primary_links", 0))))
        if "redundant_links" in stats:
            rows.append(("Redundant links", str(stats.get("redundant_links", 0))))
        if "entropy_total" in stats:
            rows.append(("Entropy (H_total)", _format_number(stats.get("entropy_total", 0), 4)))
        if "robustness_rho2" in stats:
            rows.append(("Robustness (rho2)", _format_number(stats.get("robustness_rho2", 0), 4)))
        _write_summary_csv(summary_path, rows)
        stats["summary_csv_path"] = summary_path
        print(f"  ✓ Saved vector summary CSV: {summary_path}")
        _add_summary_csv_layer(summary_path, f"TerraLink Vector Summary ({layer.name()})")
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
