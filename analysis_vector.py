"""
Linkscape Corridor Analysis - Vector Workflow (v23.2)
----------------------------------------------------
Runs the single-strategy corridor optimization workflow for polygon
patch datasets.

Updates in v23.2:
- Fixed obstacle avoidance bug where corridors were clipped/pinched.
- Added "Safety Buffer" to RasterNavigator to account for corridor width during routing.
- Implemented 'State-Aware' logic for Most Connectivity strategy.
"""

from __future__ import annotations

import heapq
import math
import os
import time
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Set, Tuple

import numpy as np
from PyQt5.QtCore import QVariant
from qgis.core import (
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


def clone_geometry(geom: QgsGeometry) -> QgsGeometry:
    if geom.isEmpty():
        return QgsGeometry()
    abstract = geom.constGet()
    return QgsGeometry(abstract.clone()) if abstract is not None else QgsGeometry()


class VectorAnalysisError(RuntimeError):
    """Raised when the vector analysis cannot be completed."""


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
    obstacle_layer_id: Optional[str] = None
    obstacle_enabled: bool = False


class UnionFind:
    """Union-Find data structure for tracking connected components."""

    def __init__(self):
        self.parent: Dict[int, int] = {}
        self.size: Dict[int, float] = {}
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

    def get_size(self, x: int) -> float:
        return self.size.get(self.find(x), 0.0)

    def get_count(self, x: int) -> int:
        return self.count.get(self.find(x), 0)


def _to_dataclass(params: Dict) -> VectorRunParams:
    output_name = params.get("output_name") or "linkscape_corridors.gpkg"
    if not output_name.lower().endswith(".gpkg"):
        output_name = f"{output_name}.gpkg"
    return VectorRunParams(
        min_corridor_width=float(params.get("min_corridor_width", 200.0)),
        max_corridor_area=(
            float(params["max_corridor_area"])
            if params.get("max_corridor_area") not in (None, 0, 0.0)
            else None
        ),
        min_patch_size=float(params.get("min_patch_size", 10.0)),
        budget_area=float(params.get("budget_area", 50.0)),
        max_search_distance=float(params.get("max_search_distance", 5000.0)),
        unit_system=str(params.get("unit_system", "metric")),
        output_name=output_name,
        grid_resolution=max(float(params.get("grid_resolution", 50.0)), 1.0),
        obstacle_layer_id=params.get("obstacle_layer_id"),
        obstacle_enabled=bool(params.get("obstacle_enabled", False)),
    )


def get_utm_crs_from_extent(layer: QgsVectorLayer) -> QgsCoordinateReferenceSystem:
    """Determine an appropriate UTM CRS from the layer extent."""
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
    """Load patches, filter by size, and build a spatial index."""
    print("  Loading patches and building spatial index...")
    source_crs = layer.crs()
    transform = QgsCoordinateTransform(source_crs, target_crs, QgsProject.instance())

    patches: Dict[int, Dict] = {}
    patch_id = 1
    filtered_count = 0
    indexed_features: List[QgsFeature] = []

    for feature in layer.getFeatures():
        geom = QgsGeometry(feature.geometry())
        if geom.isEmpty():
            continue
        geom.transform(transform)
        area_ha = geom.area() / 10000.0

        if area_ha < params.min_patch_size:
            filtered_count += 1
            continue

        feat = QgsFeature()
        feat.setGeometry(geom)
        feat.setId(patch_id)
        indexed_features.append(feat)

        patches[patch_id] = {
            "geom": clone_geometry(geom),
            "area_ha": area_ha,
        }
        patch_id += 1

    spatial_index = QgsSpatialIndex(flags=QgsSpatialIndex.FlagStoreFeatureGeometries)
    if indexed_features:
        spatial_index.addFeatures(indexed_features)

    print(f"  ✓ Loaded {len(patches)} patches (filtered {filtered_count} too small)")
    return patches, spatial_index


def _explode_to_polygons(geometry: QgsGeometry) -> List[QgsGeometry]:
    """Return a list of polygon geometries from the supplied geometry."""
    if geometry.isEmpty():
        return []
    geom_type = QgsWkbTypes.geometryType(geometry.wkbType())
    if geom_type != QgsWkbTypes.PolygonGeometry:
        return []
    if geometry.isMultipart():
        polygons = []
        for poly in geometry.asMultiPolygon():
            polygons.append(QgsGeometry.fromPolygonXY(poly))
        return polygons
    return [clone_geometry(geometry)]


def _buffer_line_segment(line_geom: QgsGeometry, width: float) -> QgsGeometry:
    """Buffer a line geometry to create corridor polygon with specified width."""
    return line_geom.buffer(width / 2.0, BUFFER_SEGMENTS)


def _corridor_passes_width(corridor_geom: QgsGeometry, min_width: float) -> bool:
    """Validate that the corridor maintains the minimum width everywhere."""
    if corridor_geom.isEmpty():
        return False

    shrink_distance = max((min_width / 2.0) - max(min_width * 0.01, 0.05), 0.05)
    try:
        shrunk = corridor_geom.buffer(-shrink_distance, BUFFER_SEGMENTS)
    except Exception:
        return False

    if shrunk.isEmpty():
        return False

    return True


def _format_no_corridor_reason(
    stage: str,
    patch_count: int,
    candidate_count: int,
    params: VectorRunParams,
) -> str:
    """Produce a descriptive error string when no corridors can be generated."""
    return (
        f"{stage}: no feasible corridors could be generated.\n"
        f"- Patches meeting criteria: {patch_count}\n"
        f"- Candidate corridors generated: {candidate_count}\n"
        f"- Max search distance: {params.max_search_distance:.2f}\n"
        f"- Min corridor width: {params.min_corridor_width:.2f}\n"
        f"- Max corridor area: {params.max_corridor_area or 'None (no limit)'}\n"
        "Try increasing the search distance, lowering the minimum corridor width/area, "
        "or simplifying the patch layer."
    )


def _create_corridor_geometry(
    waypoints: List[QgsPointXY],
    source_geom: QgsGeometry,
    target_geom: QgsGeometry,
    params: VectorRunParams,
    obstacle_geoms: Optional[List[QgsGeometry]] = None,
    smooth_iterations: int = 0,
) -> Optional[QgsGeometry]:
    if not waypoints or len(waypoints) < 2:
        return None

    corridor_line = QgsGeometry.fromPolylineXY([QgsPointXY(pt) for pt in waypoints])
    if smooth_iterations > 0:
        try:
            smoothed = corridor_line.smooth(smooth_iterations)
            if smoothed and not smoothed.isEmpty():
                corridor_line = smoothed
        except Exception:
            pass
    corridor_geom = _buffer_line_segment(corridor_line, params.min_corridor_width)
    corridor_geom = corridor_geom.difference(source_geom)
    corridor_geom = corridor_geom.difference(target_geom)
    if obstacle_geoms:
        for obstacle in obstacle_geoms:
            corridor_geom = corridor_geom.difference(obstacle)
            if corridor_geom.isEmpty():
                return None
    corridor_geom = corridor_geom.makeValid()
    if corridor_geom.isEmpty():
        return None

    if not _corridor_passes_width(corridor_geom, params.min_corridor_width):
        return None

    return corridor_geom


class RasterNavigator:
    """Hybrid raster navigator that routes corridors around obstacles."""

    def __init__(
        self,
        patches: Dict[int, Dict],
        obstacle_layer: QgsVectorLayer,
        target_crs: QgsCoordinateReferenceSystem,
        params: VectorRunParams,
    ):
        if obstacle_layer is None or QgsWkbTypes.geometryType(obstacle_layer.wkbType()) != QgsWkbTypes.PolygonGeometry:
            raise VectorAnalysisError("Select a polygon obstacle layer for vector obstacle avoidance.")

        self.resolution = max(params.grid_resolution, 1.0)
        self.obstacle_geoms: List[QgsGeometry] = []

        extent: Optional[QgsRectangle] = None
        for patch in patches.values():
            bbox = patch["geom"].boundingBox()
            if extent is None:
                extent = QgsRectangle(bbox)
            else:
                extent.combineExtentWith(bbox)

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
            raise VectorAnalysisError("Unable to determine extent for obstacle avoidance routing.")

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

        # FIX: Inflate obstacles by half the corridor width to keep the centerline safe.
        # This prevents the final corridor buffer from grazing the obstacle and being clipped.
        safety_buffer = params.min_corridor_width / 2.0

        for geom in self.obstacle_geoms:
            # We buffer the obstacle geometry solely for the mask creation
            # Using 4 segments is sufficient for raster mask accuracy
            buffered_mask = geom.buffer(safety_buffer, 4)
            self._burn_geometry(buffered_mask)

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
                    # Check center point
                    if geom.contains(QgsPointXY(x, y)):
                        self.passable[row, col] = False
                except Exception:
                    pass

    def find_path(self, start_point: QgsPointXY, end_point: QgsPointXY) -> Optional[List[QgsPointXY]]:
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
    moves = [
        (-1, 0),
        (1, 0),
        (0, -1),
        (0, 1),
        (-1, -1),
        (-1, 1),
        (1, -1),
        (1, 1),
    ]

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
    progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
    progress_start: int = 30,
    progress_end: int = 55,
    navigator: Optional[RasterNavigator] = None,
) -> List[Dict]:
    """
    Find all possible corridors between patch pairs within the distance constraint.
    Uses either buffered straight lines or the raster navigator when enabled.
    """
    print("  Finding all possible corridors...")
    all_corridors: List[Dict] = []
    processed_pairs: Set[frozenset] = set()
    total = len(patches) or 1

    for idx, (pid1, pdata1) in enumerate(patches.items(), start=1):
        if idx % 10 == 0:
            print(f"    Analyzing patch {idx}/{len(patches)}...", end="\r")
        if progress_cb is not None:
            span = max(progress_end - progress_start, 1)
            progress_value = progress_start + ((idx - 1) / total) * span
            _emit_progress(progress_cb, progress_value, "Building corridor candidates…")

        geom1 = pdata1["geom"]
        rect = geom1.boundingBox()
        rect.grow(params.max_search_distance)
        candidate_ids = spatial_index.intersects(rect)

        for pid2 in candidate_ids:
            if pid1 >= pid2:
                continue

            pair = frozenset({pid1, pid2})
            if pair in processed_pairs:
                continue

            pdata2 = patches.get(pid2)
            if not pdata2:
                continue
            geom2 = pdata2["geom"]

            distance = geom1.distance(geom2)
            if distance > params.max_search_distance:
                continue

            p1 = geom1.nearestPoint(geom2).asPoint()
            p2 = geom2.nearestPoint(geom1).asPoint()
            if p1.isEmpty() or p2.isEmpty():
                continue

            if navigator:
                path_points = navigator.find_path(QgsPointXY(p1), QgsPointXY(p2))
                if not path_points:
                    continue
                corridor_geom = _create_corridor_geometry(
                    path_points,
                    geom1,
                    geom2,
                    params,
                    obstacle_geoms=navigator.obstacle_geoms,
                    smooth_iterations=3,
                )
            else:
                corridor_geom = _create_corridor_geometry(
                    [QgsPointXY(p1), QgsPointXY(p2)],
                    geom1,
                    geom2,
                    params,
                )

            if corridor_geom is None:
                continue

            corridor_area_ha = corridor_geom.area() / 10000.0
            if corridor_area_ha <= 0:
                continue
            if params.max_corridor_area is not None and corridor_area_ha > params.max_corridor_area:
                continue

            all_corridors.append(
                {
                    "patch1": pid1,
                    "patch2": pid2,
                    "geom": clone_geometry(corridor_geom),
                    "area_ha": corridor_area_ha,
                }
            )
            processed_pairs.add(pair)

    _emit_progress(progress_cb, progress_end, "Candidate corridors ready.")
    print(f"\n  ✓ Found {len(all_corridors)} possible corridors")
    return all_corridors


def find_shortest_corridor_from_group(
    start_patches: Set[int],
    all_patches: Dict[int, Dict],
    spatial_index: QgsSpatialIndex,
    params: VectorRunParams,
    navigator: Optional[RasterNavigator] = None,
) -> Optional[Tuple]:
    """
    Find the smallest-area corridor connecting the current group to any other patch.
    """
    if not start_patches:
        return None

    start_geoms = [all_patches[pid]["geom"] for pid in start_patches]
    start_union = QgsGeometry.unaryUnion(start_geoms)

    rect = start_union.boundingBox()
    rect.grow(params.max_search_distance)
    candidate_ids = spatial_index.intersects(rect)

    best_corridor = None
    best_area = float("inf")
    best_target_id = None

    for candidate_id in candidate_ids:
        if candidate_id in start_patches:
            continue

        target = all_patches.get(candidate_id)
        if not target:
            continue

        target_geom = target["geom"]
        if start_union.distance(target_geom) > params.max_search_distance:
            continue

        p1 = start_union.nearestPoint(target_geom).asPoint()
        p2 = target_geom.nearestPoint(start_union).asPoint()
        if p1.isEmpty() or p2.isEmpty():
            continue

        if navigator:
            path_points = navigator.find_path(QgsPointXY(p1), QgsPointXY(p2))
            if not path_points:
                continue
            geom_candidate = _create_corridor_geometry(
                path_points,
                start_union,
                target_geom,
                params,
                obstacle_geoms=navigator.obstacle_geoms,
                smooth_iterations=3,
            )
        else:
            geom_candidate = _create_corridor_geometry(
                [QgsPointXY(p1), QgsPointXY(p2)],
                start_union,
                target_geom,
                params,
            )

        if geom_candidate is None:
            continue

        corridor_area_ha = geom_candidate.area() / 10000.0
        if corridor_area_ha <= 0:
            continue
        if params.max_corridor_area is not None and corridor_area_ha > params.max_corridor_area:
            continue

        if corridor_area_ha < best_area:
            best_corridor = geom_candidate
            best_area = corridor_area_ha
            best_target_id = candidate_id

    if best_corridor is None or best_target_id is None:
        return None
    return best_corridor, {best_target_id}, best_area


def optimize_most_connectivity(
    patches: Dict[int, Dict],
    candidates: List[Dict],
    params: VectorRunParams,
) -> Tuple[Dict[int, Dict], Dict]:
    """Strategy 1: Most Connectivity - maximise total connected area."""
    print("  Strategy: MOST CONNECTIVITY (Dynamic State-Aware)")
    uf = UnionFind()
    component_area: Dict[int, float] = {}
    component_active: Dict[int, bool] = {}
    for pid, pdata in patches.items():
        uf.find(pid)
        uf.size[pid] = pdata["area_ha"]
        uf.count[pid] = 1
        component_area[pid] = pdata["area_ha"]
        component_active[pid] = False

    selected: Dict[int, Dict] = {}
    remaining_budget = params.budget_area
    connections_made = 0

    def _marginal_gain(root_a: int, root_b: int) -> float:
        active_a = component_active.get(root_a, False)
        active_b = component_active.get(root_b, False)
        area_a = component_area.get(root_a, 0.0)
        area_b = component_area.get(root_b, 0.0)
        # Rescue: Both isolated -> Benefit = A + B
        if not active_a and not active_b:
            return area_a + area_b
        # Merger: Both connected -> Benefit = 0
        if active_a and active_b:
            return 0.0
        # Expansion: One connected -> Benefit = Isolated one
        return area_b if active_a else area_a

    while remaining_budget > 0:
        best_candidate: Optional[Dict] = None
        best_efficiency = 0.0
        best_gain = 0.0

        for cand in candidates:
            cost = cand["area_ha"]
            if cost > remaining_budget:
                continue
            p1, p2 = cand["patch1"], cand["patch2"]
            root1, root2 = uf.find(p1), uf.find(p2)
            if root1 == root2:
                continue
            
            # Recalculate gain based on current connectivity state
            gain = _marginal_gain(root1, root2)
            if gain <= 0:
                continue
                
            efficiency = gain / cost if cost > 0 else float("inf")
            if (
                best_candidate is None
                or efficiency > best_efficiency
                or (efficiency == best_efficiency and gain > best_gain)
            ):
                best_candidate = cand
                best_efficiency = efficiency
                best_gain = gain

        if best_candidate is None:
            break

        cost = best_candidate["area_ha"]
        p1, p2 = best_candidate["patch1"], best_candidate["patch2"]
        root1, root2 = uf.find(p1), uf.find(p2)
        area1 = component_area.get(root1, patches[p1]["area_ha"])
        area2 = component_area.get(root2, patches[p2]["area_ha"])

        new_root = uf.union(p1, p2)
        uf.size[new_root] += cost
        component_area[new_root] = area1 + area2
        component_active[new_root] = True
        if new_root != root1:
            component_area.pop(root1, None)
            component_active.pop(root1, None)
        if new_root != root2:
            component_area.pop(root2, None)
            component_active.pop(root2, None)

        connections_made += 1
        remaining_budget -= cost

        connected_area = uf.get_size(new_root)
        efficiency = best_gain / cost if cost > 0 else float("inf")
        selected[len(selected) + 1] = {
            "geom": best_candidate["geom"],
            "patch_ids": {p1, p2},
            "area_ha": cost,
            "connected_area_ha": connected_area,
            "efficiency": efficiency,
        }

    active_roots = [root for root, active in component_active.items() if active]
    root_areas = {root: uf.get_size(root) for root in active_roots}
    root_counts = {root: uf.get_count(root) for root in active_roots}
    stats = {
        "strategy": "most_connectivity",
        "corridors_used": len(selected),
        "connections_made": connections_made,
        "budget_used_ha": params.budget_area - remaining_budget,
        "total_connected_area_ha": sum(root_areas.values()),
        "groups_created": len(root_areas),
        "largest_group_area_ha": max(root_areas.values()) if root_areas else 0,
        "patches_connected": sum(root_counts.values()),
    }
    print(f"  ✓ Selected {len(selected)} efficient corridors")
    return selected, stats


def optimize_largest_patch(
    patches: Dict[int, Dict],
    spatial_index: QgsSpatialIndex,
    params: VectorRunParams,
    navigator: Optional[RasterNavigator] = None,
) -> Tuple[Dict[int, Dict], Dict]:
    """Strategy 3: Largest Patch - create a single largest connected patch."""
    print("  Strategy: LARGEST PATCH")
    print("    Testing different seed patches...")

    best_result = {"corridors": {}, "final_area": 0.0, "seed_id": None}
    sorted_patches = sorted(patches.keys(), key=lambda k: patches[k]["area_ha"], reverse=True)

    for i, seed_id in enumerate(sorted_patches[: min(20, len(sorted_patches))]):
        if (i + 1) % 5 == 0:
            print(f"      Testing seed {i + 1}...", end="\r")

        current_patches = {seed_id}
        current_area = patches[seed_id]["area_ha"]
        remaining_budget = params.budget_area
        sim_corridors: Dict[int, Dict] = {}

        while remaining_budget > 0:
            result = find_shortest_corridor_from_group(
                current_patches,
                patches,
                spatial_index,
                params,
                navigator=navigator,
            )
            if result is None:
                break

            corr_geom, targets, corr_area = result
            if corr_area > remaining_budget:
                break

            new_area = sum(patches[pid]["area_ha"] for pid in targets if pid not in current_patches)
            current_area += corr_area + new_area

            efficiency = corr_area / current_area if current_area > 0 else 0

            sim_corridors[len(sim_corridors) + 1] = {
                "geom": corr_geom,
                "patch_ids": set(current_patches | targets),
                "area_ha": corr_area,
                "connected_area_ha": current_area,
                "efficiency": efficiency,
            }

            remaining_budget -= corr_area
            current_patches.update(targets)

        if current_area > best_result["final_area"]:
            best_result = {
                "corridors": sim_corridors,
                "final_area": current_area,
                "seed_id": seed_id,
            }

    print(f"\n  ✓ Found optimal seed patch {best_result['seed_id']}")

    if not best_result["corridors"]:
        return {}, {"strategy": "largest_patch", "corridors_used": 0}

    stats = {
        "strategy": "largest_patch",
        "seed_id": best_result["seed_id"],
        "seed_area_ha": patches[best_result["seed_id"]]["area_ha"],
        "final_patch_area_ha": best_result["final_area"],
        "corridors_used": len(best_result["corridors"]),
        "budget_used_ha": sum(c["area_ha"] for c in best_result["corridors"].values()),
        "patches_merged": len(set.union(*(c["patch_ids"] for c in best_result["corridors"].values()))),
        "groups_created": 1,
    }
    return best_result["corridors"], stats


def write_corridors_layer_to_gpkg(
    corridors: Dict[int, Dict],
    output_path: str,
    layer_name: str,
    target_crs: QgsCoordinateReferenceSystem,
    original_crs: QgsCoordinateReferenceSystem,
    unit_system: str,
    overwrite_file: bool = False,
) -> bool:
    """Write corridor polygons to a specific layer within a GeoPackage."""
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
                round(cdata["area_ha"] * area_factor, 4),
                round(cdata["connected_area_ha"] * area_factor, 2),
                round(cdata["efficiency"], 6),
            ]
        )
        writer.addFeature(feat)

    del writer
    print(f"  ✓ Wrote {len(corridors)} corridors to layer '{layer_name}'")
    return True


def add_layer_to_qgis_from_gpkg(gpkg_path: str, layer_name: str) -> None:
    """Add a specific layer from a GeoPackage into the project."""
    uri = f"{gpkg_path}|layername={layer_name}"
    layer = QgsVectorLayer(uri, layer_name, "ogr")
    if layer.isValid():
        QgsProject.instance().addMapLayer(layer)
        print(f"  ✓ Added '{layer_name}' to QGIS project")
    else:
        print(f"  ✗ Could not add '{layer_name}' from {gpkg_path}")


def create_memory_layer_from_corridors(
    corridors: Dict[int, Dict],
    layer_name: str,
    target_crs: QgsCoordinateReferenceSystem,
    original_crs: QgsCoordinateReferenceSystem,
    unit_system: str,
) -> Optional[QgsVectorLayer]:
    """Build an in-memory layer for corridor geometries."""
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
                round(cdata["area_ha"] * area_factor, 4),
                round(cdata["connected_area_ha"] * area_factor, 2),
                round(cdata["efficiency"], 6),
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


def run_vector_analysis(
    layer: QgsVectorLayer,
    output_dir: str,
    raw_params: Dict,
    strategy: str = "most_connectivity",
    temporary: bool = False,
    iface=None,
    progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
) -> List[Dict]:
    """Execute the vector corridor analysis for the provided polygon layer."""
    if not isinstance(layer, QgsVectorLayer) or not layer.isValid():
        raise VectorAnalysisError("Selected layer is not a valid vector layer.")

    params = _to_dataclass(raw_params)
    unit_system = params.unit_system
    area_factor = 2.471053814 if unit_system == "imperial" else 1.0
    area_label = "ac" if unit_system == "imperial" else "ha"

    if temporary:
        output_path = ""
    else:
        out_dir = output_dir or os.path.dirname(layer.source())
        os.makedirs(out_dir, exist_ok=True)
        output_path = os.path.join(out_dir, params.output_name)

    overall_start = time.time()
    print("=" * 70)
    print("LINKSCAPE VECTOR ANALYSIS v23.2")
    print("=" * 70)
    print("\n1. Loading vector layer...")
    print(f"  ✓ Using layer: {layer.name()} ({layer.featureCount()} features)")
    _emit_progress(progress_cb, 5, "Loading vector layer…")

    original_crs = layer.crs()
    print(f"  CRS: {original_crs.authid()}")

    print("\n2. Determining analysis CRS...")
    target_crs = get_utm_crs_from_extent(layer)
    print(f"  ✓ Using {target_crs.authid()} for measurements")
    _emit_progress(progress_cb, 15, "Preparing data…")

    print("\n3. Loading patches...")
    patches, spatial_index = load_and_prepare_patches(layer, target_crs, params)
    if not patches:
        raise VectorAnalysisError("No patches found after filtering.")
    _emit_progress(progress_cb, 35, "Searching for corridor candidates…")

    navigator: Optional[RasterNavigator] = None
    if params.obstacle_enabled and params.obstacle_layer_id:
        obstacle_layer = QgsProject.instance().mapLayer(params.obstacle_layer_id)
        if not isinstance(obstacle_layer, QgsVectorLayer) or not obstacle_layer.isValid():
            print("  ⚠ Selected obstacle layer is unavailable; continuing without obstacle avoidance.")
        else:
            try:
                navigator = RasterNavigator(patches, obstacle_layer, target_crs, params)
                print(
                    f"  ✓ Raster navigator grid: {navigator.cols} × {navigator.rows} cells "
                    f"@ {navigator.resolution:.1f} units"
                )
            except VectorAnalysisError as exc:
                print(f"  ⚠ Obstacle avoidance disabled: {exc}")
                navigator = None

    print("\n4. Precomputing candidate corridors...")
    all_possible = find_all_possible_corridors(
        patches,
        spatial_index,
        params,
        progress_cb=progress_cb,
        progress_start=35,
        progress_end=60,
        navigator=navigator,
    )

    strategy = (strategy or "most_connectivity").lower()
    strategy_map = {
        "most_connectivity": (
            "candidates",
            optimize_most_connectivity,
            "Corridors (Most Connectivity)",
        ),
        "largest_patch": (
            "spatial",
            optimize_largest_patch,
            "Corridors (Largest Patch)",
        ),
    }

    if strategy not in strategy_map:
        raise VectorAnalysisError(f"Unsupported strategy '{strategy}'.")

    mode, optimize_func, layer_name = strategy_map[strategy]

    print("\n5. Running optimization...")
    print("=" * 70)
    print(f"--- {strategy.replace('_', ' ').upper()} ---")
    _emit_progress(progress_cb, 65, "Running optimization…")

    if mode == "candidates" and not all_possible:
        raise VectorAnalysisError(
            _format_no_corridor_reason(
                "Precomputation",
                len(patches),
                len(all_possible),
                params,
            )
        )

    if mode == "candidates":
        corridors, stats = optimize_func(patches, all_possible, params)
    else:
        corridors, stats = optimize_func(
            patches,
            spatial_index,
            params,
            navigator=navigator,
        )

    if not corridors:
        raise VectorAnalysisError(
            _format_no_corridor_reason(
                "Optimization",
                len(patches),
                len(all_possible),
                params,
            )
        )

    stats = _convert_stats_for_units(stats, unit_system)
    stats["budget_total_display"] = params.budget_area * area_factor

    print("  Preparing outputs...")
    _emit_progress(progress_cb, 90, "Writing outputs…")
    if temporary:
        create_memory_layer_from_corridors(corridors, layer_name, target_crs, original_crs, unit_system)
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
        add_layer_to_qgis_from_gpkg(output_path, layer_name)

    elapsed = time.time() - overall_start
    _emit_progress(progress_cb, 100, "Vector analysis complete.")

    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"Strategy:          {strategy.replace('_', ' ').title()}")
    print(f"Corridors created: {stats.get('corridors_used', 0)}")
    if strategy == "largest_patch":
        print(f"Seed patch:        {stats.get('seed_id')}")
        print(f"Seed area:         {stats.get('seed_area_display', 0):.2f} {area_label}")
        print(f"Final patch size:  {stats.get('final_patch_area_display', 0):.2f} {area_label}")
        print(f"Patches merged:    {stats.get('patches_merged', 0)}")
    else:
        print(f"Connections:       {stats.get('connections_made', 0)}")
        print(f"Total connected:   {stats.get('total_connected_area_display', 0):.2f} {area_label}")
        print(f"Largest group:     {stats.get('largest_group_area_display', 0):.2f} {area_label}")
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

    return [
        {
            "strategy": strategy,
            "stats": stats,
            "output_path": output_path if not temporary else "",
            "layer_name": layer_name,
        }
    ]