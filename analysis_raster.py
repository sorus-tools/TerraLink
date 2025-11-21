"""
Linkscape Corridor Analysis - Raster Workflow (v23.0)
-----------------------------------------------------
Runs the selected raster optimization workflow for corridor analysis.

The logic is adapted from the standalone raster script and packaged so it can
be invoked from the QGIS plugin with user-supplied parameters.
"""

from __future__ import annotations

import heapq
import math
import os
import tempfile
import time
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Set, Tuple

import numpy as np
try:
    from osgeo import gdal
except ImportError:  # pragma: no cover
    gdal = None  # type: ignore

try:
    from qgis.core import QgsProject, QgsRasterLayer
except ImportError:  # pragma: no cover
    QgsProject = None  # type: ignore
    QgsRasterLayer = None  # type: ignore

# Optional OpenCV import for faster connected component labeling
try:
    import cv2

    HAS_CV2 = True
except ImportError:
    HAS_CV2 = False

GTIFF_OPTIONS = ["COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"]

PIXEL_COUNT_WARNING_THRESHOLD = 40_000_000  # warn when raster exceeds ~40 million pixels
PIXEL_SIZE_WARNING_THRESHOLD = 10.0  # warn when pixel size < 10 map units

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
    """Create a boolean mask for obstacle pixels corridors must avoid."""
    if not params.obstacle_enabled:
        return patch_mask.astype(bool)

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
    mask |= patch_mask.astype(bool)
    return mask


def find_shortest_corridor(
    start_patches: Set[int],
    labels: np.ndarray,
    habitat: np.ndarray,
    max_width: int,
    connectivity: int,
    obstacle_mask: Optional[np.ndarray] = None,
    passable_mask: Optional[np.ndarray] = None,
) -> List[Tuple[frozenset, int, float]]:
    """
    Dijkstra search to find shortest corridors connecting start_patches to other patches.
    Returns a list of (path_pixels, target_patch, length).
    """
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

    start_positions: List[Tuple[int, int]] = []
    for i in range(rows):
        for j in range(cols):
            if obstacle_mask is not None and obstacle_mask[i, j]:
                continue
            if passable_mask is not None and not passable_mask[i, j]:
                continue
            if not habitat[i, j]:
                for di, dj in moves:
                    ni, nj = i + di, j + dj
                    if 0 <= ni < rows and 0 <= nj < cols and labels[ni, nj] in start_patches:
                        start_positions.append((i, j))
                        break

    if not start_positions:
        return []

    heap: List[Tuple[float, int, int, frozenset]] = []
    best_cost: Dict[Tuple[int, int], float] = {}
    results: List[Tuple[frozenset, int, float]] = []
    visited_targets: Set[int] = set()
    for r, c in start_positions:
        if obstacle_mask is not None and obstacle_mask[r, c]:
            continue
        if passable_mask is not None and not passable_mask[r, c]:
            continue
        path = frozenset({(r, c)})
        heapq.heappush(heap, (0.0, r, c, path))
        best_cost[(r, c)] = 0.0

    while heap:
        cost, r, c, path = heapq.heappop(heap)
        if cost > max_width:
            continue

        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                lbl = labels[nr, nc]
                if lbl > 0 and lbl not in start_patches and lbl not in visited_targets:
                    visited_targets.add(lbl)
                    results.append((path, lbl, cost))

        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                if obstacle_mask is not None and obstacle_mask[nr, nc]:
                    continue
                if passable_mask is not None and not passable_mask[nr, nc]:
                    continue
                if habitat[nr, nc]:
                    continue
                move_cost = math.sqrt(2) if dr != 0 and dc != 0 else 1.0
                new_cost = cost + move_cost
                if new_cost > max_width:
                    continue
                prev_best = best_cost.get((nr, nc))
                if prev_best is not None and prev_best <= new_cost:
                    continue
                best_cost[(nr, nc)] = new_cost
                heapq.heappush(heap, (new_cost, nr, nc, path | frozenset({(nr, nc)})))

    return results


def find_all_possible_corridors(
    labels: np.ndarray,
    habitat: np.ndarray,
    patch_sizes: Dict[int, int],
    max_width: int,
    max_area: Optional[int],
    connectivity: int,
    obstacle_mask: Optional[np.ndarray] = None,
    passable_mask: Optional[np.ndarray] = None,
    progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
    progress_start: int = 45,
    progress_end: int = 75,
) -> List[Dict]:
    """Find all possible corridors between patch pairs."""
    print("  Finding all possible corridors...")
    all_corridors: List[Dict] = []
    processed_pairs: Set[frozenset] = set()

    unique_patches = [p for p in patch_sizes.keys() if p > 0]
    total = len(unique_patches) or 1
    span = max(progress_end - progress_start, 1)
    for idx, patch_id in enumerate(unique_patches):
        if (idx + 1) % 10 == 0:
            print(f"    Analyzing patch {idx + 1}/{len(unique_patches)}...", end="\r")
        if progress_cb is not None:
            pre_value = progress_start + ((idx + 0.25) / total) * span
            _emit_progress(
                progress_cb,
                pre_value,
                f"Analyzing patch {idx + 1}/{total}…",
            )

        results = find_shortest_corridor(
            {patch_id},
            labels,
            habitat,
            max_width,
            connectivity,
            obstacle_mask=obstacle_mask,
            passable_mask=passable_mask,
        )
        if not results:
            continue

        for path_pixels, target_id, path_len in results:
            pair = frozenset({patch_id, target_id})
            if pair in processed_pairs:
                continue
            processed_pairs.add(pair)

            if max_area is not None and path_len > max_area:
                continue

            all_corridors.append(
                {
                    "patch1": patch_id,
                    "patch2": target_id,
                    "pixels": path_pixels,
                    "length": path_len,
                }
            )
        if progress_cb is not None:
            post_value = progress_start + ((idx + 1) / total) * span
            _emit_progress(
                progress_cb,
                post_value,
                f"Finished patch {idx + 1}/{total}",
            )

    _emit_progress(progress_cb, progress_end, "Corridor candidates ready.")
    print(f"\n  ✓ Found {len(all_corridors)} possible corridors")
    return all_corridors


def optimize_most_connectivity(
    candidates: List[Dict], patch_sizes: Dict[int, int], budget: int
) -> Tuple[Dict, Dict]:
    """Strategy 1: Most Connectivity - maximize total connected area."""
    print("  Strategy: MOST CONNECTIVITY")

    candidates_copy = [c.copy() for c in candidates]
    for c in candidates_copy:
        benefit = patch_sizes[c["patch1"]] + patch_sizes[c["patch2"]]
        c["efficiency"] = benefit / c["length"] if c["length"] > 0 else float("inf")

    candidates_copy.sort(key=lambda x: x["efficiency"], reverse=True)

    uf = UnionFind()
    for pid, size in patch_sizes.items():
        uf.find(pid)
        uf.size[pid] = size
        uf.count[pid] = 1

    selected: Dict[int, Dict] = {}
    remaining_budget = budget
    connections_made = 0

    for cand in candidates_copy:
        if cand["length"] > remaining_budget:
            continue

        p1, p2 = cand["patch1"], cand["patch2"]
        if uf.find(p1) == uf.find(p2):
            continue

        root = uf.union(p1, p2)
        uf.size[root] += cand["length"]
        connections_made += 1

        selected[len(selected) + 1] = {
            "pixels": cand["pixels"],
            "patch_ids": {p1, p2},
            "length": cand["length"],
            "connected_size": uf.get_size(root),
        }

        remaining_budget -= cand["length"]

    root_sizes = {uf.find(pid): uf.get_size(pid) for pid in patch_sizes.keys()}
    root_counts = {uf.find(pid): uf.get_count(pid) for pid in patch_sizes.keys()}

    return selected, {
        "strategy": "most_connectivity",
        "corridors_used": len(selected),
        "connections_made": connections_made,
        "budget_used": budget - remaining_budget,
        "total_connected_size": sum(root_sizes.values()),
        "groups_created": len(root_sizes),
        "largest_group_size": max(root_sizes.values()) if root_sizes else 0,
        "largest_group_patches": max(root_counts.values()) if root_counts else 0,
        "patches_connected": max(root_counts.values()) if root_counts else 0,
    }


def optimize_largest_patch(
    candidates: List[Dict], patch_sizes: Dict[int, int], budget: int
) -> Tuple[Dict, Dict]:
    """Strategy 3: Largest Patch - expand the largest component via priority queue."""
    print("  Strategy: LARGEST PATCH (priority growth)")
    if not candidates:
        return {}, {"strategy": "largest_patch", "corridors_used": 0}

    adjacency: Dict[int, List[Dict]] = {}
    for cand in candidates:
        adjacency.setdefault(cand["patch1"], []).append(cand)
        adjacency.setdefault(cand["patch2"], []).append(cand)

    sorted_patch_ids = sorted(patch_sizes.keys(), key=lambda k: patch_sizes[k], reverse=True)

    print(f"  Testing {min(20, len(sorted_patch_ids))} seed patches...")
    best_result = {"corridors": {}, "final_size": 0, "seed_id": None, "budget_used": 0}

    for i, seed_id in enumerate(sorted_patch_ids[:20]):
        if (i + 1) % 5 == 0:
            print(f"    Seed {i + 1}/20...", end="\r")

        component: Set[int] = {seed_id}
        component_size = patch_sizes.get(seed_id, 0)
        remaining_budget = float(budget)
        sim_corridors: Dict[int, Dict] = {}
        counter = 0
        pq: List[Tuple[float, float, int, Dict]] = []

        def enqueue_neighbors(patch_id: int) -> None:
            nonlocal counter, component_size
            for cand in adjacency.get(patch_id, []):
                p1, p2 = cand["patch1"], cand["patch2"]
                if p1 in component and p2 in component:
                    continue
                if p1 == patch_id and p2 not in component:
                    target = p2
                elif p2 == patch_id and p1 not in component:
                    target = p1
                else:
                    continue
                potential_size = component_size + patch_sizes.get(target, 0) + cand["length"]
                heapq.heappush(pq, (-potential_size, cand["length"], counter, cand))
                counter += 1

        enqueue_neighbors(seed_id)

        while pq and remaining_budget > 0:
            neg_potential, cost, _, cand = heapq.heappop(pq)
            p1, p2 = cand["patch1"], cand["patch2"]
            if p1 in component and p2 in component:
                continue
            if p1 in component:
                target = p2
            elif p2 in component:
                target = p1
            else:
                continue
            if target in component:
                continue
            if cost > remaining_budget:
                continue

            remaining_budget -= cost
            component.add(target)
            component_size += patch_sizes.get(target, 0) + cost
            sim_corridors[len(sim_corridors) + 1] = {
                "pixels": cand["pixels"],
                "patch_ids": {p1, p2},
                "length": cost,
                "connected_size": component_size,
            }
            enqueue_neighbors(target)

        if component_size > best_result["final_size"]:
            best_result = {
                "corridors": sim_corridors,
                "final_size": component_size,
                "seed_id": seed_id,
                "budget_used": budget - remaining_budget,
                "patch_count": len(component),
            }

    print(f"\n  ✓ Best: Patch {best_result['seed_id']} -> {best_result['final_size']:,} px")

    if not best_result["corridors"]:
        return {}, {"strategy": "largest_patch", "corridors_used": 0}

    selected: Dict[int, Dict] = {}
    for i, corr in best_result["corridors"].items():
        selected[i] = {
            "pixels": corr["pixels"],
            "patch_ids": corr["patch_ids"],
            "length": corr["length"],
            "connected_size": best_result["final_size"],
        }

    return selected, {
        "strategy": "largest_patch",
        "seed_id": best_result["seed_id"],
        "final_patch_size": best_result["final_size"],
        "corridors_used": len(selected),
        "budget_used": best_result["budget_used"],
        "groups_created": 1,
        "patches_connected": best_result.get("patch_count", len(selected) + 1),
        "largest_group_size": best_result["final_size"],
    }


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
        return base_passable

    offsets = _corridor_offsets(min_corridor_width)
    # Ignore habitat when evaluating clearance so corridors can still touch patches.
    true_obstacles = obstacle_bool & (~habitat_bool)
    clearance_space = ~true_obstacles
    clearance_ok = _erode_mask(clearance_space, offsets)
    return base_passable & clearance_ok


def create_output_raster(
    labels: np.ndarray,
    corridors: Dict[int, Dict],
    min_corridor_width: int,
    obstacle_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Create an output raster with corridors marked by connected size."""
    output = np.zeros_like(labels, dtype=np.int32)
    offsets = _corridor_offsets(min_corridor_width)
    rows, cols = labels.shape
    use_mask = obstacle_mask is not None

    for corridor_data in corridors.values():
        score = corridor_data["connected_size"]
        for r, c in corridor_data["pixels"]:
            for dr, dc in offsets:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    if use_mask and obstacle_mask[nr, nc]:
                        continue
                    if score > output[nr, nc]:
                        output[nr, nc] = score
    return output


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
        min_patch_size=int(params.get("min_patch_size", 0)),
        budget_pixels=int(params.get("budget_pixels", 0)),
        max_search_distance=int(params.get("max_search_distance", 50)),
        max_corridor_area=(
            int(params["max_corridor_area"]) if params.get("max_corridor_area") is not None else None
        ),
        min_corridor_width=int(params.get("min_corridor_width", 1)),
        allow_bottlenecks=bool(params.get("allow_bottlenecks", True)),
    )


def run_raster_analysis(
    layer: QgsRasterLayer,
    output_dir: str,
    raw_params: Dict,
    strategy: str = "most_connectivity",
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
    ds = gdal.Open(src_path)
    if ds is None:
        raise RasterAnalysisError(f"Cannot open raster source: {src_path}")

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

    warnings: List[str] = []
    if total_pixels >= PIXEL_COUNT_WARNING_THRESHOLD:
        warnings.append(
            "Large raster detected (>{:,} pixels).".format(
                PIXEL_COUNT_WARNING_THRESHOLD
            )
        )
    if 0 < pixel_size < PIXEL_SIZE_WARNING_THRESHOLD:
        warnings.append(
            f"High-resolution data detected (≈{pixel_size:.2f} map units per pixel). "
            "Consider resampling to a coarser resolution for faster corridor modelling."
        )

    if warnings:
        warning_text = " ".join(warnings)
        if iface and hasattr(iface, "messageBar"):
            try:
                iface.messageBar().pushWarning("Linkscape", warning_text)
            except Exception:
                print(f"WARNING: {warning_text}")
        else:
            print(f"WARNING: {warning_text}")
        raise RasterAnalysisError(
            f"Raster is too large/fine for Linkscape to process efficiently. {warning_text} "
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
    patch_pixels = int(np.sum(patch_mask))
    if patch_pixels == 0:
        raise RasterAnalysisError("No patch pixels found with the current configuration.")
    print(f"  ✓ Patch pixels: {patch_pixels:,}")
    _emit_progress(progress_cb, 25, "Processing habitat patches…")

    obstacle_mask = define_obstacles(data, nodata_mask, patch_mask, params)
    if params.obstacle_enabled:
        obstacle_pixels = int(np.sum(obstacle_mask))
        if obstacle_pixels:
            print(f"  ✓ Obstacle pixels: {obstacle_pixels:,}")
        else:
            print("  ⚠ Obstacle configuration matched no pixels; proceeding without constraints.")
            obstacle_mask = patch_mask.astype(bool)
    else:
        print("  ✓ Obstacles disabled.")
        obstacle_mask = patch_mask.astype(bool)

    passable_mask = _build_passable_mask(
        patch_mask,
        obstacle_mask,
        params.min_corridor_width,
        params.allow_bottlenecks,
    )
    _emit_progress(progress_cb, 35, "Labeling patches…")

    print("\n3. Labeling patches...")
    t0 = time.time()
    labels, n_patches = label_patches(patch_mask, params.patch_connectivity)
    print(f"  ✓ Patches: {n_patches:,} in {time.time() - t0:.2f}s")

    if params.min_patch_size > 0:
        unique_labels, counts = np.unique(labels[labels > 0], return_counts=True)
        valid_labels = unique_labels[counts >= params.min_patch_size]
        new_labels = np.zeros_like(labels)
        for new_id, old_id in enumerate(valid_labels, 1):
            new_labels[labels == old_id] = new_id
        labels = new_labels
        patch_mask = (labels > 0).astype(np.uint8)
        print(f"  ✓ After filter: {len(valid_labels):,} patches")

    debug_label_env = os.environ.get("LINKSCAPE_SAVE_PATCH_LABELS")
    if debug_label_env:
        try:
            label_out_dir = output_dir or os.path.dirname(src_path)
            os.makedirs(label_out_dir, exist_ok=True)
            label_path = os.path.join(label_out_dir, "linkscape_patch_labels.tif")
            write_raster(label_path, labels, gt, proj, nodata=0)
            print(f"  ✓ Saved patch ID raster: {label_path}")
        except Exception as label_exc:  # noqa: BLE001
            print(f"  ⚠ Could not save patch ID raster: {label_exc}")

    unique_labels, counts = np.unique(labels[labels > 0], return_counts=True)
    patch_sizes = dict(zip(unique_labels.tolist(), counts.tolist()))
    if not patch_sizes:
        raise RasterAnalysisError("No valid patches remain after filtering.")
    _emit_progress(progress_cb, 45, "Searching for corridors…")

    print("\n4. Finding possible corridors...")
    candidates = find_all_possible_corridors(
        labels,
        patch_mask,
        patch_sizes,
        params.max_search_distance,
        params.max_corridor_area,
        params.patch_connectivity,
        obstacle_mask=obstacle_mask,
        passable_mask=passable_mask,
        progress_cb=progress_cb,
        progress_start=45,
        progress_end=75,
    )

    if not candidates:
        raise RasterAnalysisError("No feasible corridors found with the current configuration.")
    _emit_progress(progress_cb, 78, "Optimizing corridor selection…")

    strategy = (strategy or "most_connectivity").lower()
    strategy_map = {
        "most_connectivity": (
            optimize_most_connectivity,
            "linkscape_most_connectivity.tif",
            "Corridors (Most Connectivity)",
        ),
        "largest_patch": (
            optimize_largest_patch,
            "linkscape_largest_patch.tif",
            "Corridors (Largest Patch)",
        ),
    }

    if strategy not in strategy_map:
        raise RasterAnalysisError(f"Unsupported strategy '{strategy}'.")

    optimize_func, default_filename, layer_name = strategy_map[strategy]

    print("\n5. Running optimization...")
    print("=" * 70)
    print(f"--- {strategy.replace('_', ' ').upper()} ---")

    corridors, stats = optimize_func(candidates, patch_sizes, params.budget_pixels)
    if not corridors:
        raise RasterAnalysisError("Selected optimization did not produce any corridors.")
    _emit_progress(progress_cb, 85, "Rendering output raster…")

    print("  Creating output raster...")
    output = create_output_raster(
        labels, corridors, params.min_corridor_width, obstacle_mask=obstacle_mask
    )

    if temporary:
        temp_file = tempfile.NamedTemporaryFile(prefix="linkscape_", suffix=".tif", delete=False)
        out_path = temp_file.name
        temp_file.close()
    else:
        out_dir = output_dir or os.path.dirname(src_path)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, default_filename)

    write_raster(out_path, output, gt, proj, nodata=0)
    print(f"  ✓ Saved: {out_path}")
    _emit_progress(progress_cb, 95, "Finishing up…")

    try:
        result_layer = QgsRasterLayer(out_path, layer_name)
        if result_layer.isValid():
            QgsProject.instance().addMapLayer(result_layer)
            print("  ✓ Added to project")
    except Exception as add_exc:  # noqa: BLE001
        print(f"  ⚠ Could not add layer to project: {add_exc}")

    elapsed = time.time() - overall_start
    stats = dict(stats)
    stats["output_path"] = out_path if not temporary else ""
    stats["layer_name"] = layer_name
    stats["budget_total"] = params.budget_pixels

    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    stats_strategy = stats.get("strategy", strategy)
    print(f"Strategy:          {stats_strategy.replace('_', ' ').title()}")
    print(f"Corridors created: {stats.get('corridors_used', 0)}")
    if stats_strategy == "most_connectivity":
        print(f"Connections:       {stats.get('connections_made', 0)}")
        print(f"Largest patch:     {stats.get('largest_group_size', 0):,} px")
    elif "connections_made" in stats:
        print(f"Connections:       {stats.get('connections_made', 0)}")
    if "seed_id" in stats:
        print(f"Seed patch:        {stats.get('seed_id')}")
    print(f"Final size:        {stats.get('final_patch_size', 0):,} px")
    print(f"Budget used:       {stats.get('budget_used', 0)}/{params.budget_pixels} px")
    print(f"Processing time:   {elapsed:.1f}s")
    if temporary:
        print("Output:            Temporary raster layer")
    else:
        print(f"Output GeoTIFF:    {out_path}")
    print("=" * 70)

    ds = None
    _emit_progress(progress_cb, 100, "Raster analysis complete.")
    return [{"strategy": strategy, "stats": stats, "output_path": out_path if not temporary else ""}]
