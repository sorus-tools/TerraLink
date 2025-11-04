"""
Linkscape Corridor Analysis - Raster Workflow (v23.0)
-----------------------------------------------------
Runs the selected raster optimization workflow for corridor analysis.

The logic is adapted from the standalone raster script and packaged so it can
be invoked from the QGIS plugin with user-supplied parameters.
"""

from __future__ import annotations

import os
import struct
import tempfile
import time
from collections import deque
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from osgeo import gdal
from qgis.core import QgsProject, QgsRasterLayer

# Optional OpenCV import for faster connected component labeling
try:
    import cv2

    HAS_CV2 = True
except ImportError:
    HAS_CV2 = False

GTIFF_OPTIONS = ["COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"]


class RasterAnalysisError(RuntimeError):
    """Raised when the raster analysis cannot be completed."""


@dataclass
class RasterRunParams:
    patch_connectivity: int
    patch_mode: str
    patch_values: List[float]
    range_lower: Optional[float]
    range_upper: Optional[float]
    value_tolerance: float
    nodata_fallback: float
    min_patch_size: int
    budget_pixels: int
    max_search_distance: int
    max_corridor_area: Optional[int]
    min_corridor_width: int  # Included for completeness even if not used directly.


class UnionFind:
    """Union-Find data structure for tracking connected components."""

    def __init__(self):
        self.parent: Dict[int, int] = {}
        self.size: Dict[int, int] = {}

    def find(self, x: int) -> int:
        if x not in self.parent:
            self.parent[x] = x
            self.size[x] = 0
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
        return ra

    def get_size(self, x: int) -> int:
        return self.size.get(self.find(x), 0)


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


def read_band(band: gdal.Band, rows: int, cols: int) -> np.ndarray:
    """Read a raster band as a numpy array."""
    buf = band.ReadRaster(0, 0, cols, rows, cols, rows, gdal.GDT_Float32)
    vals = struct.unpack("f" * (rows * cols), buf)
    return np.array(vals, dtype=np.float32).reshape(rows, cols)


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


def find_shortest_corridor(
    start_patches: Set[int],
    labels: np.ndarray,
    habitat: np.ndarray,
    max_width: int,
    connectivity: int,
) -> Optional[Tuple[frozenset, Set[int], int]]:
    """
    BFS to find shortest corridor connecting start_patches to any other patch.
    Returns (path_pixels, target_patches, length) or None
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
            if not habitat[i, j]:
                for di, dj in moves:
                    ni, nj = i + di, j + dj
                    if 0 <= ni < rows and 0 <= nj < cols and labels[ni, nj] in start_patches:
                        start_positions.append((i, j))
                        break

    if not start_positions:
        return None

    queue: deque = deque()
    visited: Set[Tuple[int, int]] = set()
    for r, c in start_positions:
        queue.append((r, c, 1, frozenset({(r, c)})))
        visited.add((r, c))

    while queue:
        r, c, depth, path = queue.popleft()

        target_patches: Set[int] = set()
        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                lbl = labels[nr, nc]
                if lbl > 0 and lbl not in start_patches:
                    target_patches.add(lbl)

        if target_patches:
            return path, target_patches, depth

        if depth < max_width:
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    if not habitat[nr, nc] and (nr, nc) not in visited:
                        visited.add((nr, nc))
                        queue.append((nr, nc, depth + 1, path | frozenset({(nr, nc)})))

    return None


def find_all_possible_corridors(
    labels: np.ndarray,
    habitat: np.ndarray,
    patch_sizes: Dict[int, int],
    max_width: int,
    max_area: Optional[int],
    connectivity: int,
) -> List[Dict]:
    """Find all possible corridors between patch pairs."""
    print("  Finding all possible corridors...")
    all_corridors: List[Dict] = []
    processed_pairs: Set[frozenset] = set()

    unique_patches = [p for p in patch_sizes.keys() if p > 0]
    for idx, patch_id in enumerate(unique_patches):
        if (idx + 1) % 10 == 0:
            print(f"    Analyzing patch {idx + 1}/{len(unique_patches)}...", end="\r")

        result = find_shortest_corridor({patch_id}, labels, habitat, max_width, connectivity)
        if result is None:
            continue

        path_pixels, target_patches, path_len = result

        for target_id in target_patches:
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

    return selected, {
        "strategy": "most_connectivity",
        "corridors_used": len(selected),
        "connections_made": connections_made,
        "budget_used": budget - remaining_budget,
        "total_connected_size": sum(root_sizes.values()),
        "groups_created": len(root_sizes),
        "largest_group_size": max(root_sizes.values()) if root_sizes else 0,
    }


def optimize_largest_patch(
    candidates: List[Dict], patch_sizes: Dict[int, int], budget: int
) -> Tuple[Dict, Dict]:
    """Strategy 2: Largest Patch - create single largest connected patch."""
    print("  Strategy: LARGEST PATCH")
    candidates_copy = [c.copy() for c in candidates]
    candidates_copy.sort(key=lambda x: x["length"])

    sorted_patch_ids = sorted(patch_sizes.keys(), key=lambda k: patch_sizes[k], reverse=True)

    print(f"  Testing {min(20, len(sorted_patch_ids))} seed patches...")
    best_result = {"corridors": {}, "final_size": 0, "seed_id": None, "budget_used": 0}

    for i, seed_id in enumerate(sorted_patch_ids[:20]):
        if (i + 1) % 5 == 0:
            print(f"    Seed {i + 1}/20...", end="\r")

        uf = UnionFind()
        for pid, size in patch_sizes.items():
            uf.parent[pid] = pid
            uf.size[pid] = size

        sim_corridors: Dict[int, Dict] = {}
        remaining_budget = float(budget)

        while True:
            best_link = None
            for cand in candidates_copy:
                p1, p2 = cand["patch1"], cand["patch2"]
                is_p1_in = uf.find(p1) == uf.find(seed_id)
                is_p2_in = uf.find(p2) == uf.find(seed_id)
                if is_p1_in != is_p2_in and cand["length"] <= remaining_budget:
                    best_link = cand
                    break

            if best_link is None:
                break

            p1, p2, length = best_link["patch1"], best_link["patch2"], best_link["length"]
            remaining_budget -= length
            root = uf.union(p1, p2)
            uf.size[root] += length
            sim_corridors[len(sim_corridors)] = best_link

        final_size = uf.get_size(seed_id)
        if final_size > best_result["final_size"]:
            best_result = {
                "corridors": sim_corridors,
                "final_size": final_size,
                "seed_id": seed_id,
                "budget_used": budget - remaining_budget,
            }

    print(f"\n  ✓ Best: Patch {best_result['seed_id']} -> {best_result['final_size']:,} px")

    if not best_result["corridors"]:
        return {}, {"strategy": "largest_patch", "corridors_used": 0}

    selected: Dict[int, Dict] = {}
    for i, corr in best_result["corridors"].items():
        selected[i + 1] = {
            "pixels": corr["pixels"],
            "patch_ids": {corr["patch1"], corr["patch2"]},
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
    }


def create_output_raster(labels: np.ndarray, corridors: Dict[int, Dict]) -> np.ndarray:
    """Create an output raster with corridors marked by connected size."""
    output = np.zeros_like(labels, dtype=np.int32)
    for corridor_data in corridors.values():
        score = corridor_data["connected_size"]
        for r, c in corridor_data["pixels"]:
            output[r, c] = max(output[r, c], score)
    return output


def _to_dataclass(params: Dict) -> RasterRunParams:
    """Convert raw parameter dict into the expected dataclass."""
    return RasterRunParams(
        patch_connectivity=int(params.get("patch_connectivity", 4)),
        patch_mode=str(params.get("patch_mode", "value")).lower(),
        patch_values=list(params.get("patch_values", [])),
        range_lower=params.get("range_lower"),
        range_upper=params.get("range_upper"),
        value_tolerance=float(params.get("value_tolerance", 1e-6)),
        nodata_fallback=float(params.get("nodata_fallback", -9999)),
        min_patch_size=int(params.get("min_patch_size", 0)),
        budget_pixels=int(params.get("budget_pixels", 0)),
        max_search_distance=int(params.get("max_search_distance", 50)),
        max_corridor_area=(
            int(params["max_corridor_area"]) if params.get("max_corridor_area") is not None else None
        ),
        min_corridor_width=int(params.get("min_corridor_width", 1)),
    )


def run_raster_analysis(
    layer: QgsRasterLayer,
    output_dir: str,
    raw_params: Dict,
    strategy: str = "most_connectivity",
    temporary: bool = False,
    iface=None,
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
    print(f"  Size: {rows:,} x {cols:,} = {rows * cols:,} pixels")

    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    if nodata is None:
        nodata = params.nodata_fallback

    print("  Reading data...")
    data = read_band(band, rows, cols)
    nodata_mask = np.abs(data - nodata) < params.value_tolerance if nodata is not None else np.zeros_like(
        data, dtype=bool
    )

    print("\n2. Identifying patch pixels...")
    patch_mask = define_habitat(data, nodata_mask, params)
    patch_pixels = int(np.sum(patch_mask))
    if patch_pixels == 0:
        raise RasterAnalysisError("No patch pixels found with the current configuration.")
    print(f"  ✓ Patch pixels: {patch_pixels:,}")

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

    unique_labels, counts = np.unique(labels[labels > 0], return_counts=True)
    patch_sizes = dict(zip(unique_labels.tolist(), counts.tolist()))
    if not patch_sizes:
        raise RasterAnalysisError("No valid patches remain after filtering.")

    print("\n4. Finding possible corridors...")
    candidates = find_all_possible_corridors(
        labels,
        patch_mask,
        patch_sizes,
        params.max_search_distance,
        params.max_corridor_area,
        params.patch_connectivity,
    )

    if not candidates:
        raise RasterAnalysisError("No feasible corridors found with the current configuration.")

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

    print("  Creating output raster...")
    output = create_output_raster(labels, corridors)

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

    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"Strategy:          {strategy.replace('_', ' ').title()}")
    print(f"Corridors created: {stats.get('corridors_used', 0)}")
    if "connections_made" in stats:
        print(f"Connections:       {stats.get('connections_made', 0)}")
        print(f"Largest patch:     {stats.get('largest_group_size', 0):,} px")
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
    return [{"strategy": strategy, "stats": stats, "output_path": out_path if not temporary else ""}]
