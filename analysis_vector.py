"""
Linkscape Corridor Analysis - Vector Workflow (v23.0)
----------------------------------------------------
Runs the single-strategy corridor optimization workflow for polygon
patch datasets (e.g. shapefiles, GeoPackage layers). Results are
written to individual layers inside a GeoPackage and added to QGIS.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

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
from PyQt5.QtCore import QVariant

BUFFER_SEGMENTS = 16

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


class UnionFind:
    """Union-Find data structure for tracking connected components."""

    def __init__(self):
        self.parent: Dict[int, int] = {}
        self.size: Dict[int, float] = {}

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

    def get_size(self, x: int) -> float:
        return self.size.get(self.find(x), 0.0)


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


def _buffer_line_segment(line_geom: QgsGeometry, width: float) -> QgsGeometry:
    return line_geom.buffer(width / 2.0, BUFFER_SEGMENTS)


def find_all_possible_corridors(
    patches: Dict[int, Dict],
    spatial_index: QgsSpatialIndex,
    params: VectorRunParams,
) -> List[Dict]:
    """
    Find all possible corridors between patch pairs within the distance constraint.
    Returns list of corridor candidates with geometry and metadata.
    """
    print("  Finding all possible corridors...")
    all_corridors: List[Dict] = []
    processed_pairs: Set[frozenset] = set()

    for idx, (pid1, pdata1) in enumerate(patches.items(), start=1):
        if idx % 10 == 0:
            print(f"    Analyzing patch {idx}/{len(patches)}...", end="\r")

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

            corridor_line = QgsGeometry.fromPolylineXY([QgsPointXY(p1), QgsPointXY(p2)])
            corridor_geom = _buffer_line_segment(corridor_line, params.min_corridor_width)
            corridor_geom = corridor_geom.difference(geom1)
            corridor_geom = corridor_geom.difference(geom2)
            if corridor_geom.isEmpty():
                continue
            corridor_geom = corridor_geom.makeValid()
            if corridor_geom.isEmpty():
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

    print(f"\n  ✓ Found {len(all_corridors)} possible corridors")
    return all_corridors


def find_shortest_corridor_from_group(
    start_patches: Set[int],
    all_patches: Dict[int, Dict],
    spatial_index: QgsSpatialIndex,
    params: VectorRunParams,
) -> Optional[Tuple]:
    """Find the smallest-area corridor connecting the current group to any other patch."""
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

        corridor_line = QgsGeometry.fromPolylineXY([QgsPointXY(p1), QgsPointXY(p2)])
        corridor_geom = _buffer_line_segment(corridor_line, params.min_corridor_width)
        corridor_geom = corridor_geom.difference(start_union)
        corridor_geom = corridor_geom.difference(target_geom)
        if corridor_geom.isEmpty():
            continue
        corridor_geom = corridor_geom.makeValid()
        if corridor_geom.isEmpty():
            continue

        corridor_area_ha = corridor_geom.area() / 10000.0
        if corridor_area_ha <= 0:
            continue
        if params.max_corridor_area is not None and corridor_area_ha > params.max_corridor_area:
            continue

        if corridor_area_ha < best_area:
            best_corridor = corridor_geom
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
    print("  Strategy: MOST CONNECTIVITY")
    print("    Calculating efficiency for each corridor...")
    candidates_copy = [c.copy() for c in candidates]
    for c in candidates_copy:
        benefit = patches[c["patch1"]]["area_ha"] + patches[c["patch2"]]["area_ha"]
        c["efficiency"] = benefit / c["area_ha"] if c["area_ha"] > 0 else float("inf")

    print("    Sorting by efficiency (benefit per hectare)...")
    candidates_copy.sort(key=lambda x: x["efficiency"], reverse=True)

    uf = UnionFind()
    for pid, pdata in patches.items():
        uf.find(pid)
        uf.size[pid] = pdata["area_ha"]

    selected: Dict[int, Dict] = {}
    remaining_budget = params.budget_area
    connections_made = 0

    for cand in candidates_copy:
        if cand["area_ha"] > remaining_budget:
            continue
        p1, p2 = cand["patch1"], cand["patch2"]
        if uf.find(p1) == uf.find(p2):
            continue

        root = uf.union(p1, p2)
        uf.size[root] += cand["area_ha"]
        connections_made += 1

        connected_area = uf.get_size(root)
        efficiency = cand["area_ha"] / connected_area if connected_area > 0 else 0

        selected[len(selected) + 1] = {
            "geom": cand["geom"],
            "patch_ids": {p1, p2},
            "area_ha": cand["area_ha"],
            "connected_area_ha": connected_area,
            "efficiency": efficiency,
        }

        remaining_budget -= cand["area_ha"]

    root_areas = {uf.find(pid): uf.get_size(pid) for pid in patches.keys()}
    stats = {
        "strategy": "most_connectivity",
        "corridors_used": len(selected),
        "connections_made": connections_made,
        "budget_used_ha": params.budget_area - remaining_budget,
        "total_connected_area_ha": sum(root_areas.values()),
        "groups_created": len(root_areas),
        "largest_group_area_ha": max(root_areas.values()) if root_areas else 0,
        "patches_connected": len(patches),
    }
    print(f"  ✓ Selected {len(selected)} efficient corridors")
    return selected, stats


def optimize_largest_patch(
    patches: Dict[int, Dict],
    spatial_index: QgsSpatialIndex,
    params: VectorRunParams,
) -> Tuple[Dict[int, Dict], Dict]:
    """Strategy 2: Largest Patch - create a single largest connected patch."""
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
            result = find_shortest_corridor_from_group(current_patches, patches, spatial_index, params)
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
    print("LINKSCAPE VECTOR ANALYSIS v23.0")
    print("=" * 70)
    print("\n1. Loading vector layer...")
    print(f"  ✓ Using layer: {layer.name()} ({layer.featureCount()} features)")

    original_crs = layer.crs()
    print(f"  CRS: {original_crs.authid()}")

    print("\n2. Determining analysis CRS...")
    target_crs = get_utm_crs_from_extent(layer)
    print(f"  ✓ Using {target_crs.authid()} for measurements")

    print("\n3. Loading patches...")
    patches, spatial_index = load_and_prepare_patches(layer, target_crs, params)
    if not patches:
        raise VectorAnalysisError("No patches found after filtering.")

    print("\n4. Precomputing candidate corridors...")
    all_possible = find_all_possible_corridors(patches, spatial_index, params)

    strategy = (strategy or "most_connectivity").lower()
    strategy_map = {
        "most_connectivity": (
            optimize_most_connectivity,
            "Corridors (Most Connectivity)",
        ),
        "largest_patch": (
            optimize_largest_patch,
            "Corridors (Largest Patch)",
        ),
    }

    if strategy not in strategy_map:
        raise VectorAnalysisError(f"Unsupported strategy '{strategy}'.")

    optimize_func, layer_name = strategy_map[strategy]

    print("\n5. Running optimization...")
    print("=" * 70)
    print(f"--- {strategy.replace('_', ' ').upper()} ---")

    if strategy == "most_connectivity" and not all_possible:
        raise VectorAnalysisError("No feasible corridors available for the selected strategy.")

    if strategy == "most_connectivity":
        corridors, stats = optimize_func(patches, all_possible, params)
    else:
        corridors, stats = optimize_func(patches, spatial_index, params)

    if not corridors:
        raise VectorAnalysisError("Selected optimization did not produce any corridors.")

    stats = _convert_stats_for_units(stats, unit_system)

    print("  Preparing outputs...")
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
