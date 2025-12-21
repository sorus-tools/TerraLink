# TerraLink QGIS Plugin v1.1.0

**Ecological Corridor Optimization for Habitat Connectivity**

TerraLink is a QGIS plugin that identifies and optimizes ecological corridors to connect fragmented habitat patches. It uses graph-based and circuit-inspired optimization to maximize landscape connectivity under realistic spatial, budget, and land-use constraints.

![Example corridors](example-optimized-corridor-locations.png)  
Image: Corridors generated (pink) to strategically connect the maximum landscape area (green) while avoiding impassable land type (orange).  

---

## What’s New in v1.1.0

- **Explicit optimization strategies**
  - Two modes: **Largest Network** and **Circuit Theory**
  - Circuit Theory prioritizes corridors with the highest marginal connectivity gain and allows ecologically meaningful redundancy.

- **Improved corridor selection**
  - Avoids redundant or overlapping corridors unless they add real connectivity value.
  - Remaining budget is preferentially used to reinforce existing corridors rather than create low-value links.

- **Vector mode: boundary-aware connectivity**
  - Patch connections use true boundary distances instead of centroids.
  - Multiple boundary terminals are sampled so narrow peninsulas or patch tips can connect realistically.

- **Impassible Landclass and Barrier routing**
  - When impassable polygon layers are enabled, TerraLink routes corridors around barriers using an internal navigator grid.
  - Optional smoothing reduces stair-step artifacts from grid-based routing.

- **Landscape metrics report**
  - Each run writes a plain-text landscape metrics summary describing patch count, connected area, fragmentation, and network structure.

## Table of Contents

1. [Overview](#overview)  
2. [Installation and Access](#installation-and-access)  
3. [Quick Start Guide](#quick-start-guide)  
4. [Raster Workflow](#raster-workflow)  
5. [Vector Workflow](#vector-workflow)  
6. [Optimization Modes](#optimization-modes)  
7. [Outputs](#outputs)  
8. [Tips and Best Practices](#tips-and-best-practices)  
9. [Troubleshooting](#troubleshooting)  
10. [Credits](#credits)  

---

## Overview

### What TerraLink Does

TerraLink analyzes habitat patches (either raster cells or vector polygons) and identifies corridors that connect them based on:  

- Available budget in pixels or area units  
- Spatial constraints like maximum search distance and minimum corridor width  
- Optimization mode (Largest Network vs Circuit Theory)  

The plugin supports both raster and vector workflows so you can work directly with land cover grids or polygon habitat maps.

### Key Features

- **Dual workflows**  
  - Raster mode for land cover or habitat grids  
  - Vector mode for polygon patch layers (one feature per patch)  

- **Flexible patch definition (raster)**  
  - Define habitat patches by one or more pixel values or value ranges  
  - Filter out tiny patches below a minimum size before corridor planning  

- **Optimization modes**  
  - Largest Network  
  - Circuit Theory  

- **Impassable Land Classes**  
  - Raster: treat one or more values or ranges as impassable land classes and control whether corridors can squeeze through narrow gaps.  
  - Vector: optionally build a raster grid over one or more impassable polygon layers for pathfinding around barriers (roads, rivers, urban footprints).  

- **Budget-constrained analysis**  
  - Limit total corridor construction cost in pixels (raster) or hectares/acres (vector). Terralink will choose the optimal configuration of the corridors. 

- **Unit flexibility (vector)**  
  - Choose metric or imperial units for corridor width, distances, and budgets.  

- **Automatic CRS handling (vector)**  
  - Geographic layers are reprojected to an appropriate UTM CRS internally so distances and areas are measured in consistent units.  

---

## Installation and Access

### Prerequisites

- QGIS 3.22 or newer.  
- Standard QGIS Python stack: `numpy`, GDAL/OSGEO, and PyQt5.  

### Installing from a ZIP

1. Download the TerraLink plugin zip from the releases page.  
2. In QGIS, go to `Plugins` → `Manage and Install Plugins`.  
3. Choose `Install from ZIP`, browse to the TerraLink zip, and install.  
4. Enable TerraLink in the plugin manager if it is not already enabled.  

Once installed, you will see **Run TerraLink** in the toolbar and in the `TerraLink` menu.  
You can also launch it from the Processing Toolbox under `TerraLink → TerraLink Corridor Generation`.

---

## Quick Start Guide

1. **Prepare your data**  
   - Raster mode: add a land cover or habitat raster with integer values that distinguish habitat from non-habitat.  
   - Vector mode: add a polygon layer where each feature is a separate habitat patch.

2. **Launch TerraLink**  
   - Click **Run TerraLink** on the toolbar or from the `TerraLink` menu.  

3. **Choose mode and layer**  
   - Select the input layer in the combo box.  
   - Pick `Raster` or `Vector` in the layer type selector if the plugin does not infer it correctly.  

4. **Pick an optimization mode**  
   - Use Largest Network to prioritize one dominant connected network.  
   - Use Circuit Theory to prioritize highest ROI corridors and meaningful redundancies.  

5. **Set basic parameters**  
   - Minimum patch size (filters tiny patches out of the habitat pool).  
   - Minimum corridor width.  
   - Budget and maximum search distance.  

6. **Set impassable land classes (optional)**  
   - Raster: choose impassable values or ranges, and decide if corridors may pass through small gaps (bottlenecks).  
   - Vector: enable impassable land classes, select one or more impassable layers, and set grid resolution for the navigator.  

7. **Run and review**  
   - Click **Run**. Watch progress in the **Log** tab.  
   - Outputs are added to your QGIS project; raster runs also add a “TerraLink Raster Summary …” table layer.  
   - After a successful run, the Run button becomes **Close** to avoid accidental re-runs.  

---

## Raster Workflow

Raster mode is designed for land cover or habitat grids, where each cell can be classified as habitat, impassable, or background.

### Patch definition

- **Patch connectivity**: choose 4- or 8-neighbor connectivity when identifying contiguous patches.  
- **Patch values**: specify one or more values that represent habitat (for example 1 for forest, 2 for wetland). TerraLink requires at least one value.  
- **Minimum patch size**: patches smaller than this threshold (in pixels) are removed from the main patch set before corridor planning.  

### Impassable land classes and bottlenecks

- Enable impassable land classes and define one or more impassable values or ranges.  
- TerraLink builds an impassable mask and then a passable mask that:  
  - Treats impassable cells as fully blocked.  
  - Treats habitat cells as non-impassable so corridors can touch patch edges.  
  - Optionally erodes the available space with the corridor width kernel so corridors cannot use gaps narrower than the configured width, unless you allow bottlenecks.

### Corridor budget and search distance

- **Budget (pixels)**: total corridor pixels you are willing to “spend.”  
- **Maximum search distance**: hard cap on how far corridors can reach when searching for candidate links.

## Vector Workflow

Vector mode is designed for polygon patch datasets, where each feature is a separate patch.

### Input requirements

- One feature per patch. If the input layer has only one feature, TerraLink will warn and refuse to run, since a single patch cannot be connected to anything.  

### Units and basic parameters

- Select metric or imperial units in the dialog.  
- Set:
  - Minimum corridor width (meters or feet)  
  - Minimum patch size (hectares or acres)  
  - Budget area (same area units)  
  - Maximum search distance (meters or feet)  

### Impassable land classes and navigator grid

- Enable impassable land classes to route around non-habitat features such as roads, rivers, or urban areas.  
- Select one or more vector impassable layers and a grid resolution for the navigator. The engine builds a raster grid over the combined impassable geometry and uses it for pathfinding.  

## Optimization Modes

TerraLink supports two optimization modes:

### Largest Network

- Goal: prioritize building one dominant connected network under the budget.

### Circuit Theory

- Goal: maximize system-wide utility (ROI) under the budget; can favor robust redundant links between major hubs.

---

## Outputs

For each run, TerraLink can produce:

### Raster mode outputs

- Corridor raster where pixel values indicate the size (in pixels) of the connected component created by including that corridor cell.  
- Contiguous areas raster where pixel values indicate the size (in pixels) of the connected component for all habitat and corridor cells.  
- A small in-project summary table layer: `TerraLink Raster Summary (<input layer>)`.
- A plain-text landscape metrics report saved alongside the other outputs (or as a temporary file when using temporary outputs).

### Vector mode outputs

- Corridor layer as a GeoPackage, with each corridor feature containing:  
  - Patch ids it connects  
  - Length and cost  
  - Connected area and strategy statistics  
- A plain-text landscape metrics report saved alongside the other outputs (or as a temporary file when using temporary outputs).

### Run summaries

- The **Log** tab provides a concise run summary suitable for copy/paste into notes.
- Raster runs also add an in-project summary table layer.

---

## Tips and Best Practices

- Start with a small area and conservative budget to verify parameter settings before scaling up.  
- Keep patch definitions simple at first. For raster mode, start with a single habitat value and add ranges later.  
- Use bottleneck controls to adjust realism rather than trying to encode everything into the land cover itself.  
- Use the Log tab summary (and the raster summary table layer) to confirm budgets and connected area before iterating on parameters.

---


## Troubleshooting

- **No corridors produced**  
  - Check that budget and max search distance are large enough.  
  - For vector mode, ensure there is more than one patch feature and that impassable land classes are not blocking all routes.

- **Corridors clip into impassable land classes**  
  - Increase minimum corridor width or tighten bottleneck rules in raster mode.  
  - For vector mode, verify that impassable layers are valid and that the navigator grid resolution is appropriate.

- **Runs are slow**  
  - Reduce extent or resample rasters to a coarser resolution.  
  - In vector mode, filter out very small patches with the minimum patch size setting.

---

## Credits

TerraLink was created by Ben Bishop (SORUS) as a practical tool for habitat connectivity and corridor planning. Previous to December 15th 2025, it was refered to as Linkscape. Deprecated versions of Linkscape are still available on the QGIS repository. 

TerraLink is built on the QGIS Python API and commonly used geospatial libraries (GDAL/OGR, NumPy).

Bug reports and feature requests are welcome on the GitHub issue tracker listed in the plugin metadata or sent to SorusConsulting@gmail.com
