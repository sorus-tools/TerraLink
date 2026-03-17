# TerraLink QGIS Plugin 1.7

**Ecological Corridor Optimization for Habitat Connectivity**

TerraLink is a QGIS plugin that identifies and optimizes ecological corridors to connect fragmented habitat patches. It uses graph-based and circuit-inspired optimization to maximize landscape connectivity under realistic spatial, budget, and land-use constraints.

![Example corridors](example-optimized-corridor-locations.png)  
Image: Corridors generated (pink) to strategically connect the maximum landscape area (green) while avoiding impassable land type (orange).  

---

## What’s New in 1.7

- **Cleaner, more predictable UI + validation**
  - Controls are organized more clearly.
  - Conflicting/invalid settings are caught earlier.
  - Clearer feedback reduces “mystery” outcomes from unnoticed options.
- **Vector mode is more robust (reference path)**
  - Terminal selection, spacing, and strategy behavior are more internally consistent.
  - Better handling of small N, weird shapes, and edge effects.
- **Raster mode now uses the vector corridor engine**
  - Raster habitat and impassable masks are polygonized, then analyzed with the vector corridor workflow.
  - Raster and vector runs now share the same selection logic, budget handling, and strategy behavior.

## Table of Contents

1. [User Guide](#user-guide)  
2. [Installation and Access](#installation-and-access)  
3. [Quick Start Guide](#quick-start-guide)  
4. [Optimization Modes](#optimization-modes)  
5. [Landscape Metrics](#landscape-metrics)  
6. [Technical Reference](#technical-reference)  
7. [Raster Workflow](#raster-workflow)  
8. [Vector Workflow](#vector-workflow)  
9. [Outputs](#outputs)  
10. [Tips and Best Practices](#tips-and-best-practices)  
11. [Troubleshooting](#troubleshooting)  
12. [Works Cited](#works-cited)  
13. [Credits](#credits)  

---

## User Guide

TerraLink is built for scenario-based corridor planning. The basic workflow is simple:

1. Load habitat data.  
2. Set approximate planning parameters that match your restoration or connectivity goal.  
3. Choose the optimization mode that best matches the ecological outcome you want.  
4. Run several scenarios.  
5. Compare the landscape metrics and keep the configuration that best fits your objective.

### What TerraLink Does

TerraLink identifies candidate corridors between habitat patches and selects the set of corridors that best meets a chosen objective under your constraints.

In practice, you tell TerraLink:

- where the habitat is
- what areas are impassable or undesirable
- how wide corridors should be
- how far TerraLink is allowed to search
- how much corridor area or cost you are willing to spend
- what kind of connectivity you want to maximize

TerraLink then evaluates corridor options and returns outputs you can compare across scenarios.

### What You Can Use It For

- Compare alternative restoration budgets.
- Test whether a barrier layer changes the best corridor layout.
- See whether your goal is better served by one large backbone or several smaller connected networks.
- Compare structural improvement against more movement-oriented or species-oriented outcomes.

### Main Idea

TerraLink is most useful when you do not expect a single “perfect” answer from one run. It works best when you run a small set of plausible scenarios and use the output metrics to decide which tradeoff is best for your landscape and planning goal.

### Main Features

- Raster workflow for land-cover or habitat grids.
- Vector workflow for polygon patch layers.
- Multiple optimization modes for different connectivity goals.
- Optional impassable areas in both raster and vector runs.
- Budget-constrained corridor selection.
- Landscape metrics for comparing scenarios.
- Automatic measurement handling for vector runs.

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

1. **Load your data**  
   Raster: use a land-cover or habitat raster.  
   Vector: use a polygon layer where each feature is a habitat patch.

2. **Open TerraLink**  
   Launch **Run TerraLink** from the toolbar, menu, or Processing Toolbox.

3. **Choose the workflow**  
   Select the input layer.  
   Choose `Raster` or `Vector` if TerraLink does not infer the layer type correctly.

4. **Choose an optimization mode**  
   Pick the mode that matches your goal:
   one backbone, broad connected habitat, smoother movement, or species-oriented reachable habitat.

5. **Set approximate planning parameters**  
   Start with values that are reasonable for your project, not perfect:
   patch size threshold, corridor width, budget, search distance, and optional impassables.

6. **Run the scenario**  
   Click **Run** and watch the **Log** tab for progress and warnings.

7. **Review the outputs**  
   TerraLink adds the corridor outputs and summary tables to QGIS.  
   It also writes a plain-text PRE/POST landscape metrics report.

8. **Run a few alternatives**  
   Change one or two key settings at a time:
   optimization mode, budget, corridor width, impassables, or search distance.

9. **Compare the metrics**  
   Use the PRE/POST landscape metrics and the mapped outputs together.  
   Keep the scenario that best matches your connectivity and restoration objective.

---

## Optimization Modes

TerraLink supports four optimization modes. They are not just different names for the same thing. Each mode rewards a different kind of connectivity improvement, so the right choice depends on what you want the landscape to do (Kupfer 2012; Gustafson 1998; Wiens 1989).

For clarity and reproducibility, TerraLink stores the selected mode as a stable internal key:

- **Largest Single Network** (`largest_single_network`)  
- **Most Connected Area** (`most_connected_habitat`)  
- **Landscape Fluidity (LF‑A)** (`landscape_fluidity`)  
- **Reachable Habitat (Advanced)** (`reachable_habitat_advanced`)  

### Largest Single Network (`largest_single_network`)

- **Plain-language summary:** build one main connected backbone.  
- **What it prioritizes:** growing one dominant connected component under budget, typically by expanding an already large patch or connected group.  
- **Use it when:** you want one strong, continuous system rather than several medium-sized connected clusters.  
- **Keep in mind:** this mode is intentionally biased toward concentration and contiguity, not toward maximizing connected habitat everywhere.

### Most Connected Area (`most_connected_habitat`)

- **Plain-language summary:** connect as much habitat as possible across the whole map.  
- **What it prioritizes:** maximizing the total amount of habitat captured by connected multi-patch networks, even if that creates several separate networks instead of one backbone.  
- **Use it when:** you want broad structural improvement across the landscape and do not need one single dominant network.  
- **Keep in mind:** this is a whole-landscape structural mode, not a movement-quality or species-specific mode.

### Landscape Fluidity (LF‑A) (`landscape_fluidity`)

- **Plain-language summary:** make movement through the network easier and less bottlenecked.  
- **What it prioritizes:** whole-network ease of movement and redundancy, not just whether patches touch. It favors corridors that relieve bottlenecks, shorten routes, and improve movement through the matrix (Kupfer 2012; Wagner & Fortin 2005).  
- **Use it when:** you care about movement quality, alternative routes, and reducing detours rather than only increasing connected area.  
- **Keep in mind:** this mode can produce multiple distinct corridors for the same patch pair and also produce corridors that connect a patch to itself if there are significant connectivity gains by doing so. For example, a meandering U-shaped patch may have a corridor generated to close it to an O-shaped patch. 

### Reachable Habitat (Advanced) (`reachable_habitat_advanced`)

- **Plain-language summary:** maximize how much habitat a focal species can realistically reach.  
- **What it prioritizes:** functional accessibility under a dispersal rule or kernel, not just structural connectedness.  
- **Use it when:** you have a focal species or dispersal assumption and want the mode that is most directly tied to that ecological process.  
- **Key inputs:** species dispersal distance, dispersal kernel, optional minimum patch area for the species, optional patch quality field (vector), and patch-area scaling.

### Simple rule of thumb

- **Most Connected Area**: maximize connected habitat across the map (structural integration).  
- **Largest Single Network**: build one main backbone (consolidation).  
- **Landscape Fluidity**: make movement easier and more redundant across the network (movement quality).  
- **Reachable Habitat (Advanced)**: maximize functionally reachable habitat for a species with a defined dispersal ability (functional accessibility).  

This “different metrics capture different processes” framing is consistent with broader guidance to choose connectivity measures according to the ecological process of interest and to interpret them at an appropriate scale (Kupfer 2012; Gustafson 1998; Wiens 1989).

---

## Landscape Metrics

TerraLink writes a plain-text **Landscape Metrics (PRE/POST)** report for each run.

The main use of these metrics is **comparison**:
compare one scenario against another on the same landscape, using the same input data, mapping scheme, and scale. They are much more useful for scenario selection than for claiming one absolute, context-free “connectivity score” (Gustafson 1998; Wiens 1989; Neel et al. 2004).

### How to interpret the main metrics

#### Connected Habitat Area (structural)
How much habitat ends up inside connected multi-patch networks after corridor selection. This is intuitive and useful, but it is mainly a *structural* metric. It tells you how much habitat has been consolidated into connected systems; it does not, by itself, tell you how easy movement is for an organism (Kupfer 2012; Gustafson 1998).

#### Largest Single Network (structural)
How much habitat is captured by the single biggest connected component. This is useful when your goal is one dominant backbone, but it can be high even when substantial habitat remains disconnected elsewhere. It is therefore a good metric for consolidation, not a complete summary of whole-landscape connectivity (Kupfer 2012).

#### Reachable Habitat Score (functional)
A functional accessibility metric. It combines habitat area and connectivity under a dispersal rule or kernel, so it reflects how much habitat is realistically reachable for an organism rather than only whether habitat patches are structurally connected (Kupfer 2012).

#### Mean Effective Resistance (whole‑network movement difficulty)
A whole-network indicator of movement difficulty or isolation. Lower values generally mean movement is easier across the network. Unlike a single least-cost path, resistance-style summaries can reflect multiple pathways and redundancy, which is why they pair well with fluidity-oriented planning. They are also more parameterized and usually harder to interpret than simple structural patch metrics (Kupfer 2012; Wagner & Fortin 2005).

### Practical interpretation notes

- **Use metrics comparatively:** compare *PRE vs POST* and compare alternative parameterizations on the same landscape; avoid over‑interpreting single runs in isolation (Gustafson 1998; Wiens 1989).  
- **Mind scale and aggregation:** metric behavior can change with mapping resolution, class aggregation, and extent, which can change what a metric “means” (Wiens 1989; Neel et al. 2004).  
- **Structural ≠ functional:** structural connectivity can be a poor proxy for movement unless the metric and parameters align with the organism/process of interest (Kupfer 2012; Gustafson 1998).  

---

## Technical Reference

This section documents implementation details and advanced behavior. If you are new to the tool, stop at the sections above and come back here only when you need to understand why certain runs differ or how the modes behave internally.

### Raster vs Vector differences (important)

Raster mode now polygonizes the selected habitat and impassable masks, then delegates corridor generation to the vector engine. In practice, raster and vector runs still differ when the starting habitat representation differs, but corridor selection and optimization now follow the same backend logic.

### Landscape Fluidity advanced behavior

- **Multiple candidates per patch pair:** LF‑A can keep more than one distinct corridor candidate between the same pair of patches, so it can reward redundancy and detour relief instead of always taking only the single cheapest option.  
- **Stepping-stone logic:** if two patches are farther apart than the **maximum search distance**, TerraLink cannot connect them directly. LF‑A can still improve movement by favoring intermediate stepping-stone links when they improve the whole network.  
- **Raster delegation detail:** raster runs first polygonize habitat and impassable masks, then pass those temporary layers into the vector corridor engine.  

## Raster Workflow

Raster mode is designed for land cover or habitat grids, where each cell can be classified as habitat, impassable, or background.

### Patch definition

- **Patch connectivity**: choose 4- or 8-neighbor connectivity when identifying contiguous patches.  
- **Patch values**: specify one or more values that represent habitat (for example 1 for forest, 2 for wetland). TerraLink requires at least one value.  
- **Minimum patch size**: patches smaller than this threshold (in pixels) are removed from the main patch set before corridor planning.  
- **Raster units**: switch between Pixels, Metric (ha/m), and Imperial (ac/ft). Metric/Imperial requires a projected CRS (meters or feet).  

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

## Outputs

For each run, TerraLink produces mapped outputs plus summary information you can use for comparison.

### Raster mode outputs

- A corridor raster.
- A contiguous-areas raster showing the connected component size for habitat and corridor cells.
- A QGIS summary table layer: `TerraLink Raster Summary (<input layer>)`.
- A plain-text PRE/POST landscape metrics report, usually named `landscape_metrics_<input layer>.txt`.
- Output filenames that include the input layer name.

If you use `Assign corridor cells`, the corridor raster can store one of three values:

- `Sum area of patches directly connected`: each corridor cell stores `area(patch1) + area(patch2)` in pixels.
- `Sum area of total network`: each corridor cell stores the total connected network area in pixels.
- `Efficiency (corridor area / connected patches)`: each corridor cell stores `corridor_area / (area(patch1) + area(patch2))`.

### Vector mode outputs

- A corridor GeoPackage layer with patch ids, corridor area, connected area, efficiency, and multipart/segment flags.
- A contiguous-areas layer in the same GeoPackage.
- A plain-text PRE/POST landscape metrics report, usually named `landscape_metrics_<input layer>.txt`.
- A summary CSV, also added to QGIS as a table layer, usually named `terralink_vector_summary_<input layer>.csv`.
- GeoPackage filenames that include the input layer name unless you already include it yourself.

### Run summaries

- The **Log** tab provides a concise run summary.
- Raster and vector runs add in-project summary tables so you can compare scenarios inside QGIS.

---

## Tips and Best Practices

- Start with a small area and conservative budget to verify parameter settings before scaling up.  
- Keep patch definitions simple at first. For raster mode, start with a single habitat value and add ranges later.  
- Use bottleneck controls to adjust realism rather than trying to encode everything into the land cover itself.  
- Use the Log tab summary (and the raster summary table layer) to confirm budgets and connected area before iterating on parameters.

---


## Troubleshooting

- **1) “No feasible corridors / no corridors produced”**  
  - Increase **budget** and/or **maximum search distance**.  
  - Decrease **minimum corridor width** (a wide corridor can block all routes when impassables/bottlenecks are enabled).  
  - Ensure you have at least **two** patches after filtering (raster: min patch size; vector: min patch size + valid polygons).

- **2) Reachable Habitat mode returns no corridors**  
  - Set a **species dispersal distance** greater than 0; the mode is disabled when dispersal distance is unset/zero.  
  - If it still selects nothing, increase budget and/or dispersal distance, or reduce corridor width.

- **3) Corridors “ignore” impassables / overlap barriers**  
  - Confirm impassables are correctly configured (raster values/ranges, or polygon layers for vector).  
  - Vector: use a finer **navigator grid resolution** if corridors are slipping through thin barriers.  
  - Raster: tighten bottleneck behavior (or increase corridor width) so corridors cannot squeeze through unrealistic gaps.

- **4) Corridors disappear or look wrong after the run**  
  - Check that your patches are valid and not self‑intersecting (vector). Fix geometries if needed.  
  - For very small patch gaps, small geometry/topology differences between pipelines can matter; compare scenarios within the same pipeline.

- **5) Vector mode says there are “too many patches”**  
  - Clip to a smaller study area, raise minimum patch size, or use raster mode for large patch counts.

- **6) Runs are slow (high CPU for a long time)**  
  - Reduce study area (clip), coarsen raster resolution, or increase vector minimum patch size.  
  - Disable impassables temporarily to confirm baseline performance.  
  - Try a smaller max search distance and a smaller budget to get a “first pass” quickly, then refine.

- **7) QGIS looks frozen / UI stops responding**  
  - TerraLink runs can be compute-heavy; give it time and watch the **Log** tab.  
  - If it’s stuck for minutes with no log progress: cancel/close QGIS, then rerun with smaller extent/coarser resolution and fewer constraints.

- **8) “Impassable configuration matched no pixels” (raster) / impassable layers skipped (vector)**  
  - Raster: verify the impassable values/ranges actually exist in the raster (check the layer’s value distribution).  
  - Vector: ensure selected impassable layers are **polygon** geometry and are still present/valid in the project.

- **9) Metric/imperial units behave strangely**  
  - Make sure your layer CRS and units are appropriate. Vector runs are reprojected internally for measurement, but inputs with bad CRS/extents can cause odd distances.

- **10) Landscape metrics report is missing or shows an error**  
  - Some environments may lack optional dependencies used for metrics. Corridors can still be valid even if a metrics report fails.  
  - Use the in‑project summary table layers and rerun on a smaller area if the metrics step is failing due to resource limits.

---

## Works Cited

- Gustafson, E. J. (1998). Quantifying landscape spatial pattern: What is the state of the art? *Ecosystems, 1*, 143–156.  
- Kupfer, J. A. (2012). Landscape ecology and biogeography: Rethinking landscape metrics in a post‑FRAGSTATS landscape. *Progress in Physical Geography, 36*(3), 400–420.  
- Neel, M. C., McGarigal, K., & Cushman, S. A. (2004). Behavior of class‑level landscape metrics across gradients of class aggregation and area. *Landscape Ecology, 19*, 435–455.  
- Wiens, J. A. (1989). Spatial scaling in ecology. *Functional Ecology, 3*, 385–397.  
- Wagner, H. H., & Fortin, M.‑J. (2005). Spatial analysis of landscapes: Concepts and statistics. *Ecology, 86*(8), 1975–1987.  

---

## Credits

TerraLink was created by Ben Bishop (SORUS Tools) as a practical tool for habitat connectivity and corridor planning. Previous to December 15th 2025, it was released under a different name (Linkscape). Deprecated versions are still available on the QGIS repository. 

TerraLink is built on the QGIS Python API and commonly used geospatial libraries (GDAL/OGR, NumPy).

Bug reports and feature requests are welcome on the GitHub issue tracker listed in the plugin metadata or sent to SorusConsulting@gmail.com
