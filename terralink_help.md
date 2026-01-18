# TerraLink at a Glance

- Pick your input layer of fragmented patches. This could be a raster LULC dataset or a vector set of polygon patches. 

- Set an output folder unless you intentionally want a temporary layer.

- If using raster data, tell TerraLink which value(s) you would like to define as patches.

- Set key criteria such as minimum patch size, optimization mode and your budget. 

- Optimization modes: **Largest Single Network** and **Most Connectivity** (utility/ROI-based selection with shortcut efficiency for redundancy).

- Enable **Impassable Land Classes** when there are areas that corridors must avoid (roads, water, cities). Set values for this in raster, or select impassable polygon layers in vector mode.

- Keep the search distance reasonable—oversized values slow things down.

- If performance is slow, resample very fine rasters or trim your AOI before running.

- The output data comes embedded with helpful info. Raster cells store the connected component size, and vector corridor attributes include patch ids, corridor area, connected area, and efficiency. 
- Vector runs also create a contiguous areas layer (patches + corridors dissolved by network).
- Raster and vector runs add summary table layers in QGIS, and vector runs also write a summary CSV.
- Landscape metrics reports are written alongside outputs as `landscape_metrics_<input layer>.txt`.

Need more help? Visit the README in this plugin folder for a deeper walkthrough.
