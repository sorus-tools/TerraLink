# Linkscape at a Glance

- Pick your input layer of fragmented patches. This could be a raster LULC dataset or a vector set of polygon patches. 

- Set an output folder unless you intentionally want a temporary layer.

- If using raster data, tell Linkscape which value(s) you would like to define as patches.

- Set key criteria such as minimum patch size, optimization mode and your budget. 

- The optimization mode "Most Connectivity" simply connects as much area as possible. The alternative, "Largest Patch" seeks to create a single connected network of patches that is the largest possible, in total area. In other words, the result of "Largest Patch" will always be contiguous. 

- Enable obstacles when there are areas that corridors must avoid (roads, water, cities). Set values for this in raster, or select obstacle layer in vector mode. 

- Keep the search distance reasonable—oversized values slow things down.

- If performance is slow, resample very fine rasters or trim your AOI before running.

- The output data comes embedded with helpful info. Resulting raster cells have the value of the total area that they connect. Vector output has info on the corridors and connected patches in its attribute table. 

Need more help? Visit the README in this plugin folder for a deeper walkthrough.
