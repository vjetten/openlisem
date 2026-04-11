## The user interface

The OpenLISEM interface has 3 pages, see the figures below:  
1. Has the main categories for running the model. On the left side are the INPUT definitions and options, to the right the OUTPUT definitions and options.  
1. A tree structure that shows the names of all input maps. You can click on each name and change names to your own names, but it is far more convenient to use default names and make sure that different scenarios are stored in different folders.  
1. Model results interface that shows summarized data, hydrographs and spatial output while running.  

The main interface of the software shows the model settings. Here, the model controls, process settings and output settings are visible. 

<p align="center">
  <img width="850" height="540" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/OpenLISEM_input_interface.png">
</p>

When opening the map database, an overview of required input maps is available. Double-clicking a map name allows you to select a file, although it is usually easier to use the default filenames. The lighter text indicates an input map is not required for the current simulation, as the relevant process has not been activated. Input maps are usually constructed using the PCRaster + Nutshell software -or the new LISEM database generator-. All input maps must be of identical rows and columns, but might be of any format supported by GDAL (.map, GeoTiff, ascii, etc..).

<p align="center">
  <img width="818" height="597" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/OpenLISEM_mapDB.png">
</p>

When the simulation is running, real-time output can be visualized in the simulation tab. Here, hydrological statistics of the simulation are provided in catchment-average mm equivalent. Additionally, a spatial view of many model variables can be shown, such as infiltration, flow heights, flow velocities and soil loss. For specified points in the simulation, such as channel outlets or bridges, a discharge curve (hydrograph) can be viewed during the simulation.

<p align="center">
  <img width="782" height="503" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/OpenLISEM_output_interface.png">
</p>
