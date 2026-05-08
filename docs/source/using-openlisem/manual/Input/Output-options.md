<p align="center">
  <img width="782" height="226" src="../../../../images/Output_options_1.png">
</p>

## Catchment and outlet totals

The total values in mm and m3, the peak values and runoff fraction for the entire flow domain and the defined are written to the outlets. Basically this is the information that can be seen at the left hand side of [the model results page](../../../quick-start/user-interface.md).  
The total values of the flow domain can also be saved per timestep.  
The hydrographs and sedigraphs for each defined outlet point are saved in separate files or as columns in one large file. Discharge, water height, and in case of erosion sediment flux and sediment concetration are saved.  The area average rainfall is added as a column. Als the total domain outflow from channels and overland flow is added.  
A satellite image (or airphoto or even an image of a topographic map) can be shown as background in the model output. This can be a geotiff image in a different resolution but it is loaded on memory and occupies memory space. It has to have the same cartographic projection and bounding coordinates as the PCRaster maps.

## Units and format
Here the output can be formatted to some extent: text files are comma delimited of space/tabs between values.
A timestamp can be added to the output filenames so that scenario results can be easily recognized.  
Outlet hydrographs can be saved as separated files (with the integer number of the outlet point), or as one large file where each outlet has several columns.   
Units for flow and erosion can be chosen for the screen and file output.  
**NOTE:** erosion values are displayed as positive (detachment, added to the flow) or negative (deposition, subtracted from the flow).

## Maps and mapseries

**Map output**  
The values of the last timestep are saved for the main hydrological, hydraulic and erosion processes. 
Hydrological maps of the totals are in mm. The maximum values display the maximum water height, velocity or momentum obtained of the entire run in a cell.  
NOTE: the these maxima are not reached at the same moment, as water travels downstream/downslope. So the maximum value maps do not display a situation at a given moment in time. Values in the soil loss map can be positive (net erosion) or negative (net deposition).

**Timeseries**  
Maps can also be saved for every timestep, or for every _X_ timesteps. The maps are saved in a PCRaster format: name0000.001, name0000.002 etc, where 'name' is the name given in the interface ('ro','int' etc.) and the extension number is the timestep number. In PCRaster/NutShell these maps can be displayed as a movie.


<p align="center">
  <img width="682" height="237" src="../../../../images/Output_options_maps.png">
</p>

