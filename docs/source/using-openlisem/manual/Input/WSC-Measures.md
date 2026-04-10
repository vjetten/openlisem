With OpenLISEM many different water retention and soil conservation measures can be simulated. These measures are often included in model scenario studies to understand and quantify the effect of the implementation of different measures on runoff, flooding and soil erosion.

In the below description we aim to simulate each measure with a process based approach. The parameters that are changed by a measure are adapted in the model input.

## Water retention measures
Water retention measures mostly influence one or more of three main processes:
1. increased infiltration
1. increased storage
1. decreased flow speed

In the section below multiple measures are described with examples and explanation how to implement these into OpenLISEM. Often other measures can be implemented by adapting the approach used for the examples below.

### Additional surface storage
A measure to reduce runoff form hill slopes is to increase the storage on the soil surface. In OpenLISEM surface storage is described with the 'random roughness' (RR - cm) which is defined as the standard deviation of roughness within a raster cell perpendicular to the slope. An example of measures to increase the surface storage can be in-furrow micro dams in potato cultivation. These micro dams increase the surface storage in a structured way. 

[[/images/micro-dams_example.png|height = 500px]]  
*Micro-dams containing water, retrieved from [fiwap.be](https://fiwap.be/documentation/du-ruissellement-au-travail-du-sol/)*

A method to simulate this in OpenLISEM is adjusting the input RR for the additional storage. The RR is used in OpenLISEM to calculate the micro-depression storage which is corrected for slope [📖](https://github.com/vjetten/openlisem/wiki/Interception#storage). To achieve this a map with the additional surface storage depth (cm) by the implemented measures is required. This can be calculated by estimating the volume of water that can be stored by the measure on a square meter when the slope of the field is 0%. A PCRaster script is developed to recalculate this to RR, this script can be found [here](files/PCRscripts/Additional_roughness_to_RR.mod).

> :bulb: **Example:**  
The micro dams in the image above can each store 8 liters of water on a flat slope. We estimate 3 micro dams per square meter. This means the additional storage is 8 * 3 = 24 liter per square meter or 2.4 cm additional storage depth.

### Swales

Swales are ditch and dike features along the contour of the hillslope which aim to increase the water storage. The soil taken from the ditch is added as a dike downstream to create a buffer area for water the flows down the field.

[[/images/swales_example.png|height = 500px]]  
*Installed swale that captured runoff, the dike can be vegetated. source to be added!*

In OpenLISEM the best implementation of swales depends on the soimulated resolution and the used numerical solution for overland flow. 
1D kinematic wave does not alow for sinks in the landscape, so the buffer capacity of the swales cannot be modelled in the DEM. With the 2D dynamic wave sinks are possible and will be filled with water before further runoff occurs.

When a high resolution raster grid is used in OpenLISEM, the dimensions of the swales roughly correspond with the cell sizes. In this case the swales can be modelled by lowering some cells to simulate the ditch and increasing the height of the next downstream cell to simulate the dike. An example PCraster script was made to implement this [option](files/PCRscripts/Swales_buffers.mod). This script has the following functionality:
- swales will be placed at a given elevation interval along the contours.
- swales can be located in a specific landuse class or over the whole area.
- the height of the dike and the depth of the ditch can be given, the width of both equals the cell size.
- the hydraulic conductivity in the ditch can be adjusted.

The script will produce a 'buffers.map' which can be added in OpenLISEM under the 'mitigation measures water' section of the input.

:construction: the information below on this page has to be updated (2025-03-20)

When the used resolution is lower and cell sizes are to large to simulate the swales directly in the DEM storage has to be simulated differently.


## Soil conservation measures

**Sediment trap**: cells that trap all sediment (to simulate e.g. hedgerows or sediment traps) but not water. This uses the map `sedretmax.map`, which has the maximum sediment retention in $m^3$. The bulk density value is used to convert $kg$ deposition to $m^3$). There is no splash and flow detachment, and the water is slowed down by the Manning's n factor given in the interface. When the trap is full the Manning’s n returns to normal. 
 
**Grass strip**: defined as the width within cells with a grass strip. Be aware of diagonal LDD connections crossing corners of the grass strip when using the kinematic wave! The Manning's n of a grass strip slows down the water, 0.2 is a high value for long grass. It is assumed that grass strips have no detachment but only sedimentation. 