## Soil properties

To prepare the soil properties, a soil type unit map and a soil properties table are needed. Similar to the land use properties, the [`lookupscalar()`](https://github.com/vjetten/openlisem/wiki/Introduction-PCRaster-&-Nutshell#lookup) function can be used to create maps with soil properties.

<p align="center">
  <img width="597" height="300" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/prepare_soil_input.png">
</p>

The soil physical values in the table come from a software package [SPAW](https://hrsl.ba.ars.usda.gov/soilwater/Index.htm) (Saxton et al., 2005).  Saxton and Rawls (2006) created pedotransfer functions for all soil samples in the USDA databases, to transform texture information to soil physical parameters. Besides their data, median grain sizes are needed. These can be estimated from the USDA soil classification (USSCS, 1975). The data for a small collection of soils, together with needed grain sizes, are furthermore provided below:

_Table 1: properties related to texture class_

Soil Type   (Saxton et al., 2005) |   | ksat | porosity | psi
-- | -- | -- | -- | --
| |   | mm/h | cm3/cm3 | cm
text class | 0 | 1 | 2 | 3
C | 1 | 2.5 | 0.5 | 50.0
CL | 2 | 4.2 | 0.5 | 50.0
L | 3 | 18.5 | 0.5 | 40.0
S | 4 | 112.0 | 0.5 | 20.0
SaCL | 5 | 7.0 | 0.4 | 35.0
SaL | 6 | 80.0 | 0.5 | 40.0
Si | 7 | 42.0 | 0.5 | 40.0
SiC | 8 | 13.0 | 0.6 | 40.0
SiCL | 9 | 25.0 | 0.5 | 40.0
SiL | 10 | 17.4 | 0.5 | 40.0
Water (W) | 20 | 0.0 | 0.0 | 0.0
Urban (A) | 21 | 0.0 | 0.2 | 40.0
Salt pans (m) | 22 | 0.0 | 0.4 | 40.0
Rock/outcrops | 23 | 0.0 | 0.1 | 40.0

## Soil water content

While OpenLISEM can be used to simulate ground water flow for a period before an event, sometimes it is easier to provide an initial soil water content, and start the simulation at the starting point of the event. Several approaches can be taken to acquire an estimate of ground water content for the initial state of a simulation.

# Soil depth

Soil depth input is a crucial parameter for the infiltration potential. OpenLISEM input uses a soil depth and volumetric relative water content. Thus, by altering the soil depth, the infiltration potential is altered. In reality, spatial soil depth patterns are unknown, and detailed measurements of this parameter are not common. There are several approaches that can be taken to predict soil depth patterns spatially for use in OpenLISEM.

**Homogeneous/large-scale estimations**  
Homogeneous soil depth can be chosen for rough estimation of the land surface processes. These homogeneous values can be based on several measurements or visible depth of the weathered material (visible depth of weathered material above bed rock can be visible in the case of landslides). The exact depth that provides good representative behavior in the OpenLISEM is difficult to estimate and must therefore be based on calibration and validation. 
In the case of large scale area, a homogeneous depth value might be improved by using large-scale estimations. Besides national databases of soil depth patterns, a global dataset is provided by [soilgrids.org](https://soilgrids.org/). For larger areas, this dataset provides a good indication of absolute soil depth. Generally, the depth is elevation, slope and accumulated area depended, but is in reality based on automated learning algorithms. However, the spatial resolution of this dataset is 250 meters. It thus does not react to individual streams, channels and slopes.

**Soil Depth Equations**  
Several equations exist based on considerations of the processes involved in rock weathering, soil production, and soil transport. Typically, weathering increases soil depth, after which transport takes place driven by gravity and hydrology. Physical laws in slope stability indicate that the possible stable depth of soil increases with gentler slopes. Then, erosive processes accumulate near channels, flood plains and lower elevation. All of these considerations have led to a variety of soil depth equations in literature. Below are several examples from Soulnier et al. (1992).

$$Hs = Hs_{min} + (Hs_{max} - Hs_{min} ) \cdot \frac{z-z_{min}}{z_{max}-z_{min}}$$

$$Hs= Hs_{min} + ( Hs_{max} - Hs_{min} ) \cdot \frac{s-s_{min}}{s_{max}-s_{min}}$$  

$$Hs = Hs_{min} + ( Hs_{max} - Hs_{min} ) \cdot \frac{dxChan - dxChan_{min}}{dxChan_{max} - dxChan_{min}}$$  


The input requirements for these equations are much easier to obtain when compared to spatial soil depth data. The minimum and maximum soil depth in the area is needed. The other parameters are all based on elevation model calculations. 

