# Interception

At the start of any rain-driven event, the majority of the available precipitation is intercepted by surfaces other than the soil (Linsley et al., 1982). Interception of rainfall can be calculated for vegetation, buildings and raindrums. Evaporation, which can take place from these surfaces, is not modelled within LISEM because of two reasons. First, the event-based nature of the model allows for the assumption that slow processes, such as evaporation, can be neglected. Evapotranspiration is furthermore minimal during the rainfall events that are typically modelled in LISEM, since cloud cover is generally high during these events. Interception is thus modelled as a fixed storage that takes from the precipitation before that reaches the soil layer. Rainfall that is not intercepted reaches the soil with the same intensity as the rainfall. While this is not the case in reality, trough fall intensities are only known for few tree types.  Actual canopy interception is given by (Aston, 1979).

$$I_c = S_{max}\ (1-e^{-k\frac{P_{cum}}{S_{max}}})$$

with  
$I_c$ the total intercepted storage at a given time $(mm)$.  
$S_{max}$ the maximum canopy storage $(mm)$.  
$P_(cum)$ the total precipitation $(mm)$.  
and  

$$k = 1 - e^{-(co\ LAI)}$$

with  
$k$ a parameter related to canopy openness $(-)$.  
$co$ the canopy openness $(-)$.  
$LAI$ the leaf area index $(-)$.  

For the maximum storage, equations for several tree types, depending on leaf area index, were found by Von Hoyningen-Huene (1981). The equations that are implemented within LISEM are:  

$S_{max}=0.935+0.498\ LAI-0.00575\ LAI^2$ (Crops)  
$S_{max}=0.2331\ LAI$                      (Pinus)  
$S_{max}=0.3165\ LAI$ (Douglas)  
$S_{max}=1.46\ LAI^{0.56\ }$ (Olive)  
$S_{max}=0.0918\ LAI^{1.04}$ (Eucalypt)  
$S_{max}=0.2856\ LAI$ (Broadleaved Forest)  
$S_{max}=0.1713\ LAI$ (Bracken)  
$S_{max}=0.59\ LAI^{0.88}$ (Clumped Grass)  

Interception by roofs and raindrums is also modelled. The fraction of the rainfall that hits an area covered by these types of surfaces is stored and does not reach the soil surface. When the maximum raindrum or roof storage is reached, any extra rainfall hits the soil surface.

# Storage

Rainfall is first stored in micro depressions in the soil surface. When the water level in these depressions increases, runoff starts. At a micro-scale, runoff is a spatial process of ponds that fill up and overflow into each other. These ponds release little runoff until they fully overflow. Because of this runoff flow does not start immediately after the first rainfall hits the soil within LISEM. To estimate the fraction of water that is stored in these depressions, and the fraction that is used for runoff, the surface roughness is used to estimate the Micro Depression Storage (MDS). The equation for the MDS was determined by Kamphorst et al. (2000) from 221 digital elevation models of various types of micro relief, in a wide variety of agricultural circumstances and soil types. The analysis is based on Digital Elevation Models (DEMs) of with a spatial resolution of roughly $1\ m^2$.

$$MDS=0.243\ RR+0.010\ RR^2+\ 0.012RR\ S$$

with  
$MDS$ the micro depressional storage $(m)$  
$S$ the slope $(m\ m^{-1})$  
$RR$ the standard deviation of the surface heights $(mm)$.  

The flow width for runoff is furthermore changed depending on the estimated ponded area. The ponded fraction of a cell is given by (Jetten and De Roo, 2001):  
$$f_{pa}=1-e^{-a\ (\ h)}$$  

with
$f_{pa}$ the fraction of the cell area that is covered by ponds $(-)$  
$h$ the average water depth $(mm)$  
and  
$$a=1.406\ RR^{-0.924}$$  

It is assumed that after the water volume in a cell reaches 10 percent of the micro depressions storage, the volume of water that is used for runoff is given by:

$${\ h}_{runoff}=max{\left(0.0,\ \left(\ h-SDS\right)\ast\ \left(1-e^{-h\ \frac{h-SDS}{MDS-SDS}}\right)\right)}$$

with  
$h_{runoff}$ the height of the water that is used for runoff $(mm)$  
$SDS$ the water height at which runoff starts $(=0.1\ MDS ) (mm)$  

When the MDS is completely filled, all remaining water is used for runoff. 
