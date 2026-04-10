Two equations are available to calculate the transport capacity of overland flow in OpenLISEM. The default is the equation by Govers (19900, besides that also a equation by Hairsine and Rose (yyyy) is available. The equation by Govers (1990) was empirically derived from measurements, and is dependent on stream power. Hessel & Jetten (2007) applied 8 different transport capacity equations to a small catchment in the Chinese Loess plateau. This extensive test showed that the transport capacity of Govers (1990) performed best due to its small sensitivity to slope and grain size. The transport capacity represents the concentration for which the sediment deposition and detachment are equal, and the concentration is thus stable, based on this [detachment or deposition](https://github.com/vjetten/openlisem/wiki/Detachment-&-Deposition) can be calculated.

***

## Govers 

$$T=\rho_s\ c\ \left(\omega\ -\ \omega_{cr}\right)^d$$  

with  
$\omega$ the stream power $(m\ s^{-1})$  
$\omega_{cr}$ the critical stream power $(m\ s^{-1})$$  
$\rho_s$ is the density of the sediment material $(kg\ m^{-3})$  
and  

$$c={\frac{{(D}_{50}+\ 5)}{0.32}}^{-0.6}$$  

$$d={\frac{{(D}_{50}+\ 5)}{300}}^{0.25}$$

with  
$D_{50}$ the median grain diameter $(m)$  

## Hairsine and Rose

$$T = \frac{D_{50}}{w_s} \cdot \frac{0.013}{g} \cdot 1.650 \cdot \frac{U\ S-0.004}{h}$$

with  
$D_{50}$ the median grain diameter $(m)$   
$w_s$ the [settling velocity](https://github.com/vjetten/openlisem/wiki/Settling-velocity) (terminal velocity of the particle) $(m\ s^{-1})$  
$g$ the gravitational force  
$U$ the flow velocity $(m s^{-1})$  
$S$ the slope of the surface $(-)$   
$h$ the water height $(m)$   


