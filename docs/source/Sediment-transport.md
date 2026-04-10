## Sediment flow equations  
Transport of sediment takes place when water with a sediment concentration is transported. To model sediment flow in the 1 dimensional kinematic wave, a simple flow advection scheme is used.

$$\frac{dS}{dt}+\frac{d(Q\ C)}{dx}=dep-det$$

with  
$S$ the sediment load $(kg)$  
$C$ the sediment concentration $(kg\ m^{-3})$  
$dep$ the deposition $(kg\ s^{-1})$  
$det$ the detachment $(kg\ s^{-1})$  

In order to implement sediment transport in the 2 dimensional kinematic wave and the saint-venaint equations for flooding, the sediment transport was similarly rewritten to 2 dimensions.

$$\frac{dS}{dt}+\frac{d{(Q}_x\ C)}{dx}+\ \frac{d(Q_y\ C)}{dy}=dep-det$$
