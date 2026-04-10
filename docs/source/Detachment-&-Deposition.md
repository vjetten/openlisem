Detachment and deposition are based on the [transport capacity](https://github.com/vjetten/openlisem/wiki/Transport-capacity), detachment efficiency and [settling velocity](https://github.com/vjetten/openlisem/wiki/Settling-Velocity) for the median grain diameter. The deposition rate is proportional to the settling velocity. 

$$dep=\ \ w_s\ B\ C$$  

with  
$B$ the flow width $(m)$  
$dep$ the deposition rate $(kg\ m^{-2}s^{-1})$  

The transport capacity represents the concentration for which the sediment deposition and detachment are equal, and the concentration is thus stable. We can use this to adapt the equation for deposition.

$$dep=\ w_sB\ min{\left(0.0,\ \left(T-C\right)\right)}$$

It can be assumed that detachment follows the same form (Rauws & and Govers, 1988), with the addition of an erosion efficiency factor when detachment takes place.

$$det=\ {\gamma\ w}_sB\ max{\left(0.0,\ \left(T-C\right)\right)}$$  

with  
$det$ the detachment rate $(kg\ m^{-2}s^{-1})$  
$\gamma$ the erosion efficiency coefficient $(-)$  

The erosion efficiency coefficient is based on soil cohesion and root strength, which provides extra soil cohesion. OpenLISEM provides three different equations to calculate detachment efficiency: Rauws & Govers (1988), the equation from the EUROSEM model and the last the equation by Morgan Morgan and Finney. The equations by Govers where originally meant to describe rill erosion. As such, all erosion within OpenLISEM either is part of splash detachment or flow detachment in the form of rills. For the simulation of intense rainfall events, which is usually the case, sheet erosion has an insignificant magnitude when compared to rill erosion (Herweg, 1996). The detachment efficiency, based on soil cohesion (in kPa),  is a value from 0 - 1 multiplied with the detachment rate. 0 gives no flow detachment, 1 gives max flow detachment. OpenLLISEM is very sensitive to the efficiency. Generally the Govers equation gives the most detachment and EUROSEM the lowest detachment.

***

## Govers detachment efficiency

$$\gamma=min{\left(1.0,\ \frac{1}{0.89+0.56\left(co+co_{veg}\right)}\right)}$$  

with  
$co$ the soil cohesion $(kPa)$  
$co_{veg}$ the extra soil cohesion due to vegetation $(kPa)$  


## EUROSEM detachment efficiency

$$\gamma = min{\left(1.0,\ 0.79 \cdot e^{-0.85\ (co+co_{veg})}\right)}$$

## MMF detachment efficiency

$$\gamma = min{\left(1.0,\ \frac{1.0}{2.0\ (co+co_{veg})}\right)}$$