Sediment transport in channels can not always assume the same transport capacity equations as overland flow. Therefor additional equations are available to calculate transport capacity for channels and flooded area's. This can done only for suspended sediment, or [bed load transport can also be simulated](https://github.com/vjetten/openlisem/wiki/Two-Layer-sediment-transport). Because the water height in channels and flooder area's is larger also [diffusion](https://github.com/vjetten/openlisem/wiki/Diffusion) can influence sediment transport. For transport capacity over suspended sediment, besides [Govers](https://github.com/vjetten/openlisem/wiki/Overland-sediment-transport-capacity#govers) three other equations are available: two versions of van Rijn's transport equations, a simplified and full version, and Engelund & Hansen. Besides that OpenLISEM also calculates the [bed load transport capacity](https://github.com/vjetten/openlisem/wiki/Bed-load-transport-capacity) when using 2 layer sediment transport.

## Channel and fooding transport capacity

### Van Rijn Simplified

Van Rijn (1984b) derived a simplified empiric equation for suspended load sediment transport.

$$Q_{s,ss}=0.008\ \rho_s\ v\ D_{50}{(D_\ast)}^{-0.6}M_e^{2.4}$$  

with  
$Q_{s,ss}$ the sediment transport rate per unit width $(kg\ m^{-1}{\ s}^{-1})$  

The simplified version of the van Rijn suspended sediment load equation, similarly to the bed load equation, uses the mobility parameter, based on the excess shear stress. 

### Van Rijn full

Van Rijn (1984b) derived an empiric equation for suspended load sediment transport based on flow velocity and height.

$$Q_{s,ss}=\rho_s\ F\ u\ h\ C_r\$$

with  
$Q_{s,ss}$ the sediment transport rate per unit width $(kg\ m^{-1}{\ s}^{-1})$
and   
$$F=\frac{{\left(\frac{h_r}{h}\right)^{Z^\prime}-\left(\frac{h_r}{h}\right)}^{1.2}}{\left(1-\frac{h_r}{h}\right)^{Z^\prime}\left(1.2-Z^\prime\right)}$$

with  
$F$ the correction factor for suspended load $(-)$  
and  
$$h_r=\ \ 0.5\Delta$$

with  
$h_r$ the reference height $(m)$  
$\Delta$ the bed form height $(m)$  
and  
$$Z^\prime=Z+\ \phi$$  
with  
$Z^\prime$ the corrected suspension number $(-)$  
$\phi$ the suspension number correction factor $(-)$
and  
$$Z=\frac{w_s}{\beta\ k\ u_\ast}$$  
with  
$Z$ the suspension number $(-)$  
$k$ the constant of Von Karman $(\approx0.41) (-)$  
and  
$$\phi=2.5\ \left(\frac{w_s}{u_\ast}\right)^{0.8}\left(\frac{C_r}{c_0}\right)^{0.4}$$

with  
$w_s$ the settling velocity for the median grain size $(m\ s^{-1})$  
$c_0$ the maximum bed concentration $(kg\ m^{-3})$  
and  
$$\beta=1+2\left(\frac{w_s}{u_\ast}\right)^2$$  
with  
$\beta$ the ratio between sediment diffusion and fluid diffusion $(-)$  
and  
$$u_\ast=\ \sqrt{g\ h\ S}$$
with  
$u_\ast$ the overall shear velocity $(m\ s^{-1})$  
$S$ the slope $(m\ m^{-1})$  
and  
$$C_r=0.015\frac{D_{50}}{h_r}\frac{T^{1.5}}{D_\ast^{0.3}}$$
with  
$C_r$ the reference concentration $(kg\ m^{-3})$  

The full equations for the van Rijn suspended sediment transport do not use a threshold value. Instead a suspension paramter is used, which represents the ratio between downwards gravitational forces and upwards turbulent forces on the sediment particles. This parameter is corrected for the influence that sediment has on the turbulence of a fluid. Then, based on a reference concentration at a reference height, the suspended load along the entire vertical profile is estimated and integrated, which is represented by the F-factor. The reference height is assumed to be the top of the bed load layer. Above this height, the suspended transport takes place, an below it, bed load transport takes place.



### Engelund & Hansen