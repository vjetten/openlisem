===INSERT FIGURE TWO LAYER===

===The bed load equations===
link to [transport capacity bed load](https://github.com/vjetten/openlisem/wiki/Bed-load-transport-capacity)

A first bed load sediment transport layer is modelled right above the soil layer. This depth of this layer depends both on depth averaged flow velocity and flow depth. Above this layer, fully suspended sediment transport takes place troughout the rest of the water depth. The definition of the bed load layer is not generally agreed upon. Most definitions however agree that the sediment concentration within the layer shuold be high enough to make grain-grain interactions important when compared the grain-fluid interactions. The transport in this layer takes the form of rolling, sliding and saltating (Bagnold, 1956). Transport capacity for these layers are differently dependent on flow velocity and depth. The bed and suspended layers show more different behaviour. The settling velocity of sediment plays, for example, and important role in suspended sediment, while this is not the case for bed load. Diffusions of sediment furthermore only exists within the suspended sediment layer. The detachment and deposition process for a 2 layer sediment system is modelled in the same way as for overland flow. 

### Bed layer thickness
The bed load layer thickness is an estimation of the maximum height that bed load particles reach. To calculate the thickness of the bed load layer,  the expression from Hu and Hui (1996) for rough bed forms is used.

$$\delta=\ D\ (\ 1.78\ \left(\frac{\rho_s}{\rho_w}\right)^{0.86}\tau_\ast^{0.69})$$

with  
$\delta$ the bed load layer thickness $(m)$  
and  
$$\tau_\ast=\frac{{u_\ast^\prime}^2}{\frac{\rho_{s\ }{-\rho}_w\ \ }{\rho_w\ }g\ D}$$

with  
$\tau_\ast$ the shear stress $(Pa)$  
$g$ the gravitational accaleration $(m\ s^2)$  
and  
$$u_\ast^{\prime}=\ \sqrt{ g\frac{u}{C_c^\prime}}$$

with  
$u_\ast^\prime$ the effective shear velocity related to grains $(m\ s^{-1})$  
$u$ the depth averaged flow velocity $(m\ s^{-1})$  

$$C_c^\prime=18\ Log\left(4\ \frac{R_b}{D_{90}}\right)$$

with  
$C^\prime$ the chezy coefficient related to grains $(m^{0.5}\ s^{-1})$  
$R_b$ the hydraulic radius related to the bed $(m)$  
$D_{90}$ the grain diameter for which 90 % of the sediment mass has a smaller grain diameter $(m)$
