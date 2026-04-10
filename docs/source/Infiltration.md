Infiltration is the process where water is transported downwards from the surface to the subsurface. Depending on the hydraulic conductivity and soil water content, surface water will seep into the soil. Several infiltration models are included in the OpenLISEM model:  
* The Green & Ampt infiltration model
* Smith & Parlange model
* SWATRE multilayered soil water model


These infiltration models use the empirical Darcy equation which describes a simple vertical soil water balance.

$$\frac{\partial\theta}{\partial t}=\ -K_s\frac{\partial h}{\partial z}$$  

with 
$\theta$ the soil moisture content $(m^3\ m^{-3})$  
$h$ the hydraulic head $(m)$  
$z$ the vertical elevation $(m)$  
$K_s$ the saturated conductivity $(m\ s^{-1})$  

## Green & Ampt
The Green & Ampt (1911) infiltration method assumes that a wetting front moves downwards into the soil layers parallel to the soil surface. Above this front, the soil is saturated, while beneath this front, the soil is completely dry. Green & Ampt stated that, when the water height above the soil surface is assumed to be zero, a simplification of the Darcy equation for vertical water flow can be used.

$$f\ =\ -K_s\ \left(\frac{h_f-h_0}{Z_f}\right)=\ K_s\ \left(\frac{\psi\ }{Z_f}+1\right)$$

with  
$f$ the infiltration rate $(m\ s^{-1})$  
$h_f$ the hydraulic head at the wetting front $(m)$  
$h_0$ the hydraulic head at the soil surface $(=0)\ (m)$  
$Z_f$ the depth of the wetting front $(m)$  
$\psi$ the matric pressure at the wetting front $(h=\psi+Z)\ (m)$  
and  
$$Z_f=\frac{F}{\theta_s-\theta_i}$$
with  
$F$ the cumulative infiltrated water $(m)$  
$\theta_s$ the porosity $(m^3\ m^{-3})$  
$\theta_i$ the initial soil moisture content $(m^3\ m^{-3})$  

The value of $\psi$ depends on the soil type. Using the Green & Ampt equations, and combining these, the final equation for infiltration rate can be acquired:
$$f\ =f_{pot}=\ -K_s\left(\psi\ \frac{\theta_s-\theta_i}{F}+1\right)$$  
with  
$f_{pot}$ the potential infiltration rate $(m\ s^{-1})$  

This method can be applied for both a 1 layer or 2 layer system. Beneath these layers, an open or closed boundary can be chosen.

### Effective $K_s$ in multilayered infiltration

When the soil is modelled as a two layer system, the effective $K_s$ has to be calculated based on the depth of the wetting front and the $K_s$ of respectively the first and second soil layer. This is done by calculating the harmonic mean over the wetting front depth, as proposed by [the GSSHA model](https://www.gsshawiki.com/Infiltration:Multi-layer_Green_and_Ampt).


### $\psi$ estimates

The value of $\psi$, described as the capillary pressure at the wetting front (Rawls 1983), has to be found emperically. If no accurate data for a specific location is available, the value can be estimate based on soil texture of $K_s$. In OpenLISEM a fit through the emperical data of Rawls et al. 1983 is used to estimate $\psi$ based on the give $K_s$:  

$$\psi = e^{-0.3382\ log(K_s) + 3.3425}$$

The $\psi$ is limited to never exceed the bubbling pressure.