When the cross section of flow is above the range of overland flow, sediment will diffuse both trough small turbulent flows, and trough brownian motion. For suspended sediment layers in both flooding water and channels, diffusion is modelled. The devision between bed load and suspended sediment will be described later. For this, the differential equation governing sediment flow was rewritten to include second order diffusive terms.

$$\frac{dS}{dt}+\frac{d{(Q}_x\ C)}{dx}+\ \frac{d(Q_y\ C)}{dy}+\frac{d^2(\varepsilon\ C)}{dx^2}+\frac{d^2(\varepsilon\ C)}{dy^2}=dep-det$$

with  
$\epsilon$ the diffusion coefficient $(m^2s^{-1})$
 
The diffusion coefficient is closely related to the turbulence of the fluid. Turbulent currents, which are not described in many physical models due to their scale, cause a sediment flux in the direction of decreasing concentration (Tsujimoto, 2010). Estimation of the diffusion coefficient is difficult process. In large-scale ocean and atmospheric models, a two variable $k-\ \epsilon$ model is often used to describe both the sources of turbulent currents, and the transport of turbulent energy. LISEM uses a simpler implementation. First we note that the diffusion coefficient is related to the turbulent viscosity by a parameter that describes the ratio between water and sediment turbulent transport.

$\epsilon=\frac{\nu_t}{\sigma}$

with  
$\nu_t$ the turbulent viscosity $(m^2s^{-1})$  
$\sigma$ the turbulent Prandtl-Smith number $(-)$  
This parameter, the turbulent Prandtl-Smith number, generally takes values between 0.5 and 1 (Wu, 2001). Furthermore, a simple description for eddy viscosity by Smagorinksy (1964) is included, based on velocity gradients.

$$\epsilon=\frac{dx\ dy}{\sigma}\sqrt{\frac{du_x}{dx}+\frac{du_y}{dy}+\frac{1}{2}\frac{du_y}{dx}\frac{du_x}{dy}}$$
