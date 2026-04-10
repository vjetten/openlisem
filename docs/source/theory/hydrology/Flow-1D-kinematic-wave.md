When precipitation is neither intercepted, stored nor infiltrated, overland flow is modelled. From conservation of mass, the kinematic wave equation for overland flow can be derived. The change in cross section of flow, should be equal to the spatial derivative of the discharge to maintain water volume. Combining this with external sources of water such as rainfall and infiltration gives equation 21. To simulate flow of runoff water over the digital elevation model, this kinematic wave equation is used, together with Mannings equation for flow velocity. In its 1 dimensional form, the kinematic wave routes the water over a local drainage network (LDD).

$$\frac{dA}{dt}+\frac{dQ}{dx}=q-i$$

with
$Q$ the discharge $(m^3\ s^{-1})$  
$A$ the cross section of the flow $(m^2)$  
$i$ the infiltration $(m^3\ s^{-1})$  
$q$ the other sources of water (rainfall and snowmelt)\ $(m^3\ s^{-1})$  

$$u=R^\frac{3}{2}\ast\frac{\sqrt S}{n}$$

with   
$u$ the flow velocity $(m\ {ms}^{-1})$  
$R$ the hydraulic radius $(m)$  
$n$ the Mannings coefficient of the surface $(s\ m^{-\frac{1}{3}})$  

The relation between the cross section of the overland flow, and the discharge is given by a power function (Chow, 1959 & Chow 1988).
$$A=\alpha Q^{\beta}$$  
with   
$\beta$ a coefficient $(= 0.6) (-)$  
and  
$$\alpha=\left(\left(\frac{n}{s^{0.5}R}\right)^\frac{2}{3}\right)^\beta$$
with   
$P$ the wetted perimeter of the flow $(m^2)$  

Earlier versions of OpenLISEM could only use a LDD map to route water through the modelled area (see figure). This allowed to implement the kinematic wave in 1 dimension. In the usage of a local drainage direction network, multiple assumptions are made:
* The directions for overland flow are restricted to 9 possible values.
* The discharge of a cell always fully flows into the single next connected cell

<p align="center">
  <img width="594" height="345" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/LDD_1Dflow.png">
</p>

These assumptions can lead to unrealistic behavior in certain situations. Diagonal flow travels greater distances than horizontal or vertical flow. Furthermore, the catchment size can cause enormous amounts of water to be routed to a single cell. With large catchments, or small cells, this causes overland runoff water heights that would in reality spread out of multiple cells. Due to this, the discharge values become much higher than is realistic. Since the soil surface area beneath the flow is also smaller, infiltration is generally lower than would be realistic. This can have an impact on especially erosion values, since the total transport capacity of overland flow does not relate linearly to water height.
