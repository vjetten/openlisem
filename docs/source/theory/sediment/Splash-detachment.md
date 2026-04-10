Splash detachment is the process were incoming precipitation detaches soil on impact with the soil surface, due to its kinetic energy. In the case of a large rainfall event, the magnitude of this form of detachment is insignificant when compared to the flow detachment. For smaller rainfall events, this form of erosion makes a distinct difference. The rainfall intensity that remains after interception is used for splash detachment. The splash detachment is calculated with:

$$dep_{splash}=\ A_s\cdot Ke\cdot e^{-1.48\ h} \cdot P_h \frac{A}{dt}$$

with  
$dep_{splash}$ the splash deposition rate $(kg\ m^{-2}s^{-1})$  
$A_s$ the aggregate strength $g\ J^{-1}$  
$Ke$ the kinetic energy of the rainfall or throughfall $(\ J\ m^{-2}\ {mm}^{-1})$  
$P_h$ the rainfall or throughfall $(mm)$  
$A$ the surface area where the splash detachment takes place $(m^2)$  

If the aggregate strength is not known it can be estimated with the following general equation, which has been derived based on Lowe splash tests (unpublished data):

$$A_s = 5.3361\ A_d^{-0.238}$$

with  
$A_d$ median number of drops to decrease the aggregate mass by 50%, based on the Lowe test $(-)$  

An increasing water height causes splash detachment to decrease since the kinetic energy of the rainfall dissipates in the layer of water. The inverse exponential dependency of splash detachment on the water height accomplishes this effect. The kinematic energy of the incoming water is calculated separately for rainfall and throughfall. For throughfall we use:

$$K_{e,t}=15.3\ \sqrt{h_{veg}}-5.87\$$

with  
$h_{veg}$ the vegetation height $(m)$  
$K_{e,r}$ the kinetic energy of direct rainfall $(\ J\ m^{-2}\ {mm}^{-1})$  
$K_{e,t}$ the kinetic energy of the throughfall from the vegetation $(\ J\ m^{-2}\ {mm}^{-1})$  

To calculate the kinetic energy of the incoming rainfall several different equations are available. The equations describe non-linear relations between rainfall intensity $(mm\ h^{-1})$ and kinetic energy $(J\ mm^{-1})$:

1. $K_{e,r} = a(1-b\ e^{-c\ P_i)})$ ; A negative exponential equation with average global parameters (van Dijk 2002). values are a= 28.3, b =0.52, c=0.042.
1. $K_{e,r}=\ a+b\ Log\left(P_i\right)$ ;  A USLE type erquation with parameters as found in the EUROSEM model (Morgan et al., 1998). a = 8.950, b=8.440
1. $K_{e,r}=\ a\ P_i^b$ ; A decreasing power law (Sanchez, 2012). a = 7.6, b = 0.22. This equation is available in an intensity based and timebased version.

with  
$P_i$ the rainfall intensity $(mm\ h^{-1})$.

The ponded area, which is calculated using [the micro depression storage](https://github.com/vjetten/openlisem/wiki/Interception), is used to calculate a dry splash detachment, and a wet splash detachment. This changes the $P$ and $A$ for both these surfaces. The effective water heights for non-ponded areas is thus 0, and this water volume is divided over the ponded area of a cell. Finally the total splash detachment is taken from the cell’s soil layer. Other surfaces such as compacted soil, crusted soil and roads are taken into account to show no splash detachment. They similarly cause a change in both $A$ and $P$.

$$P_{nonroad}=\left(1-f_{road}\right)\ P$$
$$A_{nonroad}=\left(1-f_{road}\right)\ A$$

with  
$P_{nonroad}$ the rainfall on non-road surfaces $(mm)$
$A_{nonroad}$ the area on non-road surfaces $(m^2)$
$f_{road}$ the fraction of the cell that is covered by roads $(-)$

A splash delivery fraction is used to estimate the sediment splashed from dry parts in a cell to ponded areas in a cell. The ratio is cellsize dependent since splash distances are 50 cm at most.