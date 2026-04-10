## 1D and 2D flow processes

OpenLISEM has two types of overland flow: kinematic flow and dynamic flow. The user can choose one of the following flow schemes: (1) a [kinematic scheme](https://github.com/vjetten/openlisem/wiki/Flow:-1D-kinematic-wave) for runoff without channel flow or (2) kinematic wave runoff that is connected directly to [channel flow]() or to a flooded area, whereby the flood is simulated with a dynamic scheme and comes from the overflowing channel. A lastly (3) scheme where all overland flow is simulated with a [dynamic flow scheme]() interacting with channels (there is no distinction between flooding and runoff). In the second option, the flood by dynamic flow is always coupled to the presence of a channel. Channel flow is always a kinematic wave along the channel network. Reasons to have two flow solutions in the model are partly historic (backward compatibility) partly educational (explain effects of choices) and partly because different simulation problems require different solutions.

***

Both flow types can be derived from the Navier-Stokes equations for in-compressible water flow, with simplifications. The Saint-Venant approximation for shallow flow with a depth average velocity, is set for the dynamic flow scheme. The mass and momentum balance equation from this approximation are shown below. The frictional force is estimated from the Darcy-Weisbach law (Chow, 1959). Due to the incorporation of dynamic terms, this approximation is applicable to a wide range of flow types, including flash floods and other scenarios with high pressure forces and velocities.

$$\frac{\partial h}{\partial t}+\ \frac{\partial\left(hu_x\right)}{\partial x}+\ \frac{\partial\left(hu_y\right)}{\partial y}=\ R-I$$

where 
$h$ is the flow height $(m)$  
$u$ is the flow velocity $(m\ s^{-1})$  
$R$ is the rainfall $(m)$  
$I$ is the infiltration $(m)$.

$$\frac{\partial h u_x}{\partial t}+\ \frac{\partial(hu_x^2+\frac{1}{2}gh^2)}{\partial x}+\ \frac{\partial(hu_xu_y)}{\partial y}=\ gh(S_x-n^2\frac{u_x\left|\vec{u}\right|}{h})$$

$$\frac{\partial h u_y}{\partial t}+\ \frac{\partial(hu_y^2+\frac{1}{2}gh^2)}{\partial y}+\frac{\partial(hu_xu_y)}{\partial x}=gh(S_y-n^2\frac{u_y\left|\vec{u}\right|}{h})$$

where  
$g$ is the gravitational acceleration $(m\ s^{-2})$  
$S$ is the friction term $(m\ s^{-2})$  
$S_f$ is the momentum source term $(m\ s^{-2})$  
$n$ is Manning’s n friction coefficient $(s\ m^{-\frac{1}{3}})$.

When this approximation is further simplified, and the inertial term is neglected, the diffusive wave is derived. Finally, when pressure terms are furthermore ignored, the kinematic wave remains. In this approximation, the local velocity is, at any moment, instantly determined by a balance between the gravitation force and the frictional force:

$$0=\left(S_{f_x}-S_x\right)$$

$$0=\left(S_{f_y}-S_y\right)$$

$$u=R^\frac{3}{2}\frac{\sqrt S}{n}$$

where  
$u$ is the flow velocity $(m\ s^{-1})$  
$R$ is the hydraulic radius $(m)$  
$n$ is the Mannings coefficient of the surface $(s\ m^{-\frac{1}{3}})$. 
