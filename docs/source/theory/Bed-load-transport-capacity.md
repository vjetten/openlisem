Besides suspended sediment also the transport capacity of the bed load is needed. The available equations in OpenLISEM as described below.


## Van Rijn simplified bed load
Van Rijn (1984a) derived a simplified empiric equation for bed load sediment transport. For both the bed load and suspended load equation, grain sizes between 50 and 2000 $\mu m$ were used for the derivation and calibration.

$$Q_{s,bl}=0.005\ \rho_s\ u\ h\ \left(\frac{D_{50}}{h}\right)^{1.2\ }M_e^{2.4}$$

with  
$Q_{s,bl}$ the sediment transport rate per unit width $(kg\ m^{-1}{\ s}^{-1})$  

and  
$$M_e=\frac{u-u_{cr}}{\sqrt{\left(\frac{\rho_s}{\rho_w}-1\right)g{\ D}_{50}}}$$

with  
$M_e$ the mobility parameter $(-)$  

and  
$$u_{\ast,cr}=\ \left(\frac{\rho_s}{\rho_w}-1\right)g\ D_{50}\ \theta_{cr}$$

with  
$u_{cr}$ the critical shear velocity for initiation of motion $(m\ s^{-1})$

and 
<p align="center">
$\theta_{cr}=0.24\ D_{\ast}^{-1}\ \ \ \ \ D_\ast\le4$  
</p>
<p align="center">
$\theta_{cr}=0.14\ D_{\ast}^{-0.64}\ \ \ \  {4\lt D}_\ast\le 10$ 
</p>
<p align="center">
$\theta_{cr}=0.04\ D_\ast^{-0.10}\ \ \ \   {10\lt D}_\ast\le 20$  
</p>
<p align="center">
$\theta_{cr}=0.013\ D_\ast^{-0.29}\ \ \ \   {20\lt D}_\ast\le 150$  
</p>
<p align="center">
$\theta_{cr}=0.055 \ \ \ \  {150 \lt D}_\ast$  
</p> 

with  
$\theta_{cr}$ the critical shields number $(-)$  
and   
$$D_\ast=D_{50}\sqrt[3]{\frac{\left(\frac{\rho_s}{\rho_w}-1\right)\ g}{v^2}}$$  

with  
$D_\ast$ the grain size parameter $(-)$  

The critical shear velocity, which is a dimensional version of the critical shealds number, indicates the minimal stress at the soil surface which is required for the initiation of sediment particle motion. When the shear velocity is higher than this threshold, sediment transport will occurr. The remaining energy is represented by the mobility parameter.