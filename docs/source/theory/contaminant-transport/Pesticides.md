 🚧 **The first official version of OpenLISEM-pesticide is published [(Commelin et al., 2024)](https://doi.org/10.1016/j.envsoft.2024.105960). The module is not yet fully integrated in the graphical user interface, see the [manual](https://github.com/vjetten/openlisem/wiki/Pesticide-input) for further details.**

***

The movement of pesticides by runoff consists of two main processes, the uptake of the pesticide from the soil into the runoff water and the lateral transport with runoff over the soil surface. The uptake of the pesticides from the soil into runoff is often conceptualized with a ‘mixing layer’ (Havis et al., 1992). In OpenLISEM-pesticide we use a mass transfer model to simulate the uptake of dissolved pesticides from the mixing layer. This suits the temporal scale of OpenLISEM simulations. The mass transfer coefficient describes the rate on which dissolved pesticides are taken up from the mixing layer into the runoff (Havis et al., 1992; Joyce et al., 2008). The exact physical and chemical processes in the mixing layer are highly complex and the mixing layer is simplified as a finite, steady and completely mixed reactor (Havis et al., 1992; Shao et al., 2021). The thickness of the mixing layer is defined as the depth of interaction between the soil and turbulent overland flow (Havis et al., 1992).

Sorption-desorption processes in the soil generally take several minutes to hours to reach equilibrium (Felsot and Dahm, 1979; Pignatello, 2015). At the start of a rainfall event sorption equilibrium can be assumed. With the soil-water partition coefficient (Kd) a linear, instantaneous, equilibrium concentration between dissolved and sorbed pesticides in the soil matrix can be described. This is valid if the sorption process is fast compared to the uptake of pesticides from the mixing layer (Havis et al., 1992). Because simulations for rainfall-runoff events are generally not longer than several hours, degradation of pesticides is neglected.

In addition to dissolved uptake by the runoff, sorbed pesticides can be entrained with eroding sediments. This can be either through splash detachment by raindrops (De Roo et al., 1996), or flow detachment based on the transport capacity of the runoff water (Govers et al., 1990). When deposition occurs, sorbed pesticides are added to the mixing layer together with the sediment. The concentration of sorbed pesticides in runoff sediment can be higher than the concentration in the original soil. This enrichment occurs due to preferred transport of smaller soil particles and organic matter, which also tend to have a higher sorptivity for pesticides (Ghadiri and Rose, 1993, 1991). In a study comparing the enrichment rates for different soils and rainfall events, a decreasing exponential relation with the erosion rate was found (Menzel, 1980).

When the pesticides are taken up by the runoff, either dissolved in water (dissolved phase, DP) or sorbed to entrained soil particles (particulate phase, PP), lateral flow will transport the pesticides further downstream. Combining the soil-runoff interactions and lateral transport of pesticides, at each raster cell in time and space, four concentrations in respective control volumes are relevant for the transport of pesticides:  

- The dissolved concentration in the runoff $(C_{rw},\ mg\ m^{−3})$  
- The dissolved concentration in the soil water in the mixing layer $(C_{mw},\ mg\ m^{−3})$  
- The sorbed concentration in the soil matrix in the mixing layer $(C_{ms},\ mg\ kg^{−1})$  
- The sorbed concentration in the suspended sediment in the runoff $(C_{rs},\ mg\ kg^{−1})$  

These concentrations are influenced by transfer between the runoff and the mixing layer and lateral transport. The water related pesticide fluxes are (1) infiltration through the mixing layer to deeper soil layers, (2) dissolved phase overland runoff and (3) uptake of pesticides by diffusion, convection and turbulent mixing from the mixing layer into the runoff (Joyce et al., 2008). For sediments, the pesticide fluxes include (4) enriched uptake through splash or flow detachment into the runoff, (5) deposition onto the soil surface with deposited sediment and (6) suspended sediment flow which transports the pesticides sorbed to sediment particles downstream. Finally (7) within the soil matrix equilibrium sorption redistributes pesticides between the dissolved and particulate phase. Contrary to the conceptualization of Havis et al. (1992) we do not include dissolved pesticide concentrations in the precipitation but assume that this is negligible.

[[/images/conceptual_model_pesticides.png|height = 500px]] 

This conceptual model can be described with equations for the governing processes. The flux of pesticides in runoff is described as:  

$$\frac{1}{\Delta x} \left( \frac{\partial Q_{rw}}{\partial x}+\frac{\partial C_{rw}A}{\partial t} \right)=k_{film}\cdot (C_{mw}-C_{rw} )-q_{inf}\cdot C_{rw}$$

where $\Delta x$ is the cell size $(m)$, $C_{rw}$ the dissolved concentration of pesticides in the runoff $(mg\ m^{-3})$, $Q_{rw}$ the discharge of dissolved pesticides $(mg\ sec^{-1})$, $A$ the cross-sectional area of the flow $(m^2)$, $x$ is distance in the direction of runoff $(m)$, $t$ is time $(sec)$, $k_{film}$ the transfer rate of the mixing layer $(m\ {sec}^{-1})$, $C_{mw}$ the dissolved concentration of pesticides in the mixing layer $(mg\ m^{-3})$, $q_{inf}$ the infiltration rate $(m\ {sec}^{-1})$.  

The flux of dissolved pesticides in the mixing layer is calculated as:

$$n \cdot z_m\frac{\partial C_{mw}}{\partial t}=k_{film}\cdot\left(C_{mw}-C_{rw}\right)-q_{inf}\cdot\left(C_{mw}-C_{rw}\right)-z_m \cdot S$$


where $n$ is the soil porosity $(m^3\ m^{-3})$, $z_m$ the depth of the mixing layer $(m)$ and $S$ is a sorption source or sink $(mg\ m^{-3} sec^{-1})$ either based on equilibrium sorption (Havis et al, 1992) or kinetic sorption (Joyce 2008). The equilibrium sorption source or sink is calculated as:

$$S = \left( C_{mw} - \frac{C_{mw} \cdot n + C_{ms}(1 - n) \rho_s}{n + k_d (1 - n) \rho_s} \right) \cdot \frac{n}{\Delta t}$$

with $\Delta t$ the discrete timestep in the model, $k_d$ the soil-water partitioning coefficient $(ml\ g^{-1})$, $\rho_s$ is the soil particle density $(g\ cm^{-3})$ and, $C_{ms}$ the sorbed concentration of pesticides in the mixing layer $(mg\ m^{-3})$. The flux of sorbed pesticides in the mixing layer is described with:

$$ \frac{\partial C_{ms}}{\partial t} = \frac{S}{\rho_s} - \frac{S_p \cdot C_{ms} \cdot \varepsilon + S_f \cdot C_{ms} \cdot \varepsilon + S_d \cdot C_{rs}}{z_m \cdot \rho_b}$$

where $S_p$ is the splash detachment rate $(kg\ sec^{-1}\ m^{-2})$, $S_f$ is the flow detachment rate $(kg\ sec^{-1}\ m^{-2})$, $S_d$ the deposition rate $(kg\ sec^{-1}\ m^{-2})$, $\rho_b$ the soil bulk density $(kg\ m^{-3})$ and $\varepsilon$ the enrichment ratio of the sorbed pesticide in the runoff compared to the mixing layer, calculated with:

$$\varepsilon = \alpha \cdot S_j^{\beta}$$

where $\alpha$ is a coefficient for which 7.4 is proposed as a general value by Menzel (1980), $S_j$ the specific detachment rate of either splash or flow detachment $(kg\ m^{-2}\ sec^{-1}) and $\beta$ the exponent, which was found to be -0.2 for many soil and management types (Menzel 1980). Finally, the sorbed flux in suspended sediment is calculated as:

$$\frac{1}{\Delta x} \left( \frac{\partial Q_{rs}}{\partial x} + Sc \frac{\partial C_{rs} A}{\partial t} \right) = S_p \cdot C_{ms} \cdot \varepsilon + S_f \cdot C_{ms} \cdot \varepsilon + S_d \cdot C_{rs}$$

with $Q_{rs}$ the discharge of sorbed pesticides $mg\ sec^{-1}$ and $S_c$ the suspended sediment concentration in the runoff $(kg\ m^{-3})$.
