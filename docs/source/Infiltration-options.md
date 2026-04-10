<p align="center">
  <img width="547" height="439" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/input_infiltration_options.png">
</p>

**Infiltration equations**:  
* No infiltration = catchment is impermeable
* SWATRE = 1D finite difference multilayer scheme for soil water balance based on the SWATRE model (Richards equations, multilayer Darcy equation solution), needs complex database setup. Can be used for detailed simulations in agricultural areas.
* [Green and Ampt](https://github.com/vjetten/openlisem/wiki/Infiltration#green--ampt) = 1 or 2 layer infiltration based on a simplified solution for the Darcy equation (no iteration)
* Smith and Parlange = 1 or 2 layer infiltration based on a simplified exponential solution of the Darcy equation

## Impermeable 
The lower soil boundary is impermeable and the soil profile can fill up with water and overflow to create runoff. When impermeable is witched off, there is drainage from the bottom of the soil profile and the soil is infinitely large.  
When SWATRE is used the drainage is the Darcy flux from the bottom node, when G&A or S&P are used, the drainage is the estimated unsaturated conductivity (using the Brooks-Corey equation), assuming a hydraulic gradient $\frac{\delta H}{\delta z}$ of 1.

## 2 layer (Green and Ampt or Smith and Parlange)
Simulates the soil as a two layer infiltration system, all infiltration maps are needed for layer 1 and 2. For example; ksat1.map and ksat2.map, etc. The soil depth is given in mm. Note that the soil depth of layer 1 (soildep1) is the depth of the first layer, but soildep2 is the depth from the soil surface to the bottom of the profile! So soildep2 contains soildep1.

## Compacted areas
Cells which have a fraction of compacted area (compfrac.map) use the value in the maps ksatcomp.map $(mm\ h^{-1})$ and porosity porecomp.map $(-)$ to affect infiltration in the top layer.

## Crusts
Cells which have a fraction of crusted area (crustfrac.map) use the value in the map ksatcrust.map $(mm\ h^{-1})$ and crust porosity porecrust.map $(-)$ to affect infiltration in the top layers.

## SWATRE
* Geometric average: calculate the average ksat between to layers as geometric sqrt(k1*k2), 
else calculate as plain average 0.5(k1+k2)
* Profile table: text file with soil layer definitions, pointing to text files with soil physical properties (profile.inp). 
* Table folder: directory with tables with theta, h and k for different soil types. The tables contain tabular pF curves and tabular K(h) curves. 
