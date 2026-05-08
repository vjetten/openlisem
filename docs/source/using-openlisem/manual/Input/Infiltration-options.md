<p align="center">
  <img width="547" height="439" src="../../../../images/input_infiltration_options.png">
</p>

**Infiltration equations**:  
* No infiltration = catchment is impermeable
* SWATRE = 1D finite difference multilayer scheme for soil water balance based on the SWATRE model (Richards equations, multilayer Darcy equation solution), needs complex database setup. Can be used for detailed simulations in agricultural areas.
* [Green and Ampt](https://github.com/vjetten/openlisem/wiki/Infiltration#green--ampt) = 1 or 2 layer infiltration based on a simplified solution for the Darcy equation (no iteration)
* Smith and Parlange = 1 or 2 layer infiltration based on a simplified exponential solution of the Darcy equation - not well tested.

## General options
* Impermeable lower soil boundary of free drainage   
The lower soil boundary is impermeable and the soil profile can fill up with water and overflow to create runoff. When impermeable is witched off, there is drainage from the bottom of the soil profile and the soil is infinitely large.  
When SWATRE is used the drainage is the Darcy flux from the bottom node, when G&A or S&P are used, the drainage is the estimated unsaturated conductivity (using the Brooks-Corey equation), assuming a hydraulic gradient $\frac{\delta H}{\delta z}$ of 1.

* include tile drain system  
Tiledrains can be specified with maps (see category Stormdrain/tiledrains). Tiledrains can be rectangular (specify width and height maps) or circular (specify a tilediameter map). When a tile depth map is specified (depth in mm where the drain is located) water is extracted from the soil. When using Green and Ampt, water will flow in the tiledrain when the wetting front moves beyond the tile depth (using Ksat). When Swatre is used, water will move into the tile when the matric potential at the tile depth is larger than the user pecified value (default -10 cm).
IMPORTANT: "urban" storm drains are also specified as "tiledrain" maps. Where the tiledepth.map has values larger than 0, the drains are assumed to be in agricultural fields draining the soil water. Where the tiledepth.map has values 0 and the lddtile.map exists, the drain is assumed to be a street stormdrain.

* number of layers Green and Ampt (or Smith and Parlange)
Simulates the soil as a two or three layer infiltration system, all infiltration maps are needed for layer 1 and 2. For example; ksat1.map and ksat2.map, etc. The soil depth is given in mm. Note that the soil depth of layer 1 (soildep1) is the depth of the first layer, soildep2 is the depth from the surface to the bottom of the 2nd layer (soildep2 contains soildep1). Same for 3rd laer. Soildepth is in mm.

* Compacted areas
Cells which have a fraction of compacted area (compfrac.map) use the value in the maps ksatcomp.map $(mm\ h^{-1})$ and porosity porecomp.map $(-)$ to affect infiltration in the top layer.

* Crusts
Cells which have a fraction of crusted area (crustfrac.map) use the value in the map ksatcrust.map $(mm\ h^{-1})$ and crust porosity porecrust.map $(-)$ to affect infiltration in the top layers.

## SWATRE options  

* Save h and theta for all layers (map extention is the layer number)
* Replace the input inithead maps with the matrix potential given for all layers, all cells (negative value)
* Stop calculating Swatre for cells with no runoff and no rainfall. 
  This is faster because all dry cells are skipped, but the matrix potential of the layers no longer updates. Use this for short event based runs and when not using ET. Do not use this if soil hydrology between events is important.
* Minimum internal timestep for SWATRE, use no more than 1/3 of the LISEM timestep. It is always limited to half the LISEM timestep (def. 2 sec).
* SWATRE precison for iteration timestep estimation (higher is more precise, def. 12).
* Folder name and map name of the SWATRE tables and profile definition file (profile.inp). 
