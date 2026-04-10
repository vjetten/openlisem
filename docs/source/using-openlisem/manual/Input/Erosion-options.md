<p align="center">
  <img width="633" height="661" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/input_erosion_options.png">
</p>

The erosion options are divided in splash detachment and flow detachment and transport, and conservation measures.

## Splash detachment [📖 ](https://github.com/vjetten/openlisem/wiki/Splash-detachment)

**Rainfall kinetic energy**
The equations describe non-linear relations between rainfall intensity I $(mm\ h^{-1})$ and kinetic energy $(J\ mm^{-1})$:  
* A negative exponential equation with average global parameters (van Dijk 2002).
* A USLE type equation with parameters as found in the EUROSEM model (Morgan et al., 1998). 
* A decreasing power law (Sanchez, 2012).

**Splash delivery** is the fraction of sediment splashed from dry parts in a cell to ponded areas in a cell. The division of dry and ponded parts is calculated from the surface roughness and water height. 

**Splash equation** : In openLisem aggregate stability is the median drop size to halve an aggregate (Lowe test): higher stability means less detachment. In Eurosem Aggregate stability is the opposite: it is the sediment delivery in $(g J^{-1})$. Higher values mean more detachment!

## Flow detachment and transport

For settling velocity [📖 ](https://github.com/vjetten/openlisem/wiki/Settling-Velocity) the default choice is a combination of Stokes and Zanke, the alternative is an equation by Zhiyao et al 2008.

**Detachment efficiency**  [📖 ](https://github.com/vjetten/openlisem/wiki/Detachment-&-Deposition)  
The equation used to calculate the detachment by water flow can be set for both overland flow and channel flow. Detachment efficiency is a value from 0 - 1 multiplied with the detachment rate (depends on soil cohesion in kPa). 0 gives no flow detachment, 1 gives max flow detachment. openLISEM is very sensitive to the efficiency. Generally, equation 1 gives the most detachment and equation 2 the lowest detachment.

> 📝 a negative cohesion value in `coh.map` or `chancoh.map` gives zero detachment (useful for concrete channels for instance), but deposition and transport are allowed.

**Suspended sediment transport**  [📖 - overland](https://github.com/vjetten/openlisem/wiki/Overland-sediment-transport-capacity) and [📖 - channels](https://github.com/vjetten/openlisem/wiki/Channel-&-Flooded-sediment-transport)  
For overland flow the default is Govers, which after tests shows to be valid for a wide range of textures and flow regimes. In channel flow the water depth is larger and Van Rijn's equation may be more suitable.

**Include bedload transport - channels**
Sediment is simulated in two layers: suspended and bedload. Transport capacities can be different for each layer, see below. The water depth is divided into a bedload layer depth and a suspended layer depth. Suspended transport is based on the $D_{50}$ (median) value of the grainsize distribution, bedload on the upper D90 (upper 90% grainsize). Detachment and deposition are reported as the sum of the two.

When the water hight in channels increases, diffusion [📖 ](https://github.com/vjetten/openlisem/wiki/Diffusion) can have a substantial impact on suspended sediment transport, this can be included.
