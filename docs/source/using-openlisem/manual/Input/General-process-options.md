
<p align="center">
  <img width="586" height="374" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/input_general_options.png">
</p>

## Overland flow

1. 1D uses a kinematic wave on a predefined continues network that connects all cells to the outlet of the area (ldd.map in PCRaster). The network is user defined and is an input, the DEM is not used. The kinematic wave uses a finite difference iteration. This option gives only runoff and no flood.  
1. 1D and 2D: this options uses a 1D kinematic wave for the runoff and a 2D dynamic wave for the flooding. The flooding **ONLY** takes place when there are channels defined and they overflow. Flood in this case is overflowing channel water.
1. 2D dynamic flow uses the DEM to determine where water flows. The water is distributed over the cells downstream, using a full dynamic wave with depth-average velocity. The [flow]() tab gives several numerical options. It uses a semi-explicit finite volume solution with a small adaptive timestep (non-iterative).

The distinction between runoff and flood is user defined by a flood threshold (in m). All water above this threshold is considered flood (hazardous), below is considered runoff.  
**NOTE:** channel flow is always 1D kinematic wave over a channel network (lddchan.map in PCRaster).  
**NOTE:** a kinematic wave is advised for large scale problems and less than optimal DEMs (e.g. SRTM or ASTER DEMs).  

## Erosion processes

Switch erosion processes on or off. If switched on, sediment dynamics are simulated for all flows. Default suspended matter is simulated, transport equations can be chosen. Optionally bedload can be simulated. These options can be set for overland flow and channel flow separately, in the [erosion tab]().

## Channels & Rivers
The channel system can be switched off entirely in which case only overland flow is considered. Channels are 1D networks using a kinematic wave for flow.

* Channel infiltration: assumed to be saturated infiltration, the flux equals Channel Ksat $(mm\ h^{-1})$. If not checked the channel is assumed impermeable. Mutually exclusive with baseflow.
* Channel baseflow: the baseflow discharge at the outlet(s) is given in $m3\ s^{-1}$ in the map baseflow.map. This discharge is iterated over the entire channel, so that the kinematic wave gives the stated baseflow at the outlet. Mutually exclusive with channel infiltration.
* Channel Culverts: these are cells in the channel that constrain the discharge a maximum flow in $m3\ s^{-1}$ (specified in ChanMaxQ.map). It is a simple bottleneck; no pipe flow physics are used.

## Infrastructure

Switch on or off the effect of buildings and infrastructure:  
* Buildings have interception storage, are impermeable, obstruct flow and are non-erodible. They can have storage drums for rainwater (interception from roof goes into the drum). Buildings can be smaller than the gridcell in which case the characteristics of the soil and vegetation are used.
* Roads have no interception, are impermeable, smooth, have no detachment of sediment but can have deposition. Roads can be smaller than the gridcell in which case the characteristics of the soil and vegetation are used.
* Hard Surfaces (parking lots, courtyards, airstrips), have no interception, are impermeable, smooth, have no detachment of sediment but can have deposition. They behave the same as roads but are a separate class because they are often in a different land use class.
* Storm drains (tile drains) are networks of circular subsurface pipes with regular inlets. A kinematic wave transports water along this system. The maximum flow is determined by the diameter, gradient and Manning's n. Can be used for urban drains or agricultural tile drains.
