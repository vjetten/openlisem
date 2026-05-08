<p align="center">
  <img width="593" height="408" src="../../../../images/input_calibration_options.png">
</p>

Calibration factors are directly multiplied with the input maps to provide quick and dirty calibration options. A value of 1.0 means that all input maps are used as they are. The value is a multiplication factor, larger than 1 increases the map value, smaller than 1 decreases the map value. More detailed calibration can be done by adjusting the input dataset.

 ❗ The calibration factors have different effects on the simulation: 
 
* Ksat slope and channel: saturated hydraulic conductivity, increase (> 1.0) means more infiltration, less runoff. 
* Theta: initial moisture content, increase (> 1.0) means a wetter soil and less infiltration, more runoff. 
* Psi: soil suction at the wetting front, increase (> 1.0) means a drier soil with more suction, less runoff. 
* Manning N slope and Channel: increase (> 1.0) means more friction, slower runoff, and sometimes more time for infiltration, so can result in less runoff. 
* D50 and D90 (median and 90% quantile of the grainsize distribution). Larger grainsize (>1) will cause faster/more deposition, and will affect transport capacity (depends on the transport capacity equation used).
* Cohesion and Aggregate stability: increase (>1,0) means a stronger soil and less flow detachment (cohesion) or splash detachment (aggregate stability). 

> 📝  Soil cohesion and Channel bed cohesion can be set to -1 in the input maps, in which case there is no flow detachment (but there is transport and deposition). This can be used for concrete channels and surfaces.
