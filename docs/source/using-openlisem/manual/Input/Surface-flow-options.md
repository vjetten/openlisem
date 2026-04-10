<p align="center">
  <img width="547" height="439" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/input_flow_options.png">
</p>

## Flow boundary and barriers (for 2D flow)

0 - outflow not allowed except in outflow points in `outlet.map`. Water can pool at the edges if the DEM is not perfect.  
1 - outflow to all sides possible of flood and overland flow. Outlets always allow outflow. (default)  
2 - a user defined map `flowboundary.map` allows 2D outflow where the value is 1, and a closed boundary where the value = 0. A map filled with 1 is the same as unchecked, free outflow.   

> 📝 a kinematic wave always has closed boundaries except for the outlet.

## Flood threshold
When using 2D dynamic flow there is no difference between runoff and flooding. In order to avoid spurious flood levels reported as "flood" a threshold level (in m, default 0.05) can be set to distinguish between a flood and non-hazardous surface flow or runoff.

## Initial water level
A simulation can be started with an initial water level: supply a map called `whinit.map` which has initial water levels (in m).

## Buffers
`buffers.map` is a map added to the DEM, **works only for 2D flow**. Positive values (in meters) are gridcells added to the DEM (creating dikes), negative values are values that are subtracted from the DEM (creating depressions).   

## Flow Barriers between cells
Barriers up to a given height (m) to disconnect flow between cells **(only for 2D flow)**. The barriers are assumed to be on the cell boundary and have no physical thickness, only a height. Cells can be disconnected in the 4 directions NESW. This needs `flowbarrier.map` with codes 0 (no barrier), and 1, 2, 3 ... n. The numbers correspond to the text file `flowbarriers.txt` where the height of the barriers are given. 

## 2D dynamic wave parameters
The SD Shallow Water Overland Flow (SWOF) equations are a semi-implicite finite volume solution, using progressively smaller timesteps to solve the St. Venant flood equations, using a Riemann Solver and a Flux Limiter. Fluxes and states can be estimated at the cell boundaries or from the cell centers.
The code is originally based on the [FullSWOF source code](https://www.idpoisson.fr/fullswof/) of the University of Orleans.  
* SWOF 2.0: Simplified and very fast SWOF solution, using cell center values and a flux limiter only for the hydraulic pressure differences.
* SWOF Cell centered solution: water pressure and velocities on cell boundaries based on a 1st order Taylor solution, only using information in the cell itself. 
* SWOF MUSCL solution (1): Water pressure and velocities on cell boundaries based on a 2nd order Taylor solution, using information in the cell and its neighbours in X and Y direction. Slower and more accurate. A MUSCL scheme (Monotonic Upstream-centered Scheme for Conservation Laws) estimates fluxes and states at cell boundaries (slow).
* Time average velocity: Velocity of cells is a weighted average between two iteration steps (def on), use for stabler solutions and to avoid spurious velocities (def on)
* Enable diagonal flow when blocked: water can be blocked in X and Y direction in 2D flow (a local depression, while a diagonal outlet exists. Diagonal outflows are not regarded in a 2D scheme. Enabling this uses the LDD to transport water and sediment to this diagonal cell **(default on)**. The threshold pit depth is set at 0.1m (a guide value is 1% is the cell size or less).  

Courant factor: Determines minimum timestep, lower means smaller timesteps in the iteration are used. (default = 0.2, range 0.01-1.0).

Minimum timestep: User defined minimum timestep for dynamic wave solution. In case of instability lower this value (default =  0.2sec).
