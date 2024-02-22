openLISEM
============
Date: 231216
============

This software is subject to a DISCLAIMER and released under the copyright model GPLv3

For questions contact v.g.jetten AD utwente.nl

NOTE: only a 64bit version exists, 32 bit is not supported
NOTE: The code since version 5.6 is compilable under linux (checked for Ubuntu)
NOTE: since version 6.x it is fully parallel and developed with MSYS2.0, Qt5, openmp, gdal and pcraster

version 6.94
- Improved subsurface storm drains, can be used in urban envrioenments or tile drains in fields, BETA

version 6.90-6.93
- NOT YET available: Created a new Richard Equation soil water system, multilayered
- imporved some ET calculations, may still cause mass balance errors
- small changed to boundary conditions SWOF
- fixed a bug in Rusanov Riemann solver
- Improved domain boundary flow detection when open boundary

version 6.898
- 3 different ways of groundwater flow, for different scales
- experimental: added possibility to have inflow into channel from saturated parts of the soil (in advanced options)

version 6.893 - 6.897
- KNOWN bug: activating evapotranspiration gives a mass balance error, not all evaporation losses are accounted for
- Groundwater connection to the channel is still experimental
changes and added options:
- suction below the wetting front can be supplied with maps psi1 and psi2, but also intenally calculated using pedotransfer functions (Rawls 1982)
- three ways of groundwaterflow are added for different scales
- it is now possible to add discharge manually (for instance dam spill) at userdefined points in the channel
- connection of overlandflow to the channel can be made iterative (in advanced options) for large timesteps
- new libraries, especially lbgdal-32.dll, do not mix with old dlls

version 6.891 - 6.892
- CHANGED calcyulation of unsaturated Ksat for percolation and redistribution! Closer to Brooks-Corey, Saxton and Rawls 2006
- Fixed bug where matrix potential PSI was multiplied by 0.01 twice (cm to m). Infiltration will behave differently now
- Fixed bug causing profile to become saturated instantaneously because percolation from layer 1 to layer 2 saturates subsoil
- fixed bug open and closed boundary flow, to be further tested on different catchments
- fixed small bugs for reporting to files all output to files with clear separation of total outflow and channelhydrographs fopr water and sediment
  and reporting of units m3/s and l/s

version 6.883 - 6.89
- redesigned GW flow completely, only two calib. paraters now, flow to the GW layer and GWflow to the channel
- Interface works better with low resolutions screens
- erosion overland flow uses govers and harsine and rose
- included iterative version of connection between 2D flow and channel
- fixed bug in sediment in non-iterative connection of 2D flow and channel
- fixed screen display of discharge, was cumulative

version 6.883-6.885
- fixed screen display of discharge factor 1000!

version 6.882-6.883
- fixed kin wave MB error (interception of houses)
- fixed erosion MB error (caused by cell_depositInfil(r,c))
- checked all erosion functions in Splash, 1D and 2D flow for consistency and logic
- kniown bug 1D2D flow (overflow channels) has mass balance errors in water and sediment, needs to be fixed or minimized

version 6.881
- fixed major interface bug: factor 1000 in screen display of runoff! Factor is already done in Qoutput itself, does not need to be done in interface

version 6.88
- Qmax in culverts better, Qmax does not have to appear in Kin Wave
- merged GW, pressure based and pref flow as in SWAT
- Ksat crust exponential decline from kssat1 to ksatcrust, porosity untouched for problems waterbalance(?)
- better behaviour culverts
- buildings can be added to the dem form a fraction onward (def 0.3)
- show runfile name on screen en save runfile to result dir, save all screens
- known bug: rainfall of multiple stations not working well!

version 6.873
- experimental: added 2D GW flow
- KNOWN BUG: when using ET mass balance is not correct (calculations are)
- changed behaviour of flow to culverts, gradual increase of mannings n until Q < Qmax
- fixed bug in display of avg soil moisture in Map screen
- interface: fixed bug when changing river size in map screen
- interface: fixed bug in tree display of maps

version 6.872
- fixed interface bugs that crashed LISEM
- fixed zize of output dot
- fixed map display aspect ratio

version 6.87
- Fixed some problems with culvert flow (chanmaxq.map)
- Better display of culverts in the map (white on black)
- show hard surfaces in the map

version 6.86 beta
- More rigorous checking of rainfall station ID versus map ID
- Inverse distance interpolation of rainfall, fixed bug and linked to IDgauges.map
- code: rewrite of application.cmake
- new dlll libs with the lates version of MSYS (13 aug 2022)

version 6.84.9-6.85 beta
- enabled inverse distance interpolation from station rainfall
- map draw now with real ccoordinates
- fixed boundary outflow when there is no channel
- fixed small bugs in percolation nand redistribution of water
- separate calibration for ksat2, to have better control over subsoil and GW
- experimental: channel tortuosity in calibration factors to account for non-rectangular channel

version 6.84.6-6.84.8 beta
- fixed a mass balance bug when adding water to roads
- fixed a mass balance bug when activating ET

version 6.80-6.84.6 beta
- added Engelund and Hansen for river sediment transport (after Hecras)
- added direct transport factor for river transport (skipping bed erosion)
- fixed bugs in channel concentration output to file
- fixed bug when first outlet nr > 1
- experimental dynamic diagonal flow in advanced options
- fixed crash when infiltration is set to none

version 6.77-6.80 beta
- re-evaluated Splash equations and added Eurosem method. In Eurosem the aggr stab is in fact the splash delivery in g/J
- Fixed baseflow mass balance as far as possible, Mass balance for ETa might still be off
- Fixed QSall in output which was in kg and not kg/s
- Fixed all consistency problems between output of outlets and total for water and sediment. Outlets give only outlet values, 
  Qall and Qsall (hydrograph 0) give the sum of channel outlet, overland flow, boundary flow and storm drains

version 6.7-6.77 beta
- corrected ETa, added to screen
- Thetai1 and 2 in display and reported are now average of the soil layers
- baseflow according to SWAT added with stationary baseflow
- mass balance shows error bercause water from ETa and baseflow is not from rainfall 
- added check on river cross section: when the width is > cell size the depth and widt are adjusted so that the hydraulic radius is maintained 
- added a maximum timestep of 60 sec for the river kinematic wave

version 6.69beta
- fixed output timeseries
- blocked output of hydrograoh values for now, memory leak suspected
- fixed interface errors: reset values for option tabs 

version 6.68beta
- known bug: output timeseries not working byb accident

version 6.67beta
- started continuous modelling
- enable reading of long timeseries of gauges or maps with day:minute indication
- added ET calculations => needs testing
- bug fixes, there is still a small bug in the combination kinematic wave and 2 layer Green and Ampt
- updated graph drawing internally elimination redundant datastructures
- updated the interface

version 6.61beta
- updated libraries
- fixed bug in infiltration causing sometimes to ignore impermeable option

version 6.5-6.6
- SWOF sediment changed fully parallel computing
- maximize parallel comnputing efficiency in all processes
- Diagonal (LDD) flow with SWOF when pit in X or Y dircetion using user defined pit threshold
- icon changes
- NOTE Swatre has not been checked for a while

version 6.2-6.4
- bug fixing in sediment
- corrections in SWOF 2.0 so that it behaves more as the old SWOF, better now
- interface changes
- checked infiltration and percolation
- further openMP related optimizations

version 6.1 BETA (warning: new dlls, do NOT mix with pre 6.0 versions)
- fixed bug in getting values from Riemann solver. This was solved before but reappeared!
- SAFEST choice for flood modelling is the SWOF without or with MUSCL. SWOF 2.0 is experimental 

version 6.0 BETA (warning: new dlls, do NOT mix with older versions)
- extensive rewriting of the code to use parallel processing with openMP
- changed compilation to MSYS so the newest versions of QT and MINGW are used
- new 2D flow process (very fast), still being tested 

version 5.97-5.98 beta
- Bug-fix: file runoff.map did not show flood when choosing kinematic+dynamic flow
- flow to channel with kin V using perpendicular angle and from channel with gravity
- spurious velocities as a result of SWOF limited to sum of kinematic velocity+pressure velocity
- Advanced options for SWOF stuff
- Outlet shows numer when hovering with the mouse over it
- better dealng with highDPI screens

version 5.96
- Bug fix: flood height total in mm not reported on screen
- Bug fix: did not show all outlets as dots in map display
- Bug Fix: Initial food height (whinit.map) was not activated
- Interface more consistent on 1D, 1D+ and 2D flow, flood threshold only for 2D
- Corrected help files for global options

version 5.91-5.95
- fixed: bug in 2D dyn ewave seidment: flow detachment was not working
- fixed: impermeable lower boundary was always off
- extended options and checked working of sediment traps and grass trips
- 3 flow options: 1) kin wave (no flooding); 2) kin wave and flooding from channels; 3) full dynamic wave
- river vector display can be toggled and line size adjusted

version 5.7-5.91
- fixed buggy 1d kin wave and 2d flood for sediment
- revived option kin wave only
- cosmetic and code cleaning

version 5.5-5.7
- Fixed all sediment MB errors, both simple and two phase flow (suspended and bedload) for 2D flow
- Channel displayed as vector and culverts as white dots
- colours close to QGIS palettes
Still buggy: 1D kinematic coupled to 2D flow for both water and sediment, and sediment classes

version 5.5
- Added subsurface stormdrains: kinematic wave on a subsurface LDD network (lddtile.map) with inlets,  diameter (pipe), mannings n, gradient
- Cleaned up everything so that it compiles under QT 5.13, QWT 6.1.4, mingw 73, revised LISEM_EXTERNAL (no fern)

version 5.4
- fixed a bug in SWATRE where the infiltration form SWATRE was not properly subtracted from the surface water layer
- Fixed a bug in Green and Ampt using Psi in the wrong units (cm vs. m)
- made the calibration for psi also count for SWATRE (this value is multiplied with the inithead maps)
- fixed a small bug where the screenshots were not in the right directory

version 5.31
- fixed a bug in SWATRE where the infiltration form SWATRE was not properly subtracted from the surface water layer
- made the calibration for psi also count for SWATRE (this value is multiplied with the inithead maps)
- fixed a small bug where the screenshots were not in the right directory

version 5.3
- Fixed a bug in percolation for a 2 layer Green and Ampt infiltration, where the percolation flux was not converted from mm/h to m/s
- Fixed a bug that prevented the kinematic wave from causing floods
- Added theta1 and theta2 to the output screens when using G&A

version 4.94-5.2
- fixed a bug in SWATRE, water level infiltration was not adjusted after infiltration
- courant values between dynamic and diffusive wave mixed up in run file
- fixed code on reporting cumulative erosion values (may need more work)
- added momentum and soil moisture to screen map display
- merged output flooding and water height parameters, they were the same anyway
- difference between runoff and flood is user choice on water level
- Interface: all ',' in the run file replaced with ';' and decimals are always read as '.'
- Interface: reset button per tab for default values and better help information
- Interface: added variables to map display
- Interface: Geotiff file can be used as background

NOTE: mass balance errors in erosion modules still possible

version 4.92-4.93
- fixed a bug in boudary conditions dynamic wave sometimes causing too much vertical flow
- fixed a coding error in erosion slowing down calculations
- some interface cleanup (Heun in dynamic wave is not useful)
- multicore remains more unstable than single core, and less accurate
- Possible bug: flood deposition can become negative and positive (only negative allowed)

version 4.91 (180119)
- fixed a bug in the screen output causing negative values of discharge and a large mass balance error. 
  Note: when you select dynamic wave flow, the flood height in mm reported is not part of the mass balance, as it is already 
  included in overland flow.

version 3.99-4.9 (180116)
- added dynamic wave for overland flow, three numerical solutions now: kinematic (using LDD), diffusive and dynamic (using DEM). 
- Flooding is always solved with a dynamic wave, channel flow is always kinematic.
- Random roughness and surface storage slightly changed to avoid bugs in the dynamic wave solution
- Multi CPU Core application for paralklel computing. If this gives problems, select only 1 core.
- small bugfixes 

version 3.97 - 3.99 (170308)
- BugFix: Calibration factor for Cohesion and Aggregate stability (they were reversed)
- BugFix: in Channel Cohesion, soil sohesion was used instead of channel cohesion
- BugFix: cleaned up litter interception, roof interception and effect of raindrums 
- BugFix: ensures that screen information, file and map output is all the same
- Added: Fixed bilinear interpolation for sediment, other options give mass balance errors
- Added: EXPERIMENTAL: an empirical factor (1-99) to increase the flow in the direction of the steepest resistance slope for diffuse overland flow. 
- Added: calibration factor for channel cohesion
- Added: if the cohesion of slopes or channels is negative, the detachment is assumed to be zero. Deposition will take place
- Added: total interception (roofs, canopy, litter, randrums) to the screen output (in mm)
- Added: the possibility to write GeoTIFF files for the main map outptu. GTiff is not georeferenced). Tif input is also automatically possible (experimental) 
- Added: small interface changes to deal with low resolution screens

version 3.96 (170211)
- Fixed a bug in coupling overland flow and flood water
- fixed a bug in automatic adjustment of legend when switching maps

version 3.95 (170214)
- Fixed a bug in flowboundary conditions
- Fixed a minor bug in 2D dynamic flow
- Added GeoTiff writing for main output maps (not timeseries)
- Added 3 missing fluxes (in mm) to report with totals

Version 2.03-3.94 (170212)
================================================================================================================================================
The following features are inplemented in the course of 2015 and 2016:

- Multiple catchments can be done in one simulation (limited by PC memory)
- Spatial rainfall from rainfall maps (e.g. model results or satellite images)
- 1D kinematic wave for runoff and channels based on a usr defined network
- Diffusive wave 2D runoff, uses the DEM directly (does not use the LDD which is reserved for the Kinematic wave)
- Dynamic wave shallow flooding, based on the FULLSWOF project: https://www.univ-orleans.fr/mapmo/soft/FullSWOF/
- Sediment dynamics (erosion, transport, deposition) for splash, runoff, channel flow and flooding
- Flow barriers (walls) can be defined to close areas in 1-4 directions, e.g. as flood protection walls of a given hieght
- Channels can be wider than a gridcell
- Channel stationary baseflow can be added
- Sediment transport for suspended matter only, suspended matter and bedload movement, choice from different transport equations
- Sediment transport for a median grainsize class or multiple classes (based on D50-D90 or user defined)
- All fluxes are made visible as maps on screen while running, and saved to files
- Hydrographs and sedigraphs from multiple observation points and outlets can be shown onscreen while running, and saved to files

Because of diffusive runoff the range of resolutions can now be larger, LISEM is being tested from 1 cm gridcells (2 m2) to 20m gridcells (600 km2)


1) The openLISEM model 
======================

This is the event based spatial runoff and erosion model openLISEM. Thank you 
for downloading. openLISEM simulates the spatial dynamics of surface runoff and 
erosion for catchments of 1 ha to 500 km2. It is based on the LISEM model that 
is available here: www.itc.nl/lisem. Details about the theory and dataset for 
now can be found on this website (although a bit outdated). openLISEM uses the 
freeware GIS PCRaster (http://pcraster.geo.uu.nl) for database creation and 
analysis of the results. 

2) Terms of use 
===============
This software is free and open source, hosted by  sourceforge.net. The project 
details can be found on:  http://lisem.sourceforge.net

It is distributed under the GPLv3 licence, distributed with this package

Good Luck

Victor Jetten

Chair Natural Hazards and Disaster Risk Management 
Department of Earth Systems Analysis  
Faculty ITC, Twente University, 
the Netherlands 
v.g.jetten AD utwente.nl
