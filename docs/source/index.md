# Welcome to the OpenLISEM wiki

> 🚧 This wiki is an initial version and many sections still [need attention](https://github.com/vjetten/openlisem/wiki/About-OpenLISEM#wiki-development)!

***

The Open Limburg Soil Erosion Model (OpenLISEM) is a physically based numerical model with the purpose of event based or continous runoff, flooding and erosion modelling on a catchment scale. OpenLISEM uses a square grid to solve both cell specific processes, and the differential equations governing flow. 
The main features of OpenLISEM are:
* Rainfall from different sources and interception by vegetation and buildings
* Infiltration, soil hydrology, using 1 to 3 layers or a multilayer soil water balance 
* Groundwater flow
* Overlandflow flow in 1D (kinematic wave) or 2D (dynamic wave)
* Channel flow (1D, kinematic wave)
* Channel flooding (2D dynamic wave)
in all flows
* Erosion (splash and flow detachment)
* Sediment transport and deposition

**IMPORTANT**: there are two different versions of LISEM. One is this model -openLISEM- the other is LISEMHazard maintained by Dr. Bastian van de Bout. LISEMHazard can simulate slope stability, mass movement and debrisflows, which openLISEM cannot. The other main difference is that openLISEM is both event based and continuous, while LISEMHazard is event based only. Main model principles are the same but the code is different so the behaviour may not be the same. LISEMHazard can be found here: [https://lisemmodel.com/](https://lisemmodel.com/)

This wiki is designed as a guide for using OpenLISEM. For a first look see the [general introduction](https://github.com/vjetten/openlisem/wiki/Introduction) and the [quick start](https://github.com/vjetten/openlisem/wiki/getting-started). In the [manual](https://github.com/vjetten/openlisem/wiki/Manual) detailed information is provided to setup the input for a model run. The [theory](https://github.com/vjetten/openlisem/wiki/Theory) section explains the background and concepts used in OpenLISEM. Many topics in OpenLISEM are covered in three different sections, links between these topics are depicted with an icon: Input = [🔣](https://github.com/vjetten/openlisem/wiki/Input) , Data preparation = [🔨](https://github.com/vjetten/openlisem/wiki/Prepare-Data) and Theory = [📖](https://github.com/vjetten/openlisem/wiki/Theory). We aim to provide some tutorials to apply different functionalities of the model. OpenLISEM is an open source modelling project, and the model and wiki are under continuous [construction](https://github.com/vjetten/openlisem/wiki/About-OpenLISEM)!

## OpenLISEM design  
The model aims to be applicable on a variety of scales, it has been used in projects from very detailed (1 ha in 1m gridcells) to relatively large (5000 km<sup>2</sup> in 200m gridcells). It is mostly used for research that requires 5-20 m gridcells and can be used both in rural and in urban environments.

OpenLISEM will run with practically any dataset it is given, which is however not a guarantee for good results. The model is not 'smart' and will not warn you for unrealistic parameter values or unrealistic combination of parameters. This gives the user maximum flexibility, but also assumes that you know what you are doing, and are familiar with the physical meaning of the input maps, as wellas have a basic understanding of the processes involved. 

Emphasis is put on detail: characteristic about the model is the capacity to handle sub-gridcell surface properties (Figure 1). A gridcell can contain a bare soil, crusted/compacted soil, vegetated surface, a road, a building and a channel. These surface characteristics are supplied in separate layers as fractions of the total cell area. The base layer is formed by the soil surface with its hydrological characteristics and the user supplies additional maps that trigger additional hydrological processes in the model. The presence of a vegetation will, for example, result in interception on a part of the gridcell. The presence of a building will result in roof storage and a partly impermeable surface, and a road will have sedimentation but no infiltration or erosion. 

<p align="center">
  <img width="594" height="345" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/surface%20fractions.png">
</p>

Most simulations use 2D overland flow in openLISEM, based on the SaintVenant equations for shallow flooding. The numerical solution is a'semi-implicit finite volume solution' (REF). The code is based on the [FULLSWOF](https://arxiv.org/abs/1204.3210) code (Delestre et al. REF). A Riemann solver is used to create mass and momentum continuity between cells. This means that the solution is not iterative, but the smallest timestep is found in the flow domain to solve the St Venant equations and therefore the timestep varies for the flow module. Because the implementation is fully parallel, the model is relatively fast and profits from processors with multiple cores. The flow uses directly DEM information, including terrain features, obstacles such as buildings or dykes, and depressions such as rainwater buffers and small dams. Because the model is originally an erosion model, erosion and sediment dynamics are available in all flows, which makes openLSIEM fairly unique. 




**DISCLAIMER**  
_No warranties, expressed or implied, are made that the computer programs described in
this wiki are free from errors or are consistent with any particular standard of
programming language, or that they will meet a user's requirement for any particular
application._

