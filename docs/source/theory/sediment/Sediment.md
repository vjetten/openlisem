# Contents
  * [Splash detachment](https://github.com/vjetten/openlisem/wiki/Splash-detachment)
  * [Flow detachment & deposition](https://github.com/vjetten/openlisem/wiki/Detachment-&-Deposition)
  * [Overland transport](https://github.com/vjetten/openlisem/wiki/Overland-sediment-transport-capacity)
  * [Channels and flooding](https://github.com/vjetten/openlisem/wiki/Channel-&-Flooded-sediment-transport)
      * [Diffusion](https://github.com/vjetten/openlisem/wiki/Diffusion)  
      * [Two Layer sediment transport](https://github.com/vjetten/openlisem/wiki/Two-Layer-sediment-transport)
  * [Settling Velocity](https://github.com/vjetten/openlisem/wiki/Settling-Velocity) 
  * [Sediment Transport Equations](https://github.com/vjetten/openlisem/wiki/Sediment-transport) 

***

For sediment transport and erosion there are three main processes: detachment, deposition and transport. Sediment detachment can occur by [splash detachment](https://github.com/vjetten/openlisem/wiki/Splash-detachment) by raindrops or throughfall. Besides that sediments can be [detached by water flow, or deposited](https://github.com/vjetten/openlisem/wiki/Detachment-&-Deposition) if the transport capacity is lower than the current suspended sediment concentration. Flow related sediment transport is divided into two sections in OpenLISEM: [overland flow processes](https://github.com/vjetten/openlisem/wiki/Overland-sediment-transport-capacity) and [channel and flooding](https://github.com/vjetten/openlisem/wiki/Channel-&-Flooded-sediment-transport). 

While the overland (Gover and Hairsine & Rose) transport capacity was derived from data for overland flow, the channels within OpenLISEM can also use these equations. For larger catchments, the channels do however no longer classify as overland flow, and the empirical relation does not describe the channel sediment transport capacity well. [Diffusion](https://github.com/vjetten/openlisem/wiki/Diffusion) of sediment in overland flow is furthermore assumed to be insignificant because of the small water height compared to cell sizes. Within channels or flooded areas, the water height can however be high enough to allow for substantial diffusion of sediment. In order to overcome these problems, a [2 layer sediment transport model](https://github.com/vjetten/openlisem/wiki/Two-Layer-sediment-transport) was added for the channels and flood water.  
