If a channel or river network is included in OpenLISEM, many maps with the same variables, but specific for the channels are required, for example ksat, Mannings'n and slope. Besides that also channel specific maps are needed like channel width and side angle.

## Channel LDD
To create a channel LDD, a channel mask is needed. This can be created with three methods: Digitizing, Stream order and flow accumulation.  Digitizing can be done with a vector GIS package. After the channel route is traced, this path can be converted to raster format.  
To create a mask by minimum stream order, use the following code:

> `streamorder.map = streamorder(ldd.map)`  

The resulting stream order map indicates for each cell, how many confluences have, at least, preceded the current cell within all branches. The first cell in a local drainage direction link thus always has a stream order of 1. An example of a streamorder map and a derived channel mask are given in the figure below.

<p align="center">
  <img width="548" height="251" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/prepare_ldd_streamorder.png">
</p>

Another method of creating a channel mask is by using accumulated flux. To create a mask by accumulated flux, use the following code:

> `accuflux.map = accuflux(ldd.map,1)`

This map represents the water volume that passes through a cell, when 1 mm of rainfall falls on every cell, and all water is routed through the local drainage direction. To obtain a channel mask set an appropriate threshold. This threshold should result in a channel map comparable to the
real route. This can be checked with aerial photographs. An example of an accumulated flux map and a derived channel mask are given below.

> `chanmask.map = if(accuflux.map > x, 1)`

<p align="center">
  <img width="518" height="232" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/prepare_ldd_accuflux.png">
</p>

To create a channel local drainage direction, the obtained mask must be in the correct format. All non-channel cells should contain a missing value indicator, and the channel cells should contain a 1. To do this use the following code:

> `chanmask.map = scalar(if(chanmask.map,1.0))`

Then use the lddcreate function on the masked DEM:

> `channelldd.map = lddcreate(dem*chanmask, 1e20,1e20,1e20,1e20)`