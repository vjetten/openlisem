From the digital elevation model, several topography related maps can be created. 

## The LDD map
A first necessary step is the creation of a local drainage direction map. This map links cells by providing each cell with a local
drainage direction. The value 5 is reserved for pits (catchment outlets).

<p align="center">
  <img width="229" height="224" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/prepare_topo_ldd.png">
</p>

The creation of a Local drainage direction map can be done with the following code:

> `ldd.map = lddcreate(dem*mask, 1e20,1e20,1e20,1e20)`

For more information see the [`lddcreate()`](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/op_lddcreate.html) documentation.

By choosing an arbitrary high value (such as 1e20) for the 2nd trough 5th argument, the entire DEM will be connected into one catchment. By decreasing for example the “maximum local depression outflow depth”, local depressions above a certain depth can be kept. Important to note is that the pit removal is only temporarily applied in order to create the LDD map. The function `lddcreatedem()` takes identical arguments, but applies the pit removal to the elevation map. The output of both the `lddcreate()` and `lddcreatedem()` functions is furthermore influenced by a global
variable. When `-–lddin` is specified, only a single catchment may end at the DEM boundary. When `--lddout` is specified, the subcatchment process for catchments at the boundary of the DEM is done identical to the process for normal subcatchments. For [channels a separate LDD is needed]()

## Gradient

## Outlet & Outpoints
