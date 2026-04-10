Infrastructure within openLISEM usually consists of three different types: Buildings, Roads and Bridges or culverts. Often, polygon maps of buildings are available. With a vector gis, building cover can be calculated by rasterizing the map with a very fine resolution, and then resampling this to model resolution.  
Roads usually are part of the land use map. Depending on the modeled resolution in openLISEM, roads should however be provided separately. If the cell size is larger than the road width, a cell containing a road should have the land use next to this road as a dominant land use type, and the road should be set in the roadwidth map. This way, the dominant land use type can be used for infiltration and runoff on the areas within the cell that do not contain road.

**Culverts or bridges**

To find any culverts or bridges, the overlap of roads and channels can be used. 

```
pcrcalc culvert.map = if(roads.map gt 0 and scalar(lddchan.map gt 0) scalar(1), 0);
```

To apply culverts to the channel, either the channel width can be changed, or the channelmaxq can be set.

