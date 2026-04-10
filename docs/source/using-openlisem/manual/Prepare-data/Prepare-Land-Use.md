Many variables in openLISEM depend on the landuse for their spatial distribution. This data can be obtained from existing maps, field observations or remote sensing. One method to create input maps for these variables, is by setting a value for each land use type present in the catchment. For example for crop height, a value for forest, grasslands and arable land set. By combining these values with a land use map with these categories a map of the crop height in the catchment can be made.

To prepare the land use properties, a land use type unit map and a land use properties table are needed. We combine these with a [`lookup()`](https://github.com/vjetten/openlisem/wiki/Introduction-PCRaster-&-Nutshell#lookup) command.



