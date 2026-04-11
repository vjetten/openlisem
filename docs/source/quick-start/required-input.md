## Required input

OpenLISEM needs a minimum of 24 raster maps depending on the input options selected in the interface such as special surfaces, sediment modelling, and 2-layer sediment transport. In theory these maps can be derived from 5 base sources:

1. Rainfall
1. Topography
1. Land use/cover
1. Soil type
1. Infrastructure

The number of input maps needed for a OpenLISEM run depends on the selected options. Most of these options can be set from the user interface. Configuration of a OpenLISEM model run can be saved and loaded in a .run file. This is an easy to edit ASCII file containing a line for every option, with the option name and its value. Default values are used to fill up any missing data in a .run file. The options can also be used as command line arguments when OpenLISEM is used without user interface. For more information about this, see the section about running OpenLISEM.

Besides options about functionality, many options about technical parts of the model are available to users. Minimum time steps, courant factors and approximation methods are some of these options. In general, when not sure about the influence of an option on the model, it is best to use the default value for that option. The default values have been chosen in such a way that they should ensure optimal performance in a wide range of problems. Changing these values to outside their proper range can generate errors within the model, and lead to unrealistic behavior.

The raster format causes openLISEM to need a large data set. For instance a river channel is characterized by its network, shape, width, depth, bed slope, resistance, strength and infiltration rate or stationary baseflow (8 maps). Do not be scared by this!

For information about the preparation of all specific input maps and data we refer to the [manual](../using-openlisem/manual/index.md). There we discuss the preparation of the following input data:

1. Rainfall [🔨 ](../using-openlisem/manual/Prepare-data/Prepare-Rainfall.md) [🔣 ](../using-openlisem/manual/Input/Meteorological-input.md)
1. Catchment [🔨](../using-openlisem/manual/Prepare-data/Preparing-Topography.md)
1. Land use [🔨](../using-openlisem/manual/Prepare-data/Prepare-Land-Use.md)
1. Soil properties & Infiltration [🔨](../using-openlisem/manual/Prepare-data/Prepare-Soil-Properties.md)
1. Channels [🔨](../using-openlisem/manual/Prepare-data/Prepare-Channels.md) [🔣 ](../using-openlisem/manual/Input/General-process-options.md#channels--rivers)
1. Infrastructure [🔨](../using-openlisem/manual/Prepare-data/Prepare-Infrastructure.md) [🔣 ](../using-openlisem/manual/Input/General-process-options.md#infrastructure) 