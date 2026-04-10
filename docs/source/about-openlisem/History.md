## The origin
The OpenLISEM software is a physically-based modelling tool for hydrology, erosion and flooding at various spatial scales (plots to catchments and countries). It was initially developed as LISEM (Limburg Soil Erosion Model), and was created between 1989-1992 for the Limburg Waterboard to simulate detailed soil conservation management in the Netherlands (by Utrecht University, Wageningen University and the University of Amsterdam). The project lead was by Prof Dr Ad de Roo, who also developed Lisflood the basis for the [EFAS](https://www.efas.eu/en/european-flood-awareness-system-efas) model.

LISEM was developed to predict effects of on-site conservation management (conservation tillage, mulching, grass strips and intercropping), and off-site measures. Based on this around 250 storm water buffers were designed to mitigate flooding and erosion problems in the hill area of Limburg. The LISEM model was already integrated in a GIS, it was written in a scripting language that later became [PCRaster](https://pcraster.geo.uu.nl/).

## A brief history

Over the years the LISEM model has been continuously developed based on many EU research projects and PhD researches. It is applied in many studies and projects all over the world. 
LISEM (De Roo et al.)
* 1989 - 1993 – Dutch government: first GIS integrated event-based erosion model, for soil conservation in the province of Limburg (NL).
openLISEM (Jetten et al.)
* 1993 - 2013 – Moved to windows and opensource, processes added: tillage effects, crusting and compaction, pesticides (various EU projects) 
* 2013 – 2016 – Flooding and 2D flow, urban environments: buildings, storm drains, culverts, rain harvesting (World Bank and UNHabitat projects)
* 2016 – 2020 – Sediment dynamics in all 2D flows, parallel processing
* 2020 – 2022 – Continuous modelling, added Evapotranspiration, Groundwater and baseflow, Dams, large scale. 
* 2022 - 2023 - Disolved and particulate pesticide dynamics in 1D flow (by Reindert Commelin) [🔣 ](https://github.com/vjetten/openlisem/wiki/Pesticide-input) [📖 ](https://github.com/vjetten/openlisem/wiki/Pesticides)

LISEM DBASE Generator (Jetten)
* 2022 - now - Development of a Database generator that creates a complete input dbase 

LISEMHAZARD (Van de Bout)
* 2016 - 2023 – multi-hazard model: landslides and debris/mudflows, detailed, event based, disaster oriented (World Bank projects, Italy, China)
