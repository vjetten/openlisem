## Vision

Land surface processes influence human populations throughout the globe in a variety of ways. Every-day processes such as infiltration, runoff and erosion can have enormous impact on human lives. The quality of drinking water in villages depends on the transport of pesticides in nearby farmland. Sediment budgets in large catchment influence biodiversity in rivers and lakes downstream. Erosion in mountainous regions influences food security for cities. On the other hand, extreme weather events frequently drive hazardous land surface processes such as (flash) flooding. Here, human lives are directly at danger due to exposure of populated areas to inundation. All of these processes are part of a complex whole that is
the dynamic earth surface system.

Humanity is inescapably bound to the struggle with the dynamic earth surface system on planet earth. We cannot escape the fact that the system is dynamic, and as long as we live on the earth’s surface, we must adapt to how it changes. Furthermore, we must realize and fight the dangerous influence we have in shaping the earth. To increase our capability of taking these challenges, tools are required to **understand, analyze and predict** the behavior of land surface processes, both small and extreme. Understanding can help to stop dangerous processes, such as land degradation. Analysis can inform us about our role in influencing land surface process, and teach us how to better use our resources and live
in a dynamic system. Predictive power can allow us to avoid exposure to hazardous land surface processes.

A fundamental way of development has great benefits in the predictive power of a model. When processes are described at their visible level, they will always behave as expected, limiting insight. When, instead, processes are not described directly, but instead described by the underlying
fundamental processes, the model is able to predict wider ranges of possible behavior. Furthermore, the model might show apparent behavior which was unexpected from the typical observations. Finally, we envision this model to be used by all and developed by all, allowing for greater societal
benefit. This vision influences the way in which OpenLISEM is available for use. The development and knowledge behind OpenLISEM is publicly available and Open-Source. Furthermore, usage of the model is free.

## Maintainers:

- Victor Jetten (@vjetten) : main author, editor and code maintenance.
- Bastian van den Bout (@bastianvandenbout) : numerical 2D flow.
- Meindert Commelin (@mcommelin) : pesticide transport, documentation and code maintenance.  

## Development
The openLISEM model is under continuous development. The model has been developed in c++ combined with the qt UI libraries and the PCRaster Python libraries. The source code is openly available under GNU license, and development is a public process. While the model has undergone extensive testing, it is a large project and bugs and small errors are likely to be present in the code. Any of these can be reported online, or a solution can be presented to the source code directly.

openLISEM version 6.897 (2023/05/22) is created with:

- MSYS2 with MingW64, Qt and CMake (https://www.msys2.org/,http://qt.nokia.com/,https://cmake.org/)
- Qwt technical application widgets for Qt (http://qwt.sf.net)
- Flood source code derived from fullSWOF2D (http://www.univ-orleans.fr/mapmo/soft/FullSWOF/)
- Using openMP for parallel processing (https://www.openmp.org/)
- Using GDAL for map handling (https://gdal.org/)
- PCRaster lib map functions: (http://pcraster.geo.uu.nl/)


## Documentation development

The wiki was initially developed in may 2023 with basic input from available documentation of openLISEM. It is under continuous development to stay up to date with the current development of openLISEM.

On pages that are not finished a 🚧 message is placed to indicate further improvements that are needed. Below a tasks list with the main needed improvements is given:

- [x] setup basic outline of the wiki, and fill the sections with available documentation
- [ ] add section 'Install & Compile'
- [ ] add section 'Manual/Output' based on 'Output interpretation' in documentation14/15
- [ ] add correct referencing in all sections by using with footnotes
- [ ] add section 'Manual/Running OpenLISEM' based on 'Running OpenLISEM' in documentation14/15
- [ ] add section 'Manual/Common problems' based on 'Common problems' in documentation14/15
- [ ] add additional information in 'Theory/Hydrology' based on documentation14/15
- [ ] add additional information in 'Theory/Numerical' based on documentation14/15
- [ ] add links to different section using the 🔣 & 🔨 & 📖 icons
- [ ] add description of continuous mode and groundwater flow in all sections.
- [ ] add a list of studies which use openLISEM, see documentation14/15
- [ ] ...

 

