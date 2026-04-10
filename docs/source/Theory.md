The LISEM model is based on both work by the original developers, and work by many other authors, as will be indicated throughout this section. Equations and theories come from a variety of research fields due to the extent of the LISEM model. The complexity of the hydrological cycle on a full catchment scale, and the complexity of the sediment cycle, both require many processes to be modelled.

Throughout the development of the LISEM model, some processes within both the hydrological and sediment cycle were found to have insignificant research focusing on them. This lead to the usage of equations which were originally not intended for use in that particular way. These equations were, however, found to be the best known option after testing possible varieties. For processes for which many theories are known, the developers either chose to use the theory that best fitted the general use of the model, or provided multiple options within the model.

LISEM is a discrete numerical model, which requires the subdivision of both space and time into a discrete set of locations. LISEM divides the simulation into timesteps of length $\Delta t$ and square cells with width $C_{xy}$. These values are constant throughout the simulation. Cell locations are indicated by a subscript _i_ and _j_, with _i_ the row and _j_ the column within the modelled raster. Time is notated with the superscript _t_.

$$X_{i,j}^t$$

with $X_{i,j}^t$ the value of a row _i_, column _j_ and time _t_.

However, in the theoretical section, equations will in general be given in their continuous form, without indications of discrete properties. Numerical methods are presented in more detail in a separate section of this document.
