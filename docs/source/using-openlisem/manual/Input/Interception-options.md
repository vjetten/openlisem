<p align="center">
  <img width="547" height="439" src="../../../../images/input_interception_options.png">
</p>



* Canopy openness = factor describing the open area between the leaves in the canopy that causes direct throughfall. Default openness is = 0.45. This causes the throughfall fraction for instance with LAI 6 to be = exp(-0.45*LAI) = 0.067
* Activate litter effect on interception and also splash erosion: uses a map with the fraction of litter cover, litter.map (0-1), which has a user defined storage (def 1 mm).

**[Canopy storage equations](../../../theory/hydrology/Interception.md):** links LAI (in $m^2\ m^2$) to storage capacity of the canopy (Smax in mm) for a set of researched vegetation types. This is the storage capacity of the plants, not of the pixel. So a single tree in a pixel can have $2 mm$ of storage, but the average in the pixel depends on the cover fraction of that tree in the pixel. Internally the interception storage of the pixel is calculated.  
**user defined:** in case of varying land use/vegetation types, or a vegetation type not in the equations, the user can provide a canopy storage map, smax.map $(mm)$.

📝 **this section needs to be revised with the use of NDVI timeseries and a smoother way of including interception parameters for complex landuse. Beside that the canopy openness is not well documented - needed?**

