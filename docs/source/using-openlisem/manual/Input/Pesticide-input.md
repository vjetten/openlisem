 🚧 **The pesticide extension in OpenLISEM is under initial development, it is not yet fully integrated in the main code. It is not active in the current released version of openLISEM, but can be found in the branch `lisem-pest`. Initial testing and publication are currently in progress. Currently the following features work:**
* **1D kinematic flow calculations of dissolved and particulate pesticide transport**  
* **runfile based input of pesticide choices and options**  

**The following features are planned for 2024:**  
* **intergation with 2D dynamic wave flow routing**  
* **full integration with other OpenLISEM options like channels, buildings etc.**  
* **fully functional user interface**  


***

To run the pesticide module of OpenLISEM-pesticide, the OpenLISEM model must already be setup for runoff and sediment simulation. The model can also run only for water and DP pesticide transport and neglect all sediment and PP pesticide dynamics. 

> 📝 because currently no user interface is available for the pesticide module, all settings must be applied through the runfile!

To activate the the pesticide module in OpenLISEM, the option ‘include pesticides’ must be activated in the runfile. The input for pesticides consist of five maps and four parameters (see table below). Two maps are only required if erosion and PP transport is simulated, this depends on the option ‘[include erosion](https://github.com/vjetten/openlisem/wiki/General-process-options#erosion-processes)’ in OpenLISEM. The five maps are the concentration of the pesticide in (1) DP and (2) PP of the mixing layer, the (3) concentration in PP of the deeper soil. Beside that the depth of (4) the mixing layer and of (5) the deeper soil PP profile. The four parameters are specific for the simulated pesticide; the soil-water partitioning, the mass transfer coefficient of the mixing layer, and the coefficient and exponent for the enrichment equation. Besides that, also the bulk density of the soil in the mixing layer is needed, a default value is 1500 $g\ kg^{-1}$. As default OLP simulates equilibrium sorption in the mixing layer, however also kinetic sorption can be simulated. For this a sorption rate (Kr) is required, which limits the rate of the sorption-desorption process. Besides that, also the bulk density of the soil in the mixing layer is needed, a default value is 1500 $g\ kg^{-1}$.

Name | Description | unit | range / **default value** $^1$ | format | Required
-- | -- | -- | -- | -- | --
pcmixwat | DP concentration in mixing layer | $mg\ L^{-1}$ |   | Raster map | Always
pcmixsoil | PP concentration in mixing layer | $mg\ kg^{-1}$ |   | Raster map | Always
pcsoil1 | PP concentration in deeper soil | mg kg-1 |   | Raster map | For PP
pestmixdep | Mixing layer depth | m | 0.001 – 0.025 | Raster map | Always
pestsoildep1 | Pesticide depth deeper soil | m | $z_m < z_s < soildep\ ^2$ | Raster map | For PP
Kfilm | Mixing layer transfer rate | mm sec-1 | 0.00001 – 0.1 | double | Always
Kd | Partition between DP and PP  | L kg-1 | 1 - 10000   | double | Always
rho | Bulk density of the mixin glayer | g L-1 | 500 - 2650 | double | Always
ERmax | Coefficient for enrichment ratio | - | 1-20; **7.4** | double | Optional
ERbeta | Exponent for enrichment ratio | - | -0.01 - -0.05; **-0.2** | double | Optional
Kr | Sorption rate | min-1 | 0.001 – 1; **-1** $^3$ | double | Optional
Pesticide name | The name of the modeled pesticide which   is used in the output | - | - | character | Always

$^1$ Default values are shown **bold**  
$^2$ _soildep_ = the total soildepth _(m)_ as initialized in OpenLISEM  
$^3$ The default value for Kr (-1) lets the model use equilibrium sorption

## Creating input maps for pesticides  
Below a each input map is discussed and options to obtain data for them presented.

_pcmixsoil – PP concentration in the mixing layer_   
The pesticide concentrations at the start of the event are required as input. This can be obtained form field measurements but these will often not be available. An alternative is to use a long term pesticide fate model to estimate the field and mixing layer concentrations at the start of the rainfall event (e.g. PEARL (Van den Berg et al., 2016)). For an example application of long term modelling see [paper number 3].

_pcsoil1 – PP concentration in the deeper soil_  
The concentration in the deeper soil can be set equal to the mixing layer at the start. In theory it is expected that the concentration will decrease over the soil depth [refs]. The depth of the pesticide concentration in the soil should not exceed the modelled soil depth of OpenLISEM.  

_pcmixwat – DP concentration in the mixing layer_   
If the DP concentration in the mixing layer is not known from observations or a model, the dissolved phase concentration can be estimated by assuming the soil-water partitioning is in equilibrium, equilibration times are generally between minutes to several hours (Wauchope et al., 2002; Zhou et al., 2022) so if there is no pesticide application on the day of the rainfall event, this is a safe assumption.  

_pestmixdep – depth of the mixing layer_  
The thickness of the mixing layer is currently a calibration value. The exact value is difficult to measure, several lab experiments report values between 1 and 25 mm (Ahuja and Lehman, 1983; Havis et al., 1992).  

_pestsoildep1 – depth of the deeper soil pesticide profile_  
This value can be obtained from measurements or model output. In the current version of OLP only 1 soil layer and concentration is modelled, if needed this can be increased in future versions.  

All input maps should be stored in the map directory of the OpenLISEM run.

## Pesticide specific runfile options  

**Model setup**  
The sections with input or options in the runfile are marked as: `[section name] `.

In the `[General options]` section two options are available for pesticides:

```
Include Pesticides=0	   0 = not active, 1 = active
Report Pesticides=0	   0 = not active, 1 = save map series with pesticide dynamics
```

To activate the pesticide module set `Include Pesticides=1`. To save map series of all pesticide variables set `Report Pesticides=1`. This will store maps with concentrations, fluxes and masses of dissolved and particulate pesticides for the timestep as chosen in the user interface (see figure below). Saving map series will increase computation time and takes significant storage space, so only use this when the output is required!

<p align="center">
  <img width="367" height="325" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/pesticide_report_maps.png">
</p>

**Pesticide description**  
The next relevant section in the runfile is the `[Pesticides]` section, here the pesticide name, characteristics and the soil bulk density have to be given, for value ranges see table above.

The following options must be specified in the runfile:

```
Pesticide name=foobicide
Kd pesticide=0.0
Kfilm pesticide=0.0
ERbeta pesticide=-0.2
Kr pesticide=-1.0
ERmax pesticide=7.4
Rho mixing layer=0.0
```

**The map database** 

The required input maps are loaded based on the default names. If different map names are used, the correct files can be selected using the map database tab (see figure below). The map names can also be adjusted in the second `[Pesticides]` runfile section totally at the bottom of the runfile.

<p align="center">
  <img width="556" height="253" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/pesticide_map_db.png">
</p>