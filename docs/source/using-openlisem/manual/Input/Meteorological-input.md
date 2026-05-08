<p align="center">
  <img width="588" height="385" src="../../../../images/input_meteo_options.png">
</p>

Meteorological data is given in a separate folder. Default name is 'rain'. The data is given per time interval which is indicated by "ddd:mmmm", Julian day number (1-366) followed by the minute in the day (0-1440) (note that 0 = 1440 = midnight).   

## Rainfall data
There are two ways rainfall data can be supplied:  

1. Rainfall gauging stations  
Values are given in $mm\ h^{-1}$ and the time interval should be hourly or shorter. openLISEM will work with larger time intervals (daily, weekly) but the results will be nonsense. Using rain gauges needs a map (ID.map) with areas to which the rain gauges correspond, using non-zero consecutive integer numbers (characters not allowed). Number 1 corresponds to the first gauge in the rainfall file, number 2 to the second. The rain gauge data file is given as:

```
 #some header, rainfall source
 4  
 time  
 gauge 1  
 gauge 2  
 gauge 3  
 100:0230   20.0  4.3  17.2  
 100:0240   10.0  14.3  7.2  
 100:0251   11.2  8.3  87.2  
 etc  
```

2. Rainfall maps (PCRaster format).  
Values are given in $mm\ h^{-1}$ and the map format and size must be exactly the same as all other maps in the database. The map ID.map is not used. The source can be satellite images or interpolated stations. The rainfall maps are stored in the rain folder and can have any name. If there are missing values in the model domain area, the model stops.
The rainfall file with list of maps has the following format:  

```
#some header, rainfall source  
2  
time  
maps  
100:0180   rainfall1.map  
100:0210   rainfall2.map  
100:0240   rainfall3.map  
etc  
```

## Evapotranspiration data

**NOTE** Evopotranspiration data is only required in continuous mode, while openLISEM is typically run in event mode this data is not required.

The format of the input files is exactly the same as for rainfall input, either station data with a map (ETID.map) or a series of potential Evapotranspiration ETp maps are used. This can be for instance Penman ET calculated from stations, or "latent evaporation flux" derived from satellite images.  
If the time interval is less than 1 day, the values are assumed to be in $mm\ h^{-1}$ (even if the time interval is larger than 1 hour).
If the interval is 1 day or more, the values are assumed to be in $mm\ day^{-1}$. In that case select the option "ET has daily or larger timesteps", and supply the average latitude (in degrees, e.g. $22.5^{\circ}$) of the model area. This will be used to calculate the sun angle and day length. The ETp provided in $mm\ day^{-1}$  is assumed to represent the daylight period and the instantaneous ETp per timestep is calculated according to a sine curve between sunrise and sunset (with the maximum at the solar noon).  
If the time interval is longer than a day (e.g. derived MODIS 8 daily values), you still have to provide average daily ETp in $mm\ day^{-1}$, and the daily ETp is simply repeated for the 8 days until the next timestep is found in the input file.  

## Snowmelt data
_Not yet implemented._
 
