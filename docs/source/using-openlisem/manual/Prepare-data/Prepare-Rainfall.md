Rainfall data can be prepared in three ways: uniform timeseries, Station based timeseries, map time series. For all methods, a temporal distribution of rain must first be obtained. When no high-resolution products are available, a synthetic rainfall distribution can be used. While the format for a uniform timeseries is relatively simple, the other methods require some GIS work. Station based rainfall can be prepared when both the station locations are known, and the rainfall intensity for each station is known. An example of the `ID.map` is given below.

<p align="center">
  <img width="602" height="538" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/prepare_rain_id.png">
</p>

The rainfall file has a specific format with three sections:
1. A header line starting with `#`, here you can write information like, date, location and data source.
1. Number and names of the columns, each column name is written on a new line.
1. the actual data separated by a white space: `   `  

The will look like this for the `ID.map` shown:

```
# description including: 'location', 'date' and 'data source'
11
time		
station  1
station  2
station  3
station  4
station  5
station  6
station  7
station  8
station  9
station  10

1 12	   12      12	   12      12	   12      12	   12      12	   12      
2 12       16      12       16     12       16      12       16      12       16       
3 16       20      16       20     16       20      16       20      16       20       
4 20       32      20       32     20       32      20       32      20       32       
5 32       35      32       35     32       35      32       35      32       35       
6 35       22      35       22     35       22      35       22      35       22       
7 22       33      22       33     22       33      22       33      22       33       
8 33       28      33       28     33       28      33       28      33       28       
9 28       16.57   28       16.57  28       16.57   28       16.57   28       16.57    
10 16.57   20.57    16.57   20.57   16.57   20.57    16.57   20.57    16.57   20.57 
```

The areas in the ID map can be dependent on the DEM, or each pixel can be assigned to the closest rainfall station. To create an ID map from a map with rainfall stations, the following code can be used:

> `id = spreadzone(points,0,friction)`

The [`spreadzone()`](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/op_spreadzone.html) function contains three arguments, of which two are used:
* Points - A boolean map containing the locations of the rainfall stations  
* Friction - The friction map  
The friction acts as a multiplier when the algorithm tries to determine the shortest distance from a cell to each rainfall station. With a uniform friction map of 1, each cell is assigned to the closest rainfall station. When the slope is used as friction, topography based rainfall zones can be created. From the same information (rainfall stations and rainfall intensity for each station) spatial rainfall maps can also be made by the user through, for example, inverse distance interpolation.

> `rain = inversedistance(mask, rainst,order,0,0)`

The [`inversedistance()`](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/op_inversedistance.html#index-1) function contains five arguments, of which three are used:
* Mask - A mask indicating where a value should be assigned  
* Rainst - A map containing MV (missing values) except on the locations of rainfall stations. There the rainfall intensity should be present.
* Order - The order of the interpolation (usually 2)  

