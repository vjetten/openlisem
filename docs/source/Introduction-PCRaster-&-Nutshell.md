[PCRaster](https://pcraster.geo.uu.nl/) is a software for environmental modelling. It ss a collection of software targeted at the development and deployment of spatio-temporal environmental models. PCRaster is not developed to be a full-blown raster GIS. It lacks functionality for digitising, plotting and other typical GIS tasks. A user interface for PCRaster is available in [Nutshell](https://github.com/vjetten/NutShell/releases/tag/NutShell), the adds some functionality like displaying maps and manual editing of raster maps. For installation and setup of both programs see the [quick start](https://github.com/vjetten/openlisem/wiki/Getting-started#setup--installation).

**Nutshell User interface** 

The user interface of Nutshell consists of three main parts:
* The command window  
* The model Editor (script editor)  
* The explorer window  

<p align="center">
  <img width="563" height="338" src="https://github.com/vjetten/openlisem/blob/imgs_wiki/docs/imgs/nutshell_interface2.png">
</p>

Note that the working folder of the explorer window and the command window are not the same. The explorer window is controlled from the file structure on the bottom left of the interface, while the working folder of the command window is displayed on top. The current explorer directory can easily be set as the command window directory by pressing the button left from the displayed directory.  
The explorer window provides an easy method of displaying both regular PCRaster maps and time series of maps with Aguila. Double clicking any recognized map or timeseries of maps will automatically open this in Aguila. Timeseries of maps are furthermore displayed in blue, with only the first map of the timeserie visible within the explorer structure. When two files are selected, of which one is the digital elevation model, a single button press opens both maps in 3D drape style.
The model editor shows a currently opened PCRaster script file. These script files can be run by the pcrcalc application by using the play, pause and stop buttons at the top of the model editor.

**PCRaster commands in Nutshell**

A PCRaster command can be executed in the the command window of NutShell. A command will be executed in the map directory the is selected above the command window. For preparation of openLISEM three command groups are important:
* map format transformation
* making new maps
* raster calculations

_Map format transformation_  
The format of PCRaster maps is `.map`, but often input maps like a DEM have different formats, like `GeoTIFF` of `ASCII`. In PCRaster functions exist to transform maps between `ASCII` and `map` format with [`asc2map`](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/app_asc2map.html#index-0):

> `asc2map --clone mask.map -a  input.asc result.map`

with  
`--clone` - a command to link a base map with the same extent as the transformed map, this is always needed in PCRaster!  
`mask.map` - the 'clone' map, one of the maps you need to make at the start of each PCRaster project.  
`-a` - this command assumes the first 6 rows of the input file to contain extent information of the map for other options see the documentation.  
`input.asc` - the input `ASCII` file, can also have `.txt` extension.  
`result.map` - the map in the PCRaster format.  

The function also exists the other way round; [`map2asc`](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/app_map2asc.html#index-0).

_Making new maps_  
The function to make new maps is [`mapattr`](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/app_mapattr.html#index-0). In the openLISEM workflow this is mainly usefull to produce the base 'clone' map. This can be done by knowing the extent, origin and cellsize of your input data. All maps you will make, and all input maps from other sources will need to have exactly the same extent! To make a new map:

>`mapattr -s -R 19 -C 68 -B -P yb2t -x 12 -y -14.02 -l 0.8 clone.map`

for all the options, please consult the PCRaster documentation: [`mapattr`](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/app_mapattr.html#index-0).

_Raster calculations_  
This is the main use of PCRaster and is executed by starting a command with [`pcrcalc`](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/app_pcrcalc.html). As simple map calculation can be done by e.g.:

>`pcrcalc new.map = channelwidth.map * 2`  

This command will create a new map with value 2 time the input `channelwidth.map`. Many more extensive uses are [available](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/operations.html) including operations to produce all type of input maps for openLISEM, two important operations are: [`lddcreate()`](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/op_lddcreate.html) and [`lookupscalar()`](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/op_lookup.html). `lddcreate()` is used to prepare the [catchment maps 🔨 ](https://github.com/vjetten/openlisem/wiki/Preparing-Topography#the-ldd-map).

### lookup

The `lookupscalar()` command is used to combine a base map of for example land use of soil type, with values of variables that are explained by this. To prepare the land use properties, a land use type unit map and a land use properties table are needed. Using the `lookup` command, a unit map can be used to read a table for an entire map:

> `propertymap = --matrixtable lookupscalar(property.tbl, column, unit.map)`  

with  
* `--matrixtable` - a global option that allows the use of 2-dimentional tables.
* `property.tbl` - 2D table of soil properties. Row represents land use type, Column represents variable.  
* `column` - The column which is read from the table (determines the read property).  
* `unit.map` - The soil type unit map.  

An example of a `property.tbl` is given below:

(note that the   file only    contains data, and not the text) | nr | Random    Roughness | Manning's    n | Plant   Height | Cover
-- | -- | -- | -- | -- | --
| - | - | cm | - | m | -
| - | 0 | 1 | 2 | 3 | 4
Densely Vegetated Farming | 1 | 1.0 | 0.03 | 1.0 | 0.8
Eroded Agricultural Land | 2 | 0.7 | 0.03 | 0.5 | 0.5
Flatland Intensive Farming | 3 | 1.0 | 0.03 | 1.0 | 0.8
Grasslands | 4 | 0.5 | 0.10 | 0.2 | 1.0
Grasslands and Open Wood | 5 | 0.5 | 0.10 | 0.2 | 1.0

### PCRaster script files

To decrease preparation time and repeated work, the same [PCRaster script](https://pcraster.geo.uu.nl/pcraster/4.4.0/documentation/pcraster_manual/sphinx/secdyn.html#the-script) can generally be used when preparing a database. A full example script is provided below. Note that for this script to work, the necessary base maps, which are read in the binding section, should be present in the same folder. Some values within this script are furthermore catchment dependent and should be chosen based on the users research area.

> 📝 #comments start with a hashtag

```
#! --matrixtable --lddin

################################################################
# PCRASTER script for the generation of a LISEM input database #
# Victor Jetten 17/05/10                                       #
# Data for the GANSPOEL catchment                              #
# DESIRE LISEM course                                          #
################################################################


binding
##################
### input maps ###
##################

dem = gpdem.map;
# digital elevation model, area must be <= mask
unitmap = lur.map;
# field id's

texture = soils.map;
# texture/soil map

roads = roads.map;
# road map, 20 = tarred road (8 m wide), 21 is narrow dirt road (4 m wide)

chanmask= chanmask.map;
# mask for channel maps
####################
### input tables ###
####################
unittblsoil = soil.tbl; 
# table with crop and soil parameters for each field id 

# unitbase table layout #
#-----------------------#
# 01 ksat (mm/h)
# 02 porosity (cm3/cm3)
# 03 psi initial (cm)
# 04 initial moisture content (cm3/cm3)

unittblsurface= surface.tbl; 
# table with crop and soil parameters for each field id 

# unitbase table layout #
#-----------------------#
# 05 RR (cm)
# 06 Manning's n (-)
# 07 surface cover (-)
# 08 Crop height (m)
# 09 cohesion sol (kPa)
# 10 cohesion roots (kPa)
# 11 aggregate stability (number)

#######################
### input constants ###
#######################

Soildepth = 1000;
d50 = 30; # median texture loess = 30 mu
#channel properties:
Chancoh = 10; # high cohesion, kPa
Chanman = 0.2; # high man n, grass
Chanside = 0; # rectangular
Chanwidth = 2; # 2 meter
RoadWidth = 6; # 6 meter
ChanKsat = 20; # channel ksat for infil set to 20 mm/h

###################
### output maps ###
###################

# basic topography related maps
Ldd = ldd.map; # Local Drain Direction
area = area.map; # reference map for Lisem
grad = grad.map; # max slope 
id = id.map; # pluviograph influence zones
outlet = outlet.map; # location outlets and checkpoints

# impermeable roads
roadwidth = roadwidt.map;

# crop maps
coverc= per.map;
lai= lai.map;
cropheight= ch.map;
grass= grasswid.map;
# soil maps
ksat= ksat1.map;
psi= psi1.map;
pore= thetas1.map;
thetai= thetai1.map;
soildep= soildep1.map;
# maps for G&A 2nd layer
ksat2= ksat2.map;
psi2= psi2.map;
pore2= thetas2.map;
thetai2= thetai2.map;
soildep2= soildep2.map;
# surface maps
rr= rr.map;
mann= n.map;
stone= stonefrc.map; # crusted fraction, only used when option chosen in LISEM
crust= crustfrc.map;
comp= compfrc.map;
hard=hardsurf.map;
bufferid = bufferid.map;
buffervol=buffervol.map;
# erosion maps 
cohsoil = coh.map;
cohplant = cohadd.map;
D50 = d50.map;
aggrstab = aggrstab.map;

# channel maps
lddchan = lddchan.map;
chanwidth = chanwidt.map;
changrad = changrad.map;
chanman = chanman.map; 
chanside = chanside.map; 
chancoh = chancoh.map; 
chanksat = chanksat.map;

mask=mask.map;

initial


#################
### BASE MAPS ###
#################
# correct topo for local depressions
report Ldd = lddcreate (dem*mask, 1e20,1e20,1e20,1e20);
report outlet = pit(Ldd);
# reference catchment boundaries, based on watershed from outlet 
# OBSOLETE report area = catchment(Ldd, outlet);
# LDD is reference in later LISEM eversions
# sine gradient (-), make sure slope > 0.001
report grad = max(sin(atan(slope(dem*mask))),0.001);
#########################################
### MAPS WITH RAINFALL INFLUENCE ZONE ###
#########################################

report id = nominal(mask);
# use spreadzone for thiessen polygons when more than 1 rainfall station
#####################
### LAND USE MAPS ###
#####################

# fraction soil cover (including residue)
report coverc = lookupscalar(unittblsurface, 7, unitmap) * mask;
# crop height (m)
report cropheight = lookupscalar(unittblsurface, 8, unitmap) * mask;

# LAI of plants inside gridcell (m2/m2)
coverc = min(coverc, 0.95);
lai = ln(1-coverc)/-0.4;
report lai = if(coverc gt 0, lai/coverc, 0);
# or read from table:
#lookupscalar(unittbl, 9, unitmap) * mask;
###########################################################
### INFILTRATION MAPS for option one layer GREEN & AMPT ###
###########################################################

report ksat = lookupscalar(unittblsoil, 2, texture) * mask;
report pore = lookupscalar(unittblsoil, 2, texture) * mask;
report psi = abs(lookupscalar(unittblsoil, 3, texture)) * mask;
report thetai = lookupscalar(unittblsoil, 4, texture) * mask;
report soildep = scalar(200);

report ksat2 = 0.9*lookupscalar(unittblsoil, 2, texture) * mask;
report pore2 = lookupscalar(unittblsoil, 2, texture) * mask;
report psi2 = abs(lookupscalar(unittblsoil, 3, texture)) * mask;
report thetai2 = lookupscalar(unittblsoil, 4, texture) * mask;
report soildep2 = scalar(Soildepth);

#########################
### SOIL SURFACE MAPS ###
#########################

# micro relief, random roughness (=std dev in cm)
report rr = max(lookupscalar(unittblsurface, 5, unitmap) * mask, 0.01);
# Manning's n (-) 
# take from table
#report mann = lookupscalar(unittbl, 6, unitmap) * mask;
report mann = 0.051*rr+0.104*coverc;
# or use simple regression from Limburg data: CAREFULL this is not published

report crust=mask*0;
# crust fraction map, SWATRE option 2 in LISEM. Note that this demands an extra 
# profile definition in PROFILE.INP and PROFILE.MAP
report stone = 0 * mask;
# stone fraction 
report comp = 0*mask;
#fraction compacted
report hard = 0*mask;
#hard surface cells
report roadwidth = scalar(if(roads eq 20, 4, if(roads eq 21, 8, 0)))*mask;
# road width, 21 is tarred road = 8 m, dirt rods are 4 m wide

report bufferid = 0*mask;
report buffervol = 0*mask;


####################
### EROSION MAPS ###
####################

report D50 = 30*mask;
report cohsoil = lookupscalar(unittblsurface, 9, unitmap) * mask;
report cohplant = lookupscalar(unittblsurface, 10, unitmap) * mask;
report aggrstab = lookupscalar(unittblsurface, 11, unitmap) * mask;

####################
### CHANNEL MAPS ###
####################
chanmask=chanmask/chanmask;
report lddchan=lddcreate(dem*chanmask,1e20,1e20,1e20,1e20);
report changrad=max(0.001,sin(atan(slope(chanmask*dem))));
report chancoh=chanmask*scalar(Chancoh);
report chanman=chanmask*scalar(Chanman);
report chanside=chanmask*scalar(Chanside);
report chanwidth=chanmask*scalar(Chanwidth);
report hmxinit.map=mask*0;
report chanksat = chanmask*ChanKsat;
```