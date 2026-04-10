#-------------------------------------------------------------------------------
# This script creates a map with swale along the contour of the DEM, which can 
# be used when simulating flow with the dynamic wave in OpenLISEM. This map should 
# be used with the buffers function activated in the interface. The code will give 
# you an elevated dike with additional ditch one cell upstream. Besides that the 
# infiltration capacity in the ditch can be increased. 
#
#------------------------------
# code by Victor Jetten
# v1.0 - 2025-03-16
#------------------------------

binding

# input map
DEM = dem.map;
ksat = ksat1.map;
lu = landuse.map;

# output maps
buf = buffers.map;
swale = swale.map;
dike = dike.map;
 
areamap
dem.map;

initial

# change these settings to the design of the Swales
distance = 20; # (m) elevation distance between swales
height = 0.5; # (m) height dike of the swales
ditch = -0.3; # (m) depth of the ditches swale
ditch_ksat = 100; # (mm/h) infiltration in ditch of swale, -1 will use original ksat
lu_swales = -1; # set the landuse number for which swales will be implemented, -1 will use full area

# calculate 
contour = roundoff(dem.map mod distance);
swale_loc = if(lu_swales ge 0, if(scalar(lu) eq lu_swales, 1, 0), 1);

report dike = height * if(contour eq distance and windowminimum(contour,2*celllength()) eq 0, scalar(1), 0) * swale_loc;
report swale = ditch * if(contour eq 0 and windowmaximum(contour,3*celllength()) eq distance, scalar(1), 0) * swale_loc;
report ksat=if(ditch_ksat ge 0, if(swale lt 0, ditch_ksat,ksat), ksat);
report n.map=if(dike gt 0, 0.3,n.map);
report buf=dike+swale;

