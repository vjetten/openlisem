#! --clone mask.map --lddin --matrixtable
#-------------------------------------------------------------------------------
# This script creates a new random roughness map to include additional roughness
# in the micro depression storage in OpenLISEM.
# MDS (micro depression storage) is a slope corrected volume of water that can
# be stored on the soil surface before runoff occurs. This is calculated in 
# OpenLISEM based on the RR (st.dev of micro roughness) and corrected for slope.
# A runoff decreasing measure is to make small dikes in between potato ridges,
# or barbuttes in French. This type of 'additional roughness' can be added to
# MDS with the following code.
#
# The additional roughness needs to be given as the addtional storage depth by 
# the features for a 0% slope.
#
#------------------------------
# code by Meindert Commelin
# v2.0 - 2025-03-14
#------------------------------

binding

# the following four maps are needed as input
# adjust the names of the maps after the '=' sign if you use different names.
	RRo = rr.map;			# the original random roughness map. [cm]
	slope = slope.map;		# the slope map for OpenLISEM: slope % divide by 100
    RRa = ar.map;           # the additional roughness map. [cm] 
	lu = landuse.map;		# a land unit map with one landuse for trenches
	
# the resulting new random roughness map
# adjust the name of the output file to your convenience.
	RRn = rr_adj.map;		# new random roughness map which includes addition storage by trenches
							# this should now be added to OpenLISEM as input RR map.

initial
# caluculate the additional roughness by trenches
# we use multiple helper steps for ease of reading
	sp = slope * 100;		# change slope from ratio to %
	mds = RRa * 10;		    # additional roughness to mm
	b = 0.243-0.012 * sp;
	x = (-(b)+sqrt(b**2 - 4*(0.01*-mds)))/(2*0.01); # AR in mm
	AR = x/10; 				# change x to cm

# combine both roughness maps and save
report	RRn = RRo + AR;

