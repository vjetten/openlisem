/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
\file LisUIDefaultNames.cpp
\brief Default map names and descriptions DEFmaps, Default runfile variable names namelist
*/

#include "lisemqt.h"

//--------------------------------------------------------------------
void lisemqt::DefaultMapnames()
{
    DEFmaps.clear();
    //VJ 110417 delete all at start, needed when reseatAll;

    //# interface maplist, DO NOT CHANGE if you don't know what you are doing
    //# syntax: branch level; keyword; default mapname; description; variable name
    DEFmaps.append("0;Rainfall");
    DEFmaps.append("2;ID;ID.map;Raingauge zone ID numbers, correspond to columns (1,2,...) in rainfall file;ID");
    DEFmaps.append("0;Catchment");
    DEFmaps.append("2;DEM;dem.map;Digital elevation model (m);dem");
    DEFmaps.append("2;Gradient;grad.map;Sine of slope gradient in direction of flow;grad");
    DEFmaps.append("2;LDD;ldd.map;Local surface Drainage Direction network;ldd");
    DEFmaps.append("2;Outlet;outlet.map;Main catchment outlet corresponding to LDD map;outlet");
    DEFmaps.append("2;Points;outpoint.map;Reporting points for hydrograph/sedigraph (1,2,3,...);outpoint");
    DEFmaps.append("2;Watersheds;ws.map;watersheds (1,2,3,...);watershed");
    DEFmaps.append("0;Landuse");
    DEFmaps.append("2;Units;landunit.map;Classified land unit map (integers 0-n) for output of erosion values;landunit");
    DEFmaps.append("2;Cover;per.map;Fraction surface cover by vegetation and residue (-);cover");
    DEFmaps.append("2;Litter;litter.map;Fraction of surface cover by litter/herbs under trees (-);litter");
    DEFmaps.append("2;LAI;lai.map;Leaf area index of the plant cover in a gridcell (m2/m2);lai");
    DEFmaps.append("2;Height;ch.map;Plant height (m);ch");
    DEFmaps.append("2;Road width;roadwidt.map;Width of impermeable roads (m);road");
    DEFmaps.append("2;Grass strips;grasswid.map;Width of grass strips (m);grasswidth");
    DEFmaps.append("2;Canopy storage;smax.map;Maximum canopy storage (mm);smax");
    DEFmaps.append("0;Surface");
    DEFmaps.append("2;RR;rr.map;Random Roughness (here standard deviation of heights) (cm);rr");
    DEFmaps.append("2;n;n.map;Manning's n (-);manning");
    DEFmaps.append("2;Stoniness;stonefrc.map;Fraction covered by stones (affects only splash det.) (-);stonefrc");
    DEFmaps.append("2;Crust;crustfrc.map;Fraction of gridcell covered with Crust (-) (see also ksat crust);crustfrc");
    DEFmaps.append("2;Compacted;compfrc.map;Fraction of gridcell compacted (e.g. wheeltracks)(-) (see also ksat compacted);compfrc");
    DEFmaps.append("2;Hard Surface;hardsurf.map;No interception/infiltration/detachment (fraction 0-1);hardsurf");
    DEFmaps.append("0;Erosion");
    DEFmaps.append("2;Cohesion;coh.map;Cohesion (kPa);coh");
    DEFmaps.append("2;Cohesion;cohadd.map;Extra cohesion factor by e.g. plant roots (kPa);cohadd");
    DEFmaps.append("2;Aggregates;aggrstab.map;Aggregate stability for splash erosion (-);aggrstab");
    DEFmaps.append("2;D50;d50.map;Median of the texture of the suspendeed matter (mu);d50");
    DEFmaps.append("2;D90;d90.map;Median of the texture of the suspendeed matter (mu);d90");
    DEFmaps.append("2;Material;detmat.map;Detacheable material per square meter (kg/m2) (-1 = infinite);detmat");
    DEFmaps.append("2;MixingDepth;sedmixdeth.map; Mixing depth for deposited sediment (m);sedmixdepth");

    DEFmaps.append("0;Debris Flow");
    DEFmaps.append("1;Slope Stability");
    DEFmaps.append("2;SoilDensity;soildensity.map;Soil density (m);soildensity");
    DEFmaps.append("2;SoilIFA;soilifa.map;Soil internal friction angle;soilifa");
    DEFmaps.append("2;SoilRockFraction;soilrockfraction.map;Soil rock fraction;soilrockfraction");
    DEFmaps.append("2;SoilRockSize;soilrocksize.map;Soil rock fraction;soilrocksize");
    DEFmaps.append("2;FailureMask;failuremask.map;Mask indicating possible slope failure;failuremask");

    DEFmaps.append("1;Loose Material");
    DEFmaps.append("2;DebrisMaterial;debrismaterial.map;Soil debris material;debrismaterial");
    DEFmaps.append("2;RockSize;rocksize.map;debris material rock size (m);rocksize");
    DEFmaps.append("2;RockDensity;rockdensity.map;debris material density (m);rockdensity");
    DEFmaps.append("2;RockIFA;rockifa.map;debris material internal friction angle;rockifa");



    DEFmaps.append("0;Infiltration");
    DEFmaps.append("1;Swatre");
    DEFmaps.append("2;Profile soil;profile.map;ID numbers corresponding to land units in profile table;profmap");
    DEFmaps.append("2;Prof. Crust;profcrst.map;ID numbers of crusted soils (using also profile table);profcrst");
    DEFmaps.append("2;Prof. Wheel;profwltr.map;ID numbers of compacted wheel tracks (using also profile table);profwltr");
    DEFmaps.append("2;Prof. Grass;profgras.map;ID numbers of grasstrips (using also profile table);profgras");
    DEFmaps.append("2;Initial suction;inithead;initial matrix potential (cm) of layers 001 to nnn (filename witout extension);inithead");
    DEFmaps.append("2;Repellency;repel.map;Gridcells included in water repellency (1/0);repelcell");
    DEFmaps.append("1;1st layer Green&Ampt/Smith&Parlange");
    DEFmaps.append("2;Ksat1;ksat1.map;Layer 1: Saturated Hydraulic Conductivity (mm/h);ksat1");
    DEFmaps.append("2;Psi1;psi1.map;Layer 1: Average suction at the wetting front (cm);psi1");
    DEFmaps.append("2;Thetas1;thetas1.map;Layer 1: Porosity (-);thetas1");
    DEFmaps.append("2;Thetai1;thetai1.map;Layer 1: Initial moisture content (-);thetai1");
    DEFmaps.append("2;Depth1;soildep1.map;Layer 1: Depth (mm) to bottom of layer 1;soildep1");
    DEFmaps.append("1;2nd layer Green&Ampt/Smith&Parlange");
    DEFmaps.append("2;Ksat2;ksat2.map;Layer 2: Saturated Hydraulic Conductivity (mm/h);ksat2");
    DEFmaps.append("2;Psi2;psi2.map;Layer 2: Average suction at the wetting front (cm);psi2");
    DEFmaps.append("2;Thetas2;thetas2.map;Layer 2: Porosity (-);thetas2");
    DEFmaps.append("2;Thetai2;thetai2.map;Layer 2: Initial moisture content (-);thetai2");
    DEFmaps.append("2;Depth2;soildep2.map;Layer 2: Depth (mm) to bottom of layer 2;soildep2");
    DEFmaps.append("1;Ksat subtraction");
    DEFmaps.append("2;Ksat1;ksat1.map;Saturated Hydraulic Conductivity (mm/h);ksat1simple");
    DEFmaps.append("1;Special surfaces");
    DEFmaps.append("2;Ksat Crust;ksatcrst.map;Ksat of crusts (all models except SWATRE) (mm/h);ksatcrst");
    DEFmaps.append("2;Ksat Compact;ksatcomp.map;Ksat of compacted areas (all models except SWATRE) (mm/h);ksatcomp");
    DEFmaps.append("2;Ksat Grass;ksatgras.map;Ksat of grassstrips (all models except SWATRE) (mm/h);ksatgras");
    DEFmaps.append("0;Channels");
    DEFmaps.append("1;Channel properties");
    DEFmaps.append("2;LDD;lddchan.map;LDD of main channel (must be 1 branch connected to the outlet);lddchan");
    DEFmaps.append("2;Width;chanwidt.map;Channel width (m);chanwidth");
    DEFmaps.append("2;Side angle;chanside.map;Channel side angle (tan angle  channel side and surface: 0 is rectangular);chanside");
    DEFmaps.append("2;Gradient;changrad.map;Slope gradient of channel bed (-);changrad");
    DEFmaps.append("2;N;chanman.map;Mannings n of channel bed (-);chanman");
    DEFmaps.append("2;Cohesion;chancoh.map;Cohesion of channel bed (kPa);chancoh");
    DEFmaps.append("2;Channelmaterial;chandetmat.map;Detacheable material per square meter (kg/m2) (-1 = infinite);chandetmat");
    DEFmaps.append("2;ChannelMixingDepth;chansedmixdeth.map; Mixing depth for deposited sediment in channel (m);chansedmixdepth");

    DEFmaps.append("0;Surface Flow");
    DEFmaps.append("1;Flow barriers");
    DEFmaps.append("2;DemBarriers;barriers.map;Barrier depth (is added to dem) (m);barriers");
    DEFmaps.append("2;FlowbarrierIndex;flowbarrierindex.map;Flow barrier index (reads properties from table);flowbarrierindex");
    DEFmaps.append("2;MaxVolume;maxvolume.map; Maximum flow volume;maxvol");
    DEFmaps.append("1;Flow Limiting");
    DEFmaps.append("2;ChannelDepth;chandepth.map;Channel depth, zero (0) depth is considered infinite (m);chandepth");
    DEFmaps.append("2;ChannelLevee;chanlevee.map;Height of small channel levee on both sides of the channel (m);chanlevee");
    DEFmaps.append("2;ChannelMaxVolume;channelmaxvolume.map; Maximum flow volume for Channel;channelmaxvol");
    DEFmaps.append("2;ChannelConnected;channelconnected.map; Is the channel connected to overland flow;channelconnected");
    DEFmaps.append("1;Channel Infil");
    DEFmaps.append("2;Ksat;chanksat.map;Infiltration rate of channel bed (mm/h);chanksat");
    DEFmaps.append("1;Channel BaseFLow");
    DEFmaps.append("2;BaseFlow;baseflow.map; base flow discharges (m3/s);baseflow");


    DEFmaps.append("1;Initial Volume");
    DEFmaps.append("2;InitiationTime;initiationtime.map;initial time (min);initiationtime");
    DEFmaps.append("2;InitialFVolume;initialfvolume.map;initial fluid volume (m3);initialfvolume");
    DEFmaps.append("2;InitialSVolume;initialsvolume.map;initial solid volume (m3);initialsvolume");
    DEFmaps.append("2;InitialSDensity;initialsdensity.map;initial solid density (m3);initialsdensity");
    DEFmaps.append("2;InitialSRocksize;initialsrocksize.map;initial solid rocksize (m3);initialsrocksize");
    DEFmaps.append("2;InitialSIFA;initialsifa.map;initial internal friction angle;initialsifa");

    DEFmaps.append("1;Forced Volume condition");
    DEFmaps.append("2;forcedfvolume;forcedfvolume.map;forced fluid volume (m3);forcedfvolume");
    DEFmaps.append("2;ForcedSVolume;forcedsvolume.map;forced solid volume (m3);forcedsvolume");
    DEFmaps.append("2;ForcedSDensity;forcedsdensity.map;forced solid density (m3);forcedsdensity");
    DEFmaps.append("2;ForcedSRocksize;forcedsrocksize.map;forced solid rocksize (m3);forcedsrocksize");
    DEFmaps.append("2;ForcedSIFA;forcedsifa.map;forced internal friction angle;forcedsifa");

    DEFmaps.append("0;FlowBarriers");
    DEFmaps.append("2;FlowBarrierIndex;flowbarrierindex.map;An index value, indicating which flow barrier properties will be used (-);flowbarrierindex");

    DEFmaps.append("0;Snowmelt");
    DEFmaps.append("2;Snowmelt ID;snowid.map;Snowmelt zone ID number for snowmelt file starting with 1 (0 is non-snow area);SnowID");

    //houses
    DEFmaps.append("0;Houses");
    DEFmaps.append("2;House Cover;housecover.map;Fraction of hard roof surface per cell (-);housecover");
    DEFmaps.append("2;Roof Storage;roofstore.map;Size of interception storage of rainwater on roofs (mm);roofstore");
    DEFmaps.append("2;Drum Store;drumstore.map;Size of storage of rainwater drums (m3);drumstore");

}
//---------------------------------------------------------------------------
// fill namelist with default runfile values and structure
// runfile has structure: name=value
void lisemqt::defaultRunFile()
{
    int i;
    for (i = 0; i < NUMNAMES; i++)
    {
        namelist[i].name.clear();
        namelist[i].value.clear();
    }


    //VJ CLEANED UP ORDER OF VARIABLES 160409

    i = 0;
    namelist[i++].name = QString("[openLISEM runfile version 4]");
    namelist[i++].name = QString("");

    //###
    namelist[i++].name = QString("[Input]");
    namelist[i++].name = QString("Work Directory");
    namelist[i++].name = QString("Map Directory");
    namelist[i++].name = QString("Include Rainfall");
    namelist[i++].name = QString("Rainfall Directory");
    namelist[i++].name = QString("Rainfall file");
    namelist[i++].name = QString("Include Snowmelt");
    namelist[i++].name = QString("Snowmelt Directory");
    namelist[i++].name = QString("Snowmelt file");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Output]");
    namelist[i++].name = QString("Result Directory");
    namelist[i].value = QString("totals.txt");
    namelist[i++].name = QString("Main results file");
    namelist[i].value = QString("hydrograph.csv");
    namelist[i++].name = QString("Filename point output");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Report point output separate");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Timeplot as CSV");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Timeplot as PCRaster");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Report point output for SOBEK");
    namelist[i].value = QString("10\01\01");
    namelist[i++].name = QString("SOBEK date string");
    namelist[i].value = QString("rainfall.map");
    namelist[i++].name = QString("Rainfall map");
    namelist[i].value = QString("interception.map");
    namelist[i++].name = QString("Interception map");
    namelist[i].value = QString("infiltration.map");
    namelist[i++].name = QString("Infiltration map");
    namelist[i].value = QString("runoff.map");
    namelist[i++].name = QString("Runoff map");
    namelist[i].value = QString("rofraction.map");
    namelist[i++].name = QString("Runoff fraction map");
    namelist[i].value = QString("chandism3.map");
    namelist[i++].name = QString("Channel discharge map");
    namelist[i].value = QString("eros.map");
    namelist[i++].name = QString("Erosion map");
    namelist[i].value = QString("depo.map");
    namelist[i++].name = QString("Deposition map");
    namelist[i].value = QString("soilloss.map");
    namelist[i++].name = QString("Soilloss map");
    namelist[i].value = QString("totlandunit.txt");
    namelist[i++].name = QString("Filename landunit output");
    namelist[i].value = QString("chandet.map");
    namelist[i++].name = QString("Channel detachment map");
    namelist[i].value = QString("chandep.map");
    namelist[i++].name = QString("Channel deposition map");
    namelist[i].value = QString("WHmax.map");
    namelist[i++].name = QString("WH max level map");
    namelist[i].value = QString("floodmax.map");
    namelist[i++].name = QString("Flood level map");
    namelist[i].value = QString("floodtime.map");
    namelist[i++].name = QString("Flood time map");
    namelist[i].value = QString("floodstart.map");
    namelist[i++].name = QString("Flood start time");
    namelist[i].value = QString("floodmaxv.map");
    namelist[i++].value = QString("Flood Max V");
    namelist[i].value = QString("channelmaxq.map");
    namelist[i++].name = QString("Channel Max Q");
    namelist[i].value = QString("channelmaxhw.map");
    namelist[i++].name = QString("Channel Max WH");
    namelist[i].value = QString("floodstats.csv");
    namelist[i++].name = QString("Flood stats");

    namelist[i].value = QString("maxdebrisflowheight.map");
    namelist[i++].name = QString("Maximum Debris Flow Height Map");
    namelist[i].value = QString("maxdebrisflowvelocity.map");
    namelist[i++].name = QString("Maximum Debris Flow Velocity Map");
    namelist[i].value = QString("debrisflowstart.map");
    namelist[i++].name = QString("Debris Flow Start Map");
    namelist[i].value = QString("entrainment.map");
    namelist[i++].name = QString("Entrainment Map");
    namelist[i].value = QString("slopefailure.map");
    namelist[i++].name = QString("Slope Failure Map");
    namelist[i].value = QString("minimumsafetyfactor.map");
    namelist[i++].name = QString("Minimum Safety Factor Map");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Simulation times]");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Begin time");
    namelist[i].value = QString("100");
    namelist[i++].name = QString("End time");
    namelist[i].value = QString("0.15");
    namelist[i++].name = QString("Timestep");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[General options]");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Include Rainfall");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include snowmelt");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("No Erosion simulation");  // replaced below but leave in for older runfiles
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include Erosion simulation");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Advanced sediment");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include main channels");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include channel infil");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include channel baseflow");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Hard Surfaces");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Include road system");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include house storage");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include raindrum storage");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Interception]");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Use canopy storage map");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Canopy storage equation");
    namelist[i].value = QString("0.05");
    namelist[i++].name = QString("Stemflow fraction");
    namelist[i].value = QString("0.45");
    namelist[i++].name = QString("Canopy Openess");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include litter interception");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Litter interception storage");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Infiltration]");
    namelist[i].value = QString("3");
    namelist[i++].name = QString("Infil Method");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include compacted");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include crusts");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Impermeable sublayer");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include percolation");
// obsolete
//    namelist[i].value = QString("0");
//    namelist[i++].name = QString("Subsoil drainage");
    namelist[i].value = QString("c:\\");
    namelist[i++].name = QString("Table Directory");
    namelist[i].value = QString("profile.inp");
    namelist[i++].name = QString("Table File");
    namelist[i].value = QString("0.01");
    namelist[i++].name = QString("SWATRE internal minimum timestep");
    namelist[i].value = QString("");
    namelist[i++].name = QString("Matric head files");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Geometric mean Ksat");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Use Water Repellency");
    namelist[i].value = QString("1.20");
    namelist[i++].name = QString("Water Repellency A");
    namelist[i].value = QString("0.3");
    namelist[i++].name = QString("Water Repellency B");
    namelist[i].value = QString("0.12");
    namelist[i++].name = QString("Water Repellency C");
    namelist[i].value = QString("1.00");
    namelist[i++].name = QString("Water Repellency D");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include tile drains");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Surface Flow]");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Enable Solid Phase");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Enable Entrainment");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Enable Slope Stability");

    namelist[i].value = QString("0");
    namelist[i++].name = QString("Enable Slope Failure");
    namelist[i].value = QString("0.9");
    namelist[i++].name = QString("Minimum Safety Factor");
    namelist[i].value = QString("1.3");
    namelist[i++].name = QString("Maximum Safety Factor");
    namelist[i].value = QString("3.5");
    namelist[i++].name = QString("Maximum safety factor for display");
    namelist[i].value = QString("0.000122");
    namelist[i++].name = QString("Entrainment Coefficient");
    namelist[i].value = QString("0.05");
    namelist[i++].name = QString("Minimum Entrainment Height");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Minimum Failure Height");

    namelist[i].value = QString("1");
    namelist[i++].name = QString("Spatially Dynamic Timestep");

    namelist[i].value = QString("0");
    namelist[i++].name = QString("Enable Levees");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Enable Barriers");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Enable Flow Barriers");
    namelist[i].value = QString("flowbarriers.txt");
    namelist[i++].name = QString("Flow barrier table filename");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include Initial FluidSolid Mixture");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include Forced FluidSolid Mixture");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Incldue Maximum ChannelVolume");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Incldue Maximum Volume");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Flow Minimum Timestep");
    namelist[i].value = QString("2.0");
    namelist[i++].name = QString("Kinematic Timestep Power");
    namelist[i].value = QString("0.25");
    namelist[i++].name = QString("Surface Flow Courant Factor");
    namelist[i].value = QString("2");
    namelist[i++].name = QString("Surface Flow Scheme");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Drag Power Law Coefficient");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Viscosity Alpha");
    namelist[i].value = QString("20");
    namelist[i++].name = QString("Viscosity Beta");
    namelist[i].value = QString("0.0");
    namelist[i++].name = QString("Minimal Flood Water Depth");
    namelist[i].value = QString("0.0");
    namelist[i++].name = QString("Minimum Debris Flow Volumetric Sediment Fraction");
    namelist[i].value = QString("0.0");
    namelist[i++].name = QString("Include channel flooding");
    namelist[i].value = QString("0.0");
    namelist[i++].name = QString("Use HLL2");
    
    namelist[i].value = QString("0.0");
    namelist[i++].name = QString("Suspended Viscosity");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Lax Multiplier");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Friction force correction");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Erosion Cohesion Calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Erosion Grain Size Calibration");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Kinetic Energy]");
    namelist[i].value = QString("1,28.300,0.520,0.042");
    namelist[i++].name = QString("KE parameters EQ1");
    namelist[i].value = QString("0,8.950,8.440");
    namelist[i++].name = QString("KE parameters EQ2");
    namelist[i].value = QString("0,7.600,0.220");
    namelist[i++].name = QString("KE parameters EQ3");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("KE time based");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Detachment efficiency");
    //    namelist[i].value = QString("0");
    //    namelist[i++].name = QString("Detachment stoniness");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Use material depth");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Sediment]");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Advanced sediment configuration");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("BL method");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("SS method");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Estimate grain size distribution");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Number of grain size classes (simulated)");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Read grain distribution maps");
    namelist[i].value = QString("2,20,50,125,250,500");
    namelist[i++].name = QString("Grain size class maps");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Use material depth");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Sigma diffusion");

    // not active for user!!!
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Limit TC");
    namelist[i].name = QString("0");
    namelist[i++].name = QString("Limit Deposition TC");
    // namelist[i].value = QString("1");
    // namelist[i++].name = QString("All water and sediment to outlet");
    //    namelist[i].value = QString("0");
    //    namelist[i++].name = QString("No erosion at outlet");
    //    namelist[i].value = QString("0");
    //    namelist[i++].name = QString("Alternative flow detachment");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Conservation]");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include grass strips");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Grassstrip Mannings n");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include buffers");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Buffers impermeable");
    namelist[i].value = QString("1400");
    namelist[i++].name = QString("Sediment bulk density");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include Sediment traps");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Calibration]");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Ksat calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("N calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Theta calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Psi calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Channel Ksat calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Channel N calibration");

    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Erosive Power Calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Transport Capacity Calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Settling Velocity Calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Yield Stress Calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Dynamic Viscosity Calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Drag Force Calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Solid Phase Friction Calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Deposition Criteria Calibration");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Limit Failure");


    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Soil Cohesion Calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Soil Internal Friction Angle Calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Soil Depth Calibration");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Create Stable Initial Safety Factor");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Minimum Safety Factor Calibration");

    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Splash Delivery Ratio");
    namelist[i].value = QString("0.5");
    namelist[i++].name = QString("Particle Cohesion of Deposited Layer");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[OpenGL visualization]");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Light_Ambient_R");
    namelist[i].value  = QString("1");
    namelist[i++].name = QString("Light_Ambient_G");
    namelist[i].value  = QString("1");
    namelist[i++].name = QString("Light_Ambient_B");
    namelist[i].value  = QString("0.3");
    namelist[i++].name = QString("Light_Ambient_A");
    namelist[i].value  = QString("1");
    namelist[i++].name = QString("Light_Directional_R");
    namelist[i].value  = QString("1");
    namelist[i++].name = QString("Light_Directional_G");
    namelist[i].value  = QString("1");
    namelist[i++].name = QString("Light_Directional_B");
    namelist[i].value  = QString("0.7");
    namelist[i++].name = QString("Light_Directional_A");
    namelist[i].value  = QString("-1");
    namelist[i++].name = QString("Light_Directional_X");
    namelist[i].value  = QString("-1");
    namelist[i++].name = QString("Light_Directional_Y");
    namelist[i].value  = QString("-1");
    namelist[i++].name = QString("Light_Directional_Z");
    namelist[i].value  = QString("1");
    namelist[i++].name = QString("Surface_Draw");
    namelist[i].value  = QString("100");
    namelist[i++].name = QString("Surface_Micro_Elevation_Scale");
    namelist[i].value  = QString("15000");
    namelist[i++].name = QString("Surface_Mipmap_Distance_1");
    namelist[i].value  = QString("50000");
    namelist[i++].name = QString("Surface_Mipmap_Distance_2");
    namelist[i].value  = QString("0.3");
    namelist[i++].name = QString("Surface_Vegetated_Small_Color_R");
    namelist[i].value  = QString("0.5");
    namelist[i++].name = QString("Surface_Vegetated_Small_Color_G");
    namelist[i].value  = QString("0");
    namelist[i++].name = QString("Surface_Vegetated_Small_Color_B");
    namelist[i].value = QString("0.06");
    namelist[i++].name = QString("Surface_Vegetated_Large_Color_R");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Surface_Vegetated_Large_Color_G");
    namelist[i].value = QString("0.04");
    namelist[i++].name = QString("Surface_Vegetated_Large_Color_B");
    namelist[i].value = QString("0.20");
    namelist[i++].name = QString("Surface_Bare_Color_R");
    namelist[i].value = QString("0.16");
    namelist[i++].name = QString("Surface_Bare_Color_G");
    namelist[i].value = QString("0.03");
    namelist[i++].name = QString("Surface_Bare_Color_B");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Surface_Roads_Color_R");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Surface_Roads_Color_R");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Surface_Roads_Color_R");
    namelist[i].value = QString("0.3");
    namelist[i++].name = QString("Surface_Buildings_Color_R");
    namelist[i].value = QString("0.3");
    namelist[i++].name = QString("Surface_Buildings_Color_G");
    namelist[i].value = QString("0.3");
    namelist[i++].name = QString("Surface_Buildings_Color_B");
    namelist[i].value = QString("0.15");
    namelist[i++].name = QString("Surface_Erosion_Color_R");
    namelist[i].value = QString("0.14");
    namelist[i++].name = QString("Surface_Erosion_Color_G");
    namelist[i].value = QString("0.04");
    namelist[i++].name = QString("Surface_Erosion_Color_B");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Surface_Erosion_Color_A");
    namelist[i].value = QString("0.35");
    namelist[i++].name = QString("Surface_Deposition_Color_R");
    namelist[i].value = QString("0.33");
    namelist[i++].name = QString("Surface_Deposition_Color_G");
    namelist[i].value = QString("0.13");
    namelist[i++].name = QString("Surface_Deposition_Color_B");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Surface_Deposition_Color_A");
    namelist[i].value = QString("100");
    namelist[i++].name = QString("Water_Draw");
    namelist[i].value = QString("100");
    namelist[i++].name = QString("Water_Reflectivity");
    namelist[i].value = QString("100");
    namelist[i++].name = QString("Water_Refractivity");
    namelist[i].value = QString("100");
    namelist[i++].name = QString("Water_Velocity_Scale");
    namelist[i].value = QString("100");
    namelist[i++].name = QString("Water_Micro_Elevation_Scale");
    namelist[i].value = QString("100");
    namelist[i++].name = QString("Water_Transparancy");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Water_Deep_Color_R");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Water_Deep_Color_G");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Water_Deep_Color_B");
    namelist[i].value = QString("0.8");
    namelist[i++].name = QString("Water_Deep_Color_A");
    namelist[i].value = QString("0.8");
    namelist[i++].name = QString("Water_Shallow_Color_R");
    namelist[i].value = QString("0.6");
    namelist[i++].name = QString("Water_Shallow_Color_G");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Water_Shallow_Color_B");
    namelist[i].value = QString("0.05");
    namelist[i++].name = QString("Water_Shallow_Color_A");
    namelist[i].value = QString("0.15");
    namelist[i++].name = QString("Water_Sediment_Color_R");
    namelist[i].value = QString("0.14");
    namelist[i++].name = QString("Water_Sediment_Color_G");
    namelist[i].value = QString("0.04");
    namelist[i++].name = QString("Water_Sediment_Color_B");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Water_Sediment_Color_A");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Roads_Draw");
    namelist[i].value = QString("2000");
    namelist[i++].name = QString("Roads_Distance");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Buildings_Draw");
    namelist[i].value = QString("200000");
    namelist[i++].name = QString("Buildings_Distance");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Trees_Draw");
    namelist[i].value = QString("200000");
    namelist[i++].name = QString("Trees_Distance");
    namelist[i].value = QString("1000000");
    namelist[i++].name = QString("Trees_Instances");
    namelist[i].value = QString("1500");
    namelist[i++].name = QString("Trees_Increment");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Grass_Draw");
    namelist[i].value = QString("10000");
    namelist[i++].name = QString("Grass_Distance");
    namelist[i].value = QString("200000");
    namelist[i++].name = QString("Grass_Instances");
    namelist[i].value = QString("200.0");
    namelist[i++].name = QString("Grass_Increment");
    namelist[i].value = QString("100.0");
    namelist[i++].name = QString("Grass_Vertical_Scale");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Output maps]");
    namelist[i++].name = QString("Runoff maps in l/s/m");
    namelist[i++].name = QString("Timeseries as PCRaster");
    namelist[i++].name = QString("Erosion map units (0/1/2)");
    namelist[i++].name = QString("Regular runoff output");
    namelist[i++].name = QString("Output interval");
    namelist[i++].name = QString("User defined output");
    namelist[i++].name = QString("Output times");
    namelist[i++].name = QString("CheckOutputMaps");
    namelist[i++].name = QString("CheckOutputMapsNUT");
    namelist[i++].name = QString("CheckOutputMapsMC");
    namelist[i++].name = QString("CheckOutputMapsGUL");

    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Texture classes]");
    namelist[i].name = QString("2,16,32,53,75,105");
    namelist[i++].name = QString("ClassMu");

    // output maps have standard names
    // input maps names are defined in DEFmaps
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[map names]");
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[OutputBASIC]");
    namelist[i].value = QString("ro");
    namelist[i++].name = QString("OUTRUNOFF");
    namelist[i].value = QString("conc");
    namelist[i++].name = QString("OUTCONC");
    namelist[i].value = QString("wh");
    namelist[i++].name = QString("OUTWH");
    namelist[i].value = QString("roc");
    namelist[i++].name = QString("OUTRWH");
    namelist[i].value = QString("tc");
    namelist[i++].name = QString("OUTTC");
    namelist[i].value = QString("det");
    namelist[i++].name = QString("OUTEROS");
    namelist[i].value = QString("dep");
    namelist[i++].name = QString("OUTDEPO");
    namelist[i].value = QString("sloss");
    namelist[i++].name = QString("OUTSOILLOSS");
    namelist[i].value = QString("velo");
    namelist[i++].name = QString("OUTVELO");
    namelist[i].value = QString("inf");
    namelist[i++].name = QString("OUTINF");
    namelist[i].value = QString("sstor");
    namelist[i++].name = QString("OUTSS");
    namelist[i].value = QString("chvol");
    namelist[i++].name = QString("OUTCHVOL");
    namelist[i].value = QString("Qtile");
    namelist[i++].name = QString("OUTTILED");
    namelist[i].value = QString("hmx");
    namelist[i++].name = QString("OUTHMX");
    namelist[i].value = QString("qf");
    namelist[i++].name = QString("OUTQF");
    namelist[i].value = QString("vf");
    namelist[i++].name = QString("OUTVF");
    namelist[i].value = QString("hmxwh");
    namelist[i++].name = QString("OUTHMXWH");
    namelist[i].value = QString("sed");
    namelist[i++].name = QString("OUTSED");

    namelist[i].value = QString("safa");
    namelist[i++].name = QString("OUTSAFETYFACTOR");
    namelist[i].value = QString("slfa");
    namelist[i++].name = QString("OUTSLOPEFAILURE");
    namelist[i].value = QString("dfh");
    namelist[i++].name = QString("OUTDFHEIGHT");
    namelist[i].value = QString("dfv");
    namelist[i++].name = QString("OUTDFV");
    namelist[i].value = QString("fph");
    namelist[i++].name = QString("OUTFPH");
    namelist[i].value = QString("sph");
    namelist[i++].name = QString("OUTSPH");
    namelist[i].value = QString("entr");
    namelist[i++].name = QString("OUTENTRAINMENT");
    namelist[i].value = QString("ts");
    namelist[i++].name = QString("OUTTIMESTEP");


    // input maps start here !!!
    mapstartnr = i;
    int j = mapstartnr;
    for (i = 0; i < DEFmaps.count(); i++)
    {
        QStringList SL;
        SL = DEFmaps[i].split(";",QString::SkipEmptyParts);

        if (SL[0] == "0")
        {
            namelist[j].name = QString("");
            j++;
            namelist[j].name = "[" + SL[1] + "]";
            j++;
        }
        else
            if (SL[0] == "1")
            {
                namelist[j].name = "[" + SL[1] + "]";
                j++;
            }
            else
            {
                namelist[j].name = SL[4];
                namelist[j].value = SL[2];
                j++;
            }
    }
    nrnamelist = j;

    //   for (j = 0; j < nrnamelist; j++)
    //   qDebug() << namelist[j].name << "=" << namelist[j].value;


}
//---------------------------------------------------------------------------



