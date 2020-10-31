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
    DEFmaps.append("2;Snowmelt ID;snowid.map;Snowmelt zone ID number for snowmelt file starting with 1 (0 is non-snow area);SnowID");

    DEFmaps.append("0;Catchment");
    DEFmaps.append("2;DEM;dem.map;Digital elevation model (m);dem");
    DEFmaps.append("2;Buffers;buffers.map;Flood bariers and obstacles (taluts, dikes, in m);buffers");
    DEFmaps.append("2;Gradient;grad.map;Sine of slope gradient in direction of flow;grad");
    DEFmaps.append("2;LDD;ldd.map;Local surface Drainage Direction network;ldd");
    DEFmaps.append("2;Outlet;outlet.map;Main catchment outlet corresponding to LDD map;outlet");
    DEFmaps.append("2;Points;outpoint.map;Reporting points for hydrograph/sedigraph (1,2,3,...);outpoint");

    DEFmaps.append("0;Landuse");
    DEFmaps.append("2;Units;landunit.map;Classified land unit map (integers 0-n) for output of erosion values;landunit");
    DEFmaps.append("2;Cover;per.map;Fraction surface cover by vegetation and residue (-);cover");
    DEFmaps.append("2;Litter;litter.map;Fraction of surface cover by litter/herbs under trees (-);litter");
    DEFmaps.append("2;LAI;lai.map;Leaf area index of the plant cover in a gridcell (m2/m2);lai");
    DEFmaps.append("2;Height;ch.map;Plant height (m);ch");
    DEFmaps.append("2;Road width;roadwidt.map;Width of impermeable roads (m);road");
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
    DEFmaps.append("2;D90;d90.map;90% quartile of the texture of the suspendeed matter (mu);d90");
    DEFmaps.append("2;Max material;detmat.map;Detacheable material per square meter (kg/m2) (-1 = infinite);detmat");
    DEFmaps.append("2;MixingDepth;sedmixdepth.map; Mixing depth for deposited sediment (m);sedmixdepth");
    //DEFmaps.append("2;MaxDepth;maxdetdepth.map; Maximum depth for detachment (m)(-1 = infinite);maxdet");

    DEFmaps.append("0;Infiltration");
    DEFmaps.append("1;Swatre");
    DEFmaps.append("2;Profile soil;profile.map;ID numbers corresponding to land units in profile table;profmap");
    DEFmaps.append("2;Prof. Crust;profcrst.map;ID numbers of crusted soils (using also profile table);profcrst");
    DEFmaps.append("2;Prof. Wheel;profwltr.map;ID numbers of compacted wheel tracks (using also profile table);profwltr");
    DEFmaps.append("2;Prof. Grass;profgras.map;ID numbers of grasstrips (using also profile table);profgras");
    DEFmaps.append("2;Initial suction;inithead;initial matrix potential (cm) of layers 001 to nnn (filename witout extension);inithead");
    DEFmaps.append("2;Swatre Output points;swatreoutput.map;Points for swatre profile output 1-n);swatreout");
    DEFmaps.append("2;Repellency;repel.map;Gridcells included in water repellency (1/0);repelcell");

    DEFmaps.append("1;1st layer Green&Ampt/Smith&Parlange");

  //  DEFmaps.append("2;thetasub;thetasub.map;Minimum moisture content at the bottom of the profile, below which there is no Percolation (-);thetasub");
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

    DEFmaps.append("1;Special surfaces");
    DEFmaps.append("2;Ksat Crust;ksatcrst.map;Ksat of crusts (all models except SWATRE) (mm/h);ksatcrst");
    DEFmaps.append("2;Porosity Crust;porecrust.map;Porosity of crusted areas (all models except SWATRE) (-);porecrst");
    DEFmaps.append("2;Ksat Compact;ksatcomp.map;Ksat of compacted areas (all models except SWATRE) (mm/h);ksatcomp");
    DEFmaps.append("2;Porosity Compact;porecomp.map;Porosity of compacted areas (all models except SWATRE) (-);porecomp");
//    DEFmaps.append("2;Initial moisture Crust;thetacrust.map;Initial moisture content of crusted areas (all models except SWATRE) (-);thetacrst");
//    DEFmaps.append("2;Initial moisture Compact;thetacomp.map;Initial moisture content of compacted areas (all models except SWATRE) (-);thetacomp");
//    DEFmaps.append("2;Initial moisture Grass;thetagrass.map;Initial moisture content of grasstrips (all models except SWATRE) (-);thetagras");

    DEFmaps.append("0;Channels");
    DEFmaps.append("2;LDD;lddchan.map;LDD of main channel (must be 1 branch connected to the outlet);lddchan");
    DEFmaps.append("2;Width;chanwidt.map;Channel width (m);chanwidth");
    DEFmaps.append("2;Depth;chandepth.map;Channel depth, zero (0) depth is considered infinite (m);chandepth");
    DEFmaps.append("2;Side angle;chanside.map;Channel side angle (tan angle  channel side and surface: 0 is rectangular);chanside");
    DEFmaps.append("2;Gradient;changrad.map;Slope gradient of channel bed (-);changrad");
    DEFmaps.append("2;N;chanman.map;Mannings n of channel bed (-);chanman");
    DEFmaps.append("2;Ksat;chanksat.map;Infiltration rate of channel bed (mm/h);chanksat");
    DEFmaps.append("2;BaseFlow;baseflow.map; base flow discharges (m3/s);baseflow");
    DEFmaps.append("2;Cohesion;chancoh.map;Cohesion of channel bed (kPa);chancoh");
    DEFmaps.append("2;Channelmaterial;chandetmat.map;Detacheable material per square meter (kg/m2) (-1 = infinite);chandetmat");
    DEFmaps.append("2;ChannelMixingDepth;chansedmixdepth.map; Mixing depth for deposited sediment in channel (m);chansedmixdepth");
    DEFmaps.append("2;ChannelMaxQ;chanmaxq.map;Maximum limiting channel discharge, e.g. in culverts (m3/s);chanmaxq");
//    DEFmaps.append("2;ChannelLevee;chanlevee.map;Height of small channel levee on both sides of the channel (m);chanlevee");
    DEFmaps.append("2;hmxInit;WHinit.map;Initial floodlevel (m);hmxinit");
//    DEFmaps.append("2;floodZone;floodzone.map;potential flood zone to limit calculations (1 = in, 0 = out);floodzone");

    DEFmaps.append("0;Conservation");
    DEFmaps.append("2;Grass strips;grasswid.map;Width of grass strips (m);grasswidth");
    DEFmaps.append("2;Sediment traps;sedretmax.map;Max sediment volume in m2 per cell that can be trapped;sedretmax");
    DEFmaps.append("2;Ksat Grass;ksatgras.map;Ksat of grassstrips (all models except SWATRE) (mm/h);ksatgras");
    DEFmaps.append("2;Porosity Grass;poregras.map;Porosity of grasstrips (all models except SWATRE) (-);poregras");
    DEFmaps.append("2;Cohesion Grass;cohgras.map;Porosity of grasstrips (all models except SWATRE) (-);cohgras");

//    DEFmaps.append("0;Snowmelt");
//    DEFmaps.append("2;Snowmelt ID;snowid.map;Snowmelt zone ID number for snowmelt file starting with 1 (0 is non-snow area);SnowID");

    DEFmaps.append("0;Storm/Tile drains");
    DEFmaps.append("2;LDD;lddtile.map;LDD of tile drain system (must be one system connected to the outlet);lddtile");
    DEFmaps.append("2;Sink;tilesink.map;Sink holes connecting surface to tile drain system (size in m2);tilesink");
    DEFmaps.append("2;Diameter;tilediameter.map;Tile drain pipe diameter (m);tilediameter");
    DEFmaps.append("2;Width;tilewidth.map;Tile drain pipe width, total in cell if more than one drain (m);tilewidth");
    DEFmaps.append("2;Height;tileheight.map;Tile drain pipe height (m);tileheight");
    DEFmaps.append("2;Depth;tiledepth.map;Tile drain pipe depth below surface (m);tiledepth");
    DEFmaps.append("2;Gradient;tilegrad.map;Slope gradient of the tile drains (-);tilegrad");
    DEFmaps.append("2;N;tileman.map;Mannings n of the tile drains (-);tileman");
    //houses
    DEFmaps.append("0;Houses");
    DEFmaps.append("2;House Cover;housecover.map;Fraction of hard roof surface per cell (-);housecover");
    DEFmaps.append("2;Roof Storage;roofstore.map;Size of interception storage of rainwater on roofs (mm);roofstore");
    DEFmaps.append("2;Drum Store;drumstore.map;Size of storage of rainwater drums (m3);drumstore");

    DEFmaps.append("0;FlowBarriers");
    DEFmaps.append("2;FlowBarrierIndex;flowbarrierindex.map;An index value, indicating which flow barrier properties will be used (-);flowbarrierindex");
    DEFmaps.append("2;FlowBoundary;flowboundary.map;A value of 1 at the domain boundary means free outflow, a 0 means no flow (-);flowboundary");

    // example
    //   DEFmaps.append("0;Pesticides");
    //   DEFmaps.append("2;Pest Initial;pestinit.map;Inital content bla bla;pestini");

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

    i = 0;
    namelist[i++].name = QString("[openLISEM runfile version 6.0]");
    namelist[i++].name = QString("");
    //###
    namelist[i++].name = QString("[Input]");
    namelist[i++].name = QString("Work Directory");
    namelist[i++].name = QString("Map Directory");
    namelist[i++].name = QString("Rainfall Directory");
    namelist[i++].name = QString("Rainfall file");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include Snowmelt");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Result datetime");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include Satellite Image");
    namelist[i++].name = QString("satImage Directory");
    namelist[i++].name = QString("satImage file");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Advanced Options");
    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Output]");
    namelist[i++].name = QString("Result Directory");
    namelist[i].value = QString("totals.csv");
    namelist[i++].name = QString("Main results file");
    namelist[i].value = QString("totalseries.csv");
    namelist[i++].name = QString("Total Series file");
    namelist[i].value = QString("hydrograph.csv");
    namelist[i++].name = QString("Filename point output");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Report point output separate");
//    namelist[i].value = QString("1");
//    namelist[i++].name = QString("Timeplot as CSV");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Timeplot as PCRaster");
    namelist[i].value = QString("3");
    namelist[i++].name = QString("Report digits out");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Report format GTiff");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("End run report");
    namelist[i].value = QString("rainfall.map");
    namelist[i++].name = QString("Rainfall map");
    namelist[i].value = QString("interception.map");
    namelist[i++].name = QString("Interception map");
    namelist[i].value = QString("infiltration.map");
    namelist[i++].name = QString("Infiltration map");
    namelist[i].value = QString("runoff.map");
    namelist[i++].name = QString("Runoff map");
    namelist[i].value = QString("WHmax.map");
    namelist[i++].name = QString("WH max level map");
    namelist[i].value = QString("chandism3.map");
    namelist[i++].name = QString("Channel discharge map");
    namelist[i].value = QString("chanmaxq.map");
    namelist[i++].name = QString("Channel Max Q");
    namelist[i].value = QString("chanmaxwh.map");
    namelist[i++].name = QString("Channel Max WH");
    namelist[i].value = QString("Vmax.map");
    namelist[i++].name = QString("Max Velocity");
    namelist[i].value = QString("VHmax.map");
    namelist[i++].name = QString("Max Momentum");
    namelist[i].value = QString("floodtime.map");
    namelist[i++].name = QString("Flood time map");
    namelist[i].value = QString("floodstart.map");
    namelist[i++].name = QString("Flood start time");
    namelist[i].value = QString("stormdrain.map");
    namelist[i++].name = QString("Storm Drain map");
    namelist[i].value = QString("stormdrainvol.map");
    namelist[i++].name = QString("Storm Drain Vol map");
    namelist[i].value = QString("eros.map");
    namelist[i++].name = QString("Erosion map");
    namelist[i].value = QString("depo.map");
    namelist[i++].name = QString("Deposition map");
    namelist[i].value = QString("soilloss.map");
    namelist[i++].name = QString("Soilloss map");
    namelist[i].value = QString("chandet.map");
    namelist[i++].name = QString("Channel detachment map");
    namelist[i].value = QString("chandep.map");
    namelist[i++].name = QString("Channel deposition map");
    namelist[i].value = QString("totlandunit.txt");
    namelist[i++].name = QString("Filename landunit output");
    namelist[i].value = QString("floodstats.csv");
    namelist[i++].name = QString("Flood stats");

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
    namelist[i].value = QString("2");
    namelist[i++].name = QString("Routing Kin Wave 2D");
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
    namelist[i++].name = QString("Include channel culverts");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Hard Surfaces");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Include road system");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include house storage");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include raindrum storage");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include storm drains");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Nr user Cores");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Interception]");
//    namelist[i].value = QString("0");
//    namelist[i++].name = QString("Use canopy storage map");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Canopy storage equation");
//    namelist[i].value = QString("0.05");
//    namelist[i++].name = QString("Stemflow fraction");
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
    namelist[i].value = QString("1");
    // obsolete
//    namelist[i++].name = QString("Include percolation");
//    namelist[i].value = QString("0");
//    namelist[i++].name = QString("Subsoil drainage");
    namelist[i].value = QString("c:\\");
    namelist[i++].name = QString("Table Directory");
    namelist[i].value = QString("profile.inp");
    namelist[i++].name = QString("Table File");
    //namelist[i].value = QString("0.01");
    //namelist[i++].name = QString("SWATRE internal minimum timestep");
    namelist[i].value = QString("inithead");
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
    namelist[i++].name = QString("[Flow]");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Flow Boundary 2D");
    namelist[i].value = QString("0.05");
    namelist[i++].name = QString("Minimum reported flood height");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Include buffers");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include flow barriers");
    namelist[i].value = QString("flowbarriers.txt");
    namelist[i++].name = QString("Flow barrier table filename");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Flood initial level map");
    namelist[i].value = QString("0.2");
    namelist[i++].name = QString("Flooding courant factor");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Timestep flood");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Use time avg V");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Flood solution");
    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Advanced]");
    namelist[i].value = QString("4"); //HLL2
    namelist[i++].name = QString("Flooding SWOF Reconstruction");
    namelist[i].value = QString("1"); //minmod
    namelist[i++].name = QString("Flooding SWOF flux limiter");
    namelist[i].value = QString("200");
    namelist[i++].name = QString("Flood Max Iterations");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Flooding mixing coefficient");
    namelist[i].value = QString("2.0");
    namelist[i++].name = QString("Flooding runoff partitioning");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Use gravity flow");
    namelist[i].value = QString("0.02");
    namelist[i++].name = QString("Angle flow to channel");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Use fixed angle");
    namelist[i].value = QString("10.0");
    namelist[i++].name = QString("Pit Value");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Variable Timestep");
//    namelist[i].value = QString("0");
//    namelist[i++].name = QString("Minimum Timestep Method");

//    namelist[i].value = QString("0");
//    namelist[i++].name = QString("Use Heun");

    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Erosion]");
    namelist[i].value = QString("1;28.300;0.520;0.042");
    namelist[i++].name = QString("KE parameters EQ1");
    namelist[i].value = QString("0;8.950;8.440");
    namelist[i++].name = QString("KE parameters EQ2");
    namelist[i].value = QString("0;7.600;0.220");
    namelist[i++].name = QString("KE parameters EQ3");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("KE time based");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Detachment efficiency");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Settling Velocity");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("No detachment boundary");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Splash Delivery Ratio");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Use material depth");
    namelist[i].value = QString("1500.0");
    namelist[i++].name = QString("Sediment bulk density");
    namelist[i].value = QString("0.5");
    namelist[i++].name = QString("Particle Cohesion of Deposited Layer");
    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Sediment]");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Use 2 phase flow");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("River BL method");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("River SS method");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Include River diffusion");
    namelist[i].value = QString("0.5");
    namelist[i++].name = QString("River Sigma diffusion");
    namelist[i].value = QString("0"); //govers
    namelist[i++].name = QString("Flooding SS method");
    namelist[i].value = QString("1"); //van rijn simplified
    namelist[i++].name = QString("Flooding BL method");
    namelist[i].value = QString("1");
    namelist[i++].name = QString("Include diffusion");
    namelist[i].value = QString("0.5");
    namelist[i++].name = QString("Sigma diffusion");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Use grain size distribution");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Estimate grain size distribution");
    namelist[i].value = QString("6");
    namelist[i++].name = QString("Number of grain size classes (simulated)");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Read grain distribution maps");
    namelist[i].value = QString("2;20;50;125;250;500");
    namelist[i++].name = QString("Grain size class maps");
    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Conservation]");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include grass strips");
    namelist[i].value = QString("0.1");
    namelist[i++].name = QString("Grassstrip Mannings n");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("Include Sediment traps");
    namelist[i].value = QString("0.8");
    namelist[i++].name = QString("Sediment Trap Mannings n");
    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Calibration]");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Grain Size calibration");
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
    namelist[i++].name = QString("Cohesion calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Cohesion Channel calibration");
    namelist[i].value = QString("1.0");
    namelist[i++].name = QString("Aggregate stability calibration");
    //###
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[Output maps]");
    namelist[i++].name = QString("Runoff maps in l/s/m");
  //  namelist[i++].name = QString("Timeseries as PCRaster");
    namelist[i++].name = QString("Erosion map units (0/1/2)");
    namelist[i++].name = QString("Regular runoff output");
    namelist[i++].name = QString("Output interval");
    namelist[i++].name = QString("User defined output");
    namelist[i++].name = QString("Output times");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutRunoff");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutWH");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutV");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutInterception");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutSurfStor");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutInf");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutTileDrain");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutTileV");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutDet");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutDep");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutTC");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutConc");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutSed");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutSL");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutSedSS");
    namelist[i].value = QString("0");
    namelist[i++].name = QString("OutSedBL");
    // output maps have standard names
    // input maps names are defined in DEFmaps
    namelist[i++].name = QString("");
    namelist[i++].name = QString("[map names]");
    // input maps start here !!!
    mapstartnr = i;
    int j = mapstartnr;
    for (i = 0; i < DEFmaps.count(); i++)
    {
        QStringList SL;
        SL = DEFmaps[i].split(";",Qt::SkipEmptyParts);

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
    for (int i = 0; i < nrnamelist; i++)
        namelist[i].gotit = false;
}
//---------------------------------------------------------------------------



