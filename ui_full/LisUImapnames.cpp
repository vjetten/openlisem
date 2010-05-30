/*
 * lisqtmapnames.cpp
 *
 *  Created on: Jul 10, 2009
 *      Author: jetten
 */
#include "lisemqt.h"

//--------------------------------------------------------------------
void lisemqt::DefaultMapnames()
{
    DEFmaps.append("0;Catchment;0");
    DEFmaps.append("2;Gradient;grad.map;Sine of slope gradient in direction of flow;grad;2");
    DEFmaps.append("2;LDD;ldd.map;Local surface Drainage Direction network;ldd;2");
    DEFmaps.append("2;Outlet;outlet.map;Main catchment outlet corresponding to LDD map;outlet;2");
    DEFmaps.append("2;Points;outpoint.map;Reporting points for hydrograph/sedigraph (1 to nr);outlet;2");
    DEFmaps.append("2;ID;ID.map;Raingauge zone ID numbers;id;2");

    DEFmaps.append("0;Landuse;0");
    DEFmaps.append("2;Cover;per.map;Fraction surface cover by vegetation and residue;cover;2");
    DEFmaps.append("2;LAI;lai.map;Leaf area index of the plant cover in a gridcell (m2/m2);lai;2");
    DEFmaps.append("2;Height;ch.map;Plant height (m);ch;2");
    DEFmaps.append("2;Road width;roadwidt.map;Width of impermeable roads (m);road;2");
    DEFmaps.append("2;Grass strips;grasswid.map;Width of grass strips (m);grasswidth;2");
    DEFmaps.append("2;Hard surf;hardsurf.map;Impermeable and unerodable surfaces (value 1);hardsurf;2");

    DEFmaps.append("0;Surface;0");
    DEFmaps.append("2;RR;rr.map;Random Roughness (here standard deviation of heights) (cm);rr;2");
    DEFmaps.append("2;n;n.map;Manning's n (-);manning;2");
    DEFmaps.append("2;Crust;crustfrc.map;Fraction of gridcell covered with Crust (-);crustfrc;2");
    DEFmaps.append("2;Compacted;compfrc.map;Fraction of gridcell compacted (e.g. wheeltracks)(-);compfrc;2");
    DEFmaps.append("2;Stoniness;stonefrc.map;Fraction of gridcell covered by stones (-);stonefrc;2");

    DEFmaps.append("0;Erosion;0");
    DEFmaps.append("2;Cohesion;coh.map;Cohesion (kPa);coh;2");
    DEFmaps.append("2;Cohesion;cohadd.map;Extra cohesion factor by e.g. plant roots (kPa);cohadd;2");
    DEFmaps.append("2;Aggregates;aggrstab.map;Aggregate stability for splash erosion (-);aggrstab;2");
    DEFmaps.append("2;D50;d50.map;Median of the texture of the suspendeed matter (mu);d50;2");
    DEFmaps.append("2;Hard Surface;hardsurf.map;non-erodible (1) and erodible surface (0);hardsurf;2");

    DEFmaps.append("0;Infiltration;0");
    DEFmaps.append("1;Infil Swatre;0");
    DEFmaps.append("2;Profile table;profile.inp;Text file with land unit descriptions and links to tables;profinp;2");
    DEFmaps.append("2;Table Directory; ;Directory where SWATRE tables are stored;TableDir;2");
    DEFmaps.append("2;Internal timestep;5;Timestep in sec of internal SWATRE iteration (smaller than LISEM timestep);SwatreDTSEC;2");
    DEFmaps.append("2;Profile;profile.map;ID numbers corresponding to land units in profile table;profmap;2");
    DEFmaps.append("2;Prof. Crust;profcrst.map;ID numbers of crusted soils (using also profile table);profcrst;2");
    DEFmaps.append("2;Prof. Wheel;profwltr.map;ID numbers of compacted wheel tracks (using also profile table);profwltr;2");
    DEFmaps.append("2;Prof. Grass;profgras.map;ID numbers of grasstrips (using also profile table);profgras;2");
    DEFmaps.append("2;Init. suction;inithead;initial matrix potential of layers 001 to nnn (filename witout extension)(cm);inithead;2");
    DEFmaps.append("2;Output head;headout.map;Locations to write tables of the matrix potential;headout;2");
    DEFmaps.append("2;Switch output head;1;write matrix head in 4 locations to text files;Switchheadout;2");
    DEFmaps.append("2;Swicth geometric avg.;0;Geometric or arithmetric average Ksat between nodes;SwitchGeometric;2");

    DEFmaps.append("1;Infil G&A 1Layer;0");
    DEFmaps.append("2;Ksat1;ksat1.map;Layer 1: Saturated Hydraulic Conductivity (mm/h);ksat1;2");
    DEFmaps.append("2;Psi1;psi1.map;Layer 1: Average suction at the wetting front (cm);psi1;2");
    DEFmaps.append("2;Thetas1;thetas1.map;Layer 1: Porosity (-);thetas1;2");
    DEFmaps.append("2;Thetai1;thetai1.map;Layer 1: Initial moisture content (-);thetai1;2");
    DEFmaps.append("2;Depth1;soildep1.map;Layer 1: Depth (mm) to bottom of layer 1;soildep1;2");

    DEFmaps.append("1;Infil G&A 2Layers;0");
    DEFmaps.append("2;Ksat2;ksat2.map;Layer 2: Saturated Hydraulic Conductivity (mm/h);ksat2;2");
    DEFmaps.append("2;Psi2;psi2.map;Layer 2: Average suction at the wetting front (cm);psi2;2");
    DEFmaps.append("2;Thetas2;thetas2.map;Layer 2: Porosity (-);thetas2;2");
    DEFmaps.append("2;Thetai2;thetai2.map;Layer 2: Initial moisture content (-);thetai2;2");
    DEFmaps.append("2;Depth2;soildep2.map;Layer 2: Depth (mm) to bottom of layer 2;soildep2;2");

    DEFmaps.append("1;Infil Holtan;0");
    DEFmaps.append("2;A;a.map;Parameter A in Holtan equation;A;2");
    DEFmaps.append("2;FP;fp.map;Parameter FP in Holtan equation;FP;2");
    DEFmaps.append("2;P;p.map;Parameter P in Holtan equation;P;2");

    DEFmaps.append("1;Infil Ksat;0");
    DEFmaps.append("2;Ksat1;ksat1.map;Saturated Hydraulic Conductivity (mm/h);ksat1);2");

    DEFmaps.append("1;Infil special;2");
    DEFmaps.append("2;Ksat Crust;ksatcrst.map;Ksat of crusts (all models except SWATRE) (mm/h);ksatcrst;2");
    DEFmaps.append("2;Ksat Compact;ksatcomp.map;Ksat of compacted areas (all models except SWATRE) (mm/h);ksatcomp;2");
    DEFmaps.append("2;Ksat Grass;ksatgras.map;Ksat of grassstrips (all models except SWATRE) (mm/h);ksatgras;2");
    //	DEFmaps.append("2;Drain. fact.;drfactor.map;Drainage exponent in k=ks*(moist/pore)^d (-);drfactor;2");

    DEFmaps.append("0;Channels;0");
    DEFmaps.append("1;Channel properties;0");
    DEFmaps.append("2;LDD;lddchan.map;LDD of main channel (must be 1 branch connected to the outlet);lddchan;2");
    DEFmaps.append("2;Width;chanwidt.map;Channel width (m);chanwidth;2");
    DEFmaps.append("2;Side angle;chanside.map;Channel side angle (tan angle  channel side and surface: 0 is rectangular);chanside;2");
    DEFmaps.append("2;Gradient;changrad.map;Slope gradient of channel bed (-);changrad;2");
    DEFmaps.append("2;N;chanman.map;Mannings n of channel bed (-);chanman;2");
    DEFmaps.append("2;Cohesion;chancoh.map;Cohesion of channel bed (kPa);chancoh;2");

    DEFmaps.append("1;Channelinfil;0");
    DEFmaps.append("2;Ksat;chanksat.map;Infiltration rate of channel bed (mm/h);chanksat;2");

    DEFmaps.append("1;ChannelBaseflow;0");
    DEFmaps.append("2;Inflow flux;chanbaseflux.map;Incoming flux into channel from the two sides (m3/s);chanbaseflux;2");
    DEFmaps.append("2;Increase in baseflow;chanincrease.map;Increase in basevolume during rainstorm (-);chanincrease;2");
    DEFmaps.append("2;Initial volume;chanvini.map;Initial baseflow water volume in channel (m3);chanvolini;2");

    DEFmaps.append("0;Buffers;0");
    DEFmaps.append("2;Buffer ID nr;bufferid.map;ID number for each buffer starting with 1 (0 is outside area);bufferID;2");
    DEFmaps.append("2;Buffer volume;buffervol.map;Buffer volumes at the locations of the buffers (m3);bufferVolume;2");
    DEFmaps.append("2;Buffer area;bufferarea.map;Buffer area at locations of the buffers (m2);bufferarea;2");

    DEFmaps.append("0;Snowmelt;0");
    DEFmaps.append("2;Snowmelt ID;snowid.map;Snowmelt zone ID number for snowmelt file starting with 1 (0 is non-snow area);SnowID;2");

    DEFmaps.append("0;Wheeltracks;0");
    DEFmaps.append("2;LDD;lddwheel.map;LDD of wheeltrack network (can be separate branches with pits);lddwheel;2");
    DEFmaps.append("2;Number;wheelnbr.map;Number of wheeltrack channels in a gridcell (-);wheelnbr;2");
    DEFmaps.append("2;Width;wheelwid.map;Sum of widths of wheeltracks in a gridcell (m);wheelwidth;2");
    DEFmaps.append("2;Depth;wheeldep.map;Wheel track overflow depth (cm);wheeldepth;2");
    DEFmaps.append("2;Gradient;wheelgrd.map;DEFmapsope gradient of wheel tracks (-);wheelgradient;2");
    DEFmaps.append("2;N;wheelman.map;Mannings n of Wheel tracks (-);wheelman;2");
    DEFmaps.append("2;Cohesion;wheelcoh.map;Cohesion of wheel tracks (kPa);wheelcohesion;2");
    DEFmaps.append("2;Ksat;ksatwt.map;Saturated hydraulic conductivity of wheel tracks (mm/h);ksatwt;2");

    DEFmaps.append("0;Texture classes;0");
    DEFmaps.append("2;Class 0;mu0.map;Clay fraction (MUST BE CLAY < 2mu);fractionmu0;2");
    DEFmaps.append("2;Class 1;mu1.map;Soil texture fraction for class 1 (-);fractionmu1;2");
    DEFmaps.append("2;Class 2;mu2.map;Soil texture fraction for class 2 (-);fractionmu2;2");
    DEFmaps.append("2;Class 3;mu3.map;Soil texture fraction for class 3 (-);fractionmu3;2");
    DEFmaps.append("2;Class 4;mu4.map;Soil texture fraction for class 4 (-);fractionmu4;2");
    DEFmaps.append("2;Class 5;mu5.map;Soil texture fraction for class 5 (-);fractionmu5;2");

    DEFmaps.append("0;Nutrients;0");
    DEFmaps.append("2;Bulk Dens.;bulkdens.map;Bulk density of the topsoil (kg/m3);bulk;2");
    DEFmaps.append("2;P Content;pcont.map;Phosphate (P) content of the soil (kg/kg);pcont;2");
    DEFmaps.append("2;P Solute;psolut.map;Initial solute store P in surface layer (kg/m2);psolute;2");
    DEFmaps.append("2;P Efficiency;peff.map;Extraction efficiency (s-1);pefficiency;2");
    DEFmaps.append("2;P Sorption;Psorp.map;Sorption isotherm kd (m3/kg);psorp;2");
    DEFmaps.append("2;P Conversion;Pconv.map;Conversion P from soil content to clay content(-);pconv;2");
    DEFmaps.append("2;NH4 Content;nh4cont.map;Ammonium (NH4+) content of the soil (kg/kg);nh4cont;2");
    DEFmaps.append("2;NH4 Solute;nh4solut.map;Initial solute store NH4 in surface layer (kg/m2);nh4solute;2");
    DEFmaps.append("2;NH4 Efficiency;nh4eff.map;Extraction efficiency (s-1);nh4efficiency;2");
    DEFmaps.append("2;NH4 Sorption;NH4sorp.map;Sorption isotherm kd (m3/kg);nh4sorp;2");
    DEFmaps.append("2;NH4 Conversion;NH4conv.map;Conversion NH4 from soil content to clay content(-);nh4conv;2");
    DEFmaps.append("2;NO3 Content;NO3cont.map;Nitrate (NO3-) content of the soil (kg/kg);no3cont;2");
    DEFmaps.append("2;NO3 Solute;NO3solut.map;Initial solute store NO3 in surface layer (kg/m2);no3solute;2");
    DEFmaps.append("2;NO3 Efficiency;NO3eff.map;Extraction efficiency (s-1);no3efficiency;2");
    DEFmaps.append("2;NO3 Sorption;NO3sorp.map;Sorption isotherm kd (m3/kg);no3sorp;2");
    DEFmaps.append("2;NO3 Conversion;NO3conv.map;Conversion NO3 from soil content to clay content(-);no3conv;2");

    DEFmaps.append("0;Gullies;0");
    DEFmaps.append("2;DEM;dem.map;Digital elevation model (m);dem;2");
    DEFmaps.append("2;mannings N;gullyman.map;manning's n gully bottom (-);gullyn;2");
    DEFmaps.append("2;BulkDensity;bulkdens.map;Bulkdensity of topsoil (kg/m3);bulkdens1;2");
    DEFmaps.append("2;Ksat;ksat1.map;Ksat of topsoil for gully infil (mm/h);gulksat1;2");
    DEFmaps.append("2;Depth layer 2;soilDep2.map;Depth to subsoil (cm);gullydep;2");
    DEFmaps.append("2;Cohesion layer 2;coh2.map;Cohesion of subsoil (kPa);gullycoh;2");
    DEFmaps.append("2;BulkDensity 2;bulkden2.map;Bulkdensity of subsoil (kg/m3);bulkdens2;2");
    DEFmaps.append("2;Ksat 2;gulksat2.map;Ksat of subsoil for gully infil (mm/h);gulksat2;2");
    DEFmaps.append("2;Excluded areas;noncrit.map;areas to be excluded (1) and rest (0);nonfcrit;2");
    DEFmaps.append("2;Gully initial Width;gulwinit.map; initial gully width (m);gulwinit;2");
    DEFmaps.append("2;Gully initial Depth;guldinit.map; initial gully depth (m);guldinit;2");
	/*
    DEFmaps.append("0,Catchment          ,0");
    DEFmaps.append("2,Gradient,grad.map,Sine of slope gradient in direction of flow,grad                                           ,2");
    DEFmaps.append("2,LDD,ldd.map,Local surface Drainage Direction network,ldd                                                     ,2");
    DEFmaps.append("2,Outlet,outlet.map,Main catchment outlet corresponding to LDD map,outlet                                      ,2");
    DEFmaps.append("2,Points,outpoint.map,Reporting points for hydrograph/sedigraph (1 to nr),outlet                               ,2");
    DEFmaps.append("2,ID,ID.map,Raingauge zone ID numbers,id                                                                       ,2");

    DEFmaps.append("0,Landuse          ,0");
    DEFmaps.append("2,Cover,per.map,Fraction surface cover by vegetation and residue,cover                                         ,2");
    DEFmaps.append("2,LAI,lai.map,Leaf area index of the plant cover in a gridcell (m2/m2),lai                                     ,2");
    DEFmaps.append("2,Height,ch.map,Plant height (m),ch                                                                            ,2");
    DEFmaps.append("2,Road width,roadwidt.map,Width of impermeable roads (m),road                                                  ,2");
    DEFmaps.append("2,Grass strips,grasswid.map,Width of grass strips (m),grasswidth                                               ,2");
    DEFmaps.append("2,Hard surf,hardsurf.map,Impermeable and unerodable surfaces (value 1),hardsurf                                ,2");

    DEFmaps.append("0,Surface          ,0");
    DEFmaps.append("2,RR,rr.map,Random Roughness (here standard deviation of heights) (cm),rr                                      ,2");
    DEFmaps.append("2,n,n.map,Manning's n (-),manning                                                                              ,2");
    DEFmaps.append("2,Crust,crustfrc.map,Fraction of gridcell covered with Crust (-),crustfrc                                      ,2");
    DEFmaps.append("2,Compacted,compfrc.map,Fraction of gridcell compacted (e.g. wheeltracks)(-),compfrc                           ,2");
    DEFmaps.append("2,Stoniness,stonefrc.map,Fraction of gridcell covered by stones (-),stonefrc                                   ,2");

    DEFmaps.append("0,Erosion          ,0");
    DEFmaps.append("2,Cohesion,coh.map,Cohesion (kPa),coh                                                                          ,2");
    DEFmaps.append("2,Cohesion,cohadd.map,Extra cohesion factor by e.g. plant roots (kPa),cohadd                                   ,2");
    DEFmaps.append("2,Aggregates,aggrstab.map,Aggregate stability for splash erosion (-),aggrstab                                  ,2");
    DEFmaps.append("2,D50,d50.map,Median of the texture of the suspendeed matter (mu),d50                                          ,2");
    DEFmaps.append("2,Hard Surface,hardsurf.map,non-erodible (1) and erodible surface (0),hardsurf                                 ,2");

    DEFmaps.append("0,Infiltration          ,0");
    DEFmaps.append("1,Infil Swatre          ,0");
    DEFmaps.append("2,Profile table,profile.inp,Text file with land unit descriptions and links to tables,profinp                  ,2");
    DEFmaps.append("2,Table Directory, ,Directory where SWATRE tables are stored,TableDir                                          ,2");
    DEFmaps.append("2,Internal timestep,5,Timestep in sec of internal SWATRE iteration (smaller than LISEM timestep),SwatreDTSEC           ,2");
    DEFmaps.append("2,Profile,profile.map,ID numbers corresponding to land units in profile table,profmap                          ,2");
    DEFmaps.append("2,Prof. Crust,profcrst.map,ID numbers of crusted soils (using also profile table),profcrst                     ,2");
    DEFmaps.append("2,Prof. Wheel,profwltr.map,ID numbers of compacted wheel tracks (using also profile table),profwltr            ,2");
    DEFmaps.append("2,Prof. Grass,profgras.map,ID numbers of grasstrips (using also profile table),profgras                        ,2");
    DEFmaps.append("2,Init. suction,inithead,initial matrix potential of layers 001 to nnn (filename witout extension)(cm),inithead,2");
    DEFmaps.append("2,Output head,headout.map,Locations to write tables of the matrix potential,headout                            ,2");
    DEFmaps.append("2,Switch output head,1,write matrix head in 4 locations to text files,Switchheadout                            ,2");
    DEFmaps.append("2,Swicth geometric avg.,0,Geometric or arithmetric average Ksat between nodes,SwitchGeometric                  ,2");

    DEFmaps.append("1,Infil G&A 1Layer          ,0");
    DEFmaps.append("2,Ksat1,ksat1.map,Layer 1: Saturated Hydraulic Conductivity (mm/h),ksat1                                       ,2");
    DEFmaps.append("2,Psi1,psi1.map,Layer 1: Average suction at the wetting front (cm),psi1                                        ,2");
    DEFmaps.append("2,Thetas1,thetas1.map,Layer 1: Porosity (-),thetas1                                                            ,2");
    DEFmaps.append("2,Thetai1,thetai1.map,Layer 1: Initial moisture content (-),thetai1                                            ,2");
    DEFmaps.append("2,Depth1,soildep1.map,Layer 1: Depth (mm) to bottom of layer 1,soildep1                                        ,2");

    DEFmaps.append("1,Infil G&A 2Layers          ,0");
    DEFmaps.append("2,Ksat2,ksat2.map,Layer 2: Saturated Hydraulic Conductivity (mm/h),ksat2                                       ,2");
    DEFmaps.append("2,Psi2,psi2.map,Layer 2: Average suction at the wetting front (cm),psi2                                        ,2");
    DEFmaps.append("2,Thetas2,thetas2.map,Layer 2: Porosity (-),thetas2                                                            ,2");
    DEFmaps.append("2,Thetai2,thetai2.map,Layer 2: Initial moisture content (-),thetai2                                            ,2");
    DEFmaps.append("2,Depth2,soildep2.map,Layer 2: Depth (mm) to bottom of layer 2,soildep2                                        ,2");

    DEFmaps.append("1,Infil Holtan          ,0");
    DEFmaps.append("2,A,a.map,Parameter A in Holtan equation,A                                                                     ,2");
    DEFmaps.append("2,FP,fp.map,Parameter FP in Holtan equation,FP                                                                 ,2");
    DEFmaps.append("2,P,p.map,Parameter P in Holtan equation,P                                                                     ,2");

    DEFmaps.append("1,Infil Ksat          ,0");
    DEFmaps.append("2,Ksat1,ksat1.map,Saturated Hydraulic Conductivity (mm/h),ksat1)                                               ,2");

    DEFmaps.append("1,Infil special,2");
    DEFmaps.append("2,Ksat Crust,ksatcrst.map,Ksat of crusts (all models except SWATRE) (mm/h),ksatcrst                            ,2");
    DEFmaps.append("2,Ksat Compact,ksatcomp.map,Ksat of compacted areas (all models except SWATRE) (mm/h),ksatcomp                 ,2");
    DEFmaps.append("2,Ksat Grass,ksatgras.map,Ksat of grassstrips (all models except SWATRE) (mm/h),ksatgras                       ,2");
    //	DEFmaps.append("2,Drain. fact.,drfactor.map,Drainage exponent in k=ks*(moist/pore)^d (-),drfactor                              ,2");

    DEFmaps.append("0,Channels          ,0");
    DEFmaps.append("1,Channel properties          ,0");
    DEFmaps.append("2,LDD,lddchan.map,LDD of main channel (must be 1 branch connected to the outlet),lddchan                       ,2");
    DEFmaps.append("2,Width,chanwidt.map,Channel width (m),chanwidth                                                               ,2");
    DEFmaps.append("2,Side angle,chanside.map,Channel side angle (tan angle  channel side and surface: 0 is rectangular),chanside  ,2");
    DEFmaps.append("2,Gradient,changrad.map,Slope gradient of channel bed (-),changrad                                             ,2");
    DEFmaps.append("2,N,chanman.map,Mannings n of channel bed (-),chanman                                                          ,2");
    DEFmaps.append("2,Cohesion,chancoh.map,Cohesion of channel bed (kPa),chancoh                                                   ,2");

    DEFmaps.append("1,Channelinfil          ,0");
    DEFmaps.append("2,Ksat,chanksat.map,Infiltration rate of channel bed (mm/h),chanksat                                           ,2");

    DEFmaps.append("1,ChannelBaseflow          ,0");
    DEFmaps.append("2,Inflow flux,chanbaseflux.map,Incoming flux into channel from the two sides (m3/s),chanbaseflux               ,2");
    DEFmaps.append("2,Increase in baseflow,chanincrease.map,Increase in basevolume during rainstorm (-),chanincrease               ,2");
    DEFmaps.append("2,Initial volume,chanvini.map,Initial baseflow water volume in channel (m3),chanvolini                         ,2");

    DEFmaps.append("0,Buffers          ,0");
    DEFmaps.append("2,Buffer ID nr,bufferid.map,ID number for each buffer starting with 1 (0 is outside area),bufferID             ,2");
    DEFmaps.append("2,Buffer volume,buffervol.map,Buffer volumes at the locations of the buffers (m3),bufferVolume                 ,2");
    DEFmaps.append("2,Buffer area,bufferarea.map,Buffer area at locations of the buffers (m2),bufferarea                           ,2");

    DEFmaps.append("0,Snowmelt          ,0");
    DEFmaps.append("2,Snowmelt ID,snowid.map,Snowmelt zone ID number for snowmelt file starting with 1 (0 is non-snow area),SnowID ,2");

    DEFmaps.append("0,Wheeltracks          ,0");
    DEFmaps.append("2,LDD,lddwheel.map,LDD of wheeltrack network (can be separate branches with pits),lddwheel                     ,2");
    DEFmaps.append("2,Number,wheelnbr.map,Number of wheeltrack channels in a gridcell (-),wheelnbr                                 ,2");
    DEFmaps.append("2,Width,wheelwid.map,Sum of widths of wheeltracks in a gridcell (m),wheelwidth                                 ,2");
    DEFmaps.append("2,Depth,wheeldep.map,Wheel track overflow depth (cm),wheeldepth                                                ,2");
    DEFmaps.append("2,Gradient,wheelgrd.map,DEFmapsope gradient of wheel tracks (-),wheelgradient                                  ,2");
    DEFmaps.append("2,N,wheelman.map,Mannings n of Wheel tracks (-),wheelman                                                       ,2");
    DEFmaps.append("2,Cohesion,wheelcoh.map,Cohesion of wheel tracks (kPa),wheelcohesion                                           ,2");
    DEFmaps.append("2,Ksat,ksatwt.map,Saturated hydraulic conductivity of wheel tracks (mm/h),ksatwt                               ,2");

    DEFmaps.append("0,Texture classes          ,0");
    DEFmaps.append("2,Class 0,mu0.map,Clay fraction (MUST BE CLAY < 2mu),fractionmu0                                               ,2");
    DEFmaps.append("2,Class 1,mu1.map,Soil texture fraction for class 1 (-),fractionmu1                                            ,2");
    DEFmaps.append("2,Class 2,mu2.map,Soil texture fraction for class 2 (-),fractionmu2                                            ,2");
    DEFmaps.append("2,Class 3,mu3.map,Soil texture fraction for class 3 (-),fractionmu3                                            ,2");
    DEFmaps.append("2,Class 4,mu4.map,Soil texture fraction for class 4 (-),fractionmu4                                            ,2");
    DEFmaps.append("2,Class 5,mu5.map,Soil texture fraction for class 5 (-),fractionmu5                                            ,2");

    DEFmaps.append("0,Nutrients          ,0");
    DEFmaps.append("2,Bulk Dens.,bulkdens.map,Bulk density of the topsoil (kg/m3),bulk                                             ,2");
    DEFmaps.append("2,P Content,pcont.map,Phosphate (P) content of the soil (kg/kg),pcont                                          ,2");
    DEFmaps.append("2,P Solute,psolut.map,Initial solute store P in surface layer (kg/m2),psolute                                  ,2");
    DEFmaps.append("2,P Efficiency,peff.map,Extraction efficiency (s-1),pefficiency                                                ,2");
    DEFmaps.append("2,P Sorption,Psorp.map,Sorption isotherm kd (m3/kg),psorp                                                      ,2");
    DEFmaps.append("2,P Conversion,Pconv.map,Conversion P from soil content to clay content(-),pconv                               ,2");
    DEFmaps.append("2,NH4 Content,nh4cont.map,Ammonium (NH4+) content of the soil (kg/kg),nh4cont                                  ,2");
    DEFmaps.append("2,NH4 Solute,nh4solut.map,Initial solute store NH4 in surface layer (kg/m2),nh4solute                          ,2");
    DEFmaps.append("2,NH4 Efficiency,nh4eff.map,Extraction efficiency (s-1),nh4efficiency                                          ,2");
    DEFmaps.append("2,NH4 Sorption,NH4sorp.map,Sorption isotherm kd (m3/kg),nh4sorp                                                ,2");
    DEFmaps.append("2,NH4 Conversion,NH4conv.map,Conversion NH4 from soil content to clay content(-),nh4conv                       ,2");
    DEFmaps.append("2,NO3 Content,NO3cont.map,Nitrate (NO3-) content of the soil (kg/kg),no3cont                                   ,2");
    DEFmaps.append("2,NO3 Solute,NO3solut.map,Initial solute store NO3 in surface layer (kg/m2),no3solute                          ,2");
    DEFmaps.append("2,NO3 Efficiency,NO3eff.map,Extraction efficiency (s-1),no3efficiency                                          ,2");
    DEFmaps.append("2,NO3 Sorption,NO3sorp.map,Sorption isotherm kd (m3/kg),no3sorp                                                ,2");
    DEFmaps.append("2,NO3 Conversion,NO3conv.map,Conversion NO3 from soil content to clay content(-),no3conv                       ,2");

    DEFmaps.append("0,Gullies          ,0");
    DEFmaps.append("2,DEM,dem.map,Digital elevation model (m),dem                                                                  ,2");
    DEFmaps.append("2,mannings N,gullyman.map,manning's n gully bottom (-),gullyn                                                  ,2");
    DEFmaps.append("2,BulkDensity,bulkdens.map,Bulkdensity of topsoil (kg/m3),bulkdens1                                            ,2");
    DEFmaps.append("2,Ksat,ksat1.map,Ksat of topsoil for gully infil (mm/h),gulksat1                                               ,2");
    DEFmaps.append("2,Depth layer 2,soilDep2.map,Depth to subsoil (cm),gullydep                                                    ,2");
    DEFmaps.append("2,Cohesion layer 2,coh2.map,Cohesion of subsoil (kPa),gullycoh                                                 ,2");
    DEFmaps.append("2,BulkDensity 2,bulkden2.map,Bulkdensity of subsoil (kg/m3),bulkdens2                                          ,2");
    DEFmaps.append("2,Ksat 2,gulksat2.map,Ksat of subsoil for gully infil (mm/h),gulksat2                                          ,2");
    DEFmaps.append("2,Excluded areas,noncrit.map,areas to be excluded (1) and rest (0),nonfcrit                                    ,2");
    DEFmaps.append("2,Gully initial Width,gulwinit.map, initial gully width (m),gulwinit                                           ,2");
    DEFmaps.append("2,Gully initial Depth,guldinit.map, initial gully depth (m),guldinit                                           ,2");
*/
    OUTmaps.append("0,Basic Output");
    OUTmaps.append("2,Runoff,ro, runoff (l/s),outrunoff                                                                            ");
    OUTmaps.append("2,Concentration,conc, sediment concentration (g/l),outconc                                                     ");
    OUTmaps.append("2,Water height,wh, water height on surface (mm),outwh                                                          ");
    OUTmaps.append("2,Runoff height,rwh,runoff water height on surface (mm),outrwh                                                 ");
    OUTmaps.append("2,Transport cap.,tc, transport capacity (excl. channels)(g/l),outtc                                            ");
    OUTmaps.append("2,erosion,eros, erosion (kg),outeros                                                                           ");
    OUTmaps.append("2,deposition,depo, deposition (kg),outdepo                                                                     ");
    OUTmaps.append("2,Velocity,velo, velocity (m/s),outvelo                                                                        ");
    OUTmaps.append("2,Infiltration,inf,Cumulative infil (mm),outinf                                                                ");
    OUTmaps.append("2,Surface Storage,sstor, Surface storage (mm),outss                                                            ");
    OUTmaps.append("2,Channel volume,chanvol, Channel water volume (m3),outchvol                                                   ");

    OUTMCmaps.append("0,Multiclass output ");
    OUTMCmaps.append("2,Sed. Conc mu0,smu0,suspended sediment flux class 0 (mu),outmu0                                               ");
    OUTMCmaps.append("2,Sed. Conc mu1,smu1,suspended sediment flux class 1 (mu),outmu1                                               ");
    OUTMCmaps.append("2,Sed. Conc mu2,smu2,suspended sediment flux class 2 (mu),outmu2                                               ");
    OUTMCmaps.append("2,Sed. Conc mu3,smu3,suspended sediment flux class 3 (mu),outmu3                                               ");
    OUTMCmaps.append("2,Sed. Conc mu4,smu4,suspended sediment flux class 4 (mu),outmu4                                               ");
    OUTMCmaps.append("2,Sed. Conc mu5,smu5,suspended sediment flux class 5 (mu),outmu5                                               ");
    OUTMCmaps.append("2,D50,D50s,D50 suspended sediment (mu),outD50susp                                                              ");

    OUTNUTmaps.append("0,Nutrient Output                                                                                              ");
    OUTNUTmaps.append("2,Nuts P solut.,NPsol, P in solution (kg),outpsolut                                                            ");
    OUTNUTmaps.append("2,Nuts P susp.,NPsus, P in suspension (kg),outpsus                                                             ");
    OUTNUTmaps.append("2,Nuts P inf.,NPinf, P in solution (kg),outpinf                                                                ");
    OUTNUTmaps.append("2,Nuts NH4 solut.,NNH4sol, NH4 in solution (kg),outnh4solut                                                    ");
    OUTNUTmaps.append("2,Nuts NH4 susp.,NNH4sus, NH4 in suspension (kg),outnh4sus                                                     ");
    OUTNUTmaps.append("2,Nuts NH4 inf.,NNH4inf, NH4 in solution (kg),outnh4inf                                                        ");
    OUTNUTmaps.append("2,Nuts NO3 solut.,NNO3sol, NO3 in solution (kg),outNO3solut                                                    ");
    OUTNUTmaps.append("2,Nuts NO3 susp.,NNO3sus, NO3 in suspension (kg),outNO3sus                                                     ");
    OUTNUTmaps.append("2,Nuts NO3 inf.,NNO3inf, NO3 in solution (kg),outno3inf                                                        ");
    OUTNUTmaps.append("2,P Deposition,NPdep.map, P in clay deposition (kg/m2),outpdep                                                 ");
    OUTNUTmaps.append("2,NH4 Deposition,NNH4dep.map, NH4 in clay deposition (kg/m2),outnh4dep                                         ");
    OUTNUTmaps.append("2,NO3 Deposition,NNO3Dep.map, NO3 in clay deposition (kg/m2),outno3dep                                         ");
    OUTNUTmaps.append("2,P Detachment,NPdet.map, P in clay detachment (kg/m2),outpdet                                                 ");
    OUTNUTmaps.append("2,NH4 Detachment,NNH4det.map, NH4 in clay detachment (kg/m2),outnh4det                                         ");
    OUTNUTmaps.append("2,NO3 Detachment,NNO3Det.map, NO3 in clay detachment (kg/m2),outno3det                                         ");

    OUTGULmaps.append("0,Gully Output                                                                                                 ");
    OUTGULmaps.append("2,Gully depth,guld,gully depth (m),outguld                                                                     ");
    OUTGULmaps.append("2,Gully width,gulw,gully width (m),outgulw                                                                     ");
    OUTGULmaps.append("2,Gully area,gula,gully cross section area (m2),outgula                                                        ");
    OUTGULmaps.append("2,Gully fraction,gulf,cross section - fraction of foot^2 (-),outgulf                                           ");
    OUTGULmaps.append("2,Gully DEM,guldem,digital elevation model (m),outguldem                                                       ");

}
//--------------------------------------------------------------------
void lisemqt::change_MapNameModel(int parentrow, int selrow, bool setit)
{
    if (MapNameModel)
    {
        QModelIndex indexParent = MapNameModel->index(parentrow, 0);
        // select and expand chosen main and sub-branch
        if (selrow >= 0)
        {
            MapNameModel->setflag(setit, parentrow);
            // set top level node
            if (selrow == 0)
                for (int k = 0; k < MapNameModel->rowCount(indexParent); k++)
                    MapNameModel->setflag(setit, k, indexParent);
            else
                if (selrow >= 10)
                {
                selrow -= 10;
                MapNameModel->setflag(setit, selrow, indexParent);

                QModelIndex indexChild = MapNameModel->index(selrow, 0, indexParent);

                for (int k = 0; k < MapNameModel->rowCount(indexChild); k++)
                    MapNameModel->setflag(setit, k, indexChild);
            }

            if (setit)
            {
                treeView->expand(MapNameModel->index(parentrow,0));
                if (selrow >= 0)
                    treeView->expand(MapNameModel->index(selrow, 0, indexParent));
            }
            else
            {
                if (selrow >= 0)
                    treeView->collapse(MapNameModel->index(selrow, 0, indexParent));
                if (selrow == 0)
                    treeView->collapse(MapNameModel->index(parentrow,0));
            }


        }
//        treeView->resizeColumnToContents(0);
//        treeView->resizeColumnToContents(1);
//        treeView->resizeColumnToContents(2);
    }
}
//--------------------------------------------------------------------
void lisemqt::FillMapList()
{
    QStringList headers;
    headers  <<"Variable name"<< "Map name"<<"Description"<<" ";
    /*
    QFile file("d:/prgc/simpletreemodel/try.csv");
    file.open(QIODevice::ReadOnly);
    model = new TreeModel(headers, file.readAll());
    file.close();
*/
    MapNameModel = new TreeModel(headers, DEFmaps);

    treeView->setModel(MapNameModel);
    //treeView->resizeColumnToContents(0);
    //treeView->resizeColumnToContents(1);
    treeView->setColumnWidth(0,196);
    treeView->setColumnWidth(1,108);
    treeView->setColumnWidth(2,400);
    treeView->setColumnWidth(3,0);
    treeView->setColumnWidth(4,0);
    treeView->setAlternatingRowColors(true);

    change_MapNameModel(0,0, true);
    change_MapNameModel(1,0, true);
    change_MapNameModel(2,0, true);
    change_MapNameModel(3,0, true);
    treeView->collapseAll();

}
//--------------------------------------------------------------------
