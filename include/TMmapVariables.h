/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
\file TMmapVariables.h
\brief List of maps with descriptions and units. Linked directly in the model class.
*/

QVector <cTMap*> *inith; // swatre matrix potential nodes

cTMap
*DEM,                        //!< DEM [m]
*MBm,
*ShadeBW,                      //!< Shaded relief for display [0-1]
*DX,                         //!< cell length divided by cosine slope (so corrected for terrain gradient) [m]
*CellArea,                   //!< cell area = DX * _dx [m^2]
*Grad,                       //!< sine of the DEM gradient [-]
*LDD,                        //!< local drain direction map [-]
*Outlet,                     //!< main outlet of the catchment, value 5 in LDD map [-]
*PointMap,                   //!< map with output points, values > 0 [-]
*FlowBoundary,               //!< map with open boundary fior diffusive runoff (1) or closed boundary (0)
*WaterSheds,                 //!< map with numbered siubcatchments, must be 1,2,3 ... n
*QBoundFlow,
*DomainEdge,

*IDRainPoints,                   //!< rainfall zone map (clasified map, numers corrspond to raingaug number in rainfall file) [-]
*RainZone,                   //!< rainfall zone map (clasified map, numers corrspond to raingaug number in rainfall file) [-]
*ETZone,                     //!< rainfall zone map (clasified map, numers corrspond to raingaug number in rainfall file) [-]
*Rain,                       //!< map with rain from tis time intervall [m]
*Rainc,                      //!< map with rain from tis time intervall, spread over the surface (corrected or slope) [m]
*RainCumInt,                 //!< cumulative rainfall, as spreadoutover slope [m], needed for interception
*RainCumCrust,               //!< cumulative rainfall, as spreadoutover slope [m], needed for crusting
*RainCumFlat,                //!< cumulative rainfall [m]
*RainNet,                    //!< net rainfall after interception [m]
*LeafDrain,                  //!< drainge from canopy, storage overflow [m]
*CStor,                      //!< actual canopy storage [m]
*Interc,                     //!< actual canopy storage volume, corrected for surfaces without vegetation (like roads) [m^3]
*IntercETa,                     //!< actual canopy storage volume, corrected for surfaces without vegetation (like roads) [m^3]
*LCStor,                     //!< actual Litter storage [m]
*LInterc,                    //!< actual Litter storage volume, corrected for surfaces without vegetation (like roads) [m^3]
*DStor,                      //!< actual drum storage of rainwater [m^3]
*HStor,                      //!< actual roof storage of rainwater [m]
*IntercHouse,                //!< actual roof storage volume [m^3]
*HouseCover,                 //!< fraction cover of house in pixel [-]
*RoofStore,                  //!< Max storage of roof in [mm]
*DrumStore,                  //!< Max storage of rainwter drums [m^3]
*InterceptionmmCum,
*ETa,
*ETaCum,
*ETp,
*ETpCum,

*SnowmeltZone,               //!< snowmelt zone map, values corrspond to snowmelt gauge numbers [-]
*Snowcover,                  //!< snowmelt cover map, value 1.0 if there is snowcover, 0 without [-]
*Snowmelt,                   //!< snowmelt depth in water equivalent [m]
*Snowmeltc,                  //!< snowmelt depth in water equivalent, corrected for DEM gradient [m]
*SnowmeltCum,                //!< cumulative showmelt depth [m]

*WH,                         //!< water height on the surface [m]
//*WHold,                      //!< water height on the surface before infiltration [m]
//*WHnew,                      //!< water height on the surface before infiltration [m]
*WHrunoff,                   //!< water height available for runoff [m]
*WHmax,                      //!< max runoff wh in m for reporting
*WHstore,                    //!< water heigth stored in micro depressions [m]
*MicroStoreVol,
*WaterVolall,                //!< water volume total (incl surface storage) [m^3]
*WaterVolin,                 //!< water volume total before kin wave (after tochannel) [m^3]
//*flowmask,
//*WaterVolRunoff,                //!< water volume for runoff [m^3]

*FlowWidth,                  //!< width of the flow overland, based on ponded area/roughness, +roads etc [m]
*V,                          //!< velocity of overland flow [m/s]
*Alpha,                      //!< alpha in A = alphaQ^b
*Q,                          //!< discharge of overland flow before kin wave [m^3/s]
*DischargeUserPoints,
*QuserIn,
*WHbound,
*WHboundarea,
*WHboundRain,
*Qbase,
*GWVol,
*GWWH,
*GWU,
*GWV,
*GWN,
*GWWHmax,
*GWdeep,
*GWrecharge,
*GWout,
*GWz,
*GWgrad,
*Qn,                         //!< new discharge of overland flow after kin wave [m^3/s]
*Qdiag,
*VH,
*QinKW,                      //!< new Q kinematic wave
*QKW,
*Qm3total,
*Qm3max,
*FHI,
*Qoutput,                    //!< new discharge for output purposes, sum of overland flow and channel, converted [l/s]
*Qs,                         //!< sediment discharge before kin wave [kg/s]
*Qsn,                        //!< new sediment discharge after kin wave [kg/s]
*SinKW,                      //!< New Sed flux kinematic wave
*Qsoutput,                   //!< sediment outflow for screen/file output, sum of overland flow and channel [kg/s]
//*q,                          //!< infiltration surplus going in kin wave (<= 0) [m2/s]
*R,                          //!< hydraulic radius overland flow [m]
*N,                          //!< Manning's n
*Norg,                          //!< Manning's n
*RR,                         //!< Random roughness, locally converted to m [cm]
*MDS,                        //!< Maximum depression storage [m]
//*fpa,                        //!< fraction ponded area [-]
*SoilWidthDX,                //!< width of soil surface, excluding roads and channels [m]
*RoadWidthDX,                //!< width of tarred roads [m]
*RoadWidthHSDX,
*StoneFraction,              //!< fraction of stones on the surface, affects splash [-]
*CompactFraction,            //!< fraction compacted at the surface, uses ksat compact [-]
*CrustFraction,              //!< fraction crusted at the surface, uses ksat crust [-]
*CrustFraction0,              //!< fraction crusted at the surface, uses ksat crust [-]
//*RepellencyFraction,         //!< fraction of water repellency of node 1 in Swatre [-]
//*RepellencyCell,             //!< Cell included in water repellency in Swatre [-]
*HardSurface,                //!< value 1 if 'hard' surface: no interception, infiltration, detachment [-]
*fractionImperm,            //!<// 0 is fully permeable, 1 = impermeable [-]
*runoffTotalCell,
*hSwatre,
*thetaSwatre,

*PlantHeight,                //!< height of vegetation/crops [m]
*Cover,                      //!< vegetation canopy cover fraction [-]
*Litter,                     //!< vegetation litter cover fraction [-]
*CanopyStorage,              //!< canopy storage [m]
*LAI,                        //!< leaf area index [m^2/m^2]
*LandUnit,                   //!< land unit class (> 0) [-]

*Cohesion,                   //!< total cohesion of the soil surface: coh soil *(1-cover) + coh plant (cover) [kPa]
*RootCohesion,               //!< cohesion soil [kPa]
*CohesionSoil,               //!< cohesion by plant roots [kPa]
*Y,                          //!< erosion efficiency 0-1, basd on cohesion [-]
*AggrStab,                   //!< aggregate stability, median of drops in lowe test [-]
*SplashStrength,                   //!< aggregate stability, median of drops in lowe test [-]
//*splashb,                   //!< aggregate stability, median of drops in lowe test [-]
*D50,                        //!< median of grainsize distribution [mu]
*D90,                        //!< 90 % of grainsize distribution is below this value [mu]
*D50CH,                        //!< median of grainsize distribution [mu]
*D90CH,                        //!< 90 % of grainsize distribution is below this value [mu]
*cgovers,
*dgovers,
*DETSplash,                  //!< splash detachment [kg/cell]
*DETSplashCum,
*DETFlow,                    //!< flow detachment [kg/cell]
*DETFlowCum,
*DEPCum,
*DEP,                        //!< deposition [kg/cell]
*TC,                         //!< transport capacity [kg/m^3]
*Conc,                       //!< sediment concentration in flow [kg/m^3]
*Sed,                        //!< sediment content of flow [kg]
*SettlingVelocitySS,           //!< settling velocity according to Stokes [m/s]
*SettlingVelocityBL,           //!< settling velocity according to Stokes [m/s]

// Pesticides
*PMmw,                      //!< Map with mass of pesticides in soil part of mixing zone [mg]
*PMms,                      //!< Map with mass of pesticides in soil part of mixing zone [mg]
*PMrw,                      //!< mass of pesticide in runoff water [mg]
*PMrs,                      //!< mass of pesticide in runoff sediment [mg]
*PMsoil,                    //!< mass of pesticide in the soil layer without mixing zone [mg]
*PCrw,                      //!< concentration of pesticide in runoff water [mg/L]
*PCrs,                      //!< concentration of pesticide in runoff sediment [mg/kg]
*PCms,                      //!< concentration of pesticide in soil of mixing zone [mg/kg]
*PCmw,                      //!< concentration of pesticide in water of mixing zone [mg/L]
*PQrw,                      //!< flux of pesticide in runoff water [mg/sec]
*PQrs,                      //!< flux of pesticide in runoff sediment [mg/sec]
*PMinf,                     //!< mass of pesticide in infiltrating water [mg]
*zm,                        //!< depth of the mixing layer [m]
*zs,                        //!< depth of the soil layer containing pesticides [m]
*SpinKW,                    //!< sum upstream influx Qpsn [mg/sec]
*QpinKW,                    //!< sum upstream influx Qpn [mg/sec]
*Qpw,                       //!< dissolved pesticide flux based on Qp [mg/sec]
*Qps,                       //!< pesticide sediment flux based on Qs [mg/sec]
*PCs,                       //!< concentration of pesticide in pesticide soil layer 1 [mg/kg]
*Theta_mix,                 //!< theta of the mixing layer [-]
*pmsdet,                     //!< mass of detached pesticide [mg]
*pmsdep,                     //!< mass of deposited pesticide [mg]
*pmwdep,                    //!< mass of deposited pesticide [mg]
*pmwdet,                     //!< mass of detatched pesticide [mg]
*WVji1,                     //!< water volume in cell at j, i+1 [m3]
*SedMassIn,                  //!< sediment mass in to kinematic wave [kg]
*SedAfterSplash,             //!< sediment mass in flow after splash [kg]
*PMsplash,                   //!< mass detached sorbed pesticide by splash erosion [mg]
*PMflow,                    //!< mass detached sorbed pesticide by flow detachement[mg]
*PMdep,                     //!< mass deposited sorbed pesticide [mg]
*totalPPlossmap,             //!< total loss of PP pesticide [mg/m2]
*totalDPlossmap,             //!< total loss of DP pesticide [mg/m2]
*test_map,


// infiltration
*Fcum,                       //!< cumulative infiltration [m]
//*FSurplus,                   //!< surplus infiltration for kinematic wave, calculated as actual infil - potential infil [m]
//*fact,                       //!< actual infiltration rate [m/s]
//*fpot,                       //!< potential infiltration rate [m/s]
//*InfilVolKinWave,            //!< volume infiltrated in the kin wave (slope and channel) in this timestep [m^3]
*InfilVol,                   //!< volume of water infiltrated in this timestep [m^3]
*ChannelInfilVol,                   //!< volume of water infiltrated in this timestep [m^3]

*InfilVolCum,                //!< cumulative infiltration volume for mass balance and map report [m^3]
*InfilmmCum,                 //!< cumulative infiltration volume for map report and drawing [mm]
*InfilVolFlood,

*Lw,
*Lwmm,

*ThetaS1,                    //!< porosity soil layer 1 [-]
*ThetaI1,                    //!< initial moisture content soil layer 1 [-]
*ThetaI1a,                    //!< initial moisture content soil layer 1 [-]
*Psi1,                       //!< intial suction head wetting front soil layer 1 (input map is in cm) [m]
*ThetaR1,
*ThetaFC1,
*Ksat1,                      //!< saturated hydraulic conductivity soil layer 1 (input is in mm/h) [m/s]
*SoilDepth1,                 //!< depth to end soil layer 1 (input is in mm) [m]
*SoilDepth1init,                 //!< depth to end soil layer 1 (input is in mm) [m]

*ThetaS2,                    //!< porosity soil layer 2 [-]
*ThetaI2,                    //!< initial moisture content soil layer 2 [-]
*ThetaI2a,                    //!< initial moisture content soil layer 2 [-]
*ThetaR2,
*ThetaFC2,
*Psi2,                       //!< intial suction head wetting front soil layer 2 (input map is in cm) [m]
*Ksat2,                      //!< saturated hydraulic conductivity soil layer 2 (input is in mm/h) [m/s]
*SoilDepth2,                 //!< depth to end soil layer 2 (input is in mm) [m]
*SoilDepth2init,                 //!< depth to end soil layer 2 (input is in mm) [m]

*ThetaS3,                    //!< porosity soil layer 1 [-]
*ThetaI3,                    //!< initial moisture content soil layer 1 [-]
*ThetaI3a,                    //!< initial moisture content soil layer 1 [-]
*Psi3,                       //!< intial suction head wetting front soil layer 1 (input map is in cm) [m]
*ThetaR3,
*ThetaFC3,
*Ksat3,                      //!< saturated hydraulic conductivity soil layer 1 (input is in mm/h) [m/s]
*SoilDepth3,                 //!< depth to end soil layer 1 (input is in mm) [m]
*SoilDepth3init,                 //!< depth to end soil layer 1 (input is in mm) [m]

*lambda1,
*lambda2,
*lambda3,
*vgalpha1,
*vgalpha2,
*vgalpha3,
*vgn1,
*vgn2,
*vgn3,
*psi1ae,
*psi2ae,
*psi3ae,

*KsatCrust,                  //!< saturated hydraulic conductivity crusted soil surface (input is in mm/h) [m/s]
*PoreCrust,                //!< saturated hydraulic conductivity compacted soil surface (input is in mm/h) [m/s]
*KsatCompact,                //!< saturated hydraulic conductivity compacted soil surface (input is in mm/h) [m/s]
*PoreCompact,                //!< saturated hydraulic conductivity compacted soil surface (input is in mm/h) [m/s]
*KsatGrass,                  //!< saturated hydraulic conductivity grass strip (input is in mm/h) [m/s]
*PoreGrass,                  //!< Porosity grass strip (input in cm3/cm3)
*CohGrass,                   //!< Cohesion grass strip (input in kPa)
*Ksateff,                    //!< effective saturated hydraulic conductivity (input is in mm/h) [m/s]
*Poreeff,
*Thetaeff,
*chanmask3,

*Perc,
*PercmmCum,
*GrassFraction,              //!< fraction of grasstrip in a cell [-]
*SedimentFilter,             //!< sediment deposited in the sediment trap in kg/m2
*SedMaxVolume,               //!< maxvol of sediment in that can be trapped in m3
*GrassWidthDX,               //!< width of grasstrip in [m]

//swatre
*OMcorr,
*DensFact,
*ProfileID,                  //!< SWATRE profile unit number map
*ProfileIDCrust,             //!< SWATRE profile unit number map for crusted areas
*ProfileIDCompact,           //!< SWATRE profile unit number map for compacted areas
*ProfileIDGrass,             //!< SWATRE profile unit number map for grass strips
*SwatreOutput,               //!< SWATRE cells flagged for output
//*inith,                      //!< SWATRE inithead in -cm

*LDDChannel,                 //!<
*LDDbaseflow,
*ChannelWidthO,               //!<
*ChannelWidth,               //!<
*ChannelDepth,               //!<
*ChannelSide,                //!<
*ChannelQSide,                //!<
//*ChannelQb,                   //!<
*ChannelQ,                   //!<
*ChannelQn,                  //!<
*ChannelQntot,
*ChannelQs,                  //!<
*ChannelQsn,                 //!<
*ChannelQBLs,                  //!<
*ChannelQBLsn,                 //!<
*ChannelQSSs,                  //!<
*ChannelQSSsn,                 //!<
*ChannelGrad,                //!<
*ChannelV,                   //!<
*ChannelN,                   //!<
*ChannelWH,                  //!<
*ChannelPerimeter,
*ChannelWidthB,
//*ChannelCos,
//*ChannelWHExtended,                  //!<
//*ChannelVolExtended,                  //!<
*ChannelWaterVol,            //!<
//*Channelq,                   //!<
*ChannelAlpha,               //!<
*ChannelDX,                  //!<
*ChannelKsat,                //!<
*ChannelInfM3,                //!<

*ChannelAdj,                //!<
*CHAdjDX,                //!< channel adjusted DX
*BaseflowL,

*cosGrad,
*tanGrad,
*BulkDensity,
*AngleFriction,
*FSlope,

// channel erosion
//*ChannelPerimeter,           //!<
*ChannelDetFlow,             //!<
*ChannelDep,                 //!<
//*ChannelSed,                 //!<
*ChannelBLSed,                 //!<
*ChannelSSSed,                 //!<
*ChannelBLTC,                 //!<
*ChannelSSTC,                 //!<
*ChannelBLDepth,                 //!<
*ChannelSSDepth,                 //!<
*ChannelConc,                //!<
*ChannelBLConc,                //!<
*ChannelSSConc,                //!<
*ChannelTC,                  //!<
*ChannelCohesion,            //!<
*ChannelY,                   //!<
*ChannelPAngle,               //!<
*ChannelQsr,

//baseflow
*BaseFlowDischarges,
*BaseFlowInitialVolume,
*BaseFlowInflow,

// flood maps
*floodHmxMax,                    //!<
*floodTime,                    //!<
*floodTimeStart,                //!<
*floodVMax,                    //!<
*floodVHMax,                    //!<
*maxChannelflow,                    //!<
*maxChannelWH,                    //!<
*hmx,                        //!<
*hmxWH,                        //!<
*hmxrunoff,
*hmxInit,                    //!<
*FloodDomain,                //!<
*Buffers,                    //!<
*GridRetention,                    //!<
*GridRetentionAct,
*ChanRetention,                    //!<
*ChanRetentionAct,
*ChannelDiameter,                //!<
*ChannelCulvert,                //!<
*ChannelMaxQ,                //!<
*ChannelMaxAlpha,                //!<
*ChannelMaxArea,
*FloodWaterVol,                //!<
*RunoffWaterVol,                //!<

//*FloodZonePotential,                //!<
*FloodDT,
*Uflood,*Vflood,
*hs, //*vs, *us,
*gflowx,
*gflowy,
*hllx12_0,
*hlly12_0,
*hllx21_1,
*hllx21_2,
*hlly21_1,
*hlly21_2,


//FULLSWOF2D with Sediment
*BLDepthFlood,
*SSDepthFlood,
*BLDetFlood,
*BLTCFlood,
*SSTCFlood,
*SSDetFlood,
*DepFlood,
*BLCFlood,
*BLFlood,
*SSCFlood,
*SSFlood,

*LDDTile,                    //!< LDD network of tile drains, must be connected to outlet
*TileDrainSoil,              //!< drain volume from layer
*TileDiameter,                  //!< total width of drains in cell (m)
*TileArea,                  //!< total width of drains in cell (m)
*TileWidth,                  //!< total width of drains in cell (m)
*TileHeight,                 //!< height of drain (m)
*TileDepth,                  //!< depth of tiles in soil below surface (m)
//*TileInlet,               //!< sinkhole on surface connecting to tiledrains (m2)
*TileQ,                      //!< water flux in drains m3/s
*TileMaxQ,                      //!< water flux in drains m3/s
*TileQn,                     //!< new water flux in drains m3/s
*TileGrad,                   //!< gradient of the tiledrain system
*TileN,                      //!< mannings inside the tiledrains
*TileWaterVol,               //!< water volume in the tiledrains (m3)
*TileWaterVolSoil,           //!< water volume in the tiledrains from the soil only, used for mass bal corection (m3)
*RunoffVolinToTile,          //!< can be used for shortcut of surface pits to tile system
*TileAlpha,                  //!< alpha in tile drain, in A = alpha*Q^beta
*TileMaxAlpha,                      //!< water flux in drains m3/s

*TotalChanDetMap,                //!<
*TotalChanDepMap,                //!<
*TotalSoillossMap,           //!<
*TotalSed,                   //!<
*TotalConc,                  //!<

*tm,                         //!< Auxilary map
*tma,                        //!< Auxilary map
*tmb,                        //!< Auxilary map
*tmc,                        //!< Auxilary map
*tmd,                        //!< Auxilary map
*tmshow,
//display combinations
*COMBO_V,
*COMBO_SS,
*COMBO_BL,
*COMBO_TC,
*ChannelDepthExtended,
*ChannelWidthExtended,
*ChannelNeighborsExtended,
*ChannelSourceXExtended,
*ChannelSourceYExtended,
*ChannelMaskExtended,
*ChannelBoundaryExtended,
*ChannelBoundaryLExtended,
*ChannelBoundaryRExtended,

*FlowBarrier,                //!< Flow barriers type
*FlowBarrierN,               //!< Flow barriers height North of cell
*FlowBarrierW,               //!< Flow barriers height West of cell
*FlowBarrierS,               //!< Flow barriers height South of cell
*FlowBarrierE,               //!< Flow barriers height East of cell
*FlowBarrierNT,              //!< Flow barriers end timing North of cell
*FlowBarrierWT,              //!< Flow barriers end timing West of cell
*FlowBarrierST,              //!< Flow barriers end timing South of cell
*FlowBarrierET               //!< Flow barriers end timing East of cell

;
cTRGBMap * RGB_Image;
