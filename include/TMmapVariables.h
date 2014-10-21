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
\file TMmapVariables.h
\brief List of maps with descriptions and units. Linked directly in the model class.
*/

// separte here for easier Doxygen comments

CTMap

*DEM,                        //!< DEM [m]
*Shade,                      //!< Shaded relief for display [0-1]
*DX,                         //!< cell length divided by cosine slope (so corrected for terrain gradient) [m]
*CellArea,                   //!< cell area = DX * _dx [m^2]
*Grad,                       //!< sine of the DEM gradient [-]
*LDD,                        //!< local drain direction map [-]
*Outlet,                     //!< main outlet of the catchment, value 5 in LDD map [-]
*PointMap,                   //!< map with output points, values > 0 [-]

*RainZone,                   //!< rainfall zone map (clasified map, numers corrspond to raingaug number in rainfall file) [-]
*Rain,                       //!< map with rain from tis time intervall [m]
*Rainc,                      //!< map with rain from tis time intervall, spread over the surface (corrected or slope) [m]
*RainCum,                    //!< cumulative rainfall, as spreadoutover slope [m]
*RainCumFlat,                //!< cumulative rainfall [m]
*RainNet,                    //!< net rainfall after interception [m]
*LeafDrain,                  //!< drainge from canopy, storage overflow [m]
*CStor,                      //!< actual canopy storage [m]
*Interc,                     //!< actual canopy storage volume, corrected for surfaces without vegetation (like roads) [m^3]
//houses
*HStor,                      //!< actual roof storage of rainwater [m]
*DStor,                      //!< actual drum storage of rainwater [m^3]
*IntercHouse,                //!< actual roof storage volume [m^3]
*HouseCover,                 //!< fraction cover of house in pixel [-]
*HouseWidthDX,
*RoofStore,                  //!< Max storage of roof in [mm]
*DrumStore,                  //!< Max storage of rainwter drums [m^3]

*SnowmeltZone,               //!< snowmelt zone map, values corrspond to snowmelt gauge numbers [-]
*Snowcover,                  //!< snowmelt cover map, value 1.0 if there is snowcover, 0 without [-]
*Snowmelt,                   //!< snowmelt depth in water equivalent [m]
*Snowmeltc,                  //!< snowmelt depth in water equivalent, corrected for DEM gradient [m]
*SnowmeltCum,                //!< cumulative showmelt depth [m]

*WH,                         //!< water height on the surface [m]
*WHbef,                      //!< water height on the surface before infiltration [m]
*WHroad,                     //!< water height on the roads [m]
*WHrunoff,                   //!< water height available for runoff [m]
//*WHrunoffCum,                //!< cumulative runoff water height for output to file series only [mm]
*WHstore,                    //!< water heigth stored in micro depressions [m]
*WaterVolall,                //!< water volume total (incl surface storage) [m^3]
*WaterVolin,                 //!< water volume total before kin wave (after tochannel) [m^3]

*FlowWidth,                  //!< width of the flow overland, based on ponded area/roughness, +roads etc [m]
*V,                          //!< velocity of overland flow [m/s]
*Alpha,                      //!< alpha in A = alphaQ^b
//*AlphaF,
//*QF,
//*QnF,
*Q,                          //!< discharge of overland flow before kin wave [m^3/s]
*Qn,                         //!< new discharge of overland flow after kin wave [m^3/s]
//*Qoutflow,                   //!< new discharge after kin wave at outflow point [m^3/s]
*QinKW,
*Qoutput,                    //!< new discharge for output purposes, sum of overland flow and channel, converted [l/s]
*Houtput,                    //!< new discharge for output purposes, sum of overland flow and channel, converted [l/s]
*Qs,                         //!< sediment discharge before kin wave [kg/s]
*Qsn,                        //!< new sediment discharge after kin wave [kg/s]
//*Qsoutflow,                  //!< new sediment discharge after kin wave at outflow point [kg/s]
*Qsoutput,                   //!< sediemnt outflow for screen/file output, sum of overland flow and channel [kg/s]
*q,                          //!< infiltration surplus going in kin wave (<= 0) [m2/s]
*R,                          //!< hydraulic radius overland flow [m]
*Perim,                      //!< perimeter overland flow [m]
*N,                          //!< Manning's n
*RR,                         //!< Random roughness, locally converted to m [cm]
*MDS,                        //!< Maximum depression storage [m]
*fpa,                        //!< fraction ponded area [-]
*SoilWidthDX,                //!< width of soil surface, excluding roads and channels [m]
*RoadWidthDX,                //!< width of tarred roads [m]
*StoneFraction,              //!< fraction of stones on the surface, affects splash [-]
*CompactFraction,            //!< fraction compacted at the surface, uses ksat compact [-]
*CrustFraction,              //!< fraction crusted at the surface, uses ksat crust [-]
*RepellencyFraction,         //!< fraction of water repellency of node 1 in Swatre [-]
*RepellencyCell,             //!< Cell included in water repellency in Swatre [-]
*HardSurface,                //!< value 1 if 'hard' surface: no interception, infiltration, detachment [-]
*runoffFractionCell,
*runoffTotalCell,

*PlantHeight,                //!< height of vegetation/crops [m]
*Cover,                      //!< vegetation canopy cover fraction [m]
*CanopyStorage,              //!< canopy storage [m]
*LAI,                        //!< leaf area index [m^2/m^2]
*LandUnit,                   //!< land unit class (> 0) [-]
*WheelWidth,                 //!< not used yet, width of wheel tracks [m]
*WheelWidthDX,               //!< not used yet, width of wheel tracks [m]
*GullyWidthDX,               //!< not used yet, width of gullies [m]

*Cohesion,                   //!< total cohesion of the soil surface: coh soil *(1-cover) + coh plant (cover) [kPa]
*RootCohesion,               //!< cohesion soil [kPa]
*CohesionSoil,               //!< cohesion by plant roots [kPa]
*Y,                          //!< erosion efficiency 0-1, basd on cohesion [-]
*AggrStab,                   //!< aggregate stability, median of drops in lowe test [-]
*D50,                        //!< median of grainsize distribution [mu]
*DETSplash,                  //!< splash detachment [kg/m^2]
*DETFlow,                    //!< flow detachment [kg/m^2]
*DEP,                        //!< deposition [kg/m^2]
*TC,                         //!< transport capacity [kg/m^3]
*Conc,                       //!< sediment concentration in flow [kg/m^3]
*Sed,                        //!< sediment content of flow [kg]
*CG,                         //!< parameter Govers in TC equation
*DG,                         //!< parameter Govers in TC equation
*SettlingVelocity,           //!< settling velocity according to Stokes [m/s]


*Vup,                        //!< updated runoff velocity for pesticides transport
*Vup_old,
*PCA,                        //!< applied dose [kg/m2]
*epsil,                      //!< mixing layer depth (m]
*KD,                         //!< soil water partition coefficient [m3/kg]
*kr,                         //!< rate at which solute desorb [min-1]
*rhob,                       //!< soil bulk density [kg/m3]
*C,                          //!< Pesticide concentration in dissolved form in runoff water [kg/m3]
*Cold,
*CM,                         //!< Pesticide concentration in dissolved form in the mixing zone [kg/m3]
*CS,                         //!< Pesticide concentration in sorbed form in the mixing zone [kg/m3]
*C_N,
*CM_N,
*CS_N,
*C_K,
*C_Kold,
*CM_K,
*CS_K,
*C_Kexplicit,
*CS_Kexplicit,
*CM_Kexplicit,
*CM_Kexplicitold,
*CS_Kexplicitold,
*Qp,
*Qpn,
*C_Kn,
*K1,
*Kfilm,
*pestiinf,
*pestiinfold,
*poro,
*AX,
*Fkold,
*Fk,
*Fmk,
*flagpest,
*PMassApplied,
*PRunoffSpatial,
*PDisMixing,
*PSorMixing,
*PInfilt,
*PStorage,
*PRunoffSpatialex,
*PDisMixingex,
*PSorMixingex,
*PInfiltex,
*Qin,
*Sin,
*Pest,
*Fin,
*Pdetach,
*PCinfilt,
*PCfilmexit,

*Fcum,                       //!< cumulative infiltration [m]
*FSurplus,                   //!< surplus infiltration for kinematic wave, calculated as actual infil - potential infil [m]
*hesinfil,
*FFull,                      //!< map flagging when the soil is full
*fact,                       //!< actual infiltration rate [m/s]
*fpot,                       //!< potential infiltration rate [m/s]
*InfilVolKinWave,            //!< volume infiltrated in the kin wave (slope and channel) in this timestep [m^3]
*InfilVol,                   //!< volume of water infiltrated in this timestep [m^3]
*InfilVolCum,                //!< cumulative infiltration volume for mass balance and map report [m^3]
*InfilmmCum,                 //!< cumulative infiltration volume for map report and drawing [mm]
*InfilVolFlood,

//*FfSurplus,                   //!< surplus infiltration for flooding, calculated as actual infil - potential infil [m]
//*Ffcum,                       //!< cumulative infiltration [m]
//*ffact,                       //!< actual infiltration rate [m/s]
//*ffpot,                       //!< potential infiltration rate [m/s]
//*FfFull,                      //!< map flagging when the soil is full

*Lf1,
*Lf2,

*ThetaS1,                    //!< porosity soil layer 1 [-]
*ThetaI1,                    //!< initial moisture content soil layer 1 [-]
*Psi1,                       //!< intial suction head wetting front soil layer 1 (input map is in cm) [m]
*Ksat1,                      //!< saturated hydraulic conductivity soil layer 1 (input is in mm/h) [m/s]
*SoilDepth1,                 //!< depth to end soil layer 1 (input is in mm) [m]
*L1,                         //!< depth of wetting front in layer 1 [m]
*Soilwater,                  //!< actual soil water content [-]

*ThetaS2,                    //!< porosity soil layer 2 [-]
*ThetaI2,                    //!< initial moisture content soil layer 2 [-]
*Psi2,                       //!< intial suction head wetting front soil layer 2 (input map is in cm) [m]
*Ksat2,                      //!< saturated hydraulic conductivity soil layer 2 (input is in mm/h) [m/s]
*SoilDepth2,                 //!< depth to end soil layer 2 (input is in mm) [m]
*L2,                         //!< depth of wetting front in layer 2 [m]
*Soilwater2,                  //!< actual soil water content layer 2 [-]

*KsatCrust,                  //!< saturated hydraulic conductivity crusted soil surface (input is in mm/h) [m/s]
*KsatCompact,                //!< saturated hydraulic conductivity compacted soil surface (input is in mm/h) [m/s]
*KsatGrass,                  //!< saturated hydraulic conductivity grass strip (input is in mm/h) [m/s]
*Ksateff,                    //!< effective saturated hydraulic conductivity (input is in mm/h) [m/s]
*L1gr,                       //!< depth wetting front under grass strip layer 1 [m]
*L2gr,                       //!< depth wetting front under grass strip layer 2 [m]
*factgr,                     //!< actual infiltration rate fo grassstrip [m/s]
*fpotgr,                     //!< potential infiltration rate fo grassstrip [m/s]
*Fcumgr,                     //!< cumulative infiltration under grassstrips [m]
*WHGrass,                    //!< water level on a grassstrip [m]
*GrassFraction,              //!< fraction of grasstrip in a cell [-]
*GrassWidthDX,               //!< width of grasstrip in [m]
*thetaTop,                   //!< average theta of node 0 and 1 for water repelency and nutrients

*ProfileID,                  //!< SWATRE profile unit number map
*ProfileIDCrust,             //!< SWATRE profile unit number map for crusted areas
*ProfileIDCompact,           //!< SWATRE profile unit number map for compacted areas
*ProfileIDGrass,             //!< SWATRE profile unit number map for grass strips

*LDDChannel,                 //!<
*RunoffVolinToChannel,       //!<
*ChannelWidth,               //!<
*ChannelSide,                //!<
*ChannelQ,                   //!<
*ChannelQn,                  //!<
*ChannelQntot,
*ChannelQs,                  //!<
*ChannelQsn,                 //!<
*ChannelGrad,                //!<
*ChannelV,                   //!<
*ChannelN,                   //!<
*ChannelWH,                  //!<
*ChannelWaterVol,            //!<
*Channelq,                   //!<
*ChannelAlpha,               //!<
*ChannelWidthUpDX,           //!<
*ChannelAdj,
*ChannelPerimeter,           //!<
*ChannelDX,                  //!<
*ChannelKsat,                //!<
*ChannelStore,               //!<
*ChannelDetFlow,             //!<
*ChannelDep,                 //!<
*ChannelSed,                 //!<
*ChannelConc,                //!<
*ChannelTC,                  //!<
*SedToChannel,               //!<
*ChannelCohesion,            //!<
*ChannelY,                   //!<
// flood maps
*ChannelDepth,               //!<
*UVflood,                     //!<
*Qflood,
*floodHmxMax,
*floodactive,
*timeflood,
*maxChannelflow,
*maxChannelWH,
*hmxInit,
//explicit solution
*Hx,                         //!<
*Hmx,                        //!<
*hx,                         //!<
*hmx,                        //!<
*Nx,                         //!<
*dHdLx,                      //!<
*Qxsum,                      //!<
*qx0,                        //!<
*qx1,                        //!<
*qx2,                        //!<
*qx3,                        //!<
//explicit solution
*FloodDomain,                //!<
*Barriers,                    //!<
*ChannelMaxQ,                //!<
*ChannelLevee,                //!<
*FloodWaterVol,
*FloodZonePotential,
*FloodEdge,
//*FloodVoltoChannel,
*FloodTimeStart,
*FloodTimeEnd,

// FULLSWOF2D
*z1r, *z1l, *z2r, *z2l,
*delta_z1, *delta_z2,
*h1r, *h1l, *h2r, *h2l,
*h1d, *h1g, *h2d, *h2g,
*v1r, *v1l, *v2r, *v2l,
*u1r, *u1l, *u2r, *u2l,
*delzc1, *delzc2,
*delz1, *delz2,
*f1, *f2, *f3, *cflx,
*g1, *g2, *g3, *cfly,
*hs, *vs, *us,
*hsa, *vsa, *usa,
*Uflood,*Vflood,
//*q1flood,*q2flood,
*som_z1,*som_z2,

*LDDTile,                    //!< LDD network of tile drains, must be connected to outlet
*TileDrainSoil,              //!< drain volume from layer
*TileWidth,                  //!< total width of drains in cell (m)
*TileHeight,                 //!< height of drain (m)
*TileDepth,                  //!< depth of tiles in soil below surface (m)
*TileSinkhole,               //!< sinkhole on surface connecting to tiledrains (m2)
*TileQ,                      //!< water flux in drains m3/s
*TileQn,                     //!< new water flux in drains m3/s
*TileQs,                     //!< sediment flux in drains kg/s
*TileQsn,                    //!< new sediment flux in drains kg/s
//*TileQoutflow,               //!< water outflow in outlet
*TileGrad,                   //!< gradient of the tiledrain system
*TileN,                      //!< mannings inside the tiledrains
*TileWH,                     //!< water height in the tile drains (m)
*TileWaterVol,               //!< water volume in the tiledrains (m3)
*TileWaterVolSoil,           //!< water volume in the tiledrains from the soil only, used for mass bal corection (m3)
*Tileq,                      //!< possible drainage inside tiles, not used
*RunoffVolinToTile,          //!< can be used for shortcut of surface pits to tile system
*TileAlpha,                  //!< alpha in tile drain, in A = alpha*Q^beta
*TileDX,                     //!< cell length in tile drain, dx/cos angle
*TileV,                      //!< velocity in tile drain m/s

*BufferID,                   //!<
*BufferVol,                  //!<
*BufferSed,                  //!<
*ChannelBufferSed,           //!<
*ChannelBufferVol,           //!<
*BufferVolInit,              //!<
*BufferSedInit,              //!<
*ChannelBufferSedInit,       //!<
*ChannelBufferVolInit,       //!<

*TotalDetMap,                //!<
*TotalDepMap,                //!<
*TotalSoillossMap,           //!<
*TotalSed,                   //!<
*TotalConc,                  //!<

//Mu[6],                       //!< multiclass fraction of the grainsize in the 6 classes, 6 maps sum to 1.0
//CGm[6],                      //!< multiclass TC coefficient for this texture class
//DGm[6],                       //!< multiclass TC coeficient for this texture class
*difkin,
*tm,                         //!< Auxilary map
*tma,                        //!< Auxilary map
*tmb,                        //!< Auxilary map
*tmc                        //!< Auxilary map

;

