/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact:
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
 \file lisChannelflood.cpp
 \brief Channel flood using a various solutions of St Venant equations: \n
        and more stable 1st and 2nd order st Venant following the fullSWOF2D code (univ Orleans)\n
        called before ChannelFlow(), takes old channel overflow height and spreads it out, puts new channelWH \n
        back into channel before kin wave of channel is done in ChannelFlow()

functions: \n
- void TWorld::ChannelOverflow(void) Mixing of flood and overflow in channel cells, source of overflow
- void TWorld::ChannelFlood(void) Calculate channelflood height maps (hmx, QFlood, UVFlood) and FloodDomain
*/

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

#define SQRT2G 4.42869

void TWorld::ChannelOverflow(cTMap *_h, cTMap *V)
{
    if (!SwitchIncludeChannel)
         return;

     #pragma omp parallel for num_threads(userCores)
     FOR_ROW_COL_MV_CHL {
         if (ChannelMaxQ->Drc <= 0) {
             double chdepth = ChannelDepth->Drc;
             double dH = std::max(0.0, (ChannelWH->Drc-chdepth));
             double H = _h->Drc;

             if (dH <= 1e-6 && H <= 1e-6)
                 continue;
             // no flow activity then continue

             if (fabs(dH - H) < 1e-6)
                 continue;
             // no diff in water level, no flow, continue

             // VELOCITIES
             double VtoChan = V->Drc;
             double fracA = std::min(1.0, _dt*VtoChan/(0.5*ChannelAdj->Drc));
             // fraction from _h to channel based on average flood velocity
             double VfromChan = sqrt(2*GRAV*dH); //Bernoulli

             VfromChan = 0.56*std::sqrt(2*GRAV)*std::pow(dH, 0.5)*sqrt(1-pow((dH-H)/dH,1.5));

             //see https://www.engineeringtoolbox.com/velocity-head-d_916.html
             double fracC = std::min(1.0, _dt*VfromChan/(0.5*ChannelAdj->Drc));
             // fraction from channel to surrounding

             double cwa = ChannelWidth->Drc/ChannelAdj->Drc;

             bool dosimpel = false;

             if (dH > H)   // flow from channel
             {
                 double dwh = fracC * (dH-H);
                 // amount flowing from channel
                 if (H + dwh*cwa > dH-dwh) {
                     // if flow causes situation to reverse (channel dips below _h)
                     dosimpel = true;
                 } else {

                     _h->Drc  += dwh*cwa;
                     ChannelWH->Drc -= dwh;

                     if(SwitchErosion) {
                         double sed = ChannelSSConc->Drc * dwh*ChannelWidth->Drc*ChannelDX->Drc;
                         ChannelSSSed->Drc -= sed;
                         SSFlood->Drc += sed;
                     }
                 }
             }
             else   // flow to channel
             {
                 double dwh = fracA * (H-dH);
                 // amount flowing to channel
                 if (dH + dwh/cwa > H-dwh) {
                     // if too much flow
                     dosimpel = true;
                 } else {
                     _h->Drc -= dwh;
                     ChannelWH->Drc += (dwh/cwa);
                     if(SwitchErosion) {
                         double sed = fracA*SSFlood->Drc;
                         ChannelSSSed->Drc += sed;
                         SSFlood->Drc -= sed;
                     }
                 }
             }

             // instantaneous waterlevel exquilibrium acccross channel and adjacent
             if (dosimpel)
             {
                 double fc = ChannelWidth->Drc/_dx;
                 // fraction of the channel in the gridcell, 1-fc = (dx-chw)/dx = chanadj/dx
                 double whlevel = (ChannelWH->Drc-chdepth)*fc + H*(1-fc);
                 // equilibrium water level = weighed values of channel surplus level + _h
                 // can be negative if channelwh is below channel depth and low _h level
                 if(whlevel > 0)
                 {
                     double sedch = 0;
                     double sed = 0;
                     if (SwitchErosion) {
                         sedch = ChannelSSSed->Drc;
                         sed = SSFlood->Drc;
                     }
                     double oldchwh = ChannelWH->Drc;
                     double oldwh = H;
                     ChannelWH->Drc = whlevel + chdepth;
                     _h->Drc = whlevel;

                     // new equilibrium levels
                     if (SwitchErosion) {
                         // double sed_ = SSFlood->Drc + ChannelSSSed->Drc;
                         if (oldchwh > ChannelWH->Drc) {
                             double sed = (oldchwh-ChannelWH->Drc)*ChannelWidth->Drc*ChannelDX->Drc * ChannelSSConc->Drc;
                             ChannelSSSed->Drc -= sed;
                             SSFlood->Drc += sed;
                         } else {
                             double sed = (oldwh-_h->Drc)*CHAdjDX->Drc * SSCFlood->Drc;
                             SSFlood->Drc -=sed;
                             ChannelSSSed->Drc += sed;
                         }
                     }

                 }
                 else
                 {
                     ChannelWH->Drc += _h->Drc*CHAdjDX->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
                     _h->Drc = 0;
                     //DO NOTHING
                     // this happens if there is very little flood water (< 5cm) and the channelWH is below the channeldepth
                     // we assume that there is no more flow towards the channel.
                 }
             }
         }
     }}

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (ChannelMaxQ->Drc <= 0) {
            ChannelWaterVol->Drc = ChannelWH->Drc * ChannelDX->Drc * ChannelWidth->Drc;
            WaterVolall->Drc = CHAdjDX->Drc*_h->Drc + MicroStoreVol->Drc;

            if(SwitchErosion) {
                SWOFSedimentLayerDepth(r,c,_h->Drc, V->Drc);
                SWOFSedimentSetConcentration(r,c, _h->Drc, ChannelAdj->Drc);

                RiverSedimentLayerDepth(r, c);
                RiverSedimentMaxC(r, c);
                // all concentrations, possible ChannelDep when surplus
            }
        }
    }}
}

//---------------------------------------------------------------------------
// flow to and from channel based on broad crested weirs, freeflow or drowned
// TUFLOW and other models use this
// www.brighthubengineering.com

void TWorld::ChannelOverflowAlt(cTMap *_h, cTMap *V)
{
    if (!SwitchIncludeChannel)
         return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (ChannelMaxQ->Drc == 0) {
            double dCHh = ChannelWH->Drc-ChannelDepth->Drc;
            double H = _h->Drc;

            if (H < 1e-6 && dCHh < 0)
                continue; // nthing to flow

            if (fabs(H-dCHh) < 1e-6)
                continue; // no flow, already equilibrium

            double area_channel = ChannelWidth->Drc * ChannelDX->Drc;
            double area_surface = CHAdjDX->Drc;
            double needed_volume = 0;

            bool tochannel = true;
            double transfer_volume = 0;
            double Cd = 0.56; // 2/3 * 0.86
            double factor= 2*_dt*ChannelDX->Drc;
            //factor 2 is for flow on both sides of the channel over length DX

            //double H_eq = (dCHh*area_channel + H*area_surface)/CellArea->Drc;
            double H_eq= dCHh*ChannelWidth->Drc/_dx + H*(1-ChannelWidth->Drc/_dx);
            // equilibrium level
          //  qDebug() << H_eq << dCHh << H;

//            if (H_eq < 0) {
            if (dCHh < 0) {
                needed_volume = (H-std::max(0.0, H_eq))*area_surface;
                // potentially all surface water flows into channel

                double freeflow_tochan = factor*Cd*SQRT2G*std::pow(H,1.5);
                //free flow broad crested weir:

                transfer_volume = std::min(freeflow_tochan, needed_volume);
            } else {
                //dCHh >> 0
                if (H > dCHh) {
                    // surface water higher than channel, eq of drowned broad crested weir

                    needed_volume = (H - std::max(0.0,H_eq))*area_surface;
                    // vol needed to reach equilibrium level

                    double freeflow_tochan = factor*Cd*SQRT2G*std::pow(H,1.5);
                    //free flow broad crested weir:

                    double transfer_volume_tochan = freeflow_tochan;
                    if (H > 1e-6)
                        transfer_volume_tochan *= std::sqrt(1-std::pow((H-dCHh)/H,1.5));
                    //eq. broad crested weir: drowned flow

                    transfer_volume_tochan = std::max(transfer_volume_tochan, factor*V->Drc);
                    //if surface velocity is higher take that

                    transfer_volume = qMin(transfer_volume_tochan, needed_volume);
                } else {
                    // flow from channel, drowned weir
                    needed_volume = (dCHh - std::max(0.0,H_eq))*area_channel;
                    // vol needed to reach equilibrium level

                    // broad crested weir flow if channel is leadng
                    double freeflow_fromchan = factor*Cd*SQRT2G*std::pow(dCHh, 1.5);

                    double transfer_volume_fromchan = freeflow_fromchan;
                    if (dCHh > 1e-6)
                        transfer_volume_fromchan *= sqrt(1-std::pow((dCHh-H)/dCHh,1.5));
                    //eq. broad crested weir: drowned flow

                    transfer_volume = std::min(transfer_volume_fromchan, needed_volume);
                    tochannel = false;
                }
            }

            if (tochannel) {
                WaterVolall->Drc -= transfer_volume;
                ChannelWaterVol->Drc += transfer_volume;
            } else {
                WaterVolall->Drc += transfer_volume;
                ChannelWaterVol->Drc -= transfer_volume;
            }

            // Update heights
            ChannelWH->Drc = ChannelWaterVol->Drc / area_channel;
            _h->Drc = (WaterVolall->Drc-MicroStoreVol->Drc) / area_surface;

            // new equilibrium levels
            if (SwitchErosion) {
                if (tochannel) {
                    double sed = transfer_volume * SSCFlood->Drc;
                    SSFlood->Drc -=sed;
                    ChannelSSSed->Drc += sed;
                } else {
                    double sed = transfer_volume * ChannelSSConc->Drc;
                    ChannelSSSed->Drc -= sed;
                    SSFlood->Drc += sed;
                }
                SWOFSedimentLayerDepth(r,c,_h->Drc, V->Drc);
                SWOFSedimentSetConcentration(r,c, _h->Drc, ChannelAdj->Drc);

                RiverSedimentLayerDepth(r, c);
                RiverSedimentMaxC(r, c);
            }
        }
    }}
}
//---------------------------------------------------------------------------
/**
 * @fn void TWorld::ToFlood(void)
 * @brief Calculates overland flow that flows into flooding water
 *
 * Calculates overland flow of water and sediment that flows into flooding water
 * based on the runoff partitioning factor. Depending on the parameter, water
 * is either transformed quickly or slowly. This imitates the effect that overland
 * flow would have on the velocity of the flood water.
 *
 * @return void
 * @see runoff_partitioning
 */
void TWorld::ToFlood()
{
    #pragma omp parallel for  num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (hmx->Drc > HMIN && WHrunoff->Drc > HMIN) {
            double frac = 1.0;
            double dwh = frac * WHrunoff->Drc;

            hmx->Drc += dwh;
            WH->Drc = WHstore->Drc;
            WHrunoff->Drc = 0;

            hmxWH->Drc = hmx->Drc + WH->Drc;
            WaterVolall->Drc = CHAdjDX->Drc*(WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;

            if(SwitchErosion) {
                double dsed = frac*Sed->Drc;
                SSFlood->Drc += dsed;
                Sed->Drc = 0;
                Conc->Drc = 0;

                SWOFSedimentLayerDepth(r,c,hmx->Drc, V->Drc);
                SWOFSedimentSetConcentration(r,c,hmx->Drc, ChannelAdj->Drc);
               // Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
            }
        }
    }}
}
//---------------------------------------------------------------------------
// DO NOT MAKE PARALLEL
void TWorld::FloodMaxandTiming()
{
    // floodwater volume and max flood map
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (hmxWH->Drc > minReportFloodHeight) {
            floodTime->Drc += _dt/60;
            floodHmxMax->Drc = std::max(floodHmxMax->Drc, hmxWH->Drc);
            // for output
        }

        floodVMax->Drc = std::max(floodVMax->Drc, V->Drc);
        floodVHMax->Drc = std::max(floodVHMax->Drc, V->Drc*hmxWH->Drc);
        // max velocity
        WHmax->Drc = std::max(WHmax->Drc, hmxWH->Drc);
    }}
    floodVolTotMax = 0;
    floodArea = 0;
    double area = _dx*_dx;

   // #pragma omp parallel for reduction(+:floodVolTotMax,floodArea) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (floodHmxMax->Drc > minReportFloodHeight) {
            floodVolTotMax += floodHmxMax->Drc*area;
        }
        if (hmxWH->Drc > minReportFloodHeight && floodTimeStart->Drc == 0)  {
            if (SwitchWaveUser)
                floodTimeStart->Drc = (time - BeginTime)/60.0;
            else
                floodTimeStart->Drc = (time - RainstartTime)/60.0;
            // time since first pixel received rainfall
        }
        if (hmxWH->Drc > minReportFloodHeight) {
            floodArea += area;
        }
    }}

    floodAreaMax = std::max(floodArea,floodAreaMax);
}
//---------------------------------------------------------------------------
// NOTE THIS function is only called for Kinematic+dynamic wave
void TWorld::ChannelFlood(void)
{

    if (!SwitchIncludeChannel)
        return;

    ToFlood();
    // mix HWrunoff with hmx
    // if toflood before channeloverflow then MB error in sed

    if (SwitchChannel2DflowConnect)
        ChannelOverflowAlt(hmx, V);
    else
        ChannelOverflow(hmx, V);
    // determine overflow water => hmx
    // hmx is flood water, WH is overlandflow, WHrunoff etc

    startFlood = false;
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (hmx->Drc > 0)
            startFlood = true;
    }}

    double dtflood = fullSWOF2openMUSCL(hmx, Uflood, Vflood, DEM);
    // in kindyn hmx is the channel overflow/flood part of the surface water, the rest is kinwave WHrunoff

    //new flood domain
    nrFloodedCells = 0;
    FOR_ROW_COL_MV {
        if (hmx->Drc > 0) {
            FloodDomain->Drc = 1;
            nrFloodedCells += 1.0;
        }
        else
            FloodDomain->Drc = 0;
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Qflood->Drc = 0;
        if (FloodDomain->Drc > 0) {
            V->Drc = sqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
            //Qflood->Drc = V->Drc * hmx->Drc * ChannelAdj->Drc;
            Qn->Drc = V->Drc * hmx->Drc * ChannelAdj->Drc;//0;//V->Drc * WHrunoff->Drc * ChannelAdj->Drc;
            // Qn is the runoff water, must be zero in the flooded area, becomes Qflood
        }
    }}

 #pragma omp parallel for num_threads(userCores)
 FOR_ROW_COL_MV_L {
     WH->Drc = WHrunoff->Drc+ WHstore->Drc;
     // add new average waterlevel (A/dx) to stored water

     WaterVolall->Drc = CHAdjDX->Drc*(WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;

     hmxWH->Drc = WH->Drc + hmx->Drc;
     // all water on surface

     hmxflood->Drc = std::max(0.0, WHrunoff->Drc + hmx->Drc - minReportFloodHeight);

     FloodWaterVol->Drc = hmxflood->Drc * CHAdjDX->Drc;
     double WHrunoffOutput = std::min(WHrunoff->Drc + hmx->Drc, minReportFloodHeight);
     RunoffWaterVol->Drc = WHrunoffOutput * CHAdjDX->Drc;
     // these are only used for reporting totals on screen and in file

     if(SwitchErosion) {
         Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
         if (FloodDomain->Drc  > 0) {
             double sed = SSFlood->Drc + BLFlood->Drc;
             Conc->Drc =  MaxConcentration(FloodWaterVol->Drc, sed);
             Qsn->Drc += Conc->Drc*Qn->Drc;//flood->Drc;
         }
     }
  }}
/*
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Qflood->Drc = 0;
        if (FloodDomain->Drc > 0) {
            V->Drc = sqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
            Qflood->Drc = V->Drc * hmx->Drc * ChannelAdj->Drc;
            Qn->Drc = 0;//V->Drc * WHrunoff->Drc * ChannelAdj->Drc;
            // Qn is the runoff water, must be zero in the flooded area, becomes Qflood
        }

        WH->Drc = WHrunoff->Drc+ WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        WaterVolall->Drc = CHAdjDX->Drc*(WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;

        hmxWH->Drc = WH->Drc + hmx->Drc;
        // all water on surface

        hmxflood->Drc = std::max(0.0, WHrunoff->Drc + hmx->Drc - minReportFloodHeight);

        FloodWaterVol->Drc = hmxflood->Drc * CHAdjDX->Drc;
        double WHrunoffOutput = std::min(WHrunoff->Drc + hmx->Drc, minReportFloodHeight);
        RunoffWaterVol->Drc = WHrunoffOutput * CHAdjDX->Drc;
        // these are only used for reporting totals on screen and in file
    }}

    if(SwitchErosion) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
            if (FloodDomain->Drc  > 0) {
                double sed = SSFlood->Drc + BLFlood->Drc;
                Conc->Drc = MaxConcentration(FloodWaterVol->Drc, sed);
                Qsn->Drc += Conc->Drc*Qflood->Drc;
            }
        }}
    }
*/
    FloodMaxandTiming();

    double area = nrFloodedCells*_dx*_dx;
    if (area > 0)
        debug(QString("Flooding (dt %1 sec, n %2): area %3 m2, %4 cells").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1).arg(nrFloodedCells));//.arg(K2DQOutBoun));
    // some screen error reporting

}
