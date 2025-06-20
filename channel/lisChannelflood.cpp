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
- void TWorld::ChannelFlood(void) Calculate channelflood height maps (hmx, U+VFlood) and FloodDomain
*/

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

void TWorld::ChannelOverflow(cTMap *_h, cTMap *V)
{
    if (!SwitchIncludeChannel)
         return;

   #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (!crch_[i_].culvert) {

            double dH = std::max(0.0, (ChannelWH->Drc-ChannelDepth->Drc)); // water higher than channel depth
            double H = _h->Drc;

            if (dH <= 1e-6 && H <= 1e-6)
                continue;
            // no flow activity then continue

            if (fabs(dH - H) < 1e-6)
                continue;
            // no diff in water level, no flow, continue

            double Cd = 0.56;
            double cwa = ChannelWidth->Drc/ChannelAdj->Drc;
            bool dosimpel = false;

            if (dH > H) {
                // flow from channel
                double VfromChan = Cd*SQRT2G*pow(dH-H,1.5)/H;
                double fracC = std::min(1.0, _dt*VfromChan/(0.5*ChannelAdj->Drc));
                // fraction from channel to surrounding
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
                double VtoChan = V->Drc;//=Cd*SQRT2G*pow(H-dH,1.5)/H;
                double fracA = std::min(1.0, _dt*VtoChan/(0.5*ChannelAdj->Drc));
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
            if (dosimpel) {
                double fc = ChannelWidth->Drc/_dx;
                // fraction of the channel in the gridcell, 1-fc = (dx-chw)/dx = chanadj/dx
                double whlevel = (ChannelWH->Drc-ChannelDepth->Drc)*fc + H*(1-fc);
                // equilibrium water level = weighed values of channel surplus level + _h
                // can be negative if channelwh is below channel depth and low _h level
                if(whlevel > 0) {
                    double sedch = 0;
                    double sed = 0;
                    if (SwitchErosion) {
                        sedch = ChannelSSSed->Drc;
                        sed = SSFlood->Drc;
                    }
                    double oldchwh = ChannelWH->Drc;
                    double oldwh = H;
                    ChannelWH->Drc = whlevel + ChannelDepth->Drc;
                    _h->Drc = whlevel;

                    // new equilibrium levels
                    if (SwitchErosion) {
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

                } else {
                    // NB: this gives a larger mass balance error!

                    // assume everything flows into the channel
                    ChannelWH->Drc += _h->Drc*CHAdjDX->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
                    _h->Drc = 0;

                    // this happens if there is very little flood water (< 5cm) and the channelWH is below the channeldepth
                    // we assume that there is no more flow towards the channel.
                    if (SwitchErosion) {
                        ChannelSSSed->Drc += SSFlood->Drc;
                        SSFlood->Drc = 0;
                    }
                }
            } // dosimnpel

            ChannelWaterVol->Drc = ChannelWH->Drc * ChannelDX->Drc * ChannelWidth->Drc;

        }
    }}

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (!crch_[i_].culvert && _h->Drc > 0) {
            if (SwitchKinematic2D == K2D_METHOD_KINDYN) {
                hmx->Drc = _h->Drc + WHstore->Drc;
                hmxWH->Drc = hmx->Drc;
            } else {
                WH->Drc = _h->Drc + WHstore->Drc;
                hmxWH->Drc = WH->Drc;
            }

            WaterVolall->Drc = CHAdjDX->Drc*_h->Drc + MicroStoreVol->Drc;

            if(SwitchErosion) {
                SWOFSedimentLayerDepth(r,c, _h->Drc, V->Drc);
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
// NOTE _h is WHrunoff so without microdepression storage
void TWorld::ChannelOverflowAlt(cTMap *_h, cTMap *V)
{
    if (!SwitchIncludeChannel)
         return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (!crch_[i_].culvert) {

            switch (crch_[i_].shape) {
                case SHAPERECT : chanHandPRect(r,c); break;
                case SHAPECIRC : chanHandPCirc(r,c); break; // this is always a culvert!
                case SHAPETRAP : chanHandPTrap(r,c); break;
                case SHAPETRIA : chanHandPTria(r,c); break;
            }

            double dCHh = ChannelWH->Drc-ChannelDepth->Drc;
            double dCHh0 = std::max(dCHh, 0.0);
            double H = _h->Drc; // runoff height!

            if (H < 1e-6 && dCHh0 < 1e-6)
                continue; // nothing to flow

            if (fabs(H-dCHh0) < 1e-6)
                continue; // no flow, already equilibrium

            double area_channel = ChannelWidth->Drc * ChannelDX->Drc;
            double area_surface = CHAdjDX->Drc;
            double needed_volume = 0;

            bool tochannel = true;
            double transfer_volume = 0;
            double Cd = 0.56; // 2/3 * 0.86
            double lengthfactor = 2.0*_dt*ChannelDX->Drc;
            double velocityfactor = V->Drc*V->Drc/(2*GRAV);
            //do not use factor 2 for flow on both sides


            if(dCHh < 0){
                double negvol = ChannelMaxArea->Drc*ChannelDX->Drc - ChannelWaterVol->Drc;
                double freeflow_tochan = lengthfactor*Cd*SQRT2G*std::pow(H+velocityfactor,1.5);

                needed_volume = H*area_surface;
                // if flow fills up channel create equilibrium level
                if (transfer_volume > negvol) {
                    double heq = (transfer_volume-negvol)/CellArea->Drc;
                    // equilibrium level
                    needed_volume = negvol + (H-heq)*CellArea->Drc;
                    // transfer_volume = vol needed for equilibrium level
                }
                transfer_volume = std::min(freeflow_tochan, needed_volume);
                //m3 free flow broad crested weir, water flows over edge to deeper water in channel
                tochannel = true;
            } else {
                // chhannel water is bankfull or more, channelwatervolo has already shape
                // because higher. dCHh0 always refers to rectangle above surface with channelwidth
                double H_eq = (dCHh*area_channel + H*area_surface)/CellArea->Drc;
                // equilibrium level
                if (H > dCHh0) {
                    // surface water higher than channel water, drowned broad crested weir
                    needed_volume = (H - H_eq)*area_surface;
                    // vol needed to reach equilibrium level
                    double transfer_volume_tochan = lengthfactor*Cd*SQRT2G*std::pow(H+velocityfactor - dCHh0,1.5);
                    // drowned flow to channel with velocity of approach
                    transfer_volume = std::min(transfer_volume_tochan, needed_volume);
                    tochannel = true;
                } else {
                    // flow from channel, drowned weir in the other dircetion, no added velocity
                    needed_volume = (dCHh - H_eq)*area_channel;
                    // vol needed to reach equilibrium level
                    double transfer_volume_fromchan = lengthfactor*Cd*SQRT2G*std::pow(dCHh0 - H,1.5);
                    // drowned flow from channel
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

            //Update water height from volume
            switch (crch_[i_].shape) {
                case SHAPERECT : chanHandPRect(r,c); break;
                case SHAPECIRC : chanHandPCirc(r,c); break; // this is always a culvert!
                case SHAPETRAP : chanHandPTrap(r,c); break;
                case SHAPETRIA : chanHandPTria(r,c); break;
            }
            // update surface water height
            _h->Drc = std::max(0.0, WaterVolall->Drc-MicroStoreVol->Drc) / area_surface;

            if (SwitchKinematic2D == K2D_METHOD_KINDYN) {
                hmx->Drc = WaterVolall->Drc/area_surface;
                hmxWH->Drc = hmx->Drc;
            } else {
                WH->Drc =  WaterVolall->Drc/area_surface; ///_h->Drc + WHstore->Drc;
                hmxWH->Drc = WH->Drc;
            }

            // new equilibrium levels erosion
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
// NOTE THIS function is only called for Kinematic+dynamic wave
void TWorld::ToFlood()
{
    #pragma omp parallel for  num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (hmxrunoff->Drc > 0 && WHrunoff->Drc > 0) {
            double dwh = WHrunoff->Drc;

            hmxrunoff->Drc += dwh;
            hmx->Drc = hmxrunoff->Drc + WHstore->Drc;
            WHrunoff->Drc = 0;
            WH->Drc = WHstore->Drc;

            hmxWH->Drc = hmx->Drc + WH->Drc;
            WaterVolall->Drc = CHAdjDX->Drc*hmxWH->Drc;

            if(SwitchErosion) {
                double dsed = Sed->Drc;
                SSFlood->Drc += dsed;
                Sed->Drc = 0;
                Conc->Drc = 0;

                SWOFSedimentLayerDepth(r,c,hmx->Drc, V->Drc);
                //SWOFSedimentSetConcentration(r,c,hmx->Drc, ChannelAdj->Drc);
                SSCFlood->Drc = MaxConcentration(WaterVolall->Drc, SSFlood->Drc);
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

    #pragma omp parallel for reduction(+:floodVolTotMax,floodArea) num_threads(userCores)
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
    // hmx = flood equivalent of WH; hmxrunoff of WHrunoff

    if (!SwitchIncludeChannel)
        return;

    ToFlood();

    if (SwitchChannel2DflowConnect)
        ChannelOverflowAlt(hmxrunoff, V);
    else
        ChannelOverflow(hmxrunoff, V);
    // determine overflow water => hmx
    // hmx is flood water, WH is overlandflow, WHrunoff etc

    startFlood = false;
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (hmxrunoff->Drc > 0)
            startFlood = true;
    }}

    double dtflood = 0;
    if (startFlood)
        dtflood = fullSWOF2openMUSCL(hmxrunoff, Uflood, Vflood, DEM);

    //new flood domain
    nrFloodedCells = 0;
    FOR_ROW_COL_MV {
        if (hmxrunoff->Drc > 0) {
            FloodDomain->Drc = 1;
            nrFloodedCells += 1.0;
        }
        else
            FloodDomain->Drc = 0;
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (FloodDomain->Drc > 0) {
            V->Drc = sqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
            Qn->Drc = V->Drc * hmxrunoff->Drc * ChannelAdj->Drc;
        }
    }}

    updateWHandHmx();

    FloodMaxandTiming();

    double area = nrFloodedCells*_dx*_dx;
    if (area > 0)
        debug(QString("Flooding (dt %1 sec, n %2): area %3 m2, %4 cells").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1).arg(nrFloodedCells));//.arg(K2DQOutBoun));
    // some screen error reporting

}
