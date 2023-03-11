/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
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


//---------------------------------------------------------------------------
//! Get flood level in channel from 1D kin wave channel
//! Instantaneous mixing of flood water and channel water in channel cells
//! note: ChannelDepth lets you also control which channels flood:
//! those that are 0 react as usual (infinite capacity)
//doOF is always false!!!
void TWorld::ChannelOverflow(cTMap *_h, cTMap *V)
{
    if (!SwitchIncludeChannel)
        return;

    //#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (SwitchErosion) {
            SWOFSedimentSetConcentration(r,c, _h);
            RiverSedimentMaxC(r, c);
        }

        if (ChannelWidth->Drc > 0 && ChannelMaxQ->Drc <= 0)
        {
            double chdepth = ChannelDepth->Drc;
            // double charea = ChannelWidth->Drc*ChannelDX->Drc;
            double dH = std::max(0.0, (ChannelWH->Drc-chdepth));

            if (dH <= HMIN && _h->Drc <= HMIN)
                continue;
            // no flow activity then continue

            if (dH == _h->Drc)
                continue;
            // no diff in water level, no flow, continue

            // VELOCITIES
            double VtoChan = V->Drc;
            if (SwitchGravityToChannel)
                VtoChan = std::pow(_h->Drc, 2.0/3.0)*sqrt(ChannelPAngle->Drc)/N->Drc;
            double fracA = std::min(1.0, _dt*VtoChan/(0.5*ChannelAdj->Drc));
            // fraction from _h to channel based on average flood velocity

            double VfromChan = sqrt(2*9.804*dH); //Bernoulli
            //see https://www.engineeringtoolbox.com/velocity-head-d_916.html
            double fracC = std::min(1.0, _dt*VfromChan/(0.5*ChannelAdj->Drc));
            // fraction from channel to surrounding

            double cwa = ChannelWidth->Drc/ChannelAdj->Drc;

            bool dosimpel = false;

            if (dH > _h->Drc)   // flow from channel
            {
                double dwh = fracC * dH;
                // amount flowing from channel
                if (_h->Drc + dwh*cwa > dH-dwh) {
                    // if flow causes situation to reverse (channel dips below _h)
                    dosimpel = true;
                } else {

                    _h->Drc += dwh*cwa;
                    ChannelWH->Drc -= dwh;

                    if(SwitchErosion) {
                        double sed = ChannelWH->Drc > 0 ? dwh/ChannelWH->Drc * ChannelSSSed->Drc : 0;
                        ChannelSSSed->Drc -= sed;
                        SSFlood->Drc += sed;
                    }
                }
            }
            else   // flow to channel
            {
                double dwh = fracA * _h->Drc;
                // amount flowing to channel
                if (dH + dwh/cwa > _h->Drc-dwh) {
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
                double whlevel = (ChannelWH->Drc-chdepth)*fc + _h->Drc*(1-fc);
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
                    double oldwh = _h->Drc;
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
//                        double sed_2 = SSFlood->Drc + ChannelSSSed->Drc;
//                        qDebug() << sed_ << sed_2 << sed_-sed_2;

//                        RiverSedimentLayerDepth(r,c);
//                        double _sed = ChannelSSSed->Drc + SSFlood->Drc;
//                        double volch = ChannelSSDepth->Drc*ChannelWidth->Drc*ChannelDX->Drc;
//                        double volof = _h->Drc*CHAdjDX->Drc;
//                        double _concavg = _sed/(volch+volof);

//                        ChannelSSSed->Drc = _concavg * volch;
//                        ChannelSed->Drc = (SwitchUse2Phase ? ChannelBLSed->Drc : 0.0) + ChannelSSSed->Drc;

//                        SSFlood->Drc = _concavg * volof;
//                        SSCFlood->Drc = _concavg;
                    }

                }
                else
                {
                    //DO NOTHING
                    // this happens if there is very little flood water (< 5cm) and the channelWH is below the channeldepth
                    // we assume that there is no more flow towards the channel.
                }
            }

            ChannelWaterVol->Drc = ChannelWH->Drc * ChannelDX->Drc * ChannelWidth->Drc;
            // WaterVolall->Drc = _h->Drc*ChannelAdj->Drc*DX->Drc;
            // do not recalc floodvol, MB errors

            // recalc channel water vol else big MB error
            if(SwitchErosion)
            {
                SWOFSedimentLayerDepth(r,c,_h->Drc, V->Drc);
                SWOFSedimentSetConcentration(r,c, _h);

                RiverSedimentLayerDepth(r,c);
                RiverSedimentMaxC(r, c);
                // all concentrations, possible ChannelDep when surplus

            }

        }
    }}
}
//---------------------------------------------------------------------------

// NOT USED!

//! Get flood level in channel from 1D kin wave channel
//! Instantaneous mixing of flood water and channel water in channel cells
//! note: ChannelDepth lets you also control which channels flood:
//! those that are 0 react as usual (infinite capacity)
//!
void TWorld::ChannelOverflowNew(cTMap *_h, cTMap *V, bool doOF)
{
    if (!SwitchIncludeChannel)
        return;

    cTMap *_SS;
    cTMap *_SSC;
    if(SwitchErosion) {
        if (doOF) {
            // obsolete, for when this function was used instead of tochannel
            _SS = Sed;
            _SSC = Conc;
        } else {
            _SS = SSFlood;
            _SSC = SSCFlood;
        }
    }

    FOR_ROW_COL_MV_CH {
        if(ChannelMaskExtended->data[r][c] == 1)// && !pcr::isMV(LDDChannel->data[r][c]))
        {
            int rr = r;//(int)ChannelSourceYExtended->Drc;
            int cr = c;//(int)ChannelSourceXExtended->Drc;

            if(SwitchErosion) {
                if (doOF)
                    Conc->Drcr = MaxConcentration(CHAdjDX->Drc*_h->Drcr, &Sed->Drc, &DEP->Drc);
                else
                    SWOFSedimentSetConcentration(rr,cr, _h);

                RiverSedimentMaxC(rr, cr);
            }

            if (ChannelDepth->Drcr > 0 && ChannelMaxQ->Drcr <= 0)
            {
                double chdepth = ChannelDepth->Drcr;
                // double charea = ChannelWidth->Drcr*ChannelDX->Drcr;
                double dH = std::max(0.0, (ChannelWH->Drcr-chdepth));

                if (dH <= HMIN && _h->Drcr <= HMIN)
                    continue;
                // no flow activity then continue

                if (dH == _h->Drcr)
                    continue;
                // no diff in water level, no flow, continue

                // VELOCITIES
                double VtoChan = V->Drc;
                if (F_AddGravity == 1)
                    VtoChan = std::pow(_h->Drcr, 2.0/3.0)*sqrt(ChannelPAngle->Drc)/N->Drcr; //F_Angle
                double fracA = std::min(1.0, _dt*VtoChan/(0.5*_dx));
                // fraction from _h to channel based on average flood velocity

                double VfromChan = sqrt(2*9.81*dH); //Bernoulli
                double fracC = std::min(1.0, _dt*VfromChan/(0.5*_dx));
                // fraction from channel to surrounding based on overflow height and manning

                double cwa = ChannelAdj->Drc > 0 ? ChannelWidthO->Drcr/ChannelAdj->Drc : 0;

                bool dosimpel = false;

                if (dH > _h->Drcr)   // flow from channel
                {
                    double dwh = fracC * dH;
                    if (_h->Drcr + dwh*cwa > dH-dwh) {
                        // if flow causes situation to reverse (channel dips below _h)
                        dosimpel = true;
                    } else {

                        _h->Drcr += dwh*cwa;
                        ChannelWH->Drcr -= dwh;

                        if(SwitchErosion) {
                            double sed = dwh/ChannelWH->Drc * ChannelSSSed->Drcr;//
                            ChannelSSSed->Drcr -= sed;
                            _SS->Drcr += sed;
                        }
                    }
                }
                else   // flow to channel
                {
                    // water above channel so rectangular
                //    if (dH > 0){
                        double dwh = fracA * _h->Drcr;
                        if (dH + dwh/cwa > _h->Drcr-dwh) {
                            // if too much flow
                            dosimpel = true;
                        } else {
                            _h->Drcr -= dwh;
                            ChannelWH->Drcr += (dwh/cwa);
                            if(SwitchErosion) {
                                double sed = fracA*_SS->Drcr;
                                ChannelSSSed->Drcr += sed;
                                _SS->Drcr -= sed;
                            }
                        }
//                    }
//                    else
//                    {
//                        // water below channel so can have side angle
//                        double dwh = fracA * _h->Drcr;
//                        double dvol = dwh*DX->Drcr*ChannelAdj->Drcr;
//                        double cwh = channelVoltoWH(ChannelWaterVol->Drcr+dvol,rr, cr);

//                        if (cwh-chdepth > _h->Drcr-dwh) {
//                            // if too much flow
//                            dosimpel = true;
//                        } else {
//                            _h->Drcr -= dwh;
//                            ChannelWaterVol->Drcr += dvol;
//                            fromChannelVoltoWH(rr, cr);

//                            if(SwitchErosion) {
//                                double sed = fracA*_SS->Drcr;
//                                ChannelSSSed->Drcr += sed;
//                                _SS->Drcr -= sed;
//                            }
//                        }
//                    }
                }

                // instantaneous waterlevel exquilibrium acccross channel and adjacent
                if (dosimpel)
                {
                    double fc = ChannelWidthO->Drcr/_dx;
//                    double totvol = ChannelWaterVol->Drc + _h->Drc*ChannelAdj->Drc*DX->Drc;
//                    double chvol = totvol*fc;
//                    double ofvol = totvol*(1-fc);
//                    _h->Drc = ofvol/(ChannelAdj->Drc*DX->Drc);
//                    ChannelWH->Drc = channelVoltoWH(chvol,rr, cr);

                    // fraction of the channel in the gridcell, 1-fc = (dx-chw)/dx = chanadj/dx
                    double whlevel = (ChannelWH->Drcr-chdepth)*fc + _h->Drcr*(1-fc);
                    // equilibrium water level = weighed values of channel surplus level + _h
                    // can be negative if channelwh is below channel depth and low _h level
                    if(whlevel > HMIN) {
                        ChannelWH->Drcr = whlevel + chdepth;
                        _h->Drcr = whlevel;

                        // new equilibrium levels
                        if(SwitchErosion)
                        {

                            RiverSedimentLayerDepth(rr,cr);
                            //SWOFSedimentLayerDepth(rr, cr, _h->Drcr, V->Drcr);

                            double _sed = ChannelSSSed->Drcr + _SS->Drcr;
                            double volch = ChannelSSDepth->Drcr*ChannelWidth->Drcr*ChannelDX->Drcr;
                            double volof = _h->Drcr*ChannelAdj->Drcr*DX->Drcr;
                            double _concavg = _sed/(volch+volof);

                            ChannelSSSed->Drc = _concavg * volch;
                            ChannelSed->Drc = ChannelBLSed->Drc + ChannelSSSed->Drc;

                            _SS->Drcr = _concavg * volof;
                            _SSC->Drcr = _concavg;
                        }

                    }
                    else
                    {
                        //DO NOTHING
                        // this happens if there is very little flood water (< 5cm) and the channelWH is below the channeldepth
                        // we assume that there is no more flow towards the channel.
                    }

                }

                //fromChannelWHtoVol(rr, cr);
                // do not recalc floodvol, MB errors
                ChannelWaterVol->Drcr = ChannelWidth->Drcr * ChannelWH->Drcr * ChannelDX->Drcr;

                // recalc channel water vol else big MB error
                if(SwitchErosion) {
                    if (doOF)
                        Conc->Drcr = MaxConcentration(CHAdjDX->Drc*_h->Drcr, &_SS->Drc, &DEP->Drc);
                    else {
                        SWOFSedimentLayerDepth(r,c,_h->Drcr, V->Drcr);
                        SWOFSedimentSetConcentration(rr,cr, _h);
                    }

                    RiverSedimentLayerDepth(rr,cr);
                    RiverSedimentMaxC(rr, cr);
                    // all concentrations, possible ChannelDep when surplus

                }

            }
        }
    }

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
        if(hmx->Drc > 0 && WHrunoff->Drc > 0)
        {
            double frac = 1-exp(-runoff_partitioning*hmx->Drc/(WHrunoff->Drc+0.001));

//            double V = qSqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
//            double Vkin = pow(WHrunoff->Drc,2.0/3.0)*sqrt(Grad->Drc)/N->Drc;
//            frac = Vkin/(V+0.001);
            frac = std::max(std::min(frac, 1.0),0.0);
            double dwh = frac * WHrunoff->Drc;

            hmx->Drc += dwh;
            WH->Drc -= dwh;
            WHrunoff->Drc -= dwh;
            WHroad->Drc -= dwh;
            hmxWH->Drc = hmx->Drc + WH->Drc;

            if(SwitchErosion)
            {
                double dsed = frac*Sed->Drc;
                SSFlood->Drc += dsed;
                Sed->Drc -= dsed;

                SWOFSedimentSetConcentration(r,c,hmx);
                Conc->Drc = MaxConcentration(WH->Drc*CHAdjDX->Drc, &Sed->Drc, &DEP->Drc);
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
// NOTE DEM has barriers included, done in shade map calculation !!!!
void TWorld::ChannelFlood(void)
{
    if (!SwitchIncludeChannel)
        return;

    ChannelOverflow(hmx, V);
    // determine overflow water => hmx

    ToFlood();
    // mix WHrunoff with hmx

    double dtflood = 0;

    startFlood = false;
    FOR_ROW_COL_MV {
        if (hmx->Drc > 0) {
            startFlood = true;
            break;
        }
    }

    if (SwitchSWOFopen)
        dtflood = fullSWOF2open(hmx, Uflood, Vflood, DEM);
    else
        dtflood = fullSWOF2RO(hmx, Uflood, Vflood, DEM);
    // 2D dyn flow of hmx water, when using kin wave

    //new flood domain
    nrFloodedCells = 0;
    FOR_ROW_COL_MV {
        if (hmx->Drc > 0)
        {
            FloodDomain->Drc = 1;
            nrFloodedCells += 1.0;
        }
        else
            FloodDomain->Drc = 0;
    }
qDebug() << nrFloodedCells;
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Qflood->Drc = 0;
        if (FloodDomain->Drc > 0) {
            V->Drc = qSqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
            Qflood->Drc = V->Drc * hmx->Drc * ChannelAdj->Drc;
        }
    }}

    Boundary2Ddyn();
    // boundary flow

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        hmxWH->Drc = WH->Drc + hmx->Drc;

        //InfilVolFlood->Drc += Iflood->Drc;
        // addvolume infiltrated during flood process with FSurplus

        hmxflood->Drc = std::max(0.0, WHrunoff->Drc + hmx->Drc - minReportFloodHeight);

        WaterVolall->Drc = CHAdjDX->Drc*(WHrunoff->Drc + hmx->Drc) + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        // all water on surface

        FloodWaterVol->Drc = hmxflood->Drc * CHAdjDX->Drc;
        double WHrunoffOutput = std::min(WHrunoff->Drc + hmx->Drc, minReportFloodHeight);
        RunoffWaterVol->Drc = WHrunoffOutput * CHAdjDX->Drc;

    }}

    if(SwitchErosion) {
        //calculate concentration and new sediment discharge
        //WHrunoff and Qn are adapted in case of 2D routing
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (FloodDomain->Drc  > 0) {
                double sed = SSFlood->Drc + BLFlood->Drc;
                Conc->Drc =  MaxConcentration(FloodWaterVol->Drc, &sed, &DepFlood->Drc);
                Qsn->Drc += Conc->Drc*Qflood->Drc;
            }
        }}
     }

    double area = nrFloodedCells*_dx*_dx;
    if (area > 0)
        debug(QString("Flooding (dt %1 sec, n %2): area %3 m2, %4 cells").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1).arg(nrFloodedCells));//.arg(K2DQOutBoun));
    // some screen error reporting

}
