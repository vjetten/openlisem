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
//#include "model.h"
#include "operation.h"
#include "global.h"

void TWorld::ChannelOverflowNew(cTMap *_h, cTMap *V, bool doOF)
{
    cTMap *_SS;
    cTMap *_SSC;
    if (doOF) {
        _SS = Sed;
        _SSC = Conc;
    } else {
        _SS = SSFlood;
        _SSC = SSCFlood;
    }

  //  DistributeOverExtendedChannel(ChannelWaterVol,ChannelVolExtended);

    FOR_ROW_COL_MV_CH {
        if(ChannelMaskExtended->data[r][c] == 1)
        {
            int rr = r;//(int)ChannelSourceYExtended->Drc;
            int cr = c;//(int)ChannelSourceXExtended->Drc;

          //  ChannelWHExtended->Drc = ChannelWH->Drcr;
            if (doOF)
                Conc->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*_h->Drc, &Sed->Drc, &DEP->Drc);
            else
                SWOFSedimentSetConcentration(r,c, _h);

            RiverSedimentMaxC(r, c);

            if (ChannelDepth->Drcr > 0 && ChannelMaxQ->Drcr <= 0)
            {
                double chdepth = ChannelDepth->Drcr;
                double dH = std::max(0.0, (ChannelWH->Drcr-chdepth));

                if (dH <= HMIN && _h->Drc <= HMIN)
                    continue;
                // no flow activity then continue

                if (dH == _h->Drc)
                    continue;
                // no diff in water level, no flow, continue

                // VELOCITIES
                double VtoChan = V->Drcr;//std::pow(_h->Drcr, 2.0/3.0)*sqrt(ChannelPAngle->Drc)/N->Drcr; //F_Angle
               // double VfromChan = std::pow(dH, 2.0/3.0)*sqrt(ChannelPAngle->Drc)/N->Drcr;
              //  if (F_AddGravity == 1) {
                   double VfromChan = sqrt(2*9.81*dH);
             //   }
                double fracA = std::min(1.0, _dt*VtoChan/(0.5*_dx));
                // fraction from _h to channel based on average flood velocity
                double fracC = std::min(1.0, _dt*VfromChan/(0.5*_dx));
                // fraction from channel to surrounding based on overflow height and manning

                bool dosimpel = false;//obsolete (SwitchFlood1D2DCoupling == 1);

                if (!dosimpel)
                {
                    double cwa = ChannelAdj->Drc > 0 ? ChannelWidthMax->Drcr/ChannelAdj->Drc : 0;

                    if (dH > _h->Drc)   // flow from channel
                    {
                        double dwh = fracC * dH;
                        // amount flowing from channel
                        if (_h->Drc + dwh*cwa > dH-dwh) {
                            // if flow causes situation to reverse (channel dips below _h)
                            dosimpel = true;
                        } else {
                            _h->Drc += dwh*cwa;
                            //ChannelWHExtended->Drc -= dwh;
                            ChannelWH->Drcr -= dwh;
                            // assumes dH is above channel and therefore rectangular

                            if(SwitchErosion) {
                                double sed = fracC*ChannelSSSed->Drcr;
                                ChannelSSSed->Drcr -= sed;
                                _SS->Drc += sed;
                                if(SwitchUseGrainSizeDistribution)
                                {
                                    FOR_GRAIN_CLASSES
                                    {
                                        //   SS_D.Drcd += RSSC_D.Drcd * vol;
                                        //  RSS_D.Drcd -= RSSC_D.Drcd * vol;
                                    }
                                    //CALC TOTALS HERE
                                }
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
                            ChannelWaterVol->Drcr += dwh*ChannelAdj->Drcr*DX->Drcr;
                            fromChannelVoltoWH(rr,cr);
                            //ChannelWHExtended->Drc = ChannelWH->Drcr;

                            if(SwitchErosion) {
                                double sed = fracA*_SS->Drc;
                                ChannelSSSed->Drcr += sed;
                                _SS->Drc -= sed;

                                if(SwitchUseGrainSizeDistribution)
                                {
                                    FOR_GRAIN_CLASSES
                                    {
                                        //     SS_D.Drcd += RSSC_D.Drcd * vol;
                                        //     RSS_D.Drcd -= RSSC_D.Drcd * vol;
                                    }
                                    //CALC TOTALS HERE
                                }
                            }
                        }
                    }
                }

                // instantaneous waterlevel exquilibrium acccross channel and adjacent
                if (dosimpel)
                {
                    double fc = std::min(0.95,ChannelWidthMax->Drcr/_dx);
                    // fraction of the channel in the gridcell, 1-fc = (dx-chw)/dx = chanadj/dx
                    double whlevel = (ChannelWH->Drcr-chdepth)*fc + _h->Drc*(1-fc);
                    double voltot = ChannelWaterVol->Drc + DX->Drcr*_h->Drc*ChannelAdj->Drc;

                    // equilibrium water level = weighed values of channel surplus level + _h
                    // can be negative if channelwh is below channel depth and low _h level
                    if(whlevel > HMIN)
                    {
                        //ChannelWHExtended->Drc = whlevel + chdepth;
                        _h->Drcr = voltot*(1-fc)/(DX->Drcr*ChannelAdj->Drcr);
                                //whlevel;

                        ChannelWaterVol->Drcr = voltot*fc;
                        fromChannelVoltoWH(rr,cr);
                        //ChannelWHExtended->Drc = ChannelWH->Drcr;
                        // new equilibrium levels
                        if(SwitchErosion)
                        {

                            RiverSedimentLayerDepth(rr,cr);
                            //SWOFSedimentLayerDepth(rr, cr, _h->Drcr, V->Drcr);

                            double _sed = ChannelSSSed->Drcr + _SS->Drc;
                            double volch = ChannelSSDepth->Drcr*ChannelWidthExtended->Drc*ChannelDX->Drcr;
                            double volof = _h->Drc*ChannelAdj->Drc*DX->Drc;
                            double _concavg = _sed/(volch+volof);

                            ChannelSSSed->Drcr = _concavg * volch;
                            ChannelSed->Drcr = ChannelBLSed->Drcr + ChannelSSSed->Drcr;

                            _SS->Drc = _concavg * volof;
                            _SSC->Drc = _concavg;

                            if(SwitchUseGrainSizeDistribution)
                            {
                                FOR_GRAIN_CLASSES
                                {
                                    // SS_D.Drcd += RSSC_D.Drcd * vol;
                                    // RSS_D.Drcd -= RSSC_D.Drcd * vol;
                                }
                                //CALC TOTALS HERE
                            }
                        }

                    }
                    else
                    {
                        //DO NOTHING
                        // this happens if there is very little flood water (< 5cm) and the channelWH is below the channeldepth
                        // we assume that there is no more flow towards the channel.
                    }
                }

                //ChannelVolExtended->Drc = ChannelWHExtended->Drc * ChannelDX->Drcr * ChannelWidthExtended->Drc;
                //                    // do not recalc floodvol, MB errors

                // recalc channel water vol else big MB error
                if(SwitchErosion)
                {
                    if (doOF)
                        Conc->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*_h->Drc, &_SS->Drc, &DEP->Drc);
                    else {
                        SWOFSedimentLayerDepth(r,c,_h->Drc, V->Drc);
                        SWOFSedimentSetConcentration(r,c, _h);
                    }

                    RiverSedimentLayerDepth(rr,cr);
                    RiverSedimentMaxC(rr, cr);
                    // all concentrations, possible ChannelDep when surplus

                }

            }
        }
    }

//    fill(*ChannelWaterVol, 0);
//    FOR_ROW_COL_MV {
//        if(ChannelMaskExtended->data[r][c] == 1)
//        {
//            int rr = (int)ChannelSourceYExtended->Drc;
//            int cr = (int)ChannelSourceXExtended->Drc;
//            ChannelWaterVol->Drcr += ChannelVolExtended->Drc;
//        }
//    }
    CalcVelDischChannelNT();
}
//---------------------------------------------------------------------------
//! Get flood level in channel from 1D kin wave channel
//! Instantaneous mixing of flood water and channel water in channel cells
//! note: ChannelDepth lets you also control which channels flood:
//! those that are 0 react as usual (infinite capacity)
void TWorld::ChannelOverflow(cTMap *_h, cTMap *V, bool doOF)
{
    cTMap *_SS;
    cTMap *_SSC;
    if (doOF) {
        _SS = Sed;
        _SSC = Conc;
    } else {
        _SS = SSFlood;
        _SSC = SSCFlood;
    }

    FOR_ROW_COL_MV_CH {
        if(ChannelMaskExtended->data[r][c] == 1)// && !pcr::isMV(LDDChannel->data[r][c]))
        {
            int rr = r;//(int)ChannelSourceYExtended->Drc;
            int cr = c;//(int)ChannelSourceXExtended->Drc;

            if (doOF)
                Conc->Drcr = MaxConcentration(ChannelAdj->Drc*DX->Drc*_h->Drcr, &Sed->Drc, &DEP->Drc);
            else
                SWOFSedimentSetConcentration(rr,cr, _h);

            RiverSedimentMaxC(rr, cr);

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
                double VtoChan = V->Drc;//std::pow(_h->Drcr, 2.0/3.0)*sqrt(ChannelPAngle->Drc)/N->Drcr; //F_Angle
                double VfromChan = sqrt(2*9.81*dH);
                double fracA = std::min(1.0, _dt*VtoChan/(0.5*_dx));
                // fraction from _h to channel based on average flood velocity
                double fracC = std::min(1.0, _dt*VfromChan/(0.5*_dx));
                // fraction from channel to surrounding based on overflow height and manning

                double cwa = ChannelWidth->Drcr/ChannelAdj->Drcr;

                bool dosimpel = false;//obsolete (SwitchFlood1D2DCoupling == 1);

                if (dH > _h->Drcr)   // flow from channel
                {
                    double dwh = fracC * dH;
                    // amount flowing from channel
                    if (_h->Drcr + dwh*cwa > dH-dwh) {
                        // if flow causes situation to reverse (channel dips below _h)
                        dosimpel = true;
                    } else {

                        _h->Drcr += dwh*cwa;
                        ChannelWH->Drcr -= dwh;

                        if(SwitchErosion) {
                            //  double sed = dwh*ChannelSSConc->Drcr;
                            double sed = fracC*ChannelSSSed->Drcr;
                            ChannelSSSed->Drcr -= sed;
                            _SS->Drcr += sed;
                            if(SwitchUseGrainSizeDistribution)
                            {
                                FOR_GRAIN_CLASSES
                                {
                                    //   SS_D.Drcd += RSSC_D.Drcd * vol;
                                    //  RSS_D.Drcd -= RSSC_D.Drcd * vol;
                                }
                                //CALC TOTALS HERE
                            }
                        }

                    }
                }
                else   // flow to channel
                {
                    double dwh = fracA * std::max(0.0, _h->Drcr);
                    // amount flowing to channel
                    if (dH + dwh/cwa > _h->Drcr-dwh) {
                        // if too much flow
                        dosimpel = true;
                    } else {

                        _h->Drcr -= dwh;
                        ChannelWH->Drcr += (dwh/cwa);
                        if(SwitchErosion) {
                            double sed = fracA*_SS->Drcr;
                            //double sed = dwh*_SSC->Drcr;
                            ChannelSSSed->Drcr += sed;
                            _SS->Drcr -= sed;

                            if(SwitchUseGrainSizeDistribution)
                            {
                                FOR_GRAIN_CLASSES
                                {
                                    //     SS_D.Drcd += RSSC_D.Drcd * vol;
                                    //     RSS_D.Drcd -= RSSC_D.Drcd * vol;
                                }
                                //CALC TOTALS HERE
                            }
                        }
                    }
                }

                // instantaneous waterlevel exquilibrium acccross channel and adjacent
                if (dosimpel)
                {
                    double fc = std::min(0.95,ChannelWidth->Drcr/_dx);
                    // fraction of the channel in the gridcell, 1-fc = (dx-chw)/dx = chanadj/dx
                    double whlevel = (ChannelWH->Drcr-chdepth)*fc + _h->Drcr*(1-fc);
                    // equilibrium water level = weighed values of channel surplus level + _h
                    // can be negative if channelwh is below channel depth and low _h level
                    if(whlevel > HMIN)
                    {
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

                            if(SwitchUseGrainSizeDistribution)
                            {
                                FOR_GRAIN_CLASSES
                                {
                                    // SS_D.Drcd += RSSC_D.Drcd * vol;
                                    // RSS_D.Drcd -= RSSC_D.Drcd * vol;
                                }
                                //CALC TOTALS HERE
                            }
                        }

                    }
                    else
                    {
                        //DO NOTHING
                        // this happens if there is very little flood water (< 5cm) and the channelWH is below the channeldepth
                        // we assume that there is no more flow towards the channel.
                    }
                }

                ChannelWaterVol->Drcr = ChannelWH->Drcr * ChannelDX->Drcr * ChannelWidth->Drcr;
                // do not recalc floodvol, MB errors

                // recalc channel water vol else big MB error
                if(SwitchErosion)
                {
                    if (doOF)
                        Conc->Drcr = MaxConcentration(ChannelAdj->Drc*DX->Drc*_h->Drcr, &_SS->Drc, &DEP->Drc);
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
void TWorld::ToFlood()//int thread)
{
    FOR_ROW_COL_MV {
        if(hmx->Drc > 0.0 && WHrunoff->Drc > 0.0)
        {
            double frac = 1-exp(-runoff_partitioning*hmx->Drc/WHrunoff->Drc);
            frac = std::max(std::min(frac, 1.0),0.0);
            double dwh = frac * WHrunoff->Drc;
            // if flowwidth != channeladj then deal with this here
            hmx->Drc += dwh;
            WH->Drc -= dwh;
            WHrunoff->Drc -= dwh;
            WHGrass->Drc -= dwh;
            WHroad->Drc -= dwh;

            if(SwitchErosion)
            {
                double dsed = frac*Sed->Drc;
                SSFlood->Drc += dsed;
                Sed->Drc -= dsed;

                if(SwitchUseGrainSizeDistribution)
                {
                    FOR_GRAIN_CLASSES
                    {
                        SS_D.Drcd +=  Sed_D.Drcd * frac;
                        Sed_D.Drcd = Sed_D.Drcd * (1-frac);

                    }
                }
//                SWOFSedimentSetConcentration(r,c,hmx);
//                Conc->Drc = MaxConcentration(WH->Drc*ChannelAdj->Drc*DX->Drc, &Sed->Drc, &DEP->Drc);
            }
         }
    }
}
//---------------------------------------------------------------------------
void TWorld::FloodMaxandTiming(cTMap *_h, cTMap *_UV, double threshold)
{
    // floodwater volume and max flood map
    FOR_ROW_COL_MV {
        if (_h->Drc > threshold) {
            floodTime->Drc += _dt/60;
            floodHmxMax->Drc = std::max(floodHmxMax->Drc, _h->Drc);
        // for output
            floodVMax->Drc = std::max(floodVMax->Drc, _UV->Drc);
            floodVHMax->Drc = std::max(floodVMax->Drc, _UV->Drc*_h->Drc);
        // max velocity
        }
    }

    floodVolTotMax = 0;
    floodAreaMax = 0;
    double area = _dx*_dx;
    FOR_ROW_COL_MV {
        if (floodHmxMax->Drc > threshold)
        {
            floodVolTotMax += floodHmxMax->Drc*area;
            floodAreaMax += area;
        }
    }

    FOR_ROW_COL_MV {
        if (_h->Drc > threshold && floodTimeStart->Drc == 0)
        {
            floodTimeStart->Drc = (time - RainstartTime)/60.0;
            // time since first pixel received rainfall
        }
    }

}
//---------------------------------------------------------------------------
// NOTE DEM has barriers included, done in shade map calculation !!!!
void TWorld::ChannelFlood(void)
{
    if (!SwitchIncludeChannel)
        return;

    if (SwitchKinematic2D != K2D_METHOD_KINDYN)
        return;

    ChannelOverflowNew(hmx, V, false);
    // determine overflow water => hmx

 //   ToFlood();
    // mix WHrunoff with hmx

    double dtflood = 0;

    startFlood = false;
    FOR_ROW_COL_MV {
        if (hmx->Drc > 0) {
            startFlood = true;
            break;
        }
    }

    dtflood = fullSWOF2Do2light(hmx, Uflood, Vflood, DEM, true);
    // 2D dyn flow of hmx water

    //new flood domain
    nrFloodedCells = 0;
    // used in infil and addRainfall
    FOR_ROW_COL_MV {
        if (hmx->Drc > 0)
        {
            FloodDomain->Drc = 1;
            nrFloodedCells += 1.0;
        }
        else
            FloodDomain->Drc = 0;
    }


    ToFlood();

    FOR_ROW_COL_MV {
        Qflood->Drc = 0;
        if (FloodDomain->Drc > 0) {
            V->Drc = qSqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
            Qflood->Drc = V->Drc * hmx->Drc * ChannelAdj->Drc;
        }
    }
    Boundary2Ddyn();//hmx, Qflood, Uflood, Vflood);
    // boundary flow

    FOR_ROW_COL_MV {       
        if (FloodDomain->Drc > 0) {
     //       V->Drc = qSqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
            Qflood->Drc = V->Drc * hmx->Drc * ChannelAdj->Drc;

//            double R = WHrunoff->Drc*ChannelAdj->Drc/(2*WHrunoff->Drc + ChannelAdj->Drc);
//            double vv = pow(R, 2.0/3.0) * sqrt(Grad->Drc)/N->Drc;
//            double qq  = vv * WHrunoff->Drc * ChannelAdj->Drc;
//            Qn->Drc += Qflood->Drc;

        }

        hmxWH->Drc = WH->Drc + hmx->Drc;

        //InfilVolFlood->Drc += Iflood->Drc;
        // addvolume infiltrated during flood process with FSurplus

        hmxflood->Drc = std::max(0.0, WHrunoff->Drc + hmx->Drc - minReportFloodHeight);
        WHrunoffOutput->Drc = std::min(WHrunoff->Drc + hmx->Drc, minReportFloodHeight);

        WaterVolall->Drc = DX->Drc * (WHrunoff->Drc*ChannelAdj->Drc + hmx->Drc * ChannelAdj->Drc
                                      + WHstore->Drc*SoilWidthDX->Drc);
        // all water on surface

        FloodWaterVol->Drc = hmxflood->Drc*ChannelAdj->Drc*DX->Drc;
        RunoffWaterVol->Drc = WHrunoffOutput->Drc*ChannelAdj->Drc*DX->Drc;
    }

    FloodMaxandTiming(hmxWH, V, minReportFloodHeight);

    if(SwitchErosion)
    {
        //calculate concentration and new sediment discharge
        //WHrunoff and Qn are adapted in case of 2D routing
        if(!SwitchUseGrainSizeDistribution)
        {
            FOR_ROW_COL_MV {
                if (FloodDomain->Drc  > 0) {
                    double sed = SSFlood->Drc + BLFlood->Drc;
                    Conc->Drc =  MaxConcentration(FloodWaterVol->Drc, &sed, &DepFlood->Drc);
                    Qsn->Drc += Conc->Drc*Qflood->Drc;
                }
            }
        }
        else
        {
            //calculate total sediment from induvidual grain classes,
            //and calculate concentration and new sediment discharge
//            FOR_ROW_COL_MV
//            {
//                Sed->Drc = 0;
//                Conc->Drc = 0;
//            }
//            FOR_ROW_COL_MV
//            {
//                FOR_GRAIN_CLASSES
//                {
//                    Sed->Drc += Sed_D.Drcd;
//                    Conc_D.Drcd = MaxConcentration(FloodWaterVol->Drc, Sed_D.Drcd);
//                    Conc->Drc += Conc_D.Drcd;
//                }
//                Qsn->Drc = Conc->Drc*Qn->Drc;
//            }
        }
    }

    double area = nrFloodedCells*_dx*_dx;
    if (area > 0)
    debug(QString("Flooding (dt %1 sec, n %2): area %3 m2, %4 cells").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1).arg(nrFloodedCells));//.arg(K2DQOutBoun));
    // some screen error reporting

}
