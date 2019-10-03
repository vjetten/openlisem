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


void TWorld::distributeChannelSed(int r, int c, double dh, bool fromchannel)
{
    // only susp matter, bedload cannot flow out of channel

    if (fromchannel) {
        //  qDebug() << r << c << "sed from chan";
        double vol = dh*ChannelDX->Drc*ChannelWidth->Drc;
        double sed = ChannelSSConc->Drc * vol;

        SSFlood->Drc += sed;
        ChannelSSSed->Drc -= sed;
        ChannelSed->Drc = ChannelBLSed->Drc + ChannelSSSed->Drc;

        if(SwitchUseGrainSizeDistribution)
        {
            FOR_GRAIN_CLASSES
            {
                SS_D.Drcd += RSSC_D.Drcd * vol;
                RSS_D.Drcd -= RSSC_D.Drcd * vol;
            }
            //CALC TOTALS HERE
        }
    } else {
        //   qDebug() << r << c << "sed to chan";
        double vol = dh*DX->Drc*ChannelAdj->Drc;
        double sed = SSCFlood->Drc * vol;
        SSFlood->Drc -= sed;
        ChannelSSSed->Drc += sed;
        ChannelSed->Drc = ChannelBLSed->Drc + ChannelSSSed->Drc;
        if(SwitchUseGrainSizeDistribution)
        {
            FOR_GRAIN_CLASSES
            {
                SS_D.Drcd -= SSC_D.Drcd * vol;
                RSS_D.Drcd += SSC_D.Drcd * vol;
            }
            //CALC TOTALS HERE
        }
    }
    //   ChannelSSConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSSSed->Drc);
    RiverSedimentMaxC(r,c);
    //   SSCFlood->Drc = MaxConcentration(ChannelWaterVol->Drc, SSFlood->Drc);
    SWOFSedimentMaxC(r,c);
}
//---------------------------------------------------------------------------
//! Get flood level in channel from 1D kin wave channel
//! Instantaneous mixing of flood water and channel water in channel cells
//! note: ChannelDepth lets you also control which channels flood:
//! those that are 0 react as usual (infinite capacity)
void TWorld::ChannelOverflow(cTMap *_h, cTMap *V)
{

    for (int  r = 0; r < _nrRows; r++)
    {
        for (int  c = 0; c < _nrCols; c++)
        {

            if(ChannelMaskExtended->data[r][c] == 1)// && !pcr::isMV(LDDChannel->data[r][c]))
            {
                int rr = (int)ChannelSourceYExtended->Drc;
                int cr = (int)ChannelSourceXExtended->Drc;

                double WHbef = _h->Drcr;
                double CWHbef = ChannelWH->Drcr;

                SWOFSedimentSetConcentration(rr,cr, _h);
                RiverSedimentMaxC(rr, cr);

                if (ChannelDepth->Drcr > 0 && ChannelMaxQ->Drcr <= 0)// && LDD->Drcr != 5)// && FloodZonePotential->Drc > 0)
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

                    double fracA = std::min(1.0, _dt*V->Drcr/(0.5*_dx));
                    // fraction from _h to channel based on average flood velocity
                    double fracC = std::min(1.0, _dt*(std::pow(dH, 2/3)*sqrt(std::max(Grad->Drcr,MIN_SLOPE))/N->Drcr)/(0.5*_dx));
                    // fraction from channel to surrounding based on overflow height and manning
                    double fc = std::min(0.95,ChannelWidth->Drcr/_dx);
                    // fraction of the channel in the gridcell, 1-fc = (dx-chw)/dx = chanadj/dx
                    double whlevel = (ChannelWH->Drcr - chdepth)*fc + _h->Drcr*(1-fc);
                    // equilibrium water level = weighed values of channel surplus level + _h
                    // can be negative if channelwh is below channel depth and low _h level
                    double cwa = ChannelWidth->Drcr/ChannelAdj->Drcr;

                    bool dosimpel = false;//obsolete (SwitchFlood1D2DCoupling == 1);

                    if (!dosimpel)
                    {
                        if (dH > _h->Drcr)   // flow from channel
                        {
                            double dwh = fracC * dH; // amount flowing from channel
                            if (_h->Drcr + dwh*cwa > dH-dwh) {   // if flow causes situation to reverse (channel dips below _h)
                                dosimpel = true;
                            } else {
                                _h->Drcr += dwh*cwa;
                                ChannelWH->Drcr -= dwh;

                                if(SwitchErosion) {
//                                   ChannelSSDepth->Drcr = ChannelWH->Drcr - ChannelBLDepth->Drcr;
//                                   SSDepthFlood->Drcr = _h->Drcr - BLDepthFlood->Drcr;
//                                   double sed = dwh*ChannelSSConc->Drcr;
//                                   ChannelSSSed->Drcr -= sed;
//                                   ChannelSed->Drc = ChannelBLSed->Drc + ChannelSSSed->Drc;
//                                   SSFlood->Drcr += sed;
//                                   double maxsed = MAXCONC * SSFlood->Drcr*DX->Drcr*ChannelAdj->Drcr;
//                                   if (SSFlood->Drcr > maxsed) {
//                                       BLDepFloodTot->Drcr +=  SSFlood->Drcr -maxsed ;
//                                   }
//                                   //distributeChannelSed(rr,cr,dwh, true);
                                }

                            }
                        }
                        else   // flow to channel
                        {
                            double dwh = fracA * std::max(0.0, _h->Drcr); // amount flowing to channel
                            if (dH + dwh/cwa > _h->Drcr-dwh) {   // if too much flow
                                dosimpel = true;
                            } else {
                                _h->Drcr -= dwh;
                                ChannelWH->Drcr += (dwh/cwa);
                                if(SwitchErosion) {
//                                    ChannelSSDepth->Drc = ChannelWH->Drcr - ChannelBLDepth->Drcr;
//                                    SSDepthFlood->Drcr = _h->Drcr - BLDepthFlood->Drcr;
//                                    double sed = dwh*SSCFlood->Drcr;
//                                    ChannelSSSed->Drcr += sed;
//                                    ChannelSed->Drc = ChannelBLSed->Drc + ChannelSSSed->Drc;
//                                    SSFlood->Drcr -= sed;
//                                    double maxsed = MAXCONC * ChannelSSDepth->Drcr*DX->Drcr*ChannelWidth->Drcr;
//                                    if (ChannelSSSed->Drcr > maxsed) {
//                                        ChannelDep->Drcr +=  ChannelSSSed->Drcr - maxsed ;
//                                    }

//                                   // distributeChannelSed(rr, cr, dwh,  false);
                                }
                            }
                        }
                    }

                    if (dosimpel)
                    {

                        if(whlevel > 0) // instantaneous waterlevel exquilibrium acccross channel and adjacent
                        {
                          //  qDebug() << r << c << "simpel";
                            ChannelWH->Drcr = (whlevel + chdepth);
                            _h->Drcr = whlevel;


                            // new equilibrium levels
                            if(SwitchErosion)
                            {
//                                double totsed = (ChannelSSConc->Drcr*whlevel + SSFlood->Drcr);
//                                double totvol = whlevel*DX->Drcr*_dx;
//                                double concavg = MaxConcentration(totvol, totsed);
//                                // mix SS sed top layer water above channel

//                                ChannelSSDepth->Drc = ChannelWH->Drcr - ChannelBLDepth->Drcr;
//                                SSDepthFlood->Drcr = _h->Drcr - BLDepthFlood->Drcr;

//                                ChannelSSSed->Drcr = concavg*ChannelWidth->Drcr*whlevel*ChannelDX->Drcr +
//                                                     ChannelSSConc->Drcr*chdepth*ChannelDX->Drcr*ChannelWidth->Drcr;
//                                ChannelSed->Drcr = ChannelBLSed->Drcr + ChannelSSSed->Drcr;
//                                SSFlood->Drcr = concavg*SSDepthFlood->Drcr*ChannelAdj->Drcr*DX->Drcr;

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
                    // recalc channel water vol else big MB error
                    if(SwitchErosion)
                    {
                        double WH = _h->Drcr;
                        double CWH = ChannelWH->Drcr;

                        if (CWH < CWHbef){
                //            qDebug() << "from CH"  << (CWHbef-CWH) << WH-WHbef;
                            double dsedChannel = ChannelSSConc->Drcr*(CWHbef-CWH)*ChannelWidth->Drcr*ChannelDX->Drcr;
                            SSFlood->Drcr +=  dsedChannel;
                            ChannelSSSed->Drcr -= dsedChannel;
                        }
                        if (WH < WHbef){
                     //       qDebug() << "to CH"  << CWH-CWHbef << WHbef-WH;
                            double dsedFlood = SSCFlood->Drcr*(WHbef-WH)*ChannelAdj->Drcr*DX->Drcr;
                            SSFlood->Drcr -=  dsedFlood;
                            ChannelSSSed->Drcr += dsedFlood;
                        }
                        ChannelSSDepth->Drc = ChannelWH->Drcr - ChannelBLDepth->Drcr;
                        SSDepthFlood->Drcr = _h->Drcr - BLDepthFlood->Drcr;
                        ChannelSed->Drc = ChannelBLSed->Drc + ChannelSSSed->Drc;


                        SWOFSedimentSetConcentration(rr,cr, _h);
                        RiverSedimentMaxC(rr, cr);
                    }


                   // FloodWaterVol->Drc = _h->Drc*ChannelAdj->Drc*DX->Drc;
                }
            }
        }
    }
}




//---------------------------------------------------------------------------
// correct mass balance
double TWorld::getMass(cTMap *M)
{
    double sum2 = 0;
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0)
            sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
    }
    return sum2;
}
//---------------------------------------------------------------------------
// correct mass balance
void TWorld::correctMassBalance(double sum1, cTMap *M)
{
    double sum2 = 0;
    double n = 0;
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0)
        {
            sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
            if(M->Drc > 0)
                n += 1;
        }
    }
    // total and cells active for M

    //double dh = (n > 0 ? (sum1 - sum2)/n : 0);
    double dhtot = sum2 > 0 ? (sum1 - sum2)/sum2 : 0;
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0)
        {
            M->Drc = M->Drc*(1.0 + dhtot);            // <- distribution weighted to h
            //M->Drc += dh/(DX->Drc*ChannelAdj->Drc); // <- equal distribution error
            M->Drc = std::max(M->Drc , 0.0);
        }
    }
}
//---------------------------------------------------------------------------
void TWorld::FloodBoundary()
{
    FOR_ROW_COL_MV
    {
        // NOTE : DomainEdge is a copy of  && FlowBoundary if needed
        if (FlowBoundary->Drc > 0 && hmx->Drc > 0.01)
        {
            // qDebug() << Qflood->Drc*_dt/(DX->Drc*ChannelAdj->Drc);
            double qf = Qflood->Drc;

            if( qf*_dt > hmx->Drc*DX->Drc*ChannelAdj->Drc)
            {
                qf = hmx->Drc*DX->Drc*ChannelAdj->Drc/_dt;
            }
            double hmx_old = hmx->Drc;
            hmx->Drc = std::max(0.0, hmx->Drc - qf*_dt/(DX->Drc*ChannelAdj->Drc));
            floodBoundaryTot += (hmx_old - hmx->Drc)*(DX->Drc*ChannelAdj->Drc);

            if (SwitchErosion)
            {
                double frac = (hmx_old - hmx->Drc)/hmx->Drc;
                floodBoundarySedTot += (BLFlood->Drc + SSFlood->Drc)*frac;
                BLFlood->Drc *= (1-frac);
                SSFlood->Drc *= (1-frac);
            }
        }
    }
}
//---------------------------------------------------------------------------
void TWorld::FloodMaxandTiming(cTMap *_h, cTMap *_UV, double threshold)
{
    // floodwater volume and max flood map
    FOR_CELL_IN_FLOODAREA
    {
        if (_h->Drc > threshold) {
            floodTime->Drc += _dt/60;
            floodHmxMax->Drc = std::max(floodHmxMax->Drc, _h->Drc);
//            if (floodHmxMax->Drc > threshold)
//                floodHmxMax->Drc = 0;
        // for output
            floodVMax->Drc = std::max(floodVMax->Drc, _UV->Drc);
            floodVHMax->Drc = std::max(floodVMax->Drc, _UV->Drc*_h->Drc);
        // max velocity
        }
    }}
    floodVolTotMax = 0;
    floodAreaMax = 0;
    double area = _dx*_dx;
    FOR_CELL_IN_FLOODAREA
    if (floodHmxMax->Drc > threshold)
    {
        floodVolTotMax += floodHmxMax->Drc*area;
        floodAreaMax += area;
    }}

    FOR_CELL_IN_FLOODAREA
    if (_h->Drc > threshold && floodTimeStart->Drc == 0)
    {
        floodTimeStart->Drc = (time - RainstartTime)/60.0;
      // time since first pixel received rainfall
    }}

}
//---------------------------------------------------------------------------
// NOTE DEM has barriers included, done in shade map calculation !!!!
void TWorld::ChannelFlood(void)
{
    if (!SwitchIncludeChannel)
        return;
    if (!SwitchChannelFlood)
        return;
    if (SwitchKinematic2D == K2D_METHOD_DYN)
        return;
//    if (SwitchKinematic2D == K2D_METHOD_DIFF)
//        return;
    FloodBoundary();
    // boundary flow

    ChannelOverflow(hmx, UVflood);
    // mix overflow water and flood water in channel cells

    FOR_ROW_COL_MV
    {
       FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;
    }
    double sumh_t = mapTotal(*FloodWaterVol) +mapTotal(*ChannelWaterVol);

    double dtflood = 0;

    startFlood = false;
    FOR_ROW_COL_MV {
        if (hmx->Drc > 0)
        {
            startFlood = true;
            break;
        }
    }

    if(startFlood)
        dtflood = fullSWOF2Do2light(hmx, Uflood, Vflood, DEM, true);
        //  threaded flooding
      //  dtflood = fullSWOF2RO(hmx, Uflood, Vflood, DEM, true);
       // non threaded flooding
    /*
    if (SwitchFloodSWOForder2)
    {
        dtflood = fullSWOF2Do2(hmx, Uflood, Vflood, DEM, true);
    }
    else
    if (SwitchFloodSWOForder1)
    {
        dtflood = fullSWOF2Do1(hmx, Uflood, Vflood, DEM, true);
    }
    else
    {
 // for experiments not available to user
        dtflood = fullSWOF2Do2light(hmx, Uflood, Vflood, DEM, true);
    }
*/

    //infilInWave(Iflood, hmx, _dt);

    FOR_CELL_IN_FLOODAREA
    {
        UVflood->Drc = sqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
        // U and V are vectors so can be negative, UV is positive average
        Qflood->Drc = UVflood->Drc * hmx->Drc * ChannelAdj->Drc;

        // addvolume infiltrated during flood process with FSurplus
        InfilVolFlood->Drc += Iflood->Drc;
    }}

    FOR_ROW_COL_MV
    {
       FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;
    }

    double sumh_t2 = mapTotal(*FloodWaterVol)+mapTotal(*ChannelWaterVol);

    //correctMassBalance(sumh_t, FloodWaterVol, 1e-12);
    // correct mass balance, VJ 150823: not nnecessary hhere if flow is false

    double ncells = 0;
    double floodtot = 0;
    FOR_ROW_COL_MV
    {
       if(FloodWaterVol->Drc > 0.0001)
        {
            ncells += 1.0;
            floodtot += FloodWaterVol->Drc;
        }
    }
    double diff = (sumh_t-sumh_t2);

    FOR_ROW_COL_MV
    {
       FloodWaterVol->Drc = floodtot > 0? FloodWaterVol->Drc  + diff * FloodWaterVol->Drc/floodtot : 0.0;
       hmx->Drc = FloodWaterVol->Drc / (ChannelAdj->Drc *DX->Drc);
    }

    //new flood domain
    nrFloodedCells = 0;

//    sumh_t2 = mapTotal(*FloodWaterVol)+mapTotal(*ChannelWaterVol);
//    floodBoundaryTot +=(sumh_t-sumh_t2);
    // cheat!

    // used in infil and addRainfall
    FOR_CELL_IN_FLOODAREA
        if (hmx->Drc > 0)// && FloodZonePotential->Drc == 1)
        {
            FloodDomain->Drc = 1;
            nrFloodedCells += 1.0;
        }
        else
            FloodDomain->Drc = 0;
    }


    // add RO waterheight and hmx for output, and calc flood for output
    FOR_ROW_COL_MV {
        hmxWH->Drc = hmx->Drc + WH->Drc;
        hmxflood->Drc = hmxWH->Drc < minReportFloodHeight ? 0.0 : hmxWH->Drc;
        FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;
    }

    fill(*tma, 0);
    FOR_ROW_COL_MV {
       tma->Drc = UVflood->Drc > 0 ? UVflood->Drc : V->Drc;
    }

    FloodMaxandTiming(hmxWH, tma, minReportFloodHeight);
    // flood max, start and duration

    FOR_CELL_IN_FLOODAREA
       FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;
    }

    //double avgh = (cells > 0 ? (sumh_t)/cells : 0);
    double area = nrFloodedCells*_dx*_dx;
    if (area > 0)
    debug(QString("Flooding (dt %1 sec, n %2): area %3 m2, %4 cells").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1).arg(nrFloodedCells));
    // some screen error reporting

}
