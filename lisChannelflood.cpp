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
                    double whlevel = (ChannelWH->Drcr-chdepth)*fc + _h->Drcr*(1-fc);
                    double sedlevel = 0;
                    if (SwitchErosion) {
                        sedlevel = dH/ChannelWH->Drcr*ChannelSSSed->Drcr*fc + SSFlood->Drcr*(1-fc);
                    }
                    // equilibrium water level = weighed values of channel surplus level + _h
                    // can be negative if channelwh is below channel depth and low _h level
                    double cwa = ChannelWidth->Drcr/ChannelAdj->Drcr;

                    bool dosimpel = false;//obsolete (SwitchFlood1D2DCoupling == 1);

                    if (!dosimpel)
                    {
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
                                     double sed = fracC*ChannelSSSed->Drcr;
                                     ChannelSSSed->Drcr -= sed;
                                     SSFlood->Drcr += sed;
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
                                    double sed = fracA*SSFlood->Drcr;
                                    ChannelSSSed->Drcr += sed;
                                    SSFlood->Drcr -= sed;

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
                        if(whlevel > HMIN)
                        {
                          //  qDebug() << r << c << "simpel";
                            ChannelWH->Drcr = whlevel + chdepth;
                            _h->Drcr = whlevel;
/*
                            // new equilibrium levels
                            if(SwitchErosion)
                            {
                                ChannelSSSed->Drcr = sedlevel+ ChannelSSConc->Drcr * chdepth*ChannelFlowWidth->Drcr*ChannelDX->Drcr;
                                SSFlood->Drcr = sedlevel;

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
                            */
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
//                        SWOFSedimentLayerDepth(r,c,_h->Drcr, UVflood->Drcr);//V->Drcr);
//                        RiverSedimentLayerDepth(rr,cr);

                        SWOFSedimentSetConcentration(rr,cr, _h);
                        RiverSedimentMaxC(rr, cr);
                        // all concentrations, possible ChannelDep when surplus

                    }

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

    if (SwitchKinematic2D == K2D_METHOD_DYN)
        return;

    ChannelOverflow(hmx, UVflood);
    // mix overflow water and flood water in channel cells
    // use hmx which is the 2Ddyn water

    double dtflood = 0;

    startFlood = false;
    FOR_ROW_COL_MV {
        if (hmx->Drc > 0) {
            startFlood = true;
            break;
        }
    }

    dtflood = fullSWOF2Do2light(hmx, Uflood, Vflood, DEM, true);
        //  threaded flooding

    Boundary2Ddyn(hmx, Uflood, Vflood);
    // boundary flow

    //new flood domain
    nrFloodedCells = 0;
    // used in infil and addRainfall
    FOR_ROW_COL_MV {
        if (hmx->Drc > 0)// && FloodZonePotential->Drc == 1)
        {
            FloodDomain->Drc = 1;
            nrFloodedCells += 1.0;
        }
        else
            FloodDomain->Drc = 0;
    }

    FOR_ROW_COL_MV
    {
        UVflood->Drc = sqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
        Qflood->Drc = UVflood->Drc * hmx->Drc * ChannelAdj->Drc;
        //only used for report
        V->Drc = hmx->Drc > 0 ? UVflood->Drc : V->Drc;

        // addvolume infiltrated during flood process with FSurplus
        //InfilVolFlood->Drc += Iflood->Drc;
        FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;

        // for output on screen
        hmxWH->Drc = FloodDomain->Drc  == 0 ? WH->Drc : hmx->Drc;   //hmxWH is all water
        hmxflood->Drc = hmxWH->Drc < minReportFloodHeight ? 0.0 : hmxWH->Drc;
    }
    if(SwitchErosion)
    {
        //calculate concentration and new sediment discharge
        //WHrunoff and Qn are adapted in case of 2D routing
        if(!SwitchUseGrainSizeDistribution)
        {
            FOR_ROW_COL_MV
            {
                double sed = SSFlood->Drc + BLFlood->Drc;
                Conc->Drc =  MaxConcentration(FloodWaterVol->Drc, &sed, &DepFlood->Drc);
                Qsn->Drc = Conc->Drc*Qn->Drc;
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
    FloodMaxandTiming(hmxWH, UVflood, minReportFloodHeight);

    double area = nrFloodedCells*_dx*_dx;
    if (area > 0)
    debug(QString("Flooding (dt %1 sec, n %2): area %3 m2, %4 cells %5").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1).arg(nrFloodedCells).arg(K2DQOutBoun));
    // some screen error reporting

}
