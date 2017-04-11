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
#include "model.h"
#include "operation.h"
#include "global.h"


void TWorld::distributeChannelSed(int r, int c, double dh, double charea)
{
    SSFlood->Drc += ChannelConc->Drc * dh * charea;
    if(this->SwitchUse2Layer)
    {
        ChannelSSSed->Drc -= ChannelConc->Drc * dh * charea;
    }else
    {
        ChannelBLSed->Drc -= ChannelConc->Drc * dh * charea;
    }
    if(SwitchUseGrainSizeDistribution)
    {
        FOR_GRAIN_CLASSES
        {
            SS_D.Drcd +=RSSC_D.Drcd * dh * charea;
            RSS_D.Drcd -= RSSC_D.Drcd * dh * charea;
        }
    }

}
//---------------------------------------------------------------------------
//! Get flood level in channel from 1D kin wave channel
//! Instantaneous mixing of flood water and channel water in channel cells
//! note: ChannelDepth lets you also control which channels flood:
//! those that are 0 react as usual (infinite capacity)
void TWorld::ChannelOverflow()
{
    for (int  r = 0; r < _nrRows; r++)
    {
        for (int  c = 0; c < _nrCols; c++)
        {
            if(ChannelMaskExtended->data[r][c] == 1 && !pcr::isMV(LDDChannel->data[r][c]))
            {
                int rr = (int)ChannelSourceYExtended->Drc;
                int cr = (int)ChannelSourceXExtended->Drc;

                if (ChannelDepth->Drcr > 0 && ChannelMaxQ->Drcr == 0 && LDD->Drcr != 5)// && FloodZonePotential->Drc > 0)
                {

                    double levee = 0;//ChannelLevee->Drc;
                    //double chdepth = ChannelDepth->Drc + levee; // levee always assumed on both sides channel
                    double chdepth = ChannelDepth->Drcr + levee; // levee always assumed on both sides channel
                    double chwidth = ChannelWidth->Drcr;
                    double charea = chwidth*ChannelDX->Drcr;
                    double dH = std::max(0.0, (ChannelWH->Drcr-chdepth));

                    if (dH == 0 && hmx->Drc <= levee)
                        continue;

                    // no flow activity then continue
                    if (dH == hmx->Drc)
                        continue;
                    // no diff in water level, no flow, continue

                    double fracA = std::min(1.0, _dt*UVflood->Drc/(0.5*_dx));
                    // fraction from hmx to channel based on avefrage flood velocity
                    double fracC = std::min(1.0, _dt*(std::pow(dH, 2/3)*sqrt(std::max(Grad->Drc,MIN_SLOPE))/N->Drc)/(0.5*_dx));

                    double fc = std::min(0.95,chwidth/_dx);
                    // fraction of the channel in the gridcell, 1-fc = (dx-chw)/dx = chanadj/dx
                    double whlevel = (ChannelWH->Drcr - chdepth)*fc + std::max(0.0, hmx->Drc-levee)*(1-fc);
                    // equilibrium water level = weighed values of channel surplus level + hmx, levee is counted as barrier
                    // can be negative if channelwh is below channel depth and low hmx level
                    double cwa = ChannelWidth->Drcr/ChannelAdj->Drc;
                    double widechannelfraction =(ChannelWidthExtended->Drc/ChannelWidth->Drcr);

                    bool dosimpel = (SwitchFlood1D2DCoupling == 1);

                    if (SwitchFlood1D2DCoupling == 2)
                    {
                        if (dH > hmx->Drc)   // flow from channel
                        {
                            double dwh = fracC * dH; // amount flowing from channel
                            if (hmx->Drc + dwh*cwa > dH-dwh)   // if flow causes situation to reverse (channel dips below hmx)
                            {
                                dosimpel = true;
                            }
                            else
                            {
                                //qDebug() << "from" << fracC;
                                // do the flow
                                hmx->Drc += dwh*cwa;
                                ChannelWH->Drcr -= dwh;

                                //transport sediment with water
                                if(SwitchErosion)
                                {
                                    distributeChannelSed(rr,cr,dwh, charea);
                                }
                            }
                        }
                        else   // flow to channel
                        {
                            double dwh = fracA * std::max(0.0, hmx->Drc-levee); // amount flowing to channel
                            if (dH + dwh/cwa > hmx->Drc-dwh)   // if too much flow
                            {
                                dosimpel = true;
                            }
                            else
                            {
                                //qDebug() << "to" << fracA;
                                //do flow
                                hmx->Drc -= dwh;
                                ChannelWH->Drcr += (dwh/cwa);

                                //transport sediment with water
                                if(SwitchErosion)
                                {
                                    distributeChannelSed(rr, cr, dwh, charea);
                                }
                            }
                        }
                    }

                    if (dosimpel)
                    {
                        if(whlevel > 0) // instantaneous waterlevel exquilibrium acccross channel and adjacent
                        {
                            double hmxold = hmx->Drc;
                            double chwold =  widechannelfraction *ChannelWH->Drcr;

                            ChannelWH->Drcr = (whlevel + chdepth);// + (1.0-widechannelfraction) *ChannelWH->Drcr;
                            hmx->Drc = std::min(hmx->Drc, levee) + whlevel;

                            //hmx->Drc = std::max(0.0,hmx->Drc - (ChannelWH->Drcr - chwold)/cwa);
                            // cutoff hmx at levee but can be smaller
                            //ChannelWH->Drcr = ChannelWH->Drcr  - (hmx->Drc - hmxold)*cwa;


                            //transport sediment with water
                            if(SwitchErosion)
                            {
                                double dhmx = hmx->Drc - hmxold;

                                distributeChannelSed(rr,cr, dhmx,charea);
                                // if dhmx < 0 this also works because all signs are reversed
                            }
                        }
                        else
                        {
                            // this happens if there is very little flood water (< 5cm) and the channelWH is below the channeldepth
                            // we assume that there is no more flow towards the channel.
                        }
                    }
                    ChannelWaterVol->Drcr = ChannelWH->Drcr * ChannelDX->Drcr * ChannelWidth->Drcr;
                    // recalc channel water vol else big MB error
                }
            }
        }
    }

//    FOR_ROW_COL_MV_CH
//    {
//        ChannelWaterVol->Drc = ChannelWidth->Drc * ChannelWH->Drc * ChannelDX->Drc;
//        //??? what if channelside != 0:

//    }
    // recalc channel water vol else big MB error
}
//---------------------------------------------------------------------------
// correct mass balance
double TWorld::getMass(cTMap *M)
{
    double sum2 = 0;
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0)
        {
            sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
        }
    }
    return sum2;
}
//---------------------------------------------------------------------------
// correct mass balance
double TWorld::correctMassBalance(double sum1, cTMap *M, double minV)
{
    double sum2 = 0;
    double n = 0;
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0)
        {
            sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
            if(M->Drc > minV)
                n += 1;
        }
    }
    // total and cells active for M

    double dh = (n > 0 ? (sum1 - sum2)/n : 0);
    FOR_ROW_COL_MV
    {
        if(M->Drc > minV)
        {
            M->Drc += dh/(DX->Drc*ChannelAdj->Drc);
            M->Drc = std::max(M->Drc , 0.0);
        }
    }
    return dh;
}
//---------------------------------------------------------------------------
//OBSOLETE !!!
void TWorld::FloodSpuriousValues()
{
//    fill(*tm, 0.0);
//    //calc2Maps(*tma, *DEM, *hmx, ADD);
//    FOR_ROW_COL_MV
//    {
//        if (hmx->Drc > F_extremeHeight)
//        {
//            tm->Drc = getWindowAverage(*hmx, r, c, false);
//           // tm->Drc = getWindowAverage(*tma, r, c, false);
//        }
//    }

//    FOR_ROW_COL_MV
//    {
//        if ((hmx->Drc > F_extremeHeight) && (hmx->Drc > F_extremeHeight && hmx->Drc > tm->Drc + F_extremeDiff))
//        {
//            double htmp = hmx->Drc;
//            hmx->Drc = std::min( tm->Drc, std::min(hmx->Drc, Hmx->Drc));
//            qDebug() << hmx->Drc << Hmx->Drc << tm->Drc << htmp << r << c ;
//        }
//    }
}
//---------------------------------------------------------------------------
void TWorld::FloodBoundary()
{
    FOR_ROW_COL_MV
    {
        // NOTE : DomainEdge is a copy of  && FlowBoundary if needed
        if (DomainEdge->Drc > 0 && hmx->Drc > 0)
        {
            //     qDebug() << Qflood->Drc*_dt/(DX->Drc*ChannelAdj->Drc);

            if( Qflood->Drc*_dt > hmx->Drc*DX->Drc*ChannelAdj->Drc)
            {
                Qflood->Drc = hmx->Drc*DX->Drc*ChannelAdj->Drc/_dt;
            }
            double hmx_old = hmx->Drc;
            hmx->Drc = std::max(0.0, hmx->Drc - Qflood->Drc*_dt/(DX->Drc*ChannelAdj->Drc));
            floodBoundaryTot += (hmx_old - hmx->Drc)*(DX->Drc*ChannelAdj->Drc);
        }
    }
}
//---------------------------------------------------------------------------
void TWorld::FloodMaxandTiming()
{
    // floodwater volume and max flood map
    FOR_CELL_IN_FLOODAREA
    {
        floodHmxMax->Drc = std::max(floodHmxMax->Drc, hmx->Drc);
        if (hmx->Drc > 0)//minReportFloodHeight)
            floodTime->Drc += _dt/60;
        // for output
        floodVMax->Drc = std::max(floodVMax->Drc, UVflood->Drc);
        // max velocity
    }}

    floodVolTotMax = 0;
    floodAreaMax = 0;
    double area = _dx*_dx;
    FOR_CELL_IN_FLOODAREA
        if (floodHmxMax->Drc > 0)//minReportFloodHeight)
        {
            floodVolTotMax += floodHmxMax->Drc*area;
            floodAreaMax += area;
        }
    }

    FOR_CELL_IN_FLOODAREA
        if (hmx->Drc > 0 && floodTimeStart->Drc == 0)
        {
            //            FloodTimeStart->Drc = (time - RainpeakTime)/60;
            floodTimeStart->Drc = (time - RainstartTime)/60;
            // time since first pixel received rainfall
        }
    }
}
//---------------------------------------------------------------------------
 //change flood parameters while running
void TWorld::getFloodParameters(void)
{
    SwitchFloodSWOForder2 = (op.F_solution == 2);
    SwitchFloodSWOForder1 = (op.F_solution == 1);
    F_scheme = op.F_scheme;                        //MUSCL
    F_fluxLimiter = op.F_fluxLimiter;              //HLL
    F_replaceV = op.F_replaceV;
    F_maxVelocity = op.F_maxVelocity;
    F_extremeHeight = F_extremeHeight;
    F_extremeDiff = op.F_extremeDiff;
    courant_factor = op.F_courant;
    courant_factor_diffusive = op.F_courant_diffusive;
    F_MaxIter = op.F_Maxiter;
}
//---------------------------------------------------------------------------
// NOTE DEM has barriers included, done in shade map calculation !!!!
void TWorld::ChannelFlood(void)
{
    if (!SwitchIncludeChannel)
        return;
    if (!SwitchChannelFlood)
        return;

    getFloodParameters();

    FloodBoundary();
    // boundary flow

    FOR_ROW_COL_MV
    {
       FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;
    }
    double sumh_t = mapTotal(*FloodWaterVol) +mapTotal(*ChannelWaterVol);

    ChannelOverflow();
    // mix overflow water and flood water in channel cells

    FOR_ROW_COL_MV
    {
       FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;
    }

    double dtflood = 0;

    startFlood = false;
    FOR_CELL_IN_FLOODAREA
        if (hmx->Drc > 0)
        {
            startFlood = true;
            break;
        }
    }

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
        dtflood = fullSWOF2Do2light(hmx, Uflood, Vflood, DEM, true);
    }

    //FloodSpuriousValues();
    //correct extremes

    FOR_CELL_IN_FLOODAREA
    {
        UVflood->Drc = sqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
        // U and V are vectors so can be negative, UV is positive average
        Qflood->Drc = UVflood->Drc * hmx->Drc * ChannelAdj->Drc;
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

    FOR_CELL_IN_FLOODAREA
        if (hmx->Drc > 0 && FloodZonePotential->Drc == 1)
        {
            FloodDomain->Drc = 1;
            nrFloodedCells += 1.0;
        }
        else
            FloodDomain->Drc = 0;
    }


    // copy(*Hmx, *hmx);
    // copy flood level for next dt
            // for spurious values, obsolete

    calc2Maps(*hmxWH, *hmx, *WH, ADD);
    // add RO waterheight and hmx for output

    FloodMaxandTiming();
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
