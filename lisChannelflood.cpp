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

//---------------------------------------------------------------------------
//! Get flood level in channel from 1D kin wave channel
//! Instantaneous mixing of flood water and channel water in channel cells
//! note: ChannelDepth lets you also control which channels flood:
//! those that are 0 react as usual (infinite capacity)
void TWorld::ChannelOverflow()
{
    double err = 0;
    double factor = 0.5;
    FOR_ROW_COL_MV_CH
    {
        if (ChannelDepth->Drc > 0 && ChannelMaxQ->Drc == 0 && LDD->Drc != 5)// && FloodZonePotential->Drc > 0)
        {
            double levee = ChannelLevee->Drc;
            double chdepth = ChannelDepth->Drc + levee; // levee always assumed on both sides channel
            double dH = std::max(0.0, (ChannelWH->Drc-chdepth));

            if (dH == 0 && hmx->Drc <= levee)
                continue;

            //            if (dH == 0 && hmx->Drc <= MDS->Drc)
            //                continue;

            // no flow activity then continue
            if (dH == hmx->Drc)
                continue;
            // no diff in water level, no flow, continue

            double fracA = std::min(1.0, _dt*UVflood->Drc/(0.5*_dx));//ChannelAdj->Drc));
            // fraction from hmx to channel based on avefrage flood velocity
            //double fracC = std::min(1.0, _dt*sqrt(ChannelWH->Drc*9.8067)/(0.5*_dx));//ChannelWidthUpDX->Drc));
            // fraction from channel to surrounding adj area based on gravity flow perpedicular to channel
            double fracC = std::min(1.0, _dt*(std::pow(dH, 2/3)*sqrt(Grad->Drc)/N->Drc)/(0.5*_dx));//ChannelWidthUpDX->Drc));
            //            tm->Drc = fracA;
            //            tma->Drc = fracC;
            double fc = ChannelWidthUpDX->Drc/_dx; // 1-fc = (dx-chw)/dx = chanadj/dx
            // fraction of the channel in the gridcell
            double whlevel = (ChannelWH->Drc - chdepth)*fc + std::max(0.0, hmx->Drc-levee)*(1-fc);
            // equilibrium water level = weighed values of channel surplus level + hmx, levee is counted as barrier
            // can be negative if channelwh is below channel depth and low hmx level
            double cwa = ChannelWidthUpDX->Drc/ChannelAdj->Drc;

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
                        ChannelWH->Drc -= dwh;
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
                        ChannelWH->Drc += dwh/cwa;
                    }
                }
            }

            if (dosimpel)
            {
                if(whlevel > 0)
                {
                    hmx->Drc = std::min(hmx->Drc, levee) + whlevel;
                    // cutoff hmx at levee but can be smaller
                    ChannelWH->Drc = whlevel + chdepth;
                }
                else
                {

                    err += factor*hmx->Drc*ChannelAdj->Drc*DX->Drc;
                    hmx->Drc *= (1-factor);
                }
            }
        }
    }
    qDebug() << err;
    FOR_CELL_IN_FLOODAREA
            if (hmx->Drc > 0)
    {
        hmx->Drc += (err/(nrFloodedCells+1))/(ChannelAdj->Drc*DX->Drc);
        hmx->Drc = std::max(0.0, hmx->Drc);
    }}

//FOR_CELL_IN_FLOODAREA
//{
//    if (hmx->Drc > 0)
//        hmx->Drc += (err/(nrFloodedCells+1))/(DX->Drc*ChannelAdj->Drc);
//    else
//        if (hmx->Drc == -9999)
//            hmx->Drc = (err/(nrFloodedCells+1))/(DX->Drc*ChannelAdj->Drc);
//}}



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
            sum2 += M->Drc;
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
            M->Drc += dh;
            M->Drc = std::max(M->Drc , 0.0);
        }
    }
    return dh;
}
//---------------------------------------------------------------------------
void TWorld::FloodSpuriousValues()
{
    fill(*tm, 0.0);
    //calc2Maps(*tma, *DEM, *hmx, ADD);
    FOR_ROW_COL_MV
    {
        if (hmx->Drc > F_extremeHeight)
        {
            tm->Drc = getWindowAverage(*hmx, r, c, false);
           // tm->Drc = getWindowAverage(*tma, r, c, false);
        }
    }

    FOR_ROW_COL_MV
    {
        if ((hmx->Drc > F_extremeHeight*2) || (hmx->Drc > F_extremeHeight && hmx->Drc > tm->Drc + F_extremeDiff))
        {
            double htmp = hmx->Drc;
            hmx->Drc = std::min( tm->Drc, std::min(hmx->Drc, Hmx->Drc));
            qDebug() << hmx->Drc << Hmx->Drc << tm->Drc << htmp << r << c ;
        }
    }

}
//---------------------------------------------------------------------------
void TWorld::FloodBoundary()
{
    FOR_ROW_COL_MV
    {
        if (FloodEdge->Drc > 0 && hmx->Drc > 0)
        {
            //     qDebug() << Qflood->Drc*_dt/(DX->Drc*ChannelAdj->Drc);
            if( Qflood->Drc*_dt > hmx->Drc*DX->Drc*ChannelAdj->Drc)
                Qflood->Drc = hmx->Drc*DX->Drc*ChannelAdj->Drc/_dt;
            hmx->Drc = std::max(0.0, hmx->Drc - Qflood->Drc*_dt/(DX->Drc*ChannelAdj->Drc));
            floodBoundaryTot += Qflood->Drc*_dt;
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

    FOR_CELL_IN_FLOODAREA
    {
        if (hmx->Drc > 0)
        {
            hmx->Drc += 0.5*(MBeM3/(nrFloodedCells+1))/(DX->Drc*ChannelAdj->Drc);
            hmx->Drc = std::max(0.0, hmx->Drc);
        }
    }}

    ChannelOverflow();
    // mix overflow water and flood water in channel cells

    double sumh_t = mapTotal(*hmx);
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

    FloodSpuriousValues();
    //correct extremes

   // copy(*QfloodPrev, *Qflood);
    // copy the previous flood level, at t-1

    FOR_CELL_IN_FLOODAREA
    {
        //UVflood->Drc = 0.5*(fabs(Uflood->Drc)+fabs(Vflood->Drc));
        UVflood->Drc = sqrt(Uflood->Drc*Uflood->Drc+Vflood->Drc*Vflood->Drc);
        // U and V are vectors so can be negative, UV is positive average

        Qflood->Drc = UVflood->Drc * hmx->Drc * ChannelAdj->Drc;

        //AlphaFlood->Drc = std::pow(N->Drc/Grad->Drc * std::pow(2*hmx->Drc+ChannelAdj->Drc, (2/3)),0.6);
        // calc a fake alpha for sediment routing
    }}

    ChannelOverflow();
    // mix overflow water and flood water in channel cells

    correctMassBalance(sumh_t, hmx, 1e-12);
    // correct mass balance, VJ 150823: not nnecessary hhere if flow is false

    //new flood domain
    nrFloodedCells = 0;
    sumh_t = 0;
    FOR_CELL_IN_FLOODAREA
        if (hmx->Drc > 0 && FloodZonePotential->Drc == 1)
        {
            FloodDomain->Drc = 1;
            nrFloodedCells += 1.0;
            sumh_t += hmx->Drc;
        }
        else
            FloodDomain->Drc = 0;
    }

    FloodBoundary();
    // boundary flow

    copy(*Hmx, *hmx);
    // copy flood level for next dt

    FloodMaxandTiming();
    // flood max, start and duration

  //  fill(*FloodWaterVol, 0);
    FOR_CELL_IN_FLOODAREA
       FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;
    }

    //double avgh = (cells > 0 ? (sumh_t)/cells : 0);
    double area = nrFloodedCells*_dx*_dx;
    debug(QString("Flooding (dt %1 sec, n %2): area %3 m2, %4 cells").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1).arg(nrFloodedCells));
    // some screen error reporting

}
