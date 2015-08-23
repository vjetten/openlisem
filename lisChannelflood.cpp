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

#define MINHMX 0.0
//---------------------------------------------------------------------------
//! Get flood level in channel from 1D kin wave channel
//! Instantaneous mixing of flood water and channel water in channel cells
//! note: ChannelDepth lets you also control which channels flood:
//! those that are 0 react as usual (infinite capacity)
void TWorld::ChannelOverflow(bool flow)
{
    FOR_ROW_COL_MV_CH
    {
        if (ChannelDepth->Drc > 0 && ChannelMaxQ->Drc == 0 && LDD->Drc != 5)// && FloodZonePotential->Drc > 0)
        {
            double levee = ChannelLevee->Drc;
            double chdepth = ChannelDepth->Drc + levee; // levee always assumed on both sides channel
            double dH = std::max(0.0, (ChannelWH->Drc-chdepth));

            if (dH == 0 && hmx->Drc <= levee)
                continue;
            // no flow activity then continue

            double fracA = std::min(1.0, _dt*UVflood->Drc/(0.5*ChannelAdj->Drc));
            // fraction from hmx to channel based on avefrage flood velocity
            double fracC = std::min(1.0, _dt*sqrt(ChannelWH->Drc*9.8067)/(0.5*ChannelWidthUpDX->Drc));
            // fraction from channel to surrounding adj area based on gravity flow perpedicular to channel
            double fc = ChannelWidthUpDX->Drc/_dx; // 1-fc = (dx-chw)/dx = chanadj/dx
            // fraction of the channel in the gridcell
            double whlevel = (ChannelWH->Drc - chdepth)*fc + std::max(0.0, hmx->Drc-levee)*(1-fc);
            // equilibrium water level = weighed values of channel surplus level + hmx, levee is counted as barrier
            // can be negative if channelwh is below channel depth and low hmx level

            //State 1: channel water is higher than surroundin
            if (dH > hmx->Drc)
            {
                double dHC = (dH-hmx->Drc)*fracC; // diff between channel WH and surrounding hmx
                qDebug() << dH << whlevel << hmx->Drc << dH-dHC << fracC;
                if ((dH - dHC) < hmx->Drc)   // if too much outflow, equilibrium level
                {
                    hmx->Drc = std::min(hmx->Drc, levee);
                    hmx->Drc +=  whlevel;
                    ChannelWH->Drc = whlevel + chdepth;

                }
                else
                    if (flow) // flow out of the channel
                {
                    // gives huge mass balance errors, must be false: gives huge mass balance errors!
                    ChannelWH->Drc -= dHC;
                    hmx->Drc += dHC*ChannelWidthUpDX->Drc/ChannelAdj->Drc;
                }
            }
            else  //flow in
            {
                double dhmx = (hmx->Drc-dH)*fracA; // diff between channel WH and surrounding hmx
                if (hmx->Drc - dhmx < dH) // if results in hmx dropping below channelWH
                {
                    hmx->Drc = std::min(hmx->Drc, levee);
                    hmx->Drc +=  whlevel;
                    ChannelWH->Drc = whlevel + chdepth;
                }
                else
                    if (flow)  // flow into channel, must be false: gives huge mass balance errors!
                {
                    hmx->Drc -= dhmx;
                    ChannelWH->Drc += dhmx*ChannelAdj->Drc/ChannelWidthUpDX->Drc;
                }
            }

/*
ALL IS ALWAYS IN EQUILIBRIUM

            if(whlevel > 0)
            {
                hmx->Drc = std::min(hmx->Drc, levee);
                hmx->Drc +=  whlevel;
                // cutoff hmx at levee but can be smaller
                ChannelWH->Drc = whlevel + chdepth;
            }
            else
            {
//                double frac = std::min(1.0, _dt*UVflood->Drc/std::max(0.01*_dx,0.5*ChannelAdj->Drc));

//                dH = std::min(0.0,frac*(hmx->Drc-levee));

//                ChannelWH->Drc += dH*ChannelAdj->Drc/ChannelWidthUpDX->Drc;
//                hmx->Drc = std::max(0.0, hmx->Drc - dH);

          //      qDebug() << dH << hmx->Drc << frac;

            }
*/
        }
    }
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
    // toal and cells active for M

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
    FOR_ROW_COL_MV
    {
        if (hmx->Drc > F_extremeHeight)
        {
            tm->Drc = getWindowAverage(*hmx, r, c);
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
    //  tm->copy(hmx);
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
    // tm->calcMap(hmx, SUB);
    //  tm->report("diffh");

    //volumme weg is sum diff * cell size
}
//---------------------------------------------------------------------------
void TWorld::FloodMaxandTiming()
{
    // floodwater volume and max flood map
    fill(*FloodWaterVol, 0.0);
    FOR_ROW_COL_MV
    {
       // FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;

        floodHmxMax->Drc = std::max(floodHmxMax->Drc, hmx->Drc);
        if (hmx->Drc > 0)//minReportFloodHeight)
            floodTime->Drc += _dt/60;
        // for output
        floodVMax->Drc = std::max(floodVMax->Drc, UVflood->Drc);
        // max velocity
    }

    floodVolTotMax = 0;
    floodAreaMax = 0;
    double area = _dx*_dx;
    FOR_ROW_COL_MV
    {
        if (floodHmxMax->Drc > 0)//minReportFloodHeight)
        {
            floodVolTotMax += floodHmxMax->Drc*area;
            floodAreaMax += area;
        }
    }

    FOR_ROW_COL_MV
    {
    //    if (hmx->Drc > minReportFloodHeight && floodTimeStart->Drc == 0)
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

    ChannelOverflow(false);
    // mix overflow water and flood water in channel cells

    double sumh_t = mapTotal(*hmx);
    double dtflood = 0;

    startFlood = false;
    FOR_ROW_COL_MV
    {
        if (hmx->Drc > 0)
        {
            startFlood = true;
            break;
        }
    }

    if (SwitchFloodSWOForder2)
    {
        dtflood = fullSWOF2Do2(hmx, Uflood, Vflood, DEM);
    }
    else
        if (SwitchFloodSWOForder1)
        {
            dtflood = fullSWOF2Do1(hmx, Uflood, Vflood, DEM);
        }
//        else
//            if (SwitchFloodExplicit)
//            {
//                dtflood = floodExplicit();
//            }

    FloodSpuriousValues();
    //correct extremes

    FOR_ROW_COL_MV
    {
        tmb->Drc = 0;
        UVflood->Drc = 0.5*(fabs(Uflood->Drc)+fabs(Vflood->Drc));
        // U and V are vectors so can be negative, UV is positive average
        Qflood->Drc = UVflood->Drc * hmx->Drc * ChannelAdj->Drc;
        tmb->Drc = (UVflood->Drc + sqrt(hmx->Drc * 9.8))*dtflood/_dx;
    }

    ChannelOverflow(false);
    // mix overflow water and flood water in channel cells

  //  correctMassBalance(sumh_t, hmx, 1e-12);
    // correct mass balance, VJ 150823: not nnecessary hhere if flow is false

    //new flood domain
    double cells = 0;
    sumh_t = 0;
    FOR_ROW_COL_MV
    {
        if (hmx->Drc > 0 && FloodZonePotential->Drc == 1)
        {
            FloodDomain->Drc = 1;
            cells += 1.0;
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

 //   fill(*FloodWaterVol, 0);
    FOR_ROW_COL_MV
    {
        FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;
    }

    //double avgh = (cells > 0 ? (sumh_t)/cells : 0);
    double area = cells*_dx*_dx;
    debug(QString("Flooding (dt %1 sec, n %2): area %3 m2, %4 cells").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1).arg(cells));
    // some error reporting


}
