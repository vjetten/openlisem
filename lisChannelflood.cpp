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

#include "lisemqt.h"
#include "model.h"
#include "global.h"


//---------------------------------------------------------------------------
//! Get flood level in channel from 1D kin wave channel
//! Instantaneous mixing of flood water and channel water in channel cells
//! note: ChannelDepth lets you also control which channels flood:
//! those that are 0 react as usual (infinite capacity)
void TWorld::ChannelOverflow(void)
{
    //  tmc->fill(0);

    FOR_ROW_COL_MV_CH
    {
        if (ChannelDepth->Drc > 0 && ChannelMaxQ->Drc == 0 && LDD->Drc != 5)
        {

            if (hmx->Drc == 0 && ChannelWH->Drc < ChannelDepth->Drc)
                continue;

            double fc = ChannelWidthUpDX->Drc/_dx;

            // fraction reaching the channel
            double levee = ChannelLevee->Drc;
            double chdepth = ChannelDepth->Drc + levee;

            double whlevel = (ChannelWH->Drc - chdepth)*fc + max(0, hmx->Drc-levee)*(1-fc);
            // new water level = weighed values of channel surplus level + hmx, levee is counted as barrier
            // can be negative if channelwh is below channel depth and low hmx level

            //if average water level is positive, water redistributes instantaneously and
            // hmx and channel wh are equal
            // normally this goes by a wave velocity probably because water level is flat V = sqrt(gh)
            if (whlevel > 0)
            {
                hmx->Drc = min(hmx->Drc, levee);
                // cutoff hmx at levee but can be smaller
                hmx->Drc += whlevel;
                ChannelWH->Drc = whlevel + chdepth;
            }
            else
                if (hmx->Drc > levee)
                {
                    hmx->Drc = levee;
                }

        }
        ChannelWaterVol->Drc = ChannelWH->Drc * (ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2.0 * ChannelDX->Drc;
        // recalc chjan volume for looping

    } // channel cells
}
//---------------------------------------------------------------------------
// correct mass balance
double TWorld::correctMassBalance(double sum1, TMMap *M, double minV)
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
            M->Drc = max(M->Drc , 0);
        }
    }
    return dh;
}
//---------------------------------------------------------------------------

// NOTE DEM has barriers included, done in shade map calculation !!!!
void TWorld::ChannelFlood(void)
{
    if (!SwitchIncludeChannel)
        return;
    if (!SwitchChannelFlood)
        return;

    ChannelOverflow();
    // mix overflow water and flood water in channel cells

    double sumh_t = hmx->mapTotal();
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
        dtflood = fullSWOF2Do2(hmx, Uflood, Vflood, DEM);//, q1flood, q2flood);
        FOR_ROW_COL_MV
        {
            UVflood->Drc = 0.5*(Uflood->Drc+Vflood->Drc);
            Qflood->Drc = UVflood->Drc * hmx->Drc * ChannelAdj->Drc;
        }
    }
    else
        if (SwitchFloodSWOForder1)
        {
            dtflood = fullSWOF2Do1(hmx, Uflood, Vflood, DEM);//, q1flood, q2flood);
            FOR_ROW_COL_MV
            {
                UVflood->Drc = 0.5*(Uflood->Drc+Vflood->Drc);
                Qflood->Drc = UVflood->Drc * hmx->Drc * ChannelAdj->Drc;
            }
        }
        else
            if (SwitchFloodExplicit)
            {
                dtflood = floodExplicit();
            }

    ChannelOverflow();
    // mix overflow water and flood water in channel cells

    /* double dh = */ correctMassBalance(sumh_t, hmx, 0);
    // correct mass balance

    // floodwater volume and max flood map
    FloodWaterVol->fill(0);
    FOR_ROW_COL_MV
    {
        FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;

        maxflood->Drc = max(maxflood->Drc, hmx->Drc);
        if (hmx->Drc > minReportFloodHeight)
            timeflood->Drc += _dt/60;
        // for output
    }

    floodVolTotMax = 0;
    floodAreaMax = 0;
    double area = _dx*_dx;
    FOR_ROW_COL_MV
    {
        if (maxflood->Drc > minReportFloodHeight)
        {
            floodVolTotMax += maxflood->Drc*area;
            floodAreaMax += area;
        }
    }

    //new flood domain
    double cells = 0;
    sumh_t = 0;
    FOR_ROW_COL_MV
    {
        if (hmx->Drc > 0)
        {
            FloodDomain->Drc = 1;
            cells += 1.0;
            sumh_t += hmx->Drc;
        }
        else
            FloodDomain->Drc = 0;
    }


    //double avgh = (cells > 0 ? (sumh_t)/cells : 0);
    area = cells*_dx*_dx;
    //    debug(QString("Flooding (dt %1 sec, n %2): avg h%3 m, area %4 m2").arg(dtflood,6,'f',3).arg(iter_n,4).arg(dh ,6,'e',1).arg(area,8,'f',1));
    debug(QString("Flooding (dt %1 sec, n %2): area %3 m2").arg(dtflood,6,'f',3).arg(iter_n,4).arg(area,8,'f',1));
    // some error reporting
}
