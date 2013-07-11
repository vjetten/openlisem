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
    //   tm->fill(0);
    FOR_ROW_COL_MV_CH
    {
        if (ChannelDepth->Drc > 0 && ChannelMaxQ->Drc == 0)
        {
            double fc = min(1.0, ChannelWidthUpDX->Drc/_dx);
            double levee = ChannelLevee->Drc;

            double whlevel = (ChannelWH->Drc - ChannelDepth->Drc - levee)*fc +
                    max(0, hmx->Drc-levee)*(1-fc);
            //average water level
            //if average water level is positive, water redistributes instantaneously and
            // hmx and channel wh are equal
            if (whlevel > 0)
            {
                hmx->Drc = min(hmx->Drc, levee);
                hmx->Drc += whlevel+levee;
                ChannelWH->Drc = whlevel + ChannelDepth->Drc + levee;
            }
            // if average water level is negative, channel wh < depth, but there can be hmx
            // some flood water moves into the channel accoridng to flood velocity
            // if fv = 1 and whlevel < 0 hmx always 0, see excel
            else
                if (hmx->Drc > levee)
                {
                    //                    double fv;
                    //                    double dh = hmx->Drc;
                    //                    fv = 1.0;//ChannelAdj->Drc > 0 ? min(_dt*UVflood->Drc/(0.5*ChannelAdj->Drc), 1.0) : 1.0;

                    hmx->Drc = levee;
                    //                   ChannelWH->Drc += dh*fv*ChannelAdj->Drc/ChannelWidthUpDX->Drc;

                }
        }
    }
    // tm->report("whl");
}
//---------------------------------------------------------------------------
// correct mass balance
double TWorld::correctMassBalance(double sum1, TMMap *M)
{
    double sum2 = 0;
    double n = 0;
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0.0)
        {
            n += 1;
            sum2 += M->Drc;
        }
    }
    double dh = (n > 0 ? (sum1 - sum2)/n : 0);
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0.0)
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
        dtflood = fullSWOF2Do2(hmx, Uflood, Vflood, DEM, q1flood, q2flood);
        FOR_ROW_COL_MV
        {
            UVflood->Drc = 0.5*(Uflood->Drc+Vflood->Drc);
            Qflood->Drc = UVflood->Drc * hmx->Drc * ChannelAdj->Drc;
        }
    }
    else
        if (SwitchFloodSWOForder1)
        {
            dtflood = fullSWOF2Do1(hmx, Uflood, Vflood, DEM, q1flood, q2flood);
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

    double dh = correctMassBalance(sumh_t, hmx);
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

    FOR_ROW_COL_MV_CH
    {
        maxChannelflow->Drc = max(maxChannelflow->Drc, ChannelQn->Drc);
        maxChannelWH->Drc = max(maxChannelWH->Drc, ChannelWH->Drc);
    }

    maxflood->report("maxflood.map");
    timeflood->report("timeflood.map");
    maxChannelflow->report("maxhannelq.map");
    maxChannelWH->report("maxhannelwh.map");

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

    double avgh = (cells > 0 ? (sumh_t)/cells : 0);
    area = cells*_dx*_dx;
    debug(QString("Flooding (dt %1 sec, n %2): avg h%3 m, area %4 m2").arg(dtflood,6,'f',3).arg(iter_n,4).arg(dh ,8,'e',3).arg(area,8,'f',1));
    // some error reporting
}

