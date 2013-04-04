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
 \brief Channel flood using a simple solution of St Venant equations following Bates et al. 2010.\n
        called before ChannelFlow(), takes old channel overflow height and spreads it out, puts new channelWH \n
        back into channel before kin wave of channel is done in ChannelFlow()

functions: \n
- void TWorld::ChannelFlood(void) calculate maps channelflood height (hmx) and FloodDomain
*/

#include "lisemqt.h"
#include "model.h"
#include "global.h"

//---------------------------------------------------------------------------
#define Drci Data[r+dr[i]][c+dc[i]]

//---------------------------------------------------------------------------
// flood level function consists of the following parts:
// 1) get channel overflow level
// 2) do loop with varying timestep depending on max level
// in that loop
// 3) get hydraulic head and manning for flood level cells and surrounding cells in 3x3 window
// 4) calc flux to/from central cell in 8 directions and limit to wave velocity flux
// 5) update flood level with flux and correct errors: level cannot be higher than surrounding cells outside channel
// 6) calc avg water level (mass balance) error


// NOTE DEM has barriers included, done in shade map calculation !!!!
void TWorld::ChannelFlood(void)
{
    if (!SwitchIncludeChannel)
        return;
    if (!SwitchChannelFlood)
        return;

    SwitchFloodExplicit = false;//true;
    SwitchFloodSWOForder1 = true;
    SwitchFloodSWOForder2 = false;//false;

    // get flood level in channel from 1D kin wave channel
    FOR_ROW_COL_MV_CH
    {
        if (ChannelDepth->Drc > 0 && ChannelMaxQ->Drc == 0)
        {
            double whsurp = max(0, ChannelWH->Drc - ChannelDepth->Drc);
            hmx->Drc = whsurp*ChannelWidthUpDX->Drc/_dx + hmx->Drc*(_dx-ChannelWidthUpDX->Drc)/_dx;
        }
        // note: ChannelDepth lets you also control which channels flood: those that are 0 react as usual
    }

    double sumh_t = 0, sumh_t1 = 0, diff = 0, cells = 1, dtflood;

    sumh_t = hmx->mapTotal();
    // sum all levels for mass balance


    if (SwitchFloodExplicit)
    {
        dtflood = floodExplicit();
    }
    if (SwitchFloodSWOForder2)
    {
        dtflood = fullSWOF2D(hmx, Uflood, Vflood, DEM, q1flood, q2flood);
    }

    if (SwitchFloodSWOForder1)
    {
        dtflood = fullSWOF2Do1(hmx, Uflood, Vflood, DEM, q1flood, q2flood);
    }

    FOR_ROW_COL_MV
    {
        if (hmx->Drc > 0)
            cells++;
    }

    sumh_t1 = hmx->mapTotal();
    diff = (sumh_t - sumh_t1)/cells;
    double avgh = sumh_t1/cells;
    debug(QString("Flooding (dt %3): avg h%1, avg err h %2 m").arg(avgh,8,'f',3).arg(diff,8,'e',3).arg(dtflood,6,'f',3));
    // some error reporting

    // put new flood level in channel for next 1D kin wave channel
    FOR_ROW_COL_MV_CH
    {
        if (hmx->Drc > 0)// && tmc->Drc == 1)// && ChannelMaxQ->Drc == 0)
        {
            ChannelWH->Drc = hmx->Drc + ChannelDepth->Drc;
            if (ChannelMaxQ->Drc > 0)
            {
                ChannelWH->Drc = min(ChannelDepth->Drc, ChannelWH->Drc);
                hmx->Drc = 0;
            }

            ChannelWaterVol->Drc = ChannelWH->Drc *
                    ((ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2) * ChannelDX->Drc;
        }
    }

    // floodwater volume and max flood map
    FloodWaterVol->fill(0);
    FOR_ROW_COL_MV
    {
        if (ChannelDepth->Drc == 0)
        {
            FloodWaterVol->Drc = hmx->Drc*_dx*DX->Drc;
        }
        maxflood->Drc = max(maxflood->Drc, hmx->Drc);
        // for output
    }

    maxflood->report("maxflood.map");
}

/*
    // if there is no flood skip everything
    if (startflood)
    {
        int dc[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
        int dr[9] = {1,1,1,  0, 0, 0,  -1, -1, -1};
        double courant_number = courant_factor;
        double froude_limit = 0.8;
        double gravity = 9.81;
        double h_min = 1e-6;
        double timestep = 0;
        double timestep1 = 0;
        double timesum = 0;
        // do one _dt with varying timestep based on courant condition
        do {
            // make Hydraulic head, gravity (dem+barriers) + water level
            FOR_ROW_COL_MV
            {
                Hmx->Drc = DEM->Drc + Barriers->Drc + hmx->Drc;
            }

            double maxdepth = qMax(0.01, hmx->mapMaximum());
            // find maxdepth
            timestep = courant_number*_dx/qSqrt(gravity*maxdepth);
            timestep = qMax(0.001, timestep);
            timestep = qMax(minFloodDt, timestep);
            timestep1 = timestep;
            timestep = qMin(timestep, _dt-timesum);
            // determine timestep

            Qxsum->fill(0);
            // map with fluxes to/from central cell

            // flag which cells need processing, flood domain
            // this is the cell where hmx > 0 and any cell adjacent even if dry
            tma->fill(0);
            for (int i = 0; i < 9; i++)
            {
                FOR_ROW_COL_MV
                        if (r+dr[i] > 0 && r+dr[i] < _nrRows &&
                            c+dc[i] > 0 && c+dc[i] < _nrCols &&
                            !IS_MV_REAL8(&hmx->Drci))
                {
                    if (hmx->Drci > 0)
                        tma->Drc = 1;
                    // flag
                }
            }

            // prepare maps in direction i;
            for (int i = 1; i < 9; i+=2)
            {
                double _dx2 = ((i+1) % 2 == 0 ? _dx : _dx*qSqrt(2));
                // diagonal cells have sqrt(2) dx

                // get the pressure and N in 8 directions, and fill with 0 if outside boundaries
                FOR_ROW_COL_MV
                {
                    if (r+dr[i] > 0 && r+dr[i] < _nrRows &&
                            c+dc[i] > 0 && c+dc[i] < _nrCols &&
                            !IS_MV_REAL8(&hmx->Drci))
                    {
                        Hx->Drc = DEM->Drci + hmx->Drci;
                        hx->Drc = hmx->Drci;
                        Nx->Drc = N->Drci;
                        dHdLx->Drc = (Hx->Drc - Hmx->Drc)/_dx2;
                    }
                    else
                    {
                        Hx->Drc = 0;
                        hx->Drc = 0;
                        Nx->Drc = N->Drc;
                        dHdLx->Drc = 0;
                    }
                }

                // calc and sum fluxes in 8 directions, 4 = central cell
                FOR_ROW_COL_MV
                        if (tma->Drc == 1)
                {
                    double Qx = 0, qlx = 0, qlx1 = 0, qlx2 = 0;

                    double dHdLxi = dHdLx->Drc;
                    double signx = (dHdLxi < 0 ? -1.0 : 1.0);
                    double hxi = hx->Drc;
                    double qxi = qx[i].m->Drc;
                    double NN = Nx->Drc;

                    if (hxi > h_min)
                        qlx = (qxi - (gravity*hxi*timestep*dHdLxi))/
                                (1.0 + (gravity*hxi*timestep*NN*NN*qxi)/qPow(hxi, 10.0/3.0));
                    else
                        qlx = 0;
                    // explicit solution of saint venant equation (m2/s), Bates et al. 2010
                    qlx1 = hxi*qSqrt(gravity*hxi)*froude_limit;
                    // limit max flux to h * wave velocity * froude number (m * m/s = m2/s)
                    qlx2 = (hxi - hmx->Drc)*_dx/timestep;
                    // limit max flux to width * water difference with central cell (w*dh/dt  in m2/s)
                    Qx = signx * qMin(qMin(qAbs(qlx1), qAbs(qlx)),qAbs(qlx2));
                    // Qx is the min of all possible fluxes, preserve sign
                    if (i != 4)
                        Qxsum->Drc += Qx/_dx2;//*0.5;
                    // using 4 directions
                    //                    if (i == 1 || i == 3 || i == 5 || i == 7)
                    //                          Qxsum->Drc += Qx/_dx;
                    // add fluxes of 8 directions, qsum has unit m2/s /m = m/s
                    // 0.5 accounts for the fact that the central cell has a boundary
                    // of 4 sides is 4*dx, touching only EW and NS directions. Adding
                    // four more diagonal fluxes would cause twice the flow, so all fluxes
                    // are assumed to have a width of 0.5*dx * 8 = 4 dx.
                    // using 8 instead of 4 directions seems to give a much better flow
                    // for 4 directions do i += 2 instead of i++

                    if (i == 4) // central cell
                        qx[i].m->Drc = signx * qMin(qlx, qlx1); //qlx2 is always 0 for i = 4
                    else
                        qx[i].m->Drc = Qx;
                    // save flux in direction i for next flood timestep
                } // for all flood cells
            } // for i = 1 to 9


            // add sum q*dt (= m/s * s) to h and h > 0
            cells = 1;
            FOR_ROW_COL_MV
                    if (tma->Drc == 1)
            {
                hmx->Drc += timestep*Qxsum->Drc;
                hmx->Drc = max(0, hmx->Drc);
                if (hmx->Drc > 0)
                    cells++;
            }
/*
//            FOR_ROW_COL_MV
//                    Hmx->Drc = DEM->Drc + Barriers->Drc + hmx->Drc;
            // new hydraulic head

            FOR_ROW_COL_MV
                    if (tma->Drc == 1)
            {
                double hmax = 0;
//                double Hmax = 0;
                for (int i = 0; i < 9; i++)
                    if (r+dr[i] > 0 && r+dr[i] < _nrRows &&
                            c+dc[i] > 0 && c+dc[i] < _nrCols &&
                            i != 4 &&                     // not the centre cell
                            !IS_MV_REAL8(&hmx->Drci))
                    {
//                        Hmax = qMax(hmax, Hmx->Drci);
                        hmax = qMax(hmax, hmx->Drci);
                        // find the highest water level around centre cell
                    }
                if (ChannelDepth->Drc == 0)
                {
//                    Hmx->Drc = qMin(Hmx->Drc, Hmax);
//                    hmx->Drc = qMin(hmx->Drc, hmax);
                }
                // correct hmx if not channel cell
            }
            //simple correction,outside channels flow cannot lead to water level rise above neighbours outside channel,
            //water cannot flow uphill in this simplified solution (no momentum)
            //if we use Hydraulic Head H (=h+z), instabilities occur !!!

//            FOR_ROW_COL_MV
//                    hmx->Drc = max(0, Hmx->Drc - DEM->Drc - Barriers->Drc);
*/
/*
            FOR_ROW_COL_MV
                    if (ChannelMaxQ->Drc > 0)
            {
                hmx->Drc = 0;
            }
            // no flood in culvert cells

            // find current flood domain (hmx > 0) and nr of flooded cells
            // used in the other processes, infiltration, runoff etc
            cells = 1;
            FOR_ROW_COL_MV
            {
                if (hmx->Drc > 0)
                {
                    FloodDomain->Drc = 1;
                    cells++;
                }
                else
                    FloodDomain->Drc = 0;
            }

            timesum = timesum + timestep;
            // sum to reach _dt

        } while (timesum  < _dt);
        // continue while _dt is not reached
    }

    sumh_t1 = hmx->mapTotal();
    diff = (sumh_t - sumh_t1)/cells;
    double avgh = sumh_t1/cells;
    debug(QString("Flooding (dt %3): avg h%1, avg err h %2 m").arg(avgh,8,'f',3).arg(diff,8,'e',3).arg(timestep1,6,'f',3));
    // some error reporting

    Qflood->fill(0);
    FOR_ROW_COL_MV
            if(hmx->Drc > 0)// && ChannelDepth->Drc == 0)
    {
        Vflood->Drc = qPow(hmx->Drc, 0.667)*qSqrt(Grad->Drc)/N->Drc;
        Qflood->Drc = Vflood->Drc * hmx->Drc * _dx;
        // estimate resulting flux simply by manning
    }

    // put new flood level in channel for next 1D kin wave channel
    FOR_ROW_COL_MV_CH
    {
        if (hmx->Drc > 0)// && tmc->Drc == 1)// && ChannelMaxQ->Drc == 0)
        {
            ChannelWH->Drc = hmx->Drc + ChannelDepth->Drc;
            if (ChannelMaxQ->Drc > 0)
            {
                ChannelWH->Drc = min(ChannelDepth->Drc, ChannelWH->Drc);
                hmx->Drc = 0;
            }

            ChannelWaterVol->Drc = ChannelWH->Drc *
                    ((ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2) * ChannelDX->Drc;
        }
    }

    FloodWaterVol->fill(0);
    FOR_ROW_COL_MV
    {
        if (ChannelDepth->Drc == 0)
        {
            FloodWaterVol->Drc = hmx->Drc*_dx*DX->Drc;
        }
        maxflood->Drc = max(maxflood->Drc, hmx->Drc);
        // for output
    }

    //    FOR_ROW_COL_MV
    //    {
    //        if(FloodDomain->Drc == 1)
    //        {
    //            WH->Drc = hmx->Drc;
    //            WHroad->Drc = hmx->Drc;
    //        }
    //    }
    // put flood water level back in overland flow water level in flooddomain

    //    hmx->report("hmx");
    //    Qflood->report("Qf");
    //    Vflood->report("Vf");
    maxflood->report("maxflood.map");
}
*/


double TWorld::floodExplicit()//TMMap *hmx, TMMap *Vflood, TMMap *DEM, TMMap *Qflood)
{
    int n = 1;
    double timesum = 0;

    // get flood level in channel from 1D kin wave channel
    bool startflood = false;

    FOR_ROW_COL_MV
    {
        if (hmx->Drc > 0)
        {
            startflood = true;
            break;
        }
    }

    // if there is no flood skip everything
    if (startflood)
    {
        int dc[4] = {0, -1, 1, 0};
        int dr[4] = {-1, 0,  0,1};
        double courant_number = courant_factor;
        double froude_limit = 0.8;
        double gravity = 9.81;
        double h_min = 1e-6;
        double timestep = 0;

        // do one _dt with varying timestep based on courant condition
        do {
            n++;
            // make Hydraulic head, gravity (dem+barriers) + water level
            // barriers are already in the dem
            Hmx->calc2Maps(DEM, hmx, ADD);

            double maxdepth = qMax(0.01, hmx->mapMaximum());
            // find maxdepth

            timestep = courant_number*_dx/qSqrt(gravity*maxdepth);
            // determine timestep
            timestep = qMax(0.001, timestep);
            timestep = qMin(timestep, _dt-timesum);

            Qxsum->fill(0);
            // map with fluxes to/from central cell

            // flag which cells need processing, flood domain
            // this is the cell where hmx > 0 and any cell adjacent even if dry
            tma->fill(0);
            for (int i = 0; i < 9; i++)
            {
                FOR_ROW_COL_MV
                        if (r+dr[i] > 0 && r+dr[i] < _nrRows &&
                            c+dc[i] > 0 && c+dc[i] < _nrCols &&
                            !IS_MV_REAL8(&hmx->Drci))
                {
                    if (hmx->Drci > 0)
                        tma->Drc = 1;
                    // flag
                }
            }

            // prepare maps in direction i;
            for (int i = 0; i < 4; i++)
            {

                // get the pressure and N in 8 directions, and fill with 0 if outside boundaries
                FOR_ROW_COL_MV
                {
                    if (r+dr[i] > 0 && r+dr[i] < _nrRows &&
                            c+dc[i] > 0 && c+dc[i] < _nrCols &&
                            !IS_MV_REAL8(&hmx->Drci))
                    {
                        Hx->Drc = DEM->Drci + hmx->Drci;
                        hx->Drc = hmx->Drci;
                        Nx->Drc = N->Drci;
                        dHdLx->Drc = (Hx->Drc - Hmx->Drc)/_dx;
                    }
                    else
                    {
                        Hx->Drc = 0;
                        hx->Drc = 0;
                        Nx->Drc = N->Drc;
                        dHdLx->Drc = 0;
                    }
                }

                // calc and sum fluxes in 8 directions, 4 = central cell
                FOR_ROW_COL_MV
                        if (tma->Drc == 1)
                {
                    double Qx = 0, qlx = 0, qlx1 = 0, qlx2 = 0;

                    double dHdLxi = dHdLx->Drc;
                    double signx = (dHdLxi < 0 ? -1.0 : 1.0);
                    double hxi = hx->Drc;
                    double qxi = 0;
                    double NN = Nx->Drc;

                    if (i == 0) qxi = qx0->Drc;
                    if (i == 1) qxi = qx1->Drc;
                    if (i == 2) qxi = qx2->Drc;
                    if (i == 3) qxi = qx3->Drc;

                    if (hxi > h_min)
                        qlx = (qxi - (gravity*hxi*timestep*dHdLxi))/
                                (1.0 + (gravity*hxi*timestep*NN*NN*qxi)/qPow(hxi, 10.0/3.0));
                    else
                        qlx = 0;
                    // explicit solution of saint venant equation (m2/s), Bates et al. 2010
                    qlx1 = hxi*qSqrt(gravity*hxi)*froude_limit;
                    // limit max flux to h * wave velocity * froude number (m * m/s = m2/s)
                    qlx2 = 1e6;//(hxi - hmx->Drc)*_dx/timestep;
                    // limit max flux to width * water difference with central cell (w*dh/dt  in m2/s)
                    Qx = signx * qMin(qMin(qAbs(qlx1), qAbs(qlx)),qAbs(qlx2));
                    // Qx is the min of all possible fluxes, preserve sign
                    Qxsum->Drc += Qx/_dx;

                    if (i == 0) qx0->Drc = Qx;
                    if (i == 1) qx1->Drc = Qx;
                    if (i == 2) qx2->Drc = Qx;
                    if (i == 3) qx3->Drc = Qx;
                    // save flux in direction i for next flood timestep
                } // for all flood cells
            } // for i = 1 to 9


            // add sum q*dt (= m/s * s) to h and h > 0
            FOR_ROW_COL_MV
                    if (tma->Drc == 1)
            {
                hmx->Drc += timestep*Qxsum->Drc;
                hmx->Drc = max(0, hmx->Drc);
                if (ChannelMaxQ->Drc > 0)
                    hmx->Drc = 0;
                // no flood in culvert cells
            }


            // find current flood domain (hmx > 0) and nr of flooded cells
            // used in the other processes, infiltration, runoff etc
            FOR_ROW_COL_MV
            {
                if (hmx->Drc > 0)
                {
                    FloodDomain->Drc = 1;
                }
                else
                    FloodDomain->Drc = 0;
            }

            timesum = timesum + timestep;
            // sum to reach _dt

        } while (timesum  < _dt);
        // continue while _dt is not reached
    }

    Qflood->fill(0);
    FOR_ROW_COL_MV
            if(hmx->Drc > 0)// && ChannelDepth->Drc == 0)
    {
        Vflood->Drc = qPow(hmx->Drc, 0.667)*qSqrt(Grad->Drc)/N->Drc;
        Qflood->Drc = Vflood->Drc * hmx->Drc * _dx;
        // estimate resulting flux simply by manning
    }

    return(timesum/(n-1));
}



