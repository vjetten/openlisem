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
 \file lisChannelfloodexpl.cpp
 \brief Channel flood using an explicit solution of the St Venant equations following Bates et al. 2010.\n

functions: \n
- double TWorld::floodExplicit() Do flooding with explicit solution Bates et al (old lisflood)
*/

#include "lisemqt.h"
#include "model.h"
#include "global.h"


//---------------------------------------------------------------------------
// explicit LISFLOOD solution Bates et al journal of hydrology 2010

double TWorld::floodExplicit()
{
    int n = 0;
    double timesum = 0;

    // if there is no flood skip everything
    if (startFlood)
    {
        int dc[4] = {0, -1, 1, 0};
        int dr[4] = {-1, 0,  0,1};
        double courant_number = courant_factor;
        double froude_limit = 0.8;
        double gravity = 9.81;
        double hmin = 1e-6;
        double timestep = 0;

        // do one _dt with varying timestep based on courant condition
        do {
            n++;
            // make Hydraulic head, gravity (dem+barriers) + water level
            // barriers are already in the dem
            Hmx->calc2Maps(DEM, hmx, ADD);

            double maxdepth = _max(0.01, hmx->mapMaximum());
            // find maxdepth
            double maxv = _max(0.01, Vflood->mapMaximum());

            timestep = courant_number*_dx/qSqrt(gravity*maxdepth);
            timestep = courant_number*_dx/(maxv+qSqrt(gravity*maxdepth));
            // determine timestep
            timestep = _max(0.001, timestep);
            timestep = _min(timestep, _dt-timesum);

            Qxsum->fill(0);
            // map with fluxes to/from central cell

            // flag which cells need processing, flood domain
            // this is the cell where hmx > 0 and any cell adjacent even if dry
            tma->fill(0);
            for (int i = 0; i < 4; i++)
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

                    if (hxi > hmin)
                        qlx = (qxi - (gravity*hxi*timestep*dHdLxi))/
                                (1.0 + (gravity*hxi*timestep*NN*NN*qxi)/qPow(hxi, 10.0/3.0));
                    else
                        qlx = 0;
                    // explicit solution of saint venant equation (m2/s), Bates et al. 2010
                    qlx1 = hxi*qSqrt(gravity*hxi)*froude_limit;
                    // limit max flux to h * wave velocity * froude number (m * m/s = m2/s)
                    qlx2 = 1e6;//(hxi - hmx->Drc)*_dx/timestep;
                    // limit max flux to width * water difference with central cell (w*dh/dt  in m2/s)
                    Qx = signx * _min(_min(qAbs(qlx1), qAbs(qlx)),qAbs(qlx2));
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
                hmx->Drc = _max(0.0, hmx->Drc);
                if (ChannelMaxQ->Drc > 0)
                    hmx->Drc = 0;
                // no flood in culvert cells
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
        UVflood->Drc = qPow(hmx->Drc, 0.667)*qSqrt(Grad->Drc)/N->Drc;
        Qflood->Drc = UVflood->Drc * hmx->Drc * DX->Drc;
        // estimate resulting flux simply by manning
    }

    iter_n = n;
    return(timesum/(n+1));
}

/*
8 DIRECTION SOLUTION
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

                    if (hxi > hmin)
                        qlx = (qxi - (gravity*hxi*timestep*dHdLxi))/
                                (1.0 + (gravity*hxi*timestep*NN*NN*qxi)/qPow(hxi, 10.0/3.0));
                    else
                        qlx = 0;
                    // explicit solution of saint venant equation (m2/s), Bates et al. 2010
                    qlx1 = hxi*qSqrt(gravity*hxi)*froude_limit;
                    // limit max flux to h * wave velocity * froude number (m * m/s = m2/s)
                    qlx2 = (hxi - hmx->Drc)*_dx/timestep;
                    // limit max flux to width * water difference with central cell (w*dh/dt  in m2/s)
                    Qx = signx * _min(_min(qAbs(qlx1), qAbs(qlx)),qAbs(qlx2));
                    // Qx is the min of all possible fluxes, preserve sign
                    if (i != 4)
                        Qxsum->Drc += Qx/_dx2;// *0.5;
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
                        qx[i].m->Drc = signx * _min(qlx, qlx1); //qlx2 is always 0 for i = 4
                    else
                        qx[i].m->Drc = Qx;
                    // save flux in direction i for next flood timestep
                } // for all flood cells
            } // for i = 1 to 9
*/
//---------------------------------------------------------------------------
