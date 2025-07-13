/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

#include "model.h"


//---------------------------------------------------------------------------
// just use the flux outward if the sign of u or v is pointing outward at the boundary
// nothing fancy, this is the most stable
void TWorld::Boundary2Ddyn(double dt, cTMap *h, cTMap *u, cTMap *v)
{
    QBoundary = 0;
    QsBoundary = 0;

    // TODO barriers!
//    #pragma omp parallel for num_threads(userCores)
    #pragma omp parallel for reduction(+:QBoundary, QsBoundary) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (FlowBoundary->Drc > 0) {
            int flag = 0;
            double count = 0;
            double Qbflux1 = 0;
            double Qbflux2 = 0;
            double Qbflux3 = 0;
            double Qbflux4 = 0;
            double Area = h->Drc*ChannelAdj->Drc;

            if (c > 0 && c < _nrCols-1 ) {
                if (MV(r,c-1) && !MV(r,c+1) && u->Drc < 0) {
                    flag  = 1;
                    count += 1.0;
                    Qbflux1 = -u->Drc*Area;
                }
                if (MV(r,c+1) && !MV(r,c-1) && u->Drc > 0) {
                    Qbflux2 = u->Drc*Area;
                    count += 1.0;
                    flag = 2;
                }
            }

            if (r > 0 && r < _nrRows-1) {
                if (MV(r-1,c) && !MV(r+1,c) && v->Drc < 0) {
                    Qbflux3 = -v->Drc*Area;
                    count += 1.0;
                    flag = 3;
                }
                if (MV(r+1,c) && !MV(r-1,c) && v->Drc > 0) {
                    Qbflux4 = v->Drc*Area;
                    count += 1.0;
                    flag = 4;
                }
            }


            if (flag > 0) {
                double Qbflux = sqrt(u->Drc*u->Drc + v->Drc*v->Drc)*Area;
                h->Drc = qMax(0.0, h->Drc - Qbflux*dt/CHAdjDX->Drc);
                //adjust boundary cells

                QBoundary += Qbflux;
                QBoundFlow->Drc = Qbflux; // not used anywhere !!!! use it to sum watershed boundary flow later
                if (SwitchErosion) {
                    double ds = qMin(SSFlood->Drc, SSCFlood->Drc*QBoundFlow->Drc*dt);
                    SSFlood->Drc -= ds;
                    QsBoundary += ds/dt; //in kg/s

                    if (SwitchUse2Phase) {
                        ds = qMin(BLFlood->Drc, BLCFlood->Drc*QBoundFlow->Drc*dt);
                        BLFlood->Drc -= ds;
                        QsBoundary += ds/dt;
                    }
                }
            } // flag
        } //flowboundary
    }}

//    qDebug() << "boundary flux m3/s" << QBoundary << "kg/s" << QsBoundary;
}
