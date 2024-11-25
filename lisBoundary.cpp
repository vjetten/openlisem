/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/
#include "model.h"

//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::OverlandFlow(void)
 * @brief Calls the kinematic wave or diffusive wave functions and calculates new discharge, water height and sediment presence
 *
 * Calls the kinematic, diffusive or dynamic wave functions and calculates new discharge, water height and sediment presence
 * During this process, surpluss potential infilration is subtracted from the water content.
 * Based on the options in the run file, either the 1D or 2D kinematic wave is used.
 * Sediment transport in overland flow is automatically taken into accaunt.
 */


//---------------------------------------------------------------------------
// all points that flow outward of the domain by slope and water pressure
void TWorld::dynOutflowPoints()
{
    //if boundary = 0 only outflow on pits
    if (FlowBoundaryType == 0)
        return;

    // for boundary 1 or 2, find all outflow points
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double Dhx = 0;
        double Dhy = 0;

        //DEM + water height and barriers if switched on
        double dem = DEMFB(r,c,0,0,true);

        double demx1 = DEMFB(r,c,0,1,true); //look right
        double demx2 = DEMFB(r,c,0,-1,true); // look left
        double demy1 = DEMFB(r,c,1,0,true);
        double demy2 = DEMFB(r,c,-1,0,true);

        if (FlowBoundary->Drc == 6)//  OUTORMV(r,c+1))
        {
            if(demx1 < demx2)
                K2DOutlets->Drc = 1;
        }
        if (FlowBoundary->Drc == 4)// (OUTORMV(r,c-1))
        {
            if(demx2 < demx1)
                K2DOutlets->Drc = 1;
        }

        if(FlowBoundary->Drc == 2) // OUTORMV(r+1,c))
        {
            if(demy1 < demy2)
                K2DOutlets->Drc = 1;
        }
        if(FlowBoundary->Drc == 8)   //OUTORMV(r-1,c))
        {
            if(demy2 < demy1)
                K2DOutlets->Drc = 1;
        }

        if(demx1 < demx2)
        {
            Dhx = -(demx1-dem);
        }else
        {
            Dhx = (demx2-dem);
        }

        if(demy1 < demy2)
        {
            Dhy = -(demy1-dem);
        }else
        {
            Dhy = (demy2-dem);
        }

        if (FlowBoundary->Drc == 4 && FlowBoundary->Drc == 6)// OUTORMV(r,c+1) && OUTORMV(r,c-1))
        {
            Dhx = 0;
            K2DOutlets->Drc = 1;
        }
        if (FlowBoundary->Drc == 2 && FlowBoundary->Drc == 8)//(OUTORMV(r+1,c) && OUTORMV(r-1,c))
        {
            Dhy = 0;
            K2DOutlets->Drc = 1;
        }

        //at corners, set cell as outflow cell when slope is in the direction of the boundary

        if(r == 0)
        {
            if( Dhy < 0)
            {
                K2DOutlets->Drc = 1;
            }
        }

        if(r == _nrRows-1)
        {
            if( Dhy > 0)
            {
               K2DOutlets->Drc = 1;
            }
        }

        if(c == 0)
        {
            if( Dhx < 0)
            {
                K2DOutlets->Drc = 1;
            }
        }

        if(c == _nrCols-1)
        {
            if( Dhx > 0)
            {
                K2DOutlets->Drc = 1;
            }
        }
    }}

    //flowboundary 2 use the map
//    if (FlowBoundaryType == 2) {
//        #pragma omp parallel for num_threads(userCores)
//        FOR_ROW_COL_MV_L {
//            K2DOutlets->Drc *= FlowBoundary->Drc;
//        }}
//    }
}
//---------------------------------------------------------------------------
void TWorld::Boundary2Ddyn()
{

    cTMap *Q = Qn;
    cTMap *h = WHrunoff;
    if(SwitchKinematic2D == K2D_METHOD_KINDYN) {
        Q = Qflood;
        h = hmx;
    }
    // correct the water height in the outlet(s) for a perfect WB!
//    FOR_ROW_COL_LDD5 {
//        double dh = Q->Drc*_dt/CHAdjDX->Drc;
//        h->Drc = std::max(0.0,h->Drc-dh);

//        if (SwitchErosion) {
//            double ds = std::min(SSFlood->Drc, SSCFlood->Drc*Q->Drc*_dt);
//            SSFlood->Drc -= ds;
//            if (SwitchUse2Phase) {
//                ds = std::min(BLFlood->Drc, BLCFlood->Drc*Q->Drc*_dt);
//                BLFlood->Drc -= ds;
//            }
//        }
//    }}

    if (FlowBoundaryType == 0)
        return;

    dynOutflowPoints();
    // find all points flowing to outside because of water level
    // includes effect of boundary condition 2 (user defined)

    Fill(*K2DOutlets,0);

//    FOR_ROW_COL_MV_L {
//        tma->Drc = DEM->Drc + h->Drc;
//    }}
//    // outlets already done
//    FOR_ROW_COL_LDD5 {
//        tma->Drc = 0;
//    }}

    //NOTE Uflood negative is flow to the left, positive to the right, u = x col, v = y row
    //Vflood negative is flow up, positive is flow down
    // 2,4,6,8, are ldd directions
    FOR_ROW_COL_MV_L {
        if (FlowBoundary->Drc == 4 && Uflood->Drc < 0) {
            K2DOutlets->Drc = 1;
        }
        if (FlowBoundary->Drc == 6 && Uflood->Drc > 0) {
            K2DOutlets->Drc = 1;
        }
        if (FlowBoundary->Drc == 2 && Vflood->Drc > 0) {
            K2DOutlets->Drc = 1;
        }
        if (FlowBoundary->Drc == 8 && Vflood->Drc < 0) {
            K2DOutlets->Drc = 1;
        }
    }}

    FOR_ROW_COL_LDD5 {
        K2DOutlets->Drc = 0;
    }}

    BoundaryQ = 0;
    BoundaryQs = 0;

    //#pragma omp parallel for reduction(+:BoundaryQ, BoundaryQs) num_threads(userCores)

    FOR_ROW_COL_MV_L {
        if (K2DOutlets->Drc == 1 && h->Drc > HMIN) {
            double dh = Q->Drc*_dt/CHAdjDX->Drc;

//            if (dh > h->Drc) {
//                dh = h->Drc;
//                Q->Drc = dh/_dt*CHAdjDX->Drc;
//            }

//            h->Drc -= dh;
            BoundaryQ += Q->Drc;

            if (SwitchErosion) {
                double ds = std::min(SSFlood->Drc, SSCFlood->Drc*Q->Drc*_dt);
                BoundaryQs += ds/_dt; //in kg/s
                SSFlood->Drc -= ds;
                if (SwitchUse2Phase) {
                    ds = std::min(BLFlood->Drc, BLCFlood->Drc*Q->Drc*_dt);
                    BoundaryQs += ds/_dt;
                    BLFlood->Drc -= ds;
                }
            }
        }
    }}

    //qDebug() << BoundaryQ << MB;
}
