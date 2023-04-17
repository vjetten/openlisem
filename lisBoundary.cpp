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

        if(OUTORMV(r,c+1)) // returns true if outside rows. cols or mv
        {
            if(demx1 < demx2)
                K2DOutlets->Drc = 1;
        }
        if(OUTORMV(r,c-1))
        {
            if(demx2 <demx1)
                K2DOutlets->Drc = 1;
        }

        if(OUTORMV(r+1,c))
        {
            if(demy1 < demy2)
                K2DOutlets->Drc = 1;
        }
        if(OUTORMV(r-1,c))
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

        if(OUTORMV(r,c+1) && OUTORMV(r,c-1))
        {
            Dhx = 0;
            K2DOutlets->Drc = 1;
        }
        if(OUTORMV(r+1,c) && OUTORMV(r-1,c))
        {
            Dhy = 0;
            K2DOutlets->Drc = 1;
        }

        //at boundaries, set cell as outflow cell when slope is in the direction of the boundary

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
    if (FlowBoundaryType == 2) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            K2DOutlets->Drc *= FlowBoundary->Drc;
        }}
    }
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

    dynOutflowPoints();
    // find all points flowing to outside because of water level

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (K2DOutlets->Drc == 1)
        {
            if (c > 0 && MV(r,c-1)) // U = x; V = y
                if (Uflood->Drc < 0) {
                    tma->Drc = 1;
                }
            if (c < _nrCols-1 && MV(r,c+1))
                if (Uflood->Drc > 0) {
                    tma->Drc = 1;
                }
            if (r > 0 && MV(r-1,c))
                if (Vflood->Drc < 0) {
                    tma->Drc = 1;
                }
            if (r < _nrRows-1 && MV(r+1,c))
                if (Vflood->Drc > 0) {
                    tma->Drc = 1;
                }
        }
    }}

    FOR_ROW_COL_LDD5 {
        double _q = Qout.at(i_);
        double dh = _q*_dt/CHAdjDX->Drc;
        h->Drc = std::max(0.0,h->Drc-dh);

        Q->Drc = _q;

        if (SwitchErosion) {
            double ds = std::min(SSFlood->Drc, SSCFlood->Drc*_q*_dt);
            SSFlood->Drc -= ds;
            if (SwitchUse2Phase) {
                ds = std::min(BLFlood->Drc, BLCFlood->Drc*_q*_dt);
                BLFlood->Drc -= ds;
            }
        }
    }}

    if (FlowBoundaryType == 0)
        return;


    BoundaryQ = 0;
    BoundaryQs = 0;

    #pragma omp parallel for reduction(+:BoundaryQ, BoundaryQs) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (tma->Drc == 1 && h->Drc > HMIN) {
            double _q = Q->Drc;
            double dh = _q*_dt/CHAdjDX->Drc;
            if (h->Drc-dh < 0)
                dh = h->Drc;
            h->Drc -= dh;

            BoundaryQ += _q;
            Q->Drc = _q;

            if (SwitchErosion) {
                double ds = std::min(SSFlood->Drc, SSCFlood->Drc*_q*_dt);
                BoundaryQs += ds/_dt; //in kg/s
                SSFlood->Drc -= ds;
                if (SwitchUse2Phase) {
                    ds = std::min(BLFlood->Drc, BLCFlood->Drc*_q*_dt);
                    BoundaryQs += ds/_dt;
                    BLFlood->Drc -= ds;
                }
            }
        }
    }}
}
