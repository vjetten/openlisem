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
/**
 * @fn void TWorld::DEMFB()
 * @brief Returns the digital elevation model height, with the addition of flow barriers
 *
 * @param r : row number
 * @param c : column number
 * @param rd : row direction (-1 for top, 1 for bottom)
 * @param cd : column direction (-1 for left, 1 for right)
 * @param addwh : include water height for overland flow
 * @return digital elevation model height, with the addition of flow barriers
 * @see K2DDEMA
 */
double TWorld::DEMFB(int r, int c, int rd, int cd, bool addwh)
{
    cTMap *h = WHrunoff;
    if(SwitchKinematic2D == K2D_METHOD_KINDYN) {
        h = hmx;
    }

    double wh = 0;
    double dem = 0;
    if(INSIDE(r+rd,c+cd)) {
        if(!pcr::isMV(LDD->data[r+rd][c+cd]))
        {
            if(addwh)
                wh = h->data[r + rd][c + cd];
            dem = DEM->data[r + rd][c + cd];
        } else {
            if(!pcr::isMV(LDD->data[r][c])) {
                wh = 0;
                dem = DEM->Drc;
            } else {
               return 0;
            }
        }

    } else
        if(INSIDE(r,c)) {
        if(!pcr::isMV(LDD->Drc))
        {
            wh = 0;
            dem = DEM->Drc;
        } else {
           return 0;  // returns always zero because demb r c is inside and not mv
        }

    } else {
        return 0;
    }

    if(OUTORMV(r+rd,c+cd))
    {
        return dem;
    }


    if(rd == 0 && cd == 0)
    {
        return dem + wh;
    }

    if(rd == 1 && cd == 0)
    {
        return dem + std::max(wh,(FlowBarrierS->Drc));
    }
    // else if(rd == 1 && cd == 1)
    // {
    //     return dem + std::max(wh,(std::max(std::max(FlowBarrierS->Drc,FlowBarrierE->Drc),std::max(FB(r,c +cd,0,rd),FB(r+rd,c,cd,0)))));
    // }
    else if(rd == 0 && cd == 1)
    {
        return dem + std::max(wh,(FlowBarrierE->Drc));
    }
    // else if(rd == -1 && cd == 1)
    // {
    //     return dem + std::max(wh,(std::max(std::max(FlowBarrierN->Drc,FlowBarrierE->Drc),std::max(FB(r,c  +cd,0,rd),FB(r+rd,c,cd,0)))));
    // }
    else if(rd == -1 && cd == 0)
    {
        return dem + std::max(wh,(FlowBarrierN->Drc));
    }
    // else if(rd == -1 && cd == -1)
    // {
    //     return dem + std::max(wh,(std::max(std::max(FlowBarrierN->Drc,FlowBarrierW->Drc),std::max(FB(r,c  +cd,0,rd),FB(r+rd,c,cd,0)))));
    // }
    else if(rd == 0 && cd == -1)
    {
        return dem + std::max(wh,(FlowBarrierW->Drc));
    }else
    //     if(rd == 1 && cd == -1)
    // {
    //     return dem + std::max(wh,(std::max(std::max(FlowBarrierS->Drc,FlowBarrierW->Drc),std::max(FB(r,c +cd,0,rd),FB(r+rd,c,cd,0)))));
    // }
    return 0;
}

//---------------------------------------------------------------------------
// all points that flow outward of the domain by slope and water pressure
void TWorld::dynOutflowPoints(cTMap *h)
{
    //if boundary = 0 only outflow on pits
    if (FlowBoundaryType == 0)
        return;

    // for boundary 1 or 2, find all outflow points
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        K2DOutlets->Drc = 0;
        if (DomainEdge->Drc) {
            double Dhx = 0;
            double Dhy = 0;

            //DEM + water height if true
            //double dem = DEMFB(r,c,0,0,true);

            // double demx1 = DEMFB(r,c,0,1,true); //look right
            // double demx2 = DEMFB(r,c,0,-1,true); // look left
            // double demy1 = DEMFB(r,c,1,0,true); // look up
            // double demy2 = DEMFB(r,c,-1,0,true); // look down

            double dem = DEM->Drc + h->Drc;
            // domainedge can get ws number?
            double demx1, demx2, demy1, demy2 = 0;
            // find the inland value of dem and h on the not mv side of the edge, edge = r,c
            int situation = 0;
            if (DomainEdge->Drc > 0) {
                // look left, not MV
                int rr = r;
                int cr = c-1;
                if (!pcr::isMV(LDD->Drcr)) {
                    demx1 = DEM->Drcr + h->Drcr;
                    situation = 1;
                }
                // look right not MV
                rr = r;
                cr = c+1;
                if (pcr::isMV(LDD->Drcr)) {
                    demx2 = DEM->Drcr + h->Drcr;
                    situation = 2;
                }
                // look up not MV
                rr = r-1;
                cr = c;
                if (pcr::isMV(LDD->Drcr)) {
                    demy1 = DEM->Drcr + h->Drcr;
                    situation = 3;
                }
                // look down not MV
                rr = r+1;
                cr = c;
                if (pcr::isMV(LDD->Drcr)) {
                    demy1 = DEM->Drcr + h->Drcr;
                    situation = 4;
                }
            }

            if (situation == 1 && demx1 > dem) {
                K2DOutlets->Drc = 1;
            }
            if (situation == 2 && demx1 > dem) {
                K2DOutlets->Drc = 1;
            }
            if (situation == 3 && demy1 > dem) {
                K2DOutlets->Drc = 1;
            }
            if (situation == 4 && demy2 > dem) {
                K2DOutlets->Drc = 1;
            }

            /*
            if (OUTORMV(r,c+1))
            {
                if(demx1 < demx2)
                    K2DOutlets->Drc = 1;
            }
            if ((OUTORMV(r,c-1))
            {
                if(demx2 < demx1)
                    K2DOutlets->Drc = 1;
            }

            if( OUTORMV(r+1,c))
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

            if ( OUTORMV(r,c+1) && OUTORMV(r,c-1))
            {
                Dhx = 0;
                K2DOutlets->Drc = 1;
            }
            if ((OUTORMV(r+1,c) && OUTORMV(r-1,c))
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
            */
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
void TWorld::Boundary2Ddyn(cTMap *h, cTMap *u, cTMap *v)
{

    // cTMap *Q = Qn;
    // cTMap *h = WHrunoff;
    // if(SwitchKinematic2D == K2D_METHOD_KINDYN) {
    //     Q = Qflood;
    //     h = hmx;
    // }

    if (FlowBoundaryType == 0)
        return;

    dynOutflowPoints(h);
    // find all points flowing to outside because of water level
    // includes effect of boundary condition 2 (user defined)

    FOR_ROW_COL_LDD5 {
        K2DOutlets->Drc = 0;
    }}

    BoundaryQ = 0;
    BoundaryQs = 0;

    //#pragma omp parallel for reduction(+:BoundaryQ, BoundaryQs) num_threads(userCores)

    // do not subtract outgoing flux from the volume, klike in the kin wave this is not necessary the flux is there and goes out
    // and will be calculte din the mass balance as outgoing
    FOR_ROW_COL_MV_L {
        if (K2DOutlets->Drc == 1 && h->Drc > 0) {
            BoundaryQ += Q->Drc;

            if (SwitchErosion) {
                double ds = std::min(SSFlood->Drc, SSCFlood->Drc*Q->Drc*_dt);
                BoundaryQs += ds/_dt; //in kg/s
              //  SSFlood->Drc -= ds;
                if (SwitchUse2Phase) {
                    ds = std::min(BLFlood->Drc, BLCFlood->Drc*Q->Drc*_dt);
                    BoundaryQs += ds/_dt;
                  //  BLFlood->Drc -= ds;
                }
            }
        }
    }}

    //qDebug() << BoundaryQ << MB;
}
