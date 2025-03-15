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
void TWorld::Boundary2Ddyn()
{
    QBoundary = 0;
    QsBoundary = 0;

    Fill(*tma,0);
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (FlowBoundary->Drc > 0) {
            //flow left bpoundary to the left etc
            if (c-1 >= 0 && MV(r,c-1) && !MV(r,c+1)) {
                if (Uflood->Drc < 0)
                    tma->Drc = 1;
            }
            if (c+1 <= _nrCols-1 && MV(r,c+1) && !MV(r,c-1)) {
                if (Uflood->Drc > 0)
                    tma->Drc = 1;
            }
            if (r-1 >= 0 && MV(r-1,c) && !MV(r+1,c)) {
                if (Vflood->Drc < 0)
                    tma->Drc = 1;
            }
            if (r+1 <= _nrRows-1 && MV(r+1,c) && !MV(r-1,c)) {
                if (Vflood->Drc > 0)
                    tma->Drc = 1;
            }
        }
    }}

    FOR_ROW_COL_MV_L {
        if (tma->Drc == 1) {
            double Q = Qn->Drc;
            Q = std::min(Qn->Drc, (WaterVolall->Drc-MicroStoreVol->Drc)/_dt);
            //WaterVolall->Drc -= Q*_dt;

            QBoundary += Q;
            // Qn based on vector combination Uflood and Vflood, calculated before

            if (SwitchErosion) {
                double ds = std::min(SSFlood->Drc, SSCFlood->Drc*Qn->Drc*_dt);
                // because concentrations can be spurious take the min of the two
                QsBoundary += ds/_dt; //in kg/s
                if (SwitchUse2Phase) {
                    ds = std::min(BLFlood->Drc, BLCFlood->Drc*Qn->Drc*_dt);
                    QsBoundary += ds/_dt;
                }
            }
        }
    }}
    qDebug() << "boundary flux m3/s" << QBoundary << QsBoundary;
}


//---------------------------------------------------------------------------
void TWorld::Boundary2DdynUV(cTMap * U, cTMap *V)
{
    Fill(*tma,0);
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (FlowBoundary->Drc > 0) {
            //flow left bpoundary to the left etc
            if (c-1 >= 0 && MV(r,c-1) && !MV(r,c+1)) {
                if (U->Drc < 0)
                    tma->Drc = 1;
            }
            if (c+1 <= _nrCols-1 && MV(r,c+1) && !MV(r,c-1)) {
                if (U->Drc > 0)
                    tma->Drc = 1;
            }
            if (r-1 >= 0 && MV(r-1,c) && !MV(r+1,c)) {
                if (V->Drc < 0)
                    tma->Drc = 1;
            }
            if (r+1 <= _nrRows-1 && MV(r+1,c) && !MV(r-1,c)) {
                if (V->Drc > 0)
                    tma->Drc = 1;
            }
        }
    }}

    FOR_ROW_COL_MV_L {
        if (tma->Drc == 1) {
          //  double Q = sqrt(U->Drc*U->Drc + V->Drc*V->Drc) *h->Drc * ChannelAdj->Drc;
          //  Q = std::min(Qn->Drc, (WaterVolall->Drc-MicroStoreVol->Drc)/_dt);
            //WaterVolall->Drc -= Q*_dt;

       //     QBoundary += Q;
            // Qn based on vector combination Uflood and Vflood, calculated before

            // if (SwitchErosion) {
            //     double ds = std::min(SSFlood->Drc, SSCFlood->Drc*Qn->Drc*_dt);
            //     // because concentrations can be spurious take the min of the two
            //     QsBoundary += ds/_dt; //in kg/s
            //     if (SwitchUse2Phase) {
            //         ds = std::min(BLFlood->Drc, BLCFlood->Drc*Qn->Drc*_dt);
            //         QsBoundary += ds/_dt;
            //     }
            // }
        }
    }}
}
// OBSOLETE

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
