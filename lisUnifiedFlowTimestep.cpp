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
 \file lisUnifiedFlow.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"


////Timestep functions
double TWorld::UF_InitTimeStep(cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                               cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                               cTMap * _fu1D,cTMap * _s1D,
                               cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                               cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                               cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                               cTMap * _su2D,cTMap * _sv2D,cTMap * out_t2d,cTMap * out_dt2d, cTMap * out_dtstep2d,cTMap * out_t1d,cTMap * out_dt1d, cTMap * out_dtstep1d)
{
    FOR_ROW_COL_UF2D
    {
        out_t2d->Drc = 0;
        out_dt2d->Drc = 0;
        out_dtstep2d->Drc = 0;
    }

    FOR_ROW_COL_UF1D
    {
        out_t1d->Drc = 0;
        out_dt1d->Drc = 0;
        out_dtstep1d->Drc = 0;
    }


}

double TWorld::UF_TimeStep(double t, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                           cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                           cTMap * _fu1D,cTMap * _s1D,
                           cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                           cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                           cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                           cTMap * _su2D,cTMap * _sv2D,cTMap * out_t2d,cTMap * out_dt2d, cTMap * out_dtstep2d,cTMap * out_t1d,cTMap * out_dt1d, cTMap * out_dtstep1d)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};
    UF_DTMIN = _dt;
    fill(*UF2D_CellR,-1);
    fill(*UF2D_CellC,-1);
    FOR_ROW_COL_UF2D
    {
        if(_s2D->Drc + _f2D->Drc > UF_VERY_SMALL)
        {
            double v = std::max(std::max(std::max(_fu2D->Drc,_fv2D->Drc),_su2D->Drc),_sv2D->Drc);
            out_dt2d->Drc = (UF_Courant *_dx/v);

        }else
        {
            out_dt2d->Drc = _dt;
        }

        UF_DTMIN = std::max(UF2D_MinimumDT,std::min(out_dt2d->Drc,UF_DTMIN));
        out_t2d->Drc = 0;
        out_dtstep2d->Drc = 0;

        if(!UF_OUTORMV(_ldd,r,c))
        {
            if(_s1D->Drc + _f1D->Drc > UF_VERY_SMALL)
            {
                double v = std::max(_fu1D->Drc,_su1D->Drc);
                out_dt1d->Drc = (UF_Courant *_dx/v);

            }else
            {
                out_dt1d->Drc = _dt;
            }

            UF_DTMIN = std::max(UF1D_MinimumDT,std::min(out_dt1d->Drc,UF_DTMIN));
            out_t1d->Drc = 0;
            out_dtstep1d->Drc = 0;


        }
    }

    FOR_ROW_COL_UF2D
    {
        UF_t1->Drc = std::max(1.0,std::floor(std::min(out_dt2d->Drc,_dt)/UF_DTMIN));
        if(!UF_OUTORMV(_ldd,r,c))
        {
            UF_t2->Drc = std::max(1.0,std::floor(std::min(out_dt1d->Drc,_dt)/UF_DTMIN));
        }
    }

    FOR_ROW_COL_UF2D
    {
        double dt = UF_t1->Drc;
        double dt1 = !UF_OUTORMV(_dem,r+1,c)? UF_t1->data[r+1][c]:dt;
        double dt2 = !UF_OUTORMV(_dem,r-1,c)? UF_t1->data[r-1][c]:dt;
        double dt3 = !UF_OUTORMV(_dem,r,c+1)? UF_t1->data[r][c+1]:dt;
        double dt4 = !UF_OUTORMV(_dem,r,c-1)? UF_t1->data[r][c-1]:dt;

        out_dtstep2d->Drc = std::min(dt,std::min(dt1,std::min(dt2,std::min(dt3,dt4))));

        if(!UF_OUTORMV(_ldd,r,c))
        {
            double dt = UF_t2->Drc;
            for (int i=1;i<=9;i++)
            {
                int r2, c2, ldd = 0;
                if (i==5)  // Skip current cell
                    continue;
                r2 = r+dy[i];
                c2 = c+dx[i];
                if (!UF_OUTORMV(_ldd,r2,c2))
                    ldd = (int) _ldd->data[r2][c2];
                else
                    continue;
                if (!UF_OUTORMV(_ldd,r2,c2) &&
                        FLOWS_TO(ldd, r2,c2,r,c))
                {
                    dt = std::min(dt,UF_t2->data[r2][c2]);
                }
            }
            out_dtstep1d->Drc = std::min(dt,std::min(dt1,std::min(dt2,std::min(dt3,dt4))));
        }
    }
    int rc = 0;
    int cc = 0;
    FOR_ROW_COL_UF2D
    {
        if(!(t + UF_DTMIN < _dt ))
        {
            out_dt2d->Drc = _dt - out_t2d->Drc;
        }else if(out_t2d->Drc + out_dtstep2d->Drc* UF_DTMIN < t)
        {
            out_dt2d->Drc = t - out_t2d->Drc + out_dtstep2d->Drc* UF_DTMIN;
            out_dt2d->Drc = std::min(out_dt2d->Drc,_dt - out_t2d->Drc);
        }else
        {
            out_dt2d->Drc = 0;
        }

        out_t2d->Drc += out_dt2d->Drc;

        if(out_dt2d->Drc > 0)
        {
            UF2D_CellR->data[rc][cc] = r;
            UF2D_CellC->data[rc][cc] = c;
            cc ++;
            if(cc == _nrCols)
            {
                rc ++;
                cc = 0;
            }
        }

        if(!UF_OUTORMV(_ldd,r,c))
        {

            if(!(t + UF_DTMIN < _dt ))
            {
                out_dt1d->Drc = _dt - out_t1d->Drc;
            }else if(out_t1d->Drc + out_dtstep1d->Drc* UF_DTMIN < t)
            {
                out_dt1d->Drc = t - out_t1d->Drc + out_dtstep1d->Drc* UF_DTMIN;
                out_dt1d->Drc = std::min(out_dt1d->Drc,_dt - out_t1d->Drc);
            }else
            {
                out_dt1d->Drc = 0;
            }

            out_t1d->Drc += out_dt1d->Drc;
        }
    }

    return std::min(_dt- t,UF_DTMIN);
}