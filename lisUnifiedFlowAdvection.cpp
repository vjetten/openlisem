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
 \file lisUnifiedFlowAdvection.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"

void TWorld::UF2D_Diffuse_mass(int thread,cTMap* dt, cTMap * _dem,cTMap * _m,cTMap * _f, cTMap * _fu, cTMap * _fv, cTMap * out_m)
{
    bool write_to_input = false;
    if(out_m ==0)
    {
        write_to_input = true;
        out_m = ThreadPool->UF_t2.at(thread);
    }
    FOR_ROW_COL_UF2DMT
    {
        out_m->Drc = _m->Drc;
    }}}

    FOR_ROW_COL_UF2DMT_DT
    {
        //cell sizes
        double cdx = DX->Drc;
        double cdy = _dx;
        //here it is about spacing, not flow width, so use _dx instead of ChannelAdj->Drc

        //mixing coefficient
        double dux1 = c  > 0 ? std::abs(_fu->data[r][c] - _fu->data[r][c-1]) : 0;
        double dvy1 = r  > 0 ? std::abs(_fv->data[r][c] - _fv->data[r-1][c]) : 0;
        double dvx1 = c  > 0 ? std::abs(_fv->data[r][c] - _fv->data[r][c-1]) : 0;
        double duy1 = r  > 0 ? std::abs(_fu->data[r][c] - _fu->data[r-1][c]) : 0;
        double dux2 = c  < _nrCols-1 ? std::abs(_fu->data[r][c+1] - _fu->data[r][c]) : 0;
        double dvy2 = r  < _nrRows-1 ? std::abs(_fv->data[r+1][c] - _fv->data[r][c]) : 0;
        double dvx2 = c  < _nrCols-1 ? std::abs(_fv->data[r][c+1] - _fv->data[r][c]) : 0;
        double duy2 = r  < _nrRows-1 ? std::abs(_fu->data[r+1][c] - _fu->data[r][c]) : 0;

        //drop term if at boundary
        if(std::isnan(dux1))
            dux1 = 0;
        if(std::isnan(dvx1))
            dvx1 = 0;
        if(std::isnan(duy1))
            duy1 = 0;
        if(std::isnan(dvy1))
            dvy1 = 0;

        if(std::isnan(dux2))
            dux2 = 0;
        if(std::isnan(dvx2))
            dvx2 = 0;
        if(std::isnan(duy2))
            duy2 = 0;
        if(std::isnan(dvy2))
            dvy2 = 0;

        double dux = std::max(dux1,dux2);
        double dvy = std::max(dvy1,dvy2);
        double dvx = std::max(dvx1,dvx2);
        double duy = std::max(duy1,duy2);

        double eddyvs = cdx * cdy * sqrt(dux*dux + dvy*dvy +  0.5 * (dvx +duy)*(dvx +duy));
        double eta = eddyvs/UF_SigmaDiffusion;

        //cell directions
        int dx[4] = {0, 1, -1, 0};
        int dy[4] = {1, 0, 0, -1};
        double h = _f->Drc/(cdx*cdy);

        //use the calculated weights to distribute flow
        for (int i=0; i<4; i++)
        {
            int r2, c2;

            //must multiply the cell directions by the sign of the slope vector components
            r2 = r+dy[i];
            c2 = c+dx[i];

            if(!UF_OUTORMV(_dem,r2,c2))
            {
                double h2 = _f->data[r2][c2]/(cdx*cdy);
                double coeff = std::min(dt->Drc*eta *std::min(1.0,h2/h),UF_Courant/4.0) * _m->Drc;

                out_m->data[r2][c2] += coeff;
                out_m->data[r][c] -= coeff;
            }
        }
    }}}

    if(write_to_input)
    {
        FOR_ROW_COL_UF2DMT
        {
            _m->Drc = out_m->Drc;
        }}}
    }
    return;

}


void TWorld::UF1D_Diffuse_mass(int thread,cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap * _m, cTMap * _f,cTMap * _fu, cTMap * out_m)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    bool write_to_input = false;
    if(out_m ==0)
    {
        write_to_input = true;
        out_m = ThreadPool->UF_t2.at(thread);
    }

    FOR_ROW_COL_UF1DMT
    {
        out_m->Drc = _m->Drc;
    }}}

    //diffusion of Suspended Sediment layer
    FOR_ROW_COL_UF1DMT_DT
    {
        double dux = UF1D_Derivative(_ldd,_lddw,_fu,r,c);
        double cdx = DX->Drc;
        //diffusion coefficient according to J.Smagorinski (1964)
        double eddyvs = cdx * dux;
        //and devide by turbulent prandtl-smidth number
        double eta = eddyvs/UF_SigmaDiffusion;

        double h = _f->Drc / (_lddw->Drc * DX->Drc);
        //frontward flux

        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!lddself == 5)
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){

                double h2 = _f->data[r2][c2] / (_lddw->data[r2][c2] * DX->data[r2][c2]);
                //calculate diffusive fluxes in channel cell
                double coeff = std::min(dt->Drc*eta *std::min(1.0,h2/h),UF_Courant/2.0) *_m->Drc;
                out_m->data[r2][c2] += coeff;
                out_m->data[r][c] -= coeff;
            }
        }

        //back cells
        double totalwidth = 0;
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
                totalwidth += _lddw->data[r2][c2];
            }
        }
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
                if(!UF_OUTORMV(_ldd,r2,c2)){

                    double h2 = _f->data[r2][c2] / (_lddw->data[r2][c2] * DX->data[r2][c2]);
                    //calculate diffusive fluxes in channel cell
                    double coeff = std::min(dt->Drc*eta *std::min(1.0,h2/h),UF_Courant/2.0) *_m->Drc * _lddw->data[r2][c2]/totalwidth;
                    out_m->data[r2][c2] += coeff;
                    out_m->data[r][c] -= coeff;
                }
            }
        }
    }}}

    if(write_to_input)
    {
        FOR_ROW_COL_UF1DMT
        {
            _m->Drc = out_m->Drc;
        }}}
    }


}



