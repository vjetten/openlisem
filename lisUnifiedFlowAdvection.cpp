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


void TWorld::UF2D_Advect_Momentum(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_fu, cTMap * out_fv,  cTMap * out_su, cTMap * out_sv)
{
    FOR_ROW_COL_UF2D
    {
        out_fu->Drc = _fu->Drc;
        out_fv->Drc = _fv->Drc;
        out_su->Drc = _su->Drc;
        out_sv->Drc = _sv->Drc;
        UF_t1->Drc = _f->Drc;
        UF_t2->Drc = _s->Drc;

    }

    FOR_ROW_COL_UF2D_DT
    {
        double hf = _f->Drc/(_dx*_dx);
        double hs = _s->Drc/(_dx*_dx);

        double qfx = std::min(std::fabs(dt->Drc * hf * _fu->Drc*_dx),0.5 * UF_Courant * _f->Drc);
        double qfy = std::min(std::fabs(dt->Drc * hf * _fv->Drc*_dx),0.5 * UF_Courant * _f->Drc);

        int qfxc = c + ((_fu->Drc > 0)? -1 : 1);
        int qfyr = r + ((_fv->Drc > 0)? -1 : 1);


        if(!UF_OUTORMV(_dem,r,qfxc)){
            out_fu->data[r][qfxc] = ((UF_t1->data[r][qfxc] + qfx) > 0)? (UF_t1->data[r][qfxc] * out_fu->data[r][qfxc] + qfx * _fu->Drc)/(UF_t1->data[r][qfxc] + qfx) : 0.0;
            out_fv->data[r][qfxc] = ((UF_t1->data[r][qfxc] + qfx) > 0)? (UF_t1->data[r][qfxc] * out_fv->data[r][qfxc] + qfx * _fv->Drc)/(UF_t1->data[r][qfxc] + qfx) : 0.0;
            UF_t1->data[r][qfxc] += qfx;
            UF_t1->Drc -= qfx;
        }
        if(!UF_OUTORMV(_dem,qfyr,c)){
            out_fu->data[qfyr][c] = ((UF_t1->data[qfyr][c] + qfy) > 0)? (UF_t1->data[qfyr][c] * out_fu->data[qfyr][c] + qfy * _fu->Drc)/(UF_t1->data[qfyr][c] + qfy) : 0.0;
            out_fv->data[qfyr][c] = ((UF_t1->data[qfyr][c] + qfy) > 0)? (UF_t1->data[qfyr][c] * out_fv->data[qfyr][c] + qfy * _fv->Drc)/(UF_t1->data[qfyr][c] + qfy) : 0.0;
            UF_t1->data[qfyr][c] += qfy;
            UF_t1->Drc -= qfy;
        }

        if(UF_SOLIDPHASE)
        {
            double qsx = std::min(std::fabs(dt->Drc * hs * _su->Drc*_dx),0.5 * UF_Courant * _s->Drc);
            double qsy = std::min(std::fabs(dt->Drc * hs * _sv->Drc*_dx),0.5 * UF_Courant * _s->Drc);

            int qsxc = c + ((_su->Drc > 0)? -1 : 1);
            int qsyr = r + ((_sv->Drc > 0)? -1 : 1);

            if(!UF_OUTORMV(_dem,r,qsxc)){
                out_su->data[r][qsxc] = ((UF_t2->data[r][qsxc] + qsx) > 0)? (UF_t2->data[r][qsxc] * out_su->data[r][qsxc] + qsx * _su->Drc)/(UF_t2->data[r][qsxc] + qsx) : 0.0;
                out_sv->data[r][qsxc] = ((UF_t2->data[r][qsxc] + qsx) > 0)? (UF_t2->data[r][qsxc] * out_sv->data[r][qsxc] + qsx * _sv->Drc)/(UF_t2->data[r][qsxc] + qsx) : 0.0;
                UF_t2->data[r][qsxc] += qsx;
                UF_t2->Drc -= qsx;
            }
            if(!UF_OUTORMV(_dem,qsyr,c)){
                out_su->data[qsyr][c] = ((UF_t2->data[qsyr][c] + qsy) > 0)? (UF_t2->data[qsyr][c] * out_su->data[qsyr][c] + qsy * _su->Drc)/(UF_t2->data[qsyr][c] + qsy) : 0.0;
                out_sv->data[qsyr][c] = ((UF_t2->data[qsyr][c] + qsy) > 0)? (UF_t2->data[qsyr][c] * out_sv->data[qsyr][c] + qsy * _sv->Drc)/(UF_t2->data[qsyr][c] + qsy) : 0.0;
                UF_t2->data[qsyr][c] += qsy;
                UF_t2->Drc -= qsy;
            }
        }

    }}}

}


double TWorld::UF2D_Advect_mass(cTMap * dt, cTMap * _dem,cTMap * _m, cTMap * _mu, cTMap * _mv, cTMap * out_m)
{
    double outflow = 0;
    bool write_to_input = false;
    if(out_m ==0)
    {
        write_to_input = true;
        out_m = UF_t2;
    }

    FOR_ROW_COL_UF2D
    {
        out_m->Drc = _m->Drc;
    }

    FOR_ROW_COL_UF2D_DT
    {
        double hm = _m->Drc/(_dx*_dx);

        double qmx = std::min(std::fabs(dt->Drc * hm * _mu->Drc*_dx), 0.5 * UF_Courant * out_m->Drc);
        double qmy = std::min(std::fabs(dt->Drc * hm * _mv->Drc*_dx), 0.5 * UF_Courant * out_m->Drc);

        int qmxc = c + ((_mu->Drc > 0)? -1 : 1);
        int qmyr = r + ((_mv->Drc > 0)? -1 : 1);

        if(!UF_OUTORMV(_dem,r,qmxc)){
            out_m->data[r][qmxc] += qmx;
            out_m->data[r][c] -= qmx;
        }else{
            outflow += qmx;
            out_m->data[r][c] -= qmx;
        }
        if(!UF_OUTORMV(_dem,qmyr,c)){
            out_m->data[qmyr][c] += qmy;
            out_m->data[r][c] -= qmy;
        }else{
            outflow += qmy;
            out_m->data[r][c] -= qmy;
        }
    }}}
    if(write_to_input)
    {
        FOR_ROW_COL_UF2D
        {
            _m->Drc = out_m->Drc;
        }
    }
    return outflow;
}

void TWorld::UF2D_Advect_prop(cTMap * dt, cTMap * _dem,cTMap * _m, cTMap * _mu, cTMap * _mv,cTMap *_prop, cTMap * out_prop)
{
    bool write_to_input = false;
    if(out_prop ==0)
    {
        write_to_input = true;
        out_prop = UF_t2;
    }
    FOR_ROW_COL_UF2D
    {
        out_prop->Drc = _prop->Drc;
        UF_t1->Drc = _m->Drc;
    }

    FOR_ROW_COL_UF2D_DT
    {
        if(_m->Drc < UF_VERY_SMALL)
        {
            continue;
        }

        double hm = _m->Drc/(_dx*_dx);

        double qmx = std::min(std::fabs(dt->Drc * hm * _mu->Drc*_dx), 0.5 * UF_Courant * UF_t1->Drc);
        double qmy = std::min(std::fabs(dt->Drc * hm * _mv->Drc*_dx), 0.5 * UF_Courant * UF_t1->Drc);

        int qmxc = c + ((_mu->Drc > 0)? -1 : 1);
        int qmyr = r + ((_mv->Drc > 0)? -1 : 1);


        if(!UF_OUTORMV(_dem,r,qmxc)){
            out_prop->data[r][qmxc] = (UF_t1->data[r][qmxc] +  qmx > UF_VERY_SMALL)? (out_prop->data[r][qmxc]  * UF_t1->data[r][qmxc] + _prop->Drc * qmx)/(UF_t1->data[r][qmxc] +  qmx) : 0.0;
        }
        if(!UF_OUTORMV(_dem,qmyr,c)){
            out_prop->data[qmyr][c] = (UF_t1->data[qmyr][c] +  qmy > UF_VERY_SMALL)? (out_prop->data[qmyr][c]  * UF_t1->data[qmyr][c] + _prop->Drc * qmy)/(UF_t1->data[qmyr][c] +  qmy) : 0.0;
        }

        if(!UF_OUTORMV(_dem,r,qmxc)){
            UF_t1->data[r][qmxc] += qmx;
            UF_t1->data[r][c] -= qmx;
        }
        if(!UF_OUTORMV(_dem,qmyr,c)){
            UF_t1->data[qmyr][c] += qmy;
            UF_t1->data[r][c] -= qmy;
        }


    }}}

    if(write_to_input)
    {
        FOR_ROW_COL_UF2D
        {
            _prop->Drc = out_prop->Drc;
        }
    }
}


void TWorld::UF2D_Diffuse_mass(cTMap* dt, cTMap * _dem,cTMap * _m,cTMap * _f, cTMap * _s, cTMap * _fu, cTMap * _fv, cTMap * _su, cTMap * _sv, cTMap * out_m)
{
    bool write_to_input = false;
    if(out_m ==0)
    {
        write_to_input = true;
        out_m = UF_t2;
    }
    FOR_ROW_COL_UF2D
    {
        out_m->Drc = _m->Drc;
    }

    FOR_ROW_COL_UF2D_DT
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
        FOR_ROW_COL_UF2D
        {
            _m->Drc = out_m->Drc;
        }
    }
    return;

}


//1D version
void TWorld::UF1D_Advect_Momentum(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_fu,  cTMap * out_su )
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    FOR_ROW_COL_UF1D
    {
        out_fu->Drc = _fu->Drc;
        out_su->Drc = _su->Drc;
        UF_t1->Drc = _f->Drc;
        UF_t2->Drc = _s->Drc;
    }

    FOR_ROW_COL_UF1D_DT
    {
        double hf = _f->Drc/(_dx*_lddw->Drc);
        double hs = _s->Drc/(_dx*_lddw->Drc);

        double qfx = std::min(std::fabs(dt->Drc * hf * _fu->Drc*_lddw->Drc), 0.5 * UF_Courant * _f->Drc);
        double qsx = std::min(std::fabs(dt->Drc * hs * _su->Drc*_lddw->Drc), 0.5 * UF_Courant * _s->Drc);

        //frontward flux
        if(_fu->Drc > 0)
        {
            //front cell
            int lddself = (int) _ldd->data[r][c];
            if(!lddself == 5)
            {
                int r2 = r+dy[lddself];
                int c2 = c+dx[lddself];

                if(!UF_OUTORMV(_ldd,r2,c2)){
                    out_fu->data[r2][c2] = ((UF_t1->data[r2][c2] + qfx) > 0)? (UF_t1->data[r2][c2] * out_fu->data[r2][c2] + qfx * _fu->Drc)/(UF_t1->data[r2][c2] + qfx) : 0.0;
                    UF_t1->Drc -= qfx;
                    UF_t1->data[r2][c2] += qfx;
                }
            }
        }else
        {
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
                        double flux = qfx*_lddw->data[r2][c2]/totalwidth;

                        out_fu->data[r2][c2] = ((UF_t1->data[r2][c2] + flux) > 0)? (UF_t1->data[r2][c2] * out_fu->data[r2][c2] + flux * _fu->Drc)/(UF_t1->data[r2][c2] + flux) : 0.0;
                        UF_t1->Drc -= qfx;
                        UF_t1->data[r2][c2] += qfx;
                    }
                }
            }
        }
        //frontward flux
        if(_su->Drc > 0)
        {
            //front cell
            int lddself = (int) _ldd->data[r][c];
            if(!lddself == 5)
            {
                int r2 = r+dy[lddself];
                int c2 = c+dx[lddself];

                if(!UF_OUTORMV(_ldd,r2,c2)){
                    out_su->data[r2][c2] = ((UF_t2->data[r2][c2] + qsx) > 0)? (UF_t2->data[r2][c2] * out_su->data[r2][c2] + qsx * _su->Drc)/(UF_t2->data[r2][c2] + qsx) : 0.0;
                    UF_t2->Drc -= qsx;
                    UF_t2->data[r2][c2] += qsx;
                }
            }
        }else
        {
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
                        double flux = qsx*_lddw->data[r2][c2]/totalwidth;

                        out_su->data[r2][c2] = ((UF_t2->data[r2][c2] + flux) > 0)? (UF_t2->data[r2][c2] * out_su->data[r2][c2] + flux * _su->Drc)/(UF_t2->data[r2][c2] + flux) : 0.0;
                        UF_t2->Drc -= qsx;
                        UF_t2->data[r2][c2] += qsx;
                    }
                }
            }
        }
    }
}

double TWorld::UF1D_Advect_mass(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _m, cTMap * _mu, cTMap * out_m)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};
    double outflow = 0;

    bool write_to_input = false;
    if(out_m ==0)
    {
        write_to_input = true;
        out_m = UF_t2;
    }

    FOR_ROW_COL_UF1D
    {
        out_m->Drc = _m->Drc;
    }

    FOR_ROW_COL_UF1D_DT
    {
        double hm = _m->Drc/(_dx*_lddw->Drc);

        double qmx = std::min(std::fabs(dt->Drc * hm * _mu->Drc*_lddw->Drc), 0.5 * UF_Courant * out_m->Drc);

        //frontward flux
        if(_mu->Drc > 0)
        {
            //front cell
            int lddself = (int) _ldd->data[r][c];
            if(!lddself == 5)
            {
                int r2 = r+dy[lddself];
                int c2 = c+dx[lddself];

                if(!UF_OUTORMV(_ldd,r2,c2)){
                    out_m->data[r2][c2] += qmx;
                    out_m->data[r][c] -= qmx;
                }
            }else
            {
                outflow += qmx;
                out_m->data[r][c] -= qmx;
            }
        }else
        {
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
                        double flux = qmx*_lddw->data[r2][c2]/totalwidth;
                        out_m->data[r2][c2] += flux;
                        out_m->data[r][c] -= flux;
                    }else
                    {
                        double flux = qmx*_lddw->data[r2][c2]/totalwidth;
                        outflow += flux;
                        out_m->data[r][c] -= flux;
                    }
                }
            }
        }
    }
    if(write_to_input)
    {
        FOR_ROW_COL_UF1D
        {
            _m->Drc = out_m->Drc;
        }
    }
    return outflow;
}
void TWorld::UF1D_Advect_prop(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _m, cTMap * _mu,cTMap *_prop, cTMap * out_prop)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    bool write_to_input = false;
    if(out_prop ==0)
    {
        write_to_input = true;
        out_prop = UF_t2;
    }

    FOR_ROW_COL_UF1D
    {
        out_prop->Drc = _prop->Drc;
        UF_t1->Drc = _m->Drc;
    }

    FOR_ROW_COL_UF1D_DT
    {
        double hm = _m->Drc/(_dx*_lddw->Drc);

        double qmx = std::min(std::fabs(dt->Drc * hm * _mu->Drc*_lddw->Drc), 0.5 * UF_Courant * UF_t1->Drc);

        //frontward flux
        if(_mu->Drc > 0)
        {
            //front cell
            int lddself = (int) _ldd->data[r][c];
            if(!lddself == 5)
            {
                int r2 = r+dy[lddself];
                int c2 = c+dx[lddself];

                if(!UF_OUTORMV(_ldd,r2,c2)){
                    out_prop->data[r2][c2] = (UF_t1->data[r2][c2] +  qmx) > 0?(out_prop->data[r2][c2]  * UF_t1->data[r2][c2] + _prop->Drc * qmx)/(UF_t1->data[r2][c2] +  qmx) : _prop->Drc;

                    UF_t1->data[r2][c2] += qmx;
                    UF_t1->data[r][c] -= qmx;
                }
            }
        }else
        {
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
                        double flux = qmx*_lddw->data[r2][c2]/totalwidth;
                        out_prop->data[r2][c2] = (UF_t1->data[r2][c2] +  flux) > 0 ? (out_prop->data[r2][c2]  * UF_t1->data[r2][c2] + _prop->Drc * flux)/(UF_t1->data[r2][c2] +  flux): _prop->Drc;

                        UF_t1->data[r2][c2] += flux;
                        UF_t1->data[r][c] -= flux;
                    }
                }
            }
        }
    }
    if(write_to_input)
    {
        FOR_ROW_COL_UF1D
        {
            _prop->Drc = out_prop->Drc;
        }
    }
}




void TWorld::UF1D_Diffuse_mass(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap * _m, cTMap * _f,cTMap * _fu,cTMap * _s,cTMap * _su, cTMap * out_m)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    bool write_to_input = false;
    if(out_m ==0)
    {
        write_to_input = true;
        out_m = UF_t2;
    }

    FOR_ROW_COL_UF1D
    {
        out_m->Drc = _m->Drc;
    }

    //diffusion of Suspended Sediment layer
    FOR_ROW_COL_UF1D_DT
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
    }

    if(write_to_input)
    {
        FOR_ROW_COL_UF1D
        {
            _m->Drc = out_m->Drc;
        }
    }


}



