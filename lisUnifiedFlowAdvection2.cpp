
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

void TWorld::UF2D_Advect2_Momentum(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_fu, cTMap * out_fv,  cTMap * out_su, cTMap * out_sv,cTMap * out_qfx1,cTMap *out_qfx2,cTMap * out_qfy1,cTMap *out_qfy2,cTMap * out_qsx1,cTMap *out_qsx2,cTMap * out_qsy1,cTMap *out_qsy2)
{
    FOR_ROW_COL_UF2D
    {
        out_fu->Drc = _fu->Drc;
        out_fv->Drc = _fv->Drc;
        out_su->Drc = _su->Drc;
        out_sv->Drc = _sv->Drc;
        UF_t3->Drc = _f->Drc;
        UF_t4->Drc = _s->Drc;
    }

    FOR_ROW_COL_UF2D_DT
    {
        double qx1 = out_qfx1->Drc;
        double qx2 = out_qfx2->Drc;
        double qy1 = out_qfy1->Drc;
        double qy2 = out_qfy2->Drc;

        if(!UF_OUTORMV(_dem,r,c+1))
        {
            if(qx1 > 0)
            {
                out_fu->data[r][c+1] = ((UF_t3->data[r][c+1] + qx1) > UF_VERY_SMALL)? (UF_t3->data[r][c+1] * out_fu->data[r][c+1] + qx1 * _fu->Drc)/(UF_t3->data[r][c+1] + qx1) : 0.0;
                out_fv->data[r][c+1] = ((UF_t3->data[r][c+1] + qx1) > UF_VERY_SMALL)? (UF_t3->data[r][c+1] * out_fv->data[r][c+1] + qx1 * _fv->Drc)/(UF_t3->data[r][c+1] + qx1) : 0.0;
                UF_t3->data[r][c+1] += qx1;
                UF_t3->Drc -= qx1;
            }else if(qx1 != 0)
            {
                qx1 = std::fabs(qx1);
                out_fu->Drc = ((UF_t3->Drc + qx1) > UF_VERY_SMALL)? (UF_t3->Drc * out_fu->Drc + qx1 * _fu->data[r][c+1])/(UF_t3->Drc + qx1) : 0.0;
                out_fv->Drc = ((UF_t3->Drc + qx1) > UF_VERY_SMALL)? (UF_t3->Drc * out_fv->Drc + qx1 * _fv->data[r][c+1])/(UF_t3->Drc + qx1) : 0.0;
                UF_t3->data[r][c+1] -= qx1;
                UF_t3->Drc += qx1;
            }
        }
        if(!UF_OUTORMV(_dem,r,c-1)){

            if(qx2 < 0)
            {
                qx2 = std::fabs(qx2);
                out_fu->data[r][c-1] = ((UF_t3->data[r][c-1] + std::fabs(qx2)) > UF_VERY_SMALL)? (UF_t3->data[r][c-1] * out_fu->data[r][c-1] + std::fabs(qx2) * _fu->Drc)/(UF_t3->data[r][c-1] + std::fabs(qx2)) : 0.0;
                out_fv->data[r][c-1] = ((UF_t3->data[r][c-1] + std::fabs(qx2)) > UF_VERY_SMALL)? (UF_t3->data[r][c-1] * out_fv->data[r][c-1] + std::fabs(qx2) * _fv->Drc)/(UF_t3->data[r][c-1] + std::fabs(qx2)) : 0.0;
                UF_t3->data[r][c-1] += qx2;
                UF_t3->Drc -= qx2;
            }else if(qx2 != 0)
            {

                out_fu->Drc = ((UF_t3->Drc + std::fabs(qx2)) > UF_VERY_SMALL)? (UF_t3->Drc * out_fu->Drc + std::fabs(qx2) * _fu->data[r][c-1])/(UF_t3->Drc + std::fabs(qx2)) : 0.0;
                out_fv->Drc = ((UF_t3->Drc + std::fabs(qx2)) > UF_VERY_SMALL)? (UF_t3->Drc * out_fv->Drc + std::fabs(qx2) * _fv->data[r][c-1])/(UF_t3->Drc + std::fabs(qx2)) : 0.0;
                UF_t3->data[r][c-1] -= qx2;
                UF_t3->Drc += qx2;
            }
        }
        if(!UF_OUTORMV(_dem,r+1,c)){
            if(qy1 > 0)
            {

                out_fu->data[r+1][c] = ((UF_t3->data[r+1][c] + qy1) > UF_VERY_SMALL)? (UF_t3->data[r+1][c] * out_fu->data[r+1][c] + qy1 * _fu->Drc)/(UF_t3->data[r+1][c] + qy1) : 0.0;
                out_fv->data[r+1][c] = ((UF_t3->data[r+1][c] + qy1) > UF_VERY_SMALL)? (UF_t3->data[r+1][c] * out_fv->data[r+1][c] + qy1 * _fv->Drc)/(UF_t3->data[r+1][c] + qy1) : 0.0;
                UF_t3->data[r+1][c] += qy1;
                UF_t3->Drc -= qy1;
            }else if(qy1 != 0)
            {
                qy1 = std::fabs(qy1);
                out_fu->Drc = ((UF_t3->Drc + qy1) > UF_VERY_SMALL)? (UF_t3->Drc * out_fu->Drc + qy1 * _fu->data[r+1][c])/(UF_t3->Drc + qy1) : 0.0;
                out_fv->Drc = ((UF_t3->Drc + qy1) > UF_VERY_SMALL)? (UF_t3->Drc * out_fv->Drc + qy1 * _fv->data[r+1][c])/(UF_t3->Drc + qy1) : 0.0;
                UF_t3->data[r+1][c] -= qy1;
                UF_t3->Drc += qy1;
            }
        }
        if(!UF_OUTORMV(_dem,r-1,c)){
            if(qy2 < 0)
            {
                qy2 = std::fabs(qy2);
                out_fu->data[r-1][c] = ((UF_t3->data[r-1][c] + std::fabs(qy2)) > UF_VERY_SMALL)? (UF_t3->data[r-1][c] * out_fu->data[r-1][c] + std::fabs(qy2) * _fu->Drc)/(UF_t3->data[r-1][c] + std::fabs(qy2)) : 0.0;
                out_fv->data[r-1][c] = ((UF_t3->data[r-1][c] + std::fabs(qy2)) > UF_VERY_SMALL)? (UF_t3->data[r-1][c] * out_fv->data[r-1][c] + std::fabs(qy2) * _fv->Drc)/(UF_t3->data[r-1][c] + std::fabs(qy2)) : 0.0;
                UF_t3->data[r-1][c] += qy2;
                UF_t3->Drc -= qy2;
            }else if(qy2 != 0)
            {

                out_fu->Drc = ((UF_t3->Drc + std::fabs(qy2)) > UF_VERY_SMALL)? (UF_t3->Drc * out_fu->Drc + std::fabs(qy2) * _fu->data[r-1][c])/(UF_t3->Drc + std::fabs(qy2)) : 0.0;
                out_fv->Drc = ((UF_t3->Drc + std::fabs(qy2)) > UF_VERY_SMALL)? (UF_t3->Drc * out_fv->Drc + std::fabs(qy2) * _fv->data[r-1][c])/(UF_t3->Drc + std::fabs(qy2)) : 0.0;
                UF_t3->data[r-1][c] -= qy2;
                UF_t3->Drc += qy2;
            }
        }

    }}}

    if(UF_SOLIDPHASE)
    {
        FOR_ROW_COL_UF2D_DT
        {
            double qx1 = out_qsx1->Drc;
            double qx2 = out_qsx2->Drc;
            double qy1 = out_qsy1->Drc;
            double qy2 = out_qsy2->Drc;

            if(!UF_OUTORMV(_dem,r,c+1))
            {
                if(qx1 > 0)
                {
                    out_su->data[r][c+1] = ((UF_t4->data[r][c+1] + qx1) > UF_VERY_SMALL)? (UF_t4->data[r][c+1] * out_su->data[r][c+1] + qx1 * _su->Drc)/(UF_t4->data[r][c+1] + qx1) : 0.0;
                    out_sv->data[r][c+1] = ((UF_t4->data[r][c+1] + qx1) > UF_VERY_SMALL)? (UF_t4->data[r][c+1] * out_sv->data[r][c+1] + qx1 * _sv->Drc)/(UF_t4->data[r][c+1] + qx1) : 0.0;
                    UF_t4->data[r][c+1] += qx1;
                    UF_t4->Drc -= qx1;
                }else if(qx1 != 0)
                {
                    qx1 = std::fabs(qx1);
                    out_su->Drc = ((UF_t4->Drc + qx1) > UF_VERY_SMALL)? (UF_t4->Drc * out_su->Drc + qx1 * _su->data[r][c+1])/(UF_t4->Drc + qx1) : 0.0;
                    out_sv->Drc = ((UF_t4->Drc + qx1) > UF_VERY_SMALL)? (UF_t4->Drc * out_sv->Drc + qx1 * _sv->data[r][c+1])/(UF_t4->Drc + qx1) : 0.0;
                    UF_t4->data[r][c+1] -= qx1;
                    UF_t4->Drc += qx1;
                }
            }
            if(!UF_OUTORMV(_dem,r,c-1)){

                if(qx2 < 0)
                {
                    qx2 = std::fabs(qx2);
                    out_su->data[r][c-1] = ((UF_t4->data[r][c-1] + std::fabs(qx2)) > UF_VERY_SMALL)? (UF_t4->data[r][c-1] * out_su->data[r][c-1] + std::fabs(qx2) * _su->Drc)/(UF_t4->data[r][c-1] + std::fabs(qx2)) : 0.0;
                    out_sv->data[r][c-1] = ((UF_t4->data[r][c-1] + std::fabs(qx2)) > UF_VERY_SMALL)? (UF_t4->data[r][c-1] * out_sv->data[r][c-1] + std::fabs(qx2) * _sv->Drc)/(UF_t4->data[r][c-1] + std::fabs(qx2)) : 0.0;
                    UF_t4->data[r][c-1] += qx2;
                    UF_t4->Drc -= qx2;
                }else if(qx2 != 0)
                {
                    out_su->Drc = ((UF_t4->Drc + std::fabs(qx2)) > UF_VERY_SMALL)? (UF_t4->Drc * out_su->Drc + std::fabs(qx2) * _su->data[r][c-1])/(UF_t4->Drc + std::fabs(qx2)) : 0.0;
                    out_sv->Drc = ((UF_t4->Drc + std::fabs(qx2)) > UF_VERY_SMALL)? (UF_t4->Drc * out_sv->Drc + std::fabs(qx2) * _sv->data[r][c-1])/(UF_t4->Drc + std::fabs(qx2)) : 0.0;
                    UF_t4->data[r][c-1] -= qx2;
                    UF_t4->Drc += qx2;
                }
            }
            if(!UF_OUTORMV(_dem,r+1,c)){
                if(qy1 > 0)
                {
                    out_su->data[r+1][c] = ((UF_t4->data[r][c+1] + qy1) > UF_VERY_SMALL)? (UF_t4->data[r+1][c] * out_su->data[r+1][c] + qy1 * _su->Drc)/(UF_t4->data[r+1][c] + qy1) : 0.0;
                    out_sv->data[r+1][c] = ((UF_t4->data[r][c+1] + qy1) > UF_VERY_SMALL)? (UF_t4->data[r+1][c] * out_sv->data[r+1][c] + qy1 * _sv->Drc)/(UF_t4->data[r+1][c] + qy1) : 0.0;
                    UF_t4->data[r+1][c] += qy1;
                    UF_t4->Drc -= qy1;
                }else  if(qy1 != 0)
                {
                    qy1 = std::fabs(qy1);
                    out_su->Drc = ((UF_t4->Drc + qy1) > UF_VERY_SMALL)? (UF_t4->Drc * out_su->Drc + qy1 * _su->data[r+1][c])/(UF_t4->Drc + qy1) : 0.0;
                    out_sv->Drc = ((UF_t4->Drc + qy1) > UF_VERY_SMALL)? (UF_t4->Drc * out_sv->Drc + qy1 * _sv->data[r+1][c])/(UF_t4->Drc + qy1) : 0.0;
                    UF_t4->data[r+1][c] -= qy1;
                    UF_t4->Drc += qy1;
                }
            }
            if(!UF_OUTORMV(_dem,r-1,c)){
                if(qy2 < 0)
                {
                    qy2 = std::fabs(qy2);
                    out_su->data[r-1][c] = ((UF_t4->data[r-1][c] + std::fabs(qy2)) > UF_VERY_SMALL)? (UF_t4->data[r-1][c] * out_su->data[r-1][c] + std::fabs(qy2) * _su->Drc)/(UF_t4->data[r-1][c] + std::fabs(qy2)) : 0.0;
                    out_sv->data[r-1][c] = ((UF_t4->data[r-1][c] + std::fabs(qy2)) > UF_VERY_SMALL)? (UF_t4->data[r-1][c] * out_sv->data[r-1][c] + std::fabs(qy2) * _sv->Drc)/(UF_t4->data[r-1][c] + std::fabs(qy2)) : 0.0;
                    UF_t4->data[r-1][c] += qy2;
                    UF_t4->Drc -= qy2;
                }else  if(qy2 != 0)
                {
                    out_su->Drc = ((UF_t4->Drc + std::fabs(qy2)) > UF_VERY_SMALL)? (UF_t4->Drc * out_su->Drc + std::fabs(qy2) * _su->data[r-1][c])/(UF_t4->Drc + std::fabs(qy2)) : 0.0;
                    out_sv->Drc = ((UF_t4->Drc + std::fabs(qy2)) > UF_VERY_SMALL)? (UF_t4->Drc * out_sv->Drc + std::fabs(qy2) * _sv->data[r-1][c])/(UF_t4->Drc + std::fabs(qy2)) : 0.0;
                    UF_t4->data[r-1][c] -= qy2;
                    UF_t4->Drc += qy2;
                }
            }

        }}}
    }

}

double TWorld::UF2D_Advect2_mass(cTMap* dt, cTMap * _dem,cTMap * _m, cTMap * f,cTMap * _qx1, cTMap * _qx2,cTMap * _qy1, cTMap * _qy2, cTMap * out_m)
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
        if(f->Drc < UF_VERY_SMALL)
        {
            continue;
        }


        double conc =(f->Drc > UF_VERY_SMALL)? (_m->Drc/f->Drc) :0.0;

        if(conc > 0)
        {
            double qx1 = std::min(_qx1->Drc,0.25 * _m->Drc);
            double qx2 = std::min(_qx2->Drc,0.25 * _m->Drc);
            double qy1 = std::min(_qy1->Drc,0.25 * _m->Drc);
            double qy2 = std::min(_qy2->Drc,0.25 * _m->Drc);

            if(!UF_OUTORMV(_dem,r,c+1))
            {
                if(qx1 > 0)
                {
                    out_m->data[r][c+1] += conc * qx1;
                    out_m->Drc -= conc *qx1;
                }else
                {
                    /*qx1 = std::fabs(qx1);
                    double c2 = (f->data[r][c+1] > UF_VERY_SMALL)? _m->data[r][c+1]/f->data[r][c+1] : 0.0;
                    out_m->data[r][c+1] -= c2 * qx1;
                    out_m->Drc += c2 * qx1;*/
                }
            }else
            {
                if(qx1 > 0)
                {
                    outflow += conc*qx1;
                    out_m->Drc -= conc*qx1;
                }
            }
            if(!UF_OUTORMV(_dem,r,c-1)){

                if(qx2 < 0)
                {
                    qx2 = std::fabs(qx2);
                    out_m->data[r][c-1] += conc * qx2;
                    out_m->Drc -= conc * qx2;
                }else
                {
                    /*double c2 = (f->data[r][c-1] > UF_VERY_SMALL)? _m->data[r][c-1]/f->data[r][c-1] : 0.0;
                    out_m->data[r][c-1] -= c2 * qx2;
                    out_m->Drc += c2 * qx2;*/
                }
            }else
            {
                if(qx2 < 0)
                {
                    qx2 = std::fabs(qx2);
                    outflow += conc*qx2;
                    out_m->Drc -= conc*qx2;
                }
            }
            if(!UF_OUTORMV(_dem,r+1,c)){
                if(qy1 > 0)
                {
                    out_m->data[r+1][c] += conc * qy1;
                    out_m->Drc -= conc * qy1;
                }else
                {
                    /*qy1 = std::fabs(qy1);
                    double c2 = (f->data[r+1][c] > UF_VERY_SMALL)? _m->data[r+1][c]/f->data[r+1][c] : 0.0;
                    out_m->data[r+1][c] -= c2 * qy1;
                    out_m->Drc += c2 * qy1;*/
                }
            }else
            {
                if(qy1 > 0)
                {
                    outflow += conc*qy1;
                    out_m->Drc -= conc*qy1;
                }
            }
            if(!UF_OUTORMV(_dem,r-1,c)){
                if(qy2 < 0)
                {
                    qy2 = std::fabs(qy2);
                    out_m->data[r-1][c] += conc * qy2;
                    out_m->Drc -= conc * qy2;
                }else
                {
                    /*double c2 = (f->data[r-1][c] > UF_VERY_SMALL)?_m->data[r-1][c]/f->data[r-1][c] : 0.0;
                    out_m->data[r-1][c] -= c2 * qy2;
                    out_m->Drc += c2 * qy2;*/
                }
            }else
            {
                if(qy2 < 0)
                {
                    qy2 = std::fabs(qy2);
                    outflow += conc*qy2;
                    out_m->Drc -= conc*qy2;
                }
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
    return outflow;

}

void TWorld::UF2D_Advect2_prop(cTMap* dt, cTMap * _dem,cTMap * _m, cTMap * f,cTMap * _qx1, cTMap * _qx2,cTMap * _qy1, cTMap * _qy2,cTMap *_prop, cTMap * out_prop)
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
        if(f->Drc < UF_VERY_SMALL)
        {
            continue;
        }
        double qx1 = _qx1->Drc;
        double qx2 = _qx2->Drc;
        double qy1 = _qy1->Drc;
        double qy2 = _qy2->Drc;

        double conc = (f->Drc > UF_VERY_SMALL)? (_m->Drc/f->Drc) : 0.0;
        if(!UF_OUTORMV(_dem,r,c+1))
        {
            if(qx1 > 0)
            {
                UF_t1->data[r][c+1] += conc * qx1;
                UF_t1->Drc -= conc *qx1;
                out_prop->data[r][c+1] = (UF_t1->data[r][c+1] +  conc*qx1 > UF_VERY_SMALL)?(out_prop->data[r][c+1]  * UF_t1->data[r][c+1] + _prop->Drc * conc*qx1)/(UF_t1->data[r][c+1] +  conc*qx1) : 0.0;
            }else
            {
                qx1 = std::fabs(qx1);
                double c2 = (f->data[r][c+1] > UF_VERY_SMALL)? _m->data[r][c+1]/f->data[r][c+1] : 0.0;
                UF_t1->data[r][c+1] -= c2 * qx1;
                UF_t1->Drc += c2 * qx1;
                out_prop->Drc = (UF_t1->Drc +  c2 * qx1 > UF_VERY_SMALL)? (out_prop->Drc  * UF_t1->Drc + _prop->data[r][c+1] * c2 * qx1)/(UF_t1->Drc +  c2 * qx1) : 0.0;
            }
        }
        if(!UF_OUTORMV(_dem,r,c-1)){

            if(qx2 < 0)
            {
                qx2 = std::fabs(qx2);
                UF_t1->data[r][c-1] += conc * qx2;
                UF_t1->Drc -= conc * qx2;
                out_prop->data[r][c-1] = (UF_t1->data[r][c-1] +  conc*qx2 > UF_VERY_SMALL)?(out_prop->data[r][c-1]  * UF_t1->data[r][c-1] + _prop->Drc * conc*qx2)/(UF_t1->data[r][c-1] +  conc*qx2) : 0.0;
            }else
            {
                double c2 = (f->data[r][c-1] > UF_VERY_SMALL)? _m->data[r][c-1]/f->data[r][c-1] : 0.0;
                UF_t1->data[r][c-1] -= c2 * qx2;
                UF_t1->Drc += c2 * qx2;
                out_prop->Drc = (UF_t1->Drc +  c2 * qx2 > UF_VERY_SMALL)? (out_prop->Drc  * UF_t1->Drc + _prop->data[r][c-1] * c2 * qx2)/(UF_t1->Drc +  c2 * qx2) : 0.0;
            }
        }
        if(!UF_OUTORMV(_dem,r+1,c)){
            if(qy1 > 0)
            {
                UF_t1->data[r+1][c] += conc * qy1;
                UF_t1->Drc -= conc * qy1;
                out_prop->data[r+1][c] = (UF_t1->data[r+1][c] +  conc*qy1 > UF_VERY_SMALL)?(out_prop->data[r+1][c]  * UF_t1->data[r+1][c] + _prop->Drc * conc*qy1)/(UF_t1->data[r+1][c] +  conc*qy1) : 0.0;
            }else
            {
                qy1 = std::fabs(qy1);
                double c2 = (f->data[r+1][c] > UF_VERY_SMALL)? _m->data[r+1][c]/f->data[r+1][c] : 0.0;
                UF_t1->data[r+1][c] -= c2 * qy1;
                UF_t1->Drc += c2 * qy1;
                out_prop->Drc = (UF_t1->Drc +  c2 * qy1 > UF_VERY_SMALL)? (out_prop->Drc  * UF_t1->Drc + _prop->data[r+1][c] * c2 * qy1)/(UF_t1->Drc +  c2 * qy1) : 0.0;
            }
        }
        if(!UF_OUTORMV(_dem,r-1,c)){
            if(qy2 < 0)
            {
                qy2 = std::fabs(qy2);
                UF_t1->data[r-1][c] += conc * qy2;
                UF_t1->Drc -= conc * qy2;
                out_prop->data[r-1][c] = (UF_t1->data[r-1][c] +  conc*qy2 > UF_VERY_SMALL)?(out_prop->data[r-1][c]  * UF_t1->data[r-1][c] + _prop->Drc * conc*qy2)/(UF_t1->data[r-1][c] +  conc*qy2) : 0.0;
            }else
            {
                double c2 = (f->data[r-1][c] > UF_VERY_SMALL)? _m->data[r-1][c]/f->data[r-1][c] : 0.0;
                UF_t1->data[r-1][c] -= c2 * qy2;
                UF_t1->Drc += c2 * qy2;
                out_prop->Drc = (UF_t1->Drc +  c2 * qy2 > UF_VERY_SMALL)? (out_prop->Drc  * UF_t1->Drc + _prop->data[r-1][c] * c2 * qy2)/(UF_t1->Drc +  c2 * qy2) : 0.0;
            }
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


void TWorld::UF1D_Advect2_Momentum(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_fu,  cTMap * out_su, cTMap * out_qf1, cTMap * out_qf2,cTMap * out_qs1, cTMap * out_qs2)
{

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    FOR_ROW_COL_UF1D
    {
        out_fu->Drc = _fu->Drc;
        out_su->Drc = _su->Drc;
        UF_t3->Drc = _f->Drc;
        UF_t4->Drc = _f->Drc;
        UF1D_fq1->Drc = UF1D_fq1->Drc;
        UF1D_fq2->Drc = UF1D_fq2->Drc;
        UF1D_sq1->Drc = UF1D_sq1->Drc;
        UF1D_sq2->Drc = UF1D_sq2->Drc;
    }

    FOR_ROW_COL_UF1D_DT
    {
        double q1 = UF1D_fq1->Drc;
        double q2 = UF1D_fq2->Drc;

        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!(lddself == 5))
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){

                if( q1 >0)
                {
                    out_fu->data[r2][c2] = ((UF_t3->data[r2][c2] + q1) > 0)? (UF_t3->data[r2][c2] * out_fu->data[r2][c2] + q1 * _fu->Drc)/(UF_t3->data[r2][c2] + q1) : 0.0;
                }else
                {
                    out_fu->Drc = ((UF_t3->Drc + std::fabs(q1)) > 0)? (UF_t3->Drc * out_fu->Drc + std::fabs(q1) * _fu->data[r2][c2])/(UF_t3->Drc + std::fabs(q1)) : 0.0;
                }

                UF_t3->Drc -= q1;
                UF_t3->data[r2][c2] += q1;
            }else if(q1 > 0)
            {
                UF_t3->Drc -= std::fabs(q1);
            }
        }else if(q1 > 0)
        {
            UF_t3->Drc -= std::fabs(q1);
        }

        //backward flux
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
        if(totalwidth > 0)
        {
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

                        if( q2 < 0)
                        {
                            out_fu->data[r2][c2] = ((UF_t3->data[r2][c2] + std::fabs(q2)) > 0)? (UF_t3->data[r2][c2] * out_fu->data[r2][c2] + std::fabs(q2) * _fu->Drc)/(UF_t3->data[r2][c2] + std::fabs(q2)) : 0.0;
                        }else
                        {
                            out_fu->Drc = ((UF_t3->Drc + q2) > 0)? (UF_t3->Drc * out_fu->Drc + q2 * _fu->data[r2][c2])/(UF_t3->Drc + q2) : 0.0;
                        }

                        UF_t3->Drc += q2;
                        UF_t3->data[r2][c2] -= q2;
                    }else if(q2 < 0)
                    {
                        UF_t3->Drc -= std::fabs(q2);
                    }
                }
            }
        }else if(q2 < 0)
        {
            UF_t3->Drc -= std::fabs(q2);
        }
    }

    FOR_ROW_COL_UF1D_DT
    {
        double q1 = UF1D_sq1->Drc;
        double q2 = UF1D_sq2->Drc;

        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!(lddself == 5))
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){
                if( q1 >0)
                {
                    out_su->data[r2][c2] = ((UF_t4->data[r2][c2] + q1) > 0)? (UF_t4->data[r2][c2] * out_su->data[r2][c2] + q1 * _su->Drc)/(UF_t4->data[r2][c2] + q1) : 0.0;
                }else
                {
                    out_su->Drc = ((UF_t4->Drc + std::fabs(q1)) > 0)? (UF_t4->Drc * out_su->Drc + std::fabs(q1) * _su->data[r2][c2])/(UF_t4->Drc + std::fabs(q1)) : 0.0;
                }

                UF_t4->Drc -= q1;
                UF_t4->data[r2][c2] += q1;
            }else if(q1 > 0)
            {
                UF_t4->Drc -= std::fabs(q1);
            }
        }else if(q1 > 0)
        {
            UF_t4->Drc -= std::fabs(q1);
        }
        //backward flux

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
        if(totalwidth > 0)
        {
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

                        if( q2 < 0)
                        {
                            out_su->data[r2][c2] = ((UF_t4->data[r2][c2] + std::fabs(q2)) > 0)? (UF_t4->data[r2][c2] * out_su->data[r2][c2] + std::fabs(q2) * _su->Drc)/(UF_t4->data[r2][c2] + std::fabs(q2)) : 0.0;
                        }else
                        {
                            out_su->Drc = ((UF_t4->Drc + q2) > 0)? (UF_t4->Drc * out_su->Drc + q2 * _su->data[r2][c2])/(UF_t4->Drc + q2) : 0.0;
                        }

                        UF_t4->Drc += q2;
                        UF_t4->data[r2][c2] -= q2;
                    }else if(q2 < 0)
                    {
                        UF_t4->Drc -= std::fabs(q2);
                    }
                }
            }
        }else if(q2 < 0)
        {
            UF_t4->Drc -= std::fabs(q2);
        }
    }
}

double TWorld::UF1D_Advect2_mass(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _m, cTMap * _f, cTMap * _q1,cTMap * _q2, cTMap * out_m)
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
        double conc = (_f->Drc > 0)? _m->Drc/_f->Drc: 0.0;

        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!(lddself == 5))
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){

                double con2 = _q1->Drc > 0? conc : (_f->data[r2][c2]>0?(_m->data[r2][c2]/_f->data[r2][c2]):0.0);
                out_m->Drc -= conc *_q1->Drc;
                out_m->data[r2][c2] += conc *_q1->Drc;

            }else if(_q1->Drc > 0)
            {
                outflow += conc *_q1->Drc;
                out_m->Drc -= conc *_q1->Drc;
            }
        }else if(_q1->Drc > 0)
        {
                outflow += conc *_q1->Drc;
                out_m->Drc -= conc *_q1->Drc;
        }

        //backward flux
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
        if(totalwidth)
        {
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
                        double q2 = (_lddw->data[r2][c2]/totalwidth)*_q2->Drc;

                        double con2 = q2 < 0? conc : (_f->data[r2][c2]>0?(_m->data[r2][c2]/_f->data[r2][c2]):0.0);
                        out_m->Drc += con2 *q2;
                        out_m->data[r2][c2] -= con2 *q2;

                    }else if(_q2->Drc < 0)
                    {
                        outflow += conc * std::fabs(_q2->Drc);
                        out_m->Drc -= conc *std::fabs(_q2->Drc);
                    }
                }
            }
        }else if(_q2->Drc < 0)
        {
            outflow += conc * std::fabs(_q2->Drc);
            out_m->Drc -= conc * std::fabs(_q2->Drc);
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

void TWorld::UF1D_Advect2_prop(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _m, cTMap * _f, cTMap * _q1,cTMap * _q2,cTMap *_prop, cTMap * out_prop)
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
        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!lddself == 5)
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){

                if(_q1->Drc > 0)
                {
                    double flux = _q1->Drc;
                    out_prop->data[r2][c2] = (UF_t1->data[r2][c2] +  flux) > 0? (out_prop->data[r2][c2]  * UF_t1->data[r2][c2] + _prop->Drc * flux)/(UF_t1->data[r2][c2] +  flux) : 0.0;
                }else
                {
                    double flux = std::fabs(_q1->Drc);
                    out_prop->Drc = (UF_t1->Drc +  flux) > 0? (out_prop->Drc  * UF_t1->Drc + _prop->data[r2][c2] * flux)/(UF_t1->Drc +  flux) : 0.0;
                }
                UF_t1->Drc -= _q1->Drc;
                UF_t1->data[r2][c2] += _q1->Drc;

            }
        }else if( _q1->Drc > 0)
        {
            UF_t1->Drc -= _q1->Drc;
        }

        //backward flux
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
        if(totalwidth > 0)
        {
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
                        double q2 = (_lddw->data[r2][c2]/totalwidth)*_q2->Drc;

                        if(_q1->Drc > 0)
                        {
                            double flux = std::fabs(_q1->Drc);
                            out_prop->Drc = (UF_t1->Drc +  flux) > 0? (out_prop->Drc  * UF_t1->Drc + _prop->data[r2][c2] * flux)/(UF_t1->Drc +  flux) : 0.0;
                        }else
                        {
                            double flux = _q1->Drc;
                            out_prop->data[r2][c2] = (UF_t1->data[r2][c2] +  flux) > 0? (out_prop->data[r2][c2]  * UF_t1->data[r2][c2] + _prop->Drc * flux)/(UF_t1->data[r2][c2] +  flux) :0.0;
                        }

                        UF_t1->Drc += q2;
                        UF_t1->data[r2][c2] -= q2;
                    }
                }
            }
        }else if( _q2->Drc < 0)
        {
            UF_t1->Drc -= std::fabs(_q2->Drc);
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
