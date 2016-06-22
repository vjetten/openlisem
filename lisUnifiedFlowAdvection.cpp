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

void TWorld::UF2D_Advect_Momentum(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_fu, cTMap * out_fv,  cTMap * out_su, cTMap * out_sv )
{
    FOR_ROW_COL_UF2D
    {
        out_fu->Drc = _fu->Drc;
        out_fv->Drc = _fv->Drc;
        out_su->Drc = _su->Drc;
        out_sv->Drc = _sv->Drc;
    }

    FOR_ROW_COL_UF2D_DT
    {
        double hf = _f->Drc/(_dx*_dx);
        double hs = _s->Drc/(_dx*_dx);

        double qfx = std::min(std::fabs(dt->Drc * hf * _fu->Drc*_dx),0.5 * UF_Courant * _f->Drc);
        double qfy = std::min(std::fabs(dt->Drc * hf * _fv->Drc*_dx),0.5 * UF_Courant * _f->Drc);
        double qsx = std::min(std::fabs(dt->Drc * hs * _su->Drc*_dx),0.5 * UF_Courant * _s->Drc);
        double qsy = std::min(std::fabs(dt->Drc * hs * _sv->Drc*_dx),0.5 * UF_Courant * _s->Drc);

        int qfxc = c + ((_fu->Drc > 0)? -1 : 1);
        int qfyr = r + ((_fv->Drc > 0)? -1 : 1);
        int qsxc = c + ((_su->Drc > 0)? -1 : 1);
        int qsyr = r + ((_sv->Drc > 0)? -1 : 1);

        if(!UF_OUTORMV(_dem,r,qfxc)){
            out_fu->data[r][qfxc] = ((_f->data[r][qfxc] + qfx) > 0)? (_f->data[r][qfxc] * out_fu->data[r][qfxc] + qfx * _fu->Drc)/(_f->data[r][qfxc] + qfx) : 0.0;
            out_fv->data[r][qfxc] = ((_f->data[r][qfxc] + qfx) > 0)? (_f->data[r][qfxc] * out_fv->data[r][qfxc] + qfx * _fv->Drc)/(_f->data[r][qfxc] + qfx) : 0.0;
        }
        if(!UF_OUTORMV(_dem,qfyr,c)){
            out_fu->data[qfyr][c] = ((_f->data[qfyr][c] + qfy) > 0)? (_f->data[qfyr][c] * out_fu->data[qfyr][c] + qfy * _fu->Drc)/(_f->data[qfyr][c] + qfy) : 0.0;
            out_fv->data[qfyr][c] = ((_f->data[qfyr][c] + qfy) > 0)? (_f->data[qfyr][c] * out_fv->data[qfyr][c] + qfy * _fv->Drc)/(_f->data[qfyr][c] + qfy) : 0.0;
        }
        if(!UF_OUTORMV(_dem,r,qsxc)){
            out_su->data[r][qsxc] = ((_s->data[r][qsxc] + qsx) > 0)? (_s->data[r][qsxc] * out_su->data[r][qsxc] + qsx * _su->Drc)/(_s->data[r][qsxc] + qsx) : 0.0;
            out_sv->data[r][qsxc] = ((_s->data[r][qsxc] + qsx) > 0)? (_s->data[r][qsxc] * out_sv->data[r][qsxc] + qsx * _sv->Drc)/(_s->data[r][qsxc] + qsx) : 0.0;
        }
        if(!UF_OUTORMV(_dem,qsyr,c)){
            out_su->data[qsyr][c] = ((_s->data[qsyr][c] + qsy) > 0)? (_s->data[qsyr][c] * out_su->data[qsyr][c] + qsy * _su->Drc)/(_s->data[qsyr][c] + qsy) : 0.0;
            out_sv->data[qsyr][c] = ((_s->data[qsyr][c] + qsy) > 0)? (_s->data[qsyr][c] * out_sv->data[qsyr][c] + qsy * _sv->Drc)/(_s->data[qsyr][c] + qsy) : 0.0;
        }

    }}}

}


double TWorld::UF2D_Advect_mass(cTMap * dt, cTMap * _dem,cTMap * _m, cTMap * _mu, cTMap * _mv, cTMap * out_m)
{
    double outflow = 0;
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
        double hm = _m->Drc/(_dx*_dx);

        double qmx = std::min(std::fabs(dt->Drc * hm * _mu->Drc*_dx), 0.5 * UF_Courant * UF_t1->Drc);
        double qmy = std::min(std::fabs(dt->Drc * hm * _mv->Drc*_dx), 0.5 * UF_Courant * UF_t1->Drc);

        int qmxc = c + ((_mu->Drc > 0)? 1 : -1);
        int qmyr = r + ((_mv->Drc > 0)? 1 : -1);


        if(!UF_OUTORMV(_dem,r,qmxc)){
            out_prop->data[r][qmxc] = (out_prop->data[r][qmxc]  * UF_t1->data[r][qmxc] + _prop->Drc * qmx)/(UF_t1->data[r][qmxc] +  qmx);
        }
        if(!UF_OUTORMV(_dem,qmyr,c)){
            out_prop->data[qmyr][c] = (out_prop->data[qmyr][c]  * UF_t1->data[qmyr][c] + _prop->Drc * qmy)/(UF_t1->data[qmyr][c] +  qmy);
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




void TWorld::UF2D_FluidMomentumSource(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv)
{
    FOR_ROW_COL_UF2D
    {

        out_fu->Drc = _fu->Drc;
        out_fv->Drc = _fv->Drc;
        double h = (_f->Drc + _s->Drc)/(_dx * _dx);
        UF_t1->Drc = h*h*(-UF_Gravity*h)/2.0;

        if(_f->Drc + _s->Drc > UF_VERY_SMALL)
        {
            UF_t2->Drc = _f->Drc/(_f->Drc + _s->Drc);
            UF_t3->Drc = _s->Drc/(_f->Drc + _s->Drc);
        }else
        {
            UF_t2->Drc = 0;
            UF_t3->Drc = 0;
        }

    }

    FOR_ROW_COL_UF2D_DT
    {
        if((_f->Drc < UF_VERY_SMALL))
        {
            out_fu->Drc = 0;
            out_fv->Drc = 0;
            continue;
        }

        double h = (_f->Drc + _s->Drc)/(_dx * _dx);
        double ff = UF_t2->Drc;
        double sf = UF_t3->Drc;
        double Nr = 1500;//UF_Reynolds(1000.0,_visc->Drc,ff);
        double Nra = 15000;//UF_Reynolds(1000.0,_visc->Drc,ff);
        double gamma = _d->Drc == 0? 0.5 : 1000.0/_d->Drc;
        double dc = 0;//UF_DragCoefficient(ff,sf,UF_TerminalVelocity(UF2D_rocksize->Drc,ff,_visc->Drc,sf,_d->Drc),gamma);
        double pbf = -UF_Gravity*h;

        double dh2pbdx = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_X);
        double dh2pbdy = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_Y);

        double dsfdx = UF2D_Derivative(_dem,UF_t3,r,c,UF_DIRECTION_X);
        double dsfdy = UF2D_Derivative(_dem,UF_t3,r,c,UF_DIRECTION_Y);
        double ddsfdxx = UF2D_Derivative2(_dem,UF_t3,r,c,UF_DIRECTION_X);
        double ddsfdyy = UF2D_Derivative2(_dem,UF_t3,r,c,UF_DIRECTION_Y);
        double ddsfdxy = UF2D_Derivative2(_dem,UF_t3,r,c,UF_DIRECTION_XY);

        double dfudx = UF2D_Derivative(_dem,_fu,r,c,UF_DIRECTION_X);
        double dfudy = UF2D_Derivative(_dem,_fu,r,c,UF_DIRECTION_Y);
        double dfvdx = UF2D_Derivative(_dem,_fv,r,c,UF_DIRECTION_X);
        double dfvdy = UF2D_Derivative(_dem,_fv,r,c,UF_DIRECTION_Y);
        double ddfudxx = UF2D_Derivative2(_dem,_fu,r,c,UF_DIRECTION_X);
        double ddfudyy = UF2D_Derivative2(_dem,_fu,r,c,UF_DIRECTION_X);
        double ddfvdxy = UF2D_Derivative2(_dem,_fv,r,c,UF_DIRECTION_XY);

        double dsudx = UF2D_Derivative(_dem,_su,r,c,UF_DIRECTION_X);
        double dsudy = UF2D_Derivative(_dem,_su,r,c,UF_DIRECTION_Y);
        double dsvdx = UF2D_Derivative(_dem,_sv,r,c,UF_DIRECTION_X);
        double dsvdy = UF2D_Derivative(_dem,_sv,r,c,UF_DIRECTION_Y);

        if(ff > UF_VERY_SMALL)
        {
            out_fu->Drc += dt->Drc *
                (UF_Gravity * sin(UF2D_SlopeX->Drc) -
                 UF_Aspect *(
                     (dh2pbdx)/h
                     +(pbf * UF2D_SlopeX->Drc)
                     -(1.0/(ff * Nr))*(
                         2.0*ddfudxx + ddfvdxy + ddfudyy + UF_Chi * _fu->Drc/(UF_Aspect*UF_Aspect * h* h)
                         )
                     +(1.0/(ff * Nra))*(
                         2.0 *(dsfdx*(dfudx - dsudx) + ddsfdxx*(_fu->Drc - _su->Drc))
                         +(dsfdx*(dfvdy - dsvdy) + ddsfdxy*(_fv->Drc - _sv->Drc))
                         +(dsfdy*(dfudy - dsudy) + ddsfdyy*(_fu->Drc - _su->Drc))
                         + UF_Chi * _fu->Drc/(UF_Aspect*UF_Aspect*ff*Nra*h*h)
                         )
                     -(UF_Ksi*sf*(_fu->Drc - _su->Drc)/(UF_Aspect*UF_Aspect*ff*Nra*h*h))
                     )
                 );

            out_fv->Drc += dt->Drc *
                    (UF_Gravity * sin(UF2D_SlopeY->Drc) -
                     UF_Aspect * (
                         (dh2pbdy)/h
                         +(pbf * UF2D_SlopeY->Drc)
                         -(1.0/(ff * Nr))*(
                             2.0*ddfudyy + ddfvdxy + ddfudxx + UF_Chi * _fv->Drc/(UF_Aspect*UF_Aspect * h* h)
                             )
                         +(1.0/(ff * Nra))*(
                             2.0 *(dsfdy*(dfvdy - dsvdy) + ddsfdyy*(_fv->Drc - _sv->Drc))
                             +(dsfdy*(dfudx - dsudx) + ddsfdxy*(_fu->Drc - _su->Drc))
                             +(dsfdx*(dfvdx - dsvdx) + ddsfdxx*(_fv->Drc - _sv->Drc))
                             + UF_Chi * _fv->Drc/(UF_Aspect*UF_Aspect*ff*Nra*h*h)
                             )
                         -(UF_Ksi*sf*(_fv->Drc - _sv->Drc)/(UF_Aspect*UF_Aspect*ff*Nra*h*h))
                         )
                     );

            if(sf > UF_VERY_SMALL)
            {
                double vpow = pow((_fu->Drc-_su->Drc)*(_fu->Drc-_su->Drc)+(_fv->Drc-_sv->Drc)*(_fv->Drc-_sv->Drc),0.5*(UF_j-1.0));
                double facu = std::min(1.0,std::max(0.0,std::fabs(dt->Drc * dc * sf *(_fu->Drc - _su->Drc) *vpow)));
                double facv = std::min(1.0,std::max(0.0,std::fabs(dt->Drc * dc * sf *(_fv->Drc - _sv->Drc) *vpow)));
                out_fu->Drc = (1.0-facu) * out_fu->Drc + (facu) * _su->Drc;
                out_fv->Drc = (1.0-facv) * out_fv->Drc + (facv) * _sv->Drc;
            }
        }


        double nsq = (0.1 + N->Drc)*(0.1 + N->Drc)*UF_Gravity*sqrt(_fu->Drc*_fu->Drc+_fv->Drc*_fv->Drc)*dt->Drc/pow(h,4.0/3.0) ;

        double fac = h > 0? 1.0/std::min(h* 3.0,1.0) : 0.0;
        out_fu->Drc = (out_fu->Drc/(1.0+nsq * fac));
        out_fv->Drc = (out_fv->Drc/(1.0+nsq * fac));
    }}}


}

void TWorld::UF2D_SolidMomentumSource(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv)
{
    FOR_ROW_COL_UF2D
    {
        out_su->Drc = _su->Drc;
        out_sv->Drc = _sv->Drc;
        double h = (_f->Drc + _s->Drc)/(_dx * _dx);
        double hs = (_s->Drc)/(_dx * _dx);
        UF_t1->Drc = h;

        double ifa = 0.3;

        double ifasin = sin(ifa);
        double ifatan = tan(ifa);

        double dudxx = UF2D_Derivative2(_dem,_su,r,c,UF_DIRECTION_X);
        double dudyy = UF2D_Derivative2(_dem,_sv,r,c,UF_DIRECTION_Y);

        double k_act = (1- ifasin)/(1+ifasin);
        double k_pass = (1+ifasin)/(1-ifasin);
        double k_x = (dudxx) > 0? k_act: k_pass;
        double k_y = (dudyy) > 0? k_act: k_pass;

        UF_t2->Drc = UF_Gravity * k_x *hs * hs * UF_Aspect/2.0;
        UF_t3->Drc = UF_Gravity * k_x *hs * hs * UF_Aspect/2.0;
    }

    FOR_ROW_COL_UF2D_DT

        if(!(_s->Drc > UF_VERY_SMALL))
        {
            out_su->Drc = 0;
            out_sv->Drc = 0;
            continue;
        }

        double h = (_f->Drc + _s->Drc)/(_dx * _dx);
        double ff = _f->Drc / (_f->Drc + _s->Drc);
        double sf = _s->Drc / (_f->Drc + _s->Drc);
        double gamma = 0.5;//1000.0/_d->Drc;
        double dc = 0.5;//UF_DragCoefficient(ff,sf,gamma,_visc->Drc, _rocksize->Drc, _d->Drc);
        double pbf = -UF_Gravity*h;
        double pbs = (1-gamma)*pbf;
        double ifa = 0.3;

        double dhdx = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_X);
        double dhdy = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_Y);
        double dbdx = UF2D_Derivative(_dem,UF_t2,r,c,UF_DIRECTION_X);
        double dbdy = UF2D_Derivative(_dem,UF_t3,r,c,UF_DIRECTION_Y);


        if(sf > UF_VERY_SMALL)
        {
            out_su->Drc += dt->Drc *
                (UF_Gravity * sin(UF2D_SlopeX->Drc) - (_su->Drc > 0? 1.0 : -1.0)*std::tan(ifa)*pbs-UF_Aspect*pbs*UF2D_SlopeX->Drc
                -UF_Aspect * gamma * pbf * ( dhdx +  UF2D_SlopeX->Drc ) -
                - dbdy
                 );

            out_sv->Drc +=  dt->Drc *
                (UF_Gravity * sin(UF2D_SlopeY->Drc) - (_su->Drc > 0? 1.0 : -1.0)*std::tan(ifa)*pbs-UF_Aspect*pbs*UF2D_SlopeY->Drc
                 -UF_Aspect * gamma * pbf * ( dhdy +  UF2D_SlopeY->Drc )
                 - dbdx
                );

            if(ff > UF_VERY_SMALL)
            {
                double vpow = pow((_fu->Drc-_su->Drc)*(_fu->Drc-_su->Drc)+(_fv->Drc-_sv->Drc)*(_fv->Drc-_sv->Drc),0.5*(UF_j-1.0));
                double facu = std::min(1.0,std::max(-1.0,std::fabs(dt->Drc * dc * ff *(_su->Drc - _fu->Drc) *vpow)));
                double facv = std::min(1.0,std::max(-1.0,std::fabs(dt->Drc * dc * ff *(_sv->Drc - _fv->Drc) *vpow)));
                out_su->Drc = (1.0-facu) * out_su->Drc + (facu) * _fu->Drc;
                out_sv->Drc = (1.0-facv) * out_sv->Drc + (facv) * _fv->Drc;
            }
        }
    }}

}

void TWorld::UF2D_Diffuse_mass(cTMap* dt, cTMap * _dem,cTMap * _m,cTMap * _f, cTMap * _s, cTMap * _fu, cTMap * _fv, cTMap * _su, cTMap * _sv, cTMap * out_m)
{
    FOR_ROW_COL_UF2D
    {
        out_m->Drc = _m->Drc;
    }

    FOR_ROW_COL_UF2D_DT
    {
        /*double hf = ;
        double hs = ;
        double mc = _m->Drc/(hf*_dx*_dx);

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
        }*/
    }}}
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
    }

    FOR_ROW_COL_UF1D_DT
    {
        double hf = _f->Drc/(_dx*_lddw->Drc);
        double hs = _s->Drc/(_dx*_lddw->Drc);

        double qfx = std::min(std::fabs(dt->Drc * hf * _fu->Drc*_lddw->Drc), 0.5 * UF_Courant * _f->Drc);
        double qsx = std::min(std::fabs(dt->Drc * hs * _su->Drc*_lddw->Drc), 0.5 * UF_Courant * _s->Drc);

        //frontward flux
        if(_fu->Drc < 0)
        {
            //front cell
            int lddself = (int) _ldd->data[r][c];
            if(!lddself == 5)
            {
                int r2 = r+dy[lddself];
                int c2 = c+dx[lddself];

                if(!UF_OUTORMV(_ldd,r2,c2)){
                    out_fu->data[r2][c2] = ((_f->data[r2][c2] + qfx) > 0)? (_f->data[r2][c2] * out_fu->data[r2][c2] + qfx * _fu->Drc)/(_f->data[r2][c2] + qfx) : 0.0;
                }
            }
        }else
        //backward flux
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

                        out_fu->data[r2][c2] = ((_f->data[r2][c2] + flux) > 0)? (_f->data[r2][c2] * out_fu->data[r2][c2] + flux * _fu->Drc)/(_f->data[r2][c2] + flux) : 0.0;
                    }
                }
            }
        }
        //frontward flux
        if(_su->Drc < 0)
        {
            //front cell
            int lddself = (int) _ldd->data[r][c];
            if(!lddself == 5)
            {
                int r2 = r+dy[lddself];
                int c2 = c+dx[lddself];

                if(!UF_OUTORMV(_ldd,r2,c2)){
                    out_su->data[r2][c2] = ((_s->data[r2][c2] + qsx) > 0)? (_s->data[r2][c2] * out_su->data[r2][c2] + qsx * _su->Drc)/(_s->data[r2][c2] + qsx) : 0.0;
                }
            }
        }else
        //backward flux
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

                        out_su->data[r2][c2] = ((_s->data[r2][c2] + flux) > 0)? (_s->data[r2][c2] * out_su->data[r2][c2] + flux * _su->Drc)/(_s->data[r2][c2] + flux) : 0.0;
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
    FOR_ROW_COL_UF1D
    {
        out_m->Drc = _m->Drc;
    }

    FOR_ROW_COL_UF1D_DT
    {
        double hm = _m->Drc/(_dx*_lddw->Drc);

        double qmx = std::min(std::fabs(dt->Drc * hm * _mu->Drc*_lddw->Drc), 0.5 * UF_Courant * out_m->Drc);

        //frontward flux
        if(_mu->Drc < 0)
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
        //backward flux
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
                    }
                }
            }
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
        if(_mu->Drc < 0)
        {
            //front cell
            int lddself = (int) _ldd->data[r][c];
            if(!lddself == 5)
            {
                int r2 = r+dy[lddself];
                int c2 = c+dx[lddself];

                if(!UF_OUTORMV(_ldd,r2,c2)){
                    out_prop->data[r2][c2] = (out_prop->data[r2][c2]  * UF_t1->data[r2][c2] + _prop->Drc * qmx)/(UF_t1->data[r2][c2] +  qmx);

                    UF_t1->data[r2][c2] += qmx;
                    UF_t1->data[r][c] -= qmx;
                }
            }
        }else
        //backward flux
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
                        out_prop->data[r2][c2] = (out_prop->data[r2][c2]  * UF_t1->data[r2][c2] + _prop->Drc * flux)/(UF_t1->data[r2][c2] +  flux);

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


void TWorld::UF1D_SolidMomentumSource(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_su)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    FOR_ROW_COL_UF1D
    {
        out_su->Drc = _su->Drc;
        double h = (_f->Drc + _s->Drc)/(_dx * _lddw->Drc);
        double hs = (_s->Drc)/(_dx * _lddw->Drc);
        UF_t1->Drc = h;

        double ifa = 0.3;

        double ifasin = sin(ifa);
        double ifatan = tan(ifa);

        double dudxx = UF1D_Derivative2(_ldd,_lddw,_su,r,c);

        double k_act = (1- ifasin)/(1+ifasin);
        double k_pass = (1+ifasin)/(1-ifasin);
        double k_x = (dudxx) > 0? k_act: k_pass;

        UF_t2->Drc = UF_Gravity * k_x *hs * hs * UF_Aspect/2.0;
    }

    FOR_ROW_COL_UF1D_DT
    {
        if(!(_s->Drc > UF_VERY_SMALL))
        {
            out_su->Drc = 0;
            continue;
        }

        double h = (_f->Drc + _s->Drc)/(_dx * _lddw->Drc);
        double ff = UF_t2->Drc;
        double sf = UF_t3->Drc;
        double gamma = 0.5;//1000.0/_d->Drc;
        double dc = 0.5;//UF_DragCoefficient(ff,sf,gamma,_visc->Drc, _rocksize->Drc, _d->Drc);
        double pbf = -UF_Gravity*h;
        double pbs = 1-gamma*pbf;
        double ifa = 0.3;

        double dhdx = UF1D_Derivative(_ldd,_lddw,UF_t1,r,c);
        double dbdx = UF1D_Derivative(_ldd,_lddw,UF_t2,r,c);

        if(sf > UF_VERY_SMALL)
        {
            out_su->Drc += dt->Drc *
                (UF_Gravity * sin(UF1D_Slope->Drc) - (_su->Drc > 0? 1.0 : -1.0)*std::tan(ifa)*pbs-UF_Aspect*pbs*UF1D_Slope->Drc
                -UF_Aspect * gamma * pbf * ( dhdx +  UF1D_Slope->Drc) -
                - dbdx
                 );

            double vpow = pow(std::fabs((_fu->Drc-_su->Drc)),0.5*(UF_j-1.0));
            double facu = std::min(1.0,std::max(-1.0,std::fabs(dt->Drc * dc * ff *(_su->Drc - _fu->Drc) *vpow)));
            out_su->Drc = (1.0-facu) * out_su->Drc + (facu) * _fu->Drc;
        }
    }

}
void TWorld::UF1D_FluidMomentumSource(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu)
{
    FOR_ROW_COL_UF1D
    {
        out_fu->Drc = _fu->Drc;
        double h = (_f->Drc + _s->Drc)/(_dx * _lddw->Drc);
        UF_t1->Drc = h*h*(-UF_Gravity*h)/2.0;
        if(_f->Drc + _s->Drc > UF_VERY_SMALL)
        {
            UF_t2->Drc = _f->Drc/(_f->Drc + _s->Drc);
            UF_t3->Drc = _s->Drc/(_f->Drc + _s->Drc);
        }else
        {
            UF_t2->Drc = 0;
            UF_t3->Drc = 0;
        }
    }

    FOR_ROW_COL_UF1D_DT
    {
        if((_f->Drc < UF_VERY_SMALL))
        {
            out_fu->Drc = 0;
            continue;
        }

        double h = (_f->Drc + _s->Drc)/(_dx *_lddw->Drc);
        double ff = UF_t2->Drc;
        double sf = UF_t3->Drc;
        double Nr = 1500;//UF_Reynolds(1000.0,_visc->Drc,ff);
        double Nra = 15000;//UF_Reynolds(1000.0,_visc->Drc,ff);
        double gamma = 0.5;//_d->Drc == 0? 0.5 : 1000.0/_d->Drc;
        double dc = 0.5;//UF_DragCoefficient(ff,sf,UF_TerminalVelocity(UF2D_rocksize->Drc,ff,_visc->Drc,sf,_d->Drc),gamma);
        double pbf = -UF_Gravity*h;

        double dh2pbdx = UF1D_Derivative(_ldd,_lddw,UF_t1,r,c);
        double dsfdx = UF1D_Derivative(_ldd,_lddw,UF_t3,r,c);
        double ddsfdxx = UF1D_Derivative2(_ldd,_lddw,UF_t3,r,c);
        double dfudx = UF1D_Derivative(_ldd,_lddw,_fu,r,c);
        double ddfudxx = UF1D_Derivative2(_ldd,_lddw,_fu,r,c);
        double dsudx = UF1D_Derivative(_ldd,_lddw,_su,r,c);

        if(ff > UF_VERY_SMALL)
        {
            out_fu->Drc += dt->Drc *
                (UF_Gravity * sin(UF1D_Slope->Drc) -
                 UF_Aspect *(
                     (dh2pbdx)/h
                     +(pbf * UF1D_Slope->Drc)
                     -(1.0/(ff * Nr))*(
                         2.0*ddfudxx + UF_Chi * _fu->Drc/(UF_Aspect*UF_Aspect * h* h)
                         )
                     +(1.0/(ff * Nra))*(
                         2.0 *(dsfdx*(dfudx - dsudx) + ddsfdxx*(_fu->Drc - _su->Drc))
                         + UF_Chi * _fu->Drc/(UF_Aspect*UF_Aspect*ff*Nra*h*h)
                         )
                     -(UF_Ksi*sf*(_fu->Drc - _su->Drc)/(UF_Aspect*UF_Aspect*ff*Nra*h*h))
                     )
                  );

            double vpow = pow(std::fabs(_fu->Drc-_su->Drc),0.5*(UF_j-1.0));
            double facu = std::min(1.0,std::max(0.0,std::fabs(dt->Drc * dc * sf *(_fu->Drc - _su->Drc) *vpow)));
            out_fu->Drc = (1.0-facu) * out_fu->Drc + (facu) * _su->Drc;
        }

        double nsq = (0.1 + N->Drc)*(0.1 + N->Drc)*UF_Gravity*fabs(_fu->Drc)*dt->Drc/pow(h,4.0/3.0) ;
        double fac = h > 0? 1.0/std::min(h* 3.0,1.0) : 0.0;
        out_fu->Drc = (out_fu->Drc/(1.0+nsq * fac));

    }

}

void TWorld::UF1D_Diffuse_mass(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap * _m, cTMap * _f,cTMap * _fu,cTMap * _s,cTMap * _su, cTMap * out_m)
{





}
