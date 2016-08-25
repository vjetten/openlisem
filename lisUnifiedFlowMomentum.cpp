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
 \file lisUnifiedFlowMomentum.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"

void TWorld::UF2D_FluidApplyMomentum(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv)
{
    FOR_ROW_COL_UF2D
    {
        out_fu->Drc = _fu->Drc;
        out_fv->Drc = _fv->Drc;
    }

    FOR_ROW_COL_UF2D_DT
    {

        if(!(_f->Drc > UF_VERY_SMALL))
        {
            out_fu->Drc = 0;
            out_fv->Drc = 0;
            continue;
        }

        out_fu->Drc += dt->Drc * UF2D_fax->Drc;
        out_fv->Drc += dt->Drc * UF2D_fay->Drc;

    }}}

}

void TWorld::UF2D_SolidApplyMomentum(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv)
{
    FOR_ROW_COL_UF2D
    {
        out_su->Drc = _su->Drc;
        out_sv->Drc = _sv->Drc;
    }

    FOR_ROW_COL_UF2D_DT
    {
        if(!(_s->Drc > UF_VERY_SMALL))
        {
            out_su->Drc = 0;
            out_sv->Drc = 0;
            continue;
        }


        out_su->Drc += dt->Drc * UF2D_sax->Drc;
        out_sv->Drc += dt->Drc * UF2D_say->Drc;

    }}}
}

void TWorld::UF1D_FluidApplyMomentum(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu)
{
    FOR_ROW_COL_UF1D
    {
        out_fu->Drc = _fu->Drc;
    }

    FOR_ROW_COL_UF1D_DT
    {
        if(!(_f->Drc > UF_VERY_SMALL))
        {
            out_fu->Drc = 0;
            continue;
        }
        out_fu->Drc += dt->Drc * UF1D_fa->Drc;
    }
}


void TWorld::UF1D_SolidApplyMomentum(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_su)
{
    FOR_ROW_COL_UF1D
    {
        out_su->Drc = _su->Drc;
    }

    FOR_ROW_COL_UF1D_DT
    {
        if(!(_s->Drc > UF_VERY_SMALL))
        {
            out_su->Drc = 0;
            continue;
        }
        out_su->Drc += dt->Drc * UF1D_sa->Drc;
    }
}


void TWorld::UF2D_FluidMomentumSource(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv)
{
    FOR_ROW_COL_UF2D
    {
        double h = (_f->Drc + _s->Drc)/(_dx * _dx);
        UF_t1->Drc = h*h*(UF_Gravity*h)/2.0;

        if(_f->Drc + _s->Drc > UF_VERY_SMALL)
        {
            UF_t2->Drc = _f->Drc/(_f->Drc + _s->Drc);
            UF_t3->Drc = _s->Drc/(_f->Drc + _s->Drc);
        }else
        {
            UF_t2->Drc = 0;
            UF_t3->Drc = 0;
        }

        UF_t4->Drc = _f->Drc/(_dx * _dx);
        UF_t5->Drc = _s->Drc/(_dx * _dx);
        UF2D_fax->Drc = 0;
        UF2D_fay->Drc = 0;
    }

    FOR_ROW_COL_UF2D_DT
    {
        double h = (_f->Drc + _s->Drc)/(_dx * _dx);
        double ff = UF_t2->Drc;
        double sf = UF_t3->Drc;
        double visc = UF_DynamicViscosity(sf);
        double Nr = std::max(0.5,UF_Reynolds(_d->Drc,visc,ff,sf, _rocksize->Drc));
        double Nra = 150000;
        double gamma = _d->Drc > UF_VERY_SMALL? 1000.0/_d->Drc : 0.5;
        double dc = UF_DragCoefficient(ff, sf, gamma,visc, _rocksize->Drc, _d->Drc);
        double pbf = UF_Gravity*h;

        double dhfdx = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_X);
        double dhfdy = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_Y);

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
        double ddfudyy = UF2D_Derivative2(_dem,_fu,r,c,UF_DIRECTION_Y);
        double ddfvdxy = UF2D_Derivative2(_dem,_fv,r,c,UF_DIRECTION_XY);
        double ddfvdxx = UF2D_Derivative2(_dem,_fv,r,c,UF_DIRECTION_X);
        double ddfvdyy = UF2D_Derivative2(_dem,_fv,r,c,UF_DIRECTION_Y);
        double ddfudxy = UF2D_Derivative2(_dem,_fu,r,c,UF_DIRECTION_XY);

        double dsudx = UF2D_Derivative(_dem,_su,r,c,UF_DIRECTION_X);
        double dsudy = UF2D_Derivative(_dem,_su,r,c,UF_DIRECTION_Y);
        double dsvdx = UF2D_Derivative(_dem,_sv,r,c,UF_DIRECTION_X);
        double dsvdy = UF2D_Derivative(_dem,_sv,r,c,UF_DIRECTION_Y);

        if(!(_f->Drc > UF_VERY_SMALL))
        {

            UF2D_fax->Drc = 0;
            UF2D_fay->Drc = 0;
            continue;
        }

        UF2D_fax->Drc = UF2D_MomentumBalanceFluid(true,_f->Drc,_s->Drc,_fu->Drc, _fv->Drc, _su->Drc, _sv->Drc, ff, sf, Nr, Nra, _ifa->Drc, gamma, visc, pbf, UF2D_SlopeX->Drc, UF2D_SlopeY->Drc,
                                                 dhfdx,  dhfdy,   dh2pbdx,   dh2pbdy,   dsfdx,   dsfdy,   ddsfdxx,   ddsfdyy,   ddsfdxy,
                                                 dfudx,   dfudy,   dfvdx,   dfvdy,   ddfudxx,   ddfudyy,   ddfvdxy,   ddfvdxx,   ddfvdyy,   ddfudxy,
                                                   dsudx,   dsudy,   dsvdx,  dsvdy);

        UF2D_fay->Drc = UF2D_MomentumBalanceFluid(false, _f->Drc,_s->Drc,_fu->Drc, _fv->Drc, _su->Drc, _sv->Drc, ff, sf, Nr, Nra, _ifa->Drc, gamma, visc, pbf, UF2D_SlopeX->Drc, UF2D_SlopeY->Drc,
                                                   dhfdx,  dhfdy,   dh2pbdx,   dh2pbdy,   dsfdx,   dsfdy,   ddsfdxx,   ddsfdyy,   ddsfdxy,
                                                   dfudx,   dfudy,   dfvdx,   dfvdy,   ddfudxx,   ddfudyy,   ddfvdxy,   ddfvdxx,   ddfvdyy,   ddfudxy,
                                                   dsudx,   dsudy,   dsvdx,  dsvdy);

        UF_t6->Drc = _fu->Drc;
        UF_t7->Drc = _fv->Drc;
        double dttemp = dt->Drc/((double)UF_FrictionIterations);
        for(int i = 0; i < UF_FrictionIterations; i++)
        {
            UF_t6->Drc += dttemp* UF2D_fax->Drc;
            UF_t7->Drc += dttemp* UF2D_fay->Drc;

            double tempu = UF_t6->Drc;
            UF_t6->Drc = UF_Friction(dttemp,UF_t6->Drc,UF_t7->Drc,N->Drc,h,UF2D_SlopeX->Drc);
            UF_t7->Drc = UF_Friction(dttemp,UF_t7->Drc,tempu,N->Drc,h,UF2D_SlopeY->Drc);
        }

        if(sf > UF_VERY_SMALL)
        {
            double vpow = pow((UF_t6->Drc-_su->Drc)*(UF_t6->Drc-_su->Drc)+(UF_t7->Drc-_sv->Drc)*(UF_t7->Drc-_sv->Drc),0.5*(UF_j-1.0));
            double facu = std::min(1.0,std::max(0.0,std::fabs(dt->Drc * dc * sf *(UF_t6->Drc - _su->Drc) *vpow)));
            double facv = std::min(1.0,std::max(0.0,std::fabs(dt->Drc * dc * sf *(UF_t7->Drc - _sv->Drc) *vpow)));
            UF_t6->Drc = (1.0-facu) * UF_t6->Drc + (facu) * _su->Drc;
            UF_t7->Drc = (1.0-facv) * UF_t7->Drc + (facv) * _sv->Drc;
        }

        UF2D_fax->Drc = (UF_t6->Drc - _fu->Drc)/dt->Drc;
        UF2D_fay->Drc = (UF_t7->Drc - _fv->Drc)/dt->Drc;

    }}}

}

void TWorld::UF2D_SolidMomentumSource(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv)
{
    FOR_ROW_COL_UF2D
    {

        double h = (_f->Drc + _s->Drc)/(_dx * _dx);
        double hs = (_s->Drc)/(_dx * _dx);
        UF_t1->Drc = h;

        double ifa = _ifa->Drc;

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

        UF_t4->Drc = _f->Drc/(_dx * _dx);
        UF_t5->Drc = _s->Drc/(_dx * _dx);
        UF2D_sax->Drc = 0;
        UF2D_say->Drc = 0;
    }

    FOR_ROW_COL_UF2D_DT
    {
        if(!(_s->Drc > UF_VERY_SMALL))
        {
            UF2D_sax->Drc = 0;
            UF2D_say->Drc = 0;
            continue;
        }

        double h = (_f->Drc + _s->Drc)/(_dx * _dx);
        double ff = _f->Drc / (_f->Drc + _s->Drc);
        double sf = _s->Drc / (_f->Drc + _s->Drc);
        double gamma = _d->Drc > UF_VERY_SMALL? 1000.0/_d->Drc : 0.5;
        double dc = UF_DragCoefficient(ff,sf,gamma,_visc->Drc, _rocksize->Drc, _d->Drc);
        double pbf = -UF_Gravity*h;
        double pbs = (1-gamma)*pbf;
        double ifa = 0.3;

        double dhsdx = UF2D_Derivative(_dem,UF_t5,r,c,UF_DIRECTION_X);
        double dhsdy = UF2D_Derivative(_dem,UF_t5,r,c,UF_DIRECTION_Y);

        double dhdx = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_X);
        double dhdy = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_Y);
        double dbdx = UF2D_Derivative(_dem,UF_t2,r,c,UF_DIRECTION_X);
        double dbdy = UF2D_Derivative(_dem,UF_t3,r,c,UF_DIRECTION_Y);

        UF2D_sax->Drc = UF2D_MomentumBalanceSolid(true,_f->Drc,_s->Drc,_fu->Drc, _fv->Drc, _su->Drc, _sv->Drc, ff, sf, 0, 0, ifa, gamma, _visc->Drc,pbs, pbf, UF2D_SlopeX->Drc, UF2D_SlopeY->Drc,
                                                  dhsdx, dhsdy, dhdx, dhdy, dbdx, dbdy);

        UF2D_say->Drc = UF2D_MomentumBalanceSolid(false,_f->Drc,_s->Drc,_fu->Drc, _fv->Drc, _su->Drc, _sv->Drc, ff, sf, 0, 0, ifa, gamma, _visc->Drc,pbs, pbf, UF2D_SlopeX->Drc, UF2D_SlopeY->Drc,
                                                 dhsdx, dhsdy, dhdx, dhdy, dbdx, dbdy);



        UF_t6->Drc = _su->Drc;
        UF_t7->Drc = _sv->Drc;

        double dttemp = dt->Drc/((double)UF_FrictionIterations);
        for(int i = 0; i < UF_FrictionIterations; i++)
        {
            UF_t6->Drc += dttemp * UF2D_sax->Drc;
            UF_t7->Drc += dttemp * UF2D_say->Drc;

            double tempu = UF_t6->Drc;
            UF_t6->Drc = UF_Friction(dttemp,UF_t6->Drc,UF_t7->Drc,N->Drc,h,UF2D_SlopeX->Drc);
            UF_t7->Drc = UF_Friction(dttemp,UF_t7->Drc,tempu,N->Drc,h,UF2D_SlopeY->Drc);
        }

        if(ff > UF_VERY_SMALL)
        {
            double vpow = pow((_fu->Drc-UF_t6->Drc)*(_fu->Drc-UF_t6->Drc)+(_fv->Drc-UF_t7->Drc)*(_fv->Drc-UF_t7->Drc),0.5*(UF_j-1.0));
            double facu = std::min(1.0,std::max(0.0,std::fabs(dt->Drc * dc * ff *(UF_t6->Drc - _fu->Drc) *vpow)));
            double facv = std::min(1.0,std::max(0.0,std::fabs(dt->Drc * dc * ff *(UF_t7->Drc - _fv->Drc) *vpow)));
            UF_t6->Drc = (1.0-facu) * UF_t6->Drc + (facu) * _fu->Drc;
            UF_t7->Drc = (1.0-facv) * UF_t7->Drc + (facv) * _fv->Drc;
        }

        UF2D_sax->Drc = (UF_t6->Drc - _su->Drc)/dt->Drc;
        UF2D_say->Drc = (UF_t7->Drc - _sv->Drc)/dt->Drc;

    }}}

}


void TWorld::UF1D_FluidMomentumSource(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu)
{
    FOR_ROW_COL_UF1D
    {
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

        UF_t4->Drc = _f->Drc/(_dx * _dx);
        UF_t5->Drc = _s->Drc/(_dx * _dx);
        UF1D_fa->Drc = 0;
    }

    FOR_ROW_COL_UF1D_DT
    {
        if(!(_f->Drc > UF_VERY_SMALL))
        {
            UF1D_fa->Drc = 0;
            continue;
        }
        double h = (_f->Drc + _s->Drc)/(_dx *_lddw->Drc);
        double ff = UF_t2->Drc;
        double sf = UF_t3->Drc;
        _visc->Drc = UF_DynamicViscosity(sf);
        if(_d->Drc == 0)
        {
            qDebug() << r << c << 0;
        }
        _d->Drc = 2000;
        double Nr = UF_Reynolds(_d->Drc,_visc->Drc,ff,sf, _rocksize->Drc);
        double Nra = 15000;
        double gamma =(!(_d->Drc > UF_VERY_SMALL))? 0.5 : 1000.0/_d->Drc;
        double dc = UF_DragCoefficient(ff, sf, gamma,_visc->Drc, _rocksize->Drc, _d->Drc);
        double pbf = -UF_Gravity*h;

        double dhfdx = UF1D_Derivative(_ldd,_lddw,UF_t4,r,c);

        double dh2pbdx = UF1D_Derivative(_ldd,_lddw,UF_t1,r,c);
        double dsfdx = UF1D_Derivative(_ldd,_lddw,UF_t3,r,c);
        double ddsfdxx = UF1D_Derivative2(_ldd,_lddw,UF_t3,r,c);
        double dfudx = UF1D_Derivative(_ldd,_lddw,_fu,r,c);
        double ddfudxx = UF1D_Derivative2(_ldd,_lddw,_fu,r,c);
        double dsudx = UF1D_Derivative(_ldd,_lddw,_su,r,c);

        UF1D_fa->Drc = UF1D_MomentumBalanceFluid( _f->Drc,_s->Drc,_fu->Drc, _su->Drc, ff, sf, Nr, Nra, _ifa->Drc, gamma, _visc->Drc, pbf, UF1D_Slope->Drc,
                                                 dhfdx, dh2pbdx, dsfdx, ddsfdxx, dfudx, ddfudxx, dsudx);

        UF_t6->Drc = _fu->Drc;

        double dttemp = dt->Drc/((double)UF_FrictionIterations);
        for(int i = 0; i < UF_FrictionIterations; i++)
        {

            UF_t6->Drc += dttemp * UF1D_fa->Drc;
            UF_t6->Drc = UF_Friction(dttemp,UF_t6->Drc,0,N->Drc,h,UF1D_Slope->Drc);
        }

        if(sf > UF_VERY_SMALL)
        {
            double vpow = pow(std::fabs(UF_t6->Drc-_su->Drc),0.5*(UF_j-1.0));
            double facu = std::min(1.0,std::max(0.0,std::fabs(dt->Drc * dc * sf *(UF_t6->Drc - _su->Drc) *vpow)));
            UF_t6->Drc = (1.0-facu) * UF_t6->Drc + (facu) * _su->Drc;
        }

        UF1D_fa->Drc = (UF_t6->Drc - _fu->Drc)/dt->Drc;
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

        double ifa = _ifa->Drc;

        double ifasin = sin(ifa);
        double ifatan = tan(ifa);

        double dudxx = UF1D_Derivative2(_ldd,_lddw,_su,r,c);

        double k_act = (1- ifasin)/(1+ifasin);
        double k_pass = (1+ifasin)/(1-ifasin);
        double k_x = (dudxx) > 0? k_act: k_pass;

        UF_t2->Drc = UF_Gravity * k_x *hs * hs * UF_Aspect/2.0;

        UF_t4->Drc = _f->Drc/(_dx * _dx);
        UF_t5->Drc = _s->Drc/(_dx * _dx);
        UF1D_sa->Drc = 0;
    }

    FOR_ROW_COL_UF1D_DT
    {
        if(!(_s->Drc > UF_VERY_SMALL))
        {
            UF1D_sa->Drc = 0;
            continue;
        }

        double h = (_f->Drc + _s->Drc)/(_dx * _lddw->Drc);
        double ff = UF_t2->Drc;
        double sf = UF_t3->Drc;
        double gamma = _d->Drc > UF_VERY_SMALL? 1000.0/_d->Drc : 0.5;
        double dc = UF_DragCoefficient(ff,sf,gamma,_visc->Drc, _rocksize->Drc, _d->Drc);
        double pbf = -UF_Gravity*h;
        double pbs = 1-gamma*pbf;
        double ifa = _ifa->Drc;

        double dhsdx = UF1D_Derivative(_ldd,_lddw,UF_t5,r,c);

        double dhdx = UF1D_Derivative(_ldd,_lddw,UF_t1,r,c);
        double dbdx = UF1D_Derivative(_ldd,_lddw,UF_t2,r,c);

        UF1D_sa->Drc = UF1D_MomentumBalanceSolid(_f->Drc,_s->Drc,_fu->Drc, _su->Drc, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs,pbf, UF1D_Slope->Drc,
                                                 dhsdx, dhdx, dbdx);

        UF_t6->Drc = _su->Drc;

        double dttemp = dt->Drc/((double)UF_FrictionIterations);
        for(int i = 0; i < UF_FrictionIterations; i++)
        {
            UF_t6->Drc += dttemp * UF1D_sa->Drc;
            UF_t6->Drc = UF_Friction(dttemp,UF_t6->Drc,0,N->Drc,h,UF1D_Slope->Drc);
        }

        if(ff > UF_VERY_SMALL)
        {
            double vpow = pow(std::fabs((_fu->Drc-UF_t6->Drc)),0.5*(UF_j-1.0));
            double facu = std::min(1.0,std::max(0.0,std::fabs(dt->Drc * dc * ff *(UF_t6->Drc - _fu->Drc) *vpow)));
            UF_t6->Drc = (1.0-facu) * UF_t6->Drc + (facu) * _fu->Drc;
        }

        UF1D_sa->Drc = (UF_t6->Drc - _su->Drc)/dt->Drc;

    }

}


