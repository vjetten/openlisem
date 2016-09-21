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

void TWorld::UF2D_FluidApplyMomentum2(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv)
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

void TWorld::UF2D_SolidApplyMomentum2(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv)
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

void TWorld::UF1D_FluidApplyMomentum2(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu)
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


void TWorld::UF1D_SolidApplyMomentum2(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_su)
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


void TWorld::UF2D_FluidMomentum2Source(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv)
{
    FOR_ROW_COL_UF2D
    {
        UF_t1->Drc = (_f->Drc)/(_dx * _dx);
    }

    UF2D_MUSCLE(_dem,dt,UF_t1,UF_MUSCLE_TARGET_IN1);
    UF2D_MUSCLE(_dem,dt,_fu,UF_MUSCLE_TARGET_IN2, UF_DIRECTION_X);
    UF2D_MUSCLE(_dem,dt,_fv,UF_MUSCLE_TARGET_IN3, UF_DIRECTION_Y);

    FOR_ROW_COL_UF2D
    {
        double h = (_s->Drc + _f->Drc)/(_dx * _dx);
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

        if(!(_f->Drc > UF_VERY_SMALL) ||!(_s->Drc + _f->Drc > UF_VERY_SMALL) )
        {

            UF2D_fax->Drc = 0;
            UF2D_fay->Drc = 0;
            UF2D_fqx1->Drc = 0;
            UF2D_fqx2->Drc = 0;
            UF2D_fqy1->Drc = 0;
            UF2D_fqy2->Drc = 0;
            UF2D_fax1->Drc = 0;
            UF2D_fax2->Drc = 0;
            UF2D_fay1->Drc = 0;
            UF2D_fay2->Drc = 0;
            continue;
        }

        ////////////input parameters

        double lfu = UF2D_MUSCLE_2_x2->Drc;
        double lfv = UF2D_MUSCLE_3_y2->Drc;
        double rfu = UF2D_MUSCLE_2_x1->Drc;
        double rfv = UF2D_MUSCLE_3_y1->Drc;
        double lxh = std::max(0.0,UF2D_MUSCLE_1_x2->Drc);
        double rxh = std::max(0.0,UF2D_MUSCLE_1_x1->Drc);
        double lyh = std::max(0.0,UF2D_MUSCLE_1_y2->Drc);
        double ryh = std::max(0.0,UF2D_MUSCLE_1_y1->Drc);

        double lxslope = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double lyslope = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rxslope = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double ryslope = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double ff = UF_t2->Drc;
        double sf = UF_t3->Drc;
        _visc->Drc = UF_DynamicViscosity(sf + (SwitchErosion? (((UF2D_blm->Drc + UF2D_ssm->Drc)/2000.0 + ff) > UF_VERY_SMALL ?((UF2D_blm->Drc + UF2D_ssm->Drc)/2000.0)/((UF2D_blm->Drc + UF2D_ssm->Drc)/2000.0 + ff):0.0):0.0));
        double Nr = std::max(0.5,UF_Reynolds(_d->Drc,_visc->Drc,ff,sf, _rocksize->Drc));
        UF2D_Nr->Drc = Nr;
        double Nra = UF_NRA;
        double gamma = _d->Drc > UF_VERY_SMALL? 1000.0/_d->Drc : 0.5;
        double dc = UF_DragCoefficient(ff,sf,gamma,_visc->Drc, _rocksize->Drc, _d->Drc);

        double lxpbf = UF_Gravity*lxh;
        double rxpbf = UF_Gravity*rxh;
        double lypbf = UF_Gravity*lyh;
        double rypbf = UF_Gravity*ryh;

        double dhfdx = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_X);
        double dhfdy = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_Y);
        double ldhfdx = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double ldhfdy = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rdhfdx = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rdhfdy = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        /*double dhfdx = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_X, UF_DERIVATIVE_LR, true);
        double dhfdy = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_Y, UF_DERIVATIVE_LR,true);
        double ldhfdx = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_X,UF_DERIVATIVE_L, true);
        double ldhfdy = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L, true);
        double rdhfdx = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_X,UF_DERIVATIVE_R, true);
        double rdhfdy = UF2D_Derivative(_dem,UF_t4,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R, true);*/

        double dh2pbdx = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_X);
        double dh2pbdy = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_Y);
        double ldh2pbdx = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double ldh2pbdy = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rdh2pbdx = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rdh2pbdy = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

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

        ////////////momentum balance
        double lfax = UF2D_MomentumBalanceFluid(true,lxh*_dx*_dx,_s->Drc,lfu, _fv->Drc, _su->Drc, _sv->Drc, ff, sf, Nr, Nra, _ifa->Drc, gamma, _visc->Drc, lxpbf, lxslope, 0,
                                                 ldhfdx,  dhfdy,   ldh2pbdx,   dh2pbdy,   dsfdx,   dsfdy,   ddsfdxx,   ddsfdyy,   ddsfdxy,
                                                 dfudx,   dfudy,   dfvdx,   dfvdy,   ddfudxx,   ddfudyy,   ddfvdxy,   ddfvdxx,   ddfvdyy,   ddfudxy,
                                                   dsudx,   dsudy,   dsvdx,  dsvdy);

        double rfax = UF2D_MomentumBalanceFluid(true, rxh*_dx*_dx,_s->Drc,rfu, _fv->Drc, _su->Drc, _sv->Drc, ff, sf, Nr, Nra, _ifa->Drc, gamma, _visc->Drc, rxpbf, rxslope, 0,
                                                   rdhfdx,  dhfdy,   rdh2pbdx,   dh2pbdy,   dsfdx,   dsfdy,   ddsfdxx,   ddsfdyy,   ddsfdxy,
                                                   dfudx,   dfudy,   dfvdx,   dfvdy,   ddfudxx,   ddfudyy,   ddfvdxy,   ddfvdxx,   ddfvdyy,   ddfudxy,
                                                   dsudx,   dsudy,   dsvdx,  dsvdy);

        double lfay = UF2D_MomentumBalanceFluid(false,lyh*_dx*_dx,_s->Drc,_fu->Drc, lfv, _su->Drc, _sv->Drc, ff, sf, Nr, Nra, _ifa->Drc, gamma, _visc->Drc, lypbf, 0, lyslope,
                                                 dhfdx,  ldhfdy,   dh2pbdx,   ldh2pbdy,   dsfdx,   dsfdy,   ddsfdxx,   ddsfdyy,   ddsfdxy,
                                                 dfudx,   dfudy,   dfvdx,   dfvdy,   ddfudxx,   ddfudyy,   ddfvdxy,   ddfvdxx,   ddfvdyy,   ddfudxy,
                                                   dsudx,   dsudy,   dsvdx,  dsvdy);

        double rfay = UF2D_MomentumBalanceFluid(false,ryh*_dx*_dx,_s->Drc,_fu->Drc, rfv, _su->Drc, _sv->Drc, ff, sf, Nr, Nra, _ifa->Drc, gamma, _visc->Drc, rypbf, 0, ryslope,
                                                   dhfdx,  rdhfdy,   dh2pbdx,   rdh2pbdy,   dsfdx,   dsfdy,   ddsfdxx,   ddsfdyy,   ddsfdxy,
                                                   dfudx,   dfudy,   dfvdx,   dfvdy,   ddfudxx,   ddfudyy,   ddfvdxy,   ddfvdxx,   ddfvdyy,   ddfudxy,
                                                   dsudx,   dsudy,   dsvdx,  dsvdy);

        UF2D_fax1->Drc = rfax;
        UF2D_fax2->Drc = lfax;
        UF2D_fay1->Drc = rfay;
        UF2D_fay2->Drc = lfay;

        ////////////friction and actual accaleration

            UF_t6->Drc = lfu;
            UF_t7->Drc = rfu;

            UF_t6->Drc += dt->Drc*(2.0/3.0)* lfax;
            UF_t7->Drc += dt->Drc*(2.0/3.0)* rfax;
            UF_t6->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t6->Drc,0,N->Drc,lxh,lxslope,false);
            UF_t7->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t7->Drc,0,N->Drc,rxh,rxslope,false);
            double ut1 = UF_t6->Drc + dt->Drc*(2.0/3.0)* lfax;
            double ut2 = UF_t7->Drc + dt->Drc*(2.0/3.0)* rfax;
            UF_t6->Drc = (UF_t6->Drc + UF_Friction(dt->Drc*(2.0/3.0),ut1,0,N->Drc,lxh,lxslope,false))/2.0;
            UF_t7->Drc = (UF_t7->Drc + UF_Friction(dt->Drc*(2.0/3.0),ut2,0,N->Drc,rxh,rxslope,false))/2.0;

            if(sf > UF_VERY_SMALL)
            {
                double lvpow = pow((UF_t6->Drc-_su->Drc)*(UF_t6->Drc-_su->Drc)+(_fv->Drc-_sv->Drc)*(_fv->Drc-_sv->Drc),0.5*(UF_j-1.0));
                double rvpow = pow((UF_t7->Drc-_su->Drc)*(UF_t7->Drc-_su->Drc)+(_fv->Drc-_sv->Drc)*(_fv->Drc-_sv->Drc),0.5*(UF_j-1.0));
                double lfacu = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * sf*lvpow)));
                double rfacu = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * sf *rvpow)));
                UF_t6->Drc = (1.0-lfacu) * UF_t6->Drc + (lfacu) * _su->Drc;
                UF_t7->Drc = (1.0-rfacu) * UF_t7->Drc + (rfacu) * _su->Drc;
            }

            lfax = (UF_t6->Drc - lfu)/dt->Drc;
            rfax = (UF_t7->Drc - rfu)/dt->Drc;

            UF_t6->Drc = lfv;
            UF_t7->Drc = rfv;

            UF_t6->Drc += dt->Drc*(2.0/3.0)* lfay;
            UF_t7->Drc += dt->Drc*(2.0/3.0)* rfay;
            UF_t6->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t6->Drc,0,N->Drc,lyh,lyslope,false);
            UF_t7->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t7->Drc,0,N->Drc,ryh,ryslope,false);
            double vt1 = UF_t6->Drc + dt->Drc*(2.0/3.0)* lfay;
            double vt2 = UF_t7->Drc + dt->Drc*(2.0/3.0)* rfay;
            UF_t6->Drc = (UF_t6->Drc + UF_Friction(dt->Drc*(2.0/3.0),vt1,0,N->Drc,lyh,lyslope,false))/2.0;
            UF_t7->Drc = (UF_t7->Drc + UF_Friction(dt->Drc*(2.0/3.0),vt2,0,N->Drc,ryh,ryslope,false))/2.0;

            if(sf > UF_VERY_SMALL)
            {
                double lvpow = pow((UF_t6->Drc-_sv->Drc)*(UF_t6->Drc-_sv->Drc)+(_fu->Drc-_su->Drc)*(_fu->Drc-_su->Drc),0.5*(UF_j-1.0));
                double rvpow = pow((UF_t7->Drc-_sv->Drc)*(UF_t7->Drc-_sv->Drc)+(_fu->Drc-_su->Drc)*(_fu->Drc-_su->Drc),0.5*(UF_j-1.0));
                double lfacv = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * sf *lvpow)));
                double rfacv = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * sf *rvpow)));
                UF_t6->Drc = (1.0-lfacv) * UF_t6->Drc + (lfacv) * _sv->Drc;
                UF_t7->Drc = (1.0-rfacv) * UF_t7->Drc + (rfacv) * _sv->Drc;
            }

            lfay = (UF_t6->Drc - lfv)/dt->Drc;
            rfay = (UF_t7->Drc - rfv)/dt->Drc;

        ////////////average accalerations for cell centers
        UF2D_fax->Drc = (lfax + rfax)/2.0;
        UF2D_fay->Drc = (lfay + rfay)/2.0;

        ////////////fluxes
        double hdem = (_f->Drc/(_dx*_dx) + _dem->Drc);
        double volx1 = UF_OUTORMV(_dem,r,c+1)?0.0:std::max(0.0,(_dx*_dx) * (hdem - (_f->data[r][c+1]/(_dx*_dx) + _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1))));
        double volx2 = UF_OUTORMV(_dem,r,c-1)?0.0:std::max(0.0,(_dx*_dx) * (hdem - (_f->data[r][c-1]/(_dx*_dx) + _dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1))));
        double voly1 = UF_OUTORMV(_dem,r+1,c)?0.0:std::max(0.0,(_dx*_dx) * (hdem - (_f->data[r+1][c]/(_dx*_dx) + _dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0))));
        double voly2 = UF_OUTORMV(_dem,r-1,c)?0.0:std::max(0.0,(_dx*_dx) * (hdem - (_f->data[r-1][c]/(_dx*_dx) + _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0))));

        double dtx1 = dt->Drc;//(UF2D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
        double dtx2 = dt->Drc;//(UF2D_MUSCLE_OUT_x2->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c-1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c-1)? 0.0: dt->data[r][c-1]));
        double dty1 = dt->Drc;//(UF2D_MUSCLE_OUT_y1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r+1,c)? dt->Drc : (UF_NOTIME(_dem,dt,r+1,c)? 0.0: dt->data[r+1][c]));
        double dty2 = dt->Drc;//(UF2D_MUSCLE_OUT_y2->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r-1,c)? dt->Drc : (UF_NOTIME(_dem,dt,r-1,c)? 0.0: dt->data[r-1][c]));

        double vx1 = rfu + dtx1 * rfax;
        double vx2 = lfu + dtx2 * lfax;
        double vy1 = rfv + dty1 * rfay;
        double vy2 = lfv + dty2 * lfay;

        rxh = UF_OUTORMV(_dem,r,c+1)?0.0:std::min(rxh,std::max(0.0,(vx1 > 0? 1.0:0.0)*((_dem->Drc + rxh) - (_f->data[r][c+1]/(_dx*_dx) + _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1)))));
        lxh = UF_OUTORMV(_dem,r,c-1)?0.0:std::min(lxh,std::max(0.0,(vx2 > 0? 0.0:1.0)*((_dem->Drc + lxh) - (_f->data[r][c-1]/(_dx*_dx) + _dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1)))));
        ryh = UF_OUTORMV(_dem,r+1,c)?0.0:std::min(ryh,std::max(0.0,(vy1 > 0? 1.0:0.0)*((_dem->Drc + ryh) - (_f->data[r+1][c]/(_dx*_dx) + _dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0)))));
        lyh = UF_OUTORMV(_dem,r-1,c)?0.0:std::min(lyh,std::max(0.0,(vy2 > 0? 0.0:1.0)*((_dem->Drc + lyh) - (_f->data[r-1][c]/(_dx*_dx) + _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0)))));

        double cq = 0.25 * UF_Courant * _f->Drc;
        double cqx1 = 0.25 * UF_Courant * (UF_OUTORMV(_dem,r,c+1)? 0.0 : volx1);
        double cqx2 = 0.25 * UF_Courant * (UF_OUTORMV(_dem,r,c-1)? 0.0 : volx2);
        double cqy1 = 0.25 * UF_Courant * (UF_OUTORMV(_dem,r+1,c)? 0.0 : voly1);
        double cqy2 = 0.25 * UF_Courant * (UF_OUTORMV(_dem,r-1,c)? 0.0 : voly2);

        double qx1 = 0.5 *dtx1 * (vx1) * rxh *_dx;
        double qx2 = 0.5 *dtx2 * (vx2) * lxh *_dx;
        double qy1 = 0.5 *dty1 * (vy1) * ryh *_dx;
        double qy2 = 0.5 *dty2 * (vy2) * lyh *_dx;

        qx1 = ((qx1 > 0)? 1.0 : -1.0) * std::min(std::fabs(qx1),(qx1 > 0)? cq : cqx1);
        qx2 = ((qx2 > 0)? 1.0 : -1.0) * std::min(std::fabs(qx2),(qx2 > 0)? cqx2 : cq);
        qy1 = ((qy1 > 0)? 1.0 : -1.0) * std::min(std::fabs(qy1),(qy1 > 0)? cq : cqy1);
        qy2 = ((qy2 > 0)? 1.0 : -1.0) * std::min(std::fabs(qy2),(qy2 > 0)? cqy2 : cq);

        qx1 = UF_OUTORMV(_dem,r,c+1)? UF_BoundaryFlux2D(dtx1,_dx,_dx,_f->Drc,0,_fu->Drc,_fv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc,0.1 + N->Drc, 0,1) : qx1;
        qx2 = UF_OUTORMV(_dem,r,c-1)? -UF_BoundaryFlux2D(dtx2,_dx,_dx,_f->Drc,0,_fu->Drc,_fv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc,0.1 + N->Drc, 0,-1) : qx2;
        qy1 = UF_OUTORMV(_dem,r+1,c)? UF_BoundaryFlux2D(dty1,_dx,_dx,_f->Drc,0,_fu->Drc,_fv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc,0.1 + N->Drc, 1,0) : qy1;
        qy2 = UF_OUTORMV(_dem,r-1,c)? -UF_BoundaryFlux2D(dty2,_dx,_dx,_f->Drc,0,_fu->Drc,_fv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc,0.1 + N->Drc, -1,0) : qy2;

        UF2D_fqx1->Drc = qx1;
        UF2D_fqx2->Drc = qx2;
        UF2D_fqy1->Drc = qy1;
        UF2D_fqy2->Drc = qy2;

    }}}

}

void TWorld::UF2D_SolidMomentum2Source(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv)
{

    FOR_ROW_COL_UF2D
    {
        UF_t1->Drc = (_s->Drc)/(_dx * _dx);
    }

    UF2D_MUSCLE(_dem,dt,UF_t1,UF_MUSCLE_TARGET_IN1);
    UF2D_MUSCLE(_dem,dt,_su,UF_MUSCLE_TARGET_IN2, UF_DIRECTION_X);
    UF2D_MUSCLE(_dem,dt,_sv,UF_MUSCLE_TARGET_IN3, UF_DIRECTION_Y);

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

        UF_t2->Drc = UF_Gravity * 1 *hs * hs * UF_Aspect/2.0;
        UF_t3->Drc = UF_Gravity * 1 *hs * hs * UF_Aspect/2.0;

        UF_t4->Drc = _f->Drc/(_dx * _dx);
        UF_t5->Drc = _s->Drc/(_dx * _dx);
        UF2D_sax->Drc = 0;
        UF2D_say->Drc = 0;
    }

    FOR_ROW_COL_UF2D_DT
    {
        if(!(_s->Drc > UF_VERY_SMALL) ||!(_s->Drc + _f->Drc > UF_VERY_SMALL) )
        {
            UF2D_sax->Drc = 0;
            UF2D_say->Drc = 0;
            UF2D_sqx1->Drc = 0;
            UF2D_sqx2->Drc = 0;
            UF2D_sqy1->Drc = 0;
            UF2D_sqy2->Drc = 0;
            UF2D_sax1->Drc = 0;
            UF2D_sax2->Drc = 0;
            UF2D_say1->Drc = 0;
            UF2D_say2->Drc = 0;
            continue;
        }

        ///////////input parameters
        double lsu = UF2D_MUSCLE_2_x2->Drc;
        double lsv = UF2D_MUSCLE_3_y2->Drc;
        double rsu = UF2D_MUSCLE_2_x1->Drc;
        double rsv = UF2D_MUSCLE_3_y1->Drc;
        double lxh = std::max(0.0,UF2D_MUSCLE_1_x2->Drc);
        double rxh = std::max(0.0,UF2D_MUSCLE_1_x1->Drc);
        double lyh = std::max(0.0,UF2D_MUSCLE_1_y2->Drc);
        double ryh = std::max(0.0,UF2D_MUSCLE_1_y1->Drc);

        double h = (_s->Drc)/(_dx * _dx);
        double ff = _f->Drc / (_f->Drc + _s->Drc);
        double sf = _s->Drc / (_f->Drc + _s->Drc);
        double gamma = _d->Drc > UF_VERY_SMALL? 1000.0/_d->Drc : 0.5;
        double dc = UF_DragCoefficient(ff,sf,gamma,_visc->Drc, _rocksize->Drc, _d->Drc);
        double pbf = -UF_Gravity*h;
        double lxpbf = -UF_Gravity*lxh;
        double rxpbf = -UF_Gravity*rxh;
        double lypbf = -UF_Gravity*lyh;
        double rypbf = -UF_Gravity*ryh;

        double pbs = (1-gamma)*pbf;
        double ifa = 0.3;

        double lxslope = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double lyslope = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rxslope = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double ryslope = UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double dhsdx = UF2D_Derivative(_dem,UF_t5,r,c,UF_DIRECTION_X);
        double dhsdy = UF2D_Derivative(_dem,UF_t5,r,c,UF_DIRECTION_Y);
        double ldhsdx = UF2D_Derivative(_dem,UF_t5,r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double ldhsdy = UF2D_Derivative(_dem,UF_t5,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rdhsdx = UF2D_Derivative(_dem,UF_t5,r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rdhsdy = UF2D_Derivative(_dem,UF_t5,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double dhdx = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_X);
        double dhdy = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_Y);
        double ldhdx = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double ldhdy = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rdhdx = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rdhdy = UF2D_Derivative(_dem,UF_t1,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);
        double dbdx = UF2D_Derivative(_dem,UF_t2,r,c,UF_DIRECTION_X);
        double dbdy = UF2D_Derivative(_dem,UF_t3,r,c,UF_DIRECTION_Y);
        double ldbdx = UF2D_Derivative(_dem,UF_t2,r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double ldbdy = UF2D_Derivative(_dem,UF_t3,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rdbdx = UF2D_Derivative(_dem,UF_t2,r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rdbdy = UF2D_Derivative(_dem,UF_t3,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double lsax = UF2D_MomentumBalanceSolid(true,_f->Drc,lxh*(_dx*_dx),_fu->Drc, _fv->Drc, lsu, _sv->Drc, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs,lxpbf, lxslope, UF2D_SlopeY->Drc,
                                                  ldhsdx, dhsdy, ldhdx, dhdy, ldbdx, dbdy);

        double rsax = UF2D_MomentumBalanceSolid(true,_f->Drc,rxh*(_dx*_dx),_fu->Drc, _fv->Drc, rsu, _sv->Drc, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs,rxpbf, rxslope, UF2D_SlopeY->Drc,
                                                 rdhsdx, dhsdy, rdhdx, dhdy, rdbdx, dbdy);

        double lsay = UF2D_MomentumBalanceSolid(false,_f->Drc,lyh*(_dx*_dx),_fu->Drc, _fv->Drc, _su->Drc, lsv, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs,lypbf, UF2D_SlopeX->Drc, lyslope,
                                                  dhsdx, ldhsdy, dhdx, ldhdy, dbdx, ldbdy);

        double rsay = UF2D_MomentumBalanceSolid(false,_f->Drc,ryh*(_dx*_dx),_fu->Drc, _fv->Drc, _su->Drc, rsv, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs,rypbf, UF2D_SlopeX->Drc, ryslope,
                                                 dhsdx, rdhsdy, dhdx, rdhdy, dbdx, rdbdy);

        UF2D_sax1->Drc = rsax;
        UF2D_sax2->Drc = lsax;
        UF2D_say1->Drc = rsay;
        UF2D_say2->Drc = lsay;

        ////////////friction and actual accaleration

            UF_t6->Drc = lsu;
            UF_t7->Drc = rsu;

            UF_t6->Drc += dt->Drc*(2.0/3.0)* lsax;
            UF_t7->Drc += dt->Drc*(2.0/3.0)* rsax;
            UF_t6->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t6->Drc,0,N->Drc,lxh,lxslope,true);
            UF_t7->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t7->Drc,0,N->Drc,rxh,rxslope,true);
            double ut1 = UF_t6->Drc + dt->Drc*(2.0/3.0)* lsax;
            double ut2 = UF_t7->Drc + dt->Drc*(2.0/3.0)* rsax;
            UF_t6->Drc = (UF_t6->Drc + UF_Friction(dt->Drc*(2.0/3.0),ut1,0,N->Drc,lxh,lxslope,true))/2.0;
            UF_t7->Drc = (UF_t7->Drc + UF_Friction(dt->Drc*(2.0/3.0),ut2,0,N->Drc,rxh,rxslope,true))/2.0;

            if(ff > UF_VERY_SMALL)
            {
                double lvpow = pow((UF_t6->Drc-_fu->Drc)*(UF_t6->Drc-_fu->Drc)+(_fv->Drc-_sv->Drc)*(_fv->Drc-_sv->Drc),0.5*(UF_j-1.0));
                double rvpow = pow((UF_t7->Drc-_fu->Drc)*(UF_t7->Drc-_fu->Drc)+(_fv->Drc-_sv->Drc)*(_fv->Drc-_sv->Drc),0.5*(UF_j-1.0));
                double lsacu = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * ff *lvpow)));
                double rsacu = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * ff *rvpow)));
                UF_t6->Drc = (1.0-lsacu) * UF_t6->Drc + (lsacu) * _fu->Drc;
                UF_t7->Drc = (1.0-rsacu) * UF_t7->Drc + (rsacu) * _fu->Drc;
            }
            lsax = (UF_t6->Drc - lsu)/dt->Drc;
            rsax = (UF_t7->Drc - rsu)/dt->Drc;

            UF_t6->Drc = lsv;
            UF_t7->Drc = rsv;

            UF_t6->Drc += dt->Drc*(2.0/3.0)* lsay;
            UF_t7->Drc += dt->Drc*(2.0/3.0)* rsay;
            UF_t6->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t6->Drc,0,N->Drc,lyh,lyslope,true);
            UF_t7->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t7->Drc,0,N->Drc,ryh,ryslope,true);
            double vt1 = UF_t6->Drc + dt->Drc*(2.0/3.0)* lsay;
            double vt2 = UF_t7->Drc + dt->Drc*(2.0/3.0)* rsay;
            UF_t6->Drc = (UF_t6->Drc + UF_Friction(dt->Drc*(2.0/3.0),vt1,0,N->Drc,lyh,lyslope,true))/2.0;
            UF_t7->Drc = (UF_t7->Drc + UF_Friction(dt->Drc*(2.0/3.0),vt2,0,N->Drc,ryh,ryslope,true))/2.0;

            if(ff > UF_VERY_SMALL)
            {

                double lvpow = pow((UF_t6->Drc-_fv->Drc)*(UF_t6->Drc-_fv->Drc)+(_fu->Drc-_su->Drc)*(_fu->Drc-_su->Drc),0.5*(UF_j-1.0));
                double rvpow = pow((UF_t7->Drc-_fv->Drc)*(UF_t7->Drc-_fv->Drc)+(_fu->Drc-_su->Drc)*(_fu->Drc-_su->Drc),0.5*(UF_j-1.0));
                double lsacv = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * ff *lvpow)));
                double rsacv = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * ff *rvpow)));

                UF_t6->Drc = (1.0-lsacv) * UF_t6->Drc + (lsacv) * _fv->Drc;
                UF_t7->Drc = (1.0-rsacv) * UF_t7->Drc + (rsacv) * _fv->Drc;
            }
            lsay = (UF_t6->Drc - lsv)/dt->Drc;
            rsay = (UF_t7->Drc - rsv)/dt->Drc;

        /*UF2D_sax1->Drc = rsax;
        UF2D_sax2->Drc = lsax;
        UF2D_say1->Drc = rsay;
        UF2D_say2->Drc = lsay;*/

        ////////////average accalerations for cell centers
        UF2D_sax->Drc = (lsax + rsax)/2.0;
        UF2D_say->Drc = (lsay + rsay)/2.0;


        ////////////fluxes
        double hdem = (_s->Drc/(_dx*_dx) + _dem->Drc);
        double volx1 = UF_OUTORMV(_dem,r,c+1)?0.0:std::max(0.0,(_dx*_dx) * (hdem - (_s->data[r][c+1]/(_dx*_dx) + _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1))));
        double volx2 = UF_OUTORMV(_dem,r,c-1)?0.0:std::max(0.0,(_dx*_dx) * (hdem - (_s->data[r][c-1]/(_dx*_dx) + _dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1))));
        double voly1 = UF_OUTORMV(_dem,r+1,c)?0.0:std::max(0.0,(_dx*_dx) * (hdem - (_s->data[r+1][c]/(_dx*_dx) + _dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0))));
        double voly2 = UF_OUTORMV(_dem,r-1,c)?0.0:std::max(0.0,(_dx*_dx) * (hdem - (_s->data[r-1][c]/(_dx*_dx) + _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0))));

        double dtx1 = dt->Drc;//(UF2D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
        double dtx2 = dt->Drc;//(UF2D_MUSCLE_OUT_x2->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c-1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c-1)? 0.0: dt->data[r][c-1]));
        double dty1 = dt->Drc;//(UF2D_MUSCLE_OUT_y1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r+1,c)? dt->Drc : (UF_NOTIME(_dem,dt,r+1,c)? 0.0: dt->data[r+1][c]));
        double dty2 = dt->Drc;//(UF2D_MUSCLE_OUT_y2->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r-1,c)? dt->Drc : (UF_NOTIME(_dem,dt,r-1,c)? 0.0: dt->data[r-1][c]));

        double vx1 = rsu + dtx1 * rsax;
        double vx2 = lsu + dtx2 * lsax;
        double vy1 = rsv + dty1 * rsay;
        double vy2 = lsv + dty2 * lsay;

        rxh = UF_OUTORMV(_dem,r,c+1)?0.0:std::min(rxh,std::max(0.0,(vx1 > 0? 1.0:0.0)*((_dem->Drc + rxh) - (_s->data[r][c+1]/(_dx*_dx) + _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1)))));
        lxh = UF_OUTORMV(_dem,r,c-1)?0.0:std::min(lxh,std::max(0.0,(vx2 > 0? 0.0:1.0)*((_dem->Drc + lxh) - (_s->data[r][c-1]/(_dx*_dx) + _dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1)))));
        ryh = UF_OUTORMV(_dem,r+1,c)?0.0:std::min(ryh,std::max(0.0,(vy1 > 0? 1.0:0.0)*((_dem->Drc + ryh) - (_s->data[r+1][c]/(_dx*_dx) + _dem->data[r+1][c] +GetFlowBarrierHeight(r,c,-1,0)))));
        lyh = UF_OUTORMV(_dem,r-1,c)?0.0:std::min(lyh,std::max(0.0,(vy2 > 0? 0.0:1.0)*((_dem->Drc + lyh) - (_s->data[r-1][c]/(_dx*_dx) + _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0)))));

        double cq = 0.25 * UF_Courant * _s->Drc;
        double cqx1 = 0.25 * UF_Courant * (UF_OUTORMV(_dem,r,c+1)? 0.0 : volx1);
        double cqx2 = 0.25 * UF_Courant * (UF_OUTORMV(_dem,r,c-1)? 0.0 : volx2);
        double cqy1 = 0.25 * UF_Courant * (UF_OUTORMV(_dem,r+1,c)? 0.0 : voly1);
        double cqy2 = 0.25 * UF_Courant * (UF_OUTORMV(_dem,r-1,c)? 0.0 : voly2);

        double qx1 = 0.5 *dtx1 * (vx1) * rxh *_dx;
        double qx2 = 0.5 *dtx2 * (vx2) * lxh *_dx;
        double qy1 = 0.5 *dty1 * (vy1) * ryh *_dx;
        double qy2 = 0.5 *dty2 * (vy2) * lyh *_dx;

        qx1 = ((qx1 > 0)? 1.0 : -1.0) * std::min(std::fabs(qx1),(qx1 > 0)? cq : cqx1);
        qx2 = ((qx2 > 0)? 1.0 : -1.0) * std::min(std::fabs(qx2),(qx2 > 0)? cqx2 : cq);
        qy1 = ((qy1 > 0)? 1.0 : -1.0) * std::min(std::fabs(qy1),(qy1 > 0)? cq : cqy1);
        qy2 = ((qy2 > 0)? 1.0 : -1.0) * std::min(std::fabs(qy2),(qy2 > 0)? cqy2 : cq);

        qx1 = UF_OUTORMV(_dem,r,c+1)? UF_BoundaryFlux2D(dtx1,_dx,_dx,0,_s->Drc,_fu->Drc,_fv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc,0.1 + N->Drc, 0,1) : qx1;
        qx2 = UF_OUTORMV(_dem,r,c-1)? -UF_BoundaryFlux2D(dtx2,_dx,_dx,0,_s->Drc,_fu->Drc,_fv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc,0.1 + N->Drc, 0,-1) : qx2;
        qy1 = UF_OUTORMV(_dem,r+1,c)? UF_BoundaryFlux2D(dty1,_dx,_dx,0,_s->Drc,_fu->Drc,_fv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc,0.1 + N->Drc, 1,0) : qy1;
        qy2 = UF_OUTORMV(_dem,r-1,c)? -UF_BoundaryFlux2D(dty2,_dx,_dx,0,_s->Drc,_fu->Drc,_fv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc,0.1 + N->Drc, -1,0) : qy2;

        UF2D_sqx1->Drc = qx1;
        UF2D_sqx2->Drc = qx2;
        UF2D_sqy1->Drc = qy1;
        UF2D_sqy2->Drc = qy2;
    }}}

}


void TWorld::UF1D_FluidMomentum2Source(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu)
{

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

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

        UF_t4->Drc = _f->Drc/(_dx * _lddw->Drc);
        UF_t5->Drc = _s->Drc/(_dx * _lddw->Drc);
        UF1D_fq1->Drc = 0;
        UF1D_fq2->Drc = 0;
        UF1D_fa->Drc = 0;
        UF1D_fa1->Drc = 0;
        UF1D_fa2->Drc = 0;
    }

    UF1D_MUSCLE(_ldd,_lddw,dt,UF_t4,UF_MUSCLE_TARGET_IN1);
    UF1D_MUSCLE(_ldd,_lddw,dt,_fu,UF_MUSCLE_TARGET_IN2);

    FOR_ROW_COL_UF1D_DT
    {
        if(!(_f->Drc > UF_VERY_SMALL))
        {
            UF1D_fa->Drc = 0;
            continue;
        }
        double h = (_f->Drc + _s->Drc)/(_dx *_lddw->Drc);
        double rhf = std::max(0.0,(UF1D_MUSCLE_1_x1->Drc));
        double lhf = std::max(0.0,(UF1D_MUSCLE_1_x2->Drc));

        double rfu = (UF1D_MUSCLE_2_x1->Drc);
        double lfu = (UF1D_MUSCLE_2_x2->Drc);
        double ff = UF_t2->Drc;
        double sf = UF_t3->Drc;
        _visc->Drc = UF_DynamicViscosity(sf + (SwitchErosion? (((UF1D_blm->Drc + UF1D_ssm->Drc)/2000.0 + ff) > UF_VERY_SMALL ?((UF1D_blm->Drc + UF1D_ssm->Drc)/2000.0)/((UF1D_blm->Drc + UF1D_ssm->Drc)/2000.0 + ff) :0.0):0.0));
        _d->Drc = 2000.0;
        double Nr = UF_Reynolds(_d->Drc,_visc->Drc,ff,sf, _rocksize->Drc);
        UF1D_Nr->Drc = Nr;
        double Nra = 15000.0;
        double gamma =(!(_d->Drc > UF_VERY_SMALL))? 0.5 : 1000.0/_d->Drc;
        double dc = UF_DragCoefficient(ff, sf, gamma,_visc->Drc, _rocksize->Drc, _d->Drc);
        double pbf = UF_Gravity*h;
        double lpbf = UF_Gravity*lhf;
        double rpbf = UF_Gravity*rhf;

        double dhfdx = UF1D_Derivative(_ldd,_lddw,UF_t4,r,c);
        double ldhfdx = UF1D_Derivative(_ldd,_lddw,UF_t4,r,c,false,UF_DERIVATIVE_L);
        double rdhfdx = UF1D_Derivative(_ldd,_lddw,UF_t4,r,c,false,UF_DERIVATIVE_R);

        double dh2pbdx = UF1D_Derivative(_ldd,_lddw,UF_t1,r,c);
        double ldh2pbdx = UF1D_Derivative(_ldd,_lddw,UF_t1,r,c,false,UF_DERIVATIVE_L);
        double rdh2pbdx = UF1D_Derivative(_ldd,_lddw,UF_t1,r,c,false,UF_DERIVATIVE_R);

        double dsfdx = UF1D_Derivative(_ldd,_lddw,UF_t3,r,c);
        double ddsfdxx = UF1D_Derivative2(_ldd,_lddw,UF_t3,r,c);
        double dfudx = UF1D_Derivative(_ldd,_lddw,_fu,r,c);
        double ddfudxx = UF1D_Derivative2(_ldd,_lddw,_fu,r,c);
        double dsudx = UF1D_Derivative(_ldd,_lddw,_su,r,c);

        double lfa = UF1D_MomentumBalanceFluid( lhf*(_dx*_lddw->Drc),_s->Drc,lfu, _su->Drc, ff, sf, Nr, Nra, _ifa->Drc, gamma, _visc->Drc, lpbf, UF1D_Slope->Drc,
                                                 ldhfdx, ldh2pbdx, dsfdx, ddsfdxx, dfudx, ddfudxx, dsudx);
        double rfa = UF1D_MomentumBalanceFluid(  rhf*(_dx*_lddw->Drc),_s->Drc,rfu, _su->Drc, ff, sf, Nr, Nra, _ifa->Drc, gamma, _visc->Drc, rpbf, UF1D_Slope->Drc,
                                                 rdhfdx, rdh2pbdx, dsfdx, ddsfdxx, dfudx, ddfudxx, dsudx);

        UF1D_fa2->Drc = lfa;
        UF1D_fa1->Drc = rfa;

        UF_t6->Drc = lfu;
        UF_t7->Drc = rfu;

        UF_t6->Drc += dt->Drc*(2.0/3.0)* lfa;
        UF_t7->Drc += dt->Drc*(2.0/3.0)* rfa;
        UF_t6->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t6->Drc,0,N->Drc,lhf,UF1D_Slope->Drc, false);
        UF_t7->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t7->Drc,0,N->Drc,rhf,UF1D_Slope->Drc, false);
        double vt1 = UF_t6->Drc + dt->Drc*(2.0/3.0)* lfa;
        double vt2 = UF_t7->Drc + dt->Drc*(2.0/3.0)* rfa;
        UF_t6->Drc = (UF_t6->Drc + UF_Friction(dt->Drc*(2.0/3.0),vt1,0,N->Drc,lhf,UF1D_Slope->Drc, false))/2.0;
        UF_t7->Drc = (UF_t7->Drc + UF_Friction(dt->Drc*(2.0/3.0),vt2,0,N->Drc,rhf,UF1D_Slope->Drc, false))/2.0;

        if(sf > UF_VERY_SMALL)
        {
            double lupow = pow((UF_t6->Drc-_su->Drc)*(UF_t6->Drc-_su->Drc),0.5*(UF_j-1.0));
            double rupow = pow((UF_t7->Drc-_su->Drc)*(UF_t7->Drc-_su->Drc),0.5*(UF_j-1.0));
            double lfacu = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * sf *lupow)));
            double rfacu = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * sf *rupow)));
            UF_t6->Drc = (1.0-lfacu) * UF_t6->Drc + (lfacu) * _su->Drc;
            UF_t7->Drc = (1.0-rfacu) * UF_t7->Drc + (rfacu) * _su->Drc;
        }

        lfa = (UF_t6->Drc - lfu)/dt->Drc;
        rfa = (UF_t7->Drc - rfu)/dt->Drc;

        UF1D_fa->Drc = (lfa + rfa) / 2.0;

        double cq = 0.25 * UF_Courant * _f->Drc;
        double hself = _f->Drc/(_dx*_lddw->Drc);

        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!(lddself == 5))
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){

                ////////////fluxes
                double volxr = UF_OUTORMV(_ldd,r,c+1)?0.0:std::max(0.0,(_dx*_lddw->Drc) * (hself - (_f->data[r2][c2]/(_dx*_lddw->data[r2][c2]) - UF1D_Slope->Drc * _dx)));
                double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
                double vxr = rfu + dtx1 * rfa;
                rhf = std::min(rhf,std::max(0.0,(vxr > 0? 1.0:-1.0)*((rhf) - (_f->data[r2][c2]/(_dx*_dx) + UF1D_Slope->Drc * _dx))));

                double cqt = 0.5 * UF_Courant * (volxr);

                double q1 = 0.5 *dtx1 * (vxr) * rhf *_lddw->Drc;
                q1 = ((q1 > 0)? 1.0 : 0.0) * std::min(std::fabs(q1),(q1 > 0)? cq : cqt);

                UF1D_fq1->Drc = q1;
            }else if(rfu > 0)
            {
                double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
                double q1 = UF_BoundaryFlux1D(dtx1,_lddw->Drc,_f->Drc/(_dx*_lddw->Drc),0,rfu,_su->Drc,UF1D_Slope->Drc,ChannelN->Drc, true);
                q1 = ((q1 > 0)? 1.0 : 0.0) * std::min(std::fabs(q1),(q1 > 0)? cq : 0.0);

                UF1D_fq1->Drc = q1;
            }
        }else if (rfu > 0)
        {
            double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
            double q1 = UF_BoundaryFlux1D(dtx1,_lddw->Drc,_f->Drc,0,rfu,_su->Drc,UF1D_Slope->Drc,ChannelN->Drc, true);
            q1 = ((q1 > 0)? 1.0 : 0.0) * std::min(std::fabs(q1),(q1 > 0)? cq : 0.0);

            UF1D_fq1->Drc = q1;
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
            UF1D_fq2->Drc = 0;
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

                        ////////////fluxes
                        double volxl = UF_OUTORMV(_ldd,r,c+1)?0.0:std::max(0.0,(_dx*_lddw->Drc) * (hself - (_f->data[r2][c2]/(_dx*_lddw->data[r2][c2]) + UF1D_Slope->Drc * _dx)));
                        double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
                        double vxl = lfu + dtx1 * lfa;
                        lhf = std::min(lhf,std::max(0.0,(vxl > 0? 1.0:-1.0)*((lhf) - (_f->data[r2][c2]/(_dx*_dx) - UF1D_Slope->Drc * _dx))));
                        double cqt = 0.5 * UF_Courant * (volxl);
                        double q2 = 0.5 *dtx1 * (vxl) * lhf *_lddw->Drc* (_lddw->data[r2][c2]/totalwidth);

                        q2 = ((q2 > 0)? 0.0 : -1.0) * std::min(std::fabs(q2),(q2 > 0)? cqt : cq);

                        UF1D_fq2->Drc += q2;
                    }
                }
            }
        }else if(lfu < 0)
        {
            double dtx2 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
            double q2 = UF_BoundaryFlux1D(dtx2,_lddw->Drc,_f->Drc,0,lfu,_su->Drc,UF1D_Slope->Drc,ChannelN->Drc, true);
            q2 = ((q2 < 0)? -1.0 : 0.0) * std::min(std::fabs(q2),(q2 > 0)? cq : 0.0);

            UF1D_fq2->Drc = q2;
        }
    }
}

void TWorld::UF1D_SolidMomentum2Source(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_su)
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

        UF_t4->Drc = _f->Drc/(_dx * _lddw->Drc);
        UF_t5->Drc = _s->Drc/(_dx * _lddw->Drc);
        UF1D_sq1->Drc = 0;
        UF1D_sq2->Drc = 0;
        UF1D_sa->Drc = 0;
        UF1D_sa1->Drc = 0;
        UF1D_sa2->Drc = 0;
    }

    UF1D_MUSCLE(_ldd,_lddw,dt,UF_t5,UF_MUSCLE_TARGET_IN1);
    UF1D_MUSCLE(_ldd,_lddw,dt,_su,UF_MUSCLE_TARGET_IN2);

    FOR_ROW_COL_UF1D_DT
    {
        if(!(_s->Drc > UF_VERY_SMALL))
        {
            UF1D_sa->Drc = 0;
            continue;
        }

        double h = (_f->Drc + _s->Drc)/(_dx * _lddw->Drc);
        double ff = _f->Drc/(_f->Drc + _s->Drc);
        double sf = _s->Drc/(_f->Drc + _s->Drc);
        double gamma = _d->Drc > UF_VERY_SMALL? 1000.0/_d->Drc : 0.5;
        double dc = UF_DragCoefficient(ff,sf,gamma,_visc->Drc, _rocksize->Drc, _d->Drc);

        double pbf = -UF_Gravity*h * ff;
        double pbs = (1.0-gamma)*pbf;
        double ifa = _ifa->Drc;

        double dhsdx = UF1D_Derivative(_ldd,_lddw,UF_t5,r,c);
        double ldhsdx = UF1D_Derivative(_ldd,_lddw,UF_t5,r,c, false,UF_DERIVATIVE_L);
        double rdhsdx = UF1D_Derivative(_ldd,_lddw,UF_t5,r,c, false, UF_DERIVATIVE_R);

        double rsu = (UF1D_MUSCLE_2_x1->Drc);
        double lsu = (UF1D_MUSCLE_2_x2->Drc);
        double rhs = std::max(0.0,(UF1D_MUSCLE_1_x1->Drc));
        double lhs = std::max(0.0,(UF1D_MUSCLE_1_x2->Drc));

        double dhdx = UF1D_Derivative(_ldd,_lddw,UF_t1,r,c);
        double dbdx = UF1D_Derivative(_ldd,_lddw,UF_t2,r,c);
        double ldhdx = UF1D_Derivative(_ldd,_lddw,UF_t1,r,c, false,UF_DERIVATIVE_L);
        double ldbdx = UF1D_Derivative(_ldd,_lddw,UF_t2,r,c, false,UF_DERIVATIVE_L);
        double rdhdx = UF1D_Derivative(_ldd,_lddw,UF_t1,r,c, false,UF_DERIVATIVE_R);
        double rdbdx = UF1D_Derivative(_ldd,_lddw,UF_t2,r,c, false,UF_DERIVATIVE_R);

        double lsa = UF1D_MomentumBalanceSolid(_f->Drc,lhs*(_dx*_lddw->Drc),_fu->Drc, lsu, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs, pbf, UF1D_Slope->Drc,
                                                 ldhsdx, ldhdx, ldbdx);
        double rsa = UF1D_MomentumBalanceSolid(_f->Drc,rhs*(_dx*_lddw->Drc),_fu->Drc, rsu, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs, pbf, UF1D_Slope->Drc,
                                                 rdhsdx, rdhdx, rdbdx);

        UF1D_sa2->Drc = lsa;
        UF1D_sa1->Drc = rsa;

        UF_t6->Drc = lsu;
        UF_t7->Drc = rsu;

        UF_t6->Drc += dt->Drc*(2.0/3.0)* lsa;
        UF_t7->Drc += dt->Drc*(2.0/3.0)* rsa;
        UF_t6->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t6->Drc,0,N->Drc,lhs,UF1D_Slope->Drc, true);
        UF_t7->Drc = UF_Friction(dt->Drc*(2.0/3.0),UF_t7->Drc,0,N->Drc,rhs,UF1D_Slope->Drc, true);
        double vt1 = UF_t6->Drc + dt->Drc*(2.0/3.0)* lsa;
        double vt2 = UF_t7->Drc + dt->Drc*(2.0/3.0)* rsa;
        UF_t6->Drc = (UF_t6->Drc + UF_Friction(dt->Drc*(2.0/3.0),vt1,0,N->Drc,lhs,UF1D_Slope->Drc, true))/2.0;
        UF_t7->Drc = (UF_t7->Drc + UF_Friction(dt->Drc*(2.0/3.0),vt2,0,N->Drc,rhs,UF1D_Slope->Drc, true))/2.0;

        if(sf > UF_VERY_SMALL)
        {
            double lupow = pow((UF_t6->Drc-_fu->Drc)*(UF_t6->Drc-_fu->Drc),0.5*(UF_j-1.0));
            double rupow = pow((UF_t7->Drc-_fu->Drc)*(UF_t7->Drc-_fu->Drc),0.5*(UF_j-1.0));
            double lsacu = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * ff *lupow)));
            double rsacu = std::min(0.25,std::max(0.0,std::fabs(dt->Drc * dc * ff *rupow)));
            UF_t6->Drc = (1.0-lsacu) * UF_t6->Drc + (lsacu) * _fu->Drc;
            UF_t7->Drc = (1.0-rsacu) * UF_t7->Drc + (rsacu) * _fu->Drc;
        }

        lsa = (UF_t6->Drc - lsu)/dt->Drc;
        rsa = (UF_t7->Drc - rsu)/dt->Drc;

        UF1D_sa->Drc = (lsa + rsa) / 2.0;

        double cq = 0.25 * UF_Courant * _s->Drc;

        double hself = _s->Drc/(_dx * _lddw->Drc);
        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!(lddself == 5))
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){

                ////////////fluxes
                double volxr = UF_OUTORMV(_ldd,r,c+1)?0.0:std::max(0.0,(_dx*_lddw->Drc) * (hself - (_s->data[r2][c2]/(_dx*_lddw->data[r2][c2]) - UF1D_Slope->Drc * _dx)));
                double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
                double vxr = rsu + dtx1 * rsa;
                rhs = std::min(rhs,std::max(0.0,(vxr > 0? 1.0:-1.0)*((rhs) - (_s->data[r2][c2]/(_dx*_dx) + UF1D_Slope->Drc * _dx))));

                double cqt = 0.5 * UF_Courant * (volxr);

                double q1 = 0.5 *dtx1 * (vxr) * rhs *_lddw->Drc;
                q1 = ((q1 > 0)? 1.0 : 0.0) * std::min(std::fabs(q1),(q1 > 0)? cq : cqt);

                UF1D_sq1->Drc = q1;
            }else if(_su->Drc > 0)
            {
                double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
                double q1 = UF_BoundaryFlux1D(dtx1,_lddw->Drc,_s->Drc,0,rsu,_fu->Drc,UF1D_Slope->Drc,0.1 + N->Drc, true);
                q1 = ((q1 > 0)? 1.0 : 0.0) * std::min(std::fabs(q1),(q1 > 0)? cq : 0.0);

                UF1D_sq1->Drc = q1;
            }
        }else
        {
            if(_su->Drc > 0)
            {
                double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
                double q1 = UF_BoundaryFlux1D(dtx1,_lddw->Drc,_s->Drc,0,rsu,_fu->Drc,UF1D_Slope->Drc,0.1 + N->Drc, true);
                q1 = ((q1 > 0)? 1.0 : 0.0) * std::min(std::fabs(q1),(q1 > 0)? cq : 0.0);


                UF1D_sq1->Drc = q1;
            }
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

                        ////////////fluxes
                        double volxl = UF_OUTORMV(_ldd,r,c+1)?0.0:std::max(0.0,(_dx*_lddw->Drc) * (hself - (_s->data[r2][c2]/(_dx*_lddw->data[r2][c2]) + UF1D_Slope->Drc * _dx)));
                        double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
                        double vxl = lsu + dtx1 * lsa;
                        lhs = std::min(lhs,std::max(0.0,(vxl > 0? 1.0:-1.0)*((lhs) - (_s->data[r2][c2]/(_dx*_dx) - UF1D_Slope->Drc * _dx))));
                        double cqt = 0.5 * UF_Courant * (volxl);
                        double q2 = 0.5 *dtx1 * (vxl) * lhs *_lddw->Drc* (_lddw->data[r2][c2]/totalwidth);

                        q2 = ((q2 > 0)? 0.0 :-1.0) * std::min(std::fabs(q2),(q2 > 0)? cqt : cq);

                        UF1D_sq2->Drc += q2;
                    }
                }
            }
        }else if(_su->Drc < 0 )
        {
            double dtx2 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
            double q2 = UF_BoundaryFlux1D(dtx2,_lddw->Drc,_s->Drc,0,lsu,_su->Drc,UF1D_Slope->Drc,0.1 + N->Drc, true);
            q2 = ((q2 < 0)? -1.0 : 0.0) * std::min(std::fabs(q2),(q2 < 0)? cq : 0.0);

            UF1D_sq2->Drc = q2;
        }
    }

}


