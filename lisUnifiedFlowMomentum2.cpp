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

#define UF_TIMERATIO (2.0/3.0)

void TWorld::UF2D_FluidApplyMomentum2(int thread,cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv)
{
    FOR_ROW_COL_UF2DMTDER
    {
        out_fu->Drc = _fu->Drc;
        out_fv->Drc = _fv->Drc;
    }}}

    FOR_ROW_COL_UF2DMT_DT
    {

        if(!(_f->Drc > UF_VERY_SMALL))
        {
            out_fu->Drc = 0;
            out_fv->Drc = 0;
            continue;
        }

        out_fu->Drc =  out_fu->Drc+dt->Drc * UF2D_fax->Drc;
        out_fv->Drc =  out_fv->Drc+dt->Drc * UF2D_fay->Drc;

    }}}

}

void TWorld::UF2D_SolidApplyMomentum2(int thread,cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv)
{
    FOR_ROW_COL_UF2DMTDER
    {
        out_su->Drc = _su->Drc;
        out_sv->Drc = _sv->Drc;
    }}}

    FOR_ROW_COL_UF2DMT_DT
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

void TWorld::UF1D_FluidApplyMomentum2(int thread,cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu)
{
    FOR_ROW_COL_UF1DMTDER
    {
        out_fu->Drc = _fu->Drc;
    }}}

    FOR_ROW_COL_UF1DMT_DT
    {
        if(!(_f->Drc > UF_VERY_SMALL))
        {
            out_fu->Drc = 0;
            continue;
        }
        out_fu->Drc += dt->Drc * UF1D_fa->Drc;


    }}}
}


void TWorld::UF1D_SolidApplyMomentum2(int thread,cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_su)
{
    FOR_ROW_COL_UF1DMTDER
    {
        out_su->Drc = _su->Drc;
    }}}

    FOR_ROW_COL_UF1DMT_DT
    {
        if(!(_s->Drc > UF_VERY_SMALL))
        {
            out_su->Drc = 0;
            continue;
        }
        out_su->Drc += dt->Drc * UF1D_sa->Drc;
    }}}
}


void TWorld::UF2D_FluidMomentum2Source(int thread,cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv)
{
    FOR_ROW_COL_UF2DMTDER
    {
        ThreadPool->UF_t1.at(thread)->Drc = (_f->Drc)/(_dx * _dx);
    }}}

    UF2D_MUSCLE(thread,_dem,dt,ThreadPool->UF_t1.at(thread),UF_MUSCLE_TARGET_IN1);
    UF2D_MUSCLE(thread,_dem,dt,_fu,UF_MUSCLE_TARGET_IN2, UF_DIRECTION_X);
    UF2D_MUSCLE(thread,_dem,dt,_fv,UF_MUSCLE_TARGET_IN3, UF_DIRECTION_Y);

    FOR_ROW_COL_UF2DMTDER
    {
        double s = 0;
        if(UF_SOLIDPHASE)
        {
            s = _s->Drc;
        }

        double h = (s + _f->Drc)/(_dx * _dx);
        ThreadPool->UF_t1.at(thread)->Drc = h*(UF_Gravity*h)/2.0;

        if(_f->Drc + s> UF_VERY_SMALL)
        {
            ThreadPool->UF_t2.at(thread)->Drc = _f->Drc/(_f->Drc + s);
            ThreadPool->UF_t3.at(thread)->Drc = s/(_f->Drc + s);
        }else
        {
            ThreadPool->UF_t2.at(thread)->Drc = 0;
            ThreadPool->UF_t3.at(thread)->Drc = 0;
        }

        ThreadPool->UF_t4.at(thread)->Drc = _f->Drc/(_dx * _dx);
        ThreadPool->UF_t5.at(thread)->Drc = s/(_dx * _dx);

        ThreadPool->UF_t6.at(thread)->Drc = (s+_f->Drc)/(_dx * _dx);

        ThreadPool->UF_t8.at(thread)->Drc = UF2D_DEM->Drc + (s+_f->Drc)/(_dx * _dx);

        UF2D_fax->Drc = 0;
        UF2D_fay->Drc = 0;
    }}}

    FOR_ROW_COL_UF2DMT_DT
    {

        double s = 0;
        if(UF_SOLIDPHASE)
        {
            s = _s->Drc;
        }

        if(!(_f->Drc > UF_VERY_SMALL) ||!(s + _f->Drc > UF_VERY_SMALL) )
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

        if(!(dt->Drc > 1e-12))
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


        double ff = _f->Drc/(_f->Drc + s);
        double sf = s/(_f->Drc + s);
        _visc->Drc = std::max(1.0,UF_DynamicViscosity(sf + ((SwitchErosion && UF_SUSPENDEDVISCOSITY)? (((UF2D_blm->Drc + UF2D_ssm->Drc)/2000.0 + ff) > UF_VERY_SMALL ?((UF2D_blm->Drc + UF2D_ssm->Drc)/2000.0)/((UF2D_blm->Drc + UF2D_ssm->Drc)/2000.0 + ff):0.0):0.0)));
        double Nr = 1.0;
        if(UF_SOLIDPHASE)
        {
            Nr = std::min(10000.0,std::max(1.0,UF_Reynolds(_d->Drc,_visc->Drc,ff,s, _rocksize->Drc)));

            UF2D_Nr->Drc = Nr;
        }

        double Nra = UF_NRA;
        double gamma = 1.0;
        if(UF_SOLIDPHASE)
        {
            gamma = _d->Drc > UF_VERY_SMALL? 1000.0/_d->Drc : 0.5;
        }
        double dc = 0.0;

        if(UF_SOLIDPHASE)
        {
            dc = std::min(1000.0,std::max(0.0,UF_DragCoefficient(ff,s,gamma,_visc->Drc, _rocksize->Drc, _d->Drc)));
            UF2D_DC->Drc = dc;
        }


        double lxpbf = UF_Gravity*lxh;
        double rxpbf = UF_Gravity*rxh;
        double lypbf = UF_Gravity*lyh;
        double rypbf = UF_Gravity*ryh;

        double dalx = UF_DEMACCES(_dem,ThreadPool->UF_t6.at(thread),r,c,0,-1);
        double darx = UF_DEMACCES(_dem,ThreadPool->UF_t6.at(thread),r,c,0,1);
        double daly = UF_DEMACCES(_dem,ThreadPool->UF_t6.at(thread),r,c,-1,0);
        double dary = UF_DEMACCES(_dem,ThreadPool->UF_t6.at(thread),r,c,1,0);
        double dax = (dalx + darx ) * 0.5;
        double day = (daly + dary ) * 0.5;


        double dhfdx = dax * UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_X);
        double dhfdy = day * UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_Y);
        double ldhfdx = dalx * UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double ldhfdy = daly * UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rdhfdx = darx * UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rdhfdy = dary * UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double lddemhdx = dalx * UF2D_Derivative(_dem,ThreadPool->UF_t8.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double lddemhdy = daly * UF2D_Derivative(_dem,ThreadPool->UF_t8.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rddemhdx = darx * UF2D_Derivative(_dem,ThreadPool->UF_t8.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rddemhdy = dary * UF2D_Derivative(_dem,ThreadPool->UF_t8.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);


        /*double dhfdx = UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_X, UF_DERIVATIVE_LR, true);
        double dhfdy = UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_Y, UF_DERIVATIVE_LR,true);
        double ldhfdx = UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_L, true);
        double ldhfdy = UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L, true);
        double rdhfdx = UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_R, true);
        double rdhfdy = UF2D_Derivative(_dem,ThreadPool->UF_t4.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R, true);*/

        double dh2pbdx = dax * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_X);
        double dh2pbdy = day * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_Y);
        double ldh2pbdx = dalx * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double ldh2pbdy = daly * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rdh2pbdx = darx * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rdh2pbdy = dary * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double dsfdx = dax * UF2D_Derivative(_dem,ThreadPool->UF_t3.at(thread),r,c,UF_DIRECTION_X);
        double dsfdy = day * UF2D_Derivative(_dem,ThreadPool->UF_t3.at(thread),r,c,UF_DIRECTION_Y);
        double ddsfdxx = dax * UF2D_Derivative2(_dem,ThreadPool->UF_t3.at(thread),r,c,UF_DIRECTION_X);
        double ddsfdyy = day * UF2D_Derivative2(_dem,ThreadPool->UF_t3.at(thread),r,c,UF_DIRECTION_Y);
        double ddsfdxy = dax * day * UF2D_Derivative2(_dem,ThreadPool->UF_t3.at(thread),r,c,UF_DIRECTION_XY);

        double dfudx = dax * UF2D_Derivative(_dem,_fu,r,c,UF_DIRECTION_X,UF_DERIVATIVE_LR,std::max(_dx,1.0));
        double dfudy = day * UF2D_Derivative(_dem,_fu,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_LR,std::max(_dx,1.0));
        double dfvdx = dax * UF2D_Derivative(_dem,_fv,r,c,UF_DIRECTION_X,UF_DERIVATIVE_LR,std::max(_dx,1.0));
        double dfvdy = day * UF2D_Derivative(_dem,_fv,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_LR,std::max(_dx,1.0));
        double ddfudxx = dax * UF2D_Derivative2(_dem,_fu,r,c,UF_DIRECTION_X,UF_DERIVATIVE_LR,std::max(_dx,1.0));
        double ddfudyy = day * UF2D_Derivative2(_dem,_fu,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_LR,std::max(_dx,1.0));
        double ddfvdxy = dax * day * UF2D_Derivative2(_dem,_fv,r,c,UF_DIRECTION_XY,UF_DERIVATIVE_LR,std::max(_dx,1.0));
        double ddfvdxx = dax * UF2D_Derivative2(_dem,_fv,r,c,UF_DIRECTION_X,UF_DERIVATIVE_LR,std::max(_dx,1.0));
        double ddfvdyy = day * UF2D_Derivative2(_dem,_fv,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_LR,std::max(_dx,1.0));
        double ddfudxy = dax * day * UF2D_Derivative2(_dem,_fu,r,c,UF_DIRECTION_XY,UF_DERIVATIVE_LR,std::max(_dx,1.0));

        double dsudx = 0;
        double dsudy = 0;
        double dsvdx = 0;
        double dsvdy = 0;

        if(UF_SOLIDPHASE)
        {
            dsudx = dax * UF2D_Derivative(_dem,_su,r,c,UF_DIRECTION_X);
            dsudy = day * UF2D_Derivative(_dem,_su,r,c,UF_DIRECTION_Y);
            dsvdx = dax * UF2D_Derivative(_dem,_sv,r,c,UF_DIRECTION_X);
            dsvdy = day * UF2D_Derivative(_dem,_sv,r,c,UF_DIRECTION_Y);

        }

        UF2D_Test->Drc =  std::fabs(lddemhdx) + std::fabs(rddemhdx) + std::fabs(lddemhdy) + std::fabs(rddemhdy);


        double su = 0;
        double sv = 0;
        double ifa = 0;
        if(UF_SOLIDPHASE)
        {
            su = _su->Drc;
            sv = _sv->Drc;
            ifa = _ifa->Drc;
        }


        ////////////momentum balance
        double lfax = UF2D_MomentumBalanceFluid(true,lxh*_dx*_dx,s,lfu, _fv->Drc,su, sv, ff, sf, Nr, Nra, ifa, gamma, _visc->Drc, lxpbf, lxslope, 0,
                                                 ldhfdx,  dhfdy,   ldh2pbdx,   dh2pbdy,   dsfdx,   dsfdy,   ddsfdxx,   ddsfdyy,   ddsfdxy,
                                                 dfudx,   dfudy,   dfvdx,   dfvdy,   ddfudxx,   ddfudyy,   ddfvdxy,   ddfvdxx,   ddfvdyy,   ddfudxy,
                                                 dsudx,   dsudy,   dsvdx,  dsvdy, lddemhdx,0);

        double rfax = UF2D_MomentumBalanceFluid(true, rxh*_dx*_dx,s,rfu, _fv->Drc,su, sv, ff, sf, Nr, Nra, ifa, gamma, _visc->Drc, rxpbf, rxslope, 0,
                                                 rdhfdx,  dhfdy,   rdh2pbdx,   dh2pbdy,   dsfdx,   dsfdy,   ddsfdxx,   ddsfdyy,   ddsfdxy,
                                                 dfudx,   dfudy,   dfvdx,   dfvdy,   ddfudxx,   ddfudyy,   ddfvdxy,   ddfvdxx,   ddfvdyy,   ddfudxy,
                                                 dsudx,   dsudy,   dsvdx,  dsvdy,rddemhdx,0);

        double lfay = UF2D_MomentumBalanceFluid(false,lyh*_dx*_dx,s,_fu->Drc, lfv,su, sv, ff, sf, Nr, Nra, ifa, gamma, _visc->Drc, lypbf, 0, lyslope,
                                                 dhfdx,  ldhfdy,   dh2pbdx,   ldh2pbdy,   dsfdx,   dsfdy,   ddsfdxx,   ddsfdyy,   ddsfdxy,
                                                 dfudx,   dfudy,   dfvdx,   dfvdy,   ddfudxx,   ddfudyy,   ddfvdxy,   ddfvdxx,   ddfvdyy,   ddfudxy,
                                                 dsudx,   dsudy,   dsvdx,  dsvdy,0,lddemhdy);

        double rfay = UF2D_MomentumBalanceFluid(false,ryh*_dx*_dx,s,_fu->Drc, rfv,su, sv, ff, sf, Nr, Nra, ifa, gamma, _visc->Drc, rypbf, 0, ryslope,
                                                 dhfdx,  rdhfdy,   dh2pbdx,   rdh2pbdy,   dsfdx,   dsfdy,   ddsfdxx,   ddsfdyy,   ddsfdxy,
                                                 dfudx,   dfudy,   dfvdx,   dfvdy,   ddfudxx,   ddfudyy,   ddfvdxy,   ddfvdxx,   ddfvdyy,   ddfudxy,
                                                 dsudx,   dsudy,   dsvdx,  dsvdy,0,rddemhdy);


        /*rfax = -lddemhdx * 9.81;
        lfax = -rddemhdx * 9.81;
        rfay = -lddemhdy * 9.81;
        lfay = -rddemhdy * 9.81;*/

        /*rfax = 0.0;
        lfax = 0.0;
        rfay = 10.0;
        lfay = 10.0;*/

        /*rfax = -UF2D_SlopeX->Drc * 9.81;
        lfax = -UF2D_SlopeX->Drc * 9.81;
        rfay = -UF2D_SlopeY->Drc * 9.81;
        lfay = -UF2D_SlopeY->Drc * 9.81;*/

        UF2D_fax1->Drc = rfax;
        UF2D_fax2->Drc = lfax;
        UF2D_fay1->Drc = rfay;
        UF2D_fay2->Drc = lfay;

        UF2D_fax->Drc = 0.5*(lfax+rfax);
        UF2D_fay->Drc = 0.5*(lfay+rfay);


        ////////////friction and actual accaleration

            ThreadPool->UF_t6.at(thread)->Drc = lfu;
            ThreadPool->UF_t7.at(thread)->Drc = rfu;

            ThreadPool->UF_t6.at(thread)->Drc =  UF_Friction(dt->Drc*UF_TIMERATIO* lfax,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lxh,lxslope,false,false,0.0,0,ff,sf,Nr);
            ThreadPool->UF_t7.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* rfax,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,rxh,rxslope,false,false,0.0,0,ff,sf,Nr);


            ThreadPool->UF_t6.at(thread)->Drc = (ThreadPool->UF_t6.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* lfax,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lxh,lxslope,false))/2.0;
            ThreadPool->UF_t7.at(thread)->Drc = (ThreadPool->UF_t7.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* rfax,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,rxh,rxslope,false))/2.0;


            if(sf > UF_VERY_SMALL)
            {
                double lvpow = pow((ThreadPool->UF_t6.at(thread)->Drc-su)*(ThreadPool->UF_t6.at(thread)->Drc-su)+(_fv->Drc-sv)*(_fv->Drc-sv),0.5*(UF_j-1.0));
                double rvpow = pow((ThreadPool->UF_t7.at(thread)->Drc-su)*(ThreadPool->UF_t7.at(thread)->Drc-su)+(_fv->Drc-sv)*(_fv->Drc-sv),0.5*(UF_j-1.0));
                double lfacu = std::max(0.0,std::fabs(dt->Drc * dc * sf));
                double rfacu = std::max(0.0,std::fabs( dt->Drc * dc * sf));
                double ul_balance = ff * ThreadPool->UF_t6.at(thread)->Drc + sf * su;
                double ur_balance = ff * ThreadPool->UF_t7.at(thread)->Drc + sf * su;
                ThreadPool->UF_t6.at(thread)->Drc = ul_balance + (ThreadPool->UF_t6.at(thread)->Drc - ul_balance)*std::exp(-lfacu);
                ThreadPool->UF_t7.at(thread)->Drc = ur_balance + (ThreadPool->UF_t7.at(thread)->Drc - ur_balance)*std::exp(-rfacu);
            }

            lfax = (ThreadPool->UF_t6.at(thread)->Drc - lfu)/dt->Drc;
            rfax = (ThreadPool->UF_t7.at(thread)->Drc - rfu)/dt->Drc;

            ThreadPool->UF_t6.at(thread)->Drc = lfv;
            ThreadPool->UF_t7.at(thread)->Drc = rfv;

            ThreadPool->UF_t6.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* lfay,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lyh,lyslope,false,false,0.0,0,ff,sf,Nr);
            ThreadPool->UF_t7.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* rfay,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,ryh,ryslope,false,false,0.0,0,ff,sf,Nr);
            ThreadPool->UF_t6.at(thread)->Drc = (ThreadPool->UF_t6.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* lfay,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lyh,lyslope,false))/2.0;
            ThreadPool->UF_t7.at(thread)->Drc = (ThreadPool->UF_t7.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* rfay,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,ryh,ryslope,false))/2.0;


            if(sf > UF_VERY_SMALL)
            {
                double lvpow = pow((ThreadPool->UF_t6.at(thread)->Drc-sv)*(ThreadPool->UF_t6.at(thread)->Drc-sv)+(_fu->Drc-su)*(_fu->Drc-su),0.5*(UF_j-1.0));
                double rvpow = pow((ThreadPool->UF_t7.at(thread)->Drc-sv)*(ThreadPool->UF_t7.at(thread)->Drc-sv)+(_fu->Drc-su)*(_fu->Drc-su),0.5*(UF_j-1.0));
                double lfacu = std::max(0.0,std::fabs(dt->Drc * dc * sf));
                double rfacu = std::max(0.0,std::fabs(dt->Drc * dc * sf));
                double ul_balance = ff * ThreadPool->UF_t6.at(thread)->Drc + sf * sv;
                double ur_balance = ff * ThreadPool->UF_t7.at(thread)->Drc + sf * sv;
                ThreadPool->UF_t6.at(thread)->Drc = ul_balance + (ThreadPool->UF_t6.at(thread)->Drc - ul_balance)*std::exp(-lfacu);
                ThreadPool->UF_t7.at(thread)->Drc = ur_balance + (ThreadPool->UF_t7.at(thread)->Drc - ur_balance)*std::exp(-rfacu);
            }

            lfay = (ThreadPool->UF_t6.at(thread)->Drc - lfv)/dt->Drc;
            rfay = (ThreadPool->UF_t7.at(thread)->Drc - rfv)/dt->Drc;



        ////////////average accalerations for cell centers
        UF2D_fax->Drc = 0.5 * (lfax+rfax);
        UF2D_fay->Drc = 0.5 * (lfay+rfay);

        UF2D_fax1->Drc = rfax;//(std::fabs(UF2D_fax->Drc)/std::max(std::fabs(lfax),std::fabs(rfax)))*(rfax);
        UF2D_fax2->Drc = lfax;//(std::fabs(UF2D_fax->Drc)/std::max(std::fabs(lfax),std::fabs(rfax)))*(lfax);
        UF2D_fay1->Drc = rfay;//(std::fabs(UF2D_fay->Drc)/std::max(std::fabs(lfay),std::fabs(rfay)))*(rfay);
        UF2D_fay2->Drc = lfay;//(std::fabs(UF2D_fay->Drc)/std::max(std::fabs(lfay),std::fabs(rfay)))*(lfay);

        ////////////fluxes
        double hdem = ((s + _f->Drc)/(_dx*_dx) + _dem->Drc);
        double volx1;
        double volx2;
        double voly1;
        double voly2;

        double dtx1 = (UF_OUTORMV(_dem,r,c+1)? dt->Drc : std::max(dt->Drc ,dt->data[r][c+1]));
        double dtx2 = (UF_OUTORMV(_dem,r,c-1)? dt->Drc : std::max(dt->Drc ,dt->data[r][c-1]));
        double dty1 = (UF_OUTORMV(_dem,r+1,c)? dt->Drc : std::max(dt->Drc ,dt->data[r+1][c]));
        double dty2 = (UF_OUTORMV(_dem,r-1,c)? dt->Drc : std::max(dt->Drc ,dt->data[r-1][c]));

        double vx1 = rfu + dtx1 * rfax;
        double vx2 = lfu + dtx2 * lfax;
        double vy1 = rfv + dty1 * rfay;
        double vy2 = lfv + dty2 * lfay;

        double phi_rqx = (rfu + dt->Drc *  rfax) *(rxh);
        double phi_qx = (_fu->Drc + dt->Drc * UF2D_fax->Drc) *(_f->Drc)/(_dx*_dx);
        double phi_lqx = (lfu + dt->Drc *  lfax) *(lxh);

        double phi_rqy = (rfv + dt->Drc *  rfay) *(ryh);
        double phi_qy = (_fv->Drc + dt->Drc * UF2D_fay->Drc) *(_f->Drc)/(_dx*_dx);
        double phi_lqy = (lfv + dt->Drc *  lfay) *(lyh);

        double phi_rx = std::fabs(phi_qx - phi_lqx) > 0? (phi_rqx - phi_qx)/(phi_qx - phi_lqx) : 1.0;
        double phi_ry = std::fabs(phi_qy - phi_lqy) > 0? (phi_rqy - phi_qy)/(phi_qy - phi_lqy) : 1.0;

        double phi_rx_i = std::fabs(phi_rqx - phi_qx) > 0? (phi_qx - phi_lqx)/(phi_rqx - phi_qx) : 1.0;
        double phi_ry_i = std::fabs(phi_rqy - phi_qy) > 0? (phi_qy - phi_lqy)/(phi_rqy - phi_qy) : 1.0;

        double phi_x = std::max(0.0,std::min(1.0, std::min(phi_rx,phi_rx)));
        double phi_y = std::max(0.0,std::min(1.0, std::min(phi_ry,phi_ry)));

        //double phi_x = std::min(1.0,std::max(0.0,std::max(std::min(1.0, 2.0*std::max(phi_rx,phi_rx_i)),std::min(2.0,std::max(phi_rx,phi_rx_i)) )));
        //double phi_y = std::min(1.0,std::max(0.0,std::max(std::min(1.0, 2.0*std::max(phi_ry,phi_ry_i)),std::min(2.0,std::max(phi_ry,phi_ry_i)) )));

        if(!UF_USE_HLL2)
        {
            phi_x = 0;
            phi_y = 0;
        }


        //if(!UF_USE_HLL2)
        if(UF_SOLIDPHASE)
        {
            volx1 = UF_OUTORMV(_dem,r,c+1)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_f->data[r][c+1] + _s->data[r][c+1])/(_dx*_dx) + _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1))));
            volx2 = UF_OUTORMV(_dem,r,c-1)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_f->data[r][c-1] + _s->data[r][c-1])/(_dx*_dx) + _dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1))));
            voly1 = UF_OUTORMV(_dem,r+1,c)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_f->data[r+1][c] + _s->data[r+1][c])/(_dx*_dx) + _dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0))));
            voly2 = UF_OUTORMV(_dem,r-1,c)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_f->data[r-1][c] + _s->data[r-1][c])/(_dx*_dx) + _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0))));
        }else
        {
            volx1 = UF_OUTORMV(_dem,r,c+1)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_f->data[r][c+1])/(_dx*_dx) + _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1))));
            volx2 = UF_OUTORMV(_dem,r,c-1)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_f->data[r][c-1])/(_dx*_dx) + _dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1))));
            voly1 = UF_OUTORMV(_dem,r+1,c)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_f->data[r+1][c])/(_dx*_dx) + _dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0))));
            voly2 = UF_OUTORMV(_dem,r-1,c)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_f->data[r-1][c])/(_dx*_dx) + _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0))));
        }

        {
            volx1 = (1.0-phi_x) * volx1 + (phi_x)*_f->data[r][c];
            volx2 = (1.0-phi_x) * volx2 + (phi_x)*_f->data[r][c];
            voly1 = (1.0-phi_y) * voly1 + (phi_y)*_f->data[r][c];
            voly2 = (1.0-phi_y) * voly2 + (phi_y)*_f->data[r][c];
        }

        double rxh_lr = 0;
        double lxh_lr = 0;
        double ryh_lr = 0;
        double lyh_lr = 0;

        if(!UF_SOLIDPHASE)
        {
            rxh_lr = UF_OUTORMV(_dem,r,c+1)?0.0:std::min(rxh,std::max(0.0,(vx1 > 0? 1.0:1.0)*((_dem->Drc + rxh ) - ( _f->data[r][c+1]/(_dx*_dx) + _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1)))));
            lxh_lr = UF_OUTORMV(_dem,r,c-1)?0.0:std::min(lxh,std::max(0.0,(vx2 > 0? 1.0:1.0)*((_dem->Drc + lxh ) - (_f->data[r][c-1]/(_dx*_dx) + _dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1)))));
            ryh_lr = UF_OUTORMV(_dem,r+1,c)?0.0:std::min(ryh,std::max(0.0,(vy1 > 0? 1.0:1.0)*((_dem->Drc + ryh ) - (_f->data[r+1][c]/(_dx*_dx) + _dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0)))));
            lyh_lr = UF_OUTORMV(_dem,r-1,c)?0.0:std::min(lyh,std::max(0.0,(vy2 > 0? 1.0:1.0)*((_dem->Drc + lyh ) - (_f->data[r-1][c]/(_dx*_dx) + _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0)))));
        }else
        {
            rxh_lr = UF_OUTORMV(_dem,r,c+1)?0.0:std::min(rxh,std::max(0.0,(vx1 > 0? 1.0:1.0)*((_dem->Drc + rxh + s/(_dx*_dx)) - (_s->data[r][c+1]/(_dx*_dx) + _f->data[r][c+1]/(_dx*_dx) + _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1)))));
            lxh_lr = UF_OUTORMV(_dem,r,c-1)?0.0:std::min(lxh,std::max(0.0,(vx2 > 0? 1.0:1.0)*((_dem->Drc + lxh + s/(_dx*_dx)) - (_s->data[r][c-1]/(_dx*_dx) + _f->data[r][c-1]/(_dx*_dx) + _dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1)))));
            ryh_lr = UF_OUTORMV(_dem,r+1,c)?0.0:std::min(ryh,std::max(0.0,(vy1 > 0? 1.0:1.0)*((_dem->Drc + ryh + s/(_dx*_dx)) - (_s->data[r+1][c]/(_dx*_dx) + _f->data[r+1][c]/(_dx*_dx) + _dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0)))));
            lyh_lr = UF_OUTORMV(_dem,r-1,c)?0.0:std::min(lyh,std::max(0.0,(vy2 > 0? 1.0:1.0)*((_dem->Drc + lyh + s/(_dx*_dx)) - (_s->data[r-1][c]/(_dx*_dx) + _f->data[r-1][c]/(_dx*_dx) + _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0)))));
        }
            //new flux limiter test seems to work for both deep and shallow flow on flat and sloped surfaces. Note: Might require lax for stability
            double rxh_hr = UF_OUTORMV(_dem,r,c+1)?0.0:std::min(rxh,std::max(0.0,rxh + ((_dem->Drc) - ( _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1)))));
            double lxh_hr = UF_OUTORMV(_dem,r,c-1)?0.0:std::min(lxh,std::max(0.0,lxh + ((_dem->Drc) - (_dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1)))));
            double ryh_hr = UF_OUTORMV(_dem,r+1,c)?0.0:std::min(ryh,std::max(0.0,ryh + ((_dem->Drc) - (_dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0)))));
            double lyh_hr = UF_OUTORMV(_dem,r-1,c)?0.0:std::min(lyh,std::max(0.0,lyh + ((_dem->Drc) - ( _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0)))));

            rxh = (1.0-phi_x) * rxh_lr + (phi_x)*rxh_hr;
            lxh = (1.0-phi_x) * lxh_lr + (phi_x)*lxh_hr;
            ryh = (1.0-phi_y) * ryh_lr + (phi_y)*ryh_hr;
            lyh = (1.0-phi_y) * lyh_lr + (phi_y)*lyh_hr;


        double outlet = 0.0;
        outlet = UF_OUTORMV(_dem,r,c+1)? 1.0 : outlet;
        outlet = UF_OUTORMV(_dem,r,c-1)? 1.0 : outlet;
        outlet = UF_OUTORMV(_dem,r+1,c)? 1.0 : outlet;
        outlet = UF_OUTORMV(_dem,r-1,c)? 1.0 : outlet;

        double cq =  UF2D_COURANTSCHEMEFACTOR *UF_Courant * _f->Drc;
        double cqx1 =  UF2D_COURANTSCHEMEFACTOR *UF_Courant * (UF_OUTORMV(_dem,r,c+1)? 0.0 : volx1);
        double cqx2 =  UF2D_COURANTSCHEMEFACTOR *UF_Courant * (UF_OUTORMV(_dem,r,c-1)? 0.0 : volx2);
        double cqy1 =  UF2D_COURANTSCHEMEFACTOR *UF_Courant * (UF_OUTORMV(_dem,r+1,c)? 0.0 : voly1);
        double cqy2 =  UF2D_COURANTSCHEMEFACTOR *UF_Courant * (UF_OUTORMV(_dem,r-1,c)? 0.0 : voly2);

        double qx1 = UF_AVERAGEFACTOR * dtx1 * (vx1) * rxh *_dx;
        double qx2 = UF_AVERAGEFACTOR * dtx2 * (vx2) * lxh *_dx;
        double qy1 = UF_AVERAGEFACTOR * dty1 * (vy1) * ryh *_dx;
        double qy2 = UF_AVERAGEFACTOR * dty2 * (vy2) * lyh *_dx;

        double qx1old = qx1;
        double qx2old = qx2;
        double qy1old = qy1;
        double qy2old = qy2;

        double win = 1.0;
        double wout = 1.0;

        qx1 = ((qx1 > 0)? win : -wout) * std::min(std::fabs(qx1),(qx1 > 0)? cq : cqx1);
        qx2 = ((qx2 > 0)? wout : -win) * std::min(std::fabs(qx2),(qx2 > 0)? cqx2 : cq);
        qy1 = ((qy1 > 0)? win : -wout) * std::min(std::fabs(qy1),(qy1 > 0)? cq : cqy1);
        qy2 = ((qy2 > 0)? wout : -win) * std::min(std::fabs(qy2),(qy2 > 0)? cqy2 : cq);

        double qnextx1 = UF_OUTORMV(_dem,r,c+1)? 0.0: dtx1 * (_fu->data[r][c+1]/_dx) *_f->data[r][c+1];
        double qnextx2 = UF_OUTORMV(_dem,r,c-1)? 0.0: dtx2 * (_fu->data[r][c-1]/_dx) *_f->data[r][c-1];
        double qnexty1 = UF_OUTORMV(_dem,r+1,c)? 0.0: dty1 * (_fu->data[r+1][c]/_dx) *_f->data[r+1][c];
        double qnexty2 = UF_OUTORMV(_dem,r-1,c)? 0.0: dty2 * (_fu->data[r-1][c]/_dx) *_f->data[r-1][c];

        qx1 = qx1 + (1.0 - phi_x) *UF_MinMod(qx1old - qx1,qnextx1);
        qx2 = qx2 + (1.0 - phi_x) *UF_MinMod(qx2old - qx2,qnextx2);
        qy1 = qy1 + (1.0 - phi_y) *UF_MinMod(qy1old - qy1,qnexty1);
        qy2 = qy2 + (1.0 - phi_y) *UF_MinMod(qy2old - qy2,qnexty2);

        qx1 = UF_OUTORMV(_dem,r,c+1)? UF_BoundaryFlux2D(dtx1,_dx,_dx,_f->Drc,0,_fu->Drc,_fv->Drc,su,sv,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc,N->Drc, 0,1) : qx1;
        qx2 = UF_OUTORMV(_dem,r,c-1)? -UF_BoundaryFlux2D(dtx2,_dx,_dx,_f->Drc,0,_fu->Drc,_fv->Drc,su,sv,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc, N->Drc, 0,-1) : qx2;
        qy1 = UF_OUTORMV(_dem,r+1,c)? UF_BoundaryFlux2D(dty1,_dx,_dx,_f->Drc,0,_fu->Drc,_fv->Drc,su,sv,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc, N->Drc, 1,0) : qy1;
        qy2 = UF_OUTORMV(_dem,r-1,c)? -UF_BoundaryFlux2D(dty2,_dx,_dx,_f->Drc,0,_fu->Drc,_fv->Drc,su,sv,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc, N->Drc, -1,0) : qy2;

        qx1 = ((qx1 > 0)? win : -wout) * std::min(std::fabs(qx1),(qx1 > 0)? cq : cqx1);
        qx2 = ((qx2 > 0)? wout : -win) * std::min(std::fabs(qx2),(qx2 > 0)? cqx2 : cq);
        qy1 = ((qy1 > 0)? win : -wout) * std::min(std::fabs(qy1),(qy1 > 0)? cq : cqy1);
        qy2 = ((qy2 > 0)? wout :-win) * std::min(std::fabs(qy2),(qy2 > 0)? cqy2 : cq);

        UF2D_fqx1->Drc = qx2;
        UF2D_fqx2->Drc = qx1;
        UF2D_fqy1->Drc = qy2;
        UF2D_fqy2->Drc = qy1;


    }}}

}

void TWorld::UF2D_SolidMomentum2Source(int thread,cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv)
{

    FOR_ROW_COL_UF2DMTDER
    {
        ThreadPool->UF_t1.at(thread)->Drc = (_s->Drc)/(_dx * _dx);
    }}}

    UF2D_MUSCLE(thread,_dem,dt,ThreadPool->UF_t1.at(thread),UF_MUSCLE_TARGET_IN1);
    UF2D_MUSCLE(thread,_dem,dt,_su,UF_MUSCLE_TARGET_IN2, UF_DIRECTION_X);
    UF2D_MUSCLE(thread,_dem,dt,_sv,UF_MUSCLE_TARGET_IN3, UF_DIRECTION_Y);

    FOR_ROW_COL_UF2DMTDER
    {
        double h = (_f->Drc + _s->Drc)/(_dx * _dx);
        double hs = (_s->Drc)/(_dx * _dx);
        ThreadPool->UF_t1.at(thread)->Drc = h;

        double ifa = _ifa->Drc;

        double ifasin = sin(ifa);
        double ifatan = tan(ifa);

        double dudxx = UF2D_Derivative2(_dem,_su,r,c,UF_DIRECTION_X);
        double dudyy = UF2D_Derivative2(_dem,_sv,r,c,UF_DIRECTION_Y);

        double k_act = (1.0- ifasin)/(1.0+ifasin);
        double k_pass = (1.0+ifasin)/(1.0-ifasin);
        double k_x = (dudxx) > 0? k_act: k_pass;
        double k_y = (dudyy) > 0? k_act: k_pass;

        ThreadPool->UF_t2.at(thread)->Drc = UF_Gravity * k_x *hs * hs * UF_Aspect/2.0;
        ThreadPool->UF_t3.at(thread)->Drc = UF_Gravity * k_y *hs * hs * UF_Aspect/2.0;

        ThreadPool->UF_t4.at(thread)->Drc = _f->Drc/(_dx * _dx);
        ThreadPool->UF_t5.at(thread)->Drc = _s->Drc/(_dx * _dx);

        ThreadPool->UF_t6.at(thread)->Drc = (_s->Drc+_f->Drc)/(_dx * _dx);

        ThreadPool->UF_t8.at(thread)->Drc = _dem->Drc + (_s->Drc)/(_dx * _dx);

        UF2D_sax->Drc = 0;
        UF2D_say->Drc = 0;
    }}}

    FOR_ROW_COL_UF2DMT_DT
    {

        if(!(dt->Drc > 1e-12))
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

        double hs = (_s->Drc)/(_dx * _dx);

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
        double dc = std::max(0.0,std::min(1000.0,UF_DragCoefficient(ff,sf,gamma,_visc->Drc, _rocksize->Drc, _d->Drc)));
        double pbf = -UF_Gravity;
        double lxpbf = -UF_Gravity;
        double rxpbf = -UF_Gravity;
        double lypbf = -UF_Gravity;
        double rypbf = -UF_Gravity;

        double pbs = (1.0-gamma)*pbf;
        double ifa = _ifa->Drc < 0.01? 0.3:_ifa->Drc;

        double dalx = UF_DEMACCES(_dem,ThreadPool->UF_t6.at(thread),r,c,0,-1);
        double darx = UF_DEMACCES(_dem,ThreadPool->UF_t6.at(thread),r,c,0,1);
        double daly = UF_DEMACCES(_dem,ThreadPool->UF_t6.at(thread),r,c,-1,0);
        double dary = UF_DEMACCES(_dem,ThreadPool->UF_t6.at(thread),r,c,1,0);
        double dax = (dalx + darx ) * 0.5;
        double day = (daly + dary ) * 0.5;

        double lxslope = dalx * UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double lyslope = daly * UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rxslope = darx * UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double ryslope = dary * UF2D_Derivative(_dem,_dem,r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double dhsdx = dax * UF2D_Derivative(_dem,ThreadPool->UF_t5.at(thread),r,c,UF_DIRECTION_X);
        double dhsdy = day * UF2D_Derivative(_dem,ThreadPool->UF_t5.at(thread),r,c,UF_DIRECTION_Y);
        double ldhsdx = dalx * UF2D_Derivative(_dem,ThreadPool->UF_t5.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double ldhsdy = daly * UF2D_Derivative(_dem,ThreadPool->UF_t5.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rdhsdx = darx * UF2D_Derivative(_dem,ThreadPool->UF_t5.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rdhsdy = dary * UF2D_Derivative(_dem,ThreadPool->UF_t5.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double dhdx = dax * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_X);
        double dhdy = day * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_Y);
        double ldhdx = day * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double ldhdy = day * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rdhdx = day * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rdhdy = day * UF2D_Derivative(_dem,ThreadPool->UF_t1.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double dbdx = dax * UF2D_Derivative(_dem,ThreadPool->UF_t2.at(thread),r,c,UF_DIRECTION_X);
        double dbdy = day * UF2D_Derivative(_dem,ThreadPool->UF_t3.at(thread),r,c,UF_DIRECTION_Y);
        double ldbdx = dalx * UF2D_Derivative(_dem,ThreadPool->UF_t2.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double ldbdy = daly * UF2D_Derivative(_dem,ThreadPool->UF_t3.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rdbdx = darx * UF2D_Derivative(_dem,ThreadPool->UF_t2.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rdbdy = dary * UF2D_Derivative(_dem,ThreadPool->UF_t3.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double lddemhdx = dalx * UF2D_Derivative(_dem,ThreadPool->UF_t8.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_L);
        double lddemhdy = daly * UF2D_Derivative(_dem,ThreadPool->UF_t8.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_L);
        double rddemhdx = darx * UF2D_Derivative(_dem,ThreadPool->UF_t8.at(thread),r,c,UF_DIRECTION_X,UF_DERIVATIVE_R);
        double rddemhdy = dary * UF2D_Derivative(_dem,ThreadPool->UF_t8.at(thread),r,c,UF_DIRECTION_Y,UF_DERIVATIVE_R);

        double lsax = UF2D_MomentumBalanceSolid(true,_f->Drc,lxh*(_dx*_dx),_fu->Drc, _fv->Drc, lsu, _sv->Drc, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs,lxpbf, lxslope, UF2D_SlopeY->Drc,
                                                 ldhsdx, dhsdy, ldhdx, dhdy, ldbdx, dbdy,lddemhdx,lddemhdy);

        double rsax = UF2D_MomentumBalanceSolid(true,_f->Drc,rxh*(_dx*_dx),_fu->Drc, _fv->Drc, rsu, _sv->Drc, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs,rxpbf, rxslope, UF2D_SlopeY->Drc,
                                                 rdhsdx, dhsdy, rdhdx, dhdy, rdbdx, dbdy,rddemhdx,rddemhdy);

        double lsay = UF2D_MomentumBalanceSolid(false,_f->Drc,lyh*(_dx*_dx),_fu->Drc, _fv->Drc, _su->Drc, lsv, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs,lypbf, UF2D_SlopeX->Drc, lyslope,
                                                  dhsdx, ldhsdy, dhdx, ldhdy, dbdx, ldbdy,lddemhdx,lddemhdy);

        double rsay = UF2D_MomentumBalanceSolid(false,_f->Drc,ryh*(_dx*_dx),_fu->Drc, _fv->Drc, _su->Drc, rsv, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs,rypbf, UF2D_SlopeX->Drc, ryslope,
                                                 dhsdx, rdhsdy, dhdx, rdhdy, dbdx, rdbdy,rddemhdx,rddemhdy);

        UF2D_sax1->Drc = rsax;
        UF2D_sax2->Drc = lsax;
        UF2D_say1->Drc = rsay;
        UF2D_say2->Drc = lsay;

        UF2D_sax->Drc = 0.5*(lsax+rsax);
        UF2D_say->Drc = 0.5*(lsay+rsay);


        ////////////friction and actual accaleration

            ThreadPool->UF_t6.at(thread)->Drc = lsu;// + dt->Drc * lsax;
            ThreadPool->UF_t7.at(thread)->Drc = rsu;// + dt->Drc * rsax;
            ThreadPool->UF_t6.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* lsax,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lxh,lxslope,true);
            ThreadPool->UF_t7.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* rsax,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,rxh,rxslope,true);
            ThreadPool->UF_t6.at(thread)->Drc = (ThreadPool->UF_t6.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* lsax,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lxh,lxslope,true))/2.0;
            ThreadPool->UF_t7.at(thread)->Drc = (ThreadPool->UF_t7.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* rsax,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,rxh,rxslope,true))/2.0;

            double left_l = 0;
            double left_r = 0;


            ThreadPool->UF_t6.at(thread)->Drc = UF_SFriction(sf*0.5*UF_Gravity *std::tan(ifa),ThreadPool->UF_t6.at(thread)->Drc,dt->Drc, left_l);
            ThreadPool->UF_t7.at(thread)->Drc = UF_SFriction(sf*0.5*UF_Gravity *std::tan(ifa),ThreadPool->UF_t7.at(thread)->Drc,dt->Drc, left_r);

            if(ff > UF_VERY_SMALL)
            {
                double lvpow = pow((ThreadPool->UF_t6.at(thread)->Drc-_su->Drc)*(ThreadPool->UF_t6.at(thread)->Drc-_su->Drc)+(_fv->Drc-_sv->Drc)*(_fv->Drc-_sv->Drc),0.5*(UF_j-1.0));
                double rvpow = pow((ThreadPool->UF_t7.at(thread)->Drc-_su->Drc)*(ThreadPool->UF_t7.at(thread)->Drc-_su->Drc)+(_fv->Drc-_sv->Drc)*(_fv->Drc-_sv->Drc),0.5*(UF_j-1.0));
                double lfacu = std::max(0.0,std::fabs(dt->Drc * std::max(0.0,dc*ff-left_l)));//std::max(0.0,dc*ff - left_l)));
                double rfacu = std::max(0.0,std::fabs(dt->Drc * std::max(0.0,dc*ff - left_r)));// std::max(0.0,dc*ff - left_r)));
                double ul_balance = sf * ThreadPool->UF_t6.at(thread)->Drc + ff * _fu->Drc;
                double ur_balance = sf * ThreadPool->UF_t7.at(thread)->Drc + ff * _fu->Drc;
                ThreadPool->UF_t6.at(thread)->Drc = ul_balance + (ThreadPool->UF_t6.at(thread)->Drc - ul_balance)*std::exp(-lfacu);
                ThreadPool->UF_t7.at(thread)->Drc = ur_balance + (ThreadPool->UF_t7.at(thread)->Drc - ur_balance)*std::exp(-rfacu);
            }

            lsax = (ThreadPool->UF_t6.at(thread)->Drc - lsu)/dt->Drc;
            rsax = (ThreadPool->UF_t7.at(thread)->Drc - rsu)/dt->Drc;

            ThreadPool->UF_t6.at(thread)->Drc = lsv;// + dt->Drc * lsay;
            ThreadPool->UF_t7.at(thread)->Drc = rsv;// + dt->Drc * rsay;

            ThreadPool->UF_t6.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* lsay,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lyh,lyslope,true);
            ThreadPool->UF_t7.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* rsay,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,ryh,ryslope,true);
            ThreadPool->UF_t6.at(thread)->Drc = (ThreadPool->UF_t6.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* lsay,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lyh,lyslope,true))/2.0;
            ThreadPool->UF_t7.at(thread)->Drc = (ThreadPool->UF_t7.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* rsay,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,ryh,ryslope,true))/2.0;

            left_l = 0;
            left_r = 0;

            ThreadPool->UF_t6.at(thread)->Drc = UF_SFriction(sf*0.5*UF_Gravity *std::tan(ifa),ThreadPool->UF_t6.at(thread)->Drc,dt->Drc, left_l);
            ThreadPool->UF_t7.at(thread)->Drc = UF_SFriction(sf*0.5*UF_Gravity *std::tan(ifa),ThreadPool->UF_t7.at(thread)->Drc,dt->Drc, left_r);

            if(ff > UF_VERY_SMALL)
            {
                double lvpow = pow((ThreadPool->UF_t6.at(thread)->Drc-_sv->Drc)*(ThreadPool->UF_t6.at(thread)->Drc-_su->Drc)+(_fu->Drc-_su->Drc)*(_fv->Drc-_su->Drc),0.5*(UF_j-1.0));
                double rvpow = pow((ThreadPool->UF_t7.at(thread)->Drc-_sv->Drc)*(ThreadPool->UF_t7.at(thread)->Drc-_su->Drc)+(_fu->Drc-_su->Drc)*(_fv->Drc-_su->Drc),0.5*(UF_j-1.0));
                double lfacu = std::max(0.0,std::fabs(dt->Drc  * std::max(0.0,dc*ff-left_l)));// std::max(0.0,dc*ff - left_l)));
                double rfacu = std::max(0.0,std::fabs(dt->Drc  * std::max(0.0,dc*ff-left_r)));// std::max(0.0,dc * ff - left_r)));
                double ul_balance = sf * ThreadPool->UF_t6.at(thread)->Drc + ff * _fv->Drc;
                double ur_balance = sf * ThreadPool->UF_t7.at(thread)->Drc + ff * _fv->Drc;
                ThreadPool->UF_t6.at(thread)->Drc = ul_balance + (ThreadPool->UF_t6.at(thread)->Drc - ul_balance)*std::exp(-lfacu);
                ThreadPool->UF_t7.at(thread)->Drc = ur_balance + (ThreadPool->UF_t7.at(thread)->Drc - ur_balance)*std::exp(-rfacu);
            }

            lsay = (ThreadPool->UF_t6.at(thread)->Drc - lsv)/dt->Drc;
            rsay = (ThreadPool->UF_t7.at(thread)->Drc - rsv)/dt->Drc;

            /*rsax = -lddemhdx * 9.81;
            lsax = -rddemhdx * 9.81;
            rsay = -lddemhdy * 9.81;
            lsay = -rddemhdy * 9.81;*/

            UF2D_sax1->Drc = rsax;
            UF2D_sax2->Drc = lsax;
            UF2D_say1->Drc = rsay;
            UF2D_say2->Drc = lsay;

        ////////////average accalerations for cell centers
        UF2D_sax->Drc = 0.5*(lsax+rsax);
        UF2D_say->Drc = 0.5*(lsay+rsay);


        ////////////fluxes
        double hdem = ((_s->Drc + _f->Drc)/(_dx*_dx) + _dem->Drc);
        double volx1;
        double volx2;
        double voly1;
        double voly2;

        double dtx1 = (UF_OUTORMV(_dem,r,c+1)? dt->Drc : std::max(dt->Drc ,dt->data[r][c+1]));
        double dtx2 = (UF_OUTORMV(_dem,r,c-1)? dt->Drc : std::max(dt->Drc ,dt->data[r][c-1]));
        double dty1 = (UF_OUTORMV(_dem,r+1,c)? dt->Drc : std::max(dt->Drc ,dt->data[r+1][c]));
        double dty2 = (UF_OUTORMV(_dem,r-1,c)? dt->Drc : std::max(dt->Drc ,dt->data[r-1][c]));

        double vx1 = rsu + dtx1 * rsax;
        double vx2 = lsu + dtx2 * lsax;
        double vy1 = rsv + dty1 * rsay;
        double vy2 = lsv + dty2 * lsay;

        double phi_rqx = (rsu + dt->Drc *  rsax) *(rxh);
        double phi_qx = (_su->Drc + dt->Drc * UF2D_sax->Drc) *(_s->Drc)/(_dx*_dx);
        double phi_lqx = (lsu + dt->Drc *  lsax) *(lxh);

        double phi_rqy = (rsv + dt->Drc *  rsay) *(ryh);
        double phi_qy = (_sv->Drc + dt->Drc * UF2D_say->Drc) *(_s->Drc)/(_dx*_dx);
        double phi_lqy = (lsv + dt->Drc *  lsay) *(lyh);

        double phi_rx = std::fabs(phi_qx - phi_lqx) > 0? (phi_rqx - phi_qx)/(phi_qx - phi_lqx) : 1.0;
        double phi_ry = std::fabs(phi_qy - phi_lqy) > 0? (phi_rqy - phi_qy)/(phi_qy - phi_lqy) : 1.0;

        double phi_rx_i = std::fabs(phi_rqx - phi_qx) > 0? (phi_qx - phi_lqx)/(phi_rqx - phi_qx) : 1.0;
        double phi_ry_i = std::fabs(phi_rqy - phi_qy) > 0? (phi_qy - phi_lqy)/(phi_rqy - phi_qy) : 1.0;

        double phi_x = std::max(0.0,std::min(1.0, std::min(phi_rx,phi_rx)));
        double phi_y = std::max(0.0,std::min(1.0, std::min(phi_ry,phi_ry)));

        //double phi_x = std::min(1.0,std::max(0.0,std::max(std::min(1.0, 2.0*std::max(phi_rx,phi_rx_i)),std::min(2.0,std::max(phi_rx,phi_rx_i)) )));
        //double phi_y = std::min(1.0,std::max(0.0,std::max(std::min(1.0, 2.0*std::max(phi_ry,phi_ry_i)),std::min(2.0,std::max(phi_ry,phi_ry_i)) )));

        if(!UF_USE_HLL2)
        {
            phi_x = 0;
            phi_y = 0;
        }


        //if(!UF_USE_HLL2)
        {
            volx1 = UF_OUTORMV(_dem,r,c+1)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_s->data[r][c+1])/(_dx*_dx) + _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1))));
            volx2 = UF_OUTORMV(_dem,r,c-1)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_s->data[r][c-1])/(_dx*_dx) + _dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1))));
            voly1 = UF_OUTORMV(_dem,r+1,c)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_s->data[r+1][c])/(_dx*_dx) + _dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0))));
            voly2 = UF_OUTORMV(_dem,r-1,c)?0.0:std::max(0.0,(_dx*_dx) * (hdem - ((_s->data[r-1][c])/(_dx*_dx) + _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0))));

        }
        {
            volx1 = (1.0-phi_x) * volx1 + (phi_x)*_s->data[r][c];
            volx2 = (1.0-phi_x) * volx2 + (phi_x)*_s->data[r][c];
            voly1 = (1.0-phi_y) * voly1 + (phi_y)*_s->data[r][c];
            voly2 = (1.0-phi_y) * voly2 + (phi_y)*_s->data[r][c];

        }



            double rxh_lr = UF_OUTORMV(_dem,r,c+1)?0.0:std::min(rxh,std::max(0.0,(vx1 > 0? 1.0:1.0)*((_dem->Drc + rxh + _s->Drc/(_dx*_dx)) - (_s->data[r][c+1]/(_dx*_dx) + _s->data[r][c+1]/(_dx*_dx) + _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1)))));
            double lxh_lr = UF_OUTORMV(_dem,r,c-1)?0.0:std::min(lxh,std::max(0.0,(vx2 > 0? 1.0:1.0)*((_dem->Drc + lxh + _s->Drc/(_dx*_dx)) - (_s->data[r][c-1]/(_dx*_dx) + _s->data[r][c-1]/(_dx*_dx) + _dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1)))));
            double ryh_lr = UF_OUTORMV(_dem,r+1,c)?0.0:std::min(ryh,std::max(0.0,(vy1 > 0? 1.0:1.0)*((_dem->Drc + ryh + _s->Drc/(_dx*_dx)) - (_s->data[r+1][c]/(_dx*_dx) + _s->data[r+1][c]/(_dx*_dx) + _dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0)))));
            double lyh_lr = UF_OUTORMV(_dem,r-1,c)?0.0:std::min(lyh,std::max(0.0,(vy2 > 0? 1.0:1.0)*((_dem->Drc + lyh + _s->Drc/(_dx*_dx)) - (_s->data[r-1][c]/(_dx*_dx) + _s->data[r-1][c]/(_dx*_dx) + _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0)))));


            //new slux limiter test seems to work sor both deep and shallow slow on slat and sloped sursaces. Note: Might require lax sor stability
            double rxh_hr = UF_OUTORMV(_dem,r,c+1)?0.0:std::min(rxh,std::max(0.0,rxh + ((_dem->Drc) - ( _dem->data[r][c+1] +GetFlowBarrierHeight(r,c,0,1)))));
            double lxh_hr = UF_OUTORMV(_dem,r,c-1)?0.0:std::min(lxh,std::max(0.0,lxh + ((_dem->Drc) - (_dem->data[r][c-1] +GetFlowBarrierHeight(r,c,0,-1)))));
            double ryh_hr = UF_OUTORMV(_dem,r+1,c)?0.0:std::min(ryh,std::max(0.0,ryh + ((_dem->Drc) - (_dem->data[r+1][c] +GetFlowBarrierHeight(r,c,1,0)))));
            double lyh_hr = UF_OUTORMV(_dem,r-1,c)?0.0:std::min(lyh,std::max(0.0,lyh + ((_dem->Drc) - ( _dem->data[r-1][c] +GetFlowBarrierHeight(r,c,-1,0)))));

            rxh = (1.0-phi_x) * rxh_lr + (phi_x)*rxh_hr;
            lxh = (1.0-phi_x) * lxh_lr + (phi_x)*lxh_hr;
            ryh = (1.0-phi_y) * ryh_lr + (phi_y)*ryh_hr;
            lyh = (1.0-phi_y) * lyh_lr + (phi_y)*lyh_hr;


        double outlet = 0.0;
        outlet = UF_OUTORMV(_dem,r,c+1)? 1.0 : outlet;
        outlet = UF_OUTORMV(_dem,r,c-1)? 1.0 : outlet;
        outlet = UF_OUTORMV(_dem,r+1,c)? 1.0 : outlet;
        outlet = UF_OUTORMV(_dem,r-1,c)? 1.0 : outlet;

        double cq =  UF2D_COURANTSCHEMEFACTOR *UF_Courant * _s->Drc;
        double cqx1 =  UF2D_COURANTSCHEMEFACTOR *UF_Courant * (UF_OUTORMV(_dem,r,c+1)? 0.0 : volx1);
        double cqx2 =  UF2D_COURANTSCHEMEFACTOR *UF_Courant * (UF_OUTORMV(_dem,r,c-1)? 0.0 : volx2);
        double cqy1 =  UF2D_COURANTSCHEMEFACTOR *UF_Courant * (UF_OUTORMV(_dem,r+1,c)? 0.0 : voly1);
        double cqy2 =  UF2D_COURANTSCHEMEFACTOR *UF_Courant * (UF_OUTORMV(_dem,r-1,c)? 0.0 : voly2);

        double qx1 = UF_AVERAGEFACTOR * dtx1 * (vx1) * rxh *_dx;
        double qx2 = UF_AVERAGEFACTOR * dtx2 * (vx2) * lxh *_dx;
        double qy1 = UF_AVERAGEFACTOR * dty1 * (vy1) * ryh *_dx;
        double qy2 = UF_AVERAGEFACTOR * dty2 * (vy2) * lyh *_dx;

        double qx1old = qx1;
        double qx2old = qx2;
        double qy1old = qy1;
        double qy2old = qy2;

        double win = 1.0;
        double wout = 1.0;

        qx1 = ((qx1 > 0)? win : -wout) * std::min(std::fabs(qx1),(qx1 > 0)? cq : cqx1);
        qx2 = ((qx2 > 0)? wout : -win) * std::min(std::fabs(qx2),(qx2 > 0)? cqx2 : cq);
        qy1 = ((qy1 > 0)? win : -wout) * std::min(std::fabs(qy1),(qy1 > 0)? cq : cqy1);
        qy2 = ((qy2 > 0)? wout : -win) * std::min(std::fabs(qy2),(qy2 > 0)? cqy2 : cq);

        double qnextx1 = UF_OUTORMV(_dem,r,c+1)? 0.0: dtx1 * (_su->data[r][c+1]/_dx) *_s->data[r][c+1];
        double qnextx2 = UF_OUTORMV(_dem,r,c-1)? 0.0: dtx2 * (_su->data[r][c-1]/_dx) *_s->data[r][c-1];
        double qnexty1 = UF_OUTORMV(_dem,r+1,c)? 0.0: dty1 * (_su->data[r+1][c]/_dx) *_s->data[r+1][c];
        double qnexty2 = UF_OUTORMV(_dem,r-1,c)? 0.0: dty2 * (_su->data[r-1][c]/_dx) *_s->data[r-1][c];

        qx1 = qx1 + (1.0 - phi_x) *UF_MinMod(qx1old - qx1,qnextx1);
        qx2 = qx2 + (1.0 - phi_x) *UF_MinMod(qx2old - qx2,qnextx2);
        qy1 = qy1 + (1.0 - phi_y) *UF_MinMod(qy1old - qy1,qnexty1);
        qy2 = qy2 + (1.0 - phi_y) *UF_MinMod(qy2old - qy2,qnexty2);

        qx1 = UF_OUTORMV(_dem,r,c+1)? UF_BoundaryFlux2D(dtx1,_dx,_dx,_s->Drc,0,_su->Drc,_sv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc,N->Drc, 0,1) : qx1;
        qx2 = UF_OUTORMV(_dem,r,c-1)? -UF_BoundaryFlux2D(dtx2,_dx,_dx,_s->Drc,0,_su->Drc,_sv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc, N->Drc, 0,-1) : qx2;
        qy1 = UF_OUTORMV(_dem,r+1,c)? UF_BoundaryFlux2D(dty1,_dx,_dx,_s->Drc,0,_su->Drc,_sv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc, N->Drc, 1,0) : qy1;
        qy2 = UF_OUTORMV(_dem,r-1,c)? -UF_BoundaryFlux2D(dty2,_dx,_dx,_s->Drc,0,_su->Drc,_sv->Drc,_su->Drc,_sv->Drc,UF2D_SlopeX->Drc,UF2D_SlopeY->Drc, N->Drc, -1,0) : qy2;

        qx1 = ((qx1 > 0)? win : -wout) * std::min(std::fabs(qx1),(qx1 > 0)? cq : cqx1);
        qx2 = ((qx2 > 0)? wout : -win) * std::min(std::fabs(qx2),(qx2 > 0)? cqx2 : cq);
        qy1 = ((qy1 > 0)? win : -wout) * std::min(std::fabs(qy1),(qy1 > 0)? cq : cqy1);
        qy2 = ((qy2 > 0)? wout :-win) * std::min(std::fabs(qy2),(qy2 > 0)? cqy2 : cq);


        UF2D_sqx1->Drc = qx2;
        UF2D_sqx2->Drc = qx1;
        UF2D_sqy1->Drc = qy2;
        UF2D_sqy2->Drc = qy1;

        UF2D_sqx->Drc = (std::fabs(qx1) + std::fabs(qx2))/(_dx*(hs));
        UF2D_sqy->Drc = (std::fabs(qy1) + std::fabs(qy2))/(_dx*(hs));


    }}}

}


void TWorld::UF1D_FluidMomentum2Source(int thread,cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu)
{

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    FOR_ROW_COL_UF1DMTDER
    {
        double s = 0;
        if(UF_SOLIDPHASE)
        {
            s = _s->Drc;
        }

        double h = (_f->Drc + s)/(_dx * _lddw->Drc);
        ThreadPool->UF_t1.at(thread)->Drc = h*h*(UF_Gravity*h)/2.0;
        if(_f->Drc + s > UF_VERY_SMALL)
        {
            ThreadPool->UF_t2.at(thread)->Drc = _f->Drc/(_f->Drc + s);
            ThreadPool->UF_t3.at(thread)->Drc = s/(_f->Drc + s);
        }else
        {
            ThreadPool->UF_t2.at(thread)->Drc = 0;
            ThreadPool->UF_t3.at(thread)->Drc = 0;
        }

        ThreadPool->UF_t4.at(thread)->Drc = _f->Drc/(_dx * _lddw->Drc);
        ThreadPool->UF_t5.at(thread)->Drc = s/(_dx * _lddw->Drc);
        UF1D_fq1->Drc = 0;
        UF1D_fq2->Drc = 0;
        UF1D_fa->Drc = 0;
        UF1D_fa1->Drc = 0;
        UF1D_fa2->Drc = 0;
    }}}
    UF1D_MUSCLE(thread,_ldd,_lddw,dt,ThreadPool->UF_t4.at(thread),UF_MUSCLE_TARGET_IN1);
    UF1D_MUSCLE(thread,_ldd,_lddw,dt,_fu,UF_MUSCLE_TARGET_IN2);

    FOR_ROW_COL_UF1DMT_DT
    {
        if(!(_f->Drc > UF_VERY_SMALL))
        {
            UF1D_fa->Drc = 0;
            continue;
        }

        double s = 0;
        double su = 0;
        if(UF_SOLIDPHASE)
        {
            s = _s->Drc;
            su = _su->Drc;
        }

        double h = (_f->Drc + s)/(_dx *_lddw->Drc);
        double rhf = std::max(0.0,(UF1D_MUSCLE_1_x1->Drc));
        double lhf = std::max(0.0,(UF1D_MUSCLE_1_x2->Drc));

        double rfu = (UF1D_MUSCLE_2_x1->Drc);
        double lfu = (UF1D_MUSCLE_2_x2->Drc);
        double ff = ThreadPool->UF_t2.at(thread)->Drc;
        double sf = ThreadPool->UF_t3.at(thread)->Drc;
        _visc->Drc = 0;
        if(UF_SOLIDPHASE)
        {
            _visc->Drc = UF_DynamicViscosity(sf + ((SwitchErosion && UF_SUSPENDEDVISCOSITY)? (((UF1D_blm->Drc + UF1D_ssm->Drc)/2000.0 + ff) > UF_VERY_SMALL ?((UF1D_blm->Drc + UF1D_ssm->Drc)/2000.0)/((UF1D_blm->Drc + UF1D_ssm->Drc)/2000.0 + ff) :0.0):0.0));
        }
        if(UF_SOLIDPHASE)
        {
            //_d->Drc = 2000.0;
        }
        double Nr = 1.0;
        if(UF_SOLIDPHASE)
        {
            Nr = UF_Reynolds(_d->Drc,_visc->Drc,ff,sf, _rocksize->Drc);
            UF1D_Nr->Drc = Nr;
        }

        double Nra = 15000.0;
        double gamma = 1.0;
        if(UF_SOLIDPHASE)
        {
            gamma = (!(_d->Drc > UF_VERY_SMALL))? 0.5 : 1000.0/_d->Drc;
        }
        double dc = 0.0;
        if(UF_SOLIDPHASE)
        {
            dc =UF_DragCoefficient(ff, sf, gamma,_visc->Drc, _rocksize->Drc, _d->Drc);
        }

        double pbf = UF_Gravity*h;
        double lpbf = UF_Gravity*lhf;
        double rpbf = UF_Gravity*rhf;

        double dhfdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t4.at(thread),r,c);
        double ldhfdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t4.at(thread),r,c,false,UF_DERIVATIVE_L);
        double rdhfdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t4.at(thread),r,c,false,UF_DERIVATIVE_R);

        double dh2pbdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t1.at(thread),r,c);
        double ldh2pbdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t1.at(thread),r,c,false,UF_DERIVATIVE_L);
        double rdh2pbdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t1.at(thread),r,c,false,UF_DERIVATIVE_R);

        double dsfdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t3.at(thread),r,c);
        double ddsfdxx = UF1D_Derivative2(_ldd,_lddw,ThreadPool->UF_t3.at(thread),r,c);
        double dfudx = UF1D_Derivative(_ldd,_lddw,_fu,r,c);
        double ddfudxx = UF1D_Derivative2(_ldd,_lddw,_fu,r,c);
        double dsudx = 0;
        if(UF_SOLIDPHASE)
        {
            dsudx = UF1D_Derivative(_ldd,_lddw,_su,r,c);
        }

        double ifa = 0;
        if(UF_SOLIDPHASE)
        {
             ifa = _ifa->Drc;
        }
        double lfa = UF1D_MomentumBalanceFluid( lhf*(_dx*_lddw->Drc),s,lfu, su, ff, sf, Nr, Nra,ifa, gamma, _visc->Drc, lpbf, UF1D_Slope->Drc,
                                                 ldhfdx, ldh2pbdx, dsfdx, ddsfdxx, dfudx, ddfudxx, dsudx);
        double rfa = UF1D_MomentumBalanceFluid(  rhf*(_dx*_lddw->Drc),s,rfu,su, ff, sf, Nr, Nra, ifa, gamma, _visc->Drc, rpbf, UF1D_Slope->Drc,
                                                 rdhfdx, rdh2pbdx, dsfdx, ddsfdxx, dfudx, ddfudxx, dsudx);


        UF1D_fa2->Drc = lfa;
        UF1D_fa1->Drc = rfa;

        ThreadPool->UF_t6.at(thread)->Drc = lfu;
        ThreadPool->UF_t7.at(thread)->Drc = rfu;

        ThreadPool->UF_t6.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* lfa,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lhf,UF1D_Slope->Drc, false,true,_lddw->Drc);
        ThreadPool->UF_t7.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* rfa,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,rhf,UF1D_Slope->Drc, false,true,_lddw->Drc);
        ThreadPool->UF_t6.at(thread)->Drc = (ThreadPool->UF_t6.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* lfa,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lhf,UF1D_Slope->Drc, false,true,_lddw->Drc))/2.0;
        ThreadPool->UF_t7.at(thread)->Drc = (ThreadPool->UF_t7.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* rfa,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,rhf,UF1D_Slope->Drc, false,true,_lddw->Drc))/2.0;

        if(sf > UF_VERY_SMALL && dt->Drc > 0)
        {
            double lvpow = pow((ThreadPool->UF_t6.at(thread)->Drc-su)*(ThreadPool->UF_t6.at(thread)->Drc-su),0.5*(UF_j-1.0));
            double rvpow = pow((ThreadPool->UF_t7.at(thread)->Drc-su)*(ThreadPool->UF_t7.at(thread)->Drc-su),0.5*(UF_j-1.0));
            double lfacu = std::max(0.0,std::fabs(dt->Drc * dc * sf*lvpow));
            double rfacu = std::max(0.0,std::fabs(dt->Drc * dc * sf*rvpow));
            double ul_balance = ff * ThreadPool->UF_t6.at(thread)->Drc + sf * su;
            double ur_balance = ff * ThreadPool->UF_t7.at(thread)->Drc + sf * su;
            ThreadPool->UF_t6.at(thread)->Drc = ul_balance + (ThreadPool->UF_t6.at(thread)->Drc - ul_balance)*std::exp(-lfacu);
            ThreadPool->UF_t7.at(thread)->Drc = ur_balance + (ThreadPool->UF_t7.at(thread)->Drc - ur_balance)*std::exp(-rfacu);
        }

        lfa = (ThreadPool->UF_t6.at(thread)->Drc - lfu)/dt->Drc;
        rfa = (ThreadPool->UF_t7.at(thread)->Drc - rfu)/dt->Drc;

        UF1D_fa->Drc = (lfa + rfa) / 2.0;

        double cqself = UF1D_COURANTSCHEMEFACTOR * UF_Courant * _f->Drc;
        double hself = _f->Drc/(_dx*_lddw->Drc);


        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!(lddself == 5))
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){

                ////////////fluxes
                double volx = std::fabs(UF_LDDOUT(_ldd,r,c,true)?0.0:std::min(0.0,(_dx*_lddw->Drc) * (-(hself - ((_f->data[r2][c2])/(_dx*_lddw->data[r2][c2]) + (UF1D_Slope->Drc + UF1D_Slope->data[r2][c2])*0.5 * _dx)))));
                double cq = UF1D_COURANTSCHEMEFACTOR * volx;
                double volxr = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * (hself - (_f->data[r2][c2]/(_dx*_lddw->data[r2][c2]) + (UF1D_Slope->Drc + UF1D_Slope->data[r2][c2])*0.5  * _dx)));

                double dtx1 = dt->Drc;
                double vxr = (rfu + dtx1 * rfa);

                double cqt = UF1D_COURANTSCHEMEFACTOR * (volxr);

                double q1 = UF_AVERAGEFACTOR * dtx1 * (vxr) *hself*_lddw->Drc;

                if(SwitchChannelMaxCS)
                {
                    double h_max = UF1D_ChannelMaxCS->Drc;
                    if(h_max > 0.0)
                    {
                         q1 = UF_AVERAGEFACTOR * dtx1 * (vxr) *std::min(ff*h_max,hself)*_lddw->Drc;
                    }
                }

                double qold = q1;
                q1 = ((q1 > 0)? 1.0 : -1.0) * std::min(std::fabs(q1),(q1 > 0)? std::fabs(cq) : cqt);

                if(SwitchChannelConnection)
                {
                    double connection = UF1D_ChannelConnected->Drc;
                    if(connection == 0.0)
                    {
                        double volx_r = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * hself);
                        double volxr_r = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->data[r2][c2]) * (_f->data[r2][c2]/(_dx*_lddw->data[r2][c2])));
                        double volx_max = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * _lddh->Drc);
                        double volxr_max = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->data[r2][c2]) *  _lddh->data[r2][c2]);
                        q1 = ((q1 > 0)? 1.0 : -1.0) * std::min(std::fabs(q1),(q1 > 0)? ff*std::max(0.0,volxr_max-volxr_r) : ff*std::max(0.0,volx_max-volx_r));
                    }
                }



                double volxr_max = 0;

                double qleft = qold - q1;
                double qnext = dtx1 * (_fu->data[r2][c2]/_dx) *_f->data[r2][c2];

                UF1D_fq1->Drc = q1;// + UF_MinMod(qleft,qnext);

                if(_ldd->data[r2][c2] == 5)
                {

                    double dtx1 = dt->Drc;
                    double q1b = UF_BoundaryFlux1D(dtx1,_lddw->Drc,_f->Drc/(_dx*_lddw->Drc),0,rfu,su,UF1D_Slope->Drc,ChannelN->Drc, true);
                    q1b = ((q1b > 0)? 1.0 : 1.0) * std::min(std::fabs(q1b),(q1b > 0)? cqself : cqself);
                    UF1D_fq1->Drc =  (UF1D_fq1->Drc + q1b)/2.0;
                }

            }else// if(rfu > 0)
            {
                double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
                double q1 = UF_BoundaryFlux1D(dtx1,_lddw->Drc,_f->Drc/(_dx*_lddw->Drc),0,rfu,su,UF1D_Slope->Drc,ChannelN->Drc, true);
                q1 = ((q1 > 0)? 1.0 : 1.0) * std::min(std::fabs(q1),(q1 > 0)? cqself : cqself);

                UF1D_fq1->Drc = 0;//q1;
            }
        }else// if (rfu > 0)
        {
            double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
            double q1 = UF_BoundaryFlux1D(dtx1,_lddw->Drc,_f->Drc,0,rfu,su,UF1D_Slope->Drc,ChannelN->Drc, true);

            q1 = ((q1 > 0)? 1.0 : 1.0) * std::min(std::fabs(q1),(q1 > 0)? cqself : cqself);

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
                        double volx = std::fabs(UF_LDDOUT(_ldd,r,c,false)?0.0:std::min(0.0,(_dx*_lddw->Drc) * (-(hself - (_f->data[r2][c2]/(_dx*_lddw->data[r2][c2]) - (UF1D_Slope->Drc + UF1D_Slope->data[r2][c2])*0.5  * _dx)))));
                        double cq = UF1D_COURANTSCHEMEFACTOR *UF_Courant *  volx;

                        double volxl = UF_LDDOUT(_ldd,r,c,false)?0.0:std::max(0.0,(_dx*_lddw->Drc) * (hself - (_f->data[r2][c2]/(_dx*_lddw->data[r2][c2]) - (UF1D_Slope->Drc + UF1D_Slope->data[r2][c2])*0.5  * _dx)));
                        double dtx1 = dt->Drc;//(UF1D_MUSCLE_OUT_x1->Drc > 0)? dt->Drc : (UF_OUTORMV(_dem,r,c+1)? dt->Drc : (UF_NOTIME(_dem,dt,r,c+1)? 0.0: dt->data[r][c+1]));
                        double vxl = (lfu + dtx1 * lfa);
                        //lhf = std::min(lhf,std::max(0.0,(vxl > 0? 1.0:1.0)*((lhf) - (_f->data[r2][c2]/(_dx*_dx) - UF1D_Slope->Drc * _dx))));

                        double cqt = UF1D_COURANTSCHEMEFACTOR *(volxl);
                        double q2 = UF_AVERAGEFACTOR * dtx1 * (vxl/_dx) * _f->Drc* (_lddw->data[r2][c2]/totalwidth);
                        double qold = q2;

                        if(SwitchChannelMaxCS)
                        {
                            double h_max = UF1D_ChannelMaxCS->Drc;
                            if(h_max > 0.0)
                            {
                                 q2 = UF_AVERAGEFACTOR * dtx1 * (vxl) *std::min(ff * h_max,hself)*_lddw->Drc;
                            }
                        }

                        q2 = ((q2 > 0)? 1.0 : -1.0) *std::min(std::fabs(q2),(q2 > 0)? cqt : std::fabs(cq));


                        if(SwitchChannelConnection)
                        {
                            double connection = UF1D_ChannelConnected->Drc;
                            if(connection == 0.0)
                            {
                                double volx_r = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * hself);
                                double volxr_r = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->data[r2][c2]) * (_f->data[r2][c2]/(_dx*_lddw->data[r2][c2])));
                                double volx_max = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * _lddh->Drc);
                                double volxr_max = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->data[r2][c2]) *  _lddh->data[r2][c2]);
                                q2 = ((q2 > 0)? 1.0 : -1.0) * std::min(std::fabs(q2),(q2 > 0)? ff*std::max(0.0,volxr_max-volxr_r) : ff*std::max(0.0,volx_max-volx_r));
                            }
                        }

                        double qleft = qold - q2;
                        double qnext = dtx1 * (_fu->data[r2][c2]/_dx) *_f->data[r2][c2];


                        UF1D_fq2->Drc += q2 + UF_MinMod(qleft,qnext);
                    }
                }
            }
        }else// if(lfu < 0)
        {
            double dtx2 = dt->Drc;
            double q2 = UF_BoundaryFlux1D(dtx2,_lddw->Drc,_f->Drc,0,lfu,su,UF1D_Slope->Drc,ChannelN->Drc, true);
            q2 = ((q2 < 0)? 1.0 : 1.0) * std::min(std::fabs(q2),(q2 > 0)? cqself : cqself);

            UF1D_fq2->Drc = 0;//q2;
        }
    }}}




}

void TWorld::UF1D_SolidMomentum2Source(int thread,cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_su)
{
    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    FOR_ROW_COL_UF1DMTDER
    {
        out_su->Drc = _su->Drc;
        double h = (_f->Drc + _s->Drc)/(_dx * _lddw->Drc);
        double hs = (_s->Drc)/(_dx * _lddw->Drc);
        ThreadPool->UF_t1.at(thread)->Drc = h;

        double ifa = _ifa->Drc;

        double ifasin = sin(ifa);
        double ifatan = tan(ifa);

        double dudxx = UF1D_Derivative2(_ldd,_lddw,_su,r,c);

        double k_act = (1- ifasin)/(1+ifasin);
        double k_pass = (1+ifasin)/(1-ifasin);
        double k_x = (dudxx) > 0? k_act: k_pass;

        ThreadPool->UF_t2.at(thread)->Drc = UF_Gravity * k_x *hs * hs * UF_Aspect/2.0;

        ThreadPool->UF_t4.at(thread)->Drc = _f->Drc/(_dx * _lddw->Drc);
        ThreadPool->UF_t5.at(thread)->Drc = _s->Drc/(_dx * _lddw->Drc);
        UF1D_sq1->Drc = 0;
        UF1D_sq2->Drc = 0;
        UF1D_sa->Drc = 0;
        UF1D_sa1->Drc = 0;
        UF1D_sa2->Drc = 0;
    }}}

    UF1D_MUSCLE(thread,_ldd,_lddw,dt,ThreadPool->UF_t5.at(thread),UF_MUSCLE_TARGET_IN1);
    UF1D_MUSCLE(thread,_ldd,_lddw,dt,_su,UF_MUSCLE_TARGET_IN2);

    FOR_ROW_COL_UF1DMT_DT
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

        double dhsdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t5.at(thread),r,c);
        double ldhsdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t5.at(thread),r,c, false,UF_DERIVATIVE_L);
        double rdhsdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t5.at(thread),r,c, false, UF_DERIVATIVE_R);

        double rsu = (UF1D_MUSCLE_2_x1->Drc);
        double lsu = (UF1D_MUSCLE_2_x2->Drc);
        double rhs = std::max(0.0,(UF1D_MUSCLE_1_x1->Drc));
        double lhs = std::max(0.0,(UF1D_MUSCLE_1_x2->Drc));

        double dhdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t1.at(thread),r,c);
        double dbdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t2.at(thread),r,c);
        double ldhdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t1.at(thread),r,c, false,UF_DERIVATIVE_L);
        double ldbdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t2.at(thread),r,c, false,UF_DERIVATIVE_L);
        double rdhdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t1.at(thread),r,c, false,UF_DERIVATIVE_R);
        double rdbdx = UF1D_Derivative(_ldd,_lddw,ThreadPool->UF_t2.at(thread),r,c, false,UF_DERIVATIVE_R);

        double lsa = UF1D_MomentumBalanceSolid(_f->Drc,lhs*(_dx*_lddw->Drc),_fu->Drc, lsu, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs, pbf, UF1D_Slope->Drc,
                                                 ldhsdx, ldhdx, ldbdx);
        double rsa = UF1D_MomentumBalanceSolid(_f->Drc,rhs*(_dx*_lddw->Drc),_fu->Drc, rsu, ff, sf, 0, 0, ifa, gamma, _visc->Drc, pbs, pbf, UF1D_Slope->Drc,
                                                 rdhsdx, rdhdx, rdbdx);

        UF1D_sa2->Drc = lsa;
        UF1D_sa1->Drc = rsa;

        ThreadPool->UF_t6.at(thread)->Drc = lsu;
        ThreadPool->UF_t7.at(thread)->Drc = rsu;
        ThreadPool->UF_t6.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* lsa,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lhs,UF1D_Slope->Drc, true,true,_lddw->Drc);
        ThreadPool->UF_t7.at(thread)->Drc = UF_Friction(dt->Drc*UF_TIMERATIO* rsa,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,rhs,UF1D_Slope->Drc, true,true,_lddw->Drc);
        ThreadPool->UF_t6.at(thread)->Drc = (ThreadPool->UF_t6.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* lsa,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t6.at(thread)->Drc,0,N->Drc,lhs,UF1D_Slope->Drc, true,true,_lddw->Drc))/2.0;
        ThreadPool->UF_t7.at(thread)->Drc = (ThreadPool->UF_t7.at(thread)->Drc + UF_Friction(dt->Drc*UF_TIMERATIO* rsa,dt->Drc*UF_TIMERATIO,ThreadPool->UF_t7.at(thread)->Drc,0,N->Drc,rhs,UF1D_Slope->Drc, true,true,_lddw->Drc))/2.0;

        if(sf > UF_VERY_SMALL && dt->Drc > 0)
        {
            double lvpow = pow((ThreadPool->UF_t6.at(thread)->Drc-_su->Drc)*(ThreadPool->UF_t6.at(thread)->Drc-_su->Drc),0.5*(UF_j-1.0));
            double rvpow = pow((ThreadPool->UF_t7.at(thread)->Drc-_su->Drc)*(ThreadPool->UF_t7.at(thread)->Drc-_su->Drc),0.5*(UF_j-1.0));
            double lfacu = std::max(0.0,std::fabs(dt->Drc * dc * ff*lvpow));
            double rfacu = std::max(0.0,std::fabs(dt->Drc * dc * ff*rvpow));
            double ul_balance = sf * ThreadPool->UF_t6.at(thread)->Drc + ff * _su->Drc;
            double ur_balance = sf * ThreadPool->UF_t7.at(thread)->Drc + ff * _su->Drc;
            ThreadPool->UF_t6.at(thread)->Drc = ul_balance + (ThreadPool->UF_t6.at(thread)->Drc - ul_balance)*std::exp(-lfacu);
            ThreadPool->UF_t7.at(thread)->Drc = ur_balance + (ThreadPool->UF_t7.at(thread)->Drc - ur_balance)*std::exp(-rfacu);
        }

        lsa = (ThreadPool->UF_t6.at(thread)->Drc - lsu)/dt->Drc;
        rsa = (ThreadPool->UF_t7.at(thread)->Drc - rsu)/dt->Drc;

        UF1D_sa->Drc = (lsa + rsa) / 2.0;

        double cqself = UF1D_COURANTSCHEMEFACTOR * UF_Courant * _s->Drc;

        double hself = _s->Drc/(_dx * _lddw->Drc);
        //front cell
        int lddself = (int) _ldd->data[r][c];
        if(!(lddself == 5))
        {
            int r2 = r+dy[lddself];
            int c2 = c+dx[lddself];

            if(!UF_OUTORMV(_ldd,r2,c2)){

                ////////////fluxes
                double volx = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * (-(hself - (_s->data[r2][c2]/(_dx*_lddw->data[r2][c2]) + UF1D_Slope->Drc * _dx))));
                double cq = UF1D_COURANTSCHEMEFACTOR * UF_Courant * volx;
                double volxr = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * (hself - (_s->data[r2][c2]/(_dx*_lddw->data[r2][c2]) - UF1D_Slope->Drc * _dx)));
                double dtx1 = dt->Drc;
                double vxr = rsu + dtx1 * rsa;
                rhs = std::min(rhs,std::max(0.0,(vxr > 0? 1.0:1.0)*((rhs) - (_s->data[r2][c2]/(_dx*_dx) + UF1D_Slope->Drc * _dx))));

                double cqt = UF1D_COURANTSCHEMEFACTOR * UF_Courant * (volxr);

                double q1 = UF_AVERAGEFACTOR * dtx1 * (vxr) * _s->Drc/_dx;
                q1 = ((q1 > 0)? 1.0 : 0.0) * std::min(std::fabs(q1),(q1 > 0)? cqt : cqt);

                if(SwitchChannelMaxCS)
                {
                    double h_max = UF1D_ChannelMaxCS->Drc;
                    if(h_max > 0.0)
                    {
                         q1 = UF_AVERAGEFACTOR * dtx1 * (vxr) *std::min(sf * h_max,hself)*_lddw->Drc;
                    }
                }

                if(SwitchChannelConnection)
                {
                    double connection = UF1D_ChannelConnected->Drc;
                    if(connection == 0.0)
                    {
                        double volx_r = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * hself);
                        double volxl_r = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->data[r2][c2]) * (_s->data[r2][c2]/(_dx*_lddw->data[r2][c2])));
                        double volx_max = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * _lddh->Drc);
                        double volxl_max = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->data[r2][c2]) *  _lddh->data[r2][c2]);
                        q1 = ((q1 > 0)? 1.0 : -1.0) * std::min(std::fabs(q1),(q1 > 0)? sf * std::max(0.0,volxl_max-volxl_r) : sf * std::max(0.0,volx_max-volx_r));
                    }
                }

                UF1D_sq1->Drc = q1;
            }//else if(_su->Drc > 0)
            {
                double dtx1 = dt->Drc;
                double q1 = UF_BoundaryFlux1D(dtx1,_lddw->Drc,_s->Drc,0,rsu,_fu->Drc,UF1D_Slope->Drc,0.1 + N->Drc, true);
                q1 = ((q1 > 0)? 1.0 : -1.0) * std::min(std::fabs(q1),(q1 > 0)? cqself : cqself);

                UF1D_sq1->Drc = 0;//q1;
            }
        }else
        {
            //if(_su->Drc > 0)
            {
                double dtx1 = dt->Drc;
                double q1 = UF_BoundaryFlux1D(dtx1,_lddw->Drc,_s->Drc,0,rsu,_fu->Drc,UF1D_Slope->Drc,0.1 + N->Drc, true);
                q1 = ((q1 > 0)? 1.0 : 1.0) * std::min(std::fabs(q1),(q1 > 0)? cqself : cqself);


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
                        double volx = UF_LDDOUT(_ldd,r,c,false)?0.0:std::max(0.0,(_dx*_lddw->Drc) * (-(hself - (_s->data[r2][c2]/(_dx*_lddw->data[r2][c2]) + UF1D_Slope->Drc * _dx))));
                        double cq = UF1D_COURANTSCHEMEFACTOR * UF_Courant * volx;
                        double volxl = UF_LDDOUT(_ldd,r,c,false)?0.0:std::max(0.0,(_dx*_lddw->Drc) * (hself - (_s->data[r2][c2]/(_dx*_lddw->data[r2][c2]) + UF1D_Slope->Drc * _dx)));
                        double dtx1 = dt->Drc;
                        double vxl = lsu + dtx1 * lsa;
                        lhs = std::min(lhs,std::max(0.0,(vxl > 0? 1.0:1.0)*((lhs) - (_s->data[r2][c2]/(_dx*_dx) - UF1D_Slope->Drc * _dx))));
                        double cqt = UF1D_COURANTSCHEMEFACTOR * UF_Courant * (volxl);
                        double q2 = UF_AVERAGEFACTOR * dtx1 * (vxl) * lhs *_lddw->Drc* (_lddw->data[r2][c2]/totalwidth);

                        if(SwitchChannelMaxCS)
                        {
                            double h_max = UF1D_ChannelMaxCS->Drc;
                            if(h_max > 0.0)
                            {
                                 q2 = UF_AVERAGEFACTOR * dtx1 * (vxl) *std::min(sf * h_max,hself)*_lddw->Drc;
                            }
                        }

                        q2 = ((q2 > 0)? 0.0 :-1.0) * std::min(std::fabs(q2),(q2 > 0)? cqt : cqt);


                        if(SwitchChannelConnection)
                        {
                            double connection = UF1D_ChannelConnected->Drc;
                            if(connection == 0.0)
                            {
                                double volx_r = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * hself);
                                double volxr_r = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->data[r2][c2]) * (_s->data[r2][c2]/(_dx*_lddw->data[r2][c2])));
                                double volx_max = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->Drc) * _lddh->Drc);
                                double volxr_max = UF_LDDOUT(_ldd,r,c,true)?0.0:std::max(0.0,(_dx*_lddw->data[r2][c2]) *  _lddh->data[r2][c2]);
                                q2 = ((q2 > 0)? 1.0 : -1.0) * std::min(std::fabs(q2),(q2 > 0)? sf* std::max(0.0,volxr_max-volxr_r) :  sf*std::max(0.0,volx_max-volx_r));
                            }
                        }


                        UF1D_sq2->Drc += q2;
                    }
                }
            }
        }else// if(_su->Drc < 0 )
        {
            double dtx2 = dt->Drc;
            double q2 = UF_BoundaryFlux1D(dtx2,_lddw->Drc,_s->Drc,0,lsu,_su->Drc,UF1D_Slope->Drc,0.1 + N->Drc, true);
            q2 = ((q2 < 0)? 1.0 : 1.0) * std::min(std::fabs(q2),(q2 < 0)? cqself : cqself);

            UF1D_sq2->Drc = 0;//q2;
        }
    }}}

}


