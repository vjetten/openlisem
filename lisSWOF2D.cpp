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
 \file lisSWOF2d.cpp
 \brief Channel flood using the opensource code from the FullSWOF_2D project.
        St Vennat equations
        authors: Olivier Delestre, Christian Laguerre, Carine Lucas, Ulrich Razafison, Marie Rousseau.
        website: http://www.univ-orleans.fr/mapmo/soft/FullSWOF/

functions: \n
- void TWorld::ChannelFlood(void) calculate maps channelflood height (hmx) and FloodDomain
*/

#include "lisemqt.h"
#include "model.h"
#include "global.h"

#define he_ca 1e-12
#define ve_ca 1e-12
#define cfl_fix 0.4
#define grav 9.81
#define scheme_type 1


//---------------------------------------------------------------------------
// St Venant equations: conservation of mass and momentum:
// h u et v go in, new hs, vs, us go out
void TWorld::bloc2(double dt, TMMap *he, TMMap *ve1, TMMap *ve2,
                   TMMap *hes, TMMap *ves1, TMMap *ves2)
{
    tx=dt/dx;
    ty=dt/dy;

    FOR_ROW_COL_MV_MV
            if (FloodDomain->Drc == 1)
    {
        // Solution of the equation of mass conservation (First equation of Saint venant)
        // f1 comes from MUSCL calculations
        hes->Drc = he->Drc-tx*(f1->Data[r][c+1]-f1->Drc)-ty*(g1->Data[r+1][c]-g1->Drc);

        if (hes->Drc > he_ca)
        {
            //Solution of the equation of momentum (Second and third equation of Saint-venant)
            double qes1;
            double qes2;

            qes1 = he->Drc*ve1->Drc - tx*(f2->Data[r][c+1]-f2->Drc
                    + grav*0.5*(h1g->Drc*h1g->Drc -
                                h1l->Drc*h1l->Drc + h1r->Drc*h1r->Drc - h1d->Drc*h1d->Drc
                                + (h1l->Drc+h1r->Drc)*delzc1->Drc))
                    - ty*(g2->Data[r+1][c]-g2->Drc);

            qes2 = he->Drc*ve2->Drc -
                    tx*(f3->Data[r][c+1]-f3->Drc) -
                    ty*(g3->Data[r+1][c]-g3->Drc + grav*0.5*(h2g->Drc*h2g->Drc - h2l->Drc*h2l->Drc + h2r->Drc*h2r->Drc - h2d->Drc*h2d->Drc
                                                             + (h2l->Drc+h2r->Drc)*delzc2->Drc));
            //Calcul friction in semi-implicit.
            Fr_Manning(ve1->Drc, ve2->Drc, hes->Drc, qes1, qes2, dt, N->Drc);
            ves1->Drc = q1mod/hes->Drc;
            ves2->Drc = q2mod/hes->Drc;
        }
        else
        {
            // Case of height of water is zero.
            ves1->Drc = 0;
            ves2->Drc = 0;
        }
    }
}

//---------------------------------------------------------------------------
/*Construction varibles for hydrostatic reconstruction.
 Flux with x and y.
 Calculaton of the time steps in relation to fixed cfl.*/
double TWorld::bloc1(double dt, double dt_max)
{
    double dt_tmp, dtx, dty;
    double velocity_max_x, velocity_max_y;
    //tempory velocity to verify if clf > cflfix
    dtx = dt_max;
    dty = dt_max;
    velocity_max_x = -ve_ca;
    velocity_max_y = -ve_ca;

    FOR_ROW_COL_MV_MV
            if (FloodDomain->Drc == 1)
    {
        h1d->Data[r][c-1] = max(0, h1r->Data[r][c-1] - max(0, delz1->Data[r][c-1]));
        h1g->Drc = max(0, h1l->Drc - max(0, -delz1->Data[r][c-1]));

        // this rusanov, HLL or HLL2
        //        F_Rusanov(h1d->Data[r][c-1], u1r->Data[r][c-1], v1r->Data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
        //        F_HLL(h1d->Data[r][c-1], u1r->Data[r][c-1], v1r->Data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
        F_HLL2(h1d->Data[r][c-1], u1r->Data[r][c-1], v1r->Data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
        f1->Drc = HLL2_f1;
        f2->Drc = HLL2_f2;
        f3->Drc = HLL2_f3;

        if (HLL2_cfl*dt/dx == 0)
            dt_tmp = dt_max;
        else
            dt_tmp = cfl_fix*dx/HLL2_cfl;
        dtx = min(min(dt, dt_tmp), dtx);
        velocity_max_x = max(velocity_max_x, HLL2_cfl);
    }

    FOR_ROW_COL_MV_MV
            if (FloodDomain->Drc == 1)
    {
        h2d->Data[r-1][c] = max(0, h2r->Data[r-1][c] - max(0, delz2->Data[r-1][c]));
        h2g->Drc = max(0, h2l->Drc - max(0, -delz2->Data[r-1][c]));

        //F_HLL2
        //         F_Rusanov(h2d->Data[r-1][c],v2r->Data[r-1][c],u2r->Data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);
        //        F_HLL(h1d->Data[r][c-1], u1r->Data[r][c-1], v1r->Data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
        F_HLL2(h1d->Data[r][c-1], u1r->Data[r][c-1], v1r->Data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
        g1->Drc = HLL2_f1;
        g2->Drc = HLL2_f2;   //LET OP !!!!!!!!!!!!!!!! is omgekeerd
        g3->Drc = HLL2_f3;

        if (HLL2_cfl*dt/dy == 0)
            dt_tmp = dt_max;
        else
            dt_tmp = cfl_fix*dy/HLL2_cfl;
        dty = min(min(dt, dt_tmp), dty);
        velocity_max_y = max(velocity_max_y, HLL2_cfl);
    } //end for l

    qDebug() << dtx << dty << dt << dt_tmp << HLL2_cfl << velocity_max_x << velocity_max_y;
    if (scheme_type == 1)
        return(min(dtx,dty));
    else
    {
        //        if ((velocity_max_x*dt_fix/dx > cflfix)||(velocity_max_y*dt_fix/dy > cflfix)){
        //            qDebug() << "the CFL condition is not satisfied: CFL >"<<cflfix << endl;
        //            exit(1);
        //        }
        //        return (dt_fix);
    }

}
//---------------------------------------------------------------------------
// 1st order solution
double TWorld::fullSWOF2Do1(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2)
{
    double timesum = 0;
    int n = 0;

    dx = _dx;
    dy = _dx; //TODO could be DX->Drc ? later !!!!

    dt_max = min(dx/2,dy/2);
    dt_max = 0.8*_dx/qSqrt(grav*h->mapMaximum());

    dt1 = dt_max;

    bool startflood = false;
    tma->copy(hmx);
    FOR_ROW_COL_MV
    {
        if (hmx->Drc > 0)
        {
            startflood = true;
            break;
        }
    }

    // if there is no flood skip everything
    if (startflood)
    {
        FOR_ROW_COL_MV_MV
                delz1->Data[r][c-1] = z->Drc - z->Data[r][c-1];
        FOR_ROW_COL_MV_MV
                delz2->Data[r-1][c] = z->Drc - z->Data[r-1][c];

        do {
            n++;
            //boundary conditions
            //boundary(h,u,v,tps,Nxcell,Nycell);

            setZero(h, u, v, q1, q2);

            FOR_ROW_COL_MV_MV
                    if (FloodDomain->Drc == 1)
            {
                h1r->Drc = h->Drc;
                u1r->Drc = u->Drc;
                v1r->Drc = v->Drc;
                h1l->Data[r][c+1] = h->Data[r][c+1];
                u1l->Data[r][c+1] = u->Data[r][c+1];
                v1l->Data[r][c+1] = v->Data[r][c+1];
            }
            FOR_ROW_COL_MV_MV
                    if (FloodDomain->Drc == 1)
            {
                h2r->Drc = h->Drc;
                u2r->Drc = u->Drc;
                v2r->Drc = v->Drc;
                h2l->Data[r+1][c] = h->Data[r+1][c];
                u2l->Data[r+1][c] = u->Data[r+1][c];
                v2l->Data[r+1][c] = v->Data[r+1][c];
            }

            double dta = bloc1(dt1, dt_max);

            dt1 = qMin(dta, _dt-timesum);
            qDebug() << "bloc1" << dt1 << dta << dt_max << timesum;

            bloc2(dt1, h,u,v, hs,us,vs);

            if (dt1 <= 0)
                exit(0);

            FOR_ROW_COL_MV_MV
            {
                h->Drc = hs->Drc;
                u->Drc = us->Drc;
                v->Drc = vs->Drc;
                //                u->Drc = min(20, us->Drc);
                //                v->Drc = min(20, vs->Drc);
                q1->Drc = h->Drc*u->Drc;
                q2->Drc = h->Drc*v->Drc;
            }

            setZero(h, u, v, q1, q2);

            // Calcul_Volume_inf_ruis_rain(dt1,h,Rain,Nxcell+1,Nycell+1);

            // out->check_vol(tps,dt1,h,Vin_tot,Volrain_Tot,Vol_inf,Vol_of,Total_discharge_outflow);

            timesum = timesum + dt1;
            // sum to reach _dt

        } while (timesum  < _dt);

        FOR_ROW_COL_MV
        {
            if (hmx->Drc > 0)
                FloodDomain->Drc = 1;
            else
                FloodDomain->Drc = 0;
        }
    }

    // continue while _dt is not reached

    //        //to give the ultimate situation
    //        out->write(h,u,v,z,tps);

    //        //to give the computing time
    //        time(&end);
    //        timecomputation=difftime(end,start);

    //        //to inform about the froude number
    //        Fr=froude_number(hs,us,vs);

    //        // The quantity of water outflow is the sum of flux at the boundary multiply by one cell area, so
    //        //Outflow volum = (fluxy0_cum_T+fluxNycell_cum_T+fluxx0_cum_T+fluxNxcell_cum_T)*dx*dy
    //        //In this case n1=1 and n2=-1 because we consider the direction of the flow
    //        // is from left to the right (x=0 to x=Nxcell) and from bottom to the top (y=0 to y=Nycelle)
    //        out->result(timecomputation,Volrain_Tot, Vol_inf, Vol_of,Fr,n,Total_discharge_outflow);

    //        //storage of h u v value in the final time
    //        out->final(dx, dy,z, h, u,v);
    tm->calc2Maps(v,u,ADD);
    tm->report("vv");
    v->report("v");
    u->report("u");
    q1->report("qa");
    q2->report("qb");

    return(timesum/(n+1));
}
//---------------------------------------------------------------------------
// 2nd order solution
double TWorld::fullSWOF2D(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2)
{
    double dt2;
    double timesum = 0;
    int n = 1;

    if (prepareFlood)
    {
        verif =1;
    }

    dx = _dx;
    dy = _dx; //TODO could be DX->Drc ? later !!!!

    dt_max = min(dx/2,dy/2);

    dt_max = 0.8*_dx/qSqrt(grav*h->mapMaximum());

    bool startflood = false;
    FOR_ROW_COL_MV
    {
        if (hmx->Drc > 0)
        {
            startflood = true;
            break;
        }
    }

    // if there is no flood skip everything
    if (startflood)
    {
        do {
            n++;
            if (verif == 1)
            {
                dt1 = dt_max;

                FOR_ROW_COL_MV
                {
                    q1->Drc = u->Drc*h->Drc;
                    q2->Drc = v->Drc*h->Drc;
                }

                setZero(h, u, v, q1, q2);

                // Reconstruction for order 2
                ENO(h,u,v,z,delzc1,delzc2,delz1,delz2,h1r,u1r,v1r,h1l,u1l,v1l,h2r,u2r,v2r,h2l,u2l,v2l);
            }
            dt1 = bloc1(dt1, dt_max);

            tx=dt1/dx;
            ty=dt1/dy;
            //h, u, v go in hs, vs, us come out
            bloc2(dt1, h,u,v, hs,us,vs);
            dt2 = dt1;

            //VJ set all < 1e-12 to zero
            setZero(hs, us, vs, qs1, qs2);

            //Reconstruction for order 2
            ENO(hs,us,vs,z,delzc1,delzc2,delz1,delz2,h1r,u1r,v1r,h1l,u1l,v1l,h2r,u2r,v2r,h2l,u2l,v2l);

            dt2 = bloc1(dt2, dt_max);

            if (dt2 < dt1)
            {
                dt1 = dt2;
                tx=dt1/dx;
                ty=dt1/dy;
                verif = 0;
            }
            //            else
            //          {
            verif = 1;
            //hs, us, vs go in hsa, vsa, usa come out
            bloc2(dt1, hs,us,vs, hsa,usa,vsa);

            //Heun method (order 2 in time)
            FOR_ROW_COL_MV
            {
                if (hsa->Drc <= he_ca)
                    hsa->Drc = 0;
                //if (qAbs(usa->Drc) <= ve_ca)
                if (usa->Drc <= ve_ca)
                    usa->Drc = 0;
                //if (qAbs(vsa->Drc) <= ve_ca)
                if (vsa->Drc <= ve_ca)
                    vsa->Drc = 0;

                double tmp = 0.5*(h->Drc+hsa->Drc);
                if (tmp >= he_ca)
                {
                    q1->Drc = 0.5*(h->Drc*u->Drc+hsa->Drc*usa->Drc);
                    u->Drc = q1->Drc/tmp;
                    q2->Drc = 0.5*(h->Drc*v->Drc+hsa->Drc*vsa->Drc);
                    v->Drc = q2->Drc/tmp;
                    h->Drc = tmp;
                }
                else
                {
                    u->Drc = 0;
                    q1->Drc = 0;
                    v->Drc = 0;
                    q2->Drc = 0;
                    h->Drc = 0;
                }
            }
            //        }//end for else dt2<dt1

            dt1 = qMin(dt1, _dt-timesum);

            timesum = timesum + dt1;

        } while (timesum  < _dt);

        FOR_ROW_COL_MV
        {
            if (hmx->Drc > 0)
                FloodDomain->Drc = 1;
            else
                FloodDomain->Drc = 0;
        }
    } // if floodstart

    u->report("u");
    v->report("v");
    q1->report("qa");
    q2->report("qb");

    return(timesum/(n+1));

}
