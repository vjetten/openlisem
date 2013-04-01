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
#define scheme_type 2   //fixed timetep
#define grav_dem 4.905
#define const_cfl_x 0.5
#define const_cfl_y 0.5


//---------------------------------------------------------------------------
//friction slope
void TWorld::Fr_Manning(double uold, double vold, double hnew, double q1new, double q2new, double dt, double cf)
{
    q1mod = q1new/(1.+cf*cf*grav*sqrt(uold*uold+vold*vold)*dt/pow(hnew,4./3.));
    q2mod = q2new/(1.+cf*cf*grav*sqrt(uold*uold+vold*vold)*dt/pow(hnew,4./3.));
}

//Sf = Manning = v|v|/c^2*h^{4/3}
void TWorld::Fr_ManningSf(double h, double u, double v, double cf)
{
    if (h>he_ca)
    {
        Sf1 = cf*cf*u*sqrt(u*u+v*v)/(pow(h,4./3.));
        Sf2 = cf*cf*v*sqrt(u*u+v*v)/(pow(h,4./3.));
    }
    else
    {
        Sf1 = 0;
        Sf2 = 0;
    }
}
//---------------------------------------------------------------------------
// St Venant equations: conservation of mass and momentum: h u et v
void TWorld::bloc2(TMMap *he, TMMap *ve1, TMMap *ve2, TMMap *qe1, TMMap *qe2,
                   TMMap *hes, TMMap *ves1, TMMap *ves2, TMMap *qes1, TMMap *qes2)
//, TMMap *Vin,double tps, double, int n, double dtheta )
{

    //In case of periodic boundary condition  it's necessary that the flux at the boundary (Left and Right, Bottom and Top) are the same.
    // are the same.
    //moreover we need to consider the direction of the discharge in order to exchange the flux.
    //    if ((nchoice_Tbound==4) && (nchoice_Bbound==4))
    //    {
    //        for (int i=1 ; i<Nxcell+1 ; i++)
    //        {
    //            if ((ve2[i][1] > 0.) && (ve2[i][Nycell]> 0.)){ //the direction of flow is Bottom to the Top
    //                g1[i][1]= g1[i][Nycell+1];
    //            }
    //            if ((ve2[i][1] < 0.) && (ve2[i][Nycell]< 0.)){ //the direction of flow is Top to the Bottom
    //                g1[i][Nycell+1]	= g1[i][1];
    //            }
    //        }
    //    }

    //    if ((nchoice_Rbound==4) && (nchoice_Lbound==4))
    //    {
    //        for (int l=1 ; l<Nycell+1 ; l++){
    //            if ((ve1[1][l] > 0.) && (ve1[Nxcell][l]> 0.)){ //the direction of flow is Left to the Right
    //                f1[1][l]= f1[Nxcell+1][l];
    //            }
    //        }
    //        for (int l=1 ; l<Nycell+1 ; l++){
    //            if ((ve1[1][l] < 0.) && (ve1[Nxcell][l]< 0.)){ //the direction of flow is Right to the Left
    //                f1[Nxcell+1][l] = f1[1][l];
    //            }
    //        }
    //    }
    /*---------------------------------------------------------------------------------------------------------*/
    //Rainfall and infiltration calculated in Saint-Venant system
    //  rain(tps,P);
    //Prain->Rain_func(tps,Rain);
    /*---------------------------------------------------------------------------------------------------------*/

    tx=dt/dx;
    ty=dt/dy;

    // flux' s computation  at each boundaries-----------------------------------------------
    //    if (verif == 1)
    //    {
    //        dt_first=dt;
    //    }

    // --------------------------------------------------------------------------------------

    //    for (int i=1 ; i<Nxcell+1 ; i++)
    //    {
    //        for (int l=1 ; l<Nycell+1 ; l++)
    FOR_ROW_COL_MV_MV
    {
        // Solution of the equation of mass conservation (First equation of Saint venant)
        hes->Drc = he->Drc-tx*(f1->Data[r][c+1]-f1->Drc)-ty*(g1->Data[r+1][c]-g1->Drc)+Rain->Drc*dt;
        //------------------------------------------
        //Infiltration case  ------------------------
        //            if (Kc!=0.)
        //            {

        //                I->calcul(hes->Drc,Vin->Drc,dt,Kc,Ks,dtheta,Psi,zcrust);
        //                hes->Drc = I->get_hmod();
        //                Vin->Drc += I->get_Vin();
        //            }

        //------------------------------------------

        if (hes->Drc > he_ca)
        {
            //Solution of the equation of momentum (Second and third equation of Saint-venant)

            qes1->Drc = he->Drc*ve1->Drc - tx*(f2->Data[r][c+1]-f2->Drc
                    + grav_dem*(h1g->Drc*h1g->Drc -
                                h1l->Drc*h1l->Drc + h1r->Drc*h1r->Drc - h1d->Drc*h1d->Drc
                                + (h1l->Drc+h1r->Drc)*delzc1->Drc))
                    - ty*(g2->Data[r+1][c]-g2->Drc);

            qes2->Drc = he->Drc*ve2->Drc-tx*(f3->Data[r][c+1]-f3->Drc)-ty*(g3->Data[r+1][c]-g3->Drc+grav_dem*(h2g->Drc*h2g->Drc-h2l->Drc*h2l->Drc+h2r->Drc*h2r->Drc-h2d->Drc*h2d->Drc+(h2l->Drc+h2r->Drc)*delzc2->Drc));

            //Calcul friction in semi-implicit.
            Fr_Manning(ve1->Drc,ve2->Drc,hes->Drc,qes1->Drc,qes2->Drc,dt,N->Drc);
            qes1->Drc = q1mod;
            qes2->Drc = q2mod;
            //}
            ves1->Drc = qes1->Drc/hes->Drc;
            ves2->Drc = qes2->Drc/hes->Drc;
        }
        else
        { // Case of height of water is zero.
            ves1->Drc = 0;
            ves2->Drc = 0;

        }
    } //end for l
    //} //end for i

    //  Total_discharge_outflow =  out->boundaries_flux(tps,f1,g1, dt, dt_first,ORDER,verif,Volrain_Tot,Rain);
    //  out->boundaries_flux_LR(tps,f1);
    //  out->boundaries_flux_BT(tps,g1);


}
//---------------------------------------------------------------------------
void TWorld::F_HLL2(double hg, double ug, double vg, double hd, double ud, double vd)
{
    if (hg <= 0 && hd <= 0)
    {
        HLL2_f1 = 0;
        HLL2_f2 = 0;
        HLL2_f3 = 0;
        HLL2_cfl = 0;
    }
    else
    {
        double grav_hg = grav*hg;
        double grav_hd = grav*hd;
        double sqrt_grav_hg = qSqrt(grav_hg);
        double sqrt_grav_hd = qSqrt(grav_hd);
        double qd = ud*hd;
        double qg = ug*hg;
        double c1 = min(ug-sqrt_grav_hg,ud-sqrt_grav_hd); //we already have ug-sqrt_grav_hg<ug+sqrt_grav_hg and ud-sqrt_grav_hd<ud+sqrt_grav_hd
        double c2 = max(ug+sqrt_grav_hg,ud+sqrt_grav_hd); //so we do not need all the eigenvalues to get c1 and c2
        double tmp = 1/(c2-c1);
        double t1 = (min(c2,0)-min(c1,0))*tmp;
        double t2 = 1-t1;
        double t3 = (c2*fabs(c1)-c1*fabs(c2))*0.5*tmp;

        HLL2_f1 = t1*qd+t2*qg-t3*(hd-hg);
        HLL2_f2 = t1*(qd*ud+grav_hd*hd*0.5)+t2*(qg*ug+grav_hg*hg*0.5)-t3*(qd-qg);
        HLL2_f3 = t1*qd*vd+t2*qg*vg-t3*(hd*vd-hg*vg);
        HLL2_cfl = max(fabs(c1),fabs(c2)); //cfl is the velocity to compute the cfl condition max(fabs(c1),fabs(c2))*tx with tx=dt/dx
    }
}
//---------------------------------------------------------------------------
/*Construction varibles for hydrostatic reconstruction.
 Flux with x and y.
 Calculaton of the time steps in relation to fixed cfl.*/
void TWorld::bloc1(double cflfix, double T, double tps,
                   double dt_max, double dt, double &dt_cal, double cfl_new)
{
    double dt_tmp, dtx, dty;
    double velocity_max_x, velocity_max_y;
    //tempory velocity to verify if clf > cflfix
    dtx = dt_max;
    dty = dt_max;
    velocity_max_x = -ve_ca;
    velocity_max_y = -ve_ca;

    FOR_ROW_COL_MV_MV
    {
        h1d->Data[r][c-1] = max(0, h1r->Data[r][c-1] - max(0, delz1->Data[r][c-1]));
        h1g->Drc = max(0, h1l->Drc - max(0, -delz1->Data[r][c-1]));

        F_HLL2(h1d->Data[r][c-1],u1r->Data[r][c-1],v1r->Data[r][c-1],
                h1g->Drc,u1l->Drc,v1l->Drc);
        f1->Drc = HLL2_f1;
        f2->Drc = HLL2_f2;
        f3->Drc = HLL2_f3;

        if (HLL2_cfl*dt/dx == 0)
            dt_tmp = dt_max;
        else
            dt_tmp = cflfix*dx/HLL2_cfl;
        dtx = min(min(dt, dt_tmp), dtx);
        velocity_max_x = max(velocity_max_x, HLL2_cfl);
    }

    FOR_ROW_COL_MV_MV
    {
        h2d->Data[r-1][c] = max(0, h2r->Data[r-1][c] - max(0, delz2->Data[r-1][c]));
        h2g->Drc = max(0, h2l->Drc - max(0, -delz2->Data[r-1][c]));

        F_HLL2(h2d->Data[r-1][c],v2r->Data[r-1][c],u2r->Data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);
        g1->Drc = HLL2_f1;
        g2->Drc = HLL2_f2;   //LET OP !!!!!!!!!!!!!!!! is omgekeerd
        g3->Drc = HLL2_f3;

        if (HLL2_cfl*dt/dy == 0)
            dt_tmp = dt_max;
        else
            dt_tmp = cflfix*dy/HLL2_cfl;
        dty = min(min(dt, dt_tmp), dty);
        velocity_max_y = max(velocity_max_y, HLL2_cfl);
    } //end for l


    if (scheme_type == 1)
        dt_cal=min(dtx,dty);
    else
    {
        if ((velocity_max_x*dt_fix/dx > cflfix)||(velocity_max_y*dt_fix/dy>cflfix)){
            qDebug() << "the CFL condition is not satisfied: CFL >"<<cflfix << endl;
            exit(1);
        }
        dt_cal=dt_fix;
    }

}
//---------------------------------------------------------------------------
void TWorld::initMUSCL(TMMap *z)
{
    //	this->Nxcell = z.size()-2;
    //	this->Nycell = z[1].size()-2;

    //	z1r.resize(Nxcell+1); // i : 0->Nxcell
    //	z1l.resize(Nxcell+2); // i : 1->Nxcell+1
    //	z2r.resize(Nxcell+1); // i : 1->Nxcell
    //	z2l.resize(Nxcell+1); // i : 1->Nxcell
    //	delta_z1.resize(Nxcell+1); // i : 0->Nxcell
    //	delta_z2.resize(Nxcell+1); // i : 1->Nxcell

    //	z1r[0].resize(Nycell+1); // l : 1->Nycell
    //	delta_z1[0].resize(Nycell+1); // l : 1->Nycell

    //	const int nodex=Nxcell;
    //	const int nodey=Nycell;

    //	for (int i=1 ; i<=nodex ; i++)
    //    {
    //		z1r[i].resize(Nycell+1); // l : 1->Nycell
    //		z1l[i].resize(Nycell+1); // l : 1->Nycell
    //		z2r[i].resize(Nycell+1); // l : 0->Nycell
    //		z2l[i].resize(Nycell+2); // l : 1->Nycell+1
    //		delta_z1[i].resize(Nycell+1); // l : 1->Nycell
    //		delta_z2[i].resize(Nycell+1); // l : 0->Nycell
    //	}

    //	z1l[Nxcell+1].resize(Nycell+1);

    //	for (int l=1 ; l<=nodey ; l++)
    //    {
    //	  z1r[0][l] = z[0][l];
    //	  z1l[Nxcell+1][l] = z[Nxcell+1][l];
    //	  delta_z1[0][l] = z[1][l]-z[0][l];
    //	} //end for l

    //	for (int i=1 ; i<=nodex ; i++)
    //    {
    //	  z2r[i][0] = z[i][0];
    //	  z2l[i][Nycell+1] = z[i][Nycell+1];
    //	  delta_z2[i][0] = z[i][1]-z[i][0];
    //	  for (int l=1 ; l<=nodey ; l++)
    //      {
    //	    delta_z1->Drc = z->Data[r][c+1]-z->Drc;
    //	    delta_z2->Drc = z->Data[r+1][c]-z->Drc;
    //	  } //end for l
    //	} //end for i
}
//---------------------------------------------------------------------------
double TWorld::limiter(double a, double b)
{
    if (a>=0. && b>=0)
        return(min(a,b));
    else
        if (a<=0. && b<=0)
            return(max(a,b));
        else
            return(0);
}
//---------------------------------------------------------------------------
//reconstruction scheme
//original matrix is col, row orientation (y,x), NOT row, col (x,y)
//1 is in the y direction, 2 = x direction
//TO DO: only in flood domain
//called as
//rec->calcul(h,u,v,z,delzc1,delzc2,delz1,delz2,h1r,u1r,v1r,h1l,u1l,v1l,h2r,u2r,v2r,h2l,u2l,v2l);

void TWorld::MUSCL(TMMap *h,TMMap *u,TMMap *v,TMMap *z,
                   TMMap *delzc1,TMMap *delzc2,TMMap *delz1,TMMap *delz2,
                   TMMap *h1r,TMMap *u1r,TMMap *v1r,TMMap *h1l,TMMap *u1l,TMMap *v1l,
                   TMMap *h2r,TMMap *u2r,TMMap *v2r,TMMap *h2l,TMMap *u2l,TMMap *v2l)
{
    double delta_h1, delta_u1, delta_v1;
    double delta_h2, delta_u2, delta_v2;
    double dh, du, dv, dz_h;
    delta_h2 = 0;
    delta_u2 = 0;
    delta_v2 = 0;
    delta_h1 = 0;
    delta_u1 = 0;
    delta_v1 = 0;

    //    h1->shift(h, 0, 1);  //row, col
    //    u1->shift(u, 0, 1);  //row, col
    //    v1->shift(v, 0, 1);  //row, col
    //    z1->shift(z, 0, 1);
    //    z2->shift(z, 1, 0);

    FOR_ROW_COL_MV_MV
    {
        //        delta_z1->Drc = z->Data[r][c+1] - z->Drc;
        //        delta_z2->Drc = z->Data[r+1][c] - z->Drc;
        delta_z1->Drc = z->Data[r+1][c] - z->Drc;
        delta_z2->Drc = z->Data[r][c+1] - z->Drc;
    }


    // first do 1: x direction (col+1 - col)
    FOR_ROW_COL_MV_MV
    {

        delta_h2 = h->Data[r+1][c] - h->Drc;
        delta_u2 = u->Data[r+1][c] - u->Drc;
        delta_v2 = v->Data[r+1][c] - v->Drc;

        if (IS_MV_REAL8(&LDD->Data[r][c+1]))
        {
            delta_h1 = delta_h2;
            delta_u1 = delta_u2;
            delta_v1 = delta_v2;
        }

        dh = limiter(delta_h1, delta_h2);
        dz_h = limiter(delta_h1+delta_z1->Data[r][c-1], delta_h2+delta_z1->Drc);

        du = limiter(delta_u1, delta_u2);
        dv = limiter(delta_v1, delta_v2);

        h1r->Drc = h->Drc+dh*0.5;
        h1l->Drc = h->Drc-dh*0.5;

        z1r->Drc = z->Drc+0.5*(dz_h-dh);
        z1l->Drc = z->Drc+0.5*(dh-dz_h);
        delzc1->Drc = z1r->Drc-z1l->Drc;
        delz1->Data[r][c-1] = z1l->Drc-z1r->Data[r][c-1];

        if (h->Drc > 0)
        {
            u1r->Drc = u->Drc+h1l->Drc*du*0.5/h->Drc;
            u1l->Drc = u->Drc-h1r->Drc*du*0.5/h->Drc;
            v1r->Drc = v->Drc+h1l->Drc*dv*0.5/h->Drc;
            v1l->Drc = v->Drc-h1r->Drc*dv*0.5/h->Drc;
        }
        else
        {
            u1r->Drc = u->Drc+du*0.5;
            u1l->Drc = u->Drc-du*0.5;
            v1r->Drc = v->Drc+dv*0.5;
            v1l->Drc = v->Drc-dv*0.5;
        }
        //rajout de l'iteration des delta
        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;

    }


    //    h1->shift(h, 1, 0);  //row, col
    //    u1->shift(u, 1, 0);  //row, col
    //    v1->shift(v, 1, 0);  //row, col

    // then do 2: y direction (row+1 - row)
    FOR_ROW_COL_MV_MV
    {
        delta_h2 = h->Data[r][c+1] - h->Drc;
        delta_u2 = u->Data[r][c+1] - u->Drc;
        delta_v2 = v->Data[r][c+1] - v->Drc;

        if (IS_MV_REAL8(&LDD->Data[r+1][c]))
        {
            delta_h1 = delta_h2;
            delta_u1 = delta_u2;
            delta_v1 = delta_v2;
        }

        dh = limiter(delta_h1, delta_h2);
        dz_h = limiter(delta_h1+delta_z2->Data[r-1][c], delta_h2+delta_z2->Drc);

        du = limiter(delta_u1, delta_u2);
        dv = limiter(delta_v1, delta_v2);

        h2r->Drc = h->Drc+dh*0.5;
        h2l->Drc = h->Drc-dh*0.5;

        z2r->Drc = z->Drc+0.5*(dz_h-dh);
        z2l->Drc = z->Drc+0.5*(dh-dz_h);
        delzc2->Drc = z2r->Drc-z2l->Drc;
        delz2->Data[r-1][c] = z2l->Drc - z2r->Data[r-1][c];

        if (h->Drc > 0)
        {
            u2r->Drc = u->Drc+h2l->Drc*du*0.5/h->Drc;
            u2l->Drc = u->Drc-h2r->Drc*du*0.5/h->Drc;
            v2r->Drc = v->Drc+h2l->Drc*dv*0.5/h->Drc;
            v2l->Drc = v->Drc-h2r->Drc*dv*0.5/h->Drc;
        }
        else
        {
            u2r->Drc = u->Drc+du*0.5;
            u2l->Drc = u->Drc-du*0.5;
            v2r->Drc = v->Drc+dv*0.5;
            v2l->Drc = v->Drc-dv*0.5;
        } //end if

        //rajout de l'iteration des delta
        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;

    }
}

//---------------------------------------------------------------------------
void TWorld::boundary()
{

}
//---------------------------------------------------------------------------
void TWorld::fullSWOF2D(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2)
{
    //initialization

    double dt2;

    double cfl_new = cfl_fix;

    bool start = true;
    if (start)
    {
        //dt1=dt_max;

        //boundary conditions
        //     boundary(h,u,v,tps,Nxcell,Nycell);


        FOR_ROW_COL_MV
        {
            if (h->Drc <= he_ca)
            {
                h->Drc = 0;
                u->Drc = 0;
                v->Drc = 0;
            }
            if (fabs(u->Drc) <= ve_ca)
            {
                u->Drc = 0;
                q1->Drc = 0;
            }
            if (fabs(v->Drc) <= ve_ca)
            {
                v->Drc = 0;
                q2->Drc = 0;
            }
        }

        // Reconstruction for order 2
        MUSCL(h,u,v,z,delzc1,delzc2,delz1,delz2,h1r,u1r,v1r,h1l,u1l,v1l,h2r,u2r,v2r,h2l,u2l,v2l);

    }


    bloc1(cfl_fix, T, tps, dt_max, dt1, dt1, cfl_new);

    dt1=min(T-tps,dt1);

    tx=dt1/dx;
    ty=dt1/dy;

    bloc2(h,u,v,q1,q2,hs,us,vs,qs1,qs2);//Vin1,tps,dt1);//,n,dtheta);
    //    dt2=dt_max;
    dt2=dt1;

    //boundary conditions
    //boundary(hs,us,vs,tps+dt1,Nxcell,Nycell);

    //VJ set all < 1e-12 to zero
    FOR_ROW_COL_MV
    {
        if (hs->Drc<=he_ca)
        {
            hs->Drc=0;
            us->Drc = 0;
            vs->Drc = 0;
            qs1->Drc = 0;
            qs2->Drc = 0;
        }
        if (fabs(us->Drc)<=ve_ca)
        {
            us->Drc = 0;
            qs1->Drc = 0;
        }
        if (fabs(vs->Drc)<=ve_ca)
        {
            vs->Drc = 0;
            qs2->Drc = 0;
        }
    }


    //Reconstruction for order 2
    MUSCL(hs,us,vs,z,delzc1,delzc2,delz1,delz2,h1r,u1r,v1r,h1l,u1l,v1l,h2r,u2r,v2r,h2l,u2l,v2l);


    //commun bloc

    bloc1(cfl_fix,T,tps+dt1,dt_max,dt2,dt2,cfl_new);

    if (dt2<dt1)
    {
        dt1=dt2;
        tx=dt1/dx;
        ty=dt1/dy;
        start = false;
    }
    else
    {

        //Added to do calculus at the beginning
        start = true;
        bloc2(hs,us,vs,qs1,qs2,hsa,usa,vsa,qsa1,qsa2,Vin2,tps+dt1,dt1,n+1,dtheta);



        //Heun method (order 2 in time)
        FOR_ROW_COL_MV
        {
            if (hsa->Drc <= he_ca)
                hsa->Drc=0;
            if (abs(usa->Drc) <= ve_ca)
                usa->Drc=0;
            if (abs(vsa->Drc) <= ve_ca)
                vsa->Drc=0;
            tmp = 0.5*(h->Drc+hsa->Drc);
            if (tmp>=he_ca)
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
            Vin_tot->Drc = (Vin1->Drc + Vin2->Drc)*0.5;

            Vin1->Drc=Vin_tot->Drc;
            Vin2->Drc=Vin_tot->Drc;

        } //end for i

        tps=tps+dt1;
        n=n+1;

        out->check_vol(tps,dt1,h,Vin_tot,Volrain_Tot,Vol_inf,Vol_of,Total_discharge_outflow);

        //Display the percentage of elapsed time
        cout  << '\r' << '\t' << "[" << int((tps/T)*100) << "%] done" ;

    }//end for else dt2<dt1



} //end for while : loop in time


//to give the ultimate situation
out->write(h,u,v,z,tps);


//to give the computing time
time(&end);
timecomputation=difftime(end,start);

//to inform about the froude number
Fr=froude_number(hs,us,vs);


// The quantity of water outflow is the sum of flux at the boundary multiply by one cell area, so
//Outflow volum = (fluxNxcell_cum_T+fluxNycell_cum_T)*dx*dy*n1+(fluxx0_cum_T+fluxy0_cum_T)*dx*dy*n2
//In this case n1=1 and n2=-1 because we consider the direction of the flow
// is from left to the right (x=0 to x=Nxcell) and from bottom to the top (y=0 to y=Nycelle)
out->result(timecomputation,Volrain_Tot, Vol_inf, Vol_of,Fr,n,Total_discharge_outflow);
//storage of h u v value in the final time
out->final(dx, dy,z, h, u,v);
*/
}
//---------------------------------------------------------------------------

