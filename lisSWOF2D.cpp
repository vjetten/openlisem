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
#define dt_ca 1e-6
#define grav 9.81
#define scheme_type 1

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
void TWorld::setZero(TMMap *h, TMMap *u, TMMap *v)
{
    FOR_ROW_COL_MV
    {
        if (h->Drc <= he_ca)
        {
            h->Drc = 0;
            u->Drc = 0;
            v->Drc = 0;
        }
        //       if (fabs(u->Drc) <= ve_ca)
        if (u->Drc <= ve_ca)
        {
            u->Drc = 0;
            //q1->Drc = 0;
        }
        //if (fabs(v->Drc) <= ve_ca)
        if (v->Drc <= ve_ca)
        {
            v->Drc = 0;
            //q2->Drc = 0;
        }
    }
    h->DrcOutlet = 0;
    v->DrcOutlet = 0;
    u->DrcOutlet = 0;
}
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
void TWorld::F_HLL2(double hg, double ug, double vg, double hd, double ud, double vd)
{
    double f1, f2, f3, cfl;
    if (hg<=0. && hd<=0)
    {
        f1 = 0;
        f2 = 0;
        f3 = 0;
        cfl = 0;
    }
    else
    {
        double grav_hg = grav*hg;
        double grav_hd = grav*hd;
        double sqrt_grav_hg = sqrt(grav_hg);  //wave velocity
        double sqrt_grav_hd = sqrt(grav_hd);
        double qd = ud*hd;
        double qg = ug*hg;
        double c1 = min(ug-sqrt_grav_hg,ud-sqrt_grav_hd); //we already have ug-sqrt_grav_hg<ug+sqrt_grav_hg and ud-sqrt_grav_hd<ud+sqrt_grav_hd
        double c2 = max(ug+sqrt_grav_hg,ud+sqrt_grav_hd); //so we do not need all the eigenvalues to get c1 and c2
        double tmp = 1./(c2-c1);
        double t1 = (min(c2,0)-min(c1,0))*tmp;
        double t2 = 1.-t1;
        double t3 = (c2*fabs(c1)-c1*fabs(c2))*0.5*tmp;

        f1 = t1*qd+t2*qg-t3*(hd-hg);
        f2 = t1*(qd*ud+grav_hd*hd*0.5)+t2*(qg*ug+grav_hg*hg*0.5)-t3*(qd-qg);
        f3 = t1*qd*vd+t2*qg*vg-t3*(hd*vd-hg*vg);
        cfl = max(fabs(c1),fabs(c2)); //cfl is the velocity to compute the cfl condition max(fabs(c1),fabs(c2))*tx with tx=dt/dx
    }
    HLL2_cfl = cfl;
    HLL2_f1 = f1;
    HLL2_f2 = f2;
    HLL2_f3 = f3;
}

void TWorld::F_HLL(double hg,double ug,double vg,double hd,double ud,double vd)
{
    double f1, f2, f3, cfl;
    if (hg<=0. && hd<=0)
    {
        f1 = 0;
        f2 = 0;
        f3 = 0;
        cfl = 0;
    }
    else
    {
        double grav_hg = grav*hg;
        double grav_hd = grav*hd;
        double qd = ud*hd;
        double qg = ug*hg;
        double c1 = min(ug-sqrt(grav_hg),ud-sqrt(grav_hd));
        double c2 = max(ug+sqrt(grav_hg),ud+sqrt(grav_hd));

        //cfl is the velocity to calculate the real cfl=max(fabs(c1),fabs(c2))*tx with tx=dt/dx
        if (c1==0. && c2==0)
        {
            //dry
            f1=0;
            f2=0;
            f3=0;
            cfl=0; //max(fabs(c1),fabs(c2))=0
        }
        else
            if (c1>=0){
                //supercritical flow, from left to right : we have max(abs(c1),abs(c2))=c2>0
                f1=qg;
                f2=qg*ug+grav*hg*hg*0.5;
                f3=qg*vg;
                cfl=c2; //max(fabs(c1),fabs(c2))=c2>0
            }
            else
                if (c2<=0)
                {
                    //supercritical flow, from right to left : we have max(abs(c1),abs(c2))=-c1>0
                    f1=qd;
                    f2=qd*ud+grav*hd*hd*0.5;
                    f3=qd*vd;
                    cfl=fabs(c1); //max(fabs(c1),fabs(c2))=fabs(c1)
                }
                else
                { //subcritical flow
                    double tmp = 1./(c2-c1);
                    f1=(c2*qg-c1*qd)*tmp+c1*c2*(hd-hg)*tmp;
                    f2=(c2*(qg*ug+grav*hg*hg*0.5)-c1*(qd*ud+grav*hd*hd*0.5))*tmp+c1*c2*(qd-qg)*tmp;
                    f3=(c2*(qg*vg)-c1*(qd*vd))*tmp+c1*c2*(hd*vd-hg*vg)*tmp;
                    cfl = max(fabs(c1),fabs(c2));
                }
    }
    HLL2_cfl = cfl;
    HLL2_f1 = f1;
    HLL2_f2 = f2;
    HLL2_f3 = f3;
}


void TWorld::F_Rusanov(double hg,double ug,double vg,double hd,double ud,double vd)
{
    double c;
    double f1, f2, f3, cfl;
    if (hg<=0. && hd<=0)
    {
        c = 0;
        f1 = 0;
        f2 = 0;
        f3 = 0;
        cfl = 0;
    }
    else
    {
        c = max(fabs(ug)+sqrt(grav*hg),fabs(ud)+sqrt(grav*hd));
        double cd = c*0.5;
        double qd = ud*hd;
        double qg = ug*hg;
        f1 = (qg+qd)*0.5-cd*(hd-hg);
        f2 = ((ug*qg)+(grav*0.5*hg*hg)+(ud*qd)+(grav*0.5*hd*hd))*0.5-cd*(qd-qg);
        f3 = (qg*vg+qd*vd)*0.5-cd*(hd*vd-hg*vg);
        cfl = c;
    }
    HLL2_cfl = cfl;
    HLL2_f1 = f1;
    HLL2_f2 = f2;
    HLL2_f3 = f3;
}
//---------------------------------------------------------------------------

void TWorld::ENO(TMMap *h,TMMap *u,TMMap *v,TMMap *z)
//,TMMap *delzc1,TMMap *delzc2,TMMap *delz1,TMMap *delz2,TMMap *h1r,TMMap *u1r,TMMap *v1r,TMMap *h1l,TMMap *u1l,TMMap *v1l,TMMap *h2r,TMMap *u2r,TMMap *v2r,TMMap *h2l,TMMap *u2l,TMMap *v2l)
{
    double ddh1 = 0;
    double ddz1 = 0;
    double ddu1 = 0;
    double ddv1 = 0;
    double ddh2, ddv2, ddu2, ddz2;
    double uu1 = 0, vv1 = 0, hh1 = 0;
    double uu2, vv2, hh2;
    double delta_h1, delta_v1, delta_u1;
    double delta_h2, delta_v2, delta_u2;
    double dh, du, dv, dz_h;
    bool startloop = true;
    bool startloop1 = true;
    double amortENO = 0.25;
    delta_h1 = 0;
    delta_u1 = 0;
    delta_v1 = 0;

    if (prepareFlood)
    {
        prepareFlood = false;
        FOR_ROW_COL_MV_MV
        {
            delta_z1->Drc = z->Data[r+1][c]-z->Drc;
            delta_z2->Drc = z->Data[r+1][c]-z->Drc;
            som_z1->Drc = z->Data[r][c-1]-2*z->Drc+z->Data[r][c+1];
            som_z2->Drc = z->Data[r-1][c]-2*z->Drc+z->Data[r+1][c];
        }
    }


    FOR_ROW_COL_MV_MV
    {
        //        hh2 = h->Drc-2.*h->Data[r][c+1]+h[i+2][l];
        //        uu2 = u->Drc-2.*u->Data[r][c+1]+u[i+2][l];
        //        vv2 = v->Drc-2.*v->Data[r][c+1]+v[i+2][l];
        hh2 = h->Data[r][c-1]-2.*h->Drc+h->Data[r][c+1];
        uu2 = u->Data[r][c-1]-2.*u->Drc+u->Data[r][c+1];
        vv2 = v->Data[r][c-1]-2.*v->Drc+v->Data[r][c+1];

        ddh2 = amortENO*limiter(hh1,hh2);
        ddu2 = amortENO*limiter(uu1,uu2);
        ddz2 = amortENO*limiter(hh1+som_z1->Drc,hh2+som_z1->Data[r][c+1]);
        ddv2 = amortENO*limiter(vv1,vv2);

        delta_h2 = h->Data[r][c+1]-h->Drc;
        delta_u2 = u->Data[r][c+1]-u->Drc;
        delta_v2 = v->Data[r][c+1]-v->Drc;

        if (startloop)
        {
            ddh1 = ddh2;
            ddz1 = ddz2;
            ddu1 = ddu2;
            ddv1 = ddv2;
            delta_h1 = delta_h2;
            delta_u1 = delta_u2;
            delta_v1 = delta_v2;
            startloop = false;
        }

        dh = limiter(delta_h1+ddh1*0.5,delta_h2-ddh2*0.5);
        dz_h = limiter(delta_h1+delta_z1->Data[r+1][c]+ddz1*0.5,delta_h2+delta_z1->Drc-ddz2*0.5);

        du = limiter(delta_u1+ddu1*0.5,delta_u2-ddu2*0.5);
        dv = limiter(delta_v1+ddv1*0.5,delta_v2-ddv2*0.5);

        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;

        h1r->Drc=h->Drc+dh*0.5;
        h1l->Drc=h->Drc-dh*0.5;

        z1r->Drc=z->Drc+0.5*(dz_h-dh);
        z1l->Drc=z->Drc+0.5*(dh-dz_h);
        delzc1->Drc = z1r->Drc-z1l->Drc;
        delz1->Data[r+1][c] = z1l->Drc-z1r->Data[r+1][c];

        if (h->Drc>0)
        {
            u1r->Drc=u->Drc+h1l->Drc*du*0.5/h->Drc;
            u1l->Drc=u->Drc-h1r->Drc*du*0.5/h->Drc;
            v1r->Drc=v->Drc+h1l->Drc*dv*0.5/h->Drc;
            v1l->Drc=v->Drc-h1r->Drc*dv*0.5/h->Drc;
        }
        else
        {
            u1r->Drc=u->Drc+du*0.5;
            u1l->Drc=u->Drc-du*0.5;
            v1r->Drc=v->Drc+dv*0.5;
            v1l->Drc=v->Drc-dv*0.5;
        } //end if

        ddh1=ddh2;
        ddz1=ddz2;
        ddu1=ddu2;
        ddv1=ddv2;

    }

    FOR_ROW_COL_MV_MV
    {
        //        hh2 = h->Drc-2.*h->Data[r+1][c]+h[i][l+2];
        //        uu2 = u->Drc-2.*u->Data[r+1][c]+u[i][l+2];
        //        vv2 = v->Drc-2.*v->Data[r+1][c]+v[i][l+2];
        hh2 = h->Data[r-1][c]-2.*h->Drc+h->Data[r+1][c];
        uu2 = u->Data[r-1][c]-2.*u->Drc+u->Data[r+1][c];
        vv2 = v->Data[r-1][c]-2.*v->Drc+v->Data[r+1][c];

        ddh2 = amortENO*limiter(hh1,hh2);
        ddu2 = amortENO*limiter(uu1,uu2);
        ddz2 = amortENO*limiter(hh1+som_z2->Drc,hh2+som_z2->Data[r+1][c]);
        ddv2 = amortENO*limiter(vv1,vv2);

        delta_h2 = h->Data[r+1][c]-h->Drc;
        delta_u2 = u->Data[r+1][c]-u->Drc;
        delta_v2 = v->Data[r+1][c]-v->Drc;

        if (startloop1)
        {
            ddh1 = ddh2;
            ddz1 = ddz2;
            ddu1 = ddu2;
            ddv1 = ddv2;
            delta_h1 = delta_h2;
            delta_u1 = delta_u2;
            delta_v1 = delta_v2;
            startloop1 = false;
        }


        dh = limiter(delta_h1+ddh1*0.5,delta_h2-ddh2*0.5);
        dz_h = limiter(delta_h1+delta_z2->Data[r-1][c]+ddz1*0.5,delta_h2+delta_z2->Drc-ddz2*0.5);

        du = limiter(delta_u1+ddu1*0.5,delta_u2-ddu2*0.5);
        dv = limiter(delta_v1+ddv1*0.5,delta_v2-ddv2*0.5);

        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;

        h2r->Drc = h->Drc+dh*0.5;
        h2l->Drc = h->Drc-dh*0.5;

        z2r->Drc = z->Drc+0.5*(dz_h-dh);
        z2l->Drc = z->Drc+0.5*(dh-dz_h);
        delzc2->Drc = z2r->Drc-z2l->Drc;
        delz2->Data[r-1][c] = z2l->Drc-z2r->Data[r-1][c];

        if (h->Drc>0)
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

        ddh1=ddh2;
        ddz1=ddz2;
        ddu1=ddu2;
        ddv1=ddv2;

    } //end for


}
//---------------------------------------------------------------------------
//reconstruction scheme
//original matrix is col, row orientation (y,x), NOT row, col (x,y)
//1 is in the y direction, 2 = x direction
//TO DO: only in flood domain
//called as
//rec->calcul(h,u,v,z,delzc1,delzc2,delz1,delz2,h1r,u1r,v1r,h1l,u1l,v1l,h2r,u2r,v2r,h2l,u2l,v2l);

void TWorld::MUSCL(TMMap *h,TMMap *u,TMMap *v,TMMap *z)
/*,
                   TMMap *delzc1,TMMap *delzc2,TMMap *delz1,TMMap *delz2,
                   TMMap *h1r,TMMap *u1r,TMMap *v1r,TMMap *h1l,TMMap *u1l,TMMap *v1l,
                   TMMap *h2r,TMMap *u2r,TMMap *v2r,TMMap *h2l,TMMap *u2l,TMMap *v2l)
                   */
{
    double delta_h1, delta_u1, delta_v1;
    double delta_h2, delta_u2, delta_v2;
    double dh, du, dv, dz_h;
    bool startloop = true;
    bool startloop1 = true;


    if (prepareFlood)
    {
        prepareFlood = false;
        FOR_ROW_COL_MV_MV
        {
            delta_z1->Drc = z->Data[r+1][c] - z->Drc;
            delta_z2->Drc = z->Data[r][c+1] - z->Drc;
        }
    }

    // first do 1: x direction (col+1 - col)
    FOR_ROW_COL_MV_MV
    {
        delta_h2 = h->Data[r+1][c] - h->Drc;
        delta_u2 = u->Data[r+1][c] - u->Drc;
        delta_v2 = v->Data[r+1][c] - v->Drc;

        if (startloop)
        {
            delta_h1 = delta_h2;
            delta_u1 = delta_u2;
            delta_v1 = delta_v2;
            startloop = false;
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
        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;
    }

    // then do 2: y direction (row+1 - row)
    FOR_ROW_COL_MV_MV
    {
        delta_h2 = h->Data[r][c+1] - h->Drc;
        delta_u2 = u->Data[r][c+1] - u->Drc;
        delta_v2 = v->Data[r][c+1] - v->Drc;

        if (startloop1)
        {
            delta_h1 = delta_h2;
            delta_u1 = delta_u2;
            delta_v1 = delta_v2;
            startloop1 = false;
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

        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;
    }
}
void TWorld::MUSCL2(TMMap *h,TMMap *u,TMMap *v, TMMap *z)
/*
                    TMMap *delzc1,TMMap *delzc2,TMMap *delz1,TMMap *delz2,
                    TMMap *h1r,TMMap *u1r,TMMap *v1r,TMMap *h1l,TMMap *u1l,TMMap *v1l,
                    TMMap *h2r,TMMap *u2r,TMMap *v2r,TMMap *h2l,TMMap *u2l,TMMap *v2l)*/
{
    double delta_h1, delta_u1, delta_v1;
    double delta_h2, delta_u2, delta_v2;
    double dh, du, dv, dz_h;
    bool startloop = true;
    bool startloop1 = true;


    if (prepareFlood)
    {
        prepareFlood = false;
        FOR_ROW_COL_MV_MV
        {
            delta_z1->Drc = z->Data[r][c+1] - z->Drc;
            delta_z2->Drc = z->Data[r+1][c] - z->Drc;
        }
    }

    FOR_ROW_COL_MV_MV
    {
        delta_h2 = h->Data[r][c+1]-h->Data[r][c];
        delta_u2 = u->Data[r][c+1]-u->Data[r][c];
        delta_v2 = v->Data[r][c+1]-v->Data[r][c];

        if (startloop)
        {
            delta_h1 = delta_h2;
            delta_u1 = delta_u2;
            delta_v1 = delta_v2;
            startloop = false;
        }


        dh = limiter(delta_h1,delta_h2);
        dz_h = limiter(delta_h1+delta_z1->Data[r][c-1],delta_h2+delta_z1->Data[r][c]);

        du = limiter(delta_u1,delta_u2);
        dv = limiter(delta_v1,delta_v2);

        dh = delta_h2;
        du = delta_u2;
        dv = delta_v2;
        dz_h = delta_h2+delta_z2->Data[r][c];

        h1r->Data[r][c]=h->Data[r][c]+dh*0.5;
        h1l->Data[r][c]=h->Data[r][c]-dh*0.5;

        z1r->Data[r][c]=z->Data[r][c]+0.5*(dz_h-dh);
        z1l->Data[r][c]=z->Data[r][c]+0.5*(dh-dz_h);
        delzc1->Data[r][c] = z1r->Data[r][c]-z1l->Data[r][c];
        delz1->Data[r][c-1] = z1l->Data[r][c]-z1r->Data[r][c-1];

        if (h->Data[r][c]>0.)
        {
            u1r->Data[r][c]=u->Data[r][c]+h1l->Data[r][c]*du*0.5/h->Data[r][c];
            u1l->Data[r][c]=u->Data[r][c]-h1r->Data[r][c]*du*0.5/h->Data[r][c];
            v1r->Data[r][c]=v->Data[r][c]+h1l->Data[r][c]*dv*0.5/h->Data[r][c];
            v1l->Data[r][c]=v->Data[r][c]-h1r->Data[r][c]*dv*0.5/h->Data[r][c];
        }
        else
        {
            u1r->Data[r][c]=u->Data[r][c]+du*0.5;
            u1l->Data[r][c]=u->Data[r][c]-du*0.5;
            v1r->Data[r][c]=v->Data[r][c]+dv*0.5;
            v1l->Data[r][c]=v->Data[r][c]-dv*0.5;
        } //end if
        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;

    }

    FOR_ROW_COL_MV_MV
    {

        delta_h2 = h->Data[r+1][c]-h->Data[r][c];
        delta_u2 = u->Data[r+1][c]-u->Data[r][c];
        delta_v2 = v->Data[r+1][c]-v->Data[r][c];

        if (startloop1)
        {
            delta_h1 = delta_h2;
            delta_u1 = delta_u2;
            delta_v1 = delta_v2;
            startloop1 = false;
        }


        dh = limiter(delta_h1,delta_h2);
        dz_h = limiter(delta_h1+delta_z2->Data[r-1][c],delta_h2+delta_z2->Data[r][c]);

        du = limiter(delta_u1,delta_u2);
        dv = limiter(delta_v1,delta_v2);

        dh = delta_h2;
        du = delta_u2;
        dv = delta_v2;
        dz_h = delta_h2+delta_z2->Data[r][c];

        h2r->Data[r][c] = h->Data[r][c]+dh*0.5;
        h2l->Data[r][c] = h->Data[r][c]-dh*0.5;

        z2r->Data[r][c] = z->Data[r][c]+0.5*(dz_h-dh);
        z2l->Data[r][c] = z->Data[r][c]+0.5*(dh-dz_h);
        delzc2->Data[r][c] = z2r->Data[r][c]-z2l->Data[r][c];
        delz2->Data[r-1][c] = z2l->Data[r][c]-z2r->Data[r-1][c];

        if (h->Data[r][c]>0.){
            u2r->Data[r][c] = u->Data[r][c]+h2l->Data[r][c]*du*0.5/h->Data[r][c];
            u2l->Data[r][c] = u->Data[r][c]-h2r->Data[r][c]*du*0.5/h->Data[r][c];
            v2r->Data[r][c] = v->Data[r][c]+h2l->Data[r][c]*dv*0.5/h->Data[r][c];
            v2l->Data[r][c] = v->Data[r][c]-h2r->Data[r][c]*dv*0.5/h->Data[r][c];
        }else{
            u2r->Data[r][c] = u->Data[r][c]+du*0.5;
            u2l->Data[r][c] = u->Data[r][c]-du*0.5;
            v2r->Data[r][c] = v->Data[r][c]+dv*0.5;
            v2l->Data[r][c] = v->Data[r][c]-dv*0.5;
        } //end if

        //rajout de l'iteration des delta
        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;

    } //end for l

}

//---------------------------------------------------------------------------
// St Venant equations: conservation of mass and momentum: h u et v
void TWorld::bloc2(double dt, TMMap *he, TMMap *ve1, TMMap *ve2, /*TMMap *qe1, TMMap *qe2, */
                   TMMap *hes, TMMap *ves1, TMMap *ves2)
//, TMMap *qes1, TMMap *qes2, TMMap *Vin, double tps, double, int n, double dtheta )
{
    tx=dt/dx;
    ty=dt/dy;

    FOR_ROW_COL_MV_MV
    {
        // Solution of the equation of mass conservation (First equation of Saint venant)
        // f1 comes from MUSCL calculations
        hes->Drc = he->Drc-tx*(f1->Data[r][c+1]-f1->Drc)-ty*(g1->Data[r+1][c]-g1->Drc);

     //   hes->Drc += FfSurplus->Drc*dt/_dt;
        // add infiltration here

        if (hes->Drc > he_ca)
        {
            //Solution of the equation of momentum (Second and third equation of Saint-venant)
            double qes1;
            double qes2;

            qes1 = he->Drc*ve1->Drc - tx*(f2->Data[r][c+1]-f2->Drc
                    + grav*0.5*(h1g->Drc*h1g->Drc -
                                h1l->Drc*h1l->Drc + h1r->Drc*h1r->Drc - h1d->Drc*h1d->Drc
                                + (h1l->Drc+h1r->Drc)*delzc1->Drc))
                    - ty*(g2->Data[r][c+1]-g2->Drc);

            qes2 = he->Drc*ve2->Drc -
                    tx*(f3->Data[r][c+1]-f3->Drc) -
                    ty*(g3->Data[r+1][c]-g3->Drc +
                    grav*0.5*(h2g->Drc*h2g->Drc - h2l->Drc*h2l->Drc +
                              h2r->Drc*h2r->Drc - h2d->Drc*h2d->Drc +
                              (h2l->Drc+h2r->Drc)*delzc2->Drc));
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
            hes->Drc = 0;
        }
    }
}
//---------------------------------------------------------------------------
/*Construction varibles for hydrostatic reconstruction.
 Flux with x and y.
 Calculaton of the time steps in relation to fixed cfl.*/
double TWorld::bloc1(double dt, double dt_max)
{
    // double cflfix = cfl_fix;
    double dt_tmp, dtx, dty;
    double velocity_max_x, velocity_max_y;
    dtx = dt_max;
    dty = dt_max;
    velocity_max_x = -ve_ca;
    velocity_max_y = -ve_ca;
    double cfl;

    FOR_ROW_COL_MV_MV
    {
        h1d->Data[r][c-1] = max(0, h1r->Data[r][c-1] - max(0,  delz1->Data[r][c-1]));
        h1g->Drc          = max(0, h1l->Drc          - max(0, -delz1->Data[r][c-1]));

        if (F_scheme == 1)
            F_Rusanov(h1d->Data[r][c-1], u1r->Data[r][c-1], v1r->Data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
        else
            if (F_scheme == 2)
                F_HLL(h1d->Data[r][c-1], u1r->Data[r][c-1], v1r->Data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
            else
                F_HLL2(h1d->Data[r][c-1], u1r->Data[r][c-1], v1r->Data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
        f1->Drc = HLL2_f1;
        f2->Drc = HLL2_f2;
        f3->Drc = HLL2_f3;
        cfl = HLL2_cfl;
        //      cfl = min(1000,HLL2_cfl);

        if (HLL2_cfl*dt/dx == 0)
            dt_tmp = dt_max;
        else
            dt_tmp = cfl_fix*dx/cfl;

        dtx = min(min(dt, dt_tmp), dtx);
        velocity_max_x = max(velocity_max_x, cfl);
    }

    qDebug() << "bloc1a" << dt << dtx << dt_tmp << HLL2_cfl << velocity_max_x;

    FOR_ROW_COL_MV_MV
    {
        h2d->Data[r-1][c] = max(0, h2r->Data[r-1][c] - max(0, delz2->Data[r-1][c]));
        h2g->Drc = max(0, h2l->Drc - max(0, -delz2->Data[r-1][c]));
        if (F_scheme == 1)
            F_Rusanov(h2d->Data[r-1][c],v2r->Data[r-1][c],u2r->Data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);
        else
            if (F_scheme == 2)
                F_HLL(h2d->Data[r-1][c],v2r->Data[r-1][c],u2r->Data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);
            else
                F_HLL2(h2d->Data[r-1][c],v2r->Data[r-1][c],u2r->Data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);
        g1->Drc = HLL2_f1;
        g2->Drc = HLL2_f3;
        g3->Drc = HLL2_f2;
        //    cfl = min(1000, HLL2_cfl);
        cfl = HLL2_cfl;

        if (HLL2_cfl*dt/dy == 0)
            dt_tmp = dt_max;
        else
            dt_tmp = cfl_fix*dy/cfl;
        dty = min(min(dt, dt_tmp), dty);
        velocity_max_y = max(velocity_max_y, cfl);
    }

    qDebug() << "bloc1b" << dt << dty << dt_tmp << HLL2_cfl << velocity_max_y;

    if (scheme_type == 1)
        return(max(dt_ca,min(dtx,dty)));
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
// first order solution
double TWorld::fullSWOF2Do1(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2)
{
    double timesum = 0;
    int n = 0;

    dx = _dx;
    dy = _dx; //TODO could be DX->Drc ? later !!!!

    dt_max = min(dx/2,dy/2);
    dt1 = dt_max;

    if (prepareFlood)
    {
        prepareFlood = false;
        FOR_ROW_COL_MV_MV
        {
            delz1->Data[r][c-1] = z->Drc - z->Data[r][c-1];
            delz2->Data[r-1][c] = z->Drc - z->Data[r-1][c];
        }
    }

    // if there is no flood skip everything
    if (startFlood)
    {

        do {
            n++;

            dt1 = min(dt1*qSqrt(double(n)), dt_max);

            setZero(h, u, v);//, q1, q2);

            FOR_ROW_COL_MV_MV
            {
                h1r->Drc = h->Drc;
                u1r->Drc = u->Drc;
                v1r->Drc = v->Drc;
                h1l->Data[r][c+1] = h->Data[r][c+1];
                u1l->Data[r][c+1] = u->Data[r][c+1];
                v1l->Data[r][c+1] = v->Data[r][c+1];
            }
            FOR_ROW_COL_MV_MV
            {
                h2r->Drc = h->Drc;
                u2r->Drc = u->Drc;
                v2r->Drc = v->Drc;
                h2l->Data[r+1][c] = h->Data[r+1][c];
                u2l->Data[r+1][c] = u->Data[r+1][c];
                v2l->Data[r+1][c] = v->Data[r+1][c];
            }

            double dta = bloc1(dt1, dt_max);
            dt1 = dta;
            dt1 = qMin(dt1, _dt-timesum);

            bloc2(dt1, h,u,v, hs,us,vs);

            setZero(hs, us, vs);//, q1, q2);

            FOR_ROW_COL_MV_MV
            {
                h->Drc = hs->Drc;
                u->Drc = us->Drc;
                v->Drc = vs->Drc;
                q1->Drc = h->Drc*u->Drc;
                q2->Drc = h->Drc*v->Drc;
            }

            timesum = timesum + dt1;
            // sum to reach _dt

        } while (timesum  < _dt);

    }
    //        Fr=froude_number(hs,us,vs);
    iter_n = n;
    return(timesum/(n+1));
}
//---------------------------------------------------------------------------
double TWorld::fullSWOF2D(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2)
{
    //initialization
    double dt2;
    double timesum = 0;

    if (prepareFlood)
        verif =1;

    dx = _dx;
    dy = _dx; //TODO could be DX->Drc ? later !!!!

    dt_max = min(dx/2,dy/2);
    dt1 = dt_max;

    // if there is no flood skip everything
    int n = 1;
    if (startFlood)
    {
        do {
            n++;
            dt1 = min(dt1*qSqrt(double(n)), dt_max);
            dt1 = min(dt1*(double(n)), dt_max);
//            dt1 = dt_max;

 //           if (verif == 1)
 //           {

                setZero(h, u, v);//, q1, q2);

                // Reconstruction for order 2
                MUSCL(h,u,v,z);
     //       }
            dt1 = bloc1(dt1, dt_max);

            dt1 = qMin(dt1, _dt-timesum);

            tx=dt1/dx;
            ty=dt1/dy;

            //h, u, v go in hs, vs, us come out
            bloc2(dt1, h,u,v, hs,us,vs);
            dt2 = dt1;

            setZero(hs, us, vs);//, qs1, qs2);

            //Reconstruction for order 2
            MUSCL(hs,us,vs,z);

            dt2 = bloc1(dt2, dt_max);

            if (dt2 < dt1)
            {
                dt1 = dt2;
                tx=dt1/dx;
                ty=dt1/dy;
                verif = 0;
            }
//            else
//            {
                verif = 1;
                //hs, us, vs go in hsa, vsa, usa come out
                bloc2(dt1, hs,us,vs, hsa, usa, vsa);

                setZero(hsa, usa, vsa);//, qs1, qs2);

                //Heun method (order 2 in time)
                FOR_ROW_COL_MV
                {
                    double tmp = 0.5*(h->Drc + hsa->Drc);
                    if (tmp >= he_ca)
                    {
                        q1->Drc = 0.5*(h->Drc*u->Drc + hsa->Drc*usa->Drc);
                        u->Drc = q1->Drc/tmp;
                        q2->Drc = 0.5*(h->Drc*v->Drc + hsa->Drc*vsa->Drc);
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
                }//Heun
   //         }//end for else dt2<dt1

            timesum = timesum + dt1;
            qDebug() << n;
        } while (timesum  < _dt);
    } // if floodstart

    iter_n = n;

    return(timesum/(n+1));

}
//---------------------------------------------------------------------------
double TWorld::fullSWOF2Do1a(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2)
{
    double timesum = 0;
    int n = 0;

    dx = _dx;
    dy = _dx; //TODO could be DX->Drc ? later !!!!

    dt_max = min(dx/2,dy/2);
    dt1 = dt_max;

    if (prepareFlood)
    {
        prepareFlood = false;
        FOR_ROW_COL_MV_MV
        {
            delz1->Data[r][c-1] = z->Drc - z->Data[r][c-1];
            delz2->Data[r-1][c] = z->Drc - z->Data[r-1][c];
        }
    }

    // if there is no flood skip everything
    if (startFlood)
    {

        do {
            n++;

            dt1 = min(dt1*qSqrt(double(n)), dt_max);

            setZero(h, u, v);//, q1, q2);

            FOR_ROW_COL_MV_MV
            {
                h1r->Drc = h->Drc;
                u1r->Drc = u->Drc;
                v1r->Drc = v->Drc;
                h1l->Data[r][c+1] = h->Data[r][c+1];
                u1l->Data[r][c+1] = u->Data[r][c+1];
                v1l->Data[r][c+1] = v->Data[r][c+1];
            }
            FOR_ROW_COL_MV_MV
            {
                h2r->Drc = h->Drc;
                u2r->Drc = u->Drc;
                v2r->Drc = v->Drc;
                h2l->Data[r+1][c] = h->Data[r+1][c];
                u2l->Data[r+1][c] = u->Data[r+1][c];
                v2l->Data[r+1][c] = v->Data[r+1][c];
            }

            dt1 = bloc1(dt1, dt_max);
            dt1 = qMin(dt1, _dt-timesum);

            bloc2(dt1, h,u,v, hs,us,vs);

            setZero(hs, us, vs);

            FOR_ROW_COL_MV_MV
            {
                h1r->Drc = hs->Drc;
                u1r->Drc = us->Drc;
                v1r->Drc = vs->Drc;
                h1l->Data[r][c+1] = hs->Data[r][c+1];
                u1l->Data[r][c+1] = us->Data[r][c+1];
                v1l->Data[r][c+1] = vs->Data[r][c+1];
            }
            FOR_ROW_COL_MV_MV
            {
                h2r->Drc = hs->Drc;
                u2r->Drc = us->Drc;
                v2r->Drc = vs->Drc;
                h2l->Data[r+1][c] = hs->Data[r+1][c];
                u2l->Data[r+1][c] = us->Data[r+1][c];
                v2l->Data[r+1][c] = vs->Data[r+1][c];
            }

            dt1 = bloc1(dt1, dt_max);
            dt1 = qMin(dt1, _dt-timesum);

            bloc2(dt1, hs,us,vs, hsa,usa,vsa);

            setZero(hsa, usa, vsa);//, q1, q2);

            FOR_ROW_COL_MV
            {
                double tmp = 0.5*(h->Drc+hsa->Drc);
                if (tmp >= he_ca)
                {
                    q1->Drc = 0.5*(h->Drc*u->Drc + hsa->Drc*usa->Drc);
                    u->Drc = q1->Drc/tmp;
                    q2->Drc = 0.5*(h->Drc*v->Drc + hsa->Drc*vsa->Drc);
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
            }//Heun

            timesum = timesum + dt1;

        } while (timesum  < _dt);
    }

    iter_n = n;
    return(timesum/(n+1));
}
