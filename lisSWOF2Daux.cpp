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

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "global.h"

#define he_ca 1e-12
#define ve_ca 1e-12
#define cfl_fix 0.4
#define grav 9.81
#define scheme_type 1
#define eps 1e-11

#define diff(a, b) ((a>0 && b>0) || (a<0 && b<0) ? (2*a*b/(a+b)) : 0)
#define diffA(a, b) ((a*b<0) ? 0 : ((a*(b*b+eps)+b*(a*a+eps))/(a*a+b*b+2*eps)))



//---------------------------------------------------------------------------
double TWorld::limiter(double a, double b)
{
    //return(diff(a,b));
    //   return(diffA(a,b));
    //    if (a*b<0)
    //        return(0);
    //    else
    //        return( (a*(b*b+1e-11)+b*(a*a+1e-11))/(a*a+b*b+2.*1e-11));


    if (a>=0. && b>=0)
        return(std::min(a,b));
    else
        if (a<=0. && b<=0)
            return(std::max(a,b));
        else
            return(0);
}
//---------------------------------------------------------------------------
void TWorld::setZero(CTMap *h, CTMap *u, CTMap *v, CTMap *q1, CTMap *q2)
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
            q1->Drc = 0;
        }
        //if (fabs(v->Drc) <= ve_ca)
        if (v->Drc <= ve_ca)
        {
            v->Drc = 0;
            q2->Drc = 0;
        }
    }
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
    if (hg<=0. && hd<=0.)
    {
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }
    else
    {
        double grav_hg = grav*hg;
        double grav_hd = grav*hd;
        double sqrt_grav_hg = sqrt(grav_hg);  //wave velocity
        double sqrt_grav_hd = sqrt(grav_hd);
        double qd = ud*hd;
        double qg = ug*hg;
        double c1 = std::min(ug-sqrt_grav_hg,ud-sqrt_grav_hd); //we already have ug-sqrt_grav_hg<ug+sqrt_grav_hg and ud-sqrt_grav_hd<ud+sqrt_grav_hd
        double c2 = std::max(ug+sqrt_grav_hg,ud+sqrt_grav_hd); //so we do not need all the eigenvalues to get c1 and c2
        double tmp = 1./(c2-c1);
        double t1 = (std::min(c2,0.)-std::min(c1,0.))*tmp;
        double t2 = 1.-t1;
        double t3 = (c2*fabs(c1)-c1*fabs(c2))*0.5*tmp;

        f1 = t1*qd+t2*qg-t3*(hd-hg);
        f2 = t1*(qd*ud+grav_hd*hd*0.5)+t2*(qg*ug+grav_hg*hg*0.5)-t3*(qd-qg);
        f3 = t1*qd*vd+t2*qg*vg-t3*(hd*vd-hg*vg);
        cfl = std::max(fabs(c1),fabs(c2)); //cfl is the velocity to compute the cfl condition std::max(fabs(c1),fabs(c2))*tx with tx=dt/dx
    }
    HLL2_cfl = cfl;
    HLL2_f1 = f1;
    HLL2_f2 = f2;
    HLL2_f3 = f3;
}

void TWorld::F_HLL(double hg,double ug,double vg,double hd,double ud,double vd)
{
    double f1, f2, f3, cfl;
    if (hg<=0. && hd<=0.)
    {
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }
    else
    {
        double grav_hg = grav*hg;
        double grav_hd = grav*hd;
        double qd = ud*hd;
        double qg = ug*hg;
        double c1 = std::min(ug-sqrt(grav_hg),ud-sqrt(grav_hd));
        double c2 = std::max(ug+sqrt(grav_hg),ud+sqrt(grav_hd));

        //cfl is the velocity to calculate the real cfl=std::max(fabs(c1),fabs(c2))*tx with tx=dt/dx
        if (c1==0. && c2==0.)
        {
            //dry
            f1=0.;
            f2=0.;
            f3=0.;
            cfl=0.; //std::max(fabs(c1),fabs(c2))=0
        }
        else
            if (c1>=0.){
                //supercritical flow, from left to right : we have std::max(abs(c1),abs(c2))=c2>0
                f1=qg;
                f2=qg*ug+grav*hg*hg*0.5;
                f3=qg*vg;
                cfl=c2; //std::max(fabs(c1),fabs(c2))=c2>0
            }
            else
                if (c2<=0.)
                {
                    //supercritical flow, from right to left : we have std::max(abs(c1),abs(c2))=-c1>0
                    f1=qd;
                    f2=qd*ud+grav*hd*hd*0.5;
                    f3=qd*vd;
                    cfl=fabs(c1); //std::max(fabs(c1),fabs(c2))=fabs(c1)
                }
                else
                { //subcritical flow
                    double tmp = 1./(c2-c1);
                    f1=(c2*qg-c1*qd)*tmp+c1*c2*(hd-hg)*tmp;
                    f2=(c2*(qg*ug+grav*hg*hg*0.5)-c1*(qd*ud+grav*hd*hd*0.5))*tmp+c1*c2*(qd-qg)*tmp;
                    f3=(c2*(qg*vg)-c1*(qd*vd))*tmp+c1*c2*(hd*vd-hg*vg)*tmp;
                    cfl = std::max(fabs(c1),fabs(c2));
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
    if (hg<=0. && hd<=0.)
    {
        c = 0.;
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }
    else
    {
        c = std::max(fabs(ug)+sqrt(grav*hg),fabs(ud)+sqrt(grav*hd));
        double cd = c*0.5;
        double qd = ud*hd;
        double qg = ug*hg;
        f1 = (qg+qd)*0.5-cd*(hd-hg);
        f2 = ((ug*qg)+(grav*0.5*hg*hg)+(ud*qd)+(grav*0.5*hd*hd))*0.5-cd*(qd-qg);
        f3 = (qg*vg+qd*vd)*0.5-cd*(hd*vd-hg*vg);
        cfl = c;//*tx;
    }
    HLL2_cfl = cfl;
    HLL2_f1 = f1;
    HLL2_f2 = f2;
    HLL2_f3 = f3;
}
//---------------------------------------------------------------------------

void TWorld::MUSCL2(CTMap *h,CTMap *u,CTMap *v,CTMap *z,
                    CTMap *delzc1,CTMap *delzc2,CTMap *delz1,CTMap *delz2,
                    CTMap *h1r,CTMap *u1r,CTMap *v1r,CTMap *h1l,CTMap *u1l,CTMap *v1l,
                    CTMap *h2r,CTMap *u2r,CTMap *v2r,CTMap *h2l,CTMap *u2l,CTMap *v2l)
{
    double delta_h1, delta_u1, delta_v1;
    double delta_h2, delta_u2, delta_v2;
    double dh, du, dv, dz_h;

    delta_h1 = 0;
    delta_u1 = 0;
    delta_v1 = 0;


    FOR_ROW_COL_MV_MV
    //{
        delta_h2 = h->Data[r][c+1] - h->Drc;
        delta_u2 = u->Data[r][c+1] - u->Drc;
        delta_v2 = v->Data[r][c+1] - v->Drc;

        dh = limiter(delta_h1,delta_h2);
        dz_h =limiter(delta_h1+delta_z1->Data[r][c-1],delta_h2+delta_z1->Drc);

        du = limiter(delta_u1,delta_u2);
        dv = limiter(delta_v1,delta_v2);

        h1r->Drc=h->Drc+dh*0.5;
        h1l->Drc=h->Drc-dh*0.5;

        z1r->Drc=z->Drc+0.5*(dz_h-dh);
        z1l->Drc=z->Drc+0.5*(dh-dz_h);
        delzc1->Drc = z1r->Drc-z1l->Drc;
        delz1->Data[r][c-1] = z1l->Drc-z1r->Data[r][c-1];

        if (h->Drc>0.)
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
        //rajout de l'iteration des delta
        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;

    }
    delta_h1 = 0;
    delta_u1 = 0;
    delta_v1 = 0;

    FOR_ROW_COL_MV_MV
    //{
        delta_h2 = h->Data[r+1][c]-h->Drc;
        delta_u2 = u->Data[r+1][c]-u->Drc;
        delta_v2 = v->Data[r+1][c]-v->Drc;

        dh = limiter(delta_h1,delta_h2);
        dz_h = limiter(delta_h1+delta_z2->Data[r-1][c], delta_h2+delta_z2->Drc);

        du = limiter(delta_u1,delta_u2);
        dv = limiter(delta_v1,delta_v2);

        h2r->Drc = h->Drc+dh*0.5;
        h2l->Drc = h->Drc-dh*0.5;

        z2r->Drc = z->Drc+0.5*(dz_h-dh);
        z2l->Drc = z->Drc+0.5*(dh-dz_h);
        delzc2->Drc = z2r->Drc-z2l->Drc;
        delz2->Data[r-1][c] = z2l->Drc-z2r->Data[r-1][c];

        if (h->Drc>0.)
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


void TWorld::ENO(CTMap *h,CTMap *u,CTMap *v,CTMap *z,CTMap *delzc1,CTMap *delzc2,CTMap *delz1,CTMap *delz2,CTMap *h1r,CTMap *u1r,CTMap *v1r,CTMap *h1l,CTMap *u1l,CTMap *v1l,CTMap *h2r,CTMap *u2r,CTMap *v2r,CTMap *h2l,CTMap *u2l,CTMap *v2l)
{
    double ddh1=0.;
    double ddz1=0.;
    double ddu1=0.;
    double ddv1=0.;
    double ddh2=0.;
    double ddz2=0.;
    double ddu2=0.;
    double ddv2=0.;
    double uu1 = 0, vv1 = 0, hh1 = 0;
    double uu2, vv2, hh2;
    double delta_h1, delta_v1, delta_u1;
    double delta_h2, delta_v2, delta_u2;
    double dh, du, dv, dz_h;
    bool hoi = true;
    bool hoi1 = true;
    double amortENO = 0.25;
    delta_h1 = 0;
    delta_u1 = 0;
    delta_v1 = 0;

    FOR_ROW_COL_MV_MV
    //{
        delta_z1->Drc = z->Data[r+1][c]-z->Drc;
        delta_z2->Drc = z->Data[r+1][c]-z->Drc;
        som_z1->Drc = z->Data[r][c-1]-2*z->Drc+z->Data[r][c+1];
        som_z2->Drc = z->Data[r-1][c]-2*z->Drc+z->Data[c+1][r];
    }


    FOR_ROW_COL_MV_MV
    //{
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

        if (hoi)
        {
            ddh1=ddh2;
            ddz1=ddz2;
            ddu1=ddu2;
            ddv1=ddv2;
            delta_h1 = delta_h2;
            delta_u1 = delta_u2;
            delta_v1 = delta_v2;
            hoi = false;
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

        if (h->Drc>0.)
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


    ddh1=0.;
    ddz1=0.;
    ddu1=0.;
    ddv1=0.;

    FOR_ROW_COL_MV_MV
    //{
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

        if (hoi1)
        {
            ddh1=ddh2;
            ddz1=ddz2;
            ddu1=ddu2;
            ddv1=ddv2;
            delta_h1 = delta_h2;
            delta_u1 = delta_u2;
            delta_v1 = delta_v2;
            hoi1 = false;
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

        if (h->Drc>0.){
            u2r->Drc = u->Drc+h2l->Drc*du*0.5/h->Drc;
            u2l->Drc = u->Drc-h2r->Drc*du*0.5/h->Drc;
            v2r->Drc = v->Drc+h2l->Drc*dv*0.5/h->Drc;
            v2l->Drc = v->Drc-h2r->Drc*dv*0.5/h->Drc;
        }else{
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

void TWorld::MUSCL(CTMap *h,CTMap *u,CTMap *v,CTMap *z,
                   CTMap *delzc1,CTMap *delzc2,CTMap *delz1,CTMap *delz2,
                   CTMap *h1r,CTMap *u1r,CTMap *v1r,CTMap *h1l,CTMap *u1l,CTMap *v1l,
                   CTMap *h2r,CTMap *u2r,CTMap *v2r,CTMap *h2l,CTMap *u2l,CTMap *v2l)
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

    if (prepareFlood)
    {
        prepareFlood = false;
        FOR_ROW_COL_MV_MV
        //{
            delta_z1->Drc = z->Data[r+1][c] - z->Drc;
            delta_z2->Drc = z->Data[r][c+1] - z->Drc;
        }
    }
    bool hoi = true;
    // first do 1: x direction (col+1 - col)
    FOR_ROW_COL_MV_MV
    //{

        delta_h2 = h->Data[r+1][c] - h->Drc;
        delta_u2 = u->Data[r+1][c] - u->Drc;
        delta_v2 = v->Data[r+1][c] - v->Drc;

        if (hoi)
        {
            hoi = false;
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
        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;
    }

    // then do 2: y direction (row+1 - row)
    FOR_ROW_COL_MV_MV
    //{
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

        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;
    }
}

//---------------------------------------------------------------------------
