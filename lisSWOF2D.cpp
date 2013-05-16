/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010 - 2013  Victor Jetten
**  contact: v.g.jetten@utwente.nl
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
**  Developed in: VC2010/Qt
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
 \file lisSWOF2d.cpp
 \brief Channel flood using the opensource code from the FullSWOF_2D project.
        St Vennant equations
        authors: Olivier Delestre, Christian Laguerre, Carine Lucas, Ulrich Razafison, Marie Rousseau.
        website: http://www.univ-orleans.fr/mapmo/soft/FullSWOF/

functions: \n
-   double fullSWOF2D(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2);
-   double fullSWOF2Do1(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2);
-   void MUSCL(TMMap *h,TMMap *u,TMMap *v,TMMap *z);
-   void ENO(TMMap *h,TMMap *u,TMMap *v,TMMap *z);
-   double bloc1(double dt, double dt_max);
-   void bloc2(double dt, TMMap *he, TMMap *ve1, TMMap *ve2,TMMap *hes, TMMap *ves1, TMMap *ves2);
-   void F_HLL2(double hg,double ug,double vg,double hd,double ud,double vd);
-   void F_HLL(double hg,double ug,double vg,double hd,double ud,double vd);
-   void F_Rusanov(double hg,double ug,double vg,double hd,double ud,double vd);
-   void Fr_Manning(double uold, double vold, double hnew, double q1new, double q2new, double dt, double cf);
-   void Fr_ManningSf(double h, double u, double v, double cf);
-   void setZero(TMMap *h, TMMap *u, TMMap *v);//, TMMap *q1, TMMap *q2);
-   double limiter(double a, double b);
*/

#include "lisemqt.h"
#include "model.h"
#include "global.h"

#define he_ca 1e-12
#define ve_ca 1e-12

#define dt_ca 1e-4

#define GRAV 9.8067   // to copy swof code directly
#define EPSILON 1e-6
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
void TWorld::setZero(TMMap *_h, TMMap *_u, TMMap *_v)
{
    FOR_ROW_COL_MV
    {
        if (_h->Drc <= he_ca)
        {
            _h->Drc = 0;
            _u->Drc = 0;
            _v->Drc = 0;
        }
        //if (fabs(_u->Drc) <= ve_ca)
        if (_u->Drc <= ve_ca)
        {
            _u->Drc = 0;
            //q1->Drc = 0;
        }
        //if (fabs(_v->Drc) <= ve_ca)
        if (_v->Drc <= ve_ca)
        {
            _v->Drc = 0;
            //q2->Drc = 0;
        }
    }
}
//---------------------------------------------------------------------------
//friction slope
void TWorld::Fr_Manning(double uold, double vold, double hnew, double q1new, double q2new, double dt, double N)
{
    q1mod = q1new/(1.0+N*N*GRAV*sqrt(uold*uold+vold*vold)*dt/pow(hnew,4.0/3.0));
    q2mod = q2new/(1.0+N*N*GRAV*sqrt(uold*uold+vold*vold)*dt/pow(hnew,4.0/3.0));
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
void TWorld::F_HLL2(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    double f1, f2, f3, cfl;
    if (h_L<=0. && h_R<=0.)
    {
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }
    else
    {
        double grav_h_L = GRAV*h_L;
        double grav_h_R = GRAV*h_R;
        double sqrt_grav_h_L = sqrt(grav_h_L);
        double sqrt_grav_h_R = sqrt(grav_h_R);
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;

        double c1 = min(u_L - sqrt_grav_h_L,u_R - sqrt_grav_h_R); //we already have u_L - sqrt_grav_h_L<u_L + sqrt_grav_h_L and u_R - sqrt_grav_h_R<u_R + sqrt_grav_h_R
        double c2 = max(u_L + sqrt_grav_h_L,u_R + sqrt_grav_h_R); //so we do not need all the eigenvalues to get c1 and c2
        double tmp = 1./(c2-c1);
        double t1 = (min(c2,0.) - min(c1,0.))*tmp;
        double t2 = 1. - t1;
        double t3 = (c2*fabs(c1) - c1*fabs(c2))*0.5*tmp;

        f1 = t1*q_R + t2*q_L - t3*(h_R - h_L);
        f2 = t1*(q_R*u_R + grav_h_R*h_R*0.5) + t2*(q_L*u_L + grav_h_L*h_L*0.5) - t3*(q_R - q_L);
        f3 = t1*q_R*v_R + t2*q_L*v_L - t3*(h_R*v_R - h_L*v_L);
        cfl = max(fabs(c1),fabs(c2)); //cfl is the velocity to compute the cfl condition max(fabs(c1),fabs(c2))*tx with tx=dt/dx
    }
    HLL2_cfl = cfl;
    HLL2_f1 = f1;
    HLL2_f2 = f2;
    HLL2_f3 = f3;
}

void TWorld::F_HLL(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    double f1, f2, f3, cfl;
    if (h_L<=0. && h_R<=0.){
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }else{
        double grav_h_L = GRAV*h_L;
        double grav_h_R = GRAV*h_R;
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;
        double c1 = min(u_L-sqrt(grav_h_L),u_R-sqrt(grav_h_R));
        double c2 = max(u_L+sqrt(grav_h_L),u_R+sqrt(grav_h_R));

        //cfl is the velocity to calculate the real cfl=max(fabs(c1),fabs(c2))*tx with tx=dt/dx
        if (fabs(c1)<EPSILON && fabs(c2)<EPSILON){              //dry state
            f1=0.;
            f2=0.;
            f3=0.;
            cfl=0.; //max(fabs(c1),fabs(c2))=0
        }else if (c1>=EPSILON){ //supercritical flow, from left to right : we have max(abs(c1),abs(c2))=c2>0
            f1=q_L;
            f2=q_L*u_L+GRAV*h_L*h_L*0.5;
            f3=q_L*v_L;
            cfl=c2; //max(fabs(c1),fabs(c2))=c2>0
        }else if (c2<=-EPSILON){ //supercritical flow, from right to left : we have max(abs(c1),abs(c2))=-c1>0
            f1=q_R;
            f2=q_R*u_R+GRAV*h_R*h_R*0.5;
            f3=q_R*v_R;
            cfl=fabs(c1); //max(fabs(c1),fabs(c2))=fabs(c1)
        }else{ //subcritical flow
            double tmp = 1./(c2-c1);
            f1=(c2*q_L-c1*q_R)*tmp+c1*c2*(h_R-h_L)*tmp;
            f2=(c2*(q_L*u_L+GRAV*h_L*h_L*0.5)-c1*(q_R*u_R+GRAV*h_R*h_R*0.5))*tmp+c1*c2*(q_R-q_L)*tmp;
            f3=(c2*(q_L*v_L)-c1*(q_R*v_R))*tmp+c1*c2*(h_R*v_R-h_L*v_L)*tmp;
            cfl=max(fabs(c1),fabs(c2));
        }
    }
    HLL2_cfl = cfl;
    HLL2_f1 = f1;
    HLL2_f2 = f2;
    HLL2_f3 = f3;
}

void TWorld::F_Rusanov(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    double f1, f2, f3, cfl;
    double c;
    if (h_L<=0. && h_R<=0.){
        c = 0.;
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }else{
        c = max(fabs(u_L)+sqrt(GRAV*h_L),fabs(u_R)+sqrt(GRAV*h_R));
        double cd = c*0.5;
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;
        f1 = (q_L+q_R)*0.5-cd*(h_R-h_L);
        f2 = ((u_L*q_L)+(GRAV*0.5*h_L*h_L)+(u_R*q_R)+(GRAV*0.5*h_R*h_R))*0.5-cd*(q_R-q_L);
        f3 = (q_L*v_L+q_R*v_R)*0.5-cd*(h_R*v_R-h_L*v_L);
        cfl = c;//*tx;
    }
    HLL2_cfl = cfl;
    HLL2_f1 = f1;
    HLL2_f2 = f2;
    HLL2_f3 = f3;
}
//---------------------------------------------------------------------------

void TWorld::ENO(TMMap *h,TMMap *u,TMMap *v,TMMap *z)
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
    double amortENO = 0.25;
    delta_h1 = 0;
    delta_u1 = 0;
    delta_v1 = 0;

    FOR_ROW_COL_MV_MV
            if(!IS_MV_REAL8(&LDD->Data[r][c+2]))
    {
        hh2 = h->Drc-2.*h->Data[r][c+1]+h->Data[r][c+2];
        uu2 = u->Drc-2.*u->Data[r][c+1]+u->Data[r][c+2];
        vv2 = v->Drc-2.*v->Data[r][c+1]+v->Data[r][c+2];

        ddh2 = amortENO*limiter(hh1,hh2);
        ddu2 = amortENO*limiter(uu1,uu2);
        ddz2 = amortENO*limiter(hh1+som_z1->Drc,hh2+som_z1->Data[r][c+1]);
        ddv2 = amortENO*limiter(vv1,vv2);

        delta_h2 = h->Data[r][c+1]-h->Drc;
        delta_u2 = u->Data[r][c+1]-u->Drc;
        delta_v2 = v->Data[r][c+1]-v->Drc;

        dh = limiter(delta_h1+ddh1*0.5,delta_h2-ddh2*0.5);
        dz_h = limiter(delta_h1+delta_z1->Data[r][c-1]+ddz1*0.5,
                delta_h2+delta_z1->Drc-ddz2*0.5);

        du = limiter(delta_u1+ddu1*0.5,delta_u2-ddu2*0.5);
        dv = limiter(delta_v1+ddv1*0.5,delta_v2-ddv2*0.5);

        h1r->Drc=h->Drc+dh*0.5;
        h1l->Drc=h->Drc-dh*0.5;

        z1r->Drc=z->Drc+0.5*(dz_h-dh);
        z1l->Drc=z->Drc+0.5*(dh-dz_h);

        delzc1->Drc = z1r->Drc - z1l->Drc;
        delz1->Data[r][c-1] = z1l->Drc - z1r->Data[r][c-1];

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

        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;

        ddh1=ddh2;
        ddz1=ddz2;
        ddu1=ddu2;
        ddv1=ddv2;

    }

    FOR_ROW_COL_MV_MV
            if(!IS_MV_REAL8(&LDD->Data[r+2][c]))
    {
        hh2 = h->Drc-2.*h->Data[r+1][c]+h->Data[r+2][c];
        uu2 = u->Drc-2.*u->Data[r+1][c]+u->Data[r+2][c];
        vv2 = v->Drc-2.*v->Data[r+1][c]+v->Data[r+2][c];

        ddh2 = amortENO*limiter(hh1,hh2);
        ddu2 = amortENO*limiter(uu1,uu2);
        ddz2 = amortENO*limiter(hh1+som_z2->Drc,hh2+som_z2->Data[r+1][c]);
        ddv2 = amortENO*limiter(vv1,vv2);

        delta_h2 = h->Data[r+1][c]-h->Drc;
        delta_u2 = u->Data[r+1][c]-u->Drc;
        delta_v2 = v->Data[r+1][c]-v->Drc;

        dh = limiter(delta_h1+ddh1*0.5,delta_h2-ddh2*0.5);
        dz_h = limiter(delta_h1+delta_z2->Data[r-1][c]+ddz1*0.5,delta_h2+delta_z2->Drc-ddz2*0.5);

        du = limiter(delta_u1+ddu1*0.5,delta_u2-ddu2*0.5);
        dv = limiter(delta_v1+ddv1*0.5,delta_v2-ddv2*0.5);

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

        delta_h1 = delta_h2;
        delta_u1 = delta_u2;
        delta_v1 = delta_v2;
        ddh1=ddh2;
        ddz1=ddz2;
        ddu1=ddu2;
        ddv1=ddv2;

    } //end for


}
//---------------------------------------------------------------------------
//reconstruction scheme
//original matrix is col, row orientation (y,x), NOT row, col (x,y)
//TO DO: only in flood domain
//rec->calcul(h,u,v,z,delzc1,delzc2,delz1,delz2,h1r,u1r,v1r,h1l,u1l,v1l,h2r,u2r,v2r,h2l,u2l,v2l);
// all maps are inherited so no need to put them in the header

void TWorld::MUSCL(TMMap *ah, TMMap *au, TMMap *av, TMMap *az)
{
    double delta_h1, delta_u1, delta_v1;
    double delta_h2, delta_u2, delta_v2;
    double dh, du, dv, dz_h;

    FOR_ROW_COL_MV
    {
        tm->Drc = 0;
        tma->Drc = 0;
        tmb->Drc = 0;
    }

    FOR_ROW_COL_MV_MV
    {
        delta_h1 = tm->Drc;
        delta_u1 = tma->Drc;
        delta_v1 = tmb->Drc;

        delta_h2 = ah->Data[r][c+1] - ah->Drc;
        delta_u2 = au->Data[r][c+1] - au->Drc;
        delta_v2 = av->Data[r][c+1] - av->Drc;

        dh = limiter(delta_h1, delta_h2);
        dz_h = limiter(delta_h1 + delta_z1->Data[r][c-1],
                delta_h2 + delta_z1->Drc);
        du = limiter(delta_u1, delta_u2);
        dv = limiter(delta_v1, delta_v2);

        h1r->Drc = ah->Drc+dh*0.5;
        h1l->Drc = ah->Drc-dh*0.5;

        z1r->Drc = az->Drc+0.5*(dz_h-dh);
        z1l->Drc = az->Drc+0.5*(dh-dz_h);

        delzc1->Drc = z1r->Drc-z1l->Drc;
        delz1->Data[r][c-1] = z1l->Drc-z1r->Data[r][c-1];

        if (ah->Drc > 0.)
        {
            u1r->Drc = au->Drc + h1l->Drc*du*0.5/ah->Drc;
            u1l->Drc = au->Drc - h1r->Drc*du*0.5/ah->Drc;
            v1r->Drc = av->Drc + h1l->Drc*dv*0.5/ah->Drc;
            v1l->Drc = av->Drc - h1r->Drc*dv*0.5/ah->Drc;
        }
        else
        {
            u1r->Drc = au->Drc + du*0.5;
            u1l->Drc = au->Drc - du*0.5;
            v1r->Drc = av->Drc + dv*0.5;
            v1l->Drc = av->Drc - dv*0.5;
        }
        tm->Drc = delta_h2;
        tma->Drc = delta_u2;
        tmb->Drc = delta_v2;
    }

    FOR_ROW_COL_MV
    {
        tm->Drc = 0;
        tma->Drc = 0;
        tmb->Drc = 0;
    }

    //        for (int r = _nrRows-2; r > 0; r--)
    for (int c = 0; c < _nrCols; c++)
        for (int r = 1; r < _nrRows-1; r++)
            if(!IS_MV_REAL8(&LDD->Data[r][c]) &&
                    !IS_MV_REAL8(&LDD->Data[r-1][c]) &&
                    !IS_MV_REAL8(&LDD->Data[r+1][c]))
                //FOR_ROW_COL_MV_MV
            {
                delta_h1 = tm->Drc;
                delta_u1 = tma->Drc;
                delta_v1 = tmb->Drc;

                delta_h2 = ah->Data[r+1][c] - ah->Drc;
                delta_u2 = au->Data[r+1][c] - au->Drc;
                delta_v2 = av->Data[r+1][c] - av->Drc;

                dh = limiter(delta_h1, delta_h2);
                dz_h = limiter(delta_h1+delta_z2->Data[r-1][c],
                        delta_h2+delta_z2->Drc);

                du = limiter(delta_u1, delta_u2);
                dv = limiter(delta_v1, delta_v2);

                h2r->Drc = ah->Drc+dh*0.5;
                h2l->Drc = ah->Drc-dh*0.5;

                z2r->Drc = az->Drc+0.5*(dz_h-dh);
                z2l->Drc = az->Drc+0.5*(dh-dz_h);
                delzc2->Drc = z2r->Drc - z2l->Drc;
                delz2->Data[r-1][c] = z2l->Drc - z2r->Data[r-1][c];

                if (ah->Drc > he_ca)
                {
                    u2r->Drc = au->Drc + h2l->Drc*du*0.5/ah->Drc;
                    u2l->Drc = au->Drc - h2r->Drc*du*0.5/ah->Drc;
                    v2r->Drc = av->Drc + h2l->Drc*dv*0.5/ah->Drc;
                    v2l->Drc = av->Drc - h2r->Drc*dv*0.5/ah->Drc;
                }
                else
                {
                    u2r->Drc = au->Drc + du*0.5;
                    u2l->Drc = au->Drc - du*0.5;
                    v2r->Drc = av->Drc + dv*0.5;
                    v2l->Drc = av->Drc - dv*0.5;
                }

                tm->Drc = delta_h2;
                tma->Drc = delta_u2;
                tmb->Drc = delta_v2;
            }
}

//---------------------------------------------------------------------------
// St Venant equations: conservation of mass and momentum: h u et v
void TWorld::maincalcscheme(double dt, TMMap *he, TMMap *ve1, TMMap *ve2,
                            TMMap *hes, TMMap *ves1, TMMap *ves2)
{
    FOR_ROW_COL_MV_MV
    {
        double dx = _dx;//-ChannelWidthUpDX->Drc;
        double dy = _dx;//DX->Drc; //_dx;
        double tx = dt/dx;
        double ty = dt/dy;

        // Solution of the equation of mass conservation (First equation of Saint venant)
        // f1 comes from MUSCL calculations
        hes->Drc = he->Drc - tx*(f1->Data[r][c+1]-f1->Drc) - ty*(g1->Data[r+1][c]-g1->Drc);

        //        double hhh = hes->Drc;
        //        double inf = ChannelDepth->Drc == 0 ? Ksat1->Drc*dt/3600000.0 : 0;
        //        hhh = hhh - inf;
        //        if (hhh < 0)
        //            hhh = 0;
        //        hes->Drc = hhh;

        if (hes->Drc > he_ca)
        {
            //Solution of the equation of momentum (Second and third equation of Saint-venant)
            double qes1;
            double qes2;

            // fullswof version 1.04
            qes1 = he->Drc*ve1->Drc-tx*(f2->Data[r][c+1]-f2->Drc +
                    GRAV*0.5*((h1g->Drc-h1l->Drc)*(h1g->Drc+h1l->Drc) +
                              (h1r->Drc-h1d->Drc)*(h1r->Drc+h1d->Drc) +
                              (h1l->Drc+h1r->Drc)*delzc1->Drc)) -
                    ty*(g2->Data[r][c+1]-g2->Drc);

            // fullswof version 1.04
            qes2 = he->Drc*ve2->Drc - tx*(f3->Data[r][c+1]-f3->Drc) -
                    ty*(g3->Data[r+1][c]-g3->Drc +
                    GRAV*0.5*((h2g->Drc-h2l->Drc)*(h2g->Drc+h2l->Drc) +
                              (h2r->Drc-h2d->Drc)*(h2r->Drc+h2d->Drc) +
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
double TWorld::maincalcflux(double dt, double dt_max)
{
    // double cflfix = cfl_fix;
    double dt_tmp, dtx, dty;
    double velocity_max_x, velocity_max_y;
    dtx = dt_max;
    dty = dt_max;
    velocity_max_x = -ve_ca;
    velocity_max_y = -ve_ca;
    double cfl = 0;
    double ff1 = 0, ff2 = 0, ff3 =0, cflo = 0;
    double dx, dy;
    FOR_ROW_COL_MV_MV
    {
        dx = _dx;
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

        // correct sudden extreme alues, replace with previous
        if (cfl > 100*(cflo+1))
        {
            f1->Drc = ff1;
            f2->Drc = ff2;
            f3->Drc = ff3;
            cfl = cflo;
        }

        ff1 = HLL2_f1;
        ff2 = HLL2_f2;
        ff3 = HLL2_f3;
        cflo = cfl;

        if (qFabs(HLL2_cfl*dt/dx) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = cfl_fix*dx/cfl;
        dtx = min(min(dt, dt_tmp), dtx);
        velocity_max_x = max(velocity_max_x, cfl);

    }
    qDebug() << "mainflux x" << dt << dtx << dt_tmp << velocity_max_x;



    FOR_ROW_COL_MV_MV
    {
        dy = _dx;// DX->Drc;

        h2d->Data[r-1][c] = max(0, h2r->Data[r-1][c] - max(0,  delz2->Data[r-1][c]));
        h2g->Drc          = max(0, h2l->Drc          - max(0, -delz2->Data[r-1][c]));

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
        cfl = HLL2_cfl;

        // correct sudden extreme alues, replace with previous
        if (cfl > 100*(cflo+1))
        {
            g1->Drc = ff1;
            g2->Drc = ff2;
            g3->Drc = ff3;
            cfl = cflo;
        }

        ff1 = HLL2_f1;
        ff2 = HLL2_f2;
        ff3 = HLL2_f3;
        cflo = cfl;

        if (qFabs(HLL2_cfl*dt/dy) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = cfl_fix*dy/cfl;
        dty = min(min(dt, dt_tmp), dty);
        velocity_max_y = max(velocity_max_y, cfl);
    }

    qDebug() << "mainflux y" << dt << dty << dt_tmp << velocity_max_y;



    if (scheme_type == 1)
    {
        return(max(dt_ca, min(dtx,dty)));
    }
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
void TWorld::simpleScheme(TMMap *_h,TMMap *_u,TMMap *_v)
{
    FOR_ROW_COL_MV_MV
    {
        h1r->Drc = _h->Drc;
        u1r->Drc = _u->Drc;
        v1r->Drc = _v->Drc;
        h1l->Data[r][c+1] = _h->Data[r][c+1];
        u1l->Data[r][c+1] = _u->Data[r][c+1];
        v1l->Data[r][c+1] = _v->Data[r][c+1];
    }
    FOR_ROW_COL_MV_MV
    {
        h2r->Drc = _h->Drc;
        u2r->Drc = _u->Drc;
        v2r->Drc = _v->Drc;
        h2l->Data[r+1][c] = _h->Data[r+1][c];
        u2l->Data[r+1][c] = _u->Data[r+1][c];
        v2l->Data[r+1][c] = _v->Data[r+1][c];
    }
}
//---------------------------------------------------------------------------
// first order solution
double TWorld::fullSWOF2Do1(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2)
{
    double timesum = 0;
    int n = 0;
    double dt_max = _dx/2;
    double dt1 = dt_max;

    // do one tmime only at the start of simulation
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
            // not faster, dt_max is fastest with the same error:
            //dt1 = min(dt1*qSqrt(double(n)), dt_max);
            //dt1 = min(dt1*(double(n)), dt_max);
            dt1 = dt_max;

            setZero(h, u, v);

            simpleScheme(h, u, v);

            dt1 = maincalcflux(dt1, dt_max);
            dt1 = min(dt1, _dt-timesum);

            maincalcscheme(dt1, h,u,v, hs,us,vs);

            setZero(hs, us, vs);

            FOR_ROW_COL_MV_MV
            {
                h->Drc = hs->Drc;
                u->Drc = us->Drc;
                v->Drc = vs->Drc;
                q1->Drc = h->Drc*u->Drc;
                q2->Drc = h->Drc*v->Drc;
            }

            timesum = timesum + dt1;
            n++;

        } while (timesum  < _dt);
    }
    //        Fr=froude_number(hs,us,vs);
    // todo

    iter_n = n;
    return(timesum/(n+1));
}
//---------------------------------------------------------------------------
double TWorld::fullSWOF2Do2(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2)
{
    double dt1, dt2, timesum = 0;
    double dt_max = _dx/2;

    if (prepareFlood)
        verif = 1;
    if (prepareFlood)
    {
        prepareFlood = false;
        FOR_ROW_COL_MV_MV
        {
            delta_z1->Drc = z->Data[r][c+1] - z->Drc;
            delta_z2->Drc = z->Data[r+1][c] - z->Drc;
            som_z1->Drc = z->Data[r][c-1]-2*z->Drc+z->Data[r][c+1];
            som_z2->Drc = z->Data[r-1][c]-2*z->Drc+z->Data[r+1][c];
        }
    }

    // if there is no flood skip everything
    int n = 1;
    if (startFlood)
    {
        do {
            dt1 = dt_max;
            dt2 = dt_max*0.5;

            setZero(h, u, v);
            // Reconstruction for order 2
            // makes h1r, h1l, u1r, u1l, v1r, v1l
            // makes h2r, h2l, u2r, u2l, v2r, v2l
            // makes delzc1, delzc2, delz1, delz2
            if (SwitchMUSCL)
                MUSCL(hs,us,vs,z);
            else
                ENO(hs,us,vs,z);

            // semi-iteration: optimize the timestep
            do {
                dt1 = maincalcflux(dt2, dt_max);
                dt1 = min(dt1, _dt-timesum);

                maincalcscheme(dt1, h,u,v, hs,us,vs);
                setZero(hs, us, vs);

                if (SwitchMUSCL)
                    MUSCL(hs,us,vs,z);
                else
                    ENO(hs,us,vs,z);

                dt2 = maincalcflux(dt1, dt_max);

            } while (dt2 < dt1);


            maincalcscheme(dt1, hs,us,vs, hsa, usa, vsa);
            setZero(hsa, usa, vsa);

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

            timesum = timesum + dt1;
            n++;


        } while (timesum  < _dt);

        //        FOR_ROW_COL_MV
        //        {
        //            if(u->Drc > 1000)
        //                u->Drc = v->Drc;
        //        }

        //        do {
        //            if (verif == 1)
        //            {
        //                dt1 = dt_max;

        //                setZero(h, u, v);
        //                // Reconstruction for order 2
        //                // makes h1r, h1l, u1r, u1l, v1r, v1l
        //                // makes h2r, h2l, u2r, u2l, v2r, v2l
        //                // makes delzc1, delzc2, delz1, delz2
        //                MUSCL(h,u,v,z);
        //                //      simpleScheme(h,u,v);
        //            }
        //            //            else
        //            //               ; //reset infil!!!

        //            dt1 = maincalcflux(dt1, dt_max);
        //            dt1 = min(dt1, _dt-timesum);

        //            //h, u, v go in hs, vs, us come out
        //            maincalcscheme(dt1, h,u,v, hs,us,vs);
        //            dt2 = dt1;

        //            setZero(hs, us, vs);

        //            //Reconstruction for order 2
        //            MUSCL(hs,us,vs,z);
        //            //  simpleScheme(hs,us,vs);
        //            //    hs->report("hs");

        //            dt2 = maincalcflux(dt2, dt_max);

        //            if (dt2 < dt1)
        //            {
        //                dt1 = dt2;
        //                verif = 0;
        //            }
        //            else
        //            {
        //                verif = 1;
        //                //hs, us, vs go in hsa, vsa, usa come out
        //                maincalcscheme(dt1, hs,us,vs, hsa, usa, vsa);

        //                setZero(hsa, usa, vsa);

        //                //Heun method (order 2 in time)
        //                FOR_ROW_COL_MV
        //                {
        //                    double tmp = 0.5*(h->Drc + hsa->Drc);
        //                    if (tmp >= he_ca)
        //                    {
        //                        q1->Drc = 0.5*(h->Drc*u->Drc + hsa->Drc*usa->Drc);
        //                        u->Drc = q1->Drc/tmp;
        //                        q2->Drc = 0.5*(h->Drc*v->Drc + hsa->Drc*vsa->Drc);
        //                        v->Drc = q2->Drc/tmp;
        //                        h->Drc = tmp;
        //                    }
        //                    else
        //                    {
        //                        u->Drc = 0;
        //                        q1->Drc = 0;
        //                        v->Drc = 0;
        //                        q2->Drc = 0;
        //                        h->Drc = 0;
        //                    }
        //                }//Heun

        //                timesum = timesum + dt1;
        //                n++;

        //            }//end for else dt2<dt1

        //        } while (timesum  < _dt);

    } // if floodstart

    iter_n = n;

    return(timesum/(n+1));

}
//---------------------------------------------------------------------------
