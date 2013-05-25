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
 \file lisSWOF2D.cpp
 \brief Channel flood using the opensource code from the FullSWOF_2D project.
        St Vennant equations, 1st order and 2n order solutions
        authors: Olivier Delestre, Christian Laguerre, Carine Lucas, Ulrich Razafison, Marie Rousseau.
        website: http://www.univ-orleans.fr/mapmo/soft/FullSWOF/

functions: \n
-   double fullSWOF2Do2(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2);
-   double fullSWOF2Do1(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2);
-   double maincalcflux(double dt, double dt_max);
-   void maincalcscheme(double dt, TMMap *he, TMMap *ve1, TMMap *ve2,TMMap *hes, TMMap *ves1, TMMap *ves2);
-   void MUSCL(TMMap *h,TMMap *u,TMMap *v,TMMap *z);
-   void ENO(TMMap *h,TMMap *u,TMMap *v,TMMap *z);
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
#define dt_fix 0.125

#define GRAV 9.8067
#define EPSILON 1e-6
#define scheme_type 1   //return calculated or fixed dt

//---------------------------------------------------------------------------
/**
 * @brief TWorld::limiter: Flux limiters are used in high resolution schemes to avoid occilations
 * @param a slope on one side
 * @param b slope on oposite side
 * @return rec
 *
 * WIKI: Flux limiters are used in high resolution schemes, such as the MUSCL scheme, to avoid
 * the spurious oscillations (wiggles) that would otherwise occur with high order spatial
 * discretisation schemes due to shocks, discontinuities or sharp changes in the solution domain.
 * Use of flux limiters, together with an appropriate high resolution scheme, make the solutions
 * total variation diminishing (TVD).
 */
double TWorld::limiter(double a, double b)
{
    double eps = 1.e-15;
    double rec = 0.;

    if (a == b)
        return(b);
    if (a*b == 0)
        return(0.);
    // save time

    if (SwitchLimiter == VANLEER)
    {
        if (a*b > 0.)
            return (2.*a*b/(a+b));
    }
    else
        if (SwitchLimiter == MINMOD)
        {
            if (a >= 0. && b >= 0.)
                rec = qMin(a, b);
            else
                if (a<=0. && b<=0)
                    rec = qMax(a, b);
        }
        else
            if (SwitchLimiter == VANALBEDA)
            {
                if (a*b >= 0.)
                    rec=(a*(b*b+eps)+b*(a*a+eps))/(a*a+b*b+2.*eps);
            }

    return(rec);
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
void TWorld::Fr_Manning(double u, double v, double hnew, double q1new, double q2new, double dt, double N)
{
    double nsq = N*N*GRAV*qSqrt(u*u+v*v)*dt/qPow(hnew,4./3.);
    q1man = q1new/(1.0+nsq);
    q2man = q2new/(1.0+nsq);
}

//NOT USED!
//Sf = Manning = v|v|/c^2*h^{4/3}
void TWorld::Fr_ManningSf(double h, double u, double v, double N)
{
    double nsq = N*N*qSqrt(u*u+v*v)/qPow(h,4./3.);
    if (h>he_ca)
    {
        q1man = u*nsq;
        q1man = v*nsq;
    }
    else
    {
        q1man = 0;
        q1man = 0;
    }
}
//---------------------------------------------------------------------------
/// Numerical flux calculation on which the new velocity is based
/// U_n+1 = U_n + dt/dx* [flux]  when flux is calculated by HLL, HLL2, Rusanov
/// HLL = Harten, Lax, van Leer numerical solution
/// TODO: 1/(c1-c2) can become 0?
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
/// Essentially Non-Oscillatory schemes (ENO)
/// second order scheme: based on [r, r+1, r+2] and [c, c+1, c+2]
/// http://en.wikipedia.org/wiki/Shock_capturing_method
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

    //x-direction
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

    //y-direction
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
        }

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
/// MUSCL: Monotone Upstream-centered Schemes for Conservation Laws
/// see http://en.wikipedia.org/wiki/MUSCL_scheme
void TWorld::MUSCL(TMMap *_h, TMMap *_u, TMMap *_v, TMMap *_z)
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

        delta_h2 = _h->Data[r][c+1] - _h->Drc;
        delta_u2 = _u->Data[r][c+1] - _u->Drc;
        delta_v2 = _v->Data[r][c+1] - _v->Drc;

        dh = limiter(delta_h1, delta_h2);
        dz_h = limiter(delta_h1 + delta_z1->Data[r][c-1],
                delta_h2 + delta_z1->Drc);
        du = limiter(delta_u1, delta_u2);
        dv = limiter(delta_v1, delta_v2);

        h1r->Drc = _h->Drc+dh*0.5;
        h1l->Drc = _h->Drc-dh*0.5;

        z1r->Drc = _z->Drc+0.5*(dz_h-dh);
        z1l->Drc = _z->Drc+0.5*(dh-dz_h);

        delzc1->Drc = z1r->Drc-z1l->Drc;
        delz1->Data[r][c-1] = z1l->Drc-z1r->Data[r][c-1];

        if (_h->Drc > he_ca)
        {
            u1r->Drc = _u->Drc + h1l->Drc*du*0.5/_h->Drc;
            u1l->Drc = _u->Drc - h1r->Drc*du*0.5/_h->Drc;
            v1r->Drc = _v->Drc + h1l->Drc*dv*0.5/_h->Drc;
            v1l->Drc = _v->Drc - h1r->Drc*dv*0.5/_h->Drc;
        }
        else
        {
            u1r->Drc = _u->Drc + du*0.5;
            u1l->Drc = _u->Drc - du*0.5;
            v1r->Drc = _v->Drc + dv*0.5;
            v1l->Drc = _v->Drc - dv*0.5;
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
    //    for (int c = 0; c < _nrCols; c++)
    //        for (int r = 1; r < _nrRows-1; r++)
    //            if(!IS_MV_REAL8(&LDD->Data[r][c]) &&
    //                    !IS_MV_REAL8(&LDD->Data[r-1][c]) &&
    //                    !IS_MV_REAL8(&LDD->Data[r+1][c]))
    FOR_ROW_COL_MV_MV
    {
        delta_h1 = tm->Drc;
        delta_u1 = tma->Drc;
        delta_v1 = tmb->Drc;

        delta_h2 = _h->Data[r+1][c] - _h->Drc;
        delta_u2 = _u->Data[r+1][c] - _u->Drc;
        delta_v2 = _v->Data[r+1][c] - _v->Drc;

        dh = limiter(delta_h1, delta_h2);
        dz_h = limiter(delta_h1+delta_z2->Data[r-1][c],
                delta_h2+delta_z2->Drc);

        du = limiter(delta_u1, delta_u2);
        dv = limiter(delta_v1, delta_v2);

        h2r->Drc = _h->Drc+dh*0.5;
        h2l->Drc = _h->Drc-dh*0.5;

        z2r->Drc = _z->Drc+0.5*(dz_h-dh);
        z2l->Drc = _z->Drc+0.5*(dh-dz_h);
        delzc2->Drc = z2r->Drc - z2l->Drc;
        delz2->Data[r-1][c] = z2l->Drc - z2r->Data[r-1][c];

        if (_h->Drc > he_ca)
        {
            u2r->Drc = _u->Drc + h2l->Drc*du*0.5/_h->Drc;
            u2l->Drc = _u->Drc - h2r->Drc*du*0.5/_h->Drc;
            v2r->Drc = _v->Drc + h2l->Drc*dv*0.5/_h->Drc;
            v2l->Drc = _v->Drc - h2r->Drc*dv*0.5/_h->Drc;
        }
        else
        {
            u2r->Drc = _u->Drc + du*0.5;
            u2l->Drc = _u->Drc - du*0.5;
            v2r->Drc = _v->Drc + dv*0.5;
            v2l->Drc = _v->Drc - dv*0.5;
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
        double dy = DX->Drc; //_dx;
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
            qes1 = he->Drc*ve1->Drc - ty*(g2->Data[r+1][c]-g2->Drc) -
                    tx*(f2->Data[r][c+1]-f2->Drc +
                    GRAV*0.5*((h1g->Drc-h1l->Drc)*(h1g->Drc+h1l->Drc) +
                              (h1r->Drc-h1d->Drc)*(h1r->Drc+h1d->Drc) +
                              (h1l->Drc+h1r->Drc)*delzc1->Drc)) ;

            // fullswof version 1.04
            qes2 = he->Drc*ve2->Drc - tx*(f3->Data[r][c+1]-f3->Drc) -
                    ty*(g3->Data[r+1][c]-g3->Drc +
                    GRAV*0.5*((h2g->Drc-h2l->Drc)*(h2g->Drc+h2l->Drc) +
                              (h2r->Drc-h2d->Drc)*(h2r->Drc+h2d->Drc) +
                              (h2l->Drc+h2r->Drc)*delzc2->Drc));

            //Calcul friction in semi-implicit with old v and u and new h
            Fr_Manning(ve1->Drc, ve2->Drc, hes->Drc, qes1, qes2, dt, N->Drc);

            ves1->Drc = q1man/hes->Drc;
            ves2->Drc = q2man/hes->Drc;
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
/**
 Construction of variables for hydrostatic reconstruction.
 Flux in x and y direction.
 Calculaton of the time steps in relation to cfl.
*/
double TWorld::maincalcflux(double dt, double dt_max)
{
    double dt_tmp, dtx, dty;
    double velocity_max_x, velocity_max_y;
    dtx = dt_max;
    dty = dt_max;
    velocity_max_x = -ve_ca;
    velocity_max_y = -ve_ca;
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
        cflx->Drc = HLL2_cfl;
    }


    FOR_ROW_COL_MV_MV
    {
        dy = DX->Drc;

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
        cfly->Drc = HLL2_cfl;
    }

    // VJ 130517: not in the original code!
    // correct sudden extreme alues, swap x or y direction
    // cfl = v+sqrt(v), cannot be extremely large
    FOR_ROW_COL_MV_MV
    {
        if (cflx->Drc > 100 || cfly->Drc > 100)
        {
            qDebug() << "oh oh";
            if (cflx->Drc > 100)
            {
                qDebug() << "mainflux x" << dt << dtx << cflx->Drc << cfly->Drc;
                cflx->Drc = cfly->Drc;
                f1->Drc = g1->Drc;
                f2->Drc = g2->Drc;
                f3->Drc = g3->Drc;
                //   qDebug() << "mainflux y" << dt << dty << velocity_max_y;

            }
            else
            {
                qDebug() << "mainflux y" << dt << dty << cflx->Drc << cfly->Drc;
                cfly->Drc = cflx->Drc;
                g1->Drc = f1->Drc;
                g2->Drc = f2->Drc;
                g3->Drc = f3->Drc;
            }

        }
    }

    // find largest velocity and determine dt
    FOR_ROW_COL_MV_MV
    {
        if (qFabs(cflx->Drc*dt/dx) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = cfl_fix*dx/cflx->Drc;

        dtx = min(min(dt, dt_tmp), dtx);
        velocity_max_x = max(velocity_max_x, cflx->Drc);
    }

    // find largest velocity and determine dt
    FOR_ROW_COL_MV_MV
    {
        if (qFabs(cfly->Drc*dt/dy) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = cfl_fix*dy/cfly->Drc;

        dty = min(min(dt, dt_tmp), dty);
        velocity_max_y = max(velocity_max_y, cfly->Drc);
    }

 //   qDebug() << "mainflux x" << dt << dtx << velocity_max_x;
 //   qDebug() << "mainflux y" << dt << dty << velocity_max_y;

    if (scheme_type == 1)
        return(max(dt_ca, min(dtx,dty)));
    else
    {
        //        if ((velocity_max_x*dt_fix/dx > cflfix)||(velocity_max_y*dt_fix/dy > cflfix)){
        //            qDebug() << "the CFL condition is not satisfied: CFL >"<<cflfix << endl;
        //        }
        return (dt_fix);
    }

}
//---------------------------------------------------------------------------
void TWorld::simpleScheme(TMMap *_h,TMMap *_u,TMMap *_v, TMMap *_z)
{

    FOR_ROW_COL_MV_MV
    {
        h1r->Drc = _h->Drc;
        u1r->Drc = _u->Drc;
        v1r->Drc = _v->Drc;
        h1l->Data[r][c+1] = _h->Data[r][c+1];
        u1l->Data[r][c+1] = _u->Data[r][c+1];
        v1l->Data[r][c+1] = _v->Data[r][c+1];

//        h1r->Data[r][c] = (3.*_h->Data[r][c] + _h->Data[r][c-1])/4.;
//        u1r->Data[r][c] = (3.*_u->Data[r][c] + _u->Data[r][c-1])/4.;
//        v1r->Data[r][c] = (3.*_v->Data[r][c] + _v->Data[r][c-1])/4.;
//        h1l->Data[r][c+1] = (3.*_h->Data[r][c+1] + _h->Data[r][c])/4.;
//        u1l->Data[r][c+1] = (3.*_u->Data[r][c+1] + _u->Data[r][c])/4.;
//        v1l->Data[r][c+1] = (3.*_v->Data[r][c+1] + _v->Data[r][c])/4.;

    }
    FOR_ROW_COL_MV_MV
    {
        h2r->Drc = _h->Drc;
        u2r->Drc = _u->Drc;
        v2r->Drc = _v->Drc;
        h2l->Data[r+1][c] = _h->Data[r+1][c];
        u2l->Data[r+1][c] = _u->Data[r+1][c];
        v2l->Data[r+1][c] = _v->Data[r+1][c];

//        h2l->Data[r][c] = (3*_h->Data[r][c] + _h->Data[r-1][c])/4.0;
//        u2l->Data[r][c] = (3*_u->Data[r][c] + _u->Data[r-1][c])/4.0;
//        v2l->Data[r][c] = (3*_v->Data[r][c] + _v->Data[r-1][c])/4.0;
//        h2l->Data[r+1][c] = (3*_h->Data[r+1][c] + _h->Drc)/4.0;
//        u2l->Data[r+1][c] = (3*_u->Data[r+1][c] + _u->Drc)/4.0;
//        v2l->Data[r+1][c] = (3*_v->Data[r+1][c] + _v->Drc)/4.0;

    }
}
//---------------------------------------------------------------------------
/**
 * @brief TWorld::fullSWOF2Do1: first order solution for the st Venant equations
 * @param h : flood water level (m)
 * @param u : velocity in x-direction(m/s)
 * @param v : velocity in y-direction(m/s)
 * @param z : DTM = DEM and obstacles
 * @param q1: flux in the x-direction(m2/s)
 * @param q2: flux in the y-direction(m2/s)
 * @return average dt in flood loop
 */
double TWorld::fullSWOF2Do1(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2)
{
    double timesum = 0;
    int n = 0;
    double dt_max = min(_dt, _dx/2);
    double dt1 = dt_max;

    // do one tmime only at the start of simulation
    if (prepareFlood)
    {
        prepareFlood = false;
        FOR_ROW_COL_MV_MV
        {
            delta_z1->Drc = z->Data[r][c+1] - z->Drc;
            delta_z2->Drc = z->Data[r+1][c] - z->Drc;

            delz1->Data[r][c-1] = z->Drc - z->Data[r][c-1];
            delz2->Data[r-1][c] = z->Drc - z->Data[r-1][c];
        }
    }

    // if there is no flood skip everything
    if (startFlood)
    {

        do {
            // not faster, dt_max is fastest with the same error:
           // dt1 = min(dt1*qSqrt(double(n)), dt_max);
            //dt1 = min(dt1*(double(n)), dt_max);
            dt1 = dt_max;

            setZero(h, u, v);

            //MUSCL(h,u,v,z);
            simpleScheme(h, u, v, z);

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
/**
 * @brief TWorld::fullSWOF2Do2: second order solution for the st Venant equations
 * @param h : flood water level (m)
 * @param u : velocity in x-direction(m/s)
 * @param v : velocity in y-direction(m/s)
 * @param z : DTM = DEM and obstacles
 * @param q1: flux in the x-direction(m2/s)
 * @param q2: flux in the y-direction(m2/s)
 * @return average dt in flood loop
 */
double TWorld::fullSWOF2Do2(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2)
{
    double dt1, dt2, timesum = 0;
    double dt_max = min(_dt, _dx/2);

    if (prepareFlood)
    {
        //verif = 1;
        prepareFlood = false;
        FOR_ROW_COL_MV_MV
        {
            delz1->Data[r][c-1] = z->Drc - z->Data[r][c-1];
            delz2->Data[r-1][c] = z->Drc - z->Data[r-1][c];
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
            dt2 = dt_max;

            setZero(h, u, v);
            // Reconstruction for order 2
            // makes h1r, h1l, u1r, u1l, v1r, v1l
            // makes h2r, h2l, u2r, u2l, v2r, v2l
            // makes delzc1, delzc2, delz1, delz2
            if (SwitchMUSCL)
                MUSCL(h,u,v,z);
            else
                ENO(h,u,v,z);
//            simpleScheme(h, u, v, z);
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
//                simpleScheme(hs,us,vs, z);
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
