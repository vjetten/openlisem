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
-   double fullSWOF2Do2(cTMap *h, cTMap *u, cTMap *v, cTMap *z, cTMap *q1, cTMap *q2);
-   double fullSWOF2Do1(cTMap *h, cTMap *u, cTMap *v, cTMap *z, cTMap *q1, cTMap *q2);
-   double maincalcflux(double dt, double dt_max);
-   void maincalcscheme(double dt, cTMap *he, cTMap *ve1, cTMap *ve2,cTMap *hes, cTMap *ves1, cTMap *ves2);
-   void MUSCL(cTMap *h,cTMap *u,cTMap *v,cTMap *z);
-   void ENO(cTMap *h,cTMap *u,cTMap *v,cTMap *z);
-   void F_HLL2(double hg,double ug,double vg,double hd,double ud,double vd);
-   void F_HLL(double hg,double ug,double vg,double hd,double ud,double vd);
-   void F_Rusanov(double hg,double ug,double vg,double hd,double ud,double vd);
-   void Fr_Manning(double uold, double vold, double hnew, double q1new, double q2new, double dt, double cf);
-   void Fr_ManningSf(double h, double u, double v, double cf);
-   void setZero(cTMap *h, cTMap *u, cTMap *v);//, cTMap *q1, cTMap *q2);
-   double limiter(double a, double b);
*/

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-10
#define ve_ca 1e-10

#define dt_ca 0.005

#define GRAV 9.8067
#define EPSILON 1e-6

//---------------------------------------------------------------------------
/**
 * @brief TWorld::limiter: Flux limiters are used in high resolution schemes to avoid occilations
 * @param a slope on one side
 * @param b slope on oposite side
 * @return rec
 *
 * ONLY used in MUSCL or ENO
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

  if (F_fluxLimiter == (int)MINMOD)
    {
      if (a >= 0 && b >= 0)
        rec = std::min(a, b);
      else
        if (a <= 0 && b <= 0)
          rec = std::max(a, b);
    }
  else
    {
      double ab = a*b;

      if (F_fluxLimiter == (int)VANLEER)
        {
          if (ab > 0)
            return (2*ab/(a+b));
        }
      else
        if (F_fluxLimiter == (int)VANALBEDA)
          {
            double aa = a*a;
            double bb = b*b;
            if (ab > 0)
              rec=(a*(bb+eps)+b*(aa+eps))/(aa+bb+2*eps);
          }
    }
  return(rec);
}
//---------------------------------------------------------------------------
void TWorld::setZero(cTMap *_h, cTMap *_u, cTMap *_v)
{
  FOR_CELL_IN_FLOODAREA  {
    if (_h->Drc <= he_ca)
      {
        _h->Drc = 0;
        _u->Drc = 0;
        _v->Drc = 0;
      }

    if (fabs(_u->Drc) <= ve_ca)
      {
        _u->Drc = 0;
      }
    if (fabs(_v->Drc) <= ve_ca)
      {
        _v->Drc = 0;
      }
  }}
}
//---------------------------------------------------------------------------
//friction slope
void TWorld::Fr_Manning(double u, double v, double hnew, double q1new, double q2new, double dt, double N)
{
  double nsq = N*N*GRAV*sqrt(u*u+v*v)*dt/qPow(hnew,4./3.);
  q1man = q1new/(1.0+nsq);
  q2man = q2new/(1.0+nsq);
}

//NOT USED!
//Sf = Manning = v|v|/c^2*h^{4/3}
void TWorld::Fr_ManningSf(double h, double u, double v, double N)
{
  double nsq = N*N*sqrt(u*u+v*v)/qPow(h,4./3.);
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
    double f1, f2, f3, cfl, tmp = 0;
    //if (h_L<=0. && h_R<=0.){
    if (h_L<=he_ca && h_R<=he_ca){
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }
    else
    {
        double grav_h_L = GRAV*h_L;
        double grav_h_R = GRAV*h_R;
        double sqrt_grav_h_L = sqrt(grav_h_L);  // wave velocity
        double sqrt_grav_h_R = sqrt(grav_h_R);
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;

        double c1 = std::min(u_L - sqrt_grav_h_L,u_R - sqrt_grav_h_R); //we already have u_L - sqrt_grav_h_L<u_L + sqrt_grav_h_L and u_R - sqrt_grav_h_R<u_R + sqrt_grav_h_R
        double c2 = std::max(u_L + sqrt_grav_h_L,u_R + sqrt_grav_h_R); //so we do not need all the eigenvalues to get c1 and c2
        tmp = 1./(c2-c1);
        double t1 = (std::min(c2,0.) - std::min(c1,0.))*tmp;
        double t2 = 1. - t1;
        double t3 = (c2*fabs(c1) - c1*fabs(c2))*0.5*tmp;

        f1 = t1*q_R + t2*q_L - t3*(h_R - h_L);
        f2 = t1*(q_R*u_R + grav_h_R*h_R*0.5) + t2*(q_L*u_L + grav_h_L*h_L*0.5) - t3*(q_R - q_L);
        f3 = t1*q_R*v_R + t2*q_L*v_L - t3*(h_R*v_R - h_L*v_L);
        cfl = std::max(fabs(c1),fabs(c2)); //cfl is the velocity to compute the cfl condition std::max(fabs(c1),fabs(c2))*tx with tx=dt/dx
    }
    HLL2_cfl = cfl;
    HLL2_f1 = f1;
    HLL2_f2 = f2;
    HLL2_f3 = f3;
}

void TWorld::F_HLL(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    double f1, f2, f3, cfl, tmp = 0;
  //  if (h_L<=0. && h_R<=0.){
    if (h_L<=he_ca && h_R<=he_ca){
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }else{
        double grav_h_L = GRAV*h_L;
        double grav_h_R = GRAV*h_R;
        double q_R = u_R*h_R;
        double q_L = u_L*h_L;
        double c1 = std::min(u_L-sqrt(grav_h_L),u_R-sqrt(grav_h_R));
        double c2 = std::max(u_L+sqrt(grav_h_L),u_R+sqrt(grav_h_R));

        //cfl is the velocity to calculate the real cfl=std::max(fabs(c1),fabs(c2))*tx with tx=dt/dx
        if (fabs(c1)<EPSILON && fabs(c2)<EPSILON){              //dry state
            f1=0.;
            f2=0.;
            f3=0.;
            cfl=0.; //std::max(fabs(c1),fabs(c2))=0
        }else if (c1>=EPSILON){ //supercritical flow, from left to right : we have std::max(abs(c1),abs(c2))=c2>0
            f1=q_L;   //flux
            f2=q_L*u_L+GRAV*h_L*h_L*0.5;  //flux*velocity + 0.5*(wave velocity squared)
            f3=q_L*v_L; //flux *velocity
            cfl=c2; //std::max(fabs(c1),fabs(c2))=c2>0
        }else if (c2<=-EPSILON){ //supercritical flow, from right to left : we have std::max(abs(c1),abs(c2))=-c1>0
            f1=q_R;
            f2=q_R*u_R+GRAV*h_R*h_R*0.5;
            f3=q_R*v_R;
            cfl=fabs(c1); //std::max(fabs(c1),fabs(c2))=fabs(c1)
        }else{ //subcritical flow
            tmp = 1./(c2-c1);
            f1=(c2*q_L-c1*q_R)*tmp+c1*c2*(h_R-h_L)*tmp;
            f2=(c2*(q_L*u_L+GRAV*h_L*h_L*0.5)-c1*(q_R*u_R+GRAV*h_R*h_R*0.5))*tmp+c1*c2*(q_R-q_L)*tmp;
            f3=(c2*(q_L*v_L)-c1*(q_R*v_R))*tmp+c1*c2*(h_R*v_R-h_L*v_L)*tmp;
            cfl=std::max(fabs(c1),fabs(c2));
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
//    if (h_L<=0. && h_R<=0.){
   if (h_L<=he_ca && h_R<=he_ca){
        c = 0.;
        f1 = 0.;
        f2 = 0.;
        f3 = 0.;
        cfl = 0.;
    }else{
        c = std::max(fabs(u_L)+sqrt(GRAV*h_L),fabs(u_R)+sqrt(GRAV*h_R));
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
/// Called reconstruction schemes because they reconstruct the flux through
/// the boundary of adjacent cells frm the values at the cell centres
/// second order scheme: based on [r, r+1, r+2] and [c, c+1, c+2]
/// http://en.wikipedia.org/wiki/Shock_capturing_method
void TWorld::ENO(cTMap *h,cTMap *u,cTMap *v,cTMap *z)
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
    double MODIFENO = 0.9;
    double a1, a2, a3, a4;

  //x-direction
    FOR_CELL_IN_FLOODAREA
      if(c > 0 && c < _nrCols-2 && !MV(r,c-1) && !MV(r, c+1) && !MV(r, c+2))
    {
      hh1 = h->data[r][c-1]-2.*h->data[r][c]+h->data[r][c+1];
      uu1 = u->data[r][c-1]-2.*u->data[r][c]+u->data[r][c+1];
      vv1 = v->data[r][c-1]-2.*v->data[r][c]+v->data[r][c+1];

      hh2 = h->Drc-2.*h->data[r][c+1]+h->data[r][c+2];
      uu2 = u->Drc-2.*u->data[r][c+1]+u->data[r][c+2];
      vv2 = v->Drc-2.*v->data[r][c+1]+v->data[r][c+2];

      ddh2 = amortENO*limiter(hh1,hh2);
      ddu2 = amortENO*limiter(uu1,uu2);
      ddz2 = amortENO*limiter(hh1+som_z1->Drc,hh2+som_z1->data[r][c+1]);
      ddv2 = amortENO*limiter(vv1,vv2);

      delta_h1 = h->Drc - h->data[r][c-1];
      delta_u1 = u->Drc - u->data[r][c-1];
      delta_v1 = v->Drc - v->data[r][c-1];
      delta_h2 = h->data[r][c+1]-h->Drc;
      delta_u2 = u->data[r][c+1]-u->Drc;
      delta_v2 = v->data[r][c+1]-v->Drc;


      if (F_scheme == (int)FENO)
      {
          dh = limiter(delta_h1+ddh1*0.5,delta_h2-ddh2*0.5);
          dz_h = limiter(delta_h1+delta_z1->data[r][c-1]+ddz1*0.5, delta_h2+delta_z1->Drc-ddz2*0.5);
      }
      else
      {
          a2 = limiter(delta_h1+ddh1*0.5,delta_h2-ddh2*0.5);
          a4 = limiter(delta_h1+delta_z1->data[r][c-1]+ddz1*0.5,delta_h2+delta_z1->Drc-ddz2*0.5);

          a1 = limiter(delta_h1,delta_h2);
          dh = limiter(2*MODIFENO*a1,a2);

          a3 = limiter(delta_h1+delta_z1->data[r][c-1],delta_h2+delta_z1->Drc);
          dz_h = limiter(2*MODIFENO*a3,a4);
      }
      du = limiter(delta_u1+ddu1*0.5,delta_u2-ddu2*0.5);
      dv = limiter(delta_v1+ddv1*0.5,delta_v2-ddv2*0.5);

      h1r->Drc=h->Drc+dh*0.5;
      h1l->Drc=h->Drc-dh*0.5;

      z1r->Drc=z->Drc+0.5*(dz_h-dh);
      z1l->Drc=z->Drc+0.5*(dh-dz_h);

      delzc1->Drc = z1r->Drc - z1l->Drc;
      delz1->data[r][c-1] = z1l->Drc - z1r->data[r][c-1];

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

    }}


    //y-direction
    FOR_CELL_IN_FLOODAREA
    if(r > 0 && r < _nrRows-2 && !MV(r-1,c) && !MV(r+1, c) && !MV(r+2, c))
    {
        hh1 = h->data[r-1][c]-2.*h->data[r][c]+h->data[r+1][c];
        uu1 = u->data[r-1][c]-2.*u->data[r][c]+u->data[r+1][c];
        vv1 = v->data[r-1][c]-2.*v->data[r][c]+v->data[r+1][c];
        hh2 = h->Drc-2.*h->data[r+1][c]+h->data[r+2][c];
        uu2 = u->Drc-2.*u->data[r+1][c]+u->data[r+2][c];
        vv2 = v->Drc-2.*v->data[r+1][c]+v->data[r+2][c];

        ddh2 = amortENO*limiter(hh1,hh2);
        ddu2 = amortENO*limiter(uu1,uu2);
        ddz2 = amortENO*limiter(hh1+som_z2->Drc,hh2+som_z2->data[r+1][c]);
        ddv2 = amortENO*limiter(vv1,vv2);

        delta_h1 = h->Drc - h->data[r-1][c];
        delta_u1 = u->Drc - u->data[r-1][c];
        delta_v1 = v->Drc - v->data[r-1][c];
        delta_h2 = h->data[r+1][c]-h->Drc;
        delta_u2 = u->data[r+1][c]-u->Drc;
        delta_v2 = v->data[r+1][c]-v->Drc;

        if (F_scheme == (int)FENO)
        {
            dh = limiter(delta_h1+ddh1*0.5,delta_h2-ddh2*0.5);
            dz_h = limiter(delta_h1+delta_z2->data[r-1][c]+ddz1*0.5,delta_h2+delta_z2->Drc-ddz2*0.5);
        }
        else
        {
            a2 = limiter(delta_h1+ddh1*0.5,delta_h2-ddh2*0.5);
            a4 = limiter(delta_h1+delta_z1->data[r-1][c]+ddz1*0.5,delta_h2+delta_z1->Drc-ddz2*0.5);

            a1 = limiter(delta_h1,delta_h2);
            dh = limiter(2*MODIFENO*a1,a2);

            a3 = limiter(delta_h1+delta_z1->data[r-1][c],delta_h2+delta_z1->Drc);
            dz_h = limiter(2*MODIFENO*a3,a4);
        }

        du = limiter(delta_u1+ddu1*0.5,delta_u2-ddu2*0.5);
        dv = limiter(delta_v1+ddv1*0.5,delta_v2-ddv2*0.5);

        h2r->Drc = h->Drc+dh*0.5;
        h2l->Drc = h->Drc-dh*0.5;

        z2r->Drc = z->Drc+0.5*(dz_h-dh);
        z2l->Drc = z->Drc+0.5*(dh-dz_h);
        delzc2->Drc = z2r->Drc-z2l->Drc;
        delz2->data[r-1][c] = z2l->Drc-z2r->data[r-1][c];

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

        ddh1=ddh2;
        ddz1=ddz2;
        ddu1=ddu2;
        ddv1=ddv2;

    }} //end for
}
//---------------------------------------------------------------------------
/// MUSCL: Monotone Upstream-centered Schemes for Conservation Laws
/// see http://en.wikipedia.org/wiki/MUSCL_scheme
///
/// This MUSCL creates left and right u,v,h, arrays for in mainfluxcalc
void TWorld::MUSCL(cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z)
{
  double delta_h1, delta_u1, delta_v1;
  double delta_h2, delta_u2, delta_v2;
  double dh, du, dv, dz_h;

  // fill EW and NS conditions with cell itself, 1st order approximation used for boundary conditions
  FOR_CELL_IN_FLOODAREA {
      h1r->Drc = _h->Drc;
      u1r->Drc = _u->Drc;
      v1r->Drc = _v->Drc;
      h1l->Drc = _h->Drc;
      u1l->Drc = _u->Drc;
      v1l->Drc = _v->Drc;

      h2r->Drc = _h->Drc;
      u2r->Drc = _u->Drc;
      v2r->Drc = _v->Drc;
      h2l->Drc = _h->Drc;
      u2l->Drc = _u->Drc;
      v2l->Drc = _v->Drc;
  }}


  FOR_CELL_IN_FLOODAREA
    if(c > 0 && c < _nrCols-1 && !MV(r,c-1) &&  !MV(r,c+1))
  {
      delta_h1 = _h->Drc - _h->data[r][c-1];
      delta_u1 = _u->Drc - _u->data[r][c-1];
      delta_v1 = _v->Drc - _v->data[r][c-1];

      delta_h2 = _h->data[r][c+1] - _h->Drc;
      delta_u2 = _u->data[r][c+1] - _u->Drc;
      delta_v2 = _v->data[r][c+1] - _v->Drc;

      dh   = 0.5*limiter(delta_h1, delta_h2);
      dz_h = 0.5*limiter(delta_h1 + delta_z1->data[r][c-1], delta_h2 + delta_z1->Drc);
      du   = 0.5*limiter(delta_u1, delta_u2);
      dv   = 0.5*limiter(delta_v1, delta_v2);

      h1r->Drc = _h->Drc+dh;
      h1l->Drc = _h->Drc-dh;

      z1r->Drc = _z->Drc+(dz_h-dh);
      z1l->Drc = _z->Drc+(dh-dz_h);

      delzc1->Drc = (long double)z1r->Drc-(long double)z1l->Drc;
      delz1->data[r][c-1] = z1l->Drc - z1r->data[r][c-1];

      if (_h->Drc > 0.)//he_ca)
      {
          double h1lh = h1l->Drc/_h->Drc;
          double h1rh = h1r->Drc/_h->Drc;

          u1r->Drc = _u->Drc + h1lh * du;
          u1l->Drc = _u->Drc - h1rh * du;
          v1r->Drc = _v->Drc + h1lh * dv;
          v1l->Drc = _v->Drc - h1rh * dv;
      }
      else
      {
          u1r->Drc = _u->Drc + du;
          u1l->Drc = _u->Drc - du;
          v1r->Drc = _v->Drc + dv;
          v1l->Drc = _v->Drc - dv;
      }
  }}


    FOR_CELL_IN_FLOODAREA
      if(r > 0 && MV(r-1,c) && r < _nrRows-1 && !MV(r+1, c))
    {
        delta_h1 = _h->Drc - _h->data[r-1][c];
        delta_u1 = _u->Drc - _u->data[r-1][c];
        delta_v1 = _v->Drc - _v->data[r-1][c];

        delta_h2 = _h->data[r+1][c] - _h->Drc;
        delta_u2 = _u->data[r+1][c] - _u->Drc;
        delta_v2 = _v->data[r+1][c] - _v->Drc;

        dh   = 0.5*limiter(delta_h1, delta_h2);
        dz_h = 0.5*limiter(delta_h1+delta_z2->data[r-1][c],delta_h2+delta_z2->Drc);
        du   = 0.5*limiter(delta_u1, delta_u2);
        dv   = 0.5*limiter(delta_v1, delta_v2);

        h2r->Drc = _h->Drc+dh;
        h2l->Drc = _h->Drc-dh;

        z2r->Drc = _z->Drc+(dz_h-dh);
        z2l->Drc = _z->Drc+(dh-dz_h);
        delzc2->Drc = (long double)z2r->Drc - (long double)z2l->Drc;
        delz2->data[r-1][c] = z2l->Drc - z2r->data[r-1][c];

        if (_h->Drc > 0.)
        {
            double h2lh = h2l->Drc/_h->Drc;
            double h2rh = h2r->Drc/_h->Drc;

            u2r->Drc = _u->Drc + h2lh * du;
            u2l->Drc = _u->Drc - h2rh * du;
            v2r->Drc = _v->Drc + h2lh * dv;
            v2l->Drc = _v->Drc - h2rh * dv;
        }
        else
        {
            u2r->Drc = _u->Drc + du;
            u2l->Drc = _u->Drc - du;
            v2r->Drc = _v->Drc + dv;
            v2l->Drc = _v->Drc - dv;
        }
    }}
}
//---------------------------------------------------------------------------
// St Venant equations: conservation of mass and momentum: h u et v
void TWorld::maincalcscheme(double dt, cTMap *he, cTMap *ve1, cTMap *ve2,
                            cTMap *hes, cTMap *ves1, cTMap *ves2)
{
  FOR_CELL_IN_FLOODAREA
  {
    double dx = ChannelAdj->Drc;
    double dy = DX->Drc;
    long double tx = dt/dx;
    long double ty = dt/dy;
    // Solution of the equation of mass conservation (First equation of Saint venant)
    // f1 comes from MUSCL calculations
    if ((r > _nrRows-2 || c > _nrCols-2) || (pcr::isMV(LDD->data[r][c+1]) || pcr::isMV(LDD->data[r+1][c])))
      hes->Drc = he->Drc;
    else
      hes->Drc = he->Drc - tx*(f1->data[r][c+1]-f1->Drc) - ty*(g1->data[r+1][c]-g1->Drc);
    if (hes->Drc > he_ca)
      {
        //Solution of the equation of momentum (Second and third equation of Saint-venant)
        double qes1;
        double qes2;

        // fullswof version 1.04
        // Solution of the equation of momentum (Second and third equation of Shallow-Water)
        // This expression for the flux (instead of the differences of the squares) avoids numerical errors
        // see http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html section "Cancellation".
        qes1 = (long double)he->Drc*(long double)ve1->Drc -
            ty*((long double)g2->data[r+1][c]-(long double)g2->Drc) -
            tx*((long double)f2->data[r][c+1]-(long double)f2->Drc +
            GRAV*0.5*(((long double)h1g->Drc-(long double)h1l->Drc)*((long double)h1g->Drc+(long double)h1l->Drc) +
                      ((long double)h1r->Drc-(long double)h1d->Drc)*((long double)h1r->Drc+(long double)h1d->Drc) +
                      ((long double)h1l->Drc+(long double)h1r->Drc)*(long double)delzc1->Drc)) ;

        // fullswof version 1.04
        qes2 = (long double)he->Drc*(long double)ve2->Drc - tx*((long double)f3->data[r][c+1]-(long double)f3->Drc) -
            ty*((long double)g3->data[r+1][c]-(long double)g3->Drc +
            GRAV*0.5*(((long double)h2g->Drc-(long double)h2l->Drc)*((long double)h2g->Drc+(long double)h2l->Drc) +
                      ((long double)h2r->Drc-(long double)h2d->Drc)*((long double)h2r->Drc+(long double)h2d->Drc) +
                      ((long double)h2l->Drc+(long double)h2r->Drc)*(long double)delzc2->Drc));


        //Calc friction semi-implicit with old v and u and new h
        //Fr_Manning(ve1->Drc, ve2->Drc, hes->Drc, qes1, qes2, dt, N->Drc);

        double nsq = N->Drc*N->Drc*GRAV*sqrt(ve1->Drc*ve1->Drc+ve2->Drc*ve2->Drc)*dt/qPow(hes->Drc,4.0/3.0);

        ves1->Drc = (qes1/(1.0+nsq))/hes->Drc;
        ves2->Drc = (qes2/(1.0+nsq))/hes->Drc;
      }
    else
      {
        // Case of height of water is zero.
        ves1->Drc = 0;
        ves2->Drc = 0;
        hes->Drc = 0;
      }
  }}
}
//---------------------------------------------------------------------------
/**
 Construction of variables for hydrostatic reconstruction using HLL or Rusanov
 Flux in x and y direction. The method is based solving PDEs with an estimate of the status
 of a domain based on spatial averages of the previous timestep, this is the hydrostatic equlibrium
 Calculaton of the time steps in relation to cfl.
*/
double TWorld::maincalcflux(double dt, double dt_max)
{
  double dt_tmp, dtx, dty;
  dtx = dt_max;
  dty = dt_max;

  FOR_CELL_IN_FLOODAREA
      if (c > 0 && !pcr::isMV(LDD->data[r][c-1]))
  {
    h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1] + std::max(FlowBarrierW->Drc,FlowBarrierE->data[r][c-1])));
    h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1] + std::max(FlowBarrierW->Drc,FlowBarrierE->data[r][c-1])));
    if (F_scheme == 1)
      F_Rusanov(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
    else
      if (F_scheme == 2)
        F_HLL(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
      else
        F_HLL2(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);

    f1->Drc = (std::abs(HLL2_f1) < he_ca ? 0: HLL2_f1);
    f2->Drc = (std::abs(HLL2_f2) < he_ca ? 0: HLL2_f2);
    f3->Drc = (std::abs(HLL2_f3) < he_ca ? 0: HLL2_f3);
//    f1->Drc = HLL2_f1;
//    f2->Drc = HLL2_f2;
//    f3->Drc = HLL2_f3;
    cflx->Drc = HLL2_cfl;
  }
  else
  {
    h1d->data[r][c] = std::max(0.0, h1r->data[r][c] - std::max(0.0,  delz1->data[r][c]));
    h1g->Drc        = std::max(0.0, h1l->Drc        - std::max(0.0, -delz1->data[r][c]));
    if (F_scheme == 1)
      F_Rusanov(h1d->data[r][c], u1r->data[r][c], v1r->data[r][c],h1g->Drc, u1l->Drc, v1l->Drc);
    else
      if (F_scheme == 2)
        F_HLL(h1d->data[r][c], u1r->data[r][c], v1r->data[r][c],h1g->Drc, u1l->Drc, v1l->Drc);
      else
        F_HLL2(h1d->data[r][c], u1r->data[r][c], v1r->data[r][c],h1g->Drc, u1l->Drc, v1l->Drc);

    f1->Drc = (std::abs(HLL2_f1) < he_ca ? 0: HLL2_f1);
    f2->Drc = (std::abs(HLL2_f2) < he_ca ? 0: HLL2_f2);
    f3->Drc = (std::abs(HLL2_f3) < he_ca ? 0: HLL2_f3);
//    f1->Drc = HLL2_f1;
//    f2->Drc = HLL2_f2;
//    f3->Drc = HLL2_f3;
    cflx->Drc = HLL2_cfl;
  }}

  FOR_CELL_IN_FLOODAREA
  if(r > 0 && !pcr::isMV(LDD->data[r-1][c]))
  {

    h2d->data[r-1][c] = std::max(0.0, h2r->data[r-1][c] - std::max(0.0,  delz2->data[r-1][c] + std::max(FlowBarrierS->Drc,FlowBarrierN->data[r-1][c])));
    h2g->Drc          = std::max(0.0, h2l->Drc          - std::max(0.0, -delz2->data[r-1][c] + std::max(FlowBarrierS->Drc,FlowBarrierN->data[r-1][c])));

    if (F_scheme == 1)
      F_Rusanov(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);
    else
      if (F_scheme == 2)
        F_HLL(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);
      else
        F_HLL2(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);

//    g1->Drc = HLL2_f1;
//    g2->Drc = HLL2_f3;
//    g3->Drc = HLL2_f2;
    g1->Drc = (std::abs(HLL2_f1) < he_ca ? 0: HLL2_f1);
    g2->Drc = (std::abs(HLL2_f3) < he_ca ? 0: HLL2_f3);
    g3->Drc = (std::abs(HLL2_f2) < he_ca ? 0: HLL2_f2);
    cfly->Drc = HLL2_cfl;
  }
  else
  {
     h2d->data[r][c] = std::max(0.0, h2r->data[r][c] - std::max(0.0,  delz2->data[r][c]));
     h2g->Drc        = std::max(0.0, h2l->Drc        - std::max(0.0, -delz2->data[r][c]));

     if (F_scheme == 1)
     F_Rusanov(h2d->data[r][c],v2r->data[r][c],u2r->data[r][c],h2g->Drc,v2l->Drc,u2l->Drc);
     else
     if (F_scheme == 2)
     F_HLL(h2d->data[r][c],v2r->data[r][c],u2r->data[r][c],h2g->Drc,v2l->Drc,u2l->Drc);
     else
     F_HLL2(h2d->data[r][c],v2r->data[r][c],u2r->data[r][c],h2g->Drc,v2l->Drc,u2l->Drc);

//     g1->Drc = HLL2_f1;
//     g2->Drc = HLL2_f3;
//     g3->Drc = HLL2_f2;
     g1->Drc = (std::abs(HLL2_f1) < he_ca ? 0: HLL2_f1);
     g2->Drc = (std::abs(HLL2_f3) < he_ca ? 0: HLL2_f3);
     g3->Drc = (std::abs(HLL2_f2) < he_ca ? 0: HLL2_f2);
     cfly->Drc = HLL2_cfl;
  }}

  // VJ 130517: not in the original code!
  // correct sudden extreme values, swap x or y direction
  // cfl = v+sqrt(v), cannot be extremely large such as 100 m/s!
/*
#define AVG(a1,a2) ((a1*a2 > 0) ? std::sqrt(a1*a2) : -std::sqrt(std::abs(a1*a2)))

     if (F_replaceV > 0)
     {
         FOR_CELL_IN_FLOODAREA {
             if (cflx->Drc > F_maxVelocity+qSqrt(GRAV* h1d->Drc) || cflx->Drc > F_maxVelocity+qSqrt(GRAV*h2d->Drc))
             {
                 double tmp1 = cflx->Drc;
                 double tmp2 = cfly->Drc;

                 double e1 = AVG(g1->Drc,f1->Drc);
                 double e2 = AVG(g2->Drc,f2->Drc);
                 double e3 = AVG(g3->Drc,f3->Drc);
                 double cfle = AVG(cfly->Drc,cflx->Drc);

                 if (cflx->Drc > F_maxVelocity)
                 {
                     cflx->Drc = cfle;
                     f1->Drc = e1;
                     f2->Drc = e2;
                     f3->Drc = e3;
                 }
                 if (cfly->Drc > F_maxVelocity)
                 {
                     cfly->Drc = cfle;
                     g1->Drc = e1;
                     g2->Drc = e2;
                     g3->Drc = e3;
                 }

                 qDebug() << "swap extreme velocity with avg XY" << tmp1 << tmp2 << cflx->Drc << cfly->Drc << r << c;
             }
         }}
    }
*/
     // find largest velocity and determine dt
     FOR_CELL_IN_FLOODAREA {
         double dx = ChannelAdj->Drc;
         if (qFabs(cflx->Drc*dt/dx) < 1e-10)
             dt_tmp = dt_max;
         else
             dt_tmp = courant_factor*dx/cflx->Drc;
         // was cfl_fix = 0.4
         dtx = std::min(std::min(dt, dt_tmp), dtx);
     }}


  // find largest velocity and determine dt
     FOR_CELL_IN_FLOODAREA {
         double dy = DX->Drc;
         if (qFabs(cfly->Drc*dt/dy) < 1e-10)
             dt_tmp = dt_max;
         else
             dt_tmp = courant_factor*dy/cfly->Drc;
         dty = std::min(std::min(dt, dt_tmp), dty);
     }}

     return(std::max(dt_ca, std::min(dtx,dty)));

}
//---------------------------------------------------------------------------

// fill the left and right h,u,v arrays to solve the 1D scheme, also used as prep for the 2D MUSCL scheme for the boundary cells
void TWorld::simpleScheme(cTMap *_h,cTMap *_u,cTMap *_v)
{

    FOR_CELL_IN_FLOODAREA {
        h1r->Drc = _h->Drc;
        u1r->Drc = _u->Drc;
        v1r->Drc = _v->Drc;
        h1l->Drc = _h->Drc;
        u1l->Drc = _u->Drc;
        v1l->Drc = _v->Drc;

        h2r->Drc = _h->Drc;
        u2r->Drc = _u->Drc;
        v2r->Drc = _v->Drc;
        h2l->Drc = _h->Drc;
        u2l->Drc = _u->Drc;
        v2l->Drc = _v->Drc;
    }}
}
//---------------------------------------------------------------------------
void TWorld::prepareFloodZ(cTMap *z)
{
    prepareFlood = false;

    fill(*delz1,-9999);
    fill(*delz2,-9999);

    // diff between z cell and adjacent
    for (int r = 0; r < _nrRows; r++)
        for (int c = 1; c < _nrCols; c++)
            if(!pcr::isMV(LDD->data[r][c]) &&
                    !pcr::isMV(LDD->data[r][c-1]))
            {
                delz1->data[r][c-1] = z->Drc - z->data[r][c-1];
                // needed in maincalcflux for 1D scheme, is calculated in MUSCL for 2D scheme
            }
    for (int r = 1; r < _nrRows; r++)
        for (int c = 0; c < _nrCols; c++)
            if(!pcr::isMV(LDD->data[r][c]) &&
                    !pcr::isMV(LDD->data[r-1][c]))
            {
                delz2->data[r-1][c] = z->Drc - z->data[r-1][c];
                // needed in maincalcflux for 1D scheme, is calculated in MUSCL for 2D scheme
            }

    // edges of domain: take cell next to edge cell or 0
    FOR_ROW_COL_MV
    {
        if (delz1->Drc == -9999)
        {
            if (!pcr::isMV(LDD->data[r][c-1]))
                delz1->Drc = delz1->data[r][c-1];
            else
                delz1->Drc = 0 ;
        }
        if (delz2->Drc == -9999)
        {
            if (!pcr::isMV(LDD->data[r-1][c]))
                delz2->Drc = delz2->data[r-1][c];
            else
                delz2->Drc = 0 ;
        }
    }

    fill(*delta_z1, -9999);
    fill(*delta_z2, -9999);
    for (int r = 0; r < _nrRows; r++)
        for (int c = 0; c < _nrCols-1; c++)
            if(!pcr::isMV(LDD->data[r][c]) &&
                    !pcr::isMV(LDD->data[r][c+1]))
                delta_z1->Drc = z->data[r][c+1] - z->Drc;

    for (int r = 0; r < _nrRows-1; r++)
        for (int c = 0; c < _nrCols; c++)
            if(!pcr::isMV(LDD->data[r][c]) &&
                    !pcr::isMV(LDD->data[r+1][c]))
                delta_z2->Drc = z->data[r+1][c] - z->Drc;

    FOR_ROW_COL_MV
    {
        if (delta_z1->Drc == -9999)
        {
            if (!pcr::isMV(LDD->data[r][c-1]))
                delta_z1->Drc = delta_z1->data[r][c-1];
            else
                delta_z1->Drc = 0 ;
        }
        if (delta_z2->Drc == -9999)
        {
            if (!pcr::isMV(LDD->data[r-1][c]))
                delta_z2->Drc = delta_z2->data[r-1][c];
            else
                delta_z2->Drc = 0 ;
        }
    }

    copy(*som_z1, *z);
    copy(*som_z2, *z);

    for (int r = 1; r < _nrRows-1; r++)
        for (int c = 1; c < _nrCols-1; c++)
            if(!pcr::isMV(LDD->data[r][c]) &&
                    !pcr::isMV(LDD->data[r-1][c]) &&
                    !pcr::isMV(LDD->data[r+1][c]) &&
                    !pcr::isMV(LDD->data[r][c-1]) &&
                    !pcr::isMV(LDD->data[r][c+1]))
            {
                som_z1->Drc = z->data[r][c-1]-2*z->Drc+z->data[r][c+1];
                som_z2->Drc = z->data[r-1][c]-2*z->Drc+z->data[r+1][c];
                // needed in ENO
            }
}
//---------------------------------------------------------------------------
/**
 * @brief TWorld::fullSWOF2Do1: first order solution for the st Venant equations
 * @param h : flood water level (m)
 * @param u : velocity in x-direction(m/s)
 * @param v : velocity in y-direction(m/s)
 * @param z : DTM = DEM and obstacles
 * @return average dt in flood loop
 */

// * @param q1: flux in the x-direction(m2/s)
// * @param q2: flux in the y-direction(m2/s)
double TWorld::fullSWOF2Do1(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct)
{
  double timesum = 0;
  int n = 1;
  double dt_max = std::min(_dt, _dx*0.5);
  double dt1 = dt_max;
  double sumh = 0;

  // do one tmime only at the start of simulation
  if (prepareFlood)
      prepareFloodZ(z);

  // if there is no flood skip everything
  if (startFlood)
  {
      if (correct)
          sumh = getMass(h);

      do {

          dt1 = dt_max;

          setZero(h, u, v);

          simpleScheme(h, u, v);

          dt1 = maincalcflux(dt1, dt_max);
          dt1 = std::min(dt1, _dt-timesum);

          //sediment
          SWOFSediment(dt1,h,u,v);

          maincalcscheme(dt1, h,u,v, hs,us,vs);

          setZero(hs, us, vs);

          FOR_CELL_IN_FLOODAREA {
            h->Drc = hs->Drc;
            u->Drc = us->Drc;
            v->Drc = vs->Drc;
          }}

          timesum = timesum + dt1;
          n++;

          if (correct)
              correctMassBalance(sumh, h, 1e-12);

          if (n > F_MaxIter)
            break;

        } while (timesum  < _dt);
    }

    iter_n = n;
    dt1 = n > 0? _dt/n : dt1;

    return(dt1);
}
//---------------------------------------------------------------------------
/**
 * @brief TWorld::fullSWOF2Do2: second order solution for the st Venant equations
 * @param h : flood water level (m)
 * @param u : velocity in x-direction(m/s)
 * @param v : velocity in y-direction(m/s)
 * @param z : DTM = DEM and obstacles
 * @return average dt in flood loop
 */

// * @param q1: flux in the x-direction(m2/s)
// * @param q2: flux in the y-direction(m2/s)
double TWorld::fullSWOF2Do2(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct)//, cTMap *q1, cTMap *q2)
{
    double dt1 = 0, dt2, timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int n = 0;
    double sumh = 0;

    if (prepareFlood)
        prepareFloodZ(z);

    // if there is no flood skip everything
    if (startFlood)
    {
        if (correct)
            sumh = getMass(h);

        verif = 1;

        do {

            if (verif == 1)
            {
                dt1 = dt_max;

                setZero(h, u, v);

                // Reconstruction for order 2
                // makes h1r, h1l, u1r, u1l, v1r, v1l
                // makes h2r, h2l, u2r, u2l, v2r, v2l
                // makes delzc1, delzc2, delz1, delz2
                    MUSCL(h,u,v,z);
            }

            dt1 = maincalcflux(dt1, dt_max);
            dt1 = std::min(dt1, _dt-timesum);
            if(SwitchErosion)
            {
                //temporarily store all the values from the MUSCL or ENO, so the sediment transport model can use these
                //otherwise they will be overwritten by the second reconstruction
                FOR_ROW_COL_MV
                {
                    temp1->Drc = h1d->Drc;
                    temp2->Drc = h1g->Drc;
                    temp3->Drc = h2d->Drc;
                    temp4->Drc = h2g->Drc;
                    temp5->Drc = u1r->Drc;
                    temp6->Drc = u1l->Drc;
                    temp7->Drc = v1r->Drc;
                    temp8->Drc = v1l->Drc;
                    temp9->Drc = u2r->Drc;
                    temp10->Drc = u2l->Drc;
                    temp11->Drc = v2r->Drc;
                    temp12->Drc = v2l->Drc;

                }
            }

            //st venant equations, h, u, v go in hs, vs, us come out
            maincalcscheme(dt1, h,u,v, hs,us,vs);
            dt2 = dt1;

            setZero(hs, us, vs);

            //Reconstruction for order 2
              MUSCL(hs,us,vs,z);

            dt2 = maincalcflux(dt2, dt_max);

            if (dt2 < dt1)
            {
                dt1 = dt2;
                verif = 0;
            }
            else
            {
                verif = 1;
                //hs, us, vs go in hsa, vsa, usa come out
                maincalcscheme(dt1, hs,us,vs, hsa, usa, vsa);
                dt1 = std::min(dt1, _dt-timesum);

                setZero(hsa, usa, vsa);

                //sediment
                SWOFSediment(dt1,h,u,v );

                //Heun method (order 2 in time)
                FOR_ROW_COL_MV
                {
                    double havg = 0.5*(h->Drc + hsa->Drc);
                    if (havg >= he_ca)
                    {
                        double q1 = 0.5*(h->Drc*u->Drc + hsa->Drc*usa->Drc);
                        u->Drc = q1/havg;
                        double q2 = 0.5*(h->Drc*v->Drc + hsa->Drc*vsa->Drc);
                        v->Drc = q2/havg;
                        h->Drc = havg;
                    }
                    else
                    {
                        u->Drc = 0;
                        v->Drc = 0;
                        h->Drc = 0;
                    }
                }//Heun

                timesum = timesum + dt1;
                n++;

                if (n > F_MaxIter)
                    break;
            }//end for else dt2<dt1

            if (correct)
                correctMassBalance(sumh, h, 1e-12);

            // qDebug() << n << timesum << dt1;
          } while (timesum  < _dt);

      } // if floodstart
    iter_n = n;
    dt1 = n > 0? _dt/n : dt1;
    return(dt1);
}
//---------------------------------------------------------------------------
// 2nd order without itertaion dt1, dt2!
double TWorld::fullSWOF2Do2light(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct)//, cTMap *q1, cTMap *q2)
{
    double dt1 = 0, dt2, timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    int n = 0;
    double sumh = 0;

    if (prepareFlood)
        prepareFloodZ(z);

    // if there is no flood skip everything
    if (startFlood)
    {
        if (correct)
            sumh = getMass(h);

        do {

            dt1 = dt_max;

            setZero(h, u, v);

            // Reconstruction for order 2
            // makes h1r, h1l, u1r, u1l, v1r, v1l
            // makes h2r, h2l, u2r, u2l, v2r, v2l
            // makes delzc1, delzc2, delz1, delz2
            MUSCL(h,u,v,z);

            dt1 = maincalcflux(dt1, dt_max);
            dt1 = std::min(dt1, _dt-timesum);

            if(SwitchErosion)
            {
                //temporarily store all the values from the MUSCL or ENO, so the sediment transport model can use these
                //otherwise they will be overwritten by the second reconstruction
                FOR_ROW_COL_MV
                {
                    temp1->Drc = h1d->Drc;
                    temp2->Drc = h1g->Drc;
                    temp3->Drc = h2d->Drc;
                    temp4->Drc = h2g->Drc;
                    temp5->Drc = u1r->Drc;
                    temp6->Drc = u1l->Drc;
                    temp7->Drc = v1r->Drc;
                    temp8->Drc = v1l->Drc;
                    temp9->Drc = u2r->Drc;
                    temp10->Drc = u2l->Drc;
                    temp11->Drc = v2r->Drc;
                    temp12->Drc = v2l->Drc;

                }
            }

            //st venant equations, h, u, v go in hs, vs, us come out
            maincalcscheme(dt1, h,u,v, hs,us,vs);

            setZero(hs, us, vs);

            //Reconstruction for order 2
            MUSCL(hs,us,vs,z);
            dt2 = maincalcflux(dt1, dt_max);

            //hs, us, vs go in hsa, vsa, usa come out
            maincalcscheme(dt1, hs,us,vs, hsa, usa, vsa);
            dt1 = std::min(dt1, _dt-timesum);

            setZero(hsa, usa, vsa);

            //sediment
            SWOFSediment(dt1,h,u,v );

            //Heun method (order 2 in time)
            FOR_ROW_COL_MV
            {
                double havg = 0.5*(h->Drc + hsa->Drc);
                if (havg >= he_ca)
                {
                    double q1 = 0.5*(h->Drc*u->Drc + hsa->Drc*usa->Drc);
                    u->Drc = q1/havg;
                    double q2 = 0.5*(h->Drc*v->Drc + hsa->Drc*vsa->Drc);
                    v->Drc = q2/havg;
                    h->Drc = havg;
                }
                else
                {
                    u->Drc = 0;
                    v->Drc = 0;
                    h->Drc = 0;
                }
            }//Heun

            timesum = timesum + dt1;
            n++;

            if (n > F_MaxIter)
                break;
            if (correct)
                correctMassBalance(sumh, h, 1e-12);

            // qDebug() << n << timesum << dt1;
        } while (timesum  < _dt);

    } // if floodstart
    iter_n = n;
    dt1 = n > 0? _dt/n : dt1;
    return(dt1);
}
