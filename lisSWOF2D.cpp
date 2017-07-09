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

#define he_ca 1e-4
#define ve_ca 1e-6

#define dt_ca 0.001

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
/// MUSCL: Monotone Upstream-centered Schemes for Conservation Laws
/// see http://en.wikipedia.org/wiki/MUSCL_scheme
///
/// This MUSCL creates left and right u,v,h, arrays for in mainfluxcalc
void TWorld::MUSCL(cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z)
{
  double delta_h1, delta_u1, delta_v1;
  double delta_h2, delta_u2, delta_v2;
  double dh, du, dv, dz_h;

  simpleScheme(_h, _u, _v);

  FOR_CELL_IN_FLOODAREA
  {
      if(c > 0 && !MV(r,c-1)) {
          delta_h1 = _h->Drc - _h->data[r][c-1];
          delta_u1 = _u->Drc - _u->data[r][c-1];
          delta_v1 = _v->Drc - _v->data[r][c-1];
      } else {
          delta_h1 = 0;
          delta_u1 = 0;
          delta_v1 = 0;
      }
      if(c < _nrCols-1 &&  !MV(r,c+1)) {
          delta_h2 = _h->data[r][c+1] - _h->Drc;
          delta_u2 = _u->data[r][c+1] - _u->Drc;
          delta_v2 = _v->data[r][c+1] - _v->Drc;
      } else {
          delta_h2 = 0;
          delta_u2 = 0;
          delta_v2 = 0;
      }

      dh   = 0.5*limiter(delta_h1, delta_h2);
      dz_h = 0.5*limiter(delta_h1 + delta_z1->data[r][c-1], delta_h2 + delta_z1->Drc);
      du   = 0.5*limiter(delta_u1, delta_u2);
      dv   = 0.5*limiter(delta_v1, delta_v2);

      h1r->Drc = _h->Drc+dh;
      h1l->Drc = _h->Drc-dh;

      z1r->Drc = _z->Drc+(dz_h-dh);
      z1l->Drc = _z->Drc+(dh-dz_h);

      delzc1->Drc = (long double)z1r->Drc-(long double)z1l->Drc;
      delz1->data[r][c-1] = z1l->Drc - z1r->data[r][c-1];  //outside noundary a value exists with preparefloodz

      if (_h->Drc > he_ca) {
          double h1lh = h1l->Drc/_h->Drc;
          double h1rh = h1r->Drc/_h->Drc;

          u1r->Drc = _u->Drc + h1lh * du;
          u1l->Drc = _u->Drc - h1rh * du;
          v1r->Drc = _v->Drc + h1lh * dv;
          v1l->Drc = _v->Drc - h1rh * dv;
      } else {
          u1r->Drc = _u->Drc + du;
          u1l->Drc = _u->Drc - du;
          v1r->Drc = _v->Drc + dv;
          v1l->Drc = _v->Drc - dv;
      }
  }}

    FOR_CELL_IN_FLOODAREA {
        if(r > 0 && MV(r-1,c)) {
            delta_h1 = _h->Drc - _h->data[r-1][c];
            delta_u1 = _u->Drc - _u->data[r-1][c];
            delta_v1 = _v->Drc - _v->data[r-1][c];
        } else {
            delta_h1 = 0;
            delta_u1 = 0;
            delta_v1 = 0;
        }
        if(r < _nrRows-1 && !MV(r+1, c)) {
            delta_h2 = _h->data[r+1][c] - _h->Drc;
            delta_u2 = _u->data[r+1][c] - _u->Drc;
            delta_v2 = _v->data[r+1][c] - _v->Drc;
        } else {
            delta_h2 = 0;
            delta_u2 = 0;
            delta_v2 = 0;
        }

        dh   = 0.5*limiter(delta_h1, delta_h2);
        dz_h = 0.5*limiter(delta_h1+delta_z2->data[r-1][c],delta_h2+delta_z2->Drc);
        du   = 0.5*limiter(delta_u1, delta_u2);
        dv   = 0.5*limiter(delta_v1, delta_v2);

        h2r->Drc = _h->Drc+dh;
        h2l->Drc = _h->Drc-dh;

        z2r->Drc = _z->Drc+(dz_h-dh);
        z2l->Drc = _z->Drc+(dh-dz_h);

        delzc2->Drc = (long double)z2r->Drc - (long double)z2l->Drc;
        delz2->data[r-1][c] = z2l->Drc - z2r->data[r-1][c]; //outside noundary a value exists with preparefloodz

        if (_h->Drc > he_ca) {
            double h2lh = h2l->Drc/_h->Drc;
            double h2rh = h2r->Drc/_h->Drc;

            u2r->Drc = _u->Drc + h2lh * du;
            u2l->Drc = _u->Drc - h2rh * du;
            v2r->Drc = _v->Drc + h2lh * dv;
            v2l->Drc = _v->Drc - h2rh * dv;
        } else {
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
    FOR_CELL_IN_FLOODAREA {
      double tx = dt/ChannelAdj->Drc;
      double ty = dt/DX->Drc;

      hes->Drc = he->Drc;
      if (c < _nrCols-1 && !MV(c+1,r))
          hes->Drc -= tx*(f1->data[r][c+1]-f1->Drc);
      if (r < _nrRows-1 && !MV(c,r+1))
          hes->Drc -= ty*(g1->data[r+1][c]-g1->Drc);

      if (hes->Drc > he_ca)
      {
          //Solution of the equation of momentum (Second and third equation of Saint-venant)
          double qes1;
          double qes2;

          // fullswof version 1.04
          // Solution of the equation of momentum (Second and third equation of Shallow-Water)
          // This expression for the flux (instead of the differences of the squares) avoids numerical errors
          // see http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html section "Cancellation".
          int rr = r < _nrRows-1 && !MV(r+1,c) ? r+1 : r;
          int cc = c < _nrCols-1 && !MV(r,c+1) ? c+1 : c;
          qes1 = he->Drc*ve1->Drc -
                  ty*(g2->data[rr][c]-g2->Drc) -
                  tx*(f2->data[r][cc]-f2->Drc +
                      GRAV*0.5*((h1g->Drc-h1l->Drc)*(h1g->Drc+h1l->Drc) +
                                (h1r->Drc-h1d->Drc)*(h1r->Drc+h1d->Drc) +
                                (h1l->Drc+h1r->Drc)*delzc1->Drc)) ;

          // fullswof version 1.04
          qes2 = he->Drc*ve2->Drc - tx*(f3->data[r][cc]-f3->Drc) -
                  ty*(g3->data[rr][c]-g3->Drc +
                      GRAV*0.5*((h2g->Drc-h2l->Drc)*(h2g->Drc+h2l->Drc) +
                                (h2r->Drc-h2d->Drc)*(h2r->Drc+h2d->Drc) +
                                (h2l->Drc+h2r->Drc)*delzc2->Drc));


          //Calc friction semi-implicit with old v and u and new h
          //Fr_Manning(ve1->Drc, ve2->Drc, hes->Drc, qes1, qes2, dt, N->Drc);

          double nsq = N->Drc*N->Drc*GRAV*sqrt(ve1->Drc*ve1->Drc+ve2->Drc*ve2->Drc)*dt/qPow(hes->Drc,4.0/3.0);

          ves1->Drc = (qes1/(1.0+nsq))/hes->Drc;
          ves2->Drc = (qes2/(1.0+nsq))/hes->Drc;

          if (std::fabs(ves1->Drc) > 25 || std::fabs(ves2->Drc) > 25) {
              ves1->Drc = hes->Drc/_dt;
              ves2->Drc = hes->Drc/_dt;
          }
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

    FOR_CELL_IN_FLOODAREA {
        int cc = (c > 0 && !MV(r,c-1)) ? c-1 : c;
        h1d->data[r][cc] = std::max(0.0, h1r->data[r][cc] - std::max(0.0,  delz1->data[r][cc] + std::max(FlowBarrierW->Drc,FlowBarrierE->data[r][cc])));
        h1g->Drc         = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][cc] + std::max(FlowBarrierW->Drc,FlowBarrierE->data[r][cc])));
        if (F_scheme == 1)
            F_Rusanov(h1d->data[r][cc], u1r->data[r][cc], v1r->data[r][cc],h1g->Drc, u1l->Drc, v1l->Drc);
        else
            if (F_scheme == 2)
                F_HLL(h1d->data[r][cc], u1r->data[r][cc], v1r->data[r][cc],h1g->Drc, u1l->Drc, v1l->Drc);
            else
                F_HLL2(h1d->data[r][cc], u1r->data[r][cc], v1r->data[r][cc],h1g->Drc, u1l->Drc, v1l->Drc);

        f1->Drc = HLL2_f1;//(std::abs(HLL2_f1) < he_ca ? 0: HLL2_f1);
        f2->Drc = HLL2_f2;//(std::abs(HLL2_f2) < he_ca ? 0: HLL2_f2);
        f3->Drc = HLL2_f3;//(std::abs(HLL2_f3) < he_ca ? 0: HLL2_f3);

        cflx->Drc = HLL2_cfl;
    }}

    FOR_CELL_IN_FLOODAREA {
        int rr = (r > 0 && !MV(r-1,c)) ? r-1 : r;
        h2d->data[rr][c] = std::max(0.0, h2r->data[rr][c] - std::max(0.0,  delz2->data[rr][c] + std::max(FlowBarrierS->Drc,FlowBarrierN->data[rr][c])));
        h2g->Drc         = std::max(0.0, h2l->Drc          - std::max(0.0, -delz2->data[rr][c] + std::max(FlowBarrierS->Drc,FlowBarrierN->data[rr][c])));
        if (F_scheme == 1)
            F_Rusanov(h2d->data[rr][c],v2r->data[rr][c],u2r->data[rr][c],h2g->Drc,v2l->Drc,u2l->Drc);
        else
            if (F_scheme == 2)
                F_HLL(h2d->data[rr][c],v2r->data[rr][c],u2r->data[rr][c],h2g->Drc,v2l->Drc,u2l->Drc);
            else
                F_HLL2(h2d->data[rr][c],v2r->data[rr][c],u2r->data[rr][c],h2g->Drc,v2l->Drc,u2l->Drc);

        g1->Drc = HLL2_f1;//(std::abs(HLL2_f1) < he_ca ? 0: HLL2_f1);
        g2->Drc = HLL2_f3;//(std::abs(HLL2_f3) < he_ca ? 0: HLL2_f3);
        g3->Drc = HLL2_f2;//(std::abs(HLL2_f2) < he_ca ? 0: HLL2_f2);
        cfly->Drc = HLL2_cfl;
    }}


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

     return(std::max(TimestepfloodMin, std::min(dtx,dty)));

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
            if(!pcr::isMV(LDD->data[r][c]) && !pcr::isMV(LDD->data[r][c-1]))
            {
                delz1->data[r][c-1] = z->Drc - z->data[r][c-1];
                // needed in maincalcflux for 1D scheme, is calculated in MUSCL for 2D scheme
            }
    for (int r = 1; r < _nrRows; r++)
        for (int c = 0; c < _nrCols; c++)
            if(!pcr::isMV(LDD->data[r][c]) ||
                    !pcr::isMV(LDD->data[r-1][c]))
            {
                delz2->data[r-1][c] = z->Drc - z->data[r-1][c];
                // needed in maincalcflux for 1D scheme, is calculated in MUSCL for 2D scheme
            }
    fill(*tm,0);
    FOR_ROW_COL_MV
    {
        double v = 0;
        double cnt = 0;
        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
        {
            if(INSIDE(r+i, c+j) && delz1->data[r+i][c+j] > -9999 && !pcr::isMV(LDD->data[r+i][c+j]))
            {
                v += delz1->data[r+i][c+j]; cnt+=1.0;
            }
            tm->Drc = cnt > 0 ? v/cnt : 0;
        }
    }
    FOR_ROW_COL_MV
        if(delz1->Drc == -9999) delz1->Drc = tm->Drc;

    fill(*tm,0);
    FOR_ROW_COL_MV
    {
        double v = 0;
        double cnt = 0;
        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
        {
            if(INSIDE(r+i, c+j) && delz2->data[r+i][c+j] > -9999 && !pcr::isMV(LDD->data[r+i][c+j]))
            {
                v += delz2->data[r+i][c+j]; cnt+=1.0;
            }
            tm->Drc = cnt > 0 ? v/cnt : 0;
        }
    }
    FOR_ROW_COL_MV
        if(delz2->Drc == -9999) delz2->Drc = tm->Drc;

    fill(*delta_z1, -9999);
    fill(*delta_z2, -9999);
    for (int r = 0; r < _nrRows; r++)
        for (int c = 0; c < _nrCols-1; c++)
            if(!pcr::isMV(LDD->data[r][c]) ||
                    !pcr::isMV(LDD->data[r][c+1]))
                delta_z1->Drc = z->data[r][c+1] - z->Drc;

    for (int r = 0; r < _nrRows-1; r++)
        for (int c = 0; c < _nrCols; c++)
            if(!pcr::isMV(LDD->data[r][c]) ||
                    !pcr::isMV(LDD->data[r+1][c]))
                delta_z2->Drc = z->data[r+1][c] - z->Drc;

    fill(*tm,0);
    FOR_ROW_COL_MV
    {
        double v = 0;
        double cnt = 0;
        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
        {
            if(INSIDE(r+i, c+j) && delta_z1->data[r+i][c+j] > -9999 && !pcr::isMV(LDD->data[r+i][c+j]))
            {
                v += delta_z1->data[r+i][c+j]; cnt+=1.0;
            }
            tm->Drc = cnt > 0 ? v/cnt : 0;
        }
    }
    FOR_ROW_COL_MV
        if(delta_z1->Drc == -9999) delta_z1->Drc = tm->Drc;

    fill(*tm,0);
    FOR_ROW_COL_MV
    {
        double v = 0;
        double cnt = 0;
        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
        {
            if(INSIDE(r+i, c+j) && delta_z2->data[r+i][c+j] > -9999 && !pcr::isMV(LDD->data[r+i][c+j]))
            {
                v += delta_z2->data[r+i][c+j]; cnt+=1.0;
            }
            tm->Drc = cnt > 0 ? v/cnt : 0;
        }
    }
    FOR_ROW_COL_MV
        if(delta_z2->Drc == -9999) delta_z2->Drc = tm->Drc;

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
          //MUSCL(h,u,v,z);

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
            simpleScheme(hs, us, vs);
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
    double dt1 = 0, dt2 = 0, timesum = 0;
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
           // simpleScheme(h, u, v);
            MUSCL(h,u,v,z);

            dt1 = maincalcflux(dt1, dt_max);
            dt1 = std::min(dt1, _dt-timesum);

            //st venant equations, h, u, v go in hs, vs, us come out
            maincalcscheme(dt1, h,u,v, hs,us,vs);

            setZero(hs, us, vs);
            setZero(hs, us, vs);

            FOR_CELL_IN_FLOODAREA {
              h->Drc = hs->Drc;
              u->Drc = us->Drc;
              v->Drc = vs->Drc;
            }}

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
