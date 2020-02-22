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
-   void F_HLL2(double hg,double ug,double vg,double hd,double ud,double vd);
-   void F_HLL(double hg,double ug,double vg,double hd,double ud,double vd);
-   void F_Rusanov(double hg,double ug,double vg,double hd,double ud,double vd);
-   void setZero(cTMap *h, cTMap *u, cTMap *v);//, cTMap *q1, cTMap *q2);
-   double limiter(double a, double b);
*/

#include <algorithm>
#include "lisemqt.h"
//#include "model.h"
#include "operation.h"
#include "global.h"

#define he_ca 1e-6
#define ve_ca 1e-6

#define dt_ca 0.001

#define GRAV 9.8067
#define EPSILON 1e-6

#define dtmaxfrac 0.5

//--------------------------------------------------------------------------------------------
// correct mass balance
double TWorld::getMass(cTMap *M)
{
    double sum2 = 0;
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0)
            sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
    }
    return sum2;
}
//---------------------------------------------------------------------------
// correct mass balance
void TWorld::correctMassBalance(double sum1, cTMap *M)
{
    double sum2 = 0;
    double n = 0;
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0)
        {
            sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
            if(M->Drc > 0)
                n += 1;
        }
    }
    // total and cells active for M

    //double dh = (n > 0 ? (sum1 - sum2)/n : 0);
    double dhtot = sum2 > 0 ? (sum1 - sum2)/sum2 : 0;
    FOR_ROW_COL_MV
    {
        if(M->Drc > 0)
        {
            M->Drc = M->Drc*(1.0 + dhtot);            // <- distribution weighted to h
            //M->Drc += dh/(DX->Drc*ChannelAdj->Drc); // <- equal distribution error
            M->Drc = std::max(M->Drc , 0.0);
        }
    }
}
//---------------------------------------------------------------------------
// used in datainit, done once
void TWorld::prepareFloodZ(cTMap *z)
{

    prepareFlood = false;

    fill(*delz1,0);
    fill(*delz2,0);
    // diff between z cell and adjacent
    for (int r = 0; r < _nrRows; r++)
        for (int c = 1; c < _nrCols; c++)
            if(!pcr::isMV(LDD->data[r][c]) && !pcr::isMV(LDD->data[r][c-1]))
            {
                delz1->data[r][c-1] = (z->Drc) - z->data[r][c-1];
                // needed in maincalcflux for 1D scheme, is calculated in MUSCL for 2D scheme
            }
    for (int r = 1; r < _nrRows; r++)
        for (int c = 0; c < _nrCols; c++)
            if(!pcr::isMV(LDD->data[r][c]) && !pcr::isMV(LDD->data[r-1][c]))
            {
                delz2->data[r-1][c] = z->Drc - z->data[r-1][c];
                // needed in maincalcflux for 1D scheme, is calculated in MUSCL for 2D scheme
            }

    fill(*delta_z1, 0);
    fill(*delta_z2, 0);
    for (int r = 0; r < _nrRows; r++)
        for (int c = 0; c < _nrCols-1; c++)
            if(!pcr::isMV(LDD->data[r][c]) || !pcr::isMV(LDD->data[r][c+1]))
                delta_z1->Drc = (z->data[r][c+1] - z->Drc);
    //                delta_z1->Drc = 0.5*((z->Drc - z->data[r][c-1]) + (z->data[r][c+1] - z->Drc));

    for (int r = 1; r < _nrRows-1; r++)
        for (int c = 0; c < _nrCols; c++)
            if(!pcr::isMV(LDD->data[r][c]) || !pcr::isMV(LDD->data[r+1][c]))
                delta_z2->Drc = (z->data[r+1][c] - z->Drc);
    //                delta_z2->Drc = 0.5*((z->Drc - z->data[r-1][c]) + (z->data[r+1][c] - z->Drc));

}
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
    // F_fluxLimiter=1;
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
void TWorld::setZero(int thread, cTMap *_h, cTMap *_u, cTMap *_v)
{
  FOR_ROW_COL_UF2DMTDER
  //FOR_ROW_COL_UF2DMT_DT
  {
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
  }}}}
}
//---------------------------------------------------------------------------
/// Numerical flux calculation on which the new velocity is based
/// U_n+1 = U_n + dt/dx* [flux]  when flux is calculated by HLL, HLL2, Rusanov
/// HLL = Harten, Lax, van Leer numerical solution
void TWorld::F_HLL2(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R)
{
    double f1, f2, f3, cfl, tmp = 0;
    if (h_L<=0. && h_R<=0.){
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
    if (h_L<=0. && h_R<=0.){
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
void TWorld::MUSCL(int thread, cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z)
{
    double delta_h1, delta_u1, delta_v1, dz1;
    double delta_h2, delta_u2, delta_v2, dz2;
    double dh, du, dv, dz_h, hlh, hrh;

    FOR_ROW_COL_UF2DMT_DT {
        if(FloodHMaskDer->Drc != 0) {

            if(c > 0 && !MV(r,c-1) && c < _nrCols-1 && !MV(r,c+1)) {
                delta_h1 = _h->Drc - _h->data[r][c-1];
                delta_u1 = _u->Drc - _u->data[r][c-1];
                delta_v1 = _v->Drc - _v->data[r][c-1];
                dz1 = _z->Drc - _z->data[r][c-1];
                delta_h2 = _h->data[r][c+1] - _h->Drc;
                delta_u2 = _u->data[r][c+1] - _u->Drc;
                delta_v2 = _v->data[r][c+1] - _v->Drc;
                dz2 = _z->data[r][c+1] - _z->Drc;
            } else {
                delta_h1 = 0;
                delta_u1 = 0;
                delta_v1 = 0;
                delta_h2 = 0;
                delta_u2 = 0;
                delta_v2 = 0;
                dz1 = 0;
                dz2 = 0;
            }

            dh   = 0.5*limiter(delta_h1, delta_h2);
            dz_h = 0.5*limiter(delta_h1 + dz1, delta_h2 + dz2);

            du   = 0.5*limiter(delta_u1, delta_u2);
            dv   = 0.5*limiter(delta_v1, delta_v2);

            h1r->Drc = _h->Drc+dh;
            h1l->Drc = _h->Drc-dh;

            z1r->Drc = _z->Drc+(dz_h-dh);
            z1l->Drc = _z->Drc+(dh-dz_h);

            delzc1->Drc = z1r->Drc-z1l->Drc;
            if (c > 0 && !MV(r,c-1))
                delz1->data[r][c-1] = z1l->Drc - z1r->data[r][c-1];

            if (_h->Drc > he_ca) {
                hlh = h1l->Drc/_h->Drc;
                hrh = h1r->Drc/_h->Drc;
            }
            else {
                hlh = 1.0;
                hrh = 1.0;
            }

            u1r->Drc = _u->Drc + hlh * du;
            u1l->Drc = _u->Drc - hrh * du;
            v1r->Drc = _v->Drc + hlh * dv;
            v1l->Drc = _v->Drc - hrh * dv;
        }
    }}}}

    FOR_ROW_COL_UF2DMT_DT {
        if(FloodHMaskDer->Drc != 0)
        {
            if(r > 0 && MV(r-1,c) && r < _nrRows-1 && !MV(r+1, c)) {
                delta_h1 = _h->Drc - _h->data[r-1][c];
                delta_u1 = _u->Drc - _u->data[r-1][c];
                delta_v1 = _v->Drc - _v->data[r-1][c];
                dz1 = _z->Drc - _z->data[r-1][c];
                delta_h2 = _h->data[r+1][c] - _h->Drc;
                delta_u2 = _u->data[r+1][c] - _u->Drc;
                delta_v2 = _v->data[r+1][c] - _v->Drc;
                dz2 = _z->data[r+1][c] - _z->Drc;
            } else {
                delta_h1 = 0;
                delta_u1 = 0;
                delta_v1 = 0;
                delta_h2 = 0;
                delta_u2 = 0;
                delta_v2 = 0;
                dz1 = 0;
                dz2 = 0;
            }
            dh   = 0.5*limiter(delta_h1, delta_h2);
            dz_h = 0.5*limiter(delta_h1+dz2,delta_h2+dz2);
            du   = 0.5*limiter(delta_u1, delta_u2);
            dv   = 0.5*limiter(delta_v1, delta_v2);

            h2r->Drc = _h->Drc+dh;
            h2l->Drc = _h->Drc-dh;

            z2r->Drc = _z->Drc+(dz_h-dh);
            z2l->Drc = _z->Drc+(dh-dz_h);

            delzc2->Drc = z2r->Drc - z2l->Drc;
            if(r > 0 && MV(r-1,c))
                delz2->data[r-1][c] = z2l->Drc - z2r->data[r-1][c];

            if (_h->Drc > he_ca) {
                hlh = h2l->Drc/_h->Drc;
                hrh = h2r->Drc/_h->Drc;
            }
            else {
                hlh = 1.0;
                hrh = 1.0;
            }

            u2r->Drc = _u->Drc + hlh * du;
            u2l->Drc = _u->Drc - hrh * du;
            v2r->Drc = _v->Drc + hlh * dv;
            v2l->Drc = _v->Drc - hrh * dv;
        }
    }}}}
}
//---------------------------------------------------------------------------

    /**
     Construction of variables for hydrostatic reconstruction using HLL or Rusanov
     Flux in x and y direction. The method is based solving PDEs with an estimate of the status
     of a domain based on spatial averages of the previous timestep, this is the hydrostatic equlibrium
     Calculaton of the time steps in relation to cfl.
    */

void TWorld::maincalcflux(int thread, double dt, double dt_max)
{
    double dtx, dty;

    cTMap *fbw = FlowBarrierW;
    cTMap *fbe = FlowBarrierE;
    cTMap *fbn = FlowBarrierN;
    cTMap *fbs = FlowBarrierS;

    FOR_ROW_COL_UF2DMT_DT {
        if(FloodHMaskDer->Drc != 0) {
            f1->Drc = 0;
            f2->Drc = 0;
            f3->Drc = 0;
            f1o->Drc = 0;
            f2o->Drc = 0;
            f3o->Drc = 0;
//            double bh = 0.1;
//            fbn->Drc = 0;
//            fbs->Drc = 0;
//            fbe->Drc = 0;
//            fbw->Drc = 0;
//            if(LDD->Drc == 1) { fbn->Drc = bh; fbw->Drc = bh;}
//            if(LDD->Drc == 2) { fbe->Drc = bh; fbw->Drc = bh;}
//            if(LDD->Drc == 3) { fbn->Drc = bh; fbe->Drc = bh;}
//            if(LDD->Drc == 4) { fbn->Drc = bh; fbs->Drc = bh;}
//            if(LDD->Drc == 6) { fbn->Drc = bh; fbs->Drc = bh;}
//            if(LDD->Drc == 7) { fbs->Drc = bh; fbw->Drc = bh;}
//            if(LDD->Drc == 8) { fbe->Drc = bh; fbw->Drc = bh;}
//            if(LDD->Drc == 9) { fbs->Drc = bh; fbe->Drc = bh;}

            if(c > 0 && !MV(r,c-1)) {
                h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
                h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
                if (F_scheme == 1)
                    F_Rusanov(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
                else
                    if (F_scheme == 2)
                        F_HLL(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
                    else
                        F_HLL2(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);

                f1->Drc = HLL2_f1;
                f2->Drc = HLL2_f2;
                f3->Drc = HLL2_f3;
                cflx->Drc = HLL2_cfl;
            } else {
                double _h1g = std::max(0.0, h1l->Drc - FlowBarrierE->Drc);
                if (F_scheme == 1)
                    F_Rusanov(0,0,0, _h1g, u1l->Drc, v1l->Drc);
                else
                    if (F_scheme == 2)
                        F_HLL(0,0,0, _h1g, u1l->Drc, v1l->Drc);
                    else
                        F_HLL2(0,0,0, _h1g, u1l->Drc, v1l->Drc);
                f1->Drc = HLL2_f1;
                f2->Drc = HLL2_f2;
                f3->Drc = HLL2_f3;
                cflx->Drc = HLL2_cfl;
            }
            // right hand side boundary
            if(c == _nrCols-1 || MV(r, c+1)){
                double _h1d = std::max(0.0, h1r->Drc - fbw->Drc);
                if (F_scheme == 1)
                    F_Rusanov(_h1d,u1r->Drc,v1r->Drc,0.,0.,0.);
                else
                    if (F_scheme == 2)
                        F_HLL(_h1d,u1r->Drc,v1r->Drc,0.,0.,0.);
                    else
                        F_HLL2(_h1d,u1r->Drc,v1r->Drc,0.,0.,0.);
                f1o->Drc = HLL2_f1;
                f2o->Drc = HLL2_f2;
                f3o->Drc = HLL2_f3;
            }

        }
    }}}}

    FOR_ROW_COL_UF2DMT_DT {
        if(FloodHMaskDer->Drc != 0){
            g1->Drc = 0;
            g2->Drc = 0;
            g3->Drc = 0;
            g1o->Drc = 0;
            g2o->Drc = 0;
            g3o->Drc = 0;

            if(r > 0 && !MV(r-1,c)) {
                h2d->data[r-1][c] = std::max(0.0, h2r->data[r-1][c] - std::max(0.0,  delz2->data[r-1][c]  + std::max(fbs->Drc,fbn->data[r-1][c])));
                h2g->Drc          = std::max(0.0, h2l->Drc          - std::max(0.0, -delz2->data[r-1][c]  + std::max(fbs->Drc,fbn->data[r-1][c])));

                if (F_scheme == 1)
                    F_Rusanov(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c], h2g->Drc,v2l->Drc,u2l->Drc);
                else
                    if (F_scheme == 2)
                        F_HLL(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c], h2g->Drc,v2l->Drc,u2l->Drc);
                    else
                        F_HLL2(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c], h2g->Drc,v2l->Drc,u2l->Drc);

                g1->Drc = HLL2_f1;
                g2->Drc = HLL2_f3;
                g3->Drc = HLL2_f2;
                cfly->Drc = HLL2_cfl;
            } else {
                double _h2g = std::max(0.0, h2l->Drc - fbn->Drc);
                if (F_scheme == 1)
                    F_Rusanov(0,0,0,_h2g,v2l->Drc,u2l->Drc);
                else
                    if (F_scheme == 2)
                        F_HLL(0,0,0,_h2g,v2l->Drc,u2l->Drc);
                    else
                        F_HLL2(0,0,0,_h2g,v2l->Drc,u2l->Drc);
                g1->Drc = HLL2_f1;
                g2->Drc = HLL2_f2;
                g3->Drc = HLL2_f3;
                cflx->Drc = HLL2_cfl;
            }
            // left hand side boundary
            if (r == _nrRows-1 || MV(r+1, c)) {
                double _h2d = std::max(0.0, h2d->Drc - fbs->Drc);
                if (F_scheme == 1)
                    F_Rusanov(_h2d,v2l->Drc,u2l->Drc,0.,0.,0.);
                else
                    if (F_scheme == 2)
                        F_HLL(_h2d,v2l->Drc,u2l->Drc,0.,0.,0.);
                    else
                        F_HLL2(_h2d,v2l->Drc,u2l->Drc,0.,0.,0.);
                g1o->Drc = HLL2_f1;
                g2o->Drc = HLL2_f3;
                g3o->Drc = HLL2_f2;
            }
        }
    }}}}

    // find largest velocity and determine dt
    dt = dt_max;
    FOR_ROW_COL_UF2DMT_DT {
        dtx = dt_max;
        dty = dt_max;
        if(FloodHMaskDer->Drc != 0){
            double dx = FlowWidth->Drc;// ChannelAdj->Drc;
            if (qFabs(cflx->Drc*dt/dx) > 1e-10)
                dtx = std::min(dt_max, courant_factor*dx/cflx->Drc);

            double dy = DX->Drc;
            if (qFabs(cfly->Drc*dt/dy) > 1e-10)
                dty = std::min(dt_max, courant_factor*dy/cfly->Drc);

            FloodDT->Drc = std::min(dtx, dty);

        }
    }}}}

}
//---------------------------------------------------------------------------
// St Venant equations: conservation of mass and momentum: h u et v (= he, ve1, ve2)
void TWorld::maincalcscheme(int thread, double dt, cTMap *he, cTMap *ve1, cTMap *ve2,
                            cTMap *hes, cTMap *ves1, cTMap *ves2)
{
    FOR_ROW_COL_UF2DMT_DT {
        if(FloodHMaskDer->Drc != 0)
        {
            if (SwitchVariableTimestep)
                dt = FloodDT->Drc;
            else
                dt = Flood_DTMIN;
            double tx = dt/FlowWidth->Drc;//ChannelAdj->Drc; //
            double ty = dt/DX->Drc;
            double _f1=0, _f2=0, _f3=0, _g1=0, _g2=0, _g3=0;

            //choose left hand boundary and normal (f1), or right hand boundary values (f1o)
            if (c < _nrCols-1 && !MV(r, c+1)){
                _f1 = f1->data[r][c+1];
                _f2 = f2->data[r][c+1];
                _f3 = f3->data[r][c+1];
            }
            if (c == _nrCols-1 || MV(r, c+1)){
                _f1 = f1o->Drc;
                _f2 = f2o->Drc;
                _f3 = f3o->Drc;
            }
            if (r < _nrRows-1 && !MV(r+1, c)) {
                _g1 = g1->data[r+1][c];
                _g2 = g2->data[r+1][c];
                _g3 = g3->data[r+1][c];
            }
            else if (r == _nrRows-1 || MV(r+1, c)) {
                _g1 = g1o->Drc;
                _g2 = g2o->Drc;
                _g3 = g3o->Drc;
            }
            hes->Drc = std::max(0.0, he->Drc - tx*(_f1 - f1->Drc) - ty*(_g1 - g1->Drc));

            if (hes->Drc > he_ca)
            {
                //Solution of the equation of momentum (Second and third equation of Saint-venant)
                double qes1;
                double qes2;

                qes1 = he->Drc*ve1->Drc -
                        ty*(_g2 - g2->Drc) -
                        tx*(_f2 - f2->Drc +
                            GRAV*0.5*((h1g->Drc-h1l->Drc)*(h1g->Drc+h1l->Drc) +
                                      (h1r->Drc-h1d->Drc)*(h1r->Drc+h1d->Drc)
                                      + (h1l->Drc+h1r->Drc)*delzc1->Drc)) ;

                qes2 = he->Drc*ve2->Drc -
                        tx*(_f3 - f3->Drc) -
                        ty*(_g3 - g3->Drc +
                            GRAV*0.5*((h2g->Drc-h2l->Drc)*(h2g->Drc+h2l->Drc) +
                                      (h2r->Drc-h2d->Drc)*(h2r->Drc+h2d->Drc)
                                      + (h2l->Drc+h2r->Drc)*delzc2->Drc));


                double sqQ = qSqrt(qes1*qes1+qes2*qes2);

                double sqUV = qSqrt(ve1->Drc*ve1->Drc+ve2->Drc*ve2->Drc);
                double nsq1 = (0.001+N->Drc)*(0.001+N->Drc)*GRAV/qPow(hes->Drc,4.0/3.0);
                double nsq = nsq1*sqUV*dt;

                ves1->Drc = (qes1/(1.0+nsq))/hes->Drc;
                ves2->Drc = (qes2/(1.0+nsq))/hes->Drc;

                double fac = 0;
                if (SwitchTimeavgV) {
                    fac = 0.5+0.5*std::min(1.0,4*hes->Drc)*std::min(1.0,4*hes->Drc);
                    fac = fac *exp(- std::max(1.0,dt) / nsq1);
                }
                ves1->Drc = fac * ve1->Drc + (1.0-fac) *ves1->Drc;
                ves2->Drc = fac * ve2->Drc + (1.0-fac) *ves2->Drc;

                double thv = 10.0;
                double dv = 5.0;
                correctSpuriousVelocities(r, c, hes, ves1, ves2,thv, dv, dt);

       //         sqUV = qSqrt(ves1->Drc*ves1->Drc+ves2->Drc*ves2->Drc);
       //         double frac = std::abs(ves1->Drc/ves2->Drc);

//                if (sqUV*dt*hes->Drc*ChannelAdj->Drc > hes->Drc*ChannelAdj->Drc*_dx) {
//                   frac = _dx/(sqUV*dt);
//                   ves1->Drc *= frac;
//                   ves2->Drc *= frac;
//                   qDebug() << frac;
//                }

            }
            else
            {
                // Case of height of water < ha.
                ves1->Drc = 0;
                ves2->Drc = 0;
            }
        }
        // dan maar even met geweld!
        if (std::isnan(ves1->Drc) || std::isnan(ves2->Drc)  )
        {
            ves1->Drc = 0;
            ves2->Drc = 0;
            hes->Drc = 0;
        }
    }}}}
}
//---------------------------------------------------------------------------
void TWorld::correctSpuriousVelocities(int r, int c, cTMap *hes, cTMap *ves1, cTMap *ves2, double thv, double dv, double dt)
{
    double vs1 = fabs(ves1->Drc);
    double vs2 = fabs(ves2->Drc);
    double Vsqrt = sqrt(vs1*vs1+vs2*vs2) ;
    double sign1 = ves1->Drc < 0? -1.0 : 1.0;
    double sign2 = ves2->Drc < 0? -1.0 : 1.0;

    if (Vsqrt < thv)
        return;

    double vu = r > 0 && !MV(r-1,c) ? fabs(ves1->data[r-1][c])+dv : vs1;
    double vd = r < _nrRows-1 && !MV(r+1,c) ? fabs(ves1->data[r+1][c])+dv : vs1;
    double vl = c > 0 && !MV(r,c-1) ? fabs(ves1->data[r][c-1])+dv : vs1;
    double vr = c < _nrCols-1 && !MV(r,c+1) ? fabs(ves1->data[r][c+1]) + dv : vs1;

    bool fv1 = (vs1 >= vu && vs1 >= vd && vs1 >= vl && vs1 >= vr);

    vu = r > 0 && !MV(r-1,c) ? fabs(ves2->data[r-1][c])+dv : vs2;
    vd = r < _nrRows-1 && !MV(r+1,c) ? fabs(ves2->data[r+1][c])+dv : vs2;
    vl = c > 0 && !MV(r,c-1) ? fabs(ves2->data[r][c-1])+dv : vs2;
    vr = c < _nrCols-1 && !MV(r,c+1) ? fabs(ves2->data[r][c+1]) + dv : vs2;

    bool fv2 = (vs2 >= vu && vs2 >= vd && vs2 >= vl && vs2 >= vr);

    if (vs1 > thv || vs2 > thv || Vsqrt > thv || fv1 || fv2) {
        double vh = hes->Drc/dt;
        double vkin = sqrt(qPow(hes->Drc, 2.0/3.0)*sqrt(Grad->Drc)/N->Drc);
        ves1->Drc = sign1*std::min(std::min(vh, vkin), vs1);
        ves2->Drc = sign2*std::min(std::min(vh, vkin), vs2);
    }
}
//---------------------------------------------------------------------------

// fill the left and right h,u,v arrays to solve the 1D scheme, also used as prep for the 2D MUSCL scheme for the boundary cells
void TWorld::simpleScheme(int thread, cTMap *_h,cTMap *_u,cTMap *_v)
{

    //FOR_ROW_COL_UF2DMT_DT
    FOR_ROW_COL_UF2DMTDER{
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
    }}}}
}
//---------------------------------------------------------------------------
// h > 0 plus surrounding cells
void TWorld::setFloodMask(cTMap * h)
{
    fill(*FloodHMaskDer, 0);
    FOR_ROW_COL_MV
    {
        if(h->data[r][c] > HMIN)
        {
            FloodHMaskDer->Drc = 1.0;

            if(!OUTORMV(r+1,c))
            {
                FloodHMaskDer->data[r+1][c] = 1.0;
            }
            if(!OUTORMV(r-1,c))
            {
                FloodHMaskDer->data[r-1][c] = 1.0;
            }
            if(!OUTORMV(r,c+1))
            {
                FloodHMaskDer->data[r][c+1] = 1.0;
            }
            if(!OUTORMV(r,c-1))
            {
                FloodHMaskDer->data[r][c-1] = 1.0;
            }
        }
    }

    int rc = 0;
    int cc = 0;
    int count = 0;

    FloodHR->data[0][0] = -1;
    FloodHC->data[0][0] = -1;

    FOR_ROW_COL_MV
    {
        if(FloodHMaskDer->Drc > 1e-12)
        {
            FloodHR->data[rc][cc] = r;
            FloodHC->data[rc][cc] = c;

            count ++;

            cc ++;
            if(cc == _nrCols)
            {
                rc ++;
                cc = 0;
            }

            if(INSIDE(rc,cc))
            {
                FloodHR->data[rc][cc] = -1;
                FloodHC->data[rc][cc] = -1;
            }
        }

    }

}
//---------------------------------------------------------------------------
void TWorld::setFloodMaskDT(cTMap * DT)
{

    FloodDTR->data[0][0] = -1;
    FloodDTC->data[0][0] = -1;

    int rc = 0;
    int cc = 0;
    int count = 0;

    for(int rc2 = 0; rc2 < _nrRows; rc2++)
    {
        bool out= false;
        for (int cc2 = 0; cc2 < _nrCols; cc2++)
        {
            int r = (int) (FloodHR->data[rc2][cc2]);
            int c = (int) (FloodHC->data[rc2][cc2]);

            if(!INSIDE(r,c)){out = true; break;}

            if(DT->Drc > 1e-12)
            {
                FloodDTR->data[rc][cc] = r;
                FloodDTC->data[rc][cc] = c;

                count ++;

                cc ++;
                if(cc == _nrCols)
                {
                    rc ++;
                    cc = 0;
                }

                if(INSIDE(rc,cc))
                {
                    FloodDTR->data[rc][cc] = -1;
                    FloodDTC->data[rc][cc] = -1;
                }
            }
        }
        if(out){break;}
    }
}

//---------------------------------------------------------------------------
// t = timesum
void TWorld::setFloodDT(double t, cTMap * h)
{
    double dt_max = std::min(_dt, _dx*0.5);
    double dtmin = dt_max;

    //find largest velocity and determine dt
    for(int rc = 0; rc < _nrRows; rc++)
    {
        bool out = false;
        for (int cc = 0; cc < _nrCols; cc++)
        {
            int r = (int) (FloodHR->data[rc][cc]);
            int c = (int) (FloodHC->data[rc][cc]);
            if(!INSIDE(r,c)){out = true; break;}

            if(h->Drc > HMIN && FloodHMaskDer->Drc != 0.0)
            {
               dtmin = std::min(FloodDT->Drc,dtmin);
            }
        }
       if(out){break;}
    }


    dtmin = std::min(_dt-t,std::max(TimestepfloodMin ,dtmin));

    Flood_DTMIN = dtmin;

    for(int rc = 0; rc < _nrRows; rc++)
    {
        bool out = false;
        for (int cc = 0; cc < _nrCols; cc++)
        {
            int r = (int) (FloodHR->data[rc][cc]);
            int c = (int) (FloodHC->data[rc][cc]);
            if(!INSIDE(r,c)){out = true; break;}

            if(h->Drc > HMIN)
            {
                //for homogeneous timestep!!

                //FloodDT->Drc = dtmin;

            }else if(FloodHMaskDer->Drc != 0.0)
            {
                //for really small h, dt becomes small. We can not ignore due to possible incoming flow, so set to minimum of other cells
                FloodDT->Drc = dtmin;//_max;
            }

            double dtm = dt_max; //_dt;
            if(FloodHMaskDer->Drc != 0.0)
            {
//                for (int i = -1; i<2; i++)
//                    for (int j = -1; j < 2; j++)
//                        if(!OUTORMV(r+i,c+j)) {
//                            double weight = 1.0;
//                            if (abs(i*j) >= 1)
//                                weight = 1.5;
//                            dtm = std::min(dtm, weight*FloodDT->data[r+i][c+j]);
//                        }

                if(!OUTORMV(r,c-1))
                    dtm = std::min(dtm, FloodDT->data[r][c-1]);
                if(!OUTORMV(r+1,c))
                    dtm = std::min(dtm, FloodDT->data[r+1][c]);
                if(!OUTORMV(r-1,c))
                    dtm = std::min(dtm, FloodDT->data[r-1][c]);
                if(!OUTORMV(r,c+1))
                    dtm = std::min(dtm, 1.5*FloodDT->data[r][c+1]);
                if(!OUTORMV(r+1,c+1))
                    dtm = std::min(dtm, 1.5*FloodDT->data[r+1][c+1]);
                if(!OUTORMV(r+1,c-1))
                    dtm = std::min(dtm, 1.5*FloodDT->data[r+1][c-1]);
                if(!OUTORMV(r-1,c+1))
                    dtm = std::min(dtm, 1.5*FloodDT->data[r-1][c+1]);
                if(!OUTORMV(r-1,c-1))
                    dtm = std::min(dtm, 1.5*FloodDT->data[r-1][c-1]);

                if(!OUTORMV(r,c+2))
                    dtm = std::min(dtm,2*FloodDT->data[r][c+2]);
                if(!OUTORMV(r,c-2))
                    dtm = std::min(dtm,2*FloodDT->data[r][c-2]);
                if(!OUTORMV(r+2,c))
                    dtm = std::min(dtm,2*FloodDT->data[r+2][c]);
                if(!OUTORMV(r-2,c))
                    dtm = std::min(dtm,2*FloodDT->data[r-2][c]);
            }
            FloodDTr->Drc = dtm;
        }
        if(out){break;}
    }

    for(int rc = 0; rc < _nrRows; rc++)
    {
        bool out  = false;
        for (int cc = 0; cc < _nrCols; cc++)
        {
            int r = (int) (FloodHR->data[rc][cc]);
            int c = (int) (FloodHC->data[rc][cc]);
            if(!INSIDE(r,c)){out = true; break;}

            if(FloodHMaskDer->Drc != 0)
            {
                if (!SwitchVariableTimestep)
                FloodDT->Drc = dtmin;
                else
                FloodDT->Drc = FloodDTr->Drc;

                //determine wether to do the timestep now or later
                if(!(t + dtmin < _dt)) // if t > _dt do the remaining timestep
                {
                    FloodDT->Drc = std::max(0.0,_dt-FloodT->Drc);

                } else
                    if(!(FloodT->Drc > t + dtmin))
                    {
                        FloodDT->Drc = std::min(std::max(dtmin, FloodDT->Drc), _dt - FloodT->Drc);
                       // FloodDT->Drc = std::min(FloodDT->Drc, _dt - FloodT->Drc);

                    } else {
                        FloodDT->Drc = 0;
                    }
            } else {
                FloodDT->Drc = 0;
            }
        }
        if(out){break;}
    }


}

//---------------------------------------------------------------------------
// 2nd order without iteration dt1, dt2!
double TWorld::fullSWOF2Do2light(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct)
{
    double dt1 = 0, timesum = 0;
    double dt_max = std::min(_dt, _dx*dtmaxfrac);
    int n = 0;
    double sumh = 0;
    bool stop;

    SwitchHeun = false;

    if (startFlood)
    {
        if (SwitchErosion) {
            FOR_ROW_COL_MV {
                SSFlood->Drc += DETSplash->Drc;
                SSCFlood->Drc = MaxConcentration(ChannelAdj->Drc * DX->Drc * h->Drc, &SSFlood->Drc, &DepFlood->Drc);
                // recalc concentration
            }
            // !!!!!!!!!!  Sed_D.Drcd += DETSplash->Drc * W_D.Drcd;
        }

        fill(*FloodT,0.0);

        if (correct)
            sumh = getMass(h);

        do {

            dt1 = dt_max;

            fill(*FloodDT,0.0);
            fill(*FloodDTr,0.0);

            //mask flow domain: h>0 and one cell more
            setFloodMask(h);

            //now set the height-dependent mask for FOR_ROW_COL_UFMT_DT
            ThreadPool->SetMask(DEM,FloodHMaskDer,FloodHR,FloodHC);

            //creation callable function object by binding TWorld oject to member function
            //this leaves only the thread id as a parameter for the compute function
            // MUSCL and maincalcflux (fluxes between cells, Rieman stuff, boundaries and smallest dt)
            flood_cellcompute = std::bind((&TWorld::fullSWOF2Do2lightWrapperCell1),this,std::placeholders::_1,h,u,v,z);
            ThreadPool->RunDynamicCompute(flood_cellcompute);
            ThreadPool->WaitForAll();
            // MUSCL and mainflux, dt1 not smallest because of threading !

            setFloodDT(timesum, h);
            dt1 = Flood_DTMIN; // smallest dt in domain for this timestep

            setFloodMaskDT(FloodDT);
            //now set the timestep-dependent mask for FOR_ROW_COL_UFMT_DT
            ThreadPool->SetMask(DEM,FloodDT,FloodDTR,FloodDTC);

            //creation callable function object by binding TWorld oject to member function
            //this leaves only the thread id as a parameter for the compute function
            //add dt1 as one of the parameters
            flood_flowcompute = std::bind((&TWorld::fullSWOF2Do2lightWrapperDynamic1),this,std::placeholders::_1,h,u,v,hs,us,vs, dt1);
            ThreadPool->RunDynamicCompute(flood_flowcompute);
            ThreadPool->WaitForAll();
            //maincalcscheme with dt1 or flooddt->Drc

           // do a second time and do Heun
            if (SwitchHeun)
            {
                // do it all a second time, but not erosion
                setFloodMask(h);
                ThreadPool->SetMask(DEM,FloodHMaskDer,FloodHR,FloodHC);

                flood_cellcompute = std::bind((&TWorld::fullSWOF2Do2lightWrapperCell1),this,std::placeholders::_1,hs,us,vs,z);
                ThreadPool->RunDynamicCompute(flood_cellcompute);
                ThreadPool->WaitForAll();

                setFloodDT(timesum, h);
                dt1 = Flood_DTMIN;

                setFloodMaskDT(FloodDT);
                ThreadPool->SetMask(DEM,FloodDT,FloodDTR,FloodDTC);

                flood_flowcompute2 = std::bind((&TWorld::fullSWOF2Do2lightWrapperDynamic2),this,std::placeholders::_1,hs,us,vs,hsa,usa,vsa, dt1);
                ThreadPool->RunDynamicCompute(flood_flowcompute2);
                ThreadPool->WaitForAll();
                // Do Heun
                FOR_ROW_COL_MV {
                    double havg = 0.5*(h->Drc + hsa->Drc);
                    if (havg >= he_ca) {
                        double q1 = 0.5*(h->Drc*u->Drc + hsa->Drc*usa->Drc);
                        u->Drc = q1/havg;
                        double q2 = 0.5*(h->Drc*v->Drc + hsa->Drc*vsa->Drc);
                        v->Drc = q2/havg;
                        h->Drc = havg;
                    }
                    else {
                        u->Drc = 0;
                        v->Drc = 0;
                        h->Drc = 0;
                    }
                }
            } //heun

            timesum = timesum + Flood_DTMIN;
            stop = timesum  > _dt-1e-6;

            FOR_ROW_COL_MV
            {
                FloodT->Drc += FloodDT->Drc;
            }

            n++;
            if (n > F_MaxIter)
                break;
            // qDebug() << n << timesum << dt1;
        } while (!stop);

        correctMassBalance(sumh, h);

        FOR_ROW_COL_MV {
            SWOFSedimentSetConcentration(r,c,h);
        }
    } // if floodstart

    iter_n = n;
    dt1 = n > 0? _dt/n : dt1;
    return(dt1);
}
//---------------------------------------------------------------------------
void TWorld::fullSWOF2Do2lightWrapperCell1(int thread, cTMap *h, cTMap *u, cTMap *v, cTMap *z)
{
    double dt_max = std::min(_dt, _dx*dtmaxfrac);
    double dt1 = dt_max;

    setZero(thread,h, u, v);
    simpleScheme(thread, h,u,v);
    if (SwitchMUSCL)
        MUSCL(thread,h,u,v,z);

    maincalcflux(thread, dt1, dt_max);
    // because of threading this does not result in smallest dt!!!

    return;
}
//---------------------------------------------------------------------------
void TWorld::fullSWOF2Do2lightWrapperDynamic1(int thread, cTMap *h, cTMap *u, cTMap *v,
                                              cTMap *hs, cTMap *us, cTMap *vs,
                                               double dt1)
{
    //st venant equations, h, u, v go in hs, vs, us come out
    maincalcscheme(thread, dt1, h,u,v, hs,us,vs);

    setZero(thread,hs, us, vs);

    //only when sediment is modelled
    if (SwitchErosion)
        SWOFSediment(thread,FloodDT,h,u,v );  //TODO why not hs, us, vs

    if (!SwitchHeun)
    {
        FOR_ROW_COL_UF2DMT_DT {
            if(FloodHMaskDer->Drc != 0) {
                h->Drc = hs->Drc;
                u->Drc = us->Drc;
                v->Drc = vs->Drc;
            }
        }}}}
    }
}
//---------------------------------------------------------------------------
/// MC Dynamic processes for Heun, do all a second time and calculate time average (Heun method)
void TWorld::fullSWOF2Do2lightWrapperDynamic2(int thread, cTMap *hs, cTMap *us, cTMap *vs,
                                              cTMap *hsa, cTMap *usa, cTMap *vsa,
                                               double dt1)
{
    //st venant equations, h, u, v go in hs, vs, us come out
    maincalcscheme(thread, dt1, hs,us,vs, hsa,usa,vsa);

    setZero(thread,hsa, usa, vsa);


}
//---------------------------------------------------------------------------
void TWorld::fullSWOF2Do2lightWrapperErosion(int thread, cTMap *h, cTMap *u, cTMap *v, double dt1)
{

    //sediment
    SWOFSediment(thread,FloodDT,h,u,v );


}
