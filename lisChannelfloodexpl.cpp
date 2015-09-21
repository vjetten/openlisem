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
 \file lisChannelfloodexpl.cpp
 \brief Channel flood using an explicit solution of the St Venant equations following Bates et al. 2010.\n

functions: \n
- double TWorld::floodExplicit() Do flooding with explicit solution Bates et al (old lisflood)
*/

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"


#define he_ca 1e-12
#define ve_ca 1e-12

#define dt_ca 0.005
#define dt_fix 0.05

#define GRAV 9.8067
#define EPSILON 1e-6
#define scheme_type 1   //return calculated or fixed dt
#define MAXITER 100

void TWorld::ENOws(int wsnr, cTMap *h ,cTMap *u, cTMap *v, cTMap *z)
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
    FOR_WATERSHED_ROW_COL(wsnr)
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
    FOR_WATERSHED_ROW_COL(wsnr)
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
/// This MUSCL creates left and right u,v,h, arrays for in mainfluxcalc
void TWorld::MUSCLws(int wsnr, cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z)
{
    double delta_h1, delta_u1, delta_v1;
    double delta_h2, delta_u2, delta_v2;
    double dh, du, dv, dz_h;

    delta_u1 = 0;
    delta_v1 = 0;
    delta_h1 = 0;
    delta_u2 = 0;
    delta_v2 = 0;
    delta_h2 = 0;

    FOR_WATERSHED_ROW_COL(wsnr)
    {
        if(c > 0 && c < _nrCols-1 && !MV(r,c-1) && !MV(r,c+1))
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
        }
    }}

    delta_u1 = 0;
    delta_v1 = 0;
    delta_h1 = 0;
    delta_u2 = 0;
    delta_v2 = 0;
    delta_h2 = 0;

    FOR_WATERSHED_ROW_COL(wsnr)
    {
      //  if(r > 0 && !MV(r-1,c))
            if(r > 0 && r < _nrRows-1 && !MV(r-1,c) && !MV(r+1,c))
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
        }
    }}
}
void TWorld::setZerows(int wsnr, cTMap *_h, cTMap *_u, cTMap *_v)
{
  FOR_WATERSHED_ROW_COL(wsnr) {
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
// St Venant equations: conservation of mass and momentum: h u et v
void TWorld::maincalcschemews(int wsnr, cTMap *he, cTMap *ve1, cTMap *ve2,
                            cTMap *hes, cTMap *ves1, cTMap *ves2)
{
    FOR_WATERSHED_ROW_COL(wsnr) {
        double dx = ChannelAdj->Drc;
        double dy = DX->Drc;
        long double tx = WS[wsnr].dt/dx; //dt->Drc/dx;
        long double ty = WS[wsnr].dt/dy;   //dt->Drc/dy;
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

            double nsq = N->Drc*N->Drc*GRAV*sqrt(ve1->Drc*ve1->Drc+ve2->Drc*ve2->Drc)*WS[wsnr].dt/qPow(hes->Drc,4.0/3.0);

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
 of a domain based on spatial averages of the previous timestep, this is the hydrostatic equilibrium
 calculaton of the time steps in relation to cfl.
*/
void TWorld::maincalcfluxws(int wsnr)
{
    FOR_WATERSHED_ROW_COL(wsnr) {
        if (c > 0 && !pcr::isMV(LDD->data[r][c-1]))
        {
            h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1]));
            h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1]));
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
        }
        else
        {
            h1d->Drc = std::max(0.0, h1r->Drc - std::max(0.0,  delz1->Drc));
            h1g->Drc        = std::max(0.0, h1l->Drc        - std::max(0.0, -delz1->Drc));
            if (F_scheme == 1)
                F_Rusanov(h1d->Drc, u1r->Drc, v1r->Drc,h1g->Drc, u1l->Drc, v1l->Drc);
            else
                if (F_scheme == 2)
                    F_HLL(h1d->Drc, u1r->Drc, v1r->Drc,h1g->Drc, u1l->Drc, v1l->Drc);
                else
                    F_HLL2(h1d->Drc, u1r->Drc, v1r->Drc,h1g->Drc, u1l->Drc, v1l->Drc);

            f1->Drc = HLL2_f1;
            f2->Drc = HLL2_f2;
            f3->Drc = HLL2_f3;
            cflx->Drc = HLL2_cfl;
        }
    }}

    FOR_WATERSHED_ROW_COL(wsnr) {
        if(r > 0 && !pcr::isMV(LDD->data[r-1][c]))
        {

            h2d->data[r-1][c] = std::max(0.0, h2r->data[r-1][c] - std::max(0.0,  delz2->data[r-1][c]));
            h2g->Drc          = std::max(0.0, h2l->Drc          - std::max(0.0, -delz2->data[r-1][c]));

            if (F_scheme == 1)
                F_Rusanov(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);
            else
                if (F_scheme == 2)
                    F_HLL(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);
                else
                    F_HLL2(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c],h2g->Drc,v2l->Drc,u2l->Drc);

            g1->Drc = HLL2_f1;
            g2->Drc = HLL2_f3;
            g3->Drc = HLL2_f2;
            cfly->Drc = HLL2_cfl;
        }
        else
        {
            h2d->Drc = std::max(0.0, h2r->Drc - std::max(0.0,  delz2->Drc));
            h2g->Drc        = std::max(0.0, h2l->Drc        - std::max(0.0, -delz2->Drc));

            if (F_scheme == 1)
                F_Rusanov(h2d->Drc,v2r->Drc,u2r->Drc,h2g->Drc,v2l->Drc,u2l->Drc);
            else
                if (F_scheme == 2)
                    F_HLL(h2d->Drc,v2r->Drc,u2r->Drc,h2g->Drc,v2l->Drc,u2l->Drc);
                else
                    F_HLL2(h2d->Drc,v2r->Drc,u2r->Drc,h2g->Drc,v2l->Drc,u2l->Drc);

            g1->Drc = HLL2_f1;
            g2->Drc = HLL2_f3;
            g3->Drc = HLL2_f2;
            cfly->Drc = HLL2_cfl;
        }
    }}

#define AVG(a1,a2) ((a1*a2 > 0) ? std::sqrt(a1*a2) : -std::sqrt(std::abs(a1*a2)))

    if (F_replaceV > 0)
    {
        FOR_WATERSHED_ROW_COL(wsnr)
                if (cflx->Drc > F_maxVelocity || cflx->Drc > F_maxVelocity)
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
        }}
    }
}
//---------------------------------------------------------------------------
/// find the smallest dt in the catchment
void TWorld::findDTws(int wsnr, bool two)
{
    double dt_tmp, dtx, dty;
    double dt_max = std::min(_dt, _dx*0.5);
    double dt = dt_max;
    dtx = dt_max;
    dty = dt_max;


    FOR_WATERSHED_ROW_COL(wsnr)
        double dx = ChannelAdj->Drc;
        if (qFabs(cflx->Drc*dt/dx) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = courant_factor*dx/cflx->Drc;
        dtx = std::min(std::min(dt, dt_tmp), dtx);

        double dy = DX->Drc;
        if (qFabs(cfly->Drc*dt/dy) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = courant_factor*dy/cfly->Drc;

        dty = std::min(std::min(dt, dt_tmp), dty);

        dt = (std::max(dt_ca, std::min(dtx,dty)));
        dt = std::min(dt, _dt-WS[wsnr].dtsum);
    }

    if (two)
        WS[wsnr].dt2 = dt;
    else
        WS[wsnr].dt = dt;
}
//---------------------------------------------------------------------------
// fill the left and right h,u,v arrays to solve the 1D scheme, also used as prep for the 2D MUSCL scheme for the boundary cells
void TWorld::simpleSchemews(int wsnr, cTMap *_h,cTMap *_u,cTMap *_v)
{

    FOR_WATERSHED_ROW_COL(wsnr) {
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
double TWorld::fullSWOF2Do1ws(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct)
{
    //double timesum = 0;
    double dt_max = std::min(_dt, _dx*0.5);
    double dt1 = dt_max;
    double sumh = 0;
    double sumn = 0;
    bool done = false;
    int n = 0;

    // do one time only at the start of simulation
    if (prepareFlood)
        prepareFloodZ(z);


    for (int l = 1; l < WS.count(); l++)
    {
        WS[l].flood = false;
        WS[l].dt = dt_max;
        WS[l].dt2 = dt_max;
        WS[l].dtsum = 0;

        if (correct)
            sumh = mapTotal(*h);

        // check if flooding in this watershed else skip
        for (long k = 0; k < WS[l].cr.count(); k++)
        {
            int c = WS[l].cr[k]._c;
            int r = WS[l].cr[k]._r;
            if (hmx->Drc > 0)
            {
                WS[l].flood = true;
                break;
            }
        }

        if (WS[l].flood)
        {
            do {
                setZerows(l, h, u, v);

                simpleSchemews(l, h, u, v);

                maincalcfluxws(l);
                findDTws(l, false);
                maincalcschemews(l, h,u,v, hs,us,vs);

                setZerows(l, hs, us, vs);

                FOR_WATERSHED_ROW_COL(l){
                    h->Drc = hs->Drc;
                    u->Drc = us->Drc;
                    v->Drc = vs->Drc;
                }}

                done = true;
                WS[l].dtsum += WS[l].dt;
                WS[l].dtsum = std::min(_dt, WS[l].dtsum);
                if (WS[l].dtsum < _dt)
                    done = false;

                n++;
                if (n > F_MaxIter)
                    done = true;

                if (correct)
                    correctMassBalance(sumh, h, 1e-12);

            } while (!done);

            sumn = sumn + n;
        }  // if flood in WS
    }  // for WS

    //calc average dt and iterations for screen output
    n = 0;
    for (int l = 1; l < WS.count(); l++)
        if (WS[l].flood)
            n++;
    dt1 = n > 0 && sumn > 0? _dt/(sumn/n) : 0;
    iter_n = n > 0 && sumn > 0? int(sumn/n) : 0;
    return(dt1);
}
//---------------------------------------------------------------------------
double TWorld::fullSWOF2Do2ws(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct)//, cTMap *q1, cTMap *q2)
{
    double dt_max = std::min(_dt, _dx*0.5);
    double dt1;
    double sumn = 0;
    int n = 0;
    double sumh = 0;
    bool done;

    // only once at the beginning of the run
    if (prepareFlood)
        prepareFloodZ(z);

    if (correct)
        sumh = mapTotal(*h);

    for (int l = 1; l < WS.count(); l++)
    {
        WS[l].flood = false;
        WS[l].dt = dt_max;
        WS[l].dt2 = dt_max;
        WS[l].dtsum = 0;

        for (long k = 0; k < WS[l].cr.count(); k++)
        {
            int c = WS[l].cr[k]._c;
            int r = WS[l].cr[k]._r;
            if (hmx->Drc > 0)
            {
                WS[l].flood = true;
                break;
            }
        }

        if (correct)
            sumh = mapTotal(*h);

        verif = 1;

        if (WS[l].flood)
        {
            do {

                if (verif == 1)
                {

                    setZerows(l, h, u, v);
                    // Reconstruction order 2
                    // makes h1r, h1l, u1r, u1l, v1r, v1l
                    // makes h2r, h2l, u2r, u2l, v2r, v2l
                    // makes delzc1, delzc2, delz1, delz2
                    simpleSchemews(l, h,u,v);
                    if (F_scheme == (int)FMUSCL)
                        MUSCLws(l, h,u,v,z);
                    else
                        ENOws(l, h,u,v,z);
                }

                maincalcfluxws(l);
                // prepare potentials
                findDTws(l, false);
                //find smallest dt
                maincalcschemews(l, h,u,v, hs,us,vs);
                //st venant equations

                setZerows(l, hs, us, vs);
                simpleSchemews(l, h,u,v);
                if (F_scheme == (int)FMUSCL)
                    MUSCLws(l, hs,us,vs,z);
                else
                    ENOws(l, hs,us,vs,z);
                //Reconstruction for order 2

                maincalcfluxws(l);
                // prepare potentials
                findDTws(l, true);
                //find smallest dt

                if (WS[l].dt2 < WS[l].dt)
                {
                    //qDebug() << WS[l].dt << WS[l].dt2;
                    WS[l].dt = WS[l].dt2;
                    verif = 0;

                }
                else
                {
                    verif = 1;

                    maincalcschemews(l, hs, us, vs, hsa, usa, vsa);
                    setZerows(l, hsa, usa, vsa);

                    //Heun method (order 2 in time)
                    FOR_WATERSHED_ROW_COL(l) {
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
                    }}

                    done = true;
                    WS[l].dtsum += WS[l].dt;
                    WS[l].dtsum = std::min(_dt, WS[l].dtsum);
                    if (WS[l].dtsum < _dt)
                        done = false;

                    n++;
                    if (n > F_MaxIter)
                        done = true;

                } //end for else dt2<dt1

                if (correct)
                    correctMassBalance(sumh, h, 1e-12);

            } while (!done);

            sumn = sumn + n;

//            ChannelOverflowWS(l, _dt, h);

        } // if WS flood
    }  // for all WS

    //calc average dt and iterations for screen output
    n = 0;
    for (int l = 1; l < WS.count(); l++)
        if (WS[l].flood)
            n++;
    dt1 = n > 0 && sumn > 0? _dt/(sumn/n) : 0;
    iter_n = n > 0 && sumn > 0? int(sumn/n) : 0;
    return(dt1);
}

//---------------------------------------------------------------------------

