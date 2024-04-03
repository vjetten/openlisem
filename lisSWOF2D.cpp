/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010 - 2013  Victor Jetten
**  contact: v.g.jetten@utwente.nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: VC2010/Qt
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

#include <algorithm>
#include "lisemqt.h"
#include "model.h"
#include "operation.h"
#include "global.h"

#define dt_ca 0.001

#define dtmaxfrac 0.5


//---------------------------------------------------------------------------
/// MUSCL: Monotone Upstream-centered Schemes for Conservation Laws
/// see http://en.wikipedia.org/wiki/MUSCL_scheme
///
/// This MUSCL creates left and right u,v,h, arrays for in mainfluxcalc
void TWorld::MUSCL(cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z)
{
    double delta_h1, delta_u1, delta_v1, dz1;
    double delta_h2, delta_u2, delta_v2, dz2;
    double dh, du, dv, dz_h, hlh, hrh;

    FOR_ROW_COL_MV {
        if(_h->Drc > he_ca) {
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
            // ?
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
    }

    FOR_ROW_COL_MV {
        if(_h->Drc > he_ca)
        {
            if(r > 0 && !MV(r-1,c) && r < _nrRows-1 && !MV(r+1, c)) {
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
    }
}
//---------------------------------------------------------------------------

    /**
     Construction of variables for hydrostatic reconstruction using HLL or Rusanov
     Flux in x and y direction. The method is based solving PDEs with an estimate of the status
     of a domain based on spatial averages of the previous timestep, this is the hydrostatic equlibrium
     Calculaton of the time steps in relation to cfl.
    */

double TWorld::maincalcflux(cTMap *_h,double dt, double dt_max)
{
    vec4 rec;
    double dtx, dty;
    double dt_tmp;

    cTMap *fbw = FlowBarrierW;
    cTMap *fbe = FlowBarrierE;
    cTMap *fbn = FlowBarrierN;
    cTMap *fbs = FlowBarrierS;

    FOR_ROW_COL_MV {
       // if(_h->Drc > he_ca) {
            f1->Drc = 0;
            f2->Drc = 0;
            f3->Drc = 0;
            f1o->Drc = 0;
            f2o->Drc = 0;
            f3o->Drc = 0;

            if(c > 0 && !MV(r,c-1)) {
                h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
                h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
                rec = F_HLL2(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
                f1->Drc =   rec.v[0];
                f2->Drc =   rec.v[1];
                f3->Drc =   rec.v[2];
                cflx->Drc = rec.v[3];
            } else {
                double _h1g = std::max(0.0, h1l->Drc - FlowBarrierE->Drc);
                rec = F_HLL2(0,0,0, _h1g, u1l->Drc, v1l->Drc);
                f1->Drc =   rec.v[0];
                f2->Drc =   rec.v[1];
                f3->Drc =   rec.v[2];
                cflx->Drc = rec.v[3];
            }
            // right hand side boundary
            if(c == _nrCols-1 || MV(r, c+1)){
                double _h1d = std::max(0.0, h1r->Drc - fbw->Drc);
                rec = F_HLL2(_h1d,u1r->Drc,v1r->Drc,0.,0.,0.);
                f1o->Drc = rec.v[0];
                f2o->Drc = rec.v[1];
                f3o->Drc = rec.v[2];
            }
      //  }
    }

    FOR_ROW_COL_MV {
   //     if(_h->Drc > he_ca){
            g1->Drc = 0;
            g2->Drc = 0;
            g3->Drc = 0;
            g1o->Drc = 0;
            g2o->Drc = 0;
            g3o->Drc = 0;

            if(r > 0 && !MV(r-1,c)) {
                h2d->data[r-1][c] = std::max(0.0, h2r->data[r-1][c] - std::max(0.0,  delz2->data[r-1][c]  + std::max(fbs->Drc,fbn->data[r-1][c])));
                h2g->Drc          = std::max(0.0, h2l->Drc          - std::max(0.0, -delz2->data[r-1][c]  + std::max(fbs->Drc,fbn->data[r-1][c])));
                rec = F_HLL2(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c], h2g->Drc,v2l->Drc,u2l->Drc);

                g1->Drc = rec.v[0];
                g2->Drc = rec.v[2];
                g3->Drc = rec.v[1];
                cfly->Drc = rec.v[3];
            } else {
                double _h2g = std::max(0.0, h2l->Drc - fbn->Drc);
                rec = F_HLL2(0,0,0,_h2g,v2l->Drc,u2l->Drc);

                g1->Drc = rec.v[0];
                g2->Drc = rec.v[2];
                g3->Drc = rec.v[1];
                cfly->Drc = rec.v[3];
            }
            // left hand side boundary
            if (r == _nrRows-1 || MV(r+1, c)) {
                double _h2d = std::max(0.0, h2d->Drc - fbs->Drc);
                rec = F_HLL2(_h2d,v2l->Drc,u2l->Drc,0.,0.,0.);
                g1o->Drc = rec.v[0];
                g2o->Drc = rec.v[2];
                g3o->Drc = rec.v[1];
            }
     //   }
    }

    // find largest velocity and determine dt
    dtx = dt_max;
    dty = dt_max;

    FOR_ROW_COL_MV
       if (_h->Drc > he_ca)
    {
        double dx = _dx;//ChannelAdj->Drc;//FlowWidth->Drc;//
        if (qFabs(cflx->Drc*dt/dx) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = courant_factor*dx/cflx->Drc;
        dtx = std::min(std::min(dt, dt_tmp), dtx);
    }

    FOR_ROW_COL_MV
        if (_h->Drc > he_ca)
    {
        double dy = _dx;//DX->Drc;
        if (qFabs(cfly->Drc*dt/dy) < 1e-10)
            dt_tmp = dt_max;
        else
            dt_tmp = courant_factor*dy/cfly->Drc;
        dty = std::min(std::min(dt, dt_tmp), dty);
    }

    return(std::max(TimestepfloodMin, std::min(dtx,dty)));
}
//---------------------------------------------------------------------------
// St Venant equations: conservation of mass and momentum: h u et v (= he, ve1, ve2)
void TWorld::maincalcscheme(double dt, cTMap *he, cTMap *ve1, cTMap *ve2,
                            cTMap *hes, cTMap *ves1, cTMap *ves2)
{
    FOR_ROW_COL_MV {
        double tx = dt/ChannelAdj->Drc;
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



        hes->Drc = std::max(0.0, he->Drc - tx*_f1 + tx*f1->Drc - ty*_g1 + ty*g1->Drc);

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


            double sqUV = qSqrt(ve1->Drc*ve1->Drc+ve2->Drc*ve2->Drc);
            double nsq1 = (0.001+N->Drc)*(0.001+N->Drc)*GRAV/std::max(0.01, qPow(hes->Drc,4.0/3.0));
            double nsq = nsq1*sqUV*dt;

            ves1->Drc = (qes1/(1.0+nsq))/std::max(0.01, hes->Drc);
            ves2->Drc = (qes2/(1.0+nsq))/std::max(0.01, hes->Drc);

            if (SwitchTimeavgV) {
                double fac = 0;
                fac = 0.5+0.5*std::min(1.0,4*hes->Drc)*std::min(1.0,4*hes->Drc);
                fac = fac *exp(- std::max(1.0,dt) / nsq1);
                ves1->Drc = fac * ve1->Drc + (1.0-fac) *ves1->Drc;
                ves2->Drc = fac * ve2->Drc + (1.0-fac) *ves2->Drc;
            }

            double vmax = 0.25*_dx/dt;
            ves1->Drc = std::max(-vmax, std::min(vmax, ves1->Drc));
            ves2->Drc = std::max(-vmax, std::min(vmax, ves2->Drc));

        }
        else
        {
            // Case of height of water < ha.
            hes->Drc = 0;
            ves1->Drc = 0;
            ves2->Drc = 0;
        }
        // dan maar even met geweld!
        if (std::isnan(ves1->Drc) || std::isnan(ves2->Drc)  )
        {
            ves1->Drc = 0;
            ves2->Drc = 0;
            hes->Drc = 0;
        }
    }
}
//---------------------------------------------------------------------------


