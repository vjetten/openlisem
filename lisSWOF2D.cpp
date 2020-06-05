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
    vec4 rec;
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

            if(c > 0 && !MV(r,c-1)) {
                h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
                h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1]  + std::max(fbw->Drc,fbe->data[r][c-1])));
                rec = F_Riemann(h1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc, u1l->Drc, v1l->Drc);
                f1->Drc =   rec.v[0];
                f2->Drc =   rec.v[1];
                f3->Drc =   rec.v[2];
                cflx->Drc = rec.v[3];
            } else {
                double _h1g = std::max(0.0, h1l->Drc - fbe->Drc);
                rec = F_Riemann(0,0,0, _h1g, u1l->Drc, v1l->Drc);
                f1->Drc = rec.v[0];
                f2->Drc = rec.v[1];
                f3->Drc = rec.v[2];
                cflx->Drc = rec.v[3];
            }

            // right hand side boundary
            if(c == _nrCols-1 || MV(r, c+1)){
                double _h1d = std::max(0.0, h1r->Drc - fbw->Drc);
                rec = F_Riemann(_h1d,u1r->Drc,v1r->Drc,0.,0.,0.);
                f1o->Drc = rec.v[0];
                f2o->Drc = rec.v[1];
                f3o->Drc = rec.v[2];
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
                rec = F_Riemann(h2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c], h2g->Drc,v2l->Drc,u2l->Drc);
                g1->Drc = rec.v[0];
                g2->Drc = rec.v[1];
                g3->Drc = rec.v[2];
                cfly->Drc = rec.v[3];
            } else {
                double _h2g = std::max(0.0, h2l->Drc - fbn->Drc);
                rec = F_Riemann(0,0,0,_h2g,v2l->Drc,u2l->Drc);
                g1->Drc = rec.v[0];
                g2->Drc = rec.v[1];
                g3->Drc = rec.v[2];
                cflx->Drc = rec.v[3];
            }
            // left hand side boundary
            if (r == _nrRows-1 || MV(r+1, c)) {
                double _h2d = std::max(0.0, h2d->Drc - fbs->Drc);
                rec = F_Riemann(_h2d,v2l->Drc,u2l->Drc,0.,0.,0.);
                g1o->Drc = rec.v[0];
                g2o->Drc = rec.v[1];
                g3o->Drc = rec.v[2];
            }
        }
    }}}}

    // find largest velocity and determine dt
    dt = dt_max;
    Flood_DTMIN = dt_max;
    FOR_ROW_COL_UF2DMT_DT {
        dtx = dt_max;
        dty = dt_max;
        if(FloodHMaskDer->Drc != 0){
            double dx = ChannelAdj->Drc;// FlowWidth->Drc;
            if (qFabs(cflx->Drc*dt/dx) > 1e-10)
                dtx = std::min(dt_max, courant_factor*dx/cflx->Drc);

            double dy = DX->Drc;
            if (qFabs(cfly->Drc*dt/dy) > 1e-10)
                dty = std::min(dt_max, courant_factor*dy/cfly->Drc);

            FloodDT->Drc = std::min(dtx, dty);
            Flood_DTMIN = std::min(Flood_DTMIN, FloodDT->Drc);

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


                double sqUV = qSqrt(ve1->Drc*ve1->Drc+ve2->Drc*ve2->Drc);
                double nsq1 = (0.001+N->Drc)*(0.001+N->Drc)*GRAV/qPow(hes->Drc,4.0/3.0);
                double nsq = nsq1*sqUV*dt;

                ves1->Drc = (qes1/(1.0+nsq))/hes->Drc;
                ves2->Drc = (qes2/(1.0+nsq))/hes->Drc;

                 correctSpuriousVelocities(r, c, hes, ves1, ves2);//,thv, dv, dt);

                double fac = 0;
                if (SwitchTimeavgV) {
                    fac = 0.5+0.5*std::min(1.0,4*hes->Drc)*std::min(1.0,4*hes->Drc);
                    fac = fac *exp(- std::max(1.0,dt) / nsq1);
                }
                ves1->Drc = fac * ve1->Drc + (1.0-fac) *ves1->Drc;
                ves2->Drc = fac * ve2->Drc + (1.0-fac) *ves2->Drc;

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
void TWorld::setFloodDT(cTMap * h)
{
    double dt_max = std::min(_dt, _dx*0.5);
    double dtmin = dt_max;

    for(int rc = 0; rc < _nrRows; rc++)
    {
        bool out = false;
        for (int cc = 0; cc < _nrCols; cc++)
        {
            int r = (int) (FloodHR->data[rc][cc]);
            int c = (int) (FloodHC->data[rc][cc]);
            //if(!INSIDE(r,c)){out = true; break;}  // no need, FloodHR/HC already is inside

            if (r < 0 || c < 0) {out = true; break;}

            if(h->Drc > HMIN)// && FloodHMaskDer->Drc != 0.0) //FloodHMaskDer = >HMIN + 1 cell
            {
               dtmin = std::min(FloodDT->Drc,dtmin);
            }
        }
       if(out){break;}
    }


    //dtmin = std::min(_dt-t,std::max(TimestepfloodMin ,dtmin));
    dtmin = std::max(TimestepfloodMin ,dtmin);
    // dtmin is smallest dt in whole map, the original dt1

  //  Flood_DTMIN = dtmin;

    if (!SwitchVariableTimestep) {
        FOR_ROW_COL_MV {
            if(FloodHMaskDer->Drc != 0)
                FloodDT->Drc = std::min(dtmin,_dt-FloodT->Drc);
            else
                FloodDT->Drc = 0;
        }
        return;
    }

    for(int rc = 0; rc < _nrRows; rc++)
    {
        bool out = false;
        for (int cc = 0; cc < _nrCols; cc++)
        {
            int r = (int) (FloodHR->data[rc][cc]);
            int c = (int) (FloodHC->data[rc][cc]);
           // if(!INSIDE(r,c)){out = true; break;}
            if(r < 0 || c < 0){out = true; break;}

            double dtm = dt_max;
            if(FloodHMaskDer->Drc != 0.0)
            {

                if(!OUTORMV(r,c-1))
                    dtm = std::min(dtm, FloodDT->data[r][c-1]);
                if(!OUTORMV(r+1,c))
                    dtm = std::min(dtm, FloodDT->data[r+1][c]);
                if(!OUTORMV(r-1,c))
                    dtm = std::min(dtm, FloodDT->data[r-1][c]);
                if(!OUTORMV(r,c+1))
                    dtm = std::min(dtm, FloodDT->data[r][c+1]);
                if(!OUTORMV(r+1,c+1))
                    dtm = std::min(dtm, 4*FloodDT->data[r+1][c+1]);
                if(!OUTORMV(r+1,c-1))
                    dtm = std::min(dtm, 4*FloodDT->data[r+1][c-1]);
                if(!OUTORMV(r-1,c+1))
                    dtm = std::min(dtm, 4*FloodDT->data[r-1][c+1]);
                if(!OUTORMV(r-1,c-1))
                    dtm = std::min(dtm, 4*FloodDT->data[r-1][c-1]);

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
         //   if(r < 0 || c < 0){out = true; break;}

            if(FloodHMaskDer->Drc != 0)
            {

                    FloodDT->Drc = std::max(dtmin, FloodDTr->Drc);

                    if (_dt-FloodT->Drc > dtmin)
                        FloodDT->Drc = std::min(FloodDT->Drc,_dt-FloodT->Drc);
                    if (FloodT->Drc >= _dt)
                        FloodDT->Drc = 0;

                //determine wether to do the timestep now or later
//                if(!(t + dtmin < _dt)) // if t > _dt do the remaining timestep
//                {
//                    FloodDT->Drc = std::max(0.0,_dt-FloodT->Drc);

//                } else
//                    if(!(FloodT->Drc > t + dtmin))
//                    {
//                        FloodDT->Drc = std::min(std::max(dtmin, FloodDT->Drc), _dt - FloodT->Drc);
//                       // FloodDT->Drc = std::min(FloodDT->Drc, _dt - FloodT->Drc);

//                    } else {
//                        FloodDT->Drc = 0;
//                    }

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
    double dt_max = std::min(_dt, _dx*0.5);
    int n = 0;
    double sumh = 0;
    bool stop;

  //  F_MaxIter = (int) std::min(1000.0, _dt/TimestepfloodMin);
fill(*tmb,0.0);
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

        //mask flow domain to start: h>0 and one cell more
        setFloodMask(h);

        do {

            dt1 = dt_max;

            fill(*FloodDT,0.0);
            fill(*FloodDTr,0.0);

            //now set the height-dependent mask for FOR_ROW_COL_UFMT_DT
            ThreadPool->SetMask(DEM,FloodHMaskDer,FloodHR,FloodHC);

            // MUSCL and maincalcflux (fluxes between cells, Rieman stuff, boundaries and smallest dt)
            flood_cellcompute = std::bind((&TWorld::fullSWOF2Do2lightWrapperCell1),this,std::placeholders::_1,h,u,v,z);
            ThreadPool->RunDynamicCompute(flood_cellcompute);
            ThreadPool->WaitForAll();
            // MUSCL and mainflux, dt1 not smallest because of threading !

            setFloodDT(h);
            dt1 = Flood_DTMIN; // smallest dt in domain for this timestep

            FOR_ROW_COL_MV {
                tma->Drc = 0;
                FloodT->Drc += FloodDT->Drc;
                if (FloodT->Drc > 0 && FloodT->Drc < _dt) {
                    tma->Drc = 1;
                    tmb->Drc += 1.0;
                }
            }
            setFloodMaskDT(tma);//FloodDT);
            //now set the timestep-dependent mask for FOR_ROW_COL_UFMT_DT
            ThreadPool->SetMask(DEM,FloodDT,FloodDTR,FloodDTC);

            // order 2 calculations, new q u and v
            flood_flowcompute = std::bind((&TWorld::fullSWOF2Do2lightWrapperDynamic1),this,std::placeholders::_1,h,u,v,hs,us,vs, dt1);
            ThreadPool->RunDynamicCompute(flood_flowcompute);
            ThreadPool->WaitForAll();
            //maincalcscheme with dt1 or floodDT->Drc

            // do a second time and do Heun
            if (SwitchHeun)
            {
                // do it all a second time, but not erosion
            //    setFloodMask(h);
                ThreadPool->SetMask(DEM,FloodHMaskDer,FloodHR,FloodHC);
                flood_cellcompute = std::bind((&TWorld::fullSWOF2Do2lightWrapperCell1),this,std::placeholders::_1,hs,us,vs,z);
                ThreadPool->RunDynamicCompute(flood_cellcompute);
                ThreadPool->WaitForAll();

                setFloodDT(h);
                dt1 = Flood_DTMIN;

                FOR_ROW_COL_MV {
                    tma->Drc = 0;
                    FloodT->Drc += FloodDT->Drc;
                    if (FloodT->Drc > 0 && FloodT->Drc < _dt) {
                        tma->Drc = 1;
                    }
                }
                setFloodMaskDT(tma);//FloodDT);
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
           // stop = timesum  > _dt-1e-6;

            double cnt=0;
            FOR_ROW_COL_MV {
                tma->Drc = 0;
                if (FloodT->Drc > 0 && FloodT->Drc < _dt) {
                    cnt+=1.0;
                    tma->Drc = 1;
                }
            }
            stop = cnt < 1;

            if (!stop)
                setFloodMask(tma);

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

    FOR_ROW_COL_MV {
        FloodDTr->Drc = _dt/std::max(1.0,tmb->Drc);
    }

    iter_n = n;
    dt1 = n > 0? _dt/n : dt1;
    return(Flood_DTMIN);//dt1);
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
