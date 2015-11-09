

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
  \file lisSWOF2DSediment.cpp
  \brief Sediment transport for the SWOF2D shallow flood model

functions: \n

- void TWorld::FloodFlowDetachment(void);
- void TWorld::FloodSedFluxReconstruction(void);

    void FloodFlowDetachment();
    void
 */

#include "model.h"
#include "operation.h"

#define signf(x)  ((x < 0)? -1.0 : 1.0)

#define he_ca 1e-12
#define ve_ca 1e-12

#define dt_ca 0.005

#define GRAV 9.8067
#define EPSILON 1e-6

void TWorld::FS_FluxWS(int l, cTMap * _sbl,cTMap * _sss,cTMap * _h1d,cTMap * _h1g,cTMap * _h2d,cTMap * _h2g,cTMap * _u1r,cTMap * _u1l,cTMap * _v1r,cTMap * _v1l,cTMap * _u2r,cTMap * _u2l,cTMap * _v2r,cTMap * _v2l)
{

    FOR_WATERSHED_ROW_COL(l) {
        if (c > 0 && !pcr::isMV(LDD->data[r][c-1]))
    {
      bl1d->data[r][c-1] = std::max(0.0, bl1r->data[r][c-1] );
      bl1g->Drc          = std::max(0.0, bl1l->Drc  );

      ss1d->data[r][c-1] = std::max(0.0, ss1r->data[r][c-1] );
      ss1g->Drc          = std::max(0.0, ss1l->Drc  );

      if (F_scheme == 1)
        FS_Rusanov(_h2d->data[r][c-1], bl1d->data[r][c-1], ss1d->data[r][c-1], _u1r->data[r][c-1], _v1r->data[r][c-1],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);
      else
        if (F_scheme == 2)
          FS_HLL(_h2d->data[r][c-1], bl1d->data[r][c-1], ss1d->data[r][c-1], _u1r->data[r][c-1], _v1r->data[r][c-1],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);
        else
            FS_HLL2(_h2d->data[r][c-1], bl1d->data[r][c-1], ss1d->data[r][c-1], _u1r->data[r][c-1], _v1r->data[r][c-1],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);

      blf1->Drc = HLL2_f1;
      ssf1->Drc = HLL2_f2;

    }
    else
    {
      bl1d->data[r][c] = std::max(0.0, bl1r->data[r][c] );
      bl1g->Drc        = std::max(0.0, bl1l->Drc);

      ss1d->data[r][c] = std::max(0.0, ss1r->data[r][c] );
      ss1g->Drc        = std::max(0.0, ss1l->Drc);

      if (F_scheme == 1)
        FS_Rusanov(_h1d->data[r][c],bl1d->data[r][c],ss1d->data[r][c], _u1r->data[r][c], _v1r->data[r][c],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);
      else
        if (F_scheme == 2)
          FS_HLL(_h1d->data[r][c],bl1d->data[r][c],ss1d->data[r][c], _u1r->data[r][c], _v1r->data[r][c],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);
        else
            FS_HLL2(_h1d->data[r][c],bl1d->data[r][c],ss1d->data[r][c], _u1r->data[r][c], _v1r->data[r][c],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);


      blf1->Drc = HLL2_f1;
      ssf1->Drc = HLL2_f2;

    }}}

    FOR_WATERSHED_ROW_COL(l) {
    if(r > 0 && !pcr::isMV(LDD->data[r-1][c]))
    {

      bl2d->data[r-1][c] = std::max(0.0, bl2r->data[r-1][c]);
      bl2g->Drc          = std::max(0.0, bl2l->Drc);

      ss2d->data[r-1][c] = std::max(0.0, ss2r->data[r-1][c]);
      ss2g->Drc          = std::max(0.0, ss2l->Drc);

      if (F_scheme == 1)
        FS_Rusanov(_h2d->data[r-1][c],bl2d->data[r-1][c],ss2d->data[r-1][c],_v2r->data[r-1][c],_u2r->data[r-1][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);
      else
        if (F_scheme == 2)
          FS_HLL(_h2d->data[r-1][c],bl2d->data[r-1][c],ss2d->data[r-1][c],_v2r->data[r-1][c],_u2r->data[r-1][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);
        else
            FS_HLL2(_h2d->data[r-1][c],bl2d->data[r-1][c],ss2d->data[r-1][c],_v2r->data[r-1][c],_u2r->data[r-1][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);


      blg1->Drc = HLL2_f1;
      ssg1->Drc = HLL2_f2;

    }
    else
    {
       bl2d->data[r][c] = std::max(0.0, bl2r->data[r][c] );
       bl2g->Drc        = std::max(0.0, bl2l->Drc);

       ss2d->data[r][c] = std::max(0.0, ss2r->data[r][c] );
       ss2g->Drc        = std::max(0.0, ss2l->Drc);

       if (F_scheme == 1)
         FS_Rusanov(_h2d->data[r][c],bl2d->data[r][c],ss2d->data[r][c],_v2r->data[r][c],_u2r->data[r][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);
       else
         if (F_scheme == 2)
           FS_HLL(_h2d->data[r][c],bl2d->data[r][c],ss2d->data[r][c],_v2r->data[r][c],_u2r->data[r][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);
         else
             FS_HLL2(_h2d->data[r][c],bl2d->data[r][c],ss2d->data[r][c],_v2r->data[r][c],_u2r->data[r][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);

       blg1->Drc = HLL2_f1;
       ssg1->Drc = HLL2_f2;

    }}}
}
void TWorld::FS_MUSCLEWS(int l, cTMap * _sbl,cTMap * _sss)
{
    double delta_sbl1, delta_sss1;
    double delta_sbl2, delta_sss2;
    double dsbl, dsss;

    // fill EW and NS conditions with cell itself, 1st order approximation used for boundary conditions
    FOR_WATERSHED_ROW_COL(l) {

        bl1r->Drc = _sss->Drc;
        bl1l->Drc = _sss->Drc;

        bl2r->Drc = _sss->Drc;
        bl2l->Drc = _sss->Drc;

        ss1r->Drc = _sss->Drc;
        ss1l->Drc = _sss->Drc;

        ss2r->Drc = _sss->Drc;
        ss2l->Drc = _sss->Drc;

    }}


    FOR_WATERSHED_ROW_COL(l) {
      if(c > 0 && c < _nrCols-1 && !MV(r,c-1) &&  !MV(r,c+1))
    {
        delta_sbl1 = _sbl->Drc - _sbl->data[r][c-1];

        delta_sbl2 = _sbl->data[r][c+1] - _sbl->Drc;

        dsbl   = 0.5*limiter(delta_sbl1, delta_sbl2);


        bl1r->Drc = _sbl->Drc+dsbl;
        bl1l->Drc = _sbl->Drc-dsbl;

        delta_sss1 = _sss->Drc - _sss->data[r][c-1];

        delta_sss2 = _sss->data[r][c+1] - _sss->Drc;

        dsss   = 0.5*limiter(delta_sss1, delta_sss2);


        ss1r->Drc = _sss->Drc+dsss;
        ss1l->Drc = _sss->Drc-dsss;
    }}}


      FOR_WATERSHED_ROW_COL(l) {
        if(r > 0 && MV(r-1,c) && r < _nrRows-1 && !MV(r+1, c))
      {
          delta_sbl1 = _sbl->Drc - _sbl->data[r-1][c];

          delta_sbl2 = _sbl->data[r+1][c] - _sbl->Drc;

          dsbl   = 0.5*limiter(delta_sbl1, delta_sbl2);
          bl2r->Drc = _sbl->Drc+dsbl;
          bl2l->Drc = _sbl->Drc-dsbl;

          delta_sss1 = _sss->Drc - _sss->data[r-1][c];

          delta_sss2 = _sss->data[r+1][c] - _sss->Drc;

          dsss   = 0.5*limiter(delta_sss1, delta_sss2);
          ss2r->Drc = _sss->Drc+dsss;
          ss2l->Drc = _sss->Drc-dsss;
      }}}
}
void TWorld::FS_ENOWS(int l, cTMap * bl,cTMap * ss)
{
    double ddbl1 = 0;
    double ddss1 = 0;
    double ddbl2, ddss2;
    double bl1,ss1;
    double bl2, ss2;
    double delta_bl1, delta_ss1;
    double delta_bl2, delta_ss2;
    double dbl, dss;
    double amortENO = 0.25;
    double MODIFENO = 0.9;
    double a1, a2;

  //x-direction
    FOR_WATERSHED_ROW_COL(l) {
      if(c > 0 && c < _nrCols-2 && !MV(r,c-1) && !MV(r, c+1) && !MV(r, c+2))
    {
      bl1 = bl->data[r][c-1]-2.*bl->data[r][c]+bl->data[r][c+1];
      ss1 = ss->data[r][c-1]-2.*ss->data[r][c]+ss->data[r][c+1];

      bl2 = bl->Drc-2.*bl->data[r][c+1]+bl->data[r][c+2];
      ss2 = ss->Drc-2.*ss->data[r][c+1]+ss->data[r][c+2];

      ddbl2 = amortENO*limiter(bl1,bl2);
      ddss2 = amortENO*limiter(ss1,ss2);

      delta_bl1 = bl->Drc - bl->data[r][c-1];
      delta_ss1 = ss->Drc - ss->data[r][c-1];

      delta_bl2 = bl->data[r][c+1]-bl->Drc;
      delta_ss2 = ss->data[r][c+1]-ss->Drc;

      if (F_scheme == (int)FENO)
      {
          dbl = limiter(delta_bl1+ddbl1*0.5,delta_bl2-ddbl2*0.5);
          dss = limiter(delta_ss1+ddss1*0.5,delta_ss2-ddss2*0.5);
      }
      else
      {
          a2 = limiter(delta_bl1+ddbl1*0.5,delta_bl2-ddbl2*0.5);
          a1 = limiter(delta_bl1,delta_bl2);
          dbl = limiter(2*MODIFENO*a1,a2);
          a2 = limiter(delta_ss1+ddss1*0.5,delta_ss2-ddss2*0.5);
          a1 = limiter(delta_ss1,delta_ss2);
          dss = limiter(2*MODIFENO*a1,a2);
      }

      bl1r->Drc=bl->Drc+dbl*0.5;
      bl1l->Drc=bl->Drc-dbl*0.5;
      ss1r->Drc=ss->Drc+dss*0.5;
      ss1l->Drc=ss->Drc-dss*0.5;


      ddbl1=ddbl2;
      ddss1=ddss2;

    }}}


    //y-direction
    FOR_WATERSHED_ROW_COL(l) {
    if(r > 0 && r < _nrRows-2 && !MV(r-1,c) && !MV(r+1, c) && !MV(r+2, c))
    {
        bl1 = bl->data[r-1][c]-2.*bl->data[r][c]+bl->data[r+1][c];
        ss1 = ss->data[r-1][c]-2.*ss->data[r][c]+ss->data[r+1][c];

        bl2 = bl->Drc-2.*bl->data[r+1][c]+bl->data[r+2][c];
        ss2 = ss->Drc-2.*ss->data[r+1][c]+ss->data[r+2][c];

        ddbl2 = amortENO*limiter(bl1,bl2);
        ddss2 = amortENO*limiter(ss1,ss2);

        delta_bl1 = bl->Drc - bl->data[r-1][c];
        delta_ss1 = ss->Drc - ss->data[r-1][c];

        delta_bl2 = bl->data[r+1][c]-bl->Drc;
        delta_ss2 = ss->data[r+1][c]-ss->Drc;

        if (F_scheme == (int)FENO)
        {
            dbl = limiter(delta_bl1+ddbl1*0.5,delta_bl2-ddbl2*0.5);
            dss = limiter(delta_ss1+ddss1*0.5,delta_ss2-ddss2*0.5);
        }
        else
        {
            a2 = limiter(delta_bl1+ddbl1*0.5,delta_bl2-ddbl2*0.5);
            a1 = limiter(delta_bl1,delta_bl2);
            dbl = limiter(2*MODIFENO*a1,a2);
            a2 = limiter(delta_ss1+ddss1*0.5,delta_ss2-ddss2*0.5);
            a1 = limiter(delta_ss1,delta_ss2);
            dss = limiter(2*MODIFENO*a1,a2);

        }

        bl2r->Drc = bl->Drc+dbl*0.5;
        bl2l->Drc = bl->Drc-dbl*0.5;
        ss2r->Drc = ss->Drc+dss*0.5;
        ss2l->Drc = ss->Drc-dss*0.5;

        ddbl1=ddbl2;
        ddss1=ddss2;

    }}}


}
void TWorld::FS_SimpleWS(int l, cTMap * _sbl,cTMap * _sss)
{
    FOR_WATERSHED_ROW_COL(l) {
        bl1r->Drc = _sss->Drc;
        bl1l->Drc = _sss->Drc;

        bl2r->Drc = _sss->Drc;
        bl2l->Drc = _sss->Drc;

        ss1r->Drc = _sss->Drc;
        ss1l->Drc = _sss->Drc;

        ss2r->Drc = _sss->Drc;
        ss2l->Drc = _sss->Drc;
    }}

}
void TWorld::FS_MainCalcWS(int l, cTMap * _h, cTMap * _sbl,cTMap * _sbln,cTMap * _sss,cTMap * _sssn, double dt)
{
    FOR_WATERSHED_ROW_COL(l) {

      double dx = ChannelAdj->Drc;
      double dy = DX->Drc;
      long double tx = dt/dx;
      long double ty = dt/dy;
      double hestemp = _h->Drc;
      // Solution of the equation of mass conservation (First equation of Saint venant)
      // f1 comes from MUSCL calculations
      if ((r > _nrRows-2 || c > _nrCols-2) || (pcr::isMV(LDD->data[r][c+1]) || pcr::isMV(LDD->data[r+1][c])))
      {
        _sbln->Drc = _sbl->Drc;
        _sssn->Drc = _sss->Drc;
      }
      else
      {
        _sbln->Drc = _sbl->Drc - tx*(blf1->data[r][c+1]-blf1->Drc) - ty*(blg1->data[r+1][c]-blg1->Drc);
        _sssn->Drc = _sss->Drc - tx*(ssf1->data[r][c+1]-ssf1->Drc) - ty*(ssg1->data[r+1][c]-ssg1->Drc);
      }
    }}
}

void TWorld::FS_Flux(cTMap * _sbl,cTMap * _sss,cTMap * _h1d,cTMap * _h1g,cTMap * _h2d,cTMap * _h2g,cTMap * _u1r,cTMap * _u1l,cTMap * _v1r,cTMap * _v1l,cTMap * _u2r,cTMap * _u2l,cTMap * _v2r,cTMap * _v2l)
{


    FOR_CELL_IN_FLOODAREA
        if (c > 0 && !pcr::isMV(LDD->data[r][c-1]))
    {
      bl1d->data[r][c-1] = std::max(0.0, bl1r->data[r][c-1] );
      bl1g->Drc          = std::max(0.0, bl1l->Drc  );

      ss1d->data[r][c-1] = std::max(0.0, ss1r->data[r][c-1] );
      ss1g->Drc          = std::max(0.0, ss1l->Drc  );

      if (F_scheme == 1)
        FS_Rusanov(_h2d->data[r][c-1], bl1d->data[r][c-1], ss1d->data[r][c-1], _u1r->data[r][c-1], _v1r->data[r][c-1],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);
      else
        if (F_scheme == 2)
          FS_HLL(_h2d->data[r][c-1], bl1d->data[r][c-1], ss1d->data[r][c-1], _u1r->data[r][c-1], _v1r->data[r][c-1],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);
        else
            FS_HLL2(_h2d->data[r][c-1], bl1d->data[r][c-1], ss1d->data[r][c-1], _u1r->data[r][c-1], _v1r->data[r][c-1],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);

      blf1->Drc = HLL2_f1;
      ssf1->Drc = HLL2_f2;

    }
    else
    {
      bl1d->data[r][c] = std::max(0.0, bl1r->data[r][c] );
      bl1g->Drc        = std::max(0.0, bl1l->Drc);

      ss1d->data[r][c] = std::max(0.0, ss1r->data[r][c] );
      ss1g->Drc        = std::max(0.0, ss1l->Drc);

      if (F_scheme == 1)
        FS_Rusanov(_h1d->data[r][c],bl1d->data[r][c],ss1d->data[r][c], _u1r->data[r][c], _v1r->data[r][c],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);
      else
        if (F_scheme == 2)
          FS_HLL(_h1d->data[r][c],bl1d->data[r][c],ss1d->data[r][c], _u1r->data[r][c], _v1r->data[r][c],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);
        else
            FS_HLL2(_h1d->data[r][c],bl1d->data[r][c],ss1d->data[r][c], _u1r->data[r][c], _v1r->data[r][c],_h1g->Drc,bl1g->Drc,ss1g->Drc, _u1l->Drc, _v1l->Drc);

      blf1->Drc = HLL2_f1;
      ssf1->Drc = HLL2_f2;

    }}

    FOR_CELL_IN_FLOODAREA
    if(r > 0 && !pcr::isMV(LDD->data[r-1][c]))
    {

      bl2d->data[r-1][c] = std::max(0.0, bl2r->data[r-1][c]);
      bl2g->Drc          = std::max(0.0, bl2l->Drc);

      ss2d->data[r-1][c] = std::max(0.0, ss2r->data[r-1][c]);
      ss2g->Drc          = std::max(0.0, ss2l->Drc);

      if (F_scheme == 1)
        FS_Rusanov(_h2d->data[r-1][c],bl2d->data[r-1][c],ss2d->data[r-1][c],_v2r->data[r-1][c],_u2r->data[r-1][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);
      else
        if (F_scheme == 2)
          FS_HLL(_h2d->data[r-1][c],bl2d->data[r-1][c],ss2d->data[r-1][c],_v2r->data[r-1][c],_u2r->data[r-1][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);
        else
            FS_HLL2(_h2d->data[r-1][c],bl2d->data[r-1][c],ss2d->data[r-1][c],_v2r->data[r-1][c],_u2r->data[r-1][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);

      blg1->Drc = HLL2_f1;
      ssg1->Drc = HLL2_f2;

    }
    else
    {
       bl2d->data[r][c] = std::max(0.0, bl2r->data[r][c] );
       bl2g->Drc        = std::max(0.0, bl2l->Drc);

       ss2d->data[r][c] = std::max(0.0, ss2r->data[r][c] );
       ss2g->Drc        = std::max(0.0, ss2l->Drc);

       if (F_scheme == 1)
         FS_Rusanov(_h2d->data[r][c],bl2d->data[r][c],ss2d->data[r][c],_v2r->data[r][c],_u2r->data[r][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);
       else
         if (F_scheme == 2)
           FS_HLL(_h2d->data[r][c],bl2d->data[r][c],ss2d->data[r][c],_v2r->data[r][c],_u2r->data[r][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);
         else
             FS_HLL2(_h2d->data[r][c],bl2d->data[r][c],ss2d->data[r][c],_v2r->data[r][c],_u2r->data[r][c],_h2g->Drc,bl2g->Drc,ss2g->Drc,_v2l->Drc,_u2l->Drc);

       blg1->Drc = HLL2_f1;
       ssg1->Drc = HLL2_f2;

    }}

}


void TWorld::FS_ENO(cTMap * bl,cTMap * ss)
{
    double ddbl1 = 0;
    double ddss1 = 0;
    double ddbl2, ddss2;
    double bl1,ss1;
    double bl2, ss2;
    double delta_bl1, delta_ss1;
    double delta_bl2, delta_ss2;
    double dbl, dss;
    double amortENO = 0.25;
    double MODIFENO = 0.9;
    double a1, a2;

   //x-direction
    FOR_CELL_IN_FLOODAREA
    if(c > 0 && c < _nrCols-2 && !MV(r,c-1) && !MV(r, c+1) && !MV(r, c+2))
    {
      bl1 = bl->data[r][c-1]-2.*bl->data[r][c]+bl->data[r][c+1];
      ss1 = ss->data[r][c-1]-2.*ss->data[r][c]+ss->data[r][c+1];

      bl2 = bl->Drc-2.*bl->data[r][c+1]+bl->data[r][c+2];
      ss2 = ss->Drc-2.*ss->data[r][c+1]+ss->data[r][c+2];

      ddbl2 = amortENO*limiter(bl1,bl2);
      ddss2 = amortENO*limiter(ss1,ss2);

      delta_bl1 = bl->Drc - bl->data[r][c-1];
      delta_ss1 = ss->Drc - ss->data[r][c-1];

      delta_bl2 = bl->data[r][c+1]-bl->Drc;
      delta_ss2 = ss->data[r][c+1]-ss->Drc;

      if (F_scheme == (int)FENO)
      {
          dbl = limiter(delta_bl1+ddbl1*0.5,delta_bl2-ddbl2*0.5);
          dss = limiter(delta_ss1+ddss1*0.5,delta_ss2-ddss2*0.5);
      }
      else
      {
          a2 = limiter(delta_bl1+ddbl1*0.5,delta_bl2-ddbl2*0.5);
          a1 = limiter(delta_bl1,delta_bl2);
          dbl = limiter(2*MODIFENO*a1,a2);
          a2 = limiter(delta_ss1+ddss1*0.5,delta_ss2-ddss2*0.5);
          a1 = limiter(delta_ss1,delta_ss2);
          dss = limiter(2*MODIFENO*a1,a2);
      }

      bl1r->Drc=bl->Drc+dbl*0.5;
      bl1l->Drc=bl->Drc-dbl*0.5;
      ss1r->Drc=ss->Drc+dss*0.5;
      ss1l->Drc=ss->Drc-dss*0.5;


      ddbl1=ddbl2;
      ddss1=ddss2;

    }}


    //y-direction
    FOR_CELL_IN_FLOODAREA
    if(r > 0 && r < _nrRows-2 && !MV(r-1,c) && !MV(r+1, c) && !MV(r+2, c))
    {
        bl1 = bl->data[r-1][c]-2.*bl->data[r][c]+bl->data[r+1][c];
        ss1 = ss->data[r-1][c]-2.*ss->data[r][c]+ss->data[r+1][c];

        bl2 = bl->Drc-2.*bl->data[r+1][c]+bl->data[r+2][c];
        ss2 = ss->Drc-2.*ss->data[r+1][c]+ss->data[r+2][c];

        ddbl2 = amortENO*limiter(bl1,bl2);
        ddss2 = amortENO*limiter(ss1,ss2);

        delta_bl1 = bl->Drc - bl->data[r-1][c];
        delta_ss1 = ss->Drc - ss->data[r-1][c];

        delta_bl2 = bl->data[r+1][c]-bl->Drc;
        delta_ss2 = ss->data[r+1][c]-ss->Drc;

        if (F_scheme == (int)FENO)
        {
            dbl = limiter(delta_bl1+ddbl1*0.5,delta_bl2-ddbl2*0.5);
            dss = limiter(delta_ss1+ddss1*0.5,delta_ss2-ddss2*0.5);
        }
        else
        {
            a2 = limiter(delta_bl1+ddbl1*0.5,delta_bl2-ddbl2*0.5);
            a1 = limiter(delta_bl1,delta_bl2);
            dbl = limiter(2*MODIFENO*a1,a2);
            a2 = limiter(delta_ss1+ddss1*0.5,delta_ss2-ddss2*0.5);
            a1 = limiter(delta_ss1,delta_ss2);
            dss = limiter(2*MODIFENO*a1,a2);

        }

        bl2r->Drc = bl->Drc+dbl*0.5;
        bl2l->Drc = bl->Drc-dbl*0.5;
        ss2r->Drc = ss->Drc+dss*0.5;
        ss2l->Drc = ss->Drc-dss*0.5;

        ddbl1=ddbl2;
        ddss1=ddss2;

    }}

}

void TWorld::FS_MUSCLE(cTMap * _sbl,cTMap * _sss)
{

    double delta_sbl1, delta_sss1;
    double delta_sbl2, delta_sss2;
    double dsbl, dsss;

    // fill EW and NS conditions with cell itself, 1st order approximation used for boundary conditions
    FOR_CELL_IN_FLOODAREA {

        bl1r->Drc = _sss->Drc;
        bl1l->Drc = _sss->Drc;

        bl2r->Drc = _sss->Drc;
        bl2l->Drc = _sss->Drc;

        ss1r->Drc = _sss->Drc;
        ss1l->Drc = _sss->Drc;

        ss2r->Drc = _sss->Drc;
        ss2l->Drc = _sss->Drc;

    }}


    FOR_CELL_IN_FLOODAREA
      if(c > 0 && c < _nrCols-1 && !MV(r,c-1) &&  !MV(r,c+1))
    {
        delta_sbl1 = _sbl->Drc - _sbl->data[r][c-1];

        delta_sbl2 = _sbl->data[r][c+1] - _sbl->Drc;

        dsbl   = 0.5*limiter(delta_sbl1, delta_sbl2);


        bl1r->Drc = _sbl->Drc+dsbl;
        bl1l->Drc = _sbl->Drc-dsbl;

        delta_sss1 = _sss->Drc - _sss->data[r][c-1];

        delta_sss2 = _sss->data[r][c+1] - _sss->Drc;

        dsss   = 0.5*limiter(delta_sss1, delta_sss2);


        ss1r->Drc = _sss->Drc+dsss;
        ss1l->Drc = _sss->Drc-dsss;
    }}


      FOR_CELL_IN_FLOODAREA
        if(r > 0 && MV(r-1,c) && r < _nrRows-1 && !MV(r+1, c))
      {
          delta_sbl1 = _sbl->Drc - _sbl->data[r-1][c];

          delta_sbl2 = _sbl->data[r+1][c] - _sbl->Drc;

          dsbl   = 0.5*limiter(delta_sbl1, delta_sbl2);
          bl2r->Drc = _sbl->Drc+dsbl;
          bl2l->Drc = _sbl->Drc-dsbl;

          delta_sss1 = _sss->Drc - _sss->data[r-1][c];

          delta_sss2 = _sss->data[r+1][c] - _sss->Drc;

          dsss   = 0.5*limiter(delta_sss1, delta_sss2);
          ss2r->Drc = _sss->Drc+dsss;
          ss2l->Drc = _sss->Drc-dsss;
      }}



}

void TWorld::FS_Simple(cTMap * _sbl,cTMap * _sss)
{
    FOR_CELL_IN_FLOODAREA {
        bl1r->Drc = _sss->Drc;
        bl1l->Drc = _sss->Drc;

        bl2r->Drc = _sss->Drc;
        bl2l->Drc = _sss->Drc;

        ss1r->Drc = _sss->Drc;
        ss1l->Drc = _sss->Drc;

        ss2r->Drc = _sss->Drc;
        ss2l->Drc = _sss->Drc;
    }}
}


void TWorld::FS_MainCalc(cTMap * _h,cTMap * _sbl, cTMap * _sbln,cTMap * _sss,cTMap * _sssn, double dt)
{
    FOR_CELL_IN_FLOODAREA
    {
      double dx = ChannelAdj->Drc;
      double dy = DX->Drc;
      long double tx = dt/dx;
      long double ty = dt/dy;
      double hestemp = _h->Drc;
      // Solution of the equation of mass conservation (First equation of Saint venant)
      // f1 comes from MUSCL calculations
      if ((r > _nrRows-2 || c > _nrCols-2) || (pcr::isMV(LDD->data[r][c+1]) || pcr::isMV(LDD->data[r+1][c])))
      {
        _sbln->Drc = _sbl->Drc;
        _sssn->Drc = _sss->Drc;
      }
      else
      {
        _sbln->Drc = _sbl->Drc - tx*(blf1->data[r][c+1]-blf1->Drc) - ty*(blg1->data[r+1][c]-blg1->Drc);
        _sssn->Drc = _sss->Drc - tx*(ssf1->data[r][c+1]-ssf1->Drc) - ty*(ssg1->data[r+1][c]-ssg1->Drc);
      }
    }}
}

void TWorld::FS_Rusanov(double h_L,double bl_L,double ss_L,double u_L,double v_L,double h_R, double bl_R,double ss_R,double u_R,double v_R)
{

    double f1, f2, f3, cfl, tmp = 0;
    double c;
    if (h_L<=0. && h_R<=0.){
        f1 = 0.;
        f2 = 0.;
    }else{
        c = std::max(fabs(u_L)+sqrt(GRAV*h_L),fabs(u_R)+sqrt(GRAV*h_R));
        double cd = c*0.5;
        double qss_R = ss_R*u_R*h_R;
        double qss_L = ss_L*u_L*h_L;
        double qbl_R = bl_R*u_R*h_R;
        double qbl_L = bl_L*u_L*h_L;
        f1 = (qbl_L+qbl_R)*0.5-cd*(h_R-h_L);
        f2 = (qss_L+qss_R)*0.5-cd*(h_R-h_L);

    }
    HLL2_f1 = f1;
    HLL2_f2 = f2;
}

void TWorld::FS_HLL2(double h_L,double bl_L,double ss_L,double u_L,double v_L,double h_R, double bl_R,double ss_R,double u_R,double v_R)
{

    double f1, f2, f3, cfl, tmp = 0;
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
        double sqrt_grav_h_L = sqrt(grav_h_L);  // wave velocity
        double sqrt_grav_h_R = sqrt(grav_h_R);
        double qss_R = ss_R*u_R*h_R;
        double qss_L = ss_L*u_L*h_L;
        double qbl_R = bl_R*u_R*h_R;
        double qbl_L = bl_L*u_L*h_L;

        double c1 = std::min(u_L - sqrt_grav_h_L,u_R - sqrt_grav_h_R); //we already have u_L - sqrt_grav_h_L<u_L + sqrt_grav_h_L and u_R - sqrt_grav_h_R<u_R + sqrt_grav_h_R
        double c2 = std::max(u_L + sqrt_grav_h_L,u_R + sqrt_grav_h_R); //so we do not need all the eigenvalues to get c1 and c2
        tmp = 1./(c2-c1);
        double t1 = (std::min(c2,0.) - std::min(c1,0.))*tmp;
        double t2 = 1. - t1;
        double t3 = (c2*fabs(c1) - c1*fabs(c2))*0.5*tmp;

        f1 = t1*qbl_R + t2*qbl_L - t3*(h_R - h_L);
        f2 = t1*qss_R + t2*qss_L - t3*(h_R - h_L);

    }
    HLL2_f1 = f1;
    HLL2_f2 = f2;
}

void TWorld::FS_HLL(double h_L,double bl_L,double ss_L,double u_L,double v_L,double h_R, double bl_R,double ss_R,double u_R,double v_R)
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
        }else if (c1>=EPSILON){ //supercritical flow, from left to right : we have std::max(abs(c1),abs(c2))=c2>0
            f1=bl_L*(q_L);
            f2=ss_L*(q_L);
        }else if (c2<=-EPSILON){ //supercritical flow, from right to left : we have std::max(abs(c1),abs(c2))=-c1>0
            f1=bl_R*(q_R);
            f2=ss_R*(q_R);
        }else{ //subcritical flow
            tmp = 1./(c2-c1);
            f1=(c2*bl_L*q_L-c1*bl_R*q_R)*tmp+c1*c2*(h_R-h_L)*tmp;

            f2=(c2*ss_L*q_L-c1*ss_R*q_R)*tmp+c1*c2*(h_R-h_L)*tmp;
        }
    }

    HLL2_f1 = f1;
    HLL2_f2 = f2;


}

void TWorld::SWOFSedimentFlowWS(int l, double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
{
    cTMap*_h = hmx;

    FOR_WATERSHED_ROW_COL(l) {

        MBLNFlood->Drc = _BL->Drc;
        MSSNFlood->Drc = _SS->Drc;
        MBLFlood->Drc = _BL->Drc;
        MSSFlood->Drc = _SS->Drc;

        //set concentration from present sediment
        MBLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MBLFlood->Drc);

        //set concentration from present sediment
        MSSCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MSSFlood->Drc);

        //set concentration from present sediment
        MBLCNFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MBLFlood->Drc);

        //set concentration from present sediment
        MSSCNFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MSSFlood->Drc);

    }}


    if (SwitchFloodSWOForder2)
    {

        //reconstruction scheme

        if (F_scheme == (int)FMUSCL)
        {
            FS_MUSCLEWS(l,MBLCFlood,MSSCFlood);
        }else
        {
            FS_ENOWS(l,MBLCFlood,MSSCFlood);
        }

        //flux calculation
        FS_FluxWS(l,MBLCFlood,MSSCFlood,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12);

        //Calculate new Sediment
        FS_MainCalcWS(l,_h,MBLFlood,bls,MSSFlood,sss,dt);

        //calculate new concentration with approximated new water height and approximated new sediment
        FOR_WATERSHED_ROW_COL(l) {

            //set concentration from present sediment
            MBLCNFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hs->Drc, bls->Drc);

            //set concentration from present sediment
            MSSCNFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hs->Drc, sss->Drc);
        }}


        //redo for heun's method
        if (F_scheme == (int)FMUSCL)
        {
            FS_MUSCLEWS(l,MBLCNFlood,MSSCNFlood);
        }else
        {
            FS_ENOWS(l,MBLCNFlood,MSSCNFlood);
        }

        //flux calculation
        FS_FluxWS(l,MBLCNFlood,MSSCNFlood,h1d,h1g,h2d,h2g,u1r,u1l,v1r,v1l,u2r,u2l,v2r,v2l);

        //Calculate new Sediment
        FS_MainCalcWS(l,hs,bls,MBLNFlood,sss,MSSNFlood,dt);

        //update variable according to heuns method
        FOR_WATERSHED_ROW_COL(l) {
          MBLNFlood->Drc = (MBLNFlood->Drc + MBLFlood->Drc)/2.0;
          MSSNFlood->Drc = (MSSNFlood->Drc + MSSFlood->Drc)/2.0 ;
        }}

    }else
    {
        //reconstruction scheme

        FS_Simple(BLCFlood,SSCFlood);

       /* if (F_scheme == (int)FMUSCL)
        {
            FS_MUSCLEWS(l,BLCFlood,SSCFlood);
        }else
        {
            FS_ENOWS(l,BLCFlood,SSCFlood);
        }*/

        //flux calculation
        FS_FluxWS(l,MBLCFlood,MSSCFlood,h1d,h1g,h2d,h2g,u1r,u1l,v1r,v1l,u2r,u2l,v2r,v2l);

        //Calculate new Sediment
        FS_MainCalcWS(l,_h,MBLFlood,bls,MSSFlood,sss,dt);

        //update variable
        FOR_WATERSHED_ROW_COL(l) {
          MBLNFlood->Drc = bls->Drc;
          MSSNFlood->Drc = sss->Drc;
        }}
    }


    FOR_WATERSHED_ROW_COL(l) {

        _BL->Drc = MBLNFlood->Drc;
        _SS->Drc = MSSNFlood->Drc;

        //set concentration from present sediment
        _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _BL->Drc);

        //set concentration from present sediment
        _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _SS->Drc);

    }}



}


void TWorld::SWOFSedimentFlow(double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
{

    cTMap*_h = hmx;

    FOR_CELL_IN_FLOODAREA

        MBLNFlood->Drc = _BL->Drc;
        MSSNFlood->Drc = _SS->Drc;
        MBLFlood->Drc = _BL->Drc;
        MSSFlood->Drc = _SS->Drc;

        //set concentration from present sediment
        MBLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MBLFlood->Drc);

        //set concentration from present sediment
        MSSCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MSSFlood->Drc);

        //set concentration from present sediment
        MBLCNFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MBLFlood->Drc);

        //set concentration from present sediment
        MSSCNFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MSSFlood->Drc);

    }


    if (SwitchFloodSWOForder2)
    {

        //reconstruction scheme

        if (F_scheme == (int)FMUSCL)
        {
            FS_MUSCLE(MBLCFlood,MSSCFlood);
        }else
        {
            FS_ENO(MBLCFlood,MSSCFlood);
        }

        //flux calculation
        FS_Flux(MBLCFlood,MSSCFlood,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12);

        //Calculate new Sediment
        FS_MainCalc(_h,MBLFlood,bls,MSSFlood,sss,dt);

        //calculate new concentration with approximated new water height and approximated new sediment
        FOR_CELL_IN_FLOODAREA

            //set concentration from present sediment
            MBLCNFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hs->Drc, bls->Drc);

            //set concentration from present sediment
            MSSCNFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hs->Drc, sss->Drc);
        }


        //redo for heun's method
        if (F_scheme == (int)FMUSCL)
        {
            FS_MUSCLE(MBLCNFlood,MSSCNFlood);
        }else
        {
            FS_ENO(MBLCNFlood,MSSCNFlood);
        }

        //flux calculation
        FS_Flux(MBLCNFlood,MSSCNFlood,h1d,h1g,h2d,h2g,u1r,u1l,v1r,v1l,u2r,u2l,v2r,v2l);

        //Calculate new Sediment
        FS_MainCalc(hs,bls,MBLNFlood,sss,MSSNFlood,dt);

        //update variable according to heuns method
        FOR_CELL_IN_FLOODAREA {
          BLFlood->Drc = (MBLNFlood->Drc + MBLFlood->Drc)/2.0;
          SSFlood->Drc = (MSSNFlood->Drc + MSSFlood->Drc)/2.0 ;
        }}

    }else
    {
        //reconstruction scheme

        FS_Simple(MBLCFlood,MSSCFlood);

        /*if (F_scheme == (int)FMUSCL)
        {
            FS_MUSCLE(BLCFlood,SSCFlood);
        }else
        {
            FS_ENO(BLCFlood,SSCFlood);
        }*/

        //flux calculation
        FS_Flux(MBLCFlood,MSSCFlood,h1d,h1g,h2d,h2g,u1r,u1l,v1r,v1l,u2r,u2l,v2r,v2l);

        //Calculate new Sediment
        FS_MainCalc(_h,MBLFlood,bls,MSSFlood,sss,dt);

        //update variable
        FOR_CELL_IN_FLOODAREA {
          MBLNFlood->Drc = bls->Drc;
          MSSNFlood->Drc = sss->Drc;
        }}
    }

    FOR_CELL_IN_FLOODAREA

        _BL->Drc = MBLNFlood->Drc;
        _SS->Drc = MSSNFlood->Drc;

        //set concentration from present sediment
        _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _BL->Drc);

        //set concentration from present sediment
        _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _SS->Drc);

    }


}


void TWorld::SWOFSedimentFlowInterpolationWS(int l, double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
{

    FOR_WATERSHED_ROW_COL(l) {

        MBLNFlood->Drc = _BL->Drc;
        MSSNFlood->Drc = _SS->Drc;
        MBLFlood->Drc = _BL->Drc;
        MSSFlood->Drc = _SS->Drc;

        //set concentration from present sediment
        MBLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MBLFlood->Drc);

        //set concentration from present sediment
        MSSCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MSSFlood->Drc);

    }}

    double courant = this->courant_factor;

    //first calculate the weights for the cells that are closest to location that flow is advected to
    FOR_WATERSHED_ROW_COL(l) {

        //no flood velocity means no flood sediment transport, so skip this cell
        if((Vflood->Drc == 0 && Uflood->Drc == 0))
        {
            continue;
        }


        //the sign of the x and y direction of flow
        double yn = signf(Vflood->Drc);
        double xn = signf(Uflood->Drc);

        double vel = sqrt(Uflood->Drc*Uflood->Drc + Vflood->Drc*Vflood->Drc);

        if(vel == 0 ||hmx->Drc == 0)
        {
            continue;
        }

        double qbl = vel*ChannelAdj->Drc *BLDepthFlood->Drc * MBLCFlood->Drc;

        if(qbl > courant * DX->Drc * ChannelAdj->Drc *BLDepthFlood->Drc * MBLCFlood->Drc)
        {
            qbl =  courant * DX->Drc * ChannelAdj->Drc *BLDepthFlood->Drc * MBLCFlood->Drc;
        }

        double qss = vel*ChannelAdj->Drc *SSDepthFlood->Drc * MSSCFlood->Drc;

        if(qss > courant * DX->Drc * ChannelAdj->Drc *SSDepthFlood->Drc * MSSCFlood->Drc)
        {
            qss =  courant * DX->Drc * ChannelAdj->Drc *SSDepthFlood->Drc * MSSCFlood->Drc;
        }

        //should not travel more distance than cell size
        double dsx = xn*std::min(fabs(Uflood->Drc)/vel,1.0);
        double dsy = yn*std::min(fabs(Vflood->Drc)/vel,1.0);

        //cell directions
        int dx[4] = {0, 1, 1, 0};
        int dy[4] = {1, 0, 1, 0};

        //weights to be saved
        double w[4] = {0.0,0.0,0.0,0.0};

        //for each cell niegbhouring the advected location of the discharge, calculate interpolation weight
        for (int i=0; i<4; i++)
        {
            int r2, c2;

            //must multiply the cell directiosby the sign of the slope vector components
            r2 = r+yn*dy[i];
            c2 = c+xn*dx[i];

            // distance we want is equal to: 1 - distance from the advected location to the neighbouring cell
            double wdx = ((double)1.0) - fabs( xn * ((double)dx[i]) - dsx);
            double wdy = ((double)1.0) - fabs( yn * ((double)dy[i]) - dsy);

            //the distribution is inverly proportional to the squared distance
            double weight = fabs(wdx) *fabs(wdy);

            if(INSIDE(r2,c2))
            {
                if( !pcr::isMV(LDD->data[r2][c2]) && hmx->data[r2][c2] > 0)
                {
                    w[i] = weight;
                }
            }
        }

        //normalize: sum of the 4 weights is equal to 1
        double wt = 0.0;
        for (int i=0; i<4; i++)
        {
            wt += w[i];
        }
        if(wt == 0)
        {
            w[3] = 1.0; wt = 1.0;
        }
        for (int i=0; i<4; i++)
        {
            w[i] = w[i]/wt;
        }



        //use the calculated weights to distribute flow
        for (int i=0; i<4; i++)
        {
            int r2, c2;

            //must multiply the cell directions by the sign of the slope vector components
            r2 = r+yn*dy[i];
            c2 = c+xn*dx[i];

            if(INSIDE(r2,c2) && !pcr::isMV(LDD->data[r2][c2]))
            {

                if(hmx->data[r2][c2] > 0)
                {
                    //weight * the flow is distributed to the ith cell that neighbours the advected flow.
                    MBLNFlood->data[r2][c2] +=  w[i]* dt * qbl;
                    MBLNFlood->data[r][c] -=  w[i]* dt* qbl;
                    MSSNFlood->data[r2][c2] +=  w[i]* dt * qss;
                    MSSNFlood->data[r][c] -=  w[i]* dt* qss;



                }

            }
        }
    }}


    FOR_WATERSHED_ROW_COL(l) {

        _BL->Drc = MBLNFlood->Drc;
        _SS->Drc = MSSNFlood->Drc;

        //set concentration from present sediment
        _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _BL->Drc);

        //set concentration from present sediment
        _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _SS->Drc);

    }}
}

void TWorld::SWOFSedimentFlowInterpolation(double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
{

    FOR_CELL_IN_FLOODAREA

        MBLNFlood->Drc = _BL->Drc;
        MSSNFlood->Drc = _SS->Drc;
        MBLFlood->Drc = _BL->Drc;
        MSSFlood->Drc = _SS->Drc;

        //set concentration from present sediment
        MBLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MBLFlood->Drc);

        //set concentration from present sediment
        MSSCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MSSFlood->Drc);

    }

   double courant = this->courant_factor;

    //first calculate the weights for the cells that are closest to location that flow is advected to
    FOR_CELL_IN_FLOODAREA

        //no flood velocity means no flood sediment transport, so skip this cell
        if((Vflood->Drc == 0 && Uflood->Drc == 0))
        {
            continue;
        }


        //the sign of the x and y direction of flow
        double yn = signf(Vflood->Drc);
        double xn = signf(Uflood->Drc);

        double vel = sqrt(Uflood->Drc*Uflood->Drc + Vflood->Drc*Vflood->Drc);

        if(vel == 0 ||hmx->Drc == 0)
        {
            continue;
        }

        double qbl = vel*ChannelAdj->Drc *BLDepthFlood->Drc * MBLCFlood->Drc;

        if(qbl > courant * DX->Drc * ChannelAdj->Drc *BLDepthFlood->Drc * MBLCFlood->Drc)
        {
            qbl =  courant * DX->Drc * ChannelAdj->Drc *BLDepthFlood->Drc * MBLCFlood->Drc;
        }

        double qss = vel*ChannelAdj->Drc *SSDepthFlood->Drc * MSSCFlood->Drc;

        if(qss > courant * DX->Drc * ChannelAdj->Drc *SSDepthFlood->Drc * MSSCFlood->Drc)
        {
            qss = courant *  DX->Drc * ChannelAdj->Drc *SSDepthFlood->Drc * MSSCFlood->Drc;
        }

        //should not travel more distance than cell size
        double dsx = xn*std::min(fabs(Uflood->Drc)/vel,1.0);
        double dsy = yn*std::min(fabs(Vflood->Drc)/vel,1.0);

        //cell directions
        int dx[4] = {0, 1, 1, 0};
        int dy[4] = {1, 0, 1, 0};

        //weights to be saved
        double w[4] = {0.0,0.0,0.0,0.0};

        //for each cell niegbhouring the advected location of the discharge, calculate interpolation weight
        for (int i=0; i<4; i++)
        {
            int r2, c2;

            //must multiply the cell directiosby the sign of the slope vector components
            r2 = r+yn*dy[i];
            c2 = c+xn*dx[i];

            // distance we want is equal to: 1 - distance from the advected location to the neighbouring cell
            double wdx = ((double)1.0) - fabs( xn * ((double)dx[i]) - dsx);
            double wdy = ((double)1.0) - fabs( yn * ((double)dy[i]) - dsy);

            //the distribution is inverly proportional to the squared distance
            double weight = fabs(wdx) *fabs(wdy);

            if(INSIDE(r2,c2))
            {
                if( !pcr::isMV(LDD->data[r2][c2]) && hmx->data[r2][c2] > 0)
                {
                    w[i] = weight;
                }
            }
        }

        //normalize: sum of the 4 weights is equal to 1
        double wt = 0.0;
        for (int i=0; i<4; i++)
        {
            wt += w[i];
        }
        if(wt == 0)
        {
            w[3] = 1.0; wt = 1.0;
        }
        for (int i=0; i<4; i++)
        {
            w[i] = w[i]/wt;
        }



        //use the calculated weights to distribute flow
        for (int i=0; i<4; i++)
        {
            int r2, c2;

            //must multiply the cell directions by the sign of the slope vector components
            r2 = r+yn*dy[i];
            c2 = c+xn*dx[i];

            if(INSIDE(r2,c2) && !pcr::isMV(LDD->data[r2][c2]))
            {

                if(hmx->data[r2][c2] > 0)
                {
                    //weight * the flow is distributed to the ith cell that neighbours the advected flow.
                    MBLNFlood->data[r2][c2] +=  w[i]* dt * qbl;
                    MBLNFlood->data[r][c] -=  w[i]* dt* qbl;
                    MSSNFlood->data[r2][c2] +=  w[i]* dt * qss;
                    MSSNFlood->data[r][c] -=  w[i]* dt* qss;



                }

            }
        }
    }


    FOR_CELL_IN_FLOODAREA

        _BL->Drc = MBLNFlood->Drc;
        _SS->Drc = MSSNFlood->Drc;

        //set concentration from present sediment
        _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _BL->Drc);

        //set concentration from present sediment
        _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _SS->Drc);

    }

}

void TWorld::SWOFSedimentCheckZero(int r, int c, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
{

        if(!(hmx->Drc > 0))
        {
            BLDepFloodT->Drc += -(_BL->Drc);
            BLDepFloodT->Drc += -(_SS->Drc);

            _BL->Drc = 0;

            //set concentration from present sediment
            _BLC->Drc = 0;

            _SS->Drc  = 0;

            //set concentration from present sediment
            _SSC->Drc = 0;

        }


}

void TWorld::SWOFSedimentSetConcentration(int r, int c, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
{

        if((hmx->Drc > 0))
        {
            //set concentration from present sediment
            _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*BLDepthFlood->Drc, _BL->Drc);

            //set concentration from present sediment
            _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*SSDepthFlood->Drc, _SS->Drc);
        }else
        {
            //set concentration from present sediment
            _BL->Drc = 0;

            //set concentration from present sediment
            _SS->Drc = 0;

        }




}
void TWorld::SWOFSedimentDiffusionWS(int l, double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
{
    FOR_WATERSHED_ROW_COL(l) {

        MBLNFlood->Drc = _BL->Drc;
        MSSNFlood->Drc = _SS->Drc;
        MBLFlood->Drc = _BL->Drc;
        MSSFlood->Drc = _SS->Drc;

        //set concentration from present sediment
        MBLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MBLFlood->Drc);

        //set concentration from present sediment
        MSSCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MSSFlood->Drc);
    }}


    //diffusion of Suspended Sediment layer
     FOR_WATERSHED_ROW_COL(l) {

         if(SSDepthFlood->Drc < MIN_HEIGHT)
         {
                continue;
         }

        //cell sizes
        double cdx = DX->Drc;
        double cdy = _dx;
        //here it is about spacing, not flow width, so use _dx instead of CHannelAdj->Drc

        //mixing coefficient
        double sigma = 1;
        double dux1 = std::abs(Uflood->data[r][c] - Uflood->data[r][c-1]);
        double dvy1 = std::abs(Vflood->data[r][c] - Vflood->data[r-1][c]);
        double dvx1 = std::abs(Vflood->data[r][c] - Vflood->data[r][c-1]);
        double duy1 = std::abs(Uflood->data[r][c] - Uflood->data[r-1][c]);
        double dux2 = std::abs(Uflood->data[r][c+1] - Uflood->data[r][c]);
        double dvy2 = std::abs(Vflood->data[r+1][c] - Vflood->data[r][c]);
        double dvx2 = std::abs(Vflood->data[r][c+1] - Vflood->data[r][c]);
        double duy2 = std::abs(Uflood->data[r+1][c] - Uflood->data[r][c]);

        //drop term if at boundary
        if(std::isnan(dux1))
            dux1 = 0;
        if(std::isnan(dvx1))
            dvx1 = 0;
        if(std::isnan(duy1))
            duy1 = 0;
        if(std::isnan(dvy1))
            dvy1 = 0;

        if(std::isnan(dux2))
            dux2 = 0;
        if(std::isnan(dvx2))
            dvx2 = 0;
        if(std::isnan(duy2))
            duy2 = 0;
        if(std::isnan(dvy2))
            dvy2 = 0;

        double dux = std::max(dux1,dux2);
        double dvy = std::max(dvy1,dvy2);
        double dvx = std::max(dvx1,dvx2);
        double duy = std::max(duy1,duy2);



        double eddyvs = cdx * cdy * sqrt(dux*dux + dvy*dvy +  0.5 * (dvx +duy)*(dvx +duy));
        double eta = eddyvs/sigma;

        if(std::isnan(eta))
        {
            eta = 0;
            qDebug() << r << c << "NAN diffusion";
        }
        //cell directions
        int dx[4] = {0, 1, -1, 0};
        int dy[4] = {1, 0, 0, -1};

        //use the calculated weights to distribute flow
        for (int i=0; i<4; i++)
        {
            int r2, c2;

            //must multiply the cell directions by the sign of the slope vector components
            r2 = r+dy[i];
            c2 = c+dx[i];

            double coeff = std::min(dt*eta *std::min(1.0,SSDepthFlood->data[r2][c2]/SSDepthFlood->data[r][c]),courant_factor_diffusive/4.0) * MSSFlood->Drc;

            if(INSIDE(r2,c2) && !pcr::isMV(LDD->data[r2][c2]))
            {
                MSSNFlood->data[r2][c2] += coeff;
                MSSNFlood->data[r][c] -= coeff;
            }
        }

    }}

     FOR_WATERSHED_ROW_COL(l) {

         _BL->Drc = MBLNFlood->Drc;
         _SS->Drc = MSSNFlood->Drc;

         //set concentration from present sediment
         _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _BL->Drc);

         //set concentration from present sediment
         _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _SS->Drc);
     }}


}

void TWorld::SWOFSedimentDiffusion(double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
{


    FOR_CELL_IN_FLOODAREA

        MBLNFlood->Drc = _BL->Drc;
        MSSNFlood->Drc = _SS->Drc;
        MBLFlood->Drc = _BL->Drc;
        MSSFlood->Drc = _SS->Drc;

        //set concentration from present sediment
        MBLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MBLFlood->Drc);

        //set concentration from present sediment
        MSSCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, MSSFlood->Drc);

    }


    //diffusion of Suspended Sediment layer
    FOR_CELL_IN_FLOODAREA


        //cell sizes
        double cdx = DX->Drc;
        double cdy = _dx;
        //here it is about spacing, not flow width, so use _dx instead of CHannelAdj->Drc

        //mixing coefficient
        double sigma = 1;
        double dux1 = std::abs(Uflood->data[r][c] - Uflood->data[r][c-1]);
        double dvy1 = std::abs(Vflood->data[r][c] - Vflood->data[r-1][c]);
        double dvx1 = std::abs(Vflood->data[r][c] - Vflood->data[r][c-1]);
        double duy1 = std::abs(Uflood->data[r][c] - Uflood->data[r-1][c]);
        double dux2 = std::abs(Uflood->data[r][c+1] - Uflood->data[r][c]);
        double dvy2 = std::abs(Vflood->data[r+1][c] - Vflood->data[r][c]);
        double dvx2 = std::abs(Vflood->data[r][c+1] - Vflood->data[r][c]);
        double duy2 = std::abs(Uflood->data[r+1][c] - Uflood->data[r][c]);

        //drop term if at boundary
        if(std::isnan(dux1))
            dux1 = 0;
        if(std::isnan(dvx1))
            dvx1 = 0;
        if(std::isnan(duy1))
            duy1 = 0;
        if(std::isnan(dvy1))
            dvy1 = 0;

        if(std::isnan(dux2))
            dux2 = 0;
        if(std::isnan(dvx2))
            dvx2 = 0;
        if(std::isnan(duy2))
            duy2 = 0;
        if(std::isnan(dvy2))
            dvy2 = 0;

        double dux = std::max(dux1,dux2);
        double dvy = std::max(dvy1,dvy2);
        double dvx = std::max(dvx1,dvx2);
        double duy = std::max(duy1,duy2);

        //diffusion coefficient according to J.Smagorinski (1964)
        double eddyvs = cdx * cdy * sqrt(dux*dux + dvy*dvy +  0.5 * (dvx +duy)*(dvx +duy));
        double eta = eddyvs/FS_SigmaDiffusion;

        if(std::isnan(eta))
        {
            eta = 0;
            qDebug() << r << c << "NAN diffusion";
        }
        //cell directions
        int dx[4] = {0, 1, -1, 0};
        int dy[4] = {1, 0, 0, -1};


        //use the calculated weights to distribute flow
        for (int i=0; i<4; i++)
        {
            int r2, c2;

            //must multiply the cell directions by the sign of the slope vector components
            r2 = r+dy[i];
            c2 = c+dx[i];

            double coeff = std::min(dt*eta *std::min(1.0,SSDepthFlood->data[r2][c2]/SSDepthFlood->data[r][c]),courant_factor_diffusive/4.0) * MSSFlood->Drc;

            if(INSIDE(r2,c2) && !pcr::isMV(LDD->data[r2][c2]))
            {
                MSSNFlood->data[r2][c2] += coeff;
                MSSNFlood->data[r][c] -= coeff;
            }
        }
    }

    FOR_CELL_IN_FLOODAREA

        _BL->Drc = MBLNFlood->Drc;
        _SS->Drc = MSSNFlood->Drc;

        //set concentration from present sediment
        _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _BL->Drc);

        //set concentration from present sediment
        _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, _SS->Drc);

    }


}

void TWorld::SWOFSedimentMaxC(int r, int c, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
{

    //maximum concentraion

        _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*SSDepthFlood->Drc, _SS->Drc);
        // limit concentration to 850 and throw rest in deposition

        double sssmax = MAXCONC * DX->Drc *ChannelAdj->Drc*SSDepthFlood->Drc;
        if(sssmax < _SS->Drc)
        {
            BLDepFloodT->Drc += -(_SS->Drc - sssmax);

            _SS->Drc = sssmax;

        }


        //set concentration from present sediment
        _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*BLDepthFlood->Drc, _BL->Drc);

        double smax = MAXCONCBL * DX->Drc *ChannelAdj->Drc*BLDepthFlood->Drc;
        if(smax < _BL->Drc)
        {
            BLDepFloodT->Drc += -(_BL->Drc - smax);
            _BL->Drc = smax;

        }




}

double TWorld::SWOFSedimentTCBL(int r, int c)
{

    if(hmx->Drc < he_ca || BLDepthFlood->Drc < he_ca)
    {
        return 0;
    }
    double v = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);
    if(v < he_ca)
    {
        return 0;
    }
    if(ChannelAdj->Drc < he_ca)
    {
        return 0;
    }

        if(FS_BL_Method == FSGOVERS)
        {
            //Govers with a maximum bed load layer depth (1980)
            double discharge = v * ChannelAdj->Drc * BLDepthFlood->Drc;
            //### Calc transport capacity
            double omega = 100.0*v*discharge;
            // V in cm/s in this formula assuming grad is SINE
            double omegacrit = 0.4;
            // critical unit streampower in cm/s
            return std::min(MAXCONCBL, 2650 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));
            // not more than 2650*0.32 = 848 kg/m3

        }else if(FS_BL_Method == FSRIJN)
        {
            //Van rijn simplified (1984)


            double ps = 2400.0;
            double pw = 1000.0;
            double ucr;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            if( d50m < 0.005)
            {
               ucr  = 0.19 * pow(d50m, 0.1) * log10(4.0* hmx->Drc/d90m);
            }else
            {
               ucr  = 0.19 * pow(d50m, 0.6) * log10(4.0* hmx->Drc/d90m);
            }
            double me = std::max((v-ucr)/(sqrt(GRAV * d50m * ((ps/pw) - 1.0))),0.0);
            double qs = 0.005 * ps*v *hmx->Drc * pow(d50m/hmx->Drc,1.2) * pow(me, 2.4);
            double tc =  qs/ (v * BLDepthFlood->Drc );
            return std::max(std::min( tc,MAXCONCBL ),0.0);
        }else
        {
            //van Rijn full (1980)

            double ps = 2400.0;
            double pw = 1000.0;
            double kinvis = 1.0;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);

            double ds = D50->Drc * pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0));
            double dh = hmx->Drc;
            double chevey = 18 * log(4 * dh/d90m);
            double us = sqrt(GRAV) * v/chevey;
            double uscr = 0.055;

            if(ds <150 && !(ds < 20))
                uscr = 0.013*pow(ds,0.29);
            if(ds < 20 && !(ds < 10))
                uscr = 0.04*pow(ds,-0.10);
            if(ds <10 && !(ds < 4))
                uscr = 0.14*pow(ds,-0.64);
            if(ds <4)
                uscr = 0.24*pow(ds,-1);

            uscr = sqrt(uscr * (ps/pw - 1)*GRAV * d50m);
            double T = std::max((us*us/(uscr*uscr) - 1),0.0);
            double qs = 0.053 * (pow(T,2.1)/pow(ds,0.3)) * sqrt((ps/pw -1)*GRAV)*d50m*sqrt(d50m);
            double tc =  ps * ChannelAdj->Drc * qs/ (v * BLDepthFlood->Drc*ChannelAdj->Drc);

            return std::max(std::min(tc,MAXCONCBL ),0.0);


        }

}

double TWorld::SWOFSedimentTCSS(int r, int c)
{

    if(hmx->Drc < he_ca || SSDepthFlood->Drc < he_ca)
    {
        return 0;
    }
    double v = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);
    if(v < he_ca)
    {
        return 0;
    }
    if(ChannelAdj->Drc < he_ca)
    {
        return 0;
    }

        if(FS_SS_Method == FSRIJN)
        {
            //Van rijn simplified (1984)


            double ucr;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            double ps = 2400.0;
            double pw = 1000.0;
            double mu = 1.0;
            if( d50m < 0.0005)
            {
               ucr  = 0.19 * pow(d50m, 0.1) * log10(4.0* SSDepthFlood->Drc/d90m);
            }else
            {
               ucr  = 8.5 * pow(d50m, 0.6) * log10(4.0* SSDepthFlood->Drc/d90m);
            }
            double me = std::max((v-ucr)/sqrt(GRAV * d50m * (ps/pw - 1)),0.0);
            double ds = d50m * GRAV * ((ps/pw)-1)/(mu*mu);
            double qs = hmx->Drc * 0.008 * ps*v * d50m * pow(me, 2.4) * pow(ds, -0.6);

            double tc =  qs/ (v * SSDepthFlood->Drc);
            return std::max(std::min(tc,MAXCONC),0.0);
        }else
        {

            //van Rijn full (1980)
            //van Rijn full (1980)

            double kinvis = 1.0;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            double ps = 2400.0;
            double pw = 1000.0;
            double ds = D50->Drc * pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0));
            double dh = hmx->Drc;
            double chevey = 18 * log10(4 * dh/d90m);
            double a = 0.1;

            double us = v* sqrt(GRAV)/chevey;
            double uscr = 0.055;

            if(ds <150 && !(ds < 20))
                uscr = 0.013*pow(ds,0.29);
            if(ds < 20 && !(ds < 10))
                uscr = 0.04*pow(ds,-0.10);
            if(ds <10 && !(ds < 4))
                uscr = 0.14*pow(ds,-0.64);
            if(ds <4)
                uscr = 0.24*pow(ds,-1);

            uscr = sqrt(uscr * (ps/pw - 1)*GRAV * d50m);

            double T = std::max((us*us/(uscr*uscr) - 1),0.0);
            double bsv = sqrt(GRAV * hmx->Drc *std::max(Grad->Drc,0.05));
            double ca = 0.015 * (d50m/a) * pow(T,1.5)/pow(ds,0.3);

            double dss = 1 + 0.011*(1.8 - 1)*(T - 25);
            double sv = 10 * (kinvis/ds) *( sqrt(1 + (ps/pw - 1) * GRAV * d50m*d50m*d50m) - 1);

            double beta = std::min(1.0 + 2.0*(sv/bsv)*(sv/bsv),5.0);
            double powcb =1;
            if(BLCFlood->Drc > 0)
            {
                powcb = 0.1;//pow(ca/BLCFlood->Drc,0.4);
            }
            double phi = 2.5 * pow(sv/bsv,0.8) * powcb;
            double Z = sv/(beta*bsv*0.40);
            double Zs = Z + phi;
            double ad = 0.1;
            double F = (pow(ad,Zs) - pow(ad,1.2))/(pow(1.0-ad,Zs)* (1.2 - Zs));
            double qs = F * v * hmx->Drc * ca;
            double tc = ps * ChannelAdj->Drc * qs/ (v * SSDepthFlood->Drc * ChannelAdj->Drc);
            return std::max(std::min(tc,MAXCONC ),0.0);
        }


}

void TWorld::SWOFSedimentLayerDepth(int r , int c)
{
    double d50m = (D50->Drc/1000000);
    double d90m = (D90->Drc/1000000);
    double ps = 2400;
    double pw = 1000;
    double velocity = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);

    //critical shear velocity for bed level motion by van rijn
    double critshearvel = velocity * sqrt(GRAV)/(18 * log10(4*(ChannelAdj->Drc * hmx->Drc/(hmx->Drc * 2 + ChannelAdj->Drc))/d90m));
    //critical shear velocity for bed level motion by van rijn
    double critsheart = (critshearvel*critshearvel)/ (((ps-pw)/pw) * GRAV*d50m);
    //rough bed bed load layer depth by Hu en Hui
    BLDepthFlood->Drc = std::min(std::min(d50m * 1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), hmx->Drc), 0.1);
    SSDepthFlood->Drc = std::max(hmx->Drc - BLDepthFlood->Drc,0.0);
}

void TWorld::SWOFSedimentDet(double dt, int r,int c)
{

    SWOFSedimentLayerDepth(r,c);



    double velocity = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);

    //discharges for both layers and watervolumes
    double bldischarge = velocity * ChannelAdj->Drc * BLDepthFlood->Drc;
    double blwatervol = ChannelAdj->Drc *DX->Drc*BLDepthFlood->Drc;

    double ssdischarge = velocity * ChannelAdj->Drc * SSDepthFlood->Drc;
    double sswatervol = ChannelAdj->Drc *DX->Drc * SSDepthFlood->Drc;

    BLTCFlood->Drc = SWOFSedimentTCBL(r,c);
    SSTCFlood->Drc = SWOFSedimentTCSS(r,c);

    double deposition;
    if(BLDepthFlood->Drc < he_ca)
    {
        BLTCFlood->Drc = 0;
    }
    if(SSDepthFlood->Drc < he_ca)
    {
        SSTCFlood->Drc = 0;
    }
    if(hmx->Drc == 0)
    {
        BLTCFlood->Drc = 0;
        BLDepFlood->Drc = 0;
        BLDetFlood->Drc = 0;
        deposition = -BLFlood->Drc;
        BLDepFloodT->Drc += deposition;
        BLFlood->Drc = 0;
        BLCFlood->Drc = 0;

        SSTCFlood->Drc = 0;
        deposition = -SSFlood->Drc;
        BLDepFloodT->Drc += deposition;
        SSFlood->Drc = 0;
        SSCFlood->Drc = 0;
    }else
    {

       //first check if sediment goes to suspended sediment layer or to bed layer
       double tobl = 0;
       double toss = 0;
        double TransportFactor;
       if (SSDepthFlood->Drc > MIN_HEIGHT)
       {
           TransportFactor = (1-exp(-dt*SettlingVelocity->Drc/SSDepthFlood->Drc)) * sswatervol;
       }else
       {
            TransportFactor =  1* sswatervol;
       }
       double maxTC = std::max(BLTCFlood->Drc - BLCFlood->Drc,0.0);
       // positive difference: TC deficit becomes detachment (ppositive)
       double minTC = std::min(BLTCFlood->Drc - BLCFlood->Drc,0.0);

       tobl = TransportFactor * minTC * BLDepthFlood->Drc * DX->Drc * ChannelAdj->Drc;

       //vertical diffusion coefficient (van rijn)
       //double emax = 0.25 * hmx->Drc * velocity * 0.41;
       //double ez = 4* BLDepthFlood->Drc/hmx->Drc *(1 - BLDepthFlood->Drc/hmx->Drc) * emax;
       //toss = std::min(dt*ez *DX->Drc * ChannelAdj->Drc * BLFlood->Drc,maxTC* SSDepthFlood->Drc * DX->Drc * ChannelAdj->Drc);



       //to bed load
       tobl = std::max(tobl,-SSFlood->Drc);
       BLFlood->Drc -= tobl;
       SSFlood->Drc += tobl;

       //to suspended sediment
       //toss = std::min(toss,BLFlood->Drc);
       //BLFlood->Drc -= toss;
       //SSFlood->Drc += toss;

       TransportFactor = dt*SettlingVelocity->Drc * DX->Drc * SoilWidthDX->Drc;
       SSDetFlood->Drc = Y->Drc * maxTC * TransportFactor;
       if (GrassFraction->Drc > 0)
          SSDetFlood->Drc = (1-GrassFraction->Drc) * SSDetFlood->Drc;
       SSDetFlood->Drc = (1-StoneFraction->Drc) * SSDetFlood->Drc ;
       if (SwitchHardsurface)
          SSDetFlood->Drc = (1-HardSurface->Drc) * SSDetFlood->Drc ;
       if (SwitchHouses)
          SSDetFlood->Drc = (1-HouseCover->Drc)*SSDetFlood->Drc;

       //### sediment balance
       SSFlood->Drc += SSDetFlood->Drc;
        SSDetFloodT->Drc += SSDetFlood->Drc;

        BLCFlood->Drc = MaxConcentration(blwatervol, BLFlood->Drc);
        // limit concentration to 850 and throw rest in deposition
        SSCFlood->Drc = MaxConcentration(sswatervol, SSFlood->Drc);
       // limit concentration to 850 and throw rest in deposition

       double sssmax = MAXCONC * DX->Drc *ChannelAdj->Drc*SSDepthFlood->Drc;
       if(sssmax < SSFlood->Drc)
       {
           BLFlood->Drc += (SSFlood->Drc - sssmax);
           SSFlood->Drc = sssmax;

       }



       //deposition and detachment
       //### calc concentration and net transport capacity
       BLDepFlood->Drc = 0;
       // init deposition for this timestep
       BLCFlood->Drc = MaxConcentration(blwatervol, BLFlood->Drc);
       // limit sed concentration to max

       maxTC = std::max(BLTCFlood->Drc - BLCFlood->Drc,0.0);
       // positive difference: TC deficit becomes detachment (ppositive)
       minTC = std::min(BLTCFlood->Drc - BLCFlood->Drc,0.0);
       // negative difference: TC surplus becomes deposition (negative)
       // unit kg/m3

       //### detachment
       TransportFactor = dt*SettlingVelocity->Drc * DX->Drc *SoilWidthDX->Drc;
       // detachment can only come from soil, not roads (so do not use flowwidth)
       // units s * m/s * m * m = m3

       BLDetFlood->Drc = Y->Drc * maxTC * TransportFactor;
       // unit = kg/m3 * m3 = kg
       //BLDetFlood->Drc = std::min(BLDetFlood->Drc, maxTC * bldischarge*dt);
       // cannot have more detachment than remaining capacity in flow
       // use discharge because standing water has no erosion

       if (GrassFraction->Drc > 0)
          BLDetFlood->Drc = (1-GrassFraction->Drc) * BLDetFlood->Drc;
       // no flow detachment on grass strips

       // Detachment edxceptions:
       BLDetFlood->Drc = (1-StoneFraction->Drc) * BLDetFlood->Drc ;
       // no flow detachment on stony surfaces

       if (SwitchHardsurface)
          BLDetFlood->Drc = (1-HardSurface->Drc) * BLDetFlood->Drc ;
       // no flow detachment on hard surfaces

       if (SwitchHouses)
          BLDetFlood->Drc = (1-HouseCover->Drc)*BLDetFlood->Drc;
       // no flow det from house roofs






       // IN KG/CELL

       //DETFlow->Drc = (1-Snowcover->Drc) * DETFlow->Drc ;
       /* TODO: CHECK THIS no flow detachment on snow */
       //is there erosion and sedimentation under the snowdeck?

       //### deposition
       if (BLDepthFlood->Drc > MIN_HEIGHT)
          TransportFactor = (1-exp(-dt*SettlingVelocity->Drc/BLDepthFlood->Drc)) * blwatervol;
       else
          TransportFactor = 1*blwatervol;
       // if settl velo is very small, transportfactor is 0 and depo is 0
       // if settl velo is very large, transportfactor is 1 and depo is max

       //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
       // deposition can occur on roads and on soil (so use flowwidth)

       double deposition = minTC * TransportFactor;
       // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

       if (SwitchLimitDepTC)
          deposition = std::max(deposition, minTC *blwatervol);
       // cannot be more than sediment above capacity
       deposition = std::max(deposition, -BLFlood->Drc);
       // cannot have more depo than sediment present
       /* TODO what about this: which one to choose */

       //deposition = (1-Snowcover->Drc) * deposition;
       /* TODO: TRUE??? no deposition on snow */
       //if (GrassFraction->Drc > 0)
       //   deposition = -Sed->Drc*GrassFraction->Drc + (1-GrassFraction->Drc)*deposition;
       // generate 100% deposition on grassstrips
       //? bit tricky, maximizes effect on grassstrips ?


       BLDepFloodT->Drc += deposition;
       BLDetFloodT->Drc += BLDetFlood->Drc;
       // IN KG/CELL

       //### sediment balance
       BLFlood->Drc += BLDetFlood->Drc;
       BLFlood->Drc += deposition;




       BLCFlood->Drc = MaxConcentration(blwatervol, BLFlood->Drc);
       // limit concentration to 850 and throw rest in deposition

    }

    SWOFSedimentMaxC(r,c,BLFlood, BLCFlood, SSFlood,SSCFlood);

}

void TWorld::SWOFSediment(double dt)
{
    if (!SwitchErosion)
       return;

    FOR_CELL_IN_FLOODAREA

        SWOFSedimentDet(dt,r,c);
    }

    FOR_CELL_IN_FLOODAREA
        SWOFSedimentCheckZero(r,c, BLFlood, BLCFlood, SSFlood,SSCFlood);
        SWOFSedimentSetConcentration(r,c, BLFlood, BLCFlood, SSFlood,SSCFlood);
    }

    if(SwitchFloodSedimentMethod)
    {
        SWOFSedimentFlowInterpolation(dt, BLFlood, BLCFlood, SSFlood,SSCFlood);
    }else
    {

        SWOFSedimentFlow(dt, BLFlood, BLCFlood, SSFlood,SSCFlood);

    }

    SWOFSedimentDiffusion(dt, BLFlood, BLCFlood, SSFlood,SSCFlood);


    FOR_CELL_IN_FLOODAREA
        SWOFSedimentMaxC(r,c, BLFlood, BLCFlood, SSFlood,SSCFlood);
    }


}

void TWorld::SWOFSedimentWS(int l, double dt)
{
    if (!SwitchErosion)
       return;

    FOR_WATERSHED_ROW_COL(l) {

        SWOFSedimentDet(dt,r,c);

    }}

    FOR_WATERSHED_ROW_COL(l) {
        SWOFSedimentCheckZero(r,c, BLFlood, BLCFlood, SSFlood,SSCFlood);
        SWOFSedimentSetConcentration(r,c, BLFlood, BLCFlood, SSFlood,SSCFlood);
    }}

    if(SwitchFloodSedimentMethod)
    {
        SWOFSedimentFlowInterpolationWS(l, dt, BLFlood, BLCFlood, SSFlood,SSCFlood);
    }else
    {

        SWOFSedimentFlowWS(l,dt, BLFlood, BLCFlood, SSFlood,SSCFlood);
    }

    FOR_WATERSHED_ROW_COL(l) {
            SWOFSedimentDiffusionWS(l,dt, BLFlood, BLCFlood, SSFlood,SSCFlood);
    }}

    FOR_WATERSHED_ROW_COL(l) {
            SWOFSedimentMaxC(r,c, BLFlood, BLCFlood, SSFlood,SSCFlood);
    }}

}

