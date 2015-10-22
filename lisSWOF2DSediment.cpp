

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

void TWorld::FS_Flux(cTMap * _s)
{


    FOR_CELL_IN_FLOODAREA
        if (c > 0 && !pcr::isMV(LDD->data[r][c-1]))
    {
      bl1d->data[r][c-1] = std::max(0.0, bl1r->data[r][c-1] );
      bl1g->Drc          = std::max(0.0, bl1l->Drc  );

      h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1]));
      h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1]));

      FS_HLL(h2d->data[r][c-1], bl1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc,bl1g->Drc, u1l->Drc, v1l->Drc);

      blf1->Drc = HLL2_f1;

    }
    else
    {
      bl1d->data[r][c] = std::max(0.0, bl1r->data[r][c] );
      bl1g->Drc        = std::max(0.0, bl1l->Drc);

      h1d->data[r][c] = std::max(0.0, h1r->data[r][c] - std::max(0.0,  delz1->data[r][c]));
      h1g->Drc        = std::max(0.0, h1l->Drc        - std::max(0.0, -delz1->data[r][c]));

     FS_HLL(h1d->data[r][c],bl1d->data[r][c], u1r->data[r][c], v1r->data[r][c],h1g->Drc,bl1g->Drc, u1l->Drc, v1l->Drc);

      blf1->Drc = HLL2_f1;

    }}

    FOR_CELL_IN_FLOODAREA
    if(r > 0 && !pcr::isMV(LDD->data[r-1][c]))
    {

      bl2d->data[r-1][c] = std::max(0.0, bl2r->data[r-1][c]);
      bl2g->Drc          = std::max(0.0, bl2l->Drc);

      h2d->data[r-1][c] = std::max(0.0, h2r->data[r-1][c] - std::max(0.0,  delz2->data[r-1][c]));
      h2g->Drc          = std::max(0.0, h2l->Drc          - std::max(0.0, -delz2->data[r-1][c]));

      FS_HLL(h2d->data[r-1][c],bl2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c],h2g->Drc,bl2g->Drc,v2l->Drc,u2l->Drc);

      blg1->Drc = HLL2_f1;

    }
    else
    {
       bl2d->data[r][c] = std::max(0.0, bl2r->data[r][c] );
       bl2g->Drc        = std::max(0.0, bl2l->Drc);

       h2d->data[r][c] = std::max(0.0, h2r->data[r][c] - std::max(0.0,  delz2->data[r][c]));
       h2g->Drc        = std::max(0.0, h2l->Drc        - std::max(0.0, -delz2->data[r][c]));

       FS_HLL(h2d->data[r][c],bl2d->data[r][c],v2r->data[r][c],u2r->data[r][c],h2g->Drc,bl2g->Drc,v2l->Drc,u2l->Drc);

       blg1->Drc = HLL2_f1;

    }}

}

void TWorld::FS_MUSCLE(cTMap * _s)
{

    double delta_s1, delta_u1, delta_v1;
    double delta_s2, delta_u2, delta_v2;
    double ds, du, dv;

    // fill EW and NS conditions with cell itself, 1st order approximation used for boundary conditions
    FOR_CELL_IN_FLOODAREA {

        bl1r->Drc = _s->Drc;
        bl1l->Drc = _s->Drc;
        bl2r->Drc = _s->Drc;
        bl2l->Drc = _s->Drc;

    }}


    FOR_CELL_IN_FLOODAREA
      if(c > 0 && c < _nrCols-1 && !MV(r,c-1) &&  !MV(r,c+1))
    {
        delta_s1 = _s->Drc - _s->data[r][c-1];

        delta_s2 = _s->data[r][c+1] - _s->Drc;

        ds   = 0.5*limiter(delta_s1, delta_s2);


        bl1r->Drc = _s->Drc+ds;
        bl1l->Drc = _s->Drc-ds;
    }}


      FOR_CELL_IN_FLOODAREA
        if(r > 0 && MV(r-1,c) && r < _nrRows-1 && !MV(r+1, c))
      {
          delta_s1 = _s->Drc - _s->data[r-1][c];

          delta_s2 = _s->data[r+1][c] - _s->Drc;

          ds   = 0.5*limiter(delta_s1, delta_s2);
          bl2r->Drc = _s->Drc+ds;
          bl2l->Drc = _s->Drc-ds;
      }}

}

void TWorld::FS_Simple(cTMap * _s)
{
    FOR_CELL_IN_FLOODAREA {
        bl1r->Drc = _s->Drc;
        bl1l->Drc = _s->Drc;

        bl2r->Drc = _s->Drc;
        bl2l->Drc = _s->Drc;
    }}
}


void TWorld::FS_MainCalc(cTMap * _h,cTMap * _s, cTMap * _ss, double dt)
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
        _ss->Drc = _s->Drc;
      else
        _ss->Drc = _s->Drc - tx*(blf1->data[r][c+1]-blf1->Drc) - ty*(blg1->data[r+1][c]-blg1->Drc);
    }}
}

void TWorld::FS_HLL(double h_L,double s_L,double u_L,double v_L,double h_R, double s_R,double u_R,double v_R)
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
        }else if (c1>=EPSILON){ //supercritical flow, from left to right : we have std::max(abs(c1),abs(c2))=c2>0
            f1=s_L*(q_L);
        }else if (c2<=-EPSILON){ //supercritical flow, from right to left : we have std::max(abs(c1),abs(c2))=-c1>0
            f1=s_R*(q_R);
        }else{ //subcritical flow
            tmp = 1./(c2-c1);
            f1=(c2*s_L*q_L-c1*s_R*q_R)*tmp+c1*c2*(h_R-h_L)*tmp;
        }
    }
    HLL2_f1 = f1;


}


void TWorld::SWOFSedimentFlow(double dt)
{

    cTMap*_h = hmx;
    //calculate concentration for transport
    FOR_ROW_COL_MV
    {
        BLTCFlood->Drc = 0;
        //set concentration from present sediment
        BLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, BLFlood->Drc);

    }


    //reconstruction scheme

    FS_Simple(BLNFlood);

    //FS_MUSCLE(SCFlood);

    //flux calculation
    FS_Flux(BLCFlood);

    //Calculate new Sediment
    FS_MainCalc(_h,BLFlood,bls,dt);

    //update variable
    FOR_CELL_IN_FLOODAREA {
      BLFlood->Drc = bls->Drc;
    }}

    //diffusion of Suspended Sediment layer
    FOR_CELL_IN_FLOODAREA


        //cell sizes
        double cdx = DX->Drc;
        double cdy = _dx;

        //mixing coefficient
        double sigma = 0.75;
        double dvx = (Uflood->data[r][c+1] - Uflood->data[r-1][c-1])/cdx;
        double dvy = (Vflood->data[r+1][c] - Vflood->data[r-1][c])/cdy;
        double dvxy = (Uflood->data[r+1][c] - Uflood->data[r-1][c])/cdy + (Vflood->data[r][c+1] - Vflood->data[r-1][c-1])/cdx;
        double eddyvs = cdx * cdy * sqrt(dvx + dvy +  0.5 * (dvxy));
        double eta = eddyvs/sigma;

        //cell directions
        int dx[4] = {0, 1, -1, 0};
        int dy[4] = {1, 0, 0, -1};

        double coeff = std::min(dt/eta,0.2) * SSFlood->Drc;

        //use the calculated weights to distribute flow
        for (int i=0; i<4; i++)
        {
            int r2, c2;

            //must multiply the cell directions by the sign of the slope vector components
            r2 = r+dy[i];
            c2 = c+dx[i];

            if(INSIDE(r2,c2) && !pcr::isMV(LDD->data[r2][c2]))
            {
                SSFlood->data[r2][c2] += coeff;
                SSFlood->data[r][c] -= coeff;
            }
        }

    }

    //recalculate concentration
    FOR_ROW_COL_MV
    {
        //set concentration from present sediment
        BLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, BLFlood->Drc);
        double smax = MAXCONC * DX->Drc *ChannelAdj->Drc*hmx->Drc;
        if(smax < BLFlood->Drc)
        {
            BLDepFloodT->Drc += -(BLFlood->Drc - smax);
            BLFlood->Drc = smax;
            qDebug() << r <<c << DX->Drc << ChannelAdj->Drc << hmx->Drc <<Sed->Drc <<BLFlood->Drc <<"High concentration";

        }
    }

}

void TWorld::SWOFSedimentFlowInterpolation(double dt)
{

    FOR_ROW_COL_MV
    {
        BLTCFlood->Drc = 0;
        BLNFlood->Drc = BLFlood->Drc;

        //set concentration from present sediment
        BLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*BLDepthFlood->Drc, BLFlood->Drc);


        SSTCFlood->Drc = 0;
        SSNFlood->Drc = SSFlood->Drc;

        //set concentration from present sediment
        SSCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*SSDepthFlood->Drc, SSFlood->Drc);

    }



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

        double qs = vel*ChannelAdj->Drc *hmx->Drc * BLCFlood->Drc;

        if(qs > DX->Drc * ChannelAdj->Drc *hmx->Drc * BLCFlood->Drc)
        {
            qs =  DX->Drc * ChannelAdj->Drc *hmx->Drc * BLCFlood->Drc;
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
                    BLNFlood->data[r2][c2] +=  w[i]* dt * qs;
                    BLNFlood->data[r][c] -=  w[i]* dt* qs;

                }

            }
        }
    }


    //diffusion of Suspended Sediment layer
    FOR_CELL_IN_FLOODAREA


        //cell sizes
        double cdx = DX->Drc;
        double cdy = _dx;

        //mixing coefficient
        double sigma = 0.75;
        double dvx = (Uflood->data[r][c+1] - Uflood->data[r-1][c-1])/cdx;
        double dvy = (Vflood->data[r+1][c] - Vflood->data[r-1][c])/cdy;
        double dvxy = (Uflood->data[r+1][c] - Uflood->data[r-1][c])/cdy + (Vflood->data[r][c+1] - Vflood->data[r-1][c-1])/cdx;
        double eddyvs = cdx * cdy * sqrt(dvx + dvy +  0.5 * (dvxy));
        double eta = eddyvs/sigma;

        //cell directions
        int dx[4] = {0, 1, -1, 0};
        int dy[4] = {1, 0, 0, -1};

        double coeff = std::min(dt/eta,0.2) * SSFlood->Drc;

        //use the calculated weights to distribute flow
        for (int i=0; i<4; i++)
        {
            int r2, c2;

            //must multiply the cell directions by the sign of the slope vector components
            r2 = r+dy[i];
            c2 = c+dx[i];

            if(INSIDE(r2,c2) && !pcr::isMV(LDD->data[r2][c2]))
            {
                SSFlood->data[r2][c2] += coeff;
                SSFlood->data[r][c] -= coeff;
            }
        }

    }

    //maximum concentraion
    FOR_ROW_COL_MV
    {

        SSFlood->Drc = SSNFlood->Drc;
        SSCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*SSDepthFlood->Drc, SSFlood->Drc);
        // limit concentration to 850 and throw rest in deposition

        double sssmax = MAXCONC * DX->Drc *ChannelAdj->Drc*SSDepthFlood->Drc;
        if(sssmax < BLFlood->Drc)
        {
            BLFlood->Drc += (BLFlood->Drc - sssmax);
            SSFlood->Drc = sssmax;

        }


        //set concentration from present sediment
        BLFlood->Drc = BLNFlood->Drc;
        BLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*BLDepthFlood->Drc, BLFlood->Drc);

        double smax = MAXCONC * DX->Drc *ChannelAdj->Drc*hmx->Drc;
        if(smax < BLFlood->Drc)
        {
            BLDepFloodT->Drc += -(BLFlood->Drc - smax);
            BLFlood->Drc = smax;
            //qDebug() << r <<c << DX->Drc << ChannelAdj->Drc << hmx->Drc <<Sed->Drc <<BLFlood->Drc <<"High concentration";

        }



    }

}

double TWorld::SWOFSedimentTCBL(int r, int c)
{
    double maxlayer = 0.1;
    if(false )
    {
        //Govers with a maximum bed load layer depth (1980)

        tm->Drc = 0;
        double velocity = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);
        double discharge = velocity * ChannelAdj->Drc * hmx->Drc;
        //### Calc transport capacity
        double omega = 100.0*velocity*discharge;
        // V in cm/s in this formula assuming grad is SINE
        double omegacrit = 0.4;
        // critical unit streampower in cm/s
        return std::min(MAXCONC, 2650 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));
        // not more than 2650*0.32 = 848 kg/m3

    }else if(true)
    {
        //Van rijn simplified (1984)

        double v = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);
        double ps = 2400;
        double pw = 1000;
        double ucr;
        double d50m = (D50->Drc/1000000);
        double d90m = (D90->Drc/1000000);
        if( d50m < 0.005)
        {
           ucr  = 0.19 * pow(d50m, 0.1) * log(4.0* hmx->Drc/d90m);
        }else
        {
           ucr  = 0.19 * pow(d50m, 0.6) * log(4.0* hmx->Drc/d90m);
        }
        double me = (v - ucr)/(sqrt(GRAV * d50m * (ps/pw) - 1));
        double qs = 0.015 * 2400*v *hmx->Drc * pow(d50m/hmx->Drc,1.2) * pow(me, 1.5);
        return std::min(MAXCONC, qs/ (v * BLDepthFlood->Drc) );
    }else
    {
        //van Rijn full (1980)


    }

}

double TWorld::SWOFSedimentTCSS(int r, int c)
{
    double maxlayer = 0.1;

    if(true)
    {
        //Van rijn simplified (1984)

        double v = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);

        double ucr;
        double d50m = (D50->Drc/1000000);
        double d90m = (D90->Drc/1000000);
        double ps = 2400;
        double pw = 1000;
        double mu = 1;
        if( d50m < 0.005)
        {
           ucr  = 0.19 * pow(d50m, 0.1) * log(4.0* hmx->Drc/d90m);
        }else
        {
           ucr  = 0.19 * pow(d50m, 0.6) * log(4.0* hmx->Drc/d90m);
        }
        double me = (v - ucr)/(GRAV * d50m * sqrt((ps/pw) - 1));
        double ds = d50m * GRAV * ((ps/pw)-1)/(mu*mu);
        double qs = 0.012 * 2400*v * d50m * pow(mu, 2.4) * pow(ds, -0.6);
        return std::min(MAXCONC, qs/ (v * SSDepthFlood->Drc) );
    }else
    {
        //van Rijn full (1980)


    }

}

void TWorld::SWOFSedimentDet(double dt, int r,int c)
{

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, -1, -1, -1, 0, 0, 0, 1, 1, 1};

    double d50m = (D50->Drc/1000000);
    double d90m = (D90->Drc/1000000);
    double ps = 2400;
    double pw = 1000;
    double velocity = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);

    //critical shear velocity for bed level motion by van rijn
    double critshearvel = velocity * sqrt(GRAV)/(18 * log(4*hmx->Drc/d90m));
    //critical shear velocity for bed level motion by van rijn
    double critsheart = (critshearvel*critshearvel)/ (((ps-pw)/pw) * GRAV*d50m);
    //rough bed bed load layer depth by Hu en Hui
    BLDepthFlood->Drc = std::min(std::min(1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), hmx->Drc), 0.25);


    double bldischarge = velocity * ChannelAdj->Drc * BLDepthFlood->Drc;
    double blwatervol = ChannelAdj->Drc *DX->Drc*BLDepthFlood->Drc;

    SSDepthFlood->Drc = std::max(hmx->Drc - BLDepthFlood->Drc,0.0);
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

       double TransportFactor = (1-exp(-dt*SettlingVelocity->Drc/SSDepthFlood->Drc)) * sswatervol;

       double maxTC = std::max(BLTCFlood->Drc - BLCFlood->Drc,0.0);
       // positive difference: TC deficit becomes detachment (ppositive)
       double minTC = std::min(BLTCFlood->Drc - BLCFlood->Drc,0.0);

       tobl = TransportFactor * minTC;

       //vertical diffusion coefficient (van rijn)
       double emax = 0.25 * hmx->Drc * velocity * 0.41;
       double ez = 4* BLDepthFlood->Drc/hmx->Drc *(1 - BLDepthFlood->Drc/hmx->Drc) * emax;
       toss = std::min(dt/ez * BLFlood->Drc,maxTC);

       //to bed load
       BLFlood->Drc += tobl;
       SSFlood->Drc -= tobl;

       //to suspended sediment
       BLFlood->Drc += tobl;
       SSFlood->Drc -= tobl;

       SSCFlood->Drc = MaxConcentration(sswatervol, SSFlood->Drc);
       // limit concentration to 850 and throw rest in deposition

       double sssmax = MAXCONC * DX->Drc *ChannelAdj->Drc*SSDepthFlood->Drc;
       if(sssmax < SSFlood->Drc)
       {
           SSFlood->Drc += (SSFlood->Drc - sssmax);
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
       TransportFactor = dt*SettlingVelocity->Drc * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
       // detachment can only come from soil, not roads (so do not use flowwidth)
       // units s * m/s * m * m = m3

       BLDetFlood->Drc = Y->Drc * maxTC * TransportFactor;
       // unit = kg/m3 * m3 = kg
       BLDetFlood->Drc = std::min(BLDetFlood->Drc, maxTC * bldischarge*dt);
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

       BLDetFloodT->Drc += BLDetFlood->Drc;

       // IN KG/CELL

       //DETFlow->Drc = (1-Snowcover->Drc) * DETFlow->Drc ;
       /* TODO: CHECK THIS no flow detachment on snow */
       //is there erosion and sedimentation under the snowdeck?

       //### deposition
       if (WH->Drc > MIN_HEIGHT)
          TransportFactor = (1-exp(-dt*SettlingVelocity->Drc/BLDepthFlood->Drc)) * blwatervol;
       else
          TransportFactor = WaterVolall->Drc;
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
       // IN KG/CELL

       //### sediment balance
       BLFlood->Drc += BLDetFlood->Drc;
       BLFlood->Drc += deposition;

       BLCFlood->Drc = MaxConcentration(blwatervol, BLFlood->Drc);
       // limit concentration to 850 and throw rest in deposition

       double blsmax = MAXCONC * DX->Drc *ChannelAdj->Drc*BLDepthFlood->Drc;
       if(blsmax < BLFlood->Drc)
       {
           BLDepFloodT->Drc += -(BLFlood->Drc - blsmax);
           BLFlood->Drc = blsmax;
           //qDebug() << r <<c << DX->Drc << ChannelAdj->Drc << hmx->Drc <<Sed->Drc <<BLFlood->Drc <<"NAN2";

       }


    }

}

void TWorld::SWOFSedimentFlowWS(int wsnr, double dt)
{


}

void TWorld::SWOFSedimentDetWS(int wsnr, double dt)
{





}


void TWorld::SWOFSediment(double dt)
{
    if (!SwitchErosion)
       return;

    FOR_CELL_IN_FLOODAREA

        SWOFSedimentDet(dt,r,c);
    }

    if(SwitchFloodSedimentMethod)
    {
        SWOFSedimentFlowInterpolation(dt);
    }else
    {
        SWOFSedimentFlow(dt);
    }


}

void TWorld::SWOFSedimentWS(int l, double dt)
{
    if (!SwitchErosion)
       return;
}

