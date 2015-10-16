

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
      s1d->data[r][c-1] = std::max(0.0, s1r->data[r][c-1] );
      s1g->Drc          = std::max(0.0, s1l->Drc  );

      h1d->data[r][c-1] = std::max(0.0, h1r->data[r][c-1] - std::max(0.0,  delz1->data[r][c-1]));
      h1g->Drc          = std::max(0.0, h1l->Drc          - std::max(0.0, -delz1->data[r][c-1]));

      FS_HLL(h2d->data[r][c-1], s1d->data[r][c-1], u1r->data[r][c-1], v1r->data[r][c-1],h1g->Drc,s1g->Drc, u1l->Drc, v1l->Drc);

      sf1->Drc = HLL2_f1;

    }
    else
    {
      h1d->data[r][c] = std::max(0.0, h1r->data[r][c] );
      h1g->Drc        = std::max(0.0, h1l->Drc);

      h1d->data[r][c] = std::max(0.0, h1r->data[r][c] - std::max(0.0,  delz1->data[r][c]));
      h1g->Drc        = std::max(0.0, h1l->Drc        - std::max(0.0, -delz1->data[r][c]));

     FS_HLL(h1d->data[r][c],s1d->data[r][c], u1r->data[r][c], v1r->data[r][c],h1g->Drc,s1g->Drc, u1l->Drc, v1l->Drc);

      sf1->Drc = HLL2_f1;

    }}

    FOR_CELL_IN_FLOODAREA
    if(r > 0 && !pcr::isMV(LDD->data[r-1][c]))
    {

      s2d->data[r-1][c] = std::max(0.0, s2r->data[r-1][c]);
      s2g->Drc          = std::max(0.0, s2l->Drc);

      h2d->data[r-1][c] = std::max(0.0, h2r->data[r-1][c] - std::max(0.0,  delz2->data[r-1][c]));
      h2g->Drc          = std::max(0.0, h2l->Drc          - std::max(0.0, -delz2->data[r-1][c]));

      FS_HLL(h2d->data[r-1][c],s2d->data[r-1][c],v2r->data[r-1][c],u2r->data[r-1][c],h2g->Drc,s2g->Drc,v2l->Drc,u2l->Drc);

      sg1->Drc = HLL2_f1;

    }
    else
    {
       s2d->data[r][c] = std::max(0.0, s2r->data[r][c] );
       s2g->Drc        = std::max(0.0, s2l->Drc);

       h2d->data[r][c] = std::max(0.0, h2r->data[r][c] - std::max(0.0,  delz2->data[r][c]));
       h2g->Drc        = std::max(0.0, h2l->Drc        - std::max(0.0, -delz2->data[r][c]));

       FS_HLL(h2d->data[r][c],s2d->data[r][c],v2r->data[r][c],u2r->data[r][c],h2g->Drc,s2g->Drc,v2l->Drc,u2l->Drc);

       sg1->Drc = HLL2_f1;

    }}

}

void TWorld::FS_MUSCLE(cTMap * _s)
{

    double delta_s1, delta_u1, delta_v1;
    double delta_s2, delta_u2, delta_v2;
    double ds, du, dv;

    // fill EW and NS conditions with cell itself, 1st order approximation used for boundary conditions
    FOR_CELL_IN_FLOODAREA {
        s1r->Drc = _s->Drc;


        s1l->Drc = _s->Drc;

        s2r->Drc = _s->Drc;

        s2l->Drc = _s->Drc;

    }}


    FOR_CELL_IN_FLOODAREA
      if(c > 0 && c < _nrCols-1 && !MV(r,c-1) &&  !MV(r,c+1))
    {
        delta_s1 = _s->Drc - _s->data[r][c-1];

        delta_s2 = _s->data[r][c+1] - _s->Drc;

        ds   = 0.5*limiter(delta_s1, delta_s2);


        s1r->Drc = _s->Drc+ds;
        s1l->Drc = _s->Drc-ds;
    }}


      FOR_CELL_IN_FLOODAREA
        if(r > 0 && MV(r-1,c) && r < _nrRows-1 && !MV(r+1, c))
      {
          delta_s1 = _s->Drc - _s->data[r-1][c];

          delta_s2 = _s->data[r+1][c] - _s->Drc;

          ds   = 0.5*limiter(delta_s1, delta_s2);
          s2r->Drc = _s->Drc+ds;
          s2l->Drc = _s->Drc-ds;
      }}

}

void TWorld::FS_Simple(cTMap * _s)
{
    FOR_CELL_IN_FLOODAREA {
        s1r->Drc = _s->Drc;
        s1l->Drc = _s->Drc;

        s2r->Drc = _s->Drc;
        s2l->Drc = _s->Drc;
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
        _ss->Drc = _s->Drc - tx*(sf1->data[r][c+1]-sf1->Drc) - ty*(sg1->data[r+1][c]-sg1->Drc);
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
            f1=(c2*s_L*q_L-c1*s_R*q_R)*tmp+c1*c2*(s_L*h_R-s_R*h_L)*tmp;
        }
    }
    HLL2_f1 = f1;


}


void TWorld::SWOFSedimentFlow(cTMap*_h,double dt)
{
    if (!SwitchErosion)
        return;

    if(SwitchFloodSedimentMethod)
    {
        return SWOFSedimentFlowInterpolation(dt);
    }

    //calculate concentration for transport
    FOR_ROW_COL_MV
    {
        TCFlood->Drc = 0;
        SCNFlood->Drc = SCFlood->Drc;
        SNFlood->Drc = SFlood->Drc;

        //set concentration from present sediment
        SCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, SFlood->Drc);

    }


    //reconstruction scheme

    //FS_Simple(SNFlood, u, v);

    FS_MUSCLE(SCFlood);

    //flux calculation
    FS_Flux(SCFlood);

    //Calculate new Sediment
    FS_MainCalc(_h,SFlood,ss,dt);

    //update variable
    FOR_CELL_IN_FLOODAREA {
      SFlood->Drc = ss->Drc;
    }}

    //recalculate concentration
    FOR_ROW_COL_MV
    {
        //set concentration from present sediment
        SCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, SFlood->Drc);

        SCFlood->Drc = SCNFlood->Drc;
        SFlood->Drc = SNFlood->Drc;
    }
}

void TWorld::SWOFSedimentFlowInterpolation(double dt)
{

    FOR_ROW_COL_MV
    {
        TCFlood->Drc = 0;
        SCNFlood->Drc = SCFlood->Drc;
        SNFlood->Drc = SFlood->Drc;

        //set concentration from present sediment
        SCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, SFlood->Drc);

    }



    //first calculate the weights for the cells that are closest to location that flow is advected to
    FOR_CELL_IN_FLOODAREA

        //no flood velocity means no flood sediment transport, so skip this cell
        if((Vflood->Drc == 0 && Uflood->Drc == 0))
        {
            continue;
        }

        //cell sizes
        double cdx = DX->Drc;
        double cdy = _dx;

        //the sign of the x and y direction of flow
        double yn = signf(Vflood->Drc);
        double xn = signf(Uflood->Drc);

        double vel = sqrt(Uflood->Drc*Uflood->Drc + Vflood->Drc*Vflood->Drc);

        //should not travel more distance than cell size
        double dsx = xn*std::min(fabs(Uflood->Drc)/vel,cdx);
        double dsy = yn*std::min(fabs(Vflood->Drc)/vel,cdy);

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

                //weight * the flow is distributed to the ith cell that neighbours the advected flow.
                    SNFlood->data[r2][c2] +=  w[i]* dt *vel*ChannelAdj->Drc *hmx->Drc * SCFlood->Drc;
                    SNFlood->data[r][c] -=  w[i]* dt*vel*ChannelAdj->Drc *hmx->Drc * SCFlood->Drc;

            }
        }
    }



    FOR_ROW_COL_MV
    {
        //set concentration from present sediment
        SCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, SFlood->Drc);

        SCFlood->Drc = SCNFlood->Drc;
        SFlood->Drc = SNFlood->Drc;
    }


}

void TWorld::SWOFSedimentDet(double dt)
{
    return;
    if (!SwitchErosion)
       return;

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, -1, -1, -1, 0, 0, 0, 1, 1, 1};

    FOR_CELL_IN_FLOODAREA

       tm->Drc = 0;
       double velocity = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);
       double discharge = velocity * ChannelAdj->Drc * hmx->Drc;
       //### Calc transport capacity
       double omega = 100*velocity*discharge;
       // V in cm/s in this formula assuming grad is SINE
       double omegacrit = 0.4;
       // critical unit streampower in cm/s
       TCFlood->Drc = std::min(MAXCONC, 2650 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));
       // not more than 2650*0.32 = 848 kg/m3
    }

    //VJ 110829 TC cannot be more than surrounding cells, this limits the spikes in deposition and erosion
    if (SwitchLimitTC)
    {
      FOR_CELL_IN_FLOODAREA

          double maxtc = 0;
          double avgtc = 0;
          int count = 0;

          for (int i = 1; i <= 9; i++)
             if(i != 5)
             {
                if ((r+dx[i] >= 0 && c+dy[i] >= 0 && r+dx[i] < _nrRows && c+dy[i] < _nrCols)
                    && !pcr::isMV(TCFlood->data[r+dx[i]][c+dy[i]]))
                {
                   avgtc = avgtc + TCFlood->data[r+dx[i]][c+dy[i]];
                   maxtc = std::max(maxtc,TCFlood->data[r+dx[i]][c+dy[i]]);
                   count++;
                }
             }
          TCFlood->Drc = std::max(TCFlood->Drc, avgtc/count);
          TCFlood->Drc = std::min(TCFlood->Drc, maxtc);
       }
    }

    FOR_CELL_IN_FLOODAREA

       double velocity = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);
       double discharge = velocity * ChannelAdj->Drc * hmx->Drc;
       double watervol = ChannelAdj->Drc *DX->Drc*hmx->Drc;

       //### Add sediment to flood water
       SFlood->Drc += 0;

       //### calc concentration and net transport capacity
       DEP->Drc = 0;
       // init deposition for this timestep
       SCFlood->Drc = MaxConcentration(watervol, SFlood->Drc);
       // limit sed concentration to max

       double maxTC = std::max(TCFlood->Drc - SCFlood->Drc,0.0);
       // positive difference: TC deficit becomes detachment (ppositive)
       double minTC = std::min(TCFlood->Drc - SCFlood->Drc,0.0);
       // negative difference: TC surplus becomes deposition (negative)
       // unit kg/m3

       //### detachment
       double TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
       // detachment can only come from soil, not roads (so do not use flowwidth)
       // units s * m/s * m * m = m3

       DetFlood->Drc = Y->Drc * maxTC * TransportFactor;
       // unit = kg/m3 * m3 = kg
       DetFlood->Drc = std::min(DetFlood->Drc, maxTC * discharge*_dt);
       // cannot have more detachment than remaining capacity in flow
       // use discharge because standing water has no erosion

       if (GrassFraction->Drc > 0)
          DetFlood->Drc = (1-GrassFraction->Drc) * DetFlood->Drc;
       // no flow detachment on grass strips

       // Detachment edxceptions:
       DetFlood->Drc = (1-StoneFraction->Drc) * DetFlood->Drc ;
       // no flow detachment on stony surfaces

       if (SwitchHardsurface)
          DetFlood->Drc = (1-HardSurface->Drc) * DetFlood->Drc ;
       // no flow detachment on hard surfaces

       if (SwitchHouses)
          DetFlood->Drc = (1-HouseCover->Drc)*DetFlood->Drc;
       // no flow det from house roofs

       // IN KG/CELL

       //DETFlow->Drc = (1-Snowcover->Drc) * DETFlow->Drc ;
       /* TODO: CHECK THIS no flow detachment on snow */
       //is there erosion and sedimentation under the snowdeck?

       //### deposition
       if (WH->Drc > MIN_HEIGHT)
          TransportFactor = (1-exp(-_dt*SettlingVelocity->Drc/hmx->Drc)) * watervol;
       else
          TransportFactor = WaterVolall->Drc;
       // if settl velo is very small, transportfactor is 0 and depo is 0
       // if settl velo is very large, transportfactor is 1 and depo is max

       //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
       // deposition can occur on roads and on soil (so use flowwidth)

       double deposition = minTC * TransportFactor;
       // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

       if (SwitchLimitDepTC)
          deposition = std::max(deposition, minTC *watervol);
       // cannot be more than sediment above capacity
       deposition = std::max(deposition, -Sed->Drc);
       // cannot have more depo than sediment present
       /* TODO what about this: which one to choose */

       //deposition = (1-Snowcover->Drc) * deposition;
       /* TODO: TRUE??? no deposition on snow */
       //if (GrassFraction->Drc > 0)
       //   deposition = -Sed->Drc*GrassFraction->Drc + (1-GrassFraction->Drc)*deposition;
       // generate 100% deposition on grassstrips
       //? bit tricky, maximizes effect on grassstrips ?

       DepFlood->Drc += deposition;
       // IN KG/CELL

       //### sediment balance
       SFlood->Drc += DETFlow->Drc;
       SFlood->Drc += deposition;

       SCFlood->Drc = MaxConcentration(watervol, SFlood->Drc);
       // limit concentration to 850 and throw rest in deposition
    }



}

void TWorld::SWOFSedimentFlowWS(int wsnr, double dt)
{
   if (!SwitchErosion)
      return;

   FOR_WATERSHED_ROW_COL(wsnr)

       TCFlood->Drc = 0;
       SCNFlood->Drc = SCFlood->Drc;
       SNFlood->Drc = SFlood->Drc;

       //set concentration from present sediment
       SCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, SFlood->Drc);

   }


   //first calculate the weights for the cells that are closest to location that flow is advected to
   FOR_WATERSHED_ROW_COL(wsnr)

       //no flood velocity means no flood sediment transport, so skip this cell
       if((Vflood->Drc == 0 && Uflood->Drc == 0))
       {
           continue;
       }

       //cell sizes
       double cdx = _dx;
       double cdy = DX->Drc;

       //the sign of the x and y direction of flow
       double yn = signf(Vflood->Drc);
       double xn = signf(Uflood->Drc);

       //should not travel more distance than cell size
       double dsx = xn*std::min(fabs(Uflood->Drc)* dt,cdx);
       double dsy = yn*std::min(fabs(Vflood->Drc)* dt,cdy);

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

           w[i] = weight;



       }
       //normalize: sum of the 4 weights is equal to 1
       double wt = 0.0;
       for (int i=0; i<4; i++)
       {
           wt += w[i];
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
               double cdx2 = _dx;
               double cdy2 = DX->data[r2][c2];

               //weight * the flow is distributed to the ith cell that neighbours the advected flow.
                   SNFlood->data[r2][c2] +=  w[i]* cdx*cdy*hmx->Drc * SCFlood->Drc;
                   SNFlood->data[r][c] -=  w[i]* cdx*cdy*hmx->Drc * SCFlood->Drc;

           }
       }
   }

   FOR_WATERSHED_ROW_COL(wsnr)

       //set concentration from present sediment
       SCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hmx->Drc, SFlood->Drc);

       SCFlood->Drc = SCNFlood->Drc;
       SFlood->Drc = SNFlood->Drc;
   }

}

void TWorld::SWOFSedimentDetWS(int wsnr, double dt)
{

   if (!SwitchErosion)
      return;

   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, -1, -1, -1, 0, 0, 0, 1, 1, 1};

   FOR_WATERSHED_ROW_COL(wsnr)

      tm->Drc = 0;
      double velocity = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);
      double discharge = velocity * _dx * hmx->Drc;
      //### Calc transport capacity
      double omega = 100*velocity*discharge;
      // V in cm/s in this formula assuming grad is SINE
      double omegacrit = 0.4;
      // critical unit streampower in cm/s
      TCFlood->Drc = std::min(MAXCONC, 2650 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));
      // not more than 2650*0.32 = 848 kg/m3
   }

   //VJ 110829 TC cannot be more than surrounding cells, this limits the spikes in deposition and erosion
   if (SwitchLimitTC)
   {
     FOR_WATERSHED_ROW_COL(wsnr)

         double maxtc = 0;
         double avgtc = 0;
         int count = 0;

         for (int i = 1; i <= 9; i++)
            if(i != 5)
            {
               if ((r+dx[i] >= 0 && c+dy[i] >= 0 && r+dx[i] < _nrRows && c+dy[i] < _nrCols)
                   && !pcr::isMV(TCFlood->data[r+dx[i]][c+dy[i]]))
               {
                  avgtc = avgtc + TCFlood->data[r+dx[i]][c+dy[i]];
                  maxtc = std::max(maxtc,TCFlood->data[r+dx[i]][c+dy[i]]);
                  count++;
               }
            }
         TCFlood->Drc = std::max(TCFlood->Drc, avgtc/count);
         TCFlood->Drc = std::min(TCFlood->Drc, maxtc);
      }
   }

   FOR_WATERSHED_ROW_COL(wsnr)

      double velocity = std::sqrt(Uflood->Drc *Uflood->Drc + Vflood->Drc * Vflood->Drc);
      double discharge = velocity * _dx * hmx->Drc;
      double watervol = _dx*DX->Drc*hmx->Drc;

      //### Add sediment to flood water
      SFlood->Drc += 0;

      //### calc concentration and net transport capacity
      DEP->Drc = 0;
      // init deposition for this timestep
      SCFlood->Drc = MaxConcentration(watervol, SFlood->Drc);
      // limit sed concentration to max

      double maxTC = std::max(TCFlood->Drc - SCFlood->Drc,0.0);
      // positive difference: TC deficit becomes detachment (ppositive)
      double minTC = std::min(TCFlood->Drc - SCFlood->Drc,0.0);
      // negative difference: TC surplus becomes deposition (negative)
      // unit kg/m3

      //### detachment
      double TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
      // detachment can only come from soil, not roads (so do not use flowwidth)
      // units s * m/s * m * m = m3

      DetFlood->Drc = Y->Drc * maxTC * TransportFactor;
      // unit = kg/m3 * m3 = kg
      DetFlood->Drc = std::min(DetFlood->Drc, maxTC * discharge*_dt);
      // cannot have more detachment than remaining capacity in flow
      // use discharge because standing water has no erosion

      if (GrassFraction->Drc > 0)
         DetFlood->Drc = (1-GrassFraction->Drc) * DetFlood->Drc;
      // no flow detachment on grass strips

      // Detachment edxceptions:
      DetFlood->Drc = (1-StoneFraction->Drc) * DetFlood->Drc ;
      // no flow detachment on stony surfaces

      if (SwitchHardsurface)
         DetFlood->Drc = (1-HardSurface->Drc) * DetFlood->Drc ;
      // no flow detachment on hard surfaces

      if (SwitchHouses)
         DetFlood->Drc = (1-HouseCover->Drc)*DetFlood->Drc;
      // no flow det from house roofs

      // IN KG/CELL

      //DETFlow->Drc = (1-Snowcover->Drc) * DETFlow->Drc ;
      /* TODO: CHECK THIS no flow detachment on snow */
      //is there erosion and sedimentation under the snowdeck?

      //### deposition
      if (WH->Drc > MIN_HEIGHT)
         TransportFactor = (1-exp(-_dt*SettlingVelocity->Drc/hmx->Drc)) * watervol;
      else
         TransportFactor = WaterVolall->Drc;
      // if settl velo is very small, transportfactor is 0 and depo is 0
      // if settl velo is very large, transportfactor is 1 and depo is max

      //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
      // deposition can occur on roads and on soil (so use flowwidth)

      double deposition = minTC * TransportFactor;
      // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

      if (SwitchLimitDepTC)
         deposition = std::max(deposition, minTC * watervol);
      // cannot be more than sediment above capacity
      deposition = std::max(deposition, -Sed->Drc);
      // cannot have more depo than sediment present
      /* TODO what about this: which one to choose */

      //deposition = (1-Snowcover->Drc) * deposition;
      /* TODO: TRUE??? no deposition on snow */
      //if (GrassFraction->Drc > 0)
      //   deposition = -Sed->Drc*GrassFraction->Drc + (1-GrassFraction->Drc)*deposition;
      // generate 100% deposition on grassstrips
      //? bit tricky, maximizes effect on grassstrips ?

      DepFlood->Drc += deposition;
      // IN KG/CELL

      //### sediment balance
      SFlood->Drc += DETFlow->Drc;
      SFlood->Drc += deposition;

      SCFlood->Drc = MaxConcentration(watervol, SFlood->Drc);
      // limit concentration to 850 and throw rest in deposition
   }



}
