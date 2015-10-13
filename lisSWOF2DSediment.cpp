

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
void TWorld::SWOFSedimentFlow(double dt)
{
    if (!SwitchErosion)
       return;

    FOR_ROW_COL_MV
    {
        TCFlood->Drc = 0;
        SCNFlood->Drc = SCFlood->Drc;
        SNFlood->Drc = SFlood->Drc;

        //set concentration from present sediment
        SCFlood->Drc = MaxConcentration(_dx*DX->Drc*hmx->Drc, SFlood->Drc);

    }


    //first calculate the weights for the cells that are closest to location that flow is advected to
    FOR_CELL_IN_FLOODAREA

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

    FOR_ROW_COL_MV
    {
        //set concentration from present sediment
        SCFlood->Drc = MaxConcentration(_dx*DX->Drc*hmx->Drc, SFlood->Drc);

        SCFlood->Drc = SCNFlood->Drc;
        SFlood->Drc = SNFlood->Drc;
    }

}

void TWorld::SWOFSedimentDet(double dt)
{

    if (!SwitchErosion)
       return;

    int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dy[10] = {0, -1, -1, -1, 0, 0, 0, 1, 1, 1};

    FOR_CELL_IN_FLOODAREA

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
       double discharge = velocity * _dx * hmx->Drc;
       double watervol = _dx*DX->Drc*hmx->Drc;

       //### Add sediment to flood water
       SFlood->Drc += 0;

       //### calc concentration and net transport capacity
       DEP->Drc = 0;
       // init deposition for this timestep
       SCFlood->Drc = MaxConcentration(_dx*DX->Drc*hmx->Drc, SFlood->Drc);
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
          TransportFactor = (1-exp(-_dt*SettlingVelocity->Drc/hmx->Drc)) * _dx*DX->Drc * hmx->Drc;
       else
          TransportFactor = WaterVolall->Drc;
       // if settl velo is very small, transportfactor is 0 and depo is 0
       // if settl velo is very large, transportfactor is 1 and depo is max

       //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
       // deposition can occur on roads and on soil (so use flowwidth)

       double deposition = minTC * TransportFactor;
       // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

       if (SwitchLimitDepTC)
          deposition = std::max(deposition, minTC * _dx*DX->Drc*hmx->Drc);
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

       SCFlood->Drc = MaxConcentration(_dx*DX->Drc*hmx->Drc, SFlood->Drc);
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
       SCFlood->Drc = MaxConcentration(_dx*DX->Drc*hmx->Drc, SFlood->Drc);

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
       SCFlood->Drc = MaxConcentration(_dx*DX->Drc*hmx->Drc, SFlood->Drc);

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
      SCFlood->Drc = MaxConcentration(_dx*DX->Drc*hmx->Drc, SFlood->Drc);
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
         TransportFactor = (1-exp(-_dt*SettlingVelocity->Drc/hmx->Drc)) * _dx*DX->Drc * hmx->Drc;
      else
         TransportFactor = WaterVolall->Drc;
      // if settl velo is very small, transportfactor is 0 and depo is 0
      // if settl velo is very large, transportfactor is 1 and depo is max

      //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
      // deposition can occur on roads and on soil (so use flowwidth)

      double deposition = minTC * TransportFactor;
      // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

      if (SwitchLimitDepTC)
         deposition = std::max(deposition, minTC * _dx*DX->Drc*hmx->Drc);
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

      SCFlood->Drc = MaxConcentration(_dx*DX->Drc*hmx->Drc, SFlood->Drc);
      // limit concentration to 850 and throw rest in deposition
   }



}
