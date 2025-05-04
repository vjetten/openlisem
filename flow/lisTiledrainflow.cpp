/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
  \file lisTiledrainflow.cpp
  \brief calculate tile drain system flow as a kinematic wave, no sediment functions

functions: \n

 */

#include <algorithm>
#include "model.h"
#include "operation.h"

//---------------------------------------------------------------------------
// flow in all road cells to tiledrain
//fraction of water and sediment flowing from the surface to the tiledrain system
void TWorld::ToTiledrainAll()
{
    if (SwitchIncludeStormDrains)  //SwitchIncludeTile ||
    {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_TILEL {
          RunoffVolinToTile->Drc = 0;
          double fractiontotile = 0;;
          double MaxVol = DX->Drc*TileArea->Drc; //(pi r^2 or heightxwidth, done in datainit

          if (TileWaterVol->Drc >= MaxVol)
            fractiontotile = 0;
          else {
            fractiontotile = 0.03/(RoadWidthDX->Drc*DX->Drc)*(DX->Drc/TileDrainDistance);
            fractiontotile = std::max(0.0, std::min(1.0,fractiontotile));
           // qDebug() << fractiontotile;
            // every tile cell has a subtraction of water, based on the inlet fraction in the street
            // assumed entry is 0.1 * 0.3 m
            // if a road is divided over more cells, this probably overewstimates the entrance

            double dh = fractiontotile*WHrunoff->Drc;
            if (FloodDomain ->Drc > 0) {
                dh = fractiontotile*hmxrunoff->Drc;
            }

            RunoffVolinToTile->Drc = dh*CHAdjDX->Drc;

            // adjust water height
            WaterVolall->Drc -= dh*CHAdjDX->Drc;
            if (FloodDomain->Drc == 0) {
                WHrunoff->Drc -= dh;
                WH->Drc -= dh;
            } else {
                hmxrunoff->Drc -= dh;
                hmx->Drc -= dh;
            }
            hmxWH->Drc = WH->Drc + hmx->Drc;
          }
        }}
    }
}
//---------------------------------------------------------------------------
// V, alpha and Q in the Tile
void TWorld::CalcVelDischRectangular()
{
    double Perim, Area;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_TILEL {

        Area = TileWaterVol->Drc/DX->Drc;
        Perim = TileWidth->Drc + Area/TileWidth->Drc; //(=w+2*h)
        //TileA->Drc = Area;
        TileMaxQ->Drc = Area*pow(Area/Perim,2.0/3.0) * sqrt(TileGrad->Drc)/TileN->Drc;;
        TileAlpha->Drc  = Area/std::pow(TileQ->Drc, BETArect);
    }}
}
//---------------------------------------------------------------------------
// V, alpha and Q in the Tile
// https://www.engineersedge.com/fluid_flow/partially_full_pipe_flow_calculation/partiallyfullpipeflow_calculation.htm
// https://www.ajdesigner.com/phphydraulicradius/hydraulic_radius_equation_pipe.php
// Neweton iteration to derive drain water height
void TWorld::CalcVelDischCircular()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_TILEL {
      double Area = TileWaterVol->Drc / DX->Drc;
      //TileA->Drc = Area;
      double a = Area/TileArea->Drc;
      double perim = 0;
      if (a < 1) {
          double theta_next;
          double theta = PI;
          double tol = 1e-6;
          // get angle theta from a
          for (int j = 0; j < 50; j++ ) {
              double f = (theta - sin(theta)) / (2 * PI) - a;
              double df = (1 - cos(theta)) / (2 * PI);
              theta_next = theta - f / df;
              if (abs(theta_next - theta) < tol)
                  break;
              theta = theta_next;
          }
          perim = TileDiameter->Drc/2.0*theta_next; // P = r*theta; A =
      } else {
          perim = TileDiameter->Drc;
      }

      if (perim < 1e-6)
          TileQ->Drc = 0;
      else
          TileQ->Drc = std::pow(Area/perim, 5.0/3.0) * sqrt(TileGrad->Drc)/TileN->Drc;

      TileAlpha->Drc = std::pow(std::pow(perim, 2.0/3.0)*TileN->Drc/sqrt(TileGrad->Drc), 0.6);
   }}
}
//---------------------------------------------------------------------------
//- calc Tileflow, Tileheight, kin wave
void TWorld::TileFlow(void)
{
   if (!SwitchIncludeTile && !SwitchIncludeStormDrains)
      return;

   // get water from surface
   if (SwitchIncludeStormDrains) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_TILEL {
            TileWaterVol->Drc += RunoffVolinToTile->Drc;
                // add water from the surface
        }}
   }

    // get water from soil
   if (SwitchIncludeTile) {
       #pragma omp parallel for num_threads(userCores)
       FOR_ROW_COL_MV_TILEL {
          TileWaterVol->Drc += TileWaterVolSoil->Drc;
       }}
   }

   if (SwitchDrainCircular)
      CalcVelDischCircular();
   else
      CalcVelDischRectangular();

  #pragma omp parallel for num_threads(userCores)
  FOR_ROW_COL_MV_L {
    TileQn->Drc = 0;
  }}

  // #pragma omp parallel for ordered num_threads(userCores)
  // parallel doesn't work here because the order of cells has to be maintained
  for(long i_ =  0; i_ < crlinkedlddtile_.size(); i_++) {
      int r = crlinkedlddtile_.at(i_).r;
      int c = crlinkedlddtile_.at(i_).c;
      double Qin = 0;

      if (crlinkedlddtile_.at(i_).nr > 0) {
          for(int j = 0; j < crlinkedlddtile_.at(i_).nr; j++) {
              int rr = crlinkedlddtile_.at(i_).inn[j].r;
              int cr = crlinkedlddtile_.at(i_).inn[j].c;
              Qin += TileQn->Drcr;
          }
      }

      Qin = std::min(Qin, TileMaxQ->Drc);
      TileQn->Drc = IterateToQnew(Qin, TileQ->Drc, TileAlpha->Drc, _dt, DX->Drc, TileMaxQ->Drc, TileMaxAlpha->Drc);
      TileQn->Drc = std::min(Qin+TileWaterVol->Drc/_dt, TileQn->Drc);
      TileQn->Drc = std::min(TileQn->Drc, TileMaxQ->Drc);

      TileWaterVol->Drc = TileWaterVol->Drc + _dt*(Qin - TileQn->Drc);
      TileWaterVol->Drc = std::max(0.0, TileWaterVol->Drc);
      TileWaterVol->Drc = std::min(TileWaterVol->Drc, TileArea->Drc * DX->Drc);
  }
}
//---------------------------------------------------------------------------
