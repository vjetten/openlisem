/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
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
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/

/*!
  \file lisTiledrainflow.cpp
  \brief calculate tile drain system flow as a kinematic wave, no sediment functions

functions: \n
- void TWorld::ToTiledrain(void) \n
- void TWorld::CalcVelDischTile() \n
- void TWorld::TileFlow(void)\n
 */

#include <algorithm>
#include "model.h"
#include "operation.h"

//TODO convert flow to linked list

//---------------------------------------------------------------------------
//fraction of water and sediment flowing from the surface to the tiledrain system
void TWorld::ToTiledrain()//int thread)
{

    if (SwitchIncludeTile || SwitchIncludeStormDrains)
    {
        //Fill(*RunoffVolinToTile,0);
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            RunoffVolinToTile->Drc = 0;
        }}

       //#pragma omp parallel for collapse(2) num_threads(userCores)
        FOR_ROW_COL_MV_TILE {
            if(TileSinkhole->Drc > 0)
            {
                double fractiontotile = std::max(1.0, std::min(0.0,TileSinkhole->Drc/(_dx*DX->Drc)));
                // fraction based on surface, simpel! Street inlet is assumed to be a hole in the street

                double MaxVol;
                if (SwitchStormDrainShape)
                    MaxVol = DX->Drc * PI*TileDiameter->Drc*TileDiameter->Drc*0.25; //(pi r^2)
                else
                    MaxVol = DX->Drc * TileDiameter->Drc; //rectangular drain

                if (TileWaterVol->Drc > MaxVol*0.99)
                    fractiontotile = 0;
                else {
                    double dh = (MaxVol - TileWaterVol->Drc)/CHAdjDX->Drc;

                    dh = std::min(dh, fractiontotile*WHrunoff->Drc);

                    RunoffVolinToTile->Drc = dh * FlowWidth->Drc * DX->Drc;

                    // water diverted to the channel
                    WHrunoff->Drc -= dh; //*= (1-fractiontotile);
                    // adjust water height

                    WH->Drc = WHrunoff->Drc + WHstore->Drc;

                    WaterVolall->Drc = WHrunoff->Drc*CHAdjDX->Drc + MicroStoreVol->Drc;
                }
            }
        }
    }
}
//---------------------------------------------------------------------------
// V, alpha and Q in the Tile
void TWorld::CalcVelDischRectangular()
{
    double Perim, Area;
    const double beta = 0.6;
    const double _23 = 2.0/3.0;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_TILE {
        double wh = TileWaterVol->Drc/TileWidth->Drc;
        Perim = TileWidth->Drc + 2*wh;
        Area = TileWidth->Drc * wh;
        double grad = sqrt(TileGrad->Drc);
        double TileV_ = powl(Area/Perim,_23) * grad/TileN->Drc;
        //TileV_ = std::min(TileV_, 2.0);
        //limit velocity to 2 m/s?
        TileQ->Drc = Area*TileV_;

        TileAlpha->Drc = qPow(TileN->Drc/grad * powl(Perim, _23),beta);

    }
}
//---------------------------------------------------------------------------
// V, alpha and Q in the Tile
//https://www.engineersedge.com/fluid_flow/partially_full_pipe_flow_calculation/partiallyfullpipeflow_calculation.htm
//https://www.ajdesigner.com/phphydraulicradius/hydraulic_radius_equation_pipe.php
// Neweton iteration to derive drain water height
void TWorld::CalcVelDischCircular()
{
#pragma omp parallel for num_threads(userCores)
   FOR_ROW_COL_MV_TILE {

      double gradN = sqrt(TileGrad->Drc)/TileN->Drc;
      double rr = TileDiameter->Drc/2;
      double Perim, K, theta;

      double Area = TileWaterVol->Drc / DX->Drc;

      theta = PI;
      double fx, Fx;
      double Ar = 2*Area/(rr*rr);
      for (int k = 0 ; k < 20; k++) {
                fx = 1-cos(theta);
                Fx = -sin(theta) + theta - Ar;
                theta = fx > 0 ? theta - Fx/fx : 0.0;
                if( Fx < 1e-6)
                    break;
      }
      // newton rapson iteration to get the water height in a circular pipe for the wet perimeter

      K = rr*rr*(theta-sin(theta))*0.5;
      double TWH = rr - cos(theta/2.0)*rr;
      Perim = Area < 0.5*rr*rr*PI ? rr*theta : PI*TileDiameter->Drc-rr*theta;

//      double minV = std::pow(TileWH->Drc,2.0/3.0) * gradN;
      double minV = std::pow(TWH,2.0/3.0) * gradN;
      double TileV_ = Perim > 1e-10 ? std::min(minV, std::pow(Area/Perim,2.0/3.0) * gradN) : 0.0;
      //TileV_ = std::min(TileV_, 2.0);
      //limit velocity to 2 m/s?
      TileQ->Drc = Area*TileV_;

      TileAlpha->Drc = std::pow(std::pow(Perim, 2.0/3.0)/gradN , 0.6);
   }
}
//---------------------------------------------------------------------------
//- calc Tileflow, Tileheight, kin wave
void TWorld::TileFlow(void)
{
   if (!SwitchIncludeTile)
      return;

   #pragma omp parallel for num_threads(userCores)
   FOR_ROW_COL_MV_TILE {
      TileQn->Drc = 0;

      TileWaterVol->Drc += TileDrainSoil->Drc * TileDiameter->Drc * DX->Drc;
      // asume water can come from all sides!
      // add inflow to Tile in m3, tiledrainsoil is in m per timestep

      TileWaterVolSoil->Drc += TileDrainSoil->Drc * TileDiameter->Drc  * DX->Drc;
      // soil only used for MB correction

      //TileWaterVol->Drc += RunoffVolinToTile->Drc;
      // add from the surface ????

     // TileWH->Drc = TileWaterVol->Drc/(TileWidth->Drc * TileDX->Drc);
      // water height in m per cell

   }

   if (SwitchStormDrainShape)
      CalcVelDischCircular();
   else
      CalcVelDischRectangular();

   TileQn->setAllMV();

   Fill(*QinKW, 0.0);
   // flag all new flux as missing value, needed in kin wave and replaced by new flux
   FOR_ROW_COL_MV_TILE {
      if (LDDTile->Drc == 5)
                Kinematic(r,c, LDDTile, TileQ, TileQn, TileAlpha, DX);
   }

   cover(*TileQn, *LDD, 0); // avoid missing values around Tile for adding to Qn for output

   #pragma omp parallel for num_threads(userCores)
   FOR_ROW_COL_MV_TILE {
        TileWaterVol->Drc = TileWaterVol->Drc + _dt*(QinKW->Drc - TileQn->Drc);
      TileWaterVol->Drc =  std::max(0.0, TileWaterVol->Drc);
        TileWaterVol->Drc =  std::min(TileWaterVol->Drc, TileDiameter->Drc*DX->Drc);
        TileQ->Drc = TileQn->Drc;
   }

}
//---------------------------------------------------------------------------
//- calc Tileflow, Tileheight, kin wave
void TWorld::StormDrainFlow(void)
{
    if (!SwitchIncludeStormDrains)
        return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV {
        TileWaterVol->Drc += RunoffVolinToTile->Drc;
        // add water from the surface
    }

    if (SwitchStormDrainShape)
        CalcVelDischCircular();
    else
        CalcVelDischRectangular();
    // calc Q, V Aplha for circular drain

    TileQn->setAllMV();
    Fill(*QinKW, 0.0);
    // flag all new flux as missing value, needed in kin wave and replaced by new flux
    FOR_ROW_COL_MV_TILE
    {
        if (LDDTile->Drc == 5)
            Kinematic(r,c, LDDTile, TileQ, TileQn, TileAlpha, DX);
    }

    cover(*TileQn, *LDD, 0); // avoid missing values around Tile for adding to Qn for output
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_TILE
    {
        TileWaterVol->Drc = TileWaterVol->Drc + _dt*(QinKW->Drc - TileQn->Drc);
        TileWaterVol->Drc =  std::max(0.0, TileWaterVol->Drc);
        TileWaterVol->Drc =  std::min(TileWaterVol->Drc, TileDiameter->Drc*DX->Drc);
        TileQ->Drc = TileQn->Drc;

//        // total water vol after kin wave in m3, going to the next timestep
//        double Area = TileWaterVol->Drc/TileDX->Drc;

//        double rr = TileDiameter->Drc/2.0;
//        double theta = PI;
//        double fx, Fx;
//        double Ar = 2*Area/(rr*rr);

//        for (int k = 0 ; k< 20; k++) {
//            fx = 1-cos(theta);
//            Fx = -sin(theta) + theta - Ar;
//            theta = fx > 0 ? theta - Fx/fx : 0.0;
//            if( Fx < 1e-6)
//                break;
//        }

//        TileWH->Drc = rr - cos(theta/2.0)*rr;
//        double minV = std::pow(TileWH->Drc, 2.0/3.0)*sqrt(TileGrad->Drc)/TileN->Drc;
//        TileV->Drc = Area < 1E-10 ?  minV : TileQn->Drc/Area;

    }



//   #pragma omp parallel for num_threads(userCores)
//   FOR_ROW_COL_MV {
//      TileWaterVol->Drc += RunoffVolinToTile->Drc;
//   }

//   CalcVelDischDrain();
//   // calc Q, V Aplha for circular drain

//   Fill(*QinKW, 0.0);
//   Fill(*TileQn, 0.0);

//   upstreamDrain(LDDTile, TileMaxQ, TileQ, QinKW);
//   //cover(*TileQn, *LDD, 0); // avoid missing values around Tile for adding to Qn for output

//   #pragma omp parallel for num_threads(userCores)
//   FOR_ROW_COL_MV_TILE {
//         double MaxVol = DX->Drc * PI*TileDiameter->Drc*TileDiameter->Drc*0.25; //(pi r^2)
//         TileWaterVol->Drc = TileWaterVol->Drc + _dt*(QinKW->Drc - TileQ->Drc);
//         TileWaterVol->Drc = std::min(MaxVol, TileWaterVol->Drc);
//         TileWaterVol->Drc = std::max(0.0, TileWaterVol->Drc);
//         // total water vol after kin wave in m3, going to the next timestep
//         TileQn->Drc = TileQ->Drc;
//   }
//   cover(*TileQn, *LDD, 0);
}

