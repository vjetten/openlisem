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

//---------------------------------------------------------------------------
//fraction of water and sediment flowing from the surface to the tiledrain system
void TWorld::ToTiledrain()//int thread)
{

    if (SwitchIncludeTile || SwitchIncludeStormDrains)
    {
        fill(*RunoffVolinToTile,0);
        FOR_ROW_COL_MV_TILE {
            if(TileSinkhole->Drc > 0)
            {
                double fractiontotile = std::max(1.0, std::min(0.0,TileSinkhole->Drc/(_dx*DX->Drc)));
                // fraction based on surface, simpel! Street inlet is assumed to be a hole in the street

                double MaxVol = DX->Drc * PI*TileDiameter->Drc*TileDiameter->Drc*0.25;
                if (TileWaterVol->Drc > MaxVol*0.95)
                    fractiontotile = 0;

                RunoffVolinToTile->Drc = fractiontotile*WHrunoff->Drc * FlowWidth->Drc * DX->Drc;
                // water diverted to the channel
                WHrunoff->Drc *= (1-fractiontotile);
                // adjust water height

                WH->Drc = WHrunoff->Drc + WHstore->Drc;
                WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);

            }
        }
    }
}
//---------------------------------------------------------------------------
// V, alpha and Q in the Tile
void TWorld::CalcVelDischTile()
{
    double Perim, Area;
    const double beta = 0.6;
    const double _23 = 2.0/3.0;
    double beta1 = 1/beta;

   FOR_ROW_COL_MV_TILE
   {
      double wh = TileWH->Drc;
      Perim = TileWidth->Drc + 2*wh;
      Area = TileWidth->Drc*wh;
      double grad = sqrt(TileGrad->Drc);

      TileAlpha->Drc = qPow(TileN->Drc/grad * powl(Perim, _23),beta);
      TileV->Drc = powl(Area/Perim,_23) * grad/TileN->Drc;
      TileQ->Drc = Area*TileV->Drc;

   }
   else
   {
      TileAlpha->Drc = 0;
      TileQ->Drc = 0;
   }
}
//---------------------------------------------------------------------------
//- calc Tileflow, Tileheight, kin wave
void TWorld::TileFlow(void)
{
   if (!SwitchIncludeTile)
      return;

   FOR_ROW_COL_MV_TILE
   {
      /*---- Water ----*/

      TileQsn->Drc = 0;
      Tileq->Drc = 0;
      //TileQoutflow->Drc = 0;

      //TileDrainSoil->Drc = std::min(TileDrainSoil->Drc, TileHeight->Drc );
      // cannot have more water than fits in size
      TileWaterVol->Drc += TileDrainSoil->Drc * TileWidth->Drc * TileDX->Drc;
      // add inflow to Tile in m3, tiledrainsoil is in m per timestep

      TileWaterVolSoil->Drc += TileDrainSoil->Drc * TileWidth->Drc * TileDX->Drc;
      // soil only used for MB correction

      TileWaterVol->Drc += RunoffVolinToTile->Drc;
      // add from the surface

      TileWH->Drc = TileWaterVol->Drc/(TileWidth->Drc * TileDX->Drc);
      // water height in m per cell

   }

   CalcVelDischTile();

   TileQn->setAllMV();

   fill(*QinKW, 0.0);
   // flag all new flux as missing value, needed in kin wave and replaced by new flux
   FOR_ROW_COL_MV_TILE
   {
      if (LDDTile->Drc == 5)
         Kinematic(r,c, LDDTile, TileQ, TileQn, Tileq, TileAlpha, TileDX, TileMaxQ);
   }

   cover(*TileQn, *LDD, 0); // avoid missing values around Tile for adding to Qn for output

   FOR_ROW_COL_MV_TILE
   {
       TileWaterVol->Drc = TileWaterVol->Drc + _dt*(QinKW->Drc - TileQn->Drc);
       TileWaterVol->Drc =  TileWaterVol->Drc < 0 ? 0 : TileWaterVol->Drc;

//      double TileArea = TileAlpha->Drc*pow(TileQn->Drc, 0.6);
//      InfilVolKinWave->Drc += QinKW->Drc*_dt + TileWaterVol->Drc - (TileArea * TileDX->Drc) - TileQn->Drc*_dt;
      // diff is a small error in this case added to kin wave infil!

//      TileWaterVol->Drc = TileArea * TileDX->Drc;
      // total water vol after kin wave in m3, going to the next timestep
       TileQmax->Drc = std::max(TileQn->Drc, TileQmax->Drc);
   }

}
//---------------------------------------------------------------------------
// V, alpha and Q in the Tile
//https://www.engineersedge.com/fluid_flow/partially_full_pipe_flow_calculation/partiallyfullpipeflow_calculation.htm
//https://www.ajdesigner.com/phphydraulicradius/hydraulic_radius_equation_pipe.php
// Neweton iteration to derive drain water height
void TWorld::CalcVelDischDrain()
{
    FOR_ROW_COL_MV_TILE {

        double gradN = sqrt(TileGrad->Drc)/TileN->Drc;
        double rr = TileDiameter->Drc/2;
        double  Perim, K, theta;

        double Area = TileWaterVol->Drc / TileDX->Drc;

        theta = PI;
        double fx, Fx;
        double Ar = 2*Area/(rr*rr);
        bool stop = false;
        for (int k = 0 ; k< 20; k++) {
            fx = 1-cos(theta);
            Fx = -sin(theta) + theta - Ar;
            theta = fx > 0 ? theta - Fx/fx : 0.0;
            if( Fx < 1e-6)
                break;
        }

        K = rr*rr*(theta-sin(theta))*0.5;
        TileWH->Drc = rr - cos(theta/2.0)*rr;

        Perim = Area < 0.5*rr*rr*PI ? rr*theta : PI*TileDiameter->Drc-rr*theta;

        double minV = std::pow(TileWH->Drc,2.0/3.0) * gradN;
        TileV->Drc = Perim > 1e-10 ? std::min(minV, std::pow(Area/Perim,2.0/3.0) * gradN) : 0.0;

        TileQ->Drc = Area*TileV->Drc;
        TileAlpha->Drc = std::pow(std::pow(Perim, 2.0/3.0)/gradN , 0.6);

//qDebug() << halffull << hoi << theta <<  Perim << Area << TileWH->Drc << TileV->Drc << TileQ->Drc   ;
    }
}
//---------------------------------------------------------------------------
//- calc Tileflow, Tileheight, kin wave
void TWorld::StormDrainFlow(void)
{
   if (!SwitchIncludeStormDrains)
      return;
   FOR_ROW_COL_MV {
         TileWaterVol->Drc += RunoffVolinToTile->Drc;
      // add water from the surface
   }

   CalcVelDischDrain();
   // calc Q, V Aplha for circular drain

   TileQn->setAllMV();
   fill(*QinKW, 0.0);
   // flag all new flux as missing value, needed in kin wave and replaced by new flux
   FOR_ROW_COL_MV_TILE
   {
      if (LDDTile->Drc == 5)
         Kinematic(r,c, LDDTile, TileQ, TileQn, Tileq, TileAlpha, TileDX, TileMaxQ);
   }

   cover(*TileQn, *LDD, 0); // avoid missing values around Tile for adding to Qn for output

   FOR_ROW_COL_MV_TILE
   {
      TileWaterVol->Drc = TileWaterVol->Drc + _dt*(QinKW->Drc - TileQn->Drc);
      TileWaterVol->Drc =  std::max(0.0, TileWaterVol->Drc);
      TileQ->Drc = TileQn->Drc;

      // total water vol after kin wave in m3, going to the next timestep
      double Area = TileWaterVol->Drc/TileDX->Drc;

      double rr = TileDiameter->Drc/2.0;
      double theta = PI;
      double fx, Fx;
      double Ar = 2*Area/(rr*rr);

      for (int k = 0 ; k< 20; k++) {
          fx = 1-cos(theta);
          Fx = -sin(theta) + theta - Ar;
          theta = fx > 0 ? theta - Fx/fx : 0.0;
          if( Fx < 1e-6)
              break;
      }

      TileWH->Drc = rr - cos(theta/2.0)*rr;
      double minV = std::pow(TileWH->Drc, 2.0/3.0)*sqrt(TileGrad->Drc)/TileN->Drc;
      TileV->Drc = Area < 1E-10 ?  minV : TileQn->Drc/Area;

   }
}

