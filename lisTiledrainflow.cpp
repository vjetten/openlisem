
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

#include "model.h"

//---------------------------------------------------------------------------
//fraction of water and sediment flowing from the surface to the tiledrain system
void TWorld::ToTiledrain(void)
{
   if (SwitchIncludeTile)
   {
      FOR_ROW_COL_MV_TILE
      {
         double fractiontotile = qMin(_dt*V->Drc/(0.5*sqrt(_dx*DX->Drc-TileSinkhole->Drc)), 1.0);
         double Volume = WHrunoff->Drc * FlowWidth->Drc * DX->Drc;

         RunoffVolinToTile->Drc = fractiontotile*Volume;
         // water diverted to the channel
         WHrunoff->Drc *= (1-fractiontotile);
         // adjust water height
      }
      CalcVelDisch();
      // recalc velocity and discharge
   }
}
//---------------------------------------------------------------------------
// V, alpha and Q in the Tile
void TWorld::CalcVelDischTile()
{
   FOR_ROW_COL_MV_TILE
	{
      double Perim, Radius, Area;
      const double beta = 0.6;
      const double _23 = 2.0/3.0;
		double beta1 = 1/beta;
		double wh = TileWH->Drc;
		double FW = TileWidth->Drc;
		double grad = sqrt(TileGrad->Drc);

		Perim = FW + 2*wh;
		Area = FW*wh;

		if (Perim > 0)
			Radius = Area/Perim;
		else
			Radius = 0;

		TileAlpha->Drc = pow(TileN->Drc/grad * powl(Perim, _23),beta);

		if (TileAlpha->Drc > 0)
			TileQ->Drc = pow(Area/TileAlpha->Drc, beta1);
		else
			TileQ->Drc = 0;

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

   tm->fill(0);
   tma->fill(0);
   tmb->fill(0);
   // calculate new Tile WH , WidthUp and Volume
	FOR_ROW_COL_MV_TILE
	{
		/*---- Water ----*/

      TileQsn->Drc = 0;
      Tileq->Drc = 0;
      TileQoutflow->Drc = 0;

      TileWaterVol->Drc += RunoffVolinToTile->Drc;
      // add from the surface

      //TileDrainSoil->Drc = min(TileDrainSoil->Drc, TileHeight->Drc );
      // cannot have more water than fits in size
      TileWaterVol->Drc += TileDrainSoil->Drc * TileWidth->Drc * TileDX->Drc;
      // add inflow to Tile, tiledrainsoil is in m per timestep

      TileWH->Drc =  TileDrainSoil->Drc;
      // water height in m per cell

   }

   CalcVelDischTile();


   TileQn->setMV();
   FOR_ROW_COL_MV_TILE
   {
      if (LDDTile->Drc == 5)
      {
         Kinematic(r,c, LDDTile, TileQ, TileQn, TileQs, TileQsn, Tileq, TileAlpha, TileDX,
                   TileWaterVol, tm, tma, tmb, SubsMaps);
      }
   }

//   TileQoutflow->DrcOutlet = TileQn->DrcOutlet * _dt;
   FOR_ROW_COL_MV_TILE
         TileQoutflow->Drc = TileQn->Drc * _dt;

   TileQn->cover(LDD, 0); // avoid missing values around Tile for adding to Qn for output
   TileQs->cover(LDD, 0);

   FOR_ROW_COL_MV_TILE
   {
      double TileArea = TileAlpha->Drc*pow(TileQn->Drc, 0.6);
      InfilVolKinWave->Drc += Tileq->Drc*_dt + TileWaterVol->Drc - (TileArea * TileDX->Drc) - TileQn->Drc*_dt;
      // diff is a small error in this case added to kin wave infil!

      TileWaterVol->Drc = TileArea * TileDX->Drc;
      // total water vol after kin wave in m3, going to the next timestep
   }   
}
//---------------------------------------------------------------------------

