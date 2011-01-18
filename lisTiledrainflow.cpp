/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt
website SVN: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
Functionality in lisTileflow.cpp:
- calc alpha and Q in the Tile
- calc Tileflow, Tileheight, do kin wave
*/

#include "model.h"

//---------------------------------------------------------------------------
// V, alpha and Q in the Tile
void TWorld::CalcVelDischTile()
{
	FOR_ROW_COL_MV_CH
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

		TileQsn->Drc =0;
		Tileq->Drc =0;
		TileQoutflow->Drc =0;
		TileWH->Drc = 0;

      TileWaterVol->Drc += TileDrainSoil->Drc;  // IS VOLUME ????
		// add inflow to Tile

		TileWH->Drc = TileWaterVol->Drc/(TileWidth->Drc*DX->Drc);
   }

   CalcVelDischTile();


   TileQn->setMV();
   FOR_ROW_COL_MV_TILE
   {
      if (LDDTile->Drc == 5)
      {
         Kinematic(r,c, LDDTile, TileQ, TileQn, TileQs, TileQsn, Tileq, TileAlpha, DX,
                   TileWaterVol, tm, tma, tmb);

         TileQoutflow->Drc = TileQn->Drc * _dt;
      }
   }
   TileQn->cover(0); // avoid missing values around Tile for adding to Qn for output
   TileQs->cover(0);

   FOR_ROW_COL_MV_TILE
   {
      double TileArea = TileAlpha->Drc*pow(TileQn->Drc, 0.6);
      TileWH->Drc = TileArea/TileWidth->Drc;
      // water height is not used except for output i.e. watervolume is cycled

      TileWaterVol->Drc = TileArea * DX->Drc;
      // total water vol after kin wave in m3, going to the next timestep
   }
}
//---------------------------------------------------------------------------

