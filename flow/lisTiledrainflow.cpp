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


#include "model.h"


//---------------------------------------------------------------------------
// flow in all road cells to tiledrain, according to fraction road and fraction openings in surface
// pressure flow through the inlet hole
// no sediment
void TWorld::ToTiledrain()
{
    if (SwitchIncludeStormDrains)  {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_TILEL {
            RunoffVolinToTile->Drc = 0;
            double MaxVol = DX->Drc*TileArea->Drc;

            if (TileWaterVol->Drc < MaxVol && WHrunoff->Drc > 1e-3) {

                //double volin = WaterVolall->Drc * TileDrainSize/CHAdjDX->Drc * (DX->Drc/TileDrainDistance);//* RoadWidthHSDX->Drc/_dx
                // fraction of volume, not used

                double volin = _dt*qSqrt(2*GRAV*WHrunoff->Drc)*TileDrainSize * (DX->Drc/TileDrainDistance) * RoadWidthHSDX->Drc/_dx;
                // Bernouilly flow s*m/s*m2=m3 through a hole * fraction of distance compared to cellsize

                // qDebug() << volin <<  WaterVolall->Drc << RunoffVolinToTile->Drc;
                RunoffVolinToTile->Drc = qMin(volin, WaterVolall->Drc-MicroStoreVol->Drc);//RunoffVolinToTile->Drc);
                RunoffVolinToTile->Drc = qMin(MaxVol - TileWaterVol->Drc, RunoffVolinToTile->Drc);
                WaterVolall->Drc -= RunoffVolinToTile->Drc;

                WH->Drc = WaterVolall->Drc/CHAdjDX->Drc;
                WHrunoff->Drc = qMax(0.0, WH->Drc-WHstore->Drc);
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
        if (Area >= TileArea->Drc) {
            TileQ->Drc = TileMaxQ->Drc;
            TileAlpha->Drc = TileMaxAlpha->Drc;
        } else {
            Perim = TileWidth->Drc + 2.0*Area/TileWidth->Drc; //(=w+2*h)
            TileQ->Drc = Area*pow(Area/Perim,2.0/3.0) * sqrt(TileGrad->Drc)/TileN->Drc;;
//            TileAlpha->Drc  = Area/std::pow(TileQ->Drc, BETArect); // gives nan when tileq = 0
            TileAlpha->Drc = std::pow(std::pow(Perim, 2.0/3.0)*TileN->Drc/sqrt(TileGrad->Drc), 0.6);
        }

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
      //  TileN->Drc = 0.1;
      double Area = TileWaterVol->Drc / DX->Drc;
      // if (Area >= 0.999*TileArea->Drc) {
      //     TileQ->Drc = TileMaxQ->Drc;
      //     TileAlpha->Drc = TileMaxAlpha->Drc;
      // } else {
          double a = Area/TileArea->Drc;
          double theta = pipeThetafroma(r,c,a);
          double Perim = 0.5*TileDiameter->Drc*theta;

          if (Perim < 1e-12)
              TileQ->Drc = 0;
          else
              TileQ->Drc = std::pow(Area/Perim, 5.0/3.0) * sqrt(TileGrad->Drc)/TileN->Drc;
          TileQ->Drc = qMin(TileQ->Drc, TileMaxQ->Drc);
          TileAlpha->Drc = std::pow(std::pow(Perim, 2.0/3.0)*TileN->Drc/sqrt(TileGrad->Drc), 0.6);
          TileAlpha->Drc = qMin(TileAlpha->Drc, TileMaxAlpha->Drc);

          //TileAlpha->Drc  = Area/std::pow(TileQ->Drc, BETAcirc);
   //  }
   }}
}
//---------------------------------------------------------------------------
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
    FOR_ROW_COL_MV_TILEL {
        TileQn->Drc =0;
    }}

    int full = 0;

    //  return;
    // wothout fluxes no MB errror anyway!

  //  double tot = MapTotal(*TileWaterVol);
  //  double totq = 0;

    Fill(*tmc, 0);
    for(long i_ =  0; i_ < crlinkedlddtile_.size(); i_++) {
        int r = crlinkedlddtile_.at(i_).r;
        int c = crlinkedlddtile_.at(i_).c;
        double Qin = 0;
        double volMax = TileArea->Drc*DX->Drc;

        int NR = crlinkedlddtile_.at(i_).nr;
        if (NR > 0) {
            for(int j = 0; j < NR; j++) {
                int rr = crlinkedlddtile_.at(i_).inn[j].r;
                int cr = crlinkedlddtile_.at(i_).inn[j].c;
                Qin += TileQn->Drcr;
                // total inflow from incoming brranches
            }

            // if total inflow causes vol > max volume, adjust inflow incoming TileQn
            if (TileWaterVol->Drc+_dt*(Qin-TileQ->Drc) >= volMax) {
                double maxq = qMin(TileMaxQ->Drc, (volMax - TileWaterVol->Drc)/_dt + TileQ->Drc);

                for(int j = 0; j < NR; j++) {
                    int rr = crlinkedlddtile_.at(i_).inn[j].r;
                    int cr = crlinkedlddtile_.at(i_).inn[j].c;
                    TileQn->Drcr = maxq * TileQn->Drcr/Qin;
                    // incoming TileQn is a fraction of maxq
                }
                Qin = maxq;
            }
        }
        tmc->Drc = Qin;

        TileQn->Drc = IterateToQnew(Qin, TileQ->Drc, TileAlpha->Drc, _dt, DX->Drc, TileMaxQ->Drc, TileMaxAlpha->Drc);
        TileQn->Drc = qMin(Qin+TileWaterVol->Drc/_dt, TileQn->Drc);
        TileQn->Drc = qMin(TileQn->Drc, TileMaxQ->Drc);
    }

    #pragma omp parallel for ordered num_threads(userCores)
    FOR_ROW_COL_MV_TILEL {
        TileWaterVol->Drc = TileWaterVol->Drc + _dt*(tmc->Drc - TileQn->Drc);
        TileWaterVol->Drc = qMax(0.0, TileWaterVol->Drc);
        if (TileWaterVol->Drc >= TileArea->Drc*DX->Drc)
            full+=1;
        //TileWaterVol->Drc = qMin(TileWaterVol->Drc, TileArea->Drc*DX->Drc);
        // gives always MB errors!
        //if (LDDTile->Drc == 5)
          //  totq = TileQn->Drc*_dt;
    }}

    // check if no MB loss in this function
    //double tot1 = MapTotal(*TileWaterVol);
    //qDebug() << tot << tot1 << totq << tot-tot1-totq << full << MB;

}

//---------------------------------------------------------------------------

