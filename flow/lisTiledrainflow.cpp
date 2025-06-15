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

//#include <algorithm>
#include "model.h"
//#include "operation.h"

/*
//---------------------------------------------------------------------------
// flow in all road cells to tiledrain
//fraction of water and sediment flowing from the surface to the tiledrain system
void TWorld::ToTiledrain()
{
    if (SwitchIncludeStormDrains)  //SwitchIncludeTile ||
    {
        if (SwitchIncludeStormDrains)  {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_TILEL {
                RunoffVolinToTile->Drc = 0;
                double MaxVol = DX->Drc*TileArea->Drc;

                if (TileWaterVol->Drc < MaxVol) {

                    RunoffVolinToTile->Drc = WaterVolall->Drc * TileDrainSize/CHAdjDX->Drc * (DX->Drc/TileDrainDistance);//* RoadWidthHSDX->Drc/_dx
                 //   double volin = _dt*std::sqrt(2*GRAV*WHrunoff->Drc)*TileDrainSize * (DX->Drc/TileDrainDistance);// Bernouilly flow through a hole
                 //   RunoffVolinToTile->Drc = std::min(volin, RunoffVolinToTile->Drc);
                    RunoffVolinToTile->Drc = std::min(WaterVolall->Drc, RunoffVolinToTile->Drc);
                    RunoffVolinToTile->Drc = std::min(MaxVol - TileWaterVol->Drc, RunoffVolinToTile->Drc);
                    WaterVolall->Drc -= RunoffVolinToTile->Drc;

                    WH->Drc  = WaterVolall->Drc/CHAdjDX->Drc;
                    WHrunoff->Drc = std::max(0.0, WH->Drc - WHstore->Drc);
                    hmxWH->Drc = hmx->Drc + WH->Drc;
                }
            }}
        }
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
  FOR_ROW_COL_MV_TILEL {
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

#pragma omp parallel for num_threads(userCores)
FOR_ROW_COL_MV_TILEL {
    if (LDDTile->Drc == 5) {
        qDebug() << TileWaterVol->Drc << TileArea->Drc * DX->Drc << MB;
        TileWaterVol->Drc = std::min(TileWaterVol->Drc, TileArea->Drc * DX->Drc);
    }
}}

}
//---------------------------------------------------------------------------





*/

//---------------------------------------------------------------------------
// flow in all road cells to tiledrain, according to frraction road and fraction openings in surface
// fraction of water flowing from the surface to the tiledrain system
// no sediment
void TWorld::ToTiledrain()
{
    if (SwitchIncludeStormDrains)  {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_TILEL {
            RunoffVolinToTile->Drc = 0;
            double MaxVol = DX->Drc*TileArea->Drc;

            if (TileWaterVol->Drc < MaxVol && WaterVolall->Drc > MicroStoreVol->Drc) {

                //RunoffVolinToTile->Drc = WaterVolall->Drc * TileDrainSize/CHAdjDX->Drc * (DX->Drc/TileDrainDistance);//* RoadWidthHSDX->Drc/_dx

                double volin = _dt*std::sqrt(2*GRAV*WHrunoff->Drc)*TileDrainSize * (DX->Drc/TileDrainDistance)* RoadWidthHSDX->Drc/_dx;
                // Bernouilly flow s*m/s*m2=m3 through a hole

               // qDebug() << volin <<  WaterVolall->Drc << RunoffVolinToTile->Drc;
                RunoffVolinToTile->Drc = std::min(volin, WaterVolall->Drc-MicroStoreVol->Drc);//RunoffVolinToTile->Drc);
                RunoffVolinToTile->Drc = std::min(MaxVol - TileWaterVol->Drc, RunoffVolinToTile->Drc);
                WaterVolall->Drc -= RunoffVolinToTile->Drc;

                WH->Drc = WaterVolall->Drc/CHAdjDX->Drc;
                WHrunoff->Drc = std::max(0.0, WH->Drc-WHstore->Drc);
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
        if (Area > 0.999*TileArea->Drc) {
            TileQ->Drc = TileMaxQ->Drc;
            TileAlpha->Drc = TileMaxAlpha->Drc;
        } else {
            Perim = TileWidth->Drc + Area/TileWidth->Drc; //(=w+2*h)
            TileMaxQ->Drc = Area*pow(Area/Perim,2.0/3.0) * sqrt(TileGrad->Drc)/TileN->Drc;;
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
      double Area = TileWaterVol->Drc / DX->Drc;
      if (Area > 0.999*TileArea->Drc) {
          TileQ->Drc = TileMaxQ->Drc;
          TileAlpha->Drc = TileMaxAlpha->Drc;
      } else {
          double a = Area/TileArea->Drc;
          double theta = pipeThetafroma(r,c,a);
          double Perim = TileDiameter->Drc/2.0*theta;

          if (Perim < 1e-6)
              TileQ->Drc = 0;
          else
              TileQ->Drc = std::pow(Area/Perim, 5.0/3.0) * sqrt(TileGrad->Drc)/TileN->Drc;
         // TileQ->Drc = std::min(TileQ->Drc, TileMaxQ->Drc);

          TileAlpha->Drc = std::pow(std::pow(Perim, 2.0/3.0)*TileN->Drc/sqrt(TileGrad->Drc), 0.6);
          //TileAlpha->Drc  = Area/std::pow(TileQ->Drc, BETAcirc);
      }
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
        TileQn->Drc = 0;//TileQ->Drc;
    }}
    int full = 0;
  // #pragma omp parallel for ordered num_threads(userCores)
  // parallel doesn't work here because the order of cells has to be maintained
    //return;
    double tot = MapTotal(*TileWaterVol);
    double totq = 0;

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

        if (TileWaterVol->Drc > 0.999*TileArea->Drc*DX->Drc) { //+_dt*(Qin-TileQ->Drc)
            TileQn->Drc = 0;// TileMaxQ->Drc;
            TileWaterVol->Drc = TileArea->Drc*DX->Drc;
            full += 1;
        } else {
            TileQn->Drc = 0;//TileQ->Drc*0.1;
            //TileQn->Drc = IterateToQnew(Qin, TileQ->Drc, TileAlpha->Drc, _dt, DX->Drc, TileMaxQ->Drc, TileMaxAlpha->Drc);
            //0.5*(TileQ->Drc+Qin);//
            // TileQn->Drc = std::min(Qin+TileWaterVol->Drc/_dt, TileQn->Drc);
            TileQn->Drc = std::min(TileQn->Drc, TileMaxQ->Drc);

            TileWaterVol->Drc = TileWaterVol->Drc + _dt*(Qin - TileQn->Drc);
            TileWaterVol->Drc = std::max(0.0, TileWaterVol->Drc);
            TileWaterVol->Drc = std::min(TileWaterVol->Drc, TileArea->Drc*DX->Drc);

        }

    }

    double tot1 = MapTotal(*TileWaterVol);
    if (tot-tot1-totq > 1e-8) {

        double dtot = fabs(tot1) > 0 ? (tot - tot1-totq)/tot1 : 0;
        if (dtot > 0) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_TILEL {
                TileWaterVol->Drc = TileWaterVol->Drc*(1.0 + dtot);            // <- distribution weighted to h
                TileWaterVol->Drc = std::max(TileWaterVol->Drc , 0.0);
                //TileWaterVol->Drc = std::min(TileWaterVol->Drc , tma->Drc);
            }}
        }
    }
    tot1 = MapTotal(*TileWaterVol);
    qDebug() << tot << tot1 << totq << tot-tot1-totq << full << MB;

}

/*

void TWorld::TileFlow(void)
{
    if (!SwitchIncludeTile && !SwitchIncludeStormDrains)
        return;

    // get water from surface
    if (SwitchIncludeStormDrains) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_TILEL {
            TileWaterVol->Drc += RunoffVolinToTile->Drc;
            TileWaterVol->Drc = std::min(TileArea->Drc*DX->Drc,TileWaterVol->Drc);
            // add water from the surface
            TileQn->Drc = 0;

        }}
    }

    // get water from soil
    if (SwitchIncludeTile) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_TILEL {
            TileWaterVol->Drc += TileWaterVolSoil->Drc;
            TileWaterVol->Drc = std::min(TileArea->Drc*DX->Drc,TileWaterVol->Drc);
            TileQn->Drc = 0;
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

    // for mass balance correction
    double tot = MapTotal(*TileWaterVol);
    double totq = 0;
    int full = 0 ;

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

        if (TileWaterVol->Drc+_dt*(Qin-TileQ->Drc) > 0.999*TileArea->Drc*DX->Drc) {
            TileQn->Drc = TileMaxQ->Drc;
            TileWaterVol->Drc = TileArea->Drc*DX->Drc;
            full += 1;
        } else {
            TileQn->Drc = IterateToQnew(Qin, TileQ->Drc, TileAlpha->Drc, _dt, DX->Drc, TileMaxQ->Drc, TileMaxAlpha->Drc);
            TileQn->Drc = std::min(Qin+TileWaterVol->Drc/_dt, TileQn->Drc);
            TileQn->Drc = std::min(TileQn->Drc, TileMaxQ->Drc);

            TileWaterVol->Drc = TileWaterVol->Drc + _dt*(Qin - TileQn->Drc);
            TileWaterVol->Drc = std::max(0.0, TileWaterVol->Drc);
            TileWaterVol->Drc = std::min(TileWaterVol->Drc, TileArea->Drc*DX->Drc);
        }

        if (crlinkedlddtile_.at(i_).ldd == 5)
            totq += TileQn->Drc*_dt;
    }

double tot1 = MapTotal(*TileWaterVol);
//qDebug() << tot << tot1 << totq << tot-tot1-totq << full;
    // if (tot-tot1-totq > 1e-8) {

    //     double dtot = fabs(tot1) > 0 ? (tot - tot1-totq)/tot1 : 0;
    //     if (dtot > 0) {
    //         #pragma omp parallel for num_threads(userCores)
    //         FOR_ROW_COL_MV_TILEL {
    //             TileWaterVol->Drc = TileWaterVol->Drc*(1.0 + dtot);            // <- distribution weighted to h
    //             TileWaterVol->Drc = std::max(TileWaterVol->Drc , 0.0);
    //             //TileWaterVol->Drc = std::min(TileWaterVol->Drc , tma->Drc);
    //         }}
    //     }
    // }

//tot1 = MapTotal(*TileWaterVol);
qDebug() << tot << tot1 << totq << tot-tot1-totq << full;
//report(*TileWaterVol,"TV");

}
//---------------------------------------------------------------------------
*/
