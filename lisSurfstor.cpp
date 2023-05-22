
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
  \file lisSurfstor.cpp
  \brief calculate surface storage and flow width

functions: \n
- void TWorld::GridCell(void) \n
- void TWorld::addRainfallWH(void) \n
- void TWorld::SurfaceStorage(void)\n
 */


#include <algorithm>
#include "model.h"
//#define tiny 1e-8


//---------------------------------------------------------------------------
void TWorld::GridCell()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double dxa = _dx;
    double HouseWidthDX_ = HouseCover->Drc*_dx;

     if(SwitchIncludeChannel && ChannelWidth->Drc > 0 && ChannelMaxQ->Drc == 0) {
          dxa = _dx - ChannelWidth->Drc;
     }
     //note: channelwidth <= _dx*0.95. ADD channelmaxq here, better MB

        ChannelAdj->Drc = dxa;
        CHAdjDX->Drc = dxa*DX->Drc;

      // adjust houses to cell with channels
      HouseWidthDX_ = std::min(dxa,  HouseWidthDX_);
        // adjust roads+hardsurf to cell with channels
      RoadWidthHSDX->Drc = std::min(dxa, RoadWidthHSDX->Drc);
        // decrease roadwidth if roads + houses > dx-channel
      HouseWidthDX_ = std::min(dxa-RoadWidthHSDX->Drc , HouseWidthDX_);
      // you cannot have houses and a road larger than a pixel
  //    SoilWidthDX->Drc = std::max(0.0,dxa - RoadWidthHSDX->Drc - HouseWidthDX_);
      SoilWidthDX->Drc = dxa - RoadWidthHSDX->Drc;
      // including houses in soilwidth gives large MB errors! WHY!!!

      HouseCover->Drc = HouseWidthDX_/_dx;
      //houses are impermeable in ksateff so do have to be done here, with high mannings n, but allow flow

      N->Drc = N->Drc + 1.0*HouseCover->Drc; // N is 1 for a house, very high resistance
      // adjust man N

      FlowWidth->Drc = ChannelAdj->Drc;//is the same as SoilWidthDX->Drc + RoadWidthHSDX->Drc;
      //FlowWidth->Drc = SoilWidthDX->Drc + RoadWidthHSDX->Drc + HouseWidthDX_;
    }}
//report(*HouseCover,"hc.map");
//report(*SoilWidthDX,"sw.map");
//report(*RoadWidthHSDX,"rw.map");
//report(*FlowWidth,"fw.map");

 /*
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double dxa = _dx; // dxa is dx minus the channel

        if(SwitchIncludeChannel && ChannelWidth->Drc > 0) {
            dxa = _dx - ChannelWidth->Drc; // channelwidth is limited to 0.95*dx in datainit
        }

        ChannelAdj->Drc = dxa; // dx besides the channel
        CHAdjDX->Drc = dxa*DX->Drc; // surface next to the channel

        double HouseWidthDX_ = HouseCover->Drc*_dx;
        HouseWidthDX_ = std::min(0.95*dxa, HouseWidthDX_);
        // adjust houses to the space adjacent to the channel with a little bit of open space, channels have precedence
        HouseCover->Drc = HouseWidthDX_/_dx;
        //houses are impermeable in ksateff so do have to be done here, with high mannings n, but allow flow

        RoadWidthHSDX->Drc = std::min(dxa-HouseWidthDX_, RoadWidthHSDX->Drc);
        // adjust roads+hardsurf to cell with channels and houses.

        SoilWidthDX->Drc = std::max(0.0, dxa-RoadWidthHSDX->Drc-HouseWidthDX_);
        //soil is dx - roads+hardsurf - houses - channels
        //water can infiltrate over soilwidth

        N->Drc = N->Drc * (1-HouseCover->Drc) + 1.0*HouseCover->Drc; // N is 1 for a house, very high resistance
        // adjust man N to houses, rougher
        N->Drc = N->Drc * (1-RoadWidthHSDX->Drc/_dx) + 0.015*(RoadWidthHSDX->Drc/_dx);
        // adjust man N to roads, smoother Mannings n of asphalt is 0.015

        FlowWidth->Drc = ChannelAdj->Drc;//is the same as SoilWidthDX->Drc + RoadWidthHSDX->Drc;
       // FlowWidth->Drc = SoilWidthDX->Drc + RoadWidthHSDX->Drc;

        // water can flow in houses but very high manning's n, or houses are part of the dem and then there is no problem anyway
       // FlowWidth->Drc = SoilWidthDX->Drc + RoadWidthHSDX->Drc;

    }}
*/
    if (SwitchFloodInitial) {
        WHinitVolTot = 0;
        FOR_ROW_COL_MV_L {
            if(SwitchKinematic2D == K2D_METHOD_DYN)
                WH->Drc = hmxInit->Drc;
            else
                hmx->Drc = hmxInit->Drc;

            WHinitVolTot += hmxInit->Drc * CHAdjDX->Drc;
        }}
    }

//    thetai1tot = 0;
//    #pragma omp parallel for reduction(+:thetai1tot) num_threads(userCores)
//    FOR_ROW_COL_MV_L {
//        thetai1tot += ThetaI1->Drc * SoilDepth1->Drc*CHAdjDX->Drc;
//    }}
//    thetai1cur = thetai1tot;

//    if (SwitchTwoLayer) {
//        thetai2tot = 0;
//        #pragma omp parallel for reduction(+:thetai2tot) num_threads(userCores)
//        FOR_ROW_COL_MV_L {
//            thetai2tot += ThetaI2->Drc * (SoilDepth2->Drc-SoilDepth1->Drc)*CHAdjDX->Drc;
//        }}
//        thetai2cur = thetai2tot;
//    }

}
//---------------------------------------------------------------------------
/// Adds new rainfall afterinterception to runoff water nheight or flood waterheight
// OBSOLETE not used
void TWorld::addRainfallWH()
{
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (FloodDomain->Drc > 0) {
                hmx->Drc += RainNet->Drc + Snowmeltc->Drc;
            } else {
                WH->Drc += RainNet->Drc + Snowmeltc->Drc;
                // add net to water rainfall on soil surface (in m)

            //    if (SwitchGrassStrip && GrassWidthDX->Drc > 0)
            //        WHGrass->Drc += RainNet->Drc + Snowmeltc->Drc;
                // net rainfall on grass strips, infil is calculated separately for grassstrips
            }
        }}

    if (SwitchRoadsystem || SwitchHardsurface) {  //???? separate hs from road here?
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (RoadWidthHSDX->Drc > 0)
                WHroad->Drc += Rainc->Drc + Snowmeltc->Drc;
        }}
    }

}
//---------------------------------------------------------------------------
// not used
void TWorld::SurfaceStorage()
{
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        cell_SurfaceStorage(r, c);
    }}
}
//---------------------------------------------------------------------------
void TWorld::cell_SurfaceStorage(int r, int c)
{    
    double wh = WH->Drc;
    double SW = SoilWidthDX->Drc;
    double RW = RoadWidthHSDX->Drc;
    double WHr = WHroad->Drc;
    double WHs = std::min(wh, MDS->Drc*(1-exp(-1.875*wh/(0.01*RR->Drc))));
    //surface storage on rough surfaces
    // non-linear release fo water from depression storage
    // resembles curves from GIS surface tests, unpublished
    // note: roads and houses are assumed to be smooth!

    WHrunoff->Drc = ((wh - WHs)*SW + WHr*RW)/(SW+RW);
    // WH of overlandflow above surface storage

    WHstore->Drc = WHs;
    // non moving microstorage
    MicroStoreVol->Drc = DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
    // microstore vol in m3

    //WaterVolall->Drc = DX->Drc*(wh*SW + WHr*RW);
    WaterVolall->Drc = WHrunoff->Drc*CHAdjDX->Drc + MicroStoreVol->Drc;
    // all water in the cell incl storage
}
//---------------------------------------------------------------------------
