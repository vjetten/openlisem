
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


//---------------------------------------------------------------------------
void TWorld::GridCell()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double dxa = _dx;
        double HouseWidthDX_ = HouseCover->Drc*_dx;

        if(SwitchIncludeChannel) {
            if (ChannelWidth->Drc > 0){
                dxa = _dx - ChannelWidth->Drc;
                if (SwitchCulverts && ChannelMaxQ->Drc > 0)
                    dxa = _dx;
            }
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
        SoilWidthDX->Drc = std::max(0.0, dxa - RoadWidthHSDX->Drc);
        // including houses in soilwidth gives large MB errors! WHY!!!

        HouseCover->Drc = HouseWidthDX_/_dx;
        //houses are impermeable in ksateff so do have to be done here, with high mannings n, but allow flow

        // adjust man N
        N->Drc = N->Drc + 1.0*HouseCover->Drc; // N is 1 for a house, very high resistance
        N->Drc = N->Drc * (1-RoadWidthHSDX->Drc/_dx) * 0.016 * (RoadWidthHSDX->Drc/_dx);
        //https://www.engineeringtoolbox.com/mannings-roughness-d_799.html

        FlowWidth->Drc = ChannelAdj->Drc;//is the same as SoilWidthDX->Drc + RoadWidthHSDX->Drc;
    }}

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
        if (SwitchFloodInitial  && hmxInit->Drc > 0)
            hmxInit->Drc += RainNet->Drc + Snowmeltc->Drc;

        if (FloodDomain->Drc > 0) {
            hmx->Drc += RainNet->Drc;// + Snowmeltc->Drc;
            if (SwitchFloodInitial && hmxInit-> Drc > 0)
                hmx->Drc = hmxInit->Drc;
        } else {
            WH->Drc += RainNet->Drc;// + Snowmeltc->Drc;
            // add net to water rainfall on soil surface (in m)
            if (SwitchFloodInitial && hmxInit-> Drc > 0)
                WH->Drc = hmxInit->Drc;
        }
    }}

    // if (SwitchRoadsystem || SwitchHardsurface) {  //???? separate hs from road here?
    //     #pragma omp parallel for num_threads(userCores)
    //     FOR_ROW_COL_MV_L {
    //         if (RoadWidthHSDX->Drc > 0)
    //             WHroad->Drc += Rainc->Drc;// + Snowmeltc->Drc;
    //     }}
    // }

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
    double WHs = std::max(0.0, std::min(wh, MDS->Drc*(1-exp(-1.875*wh/(0.01*RR->Drc)))));
    //surface storage on rough surfaces
    // non-linear release fo water from depression storage
    // resemles curves from GIS surface tests, unpublished
    // note: roads and houses are assumed to be smooth!

    WHrunoff->Drc = wh-WHs;// ((wh - WHs)*SW + WHr*RW)/(SW+RW);
    // WH of overlandflow above surface storage

    WHstore->Drc = WHs;
    // non moving microstorage
    MicroStoreVol->Drc = DX->Drc*WHstore->Drc*(FlowWidth->Drc-RoadWidthHSDX->Drc); //-RoadWidthHSDX->Drc
    // microstore vol in m3

    WaterVolall->Drc = WHrunoff->Drc*CHAdjDX->Drc + MicroStoreVol->Drc;
    // all water in the cell incl storage
}
//---------------------------------------------------------------------------
