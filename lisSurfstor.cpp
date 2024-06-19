
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
        SoilWidthDX->Drc = std::max(0.0, dxa - RoadWidthHSDX->Drc - HouseWidthDX_);
        // soilwidth is used in infil, evap and erosion

        HouseCover->Drc = HouseWidthDX_/_dx;
        //houses are impermeable in ksateff so do have to be done here, with high mannings n, but allow flow

        // adjust man N
        N->Drc = N->Drc + 1.0*HouseCover->Drc; // N is 1 for a house, very high resistance
        N->Drc = N->Drc * (1-RoadWidthHSDX->Drc/_dx) + 0.016 * (RoadWidthHSDX->Drc/_dx); // asphalt manning's n
        //https://www.engineeringtoolbox.com/mannings-roughness-d_799.html

        // adjust surface storage for pixels with roads
        RR->Drc = RR->Drc*(1-RoadWidthHSDX->Drc/_dx) + 0.2 * (RoadWidthHSDX->Drc/_dx); // assume smooth asphalt surface
        double RRmm = 10 * RR->Drc;
        MDS->Drc = std::max(0.0, 0.243*RRmm + 0.010*RRmm*RRmm - 0.012*RRmm*tan(asin(Grad->Drc))*100);
        MDS->Drc /= 1000; // convert to m

        FlowWidth->Drc = ChannelAdj->Drc;
        // water can flow everywhere, a house is permeable and a migh mannings n, roads are smooth
        // if hosues are part of the dem than the water automatically flows around it
    }}

    // starting with a water level
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
}
//---------------------------------------------------------------------------
/// Adds new rainfall after interception to runoff water height or flood waterheight
// OBSOLETE not used
void TWorld::addRainfallWH()
{    
    if (SwitchKinematic2D != K2D_METHOD_KINDYN) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            WH->Drc += RainNet->Drc;// + Snowmeltc->Drc;
            // add net to water rainfall on soil surface (in m)
        }}
        // Switch floodinitial is false if not 2D flow
        if (SwitchFloodInitial) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (hmxInit->Drc > 0) {
                    hmxInit->Drc += RainNet->Drc;// + Snowmeltc->Drc;
                    WH->Drc = hmxInit->Drc;
                }
            }}
        }
    }

    if (SwitchKinematic2D == K2D_METHOD_KINDYN) {
        // TODO: hmx is the flooded part when we have kin wave + flooding, else this is not used, floodDomain = 0
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (FloodDomain->Drc > 0) {
                hmx->Drc += RainNet->Drc;// + Snowmeltc->Drc;
                if (SwitchFloodInitial && hmxInit-> Drc > 0)
                    hmx->Drc = hmxInit->Drc;
            }
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
    double WHs = std::max(0.0, std::min(wh, MDS->Drc*(1-exp(-1.875*wh/(0.01*RR->Drc)))));
    // surface storage on rough surfaces
    // non-linear release fo water from depression storage
    // resemles curves from GIS surface tests, unpublished
    // note: roads and houses are assumed to be smooth!

    WHrunoff->Drc = wh-WHs;// ((wh - WHs)*SW + WHr*RW)/(SW+RW);
    // WH of overlandflow above surface storage

    WHstore->Drc = WHs;
    // non moving microstorage
    MicroStoreVol->Drc = DX->Drc*WHstore->Drc*FlowWidth->Drc; //RR is adjusted for roads so over entire flowwidth
    // microstore vol in m3

    WaterVolall->Drc = WHrunoff->Drc*CHAdjDX->Drc + MicroStoreVol->Drc;
    // all water in the cell incl storage
}
//---------------------------------------------------------------------------
