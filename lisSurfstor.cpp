
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
**  website, information and code: http://lisem.sourceforge.net
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
#define tiny 1e-8


//---------------------------------------------------------------------------
void TWorld::GridCell()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double dxa = _dx;

        if(SwitchIncludeChannel && ChannelWidth->Drc > 0) {
            dxa = std::max(0.05, _dx - ChannelWidth->Drc);
            // use adjusted channelwidth here to avoid negative adj
            //        dxa = std::max(0.05, _dx - ChannelWidthExtended->Drc);

            if (SwitchCulverts)
                dxa = _dx;   //als culvert dan geen channelwidth want channel ondergonds
        }

        ChannelAdj->Drc = dxa;
        CHAdjDX->Drc = dxa*DX->Drc;

        //HouseWidthDX->Drc = std::min(dxa*0.95, HouseWidthDX->Drc);

        // adjust roads+hardsurf to cell with channels
        RoadWidthHSDX->Drc = std::min(dxa-HouseWidthDX->Drc, RoadWidthHSDX->Drc);
        SoilWidthDX->Drc = dxa-RoadWidthHSDX->Drc;
        //soil is pixel - roads+hardsurf - channels but including houses


        //HouseWidthDX->Drc = std::min(dxa*0.95, HouseWidthDX->Drc);
        HouseWidthDX->Drc = std::min((dxa-RoadWidthHSDX->Drc)*0.95, HouseWidthDX->Drc);
        // you cannot have houses and a road larger than a pixel
        HouseCover->Drc = HouseWidthDX->Drc/_dx;
        //houses are impermeable in ksateff so do have to be done here, with high mannings n, but allow flow

        N->Drc = N->Drc * (1-HouseCover->Drc) + 0.5*HouseCover->Drc;
        // adjust man N

        FlowWidth->Drc = ChannelAdj->Drc;//is the same as SoilWidthDX->Drc + RoadWidthHSDX->Drc;

    }}

    if (SwitchFloodInitial) {
        WHinitVolTot = 0;
        FOR_ROW_COL_MV_L {
            if(SwitchKinematic2D == K2D_METHOD_DYN)
                WH->Drc = hmxInit->Drc;
            else
                hmx->Drc = hmxInit->Drc;

            WHinitVolTot += hmxInit->Drc * DX->Drc * ChannelAdj->Drc;
        }}
    }
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
    double WHs = std::min(wh, MDS->Drc*(1-exp(-1.875*(wh/std::max(0.005,0.01*RR->Drc)))));
    //surface storage on rough surfaces
    // non-linear release fo water from depression storage
    // resemles curves from GIS surface tests, unpublished

    WHrunoff->Drc = ((wh - WHs)*SW + WHr*RW)/(SW+RW);
    // moving overlandflow above surface storage

    WHstore->Drc = WHs;
    // non moving microstorage

    WaterVolall->Drc = DX->Drc*(wh*SW + WHr*RW);
    // all water in the cell incl storage
}
//---------------------------------------------------------------------------
