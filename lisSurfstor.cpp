
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
        if(SwitchIncludeChannel && ChannelWidthMax->Drc > 0)
            dxa = std::max(0.05, _dx - ChannelWidthMax->Drc);
//        dxa = std::max(0.05, _dx - ChannelWidthExtended->Drc);

        if (SwitchCulverts && ChannelMaxQ->Drc > 0)
            dxa = _dx;

        ChannelAdj->Drc = dxa;
        CHAdjDX->Drc = dxa*DX->Drc;

        RoadWidthHSDX->Drc = std::min(dxa, RoadWidthHSDX->Drc);
        dxa = std::max(0.0, dxa - RoadWidthHSDX->Drc);

        HouseWidthDX->Drc = std::min(dxa*0.95, HouseWidthDX->Drc);
        HouseCover->Drc = HouseWidthDX->Drc/_dx;
        if (Cover->Drc + HouseCover->Drc > 1.0)
            Cover->Drc = 1.0-HouseCover->Drc;

        SoilWidthDX->Drc = dxa;  //excluding roads, including houses, hard surface
        //houses are assumed to be permeable but with high mannings n
    }}
}
//---------------------------------------------------------------------------
/// Adds new rainfall afterinterception to runoff water nheight or flood waterheight
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
void TWorld::SurfaceStorage()
{
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double wh = WH->Drc;
        //double WaterVolrunoff;
        double SW = SoilWidthDX->Drc;
        double RW = RoadWidthHSDX->Drc;
        double WHr = WHroad->Drc;
        double WHs;

        //### surface storage on rough surfaces
        WHs = std::min(wh, MDS->Drc*(1-exp(-1.875*(wh/std::max(0.01,0.01*RR->Drc)))));
        // non-linear release fo water from depression storage
        // resemles curves from GIS surface tests, unpublished

        double FW = std::min(ChannelAdj->Drc, SW + RW);
        // calculate flowwidth by fpa*surface + road, excludes channel already

    //    if (FW > 0) {
          //  WaterVolrunoff = DX->Drc*( whflow*SoilWidthDX->Drc + WHroad->Drc*RoadWidthHSDX->Drc);
            // runoff volume available for flow, surface + road + hard surfaces
            // soil surface excludes houses and roads and channels and hard surfaces
           // WHrunoff->Drc = WaterVolrunoff/(DX->Drc*FW);
        WHrunoff->Drc = ((wh - WHs)*SW + WHr*RW)/FW;
      //  } else
        //    WHrunoff->Drc = 0;
        // average WHrunoff from soil surface + roads, because kin wave can only do one discharge
        // this now takes care of ponded area, so water height is adjusted
        FlowWidth->Drc = FW;

        WaterVolall->Drc = DX->Drc*(wh*SW + WHr*RW);
        // all water in the cell incl storage
        WHstore->Drc = WHs;
    }}
}
//---------------------------------------------------------------------------
void TWorld::doETa()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double tot = 0;
        double tmp = 0;
        double ETp = 5.0/(12.0*3600.0) * 0.001 *_dt;  //5m/day in m per dt
        ETp = Rain->Drc > 0 ? 0.0 : ETp;

        double ETa_soil = hmxWH->Drc > 0 ? 0.0: ThetaI1->Drc/ThetaS1->Drc * ETp;
        ETa_soil *= Cover->Drc;
        double moist = ThetaI1->Drc * SoilDepth1->Drc;
        tmp = moist;
        moist = std::max(0.0, moist - ETa_soil);
        tmp = tmp - moist;
        ThetaI1->Drc = moist/SoilDepth1->Drc;
        tot = tot + tmp;

        double ETa_int = Cover->Drc * ETp;
        tmp = CStor->Drc;
        CStor->Drc = std::max(0.0, CStor->Drc-ETa_int);
        RainCum->Drc = std::max(0.0, RainCum->Drc-ETa_int);
        tmp = tmp - CStor->Drc;
        Interc->Drc =  Cover->Drc * CStor->Drc * SoilWidthDX->Drc * DX->Drc;
        tot = tot + tmp;

        double ETa_pond = hmxWH->Drc > 0 ? ETp : 0.0;
        if (FloodDomain->Drc > 0) {
            tmp = hmx->Drc;
            hmx->Drc = std::max(0.0, hmx->Drc-ETa_pond );
            tmp = tmp - hmx->Drc;
//            FloodWaterVol->Drc = hmx->Drc*CHAdjDX->Drc;
//            hmxWH->Drc = hmx->Drc;   //hmxWH is all water
//            hmxflood->Drc = hmxWH->Drc < minReportFloodHeight ? 0.0 : hmxWH->Drc;
        }
        else {
            tmp = WHrunoff->Drc;
            WHrunoff->Drc = std::max(0.0, WHrunoff->Drc-ETa_pond );
            tmp = tmp - WHrunoff->Drc;
            WaterVolall->Drc = WHrunoff->Drc*CHAdjDX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
            WHroad->Drc = WHrunoff->Drc;
            //WHGrass->Drc = WHrunoff->Drc;
            WH->Drc = WHrunoff->Drc + WHstore->Drc;
//            hmxWH->Drc = WH->Drc;
//            hmx->Drc = std::max(0.0, WH->Drc - minReportFloodHeight);
//            hmxflood->Drc = hmxWH->Drc < minReportFloodHeight ? 0.0 : hmxWH->Drc;
        }
        tot = tot + tmp;

        ETa->Drc += tot;

//TODO fpa


    }}
}

