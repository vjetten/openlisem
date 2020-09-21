
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

        ChannelAdj->Drc = dxa;

        RoadWidthHSDX->Drc = std::min(dxa, RoadWidthHSDX->Drc);
        dxa = std::max(0.0, dxa - RoadWidthHSDX->Drc);

        HouseWidthDX->Drc = std::min(dxa*0.95, HouseWidthDX->Drc);
        HouseCover->Drc = HouseWidthDX->Drc/_dx;
        if (Cover->Drc + HouseCover->Drc > 1.0)
            Cover->Drc = 1.0-HouseCover->Drc;

        SoilWidthDX->Drc = dxa;  //excluding roads, including houses, hard surface
        //houses are assumed to be permeable but with high mannings n
    }}
report(*SoilWidthDX, "sw.map");
}
//---------------------------------------------------------------------------
/// Adds new rainfall afterinterception to runoff water nheight or flood waterheight
void TWorld::addRainfallWH()
{
    if(SwitchKinematic2D == K2D_METHOD_DYN) {
#pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            WH->Drc += RainNet->Drc + Snowmeltc->Drc;
            // add net to water rainfall on soil surface (in m)

            if (SwitchGrassStrip && GrassWidthDX->Drc > 0)
                WHGrass->Drc += RainNet->Drc + Snowmeltc->Drc;
            // net rainfall on grass strips, infil is calculated separately for grassstrips

            if (SwitchRoadsystem && RoadWidthHSDX->Drc > 0)
                WHroad->Drc += Rainc->Drc + Snowmeltc->Drc;
            // assume no interception and infiltration on roads, gross rainfall
        }}
    } else {
    #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (FloodDomain->Drc > 0) {
                hmx->Drc += RainNet->Drc + Snowmeltc->Drc;
            } else {
                WH->Drc += RainNet->Drc + Snowmeltc->Drc;
                // add net to water rainfall on soil surface (in m)

                if (SwitchGrassStrip && GrassWidthDX->Drc > 0)
                    WHGrass->Drc += RainNet->Drc + Snowmeltc->Drc;
                // net rainfall on grass strips, infil is calculated separately for grassstrips

                if (SwitchRoadsystem && RoadWidthHSDX->Drc > 0)
                    WHroad->Drc += Rainc->Drc + Snowmeltc->Drc;
                // assume no interception and infiltration on roads, gross rainfall
            }
        }}
    }
}
//---------------------------------------------------------------------------
void TWorld::SurfaceStorage()
{
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double RRm = 0.01*RR->Drc; // assume RR in cm convert to m
        double wh = WH->Drc, whflow = 0;
        double mds = MDS->Drc;  // mds is in meters
        double WaterVolrunoff;

        //### surface storage on rough surfaces

        if (RRm < 0.001)
            whflow = wh;
        else
            whflow = std::max(0.0, wh - mds*(1-exp(-1.875*(wh/RRm))) );
        // non-linear release fo water from depression storage
        // resemles curves from GIS surface tests, unpublished

        WHstore->Drc = wh - whflow;
        // average water stored on flowwidth and not available for flow, in m

        WaterVolrunoff = DX->Drc*( whflow*SoilWidthDX->Drc + WHroad->Drc*RoadWidthHSDX->Drc);
        // runoff volume available for flow, surface + road + hard surfaces
        // soil surface excludes houses and roads and channels and hard surfaces

        WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthHSDX->Drc);
        // all water in the cell incl storage

        //### the trick: use ponded area for flowwidth
        // for "natural" soil surfaces with a given roughness
//        if (RRm == 0)
//            fpa->Drc = 1;
//        else
//            fpa->Drc = 1-exp(-1.875*(wh/RRm));
        // fraction ponded area of a gridcell
        double Fpa = 1;

        double FW = std::min(_dx, Fpa*SoilWidthDX->Drc + RoadWidthHSDX->Drc);
        // calculate flowwidth by fpa*surface + road, excludes channel already

        FW = std::min(ChannelAdj->Drc, FW);

        if (FW > 0)
            WHrunoff->Drc = WaterVolrunoff/(DX->Drc*FW);
        else
            WHrunoff->Drc = 0;
        // average WHrunoff from soil surface + roads, because kin wave can only do one discharge
        // this now takes care of ponded area, so water height is adjusted
        FlowWidth->Drc = FW;
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
//            FloodWaterVol->Drc = hmx->Drc*ChannelAdj->Drc*DX->Drc;
//            hmxWH->Drc = hmx->Drc;   //hmxWH is all water
//            hmxflood->Drc = hmxWH->Drc < minReportFloodHeight ? 0.0 : hmxWH->Drc;
        }
        else {
            tmp = WHrunoff->Drc;
            WHrunoff->Drc = std::max(0.0, WHrunoff->Drc-ETa_pond );
            tmp = tmp - WHrunoff->Drc;
            WaterVolall->Drc = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
            WHroad->Drc = WHrunoff->Drc;
            WHGrass->Drc = WHrunoff->Drc;
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
