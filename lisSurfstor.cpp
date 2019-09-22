
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
void TWorld::GridCell(int thread)
{

    FOR_ROW_COL_2DMT
    {
        double dxa = _dx;
        if(SwitchIncludeChannel)
        {
            dxa = std::max(0.05*_dx, _dx - ChannelWidthExtended->Drc);
            //VJ 161222 was 0.1, now 0.05*_dx
        }

        ChannelAdj->Drc = dxa;

        RoadWidthDX->Drc = std::min(dxa, RoadWidthDX->Drc);
        dxa = std::max(0.0, dxa - RoadWidthDX->Drc);

        HouseWidthDX->Drc = std::min(dxa*0.95, HouseWidthDX->Drc);
        HouseCover->Drc = HouseWidthDX->Drc/_dx;
        if (Cover->Drc + HouseCover->Drc > 1.0)
            Cover->Drc = 1.0-HouseCover->Drc;

        SoilWidthDX->Drc = dxa;  //ezxcluding roads, including houses,
        //houses are assumed to be permeable but with high mannings n
    }}}}
}
//---------------------------------------------------------------------------
/// Adds new rainfall afterinterception to runoff water nheight or flood waterheight
void TWorld::addRainfallWH(int thread)
{
    FOR_ROW_COL_2DMT
    {
        if (FloodDomain->Drc > 0)
        {
            hmx->Drc += RainNet->Drc + Snowmeltc->Drc;
        }else
        {
            WH->Drc += RainNet->Drc + Snowmeltc->Drc;
            // add net to water rainfall on soil surface (in m)

            if (SwitchGrassStrip && GrassWidthDX->Drc > 0)
                WHGrass->Drc += RainNet->Drc + Snowmeltc->Drc;
            // net rainfall on grass strips, infil is calculated separately for grassstrips

            if (SwitchRoadsystem && RoadWidthDX->Drc > 0)
                WHroad->Drc += Rainc->Drc + Snowmeltc->Drc;
            // assume no interception and infiltration on roads, gross rainfall
        }
    }}}}
}
//---------------------------------------------------------------------------
void TWorld::SurfaceStorage(int thread)
{
    FOR_ROW_COL_2DMT
    {
        double RRm = 0.01*RR->Drc; // assume RR in cm convert to m
        double wh = WH->Drc, whflow = 0;
       // double SDS;
        double mds = MDS->Drc;  // mds is in meters
        double WaterVolrunoff;

        //### surface storage on rough surfaces
        //SDS = 0.1*mds;
        // arbitrary minimum depression storage is 10% of max depr storage, in m

        if (SwitchHardsurface)
            RRm = RRm * (1-HardSurface->Drc) + 0.001*HardSurface->Drc;

        if (RRm < 0.001)
            whflow = wh;
        else
            whflow = std::max(0.0, wh - mds*(1-exp(-1.875*(wh/RRm))) );
//        if (mds > 0)
//            whflow = std::max(0.0,wh-SDS) * (1-exp(-wh/mds));
        // non-linear release fo water from depression storage
        // resemles curves from GIS surface tests, unpublished


        // whflow = std::max(0.0, wh-SDS);
        // subtract surface storage and calc water available for runoff, in m
        // assumed on soilsurface because there is the roughness

        WHstore->Drc = wh - whflow;
        // average water stored on flowwidth and not available for flow, in m


        // soilwidth already includes houses
        //houses no surf storage
        if (SwitchHouses)
        {
            WHstore->Drc *= 1-HouseCover->Drc;
            whflow = wh - WHstore->Drc;
        }
        // hard surface no surface storage
//        if (SwitchHardsurface)
//        {
//            WHstore->Drc *= 1-HardSurface->Drc;
//            whflow = wh - WHstore->Drc;
//        }

        WaterVolrunoff = DX->Drc*( whflow*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
        // runoff volume available for flow, surface + road
        // soil surface excludes houses and roads and channels
        WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
        // all water in the cell incl storage

        //### the trick: use ponded area for flowwidth
        // for "natural" soil surfaces with a given roughness
//        if (RRm == 0)
//            fpa->Drc = 1;
//        else
//            fpa->Drc = 1-exp(-1.875*(wh/RRm));
        // fraction ponded area of a gridcell
//        if(SwitchKinematic2D == 3)
        fpa->Drc = 1;

//        if (SwitchHardsurface)
//            fpa->Drc = fpa->Drc *(1-HardSurface->Drc) + HardSurface->Drc;
//        // hard surface has no roughness so fpa is 1 there

        FlowWidth->Drc = fpa->Drc*SoilWidthDX->Drc + RoadWidthDX->Drc;
        // calculate flowwidth by fpa*surface + road, excludes channel already

//        if (GrassFraction->Drc > 0)
//            FlowWidth->Drc = GrassWidthDX->Drc + (1-GrassFraction->Drc)*FlowWidth->Drc;
        // assume grassstrip spreads water over entire width

        //Houses
//        if (SwitchHouses)
//            FlowWidth->Drc = (1-0.5*HouseCover->Drc)*FlowWidth->Drc;
        // assume house severely restricts flow width, 0.5 is arbitrary
        // cannot be zero flowwidth in 100% house pixel because watwer would not go anywhere

        FlowWidth->Drc = std::min(ChannelAdj->Drc, FlowWidth->Drc);

        if (FlowWidth->Drc > 0)
            WHrunoff->Drc = WaterVolrunoff/(DX->Drc*FlowWidth->Drc);
        else
            WHrunoff->Drc = 0;
        // average WHrunoff from soil surface + roads, because kin wave can only do one discharge
        // this now takes care of ponded area, so water height is adjusted
    }}}}


}
//---------------------------------------------------------------------------

