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
        double HouseWidthDX_ = HouseCover->Drc*_dx;// can be zero
        double RoadWidthHSDX_ = RoadWidthHSDX->Drc;

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
        RoadWidthHSDX_ = std::min(dxa, RoadWidthHSDX_);
        // decrease roadwidth if roads + houses > dx-channel
        if (RoadWidthHSDX_ + HouseWidthDX_ > dxa)
            RoadWidthHSDX_ = std::max(0.0, dxa-HouseWidthDX_);

        SoilWidthDX->Drc = std::max(0.0, dxa - RoadWidthHSDX->Drc - HouseWidthDX_);
        // soilwidth is used in infil, evap and erosion, NOT flow

        HouseCover->Drc = HouseWidthDX_/_dx;
        // recalc the cover if because it may have been adjusted
        //houses are impermeable in ksateff so do have to be done here, with high mannings n, but allow flow
        RoadWidthHSDX->Drc = RoadWidthHSDX_;

        // adjust man N
        N->Drc = N->Drc + HouseCover->Drc; // increase N for a house, very high resistance
        N->Drc = N->Drc * (1-RoadWidthHSDX->Drc/_dx) + 0.016 * (RoadWidthHSDX->Drc/_dx); // asphalt manning's n
        //https://www.engineeringtoolbox.com/mannings-roughness-d_799.html

        // adjust surface storage for pixels with roads
        RR->Drc = RR->Drc*(1-RoadWidthHSDX->Drc/_dx) + 0.2 * (RoadWidthHSDX->Drc/_dx); // assume smooth asphalt surface
        double RRmm = 10 * RR->Drc;
        MDS->Drc = std::max(0.0, 0.243*RRmm + 0.010*RRmm*RRmm - 0.012*RRmm*tan(asin(Grad->Drc))*100);
        MDS->Drc /= 1000; // convert to m

        FlowWidth->Drc = ChannelAdj->Drc * rillfactor;
        // water can flow everywhere, a house is permeable and a migh mannings n, roads are smooth
        // if hosues are part of the dem than the water automatically flows around it
    }}

    // starting with a water level
    if (SwitchFloodInitial) {
        WHinitVolTot = 0;
        FOR_ROW_COL_MV_L {
            if(SwitchKinematic2D == K2D_METHOD_DYN)
                WH->Drc = hmxInit->Drc;
            else {
                if(SwitchKinematic2D == K2D_METHOD_KINDYN)
                    hmx->Drc = hmxInit->Drc;
            }

            WHinitVolTot += hmxInit->Drc * CHAdjDX->Drc;
        }}
    }
}
//---------------------------------------------------------------------------
void TWorld::cell_SurfaceStorage(int r, int c)
{
    if (FloodDomain->Drc == 0) {
        double wh =  WH->Drc;
        double mds = std::max(0.0, MDS->Drc*(1-exp(-1.875*wh/(0.01*RR->Drc))));
        // surface storage on rough surfaces
        // non-linear release fo water from depression storage
        // resembles curves from GIS surface tests, unpublished

        //   additional Fayna Yuu type storage in m3 per cell
        //   if in a channel cell the store is taken from the channel flow (buffer)
        if (SwitchGridRetention && ChannelWidth->Drc == 0) {
            double dvol = std::max(0.0,GridRetention->Drc - GridRetentionAct->Drc);
            if(dvol > 0) {
                double dh = dvol/CHAdjDX->Drc;
                dh = std::min(wh, dh);
                wh -= dh;
                dvol = dh*CHAdjDX->Drc;
                GridRetentionAct->Drc += dvol;
            }
        }

        WHrunoff->Drc = std::max(0.0, wh-mds);
        // used to be ((wh - WHs)*SW + WHr*RW)/(SW+RW);
        // WH of overlandflow above surface storage

        WHstore->Drc = std::min(mds, wh);
        WH->Drc = WHrunoff->Drc + WHstore->Drc;
        MicroStoreVol->Drc = CHAdjDX->Drc*WHstore->Drc; //RR is adjusted for roads so over entire flowwidth
        // microstore vol in m3
        hmxWH->Drc = WH->Drc;
        WaterVolall->Drc = hmxWH->Drc * CHAdjDX->Drc;
        // non moving microstorage
    } else {
        double wh =  hmx->Drc;
        double mds = std::max(0.0, MDS->Drc*(1-exp(-1.875*wh/(0.01*RR->Drc))));
        // surface storage on rough surfaces
        // non-linear release fo water from depression storage
        // resembles curves from GIS surface tests, unpublished

        //   additional Fayna Yuu type storage in m3 per cell
        //   if in a channel cell the store is taken from the channel flow (buffer)
        if (SwitchGridRetention && ChannelWidth->Drc == 0) {
            double dvol = std::max(0.0,GridRetention->Drc - GridRetentionAct->Drc);
            if(dvol > 0) {
                double dh = dvol/CHAdjDX->Drc;
                dh = std::min(wh, dh);
                wh -= dh;
                dvol = dh*CHAdjDX->Drc;
                GridRetentionAct->Drc += dvol;
            }
        }

        hmxrunoff->Drc = std::max(0.0, wh-mds);
        // used to be ((wh - WHs)*SW + WHr*RW)/(SW+RW);
        // WH of overlandflow above surface storage

        WHstore->Drc = std::min(mds, wh);
        // non moving microstorage
        hmx->Drc = hmxrunoff->Drc + WHstore->Drc;
        MicroStoreVol->Drc = CHAdjDX->Drc*WHstore->Drc; //RR is adjusted for roads so over entire flowwidth
        // microstore vol in m3
        hmxWH->Drc = hmx->Drc;
        WaterVolall->Drc = hmxWH->Drc * CHAdjDX->Drc;
    }
}
//---------------------------------------------------------------------------
