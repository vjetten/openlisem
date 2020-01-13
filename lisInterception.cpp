
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
  \file lisRainintc.cpp
  \brief Get rainfall, make a rainfall map and calculate interception

functions: \n
- void TWorld::Interception(void) \n
- void TWorld::InterceptionHouses(void) \n
- void TWorld::addRainfallWH(void) \n
 */

#include <algorithm>
#include "model.h"

//---------------------------------------------------------------------------
/// Interception()
/// - interception seen as rigid storage SMax filling up and overflowing\n
/// - overflow flux is identical to rainfall flux in intensity\n
/// - SMax is the storage of the plants inside the gridcell, not the average storage of the gridcell\n
/// - so if a single tree inside a cell has an SMax of 2mm even if it covers 10%, the Smax of that cell is 2\n
/// - therefore the same goes for LAI: the LAI of the plants inside the gridcell\n
/// - this is also easier to observe. The LAI from a satellite image is the average LAI of a cell, must be divided by Cover
void TWorld::Interception(int thread)
{
    // all variables are in m
    if (!SwitchRainfall)
        return;
    //VJ 110113 bug fix, no interception when no rainfall and only snowmelt

    FOR_ROW_COL_2DMT
    {
        if (Cover->Drc > 0)// && Rainc->Drc > 0)
        {
            double CS = CStor->Drc;
            //actual canopy storage in m
            double Smax = CanopyStorage->Drc;
            //max canopy storage in m
            double LAIv;
            if (SwitchInterceptionLAI)
                LAIv = LAI->Drc;
            else
                LAIv = (log(1-Cover->Drc)/-0.4)/std::max(0.1,Cover->Drc);
            //Smax is based on LAI and LAI is the average of a gridcell, already including the cover
            // a low cover means a low LAI means little interception

            if (Smax > 0)
            {
                double k = 1-exp(-CanopyOpeness*LAIv);
                //a dense canopy has a low openess factor, so little direct throughfall and high CS

                CS = Smax*(1-exp(-k*RainCum->Drc/Smax));
                //      CS = Smax*(1-exp(-0.0653*LAIv*RainCum->Drc/Smax));
                //VJ 110209 direct use of openess, astons value quite open for eucalypt.
                //A good guess is using the cover LAI relation
                //and interpreting cover as openess: k = exp(-0.45*LAI) BUT this is 0.065!!!
            }
            else
                CS = 0;
            // 0.0653 is canopy openess, based on Aston (1979), based on Merriam (1960/1973), De Jong & Jetten 2003
            // is not the same as canopy cover. it also deals with how easy rainfall drips through the canopy
            // possible to use equation from Ahston but for very open Eucalypt

            //            CS = std::max(0.0, CS * (1-StemflowFraction));
            //VJ 110206 decrease storage with stemflow fraction!
            // but stemflowfraction doesn't go anywhere!!!!!!!!!!!!!!!!!
            // it doesn't matter, either stemflow is removed from the storage and added to RainNet
            // or it is not subtracted in the first place, so the storage is a bit more but leafdrain is earlier at maximum

            LeafDrain->Drc = std::max(0.0, Cover->Drc*(Rainc->Drc - (CS - CStor->Drc)));
            // diff between new and old strage is subtracted from rainfall
            // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!
            // note: cover already implicit in LAI and Smax, part falling on LAI is cover*rainfall

            CStor->Drc = CS;
            // put new storage back in map
            Interc->Drc =  Cover->Drc * CS * SoilWidthDX->Drc * DX->Drc;
            // only on soil surface, not channels or roads, in m3

            RainNet->Drc = LeafDrain->Drc + (1-Cover->Drc)*Rainc->Drc;
            // net rainfall is direct rainfall + drainage
            // rainfall that falls on the soil, used in infiltration
        }
        else
        {
            RainNet->Drc = Rainc->Drc;
        }
    }}}}
}
//---------------------------------------------------------------------------
void TWorld::InterceptionLitter(int thread)
{
    // all variables are in m
    if (!SwitchLitter)
        return;

    if (!SwitchRainfall)
        return;

    FOR_ROW_COL_2DMT
    {
        if (hmx->Drc == 0 && WH->Drc == 0 && Litter->Drc > 0 && RainNet->Drc > 0)
        {

            double Smax = LitterSmax/1000.0;
            // assume simply that the cover linearly scales between 0 and LtterSmax of storage

            double LCS = LCStor->Drc;
            //actual canopy storage in m

    //        if (SwitchHardsurface)
    //            Smax *= (1-HardSurface->Drc);
            //VJ 110111 no litter interception on hard surfaces

            LCS = std::min(LCS + RainNet->Drc, Smax);
            // add water to the storage, not more than max

    //        LCS = std::min(LRainCum->Drc, Smax);
    //        //assume direct simple filling of litter

            double drain = std::max(0.0, Litter->Drc*(RainNet->Drc - (LCS - LCStor->Drc)));
            // diff between new and old storage is subtracted from leafdrip

            LCStor->Drc = LCS;
            // put new storage back in map for next dt

            LInterc->Drc =  Litter->Drc * LCS * SoilWidthDX->Drc * DX->Drc;
            // only on soil surface, not channels or roads, in m3

            RainNet->Drc = drain + (1-Litter->Drc)*RainNet->Drc;// + (1-Cover->Drc)*Rainc->Drc;
            //recalc
        }
    }}}}
}
//---------------------------------------------------------------------------
void TWorld::InterceptionHouses(int thread)
{
    // all variables are in m
    if (!SwitchHouses)
        return;

    if (!SwitchRainfall)
        return;

    FOR_ROW_COL_2DMT
    {
        if (HouseCover->Drc > 0 &&  Rainc->Drc > 0)
        {
            double roofsurface = (SoilWidthDX->Drc * DX->Drc * HouseCover->Drc); // m2
            //house on cell in m2
            double HS = HStor->Drc;
            //actual roof storage in m

            double Hmax = RoofStore->Drc;
            //max roof storage in m

            HS = std::min(HS + RainNet->Drc, Hmax);

            double housedrain = std::max(0.0, HouseCover->Drc * (RainNet->Drc - (HS - HStor->Drc)));
            // overflow in m3/m2 of house

            HStor->Drc = HS;
            // put new storage back in maps in m

            IntercHouse->Drc =  roofsurface * HS; //HouseCover->Drc * HS * SoilWidthDX->Drc * DX->Drc;
            // total interception in m3,exclude roads, channels

            RainNet->Drc = housedrain + (1-HouseCover->Drc)*RainNet->Drc;
            // net rainfall is direct rainfall + drainage
            // rainfall that falls on the soil, used in infiltration

            // filling raindrums with surplus drainage from roofs
            // drum is recalculated to m based on roof surface
            double DS = 0;
            if (SwitchRaindrum && DrumStore->Drc > 0)
            {
                double Dmax = DrumStore->Drc/roofsurface;
                //max drum storage in m as if roof storage is more

                DS = DStor->Drc;
                //actual drum storage in m

                DS = std::min(DS + RainNet->Drc, Dmax);
                // fill tank to max
                double drumdrain = std::max(0.0, HouseCover->Drc * (RainNet->Drc - (DS - DStor->Drc)));

                DStor->Drc = DS;
                // put new drum storage back in maps in m3

                IntercHouse->Drc += roofsurface * DS;//(HouseCover->Drc * DS * SoilWidthDX->Drc * DX->Drc);
                // total interception in m3,exclude roads, channels

                RainNet->Drc = drumdrain + (1-HouseCover->Drc)*RainNet->Drc;
                // net rainfall is net rainfall + drainage
                // rainfall that falls on the soil, used in infiltration

            }
        }
    }}}}
}
//---------------------------------------------------------------------------
