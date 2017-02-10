
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
void TWorld::Interception(void)
{
    // all variables are in m
    if (!SwitchRainfall)
        return;
    //VJ 110113 bug fix, no interception when no rainfall and only snowmelt

    FOR_ROW_COL_MV
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
            // avoid division by 0

            if (SwitchBuffers && !SwitchSedtrap)
                if(BufferID->Drc > 0)
                    Smax = 0;
            // no interception with buffers, but sedtrap can have interception

//            if (SwitchHardsurface)
//                Smax *= (1-HardSurface->Drc);
            //VJ 110111 no interception on hard surfaces

//            if (PlantHeight->Drc < WH->Drc)
//            {
//                Smax = 0;
//                CS = 0;
//            }
            //VJ no interception when water level is heigher than plants
            //???? we cannot make interception 0 when water rises because of mass balance


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
//    else
//        {
//            RainNet->Drc = Rainc->Drc;
//        }
}
//---------------------------------------------------------------------------
void TWorld::InterceptionLitter(void)
{
    // all variables are in m

    if (!SwitchLitter)
        return;

    if (!SwitchRainfall)
        return;

    FOR_ROW_COL_MV
            if (hmx->Drc == 0 && WH->Drc == 0 && Litter->Drc > 0 && Rainc->Drc > 0)
    {
        //double LAI = (log(1-std::min(0.9,Litter->Drc))/-0.4);
        // Bracken equation, avoid log 0
        //double Smax = 0.001 * 0.1713 * LAI; // in m
        double Smax = 0.002*Litter->Drc;
        // assume simply that the cover linearly scales between 0 and 2 mm of storage

        double LCS = LCStor->Drc;
        //actual canopy storage in m

        if (SwitchHardsurface)
            Smax *= (1-HardSurface->Drc);
        //VJ 110111 no interception on hard surfaces

        LRainCum->Drc += LeafDrain->Drc;
        // cumulative leaf drainage falling on litter

        LCS = std::min(LRainCum->Drc, Smax);
        //assume direct simple filling of litter

        double drain = std::max(0.0, Litter->Drc*(LeafDrain->Drc - (LCS - LCStor->Drc)));
        // diff between new and old strage is subtracted from leafdrip

        LCStor->Drc = LCS;
        // put new storage back in map
        LInterc->Drc =  Litter->Drc * LCS * SoilWidthDX->Drc * DX->Drc;
        // only on soil surface, not channels or roads, in m3

        RainNet->Drc = drain + (1-Litter->Drc)*LeafDrain->Drc + (1-Cover->Drc)*Rainc->Drc;
        //recalc
    }
}
//---------------------------------------------------------------------------
void TWorld::InterceptionHouses(void)
{
    // all variables are in m
    if (!SwitchHouses)
        return;

    if (!SwitchRainfall)
        return;

    FOR_ROW_COL_MV
    {
        if (HouseCover->Drc > 0 &&  Rainc->Drc > 0)
        {
            //house on cell in m2
            double HS = HStor->Drc;
            //actual roof storage in m
            double DS = DStor->Drc;
            //actual drum storage in m3
            double Hmax = RoofStore->Drc;
            //max roof storage in m

            // GIVES MASS BALANCE ERRORS?
            //            if (Hmax > 0)
            //            {
            //                double k = 1.0;
            //                // speed of filling of roof storage, very quickly
            //                HS = Hmax*(1-exp(-k*RainCum->Drc/Hmax));
            //                //roof storage in m
            //            }
            //            else
            //                HS = 0;

            HS = HS + RainNet->Drc;
            if (HS > Hmax)
                HS = Hmax;

            double housedrain = std::max(0.0, HouseCover->Drc * (RainNet->Drc - (HS - HStor->Drc)));
            // overflow in m3/m2 of house
            HStor->Drc = HS;
            // put new storage back in maps in m

            double Dmax = 0;
            if (SwitchRaindrum)
                Dmax = DrumStore->Drc;
            //max drum storage in m3

            if (Dmax > 0)
            {
                // housedrain water from roof in m, cover is already included
                double dsm3 = (DS + housedrain)*SoilWidthDX->Drc*DX->Drc;
                if (dsm3 < Dmax)
                    dsm3 = Dmax;
                DS = (SoilWidthDX->Drc > 0)? dsm3/(SoilWidthDX->Drc*DX->Drc) : 0.0;
                housedrain = std::max(0.0, housedrain - (DS - DStor->Drc));
            }
            else
                DS = 0;

            DStor->Drc = DS;
            // put new storage back in maps in m

            IntercHouse->Drc = HouseCover->Drc * (HS+DS) * SoilWidthDX->Drc * DX->Drc;
            // total interception in m3,exclude roads, channels

            RainNet->Drc = housedrain + (1-HouseCover->Drc)*RainNet->Drc;
            // net rainfall is direct rainfall + drainage
            // rainfall that falls on the soil, used in infiltration
        }
    }
}
//---------------------------------------------------------------------------
