
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
    {
        double CS = CStor->Drc;
        //actual canopy storage in m
        double Smax = CanopyStorage->Drc;
        //max canopy storage in m
        double LAIv;
        if (SwitchInterceptionLAI)
            LAIv = LAI->Drc;
        else
            LAIv = (log(1-Cover->Drc)/-0.4)/max(0.9,Cover->Drc);
        //Smax is based on LAI and LAI is the average of a gridcell, already including the cover
        // a low cover means a low LAI means little interception
        // avoid division by 0

        if (SwitchBuffers && !SwitchSedtrap)
            if(BufferID->Drc > 0)
                Smax = 0;
        // no interception with buffers, but sedtrap can have interception

        if (SwitchHardsurface)
            Smax *= (1-HardSurface->Drc);
        //VJ 110111 no interception on hard surfaces

        if (PlantHeight->Drc < WH->Drc)
        {
            Smax = 0;
            CS = 0;
        }
        //VJ no interception when water level is heigher than plants


        if (Smax > 0)
        {
            double k = 1-exp(-CanopyOpeness*LAIv);
            //VJ !!!!!!2013 05 18 BUG ! was k=exp(-coLAI) MUST BE 1-exp
            CS = Smax*(1-exp(-k*RainCum->Drc/Smax));
            //      CS = Smax*(1-exp(-0.0653*LAIv*RainCum->Drc/Smax));
            //VJ 110209 direct use of openess, astons value too open.
            //A good guess is using the cover LAI relation
            //and interpreting cover as openess: k = exp(-0.45*LAI)
        }
        else
            CS = 0;
        // 0.0653 is canopy openess, based on Aston (1979), based on Merriam (1960/1973), De Jong & Jetten 2003
        // is not the same as canopy cover. it also deals with how easy rainfall drips through the canopy
        // possible to use equation from Ahston but for very open Eucalypt

        CS = max(0, CS * (1-StemflowFraction));
        //VJ 110206 decrease storage with stemflow fraction!

        LeafDrain->Drc = max(0, Cover->Drc*(Rainc->Drc - (CS - CStor->Drc)));
        // diff between new and old strage is subtracted from rainfall
        // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!
        // note: cover already implicit in LAI and Smax, part falling on LAI is cover*rainfall

        CStor->Drc = CS;
        // put new storage back in map
        Interc->Drc =  Cover->Drc * CS * SoilWidthDX->Drc * DX->Drc; //*
        // only on soil surface, not channels or roads, in m3
        // cover already implicit in CS, Smax

        RainNet->Drc = LeafDrain->Drc + (1-Cover->Drc)*Rainc->Drc;
        // net rainfall is direct rainfall + drainage
        // rainfall that falls on the soil, used in infiltration
    }
}
//---------------------------------------------------------------------------
void TWorld::InterceptionHouses(void)
{
    // all variables are in m
    if (!SwitchHouses)
        return;

    FOR_ROW_COL_MV
    {
        if (HouseCover->Drc > 0)
        {
            double HS, DS;
            //actual roof storage in m
            double Hmax = RoofStore->Drc;
            //max roof storage in m
            // HouseCover->Drc = qMin(HouseCover->Drc, 0.95);

            double Dmax =0;
            if (SwitchRaindrum)
                if (HouseCover->Drc > 0)
                    Dmax = DrumStore->Drc/(CellArea->Drc*HouseCover->Drc);
            //max drum storage in m

            double housedrain = 0;
            //overflow in m

            if (Hmax > 0)
            {
                double k = 1.0;
                // speed of filling of roof storage, very quickly
                HS = Hmax*(1-exp(-k*RainCum->Drc/Hmax));
                //roof storage in m
            }
            else
                HS = 0;

            if (Dmax > 0)
            {
                double k = 0.05;
                // speed of filling of raindrum near house, slower
                DS = Dmax*(1-exp(-k*RainCum->Drc/Dmax));
                //drum storage in m

            }
            else
                DS = 0;

            housedrain = max(0, (RainNet->Drc - (HS - HStor->Drc) - (DS - DStor->Drc)));
            // diff between new and old strage is subtracted from rainfall
            // rest reaches the soil surface. ASSUMPTION: with the same intensity as the rainfall!

            HStor->Drc = HS;
            DStor->Drc = DS;
            // put new storage back in maps in m and m

            IntercHouse->Drc = HouseCover->Drc * (HS+DS) * SoilWidthDX->Drc * DX->Drc;//_dx*_dx
            // total interception in m3

            RainNet->Drc = HouseCover->Drc*housedrain + (1-HouseCover->Drc)*RainNet->Drc;
            // net rainfall is direct rainfall + drainage
            // rainfall that falls on the soil, used in infiltration
        }
        else
            IntercHouse->Drc = 0;
    }
}
//---------------------------------------------------------------------------
