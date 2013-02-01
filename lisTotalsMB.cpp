
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
  \file lisTotalsMB.cpp
  \brief calculate water and sediment flux totals and mass balance

functions: \n
- void TWorld::Totals() calculate the totals for all water and fluxes for mss balance and output\n
- void TWorld::MassBalance() water and sediment mass balance\n
 */


#include "model.h"


//---------------------------------------------------------------------------
// totals for screen and file output and mass balance
void TWorld::Totals(void)
{
    double rainfall, snowmelt;
    double oldrainpeak, oldsnowpeak;
    double catchmentAreaFlatMM = 1000.0/(_dx*_dx*nrCells);

    /***** WATER *****/

    if (SwitchRainfall)
    {
        RainAvgmm = Rain->mapAverage()*1000.0;
        RainTotmm += RainAvgmm;
        // avg area rainfall in mm

        tm->calcMapValue(Rain, (_dx*_dx), MUL); //in m3
        rainfall = tm->mapTotal();
        RainTot += rainfall; // in m3

        oldrainpeak = Rainpeak;
        Rainpeak = max(Rainpeak, rainfall);
        if (oldrainpeak  < Rainpeak)
            RainpeakTime = time;
    }

    if (SwitchSnowmelt)
    {
        SnowAvgmm += Snowmelt->mapAverage()*1000;
        SnowTotmm += SnowAvgmm;

        tm->calcMapValue(Snowmelt, (_dx*_dx), MUL); //in m3
        snowmelt = tm->mapTotal();
        SnowTot += snowmelt; // in m3

        oldsnowpeak = Snowpeak;
        Snowpeak = max(Snowpeak, snowmelt);
        if (oldsnowpeak < Snowpeak)
            SnowpeakTime = time;
    }

    IntercTot = Interc->mapTotal();
    IntercTotmm = IntercTot*catchmentAreaFlatMM;
    // interception in mm and m3

    //houses
    IntercHouseTot = IntercHouse->mapTotal();
    IntercHouseTotmm = IntercHouseTot*catchmentAreaFlatMM;
    // interception in mm and m3

    InfilTot += InfilVol->mapTotal() + InfilVolKinWave->mapTotal(); //m3
    InfilKWTot += InfilVolKinWave->mapTotal(); // not really used, available for output when needed
    InfilTotmm = max(0,InfilTot*catchmentAreaFlatMM);
    // infiltration mm and m3

    tm->calcMapValue(WHstore, 1000, MUL); //mm
    SurfStoremm = tm->mapAverage();
    // surface storage CHECK THIS
    // does not go to MB, is already in tot water vol

    WaterVolTot = WaterVolall->mapTotal();//m3

    // replace water volume in flooded areas with flood volume for water balance
    if (SwitchChannelFlood)
    {
//        tma->fill(0);
//        FOR_ROW_COL_MV
//        {
//            if (FloodDomain->Drc == 1)
//                tma->Drc = WaterVolall->Drc;
//        }
//        WaterVolTot -= tma->mapTotal();
        //WaterVolTot += FloodWaterVol->mapTotal();
    }

    WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
    // water on the surface in runoff in m3 and mm
    //NOTE: surface storage is already in here so does not need to be accounted for in MB

    // sum outflow m3 for all timesteps for the outlet, is already mult by dt
    FOR_ROW_COL_MV
            if (LDD->Drc == 5)
            Qtot += Qn->Drc*_dt;
    // sum outflow m3 for all timesteps for the outlet, in m3
    // needed for mass balance
    //Qtotmm = Qtot*catchmentAreaFlatMM;
    // in mm for screen output

    QtotOutlet += Qn->DrcOutlet*_dt;
    // for screen output, total main outlet in m3
    QtotPlot += Qn->DrcPlot * _dt;
    //VJ 110701 for screen output, total of hydrograph point in m3

    TotalWatervol->copy(WaterVolall);
    // for sed conc calc output

    if (SwitchIncludeChannel)
    {
        WaterVolTot += ChannelWaterVol->mapTotal(); //m3
        // add channel vol to total
        WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
        // recalc in mm for screen output
        FOR_ROW_COL_MV_CH
                if (LDDChannel->Drc == 5)
                Qtot += ChannelQn->Drc*_dt;
        // add channel outflow (in m3) to total for all pits
        //Qtotmm = Qtot*catchmentAreaFlatMM;
        // recalc in mm for screen output

        QtotOutlet += ChannelQn->DrcOutlet * _dt;
        // add channel outflow (in m3) to total for main outlet
        QtotPlot += ChannelQn->DrcPlot * _dt;
        // add channel outflow (in m3) to total for main outlet
        TotalWatervol->calcMap(ChannelWaterVol,ADD);
        // add channel volume to total for sed conc calc

    }

    if (SwitchIncludeTile)
    {
        WaterVolSoilTot = TileWaterVolSoil->mapTotal();
        // input for mass balance, is the water seeping from the soil, input
        // this is the water before the kin wave
        tm->calc2Maps(TileDrainSoil, TileWidth, MUL); //in m3
        tm->calcMap(TileDX, MUL); //in m3
        // tm->calcV(_dx, MUL); //in m3 ??? or DX?
        TileVolTot += tm->mapTotal(); // in m3

        // water after kin wave
        WaterVolTot += TileWaterVol->mapTotal(); //m3
        // add tile vol to total
        WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
        // recalc in mm for screen output

        //Qtot += TileQoutflow->DrcOutlet;
        FOR_ROW_COL_MV_TILE
                if (LDDTile->Drc == 5)
                Qtot += TileQn->Drc * _dt;
        // add tile outflow (in m3) to total for all pits
        //Qtotmm = Qtot*catchmentAreaFlatMM;
        // recalc in mm for screen output

        QtotOutlet += TileQn->DrcOutlet * _dt;
        // add channel outflow (in m3) to total for main outlet
        QtotPlot += TileQn->DrcPlot * _dt;
        // add channel outflow (in m3) to total for subcatch outlet
        TotalWatervol->calcMap(TileWaterVol,ADD);
        // add channel volume to total for sed conc calc

    }

    if (SwitchBuffers)
    {
        BufferVolTot = BufferVol->mapTotal(); // in m3
        if (SwitchIncludeChannel)
            BufferVolTot += ChannelBufferVol->mapTotal();
        //sum up all volume remaining in all buffers (so the non-water!)
        BufferVolin = BufferVolTotInit - BufferVolTot;
        //subtract this from the initial volume to get the total water inflow in the buffers
    }

    // output fluxes for reporting to file and screen in l/s!
    FOR_ROW_COL_MV
    {
        Qoutput->Drc = 1000*qMax(1e-6, Qn->Drc + ChannelQn->Drc + TileQn->Drc); // in l/s
        // added minimum here to avoid strange maps
//        if (SwitchChannelFlood)
//            if (hmx->Drc > 0)
//        {
//                Qoutput->Drc = Qflood->Drc*1000;
//        }
    }

    Qtotmm = Qtot*catchmentAreaFlatMM;
    // recalc to mm for screen output

    oldrainpeak = Qpeak;
    Qpeak = max(Qpeak, Qoutput->DrcOutlet);
    if (oldrainpeak < Qpeak)
        QpeakTime = time;
    // peak flow and peak time calculation, based on sum channel and runoff

    /***** SEDIMENT *****/

    if (SwitchErosion)
    {
        DetSplashTot += DETSplash->mapTotal();
        DetFlowTot += DETFlow->mapTotal();
        DepTot += DEP->mapTotal();
        DetTot += DETSplash->mapTotal() + DETFlow->mapTotal();
        SedTot = Sed->mapTotal();
        // all in kg/cell

        //SoilLossTot += Qsoutflow->DrcOutlet;
        FOR_ROW_COL_MV
                if (LDD->Drc == 5)
                SoilLossTot += Qsn->Drc * _dt;
        // sum all sed in all pits (in kg), needed for mass balance

        SoilLossTotOutlet += Qsn->DrcOutlet * _dt;
        // for screen output, total main outlet sed loss in kg
        TotalSed->copy(Sed);
        // for sed conc

        if (SwitchIncludeChannel)
        {
            // units here in kg, conversion to ton in report functions
            ChannelDetTot += ChannelDetFlow->mapTotal();
            ChannelDepTot += ChannelDep->mapTotal();
            ChannelSedTot = ChannelSed->mapTotal();

            FOR_ROW_COL_MV_CH
                    if (LDDChannel->Drc == 5)
                    SoilLossTot += ChannelQsn->Drc * _dt;
            // add sed outflow for all pits to total soil loss

            SoilLossTotOutlet += ChannelQsn->DrcOutlet * _dt;
            // add channel outflow (in kg) to total for main outlet

            TotalSed->calcMap(ChannelSed, ADD);
            // needed for sed conc in file output
        }

        FOR_ROW_COL_MV
        {
            TotalConc->Drc = MaxConcentration(TotalWatervol->Drc, TotalSed->Drc);
        }
        // for file output

        if (SwitchBuffers || SwitchSedtrap)
        {
            BufferSedTot = BufferSed->mapTotal();
            if (SwitchIncludeChannel)
                BufferSedTot += ChannelBufferSed->mapTotal();

        }
        /** TODO add gully, wheeltracks etc */

        // spatial totals for output all in kg/cell
        FOR_ROW_COL_MV
        {
            Qsoutput->Drc = Qsn->Drc + ChannelQsn->Drc;  // sum channel and OF sed output in kg/s

            TotalDetMap->Drc += DETSplash->Drc + DETFlow->Drc;
            TotalDepMap->Drc += DEP->Drc;
            if (SwitchIncludeChannel)
            {
                TotalDetMap->Drc += ChannelDetFlow->Drc;
                TotalDepMap->Drc += ChannelDep->Drc;
            }
            TotalSoillossMap->Drc = TotalDetMap->Drc + TotalDepMap->Drc;
        }
    }
}
//---------------------------------------------------------------------------
void TWorld::MassBalance()
{
    // Mass Balance water, all in m3
    // VJ 110420 added tile volume here, this is the input volume coming from the soil after swatre
    if (RainTot + SnowTot > 0)
        MB = (RainTot + SnowTot + WaterVolSoilTot
              - IntercTot - IntercHouseTot - InfilTot - WaterVolTot - Qtot - BufferVolin)/
                (RainTot + SnowTot + WaterVolSoilTot)*100;
    //watervoltot includes channel and tile

    //  qDebug() << RainTot << IntercTot << IntercHouseTot << InfilTot << WaterVolTot << BufferVolin << Qtot<< InfilKWTot;

    // Mass Balance sediment, all in kg
    //   if (SwitchErosion && (DetTot + ChannelDetTot) > 0)
    //      MBs = (DetTot + ChannelDetTot - SoilLossTot - SedTot - ChannelSedTot +
    //             DepTot + ChannelDepTot - BufferSedTot)/(DetTot + ChannelDetTot)*100;
    //VJ 110825 forgot to include channeldettot in denominator in MBs!
    if (SwitchErosion && SoilLossTot > 1e-9)
        MBs = (DetTot + ChannelDetTot - SedTot - ChannelSedTot +
               DepTot + ChannelDepTot - BufferSedTot)/(SoilLossTot) *100;
    //VJ 121212 changed to mass balance relative to soil loss
}
//---------------------------------------------------------------------------
