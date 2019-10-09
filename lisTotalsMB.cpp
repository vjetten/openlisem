
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


#include <algorithm>
#include "model.h"
#include "operation.h"


//---------------------------------------------------------------------------
// totals for screen and file output and mass balance
void TWorld::Totals(void)
{
    double rainfall, snowmelt;
    double oldrainpeak, oldsnowpeak;
    double catchmentAreaFlatMM = 1000.0/(_dx*_dx*nrCells);  

    /***** WATER *****/

    FOR_ROW_COL_MV
    {
        if (SwitchKinematic2D < 3)
              WHmax->Drc = std::max(WHmax->Drc, WH->Drc + hmx->Drc);
        else
              WHmax->Drc = std::max(WHmax->Drc, std::max(WH->Drc, hmx->Drc));
    }


    //=== precipitation ===//
    if (SwitchRainfall)
    {
        RainAvgmm = mapAverage(*Rain)*1000.0;
        RainTotmm += RainAvgmm;
        // spatial avg area rainfall in mm

        calcMapValue(*tm, *Rain, (_dx*_dx), MUL); //in m3
        rainfall = mapTotal(*tm);
        RainTot += rainfall; // in m3

        oldrainpeak  = Rainpeak;
        Rainpeak = std::max(Rainpeak, rainfall);
        if (oldrainpeak  < Rainpeak)
            RainpeakTime = time;
    }

    if (SwitchSnowmelt)
    {
        SnowAvgmm = mapAverage(*Snowmelt)*1000;
        SnowTotmm += SnowAvgmm;

        calcMapValue(*tm, *Snowmelt, (_dx*_dx), MUL); //in m3
        snowmelt = mapTotal(*tm);
        SnowTot += snowmelt; // in m3

        oldsnowpeak = Snowpeak;
        Snowpeak = std::max(Snowpeak, snowmelt);
        if (oldsnowpeak < Snowpeak)
            SnowpeakTime = time;
    }
    
    //=== interception ===//
    IntercTot = mapTotal(*Interc);
    IntercTotmm = IntercTot*catchmentAreaFlatMM;
    // interception in mm and m3

    //Litter
    IntercLitterTot = mapTotal(*LInterc);
    IntercLitterTotmm = IntercLitterTot*catchmentAreaFlatMM;
    //houses
    IntercHouseTot = mapTotal(*IntercHouse);
    IntercHouseTotmm = IntercHouseTot*catchmentAreaFlatMM;
    // interception in mm and m3

    FOR_ROW_COL_MV
    {
        InterceptionmmCum->Drc = (Interc->Drc + IntercHouse->Drc + LInterc->Drc)*1000.0/CellArea->Drc;
        // for screen output only
    }

    //=== infiltration ===//
    InfilTot += mapTotal(*InfilVol) + mapTotal(*InfilVolKinWave);// + mapTotal(*InfilVolFlood); //m3

    InfilKWTot += mapTotal(*InfilVolKinWave); // not really used, available for output when needed
    InfilTotmm = std::max(0.0 ,(InfilTot)*catchmentAreaFlatMM);
    // infiltration mm and m3

    // flood infil
    // used for reporting only
    FOR_ROW_COL_MV
    {
        InfilVolCum->Drc += InfilVol->Drc + InfilVolKinWave->Drc + InfilVolFlood->Drc;
        // infilvolflood is 0 when 2Ddyn
        InfilmmCum->Drc = std::max(0.0, InfilVolCum->Drc*1000.0/(_dx*_dx));
        PercmmCum->Drc += Perc->Drc*1000.0;
    }

    //=== surf store ===//
    calcMapValue(*tm, *WHstore, 1000, MUL); //mm
    SurfStoremm = mapAverage(*tm);
    // surface storage CHECK THIS
    // does not go to MB, is already in tot water vol

    //=== surface flow ===//   
    WaterVolTot = mapTotal(*WaterVolall);//m3
    WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm

    // split into runoff and flood for screen reporting
    WaterVolRunoffmm = 0;
    WaterVolRunoffmm_F = 0;
    FOR_ROW_COL_MV
    {
      //  WaterVolRunoffmm += WHrunoff->Drc * ChannelAdj->Drc * DX->Drc;
        WaterVolRunoffmm += WHrunoff->Drc * FlowWidth->Drc * DX->Drc;
        WaterVolRunoffmm_F += std::min(WHrunoff->Drc,minReportFloodHeight) * FlowWidth->Drc * DX->Drc;
    }
    if (SwitchKinematic2D == K2D_METHOD_DYN || SwitchKinematic2D == K2D_METHOD_DIFF)
       WaterVolRunoffmm = WaterVolRunoffmm_F;
    WaterVolRunoffmm *= catchmentAreaFlatMM;
    // water on the surface in runoff in m3 and mm

    // runoff fraction per cell calc as in-out/rainfall, indication of sinks and sources of runoff
    // exclude channel cells
    FOR_ROW_COL_MV {
        runoffTotalCell->Drc += Qn->Drc * _dt * catchmentAreaFlatMM; // in mm !!!!
    }

    //=== storm drain flow ===//
    StormDrainVolTot = mapTotal(*TileWaterVol);
    StormDrainTotmm = StormDrainVolTot*catchmentAreaFlatMM;

    //=== channel flow ===//
    if (SwitchIncludeChannel)
    {
        ChannelVolTot = mapTotal(*ChannelWaterVol); //m3
        // add channel vol to total
        ChannelVolTotmm = ChannelVolTot*catchmentAreaFlatMM; //mm
        // recalc in mm for screen output
    }

    //=== flood  ===//
    if (SwitchChannelFlood)
    {
        floodVolTot = mapTotal(*FloodWaterVol);
        floodVolTotmm = floodVolTot * catchmentAreaFlatMM; // to mm
//
//        if (SwitchFloodInitial && runstep == 1)
//            floodVolTotInit = floodVolTot;
        // save initial flood level for mass balance if start with flood
    }

    // we have now:
    // WaterVolTot and WaterVolTotmm
    // ChannelVolTot and ChannelVolTotmm
    // floodVolTot and floodVolTotmm

    //=== all discharges ===//
    QtotT = 0;
    // sum all outflow in m3 for this timestep, Qtot is for all timesteps!

    // Add overland flow
    if(SwitchKinematic2D == K2D_METHOD_KIN)
    {
        FOR_ROW_COL_MV
        {
            if (LDD->Drc == 5)
                QtotT += Qn->Drc*_dt;
        }
    }
    else
    {
        QtotT += K2DQOut; // NOTE: K2DQOut already based on if(K2DOutlets->Drc == 1)
    }
    // sum outflow m3 for all timesteps for all outlets, in m3
    // needed for mass balance

    // add channel outflow
    if (SwitchIncludeChannel)
    {
        FOR_ROW_COL_MV_CH
        {
            if (LDDChannel->Drc == 5)
                QtotT += ChannelQn->Drc*_dt; //m3

            ChannelQntot->Drc += ChannelQn->Drc*_dt;
            //cumulative m3 spatial for .map output
        }
        // add channel outflow (in m3) to total for all pits
        // recalc in mm for screen output
    }

  if(SwitchIncludeStormDrains) {
    FOR_ROW_COL_MV_TILE
        if (LDDTile->Drc == 5)
        {
            QtotT += TileQn->Drc * _dt;
            QTiletot += TileQn->Drc * _dt;
        }
  }
    if (SwitchIncludeTile)
    {
        WaterVolSoilTot = mapTotal(*TileWaterVolSoil);
        // input for mass balance, is the water seeping from the soil, input
        // this is the water before the kin wave
        calc2Maps(*tm, *TileDrainSoil, *TileWidth, MUL); //in m3
        calcMap(*tm, *TileDX, MUL); //in m3
        // tm->calcV(_dx, MUL); //in m3 ??? or DX?
        TileVolTot += mapTotal(*tm); // in m3

        // water after kin wave
        //WaterVolTot += mapTotal(*TileWaterVol); //m3
        // add tile vol to total
        //WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
        // recalc in mm for screen output

        FOR_ROW_COL_MV_TILE
            if (LDDTile->Drc == 5)
            {
                QtotT += TileQn->Drc * _dt;
                QTiletot += TileQn->Drc * _dt;
            }

        // add tile outflow (in m3) to total for all pits

        //QtotOutlet += TileQn->DrcOutlet * _dt;  obsolete, now Qtot-
        // add Tile outflow (in m3) to total for main outlet
    }

    // output fluxes for reporting to file and screen in l/s!
    FOR_ROW_COL_MV
    {
        Qoutput->Drc = 1000*(Qn->Drc + ChannelQn->Drc);// + TileQn->Drc); // in l/s
        //NOTE: for 2D flow Qn = K2DQ, already done
    }

    Qtot += QtotT;

    // Total outflow in m3 for all timesteps
    // does NOT include flood water leaving domain (floodBoundaryTot)
    Qtotmm = Qtot*catchmentAreaFlatMM;
    // recalc to mm for screen output

    floodBoundaryTot += K2DQOutBoun;
    FloodBoundarymm = floodBoundaryTot*catchmentAreaFlatMM;
    // flood boundary losses are done separately in MB

    //=====***** SEDIMENT *****====//

    // DetSplashTot, DetFlowTot and DepTot are for output in file and screen
    // DetTot and DepTot are for MB

    SoilLossTotT = 0;

    if (SwitchErosion)
    {

        // Dep and Detflow are zero if 2Ddyn
        DetSplashTot += mapTotal(*DETSplash);
        DetFlowTot += mapTotal(*DETFlow);
        DepTot += mapTotal(*DEP);
        DetTot = DetFlowTot + DetSplashTot;
        SedTot = mapTotal(*Sed);        
        // all in kg/cell

        calcMap(*DETSplashCum, *DETSplash, ADD);
        calcMap(*DETFlowCum, *DETFlow, ADD);
        calcMap(*DEPCum, *DEP, ADD);
        // DEP is set to 0 each timestep
        // for total soil loss calculation: TotalSoillossMap

        //outflow from domain/channel
        if(SwitchKinematic2D == K2D_METHOD_KIN) {
            FOR_ROW_COL_MV
            {
                if (LDD->Drc == 5)
                    SoilLossTotT += Qsn->Drc * _dt;
            }
        } else {
            SoilLossTotT += K2DQSOut; // for dyn wave this is the boundary outflow
        }
        // sum all sed in all boundary (in kg), needed for mass balance

        if (SwitchIncludeChannel)
        {
            FOR_ROW_COL_MV_CH
            {
                if (LDDChannel->Drc == 5)
                    SoilLossTotT += ChannelQsn->Drc * _dt;
            }

            // units here in kg, conversion to ton in report functions
            ChannelDetTot += mapTotal(*ChannelDetFlow);
            ChannelDepTot += mapTotal(*ChannelDep);
            ChannelSedTot = mapTotal(*ChannelBLSed) + mapTotal(*ChannelSSSed);
        }
        // used for mass balance and screen output
        FloodDetTot += mapTotal(*BLDetFlood) + mapTotal(*SSDetFlood);
        FloodDepTot += mapTotal(*DepFlood);
        FloodSedTot = mapTotal(*BLFlood) + mapTotal(*SSFlood);

        calcMap(*DETFlowCum, *BLDetFlood, ADD);
        calcMap(*DETFlowCum, *SSDetFlood, ADD);
        calcMap(*DEPCum, *DepFlood, ADD);

        // SPATIAL totals for output overland flow all in kg/cell
        // variables are valid for both 1D and 2D flow dyn and diff
        FOR_ROW_COL_MV
        {
            Qsoutput->Drc = Qsn->Drc + ChannelQsn->Drc;  // in kg/s
            // for reporting sed discharge screen
        }

        floodBoundarySedTot += K2DQSOutBoun;
        SoilLossTotT += K2DQSOutBoun;
        //FloodBoundarySedmm = floodBoundarySedTot*catchmentAreaFlatMM;
        // flood boundary losses are done separately in MB

        // for reporting
        if (SwitchIncludeChannel)
        {
            fill(*tma,0.0);
            DistributeOverExtendedChannel(ChannelDetFlow,tma);
            fill(*tmb,0.0);
            DistributeOverExtendedChannel(ChannelDep,tmb);
            FOR_ROW_COL_MV
            {
                DETFlowCum->Drc += tma->Drc;
                DEPCum->Drc += tmb->Drc;
                TotalChanDetMap->Drc += ChannelDetFlow->Drc;
                TotalChanDepMap->Drc += ChannelDep->Drc;
            }
        }


        // with all det and dep calculate the soil loss, excl channel
        FOR_ROW_COL_MV
        {
            TotalSoillossMap->Drc = DETSplashCum->Drc + DETFlowCum->Drc + DEPCum->Drc;
                    //TotalDetMap->Drc + TotalDepMap->Drc;
                  // + TotalChanDetMap->Drc + TotalChanDepMap->Drc;
//            TotalDepMap->Drc = std::min(0.0, TotalSoillossMap->Drc);
//            TotalDetMap->Drc = std::max(0.0, TotalSoillossMap->Drc);
        }

        FOR_ROW_COL_MV
        {
            //double Q = Qoutput->Drc/1000;
        //    TotalConc->Drc = (Q > MIN_FLUX ? Qsoutput->Drc/Q : 0);
            // conc is okay in principle, Qs is conc (mass/vol)*Q so here this is reversed
            double sedall = Sed->Drc + BLFlood->Drc + SSFlood->Drc + ChannelSed->Drc;
            double waterall = WaterVolall->Drc + ChannelWaterVol->Drc;
            TotalConc->Drc = MaxConcentration(waterall ,sedall);
            // for output
        }

        fill(*DepFlood,0.0);
        fill(*BLDetFlood,0.0);
        fill(*SSDetFlood,0.0);
        // RESET flood variables (?)

        SoilLossTot += SoilLossTotT;
        // total sediment outflow from outlets and domain boundaries

    }


    if (SwitchPesticide)
    {
        FOR_ROW_COL_MV
        {
            // = WHoutavg->Drc*_dx*DX->Drc*C->Drc*1000*1000*1000; //µg
            PDisMixing->Drc = CM->Drc*epsil->Drc*poro->Drc*_dx*_dx*1000*1000*1000; //µg
            PSorMixing->Drc = CS->Drc*epsil->Drc*rhob->Drc*_dx*_dx*1000*1000*1000; //µg
            PInfilt->Drc = pestiinf->Drc*CM->Drc*_dx*_dx*_dt*1000*1000*1000; //µg
            PStorage->Drc= WHstore->Drc*_dx*_dx*C->Drc*1000*1000*1000; //µg
            PRunoffSpatial->Drc = Pest->Drc*1000*1000*1000; //µg

            //            PRunoffSpatialex->Drc= WHoutavg->Drc*_dx*DX->Drc*C_Kexplicit->Drc*1000*1000*1000; //µg
            //            PDisMixingex->Drc = CM_Kexplicit->Drc*epsil->Drc*poro->Drc*_dx*DX->Drc*1000*1000*1000; //µg
            //            PSorMixingex->Drc = CS_Kexplicit->Drc*epsil->Drc*rhob->Drc*_dx*DX->Drc*1000*1000*1000; //µg
            //            PInfiltex->Drc = pestiinf->Drc*CM_Kexplicit->Drc*_dx*DX->Drc*_dt*1000*1000*1000; //µg

        }

        Pestdetach += mapTotal(*Pdetach); //KCM
        PestCinfilt += mapTotal(*PCinfilt); //fc
        PestCfilmexit += mapTotal(*PCfilmexit); //KC
        //PestLossTotOutlet += Qn->DrcOutlet*C->DrcOutlet*_dt*1000*1000*1000; //µg  //DrcOutlet obsolete
        PestRunoffSpatial = mapTotal(*PRunoffSpatial);
        PestDisMixing = mapTotal(*PDisMixing);
        PestSorMixing = mapTotal(*PSorMixing);
        PestInfilt += mapTotal(*PInfilt);
        PestStorage = mapTotal(*PStorage);

    }

    SedimentSetMaterialDistribution();

}
//---------------------------------------------------------------------------
void TWorld::MassBalance()
{
    // Mass Balance water, all in m3
    // VJ 110420 added tile volume here, this is the input volume coming from the soil after swatre
  //  if (RainTot + SnowTot > 0)
    {

        double waterin = RainTot + SnowTot + WaterVolSoilTot + BaseFlow + WHinitVolTot; //floodVolTotInit
        double waterstore = IntercTot + IntercLitterTot + IntercHouseTot + InfilTot;
        double waterflow = WaterVolTot + ChannelVolTot + StormDrainVolTot + Qtot + floodBoundaryTot ;//

        if(SwitchKinematic2D == K2D_METHOD_DYN
           || (SwitchKinematic2D != K2D_METHOD_DYN && !SwitchIncludeChannel) ) {}
            else
                waterflow += floodVolTot;

        MB = waterin > 0 ? (waterin - waterstore - waterflow)/waterin *100 : 0;
      //  qDebug() << MB << waterin << waterstore << waterflow<< floodBoundaryTot;

    }
    //watervoltot includes channel and tile

    // Mass Balance sediment, all in kg

    //VJ 110825 forgot to include channeldettot in denominator in MBs!
    if (SwitchErosion)// && SoilLossTot > 1e-9)
    {
        double detachment = DetTot + ChannelDetTot + FloodDetTot;
        double deposition = DepTot + ChannelDepTot + FloodDepTot;
        double sediment = SedTot + ChannelSedTot + FloodSedTot + SoilLossTot;
        //already in soiloss: + floodBoundarySedTot;

    //    qDebug() << "S" << DetTot<< ChannelDetTot << FloodDetTot;
    //    qDebug() << DepTot << ChannelDepTot << FloodDepTot;
     //   qDebug() << SedTot << ChannelSedTot << FloodSedTot << SoilLossTot;



        MBs = detachment > 0 ? (detachment + deposition  - sediment)/detachment*100 : 0;
            //    qDebug() << MBs << detachment << deposition << sediment;
    }
    //VJ 121212 changed to mass balance relative to soil loss

    if (SwitchPesticide)
    {
        MBp = (PestMassApplied-PestLossTotOutlet-PestRunoffSpatial-PestDisMixing-PestSorMixing-PestInfilt-PestStorage)*100/PestMassApplied;
        //MBpex = (PestMassApplied-PestLossTotOutletex-PestRunoffSpatialex-PestDisMixingex-PestSorMixingex-PestInfiltex)*100/PestMassApplied;
        //(PestMassApplied-PestLossTotOutlet-PestRunoffSpatial-PestDisMixing-PestSorMixing-PestInfilt-PestStorage)*100/PestMassApplied
        debug(QString("mbp: %1").arg(MBp));
    }
}
//---------------------------------------------------------------------------
