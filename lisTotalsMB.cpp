
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
  \file lisTotalsMB.cpp
  \brief calculate water and sediment flux totals and mass balance

functions: \n
- void TWorld::Totals() calculate the totals for all water and fluxes for mss balance and output\n
- void TWorld::MassBalance() water and sediment mass balance\n
 */


#include <algorithm>
#include "model.h"
#include "operation.h"

double TWorld::MapTotal(cTMap &M)
{
    double total = 0;
    #pragma omp parallel for reduction(+:total) num_threads(userCores)
    FOR_ROW_COL_MV_L {
       total = total + M.Drc;
    }}
    return (total);
}
//---------------------------------------------------------------------------
// totals for screen and file output and mass balance
void TWorld::Totals(void)
{
    double rainfall, snowmelt;
    double oldrainpeak, oldsnowpeak;
    double catchmentAreaFlatMM = 1000.0/(_dx*_dx*nrCells);

    /***** WATER *****/

    //=== precipitation ===//
    if (SwitchRainfall)
    {
        RainAvgmm = MapTotal(*Rain)*1000.0/(double)nrValidCells;

        RainTotmm += RainAvgmm;
        // spatial avg area rainfall in mm

        rainfall = RainAvgmm/1000.0*_dx*_dx*nrCells;
        RainTot += rainfall; // in m3

        oldrainpeak  = Rainpeak;
        Rainpeak = std::max(Rainpeak, rainfall);
        if (oldrainpeak  < Rainpeak)
            RainpeakTime = time;
    }

    if (SwitchIncludeET) {
       // double ETtot = MapTotal(*ETa);
        ETaTot = MapTotal(*ETaCum);


    }

    if (SwitchSnowmelt)
    {
        SnowAvgmm = MapTotal(*Snowmelt)*1000.0/nrCells;

        SnowTotmm += SnowAvgmm;

        snowmelt = SnowAvgmm/1000.0*_dx*_dx*nrCells;
        SnowTot += snowmelt; // in m3

        oldsnowpeak = Snowpeak;
        Snowpeak = std::max(Snowpeak, snowmelt);
        if (oldsnowpeak < Snowpeak)
            SnowpeakTime = time;
    }

    //=== interception ===//
    IntercTot = MapTotal(*Interc);
    IntercTotmm = IntercTot*catchmentAreaFlatMM;
    // interception in mm and m3
    //Litter
    if (SwitchLitter) {
        IntercLitterTot = MapTotal(*LInterc);
        IntercLitterTotmm = IntercLitterTot*catchmentAreaFlatMM;
    }

    if (SwitchHouses) {
        IntercHouseTot = MapTotal(*IntercHouse);
        IntercHouseTotmm = IntercHouseTot*catchmentAreaFlatMM;
        // interception in mm and m3
    }
#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        InterceptionmmCum->Drc = (Interc->Drc + IntercHouse->Drc + LInterc->Drc)*1000.0/CellArea->Drc;
        // for screen and file output
    }}


    //==== ETa ==========//
    //   ETaTot = mapTotal(*ETa);
    //  ETaTotmm = ETaTot*catchmentAreaFlatMM;
    // interception in mm and m3

    //=== infiltration ===//
    InfilTot += MapTotal(*InfilVol) + MapTotal(*InfilVolKinWave);
report(*InfilVol,"iv");
    if (SwitchIncludeChannel && SwitchChannelInfil) {
        InfilTot += MapTotal(*ChannelInfilVol); //m3
    }

    InfilKWTot += MapTotal(*InfilVolKinWave); // not really used, available for output when needed
    InfilTotmm = std::max(0.0 ,(InfilTot)*catchmentAreaFlatMM);
    // infiltration mm and m3

    // flood infil
    // used for reporting only
#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        InfilVolCum->Drc += InfilVol->Drc + InfilVolKinWave->Drc;// + InfilVolFlood->Drc;
        if (SwitchIncludeChannel)
            InfilVolCum->Drc += ChannelInfilVol->Drc;
        InfilmmCum->Drc = std::max(0.0, InfilVolCum->Drc*1000.0/(_dx*_dx));
        PercmmCum->Drc += Perc->Drc*1000.0;
    }}

    //=== surf store ===//

    double SStot = 0;
    #pragma omp parallel for reduction(+:SStot) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        SStot += WHstore->Drc * SoilWidthDX->Drc*DX->Drc;
    }}
    SurfStoremm = SStot * catchmentAreaFlatMM;

    // does not go to MB, is already in tot water vol

    //TODO: check init WH
    //ONLY ONCE?
//    if (SwitchFloodInitial) {
//        WHinitVolTot = 0;
//#pragma omp parallel for reduction(+:WHinitVolTot) num_threads(userCores)
//        FOR_ROW_COL_MV_L {
//            WHinitVolTot = hmxInit->Drc * DX->Drc * ChannelAdj->Drc;
//        }}
//    }

    //=== surface flow ===//
    WaterVolTot = MapTotal(*WaterVolall);//m3
    WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; // => not used
    //All water on surface

    floodVolTot = MapTotal(*FloodWaterVol);
    floodVolTotmm = floodVolTot * catchmentAreaFlatMM; // to mm
    // flood water above user defined threshold

    if (SwitchKinematic2D == K2D_METHOD_DYN || SwitchKinematic2D == K2D_METHOD_KINDYN) {
       WaterVolRunoffmm = MapTotal(*RunoffWaterVol)* catchmentAreaFlatMM;//m3
       // runoff water below user defined threshold
    } else {
        WaterVolRunoffmm = 0;
        #pragma omp parallel for reduction(+:WaterVolRunoffmm) num_threads(userCores)
        FOR_ROW_COL_MV_L {
            WaterVolRunoffmm += WHrunoff->Drc * CHAdjDX->Drc;
        }}
        WaterVolRunoffmm *= catchmentAreaFlatMM;
    }
    // water on the surface in runoff in mm


    // runoff fraction per cell calc as in-out/rainfall, indication of sinks and sources of runoff
    // exclude channel cells
        #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        runoffTotalCell->Drc += (Qn->Drc +Qflood->Drc)* _dt * catchmentAreaFlatMM; // in mm !!!!
    }}

    //=== storm drain flow ===//
    StormDrainVolTot = MapTotal(*TileWaterVol);
    StormDrainTotmm = StormDrainVolTot*catchmentAreaFlatMM;

    //=== channel flow ===//
    if (SwitchIncludeChannel)
    {
        ChannelVolTot = MapTotal(*ChannelWaterVol); //m3
        // add channel vol to total
      //  if (SwitchChannelBaseflow)
       //     BaseFlowTot = MapTotal(*BaseFlowInflow) * _dt;
        ChannelVolTotmm = ChannelVolTot*catchmentAreaFlatMM; //mm
        // recalc in mm for screen output
    }

    //=== all discharges ===//
    QtotT = 0;
    // sum all outflow in m3 for this timestep, Qtot is for all timesteps!

    // Add overland flow, do this even for 2D dyn wave
    if(SwitchKinematic2D != K2D_METHOD_DYN) {
        FOR_ROW_COL_LDD5 {
            QtotT += Qn->Drc*_dt;
        }}

    }

    Qfloodout = 0;
    if(SwitchKinematic2D == K2D_METHOD_KINDYN)
    {
//#pragma omp parallel for reduction(+:Qfloodout) num_threads(userCores)
        FOR_ROW_COL_LDD5 {
            Qfloodout += Qflood->Drc * _dt;
        }}

        QfloodoutTot += Qfloodout;
    }

    // add channel outflow
    if (SwitchIncludeChannel)
    {
   //     FOR_ROW_COL_MV_CH
     //   {
       //     if (LDDChannel->Drc == 5) {

//#pragma omp parallel for reduction(+:QtotT) num_threads(userCores)
            FOR_ROW_COL_LDDCH5 {
                QtotT += ChannelQn->Drc*_dt; //m3
            }}
#pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_CHL {
                ChannelQntot->Drc += ChannelQn->Drc*_dt;
                //cumulative m3 spatial for .map output
            }}
      //  }
        // add channel outflow (in m3) to total for all pits
    }


    if(SwitchIncludeStormDrains) {
        FOR_ROW_COL_MV_TILE
                if (LDDTile->Drc == 5)
        {
            QtotT += TileQn->Drc * _dt;
            QTiletot += TileQn->Drc * _dt;
        }
    } else {
        if (SwitchIncludeTile)
        {
            WaterVolSoilTot = MapTotal(*TileWaterVolSoil);
            // input for mass balance, is the water seeping from the soil, input
            // this is the water before the kin wave
            calc2Maps(*tm, *TileDrainSoil, *TileWidth, MUL); //in m3
            calcMap(*tm, *TileDX, MUL); //in m3
            // tm->calcV(_dx, MUL); //in m3 ??? or DX?
            TileVolTot += MapTotal(*tm); // in m3

            FOR_ROW_COL_MV_TILE
                    if (LDDTile->Drc == 5)
            {
                QtotT += TileQn->Drc * _dt;
                QTiletot += TileQn->Drc * _dt;
            }

            // add tile outflow (in m3) to total for all pits
        }
    }

    floodBoundaryTot += K2DQOutBoun*_dt;
    FloodBoundarymm = floodBoundaryTot*catchmentAreaFlatMM;
    // 2D boundary losses, ALWAYS INCLUDES LDD=5 and channelLDD=5
    QtotT += K2DQOutBoun*_dt;

    // output fluxes for reporting to file and screen in l/s!]
#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L
    {
        Qoutput->Drc = 1000.0*(Qn->Drc + (SwitchIncludeChannel ? ChannelQn->Drc : 0.0) + Qflood->Drc);// in l/s
        Qoutput->Drc = Qoutput->Drc < 1e-6 ? 0.0 : Qoutput->Drc;
    }}

    // Total outflow in m3 for all timesteps
    // does NOT include flood water leaving domain (floodBoundaryTot)

    Qtot += QtotT;
    // add timestep total to run total
    Qtotmm = Qtot*catchmentAreaFlatMM;
    // recalc to mm for screen output


    //=====***** SEDIMENT *****====//

    // DetSplashTot, DetFlowTot and DepTot are for output in file and screen
    // DetTot and DepTot are for MB

    SoilLossTotT = 0;

    if (SwitchErosion)
    {
        SedTot = 0;
#pragma omp parallel for reduction(+:DetSplashTot,DetFlowTot,DepTot,SedTot) num_threads(userCores)
        FOR_ROW_COL_MV_L {
        // Dep and Detflow are zero if 2Ddyn
            DetSplashTot += DETSplash->Drc;
            DetFlowTot += DETFlow->Drc;
            DepTot += DEP->Drc;
            SedTot += Sed->Drc;
        }}
        // all in kg/cell

        DetTot = DetFlowTot + DetSplashTot;

#pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            DETSplashCum->Drc += DETSplash->Drc;
            DETFlowCum->Drc += DETFlow->Drc;
            DEPCum->Drc += DEP->Drc;
        }}
        // DEP is set to 0 each timestep
        // for total soil loss calculation: TotalSoillossMap

        //outflow from domain/channel
        if(SwitchKinematic2D == K2D_METHOD_KIN || SwitchKinematic2D == K2D_METHOD_KINDYN)
        {
//            FOR_ROW_COL_MV
//            {
//                if (LDD->Drc == 5)
//                    SoilLossTotT += Qsn->Drc * _dt;
//            }

           // #pragma omp parallel for reduction(+:SoilLossTotT) num_threads(userCores)
            FOR_ROW_COL_LDD5 {
                SoilLossTotT += Qsn->Drc * _dt;
            }}

        }
//            else {
//            SoilLossTotT += K2DQSOut; // for dyn wave this is the boundary outflow
//        }

        if (SwitchIncludeChannel)
        {
//            FOR_ROW_COL_MV_CH
//            {
//                if (LDDChannel->Drc == 5)
           // #pragma omp parallel for reduction(+:SoilLossTotT) num_threads(userCores)
            FOR_ROW_COL_LDDCH5 {
                    SoilLossTotT += ChannelQsn->Drc * _dt;
            }}

            // units here in kg, conversion to ton in report functions
//#pragma omp parallel for reduction(+:DetSplashTot,DetFlowTot,DepTot,SedTot) num_threads(userCores)
//        FOR_ROW_COL_MV_L {
//        // Dep and Detflow are zero if 2Ddyn
//            DetSplashTot += DETSplash->Drc;
//            DetFlowTot += DETFlow->Drc;
//            DepTot += DEP->Drc;
//            SedTot += Sed->Drc;
//        }}

            ChannelDetTot += MapTotal(*ChannelDetFlow);
            ChannelDepTot += MapTotal(*ChannelDep);
            ChannelSedTot = (SwitchUse2Phase ? MapTotal(*ChannelBLSed) : 0.0) + MapTotal(*ChannelSSSed);
        }

        floodBoundarySedTot += K2DQSOutBoun;
        SoilLossTotT += K2DQSOutBoun;
        //FloodBoundarySedmm = floodBoundarySedTot*catchmentAreaFlatMM;
        // flood boundary losses are done separately in MB


        // used for mass balance and screen output
        FloodDetTot += (SwitchUse2Phase ? MapTotal(*BLDetFlood) : 0.0) + MapTotal(*SSDetFlood);
        FloodDepTot += mapTotal(*DepFlood);
        FloodSedTot = (SwitchUse2Phase ? MapTotal(*BLFlood) : 0.0) + MapTotal(*SSFlood);

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            if (SwitchUse2Phase)
                DETFlowCum->Drc += BLDetFlood->Drc;
            DETFlowCum->Drc += SSDetFlood->Drc;
            DEPCum->Drc += DepFlood->Drc;
        }}
        // SPATIAL totals for output overland flow all in kg/cell
        // variables are valid for both 1D and 2D flow dyn and diff
        FOR_ROW_COL_MV_L
        {
            Qsoutput->Drc = Qsn->Drc + (SwitchIncludeChannel ? ChannelQsn->Drc : 0.0) + K2DQ->Drc;  // in kg/s
            // for reporting sed discharge screen
        }}

        // for reporting
        if (SwitchIncludeChannel)
        {
//            fill(*tma,0.0);
//            DistributeOverExtendedChannel(ChannelDetFlow,tma);
//            fill(*tmb,0.0);
//            DistributeOverExtendedChannel(ChannelDep,tmb);
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_CHL
            {
                DETFlowCum->Drc += ChannelDetFlow->Drc;
                DEPCum->Drc += ChannelDep->Drc;
                TotalChanDetMap->Drc += ChannelDetFlow->Drc;
                TotalChanDepMap->Drc += ChannelDep->Drc;
            }}
        }


        // with all det and dep calculate the soil loss, excl channel
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            TotalSoillossMap->Drc = DETSplashCum->Drc + DETFlowCum->Drc + DEPCum->Drc;
            TotalDepMap->Drc = std::min(0.0, TotalSoillossMap->Drc); // for table damage output per landunit
            TotalDetMap->Drc = std::max(0.0, TotalSoillossMap->Drc);
        }}

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            double sedall = Sed->Drc + (SwitchUse2Phase ? BLFlood->Drc : 0.0) + SSFlood->Drc +  (SwitchIncludeChannel ? ChannelSed->Drc : 0.0);
            double waterall = WaterVolall->Drc + (SwitchIncludeChannel ? ChannelWaterVol->Drc : 0.0);
            TotalConc->Drc = MaxConcentration(waterall ,&sedall, NULL);
            // for output

            // set to zero for next loop
            DepFlood->Drc = 0;
            BLDetFlood->Drc = 0;
            SSDetFlood->Drc = 0;

        }}

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

        double waterin = RainTot + SnowTot + WaterVolSoilTot + WHinitVolTot + BaseFlowTot;
        double waterout = ETaTot;
        double waterstore = IntercTot + IntercLitterTot + IntercHouseTot + InfilTot;
        double waterflow = WaterVolTot + ChannelVolTot + StormDrainVolTot + Qtot;
//qDebug() << ChannelVolTot << BaseFlowTot;

        MB = waterin > 0 ? (waterin - waterout - waterstore - waterflow)/waterin *100 : 0;
     //   qDebug() << MB << waterin << waterstore << waterflow;
     //   qDebug() << MB << WaterVolTot << ChannelVolTot << Qtot << floodBoundaryTot;

    }
    //watervoltot includes channel and tile

    // Mass Balance sediment, all in kg

    if (SwitchErosion)
    {
        double detachment = DetTot + ChannelDetTot + FloodDetTot;
        double deposition = DepTot + ChannelDepTot + FloodDepTot;
        double sediment = SedTot + ChannelSedTot + FloodSedTot + SoilLossTot;
        //already in soiloss: + floodBoundarySedTot;

      //  qDebug() << "S" << DetTot<< ChannelDetTot << FloodDetTot;
      //  qDebug() << DepTot << ChannelDepTot << FloodDepTot;
      //  qDebug() << SedTot << ChannelSedTot << FloodSedTot << SoilLossTot;

/*
        if( SwitchKinematic2D == K2D_METHOD_KINDYN)
        {
            // distribute sed errors over dep or det

            double dsed = detachment + deposition  - sediment;
            double count = 0;
            //subtract old totals
            DepTot -= mapTotal(*DEP);
            DetFlowTot -= mapTotal(*DETFlow);
            FloodDetTot -= mapTotal(*SSDetFlood);
            FloodDepTot -= mapTotal(*DepFlood);

            // more deposition than detachment
            if(dsed < 0){
                FOR_ROW_COL_MV {
                    if (DEP->Drc < 0)
                        count += 1.0;
                    if (DepFlood->Drc < 0)
                        count += 1.0;
                }
                FOR_ROW_COL_MV {
                    if (DEP->Drc < 0)
                        DEP->Drc -= dsed/count;
                    if (DEP->Drc > 0) {
                        DETFlow->Drc += DEP->Drc;
                        DEP->Drc = 0;
                    }
                    if (DepFlood->Drc < 0)
                        DepFlood->Drc -= dsed/count;
                    if (DepFlood->Drc > 0) {
                        SSDetFlood->Drc += DepFlood->Drc;
                        DepFlood->Drc = 0;
                    }
                }

            }

            // more detachment than deposition
            if(dsed > 0){
                count = 0;
                FOR_ROW_COL_MV {
                    if (DETFlow->Drc > 0)
                        count += 1.0;
                    if (SSDetFlood->Drc > 0)
                        count += 1.0;
                }

                FOR_ROW_COL_MV {
                    if (DETFlow->Drc > 0)
                        DETFlow->Drc -= dsed/count;
                    if (DETFlow->Drc < 0) {
                        DEP->Drc -= DETFlow->Drc;
                        DETFlow->Drc = 0;
                    }
                    if (SSDetFlood->Drc > 0)
                        SSDetFlood->Drc -= dsed/count;
                    if (SSDetFlood->Drc < 0) {
                        DepFlood->Drc -= SSDetFlood->Drc;
                        SSDetFlood->Drc = 0;
                    }
                }
            }

            // add new totals
            DetFlowTot += mapTotal(*DETFlow);
            DepTot += mapTotal(*DEP);
            FloodDetTot += mapTotal(*SSDetFlood);
            FloodDepTot += mapTotal(*DepFlood);
            DetTot = DetFlowTot + DetSplashTot;
            detachment = DetTot + ChannelDetTot + FloodDetTot;
            deposition = DepTot + ChannelDepTot + FloodDepTot;
        }
*/
        MBs = detachment > 0 ? (detachment + deposition  - sediment)/detachment*100 : 0;
    }

    if (SwitchPesticide)
    {
        MBp = (PestMassApplied-PestLossTotOutlet-PestRunoffSpatial-PestDisMixing-PestSorMixing-PestInfilt-PestStorage)*100/PestMassApplied;
        //MBpex = (PestMassApplied-PestLossTotOutletex-PestRunoffSpatialex-PestDisMixingex-PestSorMixingex-PestInfiltex)*100/PestMassApplied;
        //(PestMassApplied-PestLossTotOutlet-PestRunoffSpatial-PestDisMixing-PestSorMixing-PestInfilt-PestStorage)*100/PestMassApplied
        debug(QString("mbp: %1").arg(MBp));
    }
}
//---------------------------------------------------------------------------
