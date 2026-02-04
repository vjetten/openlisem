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
void TWorld::TotalsHydro(void)
{
    double rainfall, snowmelt;
    double oldrainpeak, oldsnowpeak;
    double catchmentAreaFlatMM = (1000.0/(_dx*_dx*nrCells));

    /***** WATER *****/

    //=== precipitation ===//
    if (SwitchRainfall)
    {
        double ptot = MapTotal(*Rain);
        RainAvgmm = ptot*1000.0/nrCells;
        RainTotmm += RainAvgmm;
        // spatial avg area rainfall in mm

        rainfall = RainAvgmm/catchmentAreaFlatMM;
        RainTot += ptot*_dx*_dx; // in m3

        oldrainpeak = Rainpeak;
        Rainpeak = qMax(Rainpeak, rainfall);
        if (oldrainpeak  < Rainpeak)
            RainpeakTime = time;
    }

    if (SwitchSnowmelt)
    {
        SnowAvgmm = MapTotal(*Snowmelt)*1000.0/nrCells;

        SnowTotmm += SnowAvgmm;

        snowmelt = SnowAvgmm/catchmentAreaFlatMM;
        SnowTot += snowmelt; // in m3

        oldsnowpeak = Snowpeak;
        Snowpeak = qMax(Snowpeak, snowmelt);
        if (oldsnowpeak < Snowpeak)
            SnowpeakTime = time;
    }

    //=== interception ===//
    IntercTot = MapTotal(*Interc);
    IntercTotmm = IntercTot*catchmentAreaFlatMM;
    // currently in canopy

    if (SwitchIncludeET) {
        ETaTot = MapTotal(*ETaCum);
        ETaTotmm = ETaTot * 1000.0/nrCells;

        ETaTotVol = (ETaTot-SoilETMBcorrection)*_dx*_dx; //m3
        // correct for soil water because that is not in the mass balance

        IntercETaTot = MapTotal(*IntercETa);
        IntercETaTotmm = IntercETaTot*catchmentAreaFlatMM;
        // cumulative evaporated from canopy
    }


    // interception in mm and m3
    //Litter
    if (SwitchLitter) {
        IntercLitterTot = MapTotal(*LInterc); // in m
        IntercLitterTotmm = IntercLitterTot*catchmentAreaFlatMM; // *1000/total cellarea
    }

    if (SwitchHouses) {
        IntercHouseTot = MapTotal(*IntercHouse);
        IntercHouseTotmm = IntercHouseTot*catchmentAreaFlatMM;
        // interception in mm and m3
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        InterceptionmmCum->Drc = Interc->Drc;
        if (SwitchIncludeET)
            InterceptionmmCum->Drc += IntercETa->Drc;
        if (SwitchHouses)
            InterceptionmmCum->Drc += IntercHouse->Drc;
        if (SwitchLitter)
            InterceptionmmCum->Drc += LInterc->Drc;
        InterceptionmmCum->Drc *= 1000.0/CellArea->Drc;
        // for screen and file output
    }}

    //==== ETa ==========//
    //   ETaTot = mapTotal(*ETa);
    //  ETaTotmm = ETaTot*catchmentAreaFlatMM;
    // interception in mm and m3

    //=== infiltration ===//
    if(SwitchInfiltration) {
        InfilTot += MapTotal(*InfilVol);   //obsolete + MapTotal(*InfilVolKinWave);
        InfilTotmm = qMax(0.0 ,(InfilTot)*catchmentAreaFlatMM);
        // used in reporting
        // infiltration mm and m3

        // flood infil
        // used for reporting only
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            InfilVolCum->Drc += InfilVol->Drc;
            InfilmmCum->Drc = InfilVolCum->Drc*1000.0/(_dx*_dx);
            PercmmCum->Drc += Perc->Drc*1000.0;
        }}

        theta1tot = MapTotal(*ThetaI1a)/nrCells;
        if (SwitchTwoLayer)
            theta2tot = MapTotal(*ThetaI2a)/nrCells;

    }

    //=== surf store ===//

    double SStot = 0;
    #pragma omp parallel for reduction(+:SStot) num_threads(userCores)
    FOR_ROW_COL_MV_L {
        SStot += MicroStoreVol->Drc;
    }}
    SurfStoremm = SStot * catchmentAreaFlatMM;

    // does not go to MB, is already in tot water vol

    //TODO: check init WH
    //ONLY ONCE?
//    if (SwitchFloodInitial) {
//        WHinitVolTot = 0;
//        #pragma omp parallel for reduction(+:WHinitVolTot) num_threads(userCores)
//        FOR_ROW_COL_MV_L {
//            WHinitVolTot = hmxInit->Drc * DX->Drc * ChannelAdj->Drc;
//        }}
//    }


    SoilMoistTot += SoilMoistDiff;// MapTotal(*SoilMB);
    SoilMoistTotmm = SoilMoistTot * catchmentAreaFlatMM;

}


void TWorld::TotalsFlow(void)
{
    double catchmentAreaFlatMM = (1000.0/(_dx*_dx*nrCells));

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
    // water on the surface in runoff in mm, used in screen output

    // runoff fraction per cell calc as in-out/rainfall, indication of sinks and sources of runoff
    // exclude channel cells
    // #pragma omp parallel for num_threads(userCores)
    // FOR_ROW_COL_MV_L {
    //     //runoffTotalCell->Drc += (Qn->Drc)* _dt * catchmentAreaFlatMM; // in mm !!!!
    //     runoffTotalCell->Drc = qMax(0.0, RainCumFlat->Drc*1000-InterceptionmmCum->Drc-InfilmmCum->Drc);
    // }}


    RetentionVolTot = 0;
    if (SwitchGridRetention) {
        FOR_ROW_COL_MV_L {
            if (GridRetentionAct->Drc > 0)
                RetentionVolTot += GridRetentionAct->Drc;
        }}
        RetentionVolTotmm = RetentionVolTot*catchmentAreaFlatMM;
        qDebug() << "OF" << RetentionVolTot << RetentionVolTotPot;
    }

    //=== channel flow ===//
    if (SwitchIncludeChannel) {
        ChannelVolTot = MapTotal(*ChannelWaterVol); //m3
        ChannelVolTotmm = ChannelVolTot*catchmentAreaFlatMM; //mm

        if (SwitchChannelInfil) {
            InfilTot += MapTotal(*ChannelInfilVol); //m3
            InfilTotmm = qMax(0.0 ,(InfilTot)*catchmentAreaFlatMM);

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                InfilVolCum->Drc += ChannelInfilVol->Drc;
                InfilmmCum->Drc = InfilVolCum->Drc*1000.0/(_dx*_dx);
            }}
        }

        if (SwitchGWflow) {
            BaseFlowTot += MapTotal(*Qbase); // total inflow in m3

            GWlevel = MapTotal(*GWWH);
            GWleveltot = GWlevel*catchmentAreaFlatMM;
            GWlevel /= (double)nrValidCells; // avg GW level
            // BaseFlowTotmm = BaseFlowTot*catchmentAreaFlatMM; //mm
            //qDebug() << BaseFlowTotmm;
        }

        if (SwitchChannelBaseflowStationary)
            BaseFlowTot += MapTotal(*BaseFlowInflow)*_dt;
            // stationary base inflow every timestep, counts as input in mass balance

        // recalc in mm for screen output
        // NOT USED
        // if (SwitchChannelWFinflow) {
        //     QSideVolTot += MapTotal(*ChannelQSide);
        //         //use baseflow for channel side inflow so that it is reported
        //         double tot = 0;
        //         FOR_ROW_COL_MV_CHL {
        //             tot += ChannelQSide->Drc; // total inflow in m3
        //         }}

        //    BaseFlowTot += tot;
        // }

        BaseFlowTotmm = BaseFlowTot*catchmentAreaFlatMM; //mm
        BaseFlowInitmm = BaseFlowInit*catchmentAreaFlatMM;

        if (SwitchGridRetention) {
            ChanRetentionVolTot = 0;
            FOR_ROW_COL_MV_CHL {
                if (ChanRetentionAct->Drc > 0)
                    ChanRetentionVolTot += ChanRetentionAct->Drc;
            }}
            RetentionVolTotmm += ChanRetentionVolTot*catchmentAreaFlatMM;
            qDebug() << "chan" << ChanRetentionVolTot << ChanRetentionVolTotPot;
        }

    }


    //=== all discharges ===//
    Qtot_dt = 0;
    // sum all outflow in m3 for this timestep, Qtot is for all timesteps!

    if (FlowBoundaryType > 0) {
        QBoundaryTot += QBoundary*_dt;
        Qboundtotmm = QBoundaryTot*catchmentAreaFlatMM;
        //Qtot_dt += QBoundary*_dt;
        // do not add boundary to total, report separately
    }

    // Add outlet overland flow, for all flow methods
    FOR_ROW_COL_LDD5 {
        Qtot_dt += Qn->Drc*_dt;
    }}

    //=== channel outflow ===//
    if (SwitchIncludeChannel)
    {
        FOR_ROW_COL_LDDCH5 {
            Qtot_dt += ChannelQn->Drc*_dt; //m3
        }}

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelQntot->Drc += ChannelQn->Drc*_dt;
            //cumulative m3 spatial for .map output
            //QuserInTot += ChannelQSide->Drc;
        }}
        // add channel outflow (in m3) to total for all pits
    }


    //=== storm drain flow
    if(SwitchIncludeStormDrains || SwitchIncludeTile) {
            // sum the tile outlets
        //    QTiletot = 0;
            QTile = 0;
            FOR_ROW_COL_MV_TILEL {
                if (LDDTile->Drc == 5) {
                  QTiletot += TileQn->Drc * _dt;
                  QTile += TileQn->Drc;
                }
            }}

        //urban volume in drains
        if (SwitchIncludeStormDrains) {
            StormDrainVolTot = MapTotal(*TileWaterVol);
        }

        // agriculture volume in tiles
        if (SwitchIncludeTile) {
           StormDrainVolTot += MapTotal(*TileWaterVolSoil);
        }
        // output
        StormDrainTotmm = StormDrainVolTot*catchmentAreaFlatMM;
        //qDebug() << StormDrainVolTot << QTiletot << tott << tilein;
    }

    // sum of all fluxes ONLY for display on screen
    double factor =  (QUnits == 1 ? 1.0 : 1000);
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        Qm3total->Drc += Qn->Drc * _dt; // ONLY OVERLAND FLOW
        Qm3max->Drc = qMax(Qm3max->Drc, Qn->Drc);
        Qoutput->Drc = Qn->Drc * factor;// in m3/s

        FHI->Drc = Qn->Drc*(V->Drc + 0.5);

        if(SwitchIncludeChannel) {
            Qoutput->Drc += ChannelQn->Drc * factor;
            //Qm3total->Drc += ChannelQn->Drc * _dt;
//            Qm3max->Drc = qMax(Qm3max->Drc, ChannelQn->Drc);
        }

        //if(FlowBoundaryType > 0) {
            //Qoutput->Drc += QBoundFlow->Drc * factor;
            //Qm3total->Drc += QBoundFlow->Drc * _dt;
            //Qm3max->Drc = qMax(Qm3max->Drc, QBoundFlow->Drc);
        //}

        Qoutput->Drc = Qoutput->Drc < 1e-10 ? 0.0 : Qoutput->Drc;
    }}
    // Total outflow in m3 for all timesteps
    // does NOT include flood water leaving domain (floodBoundaryTot)
    // which is reported separatedly (because it is a messy flux)!

   // report(*Qototal,"qtotm3.map");
    Qtot += Qtot_dt;
    // add timestep total to run total in m3
    Qtotmm = Qtot*catchmentAreaFlatMM;
    // recalc to mm for screen output
    PeakFlowTotmm = Qtotmm - BaseFlowTotmm;
}

void TWorld::TotalsSediment(void)
{
    //double catchmentAreaFlatMM = (1000.0/(_dx*_dx*nrCells));

    //=====***** SEDIMENT *****====//

    // DetSplashTot, DetFlowTot and DepTot are for output in file and screen
    // DetTot and DepTot are for MB

    SoilLossTot_dt = 0;
    // everything that flows out with channelqs, qs and Qsboundary * _dt

    if (SwitchErosion)
    {
        SedTot = 0;
        //#pragma omp parallel for reduction(+:DetSplashTot,DetFlowTot,DepTot,SedTot) num_threads(userCores)
        FOR_ROW_COL_MV_L {
             // Dep and Detflow are zero if 2Ddyn
            DetSplashTot += DETSplash->Drc;
            DetFlowTot += DETFlow->Drc;
            DepTot += DEP->Drc;
            SedTot += Sed->Drc;
        }}
        // all in kg/cell

        DetTot = DetFlowTot + DetSplashTot;

        // these maps combine kin wave OF and all 2D flow and channelflow
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
           // #pragma omp parallel for reduction(+:SoilLossTotT) num_threads(userCores)
            FOR_ROW_COL_LDD5 {
                SoilLossTot_dt += Qsn->Drc * _dt;
            }}

        }

        if (SwitchIncludeChannel)
        {
          //  #pragma omp parallel for reduction(+:SoilLossTotT) num_threads(userCores)
            FOR_ROW_COL_LDDCH5 {
                SoilLossTot_dt += ChannelQsn->Drc * _dt;
            }}

            ChannelDetTot += MapTotal(*ChannelDetFlow);
            ChannelDepTot += MapTotal(*ChannelDep);
            ChannelSedTot = (SwitchUse2Phase ? MapTotal(*ChannelBLSed) : 0.0) + MapTotal(*ChannelSSSed);

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_CHL {
                double sed = (SwitchUse2Phase ? ChannelBLSed->Drc : 0) + ChannelSSSed->Drc;
                //total concentration
                ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, sed);//ChannelSed->Drc);
            }}
        }

        if (FlowBoundaryType > 0) {
            floodBoundarySedTot += QsBoundary*_dt;        //reported
            SoilLossTot_dt += QsBoundary*_dt;
        }
        // boundary sediment losses (kg) in cells that are not outlet, if open boundary else 0
        // calc as cells with velocity U and V directed outwards

        // used for mass balance and screen output
        if (SwitchKinematic2D > K2D_METHOD_KIN) {
            FloodDetTot += (SwitchUse2Phase ? MapTotal(*BLDetFlood) : 0.0) + MapTotal(*SSDetFlood);
            FloodDepTot += MapTotal(*DepFlood);
            FloodSedTot = (SwitchUse2Phase ? MapTotal(*BLFlood) : 0.0) + MapTotal(*SSFlood);
        }

        if (SwitchUse2Phase) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                DETFlowCum->Drc += BLDetFlood->Drc;
            }}
        }
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            DETFlowCum->Drc += SSDetFlood->Drc;
            DEPCum->Drc += DepFlood->Drc;
        }}
        // SPATIAL totals for output overland flow all in kg/cell
        // variables are valid for both 1D and 2D flow dyn and diff

        FOR_ROW_COL_MV_L {
            Qsoutput->Drc = Qsn->Drc + (SwitchIncludeChannel ? ChannelQsn->Drc : 0.0) + QsBoundary/_dt;
            // for reporting sed discharge screen
            // in kg/s, sum of overland flow and channel flow
        }}

        // for reporting
        if (SwitchIncludeChannel)
        {
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
        // kg/cell
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            TotalSoillossMap->Drc = DETSplashCum->Drc + DETFlowCum->Drc + DEPCum->Drc;
            TotalSoillossMap->Drc = fabs(TotalSoillossMap->Drc) < 1e-3 ? 0.0 : TotalSoillossMap->Drc;
            // 0.001 kg/cellarea = 1/cellarea g/m2
        }}

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L
        {
            double sedall = Sed->Drc + (SwitchUse2Phase ? BLFlood->Drc : 0.0) + SSFlood->Drc;
            if (SwitchIncludeChannel) {
                sedall += ChannelSSSed->Drc;
                if (SwitchUse2Phase)
                    sedall += ChannelBLSed->Drc;
            }
            // Sed is kin wave, SSFlood and BL flood is dyn wave, chnnelsed = kin wave?
            double waterall = WaterVolall->Drc + (SwitchIncludeChannel ? ChannelWaterVol->Drc : 0.0);
            TotalConc->Drc = MaxConcentration(waterall ,sedall);
            // for output

            // set to zero for next loop
            DepFlood->Drc = 0;
            BLDetFlood->Drc = 0;

            SSDetFlood->Drc = 0;

        }}

        SoilLossTot += SoilLossTot_dt;
        // total sediment outflow from outlets and domain boundaries
        // this is the value reported in the screen for total soil loss (/1000 for ton)
        // so this is the total loos through the outlets and boundaries

    }
}
//---------------------------------------------------------------------------
void TWorld::MassBalance()
{
// in mm as displayed on screen
//    double in = RainTotmm + BaseFlowTotmm + BaseFlowInitmm;// + SoilMoistTotmm;
//    double store = IntercTotmm + IntercHouseTotmm + IntercLitterTotmm + InfilTotmm +
//           ETaTotmm +StormDrainTotmm +SurfStoremm + WaterVolRunoffmm + ChannelVolTotmm + floodVolTotmm;
//    double out = Qtotmm + Qboundtotmm;
//    MB = in > 0 ? (in - out - store)/in *100 : 0;

    // Mass Balance water, all in m3
    double waterin = RainTot + WHinitVolTot + BaseFlowTot + BaseFlowInit + QuserInTot;// - QSideVolTot;
    double waterstore = IntercTot + IntercLitterTot + IntercHouseTot + InfilTot  + WaterVolTot + ChannelVolTot + StormDrainVolTot + RetentionVolTot;
    double waterout = Qtot + IntercETaTot + QTiletot + QBoundaryTot;
    // floodBoundaryTot is already in Qtot
    MB = waterin > 0 ? (waterin - waterout - waterstore)/waterin*100  : 0;

    // qDebug() << RainTot << IntercTot << IntercHouseTot << InfilTot  << WaterVolTot << ChannelVolTot <<  Qtot ;

    Fill(*MBm, 0);

    if (SwitchCorrectMB_WH && fabs(MB) > 1e-6) {
        //qDebug() << "o " << MB;
        // correct WH
        FOR_ROW_COL_MV_L {
            tma->Drc = 0;
            if (WHrunoff->Drc > 0 || hmxrunoff->Drc > 0)
                tma->Drc = 1;
        }}
        double tot = MapTotal(*tma);
        double dV = (waterin - waterout - waterstore)/tot;
        waterstore -= WaterVolTot;
        FOR_ROW_COL_MV_L {
            double dH = dV/(CHAdjDX->Drc); // avg error in m on wet cells
            if (FloodDomain->Drc == 0 && WHrunoff->Drc > 0) {
                WHrunoff->Drc = qMax(0.0,WHrunoff->Drc + dH);
                WH->Drc = WHrunoff->Drc + WHstore->Drc;
                hmxWH->Drc = WH->Drc + hmx->Drc;
                WaterVolall->Drc = WH->Drc*CHAdjDX->Drc;//WHrunoff->Drc*CHAdjDX->Drc + MicroStoreVol->Drc;
            }
            if (FloodDomain->Drc > 0 && hmxrunoff->Drc > 0) {
                hmxrunoff->Drc = qMax(0.0,hmxrunoff->Drc + dH);
                hmx->Drc = hmxrunoff->Drc + WHstore->Drc;
                hmxWH->Drc = WH->Drc + hmx->Drc;
                WaterVolall->Drc = hmxWH->Drc*CHAdjDX->Drc;//WHrunoff->Drc*CHAdjDX->Drc + MicroStoreVol->Drc;
            }
        }}
        WaterVolTot = MapTotal(*WaterVolall);
        waterstore += WaterVolTot;
        //IntercTot + IntercLitterTot + IntercHouseTot + InfilTot  + WaterVolTot + ChannelVolTot + StormDrainVolTot + RetentionVolTot;

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            MBm->Drc = (waterin - waterout - waterstore)/nrCells/CellArea->Drc;
            //average error in m water height, added to WH in next timestep in infiltration
        }}

        MB = waterin > 0 ? (waterin - waterout - waterstore)/waterin*100  : 0;
        //qDebug() << "n " << MB;
    }
   // qDebug() << MB;
    // Mass Balance sediment, all in kg
    if (SwitchErosion)
    {
        double detachment = DetTot + ChannelDetTot + FloodDetTot;
        double deposition = DepTot + ChannelDepTot + FloodDepTot;
        double sediment = SedTot + ChannelSedTot + FloodSedTot + SoilLossTot;// + floodBoundarySedTot; //<= is already in total

      //  qDebug() << "S" << DetTot<< ChannelDetTot << FloodDetTot;
      //  qDebug() << DepTot << ChannelDepTot << FloodDepTot;
      //  qDebug() << SedTot << ChannelSedTot << FloodSedTot << SoilLossTot;

        MBs = detachment > 0 ? (detachment + deposition  - sediment)/detachment*100 : 0;
    }

}
//---------------------------------------------------------------------------
