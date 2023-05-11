
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
void TWorld::Totals(void)
{
    double rainfall, snowmelt;
    double oldrainpeak, oldsnowpeak;
    double catchmentAreaFlatMM = (1000.0/(_dx*_dx*nrCells));

    /***** WATER *****/

    //=== precipitation ===//
    if (SwitchRainfall)
    {
        RainAvgmm = MapTotal(*Rain)*1000.0/(double)nrValidCells;

        RainTotmm += RainAvgmm;
        // spatial avg area rainfall in mm

        rainfall = RainAvgmm/catchmentAreaFlatMM;
        RainTot += rainfall; // in m3

        oldrainpeak  = Rainpeak;
        Rainpeak = std::max(Rainpeak, rainfall);
        if (oldrainpeak  < Rainpeak)
            RainpeakTime = time;
    }

    if (SwitchSnowmelt)
    {
        SnowAvgmm = MapTotal(*Snowmelt)*1000.0/(double)nrValidCells;

        SnowTotmm += SnowAvgmm;

        snowmelt = SnowAvgmm/catchmentAreaFlatMM;
        SnowTot += snowmelt; // in m3

        oldsnowpeak = Snowpeak;
        Snowpeak = std::max(Snowpeak, snowmelt);
        if (oldsnowpeak < Snowpeak)
            SnowpeakTime = time;
    }

    //=== interception ===//
    IntercTot = MapTotal(*Interc);
    IntercTotmm = IntercTot*catchmentAreaFlatMM;
    // currently in canopy

    if (SwitchIncludeET) {
       // double ETtot = MapTotal(*ETa);
        ETaTot = MapTotal(*ETaCum);
        ETaTotmm = ETaTot * 1000.0/(double)nrValidCells;

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
        InterceptionmmCum->Drc = (IntercETa->Drc + Interc->Drc + IntercHouse->Drc + LInterc->Drc)*1000.0/CellArea->Drc;
        // for screen and file output
    }}


    //==== ETa ==========//
    //   ETaTot = mapTotal(*ETa);
    //  ETaTotmm = ETaTot*catchmentAreaFlatMM;
    // interception in mm and m3

    //=== infiltration ===//
    if(InfilMethod != INFIL_NONE) {
        InfilTot += MapTotal(*InfilVol);// + MapTotal(*InfilVolKinWave);

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
            if (SwitchIncludeChannel && SwitchChannelInfil)
                InfilVolCum->Drc += ChannelInfilVol->Drc;

            InfilmmCum->Drc = std::max(0.0, InfilVolCum->Drc*1000.0/(_dx*_dx));
            PercmmCum->Drc += Perc->Drc*1000.0;
        }}

        theta1tot = MapTotal(*ThetaI1a)/(double)nrCells;
        if (SwitchTwoLayer)
            theta2tot = MapTotal(*ThetaI2a)/(double)nrCells;
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
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        //runoffTotalCell->Drc += (Qn->Drc /*+ Qflood->Drc*/)* _dt * catchmentAreaFlatMM; // in mm !!!!
        runoffTotalCell->Drc = std::max(0.0, RainCumFlat->Drc*1000-InterceptionmmCum->Drc-InfilmmCum->Drc);
    }}

    //=== channel flow ===//
    if (SwitchIncludeChannel)
    {
        ChannelVolTot = MapTotal(*ChannelWaterVol); //m3
        // add channel vol to total
        if (SwitchChannelBaseflow) {
            BaseFlowTot += MapTotal(*Qbase); // total inflow in m3
            if (SwitchChannelBaseflowStationary)
                BaseFlowTot += MapTotal(*BaseFlowInflow)*_dt; // stationary base inflow

            GWlevel = MapTotal(*GWWH)/(double)nrValidCells;
            BaseFlowTotmm = BaseFlowTot*catchmentAreaFlatMM; //mm
            //qDebug() << BaseFlowTotmm;
        }

        ChannelVolTotmm = ChannelVolTot*catchmentAreaFlatMM; //mm
        // recalc in mm for screen output

    }

    //=== all discharges ===//
    Qtot_dt = 0;
    // sum all outflow in m3 for this timestep, Qtot is for all timesteps!

    floodBoundaryTot += BoundaryQ*_dt;
    FloodBoundarymm = floodBoundaryTot*catchmentAreaFlatMM;    
   // qDebug() << floodBoundaryTot;

    // Add outlet overland flow, for all flow methods
    FOR_ROW_COL_LDD5 {
        Qtot_dt += Qn->Drc*_dt;
    }}

    // add channel outflow
    if (SwitchIncludeChannel)
    {
            FOR_ROW_COL_LDDCH5 {
                Qtot_dt += ChannelQn->Drc*_dt; //m3
            }}
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_CHL {
                ChannelQntot->Drc += ChannelQn->Drc*_dt;
                //cumulative m3 spatial for .map output
            }}

        // add channel outflow (in m3) to total for all pits
    }

            //=== storm drain flow ===//
    if (SwitchIncludeTile || SwitchIncludeStormDrains) {
        StormDrainVolTot = MapTotal(*TileWaterVol);
        StormDrainTotmm = StormDrainVolTot*catchmentAreaFlatMM;
    }

    if(SwitchIncludeStormDrains) {
        FOR_ROW_COL_MV_TILE
                if (LDDTile->Drc == 5)
        {
            Qtot_dt += TileQn->Drc * _dt;
            QTiletot += TileQn->Drc * _dt;
        }
    } else {
        if (SwitchIncludeTile)
        {
            WaterVolSoilTileTot = MapTotal(*TileWaterVolSoil);
            // input for mass balance, is the water seeping from the soil, input
            // this is the water before the kin wave
            calc2Maps(*tm, *TileDrainSoil, *TileWidth, MUL); //in m3
            calcMap(*tm, *TileDX, MUL); //in m3
            // tm->calcV(_dx, MUL); //in m3 ??? or DX?
            TileVolTot += MapTotal(*tm); // in m3

            FOR_ROW_COL_MV_TILE
                    if (LDDTile->Drc == 5)
            {
                Qtot_dt += TileQn->Drc * _dt;
                QTiletot += TileQn->Drc * _dt;
            }

            // add tile outflow (in m3) to total for all pits
        }
    }

    // sum of all fluxes ONLY for display on screen
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L
    {
        Qoutput->Drc = (Qn->Drc + Qflood->Drc) * (QUnits == 1 ? 1.0 : 1000);// in m3/s

        if(SwitchIncludeChannel)
            Qoutput->Drc += ChannelQn->Drc * (QUnits == 1 ? 1.0 : 1000);

        Qoutput->Drc = Qoutput->Drc < 1e-6 ? 0.0 : Qoutput->Drc;
    }}
    // Total outflow in m3 for all timesteps
    // does NOT include flood water leaving domain (floodBoundaryTot)
    // which is reported separatedly (because it is a messy flux)!

    Qtot += Qtot_dt;
    // add timestep total to run total in m3
    Qtotmm = Qtot*catchmentAreaFlatMM;
    // recalc to mm for screen output


    //=====***** SEDIMENT *****====//

    // DetSplashTot, DetFlowTot and DepTot are for output in file and screen
    // DetTot and DepTot are for MB

    SoilLossTot_dt = 0;

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
        }

        floodBoundarySedTot += BoundaryQs*_dt; // not used
        SoilLossTot_dt += BoundaryQs*_dt;
        // boundary sediment losses (kg) in cells that are not outlet, if open boundary else 0
        // calc as cells with velocity U and V directed outwards

        // used for mass balance and screen output
        FloodDetTot += (SwitchUse2Phase ? MapTotal(*BLDetFlood) : 0.0) + MapTotal(*SSDetFlood);
        FloodDepTot += MapTotal(*DepFlood);
        FloodSedTot = (SwitchUse2Phase ? MapTotal(*BLFlood) : 0.0) + MapTotal(*SSFlood);

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
            Qsoutput->Drc = Qsn->Drc + (SwitchIncludeChannel ? ChannelQsn->Drc : 0.0);
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


    SedimentSetMaterialDistribution();

}
//---------------------------------------------------------------------------
void TWorld::MassBalance()
{
    // Mass Balance water, all in m3
    // VJ 110420 added tile volume here, this is the input volume coming from the soil after swatre
  //  if (RainTot + SnowTot > 0)
  //  {
        double waterin = RainTot + SnowTot + WaterVolSoilTileTot + WHinitVolTot + BaseFlowTot + BaseFlowInit;
                //qDebug() << RainTot << thetai1tot - thetai1cur << thetai2tot - thetai2cur;
        double waterout = ETaTotVol;
        double waterstore = IntercTot + IntercLitterTot + IntercHouseTot + InfilTot + IntercETaTot;// + (thetai1cur - thetai1tot) + (thetai2cur - thetai2tot);
        double waterflow = WaterVolTot + ChannelVolTot + StormDrainVolTot + Qtot + floodBoundaryTot;
        MB = waterin > 0 ? (waterin - waterout - waterstore - waterflow)/waterin *100 : 0;
     //   qDebug() << MB << BaseFlowTot << ChannelVolTot <<  Qtot;
     //   qDebug() << MB << WaterVolTot << ChannelVolTot << Qtot << floodBoundaryTot;

   // }

    //watervoltot includes channel and tile

    // Mass Balance sediment, all in kg

    if (SwitchErosion)
    {
        double detachment = DetTot + ChannelDetTot + FloodDetTot;
        double deposition = DepTot + ChannelDepTot + FloodDepTot;
        double sediment = SedTot + ChannelSedTot + FloodSedTot + SoilLossTot;
        //already in SoilLossTot: floodBoundarySedTot;

      //  qDebug() << "S" << DetTot<< ChannelDetTot << FloodDetTot;
      //  qDebug() << DepTot << ChannelDepTot << FloodDepTot;
      //  qDebug() << SedTot << ChannelSedTot << FloodSedTot << SoilLossTot;

        MBs = detachment > 0 ? (detachment + deposition  - sediment)/detachment*100 : 0;
    }

    if (SwitchPest) {
        MassPest(PMtotI, PMerr, PMtot, PMserr, PMwerr);
    }
}
//---------------------------------------------------------------------------
