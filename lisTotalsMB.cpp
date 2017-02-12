
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
       WHmax->Drc = std::max(WHmax->Drc, WH->Drc);
    }
    // find max water height

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
    IntercTot = mapTotal(*Interc) + mapTotal(*LInterc) ;
    // added litter interception
    IntercTotmm = IntercTot*catchmentAreaFlatMM;
    // interception in mm and m3

    //houses
    IntercHouseTot = mapTotal(*IntercHouse);
    IntercHouseTotmm = IntercHouseTot*catchmentAreaFlatMM;
    // interception in mm and m3

    //=== infiltration ===//
    InfilTot += mapTotal(*InfilVol) + mapTotal(*InfilVolKinWave) + mapTotal(*InfilVolFlood); //m3

    InfilKWTot += mapTotal(*InfilVolKinWave); // not really used, available for output when needed
    InfilTotmm = std::max(0.0 ,(InfilTot)*catchmentAreaFlatMM);
    // infiltration mm and m3

    // flood infil
    // used for reporting only
    FOR_ROW_COL_MV
    {
        InfilVolCum->Drc += InfilVol->Drc + InfilVolKinWave->Drc + InfilVolFlood->Drc;
        InfilmmCum->Drc = std::max(0.0, InfilVolCum->Drc*1000.0/(_dx*_dx));
    }

    //=== surf store ===//
    calcMapValue(*tm, *WHstore, 1000, MUL); //mm
    SurfStoremm = mapAverage(*tm);
    // surface storage CHECK THIS
    // does not go to MB, is already in tot water vol

    WaterVolTot = mapTotal(*WaterVolall);//m3
    WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm

    WaterVolRunoffmm = 0;
    FOR_ROW_COL_MV
    {
        WaterVolRunoffmm += WHrunoff->Drc * ChannelAdj->Drc * DX->Drc;
    }
    WaterVolRunoffmm = WaterVolRunoffmm*catchmentAreaFlatMM;
    // water on the surface in runoff in m3 and mm
    //NOTE: surface storage is already in here so does not need to be accounted for in MB

    // runoff fraction per cell calc as in-out/rainfall, indication of sinks and sources of runoff
    // exclude channel cells
    FOR_ROW_COL_MV
        if(ChannelWidthUpDX->Drc == 0)
    {
        runoffTotalCell->Drc += Qn->Drc * _dt; // in M3 !!!!
    }
    upstream(LDD, runoffTotalCell, tm);

    FOR_ROW_COL_MV
    {
        runoffFractionCell->Drc = RainCumFlat->Drc > 0 ? (runoffTotalCell->Drc-tm->Drc)/(RainCumFlat->Drc*_dx*_dx) : 0;
    }
    //TODO runoff fraction not a very clear map

    //=== all discharges ===//

    QtotT = 0;
    // sum all outflow in m3 for this timestep, Qtot is for all timesteps!

    // Add overland flow
    if(SwitchKinematic2D == 1)
    {
        FOR_ROW_COL_MV
        {
            if (LDD->Drc == 5) {
                QtotT += Qn->Drc*_dt;
            }
        }
    }
    else
    {
        QtotT += K2DQOut; // NOTE: K2DQOut already based on if(K2DOutlets->Drc == 1)
    }

    // sum outflow m3 for all timesteps for all outlets, in m3
    // needed for mass balance

    //QtotOutlet += Qn->DrcOutlet*_dt;  obsolete, now Qtot-
    // for file totals output, total main outlet in m3

    // add channel outflow
    if (SwitchIncludeChannel)
    {
        WaterVolTot += mapTotal(*ChannelWaterVol); //m3
        // add channel vol to total
        WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
        // recalc in mm for screen output
        FOR_ROW_COL_MV_CH
        {
            if (LDDChannel->Drc == 5)
            {
                QtotT += ChannelQn->Drc*_dt; //m3
            }
            ChannelQntot->Drc += ChannelQn->Drc*_dt;  //cumulative m3 spatial for .map output

        }
        // add channel outflow (in m3) to total for all pits
        // recalc in mm for screen output

        //QtotOutlet += ChannelQn->DrcOutlet * _dt;  obsolete, now Qtot-
        // sum: add channel outflow (in m3) to total for main outlet -> for totals file

        if (SwitchChannelFlood)
        {
            floodVolTot = mapTotal(*FloodWaterVol);
            floodTotmm = floodVolTot * catchmentAreaFlatMM; // to mm
        }
        if (runstep == 1)
            floodVolTotInit = floodVolTot;
        // save initial flood level for mass balance if start with flood
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
        WaterVolTot += mapTotal(*TileWaterVol); //m3
        // add tile vol to total
        WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
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

    if (SwitchBuffers)
    {
        BufferVolTot = mapTotal(*BufferVol); // in m3
        if (SwitchIncludeChannel)
            BufferVolTot += mapTotal(*ChannelBufferVol);
        //sum up all volume remaining in all buffers (so the non-water!)
        BufferVolin = BufferVolTotInit - BufferVolTot;
        //subtract this from the initial volume to get the total water inflow in the buffers
    }

    // output fluxes for reporting to file and screen in l/s!
    FOR_ROW_COL_MV
    {
        Qoutput->Drc = 1000*(Qn->Drc + ChannelQn->Drc + TileQn->Drc); // in l/s
        //NOTE: for 2D flow Qn = K2DQ, already done
    }

    Qtot += QtotT;
    // Total outflow in m3 for all timesteps
    Qtotmm = (Qtot)*catchmentAreaFlatMM;
    // recalc to mm for screen output

    // flood boundary losses are done separately in MB

    /***** SEDIMENT *****/
    // note DETFLOW, DETSPLASH AND DEP ARE IN KG/CELL

    // DetSplashTot, DetFlowTot and DepTot are for output in file and screen
    // DetTot and DepTot are for MB

    SoilLossTotT = 0;

    if (SwitchErosion)
    {
        DetSplashTot += mapTotal(*DETSplash);
        DetFlowTot += mapTotal(*DETFlow);
        DepTot += mapTotal(*DEP);
        DetTot = DetFlowTot + DetSplashTot; //+= mapTotal(*DETSplash) + mapTotal(*DETFlow);
        SedTot = mapTotal(*Sed);
        // all in kg/cell
        // totals for screen output

        calcMap(*DETSplashCum, *DETSplash, ADD);
        calcMap(*DETFlowCum, *DETFlow, ADD);
        // for screen map output

        if(SwitchKinematic2D == 1)
        {
            FOR_ROW_COL_MV
            {
                if (LDD->Drc == 5)
                {
                    SoilLossTotT += Qsn->Drc * _dt;
                }
            }
        }else
        {
            SoilLossTotT +=K2DQSOut;
        }
        // sum all sed in all pits (in kg), needed for mass balance

// obsolete
//        SoilLossTotOutlet += Qsn->DrcOutlet * _dt;
//        // for screen output, total main outlet sed loss in kg

        if (SwitchIncludeChannel)
        {
            // units here in kg, conversion to ton in report functions
            ChannelDetTot += mapTotal(*ChannelDetFlow);
            ChannelDepTot += mapTotal(*ChannelDep);
            ChannelSedTot = mapTotal(*ChannelBLSed) + mapTotal(*ChannelSSSed);

            fill(*tma,0.0);
            DistributeOverExtendedChannel(ChannelDetFlow,tma);
            calcMap(*DETFlowCum, *tma, ADD);

            FOR_ROW_COL_MV_CH
            {
                if (LDDChannel->Drc == 5)
                {
                    SoilLossTotT += ChannelQsn->Drc * _dt;
                }
            }

//obsolete
//            SoilLossTotOutlet += ChannelQsn->DrcOutlet * _dt;
//            // add channel outflow (in kg) to total for main outlet, for screen output

            if (SwitchChannelFlood)
            {
                FloodDetTot += mapTotal(*BLDetFloodT);  // used for screen output
                FloodDetTot += mapTotal(*SSDetFloodT);

                FloodDepTot += mapTotal(*BLDepFloodT);
                //FloodDepTot += mapTotal(*SSDepFloodT);

                FloodSedTot += mapTotal(*BLDepFloodT);
                FloodSedTot = mapTotal(*BLFlood) + mapTotal(*SSFlood);

                calcMap(*DETFlowCum, *BLDetFloodT, ADD);
                calcMap(*DETFlowCum, *SSDetFloodT, ADD);
            }
        }

        if (SwitchBuffers || SwitchSedtrap)
        {
            BufferSedTot = mapTotal(*BufferSed);
            if (SwitchIncludeChannel)
                BufferSedTot += mapTotal(*ChannelBufferSed);
        }

        // spatial totals for output overland flow all in kg/cell
        // variables are valid for both 1D and 2D flow
        FOR_ROW_COL_MV
        {
            Qsoutput->Drc = Qsn->Drc + ChannelQsn->Drc;  // in kg/s

            TotalDetMap->Drc += DETSplash->Drc + DETFlow->Drc;
            TotalDepMap->Drc += DEP->Drc;
            TotalSoillossMap->Drc = TotalDetMap->Drc + TotalDepMap->Drc;
        }

        if (SwitchIncludeChannel)
        {
            fill(*tma,0.0);
            DistributeOverExtendedChannel(ChannelDetFlow,tma);
            fill(*tmb,0.0);
            DistributeOverExtendedChannel(ChannelDep,tmb);
            // DistributeOverExtendedChannel(ChannelDep,tma); <- typo?
            FOR_ROW_COL_MV
            {

                TotalDetMap->Drc += tma->Drc;

                TotalDepMap->Drc += tmb->Drc;

                TotalChanDetMap->Drc += ChannelDetFlow->Drc;
                TotalChanDepMap->Drc += ChannelDep->Drc;

                if (SwitchChannelFlood)
                {
                    TotalDepMap->Drc += BLDepFloodT->Drc;
                    TotalDetMap->Drc += SSDetFloodT->Drc;
                    TotalDetMap->Drc += BLDetFloodT->Drc;
                }
            }
        }
        FOR_ROW_COL_MV
        {
            double Q = Qoutput->Drc/1000;
            TotalConc->Drc = (Q > MIN_FLUX ? Qsoutput->Drc/Q : 0);
            //WITHOUT FLOOD????
            // conc is okay in principle, Qs is conc (mass/vol)*Q so here this is reversed
        }
    }

    SoilLossTot += SoilLossTotT;

    if(SwitchErosion && SwitchChannelFlood)
    {
        fill(*BLDetFloodT,0.0);
        fill(*BLDepFloodT,0.0);
        fill(*SSDetFloodT,0.0);
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

        // double MBtest=0.0;

        // if (PestLossTotOutlet > 1e-9)
        // MBtest = (Pestdetach-PestCinfilt-PestRunoffSpatial-PestLossTotOutlet)*100/Pestdetach;

        //if (Pestdetach > 1e-9)
        // MBtest = Pestdetach-PestCinfilt-PestCfilmexit-PestLossTotOutlet;
        // qDebug()<< "pestdetach" << Pestdetach << "pestCinfilt"<< PestCinfilt << "pestCfilmexit"<< PestCfilmexit<< "pestlosstotoutlet"<<PestLossTotOutlet;
        // qDebug()<< "MBtest" << MBtest;
        double test=0.0;
        test += mapTotal(*InfilVolKinWave);

        //        PestLossTotOutletex += Qn->DrcOutlet*C_Kexplicit->DrcOutlet*_dt*1000*1000*1000; //µg
        //        PestRunoffSpatialex = PRunoffSpatialex->mapTotal();
        //        PestDisMixingex = PDisMixingex->mapTotal();
        //        PestSorMixingex = PSorMixingex->mapTotal();
        //        PestInfiltex += PInfiltex->mapTotal();

        // flux en µg
        //        double flux1=epsil->DrcOutlet*rhob->DrcOutlet*kr->DrcOutlet*KD->DrcOutlet*CM->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;
        //        double flux2=kr->DrcOutlet*CS->DrcOutlet*rhob->DrcOutlet*epsil->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;
        //        double flux3=pestiinf->DrcOutlet*CM->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;
        //        double flux4=Kfilm->DrcOutlet*CM->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;
        //        double flux5=(Kfilm->DrcOutlet+pestiinf->DrcOutlet)*C->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;
        //        double flux6=(Kfilm->DrcOutlet+RainNet->DrcOutlet/_dt)*C->DrcOutlet*_dx*DX->DrcOutlet*_dt*1000*1000*1000;


        //            QFile fout("massbalancenew.txt");
        //            fout.open(QIODevice::Append | QIODevice::Text);
        //            QTextStream out(&fout);
        //            out.setRealNumberPrecision(3);
        //            out.setFieldWidth(0);
        //            out.setRealNumberNotation(QTextStream::FixedNotation);

        //            out << time/60 << " " << PestMassApplied << " " << PestDisMixing << " " << PestSorMixing << " " << PestLossTotOutlet << " " << PestRunoffSpatial
        //                 << " " << PestInfilt << " " << (PestMassApplied-PestLossTotOutlet-PestRunoffSpatial-PestDisMixing-PestSorMixing-PestInfilt-PestStorage)*100/PestMassApplied << " "
        //                 << RainTot << " " << WaterVolSoilTot << " " << IntercTot << " " << InfilTot << " " << Qtot*1000*1000 << " "
        //                 << MBtest << " " << test << " "<< flux3 << " "<< flux4 << " "<< flux5 << " "<< flux6 <<" "<< pestiinf->DrcOutlet*pow(10.0,9)<< " "<<CM->DrcOutlet*pow(10.0,6)<<" "
        //                 << CS->DrcOutlet*pow(10.0,6)<<" "<< fact->DrcOutlet*1000<< " "<< InfilVol->DrcOutlet*1000*1000<<" "<<Qn->DrcOutlet*pow(10.0,6) << " "<< PDisMixing->DrcOutlet << " "<< poro->DrcOutlet
        //                 << " "<< epsil->DrcOutlet<< " "<< DX->DrcOutlet << " " << switchrunoff << " "<< K1->DrcOutlet << " "<< Q->DrcOutlet*pow(10.0,6)<< " "<< C->DrcOutlet*pow(10.0,10)
        //                 << " "<< WHoutavg->DrcOutlet << " "<< WHoutavgold->DrcOutlet<<" "<< (PestMassApplied-PestLossTotOutletex-PestRunoffSpatialex-PestDisMixingex-PestSorMixingex-PestInfiltex)*100/PestMassApplied
        //                 << " " << InfilVol->DrcOutlet*1000*1000 << " " << InfilVolold->DrcOutlet*1000*1000<< " " << Vup->DrcOutlet << " " << Vup_old->DrcOutlet << " "<< Cold->DrcOutlet*pow(10.0,10);
        //            out << "\n";
        //            out << MBp << "\n";
        //            fout.close();

    }

    SedimentSetMaterialDistribution();

}
//---------------------------------------------------------------------------
void TWorld::MassBalance()
{
    // Mass Balance water, all in m3
    // VJ 110420 added tile volume here, this is the input volume coming from the soil after swatre
    if (RainTot + SnowTot > 0)
    {

        double waterin = RainTot + SnowTot + WaterVolSoilTot + floodVolTotInit + BaseFlow;
        double waterstore = IntercTot + IntercHouseTot + InfilTot + BufferVolin;
        double waterflow = WaterVolTot + floodVolTot + Qtot + floodBoundaryTot;

//        MBeM3 = (RainTot + SnowTot + WaterVolSoilTot + floodVolTotInit + BaseFlow +
//                 - IntercTot - IntercHouseTot - InfilTot - WaterVolTot - floodVolTot - Qtot - BufferVolin - floodBoundaryTot);
//        MB = MBeM3/(RainTot + SnowTot + WaterVolSoilTot + floodVolTotInit)*100;

        MB = waterin > 0 ? (waterin - waterstore - waterflow)/waterin *100 : 0;
    }
    //watervoltot includes channel and tile
//    qDebug() << MB << RainTot << IntercTot << IntercHouseTot << InfilTot << WaterVolTot << floodVolTot << BufferVolin << Qtot<< InfilKWTot;

    // Mass Balance sediment, all in kg

    //VJ 110825 forgot to include channeldettot in denominator in MBs!
    if (SwitchErosion)// && SoilLossTot > 1e-9)
    {
        double detachment = DetTot + ChannelDetTot + FloodDetTot;
        double deposition = DepTot + ChannelDepTot + FloodDepTot;
        double sediment = SedTot + ChannelSedTot + FloodSedTot +BufferSedTot + SoilLossTot;

        MBs = detachment > 0 ? (detachment + deposition  - sediment)/detachment*100 : 0;


//        MBs = (1-(DetTot + ChannelDetTot + FloodDetTot - SedTot - ChannelSedTot - FloodSedTot +
//                  DepTot + ChannelDepTot + FloodDepTot - BufferSedTot)/(SoilLossTot))*100;
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
