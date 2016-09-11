
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
        SnowAvgmm += mapAverage(*Snowmelt)*1000;
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
    InfilTot += mapTotal(*InfilVol) + mapTotal(*InfilVolKinWave); //m3

    InfilKWTot += mapTotal(*InfilVolKinWave); // not really used, available for output when needed
    InfilTotmm = std::max(0.0 ,(InfilTot)*catchmentAreaFlatMM);
    // infiltration mm and m3

    // flood infil
    // used for reporting only
    FOR_ROW_COL_MV
    {
        InfilVolCum->Drc += InfilVol->Drc + InfilVolKinWave->Drc;
        InfilmmCum->Drc = std::max(0.0, InfilVolCum->Drc*1000.0/(_dx*_dx));
    }

    //=== surf store ===//
    calcMapValue(*tm, *WHstore, 1000, MUL); //mm
    SurfStoremm = mapAverage(*tm);
    // surface storage CHECK THIS
    // does not go to MB, is already in tot water vol

    WaterVolTot = mapTotal(*WaterVolall);//m3
    WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
    WaterVolRunoffmm =0;
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

    //=== all discharges ===//
    QtotT = 0;
    // sum outflow m3 for all timesteps for the outlet

    QtotT += UF2D_foutflow;

    // sum outflow m3 for all timesteps for all outlets, in m3
    // needed for mass balance

    QtotOutlet += UF1D_q->DrcOutlet*_dt;
    // for screen output, total main outlet in m3

    if (SwitchIncludeChannel)
    {
        WaterVolTot += mapTotal(*UF1D_f); //m3
        // add channel vol to total
        WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
        // recalc in mm for screen output

        QtotT += UF1D_foutflow;

        //FOR_ROW_COL_MV_CH
        //{
        //    ChannelQntot->Drc += UF1D_q->Drc*_dt;  //m3 spatial for output
        //}
        // add channel outflow (in m3) to total for all pits
        // recalc in mm for screen output

        QtotOutlet += UF1D_q->DrcOutlet * _dt;
        // sum: add channel outflow (in m3) to total for main outlet
    }

    // output fluxes for reporting to file and screen in l/s!
    FOR_ROW_COL_MV
    {
        Qoutput->Drc = 1000.0*(UF1D_q->Drc + UF2D_q->Drc); // in l/s
        if (Qoutput->Drc < 0.0001)
            Qoutput->Drc = 0.0001;
        // added minimum here to avoid strange maps
    }

    Qtot += QtotT;

    Qtotmm = (Qtot)*catchmentAreaFlatMM;
    // recalc to mm for screen output*/

    /***** SEDIMENT *****/
    // note DETFLOW, DETSPLASH AND DEP ARE IN KG/CELL

    // DetSplashTot, DetFlowTot and DepTot are for output in file and screen
    // DetTot and DepTot are for MB

    SoilLossTotT = 0;

    if (SwitchErosion)
    {
        DetSplashTot += mapTotal(*DETSplash);
        DetFlowTot += mapTotal(*UF2D_Det);
        DepTot += mapTotal(*UF2D_Dep);
        DetTot = DetFlowTot + DetSplashTot; //+= mapTotal(*DETSplash) + mapTotal(*DETFlow);
        SedTot = mapTotal(*UF2D_blm)+ mapTotal(*UF2D_ssm);
        // all in kg/cell
        // totals for screen output

        calcMap(*DETSplashCum, *DETSplash, ADD);
        calcMap(*DETFlowCum, *UF2D_Det, ADD);
        // for screen map output

        SoilLossTotT +=UF2D_fsoutflow;
        // sum all sed in all pits (in kg), needed for mass balance

        SoilLossTotOutlet += UF2D_qs->DrcOutlet * _dt;
        // for screen output, total main outlet sed loss in kg

        if (SwitchIncludeChannel)
        {
            // units here in kg, conversion to ton in report functions
            ChannelDetTot += mapTotal(*UF1D_Det);
            ChannelDepTot += mapTotal(*UF1D_Dep);
            ChannelSedTot = mapTotal(*UF1D_blm) + mapTotal(*UF1D_ssm);

            SoilLossTotT += UF1D_fsoutflow;

            SoilLossTotOutlet += UF1D_qs->DrcOutlet * _dt;
            // add channel outflow (in kg) to total for main outlet, for screen output
        }

        // spatial totals for output all in kg/cell
        FOR_ROW_COL_MV
        {
            Qsoutput->Drc = UF1D_qs->Drc + UF2D_qs->Drc;  // sum channel and OF sed output in kg/s

            TotalDetMap->Drc += DETSplash->Drc + UF2D_Det->Drc;
            TotalDepMap->Drc += UF2D_Dep->Drc;

            if (SwitchIncludeChannel)
            {
                TotalDetMap->Drc += UF1D_Det->Drc;
                TotalDepMap->Drc += UF1D_Dep->Drc;

                TotalChanDetMap->Drc += UF1D_Det->Drc;
                TotalChanDepMap->Drc += UF1D_Dep->Drc;

                TotalDepMap->Drc += UF1D_Dep->Drc;
                TotalDetMap->Drc += UF1D_Det->Drc;
            }
            TotalSoillossMap->Drc = TotalDetMap->Drc + TotalDepMap->Drc;
        }
        FOR_ROW_COL_MV
        {
            double Q = Qoutput->Drc/1000;
            TotalConc->Drc = (Q > 1e-6 ? Qsoutput->Drc/Q : 0);
        }
    }

    SoilLossTot += SoilLossTotT;

    if(SwitchErosion)
    {
        fill(*UF1D_Det,0.0);
        fill(*UF2D_Det,0.0);
        fill(*UF1D_Dep,0.0);
        fill(*UF2D_Dep,0.0);
    }
    if(SwitchEntrainment)
    {
        FOR_ROW_COL_MV
        {
            TotalEntrainmentDep->Drc += EntrainmentDep->Drc;
            TotalEntrainmentDet->Drc += EntrainmentDet->Drc;
            EntrainmentDep->Drc = 0;
            EntrainmentDet->Drc = 0;
        }
        if(SwitchIncludeChannel)
        {
            FOR_ROW_COL_MV_CH
            {
                ChannelTotalEntrainmentDep->Drc += ChannelEntrainmentDep->Drc;
                ChannelTotalEntrainmentDet->Drc += ChannelEntrainmentDet->Drc;
                ChannelEntrainmentDep->Drc = 0;
                ChannelEntrainmentDet->Drc = 0;
            }
        }
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
        MBeM3 = (RainTot + SnowTot + UF_InitializedF + WaterVolSoilTot + floodVolTotInit + BaseFlow +
                 - IntercTot - IntercHouseTot - InfilTot - WaterVolTot - floodVolTot - Qtot - floodBoundaryTot);
        MB = MBeM3/(RainTot + SnowTot + WaterVolSoilTot + floodVolTotInit)*100;
    }

    // Mass Balance sediment, all in kg

    //VJ 110825 forgot to include channeldettot in denominator in MBs!
    if (SwitchErosion && SoilLossTot > 1e-9)
        MBs = (1-(DetTot + ChannelDetTot - SedTot - ChannelSedTot +
                  DepTot + ChannelDepTot)/(SoilLossTot))*100;
    //VJ 121212 changed to mass balance relative to soil loss
}
//---------------------------------------------------------------------------
