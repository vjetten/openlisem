
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
      RainAvgmm = Rain->MapAverage()*1000.0;
      RainTotmm += RainAvgmm;
      // avg area rainfall in mm

      tm->calc2V(Rain, (_dx*_dx), MUL); //in m3
      rainfall = tm->MapTotal();
      RainTot += rainfall; // in m3

      oldrainpeak = Rainpeak;
      Rainpeak = max(Rainpeak, rainfall);
      if (oldrainpeak  < Rainpeak)
         RainpeakTime = time;
   }

   if (SwitchSnowmelt)
   {
      SnowAvgmm += Snowmelt->MapAverage()*1000;
      SnowTotmm += SnowAvgmm;

      tm->calc2V(Snowmelt, (_dx*_dx), MUL); //in m3
      snowmelt = tm->MapTotal();
      SnowTot += snowmelt; // in m3

      oldsnowpeak = Snowpeak;
      Snowpeak = max(Snowpeak, snowmelt);
      if (oldsnowpeak < Snowpeak)
         SnowpeakTime = time;
   }

   IntercTot = Interc->MapTotal();
   IntercTotmm = IntercTot*catchmentAreaFlatMM;
   // interception in mm and m3

   InfilTot += InfilVol->MapTotal() + InfilVolKinWave->MapTotal(); //m3
   InfilKWTot += InfilVolKinWave->MapTotal(); // not really used, available for output when needed
   InfilTotmm = max(0,InfilTot*catchmentAreaFlatMM);
   // infiltration mm and m3

   tm->calc2V(WHstore, 1000, MUL); //mm
   SurfStoremm = tm->MapAverage();
   // surface storage CHECK THIS
   // does not go to MB, is already in tot water vol

   WaterVolTot = WaterVolall->MapTotal();//m3
   WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
   // water on the surface in runoff in m3 and mm
   //NOTE: surface storage is already in here so does not need to be accounted for in MB

   //   Qtot += Qoutflow->DrcOutlet;
   // sum outflow m3 for all timesteps for the outlet, is already mult by dt
   Qtot += Qn->DrcOutlet*_dt;
   // sum outflow m3 for all timesteps for the outlet, in m3
   // needed for mass balance
   //Qtotmm = Qtot*catchmentAreaFlatMM;
   // in mm for screen output

   //   QtotOutlet += Qoutflow->DrcOutlet;
   QtotOutlet += Qn->DrcOutlet*_dt;
   // for screen output, total main outlet in m3
   //   QtotPlot += Qoutflow->DrcPlot;
   QtotPlot += Qn->DrcPlot * _dt;
   //VJ 110701 for screen output, total of hydrograph point in m3

   TotalWatervol->copy(WaterVolall);
   // for sed conc calc output


   if (SwitchIncludeChannel)
   {
      WaterVolTot += ChannelWaterVol->MapTotal(); //m3
      // add channel vol to total
      WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
      // recalc in mm for screen output

      //Qtot += ChannelQoutflow->DrcOutlet;
      Qtot += ChannelQn->DrcOutlet*_dt;
      // add channel outflow (in m3) to total for all pits
      //Qtotmm = Qtot*catchmentAreaFlatMM;
      // recalc in mm for screen output

      //QtotOutlet += ChannelQoutflow->DrcOutlet;
      QtotOutlet += ChannelQn->DrcOutlet * _dt;
      // add channel outflow (in m3) to total for main outlet
      //QtotPlot += ChannelQoutflow->DrcPlot;
      QtotPlot += ChannelQn->DrcPlot * _dt;
      // add channel outflow (in m3) to total for main outlet
      TotalWatervol->calc(ChannelWaterVol,ADD);
      // add channel volume to total for sed conc calc

   }

   if (SwitchIncludeTile)
   {
      WaterVolSoilTot = TileWaterVolSoil->MapTotal();
      // input for mass balance, is the water seeping from the soil, input
      // this is the water before the kin wave
      tm->calc2(TileDrainSoil, TileWidth, MUL); //in m3
      tm->calc(TileDX, MUL); //in m3
      // tm->calcV(_dx, MUL); //in m3 ??? or DX?
      TileVolTot += tm->MapTotal(); // in m3

      // water after kin wave
      WaterVolTot += TileWaterVol->MapTotal(); //m3
      // add tile vol to total
      WaterVolTotmm = WaterVolTot*catchmentAreaFlatMM; //mm
      // recalc in mm for screen output

      //Qtot += TileQoutflow->DrcOutlet;
      Qtot += TileQn->DrcOutlet * _dt;
      // add tile outflow (in m3) to total for all pits
      //Qtotmm = Qtot*catchmentAreaFlatMM;
      // recalc in mm for screen output

      //QtotOutlet += TileQoutflow->DrcOutlet;
      QtotOutlet += TileQn->DrcOutlet * _dt;
      // add channel outflow (in m3) to total for main outlet
      //QtotPlot += TileQoutflow->DrcPlot;
      QtotPlot += TileQn->DrcPlot * _dt;
      // add channel outflow (in m3) to total for main outlet
      TotalWatervol->calc(TileWaterVol,ADD);
      // add channel volume to total for sed conc calc

   }

   if (SwitchBuffers)
   {
      BufferVolTot = BufferVol->MapTotal(); // in m3
      if (SwitchIncludeChannel)
         BufferVolTot += ChannelBufferVol->MapTotal();
      //sum up all volume remaining in all buffers (so the non-water!)
      BufferVolin = BufferVolTotInit - BufferVolTot;
      //subtract this from the initial volume to get the total water inflow in the buffers
   }

   // output fluxes for reporting to file and screen in l/s!
   FOR_ROW_COL_MV
   {
      Qoutput->Drc = 1000*qMax(1e-6, Qn->Drc + ChannelQn->Drc + TileQn->Drc); // in l/s
      // added maximum here to avoid strange maps
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
      DetSplashTot += DETSplash->MapTotal();
      DetFlowTot += DETFlow->MapTotal();
      DepTot += DEP->MapTotal();
      DetTot += DETSplash->MapTotal() + DETFlow->MapTotal();
      SedTot = Sed->MapTotal();
      // all in kg/cell

      //SoilLossTot += Qsoutflow->DrcOutlet;
      SoilLossTot += Qsn->DrcOutlet * _dt;
      // sum all sed in all pits (in kg), needed for mass balance

      //SoilLossTotOutlet += Qsoutflow->DrcOutlet;
      SoilLossTotOutlet += Qsn->DrcOutlet * _dt;
      // for screen output, total main outlet sed loss in kg
      TotalSed->copy(Sed);
      // for sed conc

      if (SwitchIncludeChannel)
      {
         // units here in kg, conversion to ton in report functions
         ChannelDetTot += ChannelDetFlow->MapTotal();
         ChannelDepTot += ChannelDep->MapTotal();
         ChannelSedTot = ChannelSed->MapTotal();

         //SoilLossTot += ChannelQsoutflow->MapTotal();
         SoilLossTot += ChannelQsn->DrcOutlet * _dt;
         // add sed outflow for all pits to total soil loss

         //SoilLossTotOutlet += ChannelQsoutflow->DrcOutlet;
         SoilLossTotOutlet += ChannelQsn->DrcOutlet * _dt;
         // add channel outflow (in kg) to total for main outlet

         TotalSed->calc(ChannelSed, ADD);
         // needed for sed conc in file output
      }
TotalSed->report("sed");
      FOR_ROW_COL_MV
      {
         TotalConc->Drc = min(MAXCONC,(TotalWatervol->Drc > _dx*_dx*1e-6? TotalSed->Drc/TotalWatervol->Drc : 0));
      }
      // for file output

      if (SwitchBuffers || SwitchSedtrap)
      {
         BufferSedTot = BufferSed->MapTotal();
         if (SwitchIncludeChannel)
            BufferSedTot += ChannelBufferSed->MapTotal();

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
            - IntercTot - InfilTot - WaterVolTot - Qtot - BufferVolin)/
            (RainTot + SnowTot + WaterVolSoilTot)*100;
   //watervoltot includes channel and tile

   //qDebug() << RainTot << IntercTot << InfilTot << WaterVolTot << BufferVolin << Qtot<< InfilKWTot;

   // Mass Balance sediment, all in kg
   if (SwitchErosion && DetTot > 0)
      MBs = (DetTot + ChannelDetTot - SoilLossTot - SedTot - ChannelSedTot +
             DepTot + ChannelDepTot - BufferSedTot)/(DetTot + ChannelDetTot)*100;
   //VJ 110825 forgot to include channeldettot in denominator in MBs!
}
//---------------------------------------------------------------------------
