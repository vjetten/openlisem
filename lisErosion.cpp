#include "model.h"


//---------------------------------------------------------------------------
double TWorld::MaxChannelConcentration(int r, int c)
{
     double conc = (ChannelWaterVol->Drc > 0 ? ChannelSedVol->Drc/ChannelWaterVol->Drc: 0);
      if (conc > MAXCONC)
      {
         double sedvolnew = MAXCONC*ChannelWaterVol->Drc;
         ChannelDep->Drc += sedvolnew-ChannelSedVol->Drc;
         ChannelSedVol->Drc = sedvolnew;
         conc = MAXCONC;
      }
      return(conc);
}
//---------------------------------------------------------------------------
double TWorld::MaxConcentration(int r, int c)
{
     // VJ 100502 changed to watervolall
     double conc = (WaterVolall->Drc > 0 ? SedVol->Drc/WaterVolall->Drc : 1000);
//     double conc = (WaterVol->Drc > 0 ? SedVol->Drc/WaterVol->Drc : 1000);
     if (conc > MAXCONC)
     {
         double sedvolnew = MAXCONC*WaterVolall->Drc; //!!
         DEP->Drc += min(0, sedvolnew-SedVol->Drc);
         SedVol->Drc = sedvolnew;
         conc = MAXCONC;
     }
     return(conc);
}
//---------------------------------------------------------------------------
void TWorld::SplashDetachment(void)
{
   if (!SwitchErosion)
      return;

   FOR_ROW_COL_MV
   {
       double b, strength, DetDT1, DetDT2, DetLD1, DetLD2;
       double g_to_kg = 0.001;

	   double Int = Rain->Drc * 3600/_dt * DX->Drc/_dx * 1000;
       // intensity in mm/h, uncorrect for slope, Rain is in m corrected for slope

       double KE_DT = 28.3*(1-(0.52*exp(-0.042*Int)));
       // kin energy in J/m2/mm, Van DIjk general equation 2002

      //TODO: allow many different equations here, interface choice
      /* equation in LISEM, based on Eurosem, Morgan 1998
       if (Int > 1)
           KE_DT = 8.95+8.44*log10(Int);
       else
    	   KE_DT = 0;
      */

       double directrain = (1-Cover->Drc)*Rain->Drc * 1000; //
       // rainfall between plants in mm

       double KE_LD = max(15.3*sqrt(PlantHeight->Drc)-5.87, 0);
       // kin energy in J/m2/mm
       double throughfall = LeafDrain->Drc * (1-StemflowFraction) * 1000;
       // leaf drip in mm is already calculated with plant covern in interception function

       double WH0 = exp(-1.48*WH->Drc*1000);
       // water buffer effect on surface, WH in mm in this empirical equation from Torri ?


       if (AggrStab > 0)
       {
          strength = 2.82/AggrStab->Drc; b = 2.96;
       }
       else
       {
          strength = 0.1033/CohesionSoil->Drc; b = 3.58;
       }

       // Between plants, directrain is already with 1-cover
       DetDT1 = g_to_kg * fpa->Drc*(strength*KE_DT*WH0+b) * directrain;
       //ponded areas, kg/m2/mm * mm = kg/m2
       DetDT2 = g_to_kg * (1-fpa->Drc)*(strength*KE_DT+b) * directrain * SplashDelivery;
       //dry areas, kg/m2/mm * mm = kg/m2

       // Under plants, throughfall is already with cover
       DetLD1 = g_to_kg * fpa->Drc*(strength*KE_LD*WH0+b) * throughfall;
       //ponded areas, kg/m2/mm * mm = kg/m2
       DetLD2 = g_to_kg * (1-fpa->Drc)*(strength*KE_LD+b) * throughfall * SplashDelivery;
       //dry areas, kg/m2/mm * mm = kg/m2

       DETSplash->Drc = DetLD1 + DetLD2 + DetDT1 + DetDT2;
       // Total splash kg/m2

       DETSplash->Drc *= (SoilWidthDX->Drc*DX->Drc);
       // kg/cell, only splash over soilwidth, not roads and channels

       DETSplash->Drc = (1-StoneFraction->Drc) * DETSplash->Drc;
       // no splash on stone surfaces

       if (GrassPresent->Drc > 0)
          DETSplash->Drc = (1-GrassFraction->Drc) * DETSplash->Drc;
       // no splash on grass strips

   }
}
//---------------------------------------------------------------------------
void TWorld::FlowDetachment(void)
{
   if (!SwitchErosion)
      return;

   FOR_ROW_COL_MV
   {
	  //### Calc transport capacity

      double omega = 100*V->Drc*Grad->Drc;
      // V in cm/s in this formula assuming grad is SINE
      double omegacrit = 0.4;
      // critical unit streampower in cm/s
      TC->Drc = min(MAXCONC, 2650 * CG->Drc * pow(max(0, omega - omegacrit), DG->Drc));
      // not more than 2650*0.32 = 848 kg/m3

	  //### Add splash to sediment
      DEP->Drc = 0;
      // init deposition for this timestep

      SedVol->Drc += DETSplash->Drc;
      // add splash to sed volume

	  //### calc concentration and net transport capacity
      Conc->Drc = MaxConcentration(r, c);
      // sed concentration
      double maxTC = max(TC->Drc - Conc->Drc,0);
      double minTC = min(TC->Drc - Conc->Drc,0);
      // unit kg/m3

	  //### detachment
      double TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
      // detachment can only come from soil, not roads (so do not use flowwidth)
      // units s * m/s * m * m = m3

      DETFlow->Drc = Y->Drc * maxTC * TransportFactor;
      // unit = kg/m3 * m3 = kg
      DETFlow->Drc = min(DETFlow->Drc, maxTC * WaterVolall->Drc);
      // cannot have more detachment than remaining capacity in flow

      if (GrassPresent->Drc > 0)
         DETFlow->Drc = (1-GrassFraction->Drc) * DETFlow->Drc;
      // no flow detachment on grass strips

      //TODO no detachment with stoniness


	  //### deposition

      if (WH->Drc > MIN_HEIGHT)
         TransportFactor = (1-exp(-_dt*SettlingVelocity->Drc/WH->Drc)) * WaterVolall->Drc;
      else
         TransportFactor = WaterVolall->Drc;

   //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
      // deposition can occur on roads and on soil (so use flowwidth)

      double deposition = minTC * TransportFactor;
      // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0
      //deposition = max(deposition, minTC * WaterVol->Drc);
      // cannot be more than sediment above capacity
      deposition = max(deposition, -SedVol->Drc);
      // cannot have more depo than sediment present

      if (GrassPresent->Drc > 0)
    	  deposition = -SedVol->Drc*GrassFraction->Drc + (1-GrassFraction->Drc)*deposition;
      // generate deposition on grassstrips
      /*

	  //### sediment balance
      double sedvolume = SedVol->Drc + DETFlow->Drc + deposition;
      // temp sed balance
      if (sedvolume < 0) //too much sediment, increase DEP with net deposition
      {
         DEP->Drc += deposition;
         SedVol->Drc = 0;
      }
      else //Sed volume exists, increase SedVol with flow det
         SedVol->Drc += DETFlow->Drc;
*/
      DEP->Drc += deposition;
      SedVol->Drc += DETFlow->Drc;
      SedVol->Drc += deposition;

      Conc->Drc = MaxConcentration(r, c);
   }
}
//---------------------------------------------------------------------------
void TWorld::ChannelFlowDetachment(void)
{
   if (!SwitchErosion)
      return;

   FOR_ROW_COL_MV_CH
   {
      double omega = 100*ChannelV->Drc*ChannelGrad->Drc;
      double omegacrit = 0.4;


      ChannelTC->Drc = min(MAXCONC, 2650 * CG->Drc * pow(max(0, omega - omegacrit), DG->Drc));
      // Channel transport capacity, ot more than 2650*0.32 = 848 kg/m3

      ChannelDep->Drc = 0;
      ChannelDetFlow->Drc = 0;

      ChannelSedVol->Drc += SedToChannel->Drc;
      // add sed flow into channel from slope

      ChannelConc->Drc = MaxChannelConcentration(r, c);
      // set conc to max and add surplus sed to ChannelDep

      double maxTC = max(ChannelTC->Drc - ChannelConc->Drc,0);
      double minTC = min(ChannelTC->Drc - ChannelConc->Drc,0);
      // basic TC deficit, posiutive if detachment, negative if depositioin, units kg/m3

      double TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * ChannelWidthUpDX->Drc;
      // units s * m/s * m * m = m3
      ChannelDetFlow->Drc = ChannelY->Drc * maxTC * TransportFactor;
      // unit kg/m3 * m3 = kg
      ChannelDetFlow->Drc = min(ChannelDetFlow->Drc, maxTC * ChannelWaterVol->Drc);
      // cannot have more detachment than remaining capacity in flow

      //or:? ChannelDetFlow->Drc = ChannelY->Drc * maxTC * Q->Drc*_dt;
      // simple erosion formula

      double deposition = minTC * TransportFactor;
      // max deposition in kg/s  < 0
      deposition = max(deposition, -ChannelSedVol->Drc);
      // cannot be more than sediment above capacity
      //or:?     deposition = max(deposition, minTC * ChannelWaterVol->Drc);

/*
      double sedvolume = ChannelSedVol->Drc + ChannelDetFlow->Drc + deposition;
      // temp sed balance
      if (sedvolume < 0.000001) //too much sediment, increase DEP with deposition
      {
         ChannelSedVol->Drc = 0;
         ChannelDep->Drc += deposition;
     }
      else //Sed volume exists, increase SedVol with flow det
         ChannelSedVol->Drc += ChannelDetFlow->Drc;
*/

      ChannelDep->Drc += deposition;
      ChannelSedVol->Drc += deposition;
      ChannelSedVol->Drc += ChannelDetFlow->Drc;

      ChannelConc->Drc = MaxChannelConcentration(r, c);
   }
}
//---------------------------------------------------------------------------

