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
     double conc = (WaterVol->Drc > 0 ? SedVol->Drc/WaterVol->Drc: 0);
      if (conc > MAXCONC)
      {
         double sedvolnew = MAXCONC*WaterVol->Drc;
         DEP->Drc += sedvolnew-SedVol->Drc;
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
       double Int = Rain->Drc * 3600000/_dt * DX->Drc/_dx;
       // intensity in mm/h
       double KErain = 28.3*(1-(0.52*exp(-0.042*Int)));
       // kin energy in J/m2/mm
       double directrain = (1-Cover->Drc)*RainNet->Drc;
       double KE_DT = KErain;// * directrain;
       // kin energy direct throughfall, in J/m2
       double KEdrain = max(15.3*sqrt(PlantHeight->Drc)-5.87, 0);
       // kin energy in J/m2/mm
       double throughfall = Cover->Drc*RainNet->Drc;
       double KE_LD = KEdrain;// * throughfall;
       // kin energy leaf drainage, in J/m2
       double WH0 = exp(-1.48*WH->Drc);
       // water buffer effect on surface
       // ponded area fraction
      // double CohesionSoil = (1-Cover->Drc)*Cohesion->Drc + Cover->Drc*RootCohesion->Drc;
       // total soil cohesion
       double b, strength, DetDT1, DetDT2, DetLD1, DetLD2;

       if (AggrStab > 0)
       {
          strength = 2.82/AggrStab->Drc; b = 2.96;
       }
       else
       {
          strength = 0.1033/CohesionSoil->Drc; b = 3.58;
       }

       DetDT1 = fpa->Drc*(strength*KE_DT*WH0+b) * directrain;
       //ponded areas between plants, kg/m2/mm * mm = kg/m2
       DetDT2 = (1-fpa->Drc)*(strength*KE_DT+b) * directrain * SplashDelivery;
       //dry areas between plants, kg/m2/mm * mm = kg/m2

       DetLD1 = fpa->Drc*(strength*KE_LD*WH0+b) * throughfall;
       //ponded areas under plants, kg/m2/mm * mm = kg/m2
       DetLD2 = (1-fpa->Drc)*(strength*KE_LD+b) * throughfall * SplashDelivery;
       //dry areas under plants, kg/m2/mm * mm = kg/m2

       DETSplash->Drc = Cover->Drc*(DetLD1 + DetLD2) + (1-Cover->Drc)*(DetDT1 + DetDT2);
       // kg/m2
       DETSplash->Drc *= (SoilWidthDX->Drc*DX->Drc);
       // kg/cell, only splash over soilwidth, not roads and channels

       DETSplash->Drc = (1-StoneFraction->Drc) * DETSplash->Drc;
       // no splash on stone surfaces

   }
}
//---------------------------------------------------------------------------
void TWorld::FlowDetachment(void)
{
   if (!SwitchErosion)
      return;

   FOR_ROW_COL_MV
   {
      double omega = 100*V->Drc*Grad->Drc;
      double omegacrit = 0.4;

      DEP->Drc = 0;
      // init deposition for this timestep

      TC->Drc = min(MAXCONC, 2650 * CG->Drc * pow(max(0, omega - omegacrit), DG->Drc));
      // not more than 2650*0.32 = 848 kg/m3

      if (WH->Drc < 0.00001)
      {
         DETSplash->Drc = 0;
         DETFlow->Drc = 0;
         DEP->Drc = -SedVol->Drc;
         SedVol->Drc = 0;
         continue;
      }
      SedVol->Drc += DETSplash->Drc;
       // no streampower no detachment
       // add splash to sed volume

      Conc->Drc = MaxConcentration(r, c);

      double TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
      double maxTC = max(TC->Drc - Conc->Drc,0);
      double minTC = min(TC->Drc - Conc->Drc,0);

      DETFlow->Drc = Y->Drc * maxTC * TransportFactor;
      // unit = m/s * m*m = m3/s,  detachment in kg/s : kg/m3 * m3/s
      DETFlow->Drc = min(DETFlow->Drc, maxTC * WaterVol->Drc);
      // cannot have more detachment than remaining capacity in flow

      double deposition = minTC * TransportFactor;
      // max deposition in kg/s  < 0
      deposition = max(deposition, minTC * WaterVol->Drc);
      // cannot be more than sediment above capacity

      double sedvolume = SedVol->Drc + DETFlow->Drc + deposition;
      // temp sed balance
      if (sedvolume < 0.000001) //too much sediment, increase DEP with net deposition
      {
         DEP->Drc += deposition;
         SedVol->Drc = 0;
      }
      else //Sed volume exists, increase SedVol with flow det
         SedVol->Drc += DETFlow->Drc;

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

      ChannelDep->Drc = 0;
      ChannelDetFlow->Drc = 0;

      ChannelSedVol->Drc += SedToChannel->Drc;
      // add sed flow into channel from slope
     
      if (ChannelWH->Drc < 0.00001)
      {
         ChannelDep->Drc -= ChannelSedVol->Drc;
         ChannelSedVol->Drc = 0;
      }

      ChannelTC->Drc = min(MAXCONC, 2650 * CG->Drc * pow(max(0, omega - omegacrit), DG->Drc));
      // not more than 2650*0.32 = 848 kg/m3

      ChannelConc->Drc = MaxChannelConcentration(r, c);

      double TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * ChannelWidth->Drc;
      double maxTC = max(ChannelTC->Drc - ChannelConc->Drc,0);
      double minTC = min(ChannelTC->Drc - ChannelConc->Drc,0);

      ChannelDetFlow->Drc = ChannelY->Drc * maxTC * TransportFactor;
      // unit = m/s * m*m = m3/s,  detachment in kg/s : kg/m3 * m3/s
      ChannelDetFlow->Drc = min(ChannelDetFlow->Drc, maxTC * ChannelWaterVol->Drc);
//or:? ChannelDetFlow->Drc = ChannelY->Drc * maxTC * Q->Drc*_dt;
      // cannot have more detachment than remaining capacity in flow

      double deposition = minTC * TransportFactor;
      // max deposition in kg/s  < 0
      deposition = max(deposition, minTC * ChannelWaterVol->Drc);
  //or:?    deposition = max(deposition, -ChannelSedVol->Drc);
      // cannot be more than sediment above capacity

      double sedvolume = ChannelSedVol->Drc + ChannelDetFlow->Drc + deposition;
      // temp sed balance
      if (sedvolume < 0.000001) //too much sediment, increase DEP with deposition
      {
         ChannelSedVol->Drc = 0;
         ChannelDep->Drc += deposition;
     }
      else //Sed volume exists, increase SedVol with flow det
         ChannelSedVol->Drc += ChannelDetFlow->Drc;

      ChannelConc->Drc = MaxChannelConcentration(r, c);
   }
}
//---------------------------------------------------------------------------

