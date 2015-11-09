
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
  \file lisErosion.cpp
  \brief Flow and splash detachment functions for slopes and channels

functions: \n
- double TWorld::MaxConcentration(double watvol, double sedvol, double dep) \n
- void TWorld::SplashDetachment(void) \n
- void TWorld::FlowDetachment(void) \n
- void TWorld::ChannelFlowDetachment(void) \n
 */

#include <algorithm>
#include "operation.h"
#include "model.h"


//---------------------------------------------------------------------------
double TWorld::MaxConcentration(double watvol, double sedvol)
{
   double conc = MAXCONC;

   if (watvol > _dx*_dx*1e-6)
      conc = sedvol/watvol;
   // 1e-6 is 1 ml/m2
   else
      conc = 0;

   if (conc > MAXCONC)
      conc = MAXCONC;

   sedvol = conc * watvol;

   return conc;
}
//---------------------------------------------------------------------------
void TWorld::SplashDetachment(void)
{
   if (!SwitchErosion)
      return;

   if(SwitchUseGrainSizeDistribution)
   {

   }

   FOR_ROW_COL_MV
   {
      double b, strength, DetDT1 = 0, DetDT2 = 0, DetLD1, DetLD2;
      double g_to_kg = 0.001;

      double Int = Rain->Drc * 3600/_dt * 1000;
      // intensity in mm/h, Rain is in m

      double KE_DT = 0.0;

      switch (KEequationType)
      {
         case KE_EXPFUNCTION: KE_DT = KEParamater_a1*(1-(KEParamater_b1*exp(-KEParamater_c1*Int))); break;
         case KE_LOGFUNCTION: KE_DT = (Int > 1 ? KEParamater_a2 + KEParamater_b2*log10(Int) : 0); break;
         case KE_POWERFUNCTION: KE_DT = KEParamater_a3*pow(Int, KEParamater_b3); break;
            // kin energy in J/m2/mm
      }
      //VJ 110706  KE equations

      double directrain = (1-Cover->Drc)*Rainc->Drc * 1000; //
      // rainfall between plants in mm

      double KE_LD = std::max(15.3*sqrt(PlantHeight->Drc)-5.87, 0.0);
      // kin energy in J/m2/mm
      double throughfall = Cover->Drc * LeafDrain->Drc * 1000;
      // leaf drip in mm, is calculated as plant leaf drip in interception function so mult cover
      // VJ 110206 stemflow is also accounted for

      double WH0 = exp(-1.48*WH->Drc*1000);
      // water buffer effect on surface, WH in mm in this empirical equation from Torri ?

      if (AggrStab != 0)
      {
         strength = 2.82/AggrStab->Drc;
         b = 2.96;
      }
      else
      {
         strength = 0.1033/CohesionSoil->Drc;
         b = 3.58;
      }
      // empirical analysis based on Limburg data, dating 1989

      // Between plants, directrain is already with 1-cover
      DetDT1 = g_to_kg * fpa->Drc*(strength*KE_DT*WH0+b) * directrain;
      //ponded areas, kg/m2/mm * mm = kg/m2
      DetDT2 = g_to_kg * (1-fpa->Drc)*(strength*KE_DT+b) * directrain * SplashDelivery;
      //dry areas, kg/m2/mm * mm = kg/m2

      if (SwitchKETimebased)
      {
         if (directrain > 0)
         {
            DetDT1 = g_to_kg * fpa->Drc*(strength*KE_DT*WH0+b) * _dt/3600;
            //ponded areas, kg/m2/sec * sec = kg/m2
            DetDT2 = g_to_kg * (1-fpa->Drc)*(strength*KE_DT+b) * _dt/3600 * SplashDelivery;
            //dry areas, kg/m2/sec * sec = kg/m2
         }
      }
      //based on work by Juan Sanchez

      // Under plants, throughfall is already with cover
      DetLD1 = g_to_kg * fpa->Drc*(strength*KE_LD*WH0+b) * throughfall;
      //ponded areas, kg/m2/mm * mm = kg/m2
      DetLD2 = g_to_kg * (1-fpa->Drc)*(strength*KE_LD+b) * throughfall * SplashDelivery;
      //dry areas, kg/m2/mm * mm = kg/m2

      DETSplash->Drc = DetLD1 + DetLD2 + DetDT1 + DetDT2;
      // Total splash kg/m2

      // Deal with all exceptions:

      DETSplash->Drc *= (SoilWidthDX->Drc*DX->Drc);
      // kg/cell, only splash over soilwidth, not roads and channels
      // FROM KG/M2 TO KG/CELL

      DETSplash->Drc = (1-StoneFraction->Drc) * DETSplash->Drc;
      // no splash on stone surfaces

      if (GrassFraction->Drc > 0)
         DETSplash->Drc = (1-GrassFraction->Drc) * DETSplash->Drc;
      // no splash on grass strips

      if (SwitchBuffers && !SwitchSedtrap)
         if(BufferID->Drc > 0)
            DETSplash->Drc = 0;
      // no splash in buffers, but sedtrap can have splash

      if (SwitchHardsurface)
         DETSplash->Drc = (1-HardSurface->Drc)*DETSplash->Drc;
      // no splash on hard surfaces

      if (SwitchHouses)
         DETSplash->Drc = (1-HouseCover->Drc)*DETSplash->Drc;
      //is already contained in soilwidth
      // no splash from house roofs

      DETSplash->Drc = (1-Snowcover->Drc)*DETSplash->Drc;
      // no splash on snow deck

      // IN KG/CELL
   }
}
//---------------------------------------------------------------------------
double TWorld::OFTC(int r, int c, double d)
{

        CG->Drc = pow((d+5)/0.32, -0.6);
        DG->Drc = pow((d+5)/300, 0.25);
        //### Calc transport capacity
        double omega = 100* V->Drc*K2DSlope->Drc;
        // V in cm/s in this formula assuming grad is SINE
        double omegacrit = 0.4;
        // critical unit streampower in cm/s
        return std::min(MAXCONC, 2650 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));
        // not more than 2650*0.32 = 848 kg/m3



}
//---------------------------------------------------------------------------
double TWorld::GetTotalDW(int r, int c,QList<cTMap *> *M)
{
    double wtotal = 0;
    FOR_GRAIN_CLASSES
    {
        wtotal += (*M).Drcd;
    }
    return wtotal;
}
//---------------------------------------------------------------------------
double TWorld::GetDp(int r, int c,double p)
{
    return GetDpMat(r,c,p,&W_D);
}
//---------------------------------------------------------------------------
double TWorld::GetDpMat(int r, int c,double p,QList<cTMap *> *M)
{
    if(numgrainclasses == 1)
    {
        graindiameters.at(0);
    }
    double wtotal = 0;
    FOR_GRAIN_CLASSES
    {
        wtotal += (*M).Drcd;
    }
    wtotal = wtotal*p;
    double w = (*M).at(0)->Drc;;
    FOR_GRAIN_CLASSES
    {

        if( d == numgrainclasses - 1)
        {
            return graindiameters.at(numgrainclasses-1);
        }
        w += (*M).at(d+1)->Drc;
        if(w > wtotal)
        {
            double wmin = (w - (*M).at(d+1)->Drc);
            double wmax = w;
            double dw = wmax-wmin;
            double f = (wtotal- wmin)/dw;
            return f* graindiameters.at(d + 1) + (1.0-f) * graindiameters.at(d);

        }


    }
    return graindiameters.at(numgrainclasses-1);
}
//---------------------------------------------------------------------------
double TWorld::GetSV(int r, int c,double d)
{
    return 2*(2650-1000)*9.80*pow(d/2000000, 2)/(9*0.001);
}

//---------------------------------------------------------------------------
// IN KG/CELL
void TWorld::FlowDetachment(void)
{
   if (!SwitchErosion)
      return;

   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, -1, -1, -1, 0, 0, 0, 1, 1, 1};

   //transport capacity
   FOR_ROW_COL_MV
   {
        if(!SwitchUseGrainSizeDistribution)
        {
            TC->Drc = OFTC(r,c,D50->Drc);
        }else
        {
            TC->Drc = OFTC(r,c,GetDp(r,c,0.5));
        }
   }

   //VJ 110829 TC cannot be more than surrounding cells, this limits the spikes in deposition and erosion
   if (SwitchLimitTC)
   {
      FOR_ROW_COL_MV
      {
         double maxtc = 0;
         double avgtc = 0;
         int count = 0;

         for (int i = 1; i <= 9; i++)
            if(i != 5)
            {
               if ((r+dx[i] >= 0 && c+dy[i] >= 0 && r+dx[i] < _nrRows && c+dy[i] < _nrCols)
                   && !pcr::isMV(TC->data[r+dx[i]][c+dy[i]]))
               {
                       avgtc = avgtc + TC->data[r+dx[i]][c+dy[i]];
                       maxtc = std::max(maxtc,TC->data[r+dx[i]][c+dy[i]]);
                       count++;

               }
            }

         double tcold = TC->Drc;

         if(count > 0)
         {
             TC->Drc = std::max(TC->Drc, avgtc/count);
         }
         TC->Drc = std::min(TC->Drc, maxtc);


      }
   }


   FOR_ROW_COL_MV
   {
      double erosionwh = WH->Drc;
      double erosionwv = WaterVolall->Drc;
      if(SwitchKinematic2D > 1)
      {
          erosionwh = std::max(WHrunoff->Drc-K2DWHStore->Drc ,0.0);
          erosionwv = std::max(WHrunoff->Drc-K2DWHStore->Drc ,0.0)*ChannelAdj->Drc*DX->Drc;
      }

      //### Add splash to sediment
      Sed->Drc += DETSplash->Drc;

      //Find the distribution of detached splash sediment
      double wtotal = GetTotalDW(r,c,&W_D);
      FOR_GRAIN_CLASSES
      {
           Sed_D.Drcd += DETSplash->Drc * W_D.Drcd/wtotal;
           Conc_D.Drcd = MaxConcentration(erosionwv, Sed_D.Drcd);

      }

      // add splash to sed volume

      //### calc concentration and net transport capacity
      DEP->Drc = 0;
      // init deposition for this timestep
      Conc->Drc = MaxConcentration(erosionwv, Sed->Drc);
      // limit sed concentration to max

      double maxTC = std::max(TC->Drc - Conc->Drc,0.0);
      // positive difference: TC deficit becomes detachment (ppositive)
      double minTC = std::min(TC->Drc - Conc->Drc,0.0);
      // negative difference: TC surplus becomes deposition (negative)
      // unit kg/m3

      //### detachment
      double TransportFactor = 0;
      if(!SwitchUseGrainSizeDistribution)
      {
            TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
      }else
      {
          TransportFactor = _dt*GetSV(r,c,GetDpMat(r,c,0.5,&W_D)) * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
      }
      // detachment can only come from soil, not roads (so do not use flowwidth)
      // units s * m/s * m * m = m3

      DETFlow->Drc = Y->Drc * maxTC * TransportFactor;
      // unit = kg/m3 * m3 = kg
      if(SwitchKinematic2D > 1)
      {
          DETFlow->Drc = std::min(DETFlow->Drc, maxTC * erosionwv);
      }else
      {
          DETFlow->Drc = std::min(DETFlow->Drc, maxTC * Q->Drc * _dt);
      }
      // cannot have more detachment than remaining capacity in flow
      // use discharge because standing water has no erosion

      if (GrassFraction->Drc > 0)
         DETFlow->Drc = (1-GrassFraction->Drc) * DETFlow->Drc;
      // no flow detachment on grass strips

      // Detachment edxceptions:
      DETFlow->Drc = (1-StoneFraction->Drc) * DETFlow->Drc ;
      // no flow detachment on stony surfaces

      if (SwitchHardsurface)
         DETFlow->Drc = (1-HardSurface->Drc) * DETFlow->Drc ;
      // no flow detachment on hard surfaces

      if (SwitchHouses)
         DETFlow->Drc = (1-HouseCover->Drc)*DETFlow->Drc;
      // no flow det from house roofs


      //### sediment balance
      Sed->Drc += DETFlow->Drc;



      //Find the distribution of detached sediment
      wtotal = GetTotalDW(r,c,&W_D);
      FOR_GRAIN_CLASSES
      {
           Sed_D.Drcd += DETFlow->Drc * W_D.Drcd/wtotal;
           Conc_D.Drcd = MaxConcentration(erosionwv, Sed_D.Drcd);

      }

      double deposition = 0;

      //### deposition
      if (WH->Drc > MIN_HEIGHT)
      {
          if(!SwitchUseGrainSizeDistribution)
          {
               TransportFactor = (1-exp(-_dt*SettlingVelocity->Drc/erosionwh)) * erosionwv;
          }else
          {
             TransportFactor = (1-exp(-_dt*GetSV(r,c,GetDpMat(r,c,0.5,&Sed_D))/erosionwh)) * erosionwv;
          }

      }else
      {
         TransportFactor = erosionwv;
      }
      // if settl velo is very small, transportfactor is 0 and depo is 0
      // if settl velo is very large, transportfactor is 1 and depo is max

      //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
      // deposition can occur on roads and on soil (so use flowwidth)

      deposition = minTC * TransportFactor;
      // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

      if (SwitchLimitDepTC)
         deposition = std::max(deposition, minTC * erosionwv);
      // cannot be more than sediment above capacity
      deposition = std::max(deposition, -Sed->Drc);
      // cannot have more depo than sediment present
      /* TODO what about this: which one to choose */

      //deposition = (1-Snowcover->Drc) * deposition;
      /* TODO: TRUE??? no deposition on snow */
      if (GrassFraction->Drc > 0)
         deposition = -Sed->Drc*GrassFraction->Drc + (1-GrassFraction->Drc)*deposition;
      // generate 100% deposition on grassstrips
      //? bit tricky, maximizes effect on grassstrips ?


      wtotal = GetTotalDW(r,c,&W_D);
      FOR_GRAIN_CLASSES
      {
           Sed_D.Drcd += deposition * W_D.Drcd/wtotal;
           Conc_D.Drcd = MaxConcentration(erosionwv, Sed_D.Drcd);
      }

      DEP->Drc += deposition;
      // IN KG/CELL
      Sed->Drc += deposition;

      Conc->Drc = MaxConcentration(erosionwv, Sed->Drc);
      // limit concentration to 850 and throw rest in deposition
   }
}
//---------------------------------------------------------------------------
void TWorld::ChannelFlowDetachment(void)
{
   if (!SwitchErosion)
      return;

   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, -1, -1, -1, 0, 0, 0, 1, 1, 1};

   FOR_ROW_COL_MV_CH
   {
      double omega = 100*ChannelV->Drc*ChannelGrad->Drc;
      double omegacrit = 0.4;

      ChannelTC->Drc = std::min(MAXCONC, 2650.0 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));
      // Channel transport capacity, ot more than 2650*0.32 = 848 kg/m3
   }

   //VJ 110829 TC cannot be more than surrounding cells, this limits the spikes in deposition and erosion
   if (SwitchLimitTC)
   {
      FOR_ROW_COL_MV
      {
         double maxtc = 0;
         double avgtc = 0;
         int count = 0;

         for (int i = 1; i <= 9; i++)
            if(i != 5)
            {
               if ((r+dx[i] >= 0 && c+dy[i] >= 0 && r+dx[i] < _nrRows && c+dy[i] < _nrCols)
                   && !pcr::isMV(ChannelTC->data[r+dx[i]][c+dy[i]]))
               {
                  avgtc = avgtc + ChannelTC->data[r+dx[i]][c+dy[i]];
                  maxtc = std::max(maxtc,ChannelTC->data[r+dx[i]][c+dy[i]]);
                  count++;
               }
            }
         ChannelTC->Drc = std::max(ChannelTC->Drc,  avgtc/count);
         ChannelTC->Drc = std::min(ChannelTC->Drc, maxtc);
      }
   }

   FOR_ROW_COL_MV
   {
      // ChannelTC->Drc = tm->Drc;

      ChannelDep->Drc = 0;
      ChannelDetFlow->Drc = 0;

      ChannelSed->Drc += SedToChannel->Drc;
      // add sed flow into channel from slope

      ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc);
      // set conc to max and add surplus sed to ChannelDep

      double maxTC = std::max(ChannelTC->Drc - ChannelConc->Drc,0.0);
      double minTC = std::min(ChannelTC->Drc - ChannelConc->Drc,0.0);
      // basic TC deficit, posiutive if detachment, negative if depositioin, units kg/m3

      double TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * ChannelWidthUpDX->Drc;
      // units s * m/s * m * m = m3
      ChannelDetFlow->Drc = ChannelY->Drc * maxTC * TransportFactor;
      // unit kg/m3 * m3 = kg
      ChannelDetFlow->Drc = std::min(ChannelDetFlow->Drc, maxTC * ChannelWaterVol->Drc);
      // cannot have more detachment than remaining capacity in flow

      //or:? ChannelDetFlow->Drc = ChannelY->Drc * maxTC * Q->Drc*_dt;
      // simple erosion formula

      double deposition = minTC * TransportFactor;
      // max deposition in kg/s  < 0
      if (SwitchLimitDepTC)
         deposition = std::max(deposition, -minTC * ChannelWaterVol->Drc);


      deposition = std::max(deposition, -ChannelSed->Drc);
      // cannot be more than sediment above capacity
      //or:?     deposition = std::max(deposition, minTC * ChannelWaterVol->Drc);

      ChannelDep->Drc += deposition;
      ChannelSed->Drc += deposition;
      ChannelSed->Drc += ChannelDetFlow->Drc;

      ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc);

      ChannelQs->Drc = ChannelQ->Drc * ChannelConc->Drc;

   }
}
//---------------------------------------------------------------------------

void TWorld::SumSedimentClasses(void)
{
    fill(*Sed,0.0);

    for(int i = 0; i < numgrainclasses; i++)
    {
        FOR_ROW_COL_MV
        {
            //Sed->Drc += Sed_D.at(i)->Drc;
        }
    }

    FOR_ROW_COL_MV
    {
        //Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
    }

}
//---------------------------------------------------------------------------
