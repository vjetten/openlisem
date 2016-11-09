
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
- double TWorld::MaxConcentration(double watvol, double sedvol)\n
- void TWorld::SplashDetachment(void)\n
- double TWorld::OFTC(int r, int c, int d)\n
- double TWorld::GetTotalDW(int r, int c,QList<cTMap *> *M)\n
- double TWorld::GetDp(int r, int c,double p)\n
- double TWorld::GetDpMat(int r, int c,double p,QList<cTMap *> *M)\n
- double TWorld::GetMpMat(int r, int c,double p,QList<cTMap *> *M, QList<double> *V)\n
- double TWorld::GetSV(double d)\n
- void TWorld::SedimentSetMaterialDistribution(int r,int c)\n
- double TWorld::DetachMaterial(int r,int c, int d,bool channel, bool flood,bool bl,double detachment)\n
- void TWorld::FlowDetachment(void)\n
- void TWorld::ChannelFlowDetachment(int r, int c)\n
- void TWorld::RiverSedimentMaxC(int r, int c)\n
- void TWorld::RiverSedimentDiffusion(double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)\n
- void TWorld::RiverSedimentLayerDepth(int r , int c)\n
- double TWorld::RiverSedimentTCBL(int r,int c,int _d)\n
- double TWorld::RiverSedimentTCSS(int r,int c, int _d)\n
 */

#include <algorithm>
#include "operation.h"
#include "model.h"

#define he_ca 1e-12
#define ve_ca 1e-12
#define GRAV 9.81
#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

//---------------------------------------------------------------------------
/**
 * @fn double TWorld::MaxConcentration(double watvol, double sedvol)
 * @brief Calculates concentration with a maximum of MAXCONC
 *
 * @param watvol : the watervolume
 * @param sedvol : the sediment mass
 * @return sediment concentration (kg/m3)
 * @see MAXCONC
 *
 */
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
/**
 * @fn double TWorld::MaxConcentration(double watvol, double sedvol)
 * @brief Calculates splash detachment for each cell
 *
 * This function calculates splash detachment for each cell,
 * based on kinetic energy of direct rainfall and troughfall.
 * Detachment is taken from the top soil layer if possible.
 * The detached sediment is not directly added to sediment in flow,
 * this happens during flow detachment.
 *
 *
 * @see KEequationType
 *
 */
void TWorld::SplashDetachment(int thread)
{
   if (!SwitchErosion)
      return;


   FOR_ROW_COL_2DMT
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
      double throughfall = (1-Litter->Drc) * Cover->Drc * LeafDrain->Drc * 1000;
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

      if(SwitchUseMaterialDepth)
      {
          double depdepth = std::max((StorageDep->Drc / (1600.0))/(_dx * DX->Drc),0.0);
          double fac1 = std::max(0.0,1.0 - depdepth/SedimentMixingDepth->Drc);
          double fac2 = 1.0 - fac1;

          strength = strength * fac2 + (0.1033/DepositedCohesion) * fac1;
          b = b * fac2 + 3.58 * fac1;
      }


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

      if (SwitchHardsurface)
         DETSplash->Drc = (1-HardSurface->Drc)*DETSplash->Drc;
      // no splash on hard surfaces

      if (SwitchHouses)
         DETSplash->Drc = (1-HouseCover->Drc)*DETSplash->Drc;
      //is already contained in soilwidth
      // no splash from house roofs

      DETSplash->Drc = ((1-Snowcover->Drc)*DETSplash->Drc);
      // no splash on snow deck

      DETSplash->Drc = std::max(0.0,DETSplash->Drc);
      if(SwitchUseMaterialDepth)
      {
          if(SwitchUseGrainSizeDistribution)
          {
              StorageDep->Drc = GetTotalDW(r,c,&StorageDep_D);
              Storage->Drc = GetTotalDW(r,c,&Storage_D);
          }

          //check wat we can detach from the top and bottom layer of present material
          double dleft = DETSplash->Drc ;
          double deptake = 0;
          double mattake = 0;
          double detachment = 0;

          deptake = std::min(dleft,StorageDep->Drc);
          StorageDep->Drc -= deptake;

          if(SwitchUseGrainSizeDistribution)
          {
              double wtotal = 0;
              FOR_GRAIN_CLASSES
              {
                  wtotal += StorageDep_D.Drcd;
              }
              if(wtotal != 0)
              {
                  FOR_GRAIN_CLASSES
                  {
                      StorageDep_D.Drcd = StorageDep_D.Drcd * StorageDep->Drc/wtotal;
                  }
              }
          }

          detachment += deptake;


          if(!(Storage->Drc < 0))
          {
              mattake = std::min(dleft,Storage->Drc);
              Storage->Drc -= mattake;

              if(SwitchUseGrainSizeDistribution)
              {
                  double wtotal = 0;
                  FOR_GRAIN_CLASSES
                  {
                        wtotal += Storage_D.Drcd;
                  }
                  if(wtotal != 0)
                  {
                      FOR_GRAIN_CLASSES
                      {
                            Storage_D.Drcd = Storage_D.Drcd * Storage->Drc/wtotal;
                      }
                  }
              }
              detachment += mattake;
          }else
          {
              detachment += dleft;
          }

          DETSplash->Drc = detachment;
          UF_AddedSplash->Drc = detachment;
      }

      // IN KG/CELL
   }}}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/**
 * @fn double TWorld::GetTotalDW(int r, int c, QList<cTMap *> *M)
 * @brief Returns the total value of a map list on a cell
 *
 * Sums up the value of a map list
 * for a certain cell
 *
 * @param r : row nr of the cell
 * @param c : column nr of the cell
 * @param M : Material map list
 * @return The total value
 */
double TWorld::GetTotalDW(int r, int c,QList<cTMap *> *M)
{
    //simple iteration over maps
    double wtotal = 0;
    FOR_GRAIN_CLASSES
    {
        wtotal += (*M).Drcd;
    }
    return wtotal;
}
//---------------------------------------------------------------------------
/**
 * @fn double TWorld::GetDp(int r, int c,double p)
 * @brief get p percent grain size
 *
 * Uses a linear approximation to find the
 * grain size for which a certain percentage
 * of sediment mass has a lower grain size.
 *
 * @param r : row nr of the cell
 * @param c : column nr of the cell
 * @param p : the mass factor that should have a lower grain size class
 * @return The p percent grain size
 * @see TWorld::GetDpMat
 */
double TWorld::GetDp(int r, int c,double p)
{
    //use more generic function
    return GetDpMat(r,c,p,&W_D);
}
//---------------------------------------------------------------------------
/**
 * @fn double TWorld::GetDpMat(int r, int c,double p,QList<cTMap *> *M)
 * @brief get p percent grain size based on material distribution
 *
 * Uses a linear approximation to find the
 * grain size for which a certain percentage
 * of sediment mass has a lower grain size.
 * Uses a material distribution for each cell
 *
 * @param r : row nr of the cell
 * @param c : column nr of the cell
 * @param p : the mass factor that should have a lower grain size class
 * @param M : Material distribution map list
 * @return The p percent grain size
 */
double TWorld::GetDpMat(int r, int c,double p,QList<cTMap *> *M)
{
    //check if there is a single grain class
    //then we can return the single value
    if(numgrainclasses == 1)
    {
        graindiameters.at(0);
    }
    //find total material
    double wtotal = 0;
    FOR_GRAIN_CLASSES
    {
        wtotal += (*M).Drcd;
    }
    //find factor of material we should reach
    wtotal = wtotal*p;

    //iterate trough grain classes untill enough material is summed
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
            //return linearly interpolated value
            double wmin = (w - (*M).at(d+1)->Drc);
            double wmax = w;
            double dw = wmax-wmin;
            double f = (wtotal- wmin)/dw;
            return f* graindiameters.at(d + 1) + (1.0-f) * graindiameters.at(d);

        }


    }
    //return latest value, either because p is near 1
    //or because the vast majority of material is in this class.
    //more specificly done if: 100 * p > (100 - percentage in last grain class)
    return graindiameters.at(numgrainclasses-1);
}
//---------------------------------------------------------------------------
/**
 * @fn double TWorld::GetMpMat(int r, int c,double p,QList<cTMap *> *M)
 * @brief get p percent grain size based on material distribution
 *
 * Uses a linear approximation to find the
 * value of a parameter,
 * for the p percent grain size class.
 * Uses a material distribution.
 *
 * @param r : row nr of the cell
 * @param c : column nr of the cell
 * @param p : the mass factor that should have a lower grain size class
 * @param M : Material distribution map list
 * @param V : Parameter value map list
 * @return The value of parameter V for the p percent grain size.
 */
double TWorld::GetMpMat(int r, int c,double p,QList<cTMap *> *M, QList<double> *V)
{
    //check if there is a single grain class
    //then we can return the single value
    if(numgrainclasses == 1)
    {
        graindiameters.at(0);
    }

    //find total material
    double wtotal = 0;
    FOR_GRAIN_CLASSES
    {
        wtotal += (*M).Drcd;
    }

    //find factor of material we should reach
    wtotal = wtotal*p;
    //iterate trough grain classes untill enough material is summed
    double w = (*M).at(0)->Drc;;
    FOR_GRAIN_CLASSES
    {

        if( d == numgrainclasses - 1)
        {
            return (*V).at(numgrainclasses-1);
        }
        w += (*M).at(d+1)->Drc;
        if(w > wtotal)
        {
            //return linearly interpolated value
            double wmin = (w - (*M).at(d+1)->Drc);
            double wmax = w;
            double dw = wmax-wmin;
            double f = (wtotal- wmin)/dw;
            return f* (*V).at(d + 1) + (1.0-f) * (*V).at(d);

        }


    }
    //return latest value, either because p is near 1
    //or because the vast majority of material is in this class.
    //more specificly done if: 100 * p > (100 - percentage in last grain class)
    return (*V).at(numgrainclasses-1);
}
/**
 * @fn double TWorld::GetSV(double d)
 * @brief get settling velocity of sediment with grain size d
 *
 * @param d : the grain size (in micrometer)
 * @return The settling velocity
 */
double TWorld::GetSV(double d)
{
    //Stokes range settling velocity
    if(d < 100)
    {
        return 2*(2650-1000)*9.80*pow(d/2000000.0, 2)/(9*0.001);

    //Settling velocity by Zanke (1977)
    }else
    {
        double dm = d/1000.0;
        return 10.0 *(sqrt(1.0 + 0.01 *((2650.0-1000.0)/1000.0)* GRAV *dm*dm*dm )-1.0)/(dm);
    }

}
/**
 * @fn void TWorld::SedimentSetMaterialDistribution(int r,int c)
 * @brief Update grain size distribution of soil layers
 *
 * Update grain size distribution of both
 * original soil layer and deposited layer.
 * Step is skipped if the storage is negative,
 * which indicated infinite storage.
 *
 * @param r : row nr of the cell
 * @param c : coumn nr of the cell
 * @return void
 */
void TWorld::SedimentSetMaterialDistribution()//(int r,int c)
{
    if(!SwitchUseMaterialDepth)
        return;

    FOR_ROW_COL_MV
    {
        //set total mass from grain size distributed mass
        if(SwitchUseGrainSizeDistribution)
        {
            StorageDep->Drc = GetTotalDW(r,c,&StorageDep_D);
            if(!(Storage->Drc < -1))
            {
                Storage->Drc = GetTotalDW(r,c,&Storage_D);
            }

        }

        //update grain size distributed weights
        if(SwitchUseGrainSizeDistribution)
        {

            if(!(Storage->Drc < -1))
            {
                //calculated layer depth based on a bulk density of 1600 kg/m3
                double depdepth = std::max((StorageDep->Drc / (1600.0))/(_dx * DX->Drc),0.0);
                //linear soil layers mixing factor for effective soil properties
                double fac1 = std::max(0.0,1.0 - depdepth/SedimentMixingDepth->Drc);
                double fac2 = 1.0 - fac1;
                if(SedimentMixingDepth->Drc < MIN_HEIGHT)
                {
                    fac1 = 1.0;
                    fac2 = 0.0;
                }

                //use soil layer mixing factor to calculate effective soil properties
                FOR_GRAIN_CLASSES
                {
                    if(StorageDep->Drc > MIN_HEIGHT && Storage->Drc > MIN_HEIGHT)
                    {
                        W_D.Drcd = fac1 * Storage_D.Drcd/Storage->Drc + fac2 * StorageDep_D.Drcd/StorageDep->Drc;
                    }else if(Storage->Drc < MIN_HEIGHT && StorageDep->Drc > MIN_HEIGHT)
                    {
                        W_D.Drcd = StorageDep_D.Drcd/StorageDep->Drc;
                    }else if(StorageDep->Drc < MIN_HEIGHT && Storage->Drc > MIN_HEIGHT)
                    {
                        W_D.Drcd = Storage_D.Drcd/Storage->Drc;
                    }else
                    {
                        W_D.Drcd = 0;
                    }
                    W_D.Drcd =std::max(0.0,W_D.Drcd);

                }

                //normalize (total weight should be 1)
                double wtotal = 0;
                FOR_GRAIN_CLASSES
                {
                    wtotal += W_D.Drcd;
                }
                if(wtotal > 0)
                {
                    FOR_GRAIN_CLASSES
                    {
                        W_D.Drcd =  W_D.Drcd / wtotal;
                    }
                }
            }

            //identical process for channel soil layers
            if(SwitchIncludeChannel)
            {
                RStorageDep->Drc = GetTotalDW(r,c,&RStorageDep_D);
                if(!(RStorage->Drc < -1))
                {
                    RStorage->Drc = GetTotalDW(r,c,&RStorage_D);

                    //calculated layer depth based on a bulk density of 1600 kg/m3
                    double depdepth = std::max((RStorageDep->Drc / (1600.0))/(_dx * DX->Drc),0.0);
                    //linear soil layers mixing factor for effective soil properties
                    double fac1 = std::max(0.0,1.0 - depdepth/RSedimentMixingDepth->Drc);
                    double fac2 = 1 - fac1;
                    if(RSedimentMixingDepth->Drc < MIN_HEIGHT)
                    {
                        fac1 = 1;
                        fac2 = 0;
                    }
                    //use soil layer mixing factor to calculate effective soil properties
                    FOR_GRAIN_CLASSES
                    {
                        if(RStorageDep->Drc > MIN_HEIGHT && RStorage->Drc > MIN_HEIGHT)
                        {
                            RW_D.Drcd = fac1 * RStorage_D.Drcd/RStorage->Drc + fac2 * RStorageDep_D.Drcd/RStorageDep->Drc;
                        }else if(RStorage->Drc < MIN_HEIGHT && RStorageDep->Drc > MIN_HEIGHT)
                        {
                            RW_D.Drcd = RStorageDep_D.Drcd/RStorageDep->Drc;
                        }else if(RStorageDep->Drc < MIN_HEIGHT && RStorage->Drc > MIN_HEIGHT)
                        {
                            RW_D.Drcd = RStorage_D.Drcd/RStorage->Drc;
                        }else
                        {
                            RW_D.Drcd = 0;
                        }

                        RW_D.Drcd =std::max(0.0,RW_D.Drcd);

                    }

                    //normalize (total weight should be 1)
                    double wtotal = 0;
                    FOR_GRAIN_CLASSES
                    {
                        wtotal += RW_D.Drcd;
                    }
                    if(wtotal > 0)
                    {
                        FOR_GRAIN_CLASSES
                        {
                            RW_D.Drcd =  RW_D.Drcd / wtotal;
                        }
                    }
                }
            }
        }
    }
}

/**
 * @fn double TWorld::DetachMaterial(int r,int c, int d,bool channel, bool flood,bool bl,double detachment)
 * @brief Calculates real detachment from potential detachment.
 *
 * This cell uses the real time calculated effective erosion coefficient
 * to calculate actual erosion. When the soil layer has less sediment left then
 * the potential erosion, an analytical solution is used to converge
 * both sediment in the soil layer and sediment in tranport to a
 * stable value. When both the channel and flood parameter are false,
 * overland flow detachment is assumed.
 *
 * @param r : Row nr of the cell
 * @param c : Column nr of the cell
 * @param d : Grain diameter class
 * @param Channel : Channel detachment?
 * @param flood : Flood detachment?
 * @param bl : Bed Load detachment?
 * @param detachment : Potential detachment
 * @return Actual detachment
 */

double TWorld::DetachMaterial(int r,int c, int d,bool channel, bool flood,bool bl,double detachment)
{
    /*
     * NOTE: the actual detachment is immediately taken
     * from the soil layers. It is therefore assumed that when using
     * this function, the actual detachment is added to the
     * sediment in flow. THIS IS REQUIRED! to maintain mass
     * balance
     */

    // when there is no usage of material depth,
    // there is no deposited layer
    // actual erosion can then be calculated using
    // the original erosion efficiency coefficient
    if(!SwitchUseMaterialDepth)
    {
        if(channel)
        {
            return detachment * ChannelY->Drc;
        }else
        {
            return detachment * Y->Drc;
        }
    }

    //first check if it is channel detachment
    if(channel)
    {\
        //calculate depth of deposited layer
        double depdepth = std::max((RStorageDep->Drc / (1600.0))/(ChannelWidth->Drc * DX->Drc),0.0);

        //linear decrease in influence from lower soil layer
        //from 0 to MixingDepth, with fac1 for bottom layer, fac2 for top layer
        double fac1 = std::max(0.0,1.0 - depdepth/RSedimentMixingDepth->Drc);
        double fac2 = 1-fac1;
        if(RSedimentMixingDepth->Drc < MIN_HEIGHT)
        {
            fac1 = 1;
            fac2 = 0;
        }
        //new erosion coefficient bases on soil layer mixinfactors
        double newY = ChannelY->Drc * fac1 + fac2 * 1.0;

        //multiply potential detachment by erosion coefficient
        detachment = detachment *newY;

        //when a grain size distribution is used,
        //transport capacity if dependend on the
        //grain size distribution
        //this is therefore needed to converge to stable values
        if( SwitchUseGrainSizeDistribution)
        {
            //if there is not an indicator for infinite sediment
            //and there is less sediment in the combined layers
            //than will be detached, use analytical solution
            if(!((RStorage->Drc) < -1) &&  detachment > RStorage_D.Drcd + RStorageDep_D.Drcd)
            {
                double conc = (bl? UF1D_blm_D.Drcd:UF1D_ssm_D.Drcd) / ( UF1D_f->Drc);
                double tc = bl? UF1D_bltc_D.Drcd:UF1D_sstc_D.Drcd;

                double masstotal = GetTotalDW(r,c,&RStorage_D)+ GetTotalDW(r,c,&RStorageDep_D);
                double a = RStorage_D.Drcd + RStorageDep_D.Drcd;
                double c2 = masstotal;
                double b = conc;
                double d2 = DX->Drc;
                double e = ChannelWidth->Drc;
                double f = tc;
                double h = UF1D_f->Drc / (_dx * _dx);

                double t1 = c2 - b*d2*e*h + d2*e*f*h;
                double t2_1 = c2 - b*d2*e*h + d2*e*f*h;
                double t2_2 = b*c2*d2*e*h - a * d2 * e * f*h;
                double t2 = sqrt(t2_1 * t2_1 + 4 * t2_2);
                detachment = t1 - t2;

            }

        }
        //to remove small rounding errors that lead to negative values
        detachment = std::max(detachment,0.0);

        //check wat we can detache from the top and bottom layer of present material
        double dleft = detachment;
        double deptake = 0;
        double mattake = 0;
        detachment = 0;

        //when a grain size distribution is used
        if(SwitchUseGrainSizeDistribution)
        {
            //take from soil layer mass within that grain class
            deptake = std::min(dleft,RStorageDep_D.Drcd);
            dleft -= deptake;
            //and take the same from the total
            RStorageDep_D.Drcd -= deptake;
            RStorageDep->Drc -= deptake;
        }else
        {
            //take from the total storage
            deptake = std::min(dleft,RStorageDep->Drc);
            RStorageDep->Drc -= deptake;
        }
        //add to the detachment what we have taken from the first soil layer
        detachment += deptake;

        //if the deposited layer is empty
        //use erosion efficiency of bottom layer again
        if(newY > 0)
        {
            dleft *= ChannelY->Drc/newY;
        }else
        {
            dleft = 0;
        }

        //bottom soil layer can be infinite
        if(!((RStorage->Drc) < -1))
        {
            //take sediment from the soil layer storage
            if(SwitchUseGrainSizeDistribution)
            {
                //take from soil layer mass within that grain class
                mattake = std::min(dleft,RStorage_D.Drcd);
                dleft -= mattake;
                //take from the total storage
                RStorage_D.Drcd -= mattake;
                RStorage->Drc -= mattake;

            }else
            {
                //take from the total storage
                mattake = std::min(dleft,RStorage->Drc);
                RStorage->Drc -= mattake;
            }
            //add to the detachment what we have taken from the second soil layer
            detachment += mattake;
        }else
        {
            //all left potential detachment is added to detachment (infinite soil layer)
            detachment += dleft;
        }

    //finally, return the total detachment
    return std::max(0.0,detachment);

    }else
    {
        double depdepth = std::max((StorageDep->Drc / (1600.0))/(_dx * DX->Drc),0.0);

        //linear decrease in influence from lower soil layer
        //from 0 to MixingDepth, with fac1 for bottom layer, fac2 for top layer
        double fac1 = std::max(0.0,1.0 - depdepth/SedimentMixingDepth->Drc);
        double fac2 = 1 - fac1;
        if(SedimentMixingDepth->Drc < MIN_HEIGHT)
        {
            fac1 = 1;
            fac2 = 0;
        }
        //new erosion coefficient bases on soil layer mixinfactors
        double newY = Y->Drc * fac1 + fac2 * 1.0;
        //multiply potential detachment by erosion coefficient
        detachment = detachment *newY;

        //when a grain size distribution is used,
        //transport capacity if dependend on the
        //grain size distribution
        //this is therefore needed to converge to stable values
        if( SwitchUseGrainSizeDistribution)
        {
            //if there is not an indicator for infinite sediment
            //and there is less sediment in the combined layers
            //than will be detached, use analytical solution
            if(!((Storage->Drc) < -1) &&  detachment > Storage_D.Drcd + StorageDep_D.Drcd)
            {
                double conc = (bl? UF2D_blm_D.Drcd:UF2D_ssm_D.Drcd) / ( UF2D_f->Drc);
                double tc = bl? UF2D_bltc_D.Drcd:UF2D_sstc_D.Drcd;

                double masstotal = GetTotalDW(r,c,&Storage_D)+ GetTotalDW(r,c,&StorageDep_D);
                double a = Storage_D.Drcd + StorageDep_D.Drcd;
                double c2 = masstotal;
                double b = conc;
                double d2 = DX->Drc;
                double e = ChannelAdj->Drc;
                double f = tc;
                double h = UF2D_f->Drc/(_dx*_dx);

                double t1 = c2 - b*d2*e*h + d2*e*f*h;
                double t2_1 = c2 - b*d2*e*h + d2*e*f*h;
                double t2_2 = b*c2*d2*e*h - a * d2 * e * f*h;
                double t2 = sqrt(t2_1 * t2_1 + 4 * t2_2);
                detachment = t1 - t2;

            }

        }
        //to remove small rounding errors that lead to negative values
        detachment = std::max(detachment,0.0);

        //check wat we can detache from the top and bottom layer of present material
        double dleft = detachment;
        double deptake = 0;
        double mattake = 0;
        detachment = 0;

        //when a grain size distribution is used
        if(SwitchUseGrainSizeDistribution)
        {
            //take from soil layer mass within that grain class
            deptake = std::min(dleft,StorageDep_D.Drcd);
            dleft -= deptake;
            //and take the same from the total
            StorageDep_D.Drcd -= deptake;
            StorageDep->Drc -= deptake;
        }else
        {
            //take from the total storage
            deptake = std::min(dleft,StorageDep->Drc);
            StorageDep->Drc -= deptake;
        }
        //add to the detachment what we have taken from the first soil layer
        detachment += deptake;

        //if the deposited layer is empty
        //use erosion efficiency of bottom layer again
        if(newY > 0)
        {
            dleft *= Y->Drc/newY;
        }else
        {
            dleft = 0;
        }

        if(!((Storage->Drc) < -1))
        {
            if(SwitchUseGrainSizeDistribution)
            {
                //take from soil layer mass within that grain class
                mattake = std::min(dleft,Storage_D.Drcd);
                dleft -= mattake;
                //take from the total storage
                Storage_D.Drcd -= mattake;
                Storage->Drc -= mattake;
            }else
            {
                //take from the total storage
                mattake = std::min(dleft,Storage->Drc);
                Storage->Drc -= mattake;
            }
            //add to the detachment what we have taken from the second soil layer
            detachment += mattake;
        }else
        {
            //all left potential detachment is added to detachment (infinite soil layer)
            detachment += dleft;
        }

        //finally, return the total detachment
        return std::max(0.0,detachment);

    }
}

