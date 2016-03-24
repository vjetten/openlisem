
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
void TWorld::SplashDetachment(void)
{
   if (!SwitchErosion)
      return;


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

      if(SwitchUseMaterialDepth)
      {
          if(SwitchUseGrainSizeDistribution)
          {
              StorageDep->Drc = GetTotalDW(r,c,&StorageDep_D);
              Storage->Drc = GetTotalDW(r,c,&Storage_D);
          }

          //check wat we can detache from the top and bottom layer of present material
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
      }

      // IN KG/CELL
   }
}
//---------------------------------------------------------------------------
/**
 * @fn double TWorld::OFTC(int r, int c, int d)
 * @brief Overland Flow Sediment Transport Capacity
 *
 * The overland flow transport capacity for a specified cell
 * and a specified grain class.
 * Automatically chooses the used equation based on OF_Method
 *
 * @param r : row nr of the cell
 * @param c : column nr of the cell
 * @param d : grain size class(only used with UseGrainSizeDistribution)
 * @return The transport capacity
 * @see OF_Method
 */
double TWorld::OFTC(int r, int c, int d)
{
    /*
     * If some variables are insignificant,
     * the transport capacity equations
     * predict strange values.
     * The following checks return zero when this is the case.
     * Furthermore decreases computation time before runoff starts
     */

    if(WHrunoff->Drc < MIN_HEIGHT)
    {
        return 0;
    }

    if(V->Drc < MIN_FLUX)
    {
        return 0;
    }

    if(ChannelAdj->Drc < MIN_HEIGHT)
    {
        return 0;
    }

    if(K2DSlope->Drc < MIN_HEIGHT)
    {
        return 0;
    }

    //use Govers transport capacity equation
    if(OF_Method == OFGOVERS)
    {
        CG->Drc = pow((D50->Drc+5)/0.32, -0.6);
        DG->Drc = pow((D50->Drc+5)/300, 0.25);
        //### Calc transport capacity
        double omega = 0;

            omega = 100* V->Drc*Grad->Drc;

        // V in cm/s in this formula assuming grad is SINE
        double omegacrit = 0.4;
        // critical unit streampower in cm/s
        return std::min(MAXCONC, 2650 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));
        // not more than 2650*0.32 = 848 kg/m3

    //use the Hairsine and Rose transport capacity equation
    }else if(OF_Method == OFHAIRSINEROSE)
    {
        double om =  100* V->Drc*Grad->Drc;
        double omcr = 0.4;
        double tc =  W_D.Drcd * (1.0/settlingvelocities.at(d))*(1.0 * 0.013/9.81) * (2650.0/(2650.0 - 1000.0)) * ( std::max(0.0, (om - omcr))/WHrunoff->Drc) ;
        return std::min(MAXCONC,tc);
    }

    //govers with some assumpions about distribution of stream power
    //use-able??
    /*CG->Drc = pow((graindiameters.at(d)+5)/0.32, -0.6);
    DG->Drc = pow((graindiameters.at(d)+5)/300, 0.25);
    //### Calc transport capacity
    double omega = 0;

        omega = 100* V->Drc*K2DSlope->Drc;

    // V in cm/s in this formula assuming grad is SINE
    double omegacrit = 0.4;
    // critical unit streampower in cm/s
    return std::min(MAXCONC, W_D.Drcd * 2650 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));*/

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
void TWorld::SedimentSetMaterialDistribution(int r,int c)
{
    if(SwitchUseMaterialDepth)
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
                double masstotal = GetTotalDW(r,c,&RStorage_D)+ GetTotalDW(r,c,&RStorageDep_D);
                double a = RStorage_D.Drcd + RStorageDep_D.Drcd;
                double c2 = masstotal;
                double b = bl? RBLC_D.Drcd:RSSC_D.Drcd;
                double d2 = DX->Drc;
                double e = ChannelWidth->Drc;
                double f = bl? RBLTC_D.Drcd:RSSTC_D.Drcd;
                double h = bl? RBLD_D.Drcd:RSSD_D.Drcd;
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

    }else if(flood)
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
                double masstotal = GetTotalDW(r,c,&Storage_D)+ GetTotalDW(r,c,&StorageDep_D);
                double a = Storage_D.Drcd + StorageDep_D.Drcd;
                double c2 = masstotal;
                double b = bl? BLC_D.Drcd:SSC_D.Drcd;
                double d2 = DX->Drc;
                double e = ChannelAdj->Drc;
                double f = bl? BLTC_D.Drcd:SSTC_D.Drcd;
                double h = bl? BLD_D.Drcd:SSD_D.Drcd;

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

    //if it is neither flood nor channel detachment, overland flow is assumed
    }else
    {
        //calculate depth of deposited layer
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
            if(!((Storage->Drc) < -1) && detachment > (Storage_D.Drcd + StorageDep_D.Drcd))
            {
                //if there is not an indicator for infinite sediment
                //and there is less sediment in the combined layers
                //than will be detached, use analytical solution
                double masstotal = GetTotalDW(r,c,&Storage_D)+ GetTotalDW(r,c,&StorageDep_D);
                double a = Storage_D.Drcd + StorageDep_D.Drcd;
                double c2 = masstotal;
                double b = Conc_D.Drcd;
                double d2 = DX->Drc;
                double e = ChannelAdj->Drc;
                double f = TC_D.Drcd;
                double h = WHrunoff->Drc;

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
        //bottom soil layer can be infinite
        if(!((Storage->Drc) < -1))
        {
            //take sediment from the soil layer storage
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

//---------------------------------------------------------------------------
/**
 * @fn void TWorld::FlowDetachment(void)
 * @brief Calculates flow detachment for overland flow in entire catchment
 *
 * This function uses the function for overland flow transport capacity to
 * calculate the potential detachment/deposition based on the settling velocity of the sediment.
 * When potential detachment is found, the fuction for taking soil
 * from the soil layer is used to find actual detachment.
 * When deposition is found, this sediment is added to the deposited soil layer.
 *
 * @see TWorld:OFTC
 * @see TWorld:GetSV
 * @see TWorld:DetachMaterial
 *
 */
void TWorld::FlowDetachment(void)
{
   if (!SwitchErosion)
      return;

   //ldd directions, all adjecent cells and the cell itself at index 5
   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, -1, -1, -1, 0, 0, 0, 1, 1, 1};

   //transport capacity
   FOR_ROW_COL_MV
   {

        if(!SwitchUseGrainSizeDistribution)
        {
            //get the transport capacity for a single grain size
            TC->Drc = OFTC(r,c,-1);

        }else
        {
            //get the transport capacity for all induvidual grain sizes
            TC->Drc = 0;
            FOR_GRAIN_CLASSES
            {
                TC_D.Drcd = OFTC(r,c,d);
                TC->Drc += TC_D.Drcd;
            }

            //if the total transport capacity is too high, limit it to the maximum value (MAXCONC)
            if(TC->Drc > MAXCONC)
            {
                //rescale all induvidual tc's to conform to MAXCONC
                FOR_GRAIN_CLASSES
                {
                    TC_D.Drcd *= MAXCONC/TC->Drc;
                }
                TC->Drc = MAXCONC;
            }

            TC->Drc = 0;
            FOR_GRAIN_CLASSES
            {
                TC->Drc += TC_D.Drcd;
            }

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

         //check all surrounding cells for MV or present sediment
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

         //limit sediment current cell
         double tcold = TC->Drc;

         if(count > 0)
         {
             TC->Drc = std::max(TC->Drc, avgtc/count);
         }
         TC->Drc = std::min(TC->Drc, maxtc);

         //when a grain size distribution is used, scale ale induvidual transport capacity's
         if(SwitchUseGrainSizeDistribution)
         {
             double fact = TC->Drc/tcold;
             if(tcold != 0)
             {
                 FOR_GRAIN_CLASSES
                 {
                     TC_D.Drcd *= fact;
                 }
             }
         }
      }
   }

   //the iterator is either the number of grain classes, or 1 if no grain size distribution is used.
   int iterator = numgrainclasses;
   if(!SwitchUseGrainSizeDistribution)
   {

        iterator = 1;
   }

   FOR_ROW_COL_MV
   {
       DETFlow->Drc = 0;
       DEP->Drc = 0;
   }
   //for each grain class, calculate flow detachment seperately
   for(int d  = 0 ; d < iterator;d++)
   {

       FOR_ROW_COL_MV
       {
          double erosionwh = WHrunoff->Drc;
          double erosionwv = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc;
          if(SwitchKinematic2D > 1)
          {
              erosionwh = std::max(WHrunoff->Drc-K2DWHStore->Drc ,0.0);
              erosionwv = std::max(WHrunoff->Drc-K2DWHStore->Drc ,0.0)*ChannelAdj->Drc*DX->Drc;
          }
          //assume splash detachment has same grain size distribution as soil
          if(!SwitchUseGrainSizeDistribution)
          {
              Sed->Drc += DETSplash->Drc;

          }else
          {
              Sed_D.Drcd += DETSplash->Drc * W_D.Drcd;
          }

          if(!SwitchUseGrainSizeDistribution)
          {
              // init deposition for this timestep
              Conc->Drc = MaxConcentration(erosionwv, Sed->Drc);
              // limit sed concentration to max
          }else
          {
               Conc_D.Drcd = MaxConcentration(erosionwv, Sed_D.Drcd);

          }
          double maxTC = 0;
          double minTC = 0;

          if(!SwitchUseGrainSizeDistribution)
          {
              maxTC = std::max(TC->Drc - Conc->Drc,0.0);
              // positive difference: TC deficit becomes detachment (positive)
              minTC = std::min(TC->Drc - Conc->Drc,0.0);
              // negative difference: TC surplus becomes deposition (negative)
              // unit kg/m3
          }else
          {
              maxTC = std::max(TC_D.Drcd - Conc_D.Drcd,0.0);
              // positive difference: TC deficit becomes detachment (positive)
              minTC = std::min(TC_D.Drcd - Conc_D.Drcd,0.0);
              // negative difference: TC surplus becomes deposition (negative)
              // unit kg/m3
          }

          //### detachment
          double TransportFactor = 0;
          if(!SwitchUseGrainSizeDistribution)
          {
              TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
          }else
          {
              TransportFactor = _dt*settlingvelocities.at(d) * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
          }

          // detachment can only come from soil, not roads (so do not use flowwidth)
          // units s * m/s * m * m = m3

          double detachment = maxTC * TransportFactor;

          if(SwitchUseGrainSizeDistribution)
          {
              //detachment =  W_D.at(d)->Drc * detachment;
          }
          // unit = kg/m3 * m3 = kg
          detachment = std::min(detachment, maxTC * erosionwv);

          // cannot have more detachment than remaining capacity in flow
          // use discharge because standing water has no erosion

          if (GrassFraction->Drc > 0)
             detachment = (1-GrassFraction->Drc) * detachment;
          // no flow detachment on grass strips

          // Detachment edxceptions:
          detachment = (1-StoneFraction->Drc) * detachment ;
          // no flow detachment on stony surfaces

          if (SwitchHardsurface)
             detachment = (1-HardSurface->Drc) * detachment;
          // no flow detachment on hard surfaces

          if (SwitchHouses)
             detachment = (1-HouseCover->Drc)*detachment;
          // no flow det from house roofs

          detachment = DetachMaterial(r,c,d,false,false,false, detachment);

          //### sediment balance
          Sed->Drc += detachment;
          DETFlow->Drc += detachment;
          if(SwitchUseGrainSizeDistribution)
          {
            Sed_D.Drcd += detachment;
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
                  TransportFactor = (1-exp(-_dt*settlingvelocities.at(d)/erosionwh)) * erosionwv;
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
          if(!SwitchUseGrainSizeDistribution)
          {
              deposition = std::max(deposition, -Sed->Drc);
          }else
          {
              deposition = std::max(deposition, -Sed_D.Drcd);
          }
          //force deposition on grass strips
          if (GrassFraction->Drc > 0)
          {
              if(!SwitchUseGrainSizeDistribution)
              {
                  deposition = -Sed->Drc*GrassFraction->Drc + (1-GrassFraction->Drc)*deposition;
              }else
              {
                  deposition = -Sed_D.Drcd*GrassFraction->Drc + (1-GrassFraction->Drc)*deposition;
              }
          }
          //add deposition to soil layer
          if(SwitchUseMaterialDepth)
          {
              StorageDep->Drc += -deposition;
              if(SwitchUseGrainSizeDistribution)
              {
                    StorageDep_D.Drcd += -deposition;
              }
          }

          //keep track of deposition
          DEP->Drc += deposition;

          // add to sediment in flow (IN KG/CELL)
          Sed->Drc += deposition;
          //if needed, add to sediment in flow for a certain grain class
          if(SwitchUseGrainSizeDistribution)
          {
            Sed_D.Drcd += deposition;
          }
          //calculate concentration
          if(SwitchUseGrainSizeDistribution)
          {
              Conc_D.Drcd = MaxConcentration(erosionwv, Sed_D.Drcd);
          }

          Conc->Drc = MaxConcentration(erosionwv, Sed->Drc);
       }
   }
}
//---------------------------------------------------------------------------
/**
 * @fn void TWorld::ChannelFlowDetachment(int r, int c)
 * @brief Calculates flow detachment for channel flow in a specific cell
 *
 * This function uses the function for channel flow transport capacity to
 * calculate the potential detachment/deposition based on the settling velocity of the sediment.
 * This process is done for one or two transport layer (bed/suspended sediment load)
 * When potential detachment is found, the fuction for taking soil
 * from the soil layer is used to find actual detachment.
 * When deposition is found, this sediment is added to the deposited soil layer.
 *
 * @param r : the row nr of the cell
 * @param c : the column nr of the cell
 * @see TWorld:RiverSedimentTCBL
 * @see TWorld:RiverSedimentTCSS
 * @see TWorld:SwitchUse2Layer
 * @see TWorld:DetachMaterial
 *
 */
void TWorld::ChannelFlowDetachment(int r, int c)
{
   if (!SwitchErosion)
      return;

   //ldd directions, all adjecent cells and the cell itself at index 5
   int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
   int dy[10] = {0, -1, -1, -1, 0, 0, 0, 1, 1, 1};

   //calculate layer depth for this cell

   /* Dependent on SwitchUse2Layer,
    * a 2 layer system is modelled for Channels.
    * This is implemented through the RiverSedimentLayerDepth() function.
    * When SwitchUse2Layer is false, the bed load layer depth
    * is equal to the full water depth!
    * No suspended transport is then modelled.
    */

   RiverSedimentLayerDepth(r,c);

   //iterator is number of grain classes
   int iterator = numgrainclasses;
   if(!SwitchUseGrainSizeDistribution)
   {

        iterator = 1;
   }

   ChannelDetFlow->Drc = 0;
   ChannelDep->Drc = 0;

   //find transport capacity for bed and suspended layer
   for(int d  = 0 ; d < iterator;d++)
   {

       cTMap * TBLDepthFlood;
       cTMap * TSSDepthFlood;
       cTMap * TBLTCFlood;
       cTMap * TSSTCFlood;
       cTMap * TBLCFlood;
       cTMap * TSSCFlood;
       cTMap * TBLFlood;
       cTMap * TSSFlood;
       cTMap * TW;

       double TSettlingVelocity;


       if(!SwitchUseGrainSizeDistribution)
       {
           TBLDepthFlood = ChannelBLDepth;
           TSSDepthFlood = ChannelSSDepth;
           TBLTCFlood = ChannelBLTC;
           TSSTCFlood = ChannelSSTC;
           TBLCFlood = ChannelBLConc;
           TSSCFlood = ChannelSSConc;
           TBLFlood = ChannelBLSed;
           TSSFlood = ChannelSSSed;
           TW = unity;
           TSettlingVelocity = SettlingVelocity->Drc;
       }else
       {
           TBLDepthFlood = RBLD_D.at(d);
           TSSDepthFlood = RSSD_D.at(d);
           TBLTCFlood = RBLTC_D.at(d);
           TSSTCFlood = RSSTC_D.at(d);
           TBLCFlood = RBLC_D.at(d);
           TSSCFlood = RSSC_D.at(d);
           TBLFlood = RBL_D.at(d);
           TSSFlood = RSS_D.at(d);
           TW = RW_D.at(d);
           TSettlingVelocity = settlingvelocities.at(d);
       }

       //get transport capacit for bed/suspended load for a specific cell and grain size class
       TBLTCFlood->Drc = RiverSedimentTCBL(r,c,d);
       TSSTCFlood->Drc = RiverSedimentTCSS(r,c,d);

    }


   //check if the sum of transport capacities of all grain sizes is larger than MAXCONC, and rescale if nessecery
   if(SwitchUseGrainSizeDistribution)
   {
        ChannelBLTC->Drc = 0;
        ChannelSSTC->Drc = 0;
        FOR_GRAIN_CLASSES
        {
            ChannelBLTC->Drc += RBLTC_D.Drcd;
            ChannelSSTC->Drc += RSSTC_D.Drcd;
        }


        if(ChannelBLTC->Drc > MAXCONCBL)
       {
           FOR_GRAIN_CLASSES
           {
               RBLTC_D.Drcd *= MAXCONCBL/ChannelBLTC->Drc;
           }
           ChannelBLTC->Drc = MAXCONCBL;
       }
       if(ChannelSSTC->Drc > MAXCONC)
       {
           FOR_GRAIN_CLASSES
           {
               RSSTC_D.Drcd *= MAXCONC/ChannelSSTC->Drc;
           }
           ChannelSSTC->Drc = MAXCONC;
       }

       ChannelBLTC->Drc = 0;
       ChannelSSTC->Drc = 0;
       FOR_GRAIN_CLASSES
       {
           ChannelBLTC->Drc += RBLTC_D.Drcd;
           ChannelSSTC->Drc += RSSTC_D.Drcd;
       }

   }

   for(int d  = 0 ; d < iterator;d++)
   {
       //set all maps for this grain class

       cTMap * TBLDepthFlood;
       cTMap * TSSDepthFlood;
       cTMap * TBLTCFlood;
       cTMap * TSSTCFlood;
       cTMap * TBLCFlood;
       cTMap * TSSCFlood;
       cTMap * TBLFlood;
       cTMap * TSSFlood;
       cTMap * TW;

       double TSettlingVelocity;

       if(!SwitchUseGrainSizeDistribution)
       {
           TBLDepthFlood = ChannelBLDepth;
           TSSDepthFlood = ChannelSSDepth;
           TBLTCFlood = ChannelBLTC;
           TSSTCFlood = ChannelSSTC;
           TBLCFlood = ChannelBLConc;
           TSSCFlood = ChannelSSConc;
           TBLFlood = ChannelBLSed;
           TSSFlood = ChannelSSSed;
           TW = unity;
           TSettlingVelocity = SettlingVelocity->Drc;
       }else
       {
           TBLDepthFlood = RBLD_D.at(d);
           TSSDepthFlood = RSSD_D.at(d);
           TBLTCFlood = RBLTC_D.at(d);
           TSSTCFlood = RSSTC_D.at(d);
           TBLCFlood = RBLC_D.at(d);
           TSSCFlood = RSSC_D.at(d);
           TBLFlood = RBL_D.at(d);
           TSSFlood = RSS_D.at(d);
           TW = RW_D.at(d);
           TSettlingVelocity = settlingvelocities.at(d);
       }


       double velocity = ChannelV->Drc;

       double bldepth = TBLDepthFlood->Drc;
       double ssdepth = TSSDepthFlood->Drc;

       //discharges for both layers and watervolumes
       double bldischarge = velocity * ChannelWidth->Drc * bldepth;
       double blwatervol = ChannelWidth->Drc *DX->Drc*bldepth;

       double ssdischarge = velocity * ChannelWidth->Drc * ssdepth;
       double sswatervol = ChannelWidth->Drc *DX->Drc * ssdepth;

       double bltc = 0;
       double sstc = 0;

       bltc = TBLTCFlood->Drc;
       sstc = TSSTCFlood->Drc;

       //when waterheight is insignificant, deposite all remaining sediment
       double deposition;
       if(bldepth < he_ca)
       {
           bltc = 0;
       }
       if(ssdepth < he_ca)
       {
           sstc = 0;
       }
       if(ChannelWH->Drc == 0)
       {
           deposition = -TBLFlood->Drc;
           deposition += -TSSFlood->Drc;
           ChannelDep->Drc += deposition;

           if(SwitchUseMaterialDepth)
           {
               RStorageDep->Drc += -deposition;
               if(SwitchUseGrainSizeDistribution)
               {
                     RStorageDep_D.Drcd += -deposition;
                     if(std::isnan(RStorageDep_D.Drcd))
                     {
                         qDebug() << "NAN dep1" << d;
                     }
               }
           }

           TBLDepthFlood = 0;
           TSSDepthFlood = 0;
           TBLTCFlood = 0;
           TSSTCFlood = 0;
           TBLCFlood = 0;
           TSSCFlood = 0;
           TBLFlood = 0;
           TSSFlood = 0;

           if(SwitchUseGrainSizeDistribution)
           {
                   RBL_D.Drcd = 0;
                   RSS_D.Drcd = 0;
                   RBLTC_D.Drcd = 0;
                   RSSTC_D.Drcd = 0;
                   RBLC_D.Drcd = 0;
                   RSSC_D.Drcd = 0;
           }
       }else
       {

           //first check if sediment goes to suspended sediment layer or to bed layer
          double tobl = 0;
          double toss = 0;
          double TransportFactor;
          if (TSSDepthFlood->Drc > MIN_HEIGHT)
          {
               TransportFactor = (1-exp(-_dt*TSettlingVelocity/TSSDepthFlood->Drc)) * sswatervol;

          }else
          {
               TransportFactor =  1* sswatervol;
          }
          double maxTC = std::max(sstc - TSSCFlood->Drc,0.0) ;
          // positive difference: TC deficit becomes detachment (ppositive)
          double minTC = std::min(sstc - TSSCFlood->Drc,0.0) ;

          tobl = TransportFactor * minTC;

          tobl = std::max(tobl,-TSSFlood->Drc);
          TBLFlood->Drc -= tobl;
          TSSFlood->Drc += tobl;

          TransportFactor = _dt*TSettlingVelocity * DX->Drc * ChannelWidthUpDX->Drc;


          double detachment = TW->Drc * maxTC * TransportFactor;
          if (GrassFraction->Drc > 0)
             detachment = (1-GrassFraction->Drc) * detachment;
          detachment = (1-StoneFraction->Drc) * detachment ;
          if (SwitchHardsurface)
             detachment = (1-HardSurface->Drc) * detachment;
          if (SwitchHouses)
             detachment = (1-HouseCover->Drc)* detachment;


          detachment = DetachMaterial(r,c,d,true,false,false, detachment);

          //### sediment balance
          TSSFlood->Drc += detachment;
          ChannelDetFlow->Drc += detachment;


          double sssmax = MAXCONC * DX->Drc *ChannelWidth->Drc*ssdepth;
          if(sssmax < TSSFlood->Drc)
          {
              TBLFlood->Drc += (TSSFlood->Drc - sssmax);
              TSSFlood->Drc = sssmax;
          }


          TBLCFlood->Drc = MaxConcentration(blwatervol, TBLFlood->Drc);
          // limit concentration to 850 and throw rest in deposition
          TSSCFlood->Drc = MaxConcentration(sswatervol, TSSFlood->Drc);
          // limit concentration to 850 and throw rest in deposition

          //deposition and detachment
          //### calc concentration and net transport capacity

          maxTC = std::max(TBLTCFlood->Drc - TBLCFlood->Drc,0.0);
          // positive difference: TC deficit becomes detachment (ppositive)
          minTC = std::min(TBLTCFlood->Drc - TBLCFlood->Drc,0.0);
          // negative difference: TC surplus becomes deposition (negative)
          // unit kg/m3

          //### detachment
          TransportFactor = _dt*TSettlingVelocity * DX->Drc *ChannelWidthUpDX->Drc;
          // detachment can only come from soil, not roads (so do not use flowwidth)
          // units s * m/s * m * m = m3

          detachment =  TW->Drc *  maxTC * TransportFactor;
          // unit = kg/m3 * m3 = kg
          //BLDetFlood->Drc = std::min(BLDetFlood->Drc, maxTC * bldischarge*dt);
          // cannot have more detachment than remaining capacity in flow
          // use discharge because standing water has no erosion

          if (GrassFraction->Drc > 0)
             detachment = (1-GrassFraction->Drc) * detachment;
          // no flow detachment on grass strips

          // Detachment edxceptions:
          detachment = (1-StoneFraction->Drc) * detachment;
          // no flow detachment on stony surfaces

          if (SwitchHardsurface)
             detachment = (1-HardSurface->Drc) * detachment;
          // no flow detachment on hard surfaces

          if (SwitchHouses)
             detachment = (1-HouseCover->Drc)*detachment;
          // no flow det from house roofs

          detachment = DetachMaterial(r,c,d,true,false,true, detachment);

          // IN KG/CELL

          //DETFlow->Drc = (1-Snowcover->Drc) * DETFlow->Drc ;
          /* TODO: CHECK THIS no flow detachment on snow */
          //is there erosion and sedimentation under the snowdeck?

          //### deposition
          if (TBLDepthFlood->Drc > MIN_HEIGHT)
             TransportFactor = (1-exp(-_dt*TSettlingVelocity/bldepth)) * blwatervol;
          else
             TransportFactor = 1*blwatervol;
          // if settl velo is very small, transportfactor is 0 and depo is 0
          // if settl velo is very large, transportfactor is 1 and depo is max

          //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
          // deposition can occur on roads and on soil (so use flowwidth)

          deposition = minTC * TransportFactor;
          // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

          if (SwitchLimitDepTC)
             deposition = std::max(deposition, minTC *blwatervol);
          // cannot be more than sediment above capacity
          deposition = std::max(deposition, -TBLFlood->Drc);
          // cannot have more depo than sediment present

          if(SwitchUseMaterialDepth)
          {
              RStorageDep->Drc += -deposition;
              if(SwitchUseGrainSizeDistribution)
              {
                    RStorageDep_D.Drcd += -deposition;
              }
          }

          ChannelDep->Drc += deposition;
          ChannelDetFlow->Drc += detachment;
          TBLFlood->Drc += detachment;
          TBLFlood->Drc += deposition;

       }
   }

   if(SwitchUseGrainSizeDistribution)
   {
       ChannelBLSed->Drc = 0;
       ChannelSSSed->Drc = 0;
       ChannelBLTC->Drc = 0;
       ChannelSSTC->Drc = 0;
       ChannelBLConc->Drc = 0;
       ChannelSSConc->Drc = 0;

       FOR_GRAIN_CLASSES
       {
           ChannelBLSed->Drc += RBL_D.Drcd;
           ChannelSSSed->Drc += RSS_D.Drcd;
           ChannelBLTC->Drc += RBLTC_D.Drcd;
           ChannelSSTC->Drc += RSSTC_D.Drcd;
           ChannelBLConc->Drc += RBLC_D.Drcd;
           ChannelSSConc->Drc += RSSC_D.Drcd;
       }
   }

   //total transport capacity (bed load + suspended load)
   ChannelTC->Drc = ChannelBLTC->Drc + ChannelSSTC->Drc;

   //calculate concentrations
   ChannelBLConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelBLSed->Drc);
   ChannelSSConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSSSed->Drc);

   //total concentration
   ChannelConc->Drc = ChannelBLConc->Drc + ChannelSSConc->Drc;

   //check for concentration surpassing MAXCONC
   RiverSedimentMaxC(r,c);
}
//---------------------------------------------------------------------------
/**
 * @fn void TWorld::RiverSedimentMaxC(int r, int c)
 * @brief Limits sediment concentration to a maximum possible concentration
 *
 * Limits sediment concentration to a maximum possible concentration.
 * When a grain size distribution is used, seperate concentrations are scaled.
 * All surpassing sediment is deposited.
 *
 * @param r : the row nr of the cell
 * @param c : the column nr of the cell
 * @see MAXCONC
 *
 */
void TWorld::RiverSedimentMaxC(int r, int c)
{

    cTMap * _BL = ChannelBLSed;
    cTMap * _BLC = ChannelBLConc;
    cTMap * _SS = ChannelSSSed;
    cTMap * _SSC = ChannelSSConc;


    //maximum concentraion
    if(!SwitchUseGrainSizeDistribution)
    {

        _SSC->Drc = MaxConcentration(ChannelWidth->Drc*DX->Drc*ChannelSSDepth->Drc, _SS->Drc);
        // limit concentration to 850 and throw rest in deposition

        double sssmax = MAXCONC * DX->Drc *ChannelWidth->Drc*ChannelSSDepth->Drc;
        if(sssmax < _SS->Drc)
        {
            ChannelDep->Drc += -(_SS->Drc - sssmax);
            if(SwitchUseMaterialDepth)
            {
                RStorageDep->Drc += (_SS->Drc - sssmax);
            }
            _SS->Drc = sssmax;

        }


        //set concentration from present sediment
        _BLC->Drc = MaxConcentration(ChannelWidth->Drc*DX->Drc*ChannelBLDepth->Drc, _BL->Drc);

        double smax = MAXCONCBL * DX->Drc *ChannelWidth->Drc*ChannelBLDepth->Drc;
        if(smax < _BL->Drc)
        {
            ChannelDep->Drc += -(_BL->Drc - smax);
            if(SwitchUseMaterialDepth)
            {
                RStorageDep->Drc += (_BL->Drc - smax);
            }
            _BL->Drc = smax;

        }
    }else
    {
       FOR_GRAIN_CLASSES
       {
            RSSC_D.Drcd = MaxConcentration(ChannelWidth->Drc*DX->Drc*RSSD_D.Drcd, RSS_D.Drcd);
            // limit concentration to 850 and throw rest in deposition

            double sssmax = MAXCONC * DX->Drc *ChannelWidth->Drc*RSSD_D.Drcd;
            if(sssmax < RSS_D.Drcd)
            {
                ChannelDep->Drc += -(RSS_D.Drcd - sssmax);
                ChannelSSSed->Drc += -(RSS_D.Drcd - sssmax);
                if(SwitchUseMaterialDepth)
                {
                    RStorageDep->Drc += (RSS_D.Drcd - sssmax);
                    RStorageDep_D.Drcd += (RSS_D.Drcd - sssmax);
                    if(std::isnan(RStorageDep_D.Drcd))
                    {
                        qDebug() << "NAN dep3" << d;
                    }
                }
                RSS_D.Drcd = sssmax;

            }


            RBLC_D.Drcd = MaxConcentration(ChannelWidth->Drc*DX->Drc*RBLD_D.Drcd, RBL_D.Drcd);
            // limit concentration to 850 and throw rest in deposition

            sssmax = MAXCONCBL * DX->Drc *ChannelWidth->Drc*RBLD_D.Drcd;
            if(sssmax < BL_D.Drcd)
            {
                ChannelDep->Drc += -(RBL_D.Drcd - sssmax);
                ChannelBLSed->Drc += -(RBL_D.Drcd - sssmax);
                RBL_D.Drcd = sssmax;
                if(SwitchUseMaterialDepth)
                {
                    RStorageDep->Drc += (RBL_D.Drcd - sssmax);
                    RStorageDep_D.Drcd += (RBL_D.Drcd - sssmax);
                    if(std::isnan(RStorageDep_D.Drcd))
                    {
                        qDebug() << "NAN dep4" << d;
                    }
                }
            }
       }

       ChannelBLConc->Drc = 0;
       ChannelSSConc->Drc = 0;

       FOR_GRAIN_CLASSES
       {
           ChannelBLConc->Drc += RBLC_D.Drcd;
           ChannelSSConc->Drc += RSSC_D.Drcd;
       }



    }

    ChannelSed->Drc = ChannelBLSed->Drc + ChannelSSSed->Drc;
    ChannelConc->Drc = ChannelBLConc->Drc + ChannelSSConc->Drc;


}

//---------------------------------------------------------------------------
/**
 * @fn void TWorld::RiverSedimentDiffusion(double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
 * @brief Diffusion throughout the channels
 *
 * This function diffuses a material map based on a concentration map
 * for a timestep dt.
 * The diffusion is scaled according to the turbulent Prandtl-Smidth number.
 * Note that the _BL and _BLC are not used since there is no diffusion in
 * a bed load layer.
 *
 * @param dt : the timestep taken with this diffusion
 * @param _BL : Bed load material to be diffused
 * @param _BLC : Bed load material concentration
 * @param _SS : Suspended material to be diffused
 * @param _SSC : Suspended material concentration
 *
 * @see FS_SigmaDiffusion
 *
 */
void TWorld::RiverSedimentDiffusion(double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
{

    FOR_ROW_COL_MV
    {

        MBLNFlood->Drc = _BL->Drc;
        MSSNFlood->Drc = _SS->Drc;
        MBLFlood->Drc = _BL->Drc;
        MSSFlood->Drc = _SS->Drc;

        //set concentration from present sediment
        MBLCFlood->Drc = MaxConcentration(ChannelWaterVol->Drc, MBLFlood->Drc);

        //set concentration from present sediment
        MSSCFlood->Drc = MaxConcentration(ChannelWaterVol->Drc, MSSFlood->Drc);

    }


    //diffusion of Suspended Sediment layer
    FOR_ROW_COL_MV_CH
    {

        int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
        int dy[10] = {0, -1, -1, -1, 0, 0, 0, 1, 1, 1};

        int rp = r;
        int cp = c;

        bool foundp = false;
        /** put all points that have to be calculated to calculate the current point in the list,
         before the current point */
        for (int i=1; i<=9; i++)
        {
            int rt = 0, ct = 0;
            int ldd = 0;

            // this is the current cell
            if (i==5)
                continue;

            rt = r+dy[i];
            ct = c+dx[i];

            if (INSIDE(rt, ct) && !pcr::isMV(LDDChannel->data[rt][ct]))
                ldd = (int) LDDChannel->data[rt][ct];
            else
                continue;

            // check if there are more cells upstream, if not subCatchDone remains true
            if (pcr::isMV(ChannelQn->Drc) && INSIDE(rt, ct) &&
                    FLOWS_TO(ldd, rt, ct, r, c)
                    )
            {
                rp = rt;
                cp = ct;
                foundp = true;
                break;

            }
        }

        bool foundn = false;
        int rn = 0, cn = 0;
        int ldd = (int) LDDChannel->data[r][c];
        if(ldd == 5)
        {
            foundn = false;

        }else if (pcr::isMV(ChannelQn->Drc) &&
             INSIDE(r+dy[ldd], c+dx[ldd]))
        {
            foundn = true;
            rn = r+dy[ldd];
            cn = c+dx[ldd];
        }

        //cell sizes
        double cdx = _dx;
        //here it is about spacing, not flow width, so use _dx instead of CHannelAdj->Drc

        //mixing coefficient
        double dux1 = 0;
        if(foundn)
        {
           dux1 = std::abs(ChannelV->data[r][c] - ChannelV->data[rp][cp]);
        }
        double dux2 = 0;
        if(foundn)
        {
            dux2 = std::abs(ChannelV->data[r][c] - ChannelV->data[rn][cn]);
        }

        double dux = std::max(dux1,dux2);

        //diffusion coefficient according to J.Smagorinski (1964)
        double eddyvs = cdx * dux;
        //and devide by turbulent prandtl-smidth number
        double eta = eddyvs/FS_SigmaDiffusion;

        //add diffusive fluxes to previous cell in channel.
        if(foundp)
        {
            double coeff1 = std::min(dt*eta *std::min(1.0,ChannelSSDepth->data[rp][cp]/ChannelSSDepth->data[r][c]),courant_factor_diffusive/2.0) * MSSFlood->Drc;

            MSSNFlood->data[rp][cp] += coeff1;
            MSSNFlood->data[r][c] -= coeff1;
        }

        //add diffusive fluxes to next cell in channel.
        if(foundn)
        {
            double coeff2 = std::min(dt*eta *std::min(1.0,ChannelSSDepth->data[rn][cn]/ChannelSSDepth->data[r][c]),courant_factor_diffusive/2.0) * MSSFlood->Drc;

            MSSNFlood->data[rn][cn] += coeff2;
            MSSNFlood->data[r][c] -= coeff2;
        }
    }

    //recalculate concentrations
    FOR_ROW_COL_MV
    {
        _BL->Drc = MBLNFlood->Drc;
        _SS->Drc = MSSNFlood->Drc;

        //set concentration from present sediment
        _BLC->Drc = MaxConcentration(ChannelWaterVol->Drc, _BL->Drc);

        //set concentration from present sediment
        _SSC->Drc = MaxConcentration(ChannelWaterVol->Drc, _SS->Drc);
    }
}
//---------------------------------------------------------------------------
/**
 * @fn void TWorld::RiverSedimentLayerDepth(int r , int c)
 * @brief Calculates River bed load layer depth
 *
 *
 * @param r : the timestep taken with this diffusion
 * @param c : Bed load material to be diffused
 */
void TWorld::RiverSedimentLayerDepth(int r , int c)
{
    if(SwitchErosion)
    {
        ChannelBLDepth->Drc = ChannelWH->Drc;
        ChannelSSDepth->Drc = 0;
        return;
    }

    if(!SwitchUseGrainSizeDistribution)
    {
        //if a two layer system is modelled, calculate thickness of layer
        if(SwitchUse2Layer)
        {

            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            double ps = 2400;
            double pw = 1000;
            double velocity = ChannelV->Drc;

            //critical shear velocity for bed level motion by van rijn
            double critshearvel = velocity * sqrt(GRAV)/(18 * log10(4*(ChannelWidth->Drc * ChannelWH->Drc/(ChannelWH->Drc * 2 + ChannelWidth->Drc))/d90m));
            //critical shear stress for bed level motion by van rijn
            double critsheart = (critshearvel*critshearvel)/ (((ps-pw)/pw) * GRAV*d50m);
            //rough bed bed load layer depth by Hu en Hui
            ChannelBLDepth->Drc = std::min(std::min(d50m * 1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), ChannelWH->Drc), 0.1);
            ChannelSSDepth->Drc = std::max(ChannelWH->Drc - ChannelBLDepth->Drc,0.0);

        //if no two layer system is used
        }else
        {   //full water height to bed layer depth
            //note: this is autimatically combined with Govers transport Capacity for this layer.
            ChannelBLDepth->Drc = ChannelWH->Drc;
            ChannelSSDepth->Drc = 0;
        }
    }else
    {
        ChannelBLDepth->Drc = 0;
        ChannelSSDepth->Drc = 0;

        FOR_GRAIN_CLASSES
        {
            //for each grain size, use the grain class property
            double d50m = graindiameters.at(d)/1000000.0;
            double d90m = 1.5 * graindiameters.at(d)/1000000.0;

            double ps = 2400;
            double pw = 1000;
            double velocity = ChannelV->Drc;

            //critical shear velocity for bed level motion by van rijn
            double critshearvel = velocity * sqrt(GRAV)/(18 * log10(4*(ChannelWidth->Drc * ChannelWH->Drc/(ChannelWH->Drc * 2 + ChannelWidth->Drc))/(d90m)));
            //critical shear stress for bed level motion by van rijn
            double critsheart = (critshearvel*critshearvel)/ (((ps-pw)/pw) * GRAV*d50m);
            //rough bed bed load layer depth by Hu en Hui
            RBLD_D.Drcd = std::min(std::min(d50m * 1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), ChannelWH->Drc), 0.1);
            RSSD_D.Drcd = std::max(ChannelWH->Drc - RBLD_D.Drcd ,0.0);
            ChannelBLDepth->Drc += RBLD_D.Drcd * RW_D.Drcd;
            ChannelSSDepth->Drc += RSSD_D.Drcd * RW_D.Drcd;

        }
    }
}
//---------------------------------------------------------------------------
/**
 * @fn void TWorld::RiverSedimentLayerDepth(int r , int c)
 * @brief Calculates River bed load transport capacity
 *
 * Calculates river bed load sediment transport capacity.
 * Based on either Govers (for single layer transport),
 * Van Rijn (simplifie or full) or Wu wang & Jia (for multiclass sediment).
 *
 * @param r : The row nr of the cell
 * @param c : The column nr of the cell
 * @param _d : The grain class (only needed when grain size distribution is used)
 * @see R_BL_Method
 */
double TWorld::RiverSedimentTCBL(int r,int c,int _d)
{
    //when water height is insignificant, transport capacity is zero
    //this is necessary since some of the used equations have strange behaviour
    //for these water heights or velocities. (h and v outside of valid range)
    if(ChannelWH->Drc < MIN_HEIGHT || ChannelBLDepth->Drc < MIN_HEIGHT)
    {
        return 0;
    }
    double v = ChannelV->Drc;
    if(v < MIN_FLUX)
    {
        return 0;
    }
    if(ChannelWidth->Drc < MIN_FLUX)
    {
        return 0;
    }

        if(R_BL_Method == FSGOVERS)
            {

            //Govers with a maximum bed load layer depth (1980)
            double discharge = ChannelQ->Drc;
            //### Calc transport capacity
            double omega = 100.0*v*discharge;
            // V in cm/s in this formula assuming grad is SINE
            double omegacrit = 0.4;
            // critical unit streampower in cm/s
            return std::min(MAXCONCBL, 2650 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));
            // not more than 2650*0.32 = 848 kg/m3

        }else if(R_BL_Method == FSRIJN)
        {
            //Van rijn simplified (1984)

            double ps = 2400.0;
            double pw = 1000.0;
            double ucr;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            if( d50m < 0.005)
            {
               ucr  = 0.19 * pow(d50m, 0.1) * log10(4.0* ChannelWH->Drc/d90m);
            }else
            {
               ucr  = 0.19 * pow(d50m, 0.6) * log10(4.0* ChannelWH->Drc/d90m);
            }
            double me = std::max((v - ucr)/(sqrt(GRAV * d50m * ((ps/pw) - 1.0))),0.0);
            double qs = 0.005 * ps*v *ChannelWH->Drc * pow(d50m/ChannelWH->Drc,1.2) * pow(me, 2.4);
            double tc =  qs/ (v * ChannelBLDepth->Drc );
            return std::max(std::min( tc,MAXCONCBL ),0.0);
        }else if(R_BL_Method == FSRIJNFULL)
        {
            //van Rijn full (1980)

            double ps = 2400.0;
            double pw = 1000.0;
            double kinvis = 1.0;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);

            double ds = D50->Drc * pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0));
            double dh = ChannelWH->Drc;
            double chevey = 18 * log(4 * dh/d90m);
            double us = sqrt(GRAV) * v/chevey;
            double uscr = 0.055;

            if(ds <150 && !(ds < 20))
                uscr = 0.013*pow(ds,0.29);
            if(ds < 20 && !(ds < 10))
                uscr = 0.04*pow(ds,-0.10);
            if(ds <10 && !(ds < 4))
                uscr = 0.14*pow(ds,-0.64);
            if(ds <4)
                uscr = 0.24*pow(ds,-1);

            uscr = sqrt(uscr * (ps/pw - 1)*GRAV * d50m);
            double T = std::max((us*us/(uscr*uscr)) - 1,0.0);
            double qs = 0.053 * (pow(T,2.1)/pow(ds,0.3)) * sqrt((ps/pw -1)*GRAV)*d50m*sqrt(d50m);
            double tc = 0.1 *  ps * ChannelWidth->Drc * qs/ (v * ChannelBLDepth->Drc*ChannelWidth->Drc);

            return std::max(std::min(tc,MAXCONCBL ),0.0);


        }else if(R_BL_Method == FSWUWANGJIA)
        {
            if(RBLD_D.at(_d)->Drc < 0.04)
            {
                return 0;
            }
            double slope = ChannelGrad->Drc;
            double ps = 2400.0;
            double pw = 1000.0;
            double h = ChannelWH->Drc;
            double n = std::max(ChannelN->Drc,0.001);
            double na = pow(graindiameters.at(_d)/100000.0,(1.0/6.0))/20.0;
            double phk = 0;
            double pek = 0;
            FOR_GRAIN_CLASSES
            {
                phk += RW_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                pek += RW_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
            }
            double ppk = 1;
            ppk = pow(phk/pek,0.6);
            if(pek == 0 )
            {
                return 0;
            }

            double dh = (ChannelWidth->Drc *h)/(ChannelWidth->Drc + 2* h);
            double css = 0.03* (ps - pw) * (graindiameters.at(_d)/1000000.0) * ppk;

            double qs = 0.0053 *pow(std::max(pow(na/n,1.5)*((pw * dh * 9.81 * 0.1 * slope/css)) - 1.0, 0.0),2.2);
            qs = qs * 1 * sqrt((ps/pw - 1)*9.81*pow(graindiameters.at(_d)/1000000.0,3.0));

            double tc = ps * ChannelWidth->Drc * qs/ (v * ChannelBLDepth->Drc*ChannelWidth->Drc);

            return std::max(std::min(tc,MAXCONCBL ),0.0);

        }else
        {

            return 0;
        }



}
//---------------------------------------------------------------------------
/**
 * @fn double TWorld::RiverSedimentTCSS(int r,int c, int _d)
 * @brief Calculates River suspended layer transport capacity
 *
 * Calculates suspended load sediment transport capacity.
 * Based on either, Van Rijn (simplified or full)
 * or Wu wang & Jia (for multiclass sediment).
 *
 * @param r : The row nr
 * @param c : The column nr of the cell
 * @param _d : The grain class (only needed when grain size distribution is used)
 * @see R_BL_Method
 */
double TWorld::RiverSedimentTCSS(int r,int c, int _d)
{
    //when water height is insignificant, transport capacity is zero
    //this is necessary since some of the used equations have strange behaviour
    //for these water heights or velocities. (h and v outside of valid range)
    if(ChannelWH->Drc < MIN_HEIGHT || ChannelSSDepth->Drc < MIN_HEIGHT)
    {
        return 0;
    }
    double v = ChannelV->Drc;
    if(v < MIN_FLUX)
    {
        return 0;
    }
    if(ChannelWidth->Drc < MIN_FLUX)
    {
        return 0;
    }


        if(R_SS_Method == FSRIJN)
        {
            //Van rijn simplified (1984)


            double ucr;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            double ps = 2400.0;
            double pw = 1000.0;
            double mu = 1.0;
            if( d50m < 0.0005)
            {
               ucr  = 0.19 * pow(d50m, 0.1) * log10(4.0* ChannelSSDepth->Drc/d90m);
            }else
            {
               ucr  = 8.5 * pow(d50m, 0.6) * log10(4.0* ChannelSSDepth->Drc/d90m);
            }
            double me = std::max((v - ucr)/sqrt(GRAV * d50m * (ps/pw - 1)),0.0);
            double ds = d50m * GRAV * ((ps/pw)-1)/(mu*mu);
            double qs = ChannelWH->Drc * 0.008 * ps*v * d50m * pow(me, 2.4) * pow(ds, -0.6);

            double tc =  qs/ (v * ChannelSSDepth->Drc);
            return std::max(std::min(tc,MAXCONC),0.0);
        }else if(R_SS_Method == FSRIJNFULL)
        {

            //van Rijn full (1980)
            //van Rijn full (1980)

            double kinvis = 1.0;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            double ps = 2400.0;
            double pw = 1000.0;
            double ds = D50->Drc * pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0));
            double dh = ChannelWH->Drc;
            double chevey = 18 * log10(4 * dh/d90m);
            double a = 0.1;

            double us = v* sqrt(GRAV)/chevey;
            double uscr = 0.055;

            if(ds <150 && !(ds < 20))
                uscr = 0.013*pow(ds,0.29);
            if(ds < 20 && !(ds < 10))
                uscr = 0.04*pow(ds,-0.10);
            if(ds <10 && !(ds < 4))
                uscr = 0.14*pow(ds,-0.64);
            if(ds <4)
                uscr = 0.24*pow(ds,-1);

            uscr = sqrt(uscr * (ps/pw - 1)*GRAV * d50m);

            double T = std::max((us*us/(uscr*uscr) - 1),0.0);
            double bsv = sqrt(GRAV * ChannelWH->Drc *std::max(ChannelGrad->Drc,0.05));
            double ca = 0.015 * (d50m/a) * pow(T,1.5)/pow(ds,0.3);

            double dss = 1 + 0.011*(1.8 - 1)*(T - 25);
            double sv = 10 * (kinvis/ds) *( sqrt(1 + (ps/pw - 1) * GRAV * d50m*d50m*d50m) - 1);

            double beta = std::min(1.0 + 2.0*(sv/bsv)*(sv/bsv),5.0);
            double powcb =1;
            if(ChannelBLConc->Drc > 0)
            {
                powcb = 0.1;//pow(ca/BLCFlood->Drc,0.4);
            }
            double phi = 2.5 * pow(sv/bsv,0.8) * powcb;
            double Z = sv/(beta*bsv*0.40);
            double Zs = Z + phi;
            double ad = 0.1;
            double F = (pow(ad,Zs) - pow(ad,1.2))/(pow(1.0-ad,Zs)* (1.2 - Zs));
            double qs = F * v * ChannelWH->Drc * ca;
            double tc = ps * ChannelWidth->Drc * qs/ (v * ChannelSSDepth->Drc * ChannelWidth->Drc);
            return std::max(std::min(tc,MAXCONC ),0.0);
        }else if(R_SS_Method == FSWUWANGJIA)
        {

            if(RSSD_D.at(_d)->Drc <0.04  )
            {
                return 0;
            }
            double slope = ChannelGrad->Drc;
            double ps = 2400.0;
            double pw = 1000.0;
            double h = ChannelWH->Drc;
            double phk = 0;
            double pek = 0;
            double sv = settlingvelocities.at(_d);
            double gd = graindiameters.at(_d)/1000000.0;
            FOR_GRAIN_CLASSES
            {
                phk += RW_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                pek += RW_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
            }
            double ppk = 1;
            ppk = pow(phk/pek,0.6);
            if(pek == 0 )
            {
                return 0;
            }

            double dh = (ChannelWidth->Drc *h)/(ChannelWidth->Drc + 2.0* h);
            double css = 0.03* (ps - pw) * (gd) * ppk;

            double qs = 0.0000262 *pow(std::max(( pw * 0.01 * h * 9.81 * slope /css) - 1.0, 0.0)* v/(sqrt(sv)),2.2);
            qs = qs * 1 * sqrt((ps/pw - 1)*9.81*pow(gd,3.0));

            double tc = ps * ChannelWidth->Drc * qs/ (v * ChannelSSDepth->Drc * ChannelWidth->Drc);

            return std::max(std::min(tc,MAXCONC ),0.0);

        }else
        {

            return 0;
        }




}




//--------------------------------------------------------------------------
//NOT USED!!!!!!
void TWorld::SumSedimentClasses(void)
{
   //not needed!
   //sediment budget (different sediment classes vs total sediment) is kept tight troughout the model!
   //this way later coding is much easier

}
//---------------------------------------------------------------------------
