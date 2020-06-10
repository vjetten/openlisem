//
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
 * @brief Calculates concentration with a maximum of MAXCONC, changes sed vol and deposisition
 *
 * @param watvol : the watervolume
 * @param sedvol : the sediment mass
 * @return sediment concentration (kg/m3)
 * @see MAXCONC
 *
 */
double TWorld::MaxConcentration(double watvol, double *sedvol, double *dep)
{
    double conc = 0;
    // when if activate, MBs error in KINDYN !!! Bizarre
    if (watvol > 0)
        conc = std::min(*sedvol/(watvol+0.001), MAXCONC);   // 1e-6 is 1 ml/m2 !!
    //      conc = MAXCONC;
    // *dep -= *sedvol;
    //*sedvol = 0;
    //}

    //    if (conc > MAXCONC) {
    //   double maxsed = MAXCONC*watvol;
    //    *dep += maxsed-*sedvol;
    //      conc = MAXCONC;
    //    *sedvol = maxsed;
    //}
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
void TWorld::SplashDetachment()
{
    if (!SwitchErosion)
        return;

#pragma omp parallel for collapse(2)
    FOR_ROW_COL_MV_L
            if(WH->Drc > 0.0001)
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

        double directrain = (1-Litter->Drc) * (1-Cover->Drc)*Rainc->Drc * 1000;
        // Added litter also to directrain, assme it covers the entire cell, not only under the plant
        // rainfall between plants in mm

        double KE_LD = std::max(15.3*sqrt(PlantHeight->Drc)-5.87, 0.0);
        // kin energy in J/m2/mm
        double throughfall = (1-Litter->Drc) * Cover->Drc * LeafDrain->Drc * 1000;
        // leaf drip in mm, is calculated as plant leaf drip in interception function so mult cover
        // VJ 110206 stemflow is also accounted for

        double WH0 = exp(-1.48*hmxWH->Drc*1000);
        // water buffer effect on surface, WH in mm in this empirical equation from Torri ?

        if (AggrStab->Drc > 0)
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
            //          double depdepth = std::max((StorageDep->Drc / (1600.0))/(_dx * DX->Drc),0.0);
            double depdepth = std::max((StorageDep->Drc / BulkDens)/(_dx * DX->Drc),0.0);
            // VJ 170308 changed to bulkdensity input
            double fac1 = std::max(0.0,1.0 - depdepth/SedimentMixingDepth->Drc);
            double fac2 = 1.0 - fac1;

            strength = strength * fac2 + (0.1033/DepositedCohesion) * fac1;
            b = b * fac2 + 3.58 * fac1;
        }
        double FPA = 1.0;
        if (RR->Drc > 0.1)
            FPA =  1-exp(-1.875*(WH->Drc/(0.01*RR->Drc)));

        // Between plants, directrain is already with 1-cover
        DetDT1 = g_to_kg * FPA*(strength*KE_DT+b)*WH0 * directrain;
        //ponded areas, kg/m2/mm * mm = kg/m2
        DetDT2 = g_to_kg * (1-FPA)*(strength*KE_DT+b) * directrain * SplashDelivery;
        //dry areas, kg/m2/mm * mm = kg/m2

        if (SwitchKETimebased)
        {
            if (directrain > 0)
            {
                DetDT1 = g_to_kg * FPA*(strength*KE_DT+b)*WH0 * _dt/3600;
                //ponded areas, kg/m2/sec * sec = kg/m2
                DetDT2 = g_to_kg * (1-FPA)*(strength*KE_DT+b) * _dt/3600 * SplashDelivery;
                //dry areas, kg/m2/sec * sec = kg/m2
            }
        }
        //based on work by Juan Sanchez

        // Under plants, throughfall is already with cover
        DetLD1 = g_to_kg * FPA*(strength*KE_LD+b)*WH0 * throughfall;
        //ponded areas, kg/m2/mm * mm = kg/m2
        DetLD2 = g_to_kg * (1-FPA)*(strength*KE_LD+b) * throughfall * SplashDelivery;
        //dry areas, kg/m2/mm * mm = kg/m2

        DETSplash->Drc = DetLD1 + DetLD2 + DetDT1 + DetDT2;
        // Total splash kg/m2

        // Deal with all exceptions:

        DETSplash->Drc *= (SoilWidthDX->Drc*DX->Drc);
        // kg/cell, only splash over soilwidth, not roads and channels
        // FROM KG/M2 TO KG/CELL

        DETSplash->Drc = (1-StoneFraction->Drc) * DETSplash->Drc;
        // no splash on stone surfaces

        if (SwitchGrassStrip)
            DETSplash->Drc = (1-GrassFraction->Drc) * DETSplash->Drc;

        //      if(SwitchSedtrap)
        //          DETSplash->Drc = (1-SedimentFilter->Drc) * DETSplash->Drc;

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
//            if(SwitchUseGrainSizeDistribution)
//            {
//                StorageDep->Drc = GetTotalDW(r,c,&StorageDep_D);
//                Storage->Drc = GetTotalDW(r,c,&Storage_D);
//            }

            //check wat we can detach from the top and bottom layer of present material
            double dleft = DETSplash->Drc ;
            double deptake = 0;
            double mattake = 0;
            double detachment = 0;

            deptake = std::min(dleft,StorageDep->Drc);
            StorageDep->Drc -= deptake;
            // det not more than storage
            // decrease store depth

//            if(SwitchUseGrainSizeDistribution)
//            {
//                double wtotal = 0;
//                FOR_GRAIN_CLASSES
//                {
//                    wtotal += StorageDep_D.Drcd;
//                }
//                if(wtotal != 0)
//                {
//                    FOR_GRAIN_CLASSES
//                    {
//                        StorageDep_D.Drcd = StorageDep_D.Drcd * StorageDep->Drc/wtotal;
//                    }
//                }
//            }

            detachment += deptake;
            // detachment is now taken material


            if(!(Storage->Drc < 0))
            {
                mattake = std::min(dleft,Storage->Drc);
                Storage->Drc -= mattake;

//                if(SwitchUseGrainSizeDistribution)
//                {
//                    double wtotal = 0;
//                    FOR_GRAIN_CLASSES
//                    {
//                        wtotal += Storage_D.Drcd;
//                    }
//                    if(wtotal != 0)
//                    {
//                        FOR_GRAIN_CLASSES
//                        {
//                            Storage_D.Drcd = Storage_D.Drcd * Storage->Drc/wtotal;
//                        }
//                    }
//                }
                detachment += mattake;
            }else
            {
                detachment += dleft;
            }
            DETSplash->Drc = detachment;
        }


        //      if (hmx->Drc > 0) {
        //          SSFlood->Drc += DETSplash->Drc;
        //          SSCFlood->Drc = MaxConcentration(ChannelAdj->Drc * DX->Drc * hmx->Drc, &SSFlood->Drc, &DepFlood->Drc);

        //      } else {
        //          Sed->Drc += DETSplash->Drc;
        //          Conc->Drc = MaxConcentration(WaterVolall->Drc, &Sed->Drc, &DEP->Drc);
        //      }


        // IN KG/CELL
    }

}
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
    double w = (*M).at(0)->Drc;
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
    double w = (*M).at(0)->Drc;
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
    if (SwitchSV == 2) {
        double dm = d / 1e6;
        double ds = dm * pow((2650.0/1000.0 - 1.0)*GRAV/(1e-6*1e-6),(1.0/3.0));
        return 1e-6/dm*(ds*ds*ds)*pow(38.1+0.93*pow(ds,12.0/7.0), -7.0/8.0);
        //    // zhiyao et al, 2008
    } else {
        if(d < 100)
        {
            return 2*(2650.0-1000.0)*GRAV*pow(d/2000000.0, 2)/(9*0.001);
            //Stokes range settling velocity
        }else
        {
            double dm = d/1000.0;
            return 10.0 *(sqrt(1.0 + 0.01 *((2650.0-1000.0)/1000.0)* GRAV *dm*dm*dm )-1.0)/(dm);
            //Settling velocity by Zanke (1977)
        }
    }
}
/**
 * @fn void TWorld::SedimentSetMaterialDistribution(int r,int c)
 * @brief Update grain size distribution of soil layers
 *
 * Update grain size distribution of both
 * original soil layer and deposited layer in river and on slope.
 * Step is skipped if the storage is negative,
 * which indicated infinite storage.
 *
 * @param r : row nr of the cell
 * @param c : coumn nr of the cell
 * @return void
 */
void TWorld::SedimentSetMaterialDistribution()
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
                double depdepth = std::max((StorageDep->Drc / (BulkDens))/(_dx * DX->Drc),0.0);
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
                    double depdepth = std::max((RStorageDep->Drc / (BulkDens))/(_dx * DX->Drc),0.0);
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

// TODO: check what happens in this function
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
            return detachment *= Y->Drc;
        }
    }

    //first check if it is channel detachment
    if(channel)
    {
        //calculate depth of deposited layer
        double depdepth = std::max((RStorageDep->Drc / (BulkDens))/(ChannelWidth->Drc * DX->Drc),0.0);

        //linear decrease in influence from lower soil layer
        //from 0 to MixingDepth, with fac1 for bottom layer, fac2 for top layer
        double fac1 = 1.0;
        if(RSedimentMixingDepth->Drc > MIN_HEIGHT)
            fac1 = std::max(0.0,1.0 - depdepth/RSedimentMixingDepth->Drc);
        double fac2 = 1-fac1;

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

    } else
        if(flood) {
            double depdepth = std::max((StorageDep->Drc / (BulkDens))/(_dx * DX->Drc),0.0);

            //linear decrease in influence from lower soil layer
            //from 0 to MixingDepth, with fac1 for bottom layer, fac2 for top layer
            double fac1 = 1.0;
            if(RSedimentMixingDepth->Drc > MIN_HEIGHT)
                fac1 = std::max(0.0,1.0 - depdepth/RSedimentMixingDepth->Drc);
            double fac2 = 1-fac1;
            //new erosion coefficient bases on soil layer mixinfactors
            double newY = Y->Drc * fac1 + fac2 * 1.0;
            //            if (SwitchSedtrap && SedimentFilter->Drc > 0)
            //                newY *= 1.0-SedimentFilter->Drc;

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

            //check wat we can detach from the top and bottom layer of present material
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
        } else {
            //calculate depth of deposited layer
            double depdepth = std::max((StorageDep->Drc / (BulkDens))/(_dx * DX->Drc),0.0);

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

// TODO: CHECK consistency SWOF
void TWorld::FlowDetachment()
{
    if (!SwitchErosion)
        return;

    if (SwitchKinematic2D == K2D_METHOD_DYN)
        return;

    //transport capacity
    #pragma omp parallel for collapse(2)
    FOR_ROW_COL_MV_L  {
        DETFlow->Drc = 0;
        DEP->Drc = 0;
        if(!SwitchUseGrainSizeDistribution) {
            //get the transport capacity for a single grain size
            TC->Drc = calcTCSuspended(r,c,-1, FSGOVERS, WHrunoff->Drc, V->Drc, 2);

        }
//        else {
//            //get the transport capacity for all induvidual grain sizes
//            TC->Drc = 0;
//            FOR_GRAIN_CLASSES
//            {
//                TC_D.Drcd = calcTCSuspended(r,c,d, FSHAIRSINEROSE, WHrunoff->Drc, V->Drc, 2);
//                TC->Drc += TC_D.Drcd;
//            }

//            //if the total transport capacity is too high, limit it to the maximum value (MAXCONC)
//            if(TC->Drc > MAXCONC)
//            {
//                //rescale all induvidual tc's to conform to MAXCONC
//                FOR_GRAIN_CLASSES
//                {
//                    TC_D.Drcd *= MAXCONC/TC->Drc;
//                }
//                TC->Drc = MAXCONC;
//            }

//            TC->Drc = 0;
//            FOR_GRAIN_CLASSES
//            {
//                TC->Drc += TC_D.Drcd;
//            }
//        }
    }


    //the iterator is either the number of grain classes, or 1 if no grain size distribution is used.
    int iterator = numgrainclasses;
    if(!SwitchUseGrainSizeDistribution){
        iterator = 1;
    }

    //for each grain class, calculate flow detachment seperately
    for(int d  = 0 ; d < iterator;d++)
    {
        #pragma omp parallel for collapse(2)
        FOR_ROW_COL_MV_L
        {
            double erosionwh = WHrunoff->Drc;
            double erosionwv = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc;

            if (!SwitchUseGrainSizeDistribution)
                Sed->Drc += DETSplash->Drc;
            else
                Sed_D.Drcd += DETSplash->Drc * W_D.Drcd;

            if (erosionwh < HMIN) {
                DEP->Drc += -Sed->Drc;
                Sed->Drc = 0;
                Conc->Drc = 0;
                TC->Drc = 0;

                // TOSO: add grainsize distribution here

            } else {

                if(!SwitchUseGrainSizeDistribution)
                    Conc->Drc = MaxConcentration(erosionwv, &Sed->Drc, &DEP->Drc);
                else
                    Conc_D.Drcd = MaxConcentration(erosionwv, &Sed_D.Drcd, &DEP->Drc);

                double maxTC = 0;
                double minTC = 0;

                double deposition = 0;
                double detachment = 0;
                double TransportFactor = 0;

                if(!SwitchUseGrainSizeDistribution)
                {
                    maxTC = std::max(TC->Drc - Conc->Drc,0.0);
                    // positive difference: TC deficit becomes detachment (positive)
                    minTC = std::min(TC->Drc - Conc->Drc,0.0);
                    // negative difference: TC surplus becomes deposition (negative)
                    // unit kg/m3
                }
                //                else
                //                {
                //                    maxTC = std::max(TC_D.Drcd - Conc_D.Drcd,0.0);
                //                    minTC = std::min(TC_D.Drcd - Conc_D.Drcd,0.0);
                //                }

                if (minTC < 0) {

                    //### deposition

//                    if(!SwitchUseGrainSizeDistribution)
//                    {
//                        if (erosionwh > HMIN)
                            TransportFactor = (1-exp(-_dt*SettlingVelocitySS->Drc/erosionwh)) * erosionwv;
//                        else
//                            TransportFactor = erosionwv;
//                    }
                    //                    else
                    //                    {
                    //                        TransportFactor = (1-exp(-_dt*settlingvelocities.at(d)/erosionwh)) * erosionwv;
                    //                        //   TransportFactor = _dt*settlingvelocities.at(d) * DX->Drc * FlowWidth->Drc;
                    //                    }
                    // if settl velo is very small, transportfactor is 0 and depo is 0
                    // if settl velo is very large, transportfactor is 1 and depo is max

                    // deposition can occur on roads and on soil (so use flowwidth)

                    deposition = minTC * TransportFactor;
                    // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0

                    if(!SwitchUseGrainSizeDistribution)
                        deposition = std::max(deposition, -Sed->Drc);
                    else
                        deposition = std::max(deposition, -Sed_D.Drcd);

                    //force deposition on grass strips
                    if (SwitchGrassStrip) {
                        //                  if(!SwitchUseGrainSizeDistribution)
                        //                  {
                        //                      deposition = -Sed->Drc*GrassFraction->Drc + (1-GrassFraction->Drc)*deposition;
                        //                  } else {
                        //                      deposition = -Sed_D.Drcd*GrassFraction->Drc + (1-GrassFraction->Drc)*deposition;
                        //                  }
                    }


                    if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                        deposition = 0;
                    // VJ 190325 prevent any activity on the boundary!
                    if (SwitchSedtrap && SedMaxVolume->Drc == 0 && N->Drc == SedTrapN) {
                        N->Drc = Norg->Drc;
                    }
                    if (SwitchSedtrap && SedMaxVolume->Drc > 0)
                    {
                        if(!SwitchUseGrainSizeDistribution)
                        {
                            if (Sed->Drc > 0) {
                                double depvol = Sed->Drc * 1.0/BulkDens; // m3
                                if (SedMaxVolume->Drc < depvol)
                                    depvol = SedMaxVolume->Drc;
                                if (SedMaxVolume->Drc > 0){
                                    deposition = -depvol*BulkDens;
                                    maxTC = 0;
                                }
                                SedMaxVolume->Drc = SedMaxVolume->Drc - depvol;
                                SedimentFilter->Drc += depvol*BulkDens;
                            }

                        } else {
                            if (Sed_D.Drcd > 0) {
                                //                              deposition = -Sed_D.Drcd*SedimentFilter->Drc;
                                //                              Sed_D.Drcd *= (1.0-SedimentFilter->Drc);
                            }
                        }
                    }

                    //add deposition to soil layer
                    if(SwitchUseMaterialDepth)
                    {
                        StorageDep->Drc += -deposition;
                        //                        if(SwitchUseGrainSizeDistribution)
                        //                        {
                        //                            StorageDep_D.Drcd += -deposition;
                        //                        }
                    }
                }
                else
                if (maxTC > 0) {
                    //### detachment

                    //if(!SwitchUseGrainSizeDistribution)
                    //{
                        TransportFactor = _dt*SettlingVelocitySS->Drc * DX->Drc * SoilWidthDX->Drc; //fpa->Drc*
                        // soilwidth is erodible surface
                    //}
                    //                else
                    //                {
                    //                    TransportFactor = _dt*settlingvelocities.at(d) * DX->Drc * fpa->Drc*SoilWidthDX->Drc;
                    //                }
                   // TransportFactor = std::min(TransportFactor, Q->Drc*_dt);

                    // detachment can only come from soil, not roads (so do not use flowwidth)
                    // units s * m/s * m * m = m3

                    detachment = maxTC * TransportFactor;

//                    if(SwitchUseGrainSizeDistribution)
//                    {
//                        detachment =  W_D.at(d)->Drc * detachment;
//                    }
                    // unit = kg/m3 * m3 = kg (/cell)

                    // detachment = std::min(detachment, maxTC * erosionwv);
                    // cannot have more detachment than remaining capacity in flow

                    // exceptions
                    if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                        detachment = 0;
                    // VJ 190325 prevent any activity on the boundary!

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

                    detachment = (1-Snowcover->Drc) * detachment;
                    /* TODO: CHECK THIS no flow detachment on snow */
                    //is there erosion and sedimentation under the snowdeck?

                    detachment = DetachMaterial(r,c,d,false,false,false, detachment);

                    if(MAXCONC * erosionwv < Sed->Drc+detachment)
                        detachment = std::max(0.0, MAXCONC * erosionwv - Sed->Drc);
                    // not more detachment then is needed to keep below ssmax

                    if (SwitchSedtrap && SedMaxVolume->Drc > 0)
                        detachment = 0;
                } // minv > 0
                //### sediment balance
                // add to sediment in flow (IN KG/CELL)

                Sed->Drc += detachment;
                Sed->Drc += deposition;
                DETFlow->Drc += detachment;
                DEP->Drc += deposition;
                Conc->Drc = MaxConcentration(erosionwv, &Sed->Drc, &DEP->Drc);

//                if(SwitchUseGrainSizeDistribution)
//                {
//                    Sed_D.Drcd += detachment;
//                    Sed_D.Drcd += deposition;
//                    Conc_D.Drcd = MaxConcentration(erosionwv, &Sed_D.Drcd, &DEP->Drc);
//                }
            }
        } // FOR
    } // iterator
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
void TWorld::ChannelFlowDetachment()
{
    if (!SwitchIncludeChannel)
        return;

    if (!SwitchErosion)
        return;

    FOR_ROW_COL_MV_CH {

        RiverSedimentLayerDepth(r,c);
        //creates ChannelBLDepth and ChannelSSDepth, if 1 layer ChannelBLDepth = 0

        double sswatervol = ChannelSSDepth->Drc*DX->Drc*ChannelFlowWidth->Drc;
        double ssdischarge = ChannelV->Drc * ChannelFlowWidth->Drc * ChannelSSDepth->Drc;
        double blwatervol = 0;
        double bldischarge = 0;
        if (SwitchUse2Layer) {
            blwatervol = ChannelBLDepth->Drc*DX->Drc*ChannelFlowWidth->Drc;
            bldischarge = ChannelV->Drc * ChannelFlowWidth->Drc * ChannelBLDepth->Drc;
        }

        //iterator is number of grain classes
        int iterator = numgrainclasses;
        if(!SwitchUseGrainSizeDistribution) {
            iterator = 1;
        }

        ChannelDetFlow->Drc = 0;
        ChannelDep->Drc = 0;

        //find transport capacity for bed and suspended layer
        for(int d  = 0 ; d < iterator;d++)
        {
            cTMap * TBLTCtemp;
            cTMap * TSSTCtemp;

            if(!SwitchUseGrainSizeDistribution)
            {
                if (SwitchUse2Layer)
                    TBLTCtemp = ChannelBLTC;
                TSSTCtemp = ChannelSSTC;
            }
            //            else
            //            {
            //                TBLTCtemp = RBLTC_D.at(d);
            //                TSSTCtemp = RSSTC_D.at(d);
            //            }

            //get transport capacity for bed/suspended load for a specific cell and grain size class
            // TBLTCFlood->Drc = RiverSedimentTCBL(r,c,d, ChannelV->Drc, ChannelWH->Drc, ChannelBLDepth->Drc, ChannelWidth->Drc);
            // TSSTCFlood->Drc = RiverSedimentTCSS(r,c,d, ChannelV->Drc, ChannelWH->Drc, ChannelSSDepth->Drc, ChannelWidth->Drc);
            if (SwitchUse2Layer)
                TBLTCtemp->Drc = calcTCBedload(r, c, d, R_BL_Method, ChannelWH->Drc,ChannelV->Drc, 0);
            TSSTCtemp->Drc = calcTCSuspended(r, c, d, R_SS_Method, ChannelWH->Drc, ChannelV->Drc, 0);

            //check if the sum of transport capacities of all grain sizes is larger than MAXCONC, and rescale if nessecery
            //            if(SwitchUseGrainSizeDistribution)
            //            {
            //                ChannelBLTC->Drc = 0;
            //                ChannelSSTC->Drc = 0;
            //                FOR_GRAIN_CLASSES
            //                {
            //                    ChannelBLTC->Drc += RBLTC_D.Drcd;
            //                    ChannelSSTC->Drc += RSSTC_D.Drcd;
            //                }

            //                if(ChannelBLTC->Drc > MAXCONCBL)
            //                {
            //                    FOR_GRAIN_CLASSES
            //                    {
            //                        RBLTC_D.Drcd *= MAXCONCBL/ChannelBLTC->Drc;
            //                    }
            //                    ChannelBLTC->Drc = MAXCONCBL;
            //                }
            //                if(ChannelSSTC->Drc > MAXCONC)
            //                {
            //                    FOR_GRAIN_CLASSES
            //                    {
            //                        RSSTC_D.Drcd *= MAXCONC/ChannelSSTC->Drc;
            //                    }
            //                    ChannelSSTC->Drc = MAXCONC;
            //                }

            //                ChannelBLTC->Drc = 0;
            //                ChannelSSTC->Drc = 0;
            //                FOR_GRAIN_CLASSES
            //                {
            //                    ChannelBLTC->Drc += RBLTC_D.Drcd;
            //                    ChannelSSTC->Drc += RSSTC_D.Drcd;
            //                }
            //            }
        }
        for(int d  = 0 ; d < iterator;d++)
        {
            //set all maps for this grain class

            cTMap * TBLDepthtemp;
            cTMap * TSSDepthtemp;
            cTMap * TBLTCtemp;
            cTMap * TSSTCtemp;
            cTMap * TBLCtemp;
            cTMap * TSSCtemp;
            cTMap * TBLtemp;
            cTMap * TSStemp;
            cTMap * TW;

            double TSettlingVelocitySS;
            double TSettlingVelocityBL;

            if(!SwitchUseGrainSizeDistribution)
            {
                TSSDepthtemp = ChannelSSDepth;
                TSSTCtemp = ChannelSSTC;
                TSSCtemp = ChannelSSConc;
                TSStemp = ChannelSSSed;
                TSettlingVelocitySS = SettlingVelocitySS->Drc;
                if (SwitchUse2Layer) {
                    TBLDepthtemp = ChannelBLDepth;
                    TBLTCtemp = ChannelBLTC;
                    TBLCtemp = ChannelBLConc;
                    TBLtemp = ChannelBLSed;
                    TSettlingVelocityBL = SettlingVelocityBL->Drc;
                }
                TW = unity;
            }
            //            else
            //            {
            //                TBLDepthtemp = RBLD_D.at(d);
            //                TSSDepthtemp = RSSD_D.at(d);
            //                TBLTCtemp = RBLTC_D.at(d);
            //                TSSTCtemp = RSSTC_D.at(d);
            //                TBLCtemp = RBLC_D.at(d);
            //                TSSCtemp = RSSC_D.at(d);
            //                TBLtemp = RBL_D.at(d);
            //                TSStemp = RSS_D.at(d);
            //                TW = RW_D.at(d);
            //                TSettlingVelocity = settlingvelocities.at(d);
            //            }

            double deposition = 0;
            double detachment = 0;
            double TransportFactor = 0;
            double maxTC = 0;
            double minTC = 0;

            //when waterheight is insignificant, deposite all remaining sediment
            if(ChannelWH->Drc < MIN_HEIGHT)
            {
                deposition += -TSStemp->Drc;
                TSStemp->Drc = 0;
                TSSTCtemp->Drc = 0;
                TSSCtemp->Drc = 0;
                TSSDepthtemp->Drc = 0;

                if (SwitchUse2Layer) {
                    deposition = -TBLtemp->Drc;
                    TBLtemp->Drc = 0;
                    TBLTCtemp->Drc = 0;
                    TBLCtemp->Drc = 0;
                    TBLDepthtemp->Drc = 0;
                }

                ChannelDep->Drc += deposition;

                if(SwitchUseMaterialDepth)
                {
                    RStorageDep->Drc += -deposition;
                    //                    if(SwitchUseGrainSizeDistribution)
                    //                    {
                    //                        RStorageDep_D.Drcd += -deposition;
                    //                        if(std::isnan(RStorageDep_D.Drcd))
                    //                        {
                    //                            qDebug() << "NAN dep1" << d;
                    //                        }
                    //                    }
                }

                //                if(SwitchUseGrainSizeDistribution)
                //                {
                //                    RBL_D.Drcd = 0;
                //                    RSS_D.Drcd = 0;
                //                    RBLTC_D.Drcd = 0;
                //                    RSSTC_D.Drcd = 0;
                //                    RBLC_D.Drcd = 0;
                //                    RSSC_D.Drcd = 0;
                //                }
            } else {
                //### do suspended first

                //deposition
                maxTC = std::max(TSSTCtemp->Drc - TSSCtemp->Drc, 0.0);  // TC in kg/m3
                minTC = std::min(TSSTCtemp->Drc - TSSCtemp->Drc, 0.0);

                if (minTC < 0) {
                    TransportFactor = (1-exp(-_dt*TSettlingVelocitySS/ChannelWH->Drc)) * sswatervol;

                    deposition = std::max(TransportFactor * minTC,-TSStemp->Drc); // in kg
                    // not more than SS present
                } else {
                    //    if (maxTC > 0) {

                    //  detachment
                    TransportFactor = _dt*TSettlingVelocitySS * ChannelDX->Drc * ChannelFlowWidth->Drc;
                    TransportFactor = std::min(TransportFactor, ssdischarge*_dt);
                    //TransportFactor = ssdischarge*_dt;
                    // use discharge because standing water has no erosion

                    //NB ChannelFlowWidth and ChannelWidth the same woith rect channel
                    detachment = TW->Drc * maxTC * TransportFactor;
                    // cannot have more detachment than remaining capacity in flow

                    detachment = DetachMaterial(r,c,d,true,false,false, detachment);
                    // multiply by Y

                    if(MAXCONC * sswatervol < TSStemp->Drc+detachment)
                        detachment = std::max(0.0, MAXCONC * sswatervol - TSStemp->Drc);
                }

                //### sediment balance add suspended
                TSStemp->Drc += detachment;
                TSStemp->Drc += deposition;
                ChannelDep->Drc += deposition;
                ChannelDetFlow->Drc += detachment;
                ChannelTC->Drc = ChannelSSTC->Drc;
                ChannelSed->Drc = ChannelSSSed->Drc;
                //total transport capacity (bed load + suspended load), used for output


                if(SwitchUseMaterialDepth)
                {
                    RStorageDep->Drc += -deposition;
                    //                if(SwitchUseGrainSizeDistribution)
                    //                    RStorageDep_D.Drcd += -deposition;
                }

                //### do bedload
                if (SwitchUse2Layer) {
                    if(TBLDepthtemp->Drc < MIN_HEIGHT) {
                        ChannelDep->Drc += -ChannelBLSed->Drc;
                        ChannelBLTC->Drc = 0;
                        ChannelBLConc->Drc = 0;
                        ChannelBLSed->Drc = 0;

                    } else {
                        //there is BL

                        maxTC = std::max(TBLTCtemp->Drc - TBLCtemp->Drc,0.0);
                        minTC = std::min(TBLTCtemp->Drc - TBLCtemp->Drc,0.0);

                        if (maxTC > 0) {
                            //### detachment
                            TransportFactor = _dt*TSettlingVelocityBL * ChannelDX->Drc * ChannelFlowWidth->Drc;
                            //              TransportFactor = std::min(TransportFactor, bldischarge*_dt);
                            //              TransportFactor = bldischarge*_dt;

                            // units s * m/s * m * m = m3
                            detachment =  TW->Drc *  maxTC * TransportFactor;
                            // unit = kg/m3 * m3 = kg

                            detachment = DetachMaterial(r,c,d,true,false,true, detachment);
                            // mult by Y and mixingdepth
                            // IN KG/CELL

                            if(MAXCONC * blwatervol < TBLtemp->Drc+detachment)
                                detachment = std::max(0.0, MAXCONC * blwatervol - TBLtemp->Drc);
                            // limit detachment to what BLtemp can carry
                        } else
                            //  if (minTC < 0)
                        {
                            //### deposition
                            if (TBLDepthtemp->Drc > MIN_HEIGHT)
                                TransportFactor = (1-exp(-_dt*TSettlingVelocityBL/TBLDepthtemp->Drc)) * blwatervol;
                            else
                                TransportFactor =  1.0 * blwatervol;

                            // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0
                            deposition = std::max(minTC * TransportFactor, -TBLtemp->Drc);
                            // cannot have more depo than sediment present

                            if(SwitchUseMaterialDepth)
                            {
                                RStorageDep->Drc += -deposition;
                                //                        if(SwitchUseGrainSizeDistribution)
                                //                        {
                                //                            RStorageDep_D.Drcd += -deposition;
                                //                        }
                            }

                            TBLtemp->Drc += detachment; //ChannelBLSed
                            TBLtemp->Drc += deposition;
                            ChannelDep->Drc += deposition;
                            ChannelDetFlow->Drc += detachment;
                            ChannelTC->Drc += ChannelBLTC->Drc;
                            ChannelSed->Drc += ChannelBLSed->Drc;
                            //total transport capacity (bed load + suspended load), used for output
                        }
                    }
                }
            }
            //        if(SwitchUseGrainSizeDistribution)
            //        {
            //            FOR_GRAIN_CLASSES
            //            {
            //                RSSC_D.Drcd = MaxConcentration(ChannelWidth->Drc*DX->Drc*RSSD_D.Drcd, &RSS_D.Drcd, &ChannelDep->Drc);
            //                double sssmax = MAXCONC * DX->Drc *ChannelWidth->Drc*RSSD_D.Drcd;
            //                if(sssmax < RSS_D.Drcd) {
            //                    double surplus = RSS_D.Drcd - sssmax;
            //                    ChannelDep->Drc -= surplus;
            //                    RSS_D.Drcd = sssmax;
            //                    if(SwitchUseMaterialDepth)   // TODO: does not work with this maxconc !!!!!!
            //                    {
            //                        RStorageDep->Drc += surplus;
            //                        RStorageDep_D.Drcd += surplus;
            //                        if(std::isnan(RStorageDep_D.Drcd))
            //                        {
            //                            qDebug() << "NAN dep3" << d;
            //                        }
            //                    }
            //                }
            //                ChannelSSSed->Drc += RSS_D.Drcd;

            //                RBLC_D.Drcd = MaxConcentration(ChannelFlowWidth->Drc*DX->Drc*RBLD_D.Drcd, &RBL_D.Drcd, &ChannelDep->Drc);
            //                sssmax = MAXCONCBL * DX->Drc *ChannelFlowWidth->Drc*RBLD_D.Drcd;
            //                if(sssmax < BL_D.Drcd) {
            //                    ChannelDep->Drc -= (RBL_D.Drcd - sssmax);
            //                    RBL_D.Drcd = sssmax;
            //                    if(SwitchUseMaterialDepth)
            //                    {
            //                        RStorageDep->Drc += (RBL_D.Drcd - sssmax);
            //                        RStorageDep_D.Drcd += (RBL_D.Drcd - sssmax);
            //                        if(std::isnan(RStorageDep_D.Drcd))
            //                        {
            //                            qDebug() << "NAN dep4" << d;
            //                        }
            //                    }
            //                }
            //                ChannelBLSed->Drc += RBL_D.Drcd;
            //            }
        }


        RiverSedimentMaxC(r,c);
        //partial and total concentration ALL DONE
    }
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

    fromChannelWHtoVol(r, c);
    double frac = ChannelSSDepth->Drc/ChannelWH->Drc;
    //maximum concentration
    if(!SwitchUseGrainSizeDistribution)
    {
        _SSC->Drc = MaxConcentration(ChannelWaterVol->Drc*frac, &_SS->Drc, &ChannelDep->Drc);
        if (SwitchUse2Layer)
            _BLC->Drc = MaxConcentration(ChannelWaterVol->Drc*(1-frac), &_BL->Drc, &ChannelDep->Drc);
    }
//    else {
//        FOR_GRAIN_CLASSES
//        {
//            RSSC_D.Drcd = MaxConcentration(ChannelFlowWidth->Drc*DX->Drc*RSSD_D.Drcd, &RSS_D.Drcd, &ChannelDep->Drc);
//            // limit concentration to 850 and throw rest in deposition

//            double sssmax = MAXCONC * DX->Drc *ChannelFlowWidth->Drc*RSSD_D.Drcd;
//            if(sssmax < RSS_D.Drcd)
//            {
//                ChannelDep->Drc += (RSS_D.Drcd - sssmax);
//                ChannelSSSed->Drc += -(RSS_D.Drcd - sssmax);
//                if(SwitchUseMaterialDepth)
//                {
//                    RStorageDep->Drc += (RSS_D.Drcd - sssmax);
//                    RStorageDep_D.Drcd += (RSS_D.Drcd - sssmax);
//                    if(std::isnan(RStorageDep_D.Drcd))
//                    {
//                        qDebug() << "NAN dep3" << d;
//                    }
//                }
//                RSS_D.Drcd = sssmax;

//            }


//            RBLC_D.Drcd = MaxConcentration(ChannelFlowWidth->Drc*DX->Drc*RBLD_D.Drcd, &RBL_D.Drcd, &ChannelDep->Drc);
//            // limit concentration to 850 and throw rest in deposition

//            sssmax = MAXCONCBL * DX->Drc *ChannelFlowWidth->Drc*RBLD_D.Drcd;
//            if(sssmax < BL_D.Drcd)
//            {
//                ChannelDep->Drc += (RBL_D.Drcd - sssmax);
//                ChannelBLSed->Drc += -(RBL_D.Drcd - sssmax);
//                RBL_D.Drcd = sssmax;
//                if(SwitchUseMaterialDepth)
//                {
//                    RStorageDep->Drc += (RBL_D.Drcd - sssmax);
//                    RStorageDep_D.Drcd += (RBL_D.Drcd - sssmax);
//                    if(std::isnan(RStorageDep_D.Drcd))
//                    {
//                        qDebug() << "NAN dep4" << d;
//                    }
//                }
//            }
//        }

//        ChannelBLConc->Drc = 0;
//        ChannelSSConc->Drc = 0;

//        FOR_GRAIN_CLASSES
//        {
//            ChannelBLConc->Drc += RBLC_D.Drcd;
//            ChannelSSConc->Drc += RSSC_D.Drcd;
//        }
//    }

    ChannelSed->Drc = (SwitchUse2Layer ? _BL->Drc : 0) + _SS->Drc;
    //total concentration
    ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, &ChannelSed->Drc, &ChannelDep->Drc);


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
void TWorld::RiverSedimentDiffusion(double dt, cTMap *_SS, cTMap *_SSC)
{

    FOR_ROW_COL_MV_CH {
        _SSC->Drc = MaxConcentration(ChannelWaterVol->Drc, &_SS->Drc, &ChannelDep->Drc);
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
        //and devide by turbulent prandtl-smidth number, def 1.0
        double eta = eddyvs/R_SigmaDiffusion;

        //add diffusive fluxes to previous cell in channel.
        if(foundp)
        {
            double coeff1 = std::min(dt*eta *std::min(1.0,ChannelSSDepth->data[rp][cp]/ChannelSSDepth->data[r][c]),
                                     courant_factor/2.0) * _SS->Drc;
            //                                     courant_factor_diffusive/2.0) * MSSNFlood->Drc;

            _SS->data[rp][cp] += coeff1;
            _SS->data[r][c] -= coeff1;
        }

        //add diffusive fluxes to next cell in channel.
        if(foundn)
        {
            double coeff2 = std::min(dt*eta *std::min(1.0,ChannelSSDepth->data[rn][cn]/ChannelSSDepth->data[r][c]),
                                     courant_factor/2.0) * _SS->Drc;

            _SS->data[rn][cn] += coeff2;
            _SS->data[r][c] -= coeff2;
        }
    }

    //recalculate concentrations
    FOR_ROW_COL_MV_CH
    {
        //set concentration from present sediment
        _SSC->Drc = MaxConcentration(ChannelWaterVol->Drc, &_SS->Drc, &ChannelDep->Drc);
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
    if (!SwitchUse2Layer) {
        ChannelSSDepth->Drc = ChannelWH->Drc;
       // ChannelBLDepth->Drc = 0;
        return;
    }

    double ps = 2650;
    double pw = 1000;
    double factor = 0.5;
    double R = (ChannelWidth->Drc * ChannelWH->Drc)/(ChannelWH->Drc * 2 + ChannelWidth->Drc);

    if(!SwitchUseGrainSizeDistribution)
    {
        //if a two phase system is modelled, calculate thickness of layer
        double d50m = (D50->Drc/1000000.0);
        double d90m = (D90->Drc/1000000.0);

        //critical shear velocity for bed level motion by van rijn
        double critshearvel = ChannelV->Drc * sqrt(GRAV)/(18 * log10(4*R/d90m));
        //critical shear stress for bed level motion by van rijn
        double critsheart = (critshearvel*critshearvel)/ (((ps-pw)/pw) * GRAV*d50m);
        //rough bed bed load layer depth by Hu en Hui
        ChannelBLDepth->Drc = std::min(std::min(d50m * 1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), factor*ChannelWH->Drc), 0.1);
        ChannelSSDepth->Drc = std::max(ChannelWH->Drc - ChannelBLDepth->Drc,0.0);

    }
//    else {
//        ChannelBLDepth->Drc = 0;
//        ChannelSSDepth->Drc = 0;
//        FOR_GRAIN_CLASSES
//        {
//            //for each grain size, use the grain class property
//            double d50m = graindiameters.at(d)/1000000.0;
//            double d90m = /* 1.5 * */ graindiameters.at(d)/1000000.0; // why 1.5 ???
//            //critical shear velocity for bed level motion by van rijn
//            double critshearvel = ChannelV->Drc * sqrt(GRAV)/(18 * log10(4*R/d90m));
//            //critical shear stress for bed level motion by van rijn
//            double critsheart = (critshearvel*critshearvel)/ (((ps-pw)/pw) * GRAV*d50m);
//            //rough bed bed load layer depth by Hu en Hui
//            RBLD_D.Drcd = std::min(std::min(d50m * 1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), factor*ChannelWH->Drc), 0.1);
//            RSSD_D.Drcd = std::max(ChannelWH->Drc - RBLD_D.Drcd ,0.0);
//            ChannelBLDepth->Drc += RBLD_D.Drcd * RW_D.Drcd;
//            ChannelSSDepth->Drc += RSSD_D.Drcd * RW_D.Drcd;
//        }
//    }
}
//---------------------------------------------------------------------------
/**
 * @fn double TWorld::calcTCSuspended(int r,int c, int _d, int method,double U, int type)
 * @brief Calculates suspended layer transport capacity
 *
 * Calculates suspended load sediment transport capacity.
 * Based on govers, Van Rijn (simplified or full)
 * or Wu wang & Jia (for multiclass sediment).
 *
 * @param _d : The grain class (only needed when grain size distribution is used)
 * @param method : the TC method used
 * @param U : velicty can be channel of overland or flood
 * @param type : channel (0) or flood (1) or overland (2)
 */
double TWorld::calcTCSuspended(int r,int c, int _d, int method, double h, double U, int type)
{
    double R=0, hs=0, S = 0, w = 0;
    cTMap *Wd = nullptr;

    if (type == 0) {
        hs = ChannelSSDepth->Drc;
        S = ChannelGrad->Drc;
        w = ChannelWidth->Drc;
        R = (w*h)/(2*h+w);
        if (SwitchMulticlass)
            Wd = RW_D.at(_d);
    } else
        if (type == 1) {
            hs = SSDepthFlood->Drc;
            S = Grad->Drc;
            w = ChannelAdj->Drc;
            R = (w*h)/(2*h+w);
            if (SwitchMulticlass)
                Wd = W_D.at(_d);
        } else
            if (type == 2) {
                hs = WHrunoff->Drc;
                S = Grad->Drc;
                w = FlowWidth->Drc;
            }

    //when water height is insignificant, transport capacity is zero
    //this is necessary since some of the used equations have strange behaviour
    //for these water heights or velocities. (h and v outside of valid range)
    if(h < MIN_HEIGHT || hs < MIN_HEIGHT)
        return 0;
    if(U < MIN_FLUX)
        return 0;

    double ps = 2650.0; //2400.0;
    double pw = 1000.0;
    double d50m = (D50->Drc/1000000.0);
    double d90m = (D90->Drc/1000000.0);
    double tc = 0;

    if(method == FSHAIRSINEROSE)
    {
        //        tc =  d50m * 1.0/SettlingVelocity->Drc * 0.013/GRAV * ps/(ps-pw) * ( std::max(0.0, (om - omcr))/h) ;
        double om =  U*S;
        double omcr = 0.004;
        tc =  Wd->Drc * (1.0/settlingvelocities.at(_d))*(0.013/GRAV) * (2650.0/(2650.0 - 1000.0)) * ( std::max(0.0, (om - omcr))/h) ;

    } else
        if(method == FSGOVERS)
        {
            //### Calc transport capacity
            double uc = 100.0*U*S; //in cm/s  in this formula assuming grad is SINE
            double ucr = 0.4;   // critical unit streampower in cm/s
            double cg = pow((D50->Drc+5)/0.32, -0.6);
            double dg = pow((D50->Drc+5)/300, 0.25);
            tc = ps * cg * pow(std::max(0.0, uc-ucr), dg); // kg/m3

        } else
            if(method == FSRIJN)
            {
                //https://www.leovanrijn-sediment.com/papers/Formulaesandtransport.pdf
                double ucr;
                //  double kinvis = 1e-6;
                if( d50m < 0.0005)
                    ucr = 0.19 * pow(d50m, 0.1) * log10(2.0* h/d50m); //p5
                else
                    ucr = 8.5  * pow(d50m, 0.6) * log10(2.0* h/d50m);

                double me = std::max((U - ucr)/sqrt(GRAV * d50m * (ps/pw - 1)),0.0); //p15
                double ds = d50m * 25296; //pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0)); // let op: D* = 25296*D50m! R2 = 1

                //        double qs = 0.008 * ps*U*d50m * pow(me, 2.4) * pow(ds, -0.6);
                double qs = 0.03 * ps*U*d50m * me*me * pow(ds, -0.6); // kg/s/m
                // van rijn 2007?, p 17, eq 6.4

                tc =  qs/ (U * h); // kg/m3   => WH or WHs

            }else
                if(method == FSRIJNFULL)
                {
                    /*
                 * D50	d50m	D*
                10	0.00001	0.25
                30	0.00003	0.76
                50	0.00005	1.26
                100	0.0001	2.53
                300	0.0003	7.59
                500	0.0005	12.65
                1000	0.001	25.30
                2000	0.002	50.59
                */
                    //van Rijn full (1984) following page 1632
                    double kinvis = 1e-6;
                    double ds = d50m * pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0));
                    double chezy = 18 * log10(4 * R/d90m);
                    double uc = U * sqrt(GRAV)/chezy;

                    //shields functions
                    double uscr = 0.055;
                    if(ds <150 && ds >= 20)
                        uscr = 0.013*pow(ds,0.29);
                    if(ds < 20 && ds >= 10)
                        uscr = 0.04*pow(ds,-0.10);
                    if(ds < 10 && ds >= 4)
                        uscr = 0.14*pow(ds,-0.64);
                    if(ds <4)
                        uscr = 0.24*pow(ds,-1);
                    uscr = sqrt(uscr * (ps/pw - 1)*GRAV * d50m);

                    double T = std::max(((uc*uc)/(uscr*uscr) - 1),0.0);  //transport stage parameter
                    double bsv = sqrt(GRAV * h * S); // bed shear velocity
                    double a = 0.1;  // half of the bedform height in m
                    double ca = 0.015 * (d50m/a) * pow(T,1.5)/pow(ds,0.3); //eq 38 reference concentration
                    double sv = SettlingVelocitySS->Drc;//GetSV(D50->Drc);

                    double beta = std::min(1.0 + 2.0*(sv/bsv)*(sv/bsv),5.0);
                    double powcb = 0.75; // not clear, between 0.4 and 1
                    double phi = 2.5 * pow(sv/bsv,0.8) * powcb;
                    double Z = sv/(beta*bsv*0.41); //suspension parameter, to do with upward turbulent versus gravity
                    double Zs = Z + phi;
                    double ad = 0.1; // else F is not valid!
                    double F = (pow(ad,Zs) - pow(ad,1.2))/(pow(1.0-ad,Zs)* (1.2 - Zs));
                    double qs =  F * U * h * ca;
                    tc = ps * qs/ (U * h);

                }else if(method == FSWUWANGJIA)
                {
                    double phk = 0;
                    double pek = 0;
                    double sv = settlingvelocities.at(_d);
                    double gd = graindiameters.at(_d)/1000000.0;
                    if (type == 0) {
                        FOR_GRAIN_CLASSES
                        {
                            //LET OP : RW_D and W_D !!!!
                            phk += RW_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                            pek += RW_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
                        }
                    } else {
                        FOR_GRAIN_CLASSES
                        {
                            //LET OP : RW_D and W_D !!!!
                            phk += W_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                            pek += W_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
                        }
                    }

                    double ppk = 1;
                    ppk = pow(phk/pek,0.6);
                    if(pek == 0 )
                    {
                        return 0;
                    }

                    double css = 0.03* (ps - pw) * (gd) * ppk;

                    double qs = 0.0000262 *pow(std::max(( pw * 0.01 * h * GRAV * S /css) - 1.0, 0.0)* U/(sqrt(sv)),2.2);
                    qs = qs * 1 * sqrt((ps/pw - 1)*GRAV*pow(gd,3.0));

                    tc = ps * qs/ (U * h);

                }
    return std::max(std::min(tc,MAXCONC ),0.0);
}
//--------------------------------------------------------------------------
/**
 * @fn double TWorld::calcTCBedload(int r,int c, int _d, int method, bool river)
 * @brief Calculates suspended layer transport capacity
 *
 * Calculates suspended load sediment transport capacity.
 * Based on govers, Van Rijn (simplified or full)
 * or Wu wang & Jia (for multiclass sediment).
 *
 * @param _d : The grain class (only needed when grain size distribution is used)
 * @param method : the TC method used
 */
double TWorld::calcTCBedload(int r,int c, int _d, int method, double h, double U, int type)
{
    double R,  hb, n, S, w;

    if (type == 0) {
        //    h = ChannelWH->Drc;
        hb = ChannelBLDepth->Drc;
        n = std::max(0.001, ChannelN->Drc);
        S = ChannelGrad->Drc;
        w = ChannelWidth->Drc;
        R = (w*h)/(2*h+w);
    } else
        if (type == 1) {
            //   h = hmx->Drc;
            hb = BLDepthFlood->Drc;
            n = std::max(0.001, N->Drc);
            S = Grad->Drc;
            w = ChannelAdj->Drc;
            R = (w*h)/(2*h+w);
        }

    //when water height is insignificant, transport capacity is zero
    //this is necessary since some of the used equations have strange behaviour
    //for these water heights or velocities. (h and v outside of valid range)
    if(h < MIN_HEIGHT || hb < MIN_HEIGHT)
        return 0;
    if(U < MIN_FLUX)
        return 0;

    double ps = 2650.0; //2400.0;
    double pw = 1000.0;
    double d50m = (D50->Drc/1000000.0);
    double d90m = (D90->Drc/1000000.0);
    double tc = 0;

    if(method == FSRIJN)
    {
        //Van rijn simplified (2007?)
        double ucr;
        if( d50m < 0.0005)
            ucr  = 0.19 * pow(d50m, 0.1) * log10(4.0*R/d90m);
        else
            ucr  = 8.5 * pow(d50m, 0.6) * log10(4.0*R/d90m);

        double me = std::max((U - ucr)/(sqrt(GRAV * d50m * ((ps/pw) - 1.0))),0.0);
        //        double qs = 0.005 * ps * U * h * pow(d50m/h,1.2) * pow(me, 2.4);
        double qs = 0.015 * ps*U*h * pow(d50m/h,1.2) * pow(me, 1.5); //eq 6.2
        // in kg/m/s /(m2/s) = kg/m3
        tc =  qs/ (U * hb);

    }else if(method == FSRIJNFULL)
    {
        //van Rijn full (1984)  see page 1450 1984_JHE_VanRijn_a.pdf
        double kinvis = 1e-6;

        double _dm = d90m; //d50m; interpretation -> assume all bedload particles are d90? Van RIjn deals mostly with sand

        double ds = _dm * pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0));
        double chezy = 18 * log(4 * h/d90m);  // h or hb or radius?
        double us = sqrt(GRAV) * U/chezy;

        // shield equations, full
        double uscr = 0.055;
        if(ds < 150 && ds >= 20)
            uscr = 0.013*pow(ds,0.29);
        if(ds < 20 && ds >= 10)
            uscr = 0.04*pow(ds,-0.10);
        if(ds < 10 && ds >= 4)
            uscr = 0.14*pow(ds,-0.64);
        if(ds <4)
            uscr = 0.24*pow(ds,-1);
        uscr = sqrt(uscr * (ps/pw - 1)*GRAV * _dm);  // effective bed shear velocity

        double T = std::max((us*us)/(uscr*uscr) - 1,0.0); // transport stage parameter
        double qs = 0.053 * (pow(T,2.1)/pow(ds,0.3)) * sqrt((ps/pw -1)*GRAV)*_dm*sqrt(_dm); // eq 22
        tc = ps * qs/ (U * hb);

    }else if(method == FSWUWANGJIA)
    {
        double na = (pow(graindiameters.at(_d)/100000.0,(1.0/6.0))/20.0)/n;
        double phk = 0;
        double pek = 0;
        if (type == 0) {
            FOR_GRAIN_CLASSES
            {
                //LET OP : RW_D and W_D !!!!
                phk += RW_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                pek += RW_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
            }
        } else {
            FOR_GRAIN_CLASSES
            {
                //LET OP : RW_D and W_D !!!!
                phk += W_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                pek += W_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
            }
        }
        double ppk = 1;
        ppk = pow(phk/pek,0.6);

        if(pek == 0 )
            return 0;

        double R = w*h/(2*h+w);
        double css = 0.03* (ps - pw) * (graindiameters.at(_d)/1000000.0) * ppk;

        double qs = 0.0053 *pow(std::max(pow(na,1.5)*((pw * R * GRAV * 0.1 * S/css)) - 1.0, 0.0),2.2);
        qs = qs * 1 * sqrt((ps/pw - 1)*GRAV*pow(graindiameters.at(_d)/1000000.0,3.0));

        tc = ps * qs/ (U * hb);

    }

    return std::max(std::min(tc,MAXCONCBL),0.0);
}



