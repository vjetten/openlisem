/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: https://github.com/vjetten/openlisem
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

//---------------------------------------------------------------------------
// deposit all sediment still in flow when infiltration causes WH to become minimum
// DOES NOT WORK, MB errors
void TWorld::cell_depositInfil(int r, int c)
{
    if (!SwitchErosion)
        return;

    if(SwitchKinematic2D == K2D_METHOD_DYN) {
        if(WH->Drc < 1e-6) {
            DepFlood->Drc -= SSFlood->Drc;
            SSFlood->Drc = 0;
            SSCFlood->Drc = 0;
            SSTCFlood->Drc = 0;
            if (SwitchUse2Phase) {
                DepFlood->Drc -= BLFlood->Drc;
                BLFlood->Drc = 0;
                BLCFlood->Drc = 0;
                BLTCFlood->Drc = 0;
            }
            Conc->Drc = 0; // after dynwave conc is sum of SS and BL!
        }
    } else {
        if (FloodDomain->Drc > 0) {
            if(hmx->Drc < 1e-6) {
                DepFlood->Drc -= SSFlood->Drc;
                SSFlood->Drc = 0;
                SSCFlood->Drc = 0;
            }
        } else {
            if(WH->Drc < 1e-6) {
                DEP->Drc -= Sed->Drc;
                Sed->Drc = 0;
                Conc->Drc = 0;
            }
        }
    }
}
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
double TWorld::MaxConcentration(double watvol, double sedvol)
{
    double conc = 0;
    if (watvol > 1e-6) {
        conc = std::min(sedvol/watvol, MAXCONC);   // 1e-6 is 1 ml/m2 !!
    }
    return conc;
}
//---------------------------------------------------------------------------
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
        double ds = dm * pow(1.65*GRAV/1e-12,0.333333);
        return SVCHCalibration*1e-6/dm*(ds*ds*ds)*pow(38.1+0.93*pow(ds,12.0/7.0), -7.0/8.0);
        //    // zhiyao et al, 2008
    } else {
        if(d < 100) {
            return SVCHCalibration*2*(2650.0-1000.0)*GRAV*pow(d/2000000.0, 2)/(9*0.001);
            //Stokes range settling velocity
        } else {
            double dm = d/1000.0;
            return SVCHCalibration*10.0 *sqrt(1.0 + 0.01 *(1.65* GRAV *dm*dm*dm )-1.0)/dm;
            //Settling velocity by Zanke (1977)
        }
    }

}
//---------------------------------------------------------------------------
void TWorld::cell_SplashDetachment(int r, int c)
{
    double _WH = FloodDomain->Drc == 0 ? WH->Drc : hmx->Drc;

    DETSplash->Drc = 0;

    if(_WH > HMIN && SplashStrength->Drc >= 0)
    {
        double DetDT1 = 0, DetDT2 = 0, DetLD1, DetLD2;
        double g_to_kg = 0.001;
        double Lc = Litter->Drc;
        double Cv = Cover->Drc;
        double strength = SplashStrength->Drc;
        double Int = Rain->Drc * 3600/_dt * 1000; // intensity in mm/h, Rain is in m
        double KE_DT = 0.0;
        double DETSplash_;               

        switch (KEequationType)
        {
            case KE_EXPFUNCTION: KE_DT = KEParamater_a1*(1-(KEParamater_b1*exp(-KEParamater_c1*Int))); break;
            case KE_LOGFUNCTION: KE_DT = (Int > 1 ? KEParamater_a2 + KEParamater_b2*log10(Int) : 0); break;
            case KE_POWERFUNCTION: KE_DT = KEParamater_a3*pow(Int, KEParamater_b3); break;
            // kin energy in J/m2/mm
        }
        //VJ 110706  KE equations

        double directrain = (1-Lc) * (1-Cv)*Rainc->Drc * 1000;
        // Added litter also to directrain, assume it covers the entire cell, not only under the plant
        // rainfall between plants in mm

        double KE_LD = std::max(15.3*sqrt(PlantHeight->Drc)-5.87, 0.0);
        // kin energy in J/m2/mm
        double throughfall = (1-Lc) * Cv * LeafDrain->Drc * 1000;
        // leaf drip in mm, is calculated as plant leaf drip in interception function so mult cover

        double WH0 = exp(-1.48*_WH*1000);
        // water buffer effect on surface, WH in mm in this empirical equation from Torri ?

        if(SwitchUseMaterialDepth)
        {
            double depdepth = std::max((StorageDep->Drc / BulkDens)/(_dx * DX->Drc),0.0);
            double fac1 = std::max(0.0,1.0 - depdepth/(SedimentMixingDepth->Drc+0.01));
            double fac2 = 1.0 - fac1;

            strength = strength * fac2 + (0.1033/DepositedCohesion) * fac1;
            //b = b * fac2 + 3.58 * fac1;
        }

        // fraction ponded area
        double FPA = 1.0;
        if (RR->Drc > 0.1)
            FPA =  1-exp(-1.875*(_WH/(0.01*RR->Drc)));

        // Between plants, directrain is already with 1-cover
        DetDT1 = g_to_kg * FPA*strength*KE_DT*WH0 * directrain;
        //ponded areas, kg/m2/mm * mm = kg/m2
        DetDT2 = _WH > 0 ? g_to_kg * (1-FPA)*strength*KE_DT * directrain * SplashDelivery: 0.0;
        //dry areas, kg/m2/mm * mm = kg/m2


        if (SwitchKETimebased)
        {
            if (directrain > 0)
            {
                DetDT1 = g_to_kg * FPA*strength*KE_DT*WH0 * _dt/3600;
                //ponded areas, kg/m2/sec * sec = kg/m2
                DetDT2 = g_to_kg * (1-FPA)*strength*KE_DT * _dt/3600 * SplashDelivery;
                //dry areas, kg/m2/sec * sec = kg/m2
            }
        }
        //based on work by Juan Sanchez

        // Under plants, throughfall is already with cover
        DetLD1 = g_to_kg * FPA*(strength*KE_LD)*WH0 * throughfall;
        //ponded areas, kg/m2/mm * mm = kg/m2
        DetLD2 = g_to_kg * (1-FPA)*(strength*KE_LD) * throughfall * SplashDelivery;
        //dry areas, kg/m2/mm * mm = kg/m2

        DETSplash_ = DetLD1 + DetLD2 + DetDT1 + DetDT2;
        // Total splash kg/m2

        // Deal with all exceptions:

        DETSplash_ *= (SoilWidthDX->Drc*DX->Drc);
        // kg/cell, only splash over soilwidth, not roads/hardsurfaces and channels ! houses re not in soilwisth, need to be done here
        // FROM KG/M2 TO KG/CELL

        DETSplash_ = (1-StoneFraction->Drc) * DETSplash_;
        // no splash on stone surfaces

        if (SwitchGrassStrip)
            DETSplash_ = (1-GrassFraction->Drc) * DETSplash_;

        //      if(SwitchSedtrap)
        //          DETSplash->Drc = (1-SedimentFilter->Drc) * DETSplash->Drc;
        // assume sedtrap can have splash

      //  if (SwitchHardsurface)
      //      DETSplash_ = (1-HardSurface->Drc)*DETSplash_;
        // no splash on hard surfaces ALREADY taken care of by soilwidth which excludes roads and hard surfaces

        if (SwitchHouses)
            DETSplash_ = (1-HouseCover->Drc)*DETSplash_;
        //is already contained in soilwidth
        // no splash from house roofs

        if (SwitchSnowmelt)
            DETSplash_ = (1-Snowcover->Drc)*DETSplash_;
        // no splash on snow deck

        if(SwitchUseMaterialDepth)
        {
            //check wat we can detach from the top and bottom layer of present material
            double dleft = DETSplash_;
            double deptake = 0;
            double mattake = 0;
            double detachment = 0;

            deptake = std::min(dleft,StorageDep->Drc);
            StorageDep->Drc -= deptake;
            // det not more than storage
            // decrease store depth

            detachment += deptake;
            // detachment is now taken material

            if(!(Storage->Drc < 0))
            {
                mattake = std::min(dleft,Storage->Drc);
                Storage->Drc -= mattake;

                detachment += mattake;
            }else
            {
                detachment += dleft;
            }
            DETSplash_ = detachment;
        }


        if(SwitchKinematic2D == K2D_METHOD_DYN) {
            SSFlood->Drc += DETSplash_;
            SSCFlood->Drc = MaxConcentration(WaterVolall->Drc, SSFlood->Drc);
        } else {
            if (FloodDomain->Drc > 0) {
                SSFlood->Drc += DETSplash_;
                SSCFlood->Drc = MaxConcentration(CHAdjDX->Drc * hmx->Drc, SSFlood->Drc);

            } else {
                Sed->Drc += DETSplash_;
                Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
            }
        }

        DETSplash->Drc = DETSplash_;
        // IN KG/CELL
    }
}
//---------------------------------------------------------------------------
/**
 * @fn double TWorld::SplashDetachment(double watvol, double sedvol)
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
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
            cell_SplashDetachment(r,c);
     }}
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

// Overland flow erosion for 1D flow only
void TWorld::cell_FlowDetachment(int r, int c)
{
    double erosionwh = WHrunoff->Drc;
    double erosionwv = WHrunoff->Drc*CHAdjDX->Drc;

    //transport capacity
    DETFlow->Drc = 0;
    DEP->Drc = 0;
    TC->Drc = calcTCSuspended(r,c,-1, FS_SS_Method, WHrunoff->Drc, V->Drc, 2);
    // trasnport capacity. 2 = kin wave. 1 = 2d flow and 0 is river

    if (erosionwh < HMIN) {
        if(DO_SEDDEP == 1) {
            DEP->Drc += -Sed->Drc;
            Sed->Drc = 0;
            Conc->Drc = 0;
            TC->Drc = 0;
        }
    } else {
        double maxTC = 0;
        double minTC = 0;

        double deposition = 0;
        double detachment = 0;
        double TransportFactor = 0;

        maxTC = std::max(TC->Drc - Conc->Drc,0.0);
        // positive difference: TC defi  cit becomes detachment (positive)
        minTC = std::min(TC->Drc - Conc->Drc,0.0);
        // negative difference: TC surplus becomes deposition (negative)
        // unit kg/m3

        if (minTC < 0) {

            //### deposition
            TransportFactor = (1-exp(-_dt*SettlingVelocitySS->Drc/erosionwh)) * erosionwv;
            // in m3
            // deposition can occur on roads and on soil (so use flowwidth)

            deposition = minTC * TransportFactor;
            // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0
            deposition = std::max(deposition, -Sed->Drc);

            if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                deposition = 0;
            // VJ 190325 prevent any activity on the boundary!

            if (SwitchSedtrap && SedMaxVolume->Drc == 0 && N->Drc == SedTrapN) {
                N->Drc = Norg->Drc;
            }
            // mannings N becomes normal when sedtrap is full

            if (SwitchSedtrap && SedMaxVolume->Drc > 0)
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
            }

            //add deposition to soil layer
            if (SwitchUseMaterialDepth)
                StorageDep->Drc += -deposition;

        } else
            if (maxTC > 0 && Y->Drc > 0) {
            //### detachment

            TransportFactor = _dt*SettlingVelocitySS->Drc * DX->Drc * SoilWidthDX->Drc;
            // soilwidth is erodible surface
            // TransportFactor = std::min(TransportFactor, Q->Drc*_dt);
            // detachment can only come from soil, not roads (so do not use flowwidth)
            // units s * m/s * m * m = m3

            detachment = maxTC * TransportFactor;//std::min(TransportFactor, erosionwv);
            // unit = kg/m3 * m3 = kg (/cell)

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

            if (SwitchHouses)
                detachment = (1-HouseCover->Drc)*detachment;
            // no flow det from house roofs

            detachment = (1-Snowcover->Drc) * detachment;
            /* TODO: CHECK THIS no flow detachment on snow */
            //is there erosion and sedimentation under the snowdeck?

            //detachment = DetachMaterial(r,c,1,false,false,false, detachment);
            // reacctivate when materiallayer is reinstalled
            detachment *= Y->Drc;

            if(Sed->Drc+detachment > MAXCONC * erosionwv)
                detachment = MAXCONC * erosionwv - Sed->Drc;
            // not more detachment then is possible to keep below diff(max concetrantion-sediment inf low)

            if (SwitchSedtrap && SedMaxVolume->Drc > 0)
                detachment = 0;
        } // minv > 0

        //### sediment balance
        // add to sediment in flow (IN KG/CELL)
        Sed->Drc += detachment;
        Sed->Drc += deposition;
        DETFlow->Drc += detachment;
        DEP->Drc += deposition;
        Conc->Drc = MaxConcentration(erosionwv, Sed->Drc);

    }

}
//---------------------------------------------------------------------------
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
    // when there is no usage of material depth, there is no deposited layer
    // actual erosion can then be calculated using the original erosion efficiency coefficient
    if(!SwitchUseMaterialDepth)
    {
        if(channel)
            return detachment *= ChannelY->Drc;
        else
            return detachment *= Y->Drc;
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

        //to remove small rounding errors that lead to negative values
        detachment = std::max(detachment,0.0);

        //check wat we can detache from the top and bottom layer of present material
        double dleft = detachment;
        double deptake = 0;
        double mattake = 0;
        detachment = 0;

        //take from the total storage
        deptake = std::min(dleft,RStorageDep->Drc);
        RStorageDep->Drc -= deptake;

        //add to the detachment what we have taken from the first soil layer
        detachment += deptake;

        //if the deposited layer is empty
        //use erosion efficiency of bottom layer again
        if(newY > 0)
            dleft *= ChannelY->Drc/newY;
        else
            dleft = 0;

        //bottom soil layer can be infinite
        if(!((RStorage->Drc) < -1))
        {
            //take from the total storage
            mattake = std::min(dleft,RStorage->Drc);
            RStorage->Drc -= mattake;
            //add to the detachment what we have taken from the second soil layer
            detachment += mattake;
        } else {
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
           // if(SedimentMixingDepth->Drc > MIN_HEIGHT)
                fac1 = std::max(0.0,1.0 - depdepth/SedimentMixingDepth->Drc);
            double fac2 = 1-fac1;
            //new erosion coefficient bases on soil layer mixinfactors
            double newY = Y->Drc * fac1 + fac2 * 1.0;

            //multiply potential detachment by erosion coefficient
            detachment = detachment *newY;

            //to remove small rounding errors that lead to negative values
            detachment = std::max(detachment,0.0);

            //check wat we can detach from the top and bottom layer of present material
            double dleft = detachment;
            double deptake = 0;
            double mattake = 0;
            detachment = 0;

            //take from the total storage
            deptake = std::min(dleft,StorageDep->Drc);
            StorageDep->Drc -= deptake;

            //add to the detachment what we have taken from the first soil layer
            detachment += deptake;

            //if the deposited layer is empty
            //use erosion efficiency of bottom layer again
            if(newY > 0)
                dleft *= Y->Drc/newY;
            else
                dleft = 0;

            if(!((Storage->Drc) < -1))
            {
                //take from the total storage
                mattake = std::min(dleft,Storage->Drc);
                Storage->Drc -= mattake;
                //add to the detachment what we have taken from the second soil layer
                detachment += mattake;
            } else {
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

            //to remove small rounding errors that lead to negative values
            detachment = std::max(detachment,0.0);
            //check wat we can detache from the top and bottom layer of present material
            double dleft = detachment;
            double deptake = 0;
            double mattake = 0;
            detachment = 0;

            //take from the total storage
            deptake = std::min(dleft,StorageDep->Drc);
            StorageDep->Drc -= deptake;

            //add to the detachment what we have taken from the first soil layer
            detachment += deptake;
            //if the deposited layer is empty
            //use erosion efficiency of bottom layer again
            if(newY > 0)
                dleft *= Y->Drc/newY;
            else
                dleft = 0;

            //bottom soil layer can be infinite
            if(!((Storage->Drc) < -1))
            {
                //take from the total storage
                mattake = std::min(dleft,Storage->Drc);
                Storage->Drc -= mattake;
                //add to the detachment what we have taken from the second soil layer
                detachment += mattake;
            } else {
                //all left potential detachment is added to detachment (infinite soil layer)
                detachment += dleft;
            }

            //finally, return the total detachment
            return std::max(0.0,detachment);
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
    /*
    //simple iteration over maps
    double wtotal = 0;
    FOR_GRAIN_CLASSES
    {
        wtotal += (*M).Drcd;
    }
    return wtotal;
    */
        return 0;
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
//---------------------------------------------------------------------------
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
/*
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
        } */
    }
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
    double R=0, hs=0, S = 0, w = 0, man = 0.01;
    double d50m;

    if (type == 0) {
        // river
        d50m = D50CH->Drc/1000000.0;
        hs = ChannelSSDepth->Drc;
        S = ChannelGrad->Drc;
        w = ChannelWidth->Drc;
        R = (w*h)/(2*h+w);
        man = ChannelN->Drc;

    } else
        if (type == 1) {
            // flood
            d50m = D50->Drc/1000000.0;
            hs = SSDepthFlood->Drc;
            S = Grad->Drc;
            w = ChannelAdj->Drc;
            R = (w*h)/(2*h+w);
        } else
            if (type == 2) {
                // kin wave
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

    double ps = 2650.0;
    double pw = 1000.0;
    double tc = 0;

    if(method == FSHAIRSINEROSE)
    {
        double om =  U*S;
        double omcr = 0.004;
        tc =  d50m/SettlingVelocitySS->Drc* 0.013/GRAV * 1.650 * std::max(0.0, om - omcr)/h ;
        //    m/ (m/s)* m/s /m  dimensionless?

    } else
        if(method == FSGOVERS)
        {
            //### Calc transport capacity
            double uc = 100.0*U*S; //in cm/s  in this formula assuming grad is SINE
            double ucr = 0.4;   // critical unit streampower in cm/s
            double cg = cgovers->Drc;//pow((d50m+5)/0.32, -0.6);
            double dg = dgovers->Drc;//pow((d50m+5)/300, 0.25);
            tc = ps * cg * pow(std::max(0.0, uc-ucr), dg); // kg/m3

        } else
            if(method == FSRIJN)
            {
                //https://www.leovanrijn-sediment.com/papers/Formulaesandtransport.pdf
                //double kinvis = 1e-6;
                //double Ds = d50m * pow(1.65*GRAV/1e-12,1.0/3.0); //dimensionless sed size
                //cr,suspension = 0.3/(1+ D*) + 0.1 [1-exp(-0.05D*)
                // critical shields parameter
                //double cs = 0.3/(1+Ds)+0.1*(1-exp(-0.05*Ds));
                // Ucritical, suspension= 5.75 [log(12h/(6D50))] [cr,suspension (s-1) g D50]0.5

                // double ucr;
               // ucr = 5.75*(log10(12*hs/(6.0*d50m))*qSqrt(cs*1.650*GRAV*d50m); //1.650 = s-1
//                if( d50m < 0.0005) // 500 mu, so always the first one!
//                    ucr = 0.19 * pow(d50m, 0.1) * log10(2.0* h/d50m); //p5
//                else
//                    ucr = 8.5  * pow(d50m, 0.6) * log10(2.0* h/d50m);
//                // set in motion critical U
                double Ds = d50m * 25296; //pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0)); // let op: D* = 25296*D50m! R2 = 1
              //  double shields = 0.3/(1+Ds)+0.1*(1-exp(-0.05*DS));
              //  double ucr = 5.75*log10(2*h/d50m)*sqrt(shields*1.65*GRAV*d50m);

                double ucr = UcrCHCalibration*2.8*pow(h/d50m,0.1)*sqrt(1.65*GRAV*d50m);
//                // suspension critical U

               //page 5

                double me = std::max(0.0,U - ucr)/sqrt(GRAV * d50m * 1.65);
                //p15 mobility parameter

                double qs = 0.03 * ps*U*d50m * me*me * pow(Ds, -0.6); // kg/s/m
                // van rijn 2007?, p 17, eq 6.4
                ChannelQsr->Drc = qs;
                //qDebug() << qs;
                tc =  qs/ (U * h); //kg/s/m / (m2/s) =  kg/m3   => WH or WHs

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
                    double ds = d50m * pow(1.65*GRAV/(kinvis*kinvis),(1.0/3.0)); //
                    //double chezy = 18 * log10(4 * R/d90m);
                    double chezy = 1/man*pow(R,1/6);
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
                    uscr = UcrCHCalibration*sqrt(uscr * 1.65*GRAV * d50m);

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
                } else
                    if(method == FENGELUND)
                    {
                        //https://www.hec.usace.army.mil/confluence/rasdocs/rassed1d/1d-sediment-transport-technical-reference-manual/computing-transport-capacity/sediment-transport-potential/engelund-hansen
                        double U2 = U*U;
                        double C = 1/man*pow(R,1/6);  //Chezy
                        //double Tb = pw*U2*GRAV/(C*C);
                        //double shields = Tb/((ps-pw)*GRAV*d50m);
                        double shields = U2/(C*C*1.65*d50m);
                        double qs = U2*pow(shields, 1.5)*sqrt(d50m/(GRAV*1.65));

                        //http://ponce.sdsu.edu/onlineengelundhansen.php
                        //double shields1 = h*S/(1.65*d50m);
                        //double qs = 0.001*(U2/(2*GRAV*S*h))*pow(shields1,2.5)*dw*sqrt(1.65*GRAV*d50m*d50m*d50m);
                        //0.001 is ton to kg
                        // gives almost the same value
                        ChannelQsr->Drc = qs;

                        tc =  qs/ (U * h); //kg/s/m / (m2/s) =  kg/m3   => WH or WHs
                }else if(method == FSWUWANGJIA)
                {
                        // NOT USED, FOR MULTIPLE GRAINSIZES
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
    if (type == 0) {
        d90m = D90CH->Drc/1000000.0;
        d50m = D50CH->Drc/1000000.0;
    }

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



