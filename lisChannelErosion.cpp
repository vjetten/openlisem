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
  \file lisChannelErosion.cpp
  \brief Flow detachment functions for channels

functions: \n
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
//#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
//    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )


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
 * @see TWorld:SwitchUse2Phase
 * @see TWorld:DetachMaterial
 *
 */
void TWorld::ChannelFlowDetachmentNew()
{
    if (!SwitchErosion)
    return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {

        RiverSedimentLayerDepth(r,c);
        //creates ChannelBLDepth and ChannelSSDepth, if 1 layer ChannelBLDepth = 0

        double sswatervol = ChannelSSDepth->Drc*DX->Drc*ChannelWidth->Drc;
        double blwatervol = 0;
        if (SwitchUse2Phase) {
            blwatervol = ChannelBLDepth->Drc*DX->Drc*ChannelWidth->Drc;
        }

        //get transport capacity for bed/suspended load for a specific cell and grain size class
        double SSTC = 0;
        double BLTC = 0;
        if (SwitchUse2Phase)
            BLTC = calcTCBedload(r, c, 1, R_BL_Method, ChannelWH->Drc,ChannelV->Drc, 0);
        SSTC = calcTCSuspended(r, c, 1, R_SS_Method, ChannelWH->Drc, ChannelV->Drc, 0);

        //find transport capacity for bed and suspended layer

        ChannelSSTC->Drc = SSTC;
        //double SSDepth = ChannelSSDepth->Drc;
        double SSC = ChannelSSConc->Drc;
        double SS = ChannelSSSed->Drc;
        double TSettlingVelocitySS = SettlingVelocitySS->Drc;

        double BLDepth = 0;
        double BLC = 0;
        double BL = 0;
        double TSettlingVelocityBL = 0;

        if (SwitchUse2Phase) {
            BLDepth = ChannelBLDepth->Drc;
            TSettlingVelocityBL = SettlingVelocityBL->Drc;
            BLC = ChannelBLConc->Drc;
            BL = ChannelBLSed->Drc;
            ChannelBLTC->Drc = BLTC;
        }

        ChannelDetFlow->Drc = 0;
        ChannelDep->Drc = 0;

        double deposition = 0;
        double detachment = 0;
        double TransportFactor = 0;
        double maxTC = 0;
        double minTC = 0;

        //when waterheight is insignificant, deposite all remaining sediment
        if(ChannelWH->Drc < HMIN)
        {
            if(DO_SEDDEP == 1) {
                deposition += -SS;
                ChannelSSConc->Drc = 0;
                ChannelSSSed->Drc = 0;
                ChannelSSTC->Drc = 0;

                if (SwitchUse2Phase) {
                    deposition = -BL;
                    ChannelBLConc->Drc = 0;
                    ChannelBLSed->Drc = 0;
                    ChannelBLTC->Drc = 0;
                }

                ChannelSed->Drc = 0;
                ChannelDep->Drc += deposition;
            }
        } else {
            //there is water

            //### do suspended first

            //deposition
            maxTC = std::max(SSTC - SSC, 0.0);  // TC in kg/m3
            minTC = std::min(SSTC - SSC, 0.0);

            if (minTC < 0) {
                TransportFactor = (1-exp(-_dt*TSettlingVelocitySS/ChannelWH->Drc)) * sswatervol;
             //   TransportFactor = _dt*TSettlingVelocitySS * ChannelDX->Drc * ChannelWidth->Drc;
                //TransportFactor = std::min(TransportFactor, ssdischarge * _dt);

                deposition = std::max(TransportFactor * minTC,-SS); // in kg
                // not more than SS present

            } else
                //  detachment
                if(maxTC > 0 && ChannelY->Drc > 0){

                    TransportFactor = _dt*TSettlingVelocitySS * ChannelDX->Drc * ChannelWidth->Drc;
                    //TransportFactor = std::min(TransportFactor, ssdischarge*_dt);
                    // use discharge because standing water has no erosion

                    //NB ChannelWidth and ChannelWidth the same woith rect channel
                    detachment = maxTC * std::min(TransportFactor, sswatervol);
                    // cannot have more detachment than remaining capacity in flow

                    detachment *= ChannelY->Drc;//DetachMaterial(r,c,1,true,false,false, detachment);
                    // multiply by Y

                    if(SS + detachment > MAXCONC * sswatervol)
                        detachment = MAXCONC * sswatervol - SS;
                }

            //### sediment balance add suspended
            SS += detachment;
            SS += deposition;
            ChannelSSSed->Drc = SS;
            ChannelDep->Drc += deposition;
            ChannelDetFlow->Drc += detachment;
            ChannelTC->Drc = ChannelSSTC->Drc;
            ChannelSed->Drc = SS;
            //total transport capacity (bed load + suspended load), used for output

            if (SwitchUseMaterialDepth)
                RStorageDep->Drc += -deposition;

            //### do bedload
            if (SwitchUse2Phase) {

                if(BLDepth < MIN_HEIGHT) {
                    ChannelDep->Drc += -BL;
                    ChannelBLTC->Drc = 0;
                    ChannelBLConc->Drc = 0;
                    ChannelBLSed->Drc = 0;

                } else
                  if (ChannelY->Drc > 0){
                    //there is BL

                    maxTC = std::max(BLTC - BLC,0.0);
                    minTC = std::min(BLTC - BLC,0.0);

                    if (maxTC > 0 && ChannelY->Drc > 0) {
                        //### detachment
                        TransportFactor = _dt*TSettlingVelocityBL * ChannelDX->Drc * ChannelWidth->Drc;
                        //TransportFactor = std::min(TransportFactor, bldischarge*_dt);
                        // units s * m/s * m * m = m3
                        detachment = maxTC * std::min(TransportFactor, maxTC*sswatervol);
                        // unit = kg/m3 * m3 = kg

                        detachment = ChannelY->Drc;//DetachMaterial(r,c,1,true,false,true, detachment);
                        // mult by Y and mixingdepth
                        // IN KG/CELL

                        if(BL + detachment > MAXCONC * blwatervol)
                            detachment = MAXCONC * blwatervol - BL;

                    } else {
                        //### deposition
                        if (BLDepth > MIN_HEIGHT)
                            TransportFactor = (1-exp(-_dt*TSettlingVelocityBL/BLDepth)) * blwatervol;
                        else
                            TransportFactor =  1.0 * blwatervol;

                        // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0
                        deposition = std::max(minTC * TransportFactor, -BL);
                        // cannot have more depo than sediment present

                        if (SwitchUseMaterialDepth)
                            RStorageDep->Drc += -deposition;

                        BL += detachment;
                        BL += deposition;
                        ChannelBLSed->Drc = BL;
                        ChannelSed->Drc += BL;
                        ChannelDep->Drc += deposition;
                        ChannelDetFlow->Drc += detachment;
                        ChannelTC->Drc += ChannelBLTC->Drc;
                        //total transport capacity (bed load + suspended load), used for output
                    }
                }
            }
        }

        RiverSedimentMaxC(r,c);
        //partial and total concentration ALL DONE
    }}
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

    double frac = ChannelSSDepth->Drc/ChannelWH->Drc;
    //maximum concentration
    if(!SwitchUseGrainSizeDistribution)
    {
        _SSC->Drc = MaxConcentration(ChannelWaterVol->Drc*frac, _SS->Drc);
        if (SwitchUse2Phase)
            _BLC->Drc = MaxConcentration(ChannelWaterVol->Drc*(1-frac), _BL->Drc);
    }

    ChannelSed->Drc = (SwitchUse2Phase ? _BL->Drc : 0) + _SS->Drc;
    //total concentration
    ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc);
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
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        _SSC->Drc = MaxConcentration(ChannelWaterVol->Drc, _SS->Drc);
    }}

    //diffusion of Suspended Sediment layer
    //#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CH {

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
        double coeff = ChannelSSDepth->data[r][c] > 0 ? dt*eta *std::min(1.0,ChannelSSDepth->data[rp][cp]/ChannelSSDepth->data[r][c]) : 0.0;
        coeff = std::min(coeff, courant_factor);

        _SS->data[rp][cp] += coeff * _SS->Drc;
        _SS->data[r][c] -= coeff * _SS->Drc;
    }

    //add diffusive fluxes to next cell in channel.
    if(foundn)
    {
        double coeff = ChannelSSDepth->data[r][c] > 0 ? dt*eta *std::min(1.0,ChannelSSDepth->data[rn][cn]/ChannelSSDepth->data[r][c]) : 0.0;
        coeff = std::min(coeff, courant_factor);

        _SS->data[rn][cn] += coeff  * _SS->Drc;
        _SS->data[r][c] -= coeff  * _SS->Drc;
    }
}

    //recalculate concentrations
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        //set concentration from present sediment
        _SS->Drc = std::max(0.0,_SS->Drc);
        _SSC->Drc = MaxConcentration(ChannelWaterVol->Drc, _SS->Drc);
    }}
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
    if (!SwitchUse2Phase) {
        ChannelSSDepth->Drc = ChannelWH->Drc;
        // ChannelBLDepth->Drc = 0;
        return;
    }

    double ps = 2650;
    double pw = 1000;
    double factor = 0.5;
    double R = (ChannelWidth->Drc * ChannelWH->Drc)/(ChannelWH->Drc * 2 + ChannelWidth->Drc);

  //  if(!SwitchUseGrainSizeDistribution)
  //  {
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

  //  }
}
