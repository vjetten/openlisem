/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
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
void TWorld::ChannelFlowDetachment()
{
    if (!SwitchErosion)
        return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        double velocityfactor = 1.0;
        RiverSedimentLayerDepth(r,c);
        //creates ChannelBLDepth and ChannelSSDepth, if 1 layer ChannelBLDepth = 0

        double sswatervol = ChannelSSDepth->Drc*DX->Drc*ChannelWidth->Drc;
        double blwatervol = 0;
        double SS = ChannelSSSed->Drc;
        double BL = 0;
        ChannelDetFlow->Drc = 0; // reset to zero, this is the det and dep of this timestep
        ChannelDep->Drc = 0;
        double deposition = 0;
        double detachment = 0;
        double TransportFactor = 0;
        double maxTC = 0;
        double minTC = 0;

        //get transport capacity for bedload for a specific cell and grain size class
        if (SwitchUse2Phase) {
            BL = ChannelBLSed->Drc;
            blwatervol = ChannelBLDepth->Drc*DX->Drc*ChannelWidth->Drc;
            velocityfactor = ChannelBLDepth->Drc/ChannelWH->Drc;
            ChannelBLTC->Drc = calcTCBedload(r, c, R_BL_Method, ChannelWH->Drc, ChannelWidth->Drc, ChannelV->Drc*velocityfactor, SUSPchannel);
            //transport capacity for bedload based on D90
        }

        ChannelSSTC->Drc = calcTCSuspended(r, c, R_SS_Method, ChannelWH->Drc, ChannelWidth->Drc, ChannelV->Drc, SUSPchannel);
        //transport capacity for suspended matter based on D50

        //=== Do suspended matter first

        //when waterheight is insignificant, deposite all remaining sediment, HMIN = 1e-6 m
        if(ChannelWH->Drc < HMIN) {
            if(DO_SEDDEP == 1) {
                deposition += -SS;
                ChannelSSConc->Drc = 0;
                ChannelSSSed->Drc = 0;
                ChannelSSTC->Drc = 0;
                ChannelDep->Drc += deposition;
            }
        } else {
            //there is water

            //### do suspended first

            maxTC = qMax(ChannelSSTC->Drc - ChannelSSConc->Drc, 0.0);  // TC in kg/m3
            minTC = qMin(ChannelSSTC->Drc - ChannelSSConc->Drc, 0.0);

            if (minTC < 0) {
                //deposition

                if (SwitchDepositionLinear)
                    TransportFactor =  _dt*SettlingVelocitySS->Drc * ChannelDX->Drc * ChannelWidth->Drc;
                else
                    TransportFactor = (1-exp(-_dt*SettlingVelocitySS->Drc/ChannelWH->Drc)) * sswatervol;

                deposition = qMax(TransportFactor * minTC,-SS); // in kg
                // not more than SS present

            } else {
                //  detachment
                if(maxTC > 0 && ChannelCohesion->Drc >= 0) {
                    TransportFactor = _dt*SettlingVelocitySS->Drc * ChannelDX->Drc * ChannelWidth->Drc;

                    detachment = ChannelY->Drc * maxTC * TransportFactor;
                    //DetachMaterial(r,c,1,true,false,false, detachment);
                    // multiply by Y

                    if (SwitchCulverts && ChannelCulvert->Drc > 0)
                        detachment = 0;
                    // no detahcment in culverts

                    if(SS + detachment > MAXCONC * sswatervol)
                        detachment = qMax(0.0, MAXCONC * sswatervol - SS);

                } else {
                    detachment = 0;
                }
            }

            //### sediment balance add suspended
            SS += detachment;
            SS += deposition;
            ChannelSSSed->Drc = SS;
            ChannelDep->Drc += deposition;
            ChannelDetFlow->Drc += detachment;
            ChannelTC->Drc = ChannelSSTC->Drc;

            // if (SwitchUseMaterialDepth)
            //     RStorageDep->Drc += -deposition;

            //### do bedload
            if (SwitchUse2Phase) {
                // water height very low
                if(ChannelBLDepth->Drc < MIN_HEIGHT) {
                    if(DO_SEDDEP == 1) {
                        ChannelDep->Drc += -BL;
                        ChannelBLTC->Drc = 0;
                        ChannelBLConc->Drc = 0;
                        ChannelBLSed->Drc = 0;
                    }
                } else {
                    // there is water

                    maxTC = qMax(ChannelBLTC->Drc - ChannelBLConc->Drc,0.0);
                    minTC = qMin(ChannelBLTC->Drc - ChannelBLConc->Drc,0.0);

                    if (maxTC > 0 && ChannelCohesion->Drc >= 0) {
                        //### detachment
                        TransportFactor = _dt*SettlingVelocityBL->Drc * ChannelDX->Drc * ChannelWidth->Drc;
                        // units s * m/s * m * m = m3
                        detachment = maxTC * qMin(TransportFactor, maxTC*sswatervol);
                        // unit = kg/m3 * m3 = kg

                        detachment *= ChannelY->Drc;//DetachMaterial(r,c,1,true,false,true, detachment);
                        // mult by Y and mixingdepth
                        // IN KG/CELL

                        if(BL + detachment > MAXCONC * blwatervol)
                            detachment = MAXCONC * blwatervol - BL;

                    } else {
                        //### deposition
                        //if (ChannelBLDepth->Drc > MIN_HEIGHT)
                        TransportFactor = (1-exp(-_dt*SettlingVelocityBL->Drc/ChannelBLDepth->Drc)) * blwatervol;

                        // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0
                        deposition = qMax(minTC * TransportFactor, -BL);
                        // cannot have more depo than sediment present
                        BL += detachment;
                        BL += deposition;
                        ChannelBLSed->Drc = BL;
                        //ChannelSed->Drc += BL;
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
void TWorld::ChannelDetachmentContinuous()
{
    if (!SwitchErosion)
        return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        double velocityfactor = 1.0;
        RiverSedimentLayerDepth(r,c);
        //creates ChannelBLDepth and ChannelSSDepth, if 1 layer ChannelBLDepth = 0

        double sswatervol = ChannelSSDepth->Drc*DX->Drc*ChannelWidth->Drc;
        double blwatervol = 0;
        double SS = ChannelSSSed->Drc;
        double BL = 0;
        ChannelDetFlow->Drc = 0; // reset to zero, this is the det and dep of this timestep
        ChannelDep->Drc = 0;
        double deposition = 0;
        double detachment = 0;
        double TransportFactor = 0;
        double maxTC = 0;
        double minTC = 0;

        //get transport capacity for bedload for a specific cell and grain size class
        if (SwitchUse2Phase) {
            BL = ChannelBLSed->Drc;
            blwatervol = ChannelBLDepth->Drc*DX->Drc*ChannelWidth->Drc;
            velocityfactor = ChannelBLDepth->Drc/ChannelWH->Drc;
            ChannelBLTC->Drc = calcTCBedload(r, c, R_BL_Method, ChannelWH->Drc, ChannelWidth->Drc, ChannelV->Drc*velocityfactor, SUSPchannel);
            //transport capacity for bedload based on D90
        }

        ChannelSSTC->Drc = calcTCSuspended(r, c, R_SS_Method, ChannelWH->Drc, ChannelWidth->Drc, ChannelV->Drc, SUSPchannel);
        //transport capacity for suspended matter based on D50

        //=== Do suspended matter first

        //when waterheight is insignificant, deposite all remaining sediment, HMIN = 1e-6 m
        if(ChannelWH->Drc < HMIN) {
            if(DO_SEDDEP == 1) {
                deposition += -SS;
                ChannelSSConc->Drc = 0;
                ChannelSSSed->Drc = 0;
                ChannelSSTC->Drc = 0;
                ChannelDep->Drc += deposition;
            }
        } else {
            //there is water

            //### do suspended first

            maxTC = qMax(ChannelSSTC->Drc - ChannelSSConc->Drc, 0.0);  // TC in kg/m3

            //deposition

            if (SwitchDepositionLinear)
                TransportFactor =  qMin(1.0, _dt*SettlingVelocitySS->Drc * ChannelDX->Drc * ChannelWidth->Drc);
            else
                TransportFactor = (1-exp(-_dt*SettlingVelocitySS->Drc/ChannelWH->Drc)) * sswatervol;

            deposition = TransportFactor * ChannelSSConc->Drc; // in kg

            //  detachment
            if(maxTC > 0 && ChannelCohesion->Drc >= 0) {

                if (SwitchDepositionLinear)
                    TransportFactor =  qMin(1.0, ChannelY->Drc * _dt*SettlingVelocitySS->Drc * ChannelDX->Drc * ChannelWidth->Drc);
                else
                    TransportFactor = (1-exp(-ChannelY->Drc * _dt*SettlingVelocitySS->Drc/ChannelWH->Drc)) * sswatervol;

                detachment = maxTC * TransportFactor;
                //DetachMaterial(r,c,1,true,false,false, detachment);
                // multiply by Y

                if (SwitchCulverts && ChannelCulvert->Drc > 0)
                    detachment = 0;
                // no detahcment in culverts

                if(SS + detachment > MAXCONC * sswatervol)
                    detachment = qMax(0.0, MAXCONC * sswatervol - SS);

            }

            //### sediment balance add suspended
            SS += detachment;
            SS += deposition;
            ChannelSSSed->Drc = SS;
            ChannelDep->Drc += deposition;
            ChannelDetFlow->Drc += detachment;
            ChannelTC->Drc = ChannelSSTC->Drc;

            // if (SwitchUseMaterialDepth)
            //     RStorageDep->Drc += -deposition;

            //### do bedload
            if (SwitchUse2Phase) {
                // water height very low
                if(ChannelBLDepth->Drc < MIN_HEIGHT) {
                    if(DO_SEDDEP == 1) {
                        ChannelDep->Drc += -BL;
                        ChannelBLTC->Drc = 0;
                        ChannelBLConc->Drc = 0;
                        ChannelBLSed->Drc = 0;
                    }
                } else {
                    // there is water

                    //### deposition
                    TransportFactor = (1-exp(-_dt*SettlingVelocityBL->Drc/ChannelBLDepth->Drc)) * blwatervol;

                    // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0
                    deposition = ChannelBLConc->Drc * TransportFactor;
                    // cannot have more depo than sediment present

                    maxTC = qMax(ChannelBLTC->Drc - ChannelBLConc->Drc,0.0);

                    if (maxTC > 0 && ChannelCohesion->Drc >= 0) {
                        //### detachment
                        TransportFactor = _dt*SettlingVelocityBL->Drc * ChannelDX->Drc * ChannelWidth->Drc;
                        // units s * m/s * m * m = m3
                        detachment = maxTC * qMin(TransportFactor, maxTC*sswatervol);
                        // unit = kg/m3 * m3 = kg

                        detachment *= ChannelY->Drc;//DetachMaterial(r,c,1,true,false,true, detachment);
                        // mult by Y and mixingdepth
                        // IN KG/CELL

                        if(BL + detachment > MAXCONC * blwatervol)
                            detachment = MAXCONC * blwatervol - BL;

                        BL += detachment;
                        BL += deposition;
                        ChannelBLSed->Drc = BL;
                        //ChannelSed->Drc += BL;
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
    //ChannelSed->Drc = (SwitchUse2Phase ? ChannelBLSed->Drc : 0) + ChannelSSSed->Drc;
    double sed = (SwitchUse2Phase ? ChannelBLSed->Drc : 0) + ChannelSSSed->Drc;
    //total concentration
    ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, sed);//ChannelSed->Drc);
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
            ldd = static_cast <int>(LDDChannel->data[rt][ct]);
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
    int ldd = static_cast <int>(LDDChannel->Drc);
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

    double dux = qMax(dux1,dux2);

    //diffusion coefficient according to J.Smagorinski (1964)
    double eddyvs = cdx * dux;
    //and devide by turbulent prandtl-smidth number, def 1.0
    double eta = eddyvs/R_SigmaDiffusion;

    //add diffusive fluxes to previous cell in channel.
    if(foundp)
    {
        Real coeff = ChannelSSDepth->data[r][c] > 0 ? dt*eta *qMin(1.0,ChannelSSDepth->data[rp][cp]/ChannelSSDepth->data[r][c]) : 0.0;
        coeff = qMin(coeff, courant_factor);

        _SS->data[rp][cp] += coeff * _SS->Drc;
        _SS->data[r][c] -= coeff * _SS->Drc;
    }

    //add diffusive fluxes to next cell in channel.
    if(foundn)
    {
        Real coeff = ChannelSSDepth->data[r][c] > 0 ? dt*eta *qMin(1.0,ChannelSSDepth->data[rn][cn]/ChannelSSDepth->data[r][c]) : 0.0;
        coeff = qMin(coeff, courant_factor);

        _SS->data[rn][cn] += coeff  * _SS->Drc;
        _SS->data[r][c] -= coeff  * _SS->Drc;
    }
}

    //recalculate concentrations
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        //set concentration from present sediment
        _SS->Drc = qMax(0.0,_SS->Drc);
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
        ChannelBLDepth->Drc = qMin(qMin(d50m * 1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), factor*ChannelWH->Drc), 0.1);
        ChannelSSDepth->Drc = qMax(ChannelWH->Drc - ChannelBLDepth->Drc,0.0);

  //  }
}
