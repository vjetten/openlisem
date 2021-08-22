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
**  website, information and code: http://lisem.sourceforge.net
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
    if (!SwitchIncludeChannel)
        return;

    if (!SwitchErosion)
        return;

#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {

        RiverSedimentLayerDepth(r,c);
        //creates ChannelBLDepth and ChannelSSDepth, if 1 layer ChannelBLDepth = 0

        double sswatervol = ChannelSSDepth->Drc*DX->Drc*ChannelWidth->Drc;
        double ssdischarge = ChannelV->Drc * ChannelWidth->Drc * ChannelSSDepth->Drc;
        double blwatervol = 0;
        double bldischarge = 0;
        if (SwitchUse2Phase) {
            blwatervol = ChannelBLDepth->Drc*DX->Drc*ChannelWidth->Drc;
            bldischarge = ChannelV->Drc * ChannelWidth->Drc * ChannelBLDepth->Drc;
        }

        //get transport capacity for bed/suspended load for a specific cell and grain size class
        double SSTC = 0;
        double BLTC = 0;
        if (SwitchUse2Phase)
            BLTC = calcTCBedload(r, c, 1, R_BL_Method, ChannelWH->Drc,ChannelV->Drc, 0);
        SSTC = calcTCSuspended(r, c, 1, R_SS_Method, ChannelWH->Drc, ChannelV->Drc, 0);

        //find transport capacity for bed and suspended layer

        ChannelSSTC->Drc = SSTC;
        double SSDepth = ChannelSSDepth->Drc;
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
        if(ChannelWH->Drc < MIN_HEIGHT)
        {
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

        } else {
            //there is water

            //### do suspended first

            //deposition
            maxTC = std::max(SSTC - SSC, 0.0);  // TC in kg/m3
            minTC = std::min(SSTC - SSC, 0.0);

            if (minTC < 0) {
                TransportFactor = (1-exp(-_dt*TSettlingVelocitySS/ChannelWH->Drc)) * sswatervol;

                deposition = std::max(TransportFactor * minTC,-SS); // in kg
                // not more than SS present

            } else {

                //  detachment
                TransportFactor = _dt*TSettlingVelocitySS * ChannelDX->Drc * ChannelWidth->Drc;
                TransportFactor = std::min(TransportFactor, ssdischarge*_dt);
                // use discharge because standing water has no erosion

                //NB ChannelWidth and ChannelWidth the same woith rect channel
                detachment = maxTC * TransportFactor;
                // cannot have more detachment than remaining capacity in flow

                detachment = DetachMaterial(r,c,1,true,false,false, detachment);
                // multiply by Y

                if(MAXCONC * sswatervol < SS + detachment)
                    detachment = std::max(0.0, MAXCONC * sswatervol - SS);
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

                } else {
                    //there is BL

                    maxTC = std::max(BLTC - BLC,0.0);
                    minTC = std::min(BLTC - BLC,0.0);

                    if (maxTC > 0) {
                        //### detachment
                        TransportFactor = _dt*TSettlingVelocityBL * ChannelDX->Drc * ChannelWidth->Drc;
                        TransportFactor = std::min(TransportFactor, bldischarge*_dt);

                        // units s * m/s * m * m = m3
                        detachment =  maxTC * TransportFactor;
                        // unit = kg/m3 * m3 = kg

                        detachment = DetachMaterial(r,c,1,true,false,true, detachment);
                        // mult by Y and mixingdepth
                        // IN KG/CELL

                        if(MAXCONC * blwatervol < BL+detachment)
                            detachment = std::max(0.0, MAXCONC * blwatervol - BL);
                        // limit detachment to what BLtemp can carry
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

                        BL += detachment; //ChannelBLSed
                        BL += deposition;
                        ChannelBLSed->Drc = BL;
                        ChannelDep->Drc += deposition;
                        ChannelDetFlow->Drc += detachment;
                        ChannelTC->Drc += ChannelBLTC->Drc;
                        ChannelSed->Drc += BL;
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
void TWorld::ChannelFlowDetachment()
{
    if (!SwitchIncludeChannel)
        return;

    if (!SwitchErosion)
        return;

#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        RiverSedimentLayerDepth(r,c);
        //creates ChannelBLDepth and ChannelSSDepth, if 1 layer ChannelBLDepth = 0

        double sswatervol = ChannelSSDepth->Drc*DX->Drc*ChannelWidth->Drc;
        double ssdischarge = ChannelV->Drc * ChannelWidth->Drc * ChannelSSDepth->Drc;
        double blwatervol = 0;
        double bldischarge = 0;
        if (SwitchUse2Phase) {
            blwatervol = ChannelBLDepth->Drc*DX->Drc*ChannelWidth->Drc;
            bldischarge = ChannelV->Drc * ChannelWidth->Drc * ChannelBLDepth->Drc;
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
                if (SwitchUse2Phase)
                    TBLTCtemp = ChannelBLTC;
                TSSTCtemp = ChannelSSTC;
            }
            //            else
            //            {
            //                TBLTCtemp = RBLTC_D.at(d);
            //                TSSTCtemp = RSSTC_D.at(d);
            //            }

            //get transport capacity for bed/suspended load for a specific cell and grain size class
            if (SwitchUse2Phase)
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
                if (SwitchUse2Phase) {
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

                if (SwitchUse2Phase) {
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
                    TransportFactor = _dt*TSettlingVelocitySS * ChannelDX->Drc * ChannelWidth->Drc;
                    TransportFactor = std::min(TransportFactor, ssdischarge*_dt);
                    //TransportFactor = ssdischarge*_dt;
                    // use discharge because standing water has no erosion

                    //NB ChannelWidth and ChannelWidth the same woith rect channel
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
                if (SwitchUse2Phase) {
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
                            TransportFactor = _dt*TSettlingVelocityBL * ChannelDX->Drc * ChannelWidth->Drc;
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

            //                RBLC_D.Drcd = MaxConcentration(ChannelWidth->Drc*DX->Drc*RBLD_D.Drcd, &RBL_D.Drcd, &ChannelDep->Drc);
            //                sssmax = MAXCONCBL * DX->Drc *ChannelWidth->Drc*RBLD_D.Drcd;
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

    //fromChannelWHtoVol(r, c);
    ChannelWaterVol->Drc = ChannelWidth->Drc * ChannelWH->Drc * ChannelDX->Drc;

    double frac = ChannelSSDepth->Drc/ChannelWH->Drc;
    //maximum concentration
    if(!SwitchUseGrainSizeDistribution)
    {
        _SSC->Drc = MaxConcentration(ChannelWaterVol->Drc*frac, &_SS->Drc, &ChannelDep->Drc);
        if (SwitchUse2Phase)
            _BLC->Drc = MaxConcentration(ChannelWaterVol->Drc*(1-frac), &_BL->Drc, &ChannelDep->Drc);
    }

    ChannelSed->Drc = (SwitchUse2Phase ? _BL->Drc : 0) + _SS->Drc;
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
#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        _SSC->Drc = MaxConcentration(ChannelWaterVol->Drc, &_SS->Drc, &ChannelDep->Drc);
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
    double eta = eddyvs/FS_SigmaDiffusion;

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
    _SSC->Drc = MaxConcentration(ChannelWaterVol->Drc, &_SS->Drc, &ChannelDep->Drc);
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
}
