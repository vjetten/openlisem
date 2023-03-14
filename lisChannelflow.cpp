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
 \file lisChannelflow.cpp
 \brief Channel hydrology and sediment detachment and movement processes.

functions: \n
- void TWorld::CalcVelDischChannel(void) calculate Velocity, alpha and Q in the channel \n
- void TWorld::ChannelFlow(void) calculate channelflow, ChannelDepth, do kinematic wave \n
*/

//#include <algorithm>
#include "model.h"
//#include "operation.h"

//---------------------------------------------------------------------------
void TWorld::ChannelFlowandErosion()
{
    if (!SwitchIncludeChannel)
        return;

    SwitchChannelKinWave = true;    // set to false for experimental swof in channel

    ChannelRainandInfil();          // subtract infil, add rainfall

    ChannelBaseflow();              // calculate baseflow

    ChannelFlow();                  // channel kin wave for water and sediment

    ChannelFlowDetachmentNew();     // detachment, deposition for SS and BL

    ChannelSedimentFlow();

}
//---------------------------------------------------------------------------
void TWorld::ChannelBaseflow(void)
{
//    if(!SwitchChannelBaseflow)
//        return;

    // add a stationary part
    if(SwitchChannelBaseflow && SwitchChannelBaseflowStationary)
    {
        // first time
        if(!addedbaseflow) {
           #pragma omp parallel for num_threads(userCores)
           FOR_ROW_COL_MV_CHL {
                ChannelWaterVol->Drc += BaseFlowInitialVolume->Drc;
           }}

           addedbaseflow = true;
        }

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelWaterVol->Drc += BaseFlowInflow->Drc * _dt;
        }}
    }

    if (SwitchChannelBaseflow && (SwitchSWATGWflow || SwitchExplicitGWflow))
        GroundwaterFlow();
    // move groundwater and add baseflow to channel

}
//---------------------------------------------------------------------------
void TWorld::ChannelRainandInfil(void)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {

        if (SwitchCulverts && ChannelMaxQ->Drc > 0)
            ChannelWaterVol->Drc += 0;
        else
            ChannelWaterVol->Drc += Rainc->Drc*ChannelWidth->Drc*DX->Drc;
        // add rainfall to channel, assume no interception

        // subtract infiltration, no infil in culverts
        if (SwitchChannelInfil && ChannelMaxQ->Drc <= 0) {
            double inf = ChannelDX->Drc * ChannelKsat->Drc*_dt/3600000.0 * (ChannelWidth->Drc + 2.0*ChannelWH->Drc/cos(atan(ChannelSide->Drc)));
            // hsat based through entire wet cross section
            inf = std::min(ChannelWaterVol->Drc, inf);
            // cannot be more than there is
            ChannelWaterVol->Drc -= inf;
            ChannelInfilVol->Drc += inf;
        }        
    }}

}
//---------------------------------------------------------------------------
//! calc channelflow, ChannelDepth, kin wave
//! channel WH and V and Q are clculated before
void TWorld::ChannelFlow(void)
{
    if (SwitchChannelKinwaveDt) {
        if (_dt_user > _dtCHkin) {
            double n = _dt_user/_dtCHkin;
            _dt = _dt_user/n;
        }
    }

    for (double t = 0; t < _dt_user; t+=_dt)
    {

        //double sumvol = getMassCH(ChannelWaterVol);

        // velocity, alpha, Q
        #pragma omp parallel num_threads(userCores)
        FOR_ROW_COL_MV_CHL {

            // calc velocity and Q
            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);

            double MaxQ = ChannelMaxQ->Drc;
            double wh = ChannelWH->Drc;
            double ChannelQ_ = 0;
            double ChannelV_ = 0;
            double ChannelAlpha_ = 0;

            // calc channel V and Q, using original width
            if (wh > 1e-6) {
                double Perim, Radius, Area;
                double sqrtgrad = std::max(sqrt(ChannelGrad->Drc), 0.001);
                double N = ChannelN->Drc;

                double FW = ChannelWidth->Drc;
                double FWO = ChannelWidthO->Drc;

                Perim = (FW + 2.0*wh);
                Area = FW*wh;

                if (SwitchChannelAdjustCHW) {
                    double whn = wh * (FW/FWO);
                    Perim = FWO + whn*2; //original dimensions, wider than cell size
                    Area = FWO * whn;
                    // shallow width perim Area
                }
                Perim *= ChnTortuosity;
                Radius = (Perim > 0 ? Area/Perim : 0);
                ChannelV_ = std::min(_CHMaxV,std::pow(Radius, 2.0/3.0)*sqrtgrad/N);
                ChannelQ_ = ChannelV_ * Area;
                ChannelAlpha_ = Area/std::pow(ChannelQ_, 0.6);
                ChannelNcul->Drc  = ChannelN->Drc;

                if (SwitchCulverts) {
                    if (ChannelMaxQ->Drc > 0 ) {

                        ChannelNcul->Drc = (0.05+ChannelQ_/MaxQ) * 0.015; //0.015 is assumed to be the N of a concrete tube
                        //https://plainwater.com/water/circular-pipe-mannings-n/
                        // resistance increases with discharge, tube is getting fuller

                        double v2 = std::pow(Radius, 2.0/3.0)*sqrtgrad/ChannelNcul->Drc;
                        //max velocity not to exceed MaxQ, see excel
                        ChannelV_ = std::min(_CHMaxV,std::min(ChannelV_, v2));
                        ChannelNcul->Drc = std::min(ChannelNcul->Drc,ChannelN->Drc);
                        ChannelQ_ = ChannelV_ * Area;

                        if (ChannelQ_ > MaxQ){
                            ChannelV_ = MaxQ/Area;
                            ChannelQ_ = MaxQ;
                        }
                        ChannelAlpha_ = Area/std::pow(ChannelQ_, 0.6);
                    }
                }

                ChannelAlpha->Drc = ChannelAlpha_;
                ChannelQ->Drc = ChannelQ_;
                ChannelV->Drc = ChannelV_;

                ChannelQsn->Drc = 0;
                Channelq->Drc = 0;
                QinKW->Drc = 0;
            }
        }}

        // ChannelV and Q and alpha now based on original width and depth, channel vol is always the same

        if (SwitchLinkedList) {

            ChannelQn->setAllMV();

            FOR_ROW_COL_LDDCH5 {
                Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelAlpha, ChannelDX);
            }}
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (pcr::isMV(ChannelQn->Drc))
                    ChannelQn->Drc = 0;
            }}

        } else {
            // default
            KinematicExplicit(crlinkedlddch_, ChannelQ, ChannelQn, ChannelAlpha, ChannelDX);
        }

        // calc V and WH back from Qn (original width and depth)
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            //  ChannelQ->Drc = ChannelQn->Drc;  // NOT because needed in erosion!

            ChannelWaterVol->Drc += (QinKW->Drc - ChannelQn->Drc)*_dt;
            ChannelWaterVol->Drc = std::max(0.0,ChannelWaterVol->Drc);
            // vol is previous + in - out
            ChannelAlpha->Drc = ChannelQn->Drc > 1e-6 ? (ChannelWaterVol->Drc/ChannelDX->Drc)/std::pow(ChannelQn->Drc, 0.6) : ChannelAlpha->Drc;

            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
            // new channel WH, use adjusted channelWidth

            double ChannelArea = ChannelWaterVol->Drc/ChannelDX->Drc;
            double P = 2*ChannelWH->Drc+ChannelWidth->Drc;

            if (P > 0)
                ChannelV->Drc = std::pow(ChannelArea/P,2/3)*sqrtGrad->Drc/ChannelNcul->Drc;
            else
                ChannelV->Drc = 0;

            // get the maximum for output
            maxChannelflow->Drc = std::max(maxChannelflow->Drc, ChannelQn->Drc);
            maxChannelWH->Drc = std::max(maxChannelWH->Drc, ChannelWH->Drc);
        }}
       // correctMassBalanceCH(sumvol,ChannelWaterVol);
    }
    _dt=_dt_user;
}

void TWorld::ChannelSedimentFlow()
{
    if (!SwitchErosion)
        return;

    //double sumvol = getMassCH(ChannelWaterVol);

    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        double concss = MaxConcentration(ChannelWaterVol->Drc, &ChannelSSSed->Drc, &ChannelDep->Drc);
        ChannelQSSs->Drc = ChannelQ->Drc * concss; // m3/s *kg/m3 = kg/s
      //  ChannelQSSs->Drc = ChannelQsr->Drc*ChannelQ_; //kg/m/s *m
    }}

    if(SwitchUse2Phase) {
        #pragma omp parallel num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            double concbl = MaxConcentration(ChannelWaterVol->Drc, &ChannelBLSed->Drc, &ChannelDep->Drc);
            ChannelQBLs->Drc = ChannelQ->Drc * concbl;
        }}
    }


    if (SwitchLinkedList) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            pcr::setMV(ChannelQSSsn->Drc);
        }}

        FOR_ROW_COL_LDDCH5 {
              routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn, ChannelAlpha, ChannelDX, ChannelSSSed);
        }}

        if(SwitchUse2Phase) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                pcr::setMV(ChannelQBLsn->Drc);
            }}

            FOR_ROW_COL_LDDCH5 {
                routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn, ChannelAlpha, ChannelDX, ChannelBLSed);
            }}
        }

    } else {
            //NOTE: this is the new channel alpha, not good!
        KinematicSubstance(crlinkedlddch_, LDDChannel, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn, ChannelAlpha, ChannelDX, ChannelSSSed);
        if(SwitchUse2Phase) {
            KinematicSubstance(crlinkedlddch_, LDDChannel, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn, ChannelAlpha, ChannelDX, ChannelBLSed);
        }
    }


    if (SwitchIncludeRiverDiffusion) {
        RiverSedimentDiffusion(_dt, ChannelSSSed, ChannelSSConc);
        // note SSsed goes in and out, SSconc is recalculated inside
    }

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        RiverSedimentLayerDepth(r,c);
        RiverSedimentMaxC(r,c);
        ChannelQsn->Drc = ChannelQSSsn->Drc;
        ChannelSed->Drc = ChannelSSSed->Drc;
    }}

    if(SwitchUse2Phase) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelQsn->Drc += ChannelQBLsn->Drc;
            ChannelSed->Drc += ChannelBLSed->Drc;
        }}
    }
}




/* not used */
double TWorld::getMassCH(cTMap *M)
{
    double sum2 = 0;
   // #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        sum2 += M->Drc;
    }}
    return sum2;
}
/* not used */
void TWorld::correctMassBalanceCH(double sum1, cTMap *M)
{
    double sum2 = 0;
   // double n = 0;

  //  #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        sum2 += M->Drc;
  //      n += 1;
    }}
    // total and cells active for M
    double dhtot = fabs(sum2) > 0 ? (sum1 - sum2)/sum2 : 0;

    if (dhtot > 0) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            M->Drc = M->Drc*(1.0 + dhtot);            // <- distribution weighted to h
            M->Drc = std::max(M->Drc , 0.0);
        }}
    }
}
