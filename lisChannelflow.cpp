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
 \file lisChannelflow.cpp
 \brief Channel hydrology and sediment detachment and movement processes.

functions: \n
- void TWorld::CalcVelDischChannel(void) calculate Velocity, alpha and Q in the channel \n
- void TWorld::ChannelFlow(void) calculate channelflow, ChannelDepth, do kinematic wave \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"

#define cell(r,c,a,b,e,d) qDebug()<<a->data[r][c]<<b->data[r][c]<<e->data[r][c]<<d->data[r][c]

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
    if(!SwitchChannelBaseflow)
        return;

    // add a stationary part
    if(SwitchChannelBaseflowStationary)
    {
        if(!addedbaseflow) {
           #pragma omp parallel for num_threads(userCores)
           FOR_ROW_COL_MV_CHL {
                ChannelWaterVol->Drc += BaseFlowInitialVolume->Drc;
           }}

           addedbaseflow = true;
        }
    }

    // GW recharge and GW outflow
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double CellArea_ = CellArea->Drc;

        //=== GW recharge
        Perc->Drc = cell_Percolation(r, c, GW_recharge); // in m
        double GWrecharge = Perc->Drc * CellArea_; // m3
        // GW recharge same principle as percolation, in m3

         //=== bypass flow
        //obsolete
//        double bpflow = 0;
//        if (GW_bypass > 0 && Lw->Drc > GW_bypass && Lw->Drc < SoilDepth1->Drc) {
//            bpflow = Lw->Drc * GW_bypass * (Poreeff->Drc-Thetaeff->Drc);
//            Lw->Drc *= (1-GW_bypass);
//        }
//        if (GW_bypass > 0 && Lw->Drc > SoilDepth1->Drc+GW_bypass) {
//            double dL = std::min(SoilDepth1->Drc, Lw->Drc * GW_bypass);
//            bpflow = dL * (ThetaS2->Drc-ThetaI2->Drc);
//            Lw->Drc -= dL;
//        }
//        GWbp->Drc = bpflow*CellArea_;


        //=== lateral GW outflow
        // ksat is already in m per timestep
        double pore, ksat;
        if (SwitchTwoLayer) {
            pore = ThetaS2->Drc;
            ksat = Ksat2->Drc;
        } else {
            pore = Poreeff->Drc;
            ksat = Ksateff->Drc;
        }

        double GWVol_ = GWVol->Drc;//outflow m3
        double wh = GWVol_/CellArea_/pore;
        double GWout_ = GW_flow * CellArea_ * ksat * BaseflowL->Drc; // m3 volume out from every cell
        GWout_ = wh > GW_threshold ?  GWout_ * (wh - GW_threshold)*pore * (1-exp(-6*wh)) : 0.0;
        //m3:  GW_flow*ksat*dt * ((dx/L)^b) *crosssection of flow dh*dx; //*porosity
       // GWout_ *= (1+Grad->Drc);

        //  GWout_ = GWout_ * wh/pore;
        // stop outflow when some minimum GW level, 2.4.2.10 in SWAT
        // decay function exp(-6 * GW WH) for smooth transition

        // ==== update GW level
        GWout_ = std::min(GWout_, GWVol_+GWrecharge);
        // cannot be more than there is
        GWVol_ = GWVol_  + GWrecharge - GWout_; //m3
        //update GW volume

        GWout->Drc = GWout_;
        GWVol->Drc = GWVol_;
        GWWH->Drc = GWVol_/CellArea_/pore;  //for display

        tma->Drc = ChannelWidth->Drc > 0 ? Qbin->Drc : 0;// prev timestep Qbin

        Qbin->Drc = 0;

        // change soildepth2 with GW changes
        if (GWWH->Drc > 0) {
            if (SwitchTwoLayer) {
                double dh = std::max(0.1,SoilDepth2init->Drc - GWWH->Drc);
                SoilDepth2->Drc = dh;
                GWWH->Drc = SoilDepth2init->Drc - dh;
            } else {
                double dh = std::max(0.1,SoilDepth1init->Drc - GWWH->Drc);
                SoilDepth1->Drc = dh;
                GWWH->Drc = SoilDepth1init->Drc - dh;
            }
            GWVol->Drc = pore*GWWH->Drc * CellArea_;
        }
        GWWHmax->Drc = std::max(GWWHmax->Drc, GWWH->Drc);

    }}
    report(*GWWHmax,"gwmax.map");
//report(*GWWH,"gwwh");
    // new qbin
    AccufluxGW(crlinkedlddbase_, GWout, Qbin, ChannelWidth);
    // LDDbase, Qin, Qout, chanwidth used as flag, move the gw flow to the channel,
    // Qbin is inflow to the channel from the surrounding cells in m3 per timestep

    double factor = exp(-GW_lag);
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        Qbin->Drc *= ChannelWidth->Drc/_dx;
        // do this or not? for very small channel a lot of water is added but what haoppens to the rest
        Qbase->Drc = Qbin->Drc*(1-factor) + tma->Drc*factor;  //m3 added per timestep, for MB
        ChannelWaterVol->Drc += Qbin->Drc*(1-factor) + tma->Drc*factor;
        // flow according to SWAT 2009, page 174 manual, eq 2.4.2.8

        if (SwitchChannelBaseflowStationary)
            ChannelWaterVol->Drc += BaseFlowInflow->Drc * _dt;
    }}

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

        // subtract infiltration
        if (SwitchChannelInfil) {
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
                double sqrtgrad = sqrt(ChannelGrad->Drc);
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

                if (sqrtgrad > MIN_SLOPE) {
                    //ChannelAlpha->Drc = std::pow(N/sqrtgrad * std::pow(Perim, 2.0/3.0), 0.6);
                    //ChannelQ->Drc = std::pow(Area/ChannelAlpha->Drc, 1.0/0.6);
                    ChannelV_ = std::min(10.0,std::pow(Radius, 2.0/3.0)*sqrtgrad/N);
                    ChannelQ_ = ChannelV_ * Area;
                    if (SwitchCulverts) {
                        if (MaxQ > 0 && ChannelQ_ > MaxQ){
                            ChannelN->Drc = ChannelQ_/MaxQ *ChannelN->Drc;
                            ChannelV_ = MaxQ/Area;
                            ChannelQ_ = MaxQ;
                        }
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

        }}

        // ChannelV and Q and alpha now based on original width and depth, channel vol is always the same

        if (SwitchLinkedList) {

            ChannelQn->setAllMV();

            FOR_ROW_COL_LDDCH5 {
                Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, Channelq, ChannelAlpha, ChannelDX, ChannelMaxQ);
            }}
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (pcr::isMV(ChannelQn->Drc))
                    ChannelQn->Drc = 0;
            }}

        } else {
            KinematicExplicit(crlinkedlddch_, ChannelQ, ChannelQn, Channelq, ChannelAlpha, ChannelDX, ChannelMaxQ);
        }

        // calc V and WH back from Qn (original width and depth)
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {

//            if (SwitchCulverts) {
//                double MaxQ = ChannelMaxQ->Drc;
//                if (MaxQ > 0 && ChannelQn->Drc > MaxQ)
//                    ChannelQn->Drc = MaxQ;
//            }

            double chqn = ChannelQn->Drc;
            ChannelWaterVol->Drc += (QinKW->Drc - chqn)*_dt;
           // ChannelQn->Drc = std::min(ChannelQn->Drc, ChannelWaterVol->Drc/_dt);
            ChannelQ->Drc = chqn;
            ChannelAlpha->Drc = chqn > 1e-6 ? (ChannelWaterVol->Drc/ChannelDX->Drc)/std::pow(chqn, 0.6) : 0.0;
        }}
        //water vol from mass balance, includes any errors

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
            // new channel WH, use adjusted channelWidth

            double ChannelArea = ChannelWaterVol->Drc/ChannelDX->Drc;
            double P = 2*ChannelWH->Drc+ChannelWidth->Drc;

            if (P > 0)
                ChannelV->Drc = std::pow(ChannelArea/P,2/3)*sqrtGrad->Drc/ChannelN->Drc;
            else
                ChannelV->Drc = 0;

        }}


        // get the maximum for output
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL
        {
            maxChannelflow->Drc = std::max(maxChannelflow->Drc, ChannelQn->Drc);
            maxChannelWH->Drc = std::max(maxChannelWH->Drc, ChannelWH->Drc);
        }}

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
