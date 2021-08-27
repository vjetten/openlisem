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
void TWorld::ChannelAddBaseandRain(void)
{
    if (!SwitchIncludeChannel)
        return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {

        if (SwitchCulverts && ChannelMaxQ->Drc > 0)
            ChannelWaterVol->Drc += 0;
        else
            ChannelWaterVol->Drc += Rainc->Drc*ChannelWidth->Drc*DX->Drc;

        // subtract infiltration
        if (SwitchChannelInfil) {
            double inf = ChannelDX->Drc * ChannelKsat->Drc*_dt/3600000.0 * (ChannelWidth->Drc + 2.0*ChannelWH->Drc/cos(atan(ChannelSide->Drc)));
            inf = std::min(ChannelWaterVol->Drc, inf);
            ChannelWaterVol->Drc -= inf;
            ChannelInfilVol->Drc += inf;
        }

        if(SwitchChannelBaseflow)
        {
//            if(!addedbaseflow)
//            {
//                ChannelWaterVol->Drc += BaseFlowInitialVolume->Drc;
//                //BaseFlowTot += BaseFlowInitialVolume->Drc;
//            }
//            ChannelWaterVol->Drc += BaseFlowInflow->Drc * _dt;
//            //BaseFlowTot += BaseFlowInflow->Drc * _dt;

        }

        ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
// deeper wh because of adjested width?

    }}
    if (!addedbaseflow)
        addedbaseflow = true;

    if(SwitchChannelBaseflow) {

        if(SwitchTwoLayer) {
            // calculate GW recharge from every cell = Percolation
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                // double FC2 = 0.5 * 0.7867*exp(-0.012*Ksat2->Drc)*ThetaS2->Drc;
                //double FC2 = (-0.115*log(Ksat2->Drc) + 0.611)*ThetaS2->Drc; // wilting point
                //FC2 = 0.025*ThetaS2->Drc;
                // recharge if soilmoisture > FC2 (wilting point?)
                double pore = ThetaS2->Drc;
                double thetar = 0.025 * pore;
                double theta = ThetaI2->Drc;
                double SoilDep2 = SoilDepth2->Drc;
                double SoilDep1 = SoilDepth1->Drc;
                double Percolation = GW_recharge * Ksat2->Drc*_dt/3600000 * pow((theta-thetar)/(pore-thetar), bca2->Drc);
                double dL = 0;
                if (Lw->Drc > SoilDep1)
                    dL = SoilDep2 - Lw->Drc;
                else
                    dL = SoilDep2 - SoilDep1;
                double moisture = std::max(0.0, dL*(theta-thetar));
                Percolation = std::min(Percolation, moisture*0.9);
                moisture = moisture - Percolation;
                ThetaI2->Drc = std::min(pore, moisture /std::max(0.01,dL) + thetar);
                GWrec->Drc = ThetaI2->Drc > thetar ? Percolation*CellArea->Drc : 0.0;   //m3
            }}
        } else {
                #pragma omp parallel for num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    double pore = Poreeff->Drc;
                    double thetar = 0.025 * pore;
                    double theta = Thetaeff->Drc;
                    double SoilDep1 = SoilDepth1->Drc;
                    double Percolation = GW_recharge * Ksateff->Drc*_dt/3600000 * pow((theta-thetar)/(pore-thetar), bca1->Drc);
                    double dL = SoilDep1 - Lw->Drc;
                    double moisture = std::max(0.0, dL*(theta-thetar));
                    Percolation = std::min(Percolation, moisture*0.9);
                    moisture = moisture - Percolation;
                    Thetaeff->Drc = std::min(pore, moisture /std::max(0.01,dL) + thetar);
                    GWrec->Drc = Thetaeff->Drc > thetar ? Percolation*CellArea->Drc : 0.0;   //m3
                }}
            }

            // calc gw flow out and update GW vol and WH
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                WHGW->Drc = VolGW->Drc/CellArea->Drc;

                tma->Drc = GW_flow * (Ksat2->Drc*_dt/3600000.0)*BaseflowL->Drc * (WHGW->Drc*_dx*ThetaS2->Drc);
                // ksat * ((dx/L)^b) *crosssection of flow dh*dx
                tma->Drc = std::min(tma->Drc, VolGW->Drc+GWrec->Drc);

              //  if (ChannelQn->Drc < 0.001)
                //    tma->Drc *= std::min(1.0,ChannelQn->Drc/1000);
                //tma->Drc *= 1-1/(1+pow(10*hmxWH->Drc,5.0));
                //ChannelQb->Drc = ChannelQn->Drc;

                VolGW->Drc = VolGW->Drc + GWrec->Drc - tma->Drc; //m3
                WHGW->Drc = VolGW->Drc/CellArea->Drc;
            }}

        Accuflux(crlinkedlddbase_, tma, Qbin);
        //move the gw flow to the channel, Qb is inflow to the channel from the surrounding cellsm in m3 per timestep

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            Qbin->Drc *= ChannelWidth->Drc/_dx;

           // =(1-MIN(1,1/((D4+0.01)/$N$4^$N$3)))
//            ChannelWaterVol->Drc += Qbin->Drc * (1-std::min(1.0,1/pow((ChannelQn->Drc+1)/50,2)));
//            if (ChannelQn->Drc < 1)
//                Qbin->Drc = 0;
            if (ChannelQn->Drc < 1)
                Qbin->Drc *= ChannelQn->Drc/1;



            ChannelWaterVol->Drc += Qbin->Drc;


//            //double Area = (VolQb->Drc+1.0)/_dx;
//            double Area = std::max(VolQb->Drc/ChannelWidth->Drc,ChannelWidth->Drc*0.1); //assume min 0.5 m baseflow
//            tma->Drc = Qbin->Drc/_dt;
//            tmb->Drc = Area/std::pow(Qbin->Drc/_dt, 0.6); //alpha
//            //ChannelAlpha_ = Area/std::pow(ChannelQ_, 0.6);
//            tmc->Drc = 0;  //q
//            tmd->Drc = -1; //qmax
        }}

//        KinematicExplicit(crlinkedlddch_, tma, Qbase, tmc, tmb, ChannelDX, tmd);

//        #pragma omp parallel for num_threads(userCores)
//        FOR_ROW_COL_MV_CHL {
//            Qbase->Drc = std::min(Qbase->Drc, VolQb->Drc/_dt + QinKW->Drc);
//            VolQb->Drc = VolQb->Drc + QinKW->Drc*_dt - Qbase->Drc*_dt ;
//        }}

        cell(135,246,WHGW,GWrec,ChannelQn,Qbin);
        //  report(*tma,"tma");
        //  report(*Qbin,"qbin");
//    report(*tma,"gwout");
//    report(*tmb,"rech");
//   // report(*ThetaI2,"ti");
 //       report(*Qbase,"qb");
    }
}



//---------------------------------------------------------------------------
//! calc channelflow, ChannelDepth, kin wave
//! channel WH and V and Q are clculated before
void TWorld::ChannelFlow(void)
{
    if (!SwitchIncludeChannel)
        return;

   // double sumch = getMassCH(ChannelWH);

    // velocity, alpha, Q
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_CHL {

        // calc velocity and Q
        ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
//deeper width
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

            Perim = FW + 2.0*wh;
            Area = FW*wh;

            if (SwitchChannelAdjustCHW) {
                double whn = wh * (FW/FWO);
                Perim = FWO + whn*2; //original dimensions, wider than cell size
                Area = FWO * whn;
                // shallow width perim Area
            }

            Radius = (Perim > 0 ? Area/Perim : 0);

            if (sqrtgrad > MIN_SLOPE) {
                //ChannelAlpha->Drc = std::pow(N/sqrtgrad * std::pow(Perim, 2.0/3.0), 0.6);
                //ChannelQ->Drc = std::pow(Area/ChannelAlpha->Drc, 1.0/0.6);
                ChannelV_ = std::pow(Radius, 2.0/3.0)*sqrtgrad/N;
                ChannelQ_ = ChannelV_ * Area;
                if (SwitchCulverts) {
                    if (MaxQ > 0 && ChannelQ_ > MaxQ){
                        ChannelQ_ = MaxQ;
                        ChannelV_ = MaxQ/Area;
                    }
                }
                ChannelAlpha_ = Area/std::pow(ChannelQ_, 0.6);
            }
            //ChannelQb->Drc = 0.01 * ChannelQ->Drc;
        }

        ChannelAlpha->Drc = ChannelAlpha_;
        ChannelQ->Drc = ChannelQ_;
        ChannelV->Drc = ChannelV_;

        ChannelQsn->Drc = 0;
        Channelq->Drc = 0;
        QinKW->Drc = 0;

        if (SwitchErosion) {
            double concss = MaxConcentration(ChannelWaterVol->Drc, &ChannelSSSed->Drc, &ChannelDep->Drc);
            ChannelQSSs->Drc = ChannelQ_ * concss;

            if(SwitchUse2Phase) {
                double concbl = MaxConcentration(ChannelWaterVol->Drc, &ChannelBLSed->Drc, &ChannelDep->Drc);
                ChannelQBLs->Drc = ChannelQ_ * concbl;
            }
        }
    }}

    // ChannelV and Q and alpha now based on original width and depth, channel vol is always the same

    if (SwitchChannelKinWave) {
        int loop = 1;
//        if (SwitchChannelKinwaveDt) {
//            if (_dt > _dtCHkin) {
//                loop = int(_dt/_dtCHkin);
//                _dt = _dtCHkin;
//            }
//            if (SwitchChannelKinwaveAvg) {
//                fill(*tma,0);
//                fill(*tmb,0);
//            }
//         //Faverage dynamix
//          //  qDebug() << loop;
//        }

        for (int i = 0; i < loop; i++) {


            if (SwitchLinkedList) {

                ChannelQn->setAllMV();
                // route water 1D and sediment
                FOR_ROW_COL_LDDCH5 {
                    Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, Channelq, ChannelAlpha, ChannelDX, ChannelMaxQ);
                }}
                //cover(*ChannelQn, *LDD, 0);
                #pragma omp parallel for num_threads(userCores)
                FOR_ROW_COL_MV_L {
                    if (pcr::isMV(ChannelQn->Drc))
                        ChannelQn->Drc = 0;
                }}

            } else {
                KinematicExplicit(crlinkedlddch_, ChannelQ, ChannelQn, Channelq, ChannelAlpha, ChannelDX, ChannelMaxQ);
            }

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_CHL {
                //                if (SwitchChannelKinwaveDt && SwitchChannelKinwaveAvg) {
                //                   tma->Drc += ChannelQn->Drc;
                //                   tmb->Drc += QinKW->Drc;
                //                }
                ChannelQ->Drc = ChannelQn->Drc;
            }}
        } //loop
//        _dt = _dt_user;

//        if (SwitchChannelKinwaveDt && SwitchChannelKinwaveAvg) {
//            #pragma omp parallel for num_threads(userCores)
//            FOR_ROW_COL_MV_CHL {
//                ChannelQn->Drc = tma->Drc/(double) loop;
//                QinKW->Drc = tmb->Drc/(double) loop;
//            }}
//        }

        // calc V and WH back from Qn (original width and depth)
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelWaterVol->Drc += (Channelq->Drc + QinKW->Drc)*_dt;
            ChannelQn->Drc = std::min(ChannelQn->Drc, ChannelWaterVol->Drc/_dt);
            ChannelWaterVol->Drc -= ChannelQn->Drc*_dt ;
            //water vol from mass balance, includes any errors

            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
            // new channel WH, use adjusted channelWidth

            double ChannelArea = ChannelWaterVol->Drc/ChannelDX->Drc;
            ChannelV->Drc = (ChannelArea > 0 ? ChannelQn->Drc/ChannelArea : 0);

            if (SwitchCulverts) {
                //TO DO ?????????????????????
                double MaxQ = ChannelMaxQ->Drc;
                if (MaxQ > 0 && ChannelQn->Drc > MaxQ) {
                    ChannelQn->Drc = MaxQ;
                    ChannelV->Drc = MaxQ/ChannelArea;
                   // ChannelWaterVol->Drc = ChannelArea*ChannelDX->Drc;
                }
             }
        }}

    } else {

        ChannelSWOFopen();

    }


//    correctMassBalanceCH(sumch, ChannelWH);
//    #pragma omp parallel for num_threads(userCores)
//    FOR_ROW_COL_MV_CHL
//    {
//        ChannelWaterVol->Drc = ChannelWH->Drc*ChannelWidth->Drc*ChannelDX->Drc;
//    }}

    // get the maximum for output
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL
    {
        maxChannelflow->Drc = std::max(maxChannelflow->Drc, ChannelQn->Drc);
        maxChannelWH->Drc = std::max(maxChannelWH->Drc, ChannelWH->Drc);
    }}

    if (SwitchErosion)
    {
        if (SwitchLinkedList) {

            if(SwitchUse2Phase) {
                ChannelQBLsn->setAllMV();
                FOR_ROW_COL_LDDCH5 {
                    routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn,
                                   ChannelAlpha, ChannelDX, ChannelWaterVol, ChannelBLSed);
                }}
              //  cover(*ChannelQBLsn, *LDD, 0);
            }
            ChannelQSSsn->setAllMV();
            //route water 1D and sediment
            FOR_ROW_COL_LDDCH5 {
               routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn,
                              ChannelAlpha, ChannelDX, ChannelWaterVol, ChannelSSSed);
            }}
           // cover(*ChannelQSSsn, *LDD, 0);
        } else {

            KinematicSubstance(crlinkedlddch_, nrValidCellsCH, LDDChannel, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn, ChannelAlpha, ChannelDX, ChannelSSSed);

            if(SwitchUse2Phase) {
                KinematicSubstance(crlinkedlddch_, nrValidCellsCH, LDDChannel, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn, ChannelAlpha, ChannelDX, ChannelBLSed);
            }
        }
        if (SwitchIncludeDiffusion) {
            RiverSedimentDiffusion(_dt, ChannelSSSed, ChannelSSConc);
            // note SSsed goes in and out, SSconc is recalculated inside
        }
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            RiverSedimentLayerDepth(r,c);
            RiverSedimentMaxC(r,c);
            ChannelQsn->Drc = ChannelQSSsn->Drc;
            ChannelSed->Drc = ChannelSSSed->Drc;
            if(SwitchUse2Phase) {
                ChannelQsn->Drc += ChannelQBLsn->Drc;
                ChannelSed->Drc += ChannelBLSed->Drc;
            }
        }}
    }
}


/* not used */
double TWorld::getMassCH(cTMap *M)
{
    double sum2 = 0;
    #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        sum2 += M->Drc*ChannelDX->Drc*ChannelWidth->Drc;
    }}
    return sum2;
}
/* not used */
void TWorld::correctMassBalanceCH(double sum1, cTMap *M)
{
    double sum2 = 0;
    double n = 0;

    #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        sum2 += M->Drc*DX->Drc*ChannelAdj->Drc;
        n += 1;
    }}
    // total and cells active for M
    double dhtot = fabs(sum2) > 0 ? (sum1 - sum2)/sum2 : 0;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        M->Drc = M->Drc*(1.0 + dhtot);            // <- distribution weighted to h
        M->Drc = std::max(M->Drc , 0.0);
    }}
}
