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

/* not used */
double TWorld::channelVoltoWH(double vol, int r, int c)
{
    double cwh = 0;
    if (vol == 0) {
        return 0;
    } else {
        if (ChannelSide->Drc == 0) {
            return(vol/(ChannelWidth->Drc*ChannelDX->Drc));
        } else {
            double maxvol = ChannelDX->Drc * (ChannelDepth->Drc*(ChannelFlowWidth->Drc+ChannelWidth->Drc)/2.0);

            if (vol < maxvol) {
                // water below surface, WH from abc rule
                double aa = ChannelSide->Drc;  //=tan(a)
                double bb = ChannelWidth->Drc; //=w
                double cc = - vol/ChannelDX->Drc; //=area

                cwh = std::max(0.0,(-bb+sqrt(bb*bb - 4.0*cc*aa))/(2.0*aa));
                if (cwh < 0) {
                    ErrorString = QString("Channel water height is negative at row %1, col %2").arg(r).arg(c);
                    throw 1;
                }
            } else {
                // water above surface, WH is depth + part sticking out
                cwh = ChannelDepth->Drc + (ChannelWaterVol->Drc - maxvol)/(ChannelDX->Drc*ChannelFlowWidth->Drc);
            }
        }
    }
    return (cwh);
}

void TWorld::fromChannelVoltoWH(int r, int c)
{

    //    non-rectangular, ABC fornula
    //            dw      w       dw
    //         __|   |            |   |__ surface, above water becomes rectangular
    //            \  |            |  /
    //             \ |          h |a/  <= tan(a) is channelside = tan angle of side wall
    //              \|____________|/
    //            area = h*w + h*dw
    //            tan(a) = dw/h, dw = h*tan(a) = h*side
    //            area =volume/DX
    //            tan(a)h^2 + w*h - area = 0
    //            aa (h2)   +   bb(h) +  cc = 0

    if (ChannelWaterVol->Drc == 0) {
        ChannelWH->Drc = 0;
        ChannelFlowWidth->Drc = 0;
    } else {
        if (ChannelSide->Drc == 0) {
            ChannelFlowWidth->Drc = ChannelWidth->Drc;
            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
        } else {
            double maxvol = ChannelDX->Drc * (ChannelDepth->Drc*(ChannelFlowWidth->Drc+ChannelWidth->Drc)/2.0);

            if (ChannelWaterVol->Drc < maxvol) {
                // water below surface, WH from abc rule
                double aa = ChannelSide->Drc;  //=tan(a)
                double bb = ChannelWidth->Drc; //=w
                double cc = - ChannelWaterVol->Drc/ChannelDX->Drc; //=area

                ChannelWH->Drc = std::max(0.0,(-bb+sqrt(bb*bb - 4.0*cc*aa))/(2.0*aa));
                if (ChannelWH->Drc < 0) {
                    ErrorString = QString("Channel water height is negative at row %1, col %2").arg(r).arg(c);
                    throw 1;
                }
                ChannelFlowWidth->Drc = std::min(ChannelWidthMax->Drc, ChannelWidth->Drc + 2.0*ChannelSide->Drc * ChannelWH->Drc);
            } else {
                // water above surface, WH is depth + part sticking out
                ChannelFlowWidth->Drc = ChannelWidthMax->Drc;
                ChannelWH->Drc = ChannelDepth->Drc + (ChannelWaterVol->Drc - maxvol)/(ChannelDX->Drc*ChannelFlowWidth->Drc);
            }
        }
    }
}
//---------------------------------------------------------------------------
void TWorld::fromChannelWHtoVol(int r, int c)
{
    if (ChannelSide->Drc == 0) {
        ChannelWaterVol->Drc = ChannelWidth->Drc * ChannelWH->Drc * ChannelDX->Drc;
        ChannelFlowWidth->Drc = ChannelWidth->Drc;
        return;
    }
    if (ChannelWH->Drc > ChannelDepth->Drc) {
        ChannelFlowWidth->Drc = ChannelWidthMax->Drc;
        ChannelWaterVol->Drc = (ChannelWidth->Drc + ChannelFlowWidth->Drc)*0.5 * ChannelDepth->Drc * ChannelDX->Drc +
                (ChannelWH->Drc - ChannelDepth->Drc)*ChannelWidthMax->Drc;
    } else {
        ChannelFlowWidth->Drc = ChannelWidth->Drc + ChannelWH->Drc*ChannelSide->Drc*2.0;
        ChannelWaterVol->Drc = (ChannelWidth->Drc + ChannelFlowWidth->Drc)*0.5 * ChannelWH->Drc * ChannelDX->Drc;
    }

}
//---------------------------------------------------------------------------
/* NOT USED */

// V, alpha and Q in the channel, called after overland flow vol to channel
// called after flood and uses new channel flood water height
void TWorld::CalcVelDischChannel()
{
    if (!SwitchIncludeChannel)
        return;
    if(!SwitchChannelKinWave)
        return;

    /*
    dw      FW      dw
   \  |            |  /
    \ |         wh | /
     \|____________|/
  */
#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL  {
        double Perim, Radius, Area;
        double grad = sqrt(ChannelGrad->Drc);

        ChannelFlowWidth->Drc = ChannelWidth->Drc;
        ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);

        double wh = ChannelWH->Drc;
        double FW = ChannelFlowWidth->Drc;

        if (ChannelSide->Drc > 0)
        {
            double dw = ChannelSide->Drc * wh;
            Perim = FW + 2.0*wh/cos(atan(ChannelSide->Drc));
            Area = FW*wh + wh*dw;
        }
        else
        {
            Perim = FW + 2.0*wh;
            Area = FW*wh;
        }

        Radius = (Perim > 0 ? Area/Perim : 0);

        if (grad > MIN_SLOPE)
            ChannelAlpha->Drc = std::pow(ChannelN->Drc/grad * std::pow(Perim, 2.0/3.0),0.6);
        else
            ChannelAlpha->Drc = 0;

        if (ChannelAlpha->Drc > 0) {
            ChannelQ->Drc = std::pow(Area/ChannelAlpha->Drc, 1.0/0.6);
            if (SwitchCulverts) {
                if (ChannelMaxQ->Drc > 0 && ChannelQ->Drc > ChannelMaxQ->Drc){
                    ChannelAlpha->Drc = Area/std::pow(ChannelMaxQ->Drc, 0.6);
                    ChannelQ->Drc = ChannelMaxQ->Drc;
                }
            }
        }
        else
            ChannelQ->Drc = 0;

        ChannelV->Drc = std::pow(Radius, 2.0/3.0)*grad/ChannelN->Drc;
    }}
}
//---------------------------------------------------------------------------
void TWorld::ChannelAddBaseandRain(void)
{
    if (!SwitchIncludeChannel)
        return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {

        if (ChannelMaxQ->Drc <= 0)
            ChannelWaterVol->Drc += Rainc->Drc*ChannelWidthMax->Drc*DX->Drc;

        // subtract infiltration
        if (SwitchChannelInfil) {
            double inf = ChannelDX->Drc * ChannelKsat->Drc*_dt/3600000.0 * (ChannelWidth->Drc + 2.0*ChannelWH->Drc/cos(atan(ChannelSide->Drc)));
            inf = std::min(ChannelWaterVol->Drc, inf);
            ChannelWaterVol->Drc -= inf;
            ChannelInfilVol->Drc += inf;
        }

        if(SwitchChannelBaseflow)
        {
            if(!addedbaseflow)
            {
                ChannelWaterVol->Drc += BaseFlowInitialVolume->Drc;
                //BaseFlowTot += BaseFlowInitialVolume->Drc;
            }
            ChannelWaterVol->Drc += BaseFlowInflow->Drc * _dt;
            //BaseFlowTot += BaseFlowInflow->Drc * _dt;
        }

        ChannelFlowWidth->Drc = ChannelWidth->Drc;
        ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);

    }}
    if (!addedbaseflow)
        addedbaseflow = true;
}
//---------------------------------------------------------------------------
    /* NOT USED */
//! add runofftochannel and rainfall and calc channel WH from volume
void TWorld::ChannelWaterHeightFromVolume()
{

    if(!SwitchIncludeChannel)
        return;

#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(!pcr::isMV(LDDChannel->Drc)) {
            ChannelFlowWidth->Drc = ChannelWidth->Drc;
            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
//            fromChannelVoltoWH(r, c);
        }
    }}
}
//---------------------------------------------------------------------------
//! calc channelflow, ChannelDepth, kin wave
//! channel WH and V and Q are clculated before
void TWorld::ChannelFlow(void)
{
    if (!SwitchIncludeChannel)
        return;


    // velocity, alpha, Q
#pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_CHL {

        // calc velocity and Q
        ChannelFlowWidth->Drc = ChannelWidth->Drc;
        ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);

        double Perim, Radius, Area;
        double sqrtgrad = sqrt(ChannelGrad->Drc);
        double N = ChannelN->Drc;
        double MaxQ = ChannelMaxQ->Drc;
        ChannelV->Drc = 0;
        ChannelQ->Drc = 0;
        double wh = ChannelWH->Drc;
        double FW = ChannelWidth->Drc;
        if (wh > 1e-10) {
            Perim = FW + 2.0*wh;
            Area = FW*wh;
            Radius = (Perim > 0 ? Area/Perim : 0);

            if (sqrtgrad > MIN_SLOPE) {
                ChannelAlpha->Drc = std::pow(N/sqrtgrad * std::pow(Perim, 2.0/3.0), 0.6);
                ChannelQ->Drc = std::pow(Area/ChannelAlpha->Drc, 1.0/0.6);
                if (SwitchCulverts) {
                    if (MaxQ > 0 && ChannelQ->Drc > MaxQ){
                        ChannelAlpha->Drc = Area/std::pow(MaxQ, 0.6);
                        Perim = std::pow(std::pow(ChannelAlpha->Drc, 1.0/0.6)*sqrtgrad/N,3.0/2.0);
                        wh = Perim/FW;
                        Radius = (wh*FW)/Perim;
                        ChannelQ->Drc = MaxQ;
                        ChannelWH->Drc = wh;
                    }
                }
            }
            ChannelV->Drc = std::pow(Radius, 2.0/3.0)*sqrtgrad/N;
        }

        ChannelQsn->Drc = 0;
        Channelq->Drc = 0;
        QinKW->Drc = 0;

        if (SwitchErosion) {
            double concss = MaxConcentration(ChannelWaterVol->Drc, &ChannelSSSed->Drc, &ChannelDep->Drc);
            ChannelQSSs->Drc = ChannelQ->Drc * concss;

            if(SwitchUse2Phase) {
                double concbl = MaxConcentration(ChannelWaterVol->Drc, &ChannelBLSed->Drc, &ChannelDep->Drc);
                ChannelQBLs->Drc = ChannelQ->Drc * concbl;
            }
        }

    }}

    if (SwitchChannelKinWave) {

        bool do_avg = false;
        int loop = 1;
        double dtkin = 60.0;
        if (_dt > dtkin) {
            loop = int(_dt/dtkin);
            _dt = dtkin;
        }

        if (do_avg) {
            fill(*tma,0);
            fill(*tmb,0);
        }

        for (int i = 0; i < loop; i++) {


            if (SwitchLinkedList) {

                ChannelQn->setAllMV();
                // route water 1D and sediment
                FOR_ROW_COL_LDDCH5 {
                    Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, Channelq, ChannelAlpha, ChannelDX, ChannelMaxQ);
                }}
                cover(*ChannelQn, *LDD, 0);

            } else {
                KinematicExplicit(crlinkedlddch_, nrValidCellsCH, LDDChannel, ChannelQ, ChannelQn, Channelq, ChannelAlpha, ChannelDX, ChannelMaxQ);
            }

            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_CHL {
                if (do_avg) {
                   tma->Drc += ChannelQn->Drc;
                   tmb->Drc += QinKW->Drc;
                }
                ChannelQ->Drc = ChannelQn->Drc;
            }}
        }
        _dt = _dt_user;

        if (do_avg) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_CHL {
                ChannelQn->Drc = tma->Drc/(double) loop;
                QinKW->Drc = tmb->Drc/(double) loop;
            }}
        }

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelWaterVol->Drc = ChannelWaterVol->Drc + QinKW->Drc*_dt - ChannelQn->Drc*_dt ;
            //mass balance

            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
            // new channel WH

            double ChannelArea = ChannelWaterVol->Drc/ChannelDX->Drc;
            double R = ChannelArea/(ChannelWidth->Drc + 2*ChannelWH->Drc);
           // ChannelV->Drc = (ChannelArea > 0 ? ChannelQn->Drc/ChannelArea : 0);
            ChannelV->Drc = pow(R, 2.0/3.0) * sqrt(ChannelGrad->Drc)/ChannelN->Drc;
            // new velocity

            //fromChannelVoltoWH(r, c);
            ChannelFlowWidth->Drc = ChannelWidth->Drc;

        }}

    } else {

        ChannelSWOFopen();

    }
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

