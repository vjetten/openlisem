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
 \file lisChannelflow.cpp
 \brief Channel hydrology and sediment detachment and movement processes.

functions: \n
- void TWorld::CalcVelDischChannel(void) calculate Velocity, alpha and Q in the channel \n
- void TWorld::ChannelFlow(void) calculate channelflow, ChannelDepth, do kinematic wave \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"

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
// V, alpha and Q in the channel, called after overland flow vol to channel
// called after flood and uses new channel flood water height
void TWorld::CalcVelDischChannelNT()
{
    if (!SwitchIncludeChannel)
        return;
    /*
    dw      FW      dw
   \  |            |  /
    \ |         wh | /
     \|____________|/
  */

    FOR_ROW_COL_MV_CH  {
        double Perim, Radius, Area;
        const double beta = 0.6;
        double beta1 = 1/beta;
        double wh = ChannelWH->Drc;
        double FW = ChannelFlowWidth->Drc;
        double grad = sqrt(ChannelGrad->Drc);

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
            ChannelAlpha->Drc = std::pow(ChannelN->Drc/grad * std::pow(Perim, 2.0/3.0),beta);
        else
            ChannelAlpha->Drc = 0;

        if (ChannelAlpha->Drc > 0)
            ChannelQ->Drc = std::pow(Area/ChannelAlpha->Drc, beta1);
        else
            ChannelQ->Drc = 0;

        ChannelV->Drc = std::pow(Radius, 2.0/3.0)*grad/ChannelN->Drc;
    }
}
//---------------------------------------------------------------------------
void TWorld::ChannelAddBaseandRainNT(void)
{
    if (!SwitchIncludeChannel)
        return;

    FOR_ROW_COL_MV_L {
        if(ChannelMaskExtended->data[r][c] == 1)
        {
            int rr = (int)ChannelSourceYExtended->Drc;
            int cr = (int)ChannelSourceXExtended->Drc;

            if (ChannelMaxQ->Drc <= 0) {
                ChannelWaterVol->Drcr += Rainc->Drc*ChannelWidthMax->Drcr*DX->Drcr;
            }

            // subtract infiltration
            if (SwitchChannelInfil) {
                double inf = ChannelDX->Drc * ChannelKsat->Drc*_dt/3600000.0 * (ChannelWidth->Drc + 2.0*ChannelWH->Drc/cos(atan(ChannelSide->Drc)));
                inf = std::min(ChannelWaterVol->Drc, inf);
               ChannelWaterVol->Drc -= inf;
               ChannelInfilVol->Drc += inf;
            }

            //add baseflow
            if(SwitchChannelBaseflow)
            {
                if(!addedbaseflow)
                {
                    ChannelWaterVol->Drc += BaseFlowInitialVolume->Drc;
                    BaseFlowTot += BaseFlowInitialVolume->Drc;

                }
                ChannelWaterVol->Drc += BaseFlowInflow->Drc * _dt;
                BaseFlowTot += BaseFlowInflow->Drc * _dt;
            }
        }

        ChannelWaterVol->Drc = std::max(0.0, ChannelWaterVol->Drc);
        fromChannelVoltoWH(r, c);
    }

    if(SwitchChannelBaseflow && !addedbaseflow)
        addedbaseflow = true;
}
//---------------------------------------------------------------------------
//! add runofftochannel and rainfall and calc channel WH from volume
void TWorld::ChannelWaterHeightFromVolumeNT(void)
{

    if(!SwitchIncludeChannel)
        return;

    //TODO: extended channel!
    FOR_ROW_COL_MV_CH {
        fromChannelVoltoWH(r, c);
    }
}
//---------------------------------------------------------------------------
//! add runofftochannel and rainfall and calc channel WH from volume
void TWorld::ChannelWaterHeightFromVolume(int thread)
{

    if(!SwitchIncludeChannel)
        return;

    FOR_ROW_COL_2DMT
    {
        if(!pcr::isMV(LDDChannel->Drc)) {
           fromChannelVoltoWH(r, c);
        }
    }}}}
}
//---------------------------------------------------------------------------
//! calc channelflow, ChannelDepth, kin wave
//! channel WH and V and Q are clculated before
void TWorld::ChannelFlow(void)
{
    if (!SwitchIncludeChannel)
        return;

    // initialize some channel stuff
    FOR_ROW_COL_MV_CH
    {
        ChannelQsn->Drc = 0;
        Channelq->Drc = 0;
//        if (SwitchChannelInfil)
//        {
//            Channelq->Drc =  -(ChannelKsat->Drc *  ChannelPerimeter->Drc/3600000.0);
//            //mm/h / 1000 = m/h / 3600 = m/s * m = m2/s
//        }
    }

    //VJ calc concentrations and ingoing Qs
    if (SwitchErosion)
    {
        if(!SwitchUseGrainSizeDistribution)
        {
            FOR_ROW_COL_MV_CH
            {
                double concbl = MaxConcentration(ChannelWaterVol->Drc, &ChannelBLSed->Drc, &ChannelDep->Drc);
                double concss = MaxConcentration(ChannelWaterVol->Drc, &ChannelSSSed->Drc, &ChannelDep->Drc);
                //temp conc because we move everything with channelQ

                //ChannelConc->Drc = (concbl + concss); // allowed because of CH vol
                //ChannelQs->Drc =  ChannelQ->Drc * ChannelConc->Drc;
                ChannelQBLs->Drc = ChannelQ->Drc * concbl;
                ChannelQSSs->Drc = ChannelQ->Drc * concss;
            }

        } else {
            double concbl = 0;
            double concss = 0;
            FOR_GRAIN_CLASSES
            {
                FOR_ROW_COL_MV_CH
                 {
                    RBLC_D.Drcd = MaxConcentration(ChannelWaterVol->Drc, &RBL_D.Drcd, &ChannelDep->Drc);
                    RSSC_D.Drcd = MaxConcentration(ChannelWaterVol->Drc, &RSS_D.Drcd, &ChannelDep->Drc);
                    concbl += RBLC_D.Drcd;
                    concss += RSSC_D.Drcd;

                    ChannelConc->Drc += RBLC_D.Drcd + RSSC_D.Drcd;

                    Tempa_D.Drcd = ChannelQ->Drc * RBLC_D.Drcd;
                    Tempc_D.Drcd = ChannelQ->Drc * RSSC_D.Drcd;
                 }
                fill(*Tempb_D.at(d), 0.0);
                fill(*Tempd_D.at(d), 0.0);
            }
            FOR_ROW_COL_MV_CH
            {
                //ChannelQs->Drc =  ChannelQ->Drc * ChannelConc->Drc;
                ChannelQBLs->Drc = ChannelQ->Drc * concbl;
                ChannelQSSs->Drc = ChannelQ->Drc * concss;
            }
        }

        ChannelQBLsn->setAllMV();
        ChannelQSSsn->setAllMV();
        if(SwitchUseGrainSizeDistribution)
        {
            FOR_GRAIN_CLASSES
            {
                Tempb_D.at(d)->setAllMV();
                Tempd_D.at(d)->setAllMV();
            }
        }
    }

    ChannelQn->setAllMV();
    fill(*QinKW, 0.0);

    // route water 1D and sediment
    FOR_ROW_COL_MV_CH
    {
        if (LDDChannel->Drc == 5)
            Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, Channelq, ChannelAlpha, ChannelDX, ChannelMaxQ);
    }
    cover(*ChannelQn, *LDD, 0);

    FOR_ROW_COL_MV_CH
    {
        ChannelWaterVol->Drc = ChannelWaterVol->Drc + QinKW->Drc*_dt - ChannelQn->Drc*_dt ;

        double ChannelArea = ChannelWaterVol->Drc/ChannelDX->Drc;
      //          double ChannelArea =ChannelAlpha->Drc*std::pow(ChannelQn->Drc, 0.6);
        ChannelV->Drc = (ChannelArea > 0 ? ChannelQn->Drc/ChannelArea : 0);

        fromChannelVoltoWH(r, c);

    }

    // get the maximum for output
    FOR_ROW_COL_MV_CH
    {
        maxChannelflow->Drc = std::max(maxChannelflow->Drc, ChannelQn->Drc);
        maxChannelWH->Drc = std::max(maxChannelWH->Drc, ChannelWH->Drc);
    }

    if (SwitchErosion)
    {
        FOR_ROW_COL_MV_CH
        {
            RiverSedimentLayerDepth(r,c);
            RiverSedimentMaxC(r,c);
        }
        // route water 1D and sediment
        FOR_ROW_COL_MV_CH
        {
            if (LDDChannel->Drc == 5)
            {
                //explicit routing of matter using Q and new Qn
                if(!SwitchUseGrainSizeDistribution)
                {
                    if (SwitchUse2Layer) {
                        routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn,
                                       ChannelAlpha, ChannelDX, ChannelWaterVol, ChannelBLSed);
                    }
                    routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn,
                                   ChannelAlpha, ChannelDX, ChannelWaterVol, ChannelSSSed);
                    //note: channelwatervol not really used

                } else {
                    FOR_GRAIN_CLASSES
                    {
                        routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, Tempa_D.at(d), Tempb_D.at(d),
                                       ChannelAlpha, ChannelDX, ChannelWaterVol, RBL_D.at(d));
                        routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, Tempc_D.at(d), Tempd_D.at(d),
                                       ChannelAlpha, ChannelDX, ChannelWaterVol, RSS_D.at(d));
                    }
                }
            }
        }

        cover(*ChannelQBLsn, *LDD, 0);
        cover(*ChannelQSSsn, *LDD, 0);
        FOR_ROW_COL_MV_CH
        {
           ChannelSSConc->Drc = MaxConcentration(ChannelWaterVol->Drc, &ChannelSSSed->Drc, &ChannelDep->Drc);
        }
        if (SwitchIncludeRiverDiffusion) {
            if(SwitchUseGrainSizeDistribution)
            {
                FOR_GRAIN_CLASSES
                {
                    cover(*(Tempb_D.at(d)), *LDD, 0);
                    cover(*(Tempd_D.at(d)), *LDD, 0);
                }

                FOR_GRAIN_CLASSES
                {
                    RiverSedimentDiffusion(_dt, RSS_D.at(d),RSSC_D.at(d));
                }
                FOR_ROW_COL_MV_CH
                {
                    ChannelSed->Drc = 0;
                    ChannelConc->Drc = 0;
                    ChannelBLConc->Drc = 0;
                    ChannelSSConc->Drc = 0;
                    ChannelBLSed->Drc = 0;
                    ChannelSSSed->Drc = 0;
                }
                FOR_GRAIN_CLASSES
                {
                    FOR_ROW_COL_MV_CH
                    {
                        ChannelBLSed->Drc += RBL_D.Drcd;
                        ChannelSSSed->Drc += RSS_D.Drcd;
                        ChannelSed->Drc += RSS_D.Drcd + RBL_D.Drcd;
                        RBLC_D.Drcd = (ChannelQn->Drc > MIN_FLUX ? Tempb_D.Drcd/ChannelQn->Drc : 0);
                        RSSC_D.Drcd = (ChannelQn->Drc > MIN_FLUX ? Tempd_D.Drcd/ChannelQn->Drc : 0);
                        ChannelConc->Drc += RBLC_D.Drcd + RSSC_D.Drcd;
                        ChannelBLConc->Drc += RBLC_D.Drcd;
                        ChannelSSConc->Drc += RSSC_D.Drcd;
                    }
                }
            }
            else
            {
                RiverSedimentDiffusion(_dt, ChannelSSSed, ChannelSSConc);
                // note SSsed goes in and out, SSconc is recalculated inside

            }
        }

        if(!SwitchUseGrainSizeDistribution)
        {
            FOR_ROW_COL_MV_CH
            {
                RiverSedimentLayerDepth(r,c);
                RiverSedimentMaxC(r,c);
          //       ChannelQsn->Drc = ChannelConc->Drc * ChannelQn->Drc;
                ChannelQsn->Drc = ChannelQBLsn->Drc + ChannelQSSsn->Drc;
                ChannelSed->Drc = ChannelBLSed->Drc + ChannelSSSed->Drc;

            }
        }
        else
        {
            FOR_ROW_COL_MV_CH
            {
                ChannelQBLsn->Drc = 0;
                ChannelQSSsn->Drc= 0;
                //ChannelQs->Drc = 0;
                ChannelQsn->Drc = 0;
                ChannelSed->Drc = 0;
                ChannelConc->Drc = 0;
                ChannelBLConc->Drc = 0;
                ChannelSSConc->Drc = 0;
                ChannelBLSed->Drc = 0;
                ChannelSSSed->Drc = 0;
            }

            FOR_GRAIN_CLASSES
            {
                FOR_ROW_COL_MV_CH
                 {
                    ChannelQBLsn->Drc += Tempb_D.Drcd ;
                    ChannelQSSsn->Drc += Tempd_D.Drcd;
                    ChannelQsn->Drc += Tempb_D.Drcd + Tempd_D.Drcd;
                    ChannelBLSed->Drc += RBL_D.Drcd;
                    ChannelSSSed->Drc += RSS_D.Drcd;
                    ChannelSed->Drc += RSS_D.Drcd + RBL_D.Drcd;
                    RBLC_D.Drcd = (ChannelQn->Drc > MIN_FLUX ? Tempb_D.Drcd/ChannelQn->Drc : 0);
                    RSSC_D.Drcd = (ChannelQn->Drc > MIN_FLUX ? Tempd_D.Drcd/ChannelQn->Drc : 0);
                    ChannelConc->Drc += RBLC_D.Drcd + RSSC_D.Drcd;
                    ChannelBLConc->Drc += RBLC_D.Drcd;
                    ChannelSSConc->Drc += RSSC_D.Drcd;
                 }
            }
        }
    }
}
//---------------------------------------------------------------------------
