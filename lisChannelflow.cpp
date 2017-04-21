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

//---------------------------------------------------------------------------
// V, alpha and Q in the channel, called after overland flow vol to channel
// called after flood and uses new channel flood water height
void TWorld::CalcVelDischChannel(void)
{

    if (!SwitchIncludeChannel)
        return;
    /*
    dw      FW      dw
   \  |            |  /
    \ |         wh | /
     \|____________|/
  */
    FOR_ROW_COL_MV_CH
    {
        double Perim, Radius, Area;
        const double beta = 0.6;
        const double _23 = 2.0/3.0;
        double beta1 = 1/beta;
        double wh = ChannelWH->Drc;
        double FW = ChannelWidth->Drc;
        double grad = sqrt(ChannelGrad->Drc);

        if (ChannelSide->Drc > 0)  //BUG, was (ChannelSide > 0) !!!
        {
            double dw = ChannelSide->Drc * wh;
            Perim = FW + 2.0*sqrt(wh*wh + dw*dw);//*cos(atan(ChannelSide->Drc));
            Area = FW*wh + wh*dw;
        }
        else
        {
            Perim = FW + 2.0*wh;
            Area = FW*wh;
        }

        Radius = (Perim > 0 ? Area/Perim : 0);

        if (grad > MIN_SLOPE)
        ChannelAlpha->Drc = std::pow(ChannelN->Drc/grad * std::pow(Perim, _23),beta);
        else
            ChannelAlpha->Drc = 0;

        if (ChannelAlpha->Drc > 0)
            ChannelQ->Drc = std::pow(Area/ChannelAlpha->Drc, beta1);
        else
            ChannelQ->Drc = 0;

        //ChannelQ->Drc = pow(Area/Perim, 2.0/3.0)*sqrt(grad)/ChannelN->Drc * Area;

        if (SwitchChannelFlood)
        {
            if (ChannelMaxQ->Drc > 0)
            {
                ChannelQ->Drc = std::min(ChannelQ->Drc, ChannelMaxQ->Drc);
                Area = ChannelAlpha->Drc*std::pow(ChannelQ->Drc, beta);
                Perim = FW + Area/FW*2;
                ChannelWH->Drc = Area/FW;
            }

        }

        ChannelV->Drc = std::pow(Radius, _23)*grad/ChannelN->Drc;

        ChannelWaterVol->Drc = Area * ChannelDX->Drc;

        ChannelPerimeter->Drc = Perim;
        //VJ 110109 needed for channel infil
    }
}
//---------------------------------------------------------------------------
//! add runofftochannel and rainfall and calc channel WH from volume
void TWorld::ChannelWaterHeight(void)
{

    if (!SwitchIncludeChannel)
    {
        return;
    }


    for (int  r = 0; r < _nrRows; r++)
    {
        for (int  c = 0; c < _nrCols; c++)
        {
            if(ChannelMaskExtended->data[r][c] == 1)
            {
                int rr = (int)ChannelSourceYExtended->Drc;
                int cr = (int)ChannelSourceXExtended->Drc;

                ChannelWaterVol->Drcr += Rainc->Drc*(_dx - ChannelAdj->Drc)*DX->Drc;
                // add rainfall in m3, no interception, rainfall so do not use ChannelDX

            }
        }
    }



    // calculate new channel WH , WidthUp and Volume
    FOR_ROW_COL_MV_CH
    {
        ChannelWH->Drc = 0;

        ChannelWaterVol->Drc += RunoffVolinToChannel->Drc;
        // water from overland flow in channel cells

        //add baseflow
        if(SwitchChannelBaseflow)
        {
            if(!addedbaseflow)
            {
                ChannelWaterVol->Drc += BaseFlowInitialVolume->Drc;
                BaseFlow += BaseFlowInitialVolume->Drc;

            }
            ChannelWaterVol->Drc += BaseFlowInflow->Drc * _dt;
            BaseFlow += BaseFlowInflow->Drc * _dt;

        }

        if (SwitchBuffers && ChannelBufferVol->Drc > 0)
        {
            ChannelBufferVol->Drc -= ChannelWaterVol->Drc;
            ChannelWaterVol->Drc = 0;
            // add inflow from slopes and rainfall to buffer
        }
    }

    if(!addedbaseflow)
    {
        addedbaseflow = true;
    }

    ChannelWaterHeightFromVolume();


}
//---------------------------------------------------------------------------
//! add runofftochannel and rainfall and calc channel WH from volume
void TWorld::ChannelWaterHeightFromVolume(void)
{

    if(!SwitchIncludeChannel)
    {
        return;
    }

    FOR_ROW_COL_MV_CH
    {
        // calculate ChannelWH
        if (ChannelSide->Drc == 0) // rectangular channel
        {
            ChannelFlowWidth->Drc = ChannelWidth->Drc;
            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
        }
        else
        {
            /*
             non-rectangular, ABC fornula
           dw      w       dw
   \  |            |  /
    \ |          h |a/  <= tan(a) is channelside = tan angle of side wall
     \|____________|/
   area = h*w + h*dw
   tan(a) = dw/h
   dw = h*tan(a)
   area = wh+dw*h = wh+tan(a)h*h
   tan(a)h^2 + wh - area = 0
   aa (h2)   +   bb(h) +  cc = 0
*/
            double aa = ChannelSide->Drc;  //=tan(a)
            double bb = ChannelWidth->Drc; //=w
            double cc = -ChannelWaterVol->Drc/ChannelDX->Drc; //=area

            ChannelWH->Drc = (-bb + sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa);

            if (ChannelWH->Drc < 0)
            {
                ErrorString = QString("channel water height is negative at row %1, col %2").arg(r).arg(c);
                throw 1;
            }

            ChannelFlowWidth->Drc = std::min(ChannelWidthMax->Drc , ChannelWidth->Drc + 2.0*ChannelSide->Drc*ChannelWH->Drc);


        }
    }
}

//---------------------------------------------------------------------------
double TWorld::ChannelIterateWH(double _h, int r, int c)
{
    double y1 = _h;
    double y = _h;
    //double dy = y/1000;
    double GN = sqrt(ChannelGrad->Drc)/ChannelN->Drc;
    double w = ChannelFlowWidth->Drc;
    double _53 =5.0/3.0;
    double _23 = 2.0/3.0;
    double Q = ChannelQn->Drc;
    //double dQ; //
    int count = 0;

//    double Q = GN*qPow(w*_h/(2*_h+w),_23)*w*_h;
//    double Alpha = qPow(1/GN * qPow((2*_h+w), _23),0.6);
//    double Q = qPow(_h*w/Alpha, 1/0.6);

/*
    if (Q == 0)
        return (0);
    if (Q > ChannelQn->Drc)
    {
        do{
            y -= dy;
            Q = GN*pow(w*y/(2*y+w),_23)*w*y;
            count++;
        }while(Q > ChannelQn->Drc);
    }
    else
    {
        do{
            y += dy;
            Q = GN*pow(w*y/(2*y+w),_23)*w*y;
            count++;
        }while(Q < ChannelQn->Drc);
    }
*/


    if (Q == 0)
       return 0;
    //int count = 0;
    do{
        double P = w+2*y;
        double A = y*w;

        //y1 = y - (1-Q/(sqrtgrad/n * ((qPow(A,5.0/3.0))/(qPow(P, 2.0/3.0))))) / ((5*width+6*y)/(3*y*P));


        //s*w^(5/3)*x^(2/3)*(6*x+5*w)/(3*(2*x+w)^(5/3))

        Q = GN*A*pow(A/P,_23);

        double dQ = GN*(pow(w, _53)*pow(y,_23)*(6*y+5*w))/(3*pow(P,_53));

        if(dQ <= 0)
        {
            return _h;
        }
      //  if (fabs(dQ) > 1e-10)
            y1 = y - (Q/dQ); //_dt*_dx*w*
//        else
  //          count = 100;

        y = y1;//std::max(y1, 0.0);
if(y < MIN_HEIGHT)
{
    y = 0;
    break;
}
        count++;

    }while(fabs(Q-ChannelQn->Drc) > 1e-6 && count < 100);

//        double a = ChannelAlpha->Drc;
//        ChannelAlpha->Drc = qPow(n/sqrtgrad * powl((2*y1+width), 2/3),0.6);

        qDebug() << count << y << _h << Q << ChannelQn->Drc;

    return std::max(y, 0.0);
}
//---------------------------------------------------------------------------
//! calc channelflow, ChannelDepth, kin wave
//! channel WH and V and Q are clculated before
void TWorld::ChannelFlow(void)
{
    if (!SwitchIncludeChannel)
        return;

    if (SwitchErosion)
    {
        FOR_ROW_COL_MV_CH
        {
            ChannelFlowDetachment(r,c);
        }
    }

    // initialize some channel stuff
    FOR_ROW_COL_MV_CH
    {
        ChannelQsn->Drc = 0;
        Channelq->Drc = 0;
        if (SwitchChannelInfil)
        {
            // NOTE: for buffers channelksat = 0
            Channelq->Drc =  -(ChannelKsat->Drc *  ChannelPerimeter->Drc/3600000.0);
            //mm/h / 1000 = m/h / 3600 = m/s * m = m2/s
        }
    }

    ChannelQn->setAllMV();
    fill(*QinKW, 0.0);
    // flag all new flux as missing value, needed in kin wave and replaced by new flux

    //VJ calc concentrations and ingoing Qs
    if (SwitchErosion)
    {
        if(!SwitchUseGrainSizeDistribution)
        {
            FOR_ROW_COL_MV_CH
            {
                double concbl = MaxConcentration(ChannelWaterVol->Drc, ChannelBLSed->Drc);
                double concss = MaxConcentration(ChannelWaterVol->Drc, ChannelSSSed->Drc);
                ChannelConc->Drc = (concbl + concss);
                ChannelQs->Drc =  ChannelQ->Drc * ChannelConc->Drc;
                ChannelQBLs->Drc = ChannelQ->Drc * concbl;
                ChannelQSSs->Drc = ChannelQ->Drc * concss;
            }

        }else
        {
            double concbl = 0;
            double concss = 0;
            FOR_GRAIN_CLASSES
            {
                FOR_ROW_COL_MV_CH
                 {
                    concbl += RBLC_D.Drcd =MaxConcentration(ChannelWaterVol->Drc, RBL_D.Drcd);
                    concss += RSSC_D.Drcd =MaxConcentration(ChannelWaterVol->Drc, RSS_D.Drcd);

                    ChannelConc->Drc += RBLC_D.Drcd + RSSC_D.Drcd;

                    Tempa_D.Drcd = ChannelQ->Drc * RBLC_D.Drcd;
                    Tempc_D.Drcd = ChannelQ->Drc * RSSC_D.Drcd;
                 }
                fill(*Tempb_D.at(d), 0.0);
                fill(*Tempd_D.at(d), 0.0);
            }
            FOR_ROW_COL_MV_CH
            {
                ChannelQs->Drc =  ChannelQ->Drc * ChannelConc->Drc;
                ChannelQBLs->Drc = ChannelQ->Drc * concbl;
                ChannelQSSs->Drc = ChannelQ->Drc * concss;
            }
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

    ChannelQBLsn->setAllMV();
    ChannelQSSsn->setAllMV();

    if(SwitchUseGrainSizeDistribution)
    FOR_GRAIN_CLASSES
    {
        Tempb_D.at(d)->setAllMV();
        Tempd_D.at(d)->setAllMV();
    }

    // route water 1D and sediment
    FOR_ROW_COL_MV_CH
    {
        if (LDDChannel->Drc == 5)
        {
            Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, Channelq, ChannelAlpha, ChannelDX,
                      ChannelWaterVol, ChannelBufferVol);
            // kin wave on water

            if (SwitchErosion)
            {
                if(!SwitchUseGrainSizeDistribution)
                {
                    routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn,
                                   ChannelAlpha, ChannelDX, ChannelWaterVol, ChannelBLSed, ChannelBufferVol, ChannelBufferSed);
                    routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn,
                                   ChannelAlpha, ChannelDX, ChannelWaterVol, ChannelSSSed, ChannelBufferVol, ChannelBufferSed);
                    //explicit routing of matter using Q and new Qn
                }else
                {
                    FOR_GRAIN_CLASSES
                    {
                        routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, Tempa_D.at(d), Tempb_D.at(d),
                                       ChannelAlpha, ChannelDX, ChannelWaterVol, RBL_D.at(d), ChannelBufferVol, ChannelBufferSed);
                        routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, Tempc_D.at(d), Tempd_D.at(d),
                                       ChannelAlpha, ChannelDX, ChannelWaterVol, RSS_D.at(d), ChannelBufferVol, ChannelBufferSed);
                    }
                }
            }
        }
    }

    cover(*ChannelQn, *LDD, 0);
    // avoid missing values around channel for adding to Qn for output
    if (SwitchChannelFlood)
    {
        FOR_ROW_COL_MV_CH
        {
        if (ChannelMaxQ->Drc > 0)
            ChannelQn->Drc = std::min(ChannelQn->Drc, ChannelMaxQ->Drc);
        }
    }
    // limit channel Q when culverts > 0


    if (SwitchErosion)
    {
        cover(*ChannelQBLsn, *LDD, 0);
        cover(*ChannelQSSsn, *LDD, 0);

        if(SwitchUseGrainSizeDistribution)
        {
            FOR_GRAIN_CLASSES
            {
                cover(*(Tempb_D.at(d)), *LDD, 0);
                cover(*(Tempd_D.at(d)), *LDD, 0);
            }

            FOR_GRAIN_CLASSES
            {
                RiverSedimentDiffusion(_dt, RBL_D.at(d),RBLC_D.at(d), RSS_D.at(d),RSSC_D.at(d));
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
                    RBLC_D.Drcd = (ChannelQn->Drc > 1e-6 ? Tempb_D.Drcd/ChannelQn->Drc : 0);
                    RSSC_D.Drcd = (ChannelQn->Drc > 1e-6 ? Tempd_D.Drcd/ChannelQn->Drc : 0);
                    ChannelConc->Drc += RBLC_D.Drcd + RSSC_D.Drcd;
                    ChannelBLConc->Drc += RBLC_D.Drcd;
                    ChannelSSConc->Drc += RSSC_D.Drcd;
                 }
            }
        }
        else
        {
            RiverSedimentDiffusion(_dt, ChannelBLSed,ChannelBLConc, ChannelSSSed,ChannelSSConc);
        }

        if(!SwitchUseGrainSizeDistribution)
        {
            FOR_ROW_COL_MV_CH
            {
                ChannelQsn->Drc = ChannelQBLsn->Drc + ChannelQSSsn->Drc;
                ChannelSed->Drc = ChannelBLSed->Drc + ChannelSSSed->Drc;
                //ChannelConc->Drc = (ChannelQn->Drc > 1e-6 ? ChannelQsn->Drc/ChannelQn->Drc : 0);
                // everywhere vol is used
                ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc);
                //VJ standing water has also a concentration

            }
        }
        else
        {
            FOR_ROW_COL_MV_CH
            {
                ChannelQBLsn->Drc = 0;
                ChannelQSSsn->Drc= 0;
                ChannelQs->Drc = 0;
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
                    RBLC_D.Drcd = (ChannelQn->Drc > 1e-6 ? Tempb_D.Drcd/ChannelQn->Drc : 0);
                    RSSC_D.Drcd = (ChannelQn->Drc > 1e-6 ? Tempd_D.Drcd/ChannelQn->Drc : 0);
                    ChannelConc->Drc += RBLC_D.Drcd + RSSC_D.Drcd;
                    ChannelBLConc->Drc += RBLC_D.Drcd;
                    ChannelSSConc->Drc += RSSC_D.Drcd;
                 }
            }
        }
    }

    /*// recalc WH from flux
    FOR_ROW_COL_MV_CH
    {
        //ChannelWH->Drc = ChannelIterateWH(r, c);
       // double ChannelArea = ChannelWH->Drc * (ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2.0;
       double ChannelArea = ChannelAlpha->Drc*std::pow(ChannelQn->Drc, 0.6);
       // calc area from Q and alpha, but alpha is from before the kin wave, so this should be iterated!
//double hest = (QinKW->Drc*_dt + ChannelWaterVol->Drc - ChannelQn->Drc*_dt)/(ChannelDX->Drc*ChannelWidthUpDX->Drc);
       //double ChannelArea = ChannelIterateWH(hest, r, c)*ChannelDX->Drc;//

       //qDebug() << ChannelArea;
        // explicit...!!! QniKW flowng in, vol is from before changes, Qn outflow after iteration

        // water height is not used except for output! i.e. watervolume is cycled
        if(ChannelSide->Drc > 0)
        {
            double aa = ChannelSide->Drc;  //=tan(a)
            double bb = ChannelWidth->Drc; //=w
            double cc = -ChannelArea; //=area

            ChannelWH->Drc = (-bb + sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa);
            //VJ gaat niet goed als side zodanig is dat width groter wordt dan gridcell bij nieuwe chanarea

        }
        else
        {
            ChannelWH->Drc = ChannelArea/((ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2.0);
        }

        ChannelV->Drc = (ChannelArea > 0 ? ChannelQn->Drc/ChannelArea : 0);
        // recalc for output screen and maps, new V in channel

    //    tma->Drc = ChannelWaterVol->Drc;
        ChannelWaterVol->Drc = ChannelArea * ChannelDX->Drc;
        // needed for concentration, new volume in channel

        //double error = (ChannelQ->Drc*_dt + tma->Drc - ChannelWaterVol->Drc - ChannelQn->Drc * _dt);
    }*/

    double outflux = 0;
    FOR_ROW_COL_MV_CH
    {
        ChannelWaterVol->Drc += -ChannelQn->Drc*_dt + QinKW->Drc * _dt;

        double ChannelArea = ChannelAlpha->Drc*std::pow(ChannelQn->Drc, 0.6);
        ChannelV->Drc = (ChannelArea > 0 ? ChannelQn->Drc/ChannelArea : 0);

        if(LDDChannel->Drc == 5)
        {
            outflux += ChannelQn->Drc*_dt;
        }
    }
    ChannelWaterHeightFromVolume();


    cover(*ChannelWH, *LDD, 0);
    //VJ ?? necessary for display ??

    FOR_ROW_COL_MV_CH
    {
        if (SwitchErosion)
        {
            ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc);
            // correct for very high concentrations, max 850 g/l
            // NOTE: nothing is done with this concentration, only for display,
            // so the display shows the conc after the kin wave
         //   ChannelConc->Drc = (ChannelQn->Drc > 1e-6 ? ChannelQs->Drc/ChannelQn->Drc : 0);
        }
    }

    FOR_ROW_COL_MV_CH
    {
        maxChannelflow->Drc = std::max(maxChannelflow->Drc, ChannelQn->Drc);
        maxChannelWH->Drc = std::max(maxChannelWH->Drc, ChannelWH->Drc);
    }
}
//---------------------------------------------------------------------------
