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

#include "model.h"

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
        double dw = /*0.5* */(ChannelWidthUpDX->Drc - FW); // extra width when non-rectamgular

        if (dw > 0)
        {
            //			Perim = FW + 2*sqrt(wh*wh + dw*dw);
            Perim = FW + 2*wh/cos(atan(ChannelSide->Drc));
            Area = FW*wh + wh*dw;
        }
        else
        {
            Perim = FW + 2*wh;
            Area = FW*wh;
        }

        Radius = (Perim > 0 ? Area/Perim : 0);

        ChannelAlpha->Drc = qPow(ChannelN->Drc/grad * powl(Perim, _23),beta);

        if (ChannelAlpha->Drc > 0)
            ChannelQ->Drc = qPow(Area/ChannelAlpha->Drc, beta1);
        else
            ChannelQ->Drc = 0;

        ChannelV->Drc = pow(Radius, _23)*grad/ChannelN->Drc;

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
        return;

    // calculate new channel WH , WidthUp and Volume
    FOR_ROW_COL_MV_CH
    {
        ChannelWH->Drc = 0;

        ChannelWaterVol->Drc += RunoffVolinToChannel->Drc;
        // water from overland flow in channel cells

        ChannelWaterVol->Drc += Rainc->Drc*ChannelWidthUpDX->Drc*ChannelDX->Drc;
        // add rainfall in m3, no interception, rainfall so do not use ChannelDX

        if (SwitchBuffers && ChannelBufferVol->Drc > 0)
        {
            ChannelBufferVol->Drc -= ChannelWaterVol->Drc;
            ChannelWaterVol->Drc = 0;
            // add inflow from slopes and rainfall to buffer
        }

        // calculate ChannelWH
        if (ChannelSide->Drc == 0) // rectangular channel
        {
            ChannelWidthUpDX->Drc = ChannelWidth->Drc;
            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
        }
        else
        {
            /*
 * non-rectangular
   ABC fornula
    dw      w       dw
   \  |            |  /
    \ |          h |a/  <= tan(a) is channelside = tan angle of side wall
     \|____________|/
   area = h*w + h*dw
   dw = h*tan(a)
   area = wh+tan(a)h^2
   tan(a)h^2 + wh - area = 0
   aa (h2)   +   bb(h) +  cc = 0
*/
            double aa = ChannelSide->Drc;  //=tan(a)
            double bb = ChannelWidth->Drc; //=w
            double cc = -ChannelWaterVol->Drc/ChannelDX->Drc; //=area

            ChannelWH->Drc = (-bb + sqrt(bb*bb - 4*aa*cc))/(2*aa);
            if (ChannelWH->Drc < 0)
            {
                ErrorString = QString("channel water height is negative at row %1, col %2").arg(r).arg(c);
                throw 1;
            }
            ChannelWidthUpDX->Drc = ChannelWidth->Drc + 2*ChannelSide->Drc*ChannelWH->Drc;

            if (ChannelWidthUpDX->Drc > _dx)
            {
                ErrorString = QString("channel width > dx at row %1, col %2").arg(r).arg(c);
                throw 1;
            }

            if (ChannelWidthUpDX->Drc < 0)
            {
                ErrorString = QString("channel width < 0 at row %1, col %2").arg(r).arg(c);
                throw 1;
            }

            ChannelWidthUpDX->Drc = min(0.9*_dx, ChannelWidthUpDX->Drc);
            // new channel width with new WH, goniometric, side is top angle tan, 1 is 45 degr
            // cannot be more than 0.9*_dx

            ChannelAdj->Drc = max(0, _dx - ChannelWidthUpDX->Drc);
            // experimental if channelwidth > dx
        }
    }
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
        ChannelFlowDetachment();
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

    ChannelQn->setMV();
    //ChannelQsn->fill(0);
    QinKW->fill(0);
    // flag all new flux as missing value, needed in kin wave and replaced by new flux

    FOR_ROW_COL_MV_CH
    {
        if (LDDChannel->Drc == 5)
        {
            Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQs, ChannelQsn, Channelq, ChannelAlpha, ChannelDX,
                      ChannelWaterVol, ChannelSed, ChannelBufferVol, ChannelBufferSed);

            //routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQS, ChannelQSn, ChannelAlpha, ChannelDX, ChannelWaterVol, ChannelSubs);
            /*
                   routing of substances add here!
                   do after kin wave so that the new flux ChannelQn out of a cell is known
                   you need to have the ingoing substance flux ChannelQS (mass/s)
                   and it will give outgoing flux ChannelQSn (mass/s)
                   and the current amount ChannelSubs (mass) in suspension+solution
             */

        }
    }

    ChannelQn->cover(LDD, 0);
    ChannelQsn->cover(LDD, 0);
    // avoid missing values around channel for adding to Qn for output

    double mb = 0;
    double n = 0;
    tm->fill(0);

    FOR_ROW_COL_MV_CH
    {
        if (SwitchChannelFlood)
        {
            if (ChannelMaxQ->Drc > 0)
                ChannelQn->Drc = qMin(ChannelQn->Drc, ChannelMaxQ->Drc);
        }
        // limit channel Q when culverts > 0

        double ChannelArea = ChannelAlpha->Drc*pow(ChannelQn->Drc, 0.6);
        tm->Drc = ChannelArea;

        ChannelWH->Drc = ChannelArea/((ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2.0);
        // water height is not used except for output! i.e. watervolume is cycled

        //double diff = Channelq->Drc*_dt + ChannelWaterVol->Drc - (ChannelArea * ChannelDX->Drc) - ChannelQn->Drc*_dt;
        double diff = QinKW->Drc*_dt + ChannelWaterVol->Drc - (ChannelArea * ChannelDX->Drc) - ChannelQn->Drc*_dt;
        //difference between fluxes and store in and out of channel cell


        if (!SwitchChannelInfil)
        {
            difkin->Drc += 0;//diff;
            mb += diff;
            if (ChannelArea > 0)
                n+=1;
        }

        if (SwitchBuffers && ChannelBufferVol->Drc > 0)
        {
            //qDebug()<< ChannelBufferVol->Drc << Channelq->Drc*_dt << ChannelWaterVol->Drc << (ChannelArea * ChannelDX->Drc) << ChannelQn->Drc*_dt<< diff;
        }
            else
            if (SwitchChannelInfil)
                InfilVolKinWave->Drc += diff;
        //VJ 110111 add channel infil to infil for mass balance
    }

    // mass balance correction, throw error on cells with WH
    if (n > 0)
        mb = mb/n;
    FOR_ROW_COL_MV_CH
    {
        ChannelWaterVol->Drc = tm->Drc * ChannelDX->Drc;
        if (ChannelWaterVol->Drc > 0)
            ChannelWaterVol->Drc = max(0, ChannelWaterVol->Drc+mb);
        ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelDX->Drc*0.5*(ChannelWidthUpDX->Drc+ChannelWidth->Drc));
        ChannelQn->Drc = qPow((ChannelWaterVol->Drc/ChannelDX->Drc)/ChannelAlpha->Drc, (1/0.6));
    }

    FOR_ROW_COL_MV_CH
    {
        //ChannelWaterVol->Drc = ChannelWH->Drc * (ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2.0 * ChannelDX->Drc;
        ChannelWaterVol->Drc = tm->Drc * ChannelDX->Drc;
        // total water vol after kin wave in m3, going to the next timestep
        // in a buffer ChannelArea = 0 so channelvolume is also 0

        if (SwitchErosion)
        {
            //ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc);
            // correct for very high concentrations, max 850 g/l
            // NOTE: nothing is done with this concentration, only for display,
            // so the display shows the conc after the kin wave
            ChannelConc->Drc = (ChannelQ->Drc > 1e-6 ? ChannelQs->Drc/ChannelQ->Drc : 0);
            // CHANGED, MORE STABLE CONC 19/9/13
        }
    }

    ChannelWH->cover(LDD, 0);

    FOR_ROW_COL_MV_CH
    {
        maxChannelflow->Drc = max(maxChannelflow->Drc, ChannelQn->Drc);
        maxChannelWH->Drc = max(maxChannelWH->Drc, ChannelWH->Drc);
    }


}
//---------------------------------------------------------------------------
