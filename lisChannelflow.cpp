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
- void TWorld::ChannelFlow(void) calculate channelflow, channelheight, do kinematic wave \n
*/

#include "model.h"

//---------------------------------------------------------------------------
// V, alpha and Q in the channel
void TWorld::CalcVelDischChannel(void)
{
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
        double dw = 0.5*(ChannelWidthUpDX->Drc - FW); // extra width when non-rectamgular

        if (dw > 0)
        {
            //			Perim = FW + 2*sqrt(wh*wh + dw*dw);
            Perim = FW + 2*wh/cos(atan(ChannelSide->Drc));
            Area = FW*wh + wh*dw;// NOT*2 BECAUSE TRIANGLE !!;
        }
        else
        {
            Perim = FW + 2*wh;
            Area = FW*wh;
        }
        //Perim = ChannelWidth->Drc + 2*wh/cos(atan(ChannelSide->Drc));
        // cos atanb more expensive than sqrt ?
        //Area = ChannelWidth->Drc*wh + wh*(ChannelWidthUpDX->Drc - ChannelWidth->Drc);

        ChannelPerimeter->Drc = Perim;
        //VJ 110109 needed for channel infil

        if (Perim > 0)
            Radius = Area/Perim;
        else
            Radius = 0;

        ChannelAlpha->Drc = qPow(ChannelN->Drc/grad * powl(Perim, _23),beta);

        if (ChannelAlpha->Drc > 0)
            ChannelQ->Drc = qPow(Area/ChannelAlpha->Drc, beta1);
        else
            ChannelQ->Drc = 0;

        ChannelV->Drc = pow(Radius, _23)*grad/ChannelN->Drc;
    }
    //   else // WHY THIS??? RELIC !!!
    //   {
    //      ChannelAlpha->Drc = 0;
    //      ChannelQ->Drc = 0;
    //      ChannelV->Drc = 0;
    //      ChannelPerimeter->Drc = 0;
    //   }
}
//---------------------------------------------------------------------------
//! calc channelflow, channelheight, kin wave
void TWorld::ChannelFlow(void)
{
    if (!SwitchIncludeChannel)
        return;

    // calculate new channel WH , WidthUp and Volume
    FOR_ROW_COL_MV_CH
    {
        /*---- Water ----*/

        ChannelQsn->Drc = 0;
        Channelq->Drc = 0;
        ChannelWH->Drc = 0;

        ChannelWaterVol->Drc += RunoffVolinToChannel->Drc;
        // add inflow to channel
        ChannelWaterVol->Drc += Rainc->Drc*ChannelWidthUpDX->Drc*DX->Drc;
        // add rainfall in m3, no interception, rainfall so do not use ChannelDX

        if (SwitchBuffers && ChannelBufferVol->Drc > 0)
        {
            ChannelBufferVol->Drc -= ChannelWaterVol->Drc;
            ChannelWaterVol->Drc = 0;
            // add inflow from slopes and rainfall to buffer
        }

        if (ChannelSide->Drc == 0)// && ChannelWidth->Drc > 0)// rectangular channel
        {
            ChannelWidthUpDX->Drc = ChannelWidth->Drc;
            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
        }
        else  // non-rectangular
        {
            /*
   ABC fornula
    dw      w       dw
   \  |            |  /
    \ |          h |a/  <= tan(a) is channelside = tan angle of side wall
     \|____________|/
   area = vol/dx = h*w + h*dw
   dw = h*tan(a)
   vol/dx = w*h + dw*h = w*h + tan(a)*h*h
   tan(a)h^2 + wh - vol/dx = 0
     aa (h2)   +   bb(h) +  cc
*/
            double aa = ChannelSide->Drc;  //=tan(a)
            double bb = ChannelWidth->Drc; //=w
            double cc = -ChannelWaterVol->Drc/ChannelDX->Drc; //=vol/DX

            ChannelWH->Drc = (-bb + sqrt(bb*bb - 4*aa*cc))/(2*aa);
            if (ChannelWH->Drc < 0)
            {
                ErrorString = QString("channel water height is negative at row %1, col %2").arg(r).arg(c);
                throw 1;
            }
            ChannelWidthUpDX->Drc = ChannelWidth->Drc + 2*ChannelSide->Drc*ChannelWH->Drc;
        }

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

        if (RoadWidthDX->Drc > 0)
            ChannelWidthUpDX->Drc = min(0.9*_dx-RoadWidthDX->Drc, ChannelWidthUpDX->Drc);
        // channel cannot be wider than _dx-road
        /* TODO zit al in gridcell, nodig hier? */

        if (SwitchChannelInfil)
        {
            Channelq->Drc =  -(ChannelKsat->Drc *  ChannelPerimeter->Drc/3600000.0);
            //mm/h / 1000 = m/h / 3600 = m/s * m = m2/s
        }
        // NOTE: for buffers channelksat = 0

    }

    CalcVelDischChannel();
    // alpha, V and Q from Manning

    /*---- Sediment ----*/

    if (SwitchErosion)
    {
        ChannelFlowDetachment();

        FOR_ROW_COL_MV_CH
        {
            ChannelQs->Drc = ChannelQ->Drc * ChannelConc->Drc;
        }
    }

    ChannelQn->setMV();

    if (useSorted)
    {
        KinematicSorted(lddlistch, lddlistchnr, ChannelQ, ChannelQn, ChannelQs, ChannelQsn, Channelq, ChannelAlpha, ChannelDX,
                        ChannelWaterVol, ChannelSed, ChannelBufferVol, ChannelBufferSed);
    }
    else
    {
        FOR_ROW_COL_MV_CH
        {
            if (LDDChannel->Drc == 5)
            {
                Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQs, ChannelQsn, Channelq, ChannelAlpha, ChannelDX,
                          ChannelWaterVol, ChannelSed, ChannelBufferVol, ChannelBufferSed);
                // for checking:
                routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQs, ChannelQsn, ChannelAlpha, ChannelDX, ChannelWaterVol, ChannelSed);


                /*
                   routing of substances add here!
                   do after kin wave so that the new flux ChannelQn out of a cell is known
                   you need to have the ingoing substance flux ChannelQS (mass/s)
                   and it will give outgoing flux ChannelQSn (mass/s)
                   and the current amount ChannelSubs (mass) in suspension+solution
                */
                //routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQS, ChannelQSn, ChannelAlpha, ChannelDX, ChannelWaterVol, ChannelSubs);

            }
        }
    }

    ChannelQn->cover(LDD, 0); // avoid missing values around channel for adding to Qn for output
    ChannelQsn->cover(LDD, 0);

    FOR_ROW_COL_MV_CH
    {
        double ChannelArea = ChannelAlpha->Drc*pow(ChannelQn->Drc, 0.6);
        // in buffers ChannelQn = 0;
        //qDebug() << ChannelArea << ChannelWidthUpDX->Drc << ChannelWidth->Drc;

        ChannelWH->Drc = ChannelArea/((ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2);
        // water height is not used except for output i.e. watervolume is cycled

        double diff = Channelq->Drc*_dt + ChannelWaterVol->Drc - (ChannelArea * ChannelDX->Drc) - ChannelQn->Drc*_dt;
        //difference between fluxes and store in and out of channel cell
        // qDebug() << diff;
        if (SwitchBuffers && ChannelBufferVol->Drc > 0)
        {
            //qDebug()<< ChannelBufferVol->Drc << Channelq->Drc*_dt << ChannelWaterVol->Drc << (ChannelArea * ChannelDX->Drc) << ChannelQn->Drc*_dt<< diff;
        }
        else
            if (SwitchChannelInfil)
                InfilVolKinWave->Drc += diff;
        //VJ 110111 add channel infil to infil for mass balance

        ChannelWaterVol->Drc = ChannelArea * ChannelDX->Drc;
        // total water vol after kin wave in m3, going to the next timestep
        // in a buffer ChannelArea = 0 so channelvolume is also 0

        if (SwitchErosion)
        {
            ChannelConc->Drc = MaxConcentration(ChannelWaterVol->Drc, ChannelSed->Drc);
            // correct for very high concentrations, max 850 g/l
            // NOTE: nothing is done with this concentration, only for display,
            // so the display shows the conc after the kin wave
        }
    }
    ChannelWH->cover(LDD, 0);

}
//---------------------------------------------------------------------------
void TWorld::ChannelFlood(void)
{
    if (!SwitchIncludeChannel)
        return;
    if (!SwitchChannelFlood)
        return;

    tm->fill(0); //channel waterheight above water surface
    tma->fill(0); //dem+wh
    tmb->fill(0);

    ChannelHeight->cover(LDD,0);

    // channel overflow water sticking out above flood water level
    FOR_ROW_COL_MV_CH
    {
        tm->Drc = max(0, ChannelWH->Drc - ChannelHeight->Drc - WHflood->Drc);
        tm->Drc = 1.0;//tm->Drc*_dx/tmb->Drc;
    }

    FOR_ROW_COL_MV
    {
        if (ChannelHeight->Drc > 0)
            WHflood->Drc = 1.0;
        DEMflood->Drc = DEM->Drc + WHflood->Drc;// + tm->Drc;
    }
    // add surplus water to DEM

    double courant_number = 0.7;
    double froude_limit = 0.8;
    double local_time_factor = time_factor;
    double hflow = 0;
    //double edgeslope = 0.001;
    double hflow_threshold = 0.001;
    double time_factor = 1.0;
    int nrruns = 1;

    local_time_factor = time_factor;
    double maxdepth = qMin(10.0, qMax(0.1, WHflood->mapMaximum()));

    if (local_time_factor > (courant_number * _dx / qSqrt(9.8*maxdepth) ))
        local_time_factor = courant_number * _dx / qSqrt(9.8*maxdepth);

    if (local_time_factor > _dt)
        local_time_factor = _dt;
    else
    {
        nrruns = qFloor(_dt/(local_time_factor+0.1))+1;
        local_time_factor = _dt/(double)nrruns;
    }

    qDebug() << nrruns << local_time_factor;
    // calculate fluxes qx and qy
    DEBUG(QString("flood timestep %1, nr floodsteps %2, maxdepth %3.").arg(local_time_factor).arg(nrruns).arg(maxdepth));

    int ymax = _nrRows-1;
    int xmax = _nrCols-1;
    int y, x;



    for (int nr = 0; nr < nrruns; nr++)
    {
        for(y = 0; y <= ymax; y++)
            for(x = 0; x <= xmax; x++)
                tmb->Data[x][y] = 0;

        for(y =1; y <= ymax; y++)
        {
            int inc = 1;
            for(x = xmax; x >= 1; x--)
            {
                if ((!IS_MV_REAL8(&WHflood->Data[x][y]) &&
                     !IS_MV_REAL8(&WHflood->Data[x-1][y]) &&
                     !IS_MV_REAL8(&WHflood->Data[x+1][y]) &&
                     !IS_MV_REAL8(&WHflood->Data[x][y-1]) &&
                     !IS_MV_REAL8(&WHflood->Data[x][y+1]))
                        )
                {

                    if ((WHflood->Data[x][y]>0
                         || WHflood->Data[x-1][y] > 0
                         || WHflood->Data[x+1][y] > 0
                         || WHflood->Data[x][y-1] > 0
                         || WHflood->Data[x][y+1] > 0))
                    {
                        tmb->Data[y][inc] = (REAL8)x;
                        inc++;
                    }
                }
            }
        }
        tmb->report("tmb");

        //FOR_ROW_COL_MV
        for(y = 1; y <= ymax; y++)
        {
            int inc = 1;
            qDebug() << y << inc << tmb->Data[y][inc];
            while (tmb->Data[y][inc] > 0)
            {
                x = (int)tmb->Data[y][inc];
                inc++;
                // calc qx, in x direction
                if ((WHflood->Data[x][y] > 0 || WHflood->Data[x - 1][y] > 0) && !IS_MV_REAL8(&DEM->Data[x-1][y]))
                {
                    hflow = qMax(DEMflood->Data[x][y], DEMflood->Data[x-1][y]) - qMax(DEM->Data[x][y], DEM->Data[x-1][y]);

                    if (hflow > hflow_threshold)
                    {
                        double tempslope = (DEMflood->Data[x][y] - DEMflood->Data[x-1][y])/_dx;

                        //                if (x == xmax) tempslope = edgeslope;
                        //                if (x <= 2) tempslope = 0 - edgeslope;

                        qx->Data[x][y] = (qx->Data[x][y] - (9.8 * hflow * local_time_factor * tempslope)) /
                                (1 + 9.8 * hflow * local_time_factor * N->Data[x][y]*N->Data[x][y] * qAbs(qx->Data[x][y])/qPow(hflow,10/3));

                        //if (oldqx != 0) qx->Data[x][y] = (oldqx + qx->Data[x][y]) / 2;

                        // need to have these lines to stop too much water moving from one cell to another
                        // - resulting in -ve discharges
                        // which causes a large instability to develop - only in steep catchments really
                        if (qx->Data[x][y] > 0 && (qx->Data[x][y]/hflow)/qSqrt(9.8*hflow) > froude_limit )
                            qx->Data[x][y] = hflow * qSqrt(9.8*hflow) * froude_limit;
                        if (qx->Data[x][y] < 0 && qAbs(qx->Data[x][y]/hflow)/qSqrt(9.8*hflow) > froude_limit )
                            qx->Data[x][y] = 0 - hflow * qSqrt(9.8*hflow) * froude_limit;

                        if (qx->Data[x][y] > 0 && (qx->Data[x][y]*local_time_factor/_dx) > (WHflood->Data[x][y]/4.0))
                            qx->Data[x][y] = ((WHflood->Data[x][y]*_dx)/5.0)/local_time_factor;
                        if (qx->Data[x][y] < 0 && qAbs(qx->Data[x][y] * local_time_factor/_dx) > (WHflood->Data[x-1][y]/4.0))
                            qx->Data[x][y] = 0 - ((WHflood->Data[x-1][y]*_dx)/5.0)/local_time_factor;
                    }
                    else
                        qx->Data[x][y] = 0;
                }

                //                // calc qy, in y direction
                if ((WHflood->Data[x][y] > 0 || WHflood->Data[x][y-1] > 0) && !IS_MV_REAL8(&DEM->Data[x][y-1]))
                {

                    hflow = qMax(DEMflood->Data[x][y], DEMflood->Data[x][y-1]) - qMax(DEM->Data[x][y], DEM->Data[x][y-1]);

                    if (hflow > hflow_threshold)
                    {
                        double tempslope = (DEMflood->Data[x][y] - DEMflood->Data[x][y-1])/_dx;

                        //                if (x == xmax) tempslope = edgeslope;
                        //                if (x <= 2) tempslope = 0 - edgeslope;

                        //double oldqx = qx->Data[x][y];
                        qy->Data[x][y] = (qy->Data[x][y] - (9.8 * hflow * local_time_factor * tempslope)) /
                                (1 + 9.8 * hflow * local_time_factor * N->Data[x][y]*N->Data[x][y] * qAbs(qy->Data[x][y])/qPow(hflow,10/3));

                        //if (oldqy != 0) qy->Data[x][y] = (oldqy + qy->Data[x][y]) / 2;

                        // need to have these lines to stop too much water moving from one cell to another
                        // - resulting in -ve discharges
                        // which causes a large instability to develop - only in steep catchments really
                        if (qy->Data[x][y] > 0 && (qy->Data[x][y]/hflow)/qSqrt(9.8*hflow) > froude_limit )
                            qy->Data[x][y] = hflow * qSqrt(9.8*hflow) * froude_limit;
                        if (qy->Data[x][y] < 0 && qAbs(qy->Data[x][y]/hflow)/qSqrt(9.8*hflow) > froude_limit )
                            qy->Data[x][y] = 0 - hflow * qSqrt(9.8*hflow) * froude_limit;

                        if (qy->Data[x][y] > 0 && (qy->Data[x][y]*local_time_factor/_dx) > (WHflood->Data[x][y]/4.0))
                            qy->Data[x][y] = ((WHflood->Data[x][y]*_dx)/5.0)/local_time_factor;
                        if (qy->Data[x][y] < 0 && qAbs(qy->Data[x][y] * local_time_factor/_dx) > (WHflood->Data[x][y-1]/4.0))
                            qy->Data[x][y] = 0 - ((WHflood->Data[x][y-1]*_dx)/5.0)/local_time_factor;
                    }
                    else
                        qy->Data[x][y] = 0;
                }

                //            }
                //update water heights with fluxes qx and qy

                for(y = 1; y <= ymax; y++)
                {
                    int inc = 1;
                    while (tmb->Data[y][inc] > 0)
                    {
                        x = (int)tmb->Data[y][inc];
                        inc++;
                        WHflood->Data[x][y] += local_time_factor * (qx->Data[x + 1][y] - qx->Data[x][y] + qy->Data[x][y + 1] - qy->Data[x][y]) / _dx;

                    }

                }
                WHflood->report("whf");
                qx->report("qx");
                qy->report("qy");
            }
        }
    }
}
//---------------------------------------------------------------------------



//  ChannelHeight->cover(LDD,0);
//  actFloodArea->copy(floodArea);
//  // initialize floodarea

//  FOR_ROW_COL_MV
//      DEMflood->Drc = DEM->Drc + WHflood->Drc;
//  DEMflood->report("dfld");
//  tmb->copy(ChannelWidthUpDX);
//  tmb->cover(LDD, 0);
//  FOR_ROW_COL_MV_CH
//  {
//    tm->Drc = max(0, ChannelWH->Drc - ChannelHeight->Drc - WHflood->Drc);
//    tm->Drc = tm->Drc*_dx/tmb->Drc;
//  }
//  tm->cover(LDD, 0);
//  // channel water rising above current flood plain
//  tm->report("ch");

//  FOR_ROW_COL_MV
//      tma->Drc = DEMflood->Drc + tm->Drc;
//  // initial water + DEM

//  for (int i = 0; i < 12; i++)
//    {
//      tmb->copy(tma);
//      // reset to original level
//      tmb->areaAverage(actFloodArea);
//      // average dem+water level in flooded area

//      FOR_ROW_COL_MV
//      {
//        if(actFloodArea->Drc == 0)
//          tmb->Drc = DEMflood->Drc;
//      }
//      // Dem = if(Lake ne 0,areaaverage(DemPit ,Lake),DemDam);
//      // set to DEM outside flood area

//      FOR_ROW_COL_MV
//      {
//        if (tmb->Drc > DEMflood->Drc)
//          actFloodArea->Drc = floodArea->Drc;
//        else
//          actFloodArea->Drc = 0;
//      }
//      // Lake = if(Dem > DemDam, PotLake, 0);
//      // adapt flood area

//      DEMflood->copy(tmb);

//      FOR_ROW_COL_MV
//      {
//        WHflood->Drc = max(0, DEMflood->Drc - DEM->Drc);
//        Vflood->Drc = _dt*pow(WHflood->Drc,(2/3))*sqrt(Grad->Drc)/N->Drc;
//      }
//      tmb->copy(Vflood);
//      tmb->areaAverage(actFloodArea);
//      FOR_ROW_COL_MV
//      {
//        if (tmb->Drc > floodDist->Drc)
//          actFloodArea->Drc = 0;
//      }
//    }

//  WHflood->report("whf");
//  DEMflood->report("dflda");
//  actFloodArea->report("area");
//    FOR_ROW_COL_MV_CH
//    {
//      if (tm->Drc > 0)
//        ChannelWH->Drc = ChannelHeight->Drc + WHflood->Drc;
//    }
//    FOR_ROW_COL_MV
//    {
//      if (actFloodArea->Drc > 0)
//        WH->Drc = WHflood->Drc;
//    }
//    WH->report("wh");
//}

//binding
// DemDam = dem10m.map;
// PotLake = floodarea.map;
// WL = ch;

//timer
//1 200 1;

//initial
//Dem = DemDam;

//dynamic

// DemPit = timeinput(WL)+Dem;
// Lake = PotLake;

// Dem = if(Lake ne 0,areaaverage(DemPit ,Lake),DemDam);
// Lake = if(Dem > DemDam, PotLake, 0);
// Dem = if(Lake ne 0,areaaverage(DemPit ,Lake),DemDam);
// Lake = if(Dem > DemDam, PotLake, 0);
// Dem = if(Lake ne 0,areaaverage(DemPit ,Lake),DemDam);
// Lake = if(Dem > DemDam, PotLake, 0);
// Dem = if(Lake ne 0,areaaverage(DemPit ,Lake),DemDam);
// Lake = if(Dem > DemDam, PotLake, 0);
//## Dem = if(Lake ne 0,areaaverage(DemPit ,Lake),DemDam);
//## Lake = if(Dem > DemDam, PotLake, 0);
//## Dem = if(Lake ne 0,areaaverage(DemPit ,Lake),DemDam);
//## Lake = if(Dem > DemDam, PotLake, 0);
// report Dem = if(Lake ne 0,areaaverage(DemPit ,Lake),DemDam);
// report Lake = if(Dem > DemDam, PotLake,0);
