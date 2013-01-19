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

    //    ChannelQn->setMV();

    //    if (useSorted)
    //    {
    //        KinematicSorted(lddlistch, lddlistchnr, ChannelQ, ChannelQn, ChannelQs, ChannelQsn, Channelq, ChannelAlpha, ChannelDX,
    //                        ChannelWaterVol, ChannelSed, ChannelBufferVol, ChannelBufferSed);
    //    }
    //    else
    //    {
    FOR_ROW_COL_MV_CH
    {
        if (LDDChannel->Drc == 5)
        {
            Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQs, ChannelQsn, Channelq, ChannelAlpha, ChannelDX,
                      ChannelWaterVol, ChannelSed, ChannelBufferVol, ChannelBufferSed);

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
    // }

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
#define Drci Data[r+dr[i]][c+dc[i]]

void TWorld::ChannelFlood(void)
{
    if (!SwitchIncludeChannel)
        return;
    if (!SwitchChannelFlood)
        return;

    int dc[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dr[9] = {1,1,1,  0, 0, 0,  -1, -1, -1};

    FOR_ROW_COL_MV
    {
        if (ChannelHeight->Drc > 0)
            hmx->Drc = max(0, ChannelWH->Drc - ChannelHeight->Drc)*ChannelWidthUpDX->Drc/_dx;
    }

    double courant_number = 0.1;
    double froude_limit = 0.8;
    double gravity = 9.81;
    double h_min = 1e-6;
    double timestep = 0;
    double timesum = 0;
    double sumh_t = hmx->mapTotal();
    double diff;

    // do one _dt with varying timestep based on courant condition
    do {
        sumh_t = hmx->mapTotal();

        Hmx->calc2Maps(DEM, hmx, ADD);

        double maxdepth = qMax(0.01, hmx->mapMaximum());
        // find maxdepth
        timestep = courant_number*_dx/qSqrt(gravity*maxdepth);
        timestep = qMax(0.01, qMin(timestep, _dt-timesum));
        // determine timestep

        Qxsum->fill(0);
        // map with fluxes to/from central cell
        FloodDomain->fill(0);

        // prepare maps in direction i;
        for (int i = 0; i < 9; i++)
        {
            double _dx2 = (((i+1) % 2) == 0 ? _dx : _dx*qSqrt(2));

            // get the relevant maps and fill with 0 if outside
            Hx->fill(0);   // dem+wh
            hx->fill(0);   // wh
            Nx->copy(N);   // manning
            dHdLx->fill(0);
            FOR_ROW_COL_MV
                    if (r+dr[i] > 0 && r+dr[i] < _nrRows &&
                        c+dc[i] > 0 && c+dc[i] < _nrCols &&
                        !IS_MV_REAL8(&hmx->Drci))
            {
                Hx->Drc = DEM->Drci + hmx->Drci;
                hx->Drc = hmx->Drci;
                Nx->Drc = N->Drci;
                dHdLx->Drc = (Hx->Drc - Hmx->Drc)/_dx2;
            }

            // calc and sum fluxes in 8 directions, 4 = central cell
            FOR_ROW_COL_MV
            {
                double Qx = 0, qlx = 0, qlx1 = 0, qlx2 = 0;

                double dHdLxi = dHdLx->Drc;
                double signx = (dHdLxi < 0 ? -1.0 : 1.0);
                double hxi = hx->Drc;
                double qxi = qx[i].m->Drc;
                double NN = Nx->Drc;

                if (hxi > h_min)
                    qlx = (qxi - (gravity*hxi*timestep*dHdLxi))/
                            (1.0 + (gravity*hxi*timestep*NN*NN*qxi)/qPow(hxi, 10.0/3.0));
                else
                    qlx = 0;
                // explicit solution of saint venant equation (m2/s), Bates et al. 2010
                qlx1 = hxi*qSqrt(gravity*hxi)*froude_limit;
                // limit max flux to h * wave velocity * froude number (m * m/s = m2/s)
                qlx2 = (hxi - hmx->Drc)*_dx2/timestep;
                // limit max flux to width * water difference with central cell (m2/s)
                Qx = signx * qMin(qMin(qAbs(qlx1), qAbs(qlx2)), qAbs(qlx));
                // Qx is the min of all possible fluxes, preserve sign
                if (i != 4)
                    Qxsum->Drc += Qx/_dx2*0.5;
                // add fluxes of 8 directions, qsum has unit m2/s /m = m/s
                // 0.5 accounts for the fact that the central cell has a boundary
                // of 4 sides is 4*dx, touching only EW and NS directions. Adding
                // four more diagonal fluxes would cause twice the flow, so all fluxes
                // are assumed to have a width of 0.5*dx * 8 = 4 dx.
                // using 8 instead of 4 directions seems to give a much better flow
                // for 4 directions do i += 2 instead of i++

                //QV = max(QV, abs(qx*0.5));
                // velocity calculated below with maximum of fluxes to/from central cell

                if (i == 4)
                    qx[i].m->Drc = signx * qMin(qlx, qlx1);
                else
                    qx[i].m->Drc = Qx;
                // save flux in direction i, min of st venant and wave velocity flux
            }
        } // i = 1 to 9

        FOR_ROW_COL_MV
        {
            hmx->Drc += timestep*Qxsum->Drc;
            // add sum q*dt (= m/s * s) to h
            hmx->Drc = max(0, hmx->Drc);
            // not below 0
        }


        // correction, flow cannot lead to water level rise above neaghbours,
        // water cannot flow up in this simplified solution (no momentum)
        FOR_ROW_COL_MV
        {
            double hmax = 0;
            for (int i = 0; i < 9; i++)
                if (r+dr[i] > 0 && r+dr[i] < _nrRows &&
                        c+dc[i] > 0 && c+dc[i] < _nrCols &&
                        i != 4 && !IS_MV_REAL8(&hmx->Drci))
                {
                    hmax = qMax(hmax, hmx->Drci);
                    // find the hieghest water level around centre cell
                }
            hmx->Drc = qMin(hmx->Drc, hmax);
        }

        FOR_ROW_COL_MV
        {
            if (hmx->Drc > 0)
                FloodDomain->Drc = 1;
            else
                FloodDomain->Drc = 0;
        }
        double cells = qMax(1.0, FloodDomain->mapTotal());
        // wet cells

        // calculate mass balance and correct h with average error over wet cells, instead of true iteration
        double sumh_t1 = sumh_t1 = hmx->mapTotal();
        diff = (sumh_t - sumh_t1)/cells;

        qDebug() << timesum << timestep  << diff << sumh_t << sumh_t1 ;

        timesum = timesum + timestep;
        // sum to reach _dt

    } while (timesum  < _dt);
    // continue while _dt is not reached

//    FOR_ROW_COL_MV
//            if (FloodDomain->Drc > 0)
//    {
//        WH->Drc = hmx->Drc + WHstore->Drc;
//        WaterVolall->Drc = DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc );
//    }

    hmx->report("h0toc");
    FloodDomain->report("FD");
    tmb->calc2Maps(DEM, hmx, ADD);
    tmb->report("hxdem");
tma->fill(0);
    FOR_ROW_COL_MV_CH
    {
        if (ChannelHeight->Drc > 0 && hmx->Drc > ChannelHeight->Drc)
        {
            tma->Drc = ChannelWH->Drc - (hmx->Drc + ChannelHeight->Drc);
//            ChannelWH->Drc = hmx->Drc + ChannelHeight->Drc;
//            ChannelWaterVol->Drc = ChannelWH->Drc *
//                    ((ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2) * ChannelDX->Drc;
        }
    }
    tma->report("diff");
}
