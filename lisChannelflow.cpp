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
            if (ChannelMaxQ->Drc > 0)
                Perim = FW + 2*wh + FW;
            // box culvert more friction
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
}
//---------------------------------------------------------------------------
//! calc channelflow, ChannelDepth, kin wave
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

        if (SwitchChannelFlood)
        {
            if (ChannelMaxQ->Drc > 0)
                ChannelWH->Drc = min(ChannelDepth->Drc-0.01, ChannelWH->Drc);
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

    // NOTE
    // hmx is flood level is looped
    // hmx receives net rainfall and decreases with infiltration in infiltration module
    // there is no runoff into the channel if it is flooding

    // get flood level in channel from 1D kin wave channel
    FOR_ROW_COL_MV
    {
        if (ChannelDepth->Drc > 0)
            hmx->Drc = max(0, ChannelWH->Drc - ChannelDepth->Drc)*ChannelWidthUpDX->Drc/_dx;
        // note: ChannelDepth lets you also control which channels flood: those that are 0 react as usual
    }

    double courant_number = courant_factor;
    double froude_limit = 0.8;
    double gravity = 9.81;
    double h_min = 1e-6;
    double timestep = 0;
    double timesum = 0;
    double sumh_t, sumh_t1, diff, cells;

    sumh_t = hmx->mapTotal();
    // if there is no flood skip everything
    if (sumh_t > 0)
    {
        // do one _dt with varying timestep based on courant condition
        do {
            // make Hydraulic head, gravity (dem+barriers) + water level
            FOR_ROW_COL_MV
            {
                Hmx->Drc = DEM->Drc + Barriers->Drc + hmx->Drc;
            }

            double maxdepth = qMax(0.01, hmx->mapMaximum());
            // find maxdepth
            timestep = courant_number*_dx/qSqrt(gravity*maxdepth);
            timestep = qMax(0.01, qMin(timestep, _dt-timesum));
            // determine timestep

            Qxsum->fill(0);
            // map with fluxes to/from central cell

            //flag which cells need processing, flood domain
            tma->fill(0);
            for (int i = 0; i < 9; i++)
            {
                FOR_ROW_COL_MV
                        if (r+dr[i] > 0 && r+dr[i] < _nrRows &&
                            c+dc[i] > 0 && c+dc[i] < _nrCols &&
                            !IS_MV_REAL8(&hmx->Drci))
                {
                    if (hx->Drci > 0)
                        tma->Drc = 1;
                }
            }

            // prepare maps in direction i;
            for (int i = 0; i < 9; i++)
            {
                double _dx2 = (((i+1) % 2) == 0 ? _dx : _dx*qSqrt(2));
                // diagonal cells have sqrt(2) dx

                // get the relevant maps and fill with 0 if outside
                //            Hx->fill(0);   // dem+wh
                //            hx->fill(0);   // wh
                //            Nx->copy(N);   // manning
                //            dHdLx->fill(0);
                FOR_ROW_COL_MV
                {
                    if (r+dr[i] > 0 && r+dr[i] < _nrRows &&
                            c+dc[i] > 0 && c+dc[i] < _nrCols &&
                            !IS_MV_REAL8(&hmx->Drci))
                    {
                        Hx->Drc = DEM->Drci + hmx->Drci;
                        hx->Drc = hmx->Drci;
                        Nx->Drc = N->Drci;
                        dHdLx->Drc = (Hx->Drc - Hmx->Drc)/_dx2;
                    }
                    else
                    {
                        Hx->Drc = 0;
                        hx->Drc = 0;
                        Nx->Drc = N->Drc;
                        dHdLx->Drc = 0;
                    }
                }

                // calc and sum fluxes in 8 directions, 4 = central cell
                FOR_ROW_COL_MV
                        if (tma->Drc == 1)
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
                    qlx2 = (hxi - hmx->Drc)*_dx/timestep;
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
                    // save flux in direction i for next flood timestep
                } // for all flood cells
            } // for i = 1 to 9


            // add sum q*dt (= m/s * s) to h and h > 0
            FOR_ROW_COL_MV
                    if (tma->Drc == 1)
            {
                hmx->Drc += timestep*Qxsum->Drc;
                hmx->Drc = max(0, hmx->Drc);
            }

            // simple correction, flow cannot lead to water level rise above neighbours,
            // water cannot flow up in this simplified solution (no momentum)
            FOR_ROW_COL_MV
                    if (tma->Drc == 1)
            {
                double hmax = 0;
                for (int i = 0; i < 9; i++)
                    if (r+dr[i] > 0 && r+dr[i] < _nrRows &&
                            c+dc[i] > 0 && c+dc[i] < _nrCols &&
                            i != 4 &&                     // not the central cell
                            !IS_MV_REAL8(&hmx->Drci))
                    {
                        hmax = qMax(hmax, hmx->Drci);
                        // find the highest water level around centre cell
                    }
                hmx->Drc = qMin(hmx->Drc, hmax);
            }

            // find current flood domain and nr of flooded cells
            cells = 1;
            FOR_ROW_COL_MV
            {
                if (hmx->Drc > 0)
                {
                    FloodDomain->Drc = 1;
                    cells++;
                }
                else
                    FloodDomain->Drc = 0;
            }

            //        int n = 0;
            //        do {
            //            n++;
            //            sumh_t1 = hmx->mapTotal();
            //            diff = (sumh_t1 > 0 ? (sumh_t - sumh_t1)/sumh_t1 : 0);
            //            FOR_ROW_COL_MV
            //            {
            //                double fh = hmx->Drc * diff;
            //                hmx->Drc = max(hmx->Drc + fh, 0);
            //            }
            //        } while (n < 5);

            timesum = timesum + timestep;
            // sum to reach _dt

        } while (timesum  < _dt);
        // continue while _dt is not reached
        sumh_t1 = hmx->mapTotal();
        diff = (sumh_t - sumh_t1)/cells;
        qDebug() << sumh_t << sumh_t1 << diff;
        // echo to screen
        debug(QString("Flooding: %1 %2 avg err h in flooded cells %3").arg(sumh_t,8,'f',3).arg(sumh_t1,8,'f',3).arg(diff,8,'f',3));
    }

    // report flood maps
    hmx->report("hmx");

    // put new flood level in channel for next 1D kin wave channel
    FOR_ROW_COL_MV_CH
    {
        if (ChannelDepth->Drc > 0 && hmx->Drc > ChannelDepth->Drc)
        {
            ChannelWH->Drc = hmx->Drc + ChannelDepth->Drc;
            ChannelWaterVol->Drc = ChannelWH->Drc *
                    ((ChannelWidthUpDX->Drc+ChannelWidth->Drc)/2) * ChannelDX->Drc;
        }
    }
}
