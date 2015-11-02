
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
  \file lisOverlandflow.cpp
  \brief calculate fraction flowing in the channel, Q, V and call kin wave

functions: \n
- void TWorld::ToChannel(void) \n
- void TWorld::CalcVelDisch(void) \n
- void TWorld::OverlandFlow(void) \n
 */

#include <algorithm>
#include "model.h"
#include "operation.h"
#define tiny 1e-8

//---------------------------------------------------------------------------
void TWorld::RainfallToFlood(void)
{
    if (SwitchRainfallFlood)
    {
        FOR_CELL_IN_FLOODAREA
                if ( Grad->Drc <= rainFloodingGradient)
        {
            // if it rains, and there is no flood, and it is flat, and there is sufficient runoff water, then this water kan turn to
            // flood directly!
            if (RainNet->Drc > 0 && WHrunoff->Drc > 0.03 && hmx->Drc == 0 && ChannelWidthUpDX->Drc == 0)
            {
                double dwh =  WHrunoff->Drc;

                hmx->Drc = dwh * FlowWidth->Drc/_dx;
                WH->Drc -= dwh;
                WHrunoff->Drc = 0;
                WHGrass->Drc -= dwh;
                WHroad->Drc -= dwh;
            }
        }}
}

}
//---------------------------------------------------------------------------
//fraction of water and sediment flowing into the channel
void TWorld::ToFlood(void)
{
    if (!SwitchChannelFlood)
        return;

    FOR_CELL_IN_FLOODAREA
            if (WHrunoff->Drc > 0.01 && hmx->Drc > 0.01 && ChannelWidthUpDX->Drc == 0)
    {
        double frac = 1-exp(-runoff_partitioning*hmx->Drc/WHrunoff->Drc);
        //std::min(1.0, std::max(0.0, exp(-runoff_partitioning*WH->Drc/hmx->Drc)));
        double dwh = frac * WHrunoff->Drc;

        hmx->Drc += dwh * FlowWidth->Drc/_dx;
        WH->Drc -= dwh;
        WHrunoff->Drc -= dwh;

        WHGrass->Drc -= dwh;
        WHroad->Drc -= dwh;

    }}
}
//---------------------------------------------------------------------------
//fraction of water and sediment flowing into the channel
void TWorld::ToChannel(void)
{
    if (!SwitchIncludeChannel)
        return;

    FOR_ROW_COL_MV_CH
    {
        double fractiontochannel;
        double Volume = WHrunoff->Drc * FlowWidth->Drc * DX->Drc;

        if (Volume == 0)
        {
            SedToChannel->Drc = 0;
            RunoffVolinToChannel->Drc = 0;
            continue;
        }

        if (ChannelAdj->Drc == 0)
            fractiontochannel = 1.0;
        else
            fractiontochannel = std::min(1.0, _dt*V->Drc/std::max(0.01*_dx,0.5*ChannelAdj->Drc));
        // fraction to channel calc from half the adjacent area width and flow velocity

        if (SwitchBuffers)
            if (BufferID->Drc > 0)
                fractiontochannel = 1.0;
        // where there is a buffer in the channel, all goes in the channel

        // cannot flow into channel is water level in channel is higher than depth
        if (SwitchChannelFlood)
        {
            if (WHrunoff->Drc <= std::max(ChannelLevee->Drc, ChannelWH->Drc-ChannelDepth->Drc))
                fractiontochannel = 0;
            // no inflow when flooded
            if (ChannelMaxQ->Drc > 0)
                fractiontochannel = 0;
            // no surface inflow when culverts and bridges
        }
        if (SwitchAllinChannel)
            if (LDD->Drc == 5)
                fractiontochannel = 1.0;
        // in catchment outlet cell, throw everything in channel

        RunoffVolinToChannel->Drc = fractiontochannel*Volume;
        // water diverted to the channel
        WHrunoff->Drc *= (1-fractiontochannel);

        WH->Drc = WHrunoff->Drc + WHstore->Drc;
        //VJ 130425


        if (SwitchErosion)
        {
            SedToChannel->Drc = fractiontochannel*Sed->Drc;
            //sediment diverted to the channel
            Sed->Drc -= SedToChannel->Drc;
            // adjust sediment in suspension
        }
    }
}
//---------------------------------------------------------------------------
void TWorld::CalcVelDisch()
{
    if(SwitchKinematic2D > 1)
    {
        return K2DCalcVelDisch();
    }

    FOR_ROW_COL_MV
    {
        double Perim;
        const double beta = 0.6;
        const double _23 = 2.0/3.0;
        double beta1 = 1/beta;
        //double kinvisc = 1.1e-6; // 15 degrees celcius water
        double NN = N->Drc;


        if (SwitchChannelFlood)
            NN = N->Drc * qExp(mixing_coefficient*hmx->Drc);
        // slow down water in flood zone
        //    tma->Drc = hmx->Drc * UVflood->Drc/kinvisc;
        // Reynolds number ==> turbulent

        // avg WH from soil surface and roads, over width FlowWidth
        Perim = 2*WHrunoff->Drc+FlowWidth->Drc;

        if (Perim > 0)
            R->Drc = WHrunoff->Drc*FlowWidth->Drc/Perim;
        else
            R->Drc = 0;

        Alpha->Drc = pow(NN/sqrt(Grad->Drc) * pow(Perim, _23),beta);

        if (Alpha->Drc > 0)
            Q->Drc = pow((FlowWidth->Drc*WHrunoff->Drc)/Alpha->Drc, beta1);
        else
            Q->Drc = 0;

        V->Drc = pow(R->Drc, _23)*sqrt(Grad->Drc)/NN;

        //tm->Drc = V->Drc * R->Drc/kinvisc;
        //Reynolds number
    }


}
//---------------------------------------------------------------------------
void TWorld::OverlandFlowNew(void)
{
    // recalculate water vars after subtractions in "to channel"
    FOR_ROW_COL_MV
    {
        WaterVolin->Drc = DX->Drc * FlowWidth->Drc * WHrunoff->Drc;
        //volume runoff into the kin wave, needed to determine infil in kin wave

        // WaterVolin total water volume in m3 before kin wave, WHrunoff may be adjusted in tochannel
        q->Drc = FSurplus->Drc*SoilWidthDX->Drc/_dt;  //???FlowWidth
        // infil flux in kin wave (<= 0)negative value), in m2/s, in kiv wave DX is used
        // surplus related to infiltrating surfaces
    }
    
    fill(*Qsn, 0.0);
    fill(*QinKW, 0.0);
    fill(*QoutKW, 0.0);
    // flag all new flux as missing value, needed in kin wave and replaced by new flux

    if(SwitchKinematic2D == K1D_METHOD)
    {

        if (SwitchErosion)
        {
            // calc seediment flux going in kin wave as Qs = Q*C
            FOR_ROW_COL_MV
            {
                Qs->Drc =  Q->Drc * Conc->Drc;
                // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
            }
        }

        Qn->setAllMV();
        FOR_ROW_COL_MV
        {
            if (LDD->Drc == 5) // if outflow point, pit
            {

                Kinematic(r,c, LDD, Q, Qn, Qs, Qsn, q, Alpha, DX, WaterVolin, Sed, BufferVol, BufferSed);
                //VJ 110429 q contains additionally infiltrated water volume after kin wave in m3
            }
        }
        //
        //      routing of substances add here!
        //      do after kin wave so that the new flux Qn out of a cell is known
        //      you need to have the ingoing substance flux QS (mass/s)
        //      and it will give outgoing flux QSn (mass/s)
        //      and the current amount Subs (mass) in suspension+solution
        //

        if (SwitchPesticide)
        {
            // calc pesticide flux going in kin wave as Qp = Q*C
            FOR_ROW_COL_MV
            {
                Qp->Drc =  Q->Drc * C->Drc;
                // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
            }

            fill(*Qpn, 0.0);
            FOR_ROW_COL_MV
            {
                if (LDD->Drc == 5) // if outflow point, pit
                {
                    routeSubstance(r,c, LDD, Q, Qn, Qp, Qpn, Alpha, DX, WaterVolin, Pest);
                }
            }
        }
    }
    else
    {
        //Kinematic wave solution in 2 dimensions
        //also includes sediment and pesticides transport
        //Functions are not generic since they are only used here
        //Buffers are neglected in this method!

        //initial function, set to zero
        K2DInit();

        //calculate slopes based on dem, and resets variables
        // if WH is not added to the DEM, this has to be done only once.
        double dt = _dt/2;  //1.0;
        double tof = 0.0;
        //maximum time is the lisem-timestep _dt
        while(tof < _dt-0.001)
        {

            K2DDEMA();
            //calculats water height, and computes the discharges according to manning etc.. and fluxes in 2 dimensions

            //function returns the minimal needed time-step for stable advection (dt > 1.0 for computational speed)
            dt = K2DFlux(dt);  //why _dt here???

            //only move in time that is left of the Lisem-timestep
            dt = std::min(dt, _dt-tof);

            if(SwitchKinematic2D == (int)K2D_METHOD_FLUX)
                K2DSolvebyFlux(dt);
            if(SwitchKinematic2D == (int)K2D_METHOD_INTER)
                K2DSolvebyInterpolation(dt);

            //solve fluxes and go back from water height to new discharge
            K2DSolve(dt);

            //total time this lisem-timestep
            tof += dt;
        }
    }

    if(SwitchKinematic2D == K1D_METHOD)
    {

        double mb = 0;
        double n = 0;

        // convert calculate Qn back to WH and volume for next loop
        FOR_ROW_COL_MV
        {

            /*   VJ 140105
    //                NEWTOWNPAHSON TO iterate h from Q. Because else we use alpha from before iteration
    //                Does not make a difference NOT NECESSARY but interesting code!
            double h, h1;
            double w = ChannelAdj->Drc;
            h = w > 0 ? (Alpha->Drc*pow(Qn->Drc, 0.6))/w : 0;//ChannelAdj->Drc;
            // first guess new h with old alpha
            h1 = h;
            if (Qn->Drc > 0)
            {
                double _23 = 2.0/3.0;
                double F, dF;
                int count = 0;

                do{
                    h = h1;
                    if (h < 1e-10)
                        break;
                    double P = w+2*h;
                    double A = h*w;
                    double R = A/P;

                    F = std::max((0.0, 1 - Qn->Drc/(sqrt(Grad->Drc)/N->Drc*A*powl(R,_23)));
                    dF = (5*w+6*h)/(3*h*P);
                    h1 = h - F/dF;
                    // function divided by derivative
                    count++;
                }while(fabs(h1-h) > 1e-10 && count < 20);
            }
          */

            WHrunoff->Drc = (Alpha->Drc*pow(Qn->Drc, 0.6))/ChannelAdj->Drc;


            //new WH based on A/dx = alpha Q^beta / dx

            double WaterVolout = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc;
            // new volume

            double InfilKWact = QinKW->Drc*_dt + WaterVolin->Drc - WaterVolout - Qn->Drc*_dt;
            //diff volume is sum of incoming fluxes+volume before - outgoing flux - volume after
            // this is the actual infiltration in the kin wave

            double diff = InfilKWact;
            InfilKWact = std::min(InfilKWact, -FSurplus->Drc*SoilWidthDX->Drc*DX->Drc);
            // infil volume cannot be more than surplus infil


            if (FFull->Drc == 1)
                InfilKWact = 0;
            //if profile full no more infil, surplus is 0

            difkin->Drc = 0;//(diff - InfilKWact);
            // difkin is not used, only to denug possible error in kin wave

            mb += (diff - InfilKWact);
            if (WHrunoff->Drc > 0)
                n+=1;

            if (SwitchBuffers && BufferVol->Drc > 0)
            {
                //qDebug() << "slope" << BufferVol->Drc << q->Drc*_dt << WaterVolin->Drc << WaterVolall->Drc << Qn->Drc*_dt << diff;
                //NOTE: buffervolume is affected by sedimentation, this causes a water volume loss that is corrected in the
                // totals and mass balance functions
            }
            else
                InfilVolKinWave->Drc = InfilKWact;

        }


        // mass balance correction, throw error on cells with WH
        //qDebug() << mb;

        if (n > 0)
            mb = mb/n;

        FOR_ROW_COL_MV
        {
            if (WHrunoff->Drc > 0)
                WHrunoff->Drc += mb/(ChannelAdj->Drc*DX->Drc);
        }
    }else
    {
        FOR_ROW_COL_MV
        {


            double err =  -WHrunoff->Drc * ChannelAdj->Drc * DX->Drc - QoutKW->Drc + QinKW->Drc +  WaterVolin->Drc - K2DI->Drc;
            //throw calculation error in infiltration, error should be insignificant
            InfilVolKinWave->Drc = K2DI->Drc;

        }
    }

    FOR_ROW_COL_MV
            // if (hmx->Drc == 0)
    {
        WHroad->Drc = WHrunoff->Drc;
        // set road to average outflowing wh, no surface storage.

        WH->Drc = WHrunoff->Drc + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        if (ChannelAdj->Drc > 0 && WHrunoff->Drc > 0)
            V->Drc = Qn->Drc/(WHrunoff->Drc*ChannelAdj->Drc);
        else
            V->Drc = 0;
        // recalc velocity for output to map, is not used in other processes

        WaterVolall->Drc = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        // is the same as :         WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);


    }

    if(SwitchKinematic2D == K1D_METHOD)
    {

        FOR_ROW_COL_MV
        {
            if (SwitchErosion)
            {
                Conc->Drc = (Qn->Drc > 1e-6 ? Qs->Drc/Qn->Drc : 0);
                //MaxConcentration(WaterVolall->Drc, Sed->Drc);
                // CHANGED, MORE STABLE CONC 19/9/13
                // correct for very high concentrations, 850 after Govers et al
                // recalc sediment volume

                if (SwitchPesticide)
                {
                    //C->Drc = ConcentrationP(WaterVolall->Drc, Pest->Drc);
                    C->Drc = Qn->Drc > 1e-10 ? Qpn->Drc/Qn->Drc : 0;
                    C_N->Drc = C->Drc;
                    //qDebug()<< "ds overlandflow"<< C->Drc;
                    //qDebug()<< "ds overlandflow"<< Pest->Drc;
                }

            }
        }
    }

}
//---------------------------------------------------------------------------
void TWorld::OverlandFlow(void)
{
}
