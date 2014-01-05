
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

#include "model.h"
#define tiny 1e-8

#define FW ChannelAdj->Drc

//---------------------------------------------------------------------------
//fraction of water and sediment flowing into the channel
void TWorld::QToChannel(void)
{
    if (SwitchIncludeChannel)
    {
        FOR_ROW_COL_MV_CH
        {
            double Volume = WHrunoff->Drc * FlowWidth->Drc * DX->Drc;

            if (Volume == 0)
            {
                SedToChannel->Drc = 0;
                RunoffVolinToChannel->Drc = 0;
                continue;
            }

            double R = (2*WHrunoff->Drc+FlowWidth->Drc)/(WHrunoff->Drc*FlowWidth->Drc);
            double vtoc = pow(R,2/3)*sqrt(WHrunoff->Drc/(0.5*ChannelAdj->Drc))/N->Drc;
            //assume slope perpendicular to flow is based on water level only
            double qtoc = min(Volume, 2 * vtoc * WHrunoff->Drc*FlowWidth->Drc * _dt);
            // 2 because flux comes from both sides, cannot be more than volume itself

            if (SwitchBuffers)
                if (BufferID->Drc > 0)
                    qtoc = Volume;
            // where there is a buffer in the channel, all goes in the channel

            if (SwitchChannelFlood)
            {
                if (ChannelWH->Drc >= ChannelDepth->Drc+ChannelLevee->Drc)
                    qtoc = 0;
                // no inflow when flooded
                if (ChannelMaxQ->Drc > 0)
                    qtoc = 0;
                // no surface inflow when culverts and bridges
            }

            if (SwitchAllinChannel && LDD->Drc == 5)
                qtoc = Volume;
            // in catchment outlet cell, throw everything in channel

            RunoffVolinToChannel->Drc = qtoc;
            // water diverted to the channel
            WHrunoff->Drc = (Volume-qtoc)/(FlowWidth->Drc * DX->Drc);

            WH->Drc = WHrunoff->Drc + WHstore->Drc;

            if (SwitchErosion)
            {
                SedToChannel->Drc = ChannelConc->Drc*qtoc;
                //sediment diverted to the channel in kg
                Sed->Drc -= SedToChannel->Drc;
                // adjust sediment in suspension
            }
        }

        CalcVelDisch(true);
        // recalc velocity and discharge with new HWrunoff
    }
}
//---------------------------------------------------------------------------
//fraction of water and sediment flowing into the channel
void TWorld::ToChannel(void)
{
    if (SwitchIncludeChannel)
    {
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
                fractiontochannel = min(1.0, _dt*V->Drc/(0.5*ChannelAdj->Drc));

//            if (WHrunoff->Drc == 0)
//                fractiontochannel = 0;
//            //VJ 140105 to prevent sed to channel if there is no water

            if (SwitchBuffers)
                if (BufferID->Drc > 0)
                    fractiontochannel = 1.0;
            // where there is a buffer in the channel, all goes in the channel

            if (SwitchChannelFlood)
            {
                if (ChannelWH->Drc >= ChannelDepth->Drc+ChannelLevee->Drc)
                    fractiontochannel = 0;
                // no inflow when flooded
                if (ChannelMaxQ->Drc > 0)
                    fractiontochannel = 0;
                // no surface inflow when culverts and bridges
            }
            if (SwitchAllinChannel)
                if (LDD->Drc == 5)    ///////////// VJ 5-1-14, for all outlets
                    fractiontochannel = 1.0;
            // in catchment outlet cell, throw everything in channel

            RunoffVolinToChannel->Drc = fractiontochannel*Volume;
            // water diverted to the channel
            WHrunoff->Drc *= (1-fractiontochannel);

            WH->Drc = WHrunoff->Drc + WHstore->Drc;
            //VJ 130425

            // adjust water height
            if (SwitchErosion)
            {
                SedToChannel->Drc = fractiontochannel*Sed->Drc;
                //sediment diverted to the channel
                Sed->Drc -= SedToChannel->Drc;
                // adjust sediment in suspension


            }
        }
        CalcVelDisch(true);
        // recalc velocity and discharge
    }

}
//---------------------------------------------------------------------------
void TWorld::CalcVelDisch(bool onlychannel)
{
    //  tm->fill(0);
    //   tma->fill(0);
    FOR_ROW_COL_MV
    {
        double Perim;
        const double beta = 0.6;
        const double _23 = 2.0/3.0;
        double beta1 = 1/beta;
        //double kinvisc = 1.1e-6; // 15 degrees celcius water
        double NN = N->Drc;

        if(onlychannel && ChannelWidthUpDX->Drc == 0)
            continue;

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

        Alpha->Drc = pow(NN/sqrtGrad->Drc * pow(Perim, _23),beta);

        if (Alpha->Drc > 0)
            Q->Drc = pow((FlowWidth->Drc*WHrunoff->Drc)/Alpha->Drc, beta1);
        else
            Q->Drc = 0;

        V->Drc = pow(R->Drc, _23)*sqrtGrad->Drc/NN;

        //tm->Drc = V->Drc * R->Drc/kinvisc;
        //Reynolds number
    }
    // tm->report("reyn");
    //  tma->report("reynF");

}
//---------------------------------------------------------------------------
void TWorld::OverlandFlow(void)
{

    /*---- Water ----*/
    // recalculate water vars after subtractions in "to channel"
    FOR_ROW_COL_MV
    {
        //FSurplus->Drc = 0;

        WaterVolin->Drc = DX->Drc * FlowWidth->Drc * WHrunoff->Drc;
        //volume runoff into the kin wave, needed to determine infil in kin wave

        // WaterVolin total water volume in m3 before kin wave, WHrunoff may be adjusted in tochannel
        q->Drc = FSurplus->Drc*SoilWidthDX->Drc/_dt;
        // infil flux in kin wave (<= 0)negative value), in m2/s, in kiv wave DX is used
        // surplus related to infiltrating surfaces
    }

    //NOTE if buffers: all water into channel

    /*---- Sediment ----*/
    if (SwitchErosion)
    {
        // calc seediment flux going in kin wave as Qs = Q*C
        FOR_ROW_COL_MV
        {
            Qs->Drc =  Q->Drc * Conc->Drc;
            // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
        }
    }

    Qn->setMV();
    Qsn->fill(0);
    QinKW->fill(0);
    // flag all new flux as missing value, needed in kin wave and replaced by new flux
    FOR_ROW_COL_MV
    {
        if (LDD->Drc == 5) // if outflow point, pit
        {
            /* TODO: WHEN MORE PITS QPEAK IS FIRST INSTEAD OF MAIN PIT? */
            Kinematic(r,c, LDD, Q, Qn, Qs, Qsn, q, Alpha, DX, WaterVolin, Sed, BufferVol, BufferSed);
            //VJ 110429 q contains additionally infiltrated water volume after kin wave in m3
        }
    }
    /*
      routing of substances add here!
      do after kin wave so that the new flux Qn out of a cell is known
      you need to have the ingoing substance flux QS (mass/s)
      and it will give outgoing flux QSn (mass/s)
      and the current amount Subs (mass) in suspension+solution
    */

    if (SwitchPesticide)
    {
        // calc pesticide flux going in kin wave as Qp = Q*C
        FOR_ROW_COL_MV
        {
            Qp->Drc =  Q->Drc * C->Drc;
            // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
        }

        Qpn->fill(0);
        FOR_ROW_COL_MV
        {
            if (LDD->Drc == 5) // if outflow point, pit
            {
                routeSubstance(r,c, LDD, Q, Qn, Qp, Qpn, Alpha, DX, WaterVolin, Pest);
            }
        }
    }

    double mb = 0;
    double mbarea = 0;

    // convert calculate Qn back to WH and volume for next loop
    FOR_ROW_COL_MV
    {

//VJ 140105
// NEWTOWNPAHSON TO iterate h from Q. Does not make a difference
// NOT NECESSARY but interesting code!
//        double h, h1;
//        double w = ChannelAdj->Drc;
//        h = w > 0 ? (Alpha->Drc*pow(Qn->Drc, 0.6))/w : 0;//ChannelAdj->Drc;
//        // first guess new h with old alpha
//        h1 = h;
//        if (Qn->Drc > 0)
//        {
//            double _23 = 2.0/3.0;
//            double F, dF;
//            int count = 0;

//            do{
//                h = h1;
//                if (h < 1e-10)
//                    break;
//                double P = w+2*h;
//                double A = h*w;
//                double R = A/P;

//                F = max(0, 1 - Qn->Drc/(sqrtGrad->Drc/N->Drc*A*powl(R,_23)));
//                dF = (5*w+6*h)/(3*h*P);
//                h1 = h - F/dF;
//                // function divided by derivative
//                count++;
//            }while(fabs(h1-h) > 1e-10 && count < 20);
//        }

        double WHoutavg = (Alpha->Drc*pow(Qn->Drc, 0.6))/ChannelAdj->Drc;
        //new WH based on A/dx = alpha Q^beta / dx

        double WaterVolout = WHoutavg*ChannelAdj->Drc*DX->Drc;
        // new volume

        WHroad->Drc = WHoutavg;
        // set road to average outflowing wh, no surface storage.

        WH->Drc = WHoutavg + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        if (ChannelAdj->Drc > 0 && WHoutavg > 0)
            V->Drc = Qn->Drc/(WHoutavg*ChannelAdj->Drc);
        else
            V->Drc = 0;
        // recalc velocity for output to map, is not used in other processes

 //       WaterVolall->Drc = DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);

        WaterVolall->Drc = WaterVolout + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        // new water volume after kin wave, all water incl depr storage

        double InfilKWact = QinKW->Drc*_dt + WaterVolin->Drc - WaterVolout - Qn->Drc*_dt;
        //diff volume is sum of incoming fluxes+volume before - outgoing flux - volume after
        // this is the actual infiltration in the kin wave

        double diff = InfilKWact;
        InfilKWact = min(InfilKWact, -FSurplus->Drc*SoilWidthDX->Drc*DX->Drc);
        // infil volume cannot be more than surplus infil

        if (FFull->Drc == 1)
            InfilKWact = 0;
        //if profile full no more infil, surplus is 0

        difkin->Drc = 0;
        // difkin is not used, only to denug possible error in kin wave

//        mb += (diff - InfilKWact);
//        if (WH->Drc > 0)
//            mbarea += ChannelAdj->Drc * DX->Drc;

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
//    mb = mbarea > 0 ? mb/mbarea : 0;
//    FOR_ROW_COL_MV
//    {
//        if (WH->Drc > 0)
//        {
//            WH->Drc = max(0, WH->Drc + mb);
//            WHroad->Drc = max(0, WHroad->Drc + mb);
//        }
//        WaterVolall->Drc = DX->Drc*(WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
//    }

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
                C->Drc = Qn->Drc > 1e-10 ? Qpn->Drc/Qn->Drc : 0;//ConcentrationP(WaterVolall->Drc, Pest->Drc);
                C_N->Drc = C->Drc;
                //qDebug()<< "ds overlandflow"<< C->Drc;
                //qDebug()<< "ds overlandflow"<< Pest->Drc;
            }

        }
    }
}
//---------------------------------------------------------------------------

// newton raphson trials for h out of new Q
// Alpha->Drc = pow(NN/sqrtGrad->Drc * pow(w+2h, _23),0.6);
// wh = Alpha->Drc*pow(Qn->Drc, 0.6))
// wh = pow(NN/sqrtGrad->Drc * pow(w+2h, _23),0.6)*pow(Qn->Drc, 0.6)
// h = pow(N/G * pow(w+2h, 2/3)*Q, 0.6)/w;
// h = pow( N/G * pow(w+2h,2/3)*wh*G/N*pow(wh/(w+2h),2/3), 0.6)/w;
// h = pow(wh*pow(w+2h,2/3)*pow(wh/(w+2h),2/3), 0.6)/w;
// h = pow(wh*pow(w+2h,2/3)*pow(wh,2/3)/pow(w+2h,2/3)    , 0.6)/w;
// h = pow(pow(wh,5/3),0.6)/w;
// F = h - pow(pow(wh,5/3),0.6)/w;
//                F = 1 - Q/(gn*(powl(w*h,_53)/powl(P, _23)));
//                F = 1 - wh* gn * powl(R, 2/3)/(gn*(powl(w*h,_53)/powl(P, _23)));
//                F = 1 - wh* (powl(A,2/3)/powl(P,2/3)) / (powl(A,5/3)/powl(P, 2/3));
//                F = 1 - A*powl(A,2/3)/powl(A,5/3);
//                F = 1 - powl(A,5/3)/powl(A,5/3);

