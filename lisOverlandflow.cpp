
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
  \brief calculate interactions between channel flow, overland flow and flooding, calculate Q, V and call the kin wave

functions: \n
- void TWorld::RainfallToFlood(void)\n
- void TWorld::ToFlood(void)\n
- void TWorld::ToChannel(void)\n
- void TWorld::CalcVelDisch(void)\n
- void TWorld::OverlandFlowNew(void)\n
 */

#include <algorithm>
#include "model.h"
#include "operation.h"
#define tiny 1e-8

//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::RainfallToFlood(void)
 * @brief Calculates overland flow that becomes flooding
 *
 * Based on the overland flow energy gradient and the water height
 * overland flow can initiate flooding when TWorld::SwitchRainfallFlood is true.
 *
 * @return void
 * @see SwitchRainfallFlood
 * @see rainFloodingGradient
 */
// OBSOLETE with runoff 2D!
void TWorld::RainfallToFlood(void)
{
    /*
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
*/

}
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::ToFlood(void)
 * @brief Calculates overland flow that flows into flooding water
 *
 * Calculates overland flow of water and sediment that flows into flooding water
 * based on the runoff partitioning factor. Depending on the parameter, water
 * is either transformed quickly or slowly. This imitates the effect that overland
 * flow would have on the velocity of the flood water.
 *
 * @return void
 * @see runoff_partitioning
 */
void TWorld::ToFlood(void)
{
   /* if (!SwitchChannelFlood)
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

        if(SwitchErosion)
        {
            //better distribute this by ratio suspended Tc and
            SSFlood->Drc += Sed->Drc * frac;
            Sed->Drc = Sed->Drc * (1-frac);

            if(SwitchUseGrainSizeDistribution)
            {
                FOR_GRAIN_CLASSES
                {
                    SS_D.Drcd +=  Sed_D.Drcd * frac;
                    Sed_D.Drcd = Sed_D.Drcd * (1-frac);

                }
            }


            //immediately check for maximum concentration
            //if not done, too high concentration will show on display, before being deposited
            SWOFSedimentLayerDepth(r,c,hmx,Uflood,Vflood);

            SWOFSedimentMaxC(r,c);//,hmx,Uflood,Vflood);
        }



    }}*/
}
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::ToChannel(void)
 * @brief Calculates fraction of overland flow that flows into channel
 *
 * Calculates fraction of overland flow that flows into channel.
 * This fraction is based on channel width and flow velocity
 *
 * @return void
 */
void TWorld::ToChannel(void)
{
    if (!SwitchIncludeChannel)
        return;

    for (int  r = 0; r < _nrRows; r++)
    {
        for (int  c = 0; c < _nrCols; c++)
        {
            if(ChannelMaskExtended->data[r][c] == 1)
            {
                int rr = (int)ChannelSourceYExtended->Drc;
                int cr = (int)ChannelSourceXExtended->Drc;

                double fractiontochannel;
                double Volume = WHrunoff->Drc * FlowWidth->Drc * DX->Drc;

                if (Volume == 0)
                {
                    SedToChannel->Drcr = 0;
                    RunoffVolinToChannel->Drcr  = 0;
                    continue;
                }

                if (ChannelAdj->Drc == 0)
                    fractiontochannel = 1.0;
                else
                    fractiontochannel = std::min(1.0, _dt*V->Drc/std::max(0.01*_dx,0.5*ChannelAdj->Drc));
                // fraction to channel calc from half the adjacent area width and flow velocity

                if (SwitchBuffers)
                    if (BufferID->Drcr  > 0)
                        fractiontochannel = 1.0;
                // where there is a buffer in the channel, all goes in the channel

                // cannot flow into channel is water level in channel is higher than depth
                if (SwitchChannelFlood)
                {
                    if (WHrunoff->Drc <= std::max(ChannelLevee->Drcr, ChannelWH->Drcr -ChannelDepthExtended->Drc))
                    {
                        fractiontochannel = 0;
                    }
                    // no inflow when flooded
                    if (ChannelMaxQ->Drcr  > 0)
                    {
                        fractiontochannel = 0;
                    }
                    // no surface inflow when culverts and bridges
                }
                if (SwitchAllinChannel)
                    if (LDD->Drcr  == 5)
                        fractiontochannel = 1.0;
                // in catchment outlet cell, throw everything in channel

                RunoffVolinToChannel->Drcr  = fractiontochannel*Volume;
                // water diverted to the channel
                WHrunoff->Drc *= (1-fractiontochannel);

                WH->Drc = WHrunoff->Drc + WHstore->Drc;
                //VJ 130425

                if (SwitchErosion)
                {
                    if(!SwitchUse2Layer)
                    {
                        ChannelBLSed->Drcr  += fractiontochannel*Sed->Drc;
                    }else
                    {
                        ChannelSSSed->Drcr  += fractiontochannel*Sed->Drc;
                    }
                    //sediment diverted to the channel
                    Sed->Drc = Sed->Drc * (1 - fractiontochannel);

                    Conc->Drc = MaxConcentration(WHrunoff->Drc * DX->Drc * ChannelAdj->Drc, Sed->Drc);
                    // adjust sediment in suspension

                    if(SwitchUseGrainSizeDistribution)
                    {
                        Conc->Drc = 0;
                        FOR_GRAIN_CLASSES
                        {
                            if(SwitchUse2Layer)
                            {
                                RSS_D.Drcdr += fractiontochannel * Sed_D.Drcdr;
                            }else
                            {
                                RBL_D.Drcdr += fractiontochannel * Sed_D.Drcdr;
                            }
                            Sed_D.Drcd = Sed_D.Drcd * (1-fractiontochannel);
                            Conc_D.Drcd = MaxConcentration(WHrunoff->Drc * DX->Drc * ChannelAdj->Drc, Sed_D.Drcd);
                            Conc->Drc +=Conc_D.Drcd;
                        }
                    }

                   RiverSedimentLayerDepth(rr,cr);




                }
            }
        }
    }

    /*
     if (!SwitchIncludeChannel)
        return;

    FOR_ROW_COL_MV_CH
    {
        double fractiontochannel;
        double Volume = WHrunoff->Drc * FlowWidth->Drc * DX->Drc;
        double v = V->Drc;

        //K2D... maps are not initialized when 1D runoff is chosen
        //every 2D runoff specific code section must be performed within an if statement
        if(SwitchKinematic2D != K1D_METHOD)
        {
            //in case of a local depression, velocity was set to 0, leading to 0 channel inflow.
            //Once the depression filled above capacity, all water would in one timestep flow into the channel
            //this gave oscillating discharge.

            //now set the velocity temporarily according to manning's equation with the hydraulic slope as slope.
            if(K2DWHStore->Drc >0)
            {
                double hrunoff = std::max(WHrunoff->Drc,0.0);

                double Perim;
                const double beta = 0.6;
                const double _23 = 2.0/3.0;
                double beta1 = 1/beta;
                double NN = N->Drc;
                double R = 0;
                Perim = 2.0*hrunoff+FlowWidth->Drc;
                if (Perim > 0)
                    R = hrunoff*FlowWidth->Drc/Perim;
                else
                    R = 0;
                double Slope = hrunoff/_dx;
                v = pow(R, _23)*sqrt(Slope)/NN;
            }
        }

        if (Volume == 0)
        {
            SedToChannel->Drc = 0;
            RunoffVolinToChannel->Drc = 0;
            continue;
        }

        if (ChannelAdj->Drc == 0)
            fractiontochannel = 1.0;
        else
            fractiontochannel = std::min(1.0, _dt*v/std::max(0.01*_dx,0.5*ChannelAdj->Drc));
        // fraction to channel calc from half the adjacent area width and flow velocity

        if (SwitchBuffers)
            if (BufferID->Drc > 0)
                fractiontochannel = 1.0;
        // where there is a buffer in the channel, all goes in the channel

        // cannot flow into channel is water level in channel is higher than depth
        if (SwitchChannelFlood)
        {
            if (WHrunoff->Drc <= std::max(ChannelLevee->Drc, ChannelWH->Drc-ChannelDepth->Drc))
            {
                fractiontochannel = 0;
            }
            // no inflow when flooded
            if (ChannelMaxQ->Drc > 0)
            {
                fractiontochannel = 0;
            }
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
            if(!SwitchUse2Layer)
            {
                ChannelBLSed->Drc += fractiontochannel*Sed->Drc;
            }else
            {
                ChannelSSSed->Drc += fractiontochannel*Sed->Drc;
            }
            //sediment diverted to the channel
            Sed->Drc = Sed->Drc * (1 - fractiontochannel);

            Conc->Drc = MaxConcentration(WHrunoff->Drc * DX->Drc * ChannelAdj->Drc, Sed->Drc);
            // adjust sediment in suspension

            if(SwitchUseGrainSizeDistribution)
            {
                Conc->Drc = 0;
                FOR_GRAIN_CLASSES
                {
                    if(SwitchUse2Layer)
                    {
                        RSS_D.Drcd += fractiontochannel * Sed_D.Drcd;
                    }else
                    {
                        RBL_D.Drcd += fractiontochannel * Sed_D.Drcd;
                    }
                    Sed_D.Drcd = Sed_D.Drcd * (1-fractiontochannel);
                    Conc_D.Drcd = MaxConcentration(WHrunoff->Drc * DX->Drc * ChannelAdj->Drc, Sed_D.Drcd);
                    Conc->Drc +=Conc_D.Drcd;
                }
            }

            RiverSedimentLayerDepth(r,c);




        }
    }*/
}
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::CalcVelDisch()
 * @brief Calculates velocity and discharge based on water height for overland flow
 *
 * Calculates velocity and discharge based on water height for overland flow
 * Using the water height and energy gradient, mannings equation for flow velocity is used.
 * The manning's N is altered when flooding is present,
 * this slows down water while it is converted into flood water.
 *
 * @return void
 * @see mixing_coefficient
 */
void TWorld::CalcVelDisch()
{
    if(SwitchKinematic2D != 1)
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
        Perim = /*2*WHrunoff->Drc+*/ FlowWidth->Drc;

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
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::OverlandFlowNew(void)
 * @brief Calls the kinematic wave functions and calculates new discharge, water height and sediment presence
 *
 * Calls the kinematic wave functions and calculates new discharge, water height and sediment presence
 * During this process, surpluss potential infilration is subtracted from the water content.
 * Based on the options in the run file, either the 1D or 2D kinematic wave is used.
 * Sediment transport in overland flow is automatically taken into accaunt.
 *
 * @return void
 * @see SwitchKinematic2D
 * @see K1D_METHOD
 * @see K2D_METHOD_INTER
 */
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
    fill(*QinKW, 0.0);
    fill(*QoutKW, 0.0);
    // flag all new flux as missing value, needed in kin wave and replaced by new flux

    if(SwitchKinematic2D == K1D_METHOD)
    {

        if (SwitchErosion)
        {
            fill(*Qs, 0.0);
            fill(*Qsn, 0.0);
            // calc seediment flux going in kin wave as Qs = Q*C
            FOR_ROW_COL_MV
            {
                Conc->Drc = MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, Sed->Drc);
                Qs->Drc =  Q->Drc * Conc->Drc;
                // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
            }
        }

        Qn->setAllMV();
        FOR_ROW_COL_MV
        {
            if (LDD->Drc == 5) // if outflow point, pit
            {
                Kinematic(r,c, LDD, Q, Qn, q, Alpha, DX, WaterVolin, BufferVol);
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

        if (SwitchErosion)
        {
            if(!SwitchUseGrainSizeDistribution)
            {

                Qsn->setAllMV();
                FOR_ROW_COL_MV
                {
                    if (LDD->Drc == 5) // if outflow point, pit
                    {
                        routeSubstance(r,c, LDD, Q, Qn, Qs, Qsn, Alpha, DX, WaterVolin, Sed, BufferVol, BufferSed);
                    }
                }
            }else
            {
                FOR_GRAIN_CLASSES
                {
                    // calc seediment flux going in kin wave as Qs = Q*C
                    FOR_ROW_COL_MV
                    {
                        Conc_D.Drcd = MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, Sed_D.Drcd);
                        Tempa_D.Drcd =  Q->Drc * Conc_D.Drcd;
                        Qs->Drc +=  Q->Drc * Conc_D.Drcd;
                        // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
                    }
                    Tempb_D.at(d)->setAllMV();
                    FOR_ROW_COL_MV
                    {
                        if (LDD->Drc == 5) // if outflow point, pit
                        {
                            routeSubstance(r,c, LDD, Q, Qn, Tempa_D.at(d), Tempb_D.at(d), Alpha, DX, WaterVolin, Sed_D.at(d), BufferVol, BufferSed);
                        }
                    }

                    FOR_ROW_COL_MV
                    {
                        Qsn->Drc += Tempb_D.Drcd;
                        // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
                    }


                }

                fill(*Sed, 0.0);
                FOR_GRAIN_CLASSES
                {
                    FOR_ROW_COL_MV
                    {
                        Sed->Drc += Sed_D.Drcd;
                    }
                }

            }

        }

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
                    routeSubstance(r,c, LDD, Q, Qn, Qp, Qpn, Alpha, DX, WaterVolin, Pest, BufferVol, NULL);
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

        fill(*tmb, 0.0);

        //initial function, set to zero
        K2DInit();

        //calculate slopes based on dem, and resets variables
        //K2DDEMA()
        // if WH is not added to the DEM, this has to be done only once.

        double dt = _dt/2;
        double tof = 0.0;
        //maximum time is the lisem-timestep _dt
        while(tof < _dt-0.001)
        {

            K2DDEMA();
            //calculats water height, and computes the discharges according to manning etc.. and fluxes in 2 dimensions

            //function returns the minimal needed time-step for stable advection (dt > 1.0 for computational speed)
            dt = K2DFlux();

            //only move in time that is left of the Lisem-timestep
            dt = std::min(dt, _dt-tof);

            K2DSolvebyInterpolation(dt);

            /*sediment transport functions must be called before K2DSolve() and after K2DSolveBy..() */
            if(SwitchErosion)
            {
                //K2DQSOut is the boundary outflow that is returned by he K2DSolveBy....Sed() function.

                //advect total sediment
                if(!SwitchUseGrainSizeDistribution)
                {
                    K2DQSOut += K2DSolvebyInterpolationSed(dt,Sed, Conc);
                }else
                {
                    //advect each induvidual grain class
                    FOR_GRAIN_CLASSES
                    {
                        K2DQSOut += K2DSolvebyInterpolationSed(dt,Sed_D.at(d), Conc_D.at(d));
                    }
                }
            }

            //no longer needed
            K2DSolve(dt);

            //total time this lisem-timestep
            tof += dt;

            // add total volume outflow for this lisem timestep for average Qn
            FOR_ROW_COL_MV
            {
                tmb->Drc += K2DQ->Drc * dt;
            }
        }

        //VJ new average flux over lisem timestep
        FOR_ROW_COL_MV
        {
            K2DQ->Drc = tmb->Drc/_dt;
            Qn->Drc = tmb->Drc/_dt;
            Q->Drc = tmb->Drc/_dt;
        }

        if(SwitchErosion)
        {
            //calculate concentration and new sediment discharge
            if(!SwitchUseGrainSizeDistribution)
            {
                FOR_ROW_COL_MV
                {
                    Conc->Drc =  MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, Sed->Drc);
                    Qsn->Drc = Conc->Drc * Qn->Drc;
                    //Qs->Drc = Qsn->Drc;
                }
            }
            else
            {
                //calculate total sediment from induvidual grain classes,
                //and calculate concentration and new sediment discharge
                FOR_ROW_COL_MV
                {
                    Sed->Drc = 0;
                    Conc->Drc = 0;

                }
                FOR_ROW_COL_MV
                {
                    FOR_GRAIN_CLASSES
                    {
                        Sed->Drc += Sed_D.Drcd;
                        Conc_D.Drcd = MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, Sed_D.Drcd);
                        Conc->Drc += Conc_D.Drcd;
                    }
                }
            }
        }
    }

    if(SwitchKinematic2D == K1D_METHOD)
    {
        // convert calculate Qn back to WH and volume for next loop
        FOR_ROW_COL_MV
        {

            bool K1Dexplicit = true;
            double WaterVolout = 0;
            double InfilKWact = 0;
            if( K1Dexplicit)
            {
                WaterVolout = std::max(0.0, QinKW->Drc*_dt + WaterVolin->Drc  - Qn->Drc*_dt);
                // new water vol is mass bal diff
                WHrunoff->Drc = WaterVolout/(ChannelAdj->Drc*DX->Drc);
                // runoff based on water vol out
            }
            else
            {
                WHrunoff->Drc = (Alpha->Drc*pow(Qn->Drc, 0.6))/ChannelAdj->Drc;
                //new WH based on A/dx = alpha Q^beta / dx

                WaterVolout = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc;
                // new volume
            }

            InfilKWact = -FSurplus->Drc*SoilWidthDX->Drc*DX->Drc;
            InfilKWact = std::min(QinKW->Drc*_dt + WaterVolin->Drc - WaterVolout - Qn->Drc * _dt, InfilKWact);
            // infiltration is the surplus infil (pot infil), or infil is all that was there
            if (FFull->Drc == 1)
                InfilKWact = 0;
            //if profile full no more infil, surplus is 0

            if (SwitchBuffers && BufferVol->Drc > 0)
            {
                //qDebug() << "slope" << BufferVol->Drc << q->Drc*_dt << WaterVolin->Drc << WaterVolall->Drc << Qn->Drc*_dt << diff;
                //NOTE: buffervolume is affected by sedimentation, this causes a water volume loss that is corrected in the
                // totals and mass balance functions
            }
            else
                InfilVolKinWave->Drc = InfilKWact;

        }

    }
    else
    {
        FOR_ROW_COL_MV
        {
            //double err =  -WHrunoff->Drc * ChannelAdj->Drc * DX->Drc - QoutKW->Drc + QinKW->Drc +  WaterVolin->Drc - K2DI->Drc;
            //throw calculation error in infiltration, error should be insignificant
            InfilVolKinWave->Drc = K2DI->Drc;
        }
    }

    FOR_ROW_COL_MV
    {
        WHroad->Drc = WHrunoff->Drc;
        // set road to average outflowing wh, no surface storage.

        WH->Drc = WHrunoff->Drc + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        if(SwitchKinematic2D == K1D_METHOD)//if(SwitchKinematic2D != K1D_METHOD)
        {
            if(K2DSlope->Drc != 0 && K2DPits->Drc != 1)
            {
                if (ChannelAdj->Drc > 0 && WHrunoff->Drc > 0)
                    V->Drc = Qn->Drc/(WHrunoff->Drc*ChannelAdj->Drc);
                else
                    V->Drc = 0;
            }
            else
            {
                V->Drc = 0;
            }
        }else
        {
//<<<<<<< HEAD
//            if(Grad->Drc != 0)
//=======
            if(K2DSlope->Drc != 0 && K2DPits->Drc != 1)
            {
                if (ChannelAdj->Drc > 0 && WHrunoff->Drc > 0)
                    V->Drc = Qn->Drc/(WHrunoff->Drc*ChannelAdj->Drc);
                else
                    V->Drc = 0;
            }
            else
            {
                V->Drc = 0;
            }
        }

        WaterVolall->Drc = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        // is the same as :         WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);

    }

    if(SwitchKinematic2D == K1D_METHOD)
    {
        if (SwitchErosion)
        {
            FOR_ROW_COL_MV
            {

                //Conc->Drc = (Qn->Drc > 1e-6 ? Qs->Drc/Qn->Drc : 0);
                Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);

                // CHANGED, MORE STABLE CONC 19/9/13
                // correct for very high concentrations, 850 after Govers et al
                // recalc sediment volume

                if(SwitchUseGrainSizeDistribution)
                {
                    FOR_GRAIN_CLASSES
                    {
                        Conc_D.Drcd = (Qn->Drc > 1e-6 ? Tempb_D.Drcd/Qn->Drc : 0);
                    }
                }

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
