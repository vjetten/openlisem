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
**  website, information and code: https://github.com/vjetten/openlisem
**
*************************************************************************/
/*!
  \file lisOverlandflow.cpp
  \brief calculate interactions between channel flow, overland flow and flooding, calculate Q, V and call the kin wave
*/

#include <algorithm>
#include "model.h"
#include "operation.h"
#define tiny 1e-8

//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::OverlandFlow(void)
 * @brief Calls the kinematic wave or diffusive wave functions and calculates new discharge, water height and sediment presence
 *
 * Calls the kinematic, diffusive or dynamic wave functions and calculates new discharge, water height and sediment presence
 * During this process, surpluss potential infilration is subtracted from the water content.
 * Based on the options in the run file, either the 1D or 2D kinematic wave is used.
 * Sediment transport in overland flow is automatically taken into accaunt.
 */

//---------------------------------------------------------------------------
void TWorld::OverlandFlow(void)
{
    ToTiledrain();   // fraction going into tiledrain directly from surface, for 1D and 2D flow

    // kinematic wave or kin wave with overflow
    if(SwitchKinematic2D == K2D_METHOD_KIN || SwitchKinematic2D == K2D_METHOD_KINDYN) {

        CalcVelDisch();
        // overland flow velocity, discharge and alpha
        // V is needed in erosion

        if (SwitchErosion) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L  {
                cell_FlowDetachment(r, c);
                // kine wave based flow detachment
            }}
        }

        ToChannel();        // overland flow water and sed flux going into or out of channel, in channel cells
        OverlandFlow1D();   // kinematic wave of water and sediment

        if(SwitchKinematic2D == K2D_METHOD_KINDYN) {
            ChannelFlood();
            // st venant channel 2D flooding from channel, only for kyn wave
        }
    }

    if(SwitchKinematic2D == K2D_METHOD_DYN) {
        OverlandFlow2Ddyn();
        // dynamic wave overland flow
    }


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
void TWorld::ToChannel() //(int r, int c)
{
    if (!SwitchIncludeChannel)
        return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
    if (ChannelWidth->Drc > 0 && WHrunoff->Drc > HMIN && hmx->Drc < HMIN)
    {
        double fractiontochannel = std::min(1.0, _dt*V->Drc/(0.5*ChannelAdj->Drc));
        // fraction to channel calc from half the adjacent area width and flow velocity

        // cannot flow into channel if water level in channel is higher than runoff depth
        if (SwitchKinematic2D == K2D_METHOD_KINDYN &&
                WHrunoff->Drc <= std::max(0.0 , ChannelWH->Drc - ChannelDepth->Drc))
            fractiontochannel = 0;

        // no inflow on culverts
        if (SwitchCulverts && ChannelMaxQ->Drc  > 0)
            fractiontochannel = 0;

        if (fractiontochannel > 0) {

            double dwh = fractiontochannel*WHrunoff->Drc;
            double dvol = fractiontochannel*(WaterVolall->Drc - MicroStoreVol->Drc);

            // water diverted to the channel
            ChannelWaterVol->Drc += dvol;
            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);

            WHrunoff->Drc -= dwh ;
            WHroad->Drc -= dwh;
            WH->Drc -= dwh;
            WaterVolall->Drc = CHAdjDX->Drc*(WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;

            if (SwitchErosion)
            {
                double dsed = fractiontochannel*Sed->Drc;
                ChannelSSSed->Drc  += dsed;
                //sediment diverted to the channel
                Sed->Drc -= dsed;

                Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
                // adjust sediment in suspension
                RiverSedimentLayerDepth(r,c);
                RiverSedimentMaxC(r,c);
            }
        }
    }
   }}
}
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::CalcVelDisch()
 * @brief Calculates velocity and discharge based on water height for overland flow
 *
 * Calculates velocity and discharge and alpha based on water height for overland flow (WHrunoff)
 * Using the water height and energy gradient, mannings equation for flow velocity is used.
 * The manning's N is altered when flooding is present,
 * this slows down water while it is converted into flood water.
 *
 * @return void
 * @see mixing_coefficient
 */
void TWorld::CalcVelDisch()//(int r, int c)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {

        double NN = N->Drc;
        double alpha;
        double WHr = WHrunoff->Drc;
        double FW = FlowWidth->Drc;
//        double mixing_coefficient = 2.0;
//        if (SwitchKinematic2D == K2D_METHOD_KINDYN && SwitchIncludeChannel && hmx->Drc > 0.001)
//            NN = N->Drc * (2.0-qExp(-mixing_coefficient*hmx->Drc));
        // slow down water in flood zone, if hmx = 0 then factor = 1

        if (Grad->Drc > MIN_SLOPE)
            alpha = pow(NN/sqrtGrad->Drc * pow(FW, 2.0/3.0),0.6);
        // perimeter = FlowWidth
        else
            alpha = 0;

        if (alpha > 0)
            Q->Drc = pow((FW*WHr)/alpha, 5.0/3.0); // Q = (A/alpha)^1/beta and  beta = 6/10 = 3/5
        else
            Q->Drc = 0;
        //Q = (A/alpha)^5/3 => A^5/3 / alpha^5/3 =? aplha^5/3 = (N/sqrtS^3/5)^5/3 *((P^2/3)^3/5)^5/3 =
        //Q =  A^5/3 / [N/Sqrt * P^2/3] => A*A^2/3 / P^2/3 * sqrtS/n = A * R^2/3 sqrtS/N = AV

        V->Drc = pow(WHr, 2.0/3.0) * sqrtGrad->Drc/NN;
        // overlandflow, we do not use perimeter here but height
        // note: we can use tortuosity here: perimeter = R/(w*tortuosity) = hw/(w*tort) = h/tort
        // tortuosity can come from random roughness! use analysis from EU project
        Alpha->Drc = alpha;

    }}
}
//---------------------------------------------------------------------------
void TWorld::updateWHandHmx(void)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double WHR = WHrunoff->Drc;

        WHroad->Drc = WHR;
        // set road to average outflowing wh, no surface storage.

        WH->Drc = WHR + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water
        WaterVolall->Drc = WHR*CHAdjDX->Drc + MicroStoreVol->Drc;
        //LOGIC:
        // water layer in three parts: WHstore < (WHrunoff < minReportFloodHeight) < (hmx > minReportFloodHeight)
        // WH =  WHstore + WHrunoff

        hmxWH->Drc = WH->Drc;// + hmx->Drc; // in 2D flow hmx is 0, not used
        // hmxWH is used for reporting and calculation of velocity on screen. It combines aall waterheight irrespective of type of flow
        // it is also used in channeloverflow and in ponded evaporation

        hmxflood->Drc = std::max(0.0, WHR - minReportFloodHeight);
        // should be WH - minfloodheight?
        // is used for reporting all water aboove a user minimum, the rest is overland flow

        FloodWaterVol->Drc = hmxflood->Drc*CHAdjDX->Drc;
        // used in mass balance
        RunoffWaterVol->Drc = std::min(WHR, minReportFloodHeight)*CHAdjDX->Drc;
        // all water that isnot flood and not stored

        if (SwitchErosion) {
            double sed = (SSFlood->Drc + BLFlood->Drc);
            //Conc->Drc =  MaxConcentration(WHrunoff->Drc * CHAdjDX->Drc, sed);
            Conc->Drc =  MaxConcentration(WaterVolall->Drc, sed);
            SSCFlood->Drc = MaxConcentration(WaterVolall->Drc, sed);

            Qsn->Drc = Conc->Drc*Qn->Drc;
        }

    }}
}


/**
 * @fn void TWorld::OverlandFlow2Ddyn()
 * @brief Does 2D flow, calling SWOF functions, WHrunoff
 *
 * @return void
  */
void TWorld::OverlandFlow2Ddyn(void)
{
    double dtOF = 0;

    if (SwitchChannel2DflowConnect)
        ChannelOverflowIteration(WHrunoff, V);
    else
        ChannelOverflow(WHrunoff, V);
    // Mixing of 2D runoff with channel water, V is used to determine how much flows into the channel
    // after this new ChannelHW and WHrunoff, and Susp sediment values ChannelSSSed and SSFlood->Drc

    if (SwitchIncludeChannel) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            ChannelWaterVol->Drc = ChannelWH->Drc * ChannelDX->Drc * ChannelWidth->Drc;
            WaterVolall->Drc = CHAdjDX->Drc*(WHrunoff->Drc + hmx->Drc) + MicroStoreVol->Drc;
            // do not recalc floodvol, MB errors

            // recalc channel water vol else big MB error
            if(SwitchErosion)
            {
                SWOFSedimentLayerDepth(r,c,WHrunoff->Drc, V->Drc);
                SWOFSedimentSetConcentration(r,c, WHrunoff);

                RiverSedimentLayerDepth(r,c);
                RiverSedimentMaxC(r, c);
                // all concentrations, possible ChannelDep when surplus
            }
        }}
    }

    startFlood = false;
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV {
        if (WHrunoff->Drc > HMIN)
            startFlood = true;
    }

    if (SwitchSWOFopen)
        dtOF = fullSWOF2open(WHrunoff, Uflood, Vflood, DEM);
    else
        dtOF = fullSWOF2RO(WHrunoff, Uflood, Vflood, DEM);
    //VJ new average flux over lisem timestep, else last Qn is used

    //  infilInWave(WHrunoff, _dt);

    // calc discharge flux
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        V->Drc = sqrt(Uflood->Drc*Uflood->Drc + Vflood->Drc*Vflood->Drc);
        Qn->Drc = V->Drc*(WHrunoff->Drc*ChannelAdj->Drc);
        //Q->Drc = Qn->Drc; // just to be sure
    }}

    Boundary2Ddyn();  // do the domain boundaries for Q, h and sediment

    updateWHandHmx();
    // update all water levels and volumes and calculate partition flood and runoff for output

    FloodMaxandTiming();

    TIMEDB(QString("Average dynamic timestep in flooded cells (dt %1 sec, n %2)").arg(dtOF,6,'f',3).arg(iter_n,4));
    // some screen error reporting
}



//--------------------------------------------------------------------------------------------
void TWorld::OverlandFlow1D(void)
{
    // recalculate water vars after subtractions in "to channel"

    //double tot = 0;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        WaterVolin->Drc = DX->Drc * FlowWidth->Drc * WHrunoff->Drc;
        //volume runoff into the kin wave, needed to determine infil in kin wave
        // WaterVolin total water volume in m3 before kin wave, WHrunoff may be adjusted in tochannel
        q->Drc = 0;//FSurplus->Drc*SoilWidthDX->Drc/_dt;
        // OBSOLETE? has never work properly
        // infil flux in kin wave (<= 0)negative value), in m2/s, in kiv wave DX is used
        // surplus related to infiltrating surfaces

        QinKW->Drc = 0; // store for incoming water in a cell
        //tot = tot + WaterVolin->Drc;

        if (SwitchErosion) {
            // calc seediment flux going in kin wave as Qs = Q*C
            Qsn->Drc = 0.0;
            Conc->Drc = MaxConcentration(WHrunoff->Drc * CHAdjDX->Drc, Sed->Drc);
            Qs->Drc =  Q->Drc * Conc->Drc;
            // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
        }
    }}

    // route water
    if (SwitchLinkedList) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            pcr::setMV(Qn->Drc);
            QinKW->Drc = 0;
        }}

        FOR_ROW_COL_LDD5 {
            Kinematic(r,c, LDD, Q, Qn,  Alpha, DX);
            // tm is not used in overland flow, in channel flow it is the max flux of e.g. culverts
        }}
    } else {
        KinematicExplicit(crlinkedldd_, Q, Qn, Alpha,DX);
    }

    //convert calculate Qn back to WH and volume for next loop
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double InfilKWact = 0;
        double WaterVolout = std::max(0.0, QinKW->Drc*_dt + WaterVolin->Drc  - Qn->Drc*_dt);
        // mass balance, this includes now errors!

        // new water vol is mass bal diff
        WHrunoff->Drc = WaterVolout/CHAdjDX->Drc;
        // runoff based on water vol out
        // NOTE route substance is already an explicit solution                      

        double diff = QinKW->Drc*_dt + WaterVolin->Drc - WaterVolout - Qn->Drc * _dt;
        InfilKWact = diff;//std::min(-FSurplus->Drc*SoilWidthDX->Drc*DX->Drc, diff);

        InfilVolKinWave->Drc = InfilKWact;

        Alpha->Drc = Qn->Drc > 0 ? (WHrunoff->Drc*FlowWidth->Drc)/pow(Qn->Drc,0.6) : Alpha->Drc;
        // needed for erosion // A = alpha Q^0.6 => alpha = A/Q^0.6
        V->Drc = pow(WHrunoff->Drc, 2.0/3.0) * sqrtGrad->Drc/N->Drc;
        // new velocity

        WHroad->Drc = WHrunoff->Drc;
        // set road to average outflowing wh, no surface storage.

        WH->Drc = WHrunoff->Drc + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        hmxWH->Drc = WH->Drc + hmx->Drc;//???? hmx here? why not
        //needed for totals and output

        WaterVolall->Drc = WHrunoff->Drc*CHAdjDX->Drc + MicroStoreVol->Drc;

        if (SwitchErosion)
             Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
    }}

    //      routing of substances add here!
    if (SwitchErosion)
    {
        if (SwitchLinkedList) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                pcr::setMV(Qsn->Drc);//Qsn->setAllMV();
            }}
            FOR_ROW_COL_LDD5 {
                routeSubstance(r,c, LDD, Q, Qn, Qs, Qsn, Alpha, DX, Sed);
            }}
        } else {
            KinematicSubstance(crlinkedldd_,LDD, Q, Qn, Qs, Qsn, Alpha, DX, Sed);
        }
    }

    // route other stuff
    if (SwitchPesticide) {
        // calc pesticide flux going in kin wave as Qp = Q*C
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
           Qp->Drc =  Qn->Drc * C->Drc;
        }}

        KinematicSubstance(crlinkedldd_, LDD, Q, Qn, Qp, Qpn, Alpha, DX, Pest);

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            //C->Drc = ConcentrationP(WaterVolall->Drc, Pest->Drc);
            C->Drc = Qn->Drc > MIN_FLUX ? Qpn->Drc/Qn->Drc : 0;
            C_N->Drc = C->Drc;
            //qDebug()<< "ds overlandflow"<< C->Drc;
            //qDebug()<< "ds overlandflow"<< Pest->Drc;
        }}
    }   
}
