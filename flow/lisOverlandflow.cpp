/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
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
    if(SwitchKinematic2D == K2D_METHOD_DYN) {
        OverlandFlow2Ddyn();
        // dynamic wave overland flow, water and sediment and pesticides
    } else {
        // kin wave overland flow

        CalcVelDisch();
        // Q, V and Alpha Manning

        // calc all erosion
        if (SwitchErosion) {

           // cell_FlowDetachment(); // obsolete

            SedimentDetachmentSS(_dt, WHrunoff, ChannelAdj, V, Sed, Conc, TC, DETFlow, DEP, SettlingVelocitySS, SUSPrunoff);
            // same sed detachment and deposition as in 2D flow
            // full flowwidth is used, but adjusted inside for fractions for roads, houses etc

            if (SwitchPest) {
                PesticideFlowDetachmentSS(Sed);
            }
        }

        OverlandFlow1D();
        // routing: kinematic wave of water and sediment

        // move water and sed into channel
        if (SwitchIncludeChannel) {

            ToChannelBroadWeir();
            // kin wave interaction with channel (FloodDomain = 0)

            // TODO pesticide to channel

            // if 2D overflow do that
            if (SwitchKinematic2D == K2D_METHOD_KINDYN) {
                ToFlood();
                // transfer kin wave WHrunoff and sed to flood height hmx and SSFlood where both exist
                ChannelOverflowBroadWeir(hmxrunoff, V);
                // 2D flow part interact with channel (FloodDomain > 0)
                ChannelFlood();
                // dyn wave for flooded part
            }
        }

    }
}

//--------------------------------------------------------------------------------------------
void TWorld::OverlandFlow2Ddyn(void)
{
    double dtOF = 0;

    // NOTE: only broad crested weir works with different channel shapes!
    ChannelOverflowBroadWeir(WHrunoff, V);
    // this changes channel and surface water volume

    // obsolete this is only for a rectangular channel
    //    ChannelOverflow(WHrunoff, V);

    // after this new ChannelHW and WHrunoff, and Susp sediment values ChannelSSSed and SSFlood->Drc

    startFlood = false;
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (WHrunoff->Drc > 0)
            startFlood = true;
    }}

    if(startFlood) {
        dtOF = fullSWOF2openMUSCL(WHrunoff, Uflood, Vflood, DEM);
        TIMEDB(QString("Average dynamic timestep in flooded cells (dt %1 sec, n %2)").arg(dtOF,6,'f',3).arg(iter_n,4));
        // some screen reporting

        // calc discharge flux form the last flux in the loop
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            V->Drc = qSqrt(Uflood->Drc*Uflood->Drc + Vflood->Drc*Vflood->Drc);
            Qn->Drc = V->Drc*(WHrunoff->Drc*ChannelAdj->Drc);
        }}
    }

    updateWHandHmx();
    // update all water levels and volumes and calculate partition flood and runoff for output

    FloodMaxandTiming();

}
//--------------------------------------------------------------------------------------------
// ToChannel is ONLY called with KIN or KINDYN
/**
 * @fn void TWorld::ToChannel(void)
 * @brief Calculates fraction of overland flow that flows into channel
 *
 * Calculates fraction of overland flow that flows into channel.
 * This fraction is based on channel width and flow velocity
 *
 * @return void
 */

    //OBSOLETE: replaced wioth broadweir principles
void TWorld::ToChannel()
{
    if (!SwitchIncludeChannel)
        return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (WHrunoff->Drc > 0 && FloodDomain->Drc == 0 && !crch_.at(i_).culvert) {

            double fractiontochannel = qMin(1.0, _dt*V->Drc/(0.5*ChannelAdj->Drc));
            // fraction to channel calc from half the adjacent area width and flow velocity

            if (SwitchKinematic2D == K2D_METHOD_KINDYN &&
                    WHrunoff->Drc <= qMax(0.0 , ChannelWH->Drc - ChannelDepth->Drc))
                fractiontochannel = 0;
            // cannot flow into channel if water level in channel is higher than runoff depth

            if (fractiontochannel > 0) {
                double dwh = fractiontochannel*WHrunoff->Drc;
                double dvol = dwh*CHAdjDX->Drc;//fractiontochannel*(WaterVolall->Drc - MicroStoreVol->Drc);

                // water diverted to the channel
                ChannelWaterVol->Drc += dvol;
                ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);

                WHrunoff->Drc -= dwh;
                WH->Drc -= dwh;
                hmxWH->Drc = WH->Drc;
                WaterVolall->Drc = CHAdjDX->Drc*hmxWH->Drc;        //(WHrunoff->Drc) + MicroStoreVol->Drc;

                if (SwitchErosion) {
                    double dsed = fractiontochannel*Sed->Drc;
                    double maxsed = MAXCONC * ChannelWaterVol->Drc;
                    if (ChannelSSSed->Drc  + dsed > maxsed)
                        dsed = maxsed - ChannelSSSed->Drc;
                    if (dsed > 0) {
                        ChannelSSSed->Drc += dsed; //sediment diverted to the channel
                        Sed->Drc -= dsed;
                        Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
                        // adjust sediment in suspension
                        RiverSedimentLayerDepth(r,c);
                        RiverSedimentMaxC(r,c);
                    }
                }
            }
        }
   }}
}
//--------------------------------------------------------------------------------------------
// used only in kin wave without overflow
// based on broad crested weir function like 2D flow
void TWorld::ToChannelBroadWeir()
{
    if (!SwitchIncludeChannel)
        return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (WHrunoff->Drc > 0 && FloodDomain->Drc == 0 && !crch_.at(i_).culvert) {

            if (SwitchKinematic2D == K2D_METHOD_KINDYN &&
                    WHrunoff->Drc <= qMax(0.0 , ChannelWH->Drc - ChannelDepth->Drc))
                continue;
            // cannot flow into channel if water level in channel is higher than runoff depth

            // potentially all surface water flows into channel
            double Cd = 0.56;
            double lengthfactor = 2.0*_dt*ChannelDX->Drc;
            double velocityfactor = (V->Drc*V->Drc)/(2*GRAV);
            double freeflow_tochan = lengthfactor*Cd*SQRT2G*std::pow(WHrunoff->Drc+velocityfactor,1.5);
            //free flow broad crested weir, water flows over edge to deeper water in channel

            double volintochan = qMin(freeflow_tochan, CHAdjDX->Drc * WHrunoff->Drc);
            // in m3, minimum off what is there and broad crested flow

            WaterVolall->Drc -= volintochan;
            ChannelWaterVol->Drc += volintochan;
            WHrunoff->Drc = qMax(0.0,WaterVolall->Drc - MicroStoreVol->Drc)/CHAdjDX->Drc;
            WH->Drc = WaterVolall->Drc/CHAdjDX->Drc;
            hmxWH->Drc = WH->Drc + hmx->Drc;

            if (SwitchErosion) {
                double sed = volintochan * Conc->Drc;  //SSCFlood->Drc; //???????? why SSCFlood is only used when 2Dflow?
                Sed->Drc -= sed;
                Conc->Drc = MaxConcentration(WaterVolall->Drc, Sed->Drc);
                ChannelSSSed->Drc += sed;
                RiverSedimentLayerDepth(r,c);
                RiverSedimentMaxC(r, c);
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
void TWorld::CalcVelDisch()
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double Perim = SwitchPerimeterKW ? FlowWidth->Drc+2*WHrunoff->Drc : FlowWidth->Drc;
        double Area = FlowWidth->Drc*WHrunoff->Drc;
        V->Drc = pow(Area/Perim, (2.0/3.0)) * std::sqrt(Grad->Drc)/N->Drc; //WHrunoff->Drc
        Q->Drc = V->Drc * Area;//pow(Area/Alpha->Drc, (5.0/3.0)); // A = aplha*Q^beta => Q = (A/alpha)^1/beta and  beta = 6/10 = 3/5

        if (Grad->Drc > 1e-6)
            Alpha->Drc = pow(N->Drc/std::sqrt(Grad->Drc) * pow(Perim, 2.0/3.0),0.6);
        else
            Alpha->Drc = 0;
    }}
}
//---------------------------------------------------------------------------
void TWorld::updateWHandHmx(void)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double  hmxflood = 0;
        if (FloodDomain->Drc == 0) {
            WH->Drc = WHrunoff->Drc + WHstore->Drc;
            hmxWH->Drc = WH->Drc; // in 2D flow hmx is 0, not used
            WaterVolall->Drc = WH->Drc* CHAdjDX->Drc;// WHrunoff->Drc*CHAdjDX->Drc + MicroStoreVol->Drc;
            hmxflood = WHrunoff->Drc;
        } else {
            hmx->Drc = hmxrunoff->Drc + WHstore->Drc;
            hmxWH->Drc = WH->Drc + hmx->Drc; // in 2D flow hmx is 0, not used
            WaterVolall->Drc = hmxWH->Drc* CHAdjDX->Drc;// WHrunoff->Drc*CHAdjDX->Drc + MicroStoreVol->Drc;
            hmxflood = WHrunoff->Drc+hmxrunoff->Drc;
        }

        FloodWaterVol->Drc = qMax(0.0,hmxflood - minReportFloodHeight)*CHAdjDX->Drc;
        // used in mass balance
        RunoffWaterVol->Drc = qMin(hmxflood, minReportFloodHeight)*CHAdjDX->Drc;
        // all water that is not flood and not stored, so below min level

        if (SwitchErosion) {
            double sed = (SSFlood->Drc + BLFlood->Drc);
            Conc->Drc =  MaxConcentration(WaterVolall->Drc, sed);
            SSCFlood->Drc = MaxConcentration(WaterVolall->Drc, SSFlood->Drc);
            BLCFlood->Drc = MaxConcentration(WaterVolall->Drc, BLFlood->Drc);

            Qsn->Drc = Conc->Drc*Qn->Drc;
        }

    }}
}


//--------------------------------------------------------------------------------------------
void TWorld::OverlandFlow1D(void)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
       // tmd->Drc = DX->Drc * FlowWidth->Drc * WHrunoff->Drc;
        // temp volume on the move

        QinKW->Drc = 0; // store for incoming water in a cell
        tma->Drc = 0; // potentially available for limiting flow, does not have to be channel!
        //do not make -1!

        if (SwitchErosion) {
            // calc sediment flux going in kin wave as Qs = Q*C
            Qsn->Drc = 0.0;
            Conc->Drc = MaxConcentration(WHrunoff->Drc * CHAdjDX->Drc, Sed->Drc);
            Qs->Drc =  Q->Drc * Conc->Drc;
            // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
        }
    }}

    // route water
    KinematicExplicit(crlinkedldd_, Q, Qn, Alpha,DX, tma, tma);

    //convert calculated Qn back to WH and volume for next loop
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double WaterVolout = (CHAdjDX->Drc * WHrunoff->Drc) + _dt*(QinKW->Drc - Qn->Drc);
        // volume mass balance, WHrunoff is still the old one

        WHrunoff->Drc = WaterVolout/CHAdjDX->Drc;
        double Area = WaterVolout/ChannelAdj->Drc;
        V->Drc = Area > 1e-12 ? Qn->Drc/Area : 0.0;

        // Alpha->Drc = Qn->Drc > 0 ? Area/pow(Qn->Drc,0.6) : Alpha->Drc;
        // CAREFULL??? gave errors in channelalpha

        WH->Drc = WHrunoff->Drc + WHstore->Drc;

        hmxWH->Drc = WH->Drc + hmx->Drc;
        //needed for totals and output

        WaterVolall->Drc = WHrunoff->Drc*CHAdjDX->Drc + MicroStoreVol->Drc;

    }}

    if (SwitchErosion)
    {
        KinematicSubstance(crlinkedldd_,LDD, Q, Qn, Qs, Qsn, Alpha, DX, Sed, tma);

        FOR_ROW_COL_MV_L {
            if (Sed->Drc > MAXCONC * WaterVolall->Drc) {
                double ss = Sed->Drc;
                Sed->Drc = MAXCONC * WaterVolall->Drc;
                double ds = ss - Sed->Drc;
                DEP->Drc -= ds;
            }
        }}
    }

    if (SwitchPest) {
        //this function takes care of dissolved and sorbed kinematic wave
        PesticideFlow1D();
    }
}
