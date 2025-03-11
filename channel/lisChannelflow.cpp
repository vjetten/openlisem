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
 \file lisChannelflow.cpp
 \brief Channel hydrology and sediment detachment and movement processes.

functions: \n
- void TWorld::CalcVelDischChannel(void) calculate Velocity, alpha and Q in the channel \n
- void TWorld::ChannelFlow(void) calculate channelflow, ChannelDepth, do kinematic wave \n
*/

//#include <algorithm>
#include "model.h"
//#include "operation.h"



//---------------------------------------------------------------------------
void TWorld::ChannelFlowandErosion()
{
    if (!SwitchIncludeChannel)
        return;

    SwitchChannelKinWave = true;    // set to false for experimental swof in channel

    ChannelRainandInfil();          // subtract infil, add rainfall    

    ChannelBaseflow();              // add stationary and GW baseflow if selected

    ChannelVelocityandDischarge();  // maaings V Q Aplha

    ChannelFlow();                  // channel kin wave for water

    ChannelFlowDetachmentNew();     // detachment, deposition for SS and BL

    ChannelSedimentFlow();          // kin wave for sediment and substances

}
//---------------------------------------------------------------------------
void TWorld::ChannelVelocityandDischarge()
{
    // velocity, alpha, Q
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_CHL {

        double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
        double FWO = ChannelWidthO->Drc;
        //ChannelWH->Drc = Area/FWO;
        double CHWH = Area/FWO;
        double Perim = FWO+2*CHWH;
        double Radius = (Perim > 0 ? Area/Perim : 0);
        double sqrtgrad = std::max(sqrt(ChannelGrad->Drc), 0.001);
        double N = ChannelN->Drc;

        ChannelV->Drc = std::min(_CHMaxV,std::pow(Radius, 2.0/3.0)*sqrtgrad/N);
        ChannelQ->Drc = ChannelV->Drc * Area;
        ChannelAlpha->Drc = Area/std::pow(ChannelQ->Drc, 0.6);
        //ChannelAlpha->Drc = pow(N/sqrtgrad * pow(Perim, 2.0/3.0),0.6);  // no difference
    }}
}

//---------------------------------------------------------------------------
void TWorld::ChannelBaseflow(void)
{
    // add a stationary part
    if(SwitchChannelBaseflowStationary)
    {
        // first time
        if(!addedbaseflow) {
           #pragma omp parallel for num_threads(userCores)
           FOR_ROW_COL_MV_CHL {
                ChannelWaterVol->Drc += BaseFlowInitialVolume->Drc;
           }}

           addedbaseflow = true;
        }

        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelWaterVol->Drc += BaseFlowInflow->Drc * _dt;
        }}
    }

    // add the baseflow from GW
    if (SwitchGWflow) {

        GroundwaterFlow();
        // move groundwater, GWout is the flow itself between cells

        cTMap *pore = Thetaeff;
        cTMap *ksat = Ksateff;
        cTMap *SD = SoilDepth1init;
        if (SwitchTwoLayer) {
            pore = ThetaS2;
            ksat = Ksat2;
            SD = SoilDepth2init;
        }

        // in all channel cells
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            if (SwitchSWATGWflow) {
                Qbase->Drc = ChannelWidth->Drc/_dx * GWout->Drc;
            } else {
                double bedrock = DEM->Drc - SD->Drc;
                double chanbot = DEM->Drc - ChannelDepth->Drc;
                bedrock=chanbot;
                double dH = bedrock + GWWH->Drc - chanbot;
                if (dH > 0 && GWWH->Drc > 0) {
                   //Qbase->Drc = std::min(GWVol->Drc, 2.0 * dH/GWWH->Drc * GWout->Drc);
                //   Qbase->Drc = std::min(GWVol->Drc, 2.0 * fabs(GWout->Drc));
                   Qbase->Drc = 2*GWout->Drc;
                   // use the fraction of GWout flow that reaches the channel
                }

            }
           // Qbase->Drc *= 2.0;

            ChannelWaterVol->Drc += Qbase->Drc;
            //GWVol->Drc -= Qbase->Drc;
            GWVol->Drc = std::max(0.0, GWVol->Drc - Qbase->Drc);
            GWWH->Drc = GWVol->Drc/CHAdjDX->Drc/pore->Drc;
            // m3 added per timestep, adjust the volume and height

            // NOTE: flow is always added no matter the conditions! e.g. when GW is below surface - channeldepth!
            // But that would make channeldepth very sensitive

        }}
    }
}
//---------------------------------------------------------------------------
void TWorld::ChannelRainandInfil(void)
{
    // add rainfall to channel, assume no interception
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (SwitchCulverts && ChannelMaxQ->Drc > 0)
            ChannelWaterVol->Drc += 0;
        else
            ChannelWaterVol->Drc += Rainc->Drc*ChannelWidth->Drc*DX->Drc;

       // ChannelWaterVol->Drc += ChannelQSide->Drc;
        // add unsaturated side inflow

    }}

    // subtract infiltration, no infil in culverts
// TODO: no infiltration if moisture content or GW does not allow this
// TODO: infiltration has to change moisture in surrounding soil
    if (SwitchChannelInfil) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            if (ChannelMaxQ->Drc <= 0) {
                double inf = std::min(ChannelWaterVol->Drc, ChannelInfM3->Drc);
                // cannot be more than there is
                ChannelWaterVol->Drc -= inf;
                ChannelInfilVol->Drc = inf; // do not make infiltration cumulative, that is done in totals
            }
        }}
    }

    // add user channel inflow
    if (SwitchDischargeUser) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelWaterVol->Drc += QuserIn->Drc * _dt;
            // add user defined discharge
        }}
    }
}
//---------------------------------------------------------------------------
//! calc channelflow, ChannelDepth, kin wave
//! channel WH and V and Q are clculated before
void TWorld::ChannelFlow(void)
{

    // if (SwitchChannelKinwaveDt) {
    //     if (_dt_user > _dtCHkin) {
    //         double n = _dt_user/_dtCHkin;
    //         _dt = _dt_user/n;
    //     }
    // }

    // for (double t = 0; t < _dt_user; t+=_dt)
    // {
     //   double sumvol = getMassCH(ChannelWaterVol);

        #pragma omp parallel num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelQn->Drc = 0;
            QinKW->Drc = 0;
            tma->Drc = ChannelWaterVol->Drc;
        }}

        if (SwitchLinkedList) {

            ChannelQn->setAllMV();
            FOR_ROW_COL_LDDCH5 {
                Kinematic(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelAlpha, ChannelDX, ChannelMaxQ, ChannelMaxAlpha);
            }}
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                if (pcr::isMV(ChannelQn->Drc))
                    ChannelQn->Drc = 0;
            }}

        } else {
            // default
            KinematicExplicit(crlinkedlddch_, ChannelQ, ChannelQn, ChannelAlpha, ChannelDX, ChannelMaxQ, ChannelMaxAlpha);
        }


        // calc V and WH back from Qn (original width and depth)
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            //  ChannelQ->Drc = ChannelQn->Drc;  // NOT because needed in erosion!
            ChannelQn->Drc = std::min(QinKW->Drc+tma->Drc/_dt, ChannelQn->Drc);
            ChannelWaterVol->Drc = tma->Drc + (QinKW->Drc - ChannelQn->Drc)*_dt;
            ChannelWaterVol->Drc = std::max(0.0,ChannelWaterVol->Drc);
            // vol is previous + in - out

            // recalc to Qn for erosion kin wave?
            ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelWidth->Drc*ChannelDX->Drc);
            // new channel WH, use adjusted channelWidth

            double ChannelArea = ChannelWaterVol->Drc/ChannelDX->Drc;
            ChannelAlpha->Drc = ChannelQn->Drc > 1e-6 ? ChannelArea/std::pow(ChannelQn->Drc, 0.6) : ChannelAlpha->Drc;

            ChannelV->Drc = std::min(_CHMaxV, (ChannelArea > 1e-12 ? ChannelQn->Drc/ChannelArea : 0));

            // get the maximum for output
            maxChannelflow->Drc = std::max(maxChannelflow->Drc, ChannelQn->Drc);
            maxChannelWH->Drc = std::max(maxChannelWH->Drc, ChannelWH->Drc);
        }}
       // correctMassBalanceCH(sumvol,ChannelWaterVol);
//     }
//     _dt=_dt_user;
}

void TWorld::ChannelSedimentFlow()
{
    if (!SwitchErosion)
        return;

    //separate Suspended and baseload for separate transport
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        ChannelQsn->Drc = 0;
        double concss = MaxConcentration(ChannelWaterVol->Drc, ChannelSSSed->Drc);
        ChannelQSSs->Drc = ChannelQ->Drc * concss; // m3/s *kg/m3 = kg/s
    }}

    if(SwitchUse2Phase) {
        #pragma omp parallel num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            double concbl = MaxConcentration(ChannelWaterVol->Drc, ChannelBLSed->Drc);
            ChannelQBLs->Drc = ChannelQ->Drc * concbl;
        }}
    }

    if (SwitchLinkedList) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            pcr::setMV(ChannelQSSsn->Drc);
        }}
        // advection SS
        FOR_ROW_COL_LDDCH5 {
              routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn, ChannelAlpha, ChannelDX, ChannelSSSed);
        }}

        //advection BL
        if(SwitchUse2Phase) {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                pcr::setMV(ChannelQBLsn->Drc);
            }}

            FOR_ROW_COL_LDDCH5 {
                routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn, ChannelAlpha, ChannelDX, ChannelBLSed);
            }}
        }

    } else {
        KinematicSubstance(crlinkedlddch_, LDDChannel, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn, ChannelAlpha, ChannelDX, ChannelSSSed);
        if(SwitchUse2Phase) {
            KinematicSubstance(crlinkedlddch_, LDDChannel, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn, ChannelAlpha, ChannelDX, ChannelBLSed);
        }
    }

    if (SwitchIncludeRiverDiffusion) {
        RiverSedimentDiffusion(_dt, ChannelSSSed, ChannelSSConc);
        // note SSsed goes in and out, SSconc is recalculated inside
    }

    // recalc all totals fluxes and conc
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        RiverSedimentLayerDepth(r,c);
        RiverSedimentMaxC(r,c);
        ChannelQsn->Drc = ChannelQSSsn->Drc + (SwitchUse2Phase ? ChannelQBLsn->Drc : 0);
        //ChannelSed->Drc = ChannelSSSed->Drc; //????? this is done in riversedmaxC
    }}
}




/* not used */
double TWorld::getMassCH(cTMap *M)
{
    double sum2 = 0;
   // #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        sum2 += M->Drc;
    }}
    return sum2;
}
/* not used */
void TWorld::correctMassBalanceCH(double sum1, cTMap *M)
{
    double sum2 = 0;
   // double n = 0;

  //  #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        sum2 += M->Drc;
  //      n += 1;
    }}
    // total and cells active for M
    double dhtot = fabs(sum2) > 0 ? (sum1 - sum2)/sum2 : 0;

    if (dhtot > 0) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            M->Drc = M->Drc*(1.0 + dhtot);            // <- distribution weighted to h
            M->Drc = std::max(M->Drc , 0.0);
        }}
    }
}
