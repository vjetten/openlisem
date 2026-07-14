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
void TWorld:: ChannelFlowandErosion()
{
    if (!SwitchIncludeChannel)
        return;

    // moved to before overland flow
    ChannelRainandInfil();          // subtract infil, add rainfall
    ChannelBaseflow();              // add stationary and GW baseflow if selected

    // looping a smaller dt doesn't work or doesn't make difference
    // _dt_user = _dt;
    // for (double t = 0; t < _dt_user; t+=_dt)
    // {

        ChannelVelocityandDischarge();  // mannings V Q Aplha

        if (SwitchDepositionContinuous)
            ChannelDetachmentContinuous();     // detachment, deposition for SS and BL
        else
            ChannelFlowDetachment();     // detachment, deposition for SS and BL

        ChannelFlow();                  // kin wave for water

        ChannelSedimentFlow();          // kin wave for sediment and substances

        // restore _dt
    // _dt = _dt_user;

}
//---------------------------------------------------------------------------
void TWorld::ChannelVelocityandDischarge()
{
  //  int dy[10] = {0,1,1,1,0,0,0,-1,-1,-1};
  //  int dx[10] = {0,-1,0,1,-1,0,1,-1,0,1};
    // velocity, alpha, Q
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        double beta = BETArect;
        switch (crch_[i_].shape) {
            case SHAPEFREE :
            case SHAPERECT : ChannelPerimeter->Drc = ChannelWidthO->Drc+2*ChannelWH->Drc;
                // use real perimeter for velocity, not chanHandPRect(r,c,Area);
                ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelDX->Drc*ChannelWidthO->Drc);
                if (!SwitchConstantBeta)
                    beta = 1.0/(1.0+2.0/3.0*ChannelWidthO->Drc/ChannelPerimeter->Drc);
                break;
            case SHAPECIRC : beta = BETAcirc; chanHandPCirc(r,c); break; // this is always a culvert!
            case SHAPETRAP : beta = BETAtrap; chanHandPTrap(r,c); break;
            case SHAPETRIA : beta = BETAtria; chanHandPTria(r,c); break;
        }
        double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
        double Radius = (ChannelPerimeter->Drc > 1e-6 ? Area/ChannelPerimeter->Drc : 0);
        ChannelV->Drc = qMin(_CHMaxV,std::pow(Radius, 2.0/3.0)*qSqrt(ChannelGrad->Drc)/ChannelN->Drc);
        ChannelQ->Drc = ChannelV->Drc * Area;
        ChannelAlpha->Drc = pow(ChannelN->Drc/qSqrt(ChannelGrad->Drc) * pow(ChannelPerimeter->Drc, 2.0/3.0),beta);  // no difference

    }}
}

//---------------------------------------------------------------------------
void TWorld::ChannelBaseflow(void)
{
    // add a stationary part
    if(SwitchChannelBaseflowStationary) {
        // add switch for baseflow as map
        // if added as map then addedbaseflow = true;
        //if (SwitchChannelBaseflowMap)
       //     addedbaseflow = true;

        //CHECK?

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
        if(InfilMethod == INFIL_SWATRE) {
    //todo
        } else {

            GroundwaterFlow();
            // is all based on 2 layer G&A !
            //TODO 3 layer G&A

            cTMap *SD = SoilDepth1init;
            cTMap *pore = Thetaeff;
            if (SwitchTwoLayer) {
                SD = SoilDepth2init;
                pore = ThetaS2;
            }
            // in all channel cells
            // move groundwater, GWout is the flow between cells
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
                       //Qbase->Drc = qMin(GWVol->Drc, 2.0 * dH/GWWH->Drc * GWout->Drc);
                    //   Qbase->Drc = qMin(GWVol->Drc, 2.0 * fabs(GWout->Drc));
                       Qbase->Drc = 2*GWout->Drc;
                       // use the fraction of GWout flow that reaches the channel
                    }
                }
               // Qbase->Drc *= 2.0;

                if (!crch_[i_].culvert) {
                    ChannelWaterVol->Drc += Qbase->Drc;
                    GWVol->Drc = qMax(0.0, GWVol->Drc - Qbase->Drc);
                    GWWH->Drc = GWVol->Drc/CHAdjDX->Drc/pore->Drc;
                }
                // m3 added per timestep, adjust the volume and height, not in culverts

                // NOTE: flow is always added no matter the conditions! e.g. when GW is below surface - channeldepth!
                // But that would make channeldepth very sensitive

            }}
        }
    }
}
//---------------------------------------------------------------------------
// Channel volume with incoming rainfall and outgoing infil and retention
void TWorld::ChannelRainandInfil(void)
{
    // add rainfall to channel, assume no interception
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (!crch_[i_].culvert)
            ChannelWaterVol->Drc += Rainc->Drc*ChannelWidth->Drc*ChannelDX->Drc;
        // goes for all channel shapes

       // ChannelWaterVol->Drc += ChannelQSide->Drc;
       // add unsaturated side inflow
    }}

    // subtract infiltration, no infil in culverts
    if (SwitchChannelInfil) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            if (!crch_[i_].culvert) {
                switch (crch_[i_].shape) {
                    case SHAPERECT : chanHandPRect(r,c); break;
                    case SHAPECIRC : chanHandPCirc(r,c); break;
                    case SHAPETRAP : chanHandPTrap(r,c); break;
                    case SHAPETRIA : chanHandPTria(r,c); break;
                    case SHAPEFREE : chanHandPRect(r,c); break;
                }
                ChannelInfM3->Drc = ChannelPerimeter->Drc * ChannelKsat->Drc * _dt/3600000.0 * ChannelDX->Drc;
                // infiltration over entire perimeter !
                double inf = qMin(ChannelWaterVol->Drc, ChannelInfM3->Drc);
                // cannot be more than there is
                ChannelWaterVol->Drc -= inf;
                ChannelInfilVol->Drc = inf;
                // do not make infiltration cumulative, that is done in totals
                // TODO check this
            }
        }}
    }

// TODO: chan retention for sediment!
    if (SwitchGridRetention) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            if (ChanRetention->Drc > 0) {
                double dvol = qMax(0.0,ChanRetention->Drc - ChanRetentionAct->Drc);
                if(dvol > 0) {
                    dvol = qMin(dvol, ChannelWaterVol->Drc);
                    if (dvol > 0) {
                        ChanRetentionAct->Drc += dvol;
                        ChannelWaterVol->Drc -= dvol;
                    }
                }
            }
        }}
    }

    // add user channel inflow
    if (SwitchDischargeUser) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelWaterVol->Drc += QuserIn->Drc * _dt;
            // add user defined discharge
            //TODO add outlet dicharge from stromdrrains
        }}
    }
}
//---------------------------------------------------------------------------
// NOTE for shapes: https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/transform/kinematic-wave-transform-model
// in these eq, alpha is 1/alpha in lisem, beta = 1/m
void TWorld::ChannelFlow(void)
{
    int dy[10] = {0,1,1,1,0,0,0,-1,-1,-1};
    int dx[10] = {0,-1,0,1,-1,0,1,-1,0,1};

    bool extrapressure = true;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        ChannelQn->Drc = 0;
        QinKW->Drc = 0; // needed for sediment
        tma->Drc = ChannelMaxQ->Drc;
        tmb->Drc = ChannelMaxAlpha->Drc;
        tmc->Drc = 0;
    }}

    for(long i_ =  0; i_ < crlinkedlddch_.size(); i_++)
    {
        int r = crlinkedlddch_.at(i_).r;
        int c = crlinkedlddch_.at(i_).c;
        double Qin = 0;
        double volMax = ChannelMaxArea->Drc*DX->Drc;

        if (crlinkedlddch_.at(i_).nr > 0) {
            //get all inflow
            for(int j = 0; j < crlinkedlddch_.at(i_).nr; j++) {
                int rr = crlinkedlddch_.at(i_).inn[j].r;
                int cr = crlinkedlddch_.at(i_).inn[j].c;
                Qin += ChannelQn->Drcr;
            }

            // if total inflow causes vol > max volume, adjust inflow incoming Qn
            // if !switchculverts then ChannelCulvert has only 0
            if (ChannelCulvert->Drc > 0 && ChannelCulvert->Drc < 5 &&
                ChannelWaterVol->Drc+_dt*(Qin-ChannelQ->Drc) >= volMax) {
                // if volume +in-out is more than maxvol, recalc maxq = inflow
                double maxq = qMin(tma->Drc, (volMax - ChannelWaterVol->Drc)/_dt + ChannelQ->Drc);
                //double maxq = qMin(ChannelMaxQ->Drc, (volMax - ChannelWaterVol->Drc)/_dt + ChannelQ->Drc);

                for(int j = 0; j < crlinkedlddch_.at(i_).nr; j++) {
                    int rr = crlinkedlddch_.at(i_).inn[j].r;
                    int cr = crlinkedlddch_.at(i_).inn[j].c;
                    ChannelQn->Drcr = maxq * ChannelQn->Drcr/Qin;
                    // incoming Qn is a fraction of maxq
                }
                Qin = maxq;
            }
        }
        QinKW->Drc = Qin;

        double beta = BETArect;
        switch ((int)ChannelCulvert->Drc) {
            case SHAPEFREE :
            case SHAPERECT : beta = SwitchConstantBeta ? BETArect : beta = 1.0/(1.0+2.0/3.0*ChannelWidth->Drc/ChannelPerimeter->Drc);
            case SHAPECIRC : beta = BETAcirc; break;
            case SHAPETRAP : beta = BETAtrap; break;
            case SHAPETRIA : beta = BETAtria; break;

        }
        if (ChannelCulvert->Drc == 0 || ChannelCulvert->Drc == 5) //!SwitchCulverts) //
            ChannelQn->Drc = IterateToQnew(Qin, ChannelQ->Drc, ChannelAlpha->Drc, beta, _dt, DX->Drc, 0,0);
        else
            ChannelQn->Drc = IterateToQnew(Qin, ChannelQ->Drc, ChannelAlpha->Drc, beta, _dt, DX->Drc, tma->Drc, tmb->Drc);
        ChannelQn->Drc = qMin(Qin+ChannelWaterVol->Drc/_dt, ChannelQn->Drc);
        // no more outflow than there is water

        // check if there is a culvert downstream and limit outflow if necessary
        int ldd = fabs(crlinkedlddch_.at(i_).ldd);
        int cr = c+dx[ldd];
        int rr = r+dy[ldd];
        if (!pcr::isMV(LDDChannel->Drcr) && ChannelCulvert->Drcr > 0 && ChannelCulvert->Drcr < 5) {
            ChannelQn->Drc = qMin(ChannelQn->Drc, tma->Drcr);

            // adjust discharge and max discharge when pressure of water is more than diameter
            if (tmc->Drc > 0) {
                // if we are in the culvert and there is extra discharge, add it to downstream maxQ
                tma->Drcr += tmc->Drc;
                // adjust max Q for downstream cells
                tmb->Drcr = ChannelMaxArea->Drcr/tma->Drcr;
                // adjust maxalpha for downstream cells, beta is 1.0 for fully submerged, so Q^beta is not necessary
            }

            if (ChannelWH->Drc > ChannelDiameter->Drcr * 1.1) {
                double dQ = 0.67 * 2 * GRAV * ChannelWH->Drc-ChannelDiameter->Drcr;
                // simplified for sharp entry and short pipe. 0.67 = Cd
                tma->Drcr += dQ;
                // simply add this to the downstream cell MaxQ as entry
                tmb->Drcr = ChannelMaxArea->Drcr/tma->Drcr;
                // adjust maxalpha for downstream cells, beta is 1.0 for fully submerged, so Q^beta is not necessary
                tmc->Drcr = dQ;
                // save extra discharge
            }
        }

    }

    // calc V and WH back from Qn (original width and depth)
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {

        ChannelWaterVol->Drc = ChannelWaterVol->Drc + _dt*(QinKW->Drc - ChannelQn->Drc);
        ChannelWaterVol->Drc = qMax(0.0, ChannelWaterVol->Drc);

        if(ChannelWaterVol->Drc == 0 && ChannelQn->Drc > 0) {
            ChannelWaterVol->Drc = ChannelDX->Drc * ChannelAlpha->Drc*qPow(ChannelQn->Drc, BETArect);
        }

        // calc  channel WH and perimeter
        switch (crch_[i_].shape) {
            case SHAPERECT : chanHandPRect(r,c); break;
            case SHAPECIRC : chanHandPCirc(r,c); break; // this is always a culvert!
            case SHAPETRAP : chanHandPTrap(r,c); break;
            case SHAPETRIA : chanHandPTria(r,c); break;
            case SHAPEFREE : chanHandPRect(r,c); break;
        }
        double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
        ChannelV->Drc = qMin(_CHMaxV, (Area > 1e-20 ? ChannelQn->Drc/Area : 0.0));
        // erosion is calculated with new V
    //    if(ChannelV->Drc == 0 && ChannelQn->Drc > 0)
        //    qDebug() << r<<c<<"Q" << ChannelQn->Drc << Area << ChannelWaterVol->Drc << QinKW->Drc;

        // ChannelAlpha->Drc = Area > 1e-6 ? ChannelQn->Drc/std::pow(Area, 0.6) : 0.0;
        // DO NOT recalculate alpha after the kin wave because we need it in erosion kin wave

        // get the maximum for output
        maxChannelflow->Drc = qMax(maxChannelflow->Drc, ChannelQn->Drc);
        maxChannelWH->Drc = qMax(maxChannelWH->Drc, ChannelWH->Drc);

    }}
}

void TWorld::ChannelSedimentFlow()
{
    if (!SwitchErosion)
        return;

    //separate Suspended and baseload for separate transport
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        ChannelQSSsn->Drc = 0;
        double concss = MaxConcentration(ChannelWaterVol->Drc, ChannelSSSed->Drc);
        ChannelQSSs->Drc = ChannelQ->Drc * concss; // m3/s *kg/m3 = kg/s
    }}

    if(SwitchUse2Phase) {
        #pragma omp parallel num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            ChannelQBLsn->Drc = 0;
            double concbl = MaxConcentration(ChannelWaterVol->Drc, ChannelBLSed->Drc);
            ChannelQBLs->Drc = ChannelQ->Drc * concbl;
        }}
    }

    // if (SwitchLinkedList) {
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

    // } else {


        // KinematicSubstance(crlinkedlddch_, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn, ChannelAlpha, ChannelDX, ChannelSSSed, ChannelMaxQ);
        // if(SwitchUse2Phase) {
        //     KinematicSubstance(crlinkedlddch_, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn, ChannelAlpha, ChannelDX, ChannelBLSed, ChannelMaxQ);
        // }
//    }

    if (SwitchIncludeRiverDiffusion) {
        RiverSedimentDiffusion(_dt, ChannelSSSed, ChannelSSConc);
        // note SSsed goes in and out, SSconc is recalculated inside
    }

    // recalc all totals fluxes and conc
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {/*
        if (ChannelSSSed->Drc > MAXCONC * ChannelWaterVol->Drc) {
            double ss = ChannelSSSed->Drc;
            ChannelSSSed->Drc = MAXCONC * ChannelWaterVol->Drc;
            double ds = ss - ChannelSSSed->Drc;
            ChannelDep->Drc -= ds;
        }*/


        RiverSedimentLayerDepth(r,c);
        RiverSedimentMaxC(r,c);
        ChannelQsn->Drc = ChannelQSSsn->Drc + (SwitchUse2Phase ? ChannelQBLsn->Drc : 0);

    }}
}




/* not used */
double TWorld::getMassCH(cTMap *M)
{
    double sum2 = 0;
    #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (ChannelCulvert->Drc == 0)
            sum2 += M->Drc;
    }}
    return sum2;
}
/* not used */
void TWorld::correctMassBalanceCH(double sum1, cTMap *M)
{
    double sum2 = 0;

    #pragma omp parallel for reduction(+:sum2) num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (ChannelCulvert->Drc == 0)
            sum2 += M->Drc;
    }}
    // total and cells active for M
    double dhtot = fabs(sum2) > 0 ? (sum1 - sum2)/sum2 : 0;

    if (dhtot > 0) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            M->Drc = M->Drc*(1.0 + dhtot);            // <- distribution weighted to h
            M->Drc = qMax(M->Drc , 0.0);
        }}
    }
}

