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

   // ChannelRainandInfil();          // subtract infil, add rainfall

    //ChannelBaseflow();              // add stationary and GW baseflow if selected

    // _dt_user = _dt;
    // _dt = _dx/2.0;
    // for (double t = 0; t < _dt_user; t+=_dt)
    // {

    ChannelVelocityandDischarge();  // maaings V Q Aplha

    ChannelFlowDetachmentNew();     // detachment, deposition for SS and BL

    ChannelFlow();                  // channel kin wave for water

    //}

    //_dt = _dt_user;

    ChannelSedimentFlow();          // kin wave for sediment and substances

}
//---------------------------------------------------------------------------
void TWorld::ChannelVelocityandDischarge()
{
    // velocity, alpha, Q
    #pragma omp parallel num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        switch (crch_[i_].shape) {
            case SHAPERECT : ChannelPerimeter->Drc = ChannelWidthO->Drc+2*ChannelWH->Drc;
                ChannelWH->Drc = ChannelWaterVol->Drc/(ChannelDX->Drc*ChannelWidthO->Drc);
                // use real perimeter for velocity, not chanHandPRect(r,c,Area);
                break;
            case SHAPECIRC : chanHandPCirc(r,c); break; // this is always a culvert!
            case SHAPETRAP : chanHandPTrap(r,c); break;
            case SHAPETRIA : chanHandPTria(r,c); break;
        }
        double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
        double Radius = (ChannelPerimeter->Drc > 1e-6 ? Area/ChannelPerimeter->Drc : 0);
        ChannelV->Drc = std::min(_CHMaxV,std::pow(Radius, 2.0/3.0)*sqrt(ChannelGrad->Drc)/ChannelN->Drc);
        ChannelQ->Drc = ChannelV->Drc * Area;
        //ChannelAlpha->Drc = ChannelQ->Drc/std::pow(Area, 0.6);
        ChannelAlpha->Drc = pow(ChannelN->Drc/sqrt(ChannelGrad->Drc) * pow(ChannelPerimeter->Drc, 2.0/3.0),0.6);  // no difference
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

            if (!crch_[i_].culvert) {
                ChannelWaterVol->Drc += Qbase->Drc;
                GWVol->Drc = std::max(0.0, GWVol->Drc - Qbase->Drc);
                GWWH->Drc = GWVol->Drc/CHAdjDX->Drc/pore->Drc;
            }
            // m3 added per timestep, adjust the volume and height, not in culverts

            // NOTE: flow is always added no matter the conditions! e.g. when GW is below surface - channeldepth!
            // But that would make channeldepth very sensitive

        }}
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
            ChannelWaterVol->Drc += Rainc->Drc*ChannelWidth->Drc*DX->Drc;
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
                }
                ChannelInfM3->Drc = ChannelPerimeter->Drc * ChannelKsat->Drc * _dt/3600000.0 * ChannelDX->Drc;
                // infiltration over entire perimeter !
                double inf = std::min(ChannelWaterVol->Drc, ChannelInfM3->Drc);
                // cannot be more than there is
                ChannelWaterVol->Drc -= inf;
                ChannelInfilVol->Drc = inf;
                // do not make infiltration cumulative, that is done in totals
                // TODO check this
            }
        }}
    }

    if (SwitchGridRetention) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_CHL {
            if (GridRetention->Drc > 0) {
                double dvol = std::max(0.0,GridRetention->Drc - GridRetentionAct->Drc);
                if(dvol > 0) {
                    dvol = std::min(dvol, ChannelWaterVol->Drc);
                    if (dvol > 0) {
                        GridRetentionAct->Drc += dvol;
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
void TWorld::ChannelFlow(void)
{
    int dy[10] = {0,1,1,1,0,0,0,-1,-1,-1};
    int dx[10] = {0,-1,0,1,-1,0,1,-1,0,1};

  //  double sumvol = MapTotal(*ChannelWaterVol);
  //  double totq = 0;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        ChannelQn->Drc = 0;
        QinKW->Drc = 0; // needed for sediment
    }}

    for(long i_ =  0; i_ < crlinkedlddch_.size(); i_++)
    {
        int r = crlinkedlddch_.at(i_).r;
        int c = crlinkedlddch_.at(i_).c;
        double Qin = 0;
        double volMax = ChannelMaxArea->Drc*DX->Drc;

        if (crlinkedlddch_.at(i_).nr > 0) {
            for(int j = 0; j < crlinkedlddch_.at(i_).nr; j++) {
                int rr = crlinkedlddch_.at(i_).inn[j].r;
                int cr = crlinkedlddch_.at(i_).inn[j].c;
                Qin += ChannelQn->Drcr;
            }

            // if total inflow causes vol > max volume, adjust inflow incoming TileQn
            // if !switchculverts then ChannelCulvert has only 0
            if (ChannelCulvert->Drc > 0 &&
                ChannelWaterVol->Drc+_dt*(Qin-ChannelQ->Drc) >= volMax) {
                double maxq = std::min(ChannelMaxQ->Drc, (volMax - ChannelWaterVol->Drc)/_dt + ChannelQ->Drc);

                for(int j = 0; j < crlinkedlddch_.at(i_).nr; j++) {
                    int rr = crlinkedlddch_.at(i_).inn[j].r;
                    int cr = crlinkedlddch_.at(i_).inn[j].c;
                    ChannelQn->Drcr = maxq * ChannelQn->Drcr/Qin;
                    // incoming TileQn is a fraction of maxq
                }
                Qin = maxq;
            }
        }
        QinKW->Drc = Qin;

        if (!SwitchCulverts)
            ChannelQn->Drc = IterateToQnew(Qin, ChannelQ->Drc, ChannelAlpha->Drc, _dt, DX->Drc, 0,0);
        else
            ChannelQn->Drc = IterateToQnew(Qin, ChannelQ->Drc, ChannelAlpha->Drc, _dt, DX->Drc, ChannelMaxQ->Drc, ChannelMaxAlpha->Drc);
        ChannelQn->Drc = std::min(Qin+ChannelWaterVol->Drc/_dt, ChannelQn->Drc);
        // no more outflow than there is water

        // check if there is a culvert downstream and limit outflow if necessary
        int ldd = fabs(crlinkedlddch_.at(i_).ldd);
        int cr = c+dx[ldd];
        int rr = r+dy[ldd];
        if (!pcr::isMV(LDDChannel->Drcr) && ChannelCulvert->Drcr > 0)
            ChannelQn->Drc = std::min(ChannelQn->Drc, ChannelMaxQ->Drcr);

    }
    // int full = 0;

    // calc V and WH back from Qn (original width and depth)
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        ChannelWaterVol->Drc = ChannelWaterVol->Drc + _dt*(QinKW->Drc - ChannelQn->Drc);
        ChannelWaterVol->Drc = std::max(0.0, ChannelWaterVol->Drc);

     //   if (ChannelCulvert->Drc > 0 && ChannelWaterVol->Drc >= ChannelMaxArea->Drc*DX->Drc) {
     //       full+=1;
     //   }
        switch (crch_[i_].shape) {
            case SHAPERECT : chanHandPRect(r,c); break;
            case SHAPECIRC : chanHandPCirc(r,c); break; // this is always a culvert!
            case SHAPETRAP : chanHandPTrap(r,c); break;
            case SHAPETRIA : chanHandPTria(r,c); break;
        }
        double Area = ChannelWaterVol->Drc/ChannelDX->Drc;
        ChannelV->Drc = std::min(_CHMaxV, (Area > 1e-12 ? ChannelQn->Drc/Area : 0.0));
        // ChannelAlpha->Drc = Area > 1e-6 ? ChannelQn->Drc/std::pow(Area, 0.6) : 0.0;
        // DO NOT recalculate alpha becuase of erosion

        // get the maximum for output
        maxChannelflow->Drc = std::max(maxChannelflow->Drc, ChannelQn->Drc);
        maxChannelWH->Drc = std::max(maxChannelWH->Drc, ChannelWH->Drc);

        //   if (LDDChannel->Drc == 5)
     //        totq += ChannelQn->Drc*_dt;
    }}
//    double sumvol1 = MapTotal(*ChannelWaterVol);
 //   qDebug() << "MB chan (aft-bef)" << sumvol << sumvol1 << totq << sumvol - sumvol1 - totq << MB << full;

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

    // if (SwitchLinkedList) {
    //     #pragma omp parallel for num_threads(userCores)
    //     FOR_ROW_COL_MV_L {
    //         pcr::setMV(ChannelQSSsn->Drc);
    //     }}
    //     // advection SS
    //     FOR_ROW_COL_LDDCH5 {
    //           routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn, ChannelAlpha, ChannelDX, ChannelSSSed);
    //     }}

    //     //advection BL
    //     if(SwitchUse2Phase) {
    //         #pragma omp parallel for num_threads(userCores)
    //         FOR_ROW_COL_MV_L {
    //             pcr::setMV(ChannelQBLsn->Drc);
    //         }}

    //         FOR_ROW_COL_LDDCH5 {
    //             routeSubstance(r,c, LDDChannel, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn, ChannelAlpha, ChannelDX, ChannelBLSed);
    //         }}
    //     }

    // } else {
        KinematicSubstance(crlinkedlddch_, LDDChannel, ChannelQ, ChannelQn, ChannelQSSs, ChannelQSSsn, ChannelAlpha, ChannelDX, ChannelSSSed, ChannelMaxQ);
        if(SwitchUse2Phase) {
            KinematicSubstance(crlinkedlddch_, LDDChannel, ChannelQ, ChannelQn, ChannelQBLs, ChannelQBLsn, ChannelAlpha, ChannelDX, ChannelBLSed, ChannelMaxQ);
        }
//    }

    if (SwitchIncludeRiverDiffusion) {
        RiverSedimentDiffusion(_dt, ChannelSSSed, ChannelSSConc);
        // note SSsed goes in and out, SSconc is recalculated inside
    }

    // recalc all totals fluxes and conc
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_CHL {
        if (ChannelSSSed->Drc > MAXCONC * ChannelWaterVol->Drc) {
            double ss = ChannelSSSed->Drc;
            ChannelSSSed->Drc = MAXCONC * ChannelWaterVol->Drc;
            double ds = ss - ChannelSSSed->Drc;
            ChannelDep->Drc -= ds;
        }


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
            M->Drc = std::max(M->Drc , 0.0);
        }}
    }
}
