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
- void TWorld::ToChannel(void)\n
- void TWorld::CalcVelDisch(void)\n
- void TWorld::OverlandFlow(void)\n
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
    if(SwitchKinematic2D == K2D_METHOD_KIN || SwitchKinematic2D == K2D_METHOD_KINDYN) {
        ToChannel();           // overland flow water and sed flux going into or out of channel, in channel cells

        OverlandFlow1D();   // kinematic wave

        if(SwitchKinematic2D == K2D_METHOD_KINDYN)
            ChannelFlood(); // st venant channel 2D flooding from channel, only for kyn wave, partly parallel
    }

    if(SwitchKinematic2D == K2D_METHOD_DYN) {
        OverlandFlow2Ddyn();
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
void TWorld::ToChannel()
{
    if (!SwitchIncludeChannel)
        return;
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(ChannelMaskExtended->data[r][c] == 1)
        {
            int rr = r;//(int)ChannelSourceYExtended->Drc;
            int cr = c;//(int)ChannelSourceXExtended->Drc;

            double fractiontochannel;

            if (WHrunoff->Drc < HMIN)
                continue;
            if (hmx->Drc > HMIN)
                continue;



            double VtoChan = V->Drc;
            if (F_AddGravity == 1)
                VtoChan = std::pow(WHrunoff->Drcr, 2.0/3.0)*sqrt(ChannelPAngle->Drc)/N->Drcr; //F_Angle
            fractiontochannel = std::min(1.0, _dt*VtoChan/std::max(0.05*_dx,0.5*ChannelAdj->Drc));
            // fraction to channel calc from half the adjacent area width and flow velocity

            // cannot flow into channel if water level in channel is higher than runoff depth
            if (SwitchKinematic2D == K2D_METHOD_KINDYN &&
                WHrunoff->Drc <= std::max(0.0 , ChannelWH->Drcr - ChannelDepth->Drcr))
                fractiontochannel = 0;

            // no inflow on culverts
            if (SwitchCulverts && ChannelMaxQ->Drcr  > 0)
                fractiontochannel = 0;

            if (fractiontochannel == 0)
                continue;

            double dwh = fractiontochannel*WHrunoff->Drc;

            // water diverted to the channel
            ChannelWaterVol->Drcr += dwh* FlowWidth->Drc * DX->Drc;
            fromChannelVoltoWH(rr, cr);

            WHrunoff->Drc -= dwh ;
            WHroad->Drc -= dwh;
            //WHGrass->Drc -= dwh;
            WH->Drc -= dwh;
            WaterVolall->Drc = WHrunoff->Drcr*CHAdjDX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;

            if (SwitchErosion)
            {
                double dsed = fractiontochannel*Sed->Drcr;
                ChannelSSSed->Drcr  += dsed;
                //sediment diverted to the channel
                Sed->Drcr -= dsed;

                Conc->Drcr = MaxConcentration(WHrunoff->Drcr * DX->Drcr * ChannelAdj->Drcr, &Sed->Drcr, &DEP->Drcr);
                // adjust sediment in suspension

//                if(SwitchUseGrainSizeDistribution)
//                {
//                    Conc->Drcr = 0;
//                    FOR_GRAIN_CLASSES
//                    {
//                        RSS_D.Drcdr += fractiontochannel * Sed_D.Drcdr;
//                        Sed_D.Drcd = Sed_D.Drcd * (1-fractiontochannel);
//                        Conc_D.Drcd = MaxConcentration(WHrunoff->Drcr * DX->Drcr * ChannelAdj->Drcr, &Sed_D.Drcd, &DEP->Drcr);
//                        Conc->Drcr += Conc_D.Drcd;
//                    }
//                }
               RiverSedimentLayerDepth(rr,cr);
               RiverSedimentMaxC(rr,cr);
            }
        }
    }}
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
	if(SwitchKinematic2D == K2D_METHOD_DYN)
		return;

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
      //  double Perim, R;
        double NN = N->Drc;
        double alpha;
        double WHr = WHrunoff->Drc;
        double FW = FlowWidth->Drc;

        if (SwitchIncludeChannel && hmx->Drc > 0.001)
            NN = N->Drc * (2.0-qExp(-mixing_coefficient*hmx->Drc));
        // slow down water in flood zone, if hmx = 0 then factor = 1

        // avg WH from soil surface and roads, over width FlowWidth
     //   Perim = FW; // ssurface flow does not have sides!only roughness  2*WHr + FW;

//        if (Perim > 0)
//            R = WHr*FW/Perim;
//        else
//            R = 0;
       // R = WHr;

        if (Grad->Drc > MIN_SLOPE)
            alpha = pow(NN/sqrtGrad->Drc * pow(FW, 2.0/3.0),0.6);
        else
            alpha = 0;

        if (alpha > 0)
            Q->Drc = pow((FW*WHr)/alpha, 5.0/3.0);
        else
            Q->Drc = 0;

        V->Drc = pow(WHr, 2.0/3.0) * sqrtGrad->Drc/NN;
        Alpha->Drc = alpha;

    }}
}

//---------------------------------------------------------------------------
// DO NOT MAKE PARALLEL
void TWorld::Boundary2Ddyn()//cTMap* h, cTMap* Q, cTMap *_U, cTMap *_V)
{
    cTMap *h = WHrunoff;
    cTMap *Q = Qn;
    cTMap *_U = Uflood;
    cTMap *_V = Vflood;

    if(SwitchKinematic2D == K2D_METHOD_KINDYN) {
        Q = Qflood;
        h = hmx;
    }

    K2DQOutBoun = 0;
    K2DQSOutBoun = 0;
    fill(*tma, 0);
    fill(*K2DQ, 0);

    // find oulets based on DEM and WHrunoff
    dynOutflowPoints();

    if (FlowBoundaryType > 0) {
        //direction of velocity is in the direction of + and -
        // U is EW and V is NS
        // find which outlets on the boundary are directed to the outside based on sign U and V
        FOR_ROW_COL_MV_L {
            if (K2DOutlets->Drc == 1 && FlowBoundary->Drc == 1 && h->Drc > 0.0)
            {
                if (c > 0 && MV(r,c-1))
                    if (_U->Drc < 0)
                        tma->Drc = 1;
                if (c < _nrCols-1 && MV(r,c+1))
                    if (_U->Drc > 0)
                        tma->Drc = 1;
                if (r > 0 && MV(r-1,c))
                    if (_V->Drc < 0)
                        tma->Drc = 1;
                if (r < _nrRows-1 && MV(r+1,c))
                    if (_V->Drc > 0)
                        tma->Drc = 1;
            }
        }}
    }

//    FOR_ROW_COL_MV {
//        if(LDD->Drc == 5)
//            tma->Drc = 1;
//    }
    FOR_ROW_COL_LDD5 {
        tma->Drc = 1;
    }}

    if(SwitchIncludeChannel) {
        //FOR_ROW_COL_MV_CH {
        //  if(LDDChannel->Drc == 5)
        FOR_ROW_COL_LDDCH5 {
            tma->Drc = 1;
        }}
    }


    FOR_ROW_COL_MV_L {
        if (tma->Drc == 1) {
            double dy = ChannelAdj->Drc;
            double UV = qSqrt(_U->Drc * _U->Drc + _V->Drc*_V->Drc);
            double frac = std::min( std::max(0.0, UV*_dt/DX->Drc) , 0.9);
            double dh = frac*h->Drc;
            double _q = dh*DX->Drc*dy;

            K2DQOutBoun += _q/_dt;
            h->Drc -= dh;
            K2DQ->Drc = _q/_dt;
            Q->Drc -= _q/_dt;

            if (SwitchErosion) {
                double ds = frac * SSFlood->Drc;
                K2DQSOutBoun += ds;
                SSFlood->Drc -= ds;

                ds = frac * BLFlood->Drc;
                K2DQSOutBoun += ds;
                BLFlood->Drc -= ds;
            }
        }
    }}
}
//---------------------------------------------------------------------------
void TWorld::SolveDeepWH(void)
{
    if (F_pitValue < 0)
        return;
    long cnt = 0;

//#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (DEMdz->Drc == 1 && WHrunoff->Drc > F_pitValue) {
            int ldd = (int) LDD->Drc;
            int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1,  0,  1};
            int dy[10] = {0,  1, 1, 1,  0, 0, 0, -1, -1, -1};
            double f = 0.2;

            double H = WHrunoff->Drc;
            double H1 = WHrunoff->data[r+dy[ldd]][c+dx[ldd]];
            double dH = H*f;

           // double dH = 0;// _dt*(f*H*_dx*std::pow(f*H, 0.667)*sqrt(Grad->Drc)/N->Drc)/(_dx*DX->Drc);
            WHrunoff->data[r+dy[ldd]][c+dx[ldd]] += dH;
            WHrunoff->Drc -= dH;

            Vflood->Drc *= WHrunoff->Drc/H;
            Uflood->Drc *= WHrunoff->Drc/H;
            Vflood->data[r+dy[ldd]][c+dx[ldd]] *= WHrunoff->data[r+dy[ldd]][c+dx[ldd]]/H1;
            Uflood->data[r+dy[ldd]][c+dx[ldd]] *= WHrunoff->data[r+dy[ldd]][c+dx[ldd]]/H1;

          //  qDebug() << r << c << ldd << H << dH << WHrunoff->data[r+dy[ldd]][c+dx[ldd]];
            cnt++;
        }
    }}
    if (cnt > 0) qDebug() << cnt;
}
//---------------------------------------------------------------------------
void TWorld::OverlandFlow2Ddyn(void)
{
    double dtOF = 0;

    ChannelOverflow(WHrunoff, V, false);
        // false means flood sediment maps are used

    startFlood = false;
#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if (WHrunoff->Drc > HMIN){
            startFlood = true;
          //  break;
        }
    }}


  //   SolveDeepWH();

    if (SwitchSWOFopen)
        dtOF = fullSWOF2open(WHrunoff, Uflood, Vflood, DEM);
    else
        dtOF = fullSWOF2RO(WHrunoff, Uflood, Vflood, DEM);
    //VJ new average flux over lisem timestep, else last Qn is used

    //  infilInWave(WHrunoff, _dt);
#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        V->Drc = qSqrt(Uflood->Drc*Uflood->Drc + Vflood->Drc*Vflood->Drc);
        Qn->Drc = V->Drc*(WHrunoff->Drc*ChannelAdj->Drc);
        Q->Drc = Qn->Drc; // just to be sure
    }}

    Boundary2Ddyn();//WHrunoff, Qn, Uflood, Vflood);  // do the domain boundaries

#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double WHR = WHrunoff->Drc;

        Qn->Drc = V->Drc*(WHR*ChannelAdj->Drc);
        Q->Drc = Qn->Drc; // just to be sure

        WHroad->Drc = WHR;
        // set road to average outflowing wh, no surface storage.
        //WHGrass->Drc = WHR;

        WH->Drc = WHR + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water
        WaterVolall->Drc = WHR*CHAdjDX->Drc + WHstore->Drc*SoilWidthDX->Drc*DX->Drc;
        //LOGIC:
        // water layer in three parts: WHstore < (WHrunoff < minReportFloodHeight) < (hmx > minReportFloodHeight)
        // WH =  WHstore + WHrunoff

        hmxWH->Drc = WH->Drc;

        hmxflood->Drc = std::max(0.0, WHR - minReportFloodHeight);

        FloodWaterVol->Drc = hmxflood->Drc*CHAdjDX->Drc;
        //WHrunoffOutput->Drc = std::min(WHrunoff->Drc, minReportFloodHeight);
        //RunoffWaterVol->Drc = WHrunoffOutput->Drc*ChannelAdj->Drc*DX->Drc;
        RunoffWaterVol->Drc = std::min(WHR, minReportFloodHeight)*CHAdjDX->Drc;
        // used for screen output

        WHmax->Drc = std::max(WHmax->Drc, hmxWH->Drc);

    }}

    FloodMaxandTiming(hmxWH, V, minReportFloodHeight);

    if(SwitchErosion)
    {
        //calculate concentration and new sediment discharge
        //WHrunoff and Qn are adapted in case of 2D routing
//        if(!SwitchUseGrainSizeDistribution)
//        {
            #pragma omp parallel for num_threads(userCores)
            FOR_ROW_COL_MV_L {
                double sed = (SSFlood->Drc + BLFlood->Drc);
                Conc->Drc =  MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, &sed, &DepFlood->Drc);
                Qsn->Drc = Conc->Drc*Qn->Drc;
            }}
//        }
//        else
//        {
//            //calculate total sediment from induvidual grain classes,
//            //and calculate concentration and new sediment discharge
//            FOR_ROW_COL_MV {
//                Sed->Drc = 0;
//                Conc->Drc = 0;

//            }
//            FOR_ROW_COL_MV
//            {
//                FOR_GRAIN_CLASSES
//                {
//                    Sed->Drc += Sed_D.Drcd;
//                    Conc_D.Drcd = MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, &Sed_D.Drcd, &DepFlood->Drc);
//                    Conc->Drc += Conc_D.Drcd;
//                }
//                Qsn->Drc = Conc->Drc*Qn->Drc;
//            }
//        }
    }

    debug(QString("Average dynamic timestep (dt %1 sec, n %2)").arg(dtOF,6,'f',3).arg(iter_n,4));
    // some screen error reporting
}



//--------------------------------------------------------------------------------------------
void TWorld::OverlandFlow1D(void)
{
    // recalculate water vars after subtractions in "to channel"
#pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        WaterVolin->Drc = DX->Drc * FlowWidth->Drc * WHrunoff->Drc;
        //volume runoff into the kin wave, needed to determine infil in kin wave
        // WaterVolin total water volume in m3 before kin wave, WHrunoff may be adjusted in tochannel
        q->Drc = FSurplus->Drc*SoilWidthDX->Drc/_dt;  //???FlowWidth
        // infil flux in kin wave (<= 0)negative value), in m2/s, in kiv wave DX is used
        // surplus related to infiltrating surfaces
        QinKW->Drc = 0; // store for incoming water in a cell
        tm->Drc = -1;
        // flag for confined fow in channel culverts
    }}

    if (SwitchErosion)
    {
        // calc seediment flux going in kin wave as Qs = Q*C
#pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            Qsn->Drc = 0.0;
            Conc->Drc = MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, &Sed->Drc, &DEP->Drc);
            Qs->Drc =  Q->Drc * Conc->Drc;
            // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
        }}
    }

    // route water
    Qn->setAllMV();
//    FOR_ROW_COL_MV {
//        if (LDD->Drc == 5)
//        {

        FOR_ROW_COL_LDD5 {
            Kinematic(r,c, LDD, Q, Qn, q, Alpha, DX, tm);
            // tm is not used in overland flow, in channel flow it is the max flux of e.g. culverts
        }}
   // }

    //convert calculate Qn back to WH and volume for next loop
    fill(*tma, 0);
    FOR_ROW_COL_MV
    {
        bool K1Dexplicit = true;
        double WaterVolout = 0;
        double InfilKWact = 0;

        if( K1Dexplicit)
        {
//TODO what about q/surplus here????
            WaterVolout = std::max(0.0, QinKW->Drc*_dt + WaterVolin->Drc  - Qn->Drc*_dt);

            // new water vol is mass bal diff
            WHrunoff->Drc = ChannelAdj->Drc > 0 ? WaterVolout/(ChannelAdj->Drc*DX->Drc) : 0.0;
            // runoff based on water vol out
            // NOTE route substance is already an explicit solution
        }
        else
        {
            WHrunoff->Drc = ChannelAdj->Drc > 0 ?(Alpha->Drc*pow(Qn->Drc, 0.6))/ChannelAdj->Drc : 0.0;
            //new WH based on A/dx = alpha Q^beta / dx
            // apha is however determined from the old Q...

            WaterVolout = WHrunoff->Drc*CHAdjDX->Drc;
            // new volume
        }
        tma->Drc = WaterVolout;

        double diff = QinKW->Drc*_dt + WaterVolin->Drc - WaterVolout - Qn->Drc * _dt;
        InfilKWact = std::min(-FSurplus->Drc*SoilWidthDX->Drc*DX->Drc, diff);

        // infiltration is the surplus infil (pot infil), or infil is all that was there
//        if (FFull->Drc == 1)
//            InfilKWact = 0;

        InfilVolKinWave->Drc = InfilKWact;

        double Perim = FlowWidth->Drc > 0 ? 2*WHrunoff->Drc + FlowWidth->Drc : 0.0;
        double R = WHrunoff->Drc*FlowWidth->Drc/Perim;
        Alpha->Drc = pow(N->Drc/sqrtGrad->Drc * pow(Perim, 2.0/3.0),0.6); // for erosion
        V->Drc = pow(R, 2.0/3.0) * sqrtGrad->Drc/N->Drc;
    //    Qn->Drc = V->Drc * WHrunoff->Drc*FlowWidth->Drc;
   //     V->Drc = std::min(Qn->Drc/(WHrunoff->Drc*ChannelAdj->Drc), V->Drc);
//        V->Drc = Qn->Drc/(WHrunoff->Drc*ChannelAdj->Drc);

        WHroad->Drc = WHrunoff->Drc;
        // set road to average outflowing wh, no surface storage.

        //WHGrass->Drc = WHrunoff->Drc;

        WH->Drc = WHrunoff->Drc + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        hmxWH->Drc = WH->Drc;
        //needed for totals and output

        WaterVolall->Drc = WHrunoff->Drc*CHAdjDX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        // is the same as :         WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);

    }

    //      routing of substances add here!
    //      do after kin wave so that the new flux Qn out of a cell is known
    //      you need to have the ingoing substance flux QS (mass/s)
    //      and it will give outgoing flux QSn (mass/s)
    //      and the current amount Subs (mass) in suspension+solution

    if (SwitchErosion)
    {
        fill(*QinKW, 0.0); // reuse
        fill(*tm, 0.0);
        copy(*tm, *Sed);
        if(!SwitchUseGrainSizeDistribution)
        {
            Qsn->setAllMV();
       //     FOR_ROW_COL_MV
       //     {
       //         if (LDD->Drc == 5) // if outflow point, pit
       //         {

            FOR_ROW_COL_LDD5 {
                routeSubstance(r,c, LDD, Q, Qn, Qs, Qsn, Alpha, DX, WaterVolin, Sed);
                Conc->Drc = MaxConcentration(WaterVolall->Drc, &Sed->Drc, &DEP->Drc);
            }}
        //    }
        } else {
            /*
            FOR_GRAIN_CLASSES
            {
                // calc seediment flux going in kin wave as Qs = Q*C
                FOR_ROW_COL_MV
                {
                    Conc_D.Drcd = MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, &Sed_D.Drcd, &DEP->Drc);
                    Tempa_D.Drcd =  Q->Drc * Conc_D.Drcd;
                    Qs->Drc +=  Q->Drc * Conc_D.Drcd;
                    // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
                }
                Tempb_D.at(d)->setAllMV();
                FOR_ROW_COL_MV
                {
                    if (LDD->Drc == 5) // if outflow point, pit
                    {
                        routeSubstance(r,c, LDD, Q, Qn, Tempa_D.at(d), Tempb_D.at(d), Alpha, DX, WaterVolin, Sed_D.at(d));//, BufferVol, BufferSed);
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
*/
        }
    }

    // route other stuff
    if (SwitchPesticide)
    {
        // calc pesticide flux going in kin wave as Qp = Q*C
        FOR_ROW_COL_MV
        {
            Qp->Drc =  Qn->Drc * C->Drc;
            // calc sed flux as water flux * conc m3/s * kg/m3 = kg/s
        }

        fill(*Qpn, 0.0);
        FOR_ROW_COL_MV
        {
            if (LDD->Drc == 5) // if outflow point, pit
            {
                routeSubstance(r,c, LDD, Q, Qn, Qp, Qpn, Alpha, DX, WaterVolin, Pest);//, BufferVol, nullptr);
            }
        }
    }


    // new concentrations with new volume

    if (SwitchErosion)
    {
        FOR_ROW_COL_MV
        {

            Conc->Drc = MaxConcentration(WaterVolall->Drc, &Sed->Drc, &DEP->Drc);
            /*
            if(SwitchUseGrainSizeDistribution)
            {
                FOR_GRAIN_CLASSES
                {
                    Conc_D.Drcd = (Qn->Drc > MIN_FLUX ? Tempb_D.Drcd/Qn->Drc : 0);
                }
            }
*/
            if (SwitchPesticide)
            {
                //C->Drc = ConcentrationP(WaterVolall->Drc, Pest->Drc);
                C->Drc = Qn->Drc > MIN_FLUX ? Qpn->Drc/Qn->Drc : 0;
                C_N->Drc = C->Drc;
                //qDebug()<< "ds overlandflow"<< C->Drc;
                //qDebug()<< "ds overlandflow"<< Pest->Drc;
            }
        }
    }
}
//---------------------------------------------------------------------------
// all points that flow outward of the domain by slope and water pressure
void TWorld::dynOutflowPoints()
{
    FOR_ROW_COL_MV_L {
        K2DOutlets->Drc = 0;
    }}

    FOR_ROW_COL_MV_L {
        double Dhx = 0;
        double Dhy = 0;

        //DEM + water height and barriers if switched on
        double dem = DEMFB(r,c,0,0,true);

        double demx1 = DEMFB(r,c,0,1,true); //look right
        double demx2 = DEMFB(r,c,0,-1,true); // look left
        double demy1 = DEMFB(r,c,1,0,true);
        double demy2 = DEMFB(r,c,-1,0,true);

        if(OUTORMV(r,c+1)) // returns true if outside rows. cols or mv
        {
            if(demx1 < demx2)
                K2DOutlets->Drc = 1;
        }
        if(OUTORMV(r,c-1))
        {
            if(demx2 <demx1)
                K2DOutlets->Drc = 1;
        }

        if(OUTORMV(r+1,c))
        {
            if(demy1 < demy2)
                K2DOutlets->Drc = 1;
        }
        if(OUTORMV(r-1,c))
        {
            if(demy2 < demy1)
                K2DOutlets->Drc = 1;
        }

        if(demx1 < demx2)
        {
            Dhx = -(demx1-dem);
        }else
        {
            Dhx = (demx2-dem);
        }

        if(demy1 < demy2)
        {
            Dhy = -(demy1-dem);
        }else
        {
            Dhy = (demy2-dem);
        }

        if(OUTORMV(r,c+1) && OUTORMV(r,c-1))
        {
            Dhx = 0;
            K2DOutlets->Drc = 1;
        }
        if(OUTORMV(r+1,c) && OUTORMV(r-1,c))
        {
            Dhy = 0;
            K2DOutlets->Drc = 1;
        }

        //at boundaries, set cell as outflow cell when slope is in the direction of the boundary

        if(r == 0)
        {
            if( Dhy < 0)
            {
                K2DOutlets->Drc = 1;
            }
        }

        if(r == _nrRows-1)
        {
            if( Dhy > 0)
            {
               K2DOutlets->Drc = 1;
            }
        }

        if(c == 0)
        {
            if( Dhx < 0)
            {
                K2DOutlets->Drc = 1;
            }
        }

        if(c == _nrCols-1)
        {
            if( Dhx > 0)
            {
                K2DOutlets->Drc = 1;
            }
        }
    }}

    //VJ use flowboundary map, type 1 is open flow, else use the map
    if (FlowBoundaryType != 1) {
        FOR_ROW_COL_MV_L {
            K2DOutlets->Drc *= FlowBoundary->Drc;  //copy 1 is 2
        }}
    }
}
