

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
  \file lisKinematic2D.cpp
  \brief 2D kinematic wave routing functions and calculation of discharge and sed flux per cell.

  The routing functions use local variables to allow for possible re-usage
  Order of usage should be:

  K2DInit()
    (start loop)
    K2DDEMA()
    dt = K2DFlux()
    K2DSolveBy...(dt)
        K2DSolveBy....Sed(dt,S,C)
    K2DSolve(dt)
    (end loop)

@functions: \n
void TWorld::K2DInit()\n
double TWorld::K2DFlux()\n
void TWorld::K2DSolvebyInterpolation(double dt)\n
void TWorld::K2DSolvebyFlux(double dt)\n
double TWorld::K2DSolvebyInterpolationSed(double dt, cTMap *_S ,cTMap *_C)\n
double TWorld::K2DSolvebyFluxSed(double dt, cTMap *_S ,cTMap *_C)\n
void TWorld::K2DSolve(double dt)\n
void TWorld::K2DCalcVelDisch()\n
void TWorld::K2DDEMA()\n
 */

#include "model.h"
#include "operation.h"

#define MINGRAD 0.001

//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::K2DInit()
 * @brief Initializes the 2d kinematic wave solution
 *
 * Initializes all necessary variables for the 2d kinematic wave solution
 * Needs to be called every timestep.
 *
 * @return void
 */
void TWorld::K2DInit()
{
    //reset all maps for calculations
    FOR_ROW_COL_MV
    {
        if(SwitchErosion)
        {
            //K2DQS->Drc = 0;
            //K2DQSX->Drc = 0;
            //K2DQSY->Drc = 0;
//            K2DS->Drc = Sed->Drc;
//            K2DSC->Drc = Conc->Drc;
//            K2DSCN->Drc = 0;
//            K2DSFX->Drc = 0;
//            K2DSFY->Drc = 0;

        }
        if(SwitchPesticide)
        {
            K2DQP->Drc = 0;
            K2DQPX->Drc = 0;
            K2DQPY->Drc = 0;
            K2DP->Drc = Pest->Drc;
            K2DPC->Drc = C->Drc;
            K2DPCN->Drc = 0;
            //K2DPFX->Drc = 0;
            //K2DPFY->Drc = 0;
        }

        K2DHNew->Drc = 0;
        K2DQX->Drc = 0;
        K2DQY->Drc = 0;
        K2DQ->Drc = 0;
        K2DQN->Drc = 0;
        K2DFX->Drc = 0;
        K2DFY->Drc = 0;
        K2DI->Drc = 0;
    }

    K2DQOut = 0;
    K2DQSOut = 0;
    K2DQPOut = 0;

    //Recalculate water height to cover the full non-channel width of the cell!
    //this is a necessary step, since at the start of each LISEM timestep,
    //water is devided over the FlowWidth
    FOR_ROW_COL_MV
    {
        if(K2DSlope->Drc < MIN_SLOPE)
        {
            K2DQ->Drc = 0;
            Qn->Drc = 0;
            if(ChannelAdj->Drc > 0)
            {
                K2DHOld->Drc = WHrunoff->Drc*FlowWidth->Drc/ChannelAdj->Drc;
            }else
            {
                K2DHOld->Drc = 0;
            }


        }else
        {
            double dy = FlowWidth->Drc;

            double hrunoff = std::max(WHrunoff->Drc, 0.0);

            double Perim = /*2.0*hrunoff+*/ dy;

            if (Perim > 0)
                R->Drc = hrunoff;// *dy/Perim;
            else
                R->Drc = 0;

            Alpha->Drc = pow(N->Drc/sqrt(K2DSlope->Drc) * pow(Perim, (2.0/3.0)),0.6);

            if (Alpha->Drc == 0)
                K2DQ->Drc = 0;
            else
                K2DQ->Drc = pow((dy*hrunoff)/Alpha->Drc, 1.0/0.6);

            K2DQ->Drc = pow(R->Drc, 2.0/3.0)*sqrt(K2DSlope->Drc)/N->Drc *dy*hrunoff;

            //temporarily store in Qn (new)
            Qn->Drc = K2DQ->Drc;

            WHrunoff->Drc = (Alpha->Drc*pow(Qn->Drc, 0.6))/ChannelAdj->Drc;

            K2DHOld->Drc = WHrunoff->Drc;
        }
    }
}
//--------------------------------------------------------------------------------------------
/**
 * @fn double TWorld::K2DFlux()
 * @brief Calculates discharges for each cell
 *
 * Calculates the discharges in the current timestep
 * for each cell.
 * Must be preceded by K2DInit().
 * Afterwards K2DSolve....() must be called
 * (K2DSolveByInterpolation() or K2DSolveByFlux() )
 * to distrtibute flow, and finally K2DSolve() to calculate
 * new discharges.
 *
 * @return dt : the minimal needed timestep to ensure stability
 */
double TWorld::K2DFlux()
{
    double dtr = _dt;
    double fraction = CourantKin;
    FOR_ROW_COL_MV
    {

        double cdx = DX->Drc;
        double cdy = FlowWidth->Drc;
        //if a pit is present, set discharge to 0
        if( K2DPits->Drc == 1 || (K2DSlopeX->Drc == 0 && K2DSlopeY->Drc == 0))
        {
            K2DQ->Drc = 0;
            K2DQX->Drc = 0;
            K2DQY->Drc = 0;
            continue;
        }

        double hrunoff = std::max(K2DHOld->Drc - K2DWHStore->Drc, 0.0);

        //calculate discharge from water height, mannings N and slope
        double Perim =/* 2.0*hrunoff+*/  cdy;
        if (Perim > 0)
            R->Drc = hrunoff;//*cdy/Perim;
        else
            R->Drc = 0;

        // Q = V*A
       // K2DQ->Drc = sqrt(K2DSlope->Drc)/N->Drc * pow(R->Drc,2.0/3.0) * hrunoff*cdy;

        Alpha->Drc = pow(N->Drc/sqrt(K2DSlope->Drc) * pow(Perim, (2.0/3.0)),0.6);

        if(Alpha->Drc > 0)
            K2DQ->Drc = pow((cdy*hrunoff)/Alpha->Drc, 1.0/0.6);
        else
            K2DQ->Drc = 0;
        //within this timestep, only fraction of the cells available water should flow out
        if(K2DQ->Drc > 0)
        {
            double mindtr = fraction * (cdx*hrunoff*cdy)/K2DQ->Drc;
            dtr = std::min(mindtr ,dtr);
            if(!std::isnan(mindtr))
            {
                dtr = std::min(mindtr ,dtr);
            }
        }
    }

    dtr = std::max(dtr,TimestepKinMin);
    FOR_ROW_COL_MV
    {
//        double Qlim = fraction * DX->Drc * K2DHOld->Drc * ChannelAdj->Drc;
        double hrunoff = std::max(K2DHOld->Drc - K2DWHStore->Drc, 0.0);
        double Vollim = fraction * DX->Drc * hrunoff * ChannelAdj->Drc;   //why channeladj and not flowwidth
        //limit discharge to fraction of the cells water
        if(Vollim < K2DQ->Drc*dtr)
        {
            K2DQ->Drc = Vollim/dtr;
        }
        if(K2DOutlets->Drc == 1)
        {
            //K2DQ->Drc =std::min(0.5,KinematicBoundaryFraction*dtr) *  (DX->Drc*K2DHOld->Drc*ChannelAdj->Drc);
        }
    }
    return dtr;

}
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::K2DSolvebyInterpolation(double dt)
 * @brief Distributes flow trough interpolation
 *
 * Distributes flow trough bilinear interpolation.
 * Furthermore adds boundary outflow to K2DQout.
 * Must be preceded by K2DInit() and K2DFlux().
 * Afterwards K2DSolve must be called
 *
 * @param dt : timestep
 * @return void
 */
void TWorld::K2DSolvebyInterpolation(double dt)
{
    //this is the bilinear interpolated method!

    FOR_ROW_COL_MV
    {

        //start with old height and concentration
        K2DHNew->Drc = K2DHOld->Drc;

        K2DQN->Drc = 0;
    }

    //first calculate the weights for the cells that are closest to location that flow is advected to
    FOR_ROW_COL_MV
    {

        double DHL = sqrt(K2DSlopeX->Drc*K2DSlopeX->Drc + K2DSlopeY->Drc*K2DSlopeY->Drc);
        double dsx = K2DSlopeX->Drc/DHL;
        double dsy = K2DSlopeY->Drc/DHL;
        //the sign of the x and y direction of flow
        double yn = dsy/fabs(dsy);
        double xn = dsx/fabs(dsx);

        if(dsx == 0){xn = 1.0;};
        if(dsy == 0){yn = 1.0;};

        double cdx = DX->Drc;
        double cdy = ChannelAdj->Drc;

        //cell directions
        int dx[4] = {0, 1, 1, 0};
        int dy[4] = {1, 0, 1, 0};

        //weights to be saved
        double w[4] = {0.0,0.0,0.0,0.0};

        //wich directions are used?
        int end = 3;
        int start = 0;

        if(K2DPitsD->Drc == 1)
        {
            w[2] = 1.0;
        }else
        {
            //for each cell neigbhouring the advected location of the discharge, calculate interpolation weight
            int big = 0;
            for (int i=0; i<4; i++)
            {
                // distance we want is equal to: 1 - distance from the advected location to the neighbouring cell
                double wdx = ((double)1.0) - fabs( xn * ((double)dx[i]) - dsx);
                double wdy = ((double)1.0) - fabs( yn * ((double)dy[i]) - dsy);

                //the distribution is inverly proportional to the squared distance
                double weight = fabs(wdx*wdy);


                if (i > 0 && weight > w[i-1])
                    big = i;
//                // find largest weight direction

                w[i] = weight;
                if(SwitchFlowBarriers)
                {
                    w[i] = w[i] * FBW(K2DHOld->Drc,r,c,yn*dy[i],xn*dx[i]);
                }
            }
            // multiply with user weight
            w[big] *= ConcentrateKin;

            //normalize: sum of the 4 weights is equal to 1
            double wt = 0.0;

            wt = w[0] + w[1] + w[2] + w[3];
            if( wt == 0)
                wt = 1;
            w[0] /= wt;
            w[1] /= wt;
            w[2] /= wt;
            w[3] /= wt;
        }

        if(K2DOutlets->Drc == 1)
        {
           K2DQOut +=  dt*(K2DQ->Drc);
           QoutKW->data[r][c] += dt*K2DQ->Drc;
           K2DHNew->data[r][c] -=  dt*(K2DQ->Drc/(cdx*cdy));
        }else
        {
            //use the calculated weights to distribute flow
            for (int i=start; i<end+1; i++)
            {
                int r2, c2;

                //must multiply the cell directions by the sign of the slope vector components
                r2 = r+yn*dy[i];
                c2 = c+xn*dx[i];

                if(INSIDE(r2,c2) && !pcr::isMV(LDD->data[r2][c2]))
                {
                    double cdx2 = DX->data[r2][c2];
                    double cdy2 = ChannelAdj->data[r2][c2];

                    //weight * the flow is distributed to the ith cell that neighbours the advected flow.
                    K2DHNew->data[r2][c2] +=  w[i]*dt*(K2DQ->Drc/(cdx2*cdy2));
                    K2DHNew->data[r][c] -=  w[i]*dt*(K2DQ->Drc/(cdx*cdy));
                    QinKW->data[r2][c2] += w[i]*dt*K2DQ->Drc;
                    QoutKW->data[r2][c2] += w[i]*dt*K2DQ->Drc;
                }
            }
        }
    }

    //finish by substracting infiltration, and calculating discharge from new water height
//    FOR_ROW_COL_MV
//    {
//        double cdx = DX->Drc;
//        double cdy = ChannelAdj->Drc;

//        //calculate infiltration in time step
//        double infil = -1.0*FSurplus->Drc*dt/_dt;
//        if (K2DHNew->Drc < infil)
//            infil = K2DHNew->Drc;
//        K2DHNew->Drc -= infil;
//        FSurplus->Drc += infil*SoilWidthDX->Drc/cdy;
//        FSurplus->Drc = std::min(0.0, FSurplus->Drc);

// //        double infil = std::min(FSurplus->Drc *SoilWidthDX->Drc*DX->Drc * dt/_dt,0.0);
// //        if(K2DHNew->Drc < fabs(infil)/(cdx*cdy))
// //        {
// //            infil = -K2DHNew->Drc*(cdx*cdy);
// //        }

//        //keep track of infiltration
//        K2DI->Drc += (infil*cdx*cdy);
//        //K2DHNew->Drc = std::max(K2DHNew->Drc + infil/(cdx*cdy) ,0.0);

//        if(K2DHNew->Drc < 0)  // prob never occurs
//        {
//            K2DHNew->Drc = 0;
//            K2DQ->Drc = 0;
//            Qn->Drc = K2DQ->Drc;
//            WHrunoff->Drc = K2DHNew->Drc;
//            K2DHOld->Drc = K2DHNew->Drc;
//            qDebug() << r << c << K2DHNew->Drc<<  "WH negative?!!!";
//        }
//        else{

//            //Qn->Drc = K2DQ->Drc;
//            WHrunoff->Drc = K2DHNew->Drc;
//            K2DHOld->Drc = K2DHNew->Drc;
//        }
//    }
}

//--------------------------------------------------------------------------------------------
/**
 * @fn double TWorld::K2DSolvebyInterpolationSed(double dt, cTMap *_S ,cTMap *_C)
 * @brief Distributes sediment flow trough bilinear interpolation.
 *
 * Distributes sediment trough bilinear interpolation.
 * The concentration map is used to determine sediment discharges.
 * Must be called after K2DSolveBy...(), so that new cell water contents are known.
 *
 * @param dt : timestep
 * @param _S : Sediment map
 * @param _C : Concentration map
 * @return Total boundary outflow in this timestep
 */
double TWorld::K2DSolvebyInterpolationSed(double dt, cTMap *_S ,cTMap *_C)
{
    double K2DQMOut = 0;

    FOR_ROW_COL_MV
    {
        K2DM->Drc = _S->Drc;
        K2DQM->Drc = 0;
        K2DQMX->Drc = 0;
        K2DQMY->Drc = 0;
        //K2DFMX->Drc = 0;
        //K2DFMY->Drc = 0;
        K2DMN->Drc = K2DM->Drc;
        K2DMC->Drc = _C->Drc;
        K2DMC->Drc = MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, K2DM->Drc);

    }


    //double dtr = dt;
    double fraction = CourantKin;


    FOR_ROW_COL_MV
    {
       K2DQM->Drc =  K2DQ->Drc * K2DMC->Drc;
       if(K2DOutlets->Drc != 1 && K2DQM->Drc > fraction*K2DM->Drc)
       {
            K2DQM->Drc = fraction*K2DM->Drc;
       }

    }


    FOR_ROW_COL_MV
    {

        double DHL = sqrt(K2DSlopeX->Drc*K2DSlopeX->Drc + K2DSlopeY->Drc* K2DSlopeY->Drc);
        double dsx = K2DSlopeX->Drc/DHL;
        double dsy = K2DSlopeY->Drc/DHL;
        double yn = dsy/fabs(dsy);
        double xn = dsx/fabs(dsx);

        if(dsx == 0){xn = 1.0;};
        if(dsy == 0){yn = 1.0;};

        //cell directions
        int dx[4] = {0, 1, 1, 0};
        int dy[4] = {1, 0, 1, 0};

        double w[4] = {0.0,0.0,0.0,0.0};

        int end = 3;
        int start = 0;

        if(K2DPitsD->Drc == 1)
        {
            w[2] = 1.0;
        }else
        {
            //for each cell niegbhouring the advected location of the discharge, calculate interpolation weight
            int big = 0;
            for (int i=0; i<4; i++)
            {
                // distance we want is equal to: 1 - distance from the advected location to the neighbouring cell
                double wdx = ((double)1.0) - fabs( xn * ((double)dx[i]) - dsx);
                double wdy = ((double)1.0) - fabs( yn * ((double)dy[i]) - dsy);

                //the distribution is inverly proportional to the squared distance
                double weight = fabs(wdx*wdy);
                w[i]  = weight;
                if (i > 0 && weight > w[i-1])
                    big = i;

                if(SwitchFlowBarriers)
                {
                    w[i] = w[i] * FBW(K2DHOld->Drc,r,c,yn * dy[i],xn * dx[i]);
                }
            }
            // multiply with user weight
            w[big] *= ConcentrateKin;

            //normalize: sum of the 4 weights is equal to 1
            double wt = 0.0;

            wt = w[0] + w[1] + w[2] + w[3];
            if( wt == 0)
                wt = 1;
            w[0] /= wt;
            w[1] /= wt;
            w[2] /= wt;
            w[3] /= wt;
        }


        if(K2DOutlets->Drc == 1)
        {

            K2DMN->data[r][c] -=dt*(K2DQM->Drc);
            K2DQMOut +=dt*(K2DQM->Drc);
        }else
        {
            for (int i=start; i<end+1; i++)
            {
                int r2, c2;

                //must multiply the cell directions by the sign of the slope vector components
                r2 = r+yn*dy[i];
                c2 = c+xn*dx[i];

                //if the cell that flows needs to go to is out of bounds or missing value, skip
                bool inside = INSIDE(r2,c2);
                bool ismv = true;
                if(inside)
                    ismv = pcr::isMV(LDD->data[r2][c2]);

                if(inside && !ismv)
                {
                    K2DMN->data[r][c] -=w[i]*dt*(K2DQM->Drc);
                    K2DMN->data[r2][c2]+=w[i]*dt*(K2DQM->Drc);

                }
            }
        }
    }

    FOR_ROW_COL_MV
    {
        K2DQM->Drc = 0;
        K2DQMX->Drc = 0;
        K2DQMY->Drc = 0;
        //K2DFMX->Drc = 0;
       // K2DFMY->Drc = 0;
        _S->Drc = K2DMN->Drc;
        if(_C != NULL)
        {
            _C->Drc = MaxConcentration(K2DHNew->Drc * ChannelAdj->Drc * DX->Drc, K2DMN->Drc);

        }
    }


    //return outflow out of the catchment boundary
    return K2DQMOut;

}
//---------------------------------------------------------------------------
/**
 * @fn void TWorld::K2DSolve(double dt)
 * @brief Finalizes the solution for the kinematic wave
 *
 * Finalizes the solution for the kinematic wave,
 * sets variables for other LISEM code,
 * recalculates water height by substracting potential infiltration
 * and calculates new discharge.
 *
 * @param dt : timestep
 * @return void
 */
void TWorld::K2DSolve(double dt)
{
    //finish by substracting infiltration, and calculating discharge from new water height
    FOR_ROW_COL_MV
    {
//        double cdx = DX->Drc;
//        double cdy = ChannelAdj->Drc;

//        //calculate infiltartion in time step
//        double infil = std::min(FSurplus->Drc *SoilWidthDX->Drc*DX->Drc * dt/_dt,0.0);
//        if(K2DHNew->Drc < fabs(infil)/(cdx*cdy))
//        {
//            infil = -K2DHNew->Drc*(cdx*cdy);
//        }
//        //keep track of infiltration
//        K2DI->Drc -= (infil);
//        K2DHNew->Drc = std::max(K2DHNew->Drc + infil/(cdx*cdy) ,0.0);
//    }

            double cdx = DX->Drc;
            double cdy = ChannelAdj->Drc;

            //calculate infiltration in time step
            double infil = -1.0*FSurplus->Drc*dt/_dt;
            if (K2DHNew->Drc < infil)
                infil = K2DHNew->Drc;
            K2DHNew->Drc -= infil;
            FSurplus->Drc += infil*SoilWidthDX->Drc/cdy;
            FSurplus->Drc = std::min(0.0, FSurplus->Drc);

            Fcum->Drc += infil*SoilWidthDX->Drc/cdy; //VJ !!!

            //keep track of infiltration
            K2DI->Drc += (infil*cdx*cdy);
    }

    FOR_ROW_COL_MV
    {
        if(K2DHNew->Drc < 0)  // prob never occurs
        {
            K2DHNew->Drc = 0;
            K2DQ->Drc = 0;
            Qn->Drc = K2DQ->Drc;
            WHrunoff->Drc = K2DHNew->Drc;
            K2DHOld->Drc = K2DHNew->Drc;
            qDebug() << r << c << K2DHNew->Drc<<  "WH negative?!!!";
        }
        else{

            //Qn->Drc = K2DQ->Drc;
            WHrunoff->Drc = K2DHNew->Drc;
            K2DHOld->Drc = K2DHNew->Drc;
        }
    }
}
//---------------------------------------------------------------------------
/**
 * @fn void TWorld::K2DCalcVelDisch()
 * @brief calculates velocity and discharge for slopes
 *
 * calculates velocity and discharge for slopes
 * automatically called instead of CalcVelDisch()
 * when 2-D kinematic wave is used.
 * This function takes pits into account!
 *
 * @param dt : timestep
 * @return void
 * @see K2DWHStore
 */
void TWorld::K2DCalcVelDisch()
{
    FOR_ROW_COL_MV
    {
        if(K2DPits->Drc == 1 || K2DSlope->Drc < MIN_SLOPE)
        {
            Q->Drc = 0;
            V->Drc = 0;

        }else
        {
            double hrunoff = std::max(WHrunoff->Drc - K2DWHStore->Drc,0.0);

            double Perim;
            const double beta = 0.6;
            const double _23 = 2.0/3.0;
            double beta1 = 1/beta;
            double NN = N->Drc;


            if (SwitchChannelFlood)
                NN = N->Drc * qExp(mixing_coefficient*hmx->Drc);
            // slow down water in flood zone
            //    tma->Drc = hmx->Drc * UVflood->Drc/kinvisc;
            // Reynolds number ==> turbulent

            // avg WH from soil surface and roads, over width FlowWidth
            Perim = /* 2.0*hrunoff+ */ FlowWidth->Drc;

            if (Perim > 0)
                R->Drc = hrunoff*FlowWidth->Drc/Perim;
            else
                R->Drc = 0;

            Alpha->Drc = pow(NN/sqrt(K2DSlope->Drc) * pow(Perim, _23),beta);

            if (Alpha->Drc > 0)
                Q->Drc = pow((FlowWidth->Drc*hrunoff)/Alpha->Drc, beta1);
            else
                Q->Drc = 0;

            V->Drc = pow(R->Drc, _23)*sqrt(K2DSlope->Drc)/NN;

            if (K2DOutlets->Drc == 1)
            {
                V->Drc = (FlowWidth->Drc*hrunoff > 0 ? Q->Drc/(FlowWidth->Drc*hrunoff) : 0);
                V->Drc = std::min(V->Drc, hrunoff/_dt);
            }
            // limit the velocity on the outlet! can be extreme
        }


        //tm->Drc = V->Drc * R->Drc/kinvisc;
        //Reynolds number
    }
}

//---------------------------------------------------------------------------
/**
 * @fn void TWorld::OUTORMV()
 * @brief Returns wether the location is mv or outside the domain
 *
 * @param r : row number
 * @param c : column number
 * @see K2DDEMA
 */
bool TWorld::OUTORMV(int r, int c)
{
    if(INSIDE(r,c))
    {
        if(!pcr::isMV(LDD->data[r][c]))
        {
            return false;
        }
    }
    return true;
}

//---------------------------------------------------------------------------
/**
 * @fn void TWorld::K2DDEMA()
 * @brief Dem Analysis for 2-D kinematic wave
 *
 * Dem Analysis for 2-D kinematic wave.
 * Calculates slopes based on averages from both sides.
 * Calculates possible pit presence and depth.
 * Marks cells where only diagonal flow can take place.
 *
 * @return void
 * @see K2DWHStore
 * @see K2DPits
 * @see K2DPitsD
 */

void TWorld::K2DDEMA()
{
    FOR_ROW_COL_MV
    {
        K2DOutlets->Drc = 0;
        K2DWHStore->Drc = 0;
        K2DPits->Drc = 0;
        K2DPitsD->Drc = 0;
        K2DSlopeX->Drc = 0;
        K2DSlopeY->Drc = 0;
        K2DSlope->Drc = 0;
        //set heigt to dem + water heigt (this provides mannings flow equation with the gradient of water head)
        K2DDEM->Drc = DEM->Drc + WHrunoff->Drc;
    }

    FOR_ROW_COL_MV
    {
        double Dhx = 0;
        double Dhy = 0;

        //DEM
        double dem = DEMFB(r,c,0,0,true);

        double demx1 = DEMFB(r,c,0,1,true); //look right
        double demx2 = DEMFB(r,c,0,-1,true); // look left
        double demy1 = DEMFB(r,c,1,0,true);
        double demy2 = DEMFB(r,c,-1,0,true);

        if(OUTORMV(r,c+1)) // returns true if outside rows. cols or mv
        {
           // Outlet->Drc= 1;
            if(demx1 < demx2)
                K2DOutlets->Drc = 1;
        }
        if(OUTORMV(r,c-1))
        {
           // Outlet->Drc= 1;
            if(demx2 <demx1)
                K2DOutlets->Drc = 1;
        }

        if(OUTORMV(r+1,c))
        {
            //Outlet->Drc= 1;
            if(demy1 < demy2)
                K2DOutlets->Drc = 1;
        }
        if(OUTORMV(r-1,c))
        {
           // Outlet->Drc= 1;
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
            //Outlet->Drc= 1;
            K2DOutlets->Drc = 1;
        }
        if(OUTORMV(r+1,c) && OUTORMV(r-1,c))
        {
            Dhy = 0;
            //Outlet->Drc= 1;
            K2DOutlets->Drc = 1;
        }

        //at boundaries, set cell as outflow cell when slope is in the direction of the boundary

        if(r == 0)
        {
            if( Dhy < 0)
            {
               // Outlet->Drc= 1;
                K2DOutlets->Drc = 1;
            }
        }

        if(r == _nrRows-1)
        {
            if( Dhy > 0)
            {
             //  Outlet->Drc= 1;
               K2DOutlets->Drc = 1;
            }
        }

        if(c == 0)
        {
            if( Dhx < 0)
            {
                //Outlet->Drc= 1;
                K2DOutlets->Drc = 1;
            }
        }

        if(c == _nrCols-1)
        {
            if( Dhx > 0)
            {
                //Outlet->Drc= 1;
                K2DOutlets->Drc = 1;
            }
        }

        K2DSlopeX->Drc = Dhx/_dx;
        K2DSlopeY->Drc = Dhy/_dx;

        //calculate actual combined slope, with a minimum value of 0.01
        double Dh = fabs(Dhx) + fabs(Dhy);
        K2DSlope->Drc = Dh / sqrt(2*_dx*_dx);

        if(std::isnan(K2DSlope->Drc))
        {
            K2DSlope->Drc = 0.00;
            K2DSlopeX->Drc = 0.00;
            K2DSlopeY->Drc = 0.00;

        }

        //minimum value for the slope
        //K2DSlope->Drc = std::max(K2DSlope->Drc,0.01);
     //   K2DOutlets->Drc = Outlet->Drc;
    }

    //Detection of water available for outflow (because of local depressions)
    FOR_ROW_COL_MV
    {
        //cell directions
        int dx[8] = {0, 0, -1, 1,1,1,-1,-1};
        int dy[8] = {1, -1, 0, 0,1,-1,1,-1};
        bool pitxw= true;
        bool pityw= true;
        bool pitdw = true;
        int direction = 0;
        int mv = 0;
        double dem = DEMFB(r,c,0,0,false);
        double demw = DEMFB(r,c,0,0,true);
        double lowestneighbor = 9999999;
        double lowestneighborw = 9999999;
        int r2, c2;
        for (int i=0; i<8; i++)
        {
            //set row and column to neighbor
            r2 = r+dy[i];
            c2 = c+dx[i];
            if(INSIDE(r2,c2))
            {
                if(!pcr::isMV(LDD->data[r2][c2]))
                {
                    double demtnw = DEMFB(r,c,dy[i],dx[i],false);
                    if(demtnw <  lowestneighbor)
                    {
                        lowestneighbor = demtnw;
                    }

                    //if at least 1 neighboring cell is lower, it is not a pit
                    double demtw = DEMFB(r,c,dy[i],dx[i],true);
                    if(demtw < demw)
                    {
                        if( i < 2){
                            pitxw = false;
                        }else if( i < 4){
                            pityw = false;
                        }else
                        {
                            pitdw = false;
                            direction = i;
                        }
                    }

                    if(demtw <  lowestneighborw)
                    {
                        lowestneighborw = demtw;
                    }
                }else
                {
                    double tdem = dem;
                    if(tdem <  lowestneighbor)
                    {
                        lowestneighbor = tdem;
                    }

                    //if at least 1 neighboring cell is lower, it is not a pit
                    if(tdem < demw)
                    {
                        if( i < 2){
                            pitxw = false;
                        }else if( i < 4){
                            pityw = false;
                        }else
                        {
                            pitdw = false;
                            direction = i;
                        }
                    }
                    if(tdem <  lowestneighborw)
                    {
                        lowestneighborw = tdem;
                    }
                    mv++;
                }
            }
        }

        //only allow non-boundary cells, to be pits!
        if(mv == 0)
        {
            if(lowestneighborw > DEM->Drc)
            {
                double pitheight = lowestneighborw - DEM->Drc;
                K2DWHStore->Drc = pitheight;
            }

            if(pitxw &&pityw && !pitdw)
            {
                K2DPitsD->Drc = 1;
                K2DSlope->Drc = std::max((demw - lowestneighborw)/(sqrt(2)*_dx),0.01);
                K2DSlopeX->Drc = double(dx[(int)direction]) * K2DSlope->Drc/2.0;
                K2DSlopeY->Drc = double(dy[(int)direction]) * K2DSlope->Drc/2.0;

            }

            if(pitxw && pityw && pitdw)
            {
                K2DPits->Drc = 1;
                K2DSlope->Drc = 0;
                K2DSlopeX->Drc = 0;
                K2DSlopeY->Drc = 0;
            }

        }else
        {
            if(pitxw && pityw && pitdw)
            {
                K2DSlopeX->Drc = 0;
                K2DSlopeY->Drc = 0;
            }
        }
    }

    //VJ use flowboundary map, type 1 is open flow, else use the map
    if (FlowBoundaryType != 1)
        copy(*K2DOutlets, *FlowBoundary);  //copy 1 is 2

}
