

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
        K2DDTT->Drc = 0;
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
double TWorld::K2DFlux(double t, double tmax)
{
    double dtr = _dt;
    double fraction = courant_factor_diffusive;//CourantKin;
    FOR_ROW_COL_MV
    {
        K2DHNew->Drc = K2DHOld->Drc;
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
            K2DQ->Drc = /* 2.0 * */ pow((cdy*hrunoff)/Alpha->Drc, 1.0/0.6);
        else
            K2DQ->Drc = 0;


        //double v = pow(R->Drc,2.0/3.0)*sqrt(K2DSlope->Drc)/N->Drc;


        //within this timestep, only fraction of the cells available water should flow out
        if(K2DQ->Drc > 0)
        {
            //alternative courant condition with empricial stuff
            //the pow(v,0.3) is completely emperical, helps reduce oscillations due to spatially dynamic timestep
            //double mindtr =std::max(TimestepKinMin,2.0*fraction * (cdx)/(std::max(K2DQ->Drc,std::max(pow(v,0.3),std::max(v,(sqrt(hrunoff)* K2DSlope->Drc))))));

            double mindtr = std::min(_dt,std::max(_dx/(10.0 *hrunoff),std::max(TimestepKinMin,2.0 *fraction * (cdx*hrunoff*cdy)/K2DQ->Drc)));
          //  double mindtr = 2.0*fraction * (cdx*hrunoff*cdy)/K2DQ->Drc;

            //store this as the calculated timestep
            //mindtr = mindtr;

            //do this for temporal smoothing of timestep
            //std::min(std::max(TimestepKinMin,K2DDTr->Drc * 1.05),mindtr));


            K2DDTm->Drc = std::min(std::max(TimestepKinMin,mindtr), _dt);

            dtr = std::min(mindtr ,dtr);

        }
    }
    double dtr2 = 2.0 * dtr;
    dtr = std::min(std::max(dtr,TimestepKinMin),std::max(0.0,tmax - t));

    int rc = 0;
    int cc = 0;
    int count = 0;

    FOR_ROW_COL_MV
    {

        double dtm = std::max(dtr,std::max(dtr2,K2DDTm->Drc));//std::max(dtr,std::floor(K2DDTm->Drc/(2.0*dtr)) * (2.0*dtr));

        //Extended spatial averaging, so that there are less oscillation
        //Doesnt work so well, i think a Lax Anti-Oscillation fix (avering WH based on gradients) is needed


        double spatial_averaging = 1.0;
        
        if(!OUTORMV(r,c-1))
        {
            dtm = std::min(dtm,spatial_averaging *K2DDTm->data[r][c-1]);
        }
        if(!OUTORMV(r+1,c))
        {
            dtm = std::min(dtm,spatial_averaging *K2DDTm->data[r+1][c]);
        }
        if(!OUTORMV(r-1,c))
        {
            dtm = std::min(dtm,spatial_averaging *K2DDTm->data[r-1][c]);
        }

        if(!OUTORMV(r+1,c+1))
        {
            dtm = std::min(dtm,1.0*spatial_averaging *K2DDTm->data[r+1][c+1]);
        }
        if(!OUTORMV(r-1,c+1))
        {
            dtm = std::min(dtm,1.0*spatial_averaging *K2DDTm->data[r-1][c+1]);
        }
        if(!OUTORMV(r-1,c+1))
        {
            dtm = std::min(dtm,1.0*spatial_averaging *K2DDTm->data[r-1][c+1]);
        }
        if(!OUTORMV(r-1,c-1))
        {
            dtm = std::min(dtm,1.0*spatial_averaging *K2DDTm->data[r-1][c-1]);
        }

        dtm = std::max(dtm,TimestepKinMin);
        //write the real timestep
        K2DDTr->Drc = dtm;

        //write the used timestep (might be less than real timestep due to the lisem timestep
        //K2DDT->Drc =std::min( dtm,std::max(0.0,tmax-K2DDTT->Drc));


        //Set wether the timestep is done now or later
        if(!(t + dtr < tmax))
        {
            K2DDT->Drc = 2.0 * tmax-K2DDTT->Drc;
            {
                K2DDTR->data[rc][cc] = r;
                K2DDTC->data[rc][cc] = c;

                count ++;

                cc ++;
                if(cc == _nrCols)
                {
                    rc ++;
                    cc = 0;
                }

                if(INSIDE(rc,cc))
                {
                    K2DDTR->data[rc][cc] = -1;
                    K2DDTC->data[rc][cc] = -1;
                }
            }

        }else if(!(K2DDTT->Drc> t +dtr) )
        {
            K2DDT->Drc = 2.0 * std::max(K2DDT->Drc, tmax - K2DDTT->Drc);
            {
                K2DDTR->data[rc][cc] = r;
                K2DDTC->data[rc][cc] = c;

                count ++;

                cc ++;
                if(cc == _nrCols)
                {
                    rc ++;
                    cc = 0;
                }

                if(INSIDE(rc,cc))
                {
                    K2DDTR->data[rc][cc] = -1;
                    K2DDTC->data[rc][cc] = -1;
                }
            }

        }else
        {
            K2DDT->Drc = 0;
        }
    }

    return dtr;

}
//--------------------------------------------------------------------------------------------
void TWorld::K2DPreSolve(int thread)
{
    double fraction = courant_factor_diffusive;//CourantKin;
    FOR_ROW_COL_UF2DMT_DT
    {
        if(K2DDT->Drc > 1e-8)
        {
            double hrunoff = std::max(K2DHOld->Drc - K2DWHStore->Drc, 0.0);
            double Vollim = std::min(0.8,2.0 *fraction) * DX->Drc * hrunoff * ChannelAdj->Drc;   //why channeladj and not flowwidth
            //limit discharge to fraction of the cells water
            if(Vollim < K2DQ->Drc*K2DDT->Drc)
            {
                K2DQ->Drc = Vollim/K2DDT->Drc;
            }
            if(K2DOutlets->Drc == 1)
            {
                //K2DQ->Drc =std::min(0.5,KinematicBoundaryFraction*dtr) *  (DX->Drc*K2DHOld->Drc*ChannelAdj->Drc);
            }
        }else
        {
             K2DQ->Drc = 0;
        }
    }}}}
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
void TWorld::K2DSolvebyInterpolation(int thread)
{

    //this is the bilinear interpolated method!

    FOR_ROW_COL_UF2DMTDER
    {

        //start with old height and concentration
        K2DHNew->Drc = K2DHOld->Drc;

        K2DQN->Drc = 0;
    }}}}

    //first calculate the weights for the cells that are closest to location that flow is advected to
    FOR_ROW_COL_UF2DMT_DT
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
            w[big] *= /*2.5 * */ ConcentrateKin;

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
           K2DQOut +=  K2DDT->Drc*(K2DQ->Drc);
           //QoutKW->data[r][c] += K2DDT->Drc*K2DQ->Drc; // not used
           K2DHNew->data[r][c] -=  K2DDT->Drc*(K2DQ->Drc/(cdx*cdy));
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
                    K2DHNew->data[r2][c2] +=  w[i]*K2DDT->Drc*(K2DQ->Drc/(cdx2*cdy2));
                    K2DHNew->data[r][c] -=  w[i]*K2DDT->Drc*(K2DQ->Drc/(cdx*cdy));
                    QinKW->data[r2][c2] += w[i]*K2DDT->Drc*K2DQ->Drc;
                    //QoutKW->data[r2][c2] += w[i]*K2DDT->Drc*K2DQ->Drc;  //not used
                }
            }
        }
    }}}}

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
double TWorld::K2DSolvebyInterpolationSed(int thread, cTMap *_S ,cTMap *_C)
{
    double K2DQMOut = 0;

    FOR_ROW_COL_UF2DMTDER
    {
        K2DM->Drc = _S->Drc;
        K2DQM->Drc = 0;
        K2DQMX->Drc = 0;
        K2DQMY->Drc = 0;
        //K2DFMX->Drc = 0;
        //K2DFMY->Drc = 0;
        K2DMN->Drc = K2DM->Drc;
        K2DMC->Drc = _C->Drc;
        K2DMC->Drc = MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, &K2DM->Drc, &DEP->Drc);

    }}}}


    //double dtr = dt;
    double fraction = courant_factor_diffusive;//CourantKin;


    FOR_ROW_COL_UF2DMT_DT
    {
       K2DQM->Drc =  K2DQ->Drc * K2DMC->Drc;
       if(K2DOutlets->Drc != 1 && K2DQM->Drc > fraction*K2DM->Drc)
       {
            K2DQM->Drc = fraction*K2DM->Drc;
       }

    }}}}


    FOR_ROW_COL_UF2DMT_DT
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
            w[big] *= /* 2.5 * */ConcentrateKin;

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

            K2DMN->data[r][c] -=K2DDT->Drc*(K2DQM->Drc);
            K2DQMOut +=K2DDT->Drc*(K2DQM->Drc);
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
                    K2DMN->data[r][c] -=w[i]*K2DDT->Drc*(K2DQM->Drc);
                    K2DMN->data[r2][c2]+=w[i]*K2DDT->Drc*(K2DQM->Drc);

                }
            }
        }
    }}}}

    FOR_ROW_COL_UF2DMTDER
    {
        K2DQM->Drc = 0;
        K2DQMX->Drc = 0;
        K2DQMY->Drc = 0;
        //K2DFMX->Drc = 0;
       // K2DFMY->Drc = 0;
        _S->Drc = K2DMN->Drc;
        if(_C != nullptr)
        {
            _C->Drc = MaxConcentration(K2DHNew->Drc * ChannelAdj->Drc * DX->Drc, &K2DMN->Drc, &DEP->Drc);

        }
    }}}}

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
void TWorld::K2DSolve(int thread)
{
    //finish by substracting infiltration, and calculating discharge from new water height
    if (InfilMethod != INFIL_NONE)
    FOR_ROW_COL_UF2DMT_DT
    {
        double cdx = DX->Drc;
        double cdy = ChannelAdj->Drc;
        double frac = std::min(SoilWidthDX->Drc, cdy)/cdy;

        //calculate infiltration in time step
        double infil = -1.0*FSurplus->Drc*K2DDT->Drc/_dt;
        if (K2DHNew->Drc < infil)		
            infil = K2DHNew->Drc;
        K2DHNew->Drc -= infil;
        FSurplus->Drc += infil*frac;
        FSurplus->Drc = std::min(0.0, FSurplus->Drc);

        Fcum->Drc += infil*frac; //VJ !!!

        //keep track of infiltration
        K2DI->Drc += (infil*cdx*cdy);
    }}}}


    FOR_ROW_COL_UF2DMTDER
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

            Qn->Drc = K2DQ->Drc;
            WHrunoff->Drc = K2DHNew->Drc;
            K2DHOld->Drc = K2DHNew->Drc;
        }
    }}}}
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
void TWorld::K2DCalcVelDisch(int thread)
{
    FOR_ROW_COL_2DMT
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

            if(K2DOutlets->Drc == 1) //VJ  why zero at outlet???
            {
                V->Drc = std::min(V->Drc, hrunoff/_dt);									   
                Q->Drc = V->Drc*hrunoff*FlowWidth->Drc;
            }
        }
    }}}}
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



bool TWorld::OUTORMVc(int r, int c)
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

void TWorld::K2DDEMA(int thread)
{
    FOR_ROW_COL_UF2DMT
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
    }}}}

    FOR_ROW_COL_UF2DMT
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

        K2DSlopeX->Drc = Dhx/_dx;
        K2DSlopeY->Drc = Dhy/_dx;

        //calculate actual combined slope
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
    }}}}

    //Detection of water available for outflow (because of local depressions)
    FOR_ROW_COL_UF2DMT
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
    }}}}

    // adjust boundary slopes to avoid extremes
    FOR_ROW_COL_UF2DMT {
        if (DomainEdge->Drc > 0 && FlowBoundary->Drc == 0)
        {
            double Savg = getWindowAverage(*K2DSlope, r, c, false);
            K2DSlope->Drc = Savg;
        }
    }}}}
    //VJ use flowboundary map, type 1 is open flow, else use the map
    if (FlowBoundaryType != 1)
        copy(*K2DOutlets, *FlowBoundary);  //copy 1 is 2

}

//---------------------------------------------------------------------------
/**
 * @fn void TWorld::K2DDEMAInitial()
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

void TWorld::K2DDEMAInitial()
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
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::OverlandFlow2D(void)
 * @brief Calls the diffusive wave functions and calculates new discharge, water height and sediment presence
 *
 * Calls the diffusive wave functions and calculates new discharge, water height and sediment presence
 * During this process, surpluss potential infilration is subtracted from the water content.
 * Sediment transport in overland flow is automatically taken into accaunt.
 *
 * @return void
 * @see SwitchKinematic2D
 * @see K1D_METHOD
 * @see K2D_METHOD_INTER
 */
void TWorld::OverlandFlow2D(void)
{
    fill(*tmb, 0.0); // used for average Q during lisem timestep

    //initial function, set to zero
    K2DInit();

    double dt = _dt/2;
    double tof = 0.0;
    //maximum time is the lisem-timestep _dt
    while(tof < _dt-0.001)
    {

        //THIS IS WHERE THE TIMESTEPS ARE SET!
        dt = K2DFlux(tof,_dt);
        //function returns the minimal needed time-step for stable advection (dt > 1.0 for computational speed)
        ThreadPool->SetMask(K2DDEM,K2DDT,K2DDTR,K2DDTC);

        //run the created function on seperate threads
        flowcompute = std::bind((&TWorld::Wrapper_OverlandFlow2D),this,std::placeholders::_1);
        ThreadPool->RunDynamicCompute(flowcompute); //calls Wrapper_OverlandFlow2D
        ThreadPool->WaitForAll();

        FOR_ROW_COL_MV
        {
            K2DDTT->Drc += 0.5 *K2DDT->Drc;
        }

        tof += dt;
        //total time this lisem-timestep
    }
    for(int i = 0 ; i < ThreadPool->Double_Out1.length(); i++)
    {
        K2DQSOut += ThreadPool->Double_Out1.at(i);
        ThreadPool->Double_Out1.replace(i,0.0);
    }

    //VJ new average flux over lisem timestep, else last Qn is used
    FOR_ROW_COL_MV
    {
        WHrunoff->Drc = K2DHNew->Drc;

     //   K2DQ->Drc = tmb->Drc/_dt; //take the timestep average !?

        Qn->Drc = K2DQ->Drc;
        Q->Drc = K2DQ->Drc;
        InfilVolKinWave->Drc = K2DI->Drc; // K2DI is a volume
    }

    correctWH(WHrunoff);
    // correct extreme velocities ad waterheights at edge cells and spreads the surplus water over the entire wet domain

    FOR_ROW_COL_MV
    {
        WHroad->Drc = WHrunoff->Drc;
        // set road to average outflowing wh, no surface storage.

        WH->Drc = WHrunoff->Drc + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        if(K2DSlope->Drc > MIN_SLOPE && K2DPits->Drc != 1)
        {
            if(WHrunoff->Drc > 1e-3)
                V->Drc = Qn->Drc/(WHrunoff->Drc*ChannelAdj->Drc);
            else
                V->Drc = 0;
        }
        else
        {
            V->Drc = 0;
            Qn->Drc = 0;
        }

        WaterVolall->Drc = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        // is the same as :         WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);

    }

    double thv = 10;
    double dv = 5;
    FOR_ROW_COL_MV
    {

        if (V->Drc < thv)
            continue;

        double vs1 = V->Drc;
        double vu = r > 0 && !MV(r-1,c) ? V->data[r-1][c]+dv : vs1;
        double vd = r < _nrRows-1 && !MV(r+1,c) ? V->data[r+1][c]+dv : vs1;
        double vl = c > 0 && !MV(r,c-1) ? V->data[r][c-1]+dv : vs1;
        double vr = c < _nrCols-1 && !MV(r,c+1) ? V->data[r][c+1] + dv :vs1;

        bool fv1 = (vs1 >= vu && vs1 >= vd && vs1 >= vl && vs1 >= vr);

        if (vs1 > thv || fv1) {
            double vh = WHrunoff->Drc/dt;
            double vkin = sqrt(qPow(WHrunoff->Drc, 2.0/3.0)*sqrt(Grad->Drc)/N->Drc);
            V->Drc = std::min(std::min(vh, vkin), vs1);
            Q->Drc = V->Drc * WHrunoff->Drc*ChannelAdj->Drc;
        }
    }
    if(SwitchErosion)
    {
        //calculate concentration and new sediment discharge
        //WHrunoff and Qn are adapted in case of 2D routing
        if(!SwitchUseGrainSizeDistribution)
        {
            FOR_ROW_COL_MV
            {
                Conc->Drc =  MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, &Sed->Drc, &DEP->Drc);
                Qsn->Drc = Conc->Drc * Qn->Drc;
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
                    Conc_D.Drcd = MaxConcentration(WHrunoff->Drc * ChannelAdj->Drc * DX->Drc, &Sed_D.Drcd, &DEP->Drc);
                    Conc->Drc += Conc_D.Drcd;
                }
            }
        }
    }
}

//--------------------------------------------------------------------------------------------
void TWorld::correctWH(cTMap *_WH)
{
    double total = 0;
    double count = 1.0;
    double maxV = 5.0;
    // adjust boundary slopes to avoid extremes
    FOR_ROW_COL_MV {
        if (DomainEdge->Drc > 0 && FlowBoundary->Drc == 0 && V->Drc > maxV)
        {
            double Vavg = getWindowAverage(*V, r, c, false);
            if (V->Drc > maxV && V->Drc > Vavg*10.0)
            {
                double whavg = getWindowAverage(*_WH, r, c, false);
                double whtmp = _WH->Drc;
                _WH->Drc = std::min(_WH->Drc, whavg);
                total += whtmp-_WH->Drc;
            }
        }
        if (_WH->Drc > 0)
            count+=1.0;
    }

    if(fabs(total) > 0)
    {
        //qDebug() << total << total/count;

        FOR_ROW_COL_MV {
            if (_WH->Drc > 0)
                _WH->Drc += total/count;
        }

        FOR_ROW_COL_MV {
            if (_WH->Drc > 0)
            {
                V->Drc = pow(_WH->Drc, 2.0/3.0)*sqrt(Grad->Drc)/N->Drc;
                Qn->Drc = V->Drc*_WH->Drc*FlowWidth->Drc;
            }
        }

    }
}
//--------------------------------------------------------------------------------------------
void TWorld::Wrapper_OverlandFlow2D(int thread)
{
    K2DPreSolve(thread);
    K2DSolvebyInterpolation(thread);
    // bylinear interpolation solution for diffusive

    // sediment transport functions must be called before K2DSolve() and after K2DSolveBy..()
    if(SwitchErosion)
    {
        //K2DQSOut is the boundary outflow that is returned by he K2DSolveBy....Sed() function.

        //advect total sediment
        if(!SwitchUseGrainSizeDistribution)
        {
            ThreadPool->Double_Out1.replace(thread,ThreadPool->Double_Out1.at(thread) + K2DSolvebyInterpolationSed(thread,Sed, Conc));
        }else
        {
            //advect each induvidual grain class
            FOR_GRAIN_CLASSES
            {
                ThreadPool->Double_Out1.replace(thread,ThreadPool->Double_Out1.at(thread) + K2DSolvebyInterpolationSed(thread,Sed_D.at(d), Conc_D.at(d)));
            }
        }
    }

    K2DSolve(thread);
    //subtract infiltration and sets Qn->Drc = K2DQ->Drc; and WHrunoff->Drc = K2DHNew->Drc;


    FOR_ROW_COL_UF2DMT_DT
    {
        tmb->Drc += K2DQ->Drc *K2DDT->Drc;
    }}}}



    K2DDEMA(thread);

}
