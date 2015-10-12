

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
  \brief kinematic wave routing functions and calculation of discharge and sed flux per cell.

  The routing functions use local variables to allow for possible re-usage


 */

#include "model.h"
#include "operation.h"


//--------------------------------------------------------------------------------------------
/**
 * @brief TWorld::K2DInit: initializes all necessary variables for the 2d kinematic wave solution
 * @return void
 */
void TWorld::K2DInit()
{
    //reset all maps for calculations
    FOR_ROW_COL_MV
    {
        if(SwitchErosion)
        {
            K2DQS->Drc = 0;
            K2DQSX->Drc = 0;
            K2DQSY->Drc = 0;
            K2DS->Drc = Sed->Drc;
            K2DSC->Drc = Conc->Drc;
            K2DSCN->Drc = 0;
            K2DSFX->Drc = 0;
            K2DSFY->Drc = 0;
        }
        if(SwitchPesticide)
        {
            K2DQP->Drc = 0;
            K2DQPX->Drc = 0;
            K2DQPY->Drc = 0;
            K2DP->Drc = Pest->Drc;
            K2DPC->Drc = C->Drc;
            K2DPCN->Drc = 0;
            K2DPFX->Drc = 0;
            K2DPFY->Drc = 0;
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

    FOR_ROW_COL_MV
    {
        double dy = FlowWidth->Drc;

        double hrunoff = std::max(WHrunoff->Drc ,0.0);

        double Perim = 2.0*hrunoff+dy;

        if (Perim > 0)
            R->Drc = dy*hrunoff/Perim;
        else
            R->Drc = 0;

        Alpha->Drc = pow(N->Drc/sqrt(K2DSlope->Drc) * pow(Perim, (2.0/3.0)),0.6);
        //WHY aplha k2dlope and not grad?

        if (Alpha->Drc == 0)
            K2DQ->Drc = 0;
        else
            K2DQ->Drc = pow((dy*hrunoff)/Alpha->Drc, 1.0/0.6);

        Qn->Drc = K2DQ->Drc;

        WHrunoff->Drc = (Alpha->Drc*pow(Qn->Drc, 0.6))/ChannelAdj->Drc;
        K2DHOld->Drc = WHrunoff->Drc;

    }



}
//--------------------------------------------------------------------------------------------
/**
 * @brief TWorld::K2DFlux: Calculates discharges for each cell
 * @param dtmax : a hint for the calculation of the timestep
 * @return dt : the minimal needed timestep to ensure stability
 */
double TWorld::K2DFlux(double dtmax)
{
    double dtr = dtmax;
    double fraction = 0.5;
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
        double Perim = 2.0*hrunoff+cdy;
        if (Perim > 0)
            R->Drc = hrunoff*cdy/Perim;
        else
            R->Drc = 0;

        Alpha->Drc = pow(N->Drc/sqrt(K2DSlope->Drc) * pow(Perim, (2.0/3.0)),0.6);

        if(Alpha->Drc > 0)
            K2DQ->Drc = pow((cdy*hrunoff)/Alpha->Drc, 1.0/0.6);
        else
            K2DQ->Drc = 0;


        //K2DQ->Drc = std::min(K2DQ->Drc,fraction*cdx*cdy*hrunoff/dtr );

        //set sediment and pesticide transport
        if(SwitchErosion)
        {
            K2DQS->Drc =  K2DQ->Drc * K2DSC->Drc;
        }
        if(SwitchPesticide)
        {
            K2DQP->Drc =  K2DQ->Drc * K2DPC->Drc;
        }

        //within this timestep, only half of the cells available water should flow out
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

    dtr = std::max(dtr,1.0);
    FOR_ROW_COL_MV
    {
        //limit discharge to half of the cells water
        if(fraction * (DX->Drc*K2DHOld->Drc*ChannelAdj->Drc) < K2DQ->Drc*dtr)
        {
            K2DQ->Drc =fraction* (DX->Drc*K2DHOld->Drc*ChannelAdj->Drc)/dtr;
        }
    }


    //return the lowest needed timestep, with a minimum of 1.0 seconds
    return dtr;

}
//--------------------------------------------------------------------------------------------
/**
 * @brief TWorld::K2DSolvebyInterpolation: Solves kinematic wave in 2d using bilinear interpolation of advected flow
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

        if(SwitchErosion)
            K2DSCN->Drc = K2DSC->Drc;
        if(SwitchPesticide)
            K2DPCN->Drc = K2DPC->Drc;
    }

    //first calculate the weights for the cells that are closest to location that flow is advected to
    FOR_ROW_COL_MV
    {
        if(K2DPits->Drc == 1 || (K2DSlopeX->Drc == 0 && K2DSlopeY->Drc == 0))
        {
            continue;
        }

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
        int end = 2;
        int start = 0;

        if(K2DPitsD->Drc == 1)
        {
            w[2] = 1.0;
        }else
        {
            //for each cell niegbhouring the advected location of the discharge, calculate interpolation weight
            for (int i=0; i<4; i++)
            {
                // distance we want is equal to: 1 - distance from the advected location to the neighbouring cell
                double wdx = ((double)1.0) - fabs( xn * ((double)dx[i]) - dsx);
                double wdy = ((double)1.0) - fabs( yn * ((double)dy[i]) - dsy);

                //the distribution is inverly proportional to the squared distance
                double weight = fabs(wdx) *fabs(wdy);

                w[i] = weight;

            }
            //normalize: sum of the 4 weights is equal to 1
            double wt = 0.0;
            for (int i=start; i<end+1; i++)
            {
                wt += w[i];
            }
            for (int i=start; i<end+1; i++)
            {
                w[i] = w[i]/wt;
            }


        }


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
                QinKW->data[r2][c2] += w[i] *dt*K2DQ->Drc;
                QoutKW->data[r2][c2] += w[i] *dt*K2DQ->Drc;

            }
            else
                if(K2DOutlets->Drc == 1)
                {
                    K2DQOut +=  w[i]*dt*(K2DQ->Drc);
                    QoutKW->data[r][c] += w[i] *dt*K2DQ->Drc;
                    K2DHNew->data[r][c] -=  w[i]*dt*(K2DQ->Drc/(cdx*cdy));
                }

        }
    }

    //similar process as above for sediment and pesticides, sediment is distributed along with discharge using interpolation
    //new water heights are needed for this routine
    if(SwitchErosion || SwitchPesticide)
    {
        FOR_ROW_COL_MV
        {
            if(K2DPits->Drc == 1 || (K2DSlopeX->Drc == 0 && K2DSlopeY->Drc == 0))
            {
                continue;
            }

            double DHL = sqrt(K2DSlopeX->Drc*K2DSlopeX->Drc + K2DSlopeY->Drc* K2DSlopeY->Drc);
            double dsx = K2DSlopeX->Drc/DHL;
            double dsy = K2DSlopeY->Drc/DHL;
            double yn = dsy/fabs(dsy);
            double xn = dsx/fabs(dsx);

            if(dsx == 0){xn = 1.0;};
            if(dsy == 0){yn = 1.0;};

            double cdx = DX->Drc;
            double cdy = ChannelAdj->Drc;

            //cell directions
            int dx[4] = {0, 1, 1, 0};
            int dy[4] = {1, 0, 1, 0};

            double w[4] = {0.0,0.0,0.0,0.0};

            //for each cell niegbhouring the advected location of the discharge, calculate interpolation weight
            for (int i=0; i<3; i++)
            {
                int r2, c2;

                //must multiply the cell directions by the sign of the slope vector components
                r2 = r+yn*dy[i];
                c2 = c+xn*dx[i];

                // distance we want is equal to: 1 - distance from the advected location to the neighbouring cell
                double wdx = ((double)1.0) - fabs( xn * ((double)dx[i]) - dsx);
                double wdy = ((double)1.0) - fabs( yn * ((double)dy[i]) - dsy);

                //the distribution is inverly proportional to the squared distance
                double weight = fabs(wdx) *fabs(wdy);

                w[i] = weight;

                //if the cell that flows needs to go to is out of bounds or missing value, skip

                if(INSIDE(r2,c2) && !pcr::isMV(LDD->data[r2][c2]))
                {
                    double cdx2 = DX->data[r2][c2];
                    double cdy2 = ChannelAdj->data[r2][c2];

                    if(SwitchErosion)
                    {
                        K2DS->data[r][c] -=w[i]*dt*(K2DQS->Drc);
                        K2DS->data[r2][c2]+=w[i]*dt*(K2DQS->Drc);

                    }
                    if(SwitchPesticide)
                    {
                        K2DP->data[r][c] -=w[i]*dt*(K2DQP->Drc);
                        K2DP->data[r2][c2]+=w[i]*dt*(K2DQP->Drc);
                    }

                }else if(K2DOutlets->Drc == 1)
                {
                    if(SwitchErosion)
                    {
                        K2DS->data[r][c] -=w[i]*dt*(K2DQS->Drc);
                        K2DQSOut +=w[i]*dt*(K2DQS->Drc);

                    }
                    if(SwitchPesticide)
                    {
                        K2DP->data[r][c] -=w[i]*dt*(K2DQP->Drc);
                        K2DQPOut +=w[i]*dt*(K2DQS->Drc);
                    }
                }
            }
        }
    }
}
//--------------------------------------------------------------------------------------------
/**
 * @brief TWorld::K2DSolvebyFlux: Solves kinematic wave in 2d using a cell boundary flux method
 * @param dt : timestep
 * @return void
 */
void TWorld::K2DSolvebyFlux(double dt)
{
    FOR_ROW_COL_MV
    {
        //start with old height and concentration
        K2DHNew->Drc = K2DHOld->Drc;
        K2DQN->Drc = 0;

        if(SwitchErosion)
            K2DSCN->Drc = K2DSC->Drc;
        if(SwitchPesticide)
            K2DPCN->Drc = K2DPC->Drc;
    }

    fill(*K2DQX, 0.0);
    fill(*K2DQY, 0.0);
    FOR_ROW_COL_MV
    {
        if(K2DPits->Drc == 1 || (K2DSlopeX->Drc == 0 && K2DSlopeY->Drc == 0))
        {
            K2DQX->Drc = 0;
            K2DQY->Drc = 0;
            continue;
        }

        double slopexy = K2DSlopeX->Drc != 0 ? K2DSlopeY->Drc/K2DSlopeX->Drc : 0;
        double slopeyx = K2DSlopeY->Drc != 0 ? K2DSlopeX->Drc/K2DSlopeY->Drc : 0;
        double powslopexy_025 = sqrt(sqrt(1.0 + slopexy*slopexy));
        double powslopeyx_025 = sqrt(sqrt(1.0 + slopeyx*slopeyx));

        //the weights for the x, any y component of the flow
        //The sqrt(x^2 + y^2) can not be used for components of discharge, since then Qx + Qy = Qtotal would not hold!
        //this distribution of flow between the x and y components is based on G Tayfur (2001)
        double xw = (K2DSlopeX->Drc > 0? 1.0:-1.0)*sqrt(fabs(K2DSlopeX->Drc))/powslopexy_025;
        double yw = (K2DSlopeY->Drc > 0? 1.0:-1.0)*sqrt(fabs(K2DSlopeY->Drc))/powslopeyx_025;

        //if the slope in a direction is 0, then set the weight to 0, to correct for devisions by 0
        //Lim_{Sx ->0, Sy ->1} Wx(Sx,Sy) = 0 || Lim_{Sx ->0, Sy ->1} Wy(Sx,Sy) = 1|| Lim_{Sy ->0, Sx ->1} Wx(Sx,Sy) = 1 || Lim_{Sy ->0, Sx ->1} Wy(Sx,Sy) = 0
        if(K2DSlopeX->Drc == 0){xw = 0.0;yw = 1.0;};
        if(K2DSlopeY->Drc == 0){yw = 0.0;xw = 1.0;};

        //apply weights to components in x and y direction
        K2DQX->Drc = K2DQ->Drc*xw;
        K2DQY->Drc = K2DQ->Drc*yw;

        //similar for sediment and pesticide transport
        if(SwitchErosion)
        {
            K2DQSX->Drc = K2DQS->Drc*xw;
            K2DQSY->Drc = K2DQS->Drc*yw;
        }
        if(SwitchPesticide)
        {
            K2DQPX->Drc = K2DQP->Drc*xw;
            K2DQPY->Drc = K2DQP->Drc*yw;
        }

    }

    //reset maps for later calculations
    fill(*K2DFX, 0.0);
    fill(*K2DFY, 0.0);
    if(SwitchErosion)
    {
        fill(*K2DSFX, 0.0);
        fill(*K2DSFY, 0.0);
    }
    if(SwitchPesticide)
    {
        fill(*K2DPFX, 0.0);
        fill(*K2DPFY, 0.0);
    }

    //now calculate the sum of ingoing and outgoing discharges for each cell, in both the x and y direction
    FOR_ROW_COL_MV
    {


        if(r != 0)
        {
            if(!pcr::isMV(LDD->data[r-1][c]) && !(K2DPitsD->data[r-1][c] == 1))
            {
                //add discharge trough cell border to cell, subtract it from the source
                double fin = std::max(K2DQY->data[r-1][c],0.0);
                K2DFY->Drc += fin;
                K2DFY->data[r-1][c] -= fin;
                QinKW->Drc += fin;
                QoutKW->data[r-1][c] += fin;
                if(SwitchErosion)
                {
                    K2DSFY->Drc += std::max(K2DQSY->data[r-1][c],0.0);
                    K2DSFY->data[r -1][c] -= std::max(K2DQSY->data[r-1][c],0.0);
                }
                if(SwitchPesticide)
                {
                    K2DPFY->Drc += std::max(K2DQPY->data[r-1][c],0.0);
                    K2DPFY->data[r -1][c] -= std::max(K2DQPY->data[r-1][c],0.0);
                }
            }

        }

        if(r != _nrRows-1)
        {
            if(!pcr::isMV(LDD->data[r+1][c]) && !(K2DPitsD->data[r+1][c] == 1))
            {
                double fin = fabs(std::min(K2DQY->data[r+1][c],0.0));
                K2DFY->Drc += fin;
                K2DFY->data[r+1][c] -= fin;
                QinKW->Drc += fin;
                QoutKW->data[r+1][c] += fin;

                if(SwitchErosion)
                {
                    K2DSFY->Drc += fabs(std::min(K2DQSY->data[r+1][c],0.0));
                    K2DSFY->data[r+1][c] -= fabs(std::min(K2DQSY->data[r+1][c],0.0));
                }
                if(SwitchPesticide)
                {
                    K2DPFY->Drc += fabs(std::min(K2DQPY->data[r+1][c],0.0));
                    K2DPFY->data[r+1][c] -= fabs(std::min(K2DQPY->data[r+1][c],0.0));
                }
            }

        }
        if(c != 0 )
        {
            if(!pcr::isMV(LDD->data[r][c-1]) && !(K2DPitsD->data[r][c-1] == 1))
            {
                double fin = std::max(K2DQX->data[r][c-1],0.0);
                K2DFX->Drc +=  fin;
                K2DFX->data[r][c-1] -= fin;
                QinKW->Drc += fin;
                QoutKW->data[r][c-1] += fin;

                if(SwitchErosion)
                {
                    K2DSFY->Drc += std::max(K2DQSY->data[r][c-1],0.0);
                    K2DSFY->data[r][c-1] -= std::max(K2DQSY->data[r][c-1],0.0);
                }
                if(SwitchPesticide)
                {
                    K2DPFY->Drc += std::max(K2DQPY->data[r][c-1],0.0);
                    K2DPFY->data[r][c-1] -= std::max(K2DQPY->data[r][c-1],0.0);
                }
            }

        }

        if(c != _nrCols-1 )
        {
            if(!pcr::isMV(LDD->data[r][c+1]) && !(K2DPitsD->data[r][c+1] == 1))
            {
                double fin = fabs(std::min(K2DQX->data[r][c+1],0.0));
                K2DFX->Drc +=  fin;
                K2DFX->data[r][c+1] -= fin;
                QinKW->Drc += fin;
                QoutKW->data[r][c+1] += fin;

                if(SwitchErosion)
                {
                    K2DSFY->Drc += fabs(std::min(K2DQSY->data[r][c+1],0.0));
                    K2DSFY->data[r][c+1] -= fabs(std::min(K2DQSY->data[r][c+1],0.0));
                }
                if(SwitchPesticide)
                {
                    K2DPFY->Drc += fabs(std::min(K2DQPY->data[r][c+1],0.0));
                    K2DPFY->data[r][c+1] -= fabs(std::min(K2DQPY->data[r][c+1],0.0));
                }
            }
        }

        //add outflow from outlets to total outflow
        if(K2DOutlets->Drc == 1)
        {
            //calculate normalized direction of flow
            double dsx = K2DSlopeX->Drc;
            double dsy = K2DSlopeY->Drc;
            int r2 = r + (dsx > 0? 1: -1);
            int c2 = c + (dsy > 0? 1: -1);
            //is the cell in this direction either out of bounds, or missing value?
            if(!INSIDE(r2,c) || pcr::isMV(LDD->data[r2][c]) )
            {
                //then add the flow to outflow, and subtract from cell
                K2DFX->Drc -= fabs(K2DQX->data[r][c]);
                K2DQOut += dt* fabs(K2DQX->data[r][c]);
                QoutKW->Drc += dt* fabs(K2DQX->data[r][c]);

                if(SwitchErosion)
                {
                    K2DSFX->Drc -= fabs(K2DQSX->data[r][c]);
                    K2DQSOut += dt* fabs(K2DQSX->data[r][c]);
                }if(SwitchPesticide)
                {
                    K2DPFX->Drc -= fabs(K2DQPX->data[r][c]);
                    K2DQPOut += dt* fabs(K2DQPX->data[r][c]);
                }

            }
            //is the cell in this direction either out of bounds, or missing value?
            if(!INSIDE(r,c2) || pcr::isMV(LDD->data[r][c2]) )
            {
                //then add the flow to outflow, and subtract from cell
                K2DFY->Drc -= fabs(K2DQY->data[r][c]);
                K2DQOut += dt* fabs(K2DQY->data[r][c]);
                QoutKW->Drc += dt* fabs(K2DQY->data[r][c]);

                if(SwitchErosion)
                {
                    K2DSFX->Drc -= fabs(K2DQSY->data[r][c]);
                    K2DQSOut += dt* fabs(K2DQSY->data[r][c]);
                }if(SwitchPesticide)
                {
                    K2DPFY->Drc -= fabs(K2DQPY->data[r][c]);
                    K2DQPOut += dt* fabs(K2DQPY->data[r][c]);
                }
            }
        }

        //handle special cases were only diagonal flow is presen (Flow distance longer??discharge smaller?)
        if(K2DPitsD->Drc == 1)
        {
            double dsx = K2DSlopeX->Drc;
            double dsy = K2DSlopeY->Drc;
            int r2 = r + (dsx > 0? 1: -1);
            int c2 = c + (dsy > 0? 1: -1);

            if(INSIDE(r2,c2))
            {
                if(!pcr::isMV(LDD->data[r2][c2]))
                {
                    K2DFX->Drc -= K2DQX->Drc;
                    K2DFY->Drc -= K2DQY->Drc;
                    K2DFX->data[r2][c2] += K2DQX->Drc;
                    K2DFY->data[r2][c2] += K2DQY->Drc;
                }

            }

        }
    }

    FOR_ROW_COL_MV
    {
        //length of cell borders
        double cdx = DX->Drc;
        double cdy = ChannelAdj->Drc;
        //new cell height, using total flux at boundaries
        K2DHNew->Drc = K2DHOld->Drc +  dt*(K2DFX->Drc + K2DFY->Drc)/(cdy*cdx);

        if(SwitchErosion)
        {
            K2DS->Drc += dt*(K2DSFX->Drc + K2DSFY->Drc);
        }
        if(SwitchPesticide)
        {
            K2DP->Drc += dt*(K2DPFX->Drc + K2DPFY->Drc);
        }
    }
}
//---------------------------------------------------------------------------
/**
 * @brief TWorld::K2DSolve: Finalizes the solution for the kinematic wave, sets variables for lisem and recalculates discharge
 * @param dt : timestep
 * @return void
 */
void TWorld::K2DSolve(double dt)
{

    //finish by substracting infiltration, and calculating discharge from new water height
    FOR_ROW_COL_MV
    {
        double cdx = DX->Drc;
        double cdy = ChannelAdj->Drc;

        //calculate infiltartion in time step
        double infil = std::min(FSurplus->Drc *SoilWidthDX->Drc*DX->Drc * dt/_dt,0.0);
        if(K2DHNew->Drc < fabs(infil)/(cdx*cdy))
        {
            infil = -K2DHNew->Drc*(cdx*cdy);

        }
        //keep track of infiltration
        K2DI->Drc -= (infil);
        K2DHNew->Drc = std::max(K2DHNew->Drc + infil/(cdx*cdy) ,0.0);

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
        if(SwitchErosion)
        {

            K2DSCN->Drc = MaxConcentration(WHrunoff->Drc*ChannelAdj->Drc *DX->Drc, K2DS->Drc);
            K2DSC->Drc = K2DSCN->Drc;
            Qsn->Drc = K2DSCN->Drc * Qn->Drc;
            Qs->Drc = Qsn->Drc;
            Sed->Drc = K2DS->Drc;
        }
        if(SwitchPesticide)
        {
            K2DPCN->Drc = MaxConcentration(WHrunoff->Drc*ChannelAdj->Drc *DX->Drc, K2DP->Drc);
            K2DPC->Drc = K2DPCN->Drc;
            Qpn->Drc = K2DPCN->Drc * Qn->Drc;
            Pest->Drc = K2DP->Drc;
        }


    }

}
//---------------------------------------------------------------------------
/**
 * @brief TWorld::K2DCalcVelDisch: Calculation of velocity and discharge of overland flow, for use with the 2d solution
 * @param dt : timestep
 * @return void
 */
void TWorld::K2DCalcVelDisch()
{
    FOR_ROW_COL_MV
    {
        if(K2DPits->Drc == 1)
        {
            Q->Drc = 0;
            V->Drc = 0;
            continue;
        }

        double hrunoff = std::max(WHrunoff->Drc - K2DWHStore->Drc,0.0);

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
        Perim = 2.0*hrunoff+FlowWidth->Drc;

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


        //tm->Drc = V->Drc * R->Drc/kinvisc;
        //Reynolds number
    }

}
//---------------------------------------------------------------------------
/**
 * @brief TWorld::K2DSolve: calculate slope based on Dem height + runoff height
 * average between left and right slope is used
 * for the domain boundary's and the map boundary's, the one-sided slope is used.
 * When both the sides of a cell along an axis are missing values, the slope in that direction is set to 0
 * when a cell is located next to a missing value, but would flow into that cell, it is set as an outlet
 * @param dt : timestep
 * @return void
 */
void TWorld::K2DDEMA()
{



    FOR_ROW_COL_MV
    {
        //   K2DDEM->Drc = DEM->Drc
        K2DOutlets->Drc = 0;
        K2DWHStore->Drc = 0;
        K2DPits->Drc = 0;
        K2DPitsD->Drc = 0;
        K2DSlopeX->Drc = 0;
        K2DSlopeY->Drc = 0;
        K2DSlope->Drc = 0;
        //K2DAspect->Drc = 0;
        //set heigt to dem + water heigt (this provides mannings flow equation with the gradient of water head)
        K2DDEM->Drc = DEM->Drc + WHrunoff->Drc;
    }

    FOR_ROW_COL_MV
    {

        double Dhx = 0;
        double Dhy = 0;

        //DEM

        if(r != 0 && c != 0 && c != _nrCols-1 && r != _nrRows-1)
        {

            //if(!pcr::isMV(K2DDEM->data[r+1][c]) && !pcr::isMV(K2DDEM->data[r-1][c]) && !pcr::isMV(K2DDEM->data[r][c+1]) && !pcr::isMV(K2DDEM->data[r][c-1]))
            {
                double dem = K2DDEM->data[r][c];

                double demx1 = K2DDEM->data[r][c+1];
                double demx2 = K2DDEM->data[r][c-1];

                double demy1 = K2DDEM->data[r+1][c];
                double demy2 = K2DDEM->data[r-1][c];

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



            }

            //one sided slope
            /*if(!pcr::isMV(K2DDEM->data[r+1][c]))
                Dhy = -(K2DDEM->data[r+1][c]-K2DDEM->data[r][c]);
            else
                Dhy = 0;

            if(!pcr::isMV(K2DDEM->data[r][c+1]))
                Dhx = -(K2DDEM->data[r][c+1]-K2DDEM->data[r][c]);
            else
                Dhx = 0;*/

            if(pcr::isMV(K2DDEM->data[r-1][c]) && !pcr::isMV(K2DDEM->data[r+1][c]))
            {
                Dhy = -(K2DDEM->data[r+1][c]-K2DDEM->data[r][c]);
                if( Dhy < 0)
                    Outlet->Drc= 1;
            }
            if(pcr::isMV(K2DDEM->data[r+1][c]) && !pcr::isMV(K2DDEM->data[r-1][c]))
            {
                Dhy = -(K2DDEM->data[r][c]-K2DDEM->data[r-1][c]);
                if( Dhy > 0)
                    Outlet->Drc= 1;
            }
            if(pcr::isMV(K2DDEM->data[r][c-1]) && !pcr::isMV(K2DDEM->data[r][c+1]))
            {
                Dhx = -(K2DDEM->data[r][c+1]-K2DDEM->data[r][c]);
                if( Dhx < 0)
                    Outlet->Drc= 1;
            }
            if(pcr::isMV(K2DDEM->data[r][c+1]) && !pcr::isMV(K2DDEM->data[r][c-1]))
            {
                Dhx = -(K2DDEM->data[r][c]-K2DDEM->data[r][c-1]);
                if( Dhx > 0)
                    Outlet->Drc= 1;
            }
            if(pcr::isMV(K2DDEM->data[r][c+1]) && pcr::isMV(K2DDEM->data[r][c-1]))
            {
                Dhx = 0;
                    Outlet->Drc= 1;
            }
            if(pcr::isMV(K2DDEM->data[r+1][c]) && pcr::isMV(K2DDEM->data[r-1][c]))
            {
                Dhy = 0;
                    Outlet->Drc= 1;
            }
        }

        //at boundaries, use one-sided derivative of the DEM
        if(r == 0)
        {
            Dhy = -(K2DDEM->data[r+1][c]-K2DDEM->data[r][c]);

            if(pcr::isMV(K2DDEM->data[r+1][c]))
            {
                Dhy = 0;
            }
            if(pcr::isMV(K2DDEM->data[r][c-1]))
            {
                Dhx = -(K2DDEM->data[r][c+1]-K2DDEM->data[r][c]);
                if( Dhx < 0)
                    Outlet->Drc= 1;
            }
            if(pcr::isMV(K2DDEM->data[r][c+1]))
            {
                Dhx = -(K2DDEM->data[r][c]-K2DDEM->data[r][c-1]);
                if( Dhx > 0)
                    Outlet->Drc= 1;
            }
            if(pcr::isMV(K2DDEM->data[r][c+1]) && pcr::isMV(K2DDEM->data[r][c-1]))
            {
                Dhx = 0;
            }
        }

        if(r == _nrRows-1)
        {
            Dhy = -(K2DDEM->data[r][c]-K2DDEM->data[r-1][c]);

            if(pcr::isMV(K2DDEM->data[r -1][c]))
            {
                Dhy = 0;
            }

            if(pcr::isMV(K2DDEM->data[r][c-1]))
            {
                Dhx = -(K2DDEM->data[r][c +1]-K2DDEM->data[r][c]);
                if( Dhx < 0) { Outlet->Drc= 1;}
            }
            if(pcr::isMV(K2DDEM->data[r][c+1]))
            {
                Dhx = -(K2DDEM->data[r][c]-K2DDEM->data[r][c-1]);
                if( Dhx > 0) { Outlet->Drc= 1;}
            }
            if(pcr::isMV(K2DDEM->data[r][c+1]) && pcr::isMV(K2DDEM->data[r ][c-1]))
            {
                Dhx = 0;
            }
        }

        if(c == 0)
        {
            Dhx = -(K2DDEM->data[r][c +1]-K2DDEM->data[r][c]);
            if(pcr::isMV(K2DDEM->data[r -1][c]))
            {
                Dhy = -(K2DDEM->data[r][c+1]-K2DDEM->data[r][c]);
                if( Dhy < 0) { Outlet->Drc= 1;}
            }
            if(pcr::isMV(K2DDEM->data[r +1][c]))
            {
                Dhy = -(K2DDEM->data[r][c]-K2DDEM->data[r-1][c]);
                if( Dhy > 0) { Outlet->Drc= 1;}
            }
            if(pcr::isMV(K2DDEM->data[r][c+1]))
            {
                Dhx = 0;
            }
            if(pcr::isMV(K2DDEM->data[r +1][c]) && pcr::isMV(K2DDEM->data[r -1][c]))
            {
                Dhy = 0;
            }
        }

        if(c == _nrCols-1)
        {
            Dhx = -(K2DDEM->data[r][c]-K2DDEM->data[r][c-1]);
            if(pcr::isMV(K2DDEM->data[r -1][c]))
            {
                Dhy = -(K2DDEM->data[r][c+1]-K2DDEM->data[r][c]);
                if( Dhy < 0) { Outlet->Drc= 1;}
            }
            if(pcr::isMV(K2DDEM->data[r +1][c]))
            {
                Dhy = -(K2DDEM->data[r][c]-K2DDEM->data[r-1][c]);
                if( Dhy > 0) { Outlet->Drc= 1;}
            }
            if(pcr::isMV(K2DDEM->data[r][c-1]))
            {
                Dhx = 0;
            }
            if(pcr::isMV(K2DDEM->data[r +1][c]) && pcr::isMV(K2DDEM->data[r -1][c]))
            {
                Dhy = 0;
            }
        }

        //   K2DDX->Drc = -sqrt(Dhx*Dhx + _dx*_dx);
        //   K2DDY->Drc = -sqrt(Dhy*Dhy + _dx*_dx);

        K2DSlopeX->Drc = Dhx/_dx;
        K2DSlopeY->Drc = Dhy/_dx;

        //if the slope is extremely flat, the ldd direction could be used for flow calculations (not needed)
        double DHL = sqrt(K2DSlopeX->Drc*K2DSlopeX->Drc + K2DSlopeY->Drc* K2DSlopeY->Drc);

        //ldd directions
        int dxldd[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
        int dyldd[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

        //if slope is really flat, follow ldd direction
        /*if(DHL < 0.0001 || K2DOutlets->Drc == 1)
        {
            K2DSlopeX->Drc = double(dxldd[(int)LDD->Drc]);
            K2DSlopeY->Drc = double(dyldd[(int)LDD->Drc]);
            DHL = sqrt(K2DSlopeX->Drc*K2DSlopeX->Drc + K2DSlopeY->Drc* K2DSlopeY->Drc);
        }*/

        //Angle of direction of steepest slope, compared to positive x-axis !not used!
        //K2DAspect->Drc = atan2(Dhy,Dhx);

        //calculate actual combined slope, with a minimum value of 0.01
        double Dh = fabs(Dhx) + fabs(Dhy);
        K2DSlope->Drc = Dh / sqrt(2*_dx*_dx);

        if(std::isnan(K2DSlope->Drc))
        {
            K2DSlope->Drc = 0.01;
            K2DSlopeX->Drc = 0.00;
            K2DSlopeY->Drc = 0.00;

        }
        //minimum value for the slope
        K2DSlope->Drc = std::max(K2DSlope->Drc,0.01);

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
        double dem = DEM->Drc;
        double demw = K2DDEM->Drc;
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
                    if(DEM->data[r2][c2] <  lowestneighbor)
                    {
                        lowestneighbor = DEM->data[r2][c2];
                    }

                    //if at least 1 neighboring cell is lower, it is not a pit
                    if(K2DDEM->data[r2][c2] < demw)
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
                    if(K2DDEM->data[r2][c2] <  lowestneighborw)
                    {
                        lowestneighborw = K2DDEM->data[r2][c2];
                    }
                }else
                {
                    mv++;
                }
            }
        }

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
            //qDebug() << r << c << "PIT!!!!";
        }
    }

}

