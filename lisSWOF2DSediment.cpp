

/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file lisSWOF2DSediment.cpp
  \brief Sediment transport for the SWOF2D shallow flood model

functions: \n

- void TWorld::FloodFlowDetachment(void);
- void TWorld::FloodSedFluxReconstruction(void);

 */

#include "model.h"
#include "operation.h"

#define signf(x)  ((x < 0)? -1.0 : 1.0)

#define he_ca 1e-12
#define ve_ca 1e-12

#define dt_ca 0.005

#define GRAV 9.8067
#define EPSILON 1e-6

//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::SWOFSediment(double dt)
 * @brief Sediment for shallow floods
 *
 * This function calls functions for
 * sediment detachment/depositon, transport and diffusion.
 * During this process uses some variables from the flood calculations,
 * and should therefore be called right before the new velocity and water height are set.
 *
 * @param dt : the timestep to be taken, should be the SWOF timestep
 * @param h : the flood water height
 * @param u : the flood velocity in the x-direction
 * @param v : the flood velocity in the y-direction
 *
 * @return void
 *
 * @see SWOFSedimentDet
 * @see SWOFSedimentCheckZero
 * @see SWOFSedimentSetConcentration
 */

void TWorld::SWOFSediment(double dt, cTMap * h,cTMap * u,cTMap * v)
{
    //sediment detachment or deposition
    SWOFSedimentDetNew(dt, h,u,v);

    // susponded matter
    SWOFSedimentFlowInterpolation(dt, h,u,v, SSFlood, SSCFlood);

    // Include Bedload
    if (SwitchUse2Phase)
       SWOFSedimentFlowInterpolation(dt, h,u,v, BLFlood, BLCFlood);

    if (SwitchIncludeDiffusion)
        SWOFSedimentDiffusion(dt, h,u,v, SSFlood, SSCFlood);

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        SWOFSedimentSetConcentration(r,c,h);
    }}
}

//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::SWOFSedimentDiffusion(double dt, cTMap * _SS,cTMap * _SSC)
 * @brief Applies diffusion to bed load and suspended load based
 *
 * Applies diffusion to bed load and suspended load based on concentrations.
 * Based the concentration gradient and the velocity gradients, fluxes are transported to adjecent cells.
 * The concentration map is recalculated after the diffusion.
 *
 * @param dt : timestep to be taken
 * @param h : the flood water height
 * @param u : the flood velocity in the x-direction
 * @param v : the flood velocity in the y-direction
 * @param _SS : Suspended sediment
 * @param _SSC : Suspended sediment concentration
 *
 * @return void
 *
 * @see FS_SigmaDiffusion
 */

void TWorld::SWOFSedimentDiffusion(double dt, cTMap *h,cTMap *u,cTMap *v, cTMap *_SS, cTMap *_SSC)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tma->Drc = 0;
        tmb->Drc = 0;
        tmc->Drc = 0;
        tmd->Drc = 0;
    }}

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        // double courant = this->courant_factorSed;

        //cell sizes
        double cdx = DX->Drc;
        double cdy = _dx;
        //here it is about spacing, not flow width, so use _dx instead of ChannelAdj->Drc

        //mixing coefficient
        //double sigma = 1;
        bool bc1 = c > 0 && !MV(r,c-1)        ;
        bool bc2 = c < _nrCols-1 && !MV(r,c+1);
        bool br1 = r > 0 && !MV(r-1,c)        ;
        bool br2 = r < _nrRows-1 && !MV(r+1,c);

        double dux1 = bc1 ? std::abs(u->Drc - u->data[r][c-1]) : 0;
        double dvy1 = br1 ? std::abs(v->Drc - v->data[r-1][c]) : 0;
        double dvx1 = bc1 ? std::abs(v->Drc - v->data[r][c-1]) : 0;
        double duy1 = br1 ? std::abs(u->Drc - u->data[r-1][c]) : 0;
        double dux2 = bc2 ? std::abs(u->data[r][c+1] - u->Drc) : 0;
        double dvy2 = br2 ? std::abs(v->data[r+1][c] - v->Drc) : 0;
        double dvx2 = bc2 ? std::abs(v->data[r][c+1] - v->Drc) : 0;
        double duy2 = br2 ? std::abs(u->data[r+1][c] - u->Drc) : 0;

        double dux = std::max(dux1,dux2);
        double dvy = std::max(dvy1,dvy2);
        double dvx = std::max(dvx1,dvx2);
        double duy = std::max(duy1,duy2);

        //diffusion coefficient according to J.Smagorinski (1964)
        double eddyvs = cdx * cdy * sqrt(dux*dux + dvy*dvy +  0.5 * (dvx +duy)*(dvx +duy));
        double eta = eddyvs/FS_SigmaDiffusion;

        //cell directions
        int dx[4] = {0, 1, -1, 0};
        int dy[4] = {1, 0, 0, -1};
        double flux[4] = {0.0,0.0,0.0,0.0};

        //use the calculated weights to distribute flow
        for (int i=0; i<4; i++)
        {
            //must multiply the cell directions by the sign of the slope vector components
            int rr = r+dy[i];
            int cr = c+dx[i];

            //add fluxes to cells
            if(INSIDE(rr,cr) && !pcr::isMV(LDD->Drcr))
            {
                //diffusion coefficient eta
                double coeff = SSDepthFlood->Drc > 0 ? dt * eta * std::min(1.0, SSDepthFlood->Drcr/SSDepthFlood->Drc) : 0.0;
                coeff = std::min(coeff, courant_factorSed);
                flux[i] = coeff*_SS->Drc;
                if (i == 0) tma->Drcr += flux[i];
                if (i == 1) tmb->Drcr += flux[i];
                if (i == 2) tmc->Drcr += flux[i];
                if (i == 3) tmd->Drcr += flux[i];
            }
        }

        _SS->Drc -= (flux[0]+flux[1]+flux[2]+flux[3]);
    }}


    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        _SS->Drc = _SS->Drc + tma->Drc + tmb->Drc+tmc->Drc+tmd->Drc;
    }}

}

//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::SWOFSedimentFlowInterpolation(int r, int c, double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
 * @brief Distributes sediment flow in flood water trough bilinear interpolation.
 *
 * Distributes sediment flow in flood water trough bilinear interpolation.
 * The concentration map, together with water height and velocity
 * are used to determine sediment discharges.
 *
 * @param dt : timestep to be taken
 * @param h : the flood water height
 * @param u : the flood velocity in the x-direction
 * @param v : the flood velocity in the y-direction
 * @param _BL : Bed load sediment
 * @param _BLC : Bed load sediment
 * @param _SS : Suspended sediment
 * @param _SSC : Suspended sediment concentration
 *
 * @return void
 */

void TWorld::SWOFSedimentFlowInterpolation(double dt, cTMap *h, cTMap *u,cTMap *v,cTMap *_SS, cTMap *_SSC)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        tma->Drc = 0;
        tmb->Drc = 0;
        tmc->Drc = 0;
        tmd->Drc = 0;
    }}

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        //first calculate the weights for the cells that are closest to location that flow is advected to
        double u_ = u->Drc;
        double v_ = v->Drc;

        //the sign of the x and y direction of flow
        double yn = signf(v_);
        double xn = signf(u_);

        double velocity = sqrt(u_*u_ + v_*v_);

        if(velocity > he_ca && h->Drc > he_ca) {
           // double courant = this->courant_factorSed; // why *0.1

            double dss = dt*velocity*ChannelAdj->Drc *SSDepthFlood->Drc * _SSC->Drc;
            // s*m/s*m*m*kg/m3 = kg
            if(dss > courant_factorSed * _SS->Drc)
                dss = courant_factorSed *  _SS->Drc;

            //should not travel more distance than cell size
            double dsx = xn*std::min(fabs(u_)/velocity,1.0);
            double dsy = yn*std::min(fabs(v_)/velocity,1.0);

            //cell directions
            int dx[4] = {0, 1, 1, 0};
            int dy[4] = {1, 0, 1, 0};

            double w[4] = {0.0,0.0,0.0,0.0};
            for (int i=0; i<4; i++)
            {
                //must multiply the cell directions by the sign of the slope vector components
                int rr = r+(int)yn*dy[i];
                int cr = c+(int)xn*dx[i];

                // distance we want is equal to: 1 - distance from the advected location to the neighbouring cell
                double wdx = ((double)1.0) - fabs( xn * ((double)dx[i]) - dsx);
                double wdy = ((double)1.0) - fabs( yn * ((double)dy[i]) - dsy);

                //the distribution is inverly proportional to the squared distance
                double weight = fabs(wdx) * fabs(wdy);

                if(INSIDE(rr,cr) && !pcr::isMV(LDD->Drcr)) {
                    if(h->Drcr > he_ca) {
                        w[i] = weight;
                    }
                }
            }

            //normalize: sum of the 4 weights is equal to 1
            double wt = w[0];
            wt += w[1];
            wt += w[2];
            wt += w[3];

            if(wt == 0) {
                w[3] = 1.0;
                wt = 1.0;
            }

            w[0] = w[0]/wt;
            w[1] = w[1]/wt;
            w[2] = w[2]/wt;
            w[3] = w[3]/wt;

            double flux[4] = {0.0,0.0,0.0,0.0};

            for (int i=0; i<4; i++) {

                int rr = r+(int)yn*dy[i];
                int cr = c+(int)xn*dx[i];
                if(INSIDE(rr,cr) && !pcr::isMV(LDD->Drcr))
                {
                    if(h->Drcr > he_ca)
                    {
                        flux[i] = w[i]*dss;

                        if (i == 0) tma->Drcr += flux[i];
                        if (i == 1) tmb->Drcr += flux[i];
                        if (i == 2) tmc->Drcr += flux[i];
                        if (i == 3) tmd->Drcr += flux[i];
                    }
                }
            }
            /*
            int rr;
            int cr;
            rr = r+(int)yn;
            cr = c;
            if(INSIDE(rr,cr) && !pcr::isMV(LDD->Drcr))
            {
                flux[0] = w[0]*dss;
                tma->Drcr += flux[0];
            }
            rr = r;
            cr = c+(int)xn;
            if(INSIDE(rr,cr) && !pcr::isMV(LDD->Drcr))
            {
                flux[1] = w[1]*dss;
                tmb->Drcr += flux[1];
            }
            rr = r+(int)yn;
            cr = c+(int)xn;
            if(INSIDE(rr,cr) && !pcr::isMV(LDD->Drcr))
            {
                flux[2] = w[2]*dss;
                tmc->Drcr += flux[2];
            }
            //rr = r;
            //cr = c;
            flux[3] = w[3]*dss;
            tmc->Drc += flux[3];
*/
            // subtract the four fluxes from each cell
            _SS->Drc -= (flux[0]+flux[1]+flux[2]+flux[3]); // flux is in kg!
        } // v en h > ha
    }}

    // update SS with new values in 4 cells that have changed
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        _SS->Drc = _SS->Drc + tma->Drc + tmb->Drc + tmc->Drc + tmd->Drc;
    }}

}
//--------------------------------------------------------------------------------------------
        // OBSOLETE
/**
 * @fn void TWorld::SWOFSedimentSetConcentration(int r, int c)
 * @brief Deposits all sediment when the water height is zero
 *
 * @param r : the row nr of the cell
 * @param c : the column nr of the cell
 * @param h : the flood water height
 * @param u : the flood velocity in the x-direction
 * @param v : the flood velocity in the y-direction
 *
 * @return void
 *
 */
void TWorld::SWOFSedimentCheckZero(int r, int c, cTMap * h)//,cTMap * u,cTMap * v)
{
    cTMap * _BL = BLFlood;
    cTMap * _BLC = BLCFlood;
    cTMap * _SS = SSFlood;
    cTMap * _SSC = SSCFlood;

    if(h->Drc < he_ca)
    {
        //add sediment to deposition
        if (SwitchUse2Phase)
            DepFlood->Drc += -(_BL->Drc);
        DepFlood->Drc += -(_SS->Drc);

        //add to soil layer
        if(SwitchUseMaterialDepth) {
            StorageDep->Drc += _BL->Drc + _SS->Drc;
        }

        //set to zero
        _BL->Drc = 0;
        _SS->Drc = 0;
        _BLC->Drc = 0;
        _SSC->Drc = 0;
    }

}


//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::SWOFSedimentSetConcentration(int r, int c)
 * @brief Calculates concentration of sediment in a cell based on sediment in flow and water volume
 *
 * @param r : the row nr of the cell
 * @param c : the column nr of the cell
 * @param h : the flood water height
 * @param u : the flood velocity in the x-direction
 * @param v : the flood velocity in the y-direction
 *
 * @return void
 *
 * @see MaxConcentration
 */
void TWorld::SWOFSedimentSetConcentration(int r, int c, cTMap * h)
{
    if(h->Drc > he_ca)
    {
        if (SwitchUse2Phase)
            BLCFlood->Drc = MaxConcentration(CHAdjDX->Drc*BLDepthFlood->Drc, &BLFlood->Drc, &DepFlood->Drc);
        SSCFlood->Drc = MaxConcentration(CHAdjDX->Drc*SSDepthFlood->Drc, &SSFlood->Drc, &DepFlood->Drc);
    }
    else
    {
        BLCFlood->Drc = 0;
        SSCFlood->Drc = 0;
    }
}


//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::SWOFSedimentLayerDepth(double dt, int r,int c)
 * @brief Sets the depth of the bed layer and suspended layer for the indicated cell
 *
 * @param r : row nr of cell
 * @param c : column nr of cell
 * @param h : the flood water height
 * @param u : the flood velocity in the x-direction
 * @param v : the flood velocity in the y-direction
 *
 * @return void
 */
void TWorld::SWOFSedimentLayerDepth(int r , int c, double h, double velocity)
{
    if (!SwitchUse2Phase) {
        BLDepthFlood->Drc = 0;
        SSDepthFlood->Drc = h;
        return;
    }

    double ps = 2650;
    double pw = 1000;
    double factor = 0.5;

    double d50m = (D50->Drc/1000000.0);
    double d90m = (D90->Drc/1000000.0);
    //critical shear velocity for bed level motion by van rijn
    double critshearvel = velocity * sqrt(GRAV)/(18 * log10(4*(ChannelAdj->Drc * h/(h*2 + ChannelAdj->Drc))/d90m));
    //critical shear stress for bed level motion by van rijn
    double critsheart = (critshearvel*critshearvel)/ (((ps-pw)/pw) * GRAV*d50m);
    //rough bed bed load layer depth by Hu en Hui
    BLDepthFlood->Drc = std::min(std::min(d50m * 1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), factor*h), 0.1);
    SSDepthFlood->Drc = std::max(h - BLDepthFlood->Drc,0.0);
}
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::SWOFSedimentDet(double dt, int r,int c)
 * @brief Flow detachment and deposition for flood water
 *
 * Flow detachment and deposition for flood water for a single cell.
 * Based on the settling velocity of the grain classes and the
 * transport capacity, erosion and deposition are simulated.
 * for each grain class induvidually.
 * Detachment is taken from the upper soil layer when possible ,and the lower
 * soil layer afterwards. Deposition is added to the upper soil layer.
 * The sediment concentration can not
 * reach values above MAXCONC. Concentrations are rescaled to prevent this,
 * with surplus sediment being deposited.
 *
 * @param dt : timestep to be taken
 * @param r : row nr of cell
 * @param c : column nr of cell
 * @param h : the flood water height
 * @param u : the flood velocity in the x-direction
 * @param v : the flood velocity in the y-direction
 *
 * @return void
 *
 * @see SWOFSedimentLayerDepth
 * @see SWOFSedimentTCSS
 * @see SWOFSedimentTCBL
 * @see DetachMaterial
 */


void TWorld::SWOFSedimentDetNew(double dt, cTMap * h,cTMap * u,cTMap * v)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {

    double chadj = ChannelAdj->Drc;
    double chadjdx = CHAdjDX->Drc;
    double BLTC = 0;
    double SSTC = 0;
    double BLDepth = BLDepthFlood->Drc;
    double SSDepth = SSDepthFlood->Drc;
    double SSC = SSCFlood->Drc;
    double BLC = BLCFlood->Drc;
    double SS = SSFlood->Drc;
    double BL = BLFlood->Drc;
    double TSettlingVelocitySS = SettlingVelocitySS->Drc;
    double TSettlingVelocityBL = SettlingVelocityBL->Drc;
    double bldischarge = 0;
    double blwatervol = 0;
    double ssdischarge = 0;
    double sswatervol = 0;

    double velocity = std::sqrt(u->Drc *u->Drc + v->Drc * v->Drc);
    SWOFSedimentLayerDepth(r,c,h->Drc, velocity);
    //creates BLDepth and SSDepth, or if 1 layer ssdepth = h and bldepth = 0

    //calculate tranport capacity for bed load and suspended load
    // Bedload is based on D90, susp on D50
    if (SwitchUse2Phase) {
       BLTC = calcTCBedload(r, c, 1, FS_BL_Method, h->Drc, velocity, 1);
       BLTCFlood->Drc = BLTC;
       bldischarge = velocity * chadj * BLDepth;
       blwatervol = chadjdx * BLDepth;
    }

    SSTC = calcTCSuspended(r, c, 1, FS_SS_Method, h->Drc, velocity, 1);
    SSTCFlood->Drc = SSTC;
    ssdischarge = velocity * chadj * SSDepth;
    sswatervol = chadjdx * SSDepth;

    double deposition = 0;
    double detachment = 0;

    if(h->Drc < MIN_HEIGHT)
    {
        //set all to zero when the water height is zero
        if (SwitchUse2Phase) {
            DepFlood->Drc += -BLFlood->Drc;
            BLTCFlood->Drc = 0;
            BLFlood->Drc = 0;
            BLCFlood->Drc = 0;
        }

        DepFlood->Drc += -SSFlood->Drc;
        SSTCFlood->Drc = 0;
        SSFlood->Drc = 0;
        SSCFlood->Drc = 0;

        if(SwitchUseMaterialDepth) {
            StorageDep->Drc += -deposition;
        }

    } else {
        // there is water

        //####### DO SUSPENDED FIRST

        //first check if sediment goes to suspended sediment layer or to bed layer
        double TransportFactor;

        double maxTC = std::max(SSTC - SSC, 0.0) ;
        // positive difference: TC deficit becomes detachment (ppositive)
        double minTC = std::min(SSTC - SSC, 0.0) ;
        // negative diff, becomes deposition

        deposition = 0;
        detachment = 0;

        //deposition based on settling velocity
        if (minTC < 0) {
            //NOTE use entire depth h for deposition of SS
            TransportFactor = (1-exp(-dt*TSettlingVelocitySS/h->Drc)) * sswatervol;
          //  TransportFactor = dt * TSettlingVelocitySS * DX->Drc * SoilWidthDX->Drc;
            TransportFactor = std::min(TransportFactor, ssdischarge * dt);

            deposition  = std::max(TransportFactor * minTC,-SS);

            // exceptions
            if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                deposition = 0;
            // prevent any activity on the boundary!

            if (SwitchSedtrap && SedMaxVolume->Drc == 0 && N->Drc == SedTrapN) {
                N->Drc = Norg->Drc;
                // if sed trap is full, no effect of N increase
            }
            if (SwitchSedtrap && SedMaxVolume->Drc > 0)
            {
                if (SS > 0) {
                    double maxvol = SedMaxVolume->Drc;
                    double depvol = SS * 1.0/BulkDens; // m3
                    if (maxvol < depvol)
                        depvol = maxvol;
                    if (maxvol > 0){
                        deposition = -depvol*BulkDens;
                        maxTC = 0;
                    }
                    SedMaxVolume->Drc = maxvol - depvol;
                    SedimentFilter->Drc += depvol*BulkDens;
                }
            }

            if(SwitchUseMaterialDepth)
                StorageDep->Drc += -deposition;

        } else
            if (Y->Drc > 0) {
                TransportFactor = dt * TSettlingVelocitySS * DX->Drc * ChannelAdj->Drc;//SoilWidthDX->Drc;
                // only erosion on soilwidth

                detachment = maxTC * std::min(TransportFactor, sswatervol);

                // Detachment exceptions:

                if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                    detachment = 0;
                // prevent any activity on the boundary!

                if (GrassFraction->Drc > 0)
                    detachment = (1-GrassFraction->Drc) * detachment;
                // no flow detachment on grass strips

                //detachment = (1-StoneFraction->Drc) * detachment;
                // no flow detachment on stony surfaces

                if (SwitchRoadsystem)
                    detachment = (1-RoadWidthDX->Drc/_dx) * detachment;
                // no flow detachment on hard surfaces

                if (SwitchHardsurface)
                    detachment = (1-HardSurface->Drc) * detachment;
                // no flow detachment on hard surfaces

                if (SwitchHouses)
                    detachment = (1-HouseCover->Drc)*detachment;
                // no flow det where houses

                detachment = (1-Snowcover->Drc) * detachment;
                // TODO: CHECK THIS no flow detachment on snow
                //is there erosion and sedimentation under the snowdeck?

                if (SwitchSedtrap && SedMaxVolume->Drc > 0)
                    detachment = 0;
                //if there is still room in the sed trap then no detahcment on those cells

                //check how much of the potential detachment can be detached from soil layer
                detachment = DetachMaterial(r,c,1, false, true, false, detachment);
                //bool channel, bool flood,bool bl

               // if(MAXCONC * sswatervol < SS + detachment)
                 //   detachment = 0;//std::max(0.0, MAXCONC * sswatervol - SS);
                // not more detachment then is needed to keep below ssmax
            }

        //### sediment balance
        SSDetFlood->Drc += detachment;  // set to zero in mass balance
        DepFlood->Drc += deposition;
        SS += deposition;
        SS += detachment;
        SSFlood->Drc = std::max(0.0,SS);

        // ########### DO BEDLOAD
        if (SwitchUse2Phase) {
            deposition = 0;
            detachment = 0;
            if(BLDepth < MIN_HEIGHT) {

                // if the BLdepth is to small dump everything
                DepFlood->Drc += -BLFlood->Drc;
                BLTCFlood->Drc = 0;
                //BLDetFlood->Drc = 0;
                BLFlood->Drc = 0;
                BLCFlood->Drc = 0;

            } else {
                // there is BL transport

                //### calc concentration and net transport capacity
                maxTC = std::max(BLTC - BLC,0.0);
                minTC = std::min(BLTC - BLC,0.0);
                // unit kg/m3
                if (minTC < 0) {
                    // IN KG/CELL

                    //### deposition
                    TransportFactor = (1-exp(-dt*TSettlingVelocityBL/BLDepth)) * blwatervol;
                    //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
                    // deposition can occur on roads and on soil (so use flowwidth)

                    // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0
                    deposition = std::max(minTC * TransportFactor, -BL);
                    // cannot have more depo than sediment present

                    if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                        deposition = 0;
                    // prevent any activity on the boundary!

                    if (SwitchSedtrap && SedMaxVolume->Drc > 0)
                    {
                        if (BL > 0) {
                            double maxvol = SedMaxVolume->Drc;
                            double depvol = BL * 1.0/BulkDens; // m3
                            if (maxvol < depvol)
                                depvol = maxvol;
                            if (maxvol > 0){
                                deposition = -depvol*BulkDens;
                                maxTC = 0;
                            }
                            SedMaxVolume->Drc = maxvol - depvol;
                            SedimentFilter->Drc += depvol*BulkDens;
                        }
                    }

                    if(SwitchUseMaterialDepth)
                    {
                        StorageDep->Drc += -deposition;
                    }
                } else
                if (Y->Drc > 0) {
                    //### detachment
                    // detachment can only come from soil, not roads (so do not use flowwidth)
                    // units s * m/s * m * m = m3
                    TransportFactor = dt * TSettlingVelocityBL * DX->Drc * SoilWidthDX->Drc;
                    //TransportFactor = std::min(TransportFactor, bldischarge * dt);
                    //TransportFactor = bldischarge*DT->Drc;

                    detachment = maxTC * std::min(TransportFactor, blwatervol);
                    // unit = kg/m3 * m3 = kg

                    if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                        detachment = 0;
                    // VJ 190325 prevent any activity on the boundary!

                    if (GrassFraction->Drc > 0)
                        detachment = (1-GrassFraction->Drc) * detachment;
                    // no flow detachment on grass strips

                    // Detachment edxceptions:
                    detachment = (1-StoneFraction->Drc) * detachment;
                    // no flow detachment on stony surfaces

                    if (SwitchHardsurface)
                        detachment = (1-HardSurface->Drc) * detachment;
                    // no flow detachment on hard surfaces

                    if (SwitchHouses)
                        detachment = (1-HouseCover->Drc)*detachment;
                    // no flow det from house roofs

                    detachment = (1-Snowcover->Drc) * detachment;
                    /* TODO: CHECK THIS no flow detachment on snow */
                    //is there erosion and sedimentation under the snowdeck?

                    detachment = DetachMaterial(r,c,1,false,false,true, detachment);

                    if(MAXCONC * blwatervol < BL + detachment)
                        detachment = std::max(0.0, MAXCONC * blwatervol - BL);
                    // limit detachment to what BLflood can carry

                    if (SwitchSedtrap && SedMaxVolume->Drc > 0) {
                        detachment = 0;
                    }
                }
                //### sediment balance IN KG/CELL
                DepFlood->Drc += deposition;
                BLDetFlood->Drc += detachment;
                BL += detachment;
                BL += deposition;
                BLFlood->Drc = std::max(0.0,BL);
            }  // BL exist
        } // 2 phase       
    } // h > MIN_HEIGHT

    SWOFSedimentSetConcentration(r,c,h);

    }}
}
