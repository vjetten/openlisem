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

// This particular code uses much of the principles of the FullSWOF model
// https://arxiv.org/abs/1204.3210
// https://sourcesup.renater.fr/frs/?group_id=895&release_id=3901#fullswof_2d-_1.10.00-title-content
// the scheme is made suited for parallel processing
// LICENCE: http://cecill.info/licences/Licence_CeCILL_V2-en.html
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

void TWorld::SWOFSediment(double dt, cTMap * h, cTMap *w, cTMap * u,cTMap * v)
{
    // vector velocity for detachment
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        V->Drc = qSqrt(u->Drc*u->Drc + v->Drc*v->Drc);
    }}

    SedimentDetachmentSS(dt, h, w , V, SSFlood, SSCFlood, SSTCFlood, SSDetFlood, DepFlood, SettlingVelocitySS, SUSPflood);
    // suspended detachment (SS), same generic function as for 1D

    if (SwitchPest)
        PesticideFlowDetachmentSS(SSDetFlood, DepFlood, SSFlood);
    // uses the detachment and depositon for pesticide fractions

    if (SwitchUse2Phase) {
        SedimentDetachmentBL(dt, h, w, V);
        // includes SWOFSedimentLayerDepth that splits wh in ss and bl layer
    } else {
        copy(*SSDepthFlood, *h);
    }

    SWOFSedimentAdvection(dt, h,u,v, SSFlood, SSCFlood, SSDepthFlood);
    // susponded matter flow, advection

    if (SwitchIncludeDiffusion) {
        SWOFSedimentDiffusion(dt, h,u,v, SSFlood, SSCFlood);
    }

    SedimentSetConcentration(h, SSFlood, SSCFlood, SSDepthFlood);

    //bedload detachment and movement
    if (SwitchUse2Phase) {
        SWOFSedimentAdvection(dt, h,u,v, BLFlood, BLCFlood, BLDepthFlood);
        SedimentSetConcentration(h, BLFlood, BLCFlood, BLDepthFlood);
    }

    if (SwitchPest) {
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L {
            double vol = CHAdjDX->Drc*h->Drc + MicroStoreVol->Drc;
            PCrw->Drc = vol > 1e-9 ? PMrw->Drc/(vol*1000) : 0.0;
            //note: PCrs is done in flow detahcment
        }}

        SWOFSedimentAdvection(dt, h,u,v, PMrw, PCrw, SSDepthFlood);
        // dissolved pest distribution between cells
        SWOFSedimentAdvection(dt, h,u,v, PMrs, PCrs, SSDepthFlood);
        // absorbed pest distribution between cells
        if (SwitchIncludeDiffusion) {
            SWOFSedimentDiffusion(dt, h,u,v, PMrw, PCrw); //dissolved
            SWOFSedimentDiffusion(dt, h,u,v, PMrs, PCrs); //absorbed
        }

        // calculate new concentration
        #pragma omp parallel for num_threads(userCores)
        FOR_ROW_COL_MV_L{
            double volmw {0.0};         // L - volume of water in mixing layer
            double massms {0.0};        // kg - mass of sediment in mixing layer
            double vol = CHAdjDX->Drc*h->Drc + MicroStoreVol->Drc;
            if (vol > 0.0)
                PCrw->Drc = PMrw->Drc / (vol * 1000);
            else
                PCrw->Drc = 0.0;
            // L = m * m * m * -- * 1000
            volmw = zm->Drc * DX->Drc * SoilWidthDX->Drc * Theta_mix->Drc * 1000;
            PCmw->Drc = PMmw->Drc / volmw; //

            // kg = m * m * m * kg m_3 * --
            massms = zm->Drc * DX->Drc * SoilWidthDX->Drc * rhoPest;
            //mg kg-1 = mg / kg
            PCms->Drc = PMms->Drc / massms;
        }}
    }

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

        double dux = qMax(dux1,dux2);
        double dvy = qMax(dvy1,dvy2);
        double dvx = qMax(dvx1,dvx2);
        double duy = qMax(duy1,duy2);

        //diffusion coefficient according to J.Smagorinski (1964)
        double eddyvs = cdx * cdy * sqrt(dux*dux + dvy*dvy +  0.5 * (dvx +duy)*(dvx +duy));
        double eta = eddyvs/FS_SigmaDiffusion; // set to 0.5

        //cell directions
        int dx[4] = {0, 1, -1, 0};
        int dy[4] = {1, 0, 0, -1};
        double flux[4] = {0.0,0.0,0.0,0.0};

        //use the calculated weights to distribute flow
        for (int i = 0; i < 4; i++)
        {
            //must multiply the cell directions by the sign of the slope vector components
            int rr = r+dy[i];
            int cr = c+dx[i];

            //add fluxes to cells
            if(INSIDE(rr,cr) && !pcr::isMV(LDD->Drcr))
            {
                //diffusion coefficient eta
                double coeff = SSDepthFlood->Drc > 0 ? dt * eta * qMin(1.0, SSDepthFlood->Drcr/SSDepthFlood->Drc) : 0.0;
                coeff = qMin(coeff, courant_factorSed); //????
                flux[i] = coeff*_SS->Drc;
                if (i == 0) tma->Drcr += flux[i];
                if (i == 1) tmb->Drcr += flux[i];
                if (i == 2) tmc->Drcr += flux[i];
                if (i == 3) tmd->Drcr += flux[i];
            }
        }

        _SS->Drc -= (flux[0]+flux[1]+flux[2]+flux[3]);
        // subtract fluxes form the cell`
    }}

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        _SS->Drc = _SS->Drc + tma->Drc + tmb->Drc + tmc->Drc + tmd->Drc;
    }}

}

//--------------------------------------------------------------------------------------------

void TWorld::SWOFSedimentAdvection(double dt, cTMap *h, cTMap *u,cTMap *v,cTMap *_SS, cTMap *_SSC, cTMap *_SSD)
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

        double velocityfactor = qBound(0.0,  _SSD->Drc / h->Drc, 1.0);
        // if only suspended sediment the velocityfactor is 1.0, if 2 phaseflow the bedload is assumed to
        // move with a fraction of the velocity

        //first calculate the weights for the cells that are closest to location that flow is advected to
        double u_ = u->Drc * velocityfactor;
        double v_ = v->Drc * velocityfactor;

        //the sign of the x and y direction of flow
        double yn = signf(v_);
        double xn = signf(u_);

        double velocity = qSqrt(u_*u_ + v_*v_);

        if(velocity > he_ca && h->Drc > he_ca) {
            double dss = dt*velocity * ChannelAdj->Drc*_SSD->Drc * _SSC->Drc;
            // s*m/s*m*m*kg/m3 = kg

            // no real reason to do this? dt is already limited by courant so this would be double

            // if(dss > courant_factorSed * _SS->Drc)
            //    dss = courant_factorSed * _SS->Drc;
            // equals courant factor but not more than 0.2
            // courant_factorSed is the same as general courant factor, but limited to < 0.2

            //should not travel more distance than cell size
            double dsx = xn*qMin(fabs(u_)/velocity,1.0);
            double dsy = yn*qMin(fabs(v_)/velocity,1.0);

            //cell directions
            int dx[4] = {0, 1, 1, 0};
            int dy[4] = {1, 0, 1, 0};

            double w[4] = {0.0,0.0,0.0,0.0};
            for (int i=0; i<4; i++)
            {
                //multiply the cell directions by the sign of the slope vector components
                int rr = r+(int)yn*dy[i];
                int cr = c+(int)xn*dx[i];

                // distance we want is equal to: 1 - distance from the advected location to the neighbouring cell
                double wdx = 1.0 - qFabs( xn * ((double)dx[i]) - dsx);
                double wdy = 1.0 - qFabs( yn * ((double)dy[i]) - dsy);

                //the distribution is inverly proportional to the squared distance
                double weight = qFabs(wdx) * qFabs(wdy);

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
void TWorld::SWOFSedimentSetConcentration(int r, int c, double h, double w)
{
    if(h > he_ca)
    {
        double Area = w * DX->Drc;
        if (SwitchUse2Phase)
            BLCFlood->Drc = MaxConcentration(Area*BLDepthFlood->Drc, BLFlood->Drc);
        SSCFlood->Drc = MaxConcentration(Area*SSDepthFlood->Drc, SSFlood->Drc);
    }
    else
    {
        BLCFlood->Drc = 0;
        SSCFlood->Drc = 0;
    }
}

void TWorld::SedimentSetConcentration(cTMap *h, cTMap *SS_, cTMap *SSC_, cTMap *SSD_)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        if(h->Drc > he_ca) {
            double vol = CHAdjDX->Drc*SSD_->Drc;
            SSC_->Drc = MaxConcentration(vol, SS_->Drc);
        }
        else
            SSC_->Drc = 0;
    }}
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
        SSDepthFlood->Drc = h;
        return;
    }

    double ps = 2650;
    double pw = 1000;
    double factor = 0.5;

    double d50m = (D50->Drc/1000000.0);
    double d90m = (D90->Drc/1000000.0);
    //critical shear velocity for bed level motion by van rijn
    double critshearvel = velocity * sqrt(GRAV)/(18 * log10(4*(_dx * h/(h*2 + _dx))/d90m));
    //critical shear stress for bed level motion by van rijn
    double critsheart = (critshearvel*critshearvel)/ (((ps-pw)/pw) * GRAV*d50m);
    //rough bed bed load layer depth by Hu en Hui
    BLDepthFlood->Drc = qMin(qMin(d50m * 1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), factor*h), 0.1);
    SSDepthFlood->Drc = qMax(h - BLDepthFlood->Drc,0.0);
}
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::SWOFSedimentDet(double dt, int r,int c)
 * @brief Flow detachment and deposition for flood water
 *
 // * Flow detachment and deposition for flood water for a single cell.
 // * Based on the settling velocity of the grain classes and the
 // * transport capacity, erosion and deposition are simulated.
 // * for each grain class induvidually.
 // * Detachment is taken from the upper soil layer when possible ,and the lower
 // * soil layer afterwards. Deposition is added to the upper soil layer.
 // * The sediment concentration can not
 // * reach values above MAXCONC. Concentrations are rescaled to prevent this,
 // * with surplus sediment being deposited.
 *
 * @param dt : timestep to be taken
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

//THIS IS NOW THE GENERIC DETACHMENT USED IN 1D and 2D FLOW
void TWorld::SedimentDetachmentSS(double dt, cTMap *h, cTMap *w, cTMap *v,
                               cTMap *SS_, cTMap *SSC_, cTMap *SSTC_, cTMap *SSDet_, cTMap *Dep_, cTMap *SSVs_, int type)
{
    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double SS = SS_->Drc;

        double sswatervol = 0;

        double velocity = v->Drc;

        double wf = w->Drc;
        double hf = h->Drc;

        SSTC_->Drc = calcTCSuspended(r, c, FS_SS_Method, hf, wf, velocity, type);
        sswatervol = hf*wf*DX->Drc;

        double deposition = 0;
        double detachment = 0;

        if(h->Drc < HMIN) {
            if(DO_SEDDEP == 1) {
                //set all to zero when the water height is zero
                Dep_->Drc += -SS_->Drc; // gives extreme deposition 50000 ton/ha etc.
                SSTC_->Drc = 0;
                SS_->Drc = 0;
                SSC_->Drc = 0;
            }
        } else {
            // there is water

            //first check if sediment goes to suspended sediment layer or to bed layer
            double TransportFactor = 0;

            double maxTC = qMax(SSTC_->Drc - SSC_->Drc, 0.0) ;
            // positive difference: TC deficit becomes detachment (ppositive)
            double minTC = qMin(SSTC_->Drc - SSC_->Drc, 0.0) ;
            // negative diff, becomes deposition

            //deposition based on settling velocity
            if (minTC < 0) {
                //NOTE use entire depth h for deposition of SS
                if (SwitchDepositionLinear)
                    TransportFactor = dt * SSVs_->Drc * wf*DX->Drc;
                else
                    TransportFactor = (1-exp(-dt*SSVs_->Drc/hf)) * sswatervol;

                deposition  = qMax(TransportFactor*minTC, -SS);

                // exceptions
                // if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                //     deposition = 0;
                // prevent any activity on the boundary!

                if (SwitchSedtrap && SedMaxVolume->Drc == 0 && N->Drc == SedTrapN) {
                    N->Drc = Norg->Drc;
                    // if sed trap is full, no effect of N increase
                }
                if (SwitchSedtrap && SedMaxVolume->Drc > 0) {
                    if (SS > 0) {
                        double maxvol = SedMaxVolume->Drc; // decreases from max to zero
                        double depvol = SS/BulkDens; // m3
                        if (depvol > maxvol)
                            depvol = maxvol;
                        if (maxvol > 0){
                            deposition = -depvol*BulkDens;
                            maxTC = 0;
                        }
                        SedMaxVolume->Drc = qMax(0.0, maxvol - depvol);
                        SedimentFilter->Drc += depvol*BulkDens; // TODO, must become an output
                    }
                }

                if(SwitchGridRetention) {
                    //TODO what happens to pesticides
                    if (SS > 0) {
                        double depvol = SS/BulkDens; // sed in m3
                        if (GridRetention->Drc < depvol)
                            depvol = GridRetention->Drc;
                        if (GridRetention->Drc > 0){
                            deposition = -depvol*BulkDens;  // deposition is all that goes into trench
                            maxTC = 0;
                        }
                        GridRetention->Drc = GridRetention->Drc - depvol;
                    }
                }
            } else {
                if (maxTC > 0 && CohesionSoil->Drc >= 0) {
                    TransportFactor = dt * SSVs_->Drc * wf*DX->Drc;
                   // TransportFactor = dt * TSettlingVelocitySS * SoilWidthDX->Drc * DX->Drc;
                    // m3, detachment only erosion on soilwidth

                    detachment = Y->Drc * maxTC * TransportFactor;
                    //check how much of the potential detachment can be detached from soil layer
                    //detachment = DetachMaterial(r,c,1, false, true, false, detachment);

                    // Detachment exceptions:

                    if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                        detachment = 0;
                    // prevent any activity on the boundary!

                    if (GrassFraction->Drc > 0)
                        detachment = (1-GrassFraction->Drc) * detachment;
                    // no flow detachment on grass strips

                    // no flow detachment in sedtraps or gridretention
                    if (SwitchSedtrap && SedMaxVolume->Drc >= 0)
                        detachment = 0;

                    if(SwitchGridRetention && GridRetention->Drc >= 0)
                        detachment = 0;

                    detachment = (1-StoneFraction->Drc) * detachment;
                    // no flow detachment on stony surfaces

                    detachment *= qMin(1.0, qMax(0.0, 1.0 - (RoadWidthHSDX->Drc/_dx)));
                    // no flow detachment on hard surfaces, map is 0 is not selected

                    if (SwitchHouses)
                        detachment = (1-HouseCover->Drc) * detachment;
                    // no flow det where houses
                    if (SwitchSnowmelt)
                        detachment = (1-Snowcover->Drc) * detachment;
                    // TODO: CHECK THIS no flow detachment on snow
                    //is there erosion and sedimentation under the snowdeck?

                    if(SS + detachment > MAXCONC * sswatervol)
                        detachment = qMax(0.0, MAXCONC * sswatervol - SS);
                    // not more detachment then is needed to keep below ssmax
                }
            }
            //### sediment balance
            SSDet_->Drc += detachment;  // set to zero in mass balance
            Dep_->Drc += deposition;
            SS += deposition;
            SS += detachment;
            SS_->Drc = qMax(0.0,SS);
            SSC_->Drc = MaxConcentration(sswatervol, SS_->Drc);

        } // h > MIN_HEIGHT
    }}
}


void TWorld::SedimentDetachmentBL(double dt, cTMap * h, cTMap *w, cTMap * V)
{

    #pragma omp parallel for num_threads(userCores)
    FOR_ROW_COL_MV_L {
        double BL = BLFlood->Drc;

        double blwatervol = 0;

        double velocity = V->Drc;//std::sqrt(u->Drc *u->Drc + v->Drc * v->Drc);

        double wf = w->Drc;
        double hf = h->Drc;

        SWOFSedimentLayerDepth(r,c, hf, velocity);
        //creates BLDepth and SSDepth, or if 1 layer ssdepth = h and bldepth = 0

        //calculate tranport capacity for bed load and suspended load
        // Bedload is based on D90, susp on D50
        BLTCFlood->Drc = calcTCBedload(r, c, FS_BL_Method, hf, wf, velocity, SUSPflood);
        blwatervol = wf*DX->Drc * BLDepthFlood->Drc;

        double deposition = 0;
        double detachment = 0;

        if(h->Drc < HMIN)
        {
            if(DO_SEDDEP == 1) {
                //set all to zero when the water height is zero
                DepFlood->Drc += -BLFlood->Drc;
                BLTCFlood->Drc = 0;
                BLFlood->Drc = 0;
                BLCFlood->Drc = 0;
            }
        } else {
        // there is water

            deposition = 0;
            detachment = 0;
            if(BLDepthFlood->Drc < MIN_HEIGHT) {

                // if the BLdepth is to small dump everything
                DepFlood->Drc += -BLFlood->Drc;
                BLTCFlood->Drc = 0;
                BLFlood->Drc = 0;
                BLCFlood->Drc = 0;

            } else {
                // there is BL transport

                //### calc concentration and net transport capacity
                double maxTC = qMax(BLTCFlood->Drc - BLCFlood->Drc,0.0);
                double minTC = qMin(BLTCFlood->Drc - BLCFlood->Drc,0.0);
                // unit kg/m3
                if (minTC < 0) {
                    // IN KG/CELL

                    //### deposition
                    double TransportFactor = _dt*SettlingVelocityBL->Drc * wf*DX->Drc;
                    // deposition can occur on roads and on soil (so use flowwidth)

                    // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0
                    deposition = qMax(minTC * TransportFactor, -BL);
                    // cannot have more depo than sediment present

                    if (SwitchNoBoundarySed && FlowBoundary->Drc > 0)
                        deposition = 0;
                    // prevent any activity on the boundary!

                    if (SwitchSedtrap && SedMaxVolume->Drc > 0) {
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

                    if(SwitchGridRetention) {
                        if (Sed->Drc > 0) {
                            double depvol = BL/BulkDens; // sed in m3
                            if (GridRetention->Drc < depvol)
                                depvol = GridRetention->Drc;
                            if (GridRetention->Drc > 0){
                                deposition = -depvol*BulkDens;  // deposition is all that goes into trench
                                maxTC = 0;
                            }
                            GridRetention->Drc = GridRetention->Drc - depvol;
                        }
                    }
                } else {
                    if (maxTC > 0 && Y->Drc > 0) {

                        //### detachment ###

                        // detachment can only come from soil, not roads (so do not use flowwidth)
                        // units s * m/s * m * m = m3
                        //TransportFactor = dt * TSettlingVelocityBL * DX->Drc * SoilWidthDX->Drc;
                        double TransportFactor = dt * SettlingVelocityBL->Drc * wf*DX->Drc;
                        //TransportFactor = qMin(TransportFactor, bldischarge * dt);

                        detachment = maxTC * qMin(TransportFactor, blwatervol);
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

                        if (SwitchHouses)
                            detachment = (1-HouseCover->Drc)*detachment;

                        detachment *= qMin(1.0, qMax(0.0, 1.0 - (RoadWidthHSDX->Drc/_dx)));
                        // no flow detachment on hard surfaces, map is 0 is not selected

                        // no flow det from house roofs
                        if (SwitchSnowmelt)
                            detachment = (1-Snowcover->Drc) * detachment;
                        /* TODO: CHECK THIS no flow detachment on snow */
                        //is there erosion and sedimentation under the snowdeck?

                        detachment = qMax(0.0,detachment);

                        //detachment = DetachMaterial(r,c,1,false,false,true, detachment);
                        detachment *= Y->Drc;

                        if(BL + detachment > MAXCONC * blwatervol)
                            detachment = qMax(0.0, MAXCONC * blwatervol - BL);
                        // limit detachment to what BLflood can carry

                        if (SwitchSedtrap && SedMaxVolume->Drc > 0)
                            detachment = 0;

                        if (SwitchGridRetention && GridRetention->Drc > 0)
                            detachment = 0;
                    }
                }
                //### sediment balance IN KG/CELL
                DepFlood->Drc += deposition;
                BLDetFlood->Drc += detachment;
                BL += detachment;
                BL += deposition;
                BLFlood->Drc = qMax(0.0,BL);
            }
        } // h > MIN_HEIGHT
    }}
}
