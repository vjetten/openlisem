

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
 * @fn void TWorld::SWOFSedimentFlowInterpolation(double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC)
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

void TWorld::SWOFSedimentFlowInterpolation(int thread, cTMap *DT, cTMap *h, cTMap *u,cTMap *v,
                                           cTMap *_BL, cTMap *_BLC, cTMap *_SS, cTMap *_SSC)
{

    FOR_ROW_COL_UF2DMTDER {
        if(FloodHMaskDer->Drc != 0)
        {
            //M... only used in this function!
            MBLNFlood->Drc = _BL->Drc;
            MSSNFlood->Drc = _SS->Drc;

            //set concentration from present sediment
            MBLCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*h->Drc, MBLNFlood->Drc);

            //set concentration from present sediment
            MSSCFlood->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*h->Drc, MSSNFlood->Drc);

           //USING h HERE and NOT ssdepth bldepth

        }
    }}}}

    double courant = this->courant_factor;
    // flooding courant factor

    //first calculate the weights for the cells that are closest to location that flow is advected to
    FOR_ROW_COL_UF2DMT_DT {
        if(FloodHMaskDer->Drc != 0)
        {
            //no flood velocity means no flood sediment transport, so skip this cell
            if((v->Drc == 0 && u->Drc == 0))
            {
                continue;
            }

            //the sign of the x and y direction of flow
            double yn = signf(v->Drc);
            double xn = signf(u->Drc);

            double vel = sqrt(u->Drc*u->Drc + v->Drc*v->Drc);

            if(vel < he_ca || h->Drc < he_ca)
            {
                continue;
            }

            double qbl = DT->Drc*vel*ChannelAdj->Drc *BLDepthFlood->Drc * MBLCFlood->Drc;

            if(qbl > courant * MBLNFlood->Drc)
            {
                qbl =  courant * MBLNFlood->Drc;
            }

            double qss = DT->Drc*vel*ChannelAdj->Drc *SSDepthFlood->Drc * MSSCFlood->Drc;

            if(qss > courant * MSSNFlood->Drc)
            {
                qss = courant *  MSSNFlood->Drc;
            }

            //should not travel more distance than cell size
            double dsx = xn*std::min(fabs(u->Drc)/vel,1.0);
            double dsy = yn*std::min(fabs(v->Drc)/vel,1.0);

            //cell directions
            int dx[4] = {0, 1, 1, 0};
            int dy[4] = {1, 0, 1, 0};

            //weights to be saved
            double w[4] = {0.0,0.0,0.0,0.0};

            //for each cell niegbhouring the advected location of the discharge, calculate interpolation weight
            for (int i=0; i<4; i++)
            {
                int r2, c2;

                //must multiply the cell directions by the sign of the slope vector components
                r2 = r+(int)yn*dy[i];
                c2 = c+(int)xn*dx[i];

                // distance we want is equal to: 1 - distance from the advected location to the neighbouring cell
                double wdx = ((double)1.0) - fabs( xn * ((double)dx[i]) - dsx);
                double wdy = ((double)1.0) - fabs( yn * ((double)dy[i]) - dsy);

                //the distribution is inverly proportional to the squared distance
                double weight = fabs(wdx) * fabs(wdy);

                if(INSIDE(r2,c2))
                {
                    if( !pcr::isMV(LDD->data[r2][c2]) && h->data[r2][c2] > he_ca)//0)
                    {
                        w[i] = weight;
                    }
                }
            }

            //normalize: sum of the 4 weights is equal to 1
            double wt = 0.0;
            for (int i=0; i<4; i++)
            {
                wt += w[i];
            }
            if(wt == 0)
            {
                w[3] = 1.0; wt = 1.0;
            }
            for (int i=0; i<4; i++)
            {
                w[i] = w[i]/wt;
            }

            //use the calculated weights to distribute flow
            for (int i=0; i<4; i++)
            {
                int r2, c2;

                //must multiply the cell directions by the sign of the slope vector components
                r2 = r+(int)yn*dy[i];
                c2 = c+(int)xn*dx[i];

                if(INSIDE(r2,c2) && !pcr::isMV(LDD->data[r2][c2]))
                {

                    if(h->data[r2][c2] > he_ca)
                    {
                        //weight * the flow is distributed to the ith cell that neighbours the advected flow.
                        MBLNFlood->data[r2][c2] +=  w[i]* qbl;
                        MBLNFlood->data[r][c] -=  w[i]* qbl;

                        MSSNFlood->data[r2][c2] +=  w[i]* qss;
                        MSSNFlood->data[r][c] -=  w[i]* qss;
                    }
                }
            }
        }
    }}}}

    FOR_ROW_COL_UF2DMTDER {
        if(FloodHMaskDer->Drc != 0)
        {
            MBLNFlood->Drc = std::max(0.0,MBLNFlood->Drc);
            MSSNFlood->Drc = std::max(0.0,MSSNFlood->Drc);
            _BL->Drc = MBLNFlood->Drc;
            _SS->Drc = MSSNFlood->Drc;

            //set concentration from present sediment
            _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*h->Drc, _BL->Drc);

            //set concentration from present sediment
            _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*h->Drc, _SS->Drc);

            //USING h HERE and NOT ssdepth bldepth


        }
    }}}}
}
//--------------------------------------------------------------------------------------------
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

    if(!(h->Drc > he_ca))//0))
    {

        if(!SwitchUseGrainSizeDistribution)
        {
            //add sediment to deposition
            DepFlood->Drc += -(_BL->Drc);
            DepFlood->Drc += -(_SS->Drc);

            //add to soil layer
            if(SwitchUseMaterialDepth)
            {
                StorageDep->Drc += _BL->Drc + _SS->Drc;
            }

            //set to zero
            _BL->Drc = 0;
            _SS->Drc = 0;


            //set concentration from present sediment REDICULOUS, SS = 0 so SSC = 0
            _BLC->Drc = 0;//MaxConcentration(ChannelAdj->Drc*DX->Drc*h->Drc, _BL->Drc);

            //set concentration from present sediment
            _SSC->Drc = 0;//MaxConcentration(ChannelAdj->Drc*DX->Drc*h->Drc, _SS->Drc);
        }else
        {
            //set totals to zero
            _BL->Drc = 0;
            _SS->Drc = 0;
            _BLC->Drc = 0;
            _SSC->Drc = 0;

            //add all to deposition
            FOR_GRAIN_CLASSES
            {
                DepFlood->Drc += -(BL_D.Drcd);
                DepFlood->Drc += -(SS_D.Drcd);

                //add to soil layer
                if(SwitchUseMaterialDepth)
                {
                    StorageDep->Drc += BL_D.Drcd + SS_D.Drcd;
                    if(SwitchUseGrainSizeDistribution)
                    {
                        StorageDep_D.Drcd += BL_D.Drcd + SS_D.Drcd;
                    }
                }

                //set to zero
                BL_D.Drcd = 0;
                SS_D.Drcd = 0;
                BLC_D.Drcd = 0;
                SSC_D.Drcd = 0;
            }
        }
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

    cTMap * _BL = BLFlood;
    cTMap * _BLC = BLCFlood;
    cTMap * _SS = SSFlood;
    cTMap * _SSC = SSCFlood;
    double hh = SSDepthFlood->Drc;
    double hhb = BLDepthFlood->Drc;
            //h->Drc;
    if((h->Drc > he_ca))
    {
        if(!SwitchUseGrainSizeDistribution)
        {
            double maxsed = MAXCONCBL * ChannelAdj->Drc*DX->Drc*hhb;
            if (_BL->Drc > maxsed){
                DepFlood->Drc -= (_BL->Drc - maxsed);
                _BL->Drc = maxsed;
            }
            _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hhb, _BL->Drc);

            maxsed = MAXCONC * ChannelAdj->Drc*DX->Drc*hh;
            if (_SS->Drc > maxsed){
                DepFlood->Drc -= (_SS->Drc - maxsed);
                _SS->Drc = maxsed;
            }
            _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*hh, _SS->Drc);
        }else
        {
            FOR_GRAIN_CLASSES
            {
                //set concentration from present sediment
                BLC_D.Drcd = MaxConcentration(ChannelAdj->Drc*DX->Drc*hhb, BL_D.Drcd);

                //set concentration from present sediment
                SSC_D.Drcd = MaxConcentration(ChannelAdj->Drc*DX->Drc*hh, SS_D.Drcd);
            }
        }
    }
    else
    {
        if(!SwitchUseGrainSizeDistribution)
        {
            //set concentration from present sediment
            _BLC->Drc = 0;

            //set concentration from present sediment
            _SSC->Drc = 0;

        }else
        {
            FOR_GRAIN_CLASSES
            {
                //set concentration from present sediment
                BLC_D.Drcd = 0;

                //set concentration from present sediment
                SSC_D.Drcd = 0;
            }
        }
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
void TWorld::SWOFSedimentDiffusion(int thread, cTMap *DT, cTMap *h,cTMap *u,cTMap *v, cTMap *_SS, cTMap *_SSC)
{
    FOR_ROW_COL_UF2DMTDER {
        if(FloodHMaskDer->Drc != 0)
        {
            MSSNFlood->Drc = _SS->Drc;
            //set concentration from present sediment
            MSSCFlood->Drc = _SSC->Drc;
                    //MaxConcentration(ChannelAdj->Drc*DX->Drc*h->Drc, MSSNFlood->Drc);
        }
    }}}}


    //diffusion of Suspended Sediment layer
    FOR_ROW_COL_UF2DMT_DT {
        if(FloodHMaskDer->Drc != 0)
        {

            //cell sizes
            double cdx = DX->Drc;
            double cdy = _dx;
            //here it is about spacing, not flow width, so use _dx instead of ChannelAdj->Drc

            //mixing coefficient
            double sigma = 1;
            double dux1 = c  > 0 ? std::abs(u->data[r][c] - u->data[r][c-1]) : 0;
            double dvy1 = r  > 0 ? std::abs(v->data[r][c] - v->data[r-1][c]) : 0;
            double dvx1 = c  > 0 ? std::abs(v->data[r][c] - v->data[r][c-1]) : 0;
            double duy1 = r  > 0 ? std::abs(u->data[r][c] - u->data[r-1][c]) : 0;
            double dux2 = c  < _nrCols-1 ? std::abs(u->data[r][c+1] - u->data[r][c]) : 0;
            double dvy2 = r  < _nrRows-1 ? std::abs(v->data[r+1][c] - v->data[r][c]) : 0;
            double dvx2 = c  < _nrCols-1 ? std::abs(v->data[r][c+1] - v->data[r][c]) : 0;
            double duy2 = r  < _nrRows-1 ? std::abs(u->data[r+1][c] - u->data[r][c]) : 0;

            //drop term if at boundary
            if(std::isnan(dux1))
                dux1 = 0;
            if(std::isnan(dvx1))
                dvx1 = 0;
            if(std::isnan(duy1))
                duy1 = 0;
            if(std::isnan(dvy1))
                dvy1 = 0;

            if(std::isnan(dux2))
                dux2 = 0;
            if(std::isnan(dvx2))
                dvx2 = 0;
            if(std::isnan(duy2))
                duy2 = 0;
            if(std::isnan(dvy2))
                dvy2 = 0;

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


            //use the calculated weights to distribute flow
            for (int i=0; i<4; i++)
            {
                int r2, c2;

                //must multiply the cell directions by the sign of the slope vector components
                r2 = r+dy[i];
                c2 = c+dx[i];

                //add fluxes to cells
                if(INSIDE(r2,c2) && !pcr::isMV(LDD->data[r2][c2]))
                {
                    //diffusion coefficient
                    double coeff =  DT->Drc * eta * std::min(1.0, SSDepthFlood->data[r2][c2]/SSDepthFlood->data[r][c]);
                    coeff = std::min(coeff, courant_factor/4.0);
//                    std::min(1.0, SSDepthFlood->data[r2][c2]/SSDepthFlood->data[r][c]), courant_factor_diffusive/4.0);

                    MSSNFlood->data[r2][c2] += coeff*MSSNFlood->Drc;
                    MSSNFlood->data[r][c] -= coeff*MSSNFlood->Drc;
                }
            }
        }
    }}}}

    FOR_ROW_COL_UF2DMTDER {
        if(FloodHMaskDer->Drc != 0)
        {
            _SS->Drc = std::max(0.0,MSSNFlood->Drc);
            //set concentration from present sediment
            _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*h->Drc, _SS->Drc);
        }
    }}}}
}
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::SWOFSedimentMaxC(int r, int c)
 * @brief Corrects flood sediment concentration for maximum possible concentration
 *
 * Corrects flood sediment concentration for maximum possible concentration.
 * Concentrations for induvidual grain classas are scaled to the corrected total concentration
 *
 * @param r : row nr of cell
 * @param c : column nr of cell
 * @param h : the flood water height
 * @param u : the flood velocity in the x-direction
 * @param v : the flood velocity in the y-direction
 *
 * @return void
 *
 * @see MAXCONC
 */
void TWorld::SWOFSedimentMaxC(int r, int c)//, cTMap * h,cTMap * u,cTMap * v)
{

    cTMap * _BL = BLFlood;
    cTMap * _BLC = BLCFlood;
    cTMap * _SS = SSFlood;
    cTMap * _SSC = SSCFlood;

    //maximum concentraion
    if(!SwitchUseGrainSizeDistribution)
    {
        _SSC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*SSDepthFlood->Drc, _SS->Drc);
        // limit concentration to 850 and throw rest in deposition

        double sssmax = MAXCONC * DX->Drc *ChannelAdj->Drc*SSDepthFlood->Drc;
        if(sssmax < _SS->Drc)
        {
            DepFlood->Drc += -(_SS->Drc - sssmax);
            if(SwitchUseMaterialDepth)
            {
                StorageDep->Drc += (_SS->Drc - sssmax);
            }
            _SS->Drc = sssmax;
        }

        //set concentration from present sediment
        _BLC->Drc = MaxConcentration(ChannelAdj->Drc*DX->Drc*BLDepthFlood->Drc, _BL->Drc);

        double smax = MAXCONC * DX->Drc *ChannelAdj->Drc*BLDepthFlood->Drc;
        if(smax < _BL->Drc)
        {
            DepFlood->Drc += -(_BL->Drc - smax);
            if(SwitchUseMaterialDepth)
            {
                StorageDep->Drc += (_BL->Drc - smax);
            }
            _BL->Drc = smax;
        }
    }else
    {

        BLFlood->Drc = 0;
        SSFlood->Drc = 0;

        FOR_GRAIN_CLASSES
        {
            BLFlood->Drc += BL_D.Drcd;
            SSFlood->Drc += SS_D.Drcd;
        }

        FOR_GRAIN_CLASSES
        {
            SSC_D.Drcd = MaxConcentration(ChannelAdj->Drc*DX->Drc*SSD_D.Drcd, SS_D.Drcd);
            // limit concentration to 850 and throw rest in deposition

            double sssmax = MAXCONC * DX->Drc *ChannelAdj->Drc*SSD_D.Drcd;
            if(sssmax < SS_D.Drcd)
            {
                DepFlood->Drc += -(SS_D.Drcd - sssmax);
                SSFlood->Drc += -(SS_D.Drcd - sssmax);
                if(SwitchUseMaterialDepth)
                {
                    StorageDep->Drc += (SS_D.Drcd - sssmax);
                    StorageDep_D.Drcd += (SS_D.Drcd - sssmax);
                }
                SS_D.Drcd = sssmax;
            }

            BLC_D.Drcd = MaxConcentration(ChannelAdj->Drc*DX->Drc*BLD_D.Drcd, BL_D.Drcd);
            // limit concentration to 850 and throw rest in deposition

            sssmax = MAXCONC * DX->Drc *ChannelAdj->Drc*BLD_D.Drcd;
            if(sssmax < BL_D.Drcd)
            {
                DepFlood->Drc += -(BL_D.Drcd - sssmax);
                BLFlood->Drc += -(BL_D.Drcd - sssmax);
                if(SwitchUseMaterialDepth)
                {
                    StorageDep->Drc += (BL_D.Drcd - sssmax);
                    StorageDep_D.Drcd += (BL_D.Drcd - sssmax);
                }
                BL_D.Drcd = sssmax;

            }
        }

        if(SwitchUseGrainSizeDistribution)
        {
            BLCFlood->Drc = 0;
            SSCFlood->Drc = 0;

            FOR_GRAIN_CLASSES
            {
                BLCFlood->Drc += BLC_D.Drcd;
                SSCFlood->Drc += SSC_D.Drcd;
            }
        }
    }

    BLFlood->Drc = std::max(0.0,BLFlood->Drc);
    SSFlood->Drc = std::max(0.0,SSFlood->Drc);
}
//--------------------------------------------------------------------------------------------
/**
 * @fn double TWorld::SWOFSedimentTCBL(int r, int c, int _d)
 * @brief Calculates flooding bed load layer sediment transport capacity
 *
 * Calculates flooding bed load layer sediment transport capacity.
 * Based on either Govers, Van Rijn (simplified or full)
 * or Wu wang & Jia (for multiclass sediment).
 *
 * @param r : row nr of cell
 * @param c : column nr of cell
 * @param d : grain class nr
 * @param hm : the flood water height
 * @param um : the flood velocity in the x-direction
 * @param vm : the flood velocity in the y-direction
 *
 * @return void
 *
 * @see FS_BL_Method
 */
double TWorld::SWOFSedimentTCBL(int r, int c, int _d, cTMap * hm, double v)
{

    if(hm->Drc < MIN_HEIGHT || BLDepthFlood->Drc < MIN_HEIGHT)
    {
        return 0;
    }
    if(v < he_ca)
    {
        return 0;
    }

    if(FS_BL_Method == FSRIJN)
    {
        //Van rijn simplified (1984)
        double ps = 2400.0;
        double pw = 1000.0;
        double ucr;
        double d50m = (D50->Drc/1000000.0);
        double d90m = (D90->Drc/1000000.0);
        if( d50m < 0.005) {
            ucr  = 0.19 * pow(d50m, 0.1) * log10(4.0* hm->Drc/d90m);
        } else {
            ucr  = 0.19 * pow(d50m, 0.6) * log10(4.0* hm->Drc/d90m);
        }
        double me = std::max((v - ucr)/(sqrt(GRAV * d50m * ((ps/pw)- 1.0))),0.0);
        double qs = 0.005 * ps*v *hm->Drc * pow(d50m/hm->Drc,1.2) * pow(me, 2.4);
        double tc =  qs/ (v * BLDepthFlood->Drc );
        return std::max(std::min( tc,MAXCONC ),0.0);
    }else if(FS_BL_Method == FSRIJNFULL)
    {
        //van Rijn full (1980)

        double ps = 2400.0;
        double pw = 1000.0;
        double kinvis = 1.0;
        double d50m = (D50->Drc/1000000.0);
        double d90m = (D90->Drc/1000000.0);

        double ds = D50->Drc * pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0));
        double dh = hm->Drc;
        double chevey = 18 * log(4 * dh/d90m);
        double us = sqrt(GRAV) * v/chevey;
        double uscr = 0.055;

        if(ds <150 && !(ds < 20))
            uscr = 0.013*pow(ds,0.29);
        if(ds < 20 && !(ds < 10))
            uscr = 0.04*pow(ds,-0.10);
        if(ds <10 && !(ds < 4))
            uscr = 0.14*pow(ds,-0.64);
        if(ds <4)
            uscr = 0.24*pow(ds,-1);

        uscr = sqrt(uscr * (ps/pw - 1)*GRAV * d50m);
        double T = std::max((us*us/(uscr*uscr)),0.0);
        double qs = 0.053 * (pow(T,2.1)/pow(ds,0.3)) * sqrt((ps/pw -1)*GRAV)*d50m*sqrt(d50m);
        double tc =  ps * ChannelAdj->Drc * qs/ (v * BLDepthFlood->Drc*ChannelAdj->Drc);

        return std::max(std::min(tc,MAXCONC ),0.0);


    }else if(FS_BL_Method == FSWUWANGJIA)
    {

        double slope = Grad->Drc;
        double ps = 2400.0;
        double pw = 1000.0;
        double h = hm->Drc;
        double n = std::max(N->Drc,0.001);
        double na = pow(graindiameters.at(_d)/100000.0,(1.0/6.0))/20.0;
        double phk = 0;
        double pek = 0;
        FOR_GRAIN_CLASSES
        {
            phk += W_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
            pek += W_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
        }
        double ppk = 1;
        ppk = pow(phk/pek,0.6);
        if(pek  == 0)
        {
            return 0;
        }

        double dh = (ChannelAdj->Drc *h)/(ChannelAdj->Drc + 2* h);
        double css = 0.03* (ps - pw) * (graindiameters.at(_d)/1000000.0) * ppk;

        double qs = 0.0053 *pow(std::max(pow(na/n,1.5)*((pw * dh * 9.81 * 0.1 * slope/css) -1.0 ), 0.0),2.2);
        qs = qs * sqrt((ps/pw - 1)*9.81*pow(graindiameters.at(_d)/1000000.0,3.0));

        double tc = ps * ChannelAdj->Drc * qs/ (v * BLDepthFlood->Drc*ChannelAdj->Drc);
        return std::max(std::min(tc,MAXCONCBL ),0.0);

    }else
    {

        return 0;
    }

}
//--------------------------------------------------------------------------------------------
/**
 * @fn double TWorld::SWOFSedimentTCSS(int r, int c, int _d)
 * @brief Calculates flooding suspended layer sediment transport capacity
 *
 * Calculates flooding suspended layer sediment transport capacity.
 * Based on either, Van Rijn (simplified or full)
 * or Wu wang & Jia (for multiclass sediment).
 *
 * @param r : row nr of cell
 * @param c : column nr of cell
 * @param d : grain class nr
 * @param h : the flood water height
 * @param u : the flood velocity in the x-direction
 * @param v : the flood velocity in the y-direction
 *
 * @return void
 *
 * @see FS_SS_Method
 */
double TWorld::SWOFSedimentTCSS(int r, int c, int _d, cTMap * hm,double v)
{

    if(hm->Drc < MIN_HEIGHT || SSDepthFlood->Drc < MIN_HEIGHT)
        return 0;
    if(v < he_ca)
        return 0;

    if(FS_SS_Method == FSGOVERS)
    {
        //Govers with a maximum bed load layer depth (1980)
        //  double discharge = v * ChannelAdj->Drc * BLDepthFlood->Drc;
        //### Calc transport capacity
        double omega = 100.0*v*Grad->Drc;
        // V in cm/s in this formula assuming grad is SINE
        double omegacrit = 0.4;
        // critical unit streampower in cm/s
        return std::min(MAXCONC, 2650 * CG->Drc * pow(std::max(0.0, omega - omegacrit), DG->Drc));
        // not more than 2650*0.32 = 848 kg/m3

    }else
        if(FS_SS_Method == FSRIJN)
        {
            //Van rijn simplified (1984)
            double ucr;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            double ps = 2400.0;
            double pw = 1000.0;
            double mu = 1.0;
            if( d50m < 0.0005)
            {
                ucr  = 0.19 * pow(d50m, 0.1) * log10(4.0* SSDepthFlood->Drc/d90m);
            }else
            {
                ucr  = 8.5 * pow(d50m, 0.6) * log10(4.0* SSDepthFlood->Drc/d90m);
            }
            double me = std::max((v - ucr)/sqrt(GRAV * d50m * (ps/pw - 1)),0.0);
            double ds = d50m * GRAV * ((ps/pw)-1)/(mu*mu);
            double qs = hm->Drc * 0.008 * ps*v * d50m * pow(me, 2.4) * pow(ds, -0.6);

            double tc =  qs/ (v * SSDepthFlood->Drc);
            return std::max(std::min(tc,MAXCONC),0.0);
        }else if(FS_SS_Method == FSRIJNFULL)
        {
            //van Rijn full (1980)

            double kinvis = 1.0;
            double d50m = (D50->Drc/1000000.0);
            double d90m = (D90->Drc/1000000.0);
            double ps = 2400.0;
            double pw = 1000.0;
            double ds = D50->Drc * pow((ps/pw-1)*GRAV/(kinvis*kinvis),(1.0/3.0));
            double dh = hm->Drc;
            double chevey = 18 * log10(4 * dh/d90m);
            double a = 0.1;

            double us = v* sqrt(GRAV)/chevey;
            double uscr = 0.055;

            if(ds <150 && !(ds < 20))
                uscr = 0.013*pow(ds,0.29);
            if(ds < 20 && !(ds < 10))
                uscr = 0.04*pow(ds,-0.10);
            if(ds <10 && !(ds < 4))
                uscr = 0.14*pow(ds,-0.64);
            if(ds <4)
                uscr = 0.24*pow(ds,-1);

            uscr = sqrt(uscr * (ps/pw - 1)*GRAV * d50m);

            double T = std::max((us*us/(uscr*uscr)) - 1.0,0.0);
            double bsv = sqrt(GRAV * hm->Drc *std::max(Grad->Drc,0.05));
            double ca = 0.015 * (d50m/a) * pow(T,1.5)/pow(ds,0.3);

            double dss = 1 + 0.011*(1.8 - 1)*(T - 25);
            double sv = 10 * (kinvis/ds) *( sqrt(1 + (ps/pw - 1) * GRAV * d50m*d50m*d50m) - 1);

            double beta = std::min(1.0 + 2.0*(sv/bsv)*(sv/bsv),5.0);
            double powcb =1;
            if(BLCFlood->Drc > 0)
            {
                powcb = 0.1;//pow(ca/BLCFlood->Drc,0.4);
            }
            double phi = 2.5 * pow(sv/bsv,0.8) * powcb;
            double Z = sv/(beta*bsv*0.40);
            double Zs = Z + phi;
            double ad = 0.1;
            double F = (pow(ad,Zs) - pow(ad,1.2))/(pow(1.0-ad,Zs)* (1.2 - Zs));
            double qs = F * v * hm->Drc * ca;
            double tc = ps * ChannelAdj->Drc * qs/ (v * SSDepthFlood->Drc * ChannelAdj->Drc);
            return std::max(std::min(tc,MAXCONC ),0.0);
        }else if(FS_SS_Method == FSWUWANGJIA)
        {
            if(SSD_D.at(_d)->Drc < 0.004)
            {
                return 0;
            }
            double slope = Grad->Drc;
            double ps = 2400.0;
            double pw = 1000.0;
            double h = hm->Drc;
            double phk = 0;
            double pek = 0;
            double sv = settlingvelocities.at(_d);
            double gd = graindiameters.at(_d)/1000000.0;
            FOR_GRAIN_CLASSES
            {
                phk += W_D.Drcd * (graindiameters.at(d)/(graindiameters.at(_d) + graindiameters.at(d)));
                pek += W_D.Drcd * (graindiameters.at(_d)/(graindiameters.at(d) + graindiameters.at(_d)));
            }

            double ppk = 1;
            ppk = pow(phk/pek,0.6);
            if(pek == 0)
            {
                return 0;
            }

            //double dh = (ChannelAdj->Drc *h)/(ChannelAdj->Drc + 2.0* h);

            double css = 0.03* (ps - pw) * (gd) * ppk;

            double qs = 0.0000262 *pow(std::max(( pw * 0.01 * h /css) - 1.0, 0.0)* v/(sqrt(sv)),2.2);
            qs =  qs * 1 * sqrt((ps/pw - 1)*9.81*pow(gd,3.0));

            double tc = ps * ChannelAdj->Drc * qs/ (v * SSDepthFlood->Drc * ChannelAdj->Drc);
            return std::max(std::min(tc,MAXCONC ),0.0);

        }else
        {

            return 0;
        }


}
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::SWOFSedimentDet(double dt, int r,int c)
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
    if (!SwitchUse2Layer) {
        BLDepthFlood->Drc = 0;
        SSDepthFlood->Drc = h;
        return;
    }

    double ps = 2400;
    double pw = 1000;
    double factor = 0.5;
    if(!SwitchUseGrainSizeDistribution)
    {
        double d50m = (D50->Drc/1000000.0);
        double d90m = (D90->Drc/1000000.0);
        //critical shear velocity for bed level motion by van rijn
        double critshearvel = velocity * sqrt(GRAV)/(18 * log10(4*(ChannelAdj->Drc * h/(h*2 + ChannelAdj->Drc))/d90m));
        //critical shear stress for bed level motion by van rijn
        double critsheart = (critshearvel*critshearvel)/ (((ps-pw)/pw) * GRAV*d50m);
        //rough bed bed load layer depth by Hu en Hui
        BLDepthFlood->Drc = std::min(std::min(d50m * 1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), factor*h), 0.1);
        SSDepthFlood->Drc = std::max(h - BLDepthFlood->Drc,0.0);
    }else
    {
        FOR_GRAIN_CLASSES
        {
            double d50m = graindiameters.at(d)/1000000.0;
            double d90m = 1.5 * graindiameters.at(d)/1000000.0;

            //critical shear velocity for bed level motion by van rijn
            double critshearvel = velocity * sqrt(GRAV)/(18 * log10(4*(ChannelAdj->Drc * h/(h*2 + ChannelAdj->Drc))/(d90m)));
            //critical shear stress for bed level motion by van rijn
            double critsheart = (critshearvel*critshearvel)/ (((ps-pw)/pw) * GRAV*d50m);
            //rough bed bed load layer depth by Hu en Hui
            BLD_D.Drcd = std::min(std::min(d50m * 1.78 * (pow(ps/pw,0.86)*pow(critsheart,0.69)), factor*h), 0.1);
            SSD_D.Drcd = std::max(h - BLD_D.Drcd,0.0);
            BLDepthFlood->Drc += BLD_D.Drcd * W_D.Drcd;
            SSDepthFlood->Drc += SSD_D.Drcd * W_D.Drcd;

        }
    }
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
void TWorld::SWOFSedimentDet(cTMap * DT, int r,int c, cTMap * h,cTMap * u,cTMap * v)
{
    //first calculate layer depth
    double velocity = std::sqrt(u->Drc *u->Drc + v->Drc * v->Drc);
    // gives BLDepthFlood->Drc and SSDepthFlood->Drc

    SWOFSedimentLayerDepth(r,c,h->Drc, velocity);
        //creates BLDepth and SSDepth, or if 1 layer ssdepth = h and bldepth = 0

    //iterator is the number of grain classes
    int iterator = numgrainclasses;
    if(!SwitchUseGrainSizeDistribution)
        iterator = 1;

    //BLDetFlood->Drc = 0;
    //SSDetFlood->Drc = 0;
    //BLDepFlood->Drc = 0;

    // get TCs
    for(int d = 0 ; d < iterator; d++)
    {
        cTMap * TBLTCFlood;
        cTMap * TSSTCFlood;
        //set maps for this grain class
        if(!SwitchUseGrainSizeDistribution) {
            TBLTCFlood = BLTCFlood;
            TSSTCFlood = SSTCFlood;
        } else {
            TBLTCFlood = BLTC_D.at(d);
            TSSTCFlood = SSTC_D.at(d);
        }
        //calculate tranport capacity for bed load and suspended load
        TBLTCFlood->Drc = SWOFSedimentTCBL(r,c,d,h,velocity);  // van Rijn etc
        TSSTCFlood->Drc = SWOFSedimentTCSS(r,c,d,h,velocity);   // Govers, van Rijn etc
    }
    //check for concentrations above MAXCONC
    if(SwitchUseGrainSizeDistribution)
    {

        //set total load as sum of induvidual grain size classes
        BLTCFlood->Drc = 0;
        SSTCFlood->Drc = 0;
        FOR_GRAIN_CLASSES
        {
            BLTCFlood->Drc += BLTC_D.Drcd;
            SSTCFlood->Drc += SSTC_D.Drcd;
        }

        //check if bed load concentration is too high
        if(BLTCFlood->Drc > MAXCONCBL)
        {
            //rescale concentration of grain classes
            FOR_GRAIN_CLASSES
            {
                BLTC_D.Drcd *= MAXCONCBL/BLTCFlood->Drc;
            }
            BLTCFlood->Drc = MAXCONCBL;
        }
        //check if suspended load concentration is too high
        if(SSTCFlood->Drc > MAXCONC)
        {
            //rescale concentration of grain classes
            FOR_GRAIN_CLASSES
            {
                SSTC_D.Drcd *= MAXCONC/SSTCFlood->Drc;
            }
            SSTCFlood->Drc = MAXCONC;
        }

        //set total load as sum of induvidual grain size classes
        BLTCFlood->Drc = 0;
        SSTCFlood->Drc = 0;
        FOR_GRAIN_CLASSES
        {
            BLTCFlood->Drc += BLTC_D.Drcd;
            SSTCFlood->Drc += SSTC_D.Drcd;
        }

    }

    for(int d  = 0 ; d < iterator;d++)
    {
        //set maps for this grain class
        cTMap * TBLTCFlood;
        cTMap * TSSTCFlood;
        cTMap * TBLDepthFlood;
        cTMap * TSSDepthFlood;
        cTMap * TBLCFlood;
        cTMap * TSSCFlood;
        cTMap * TBLFlood;
        cTMap * TSSFlood;
        cTMap * TW;

        double TSettlingVelocity;
        if (!SwitchUseGrainSizeDistribution) {
            TBLDepthFlood = BLDepthFlood;
            TSSDepthFlood = SSDepthFlood;
            TBLTCFlood = BLTCFlood;
            TSSTCFlood = SSTCFlood;
            TBLCFlood = BLCFlood;
            TSSCFlood = SSCFlood;
            TBLFlood = BLFlood;
            TSSFlood = SSFlood;
            TSettlingVelocity = SettlingVelocity->Drc;
            TW = unity;
        } else {
            TBLDepthFlood = BLD_D.at(d);
            TSSDepthFlood = SSD_D.at(d);
            TBLTCFlood = BLTC_D.at(d);
            TSSTCFlood = SSTC_D.at(d);
            TBLCFlood = BLC_D.at(d);
            TSSCFlood = SSC_D.at(d);
            TBLFlood = BL_D.at(d);
            TSSFlood = SS_D.at(d);
            TW = W_D.at(d);
            TSettlingVelocity = settlingvelocities.at(d);
        }

        //discharges for both layers and watervolumes
        //note for non-advanced erosion BLdepth = 0
        double bldischarge = velocity * ChannelAdj->Drc * TBLDepthFlood->Drc;
        double blwatervol = ChannelAdj->Drc * DX->Drc * TBLDepthFlood->Drc;

        double ssdischarge = velocity * ChannelAdj->Drc * TSSDepthFlood->Drc;
        double sswatervol = ChannelAdj->Drc * DX->Drc * TSSDepthFlood->Drc;

        double deposition = 0;
        double detachment = 0;

        //NOTE depostion is independent for SS and BL, because they have diff grainsizes

        if(h->Drc < MIN_HEIGHT)//he_ca)
        {
            //set all to zero when the water height is zero
            DepFlood->Drc += -BLFlood->Drc;
            BLTCFlood->Drc = 0;
            //BLDetFlood->Drc = 0;
            BLFlood->Drc = 0;
            BLCFlood->Drc = 0;

            DepFlood->Drc += -SSFlood->Drc;
            SSTCFlood->Drc = 0;
            SSFlood->Drc = 0;
            SSCFlood->Drc = 0;
            //SSDetFlood->Drc = 0;

            if(SwitchUseGrainSizeDistribution)
            {
                BL_D.Drcd = 0;
                SS_D.Drcd = 0;
                BLTC_D.Drcd = 0;
                SSTC_D.Drcd = 0;
                BLC_D.Drcd = 0;
                SSC_D.Drcd = 0;
            }

            if(SwitchUseMaterialDepth)
            {
                StorageDep->Drc += -deposition;
                if(SwitchUseGrainSizeDistribution)
                {
                    StorageDep_D.Drcd += -deposition;
                }
            }

        } else {
            // there is water

            //first check if sediment goes to suspended sediment layer or to bed layer
            double TransportFactor;

            //deposition based on settling velocity
            //if there is a significant water height
            if (TSSDepthFlood->Drc > MIN_HEIGHT) {
                TransportFactor = (1-exp(-DT->Drc*TSettlingVelocity/h->Drc)) * sswatervol; //NOTE use entire depth h for deposition of SS
            } else {
                TransportFactor =  1.0*sswatervol;
                // all sed is deposited
            }
            double maxTC = std::max(TSSTCFlood->Drc - TSSCFlood->Drc,0.0) ;
            // positive difference: TC deficit becomes detachment (ppositive)
            double minTC = std::min(TSSTCFlood->Drc - TSSCFlood->Drc,0.0) ;
            // negative diff, becomes deposition

            deposition  = std::max(TransportFactor * minTC,-TSSFlood->Drc);

            if(SwitchUseMaterialDepth)
            {
                StorageDep->Drc += -deposition;
                if(SwitchUseGrainSizeDistribution)
                {
                    StorageDep_D.Drcd += -deposition;
                }
            }


            //erosion values based on discharge
            TransportFactor =// ssdischarge*DT->Drc;
                DT->Drc*TSettlingVelocity * DX->Drc * ChannelAdj->Drc;//SoilWidthDX->Drc;
            detachment = TW->Drc * maxTC * TransportFactor;  // TW is 1 or grainsize fraction
            detachment = std::min(detachment, maxTC * ssdischarge*DT->Drc); //VJ 0518 this line is new

            // exceptions
            //    if (FlowBoundary->Drc > 0)
            //       detachment = 0;
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

            //check how much of the potential detachment can be detached from soil layer
            detachment = Y->Drc * detachment;
//DetachMaterial(r,c,d, false, true, false, detachment);
            //bool channel, bool flood,bool bl

            if(MAXCONC * sswatervol < TSSFlood->Drc+detachment)
                detachment = std::max(0.0, MAXCONC * sswatervol - TSSFlood->Drc);
            // not more detachment then is needed to keep below ssmax

            //### sediment balance
            //SSDetFlood->Drc += detachment;// set to zero this timestep
            // actually redundant
            TSSFlood->Drc += deposition;
            TSSFlood->Drc += detachment;
            TSSFlood->Drc = std::max(0.0,TSSFlood->Drc);
            SSDetFlood->Drc += detachment;  // set to zero in mass balance
            DepFlood->Drc += deposition;

            // ########### do bedload
            if(TBLDepthFlood->Drc < MIN_HEIGHT) {

                // if the BLdepth is to small dump everything
                //BLDepFlood->Drc = -BLFlood->Drc;
                DepFlood->Drc += -BLFlood->Drc;
                BLTCFlood->Drc = 0;
                //BLDetFlood->Drc = 0;
                BLFlood->Drc = 0;
                BLCFlood->Drc = 0;

            } else {
                // there is BL transport

                //deposition and detachment

                //### calc concentration and net transport capacity
                maxTC = std::max(TBLTCFlood->Drc - TBLCFlood->Drc,0.0);
                // positive difference: TC deficit becomes detachment (ppositive)
                minTC = std::min(TBLTCFlood->Drc - TBLCFlood->Drc,0.0);
                // negative difference: TC surplus becomes deposition (negative)
                // unit kg/m3

                //### detachment
                TransportFactor = bldischarge*DT->Drc;
                //DT->Drc*TSettlingVelocity * DX->Drc *SoilWidthDX->Drc;
                // detachment can only come from soil, not roads (so do not use flowwidth)
                // units s * m/s * m * m = m3

                detachment = TW->Drc * maxTC * TransportFactor;
                // unit = kg/m3 * m3 = kg


                if(MAXCONC * blwatervol < TBLFlood->Drc+detachment)
                    detachment = std::max(0.0, MAXCONC * blwatervol - TBLFlood->Drc);
                // limit detachment to what BLflood can carry

                //    if (FlowBoundary->Drc > 0)
                //        detachment = 0;
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

                //detachment = DetachMaterial(r,c,d,false,true,true, detachment);
                detachment = detachment * Y->Drc;
                //DetachMaterial(r,c,d,false,false,true, detachment);

                // IN KG/CELL

                //### deposition
                TransportFactor = (1-exp(-DT->Drc*TSettlingVelocity/BLDepthFlood->Drc)) * blwatervol;

                //   TransportFactor = _dt*SettlingVelocity->Drc * DX->Drc * FlowWidth->Drc;
                // deposition can occur on roads and on soil (so use flowwidth)

                // max depo, kg/m3 * m3 = kg, where minTC is sediment surplus so < 0
                deposition = std::max(minTC * TransportFactor, -TBLFlood->Drc);
                // cannot have more depo than sediment present

                //   if (FlowBoundary->Drc > 0)
                //       deposition = 0;
                // VJ 190325 prevent any activity on the boundary!

                if(SwitchUseMaterialDepth)
                {
                    StorageDep->Drc += -deposition;
                    if(SwitchUseGrainSizeDistribution)
                    {
                        StorageDep_D.Drcd += -deposition;
                    }
                }

                //### sediment balance IN KG/CELL
                DepFlood->Drc += deposition;
                BLDetFlood->Drc += detachment;
                TBLFlood->Drc += detachment;
                TBLFlood->Drc += deposition;
                TBLFlood->Drc = std::max(0.0,TBLFlood->Drc);
            }  // BL exist
        } // h > MIN_HEIGHT
    }   // iterator

    if(SwitchUseGrainSizeDistribution)
    {
        BLFlood->Drc = 0;
        SSFlood->Drc = 0;
        BLTCFlood->Drc = 0;
        SSTCFlood->Drc = 0;

        FOR_GRAIN_CLASSES
        {
            BLFlood->Drc += BL_D.Drcd;
            SSFlood->Drc += SS_D.Drcd;
            BLTCFlood->Drc += BLTC_D.Drcd;
            SSTCFlood->Drc += SSTC_D.Drcd;
        }
    }


}
//--------------------------------------------------------------------------------------------
/**
 * @fn void TWorld::SWOFSedimentBalance()
 * @brief Calculates Bed Load and Suspended load in flood water as the sum of all grain size classes
 *
 * @return void
 */
void TWorld::SWOFSedimentBalance(int thread)
{
    if(SwitchUseGrainSizeDistribution)
    {
        //first set to zero
        FOR_ROW_COL_UF2DMT_DT {
            if(FloodHMaskDer->Drc != 0)
            {
                BLFlood->Drc = 0;
                BLCFlood->Drc = 0;
                SSFlood->Drc = 0;
                SSCFlood->Drc = 0;
            }
        }}}}

        //then sum up all induvidual grain size classes
        FOR_ROW_COL_UF2DMT_DT {
            if(FloodHMaskDer->Drc != 0)
            {
                FOR_GRAIN_CLASSES
                {
                    BLFlood->Drc += BL_D.Drcd;
                    BLCFlood->Drc += BLC_D.Drcd;
                    SSFlood->Drc += SS_D.Drcd;
                    SSCFlood->Drc += SSC_D.Drcd;
                }
            }
        }}}}
    }
}
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

void TWorld::SWOFSediment(int thread,cTMap* DT,double dt,cTMap * h,cTMap * u,cTMap * v)
{
    //sediment detachment or deposition
    FOR_ROW_COL_UF2DMT_DT {
        if(FloodHMaskDer->Drc != 0)// && DT->Drc > 1e-6)
        {
            SWOFSedimentDet(DT,r,c,h,u,v);
        }
    }}}}

    //check for cells with insignificant water height and calculate concentration
    FOR_ROW_COL_UF2DMT_DT {
        if(FloodHMaskDer->Drc != 0)
        {
            SWOFSedimentCheckZero(r,c,h);

            SWOFSedimentSetConcentration(r,c,h);
        }
    }}}}


    //transport sediment using velocities and water heights from SWOF
    if(!SwitchUseGrainSizeDistribution)
    {
        SWOFSedimentFlowInterpolation(thread,DT,h,u,v, BLFlood, BLCFlood, SSFlood, SSCFlood);

        SWOFSedimentDiffusion(thread,DT,h,u,v, SSFlood, SSCFlood);

    } else {
         //or if there are multiple grain size classes

         //calculate total sediment as sum of grain classes
        SWOFSedimentBalance(thread);

        //transport sediment using velocities and water heights from SWOF
        FOR_GRAIN_CLASSES
        {
            SWOFSedimentFlowInterpolation(thread,DT,h,u,v, BL_D.at(d), BLC_D.at(d), SS_D.at(d),SSC_D.at(d));
        }

        //calculate total sediment as sum of grain classes
        //SWOFSedimentBalance(thread);

        //diffusion
        FOR_GRAIN_CLASSES
        {
            SWOFSedimentDiffusion(thread,DT,h,u,v,  SS_D.at(d), SSC_D.at(d));
        }

        //calculate total sediment as sum of grain classes
        SWOFSedimentBalance(thread);
    }


    FOR_ROW_COL_UF2DMT_DT {
        if(FloodHMaskDer->Drc != 0)
        {
            SWOFSedimentSetConcentration(r,c,h);
        }
    }}}}


}
