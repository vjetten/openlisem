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
 \file lisUnifiedFlow.cpp
 \brief

functions: \n
*/

#include <algorithm>
#include "model.h"
#include "operation.h"


//General Function
void TWorld::UnifiedFlow()
{


    ////START ALGORITHM
    ////from now on all input and output is provided as funciton arguments
    ////This increases re-usablitiy of the code (use another DEM, Velocity map.. etc.. and everything will remain functional)
    /// Note that fluxes trough N, S, W, E boundaries and accalerations at these boundaries are temporarily stored in certain maps
    //sets up the variables for the spatially dynamic timestep
    UF_DTMIN = UF_InitTimeStep( UF2D_DEM,                                                   //dem info
                                UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                               //channel info
                                UF1D_f,UF1D_visc,UF1D_fu,                                   //1d fluid phase
                                UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,               //1d solid phase
                                UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                           //2d fluid phase
                                UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv,       //2d solid phase
                                UF2D_T,UF2D_DT,UF2D_DTStep,UF1D_T,UF1D_DT,UF1D_DTStep);     //output timesteps

    double t = 0;
    double dt = UF_DTMIN;


    ////START MAIN LOOP
    ////from now on all input and output is provided as function arguments
    ////This increases re-usablitiy of the code

    //continue while we have not made a timestep of _dt
    while(t + UF_VERY_SMALL < _dt)
    {

        cTMap * _dem = UF2D_DEMOriginal;
        FOR_ROW_COL_UF2D
        {
            UF2D_DEM->Drc = UF2D_DEMOriginal->Drc + DEMChange->Drc;
        }
        ////TOPOGRAPHY ANALYSIS
        UF_DEMLDDAnalysis(UF2D_DEM,UF1D_LDD,UF1D_LDDw,UF1D_LDDh,UF1D_f,UF1D_s,UF2D_f,UF2D_s);

        if(UF_NeedsInitial)
        {
            //UF_NeedsInitial = false;
            ////initial conditions
            UF_Initial(    UF2D_DEM,                                        //dem info
                                UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                                UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                                UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                                UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                                UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase
        }

        ////TIMESTEP ANALYSIS
        //analyzes spatially dynamic timstep that should be made
        dt = UF_TimeStep(t,     UF2D_DEM,                                                    //dem info
                                UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                                UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                                UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                                UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                                UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv,        //2d solid phase
                                UF2D_T,UF2D_DT,UF2D_DTStep,UF1D_T,UF1D_DT,UF1D_DTStep);      //output timesteps

            //set threadpool masks;
            ThreadPool->SetMask(UF2D_DEM, UF2D_DT, UF2D_CellR, UF2D_CellC,UF_1DACTIVE? UF1D_LDD : 0,UF1D_DT,UF1D_CellR,UF1D_CellC);


            //run the created function on seperate threads
            ThreadPool->RunDynamicCompute(flowcompute);

            //now wait till all threads are done!
            ThreadPool->WaitForAll();

        //increase timer by the current timstep!
        t = t + dt;

        DEBUG(QString("UF Step t: %1  dt: %2").arg(t).arg(dt));

    }

    ////STOP MAIN LOOP
    ////again uses non-functionparameter variables


    //just a test to see the effective timestep
    /*int count = 0;
    UF_DTAverage = 0;
    cTMap* _dem = UF2D_DEM;
    FOR_ROW_COL_UF2D
    {
        count ++;
        UF_DTAverage += UF2D_DT->Drc;
    }
    UF_DTAverage = UF_DTAverage/count;*/

}

void TWorld::UF_Compute(int thread)
{

    ////SOURCE TERMS
    //both material and momentum source terms are called here
    UF2D_Source(thread,UF2D_DT,UF2D_DEM,UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);

    if(UF_1DACTIVE)
    {
        UF1D_Source(thread,UF1D_DT,UF1D_LDD,UF1D_LDDw,UF1D_LDDh,UF1D_f,UF1D_visc,UF1D_fu,UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su);
    }



    ////2D SCHEME
    //the actual momentum and mass advection/iteration equations are solved here
    UF2D_Scheme(thread,UF2D_DT,UF2D_DEM,UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);

    if(UF_1DACTIVE)
    {

        ////1D SCHEME
        UF1D_Scheme(thread,UF1D_DT,UF1D_LDD,UF1D_LDDw,UF1D_LDDh,UF1D_f,UF1D_visc,UF1D_fu,UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su);

    }


    ////LAX relaxation scheme to correct for courant limitation
    //the connection
    UF2D1D_LaxNumericalCorrection(thread,UF2D_DT,     UF2D_DEM,                                        //dem info
                        UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                        UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                        UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                        UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                        UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv,         //2d solid phase
                        1,3 );


    /////INFILTRATION
    //substract any possible infiltration from the flow water volume
    UF2D1D_Infiltration(thread,UF2D_DT,     UF2D_DEM,                                        //dem info
                        UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                        UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                        UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                        UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                        UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase

    ////FORCED CONDITIONS
    //apply the forced conditions to where this has been set by the user
    //this involves solid and fluid mass, and their respective properties
    UF_ForcedConditions(thread,UF2D_DT,     UF2D_DEM,                                        //dem info
                        UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                        UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                        UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                        UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                        UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase


    FOR_ROW_COL_UF2DMT
    {
        UF2D_tsf->Drc = UF2D_f->Drc> 0? (UF2D_s->Drc + ((UF2D_ssm->Drc + UF2D_blm->Drc)/UF_DENSITY_SUSPENDED))/(UF2D_f->Drc+UF2D_s->Drc + ((UF2D_ssm->Drc + UF2D_blm->Drc)/UF_DENSITY_SUSPENDED)) : 0.0;
    }}}

    ////SEDIMENT AND OTHER SOIL INTERACTIONS
    UnifiedFlowSediment(thread);

    FOR_ROW_COL_UF2DMT
    {
        UF2D_tsf->Drc = UF2D_f->Drc> 0? (UF2D_s->Drc + ((UF2D_ssm->Drc + UF2D_blm->Drc)/UF_DENSITY_SUSPENDED))/(UF2D_f->Drc+UF2D_s->Drc + ((UF2D_ssm->Drc + UF2D_blm->Drc)/UF_DENSITY_SUSPENDED)) : 0.0;
    }}}

    ////INTERACTIONS WITH THE TERRAIN ELEVATION
    UFDEMLDD_Connection(thread);


    if(UF_1DACTIVE)
    {
        ////CONNECTIONS
        //the connection between 2d and 1d flow is solved in this function.
        //fraction of inflow determined by flow velocity and channel width
        UF2D1D_Connection(thread,UF2D_DT,     UF2D_DEM,                                        //dem info
                           UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                           UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                           UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                           UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                           UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase

        UF2D1D_ChannelWater(thread,UF2D_DT,     UF2D_DEM,                                        //dem info
                           UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                           UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                           UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                           UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                           UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase
    }

    ////update the discharge maps
    UF_UpdateDisplayMaps(thread,UF2D_DT,     UF2D_DEM,                                        //dem info
                         UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                         UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                         UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                         UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                         UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase



}

//#define UF_SWAP(x, y, T) do { MaskedRaster<double> SWAP = x->data; x->data = y->data; y->data = SWAP; } while (0)

#define UF_SWAP1D(x, y, T) UF_SWAP1D_MT(thread,x,y)

#define UF_SWAP2D(x, y, T) UF_SWAP2D_MT(thread,x,y)

void TWorld::UF_SWAP2D_MT(int thread,cTMap *x, cTMap *y)
{
    FOR_ROW_COL_UF2DMTDER
    {
        //double temp = x->Drc;
        x->Drc = y->Drc;
        //y->Drc = temp;
    }}}
}
void TWorld::UF_SWAP1D_MT(int thread,cTMap *x, cTMap *y)
{
    FOR_ROW_COL_UF1DMTDER
    {
        //double temp = x->Drc;
        x->Drc = y->Drc;
        //y->Drc = temp;
    }}}
}

double TWorld::UF2D_Source(int thread, cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv)
{

    //fluid and solid mass sources (slope stability, splash detechment, etc..)

    UF2D_FluidSource(thread,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fn);

    if(UF_SOLIDPHASE)
    {
        UF2D_SolidSource(thread,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sn);
    }

    UF_SWAP2D(_f,UF2D_fn,cTMap*);
    if(UF_SOLIDPHASE)
    {
        UF_SWAP2D(_s,UF2D_sn,cTMap*);
    }

    if(UF_SCHEME == UF_SCHEME_BOUNDARYMUSCLE)
    {
        //first recalculate for the momentum source terms
        UF2D_FluidMomentum2Source(thread,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn);

        if(UF_SOLIDPHASE)
        {
            UF2D_SolidMomentum2Source(thread,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sun,UF2D_svn);
        }
    }

}

//2D version
double TWorld::UF2D_Scheme(int thread,cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv)
{

    if(UF_SCHEME == UF_SCHEME_BOUNDARYMUSCLE)
    {
        //advect momentum
        UF2D_Advect2_Momentum(thread,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn,UF2D_sun,UF2D_svn,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2);

        //advect mass
        UF2D_foutflow += UF2D_Advect2_mass(thread,dt,_dem,_f,_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,UF2D_fn,UF2D_qout);
        if(UF_SOLIDPHASE)
        {
            UF2D_soutflow += UF2D_Advect2_mass(thread,dt,_dem,_s,_s,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2,UF2D_sn,UF2D_qsout);
        }

        if(SwitchErosion)
        {
            UF2D_fsoutflow += UF2D_Advect2_mass(thread,dt,_dem,UF2D_blm,_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,0, UF1D_qblout);
            UF2D_fsoutflow += UF2D_Advect2_mass(thread,dt,_dem,UF2D_ssm,_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,0, UF1D_qssout);
            if(SwitchUseGrainSizeDistribution)
            {
                FOR_GRAIN_CLASSES
                {
                    UF2D_fsoutflow += UF2D_Advect2_mass(thread,dt,_dem,UF2D_blm_D.at(d),_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,0);
                    UF2D_fsoutflow += UF2D_Advect2_mass(thread,dt,_dem,UF2D_ssm_D.at(d),_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,0);

                }
            }
        }

        //advect properties
        UF2D_Advect2_prop(thread,dt,_dem,_f,_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,UF2D_visc,0);
        if(UF_SOLIDPHASE)
        {
            UF2D_Advect2_prop(thread,dt,_dem,_s,_s,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2,UF2D_d,0);
            UF2D_Advect2_prop(thread,dt,_dem,_s,_s,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2,UF2D_ifa,0);
            UF2D_Advect2_prop(thread,dt,_dem,_s,_s,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2,UF2D_rocksize,0);
        }

        UF_SWAP2D(_f,UF2D_fn,cTMap*);
        UF_SWAP2D(_fu,UF2D_fun,cTMap*);
        UF_SWAP2D(_fv,UF2D_fvn,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF_SWAP2D(_s,UF2D_sn,cTMap*);
            UF_SWAP2D(_su,UF2D_sun,cTMap*);
            UF_SWAP2D(_sv,UF2D_svn,cTMap*);
        }

        //first recalculate for the momentum source terms
        UF2D_FluidApplyMomentum2(thread,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn);
        UF_SWAP2D(_fu,UF2D_fun,cTMap*);
        UF_SWAP2D(_fv,UF2D_fvn,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF2D_SolidApplyMomentum2(thread,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sun,UF2D_svn);
            UF_SWAP2D(_su,UF2D_sun,cTMap*);
            UF_SWAP2D(_sv,UF2D_svn,cTMap*);
        }

    }



}

void TWorld::UF1D_Source(int thread,cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su)
{
    //fluid and solid mass sources are added here (currently none)
    UF1D_FluidSource(thread,dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fn);

    if(UF_SOLIDPHASE)
    {
        UF1D_SolidSource(thread,dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sn);
        UF_SWAP1D(_s,UF1D_sn,cTMap*);
    }

    UF_SWAP1D(_f,UF1D_fn,cTMap*);

    if(UF_SCHEME == UF_SCHEME_BOUNDARYMUSCLE)
    {
        //first recalculate for the momentum source terms
        UF1D_FluidMomentum2Source(thread,dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun);

        if(UF_SOLIDPHASE)
        {
            UF1D_SolidMomentum2Source(thread,dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sun);
        }
    }

}
void TWorld::UF1D_Scheme(int thread,cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su)
{

    if(UF_SCHEME == UF_SCHEME_BOUNDARYMUSCLE)
    {

        //advect momentum
        UF1D_Advect2_Momentum(thread,dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun,UF1D_sun,UF1D_fq1,UF1D_fq2,UF1D_sq1,UF1D_sq2);

        //advect mass
        UF1D_foutflow += UF1D_Advect2_mass(thread,dt,_ldd,_lddw,_lddh,UF1D_f,UF1D_f,UF1D_fq1,UF1D_fq2,UF1D_fn, UF1D_qout);
        if(UF_SOLIDPHASE)
        {
            UF1D_soutflow += UF1D_Advect2_mass(thread,dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_s,UF1D_sq1,UF1D_sq2,UF1D_sn, UF1D_qsout);
        }

        if(SwitchErosion)
        {
            UF1D_fsoutflow += UF1D_Advect2_mass(thread,dt,_ldd,_lddw,_lddh,UF1D_blm,UF1D_f,UF1D_fq1,UF1D_fq2,0, UF1D_qblout);
            UF1D_fsoutflow += UF1D_Advect2_mass(thread,dt,_ldd,_lddw,_lddh,UF1D_ssm,UF1D_f,UF1D_fq1,UF1D_fq2,0, UF1D_qssout);
            if(SwitchUseGrainSizeDistribution)
            {
                FOR_GRAIN_CLASSES
                {
                    UF1D_fsoutflow += UF1D_Advect2_mass(thread,dt,_ldd,_lddw,_lddh,UF1D_blm_D.at(d),UF1D_f,UF1D_fq1,UF1D_fq2,0);
                    UF1D_fsoutflow += UF1D_Advect2_mass(thread,dt,_ldd,_lddw,_lddh,UF1D_ssm_D.at(d),UF1D_f,UF1D_fq1,UF1D_fq2,0);

                }
            }
        }

        //advect properties
        UF1D_Advect2_prop(thread,dt,_ldd,_lddw,_lddh,UF1D_f,UF1D_f,UF1D_fq1,UF1D_fq2,UF1D_visc,0);

        if(UF_SOLIDPHASE)
        {
            //UF1D_Advect2_prop(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_s,UF1D_sq1,UF1D_sq2,UF1D_d,0);
            UF1D_Advect2_prop(thread,dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_s,UF1D_sq1,UF1D_sq2,UF1D_ifa,0);
            UF1D_Advect2_prop(thread,dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_s,UF1D_sq1,UF1D_sq2,UF1D_rocksize,0);
        }

        UF_SWAP1D(_f,UF1D_fn,cTMap*);
        UF_SWAP1D(_fu,UF1D_fun,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF_SWAP1D(_s,UF1D_sn,cTMap*);
            UF_SWAP1D(_su,UF1D_sun,cTMap*);
        }

        //first recalculate for the momentum source terms
        UF1D_FluidApplyMomentum2(thread,dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun);
        UF_SWAP1D(_fu,UF1D_fun,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF1D_SolidApplyMomentum2(thread,dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sun);
            UF_SWAP1D(_su,UF1D_sun,cTMap*);
        }
    }



}


//set output
void TWorld::UF_SetInput()
{
    cTMap * _dem = UF2D_DEM;
    cTMap * _ldd = UF1D_LDD;
    FOR_ROW_COL_UF2D
    {
        UF2D_Test->Drc = 0;
        UF2D_f->Drc = WHrunoff->Drc * FlowWidth->Drc * DX->Drc;

        UF2D_q->Drc = 0;
        UF2D_qs->Drc = 0;

        UF2D_qout->Drc = 0;
        UF2D_qsout->Drc = 0;
        UF2D_qblout->Drc = 0;
        UF2D_qssout->Drc = 0;
        UF1D_qout->Drc = 0;
        UF1D_qsout->Drc = 0;
        UF1D_qblout->Drc = 0;
        UF1D_qssout->Drc = 0;
    }
    FOR_ROW_COL_UF1D
    {

        UF1D_q->Drc = 0;
        UF1D_qs->Drc = 0;
    }

    UF2D_foutflow = 0;
    UF2D_fsoutflow = 0;
    UF2D_soutflow = 0;
    UF1D_foutflow = 0;
    UF1D_fsoutflow = 0;
    UF1D_soutflow = 0;


    if(UF_Input_first)
    {
        UF_Input_first = false;

        //here is the place to add water or solid volumes for testing extreme scenarios
        FOR_ROW_COL_UF2D
        {
            if(r != 250 && c!= 250)
            {
                //UF2D_f->Drc += 10;
                //UF2D_s->Drc += 5;
            }
        }
        FOR_ROW_COL_UF1D
        {
            //UF1D_f->Drc += 100;
            //UF1D_s->Drc += 50;
        }
    }
}

//set output
void TWorld::UF_SetOutput()
{


    cTMap * _dem = UF2D_DEM;
    FOR_ROW_COL_UF2D
    {
        //just for display
        UF2D_h->Drc = (UF2D_f->Drc + UF2D_s->Drc)/(_dx*_dx);
        UF2D_fsConc->Drc = (UF2D_f->Drc + UF2D_s->Drc) > UF_VERY_SMALL? (UF2D_ssm->Drc + UF2D_blm->Drc)/(UF2D_f->Drc + UF2D_s->Drc) : 0.0;
        UF2D_sConc->Drc = (UF2D_f->Drc + UF2D_s->Drc) > UF_VERY_SMALL? (UF2D_s->Drc * UF2D_d->Drc)/(UF2D_f->Drc + UF2D_s->Drc) : 0.0 ;
        UF2D_tConc->Drc = UF2D_fsConc->Drc + UF2D_sConc->Drc;
        UF2D_u->Drc = (UF2D_f->Drc + UF2D_s->Drc) > UF_VERY_SMALL? (UF2D_f->Drc * UF2D_fu->Drc + UF2D_s->Drc * UF2D_su->Drc)/(UF2D_f->Drc + UF2D_s->Drc) : 0.0;
        UF2D_v->Drc = (UF2D_f->Drc + UF2D_s->Drc) > UF_VERY_SMALL? (UF2D_f->Drc * UF2D_fv->Drc + UF2D_s->Drc * UF2D_sv->Drc)/(UF2D_f->Drc + UF2D_s->Drc) : 0.0;
        UF2D_velocity->Drc = std::sqrt(UF2D_u->Drc * UF2D_u->Drc + UF2D_v->Drc * UF2D_v->Drc);
        UF2D_q->Drc = std::fabs(UF2D_q->Drc/_dt);
        UF2D_qs->Drc = std::fabs(UF2D_qs->Drc/_dt);


        UF1D_h->Drc = (UF1D_f->Drc + UF1D_s->Drc)/(_dx * UF1D_LDDw->Drc);
        UF1D_fsConc->Drc = (UF1D_f->Drc + UF1D_s->Drc) > UF_VERY_SMALL? (UF1D_ssm->Drc + UF1D_blm->Drc)/(UF1D_f->Drc + UF1D_s->Drc) : 0.0;
        UF1D_sConc->Drc = (UF1D_f->Drc + UF1D_s->Drc) > UF_VERY_SMALL?(UF1D_s->Drc * UF1D_d->Drc)/(UF1D_f->Drc + UF1D_s->Drc) : 0.0;
        UF1D_tConc->Drc = UF1D_fsConc->Drc + UF1D_sConc->Drc;
        UF1D_velocity->Drc = (UF1D_f->Drc + UF1D_s->Drc) > UF_VERY_SMALL? std::fabs((UF1D_f->Drc * UF1D_fu->Drc + UF1D_s->Drc * UF1D_su->Drc)/(UF1D_f->Drc + UF1D_s->Drc)) : 0.0;
        UF1D_q->Drc = std::fabs(UF1D_q->Drc/_dt);
        UF1D_qs->Drc = std::fabs(UF1D_qs->Drc/_dt);

        UF2D_qout->Drc = UF2D_qout->Drc/_dt;
        UF2D_qsout->Drc = UF2D_qsout->Drc/_dt;
        UF2D_qblout->Drc = UF2D_qblout->Drc/_dt;
        UF2D_qssout->Drc = UF2D_qssout->Drc/_dt;
        UF1D_qout->Drc = UF1D_qout->Drc/_dt;
        UF1D_qsout->Drc = UF1D_qsout->Drc/_dt;
        UF1D_qblout->Drc = UF1D_qblout->Drc/_dt;
        UF1D_qssout->Drc = UF1D_qssout->Drc/_dt;

        UF2D_TimeStep->Drc = UF_DTMIN * UF2D_DTStep->Drc;
        UF2D_FPH->Drc = UF2D_f->Drc/(_dx*_dx);
        UF2D_SPH->Drc = UF2D_s->Drc/(_dx*_dx);
        //return water height to the rest of OpenLisem
        if(ChannelAdj->Drc > 0)
        {
            WHrunoff->Drc = UF2D_f->Drc/(DX->Drc*ChannelAdj->Drc);
        }else
        {
            WHrunoff->Drc = 0;
        }
        Qn->Drc = UF2D_q->Drc + UF1D_q->Drc;
        Q->Drc = Qn->Drc;

        WHroad->Drc = WHrunoff->Drc;
        // set road to average outflowing wh, no surface storage.


        WH->Drc = WHrunoff->Drc + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        InfilVolKinWave->Drc = UF2D_Infiltration->Drc;

        if(SwitchIncludeChannel && SwitchChannelInfil)
        {
            InfilVolKinWave->Drc += UF1D_Infiltration->Drc;
        }

        UF2D_Infiltration->Drc = 0;
        UF1D_Infiltration->Drc = 0;

        V->Drc = UF2D_velocity->Drc;

        // recalc velocity for output to map, is not used in other processes
        WaterVolall->Drc = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        // is the same as :         WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);



    }


}

void TWorld::UF_UpdateDisplayMaps(int thread,cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                                 cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                                 cTMap * _fu1D,cTMap * _s1D,
                                 cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                                 cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                                 cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                                 cTMap * _su2D,cTMap * _sv2D)
{

    FOR_ROW_COL_UF2DMT_DT
    {
        double q = (UF2D_fqx1->Drc - UF2D_fqx2->Drc) + (UF2D_fqy1->Drc - UF2D_fqy2->Drc);
        double fsc = (UF2D_f->Drc) > UF_VERY_SMALL? (UF2D_ssm->Drc + UF2D_blm->Drc)/(UF2D_f->Drc) : 0.0;
        UF2D_q->Drc += q;
        UF2D_qs->Drc += (UF2D_sqx1->Drc - UF2D_sqx2->Drc) + (UF2D_sqy1->Drc - UF2D_sqy2->Drc) + q * fsc;
    }}}

    FOR_ROW_COL_UF1DMT_DT
    {
        double q = (UF1D_fq1->Drc - UF1D_fq2->Drc);
        double fsc = (UF1D_f->Drc) > UF_VERY_SMALL? (UF1D_ssm->Drc + UF1D_blm->Drc)/(UF1D_f->Drc) : 0.0;
        UF1D_q->Drc += q;
        UF1D_qs->Drc += (UF1D_sq1->Drc - UF1D_sq2->Drc) + q * fsc;
    }}}
}

void TWorld::UF2D1D_LaxNumericalCorrection(int thread,cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                                 cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                                 cTMap * _fu1D,cTMap * _s1D,
                                 cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                                 cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                                 cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                                 cTMap * _su2D,cTMap * _sv2D, int nh, int nv)
{

    //The Lax numerical correction is usefull to correct for artifacts of the courant factor limitation.
    //when velocity exceeds the courant limit, block patterning appears as a result of flow limitations. Using the Lax-scheme, values are averaged with surrounding cells while maintaining mass balance
    //in this implementation, flow barriers and the dem are taken into accaunt

    FOR_ROW_COL_UF2DMT_DT
    {
        double hf = _f2D->data[r][c]/(_dx*_dx);
        double hs = _s2D->data[r][c]/(_dx*_dx);
        double w = std::min(1.0, std::max(0.0,((std::max(std::fabs(hf * 2.5* _fu2D->Drc),std::max(std::fabs(hf *2.5 * _fv2D->Drc),std::max(std::fabs(hs *_su2D->Drc),std::fabs(hs *_sv2D->Drc)))))/(UF_Courant * _dx) - 0.3)*0.5));
        UF2D_Test->Drc = w;

        if(w > 0.25)
        {
            int count = 0;
            double h = _f2D->data[r][c]; double u = _fu2D->data[r][c]; double v = _fv2D->data[r][c];
            double h1 = h,h2 = h,h3 = h,h4 = h,u1 = u,u2 = u,u3 = u,u4 = u,v1 = v,v2 = v,v3 = v,v4 = v;

            if(!UF_OUTORMV(_dem,r-1,c))
            {
                //if(GetFlowBarrierHeight(r,c,-1,0) == 0 && GetFlowBarrierHeight(r-1,c,1,0) == 0)
                {
                    count ++;
                    h1 = _f2D->data[r-1][c];
                    u1 = _fu2D->data[r-1][c];
                    v1 = _fv2D->data[r-1][c];
                    double fn = (1.0*(5.0 - w) * _f2D->Drc + w *(h1))/5.0;
                    double dif = ((fn - _f2D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,h1/(_dx*_dx)-GetFlowBarrierHeight(r,c,-1,0)) + std::min(0.0,_dem->Drc + GetFlowBarrierHeight(r,c,-1,0) -_dem->data[r-1][c] - GetFlowBarrierHeight(r-1,c,1,0)))*_dx*_dx),std::fabs(fn - _f2D->Drc));
                    _f2D->data[r-1][c] = std::max(0.0,_f2D->data[r-1][c] - dif);
                    double dif2 = h1 - _f2D->data[r-1][c];
                    _f2D->Drc += dif2;
                }
            }
            if(!UF_OUTORMV(_dem,r+1,c))
            {
                //if(GetFlowBarrierHeight(r,c,+1,0) == 0 && GetFlowBarrierHeight(r+1,c,-1,0) == 0)
                {
                    count ++;
                    h2 = _f2D->data[r+1][c];
                    u2 = _fu2D->data[r+1][c];
                    v2 = _fv2D->data[r+1][c];
                    double fn = (1.0*(5.0 - w) * _f2D->Drc + w *(h2))/5.0;
                    double dif = ((fn - _f2D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,h2/(_dx*_dx)-GetFlowBarrierHeight(r,c,1,0)) + std::min(0.0,_dem->Drc + GetFlowBarrierHeight(r,c,1,0) -_dem->data[r+1][c] - GetFlowBarrierHeight(r+1,c,-1,0)))*_dx*_dx),std::fabs(fn - _f2D->Drc));
                    _f2D->data[r+1][c] = std::max(0.0,_f2D->data[r+1][c] - dif);
                    double dif2 = h2 - _f2D->data[r+1][c];
                    _f2D->Drc += dif2;
                }
            }
            if(!UF_OUTORMV(_dem,r,c-1))
            {
                //if(GetFlowBarrierHeight(r,c,0,-1) == 0 && GetFlowBarrierHeight(r,c-1,0,1) == 0)
                {
                    count ++;
                    h3 = _f2D->data[r][c-1];
                    u3 = _fu2D->data[r][c-1];
                    v3 = _fv2D->data[r][c-1];
                    double fn = (1.0*(5.0 - w) * _f2D->Drc + w *(h3))/5.0;
                    double dif = ((fn - _f2D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,h3/(_dx*_dx)-GetFlowBarrierHeight(r,c,0,-1)) + std::min(0.0,_dem->Drc + GetFlowBarrierHeight(r,c,0,-1) -_dem->data[r][c-1] - GetFlowBarrierHeight(r,c-1,0,1)))*_dx*_dx),std::fabs(fn - _f2D->Drc));
                    _f2D->data[r][c-1] = std::max(0.0,_f2D->data[r][c-1] - dif);
                    double dif2 = h3 - _f2D->data[r][c-1];
                    _f2D->Drc += dif2;
                }
            }
            if(!UF_OUTORMV(_dem,r,c+1))
            {
                //if(GetFlowBarrierHeight(r,c,0,1) == 0 && GetFlowBarrierHeight(r,c+1,0,-1) == 0)
                {
                    count ++;
                    h4 = _f2D->data[r][c+1];
                    u4 = _fu2D->data[r][c+1];
                    v4 = _fv2D->data[r][c+1];
                    double fn = (1.0*(5.0 - w) * _f2D->Drc + w *(h4))/5.0;
                    double dif = ((fn - _f2D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,h4/(_dx*_dx)-GetFlowBarrierHeight(r,c,0,1)) + std::min(0.0,_dem->Drc + GetFlowBarrierHeight(r,c,0,1) -_dem->data[r][c+1] - GetFlowBarrierHeight(r,c+1,0,-1)))*_dx*_dx),std::fabs(fn - _f2D->Drc));
                    _f2D->data[r][c+1] = std::max(0.0,_f2D->data[r][c+1] - dif);
                    double dif2 = h4 - _f2D->data[r][c+1];
                    _f2D->Drc += dif2;
                }
            }

            {
                _fu2D->Drc = ((6.0 - 4.0 *w) * _fu2D->Drc + w *(u1 + u2 + u3 + u4))/(6.0 );
                _fv2D->Drc = ((6.0 - 4.0 *w) * _fv2D->Drc + w *(v1 + v2 + v3 + v4))/(6.0 );
            }
        }

        _fu2D->Drc = (_fu2D->Drc > 0? 1.0:-1.0) * std::min(std::fabs(_fu2D->Drc),UF_MAX_NUM_VEL);
        _fv2D->Drc = (_fv2D->Drc > 0? 1.0:-1.0) * std::min(std::fabs(_fv2D->Drc),UF_MAX_NUM_VEL);

        if(std::isnan(_fu2D->Drc) || std::isnan(_fv2D->Drc))
        {
            _fu2D->Drc = 0;
            _fv2D->Drc = 0;
        }

        if(UF_SOLIDPHASE)
        {
            if(w > 0.25)
            {
                int count = 0;
                double sh = _s2D->data[r][c]; double su = _su2D->data[r][c]; double sv = _sv2D->data[r][c];
                double sh1 = sh, sh2 = sh, sh3 = sh, sh4 = sh,su1 =su,su2 = su,su3 = su,su4 = su,sv1 = sv,sv2 = sv,sv3 = sv,sv4 = sv;

                if(!UF_OUTORMV(_dem,r-1,c))
                {
                    //if(GetFlowBarrierHeight(r,c,-1,0) == 0 && GetFlowBarrierHeight(r-1,c,1,0) == 0)
                    {
                        count ++;
                        sh1 = _s2D->data[r-1][c];
                        su1 = _su2D->data[r-1][c];
                        sv1 = _sv2D->data[r-1][c];
                        double sn = (1.0*(5.0 - w) * _s2D->Drc + w *(sh1))/5.0;
                        double dif = ((sn - _s2D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,sh1/(_dx*_dx)-GetFlowBarrierHeight(r,c,-1,0)) + std::min(0.0,_dem->Drc + GetFlowBarrierHeight(r,c,-1,0) -_dem->data[r-1][c] - GetFlowBarrierHeight(r-1,c,1,0)))*_dx*_dx),std::fabs(sn - _s2D->Drc));
                        _s2D->data[r-1][c] = std::max(0.0,_s2D->data[r-1][c] - dif);
                        double dif2 = sh1 - _s2D->data[r-1][c];
                        _s2D->Drc += dif2;
                    }

                }
                if(!UF_OUTORMV(_dem,r+1,c))
                {
                    //if(GetFlowBarrierHeight(r,c,+1,0) == 0 && GetFlowBarrierHeight(r+1,c,-1,0) == 0)
                    {
                        count ++;
                        sh2 = _s2D->data[r+1][c];
                        su2 = _su2D->data[r+1][c];
                        sv2 = _sv2D->data[r+1][c];
                        double sn = (1.0*(5.0 - w) * _s2D->Drc + w *(sh2))/5.0;
                        double dif = ((sn - _s2D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,sh2/(_dx*_dx)-GetFlowBarrierHeight(r,c,1,0)) + std::min(0.0,_dem->Drc + GetFlowBarrierHeight(r,c,1,0) -_dem->data[r+1][c] - GetFlowBarrierHeight(r+1,c,-1,0)))*_dx*_dx),std::fabs(sn - _s2D->Drc));
                        _s2D->data[r+1][c] = std::max(0.0,_s2D->data[r+1][c] - dif);
                        double dif2 = sh2 - _s2D->data[r+1][c];
                        _s2D->Drc += dif2;
                    }
                }
                if(!UF_OUTORMV(_dem,r,c-1))
                {
                    //if(GetFlowBarrierHeight(r,c,0,-1) == 0 && GetFlowBarrierHeight(r,c-1,0,1) == 0)
                    {
                        count ++;
                        sh3 = _s2D->data[r][c-1];
                        su3 = _su2D->data[r][c-1];
                        sv3 = _sv2D->data[r][c-1];
                        double sn = (1.0*(5.0 - w) * _s2D->Drc + w *(sh3))/5.0;
                        double dif = ((sn - _s2D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,sh3/(_dx*_dx)-GetFlowBarrierHeight(r,c,0,-1)) + std::min(0.0,_dem->Drc + GetFlowBarrierHeight(r,c,0,-1) -_dem->data[r][c-1] - GetFlowBarrierHeight(r,c-1,0,1)))*_dx*_dx),std::fabs(sn - _s2D->Drc));
                        _s2D->data[r][c-1] = std::max(0.0,_s2D->data[r][c-1] - dif);
                        double dif2 = sh3 - _s2D->data[r][c-1];
                        _s2D->Drc += dif2;
                    }
                }
                if(!UF_OUTORMV(_dem,r,c+1))
                {
                    //if(GetFlowBarrierHeight(r,c,0,1) == 0 && GetFlowBarrierHeight(r,c+1,0,-1) == 0)
                    {
                        count ++;
                        sh4 = _s2D->data[r][c+1];
                        su4 = _su2D->data[r][c+1];
                        sv4 = _sv2D->data[r][c+1];
                        double sn = (1.0*(5.0 - w) * _s2D->Drc + w *(sh3))/5.0;
                        double dif = ((sn - _s2D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,sh4/(_dx*_dx)-GetFlowBarrierHeight(r,c,0,1)) + std::min(0.0,_dem->Drc + GetFlowBarrierHeight(r,c,0,1) -_dem->data[r][c+1] - GetFlowBarrierHeight(r,c+1,0,-1)))*_dx*_dx),std::fabs(sn - _s2D->Drc));
                        _s2D->data[r][c+1] = std::max(0.0,_s2D->data[r][c+1] - dif);
                        double dif2 = sh4 - _s2D->data[r][c+1];
                        _s2D->Drc += dif2;
                    }
                }


                if(count == 4)
                {
                    _su2D->Drc = ((6.0 - 4.0 *w) * _su2D->Drc + w *(su1 + su2 + su3 + su4))/(6.0);
                    _sv2D->Drc = ((6.0 - 4.0 *w) * _sv2D->Drc + w *(sv1 + sv2 + sv3 + sv4))/(6.0);
                }

            }

            _su2D->Drc = (_su2D->Drc > 0? 1.0:-1.0) * std::min(std::fabs(_su2D->Drc),UF_MAX_NUM_VEL);
            _sv2D->Drc = (_sv2D->Drc > 0? 1.0:-1.0) * std::min(std::fabs(_sv2D->Drc),UF_MAX_NUM_VEL);

            if(std::isnan(_su2D->Drc) || std::isnan(_sv2D->Drc))
            {
                _su2D->Drc = 0;
                _sv2D->Drc = 0;
            }



        }
    }}}

    FOR_ROW_COL_UF1DMT_DT
    {
        double w = std::max(0.0,std::min(1.0,(_f1D->Drc/(_dx *_lddw->Drc)) / 1.5));//std::min(1.0,std::max(0.0,((_f1D->Drc/(_dx *_lddw->Drc))*std::max(std::fabs(2.5* _fu1D->Drc),std::fabs(2.5 * _su1D->Drc))/(UF_Courant * _dx) - 0.65)*0.3));

        if(w > 0.25)
        {
            int count = 0;
            double h = _f1D->data[r][c]; double u = _fu1D->data[r][c];
            double h1 = h,h2 = h,u1 = u,u2 = u;
            double slope = UF1D_Slope->Drc;

            if(!UF_LDDOUT(_ldd,r,c,true))
            {
                count ++;
                h1 = UF1D_Value(_ldd,_lddw,r,c,true,_f1D);
                u1 = UF1D_Value(_ldd,_lddw,r,c,true,_fu1D);

                double fn = (1.0*(3.0 - w) * _f1D->Drc + w *(h1))/3.0;

                double dem2 =_dem->Drc + slope * _dx;
                double dif = ((fn - _f1D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,h1/(_dx*_lddw->Drc)) + std::min(0.0,_dem->Drc -dem2 ))*_dx*_lddw->Drc),std::fabs(fn - _f1D->Drc));

                UF1D_AddValue(_ldd,_lddw,r,c,true,_f1D,- dif,true);
                double dif2 = h1 - UF1D_Value(_ldd,_lddw,r,c,true,_f1D);;
                _f1D->Drc += dif2;
            }
            if(!UF_LDDOUT(_ldd,r,c,false))
            {
                count ++;
                h2 = UF1D_Value(_ldd,_lddw,r,c,false,_f1D);
                u2 = UF1D_Value(_ldd,_lddw,r,c,false,_fu1D);

                double fn = (1.0*(3.0 - w) * _f1D->Drc + w *(h1))/3.0;

                double dem2 =_dem->Drc - slope * _dx;

                double dif = ((fn - _f1D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,h2/(_dx*_lddw->Drc)) + std::min(0.0,_dem->Drc -dem2 ))*_dx*_lddw->Drc),std::fabs(fn - _f1D->Drc));

                UF1D_AddValue(_ldd,_lddw,r,c,false,_f1D,- dif,false);
                double dif2 = h2 - UF1D_Value(_ldd,_lddw,r,c,false,_f1D);
                _f1D->Drc += dif2;
            }

            //_fu1D->Drc = (4.0*(1.0 - w) * _fu1D->Drc + w *(u1 + u2))/(4.0 + count);

        }
        _f1D->Drc = std::max(0.0,_f1D->Drc);
        _fu1D->Drc = (_fu1D->Drc > 0? 1.0:-1.0) * std::min(std::fabs(_fu1D->Drc),UF_MAX_NUM_VEL);

        if(std::isnan(_fu1D->Drc))
        {
            _fu1D->Drc = 0;
        }

        if(UF_SOLIDPHASE)
        {
            if(w > 0.25)
            {
                int count = 0;
                double h = _s1D->data[r][c]; double u = _su1D->data[r][c];
                double h1 = h,h2 = h,u1 = u,u2 = u;
                double slope = UF1D_Slope->Drc;

                if(!UF_LDDOUT(_ldd,r,c,true))
                {
                    count ++;
                    h1 = UF1D_Value(_ldd,_lddw,r,c,true,_s1D);
                    u1 = UF1D_Value(_ldd,_lddw,r,c,true,_su1D);

                    double sn = (1.0*(3.0 - w) * _s1D->Drc + w *(h1))/3.0;

                    double dem2 =_dem->Drc + slope * _dx;
                    double dif = ((sn - _s1D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,h1/(_dx*_lddw->Drc)) + std::min(0.0,_dem->Drc -dem2 ))*_dx*_lddw->Drc),std::fabs(sn - _s1D->Drc));

                    UF1D_AddValue(_ldd,_lddw,r,c,true,_s1D,- dif,true);
                    double dif2 = h1 - UF1D_Value(_ldd,_lddw,r,c,true,_s1D);
                    _s1D->Drc += dif2;
                }
                if(!UF_LDDOUT(_ldd,r,c,false))
                {
                    count ++;
                    h2 = UF1D_Value(_ldd,_lddw,r,c,false,_s1D);
                    u2 = UF1D_Value(_ldd,_lddw,r,c,false,_su1D);

                    double sn = (1.0*(3.0 - w) * _s1D->Drc + w *(h1))/3.0;

                    double dem2 =_dem->Drc - slope * _dx;

                    double dif = ((sn - _s1D->Drc) >0? 1.0:-1.0 )*std::min(std::max(0.0,0.5 * (std::max(0.0,h2/(_dx*_lddw->Drc)) + std::min(0.0,_dem->Drc -dem2 ))*_dx*_lddw->Drc),std::fabs(sn - _s1D->Drc));

                    UF1D_AddValue(_ldd,_lddw,r,c,false,_s1D,- dif,false);
                    double dif2 = h2 - UF1D_Value(_ldd,_lddw,r,c,false,_s1D);
                    _s1D->Drc += dif2;
                }

                //_su1D->Drc = (4.0*(1.0 - w) * _su1D->Drc + w *(u1 + u2))/(4.0 + count);

            }
            _s1D->Drc = std::max(0.0,_s1D->Drc);
            _su1D->Drc = (_su1D->Drc > 0? 1.0:-1.0) * std::min(std::fabs(_su1D->Drc),UF_MAX_NUM_VEL);

            if(std::isnan(_su1D->Drc))
            {
                _su1D->Drc = 0;
            }


        }
    }}}

}
