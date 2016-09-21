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
    qDebug() << "inittimestep";
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

    qDebug() << "start loop";
    //continue while we have not made a timestep of _dt
    while(t + UF_VERY_SMALL < _dt)
    {


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

        qDebug() << "timestep";
        ////TIMESTEP ANALYSIS
        //analyzes spatially dynamic timstep that should be made
        dt = UF_TimeStep(t,     UF2D_DEM,                                                    //dem info
                                UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                                UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                                UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                                UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                                UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv,        //2d solid phase
                                UF2D_T,UF2D_DT,UF2D_DTStep,UF1D_T,UF1D_DT,UF1D_DTStep);      //output timesteps


        //qDebug() << "set threadpool masks";
        //ThreadPool->SetMask(UF2D_DEM, UF2D_DT, UF2D_CellR, UF2D_CellC, true);

        qDebug() << "2d source";
        ////SOURCE TERMS
        //both material and momentum source terms are called here
        UF2D_Source(UF2D_DT,UF2D_DEM,UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);

        if(UF_1DACTIVE)
        {
            qDebug() << "1d source";
            UF1D_Source(UF1D_DT,UF1D_LDD,UF1D_LDDw,UF1D_LDDh,UF1D_f,UF1D_visc,UF1D_fu,UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su);
        }

        qDebug() << "2d scheme";
        ////2D SCHEME
        //the actual momentum and mass advection/iteration equations are solved here
        UF2D_Scheme(UF2D_DT,UF2D_DEM,UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);

        if(UF_1DACTIVE)
        {

            qDebug() << "1D scheme";
            ////1D SCHEME
            UF1D_Scheme(UF1D_DT,UF1D_LDD,UF1D_LDDw,UF1D_LDDh,UF1D_f,UF1D_visc,UF1D_fu,UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su);

            qDebug() << "Connections";
            ////CONNECTIONS
            //the connection between 2d and 1d flow is solved in this function.
            //fraction of inflow determined by flow velocity and channel width
            UF2D1D_Connection(UF2D_DT,     UF2D_DEM,                                        //dem info
                               UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                               UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                               UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                               UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                               UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase
        }

        UF2D1D_LaxNumericalCorrection(UF2D_DT,     UF2D_DEM,                                        //dem info
                            UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                            UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                            UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                            UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                            UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase

        /////INFILTRATION
        //substract any possible infiltration from the flow water volume
        UF2D1D_Infiltration(UF2D_DT,     UF2D_DEM,                                        //dem info
                            UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                            UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                            UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                            UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                            UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase

        ////FORCED CONDITIONS
        //apply the forced conditions to where this has been set by the user
        //this involves solid and fluid mass, and their respective properties
        UF_ForcedConditions(UF2D_DT,     UF2D_DEM,                                        //dem info
                            UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                            UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                            UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                            UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                            UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase

        //increase timer by the current timstep!
        t = t + dt;

        DEBUG(QString("UF Step t: %1  dt: %2").arg(t).arg(dt));
    }

    ////STOP MAIN LOOP
    ////again uses non-functionparameter variables

    ////SOIL INTERACTIONS
    //(sediment transport and solid phase transport is done together with the fluid equations since these are completely integrated)
    if(SwitchErosion)
    {
        qDebug() << "sediment";
        UnifiedFlowSediment();
    }if(SwitchEntrainment && UF_SOLIDPHASE)
    {
        qDebug() << "entrainment";
        UnifiedFlowEntrainment();
    }

    ////STOP ALGORITHM
    ////set output maps for display etc..
    UF_SetOutput();
}

#define UF_SWAP(x, y, T) do { MaskedRaster<double> SWAP = x->data; x->data = y->data; y->data = SWAP; } while (0)

double TWorld::UF2D_Source(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv)
{

    //fluid and solid mass sources (slope stability, splash detechment, etc..)

    //Tried using a custom-made threadpool, to do parrallel computing of induvidual functions
    //Can only be done when function parameters are unrelated of course

    /*RUN_ON_THREAD(UF2D_FluidSource,TWorld,this,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fn);

    if(UF_SOLIDPHASE)
    {
        RUN_ON_THREAD(UF2D_SolidSource,TWorld,this,dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sn);
    }

    ThreadPool->WaitForAll();*/

    UF2D_FluidSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fn);

    if(UF_SOLIDPHASE)
    {
        UF2D_SolidSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sn);
    }

    UF_SWAP(_f,UF2D_fn,cTMap*);
    if(UF_SOLIDPHASE)
    {
        UF_SWAP(_s,UF2D_sn,cTMap*);
    }

    if(UF_SCHEME == UF_SCHEME_CENTRALSIMPLE)
    {
        //first recalculate for the momentum source terms
        UF2D_FluidMomentumSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn);

        if(UF_SOLIDPHASE)
        {
            UF2D_SolidMomentumSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sun,UF2D_svn);
        }

    }else if(UF_SCHEME == UF_SCHEME_BOUNDARYMUSCLE)
    {
        //first recalculate for the momentum source terms
        UF2D_FluidMomentum2Source(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn);

        if(UF_SOLIDPHASE)
        {
            UF2D_SolidMomentum2Source(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sun,UF2D_svn);
        }
    }

}

//2D version
double TWorld::UF2D_Scheme(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv)
{

    //If the less accurate but slightly faster simple central scheme is set by the user, perform this
    if(UF_SCHEME == UF_SCHEME_CENTRALSIMPLE)
    {
        //advect momentum
        UF2D_Advect_Momentum(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn,UF2D_sun,UF2D_svn);

        //advect mass
        UF2D_foutflow += UF2D_Advect_mass(dt,_dem,_f,_fu,_fv,UF2D_fn);
        if(UF_SOLIDPHASE)
        {
            UF2D_soutflow += UF2D_Advect_mass(dt,_dem,_s,_su,_sv,UF2D_sn);
        }

        if(SwitchErosion)
        {
            UF2D_fsoutflow += UF2D_Advect_mass(dt,_dem,UF2D_blm,_fu,_fv,0);
            UF2D_fsoutflow += UF2D_Advect_mass(dt,_dem,UF2D_ssm,_fu,_fv,0);
            if(SwitchUseGrainSizeDistribution)
            {
                FOR_GRAIN_CLASSES
                {
                    UF2D_fsoutflow += UF2D_Advect_mass(dt,_dem,UF2D_blm_D.at(d),_fu,_fv,0);
                    UF2D_fsoutflow += UF2D_Advect_mass(dt,_dem,UF2D_ssm_D.at(d),_fu,_fv,0);

                }
            }
        }

        //advect properties
        UF2D_Advect_prop(dt,_dem,_f,_fu,_fv,UF2D_visc,0);
        if(UF_SOLIDPHASE)
        {
            UF2D_Advect_prop(dt,_dem,_s,_su,_sv,UF2D_d,0);
            UF2D_Advect_prop(dt,_dem,_s,_su,_sv,UF2D_ifa,0);
            UF2D_Advect_prop(dt,_dem,_s,_su,_sv,UF2D_rocksize,0);
        }

        UF_SWAP(_f,UF2D_fn,cTMap*);
        UF_SWAP(_fu,UF2D_fun,cTMap*);
        UF_SWAP(_fv,UF2D_fvn,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF_SWAP(_s,UF2D_sn,cTMap*);
            UF_SWAP(_su,UF2D_sun,cTMap*);
            UF_SWAP(_sv,UF2D_svn,cTMap*);
        }


        //first recalculate for the momentum source terms
        UF2D_FluidApplyMomentum(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn);
        UF_SWAP(_fu,UF2D_fun,cTMap*);
        UF_SWAP(_fv,UF2D_fvn,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF2D_SolidApplyMomentum(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sun,UF2D_svn);
            UF_SWAP(_su,UF2D_sun,cTMap*);
            UF_SWAP(_sv,UF2D_svn,cTMap*);
        }


    //if the more accurate but slightly slower boundary muscle scheme is set by the user, use this scheme
    }else if(UF_SCHEME == UF_SCHEME_BOUNDARYMUSCLE)
    {
        //advect momentum
        UF2D_Advect2_Momentum(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn,UF2D_sun,UF2D_svn,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2);

        //ThreadPool->WaitForAll();


        //advect mass
        UF2D_foutflow += UF2D_Advect2_mass(dt,_dem,_f,_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,UF2D_fn);
        if(UF_SOLIDPHASE)
        {
            UF2D_soutflow += UF2D_Advect2_mass(dt,_dem,_s,_s,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2,UF2D_sn);
        }

        if(SwitchErosion)
        {
            UF2D_fsoutflow += UF2D_Advect2_mass(dt,_dem,UF2D_blm,_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,0);
            UF2D_fsoutflow += UF2D_Advect2_mass(dt,_dem,UF2D_ssm,_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,0);
            if(SwitchUseGrainSizeDistribution)
            {
                FOR_GRAIN_CLASSES
                {
                    UF2D_fsoutflow += UF2D_Advect2_mass(dt,_dem,UF2D_blm_D.at(d),_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,0);
                    UF2D_fsoutflow += UF2D_Advect2_mass(dt,_dem,UF2D_ssm_D.at(d),_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,0);

                }
            }
        }

        //advect properties
        UF2D_Advect2_prop(dt,_dem,_f,_f,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,UF2D_visc,0);
        if(UF_SOLIDPHASE)
        {
            UF2D_Advect2_prop(dt,_dem,_s,_s,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2,UF2D_d,0);
            UF2D_Advect2_prop(dt,_dem,_s,_s,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2,UF2D_ifa,0);
            UF2D_Advect2_prop(dt,_dem,_s,_s,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2,UF2D_rocksize,0);
        }

        UF_SWAP(_f,UF2D_fn,cTMap*);
        UF_SWAP(_fu,UF2D_fun,cTMap*);
        UF_SWAP(_fv,UF2D_fvn,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF_SWAP(_s,UF2D_sn,cTMap*);
            UF_SWAP(_su,UF2D_sun,cTMap*);
            UF_SWAP(_sv,UF2D_svn,cTMap*);
        }

        //first recalculate for the momentum source terms
        UF2D_FluidApplyMomentum2(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn);
        UF_SWAP(_fu,UF2D_fun,cTMap*);
        UF_SWAP(_fv,UF2D_fvn,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF2D_SolidApplyMomentum2(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sun,UF2D_svn);
            UF_SWAP(_su,UF2D_sun,cTMap*);
            UF_SWAP(_sv,UF2D_svn,cTMap*);
        }

    }



}

void TWorld::UF1D_Source(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su)
{
    //fluid and solid mass sources are added here (currently none)
    UF1D_FluidSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fn);

    if(UF_SOLIDPHASE)
    {
        UF1D_SolidSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sn);
        UF_SWAP(_s,UF1D_sn,cTMap*);
    }

    UF_SWAP(_f,UF1D_fn,cTMap*);

    if(UF_SCHEME == UF_SCHEME_CENTRALSIMPLE)
    {
        //first recalculate for the momentum source terms
        UF1D_FluidMomentumSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun);

        if(UF_SOLIDPHASE)
        {
            UF1D_SolidMomentumSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sun);
        }

    }else if(UF_SCHEME == UF_SCHEME_BOUNDARYMUSCLE)
    {
        //first recalculate for the momentum source terms
        UF1D_FluidMomentum2Source(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun);

        if(UF_SOLIDPHASE)
        {
            UF1D_SolidMomentum2Source(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sun);
        }
    }

}
void TWorld::UF1D_Scheme(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su)
{

    //if the less accurate but slightly faster central simple scheme is set by the user, use this
    if(UF_SCHEME == UF_SCHEME_CENTRALSIMPLE)
    {
        //advect momentum
        UF1D_Advect_Momentum(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun,UF1D_sun);

        //advect mass
        UF1D_foutflow += UF1D_Advect_mass(dt,_ldd,_lddw,_lddh,UF1D_f,UF1D_fu,UF1D_fn);
        if(UF_SOLIDPHASE)
        {
            UF1D_soutflow += UF1D_Advect_mass(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_su,UF1D_sn);
        }

        if(SwitchErosion)
        {
            UF1D_fsoutflow += UF1D_Advect_mass(dt,_ldd,_lddw,_lddh,UF1D_blm,_fu,0);
            if(UF_SOLIDPHASE)
            {
                UF1D_fsoutflow += UF1D_Advect_mass(dt,_ldd,_lddw,_lddh,UF1D_ssm,_fu,0);
            }
            if(SwitchUseGrainSizeDistribution)
            {
                FOR_GRAIN_CLASSES
                {
                    UF1D_fsoutflow += UF1D_Advect_mass(dt,_ldd,_lddw,_lddh,UF1D_blm_D.at(d),_fu,0);
                    if(UF_SOLIDPHASE)
                    {
                        UF1D_fsoutflow += UF1D_Advect_mass(dt,_ldd,_lddw,_lddh,UF1D_ssm_D.at(d),_fu,0);
                    }
                }
            }
        }

        //advect properties
        UF1D_Advect_prop(dt,_ldd,_lddw,_lddh,UF1D_f,UF1D_fu,UF1D_visc,0);
        if(UF_SOLIDPHASE)
        {
            UF1D_Advect_prop(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_su,UF1D_d,0);
            UF1D_Advect_prop(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_su,UF1D_ifa,0);
            UF1D_Advect_prop(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_su,UF1D_rocksize,0);
        }

        UF_SWAP(_f,UF1D_fn,cTMap*);
        UF_SWAP(_fu,UF1D_fun,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF_SWAP(_s,UF1D_sn,cTMap*);
            UF_SWAP(_su,UF1D_sun,cTMap*);
        }

        //first recalculate for the momentum source terms
        UF1D_FluidApplyMomentum(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun);
        UF_SWAP(_fu,UF1D_fun,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF1D_SolidApplyMomentum(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sun);
            UF_SWAP(_su,UF1D_sun,cTMap*);
        }

    //if the slightly less fast but more accurate and stable boundary muscle solution is set by the user, use this
    }else if(UF_SCHEME == UF_SCHEME_BOUNDARYMUSCLE)
    {

        //advect momentum
        UF1D_Advect2_Momentum(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun,UF1D_sun,UF1D_fq1,UF1D_fq2,UF1D_sq1,UF1D_sq2);

        //advect mass
        UF1D_foutflow += UF1D_Advect2_mass(dt,_ldd,_lddw,_lddh,UF1D_f,UF1D_f,UF1D_fq1,UF1D_fq2,UF1D_fn);
        if(UF_SOLIDPHASE)
        {
            UF1D_soutflow += UF1D_Advect2_mass(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_s,UF1D_sq1,UF1D_sq2,UF1D_sn);
        }

        if(SwitchErosion)
        {
            UF1D_fsoutflow += UF1D_Advect2_mass(dt,_ldd,_lddw,_lddh,UF1D_blm,UF1D_f,UF1D_fq1,UF1D_fq2,0);
            UF1D_fsoutflow += UF1D_Advect2_mass(dt,_ldd,_lddw,_lddh,UF1D_ssm,UF1D_f,UF1D_fq1,UF1D_fq2,0);
            if(SwitchUseGrainSizeDistribution)
            {
                FOR_GRAIN_CLASSES
                {
                    UF1D_fsoutflow += UF1D_Advect2_mass(dt,_ldd,_lddw,_lddh,UF1D_blm_D.at(d),UF1D_f,UF1D_fq1,UF1D_fq2,0);
                    UF1D_fsoutflow += UF1D_Advect2_mass(dt,_ldd,_lddw,_lddh,UF1D_ssm_D.at(d),UF1D_f,UF1D_fq1,UF1D_fq2,0);

                }
            }
        }

        //advect properties
        UF1D_Advect2_prop(dt,_ldd,_lddw,_lddh,UF1D_f,UF1D_f,UF1D_fq1,UF1D_fq2,UF1D_visc,0);

        if(UF_SOLIDPHASE)
        {
            //UF1D_Advect2_prop(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_s,UF1D_sq1,UF1D_sq2,UF1D_d,0);
            UF1D_Advect2_prop(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_s,UF1D_sq1,UF1D_sq2,UF1D_ifa,0);
            UF1D_Advect2_prop(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_s,UF1D_sq1,UF1D_sq2,UF1D_rocksize,0);
        }

        UF_SWAP(_f,UF1D_fn,cTMap*);
        UF_SWAP(_fu,UF1D_fun,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF_SWAP(_s,UF1D_sn,cTMap*);
            UF_SWAP(_su,UF1D_sun,cTMap*);
        }

        //first recalculate for the momentum source terms
        UF1D_FluidApplyMomentum2(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun);
        UF_SWAP(_fu,UF1D_fun,cTMap*);

        if(UF_SOLIDPHASE)
        {
            UF1D_SolidApplyMomentum2(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sun);
            UF_SWAP(_su,UF1D_sun,cTMap*);
        }
    }



}


//General Function
void TWorld::UnifiedFlowSediment()
{
    //add any initial sources of sediment to the flow
    UF_SedimentSource(_dt);

    //solve flow detachment and add to flow
    UF_FlowDetachment(_dt);

    //solve flow entrainment (by solid phase) and add to flow
    UF_FlowEntrainment(_dt);

    //NOTE: transport is not performed here, but during the normal scheme were fluid and solid transport is also performed.
    //if any new substance needs transport, use the generic functions that are called there!!
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
        UF2D_q->Drc = UF2D_velocity->Drc * (UF2D_h->Drc) * _dx ;
        UF2D_qs->Drc = UF2D_q->Drc * UF2D_fsConc->Drc;

        UF1D_h->Drc = (UF1D_f->Drc + UF1D_s->Drc)/(_dx * UF1D_LDDw->Drc);
        UF1D_fsConc->Drc = (UF1D_f->Drc + UF1D_s->Drc) > UF_VERY_SMALL? (UF1D_ssm->Drc + UF1D_blm->Drc)/(UF1D_f->Drc + UF1D_s->Drc) : 0.0;
        UF1D_sConc->Drc = (UF1D_f->Drc + UF1D_s->Drc) > UF_VERY_SMALL?(UF1D_s->Drc * UF1D_d->Drc)/(UF1D_f->Drc + UF1D_s->Drc) : 0.0;
        UF1D_tConc->Drc = UF1D_fsConc->Drc + UF1D_sConc->Drc;
        UF1D_velocity->Drc = (UF1D_f->Drc + UF1D_s->Drc) > UF_VERY_SMALL? std::fabs((UF1D_f->Drc * UF1D_fu->Drc + UF1D_s->Drc * UF1D_su->Drc)/(UF1D_f->Drc + UF1D_s->Drc)) : 0.0;
        UF1D_q->Drc = UF1D_velocity->Drc * (UF1D_h->Drc) * _dx;
        UF1D_qs->Drc = UF1D_q->Drc * UF1D_fsConc->Drc;

        //return water height to the rest of OpenLisem
        if(ChannelAdj->Drc > 0)
        {
            WHrunoff->Drc = UF2D_f->Drc/(DX->Drc*ChannelAdj->Drc);
        }else
        {
            WHrunoff->Drc = 0;
        }
        Qn->Drc = UF2D_f->Drc * UF2D_velocity->Drc * _dx;
        Q->Drc = Qn->Drc;

        WHroad->Drc = WHrunoff->Drc;
        // set road to average outflowing wh, no surface storage.


        WH->Drc = WHrunoff->Drc + WHstore->Drc;
        // add new average waterlevel (A/dx) to stored water

        InfilVolKinWave->Drc = UF2D_Infiltration->Drc;

        UF2D_Infiltration->Drc = 0;

        V->Drc = UF2D_velocity->Drc;

        // recalc velocity for output to map, is not used in other processes
        WaterVolall->Drc = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        // is the same as :         WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);

        if(UF2D_h->Drc > UF_DISPLAYFLOODMINIMUM)
        {
            hmx->Drc = UF2D_h->Drc;
            UVflood->Drc = UF2D_velocity->Drc;
            if(floodTime->Drc == 0)
            {
                floodTimeStart->Drc = std::max((_dt/10.0)/60.0,(this->time/60.0));
            }
            floodTime->Drc += _dt/60.0;

            floodHmxMax->Drc = std::max(floodHmxMax->Drc,hmx->Drc);
            floodVMax->Drc = std::max(floodVMax->Drc,UVflood->Drc);
        }else
        {
            hmx->Drc = 0;
            UVflood->Drc = 0;
        }


        if(UF2D_sConc->Drc > UF_DISPLAYDEBRISFLOWMINIMUM)
        {
            dfhmx->Drc = UF2D_h->Drc;
            dfUV->Drc = UF2D_velocity->Drc;
            if(dfTime->Drc == 0)
            {
                dfTimeStart->Drc = std::max((_dt/10.0)/60.0,(this->time/60.0));
            }
            dfTime->Drc += _dt/60.0;

            dfHmxMax->Drc = std::max(dfHmxMax->Drc,dfhmx->Drc);
            dfVMax->Drc = std::max(dfVMax->Drc,dfUV->Drc);
        }else
        {
            dfhmx->Drc = 0;
            dfUV->Drc = 0;
        }

    }


}


void TWorld::UF2D1D_LaxNumericalCorrection(cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                                 cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                                 cTMap * _fu1D,cTMap * _s1D,
                                 cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                                 cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                                 cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                                 cTMap * _su2D,cTMap * _sv2D)
{

    //The Lax numerical correction is usefull to correct for artifacts of the courant factor limitation.
    //when velocity exceeds the courant limit, blocking appears. Using the Lax-scheme, values are averaged with surrounding cells while maintaining mass balance
    //in this implementation, flow barriers and the dem are taken into accaunt

    FOR_ROW_COL_UF2D_DT
    {
        double w = std::min(1.0,std::max(0.0,((std::max(std::fabs(2.5* _fu2D->Drc),std::max(std::fabs(2.5 * _fv2D->Drc),std::max(std::fabs(_su2D->Drc),std::fabs(_sv2D->Drc)))))/(UF_Courant * _dx) - 0.65)*0.1));
        UF2D_Test->Drc = w;

        if(w > 0.5)
        {
            double h = _f2D->data[r][c]; double u = _fu2D->data[r-1][c]; double v = _fv2D->data[r-1][c];
            double h1 = h,h2 = h,h3 = h,h4 = h,u1 = u,u2 = u,u3 = u,u4 = u,v1 = v,v2 = v,v3 = v,v4 = v;

            if(!UF_OUTORMV(_dem,r-1,c))
            {
                h1 = _f2D->data[r-1][c];
                u1 = _fu2D->data[r-1][c];
                v1 = _fv2D->data[r-1][c];
                double fn = (1.0*(5.0 - w) * _f2D->Drc + w *(h1))/5.0;
                double dif = ((fn - _f2D->Drc) >0? 1.0:-1.0 )*std::max((fn - _f2D->Drc) <0? UF2D_MaxFlux(_dem,_f2D,r,c,0,-1):UF2D_MaxFlux(_dem,_f2D,r,c-1,0,1),std::fabs(fn - _f2D->Drc));
                _f2D->data[r-1][c] = std::max(0.0,_f2D->data[r-1][c] - dif);
                double dif2 = h1 - _f2D->data[r-1][c];
                _f2D->Drc += dif2;
            }
            if(!UF_OUTORMV(_dem,r+1,c))
            {
                h2 = _f2D->data[r+1][c];
                u2 = _fu2D->data[r+1][c];
                v2 = _fv2D->data[r+1][c];
                double fn = (1.0*(5.0 - w) * _f2D->Drc + w *(h2))/5.0;
                double dif = ((fn - _f2D->Drc) >0? 1.0:-1.0 )*std::max((fn - _f2D->Drc) <0? UF2D_MaxFlux(_dem,_f2D,r,c,0,1):UF2D_MaxFlux(_dem,_f2D,r,c+1,0,-1),std::fabs(fn - _f2D->Drc));
                _f2D->data[r+1][c] = std::max(0.0,_f2D->data[r+1][c] - dif);
                double dif2 = h2 - _f2D->data[r+1][c];
                _f2D->Drc += dif2;
            }
            if(!UF_OUTORMV(_dem,r,c-1))
            {
                h3 = _f2D->data[r][c-1];
                u3 = _fu2D->data[r][c-1];
                v3 = _fv2D->data[r][c-1];
                double fn = (1.0*(5.0 - w) * _f2D->Drc + w *(h3))/5.0;
                double dif = ((fn - _f2D->Drc) >0? 1.0:-1.0 )*std::max((fn - _f2D->Drc) <0? UF2D_MaxFlux(_dem,_f2D,r,c,-1,0):UF2D_MaxFlux(_dem,_f2D,r-1,c,1,0),std::fabs(fn - _f2D->Drc));
                _f2D->data[r][c-1] = std::max(0.0,_f2D->data[r][c-1] - dif);
                double dif2 = h3 - _f2D->data[r][c-1];
                _f2D->Drc += dif2;
            }
            if(!UF_OUTORMV(_dem,r,c+1))
            {
                h4 = _f2D->data[r][c+1];
                u4 = _fu2D->data[r][c+1];
                v4 = _fv2D->data[r][c+1];
                double fn = (1.0*(5.0 - w) * _f2D->Drc + w *(h4))/5.0;
                double dif = ((fn - _f2D->Drc) >0? 1.0:-1.0 )*std::max((fn - _f2D->Drc) <0? UF2D_MaxFlux(_dem,_f2D,r,c,1,0):UF2D_MaxFlux(_dem,_f2D,r+1,c,-1,0),std::fabs(fn - _f2D->Drc));
                _f2D->data[r][c+1] = std::max(0.0,_f2D->data[r][c+1] - dif);
                double dif2 = h4 - _f2D->data[r][c+1];
                _f2D->Drc += dif2;
            }

            //_f2D->Drc = (4.0*(1.0 + w) * _f2D->Drc + w *(h1 + h2 + h3 + h4))/8.0;
            _fu2D->Drc = (4.0*(1.0 + w) * _fu2D->Drc + w *(u1 + u2 + u3 + u4))/8.0;
            _fv2D->Drc = (4.0*(1.0 + w) * _fv2D->Drc + w *(v1 + v2 + v3 + v4))/8.0;

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
            if(w > 0.5)
            {
                double sh = _s2D->data[r][c]; double su = _su2D->data[r-1][c]; double sv = _sv2D->data[r-1][c];
                double sh1 = sh, sh2 = sh, sh3 = sh, sh4 = sh,su1 =su,su2 = su,su3 = su,su4 = su,sv1 = sv,sv2 = sv,sv3 = sv,sv4 = sv;

                if(!UF_OUTORMV(_dem,r-1,c))
                {

                    sh1 = _s2D->data[r-1][c];
                    su1 = _su2D->data[r-1][c];
                    sv1 = _sv2D->data[r-1][c];
                    double sn = (1.0*(5.0 - w) * _s2D->Drc + w *(sh1))/5.0;
                    double dif = ((sn - _s2D->Drc) >0? 1.0:-1.0 )*std::max((sn - _s2D->Drc) <0? UF2D_MaxFlux(_dem,_s2D,r,c,0,-1):UF2D_MaxFlux(_dem,_s2D,r,c-1,0,1),std::fabs(sn - _s2D->Drc));
                    _s2D->data[r-1][c] = std::max(0.0,_s2D->data[r-1][c] - dif);
                    double dif2 = sh1 - _s2D->data[r-1][c];
                    _s2D->Drc += dif2;

                }
                if(!UF_OUTORMV(_dem,r+1,c))
                {
                    sh2 = _s2D->data[r+1][c];
                    su2 = _su2D->data[r+1][c];
                    sv2 = _sv2D->data[r+1][c];
                    double sn = (1.0*(5.0 - w) * _s2D->Drc + w *(sh2))/5.0;
                    double dif = ((sn - _s2D->Drc) >0? 1.0:-1.0 )*std::max((sn - _s2D->Drc) <0? UF2D_MaxFlux(_dem,_s2D,r,c,0,1):UF2D_MaxFlux(_dem,_s2D,r,c+1,0,-1),std::fabs(sn - _s2D->Drc));
                    _s2D->data[r+1][c] = std::max(0.0,_s2D->data[r+1][c] - dif);
                    double dif2 = sh2 - _s2D->data[r+1][c];
                    _s2D->Drc += dif2;
                }
                if(!UF_OUTORMV(_dem,r,c-1))
                {
                    sh3 = _s2D->data[r][c-1];
                    su3 = _su2D->data[r][c-1];
                    sv3 = _sv2D->data[r][c-1];
                    double sn = (1.0*(5.0 - w) * _s2D->Drc + w *(sh3))/5.0;
                    double dif = ((sn - _s2D->Drc) >0? 1.0:-1.0 )*std::max((sn - _s2D->Drc) <0? UF2D_MaxFlux(_dem,_s2D,r,c,-1,0):UF2D_MaxFlux(_dem,_s2D,r-1,c,1,0),std::fabs(sn - _s2D->Drc));
                    _s2D->data[r][c-1] = std::max(0.0,_s2D->data[r][c-1] - dif);
                    double dif2 = sh3 - _s2D->data[r][c-1];
                    _s2D->Drc += dif2;
                }
                if(!UF_OUTORMV(_dem,r,c+1))
                {
                    sh4 = _s2D->data[r][c+1];
                    su4 = _su2D->data[r][c+1];
                    sv4 = _sv2D->data[r][c+1];
                    double sn = (1.0*(5.0 - w) * _s2D->Drc + w *(sh3))/5.0;
                    double dif = ((sn - _s2D->Drc) >0? 1.0:-1.0 )*std::max((sn - _s2D->Drc) <0? UF2D_MaxFlux(_dem,_s2D,r,c,1,0):UF2D_MaxFlux(_dem,_s2D,r+1,c,-1,0),std::fabs(sn - _s2D->Drc));
                    _s2D->data[r][c+1] = std::max(0.0,_s2D->data[r][c+1] - dif);
                    double dif2 = sh4 - _s2D->data[r][c+1];
                    _s2D->Drc += dif2;
                }


                _su2D->Drc = (4.0*(1.0 + w) * _su2D->Drc + w *(su1 + su2 + su3 + su4))/8.0;
                _sv2D->Drc = (4.0*(1.0 + w) * _sv2D->Drc + w *(sv1 + sv2 + sv3 + sv4))/8.0;

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

}
