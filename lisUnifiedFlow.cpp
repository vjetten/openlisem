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

void TWorld::UF_Init()
{

    //constants
    UF_Courant = getvaluedouble("Surface Flow Courant Factor");
    UF_Aspect = 0.1;
    UF_Chi = 3;
    UF_Ksi = 3;
    UF_j = 2;
    UF_Gravity = 9.81;
    UF_GravitySqrt = std::sqrt(9.81);
    UF2D_MinimumDT = getvaluedouble("Surface Flow Minimum Timestep");
    UF1D_MinimumDT = getvaluedouble("Channel Flow Minimum Timestep");
    UF_SigmaDiffusion = 1.0;
    UF_MANNINGCOEFFICIENT = 0.05;

    UF_Alpha_DV  =0.1;
    UF_Beta_DV = 20.0;

    UF_Alpha_YS = 0.1;
    UF_Beta_YS = 20.0;

    UF_Alpha_DR = 0.0538;
    UF_Beta_DR = 6.0896;

    UF_Intersect_K = 30;
    UF_Slope_K = 100000.0;

    //internal use
    UF_1DACTIVE = SwitchIncludeChannel;
    UF_SCHEME = UF_SCHEME_BOUNDARYMUSCLE;
    UF_DTMIN = 0;
    UF_SOLIDPHASE = false;

    UF2D_Test = NewMap(0.0);

    //just for display
    UF2D_h = NewMap(0.0);
    UF2D_fsConc = NewMap(0.0);
    UF2D_sConc = NewMap(0.0);
    UF2D_tConc = NewMap(0.0);
    UF2D_velocity = NewMap(0.0);
    UF2D_u = NewMap(0.0);
    UF2D_v = NewMap(0.0);
    UF2D_q = NewMap(0.0);
    UF2D_qs = NewMap(0.0);

    UF1D_h = NewMap(0.0);
    UF1D_fsConc = NewMap(0.0);
    UF1D_sConc = NewMap(0.0);
    UF1D_tConc = NewMap(0.0);
    UF1D_velocity = NewMap(0.0);
    UF1D_q = NewMap(0.0);
    UF1D_qs = NewMap(0.0);

    //internal slope functions
    UF2D_Slope = NewMap(0.0);
    UF2D_SlopeX = NewMap(0.0);
    UF2D_SlopeY = NewMap(0.0);
    UF1D_LDDs = NewMap(0.0);

    //actual calculation variables
    ////2D
    UF2D_DEM = NewMap(0.0);
    copy(*UF2D_DEM,*DEM);
    UF2D_T = NewMap(0.0);
    UF2D_DT = NewMap(0.0);
    UF2D_DTStep = NewMap(0.0);
    UF2D_CellR = NewMap(0.0);
    UF2D_CellC = NewMap(0.0);
    UF2D_Courant = NewMap(0.0);

    //fluid phase
    UF2D_f = NewMap(0.0);
    UF2D_visc = NewMap(1.0);
    UF2D_fu = NewMap(0.0);
    UF2D_fv = NewMap(0.0);
    UF2D_fax = NewMap(0.0);
    UF2D_fay = NewMap(0.0);
    UF2D_fqx1 = NewMap(0.0);
    UF2D_fqy1 = NewMap(0.0);
    UF2D_fqx2 = NewMap(0.0);
    UF2D_fqy2 = NewMap(0.0);
    UF2D_ssm = NewMap(0.0);
    UF2D_blm = NewMap(0.0);
    UF2D_sstc = NewMap(0.0);
    UF2D_bltc = NewMap(0.0);
    UF2D_fsc = NewMap(0.0);
    UF2D_fsd = NewMap(0.0);

    //solid phase
    UF2D_s = NewMap(0.0);
    UF2D_d = NewMap(2000.0);
    UF2D_ifa = NewMap(0.3);
    UF2D_rocksize = NewMap(0.1);
    UF2D_su = NewMap(0.0);
    UF2D_sv = NewMap(0.0);
    UF2D_sax = NewMap(0.0);
    UF2D_say = NewMap(0.0);
    UF2D_sqx1 = NewMap(0.0);
    UF2D_sqy1 = NewMap(0.0);
    UF2D_sqx2 = NewMap(0.0);
    UF2D_sqy2 = NewMap(0.0);


    //for new timestep
    //fluid phase
    UF2D_fn = NewMap(0.0);
    UF2D_fun = NewMap(0.0);
    UF2D_fvn = NewMap(0.0);
    //solid phase
    UF2D_sn = NewMap(0.0);
    UF2D_sun = NewMap(0.0);
    UF2D_svn = NewMap(0.0);

    ////1D
    UF1D_LDD = NewMap(0.0);
    UF1D_LDDw = NewMap(0.0);
    UF1D_LDDh = NewMap(0.0);
    UF1D_Slope = NewMap(0.0);
    UF1D_LDD->setAllMV();
    if(SwitchIncludeChannel)
    {
        FOR_ROW_COL_MV_CH
        {
            UF1D_LDD->Drc = LDDChannel->Drc;
            UF1D_LDDw->Drc = ChannelWidth->Drc;
            UF1D_Slope->Drc = ChannelGrad->Drc;
        }
    }
    UF1D_T = NewMap(0.0);
    UF1D_DT = NewMap(0.0);

    UF1D_DTStep = NewMap(0.0);
    UF1D_Courant = NewMap(0.0);
    //fluid phase
    UF1D_f = NewMap(0.0);
    UF1D_fstore = NewMap(0.0);
    UF1D_visc = NewMap(1.0);
    UF1D_fu = NewMap(0.0);
    UF1D_fa = NewMap(0.0);
    UF1D_fq1 = NewMap(0.0);
    UF1D_fq2 = NewMap(0.0);
    UF1D_ssm = NewMap(0.0);
    UF1D_blm = NewMap(0.0);
    UF1D_bltc = NewMap(0.0);
    UF1D_sstc = NewMap(0.0);
    UF1D_fsc = NewMap(0.0);
    UF1D_fsd = NewMap(0.0);

    //solid phase
    UF1D_sstore = NewMap(0.0);
    UF1D_s = NewMap(0.0);
    UF1D_d = NewMap(2000.0);
    UF1D_ifa = NewMap(0.3);
    UF1D_rocksize = NewMap(0.1);
    UF1D_su = NewMap(0.0);
    UF1D_sa = NewMap(0.0);
    UF1D_sq1 = NewMap(0.0);
    UF1D_sq2 = NewMap(0.0);

    //for new timestep
    //fluid phase
    UF1D_fn = NewMap(0.0);
    UF1D_fun = NewMap(0.0);
    //solid phase
    UF1D_sn = NewMap(0.0);
    UF1D_sun = NewMap(0.0);

    //Multiclass sediment maps
    UF2D_ssm_D.clear();
    UF1D_ssm_D.clear();
    UF2D_blm_D.clear();
    UF1D_blm_D.clear();

    //Multiclass sediment maps
    UF2D_sstc_D.clear();
    UF1D_sstc_D.clear();
    UF2D_bltc_D.clear();
    UF1D_bltc_D.clear();

    if(SwitchUseGrainSizeDistribution)
    {
        FOR_GRAIN_CLASSES
        {
            UF2D_ssm_D.append(NewMap(0.0));
            UF1D_ssm_D.append(NewMap(0.0));
            UF2D_blm_D.append(NewMap(0.0));
            UF1D_blm_D.append(NewMap(0.0));

            UF2D_sstc_D.append(NewMap(0.0));
            UF1D_sstc_D.append(NewMap(0.0));
            UF2D_bltc_D.append(NewMap(0.0));
            UF1D_bltc_D.append(NewMap(0.0));
        }
    }

    UF1D_Dep = NewMap(0.0);
    UF1D_Det = NewMap(0.0);
    UF2D_Dep = NewMap(0.0);
    UF2D_Det = NewMap(0.0);

    UF2D_Infiltration = NewMap(0.0);
    UF1D_Infiltration = NewMap(0.0);

    //temporary maps for generic advection functions
    UF_t1 = NewMap(0.0);
    UF_t2 = NewMap(0.0);
    UF_t3 = NewMap(0.0);
    UF_t4 = NewMap(0.0);
    UF_t5 = NewMap(0.0);

    UF2D_MUSCLE_1_x1 = NewMap(0.0);
    UF2D_MUSCLE_1_x2 = NewMap(0.0);
    UF2D_MUSCLE_1_y1 = NewMap(0.0);
    UF2D_MUSCLE_1_y2 = NewMap(0.0);
    UF2D_MUSCLE_2_x1 = NewMap(0.0);
    UF2D_MUSCLE_2_x2 = NewMap(0.0);
    UF2D_MUSCLE_2_y1 = NewMap(0.0);
    UF2D_MUSCLE_2_y2 = NewMap(0.0);
    UF2D_MUSCLE_3_x1 = NewMap(0.0);
    UF2D_MUSCLE_3_x2 = NewMap(0.0);
    UF2D_MUSCLE_3_y1 = NewMap(0.0);
    UF2D_MUSCLE_3_y2 = NewMap(0.0);
    UF2D_MUSCLE_OUT_x1 = NewMap(0.0);
    UF2D_MUSCLE_OUT_x2 = NewMap(0.0);
    UF2D_MUSCLE_OUT_y1 = NewMap(0.0);
    UF2D_MUSCLE_OUT_y2 = NewMap(0.0);

    UF1D_MUSCLE_1_x1 = NewMap(0.0);
    UF1D_MUSCLE_1_x2 = NewMap(0.0);
    UF1D_MUSCLE_2_x1 = NewMap(0.0);
    UF1D_MUSCLE_2_x2 = NewMap(0.0);
    UF1D_MUSCLE_3_x1 = NewMap(0.0);
    UF1D_MUSCLE_3_x2 = NewMap(0.0);
    UF1D_MUSCLE_OUT_x1 = NewMap(0.0);
    UF1D_MUSCLE_OUT_x2 = NewMap(0.0);

}

//General Function
void TWorld::UnifiedFlow()
{

    //set input from the rest of the OpenLisem model
    UF_SetInput();

    ////START ALGORITHM
    ////from now on all input and output is provided as funciton arguments
    ////This increases re-usablitiy of the code (use another DEM, Velocity map.. etc.. and everything will remain functional)

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
        ////TOPOGRAPHY ANALYSIS
        UF_DEMLDDAnalysis(UF2D_DEM,UF1D_LDD,UF1D_LDDw,UF1D_LDDh,UF1D_f,UF1D_s,UF2D_f,UF2D_s);


        ////TIMESTEP ANALYSIS
        //analyzes spatially dynamic timstep that should be made
        dt = UF_TimeStep(t,     UF2D_DEM,                                                    //dem info
                                UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                                UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                                UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                                UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                                UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv,        //2d solid phase
                                UF2D_T,UF2D_DT,UF2D_DTStep,UF1D_T,UF1D_DT,UF1D_DTStep);      //output timesteps

        UF2D_Source(UF2D_DT,UF2D_DEM,UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);

        if(UF_1DACTIVE)
        {
            UF1D_Source(UF1D_DT,UF1D_LDD,UF1D_LDDw,UF1D_LDDh,UF1D_f,UF1D_visc,UF1D_fu,UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su);
        }

        ////2D SCHEME
        UF2D_Scheme(UF2D_DT,UF2D_DEM,UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);

        if(UF_1DACTIVE)
        {
            ////1D SCHEME
            UF1D_Scheme(UF1D_DT,UF1D_LDD,UF1D_LDDw,UF1D_LDDh,UF1D_f,UF1D_visc,UF1D_fu,UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su);

            ////CONNECTIONS
            UF2D1D_Connection(UF2D_DT,     UF2D_DEM,                                        //dem info
                               UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                               UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                               UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                               UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                               UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase
        }

        /////INFILTRATION
        UF2D1D_Infiltration(UF2D_DT,     UF2D_DEM,                                        //dem info
                            UF1D_LDD,UF1D_LDDw,UF1D_LDDh,                                //channel info
                            UF1D_f,UF1D_visc,UF1D_fu,                                    //1d fluid phase
                            UF1D_s,UF1D_d,UF1D_ifa,UF1D_rocksize,UF1D_su,                //1d solid phase
                            UF2D_f,UF2D_visc,UF2D_fu,UF2D_fv,                            //2d fluid phase
                            UF2D_s,UF2D_d,UF2D_ifa,UF2D_rocksize,UF2D_su,UF2D_sv);       //2d solid phase

        //increase timer by the current timstep!
        t = t + dt;

        DEBUG(QString("UF Step t: %1  dt: %2").arg(t).arg(dt));
    }

    //again uses non-functionparameter variables

    //soil interactions
    //(sediment transport and solid phase transport is done together with the fluid equations since these are completely integrated)
    if(SwitchErosion)
    {
        UnifiedFlowSediment();
    }

    //set output maps for display etc..
    UF_SetOutput();
}

#define UF_SWAP(x, y, T) do { MaskedRaster<double> SWAP = x->data; x->data = y->data; y->data = SWAP; } while (0)

double TWorld::UF2D_Source(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv)
{
    UF2D_FluidSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fn);
    UF_SWAP(_f,UF2D_fn,cTMap*);

    if(UF_SOLIDPHASE)
    {
        UF2D_SolidSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sn);
        UF_SWAP(_s,UF2D_sn,cTMap*);
    }

    //first recalculate for the momentum source terms
    UF2D_FluidMomentumSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn);

    if(UF_SOLIDPHASE)
    {
        UF2D_SolidMomentumSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sun,UF2D_svn);
    }

}

//2D version
double TWorld::UF2D_Scheme(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv)
{

    FOR_ROW_COL_UF2D_DT
    {
        UF2D_Test->Drc += 1;
    }}}

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
    }else if(UF_SCHEME == UF_SCHEME_BOUNDARYMUSCLE)
    {
        //advect momentum
        UF2D_Advect2_Momentum(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn,UF2D_sun,UF2D_svn,UF2D_fqx1,UF2D_fqx2,UF2D_fqy1,UF2D_fqy2,UF2D_sqx1,UF2D_sqx2,UF2D_sqy1,UF2D_sqy2);

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
    }


}

void TWorld::UF1D_Source(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su)
{
    UF1D_FluidSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fn);

    if(UF_SOLIDPHASE)
    {
        UF1D_SolidSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sn);
        UF_SWAP(_s,UF1D_sn,cTMap*);
    }

    UF_SWAP(_f,UF1D_fn,cTMap*);

    //first recalculate for the momentum source terms
    UF1D_FluidMomentumSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun);

    if(UF_SOLIDPHASE)
    {
        UF1D_SolidMomentumSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sun);
    }

}
void TWorld::UF1D_Scheme(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su)
{

    //first recalculate for the momentum source terms
    UF1D_FluidApplyMomentum(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun);
    UF_SWAP(_fu,UF1D_fun,cTMap*);

    if(UF_SOLIDPHASE)
    {
        UF1D_SolidApplyMomentum(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sun);
        UF_SWAP(_su,UF1D_sun,cTMap*);
    }


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
    }


}


//General Function
void TWorld::UnifiedFlowSediment()
{
    UF_FlowDetachment(_dt);

    UF_FlowEntrainment(_dt);

}
//set output
void TWorld::UF_SetInput()
{
    cTMap * _dem = UF2D_DEM;
    cTMap * _ldd = UF1D_LDD;
    FOR_ROW_COL_UF2D
    {
        UF2D_Test->Drc = 0;
        UF2D_f->Drc = WHrunoff->Drc * FlowWidth->Drc * _dx;
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
            //UF1D_f->Drc += 1;
            //UF1D_s->Drc += 0.5;
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
            WHrunoff->Drc = UF2D_f->Drc/(_dx*ChannelAdj->Drc);
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
        UF2D_Test->Drc = UF2D_Infiltration->Drc;
        UF2D_Infiltration->Drc = 0;

        V->Drc = UF2D_velocity->Drc;

        // recalc velocity for output to map, is not used in other processes
        WaterVolall->Drc = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        // is the same as :         WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
    }


}


