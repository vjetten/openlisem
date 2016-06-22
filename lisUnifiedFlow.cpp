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
    UF_Courant = 0.25;
    UF_Aspect = 0.1;
    UF_Chi = 3;
    UF_Ksi = 3;
    UF_j = 1;
    UF_Gravity = 9.81;
    UF_GravitySqrt = std::sqrt(9.81);
    UF2D_MinimumDT = 0.5;
    UF1D_MinimumDT = 0.5;

    //internal use
    UF_1DACTIVE = false;
    UF_DTMIN = 0;

    //just for display
    UF2D_h = NewMap(0.0);
    UF2D_fsConc = NewMap(0.0);
    UF2D_sConc = NewMap(0.0);
    UF2D_tConc = NewMap(0.0);
    UF2D_velocity = NewMap(0.0);
    UF2D_u = NewMap(0.0);
    UF2D_v = NewMap(0.0);

    UF1D_h = NewMap(0.0);
    UF1D_fsConc = NewMap(0.0);
    UF1D_sConc = NewMap(0.0);
    UF1D_tConc = NewMap(0.0);
    UF1D_velocity = NewMap(0.0);

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
    cTMap * _dem = UF2D_DEM;
    FOR_ROW_COL_UF2D
    {
        if( r < 150)
        {
            UF2D_f->Drc = 0.0;
        }
    }
    UF2D_visc = NewMap(1.0);
    UF2D_fu = NewMap(0.0);
    UF2D_fv = NewMap(0.0);
    UF2D_ssm = NewMap(0.0);
    UF2D_blm = NewMap(0.0);
    UF2D_sstc = NewMap(0.0);
    UF2D_bltc = NewMap(0.0);
    UF2D_fsc = NewMap(0.0);
    UF2D_fsd = NewMap(0.0);

    //solid phase
    UF2D_sm = NewMap(0.0);
    UF2D_s = NewMap(0.0);
    UF2D_d = NewMap(2000.0);
    UF2D_ifa = NewMap(0.3);
    UF2D_rocksize = NewMap(0.1);
    UF2D_su = NewMap(0.0);
    UF2D_sv = NewMap(0.0);

    //for new timestep
    //fluid phase
    UF2D_fn = NewMap(0.0);
    UF2D_viscn = NewMap(0.0);
    UF2D_fun = NewMap(0.0);
    UF2D_fvn = NewMap(0.0);
    //solid phase
    UF2D_smn = NewMap(0.0);
    UF2D_sn = NewMap(0.0);
    UF2D_dn = NewMap(0.0);
    UF2D_ifan = NewMap(0.0);
    UF2D_rocksizen = NewMap(0.0);
    UF2D_sun = NewMap(0.0);
    UF2D_svn = NewMap(0.0);

    ////1D
    UF1D_LDD = NewMap(0.0);
    UF1D_LDDw = NewMap(0.0);
    UF1D_LDDh = NewMap(0.0);
    UF1D_LDD->setAllMV();
    FOR_ROW_COL_MV_CH
    {
        UF1D_LDD->Drc = LDD->Drc;
        UF1D_LDDw->Drc = ChannelWidth->Drc;
    }
    UF1D_T = NewMap(0.0);
    UF1D_DT = NewMap(0.0);
    UF1D_Slope = NewMap(0.1);
    UF1D_DTStep = NewMap(0.0);
    UF1D_Courant = NewMap(0.0);
    //fluid phase
    UF1D_f = NewMap(0.0);
    UF1D_visc = NewMap(0.0);
    UF1D_fu = NewMap(0.0);
    UF1D_ssm = NewMap(0.0);
    UF1D_blm = NewMap(0.0);
    UF1D_bltc = NewMap(0.0);
    UF1D_sstc = NewMap(0.0);
    UF1D_fsc = NewMap(0.0);
    UF1D_fsd = NewMap(0.0);

    //solid phase
    UF1D_sm = NewMap(0.0);
    UF1D_s = NewMap(0.0);
    UF1D_d = NewMap(0.0);
    UF1D_ifa = NewMap(0.0);
    UF1D_rocksize = NewMap(0.0);
    UF1D_su = NewMap(0.0);

    //for new timestep
    //fluid phase
    UF1D_fn = NewMap(0.0);
    UF1D_viscn = NewMap(0.0);
    UF1D_fun = NewMap(0.0);
    //solid phase
    UF1D_smn = NewMap(0.0);
    UF1D_sn = NewMap(0.0);
    UF1D_dn = NewMap(0.0);
    UF1D_ifan = NewMap(0.0);
    UF1D_rocksizen = NewMap(0.0);
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

    //temporary maps for generic advection functions
    UF_t1 = NewMap(0.0);
    UF_t2 = NewMap(0.0);
    UF_t3 = NewMap(0.0);

    UF_MUSCLE_1_N = NewMap(0.0);
    UF_MUSCLE_1_E = NewMap(0.0);
    UF_MUSCLE_1_S = NewMap(0.0);
    UF_MUSCLE_1_W = NewMap(0.0);
    UF_MUSCLE_2_N = NewMap(0.0);
    UF_MUSCLE_2_E = NewMap(0.0);
    UF_MUSCLE_2_S = NewMap(0.0);
    UF_MUSCLE_2_W = NewMap(0.0);
    UF_MUSCLE_OUT_N = NewMap(0.0);
    UF_MUSCLE_OUT_E = NewMap(0.0);
    UF_MUSCLE_OUT_S = NewMap(0.0);
    UF_MUSCLE_OUT_W = NewMap(0.0);

}

//General Function
void TWorld::UnifiedFlow()
{

    //set input from the rest of the OpenLisem model
    UF_SetInput();

    //soil interactions
    if(SwitchErosion)
    {
        UnifiedFlowSediment();
    }

    ////START ALGORITHM
    ////from now on all input and output is provided as funciton arguments
    ////This increases re-usablitiy of the code

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
    ////from now on all input and output is provided as funciton arguments
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
                                UF2D_T,UF2D_DT,UF2D_DTStep,UF1D_T,UF1D_DT,UF1D_DTStep);              //output timesteps

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
        //increase timer by the current timstep!
        t = t + dt;

        DEBUG(QString("UF Step t: %1  dt: %2").arg(t,dt));
    }

   //again uses non-functionparameter variables to set output maps
   UF_SetOutput();
}

#define UF_SWAP(x, y, T) do { T SWAP = x; x = y; y = SWAP; } while (0)

//2D version
double TWorld::UF2D_Scheme(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv)
{

    UF2D_FluidSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fn,UF2D_viscn);
    UF2D_SolidSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sn,UF2D_dn,UF2D_ifan,UF2D_rocksizen);

    copy(*_f,*UF2D_fn);
    copy(*_s,*UF2D_sn);
    copy(*_visc,*UF2D_viscn);
    copy(*_d,*UF2D_dn);
    copy(*_ifa,*UF2D_ifan);
    copy(*_rocksize,*UF2D_rocksizen);

    //first recalculate for the momentum source terms
    UF2D_FluidMomentumSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn);
    UF2D_SolidMomentumSource(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_sun,UF2D_svn);

    copy(*_fu,*UF2D_fun);
    copy(*_fv,*UF2D_fvn);
    copy(*_su,*UF2D_sun);
    copy(*_sv,*UF2D_svn);

    //advect momentum
    UF2D_Advect_Momentum(dt,_dem,_f,_visc,_fu,_fv,_s,_d,_ifa,_rocksize,_su,_sv,UF2D_fun,UF2D_fvn,UF2D_sun,UF2D_svn);

    //advect mass
    UF2D_foutflow += UF2D_Advect_mass(dt,_dem,_f,_fu,_fv,UF2D_fn);
    UF2D_soutflow += UF2D_Advect_mass(dt,_dem,_s,_su,_sv,UF2D_sn);

    //advect properties
    UF2D_Advect_prop(dt,_dem,_f,_fu,_fv,UF2D_visc,0);
    UF2D_Advect_prop(dt,_dem,_s,_su,_sv,UF2D_d,0);
    UF2D_Advect_prop(dt,_dem,_s,_su,_sv,UF2D_ifa,0);
    UF2D_Advect_prop(dt,_dem,_s,_su,_sv,UF2D_rocksize,0);

    copy(*_f,*UF2D_fn);
    copy(*_s,*UF2D_sn);

    copy(*_fu,*UF2D_fun);
    copy(*_fv,*UF2D_fvn);
    copy(*_su,*UF2D_sun);
    copy(*_sv,*UF2D_svn);

}

void TWorld::UF1D_Scheme(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su)
{
    UF1D_FluidSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fn,UF1D_viscn);
    UF1D_SolidSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sn,UF1D_dn,UF1D_ifan,UF1D_rocksizen);

    copy(*_f,*UF1D_fn);
    copy(*_s,*UF1D_sn);
    copy(*_visc,*UF1D_viscn);
    copy(*_d,*UF1D_dn);
    copy(*_ifa,*UF1D_ifan);
    copy(*_rocksize,*UF1D_rocksizen);

    //first recalculate for the momentum source terms
    UF1D_FluidMomentumSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun);
    UF1D_SolidMomentumSource(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_sun);

    copy(*_fu,*UF1D_fun);
    copy(*_su,*UF1D_sun);

    //advect momentum
    UF1D_Advect_Momentum(dt,_ldd,_lddw,_lddh,_f,_visc,_fu,_s,_d,_ifa,_rocksize,_su,UF1D_fun,UF1D_sun);

    //advect mass
    UF1D_foutflow += UF1D_Advect_mass(dt,_ldd,_lddw,_lddh,UF1D_f,UF1D_fu,UF1D_fn);
    UF1D_soutflow += UF1D_Advect_mass(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_su,UF1D_sn);

    //advect properties
    UF1D_Advect_prop(dt,_ldd,_lddw,_lddh,UF1D_f,UF1D_fu,UF1D_visc,0);
    UF1D_Advect_prop(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_su,UF1D_d,0);
    UF1D_Advect_prop(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_su,UF1D_ifa,0);
    UF1D_Advect_prop(dt,_ldd,_lddw,_lddh,UF1D_s,UF1D_su,UF1D_rocksize,0);

    copy(*_f,*UF1D_fn);
    copy(*_s,*UF1D_sn);

    copy(*_fu,*UF1D_fun);
    copy(*_su,*UF1D_sun);
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
    FOR_ROW_COL_UF2D
    {
        UF2D_f->Drc = WHrunoff->Drc * FlowWidth->Drc * _dx;
    }
    UF2D_foutflow = 0;
    UF2D_soutflow = 0;
    UF1D_foutflow = 0;
    UF1D_soutflow = 0;
}

//set output
void TWorld::UF_SetOutput()
{
    cTMap * _dem = UF2D_DEM;
    FOR_ROW_COL_UF2D
    {
        //just for display
        UF2D_h->Drc = (UF2D_f->Drc + UF2D_s->Drc)/(_dx*_dx);
        UF2D_fsConc->Drc = (UF2D_ssm->Drc + UF2D_blm->Drc)/(UF2D_f->Drc + UF2D_s->Drc);
        UF2D_sConc->Drc = (UF2D_s->Drc * UF2D_d->Drc)/(UF2D_f->Drc + UF2D_s->Drc);
        UF2D_tConc->Drc = UF2D_fsConc->Drc + UF2D_sConc->Drc;
        UF2D_u->Drc = (UF2D_f->Drc * UF2D_fu->Drc + UF2D_s->Drc * UF2D_su->Drc)/(UF2D_f->Drc + UF2D_s->Drc);
        UF2D_v->Drc = (UF2D_f->Drc * UF2D_fv->Drc + UF2D_s->Drc * UF2D_sv->Drc)/(UF2D_f->Drc + UF2D_s->Drc);
        UF2D_velocity->Drc = std::sqrt(UF2D_u->Drc * UF2D_u->Drc + UF2D_v->Drc * UF2D_v->Drc);

        UF1D_h->Drc = (UF1D_f->Drc + UF1D_s->Drc)/(_dx * UF1D_LDDw->Drc);
        UF1D_fsConc->Drc = (UF1D_ssm->Drc + UF1D_blm->Drc)/(UF1D_f->Drc + UF1D_s->Drc);
        UF1D_sConc->Drc = (UF1D_s->Drc * UF1D_d->Drc)/(UF1D_f->Drc + UF1D_s->Drc);
        UF1D_tConc->Drc = UF1D_fsConc->Drc + UF1D_sConc->Drc;
        UF1D_velocity->Drc = (UF1D_f->Drc * UF1D_fu->Drc + UF1D_s->Drc * UF1D_su->Drc)/(UF1D_f->Drc + UF1D_s->Drc);

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

        V->Drc = UF2D_velocity->Drc;

        // recalc velocity for output to map, is not used in other processes

        WaterVolall->Drc = WHrunoff->Drc*ChannelAdj->Drc*DX->Drc + DX->Drc*WHstore->Drc*SoilWidthDX->Drc;
        // is the same as :         WaterVolall->Drc = DX->Drc*( WH->Drc*SoilWidthDX->Drc + WHroad->Drc*RoadWidthDX->Drc);
    }


}


