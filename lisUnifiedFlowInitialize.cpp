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
    UF_Courant = 0.25; getvaluedouble("Surface Flow Courant Factor");
    UF_Aspect = 1.0;
    UF_Chi = 3.0;
    UF_Ksi = 3.0;
    UF_j = 2.0;
    UF_Gravity = 9.81;
    UF_GravitySqrt = std::sqrt(9.81);
    UF2D_MinimumDT = getvaluedouble("Flow Minimum Timestep");
    UF1D_MinimumDT = getvaluedouble("Flow Minimum Timestep");
    UF_DISPLAYFLOODMINIMUM = getvaluedouble("Minimal Flood Water Depth");
    UF_DISPLAYDEBRISFLOWMINIMUM = getvaluedouble("Minimum Debris Flow Volumetric Sediment Fraction");
    UF_SigmaDiffusion = 1.0;
    UF_ENTRAINMENTCONSTANT = getvaluedouble("Entrainment Coefficient");
    UF_MAXSOLIDCONCENTRATION = 0.6;
    UF_MINIMUMENTRAINMENTHEIGHT = getvaluedouble("Minimum Entrainment Height");
    UF_MANNINGCOEFFICIENT_FLUID = 0.2;
    UF_MANNINGCOEFFICIENT_SOLID = 0.05;
    UF_FrictionIterations = 1;

    UF_MAX_NUM_VEL = 50;

    UF_NRA = 150000;

    UF_Alpha_DV  =1.0;
    UF_Beta_DV = 5.0;

    UF_Alpha_YS = 0.1;
    UF_Beta_YS = 20.0;

    UF_Alpha_DR = 0.0538;
    UF_Beta_DR = 6.0896;

    UF_Intersect_K = 30;
    UF_Slope_K = 100000.0;

    //internal use
    UF_1DACTIVE = SwitchIncludeChannel;
    UF_SCHEME =/* getvaluedouble("Surface Flow Scheme") == 1? UF_SCHEME_CENTRALSIMPLE :*/ UF_SCHEME_BOUNDARYMUSCLE;
    UF_DTMIN = 0;
    UF_SOLIDPHASE = SwitchSolidPhase;
    UF_CHANNELFLOOD = SwitchChannelFlood;

    UF2D_Test = NewMap(0.0);

    //just for display
    UF2D_h = NewMap(0.0);
    UF2D_fsConc = NewMap(0.0);
    UF2D_sConc = NewMap(0.0);
    UF2D_tConc = NewMap(0.0);
    UF2D_velocity = NewMap(0.0);
    UF2D_u = NewMap(0.0);
    UF2D_v = NewMap(0.0);
    UF2D_Nr = NewMap(0.0);
    UF2D_q = NewMap(0.0);
    UF2D_qs = NewMap(0.0);

    UF1D_h = NewMap(0.0);
    UF1D_fsConc = NewMap(0.0);
    UF1D_sConc = NewMap(0.0);
    UF1D_tConc = NewMap(0.0);
    UF1D_velocity = NewMap(0.0);
    UF1D_Nr = NewMap(0.0);
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
    if(SwitchBarriers)
    {
        Barriers = ReadMap(DEM,getvaluename("barriers"));
        FOR_ROW_COL_MV
        {
            if(!pcr::isMV(Barriers->Drc))
            {
                UF2D_DEM->Drc += Barriers->Drc;
                DEM->Drc += Barriers->Drc;
            }
        }
    }
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
    UF2D_fax1 = NewMap(0.0);
    UF2D_fay1 = NewMap(0.0);
    UF2D_fax2 = NewMap(0.0);
    UF2D_fay2 = NewMap(0.0);
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
    UF2D_sax1 = NewMap(0.0);
    UF2D_say1 = NewMap(0.0);
    UF2D_sax2 = NewMap(0.0);
    UF2D_say2 = NewMap(0.0);
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
    UF1D_LDDh = NewMap(5.0);
    UF1D_Slope = NewMap(0.0);
    UF1D_LDD->setAllMV();
    if(SwitchIncludeChannel)
    {
        FOR_ROW_COL_MV_CH
        {
            UF1D_LDD->Drc = LDDChannel->Drc;
            UF1D_LDDw->Drc = ChannelWidth->Drc;
            UF1D_Slope->Drc = -ChannelGrad->Drc;
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
    UF1D_fa1 = NewMap(0.0);
    UF1D_fa2 = NewMap(0.0);
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
    UF1D_sa1 = NewMap(0.0);
    UF1D_sa2 = NewMap(0.0);
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
    UF_t6 = NewMap(0.0);
    UF_t7 = NewMap(0.0);
    UF_t8 = NewMap(0.0);
    UF_t9 = NewMap(0.0);
    UF_t10 = NewMap(0.0);
    UF_t11 = NewMap(0.0);

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

    //maps for forced and initial conditions
    UF2D_ForcedFVolume = NewMap(0.0);
    UF2D_ForcedSVolume = NewMap(0.0);
    UF2D_ForcedSDensity = NewMap(0.0);
    UF2D_ForcedSIFA = NewMap(0.0);
    UF2D_ForcedSRocksize = NewMap(0.0);

    UF2D_InitialFVolume = NewMap(0.0);
    UF2D_InitialSVolume = NewMap(0.0);
    UF2D_InitialSDensity = NewMap(0.0);
    UF2D_InitialSIFA = NewMap(0.0);
    UF2D_InitialSRocksize = NewMap(0.0);

    if(SwitchUFForced)
    {
        UF2D_ForcedFVolume = ReadMap(DEM,getvaluename("forcedfvolume"));
        UF2D_ForcedSVolume = ReadMap(DEM,getvaluename("forcedsvolume"));
        UF2D_ForcedSDensity = ReadMap(DEM,getvaluename("forcedsdensity"));
        UF2D_ForcedSIFA = ReadMap(DEM,getvaluename("forcedsifa"));
        UF2D_ForcedSRocksize = ReadMap(DEM,getvaluename("forcedsrocksize"));
    }
    UF_NeedsInitial = true;
    if(SwitchUFInitial)
    {
        qDebug() << "initial f + s";
        UF2D_Initialized = NewMap(0.0);
        UF2D_InitialTime = ReadMap(DEM,getvaluename("initiationtime"));
        UF2D_InitialFVolume = ReadMap(DEM,getvaluename("initialfvolume"));
        UF2D_InitialSVolume = ReadMap(DEM,getvaluename("initialsvolume"));
        UF2D_InitialSDensity = ReadMap(DEM,getvaluename("initialsdensity"));
        UF2D_InitialSIFA = ReadMap(DEM,getvaluename("initialsifa"));
        UF2D_InitialSRocksize = ReadMap(DEM,getvaluename("initialsrocksize"));
    }

    UF_InitializedF = 0;
    UF_InitializedS = 0;

    SourceSolid = NewMap(0.0);
    SourceSolidDensity = NewMap(0.0);
    SourceSolidRocksize = NewMap(0.0);
    SourceSolidIFA = NewMap(0.0);
    SourceFluid = NewMap(0.0);

    ChannelSourceSolid = NewMap(0.0);
    ChannelSourceSolidDensity = NewMap(0.0);
    ChannelSourceSolidRocksize = NewMap(0.0);
    ChannelSourceSolidIFA = NewMap(0.0);
    ChannelSourceFluid = NewMap(0.0);

    //multithreading library
    //ltr = new LisemThread();
    //ThreadPool = new LisemThreadPool();
    //ThreadPool->InitThreads(this);
    //ThreadPool->SetMaskInitial(UF2D_DEM);
}

void TWorld::UF_Close()
{

    //ThreadPool->clear();
    //delete ThreadPool;

}
