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

    UF_ENTRAINMENTCCONSTANT = getvaluedouble("Minimum Entrainment Height");
    UF_ENTRAINMENTCONSTANT = 1.0;
    UF_ENTRAINMENTTHRESHOLDCONSTANT = 1.0;
    UF_SUSPENDEDVISCOSITY=(getvaluedouble("Suspended Viscosity")== 1.0)? true:false;
    UF_LAXMULTIPLIER=0.1 *getvaluedouble("Lax Multiplier");
    UF_FRICTIONCORRECTION=getvaluedouble("Friction force correction");

    UF_SOLIDFLUIDDRAG =getvaluedouble("Solid-Fluid Drag Coefficient");
    UF_ENTRAINMENTCONSTANT =  getvaluedouble("Entrainment Coefficient")/1000000.0;
    UF_DEPOSITIONCONSTANT = 1;
    UF_DEPOSITIONTHRESHOLDCONSTANT = 0.6;
    UF_ENTRAINMENTROOTDEPTH = 0.25;

    UF_MAXSOLIDCONCENTRATION = 0.9;
    UF_MINIMUMENTRAINMENTHEIGHT = 0.00000001;

    UF_MANNINGCOEFFICIENT_FLUID = 1.0;
    UF_MANNINGCOEFFICIENT_SOLID = 0.5;
    UF_FrictionIterations = 1;
    UF_KINEMATIC_TIMESTEP_POWER= getvaluedouble("Kinematic Timestep Power");

    UF_USE_HLL2 = getvaluedouble("Use HLL2");

    UF_DENSITY_SUSPENDED = 2000;

    UF_FLOWCOMPACTION_MAXCONCENTRATION = 0.5;
    UF_FLOWCOMPACTION_CRITICALVELOCITY = 0.25;

    UF_FLOWCOMPACTION_DEPOSITIONPOROSITY = 0.5;

    UF_AVERAGEFACTOR = 0.5;
    UF2D_COURANTSCHEMEFACTOR = 0.25;
    UF1D_COURANTSCHEMEFACTOR = 0.5;

    UF_AddedSplash = NewMap(0.0);

    UF_DEMFEEDBACK = true;

    UF_DTAverage = 0;

    UF_MAX_NUM_VEL = 100;

    UF_NRA = 150000;

    UF_Alpha_DV  = getvaluedouble("Viscosity Alpha");
    UF_Beta_DV = getvaluedouble("Viscosity Beta");

    UF_Alpha_YS = 0.1;
    UF_Beta_YS = 20.0;

    UF_Alpha_DR = 0.0538;
    UF_Beta_DR = 6.0896;

    UF_Intersect_K = 30;
    UF_Slope_K = 100000.0;

    UF_t1 = NewMap(0.0);
    UF_t2 = NewMap(0.0);
    UF_t3 = NewMap(0.0);
    UF_t4 = NewMap(0.0);

    //internal use
    UF_1DACTIVE = SwitchIncludeChannel;

    //old scheme not supported anymore
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

    UF2D_qout = NewMap(0.0);
    UF2D_qsout = NewMap(0.0);
    UF2D_qblout = NewMap(0.0);
    UF2D_qssout = NewMap(0.0);
    UF1D_qout = NewMap(0.0);
    UF1D_qsout = NewMap(0.0);
    UF1D_qblout = NewMap(0.0);
    UF1D_qssout = NewMap(0.0);

    UF2D_TimeStep = NewMap(0.0);
    UF2D_SPH = NewMap(0.0);
    UF2D_FPH = NewMap(0.0);

    //internal slope functions
    UF2D_Slope = NewMap(0.0);
    UF2D_SlopeX = NewMap(0.0);
    UF2D_SlopeY = NewMap(0.0);
    UF1D_LDDs = NewMap(0.0);

    //actual calculation variables
    ////2D
    UF2D_DEM = NewMap(0.0);
    copy(*UF2D_DEMOriginal,*DEM);
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
    UF2D_tsf = NewMap(0.0);

    //suspended matter
    UF2D_ssqx1 = NewMap(0.0);
    UF2D_ssqy1 = NewMap(0.0);
    UF2D_ssqx2 = NewMap(0.0);
    UF2D_ssqy2 = NewMap(0.0);
    UF2D_empty = NewMap(0.0);

    //solid phase
    if(UF_SOLIDPHASE)
    {

        UF2D_Nr = NewMap(0.0);

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

        UF2D_sqx = NewMap(0.0);
        UF2D_sqy = NewMap(0.0);

        UF2D_sqx1 = NewMap(0.0);
        UF2D_sqy1 = NewMap(0.0);
        UF2D_sqx2 = NewMap(0.0);
        UF2D_sqy2 = NewMap(0.0);

        UF2D_Compaction = NewMap(0.0);

        UF2D_STL = NewMap(0.0);
        UF2D_STLA = NewMap(0.0);
        UF2D_STLH = NewMap(0.0);
        UF2D_ST = NewMap(0.0);
        UF1D_ST = NewMap(0.0);

        UF2D_SolidFrictionFraction = NewMap(0.0);
        UF1D_SolidFrictionFraction = NewMap(0.0);


        UF2D_DC = NewMap(0.0);


        //solid phase
        UF2D_sn = NewMap(0.0);
        UF2D_sun = NewMap(0.0);
        UF2D_svn = NewMap(0.0);
    }

    //for new timestep
    //fluid phase
    UF2D_fn = NewMap(0.0);
    UF2D_fun = NewMap(0.0);
    UF2D_fvn = NewMap(0.0);


    ////1D
    UF1D_LDD = NewMap(0.0);
    UF1D_LDDw = NewMap(0.0);
    UF1D_LDDh = NewMap(0.0);
    UF1D_CellR = NewMap(0.0);
    UF1D_CellC = NewMap(0.0);
    UF1D_Slope = NewMap(0.0);
    UF1D_LDD->setAllMV();

    if(UF_1DACTIVE)
    {
        FOR_ROW_COL_MV_CH
        {
            UF1D_LDD->Drc = LDDChannel->Drc;
            UF1D_LDDw->Drc = ChannelWidth->Drc;
            UF1D_Slope->Drc = -ChannelGrad->Drc;
            UF1D_LDDs->Drc = std::fabs(ChannelGrad->Drc);

            if(UF_CHANNELFLOOD)
            {
                UF1D_LDDh->Drc = ChannelDepth->Drc;
            }
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
    if(SwitchSolidPhase)
    {
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
    }

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

    UF2D_EntrainmentSF = NewMap(0.0);

    UF2D_Infiltration = NewMap(0.0);
    UF1D_Infiltration = NewMap(0.0);

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



    if(SwitchUFForced)
    {

        //maps for forced and initial conditions
        UF2D_ForcedFVolume = NewMap(0.0);
        if(SwitchSolidPhase)
        {
            UF2D_ForcedSVolume = NewMap(0.0);
            UF2D_ForcedSDensity = NewMap(0.0);
            UF2D_ForcedSIFA = NewMap(0.0);
            UF2D_ForcedSRocksize = NewMap(0.0);
        }

        UF2D_ForcedFVolume = ReadMap(DEM,getvaluename("forcedfvolume"));
        if(SwitchSolidPhase)
        {
            UF2D_ForcedSVolume = ReadMap(DEM,getvaluename("forcedsvolume"));
            UF2D_ForcedSDensity = ReadMap(DEM,getvaluename("forcedsdensity"));
            UF2D_ForcedSIFA = ReadMap(DEM,getvaluename("forcedsifa"));
            UF2D_ForcedSRocksize = ReadMap(DEM,getvaluename("forcedsrocksize"));
        }

        if(SwitchSolidPhase)
        {
            FOR_ROW_COL_MV
            {
                UF2D_ForcedSIFA->Drc = UF2D_ForcedSIFA->Drc * Calibrate_YS;

            }
        }
    }
    UF_NeedsInitial = true;
    if(SwitchUFInitial)
    {

        UF2D_InitialFVolume = NewMap(0.0);
        if(SwitchSolidPhase)
        {
            UF2D_InitialSVolume = NewMap(0.0);
            UF2D_InitialSDensity = NewMap(0.0);
            UF2D_InitialSIFA = NewMap(0.0);
            UF2D_InitialSRocksize = NewMap(0.0);
        }

        UF2D_Initialized = NewMap(0.0);
        UF2D_InitialTime = ReadMap(DEM,getvaluename("initiationtime"));
        UF2D_InitialFVolume = ReadMap(DEM,getvaluename("initialfvolume"));

        if(SwitchSolidPhase)
        {
            UF2D_InitialSVolume = ReadMap(DEM,getvaluename("initialsvolume"));
            UF2D_InitialSDensity = ReadMap(DEM,getvaluename("initialsdensity"));
            UF2D_InitialSIFA = ReadMap(DEM,getvaluename("initialsifa"));
            UF2D_InitialSRocksize = ReadMap(DEM,getvaluename("initialsrocksize"));
        }

        if(SwitchSolidPhase)
        {
            FOR_ROW_COL_MV
            {
                UF2D_InitialSIFA->Drc = UF2D_InitialSIFA->Drc * Calibrate_YS;

            }
        }
    }

    UF_InitializedF = 0;
    UF_InitializedS = 0;

    SourceFluid = NewMap(0.0);

    if(SwitchSolidPhase)
    {
        SourceSolid = NewMap(0.0);
        SourceSolidDensity = NewMap(0.0);
        SourceSolidRocksize = NewMap(0.0);
        SourceSolidIFA = NewMap(0.0);


        ChannelSourceSolid = NewMap(0.0);
        ChannelSourceSolidDensity = NewMap(0.0);
        ChannelSourceSolidRocksize = NewMap(0.0);
        ChannelSourceSolidIFA = NewMap(0.0);


    }

    ChannelSourceFluid = NewMap(0.0);

    //multithreading library
    ThreadPool = new LisemThreadPool();

    ThreadPool->InitThreads(this);

    ThreadPool->SetMaskInitial(UF2D_DEM,UF1D_LDD);

    UF1D_OutletDistance = NewMap(0.0);

    if(UF_1DACTIVE)
    {

        FOR_ROW_COL_MV_CH
        {
            pcr::setMV(tmb->Drc);
        }


        for (int  ro = 0; ro < _nrRows; ro++){
        for (int  co = 0; co < _nrCols; co++){
        if(!pcr::isMV(LDDChannel->data[ro][co]))
        {
            if(LDDChannel->data[ro][co] == 5)
            {

                int ncells = 0;

                int dx[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
                int dy[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

                /// Linked list of cells in order of LDD flow network, ordered from pit upwards
                LDD_LINKEDLIST *list = NULL, *temp = NULL;
                list = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));

                list->prev = NULL;
                /// start gridcell: outflow point of area
                list->rowNr = ro;
                list->colNr = co;

                while (list != NULL)
                {
                    int i = 0;
                    bool  subCachDone = true; // are sub-catchment cells done ?
                    int rowNr = list->rowNr;
                    int colNr = list->colNr;

                    /** put all points that have to be calculated to calculate the current point in the list,
                     before the current point */
                    for (i=1; i<=9; i++)
                    {
                        int r, c;
                        int ldd = 0;

                        // this is the current cell
                        if (i==5)
                            continue;

                        r = rowNr+dy[i];
                        c = colNr+dx[i];

                        if (INSIDE(r, c) && !pcr::isMV(LDDChannel->Drc))
                            ldd = (int) LDDChannel->Drc;
                        else
                            continue;

                        // check if there are more cells upstream, if not subCatchDone remains true
                        if (pcr::isMV(tmb->Drc) &&
                                FLOWS_TO(ldd, r, c, rowNr, colNr) &&
                                INSIDE(r, c))
                        {
                            temp = (LDD_LINKEDLIST *)malloc(sizeof(LDD_LINKEDLIST));
                            temp->prev = list;
                            UF1D_OutletDistance->Drc = ncells++;
                            list->rowNr = r;
                            list->colNr = c;
                            subCachDone = false;
                        }
                    }

                    // all cells above a cell are linked in a "sub-catchment or branch
                    // continue with water and sed calculations
                    // rowNr and colNr are the last upstreM cell linked
                    if (subCachDone)
                    {
                        int r = rowNr;
                        int c = colNr;
                        tmb->Drc = 0;

                        temp=list;
                        list=list->prev;
                        if(list != NULL)
                        {
                            ncells =UF1D_OutletDistance->data[list->rowNr][list->colNr];
                        }

                        free(temp);
                        // go to the previous cell in the list

                    }/* eof subcatchment done */
                } /* eowhile list != NULL */
            }
        }}}
    }

}

void TWorld::UF_ExtendChannel()
{

    /*ChannelDepthExtended = NewMap(0.0);
    ChannelWidthExtended = NewMap(0.0);
    ChannelMaskExtended = NewMap(0.0);
    ChannelFlowWidth = NewMap(0.0);


    if(!SwitchIncludeChannel)
    {
        return;
    }

    copy(*ChannelWidthExtended, *ChannelWidthMax);
    copy(*ChannelDepthExtended, *ChannelDepth);

    ChannelNeighborsExtended = NewMap(0.0);
    ChannelSourceXExtended = NewMap(0.0);
    ChannelSourceYExtended = NewMap(0.0);
    ChannelBoundaryExtended = NewMap(0.0);
    ChannelBoundaryLExtended = NewMap(0.0);
    ChannelBoundaryRExtended = NewMap(0.0);


    //Channel cells are part of the extended channel
    FOR_ROW_COL_MV_CH
    {
        ChannelMaskExtended->Drc = (ChannelWidthMax->Drc > 0 ? 1.0 : 0.0);
    }

    int dxl[10] = {0, -1, 0, 1, -1, 0, 1, -1, 0, 1};
    int dyl[10] = {0, 1, 1, 1, 0, 0, 0, -1, -1, -1};

    //set width etc for channel cells, for these we know they are within the extende channel
    double MaxWidth = 0;
    FOR_ROW_COL_MV_CH
    {
        if(ChannelWidthMax->Drc > MaxWidth)
        {
            MaxWidth = ChannelWidthMax->Drc;
        }
        ChannelMaskExtended->Drc = 1;
        ChannelDepthExtended->Drc = ChannelDepth->Drc;
        ChannelWidthExtended->Drc = std::min(_dx,ChannelWidthMax->Drc);
        ChannelNeighborsExtended->Drc = -1;
    }

    //for iteration
    int maxdistance = ((int) ((std::max(0.0,MaxWidth-_dx)/2.0)/_dx)) + 1;

    //for checking wheter channel lies within cell
    //double maxd = (std::max(0.0,MaxWidth-_dx)/2.0);

    FOR_ROW_COL_MV
    {

        if(!pcr::isMV(LDDChannel->Drc))
        {
            ChannelSourceXExtended->Drc = c;
            ChannelSourceYExtended->Drc = r;
            continue;
        }
        bool found = false;
        double found_distance = 9999999.9;
        int i = 1;
        while(i < maxdistance + 1 && !(found_distance < double(i) * _dx))
        {

            for(int r2 = r - i; r2 < r + i + 1; r2++)
            {
                int c2 = c - i;
                if(!OUTORMV(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            //double dsx = std::max(0.0,dx-0.5*(dx/rd));
                            //double dsy = std::max(0.0,dy-0.5*(dy/rd));
                            //double d = sqrt(dsx * dsx + dsy*dsy);

                            double width =  std::min(_dx,std::max(0.0,0.5 * ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                                ChannelBoundaryExtended->Drc = 1.0;
                                if(-dyl[(int)LDDChannel->data[r2][c2]] * double(c2-c) + dxl[(int)LDDChannel->data[r2][c2]] * double(r2-r) < 0)
                                {
                                     ChannelBoundaryLExtended->Drc = 1;
                                     ChannelBoundaryRExtended->Drc = 0;
                                }else
                                {
                                     ChannelBoundaryRExtended->Drc = 1;
                                     ChannelBoundaryLExtended->Drc = 0;
                                }
                            }
                        }
                    }
                }
            }
            for(int r2 = r - i; r2 < r + i + 1; r2++)
            {
                int c2 = c+i;

                if(!OUTORMV(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            //double dsx = std::max(0.0,dx-0.5*(dx/rd));
                            //double dsy = std::max(0.0,dy-0.5*(dy/rd));
                            //double d = sqrt(dsx * dsx + dsy*dsy);
                            double width =  std::min(_dx,std::max(0.0,0.5 * ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                                ChannelBoundaryExtended->Drc = 1.0;
                                if(-dyl[(int)LDDChannel->data[r2][c2]] * double(c2-c) + dxl[(int)LDDChannel->data[r2][c2]] * double(r2-r) < 0)
                                {
                                     ChannelBoundaryLExtended->Drc = 1;
                                     ChannelBoundaryRExtended->Drc = 0;
                                }else
                                {
                                     ChannelBoundaryRExtended->Drc = 1;
                                     ChannelBoundaryLExtended->Drc = 0;
                                }
                            }
                        }
                    }
                }
            }
            for(int c2 = c - i + 1; c2 < c + i; c2++)
            {
                int r2 = r - i;
                if(!OUTORMV(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            //double dsx = std::max(0.0,dx-0.5*(dx/rd));
                            //`double dsy = std::max(0.0,dy-0.5*(dy/rd));
                            //double d = sqrt(dsx * dsx + dsy*dsy);

                            double width =  std::min(_dx,std::max(0.0,0.5 * ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                                ChannelBoundaryExtended->Drc = 1.0;
                                if(-dyl[(int)LDDChannel->data[r2][c2]] * double(c2-c) + dxl[(int)LDDChannel->data[r2][c2]] * double(r2-r) < 0)
                                {
                                     ChannelBoundaryLExtended->Drc = 1;
                                     ChannelBoundaryRExtended->Drc = 0;
                                }else
                                {
                                     ChannelBoundaryRExtended->Drc = 1;
                                     ChannelBoundaryLExtended->Drc = 0;
                                }
                            }
                        }
                    }
                }
            }
            for(int c2 = c - i + 1; c2 < c + i; c2++)
            {
                int r2 = r+i;
                if(!OUTORMV(r2,c2))
                {
                    if(!pcr::isMV(LDDChannel->data[r2][c2]))
                    {
                        double dx = fabs(double(c2-c));
                        double dy = fabs(double(r2-r));
                        double rd = sqrt(dx*dx+dy*dy)*_dx;

                        if( rd > 0)
                        {
                            //double dsx = std::max(0.0,dx-0.5*(dx/rd));
                            //double dsy = std::max(0.0,dy-0.5*(dy/rd));
                            //double d = sqrt(dsx * dsx + dsy*dsy);

                            double width =  std::min(_dx,std::max(0.0,0.5 * ChannelWidthMax->data[r2][c2] - rd));
                            if(width > 0)
                            {
                                found_distance = rd;
                                found = true;
                                ChannelDepthExtended->Drc += std::min(_dx,std::max(0.0,ChannelWidthMax->data[r2][c2] - rd)) * ChannelDepth->data[r2][c2];
                                ChannelWidthExtended->Drc += width;
                                ChannelNeighborsExtended->Drc += 1;
                                ChannelSourceXExtended->Drc = c2;
                                ChannelSourceYExtended->Drc = r2;
                                ChannelMaskExtended->Drc =1;
                                ChannelBoundaryExtended->Drc = 1.0;
                                if(-dyl[(int)LDDChannel->data[r2][c2]] * double(c2-c) + dxl[(int)LDDChannel->data[r2][c2]] * double(r2-r) < 0)
                                {
                                     ChannelBoundaryLExtended->Drc = 1;
                                     ChannelBoundaryRExtended->Drc = 0;
                                }else
                                {
                                     ChannelBoundaryRExtended->Drc = 1;
                                     ChannelBoundaryLExtended->Drc = 0;
                                }
                            }
                        }
                    }
                }
            }

            i++;
        }

        if(found)
        {
            if(ChannelWidthExtended->Drc > 0 && ChannelNeighborsExtended->Drc > 0) //VJ!!!! ChannelWidthExtended without Drc
            {
                ChannelDepthExtended->Drc /= ChannelWidthExtended->Drc;
                ChannelWidthExtended->Drc /= ChannelNeighborsExtended->Drc;
            }
        }
    }

    FOR_ROW_COL_MV
    {
        if(ChannelMaskExtended->Drc == 1)
        {
            if(!pcr::isMV(LDDChannel->Drc) && !(ChannelWidthMax->Drc > _dx))
            {

                    ChannelBoundaryExtended->Drc = 1;
                    ChannelBoundaryRExtended->Drc = 1;
                    ChannelBoundaryLExtended->Drc = 1;
            }

            double nn = 0;

            if(IsExtendedChannel(r,c,1,0)) { nn ++;};
            if(IsExtendedChannel(r,c,-1,0)) { nn ++;};
            if(IsExtendedChannel(r,c,0,1)) { nn ++;};
            if(IsExtendedChannel(r,c,0,-1)) { nn ++;};
            if(nn == 4){
                ChannelBoundaryExtended->Drc = 0;
                ChannelBoundaryRExtended->Drc = 0;
                ChannelBoundaryLExtended->Drc = 0;
            }
//            if(!ChannelBoundaryExtended->Drc == 1 && nn < 4)
                if(ChannelBoundaryExtended->Drc != 1 && nn < 4)
            {
                ChannelBoundaryExtended->Drc = 1;
                ChannelBoundaryRExtended->Drc = 1;
                ChannelBoundaryLExtended->Drc = 1;
            }
            if(IsExtendedChannel(r,c,-1,1)) { nn ++;};
            if(IsExtendedChannel(r,c,-1,-1)) { nn ++;};
            if(IsExtendedChannel(r,c,1,1)) { nn ++;};
            if(IsExtendedChannel(r,c,1,-1)) { nn ++;};
        }
    }


    return;*/
}

bool TWorld::IsExtendedChannel(int r, int c, int dr, int dc)
{
    /*if(!OUTORMV(r+dr,c+dc))
    {
        return ChannelMaskExtended->data[r+dr][c+dc] == 1;
    }else
    {
        return true;
    }*/

}


//---------------------------------------------------------------------------
// Distributes a certain value over the actual channel width (used for display stuff)

void TWorld::DistributeOverExtendedChannel(cTMap * _In, cTMap * _Out, bool do_not_divide,bool proportional)
{
    /*double totala=0;
    double totalb=0;
    FOR_ROW_COL_MV
    {
        if(!pcr::isMV(LDDChannel->Drc))
        {
            totala += _In->Drc;
        }
        if(ChannelMaskExtended->Drc == 1)
        {
            double ow = ChannelWidthMax->data[(int)ChannelSourceYExtended->Drc][(int)ChannelSourceXExtended->Drc];
            if(ow> 0 )
            {
                double div = do_not_divide? (proportional? ChannelWidthExtended->Drc/ _dx : 1.0): (ChannelWidthExtended->Drc / ow);
                _Out->Drc = _In->data[(int)ChannelSourceYExtended->Drc][(int)ChannelSourceXExtended->Drc] * div;
            }
            totalb += _Out->Drc;
        }else
        {
            _Out->Drc  = 0.0;
        }
    }
    if(totalb > 0)
    {
        FOR_ROW_COL_MV
        {
            if(ChannelMaskExtended->Drc == 1)
            {
                _Out->Drc *= totala/totalb;
            }
        }
    }*/
}
