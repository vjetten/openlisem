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
\file UFModel.h
\brief list of class members and functions related to the unified flow solutions. Is included in the model.h class!
*/

    //////////////
    //Unified FLow
    //////////////
    //constants
    double UF_Courant;
    double UF_Aspect;
    double UF_Chi;
    double UF_Ksi;
    double UF_j;
    double UF_Gravity;
    double UF_GravitySqrt;
    double UF2D_MinimumDT;
    double UF1D_MinimumDT;
    double UF_SigmaDiffusion;
    double UF_MANNINGCOEFFICIENT_FLUID;
    double UF_MANNINGCOEFFICIENT_SOLID;

    double UF_SOLIDFLUIDDRAG;

    bool UF_SUSPENDEDVISCOSITY;
    double UF_LAXMULTIPLIER;
    double UF_FRICTIONCORRECTION;
    double UF_ENTRAINMENTCCONSTANT;
    double UF_ENTRAINMENTCONSTANT;
    double UF_ENTRAINMENTTHRESHOLDCONSTANT;
    double UF_DEPOSITIONCONSTANT;
    double UF_DEPOSITIONTHRESHOLDCONSTANT;
    double UF_ENTRAINMENTROOTDEPTH;

    double UF_MAXSOLIDCONCENTRATION;
    double UF_MINIMUMENTRAINMENTHEIGHT;

    double UF_AVERAGEFACTOR;

    cTMap* UF_AddedSplash;
    bool UF_DEMFEEDBACK;


    double UF_DENSITY_SUSPENDED;
    double UF_FLOWCOMPACTION_MAXCONCENTRATION;
    double UF_FLOWCOMPACTION_CRITICALVELOCITY;

    double UF_FLOWCOMPACTION_DEPOSITIONPOROSITY;

    double UF2D_COURANTSCHEMEFACTOR;
    double UF1D_COURANTSCHEMEFACTOR;

    double UF_KINEMATIC_TIMESTEP_POWER;
    double UF_USE_HLL2;
    double UF_DTAverage;

    int UF_FrictionIterations;

    bool SpatiallyDynamicTimestep = true;

    double UF_DISPLAYFLOODMINIMUM = 0;
    double UF_DISPLAYDEBRISFLOWMINIMUM = 0;

    double UF_MAX_NUM_VEL;

    double UF_NRA;

    double UF_Alpha_DV;
    double UF_Beta_DV;

    double UF_Alpha_YS;
    double UF_Beta_YS;

    double UF_Alpha_DR;
    double UF_Beta_DR;

    double UF_Intersect_K;
    double UF_Slope_K;

    //internal use
    double UF_1DACTIVE;
    double UF_DTMIN;
    int UF_SCHEME;
    bool UF_SOLIDPHASE;
    bool UF_CHANNELFLOOD;
    static const int UF_SCHEME_CENTRALSIMPLE =1;
    static const int UF_SCHEME_BOUNDARYMUSCLE =2;

    cTMap * UF2D_Test;

    cTMap * UF2D_dhdx;
    cTMap * UF2D_dhdy;

    cTMap * UF_t1;
    cTMap * UF_t2;
    cTMap * UF_t3;
    cTMap * UF_t4;

    cTMap * UF1D_fq1lim;
    cTMap * UF1D_fq2lim;
    cTMap * UF1D_sq1lim;
    cTMap * UF1D_sq2lim;

    //just for display
    cTMap * UF2D_h;
    cTMap * UF2D_q;
    cTMap * UF2D_qout;
    cTMap * UF2D_qs;
    cTMap * UF2D_qsout;
    cTMap * UF2D_qblout;
    cTMap * UF2D_qssout;
    cTMap * UF2D_fsConc;
    cTMap * UF2D_sConc;
    cTMap * UF2D_tConc;
    cTMap * UF2D_velocity;
    cTMap * UF2D_Nr;
    cTMap * UF2D_u;
    cTMap * UF2D_v;
    cTMap * UF2D_TimeStep;
    cTMap * UF2D_SPH;
    cTMap * UF2D_FPH;
    cTMap * UF2D_tsf;
    cTMap * UF2D_SolidFrictionFraction;
    cTMap * UF1D_SolidFrictionFraction;

    cTMap * UF1D_h;
    cTMap * UF1D_q;
    cTMap * UF1D_qout;
    cTMap * UF1D_qs;
    cTMap * UF1D_qsout;
    cTMap * UF1D_qblout;
    cTMap * UF1D_qssout;
    cTMap * UF1D_fsConc;
    cTMap * UF1D_sConc;
    cTMap * UF1D_tConc;
    cTMap * UF1D_Nr;
    cTMap * UF1D_velocity;


    //internal slope functions
    cTMap * UF2D_Slope;
    cTMap * UF2D_SlopeX;
    cTMap * UF2D_SlopeY;
    cTMap * UF1D_LDDs;

    //actual calculation variables
    ////2D
    cTMap * UF2D_DEMOriginal;
    cTMap * UF2D_DEM;
    cTMap * UF2D_T;
    cTMap * UF2D_DT;
    cTMap * UF2D_DTStep;
    cTMap * UF2D_CellR;
    cTMap * UF2D_CellC;
    cTMap * UF2D_Courant;

    //fluid phase
    cTMap * UF2D_f;
    cTMap * UF2D_visc;
    cTMap * UF2D_fu;
    cTMap * UF2D_fv;
    cTMap * UF2D_fax;
    cTMap * UF2D_fay;
    cTMap * UF2D_fax1;
    cTMap * UF2D_fay1;
    cTMap * UF2D_fax2;
    cTMap * UF2D_fay2;
    cTMap * UF2D_fqx1;
    cTMap * UF2D_fqy1;
    cTMap * UF2D_fqx2;
    cTMap * UF2D_fqy2;
    cTMap * UF2D_ssm;
    cTMap * UF2D_blm;
    cTMap * UF2D_sstc;
    cTMap * UF2D_bltc;
    cTMap * UF2D_fsc;
    cTMap * UF2D_fsd;
    //solid phase
    cTMap * UF2D_s;
    cTMap * UF2D_d;
    cTMap * UF2D_ifa;
    cTMap * UF2D_rocksize;
    cTMap * UF2D_su;
    cTMap * UF2D_sv;
    cTMap * UF2D_sax;
    cTMap * UF2D_say;
    cTMap * UF2D_sax1;
    cTMap * UF2D_say1;
    cTMap * UF2D_sax2;
    cTMap * UF2D_say2;
    cTMap * UF2D_sqx;
    cTMap * UF2D_sqy;
    cTMap * UF2D_sqx1;
    cTMap * UF2D_sqy1;
    cTMap * UF2D_sqx2;
    cTMap * UF2D_sqy2;

    cTMap * UF2D_Compaction;

    cTMap * UF2D_DC;

    cTMap * UF2D_STL;
    cTMap * UF2D_STLA;
    cTMap * UF2D_STLH;
    cTMap * UF2D_ST;
    cTMap * UF1D_ST;

    //for new timestep
    //fluid phase
    cTMap * UF2D_fn;
    cTMap * UF2D_fun;
    cTMap * UF2D_fvn;
    //solid phase
    cTMap * UF2D_sn;
    cTMap * UF2D_sun;
    cTMap * UF2D_svn;


    ////1D
    cTMap * UF1D_OutletDistance;

    cTMap * UF1D_LDD;
    cTMap * UF1D_LDDw;
    cTMap * UF1D_LDDh;
    cTMap * UF1D_T;
    cTMap * UF1D_DT;
    cTMap * UF1D_CellR;
    cTMap * UF1D_CellC;
    cTMap * UF1D_Slope;
    cTMap * UF1D_DTStep;
    cTMap * UF1D_Courant;
    //fluid phase
    cTMap * UF1D_f;
    cTMap * UF1D_fstore;
    cTMap * UF1D_visc;
    cTMap * UF1D_fu;
    cTMap * UF1D_fa;
    cTMap * UF1D_fa1;
    cTMap * UF1D_fa2;
    cTMap * UF1D_fq1;
    cTMap * UF1D_fq2;
    cTMap * UF1D_ssm;
    cTMap * UF1D_blm;
    cTMap * UF1D_sstc;
    cTMap * UF1D_bltc;
    cTMap * UF1D_fsc;
    cTMap * UF1D_fsd;

    //solid phase
    cTMap * UF1D_sstore;
    cTMap * UF1D_s;
    cTMap * UF1D_d;
    cTMap * UF1D_ifa;
    cTMap * UF1D_rocksize;
    cTMap * UF1D_su;
    cTMap * UF1D_sa;
    cTMap * UF1D_sa1;
    cTMap * UF1D_sa2;
    cTMap * UF1D_sq1;
    cTMap * UF1D_sq2;

    //for new timestep
    //fluid phase
    cTMap * UF1D_fn;
    cTMap * UF1D_fun;
    //solid phase
    cTMap * UF1D_sn;
    cTMap * UF1D_sun;

    //Multiclass sediment functions
    QList<cTMap *> UF2D_ssm_D;
    QList<cTMap *> UF1D_ssm_D;
    QList<cTMap *> UF2D_blm_D;
    QList<cTMap *> UF1D_blm_D;

    QList<cTMap *> UF2D_sstc_D;
    QList<cTMap *> UF1D_sstc_D;
    QList<cTMap *> UF2D_bltc_D;
    QList<cTMap *> UF1D_bltc_D;

    cTMap * UF1D_Dep;
    cTMap * UF1D_Det;
    cTMap * UF2D_Dep;
    cTMap * UF2D_Det;

    cTMap * UF2D_Infiltration;
    cTMap * UF1D_Infiltration;

    double UF2D_foutflow = 0;
    double UF2D_fsoutflow = 0;
    double UF2D_soutflow = 0;
    double UF1D_foutflow = 0;
    double UF1D_fsoutflow = 0;
    double UF1D_soutflow = 0;


    //General Function
    void UnifiedFlow();
    void UF_Init();
    void UF_ExtendChannel();
    bool IsExtendedChannel(int r, int c, int dr, int dc);
    void DistributeOverExtendedChannel(cTMap * _In, cTMap * _Out, bool do_not_divide,bool proportional);

    void UF_SetInput();
    void UF_SetOutput();

    void UF_Compute(int thread);
    void UF_SWAP2D_MT(int thread,cTMap *x, cTMap *y);
    void UF_SWAP1D_MT(int thread,cTMap *x, cTMap *y);

    //set display maps
    void UF_UpdateDisplayMaps(int thread,cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                                     cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                                     cTMap * _fu1D,cTMap * _s1D,
                                     cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                                     cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                                     cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                                     cTMap * _su2D,cTMap * _sv2D);

    bool UF_Input_first = true;

    ////2D version

    //main calculation scheme
    double UF2D_Scheme(int thread,cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv);

    //source function
    double UF2D_Source(int thread, cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv);
    void UF2D_FluidSource(int thread, cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_f);
    void UF2D_SolidSource(int thread, cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_s);

    //advection functions for simple scheme
    void UF2D_Diffuse_mass(int thread, cTMap* dt, cTMap * _dem,cTMap * _m,cTMap * _f,cTMap * _fu, cTMap * _fv, cTMap * out_m);

    //advection with MUSCLE scheme
    void UF2D_Advect2_Momentum(int thread, cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_fu, cTMap * out_fv,  cTMap * out_su, cTMap * out_sv,cTMap * out_qfx1,cTMap *out_qfx2,cTMap * out_qfy1,cTMap *out_qfy2,cTMap * out_qsx1,cTMap *out_qsx2,cTMap * out_qsy1,cTMap *out_qsy2);
    double UF2D_Advect2_mass(int thread, cTMap* dt, cTMap * _dem,cTMap * _m, cTMap * f,cTMap * _qx1, cTMap * _qy1,cTMap * _qx2, cTMap * _qy2, cTMap * out_m,cTMap *outflow = 0);
    void UF2D_Advect2_prop(int thread, cTMap* dt, cTMap * _dem,cTMap * _m, cTMap * f,cTMap * _qx1, cTMap * _qy1,cTMap * _qx2, cTMap * _qy2,cTMap *_prop, cTMap * out_prop = 0);

    //momentum functions
    void UF2D_FluidMomentum2Source(int thread, cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv);
    void UF2D_SolidMomentum2Source(int thread, cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv);
    void UF2D_FluidApplyMomentum2(int thread, cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv);
    void UF2D_SolidApplyMomentum2(int thread, cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv);

    //not used
    void UF2D_Stored_mass(cTMap * dt, cTMap* dem, cTMap *_f,cTMap * _s, cTMap * out_f, cTMap * out_s);

    ////1D version

    //main calculation scheme
    void UF1D_Scheme(int thread, cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su);

    //source functions
    void UF1D_Source(int thread, cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su);
    void UF1D_FluidSource(int thread, cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_f);
    void UF1D_SolidSource(int thread, cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_s);

    //advection functions for simple scheme
    void UF1D_Diffuse_mass(int thread, cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap * _m, cTMap * _f,cTMap * _fu,cTMap * out_m);

    //advection functions for MUSCLE scheme
    void UF1D_Advect2_Momentum(int thread, cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_fu,  cTMap * out_su , cTMap * out_fq1, cTMap * out_fq2,cTMap * out_sq1, cTMap * out_sq2);
    double UF1D_Advect2_mass(int thread, cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _m, cTMap * _f, cTMap * _q1,cTMap * _q2, cTMap * out_m, cTMap *outflow = 0);
    void UF1D_Advect2_prop(int thread, cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _m, cTMap * _f, cTMap * _q1,cTMap * _q2,cTMap *_prop, cTMap * out_prop = 0);

    //momentum functions
    void UF1D_SolidMomentum2Source(int thread, cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu);
    void UF1D_FluidMomentum2Source(int thread, cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_su);
    void UF1D_SolidApplyMomentum2(int thread, cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_su);
    void UF1D_FluidApplyMomentum2(int thread, cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu);


    //not used
    void UF1D_Stored_mass(cTMap * dt, cTMap * _ldd,cTMap * _lddw, cTMap *_f,cTMap * _s, cTMap * out_f, cTMap * out_s);

    ////boundary conditions
    double UF_BoundaryFlux2D(double dt, double cellx, double celly, double f, double s, double fu, double fv, double su, double sv, double slopeX, double slopeY,double NN, int dr, int dc );
    double UF_BoundaryFlux1D(double dt, double width, double f, double s, double fu, double su, double slope, double NN, bool front );

    void UF2D1D_LaxNumericalCorrection(int thread, cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                                    cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                                    cTMap * _fu1D,cTMap * _s1D,
                                    cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                                    cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                                    cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                                    cTMap * _su2D,cTMap * _sv2D, int nh = 1, int nv = 1);

    ////Timestep functions
    double UF_InitTimeStep(cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                           cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                           cTMap * _fu1D,cTMap * _s1D,
                           cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                           cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                           cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                           cTMap * _su2D,cTMap * _sv2D,cTMap * out_t2d,cTMap * out_dt2d, cTMap * out_dtstep2d,cTMap * out_t1d,cTMap * out_dt1d, cTMap * out_dtstep1d);
    double UF_TimeStep(double t, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                       cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                       cTMap * _fu1D,cTMap * _s1D,
                       cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                       cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                       cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                       cTMap * _su2D,cTMap * _sv2D,cTMap * out_t2d,cTMap * out_dt2d, cTMap * out_dtstep2d,cTMap * out_t1d,cTMap * out_dt1d, cTMap * out_dtstep1d);

    ////connections
    void UF2D1D_Connection(int thread,cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                           cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                           cTMap * _fu1D,cTMap * _s1D,
                           cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                           cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                           cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                           cTMap * _su2D,cTMap * _sv2D);

    void UF2D1D_ChannelWater(int thread,cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                                     cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                                     cTMap * _fu1D,cTMap * _s1D,
                                     cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                                     cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                                     cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                                     cTMap * _su2D,cTMap * _sv2D);

    void UF2D1D_Infiltration(int thread,cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                             cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                             cTMap * _fu1D,cTMap * _s1D,
                             cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                             cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                             cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                             cTMap * _su2D,cTMap * _sv2D);

    ////Momentum balance functions
    double UF_Friction(double a,double dt,double velx,double vely, double NN, double h, double slope,bool solid, bool channel = false, double flowwidth = 0, double solids = 0,double ff = 1.0, double sf = 0.0, double Nr = 1.0);
    double UF_SFriction(double a, double v, double dt);
    double UF_Friction2(double vel,double dt,double velx,double vely, double NN, double h, double slope,bool solid, bool channel = false, double flowwidth = 0);

    double UF2D_MomentumBalanceFluid(bool x, double _f,double _s,double fu, double fv, double su, double sv, double ff, double sf, double Nr, double Nra, double ifa, double gamma, double visc, double pbf, double SlopeX, double SlopeY,
                                     double dhfdx,double dhfdy, double dh2pbdx, double dh2pbdy, double dsfdx, double dsfdy, double ddsfdxx, double ddsfdyy, double ddsfdxy,
                                     double dfudx, double dfudy, double dfvdx, double dfvdy, double ddfudxx, double ddfudyy, double ddfvdxy, double ddfvdxx, double ddfvdyy, double ddfudxy,
                                     double dsudx, double dsudy, double dsvdx,double dsvdy);
    double UF2D_MomentumBalanceSolid(bool x, double _f,double _s,double fu, double fv, double su, double sv, double ff, double sf, double Nr, double Nra, double ifa, double gamma, double visc, double pbs,double pbf, double SlopeX, double SlopeY,
                                     double dhsdx, double dhsdy, double dhdx, double dhdy, double dbdx, double dbdy);
    double UF1D_MomentumBalanceFluid(double _f,double _s,double fu, double su, double ff, double sf, double Nr, double Nra, double ifa, double gamma, double visc, double pbf, double Slope,
                                     double dhfdx, double dh2pbdx, double dsfdsx, double ddsfdxx, double dfudx, double ddfudxx, double dsudx);
    double UF1D_MomentumBalanceSolid(double _f,double _s,double fu, double su, double ff, double sf, double Nr, double Nra, double ifa, double gamma, double visc, double pbs, double pbf, double Slope,
                                     double dhsdx, double dhdx, double dbdx);

    ////common functions
    //dependent -> makes calls to independent common functions
    double UF_DragDistribution(double ffraction, double sfraction, double fvel, double svel);
    double UF_DragPower(double ffraction, double sfraction, double fvel, double svel);
    double UF_DragCoefficient(double ffraction, double sfraction, double gamma,double viscosity, double rocksize, double density);

    //independent
    double UF_Viscosity(double fsConc);
    double UF_VirtualMassCoeff(double ffraction, double sfraction);
    double UF_TerminalVelocity(double rocksize, double ffraction, double viscosity,double sfraction, double density);
    double UF_QuasiReynolds(double density, double viscosity, double fraction);
    double UF_Reynolds(double density, double viscosity, double ffraction,double sfraction, double rocksize);
    double UF_P(double rocksize, double ffraction, double viscosity,double sfraction, double density);
    double UF_DynamicViscosity(double sfraction);
    double UF_GetYieldStress(double sfraction);
    double UF_GetFlowResistence(double n);
    double UF_GetDispersiveResistence(double n, double sfraction);

    ////common math functions
    //temporary maps for MUSCL scheme
    cTMap * UF2D_MUSCLE_1_x1;
    cTMap * UF2D_MUSCLE_1_x2;
    cTMap * UF2D_MUSCLE_1_y1;
    cTMap * UF2D_MUSCLE_1_y2;
    cTMap * UF2D_MUSCLE_2_x1;
    cTMap * UF2D_MUSCLE_2_x2;
    cTMap * UF2D_MUSCLE_2_y1;
    cTMap * UF2D_MUSCLE_2_y2;
    cTMap * UF2D_MUSCLE_3_x1;
    cTMap * UF2D_MUSCLE_3_x2;
    cTMap * UF2D_MUSCLE_3_y1;
    cTMap * UF2D_MUSCLE_3_y2;
    cTMap * UF2D_MUSCLE_OUT_x1;
    cTMap * UF2D_MUSCLE_OUT_x2;
    cTMap * UF2D_MUSCLE_OUT_y1;
    cTMap * UF2D_MUSCLE_OUT_y2;

    cTMap * UF1D_MUSCLE_1_x1;
    cTMap * UF1D_MUSCLE_1_x2;
    cTMap * UF1D_MUSCLE_2_x1;
    cTMap * UF1D_MUSCLE_2_x2;
    cTMap * UF1D_MUSCLE_3_x1;
    cTMap * UF1D_MUSCLE_3_x2;
    cTMap * UF1D_MUSCLE_OUT_x1;
    cTMap * UF1D_MUSCLE_OUT_x2;

    static const int UF_MUSCLE_TARGET_IN1 = 1;
    static const int UF_MUSCLE_TARGET_IN2 = 2;
    static const int UF_MUSCLE_TARGET_IN3 = 3;
    static const int UF_MUSCLE_TARGET_OUT = 0;

    static const int UF_MUSCLE_mult = 1;
    static const int UF_MUSCLE_div = 2;
    static const int UF_MUSCLE_add = 3;
    static const int UF_MUSCLE_sub = 4;
    static const int UF_MUSCLE_pow = 5;

    static const int UF_DIRECTION_X = 1;
    static const int UF_DIRECTION_Y = 2;
    static const int UF_DIRECTION_XY = 3;

    void UF2D_MUSCLE(int thread,cTMap * _mask,cTMap *dt, cTMap * in, int target, int dir = UF_DIRECTION_XY);
    void UF2D_MUSCLE_operate(int thread,cTMap * _mask,cTMap *dt,int in_map1, int in_map2,int operation, int target);
    void UF2D_MUSCLE_operate(int thread,cTMap * _dem,cTMap * dt,int in_map1, int in_map2x,int in_map2y,int operation, int target);
    void UF2D_MUSCLE_operate(int thread,cTMap * _mask,cTMap *dt, int in_map,double in_double,int operation, int target);
    void UF2D_MUSCLE_operate(int thread,cTMap * _mask,cTMap * dt,int in_map,double in_double, int operation, cTMap * outx1,cTMap * outx2,cTMap * outy1,cTMap * outy2);
    void UF2D_MUSCLE_operate(int thread,cTMap * _mask,cTMap * dt,int in_map, int in_map2, int operation, cTMap * outx1,cTMap * outx2,cTMap * outy1,cTMap * outy2);
    void UF2D_MUSCLE_operate(int thread,cTMap * _mask,cTMap * dt,int in_map, int in_map2x,int in_map2y, int operation, cTMap * outx1,cTMap * outx2,cTMap * outy1,cTMap * outy2);

    void UF1D_MUSCLE(int thread,cTMap * _ldd,cTMap * _lddw,cTMap *dt,cTMap * in, int target);
    void UF1D_MUSCLE_operate(int thread,cTMap * _ldd,cTMap * _lddw,cTMap *dt,int in_map1, int in_map2,int operation, int target);
    void UF1D_MUSCLE_operate(int thread,cTMap * _ldd,cTMap * _lddw,cTMap *dt,int in_map,double in_1, int operation, int target);
    void UF1D_MUSCLE_operate(int thread,cTMap * _ldd,cTMap * _lddw,cTMap *dt,int in_map,double in_1, int operation,cTMap * outx1,cTMap * outx2);
    void UF1D_MUSCLE_operate(int thread,cTMap * _ldd,cTMap * _lddw,cTMap *dt,int in_map, int in_map2, int operation, cTMap * outx1,cTMap * outx2);

    double UF2D_MaxFlux(cTMap * _dem,cTMap * _f, int r, int c, int dr, int dc);

    //slope analysis and map derivative functions
    void UF_DEMLDDAnalysis(cTMap * _dem, cTMap * _ldd,cTMap * _lddw,cTMap * _lddh,cTMap * _f1D,cTMap * _s1D,cTMap * _f2D,cTMap * _s2D);

    static const int UF_DERIVATIVE_LR = 0;
    static const int UF_DERIVATIVE_L = 1;
    static const int UF_DERIVATIVE_R = 2;

    double UF2D_Derivative(cTMap * _dem, cTMap * _in, int r, int c, int direction, int calculationside = UF_DERIVATIVE_LR, bool useflowbarriers = false);
    double UF2D_Derivative_scaled(cTMap * _dem, cTMap * _in, int r, int c, int direction, double scale, int calculationside = UF_DERIVATIVE_LR);
    double UF1D_Derivative(cTMap * _ldd,cTMap * _lddw, cTMap * _in, int r, int c, bool minmod = false, int calculationside = UF_DERIVATIVE_LR);
    double UF1D_Derivative_scaled(cTMap * _ldd,cTMap * _lddw, cTMap * _in, int r, int c, double scale, bool minmod = false, int calculationside = UF_DERIVATIVE_LR);
    double UF2D_Derivative2(cTMap * _dem, cTMap * _in, int r, int c, int direction,int calculationside = UF_DERIVATIVE_LR);
    double UF1D_Derivative2(cTMap * _ldd,cTMap * _lddw, cTMap * _in, int r, int c, int calculationside = UF_DERIVATIVE_LR);

    double UF_DEMACCES(cTMap * dem, cTMap * h,int r, int c, int dr, int dc);
    double UF_5CellAverage(cTMap * m,int r, int c);


    double FluxLimiter2D();
    double FluxLimiter1D();
    double FL_MinMod();
    double FL_HLL2();
    double UF_MinMod(double a, double b);

    bool UF_OUTORMV(cTMap * mask, int r, int c);
    bool UF_LDDOUT(cTMap * mask, int r, int c, bool front);
    bool UF_NOTIME(cTMap * mask,cTMap * dt, int r, int c);

    void UF1D_AddValue(cTMap * _ldd,cTMap * _lddw,int r, int c, bool front, cTMap * _in, double add, bool zeromin = false);
    double UF1D_Value(cTMap * _ldd,cTMap * _lddw,int r, int c,bool front, cTMap * _in);

    void upstream(cTMap *_LDD, cTMap *_M, cTMap *out);

    ////Unified Flow Soil Interactions
    //General Function
    void UnifiedFlowSediment(int thread);
    void UF_SedimentSource(int thread);
    void UF_FlowDetachment(int thread);
    void UF_FlowEntrainment(int thread);

    double UF_FlowEntrainmentST(int r, int c, bool channel);
    double UnifiedFlowActiveEntrainmentST(double slope, double _f, double _s,double area, double _fv, double _sv, double _sc, double visc, double d, double ifa,double rocksize, double d_bed, double ifa_bed, double coh_bed, double veg_coh, double manning, int r, int c);
    void UF_LateralEntrainment(int thread);
    double UF_LateralEntrainmentArea(int r, int c, int dr, int dc);

    void UF_FlowDetachment(double dt, int r, int c,int d, bool channel);
    void UF_FlowEntrainment(double dt, int r, int c, bool channel);

    void UF_FlowCompaction(int thread);
    void UF_FlowCompaction(double dt, int r, int c, bool channel);
    void UF_FlowCompactionThreshold(double dt, int r, int c, bool channel);

    double UF_SoilTake(int r, int c, int d, double potential,bool channel,bool bedload);
    void UF_SoilAdd(int r, int c, int d, double mass, bool channel);

    double UF_RockTake(int r, int c, double entrainment, bool channel);
    double UF_RockAdd(int r, int c, double entrainment, bool channel);

    void UF_SumGrainClasses();

    cTMap *UF2D_EntrainmentSF;

    //transport capacity
    double UnifiedFlowTransportCapacity(int r, int c, int d, bool channel, bool bedload);
    //active entrainment
    double UnifiedFlowActiveEntrainment(double dt,double st,double slope, double _f, double _s,double area, double _fv, double _sv, double _sc, double visc, double d, double ifa,double rocksize, double d_bed, double ifa_bed, double coh_bed, double veg_coh, double manning, int r, int c);
    double UnifiedFlowActiveEntrainmentLat(double dt,double st, double slope, double h, double _f, double _s,double area, double _fv, double _sv, double _sc, double visc, double d, double ifa,double rocksize, double d_bed, double ifa_bed, double coh_bed, double veg_coh, double manning, int r, int c);

    double UnifiedFlowActiveDeposition(double dt,double slope, double _f, double _s,double area, double _fv, double _sv, double _sc, double visc, double d, double ifa,double rocksize, double d_bed, double ifa_bed, int r, int c);

    double UnifiedFlowEntrainmentAvailableDepth(int r,int c, double vx, double vy);
    double UnifiedFlowDepositionAvailableDepth(int r,int c);

    //collapse of entrainment sides
    double UF_EntrainmentSideSlopeFailure(double dt, int r, int c);


    //connection to the digital elevation model
    void UFDEMLDD_Connection(int thread);

    ////Initial and Forcing stuff
    cTMap * UF2D_ForcedFVolume;
    cTMap * UF2D_ForcedSVolume;
    cTMap * UF2D_ForcedSDensity;
    cTMap * UF2D_ForcedSIFA;
    cTMap * UF2D_ForcedSRocksize;

    cTMap * UF2D_Initialized;
    cTMap * UF2D_InitialTime;
    cTMap * UF2D_InitialFVolume;
    cTMap * UF2D_InitialSVolume;
    cTMap * UF2D_InitialSDensity;
    cTMap * UF2D_InitialSIFA;
    cTMap * UF2D_InitialSRocksize;

    bool UF_NeedsInitial = false;
    void UF_ForcedConditions(int thread,cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                           cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                           cTMap * _fu1D,cTMap * _s1D,
                           cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                           cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                           cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                           cTMap * _su2D,cTMap * _sv2D);

    void UF_Initial( cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                           cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                           cTMap * _fu1D,cTMap * _s1D,
                           cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                           cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                           cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                           cTMap * _su2D,cTMap * _sv2D);

    double UF_InitializedF;
    double UF_InitializedS;

    cTMap * SourceSolid;
    cTMap * SourceSolidDensity;
    cTMap * SourceSolidRocksize;
    cTMap * SourceSolidIFA;
    cTMap * SourceFluid;

    cTMap * ChannelSourceSolid;
    cTMap * ChannelSourceSolidDensity;
    cTMap * ChannelSourceSolidRocksize;
    cTMap * ChannelSourceSolidIFA;
    cTMap * ChannelSourceFluid;

    void AddSource(int r, int c, double f, double s, double d, double rocksize, double ifa, bool channel);

    LisemThreadPool * ThreadPool;

////END FILE
