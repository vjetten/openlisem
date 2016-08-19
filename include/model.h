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
  \file model.h:
  \brief main class, describing the model world with all variables and processes
*
* - TWorld class that combines ALL model variables
* - all processes
* - global defines Drc, MV, FOR_ROW_COL_MV etc
* - global defines for lisem type; infiltration type etc.
*/

#ifndef modelH
#define modelH

//#include <math.h>
//#include <stdlib.h>

#include <QtGui>
#include <QMutex>

#include "CsfMap.h"
#include "io.h"
#include "error.h"
#include "swatre_p.h"
#include "swatre_g.h"


#define OLDSWATRE 1

//---------------------------------------------------------------------------
#define PI 3.14159265

#define DEBUG(s) emit debug(QString(s))

//#define mwrite(name) writeRaster(QString(resultDir+name))
#define report(raster, name) WriteMapSeries(raster, resultDir,QString(name), printstep)

// defines to make life easier

/// shortcut to access data
#define Drc     data[r][c]
#define Drcd    at(d)->data[r][c]
#define Drci data[r+dr[i]][c+dc[i]]

/// shortcut to access the outlet point data
#define DrcOutlet     data[r_outlet][c_outlet]

/// shortcut missing value in map
#define MV(r,c) pcr::isMV(LDD->data[r][c])

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_MV for(int r = 0; r < _nrRows; r++)\
    for (int c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(LDD->data[r][c]))

#define FOR_GRAIN_CLASSES for(int d  = 0 ; d < numgrainclasses;d++)

#define FOR_CELL_IN_FLOODAREA for (long _i = 0; _i < nrFloodcells ; _i++)\
{\
    int r = floodRow[_i];\
    int c = floodCol[_i];

/// shortcut for all cell in watershed with nr wsnr
#define FOR_WATERSHED_ROW_COL(wsnr) for (long k = 0; k < WS[wsnr].cr.count(); k++) {\
    int c = WS[wsnr].cr[k]._c;\
    int r = WS[wsnr].cr[k]._r;\

/// shortcut for channel row and col loop
#define FOR_ROW_COL_MV_CH for (int  r = 0; r < _nrRows; r++)\
    for (int  c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(LDDChannel->data[r][c]))

/// shortcut for tile network row and col loop.
#define FOR_ROW_COL_MV_TILE for (int  r = 0; r < _nrRows; r++)\
    for (int  c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(LDDTile->data[r][c]))

/// shortcut to check if r,c is inside map boundaries, used in kinematic and flooding
#define INSIDE(r, c) (r>=0 && r<_nrRows && c>=0 && c<_nrCols)


/// shortcut for LDD row and col loop
#define FOR_ROW_COL_UF2D for(int r = 0; r < _nrRows; r++)\
    for (int c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(_dem->data[r][c]) )

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_UF2D_DT for(int rc = 0; rc < _nrRows; rc++)\
    {for (int cc = 0; cc < _nrCols; cc++)\
    {int r = (int) (UF2D_CellR->data[rc][cc]);\
    int c = (int) (UF2D_CellC->data[rc][cc]);\
    if(!INSIDE(r,c)){break;}\
    if(!pcr::isMV(_dem->data[r][c]) && !(dt->Drc == 0))

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_UF1D for(int r = 0; r < _nrRows; r++)\
    for (int c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(_ldd->data[r][c]) )

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_UF1D_DT for(int r = 0; r < _nrRows; r++)\
    for (int c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(_ldd->data[r][c]) && !(dt->Drc == 0))



// check if cell From flows to To
#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )


#define MAX_ITERS 50

/*
  local drain direction maps have values for directions as follows:
    7  8  9
     \ | /
   4 - 5 - 6
     / | \
    1  2  3
 */

#define NUMNAMES 2000   /// \def NUMNAMES runfile namelist max
#define NUMMAPS 1000    /// \def max nr maps
#define MIN_FLUX 1e-12 /// \def minimum flux (m3/s) in kinematic wave
#define MIN_HEIGHT 1e-6 /// \def minimum water height (m) for transport of sediment
#define MAXCONC 848.0    /// \def max concentration susp. sed. in kg/m3 0.32 * 2650 = max vol conc from experiments Govers x bulk density
#define MAXCONCBL 848.0    /// \def max concentration susp. sed. in kg/m3 0.32 * 2650 = max vol conc from experiments Govers x bulk density
#define UF_VERY_SMALL 1e-8 /// \def min timestep/flux/height in unified flow equations

#define INFIL_NONE 0
#define INFIL_SWATRE 1
#define INFIL_HOLTAN 2
#define INFIL_GREENAMPT 3
#define INFIL_GREENAMPT2 4
#define INFIL_KSAT 5
#define INFIL_MOREL 21
#define INFIL_SMITH 22
#define INFIL_SMITH2 23

#define KE_EXPFUNCTION 0
#define KE_LOGFUNCTION 1
#define KE_POWERFUNCTION 2

#define MINMOD 1
#define VANALBEDA 2
#define VANLEER 3

#define FMUSCL 1
#define FENO 2
#define FENOMOD 3

#define FSGOVERS 0
#define FSRIJN 1
#define FSRIJNFULL 2
#define FSWUWANGJIA 3

#define RGOVERS 0
#define RRIJN 1
#define RRIJNFULL 2
#define RWUWANGJIA 3

#define OFGOVERS 0
#define OFHAIRSINEROSE 1

#define K1D_METHOD       1
#define K2D_METHOD_FLUX  2
#define K2D_METHOD_INTER 3


//---------------------------------------------------------------------------
/// structure containing pointers to all maps

/** structure containing pointers to all maps for automatic destruction after runs
 so memory doesn't have to be freed for each map. The functions Newmap(double) and
ReadMap(cTMap *Mask, QString name) put a map on this list
*/
typedef struct MapListStruct {
    cTMap *m;
}  MapListStruct;
//---------------------------------------------------------------------------list
/// linked list structure for network in kin wave
typedef struct LDD_LINKEDLIST {
    int rowNr;
    int colNr;
    struct LDD_LINKEDLIST *prev;
}  LDD_LINKEDLIST;
//---------------------------------------------------------------------------
/// name list structure used to read run file
typedef struct NAME_LIST {
    QString name;
    QString value;
} NAME_LIST;
//---------------------------------------------------------------------------
/// structure for output of land unit stats
typedef struct UNIT_LIST {
    long nr;
    double var0;
    double var1;
    double var2;
    double var3;
    double var4;
    double var5;
} UNIT_LIST;
//---------------------------------------------------------------------------
/// structure for output of land unit stats
typedef struct COORD {
    int _r;
    int _c;
} COORD;
//--------------------------------------------------------------------------
/// structure for watershed coordinates for flooding
typedef struct WS_LIST {
    bool flood;
   QList <COORD> cr;
   //  QVector <COORD> cr;
    int ws;
    double dt;
    double dt2;
    double dtsum;
} WS_LIST;
//---------------------------------------------------------------------------
/// Strunture to store rain station values of rainfile mapnames
typedef struct RAIN_LIST {
    double time;
    //double *intensity;
    QVector <double> intensity;
    bool isMap;
    QString name;
} RAIN_LIST;
//---------------------------------------------------------------------------

/// \class TWorld model.h contains the model 'World': constants, variables and erosion processes

/** The model 'world': the main class containing all variables, maps, options, filenames.\n
  The class contains hydrological and erosion processes which are run in a time loop.\n
  The main function is <B>void TWorld::DoModel()</B>, which has the time loop calling ll processes\n
  Every timestep the mass balance is calculated and output is reported to the UI and disk.
 */

//http://blog.exys.org/entries/2010/QThread_affinity.html
//http://thesmithfam.org/blog/2009/09/30/lock-free-multi-threading-in-qt/

class TWorld: public QThread
{
    Q_OBJECT

public:
    TWorld(QObject *parent = 0);
    ~TWorld();

    /// copy of overall rows and columns, set in initmask
    long _nrRows;
    long _nrCols;

    /// map management structure, automatic adding and deleting of all cTMap variables
    MapListStruct maplistCTMap[NUMNAMES];
    int maplistnr;

    /// variable declaration list of all maps with comments:
#include "TMmapVariables.h"

    /// SwitchXXX are boolean options that are set in interface and runfile, mainly corrsponding to checkboxes in the UI
    bool SwitchRoadsystem, SwitchHardsurface, SwatreInitialized, SwitchInfilGA2, SwitchLimitTC, SwitchLimitDepTC,
    SwitchWheelPresent, SwitchCompactPresent, SwitchIncludeChannel, SwitchChannelBaseflow,
    startbaseflowincrease, SwitchChannelInfil, SwitchAllinChannel, SwitchErosion, SwitchAltErosion,
    SwitchSimpleDepression, SwitchBuffers, SwitchBuffersImpermeable, SwitchSedtrap, SwitchSnowmelt, SwitchRainfall, SwitchRunoffPerM, SwitchInfilCompact,
    SwitchInfilCrust, SwitchGrassStrip, SwitchImpermeable, SwitchPercolation, SwitchDumphead, SwitchWaterRepellency,
    SwitchWheelAsChannel, SwitchMulticlass, SwitchNutrients, SwitchGullies, SwitchGullyEqualWD, SwitchGullyInfil,
    SwitchGullyInit, SwitchOutputTimeStep, SwitchOutputTimeUser, SwitchMapoutRunoff, SwitchMapoutConc,
    SwitchMapoutWH, SwitchMapoutWHC, SwitchMapoutTC, SwitchMapoutEros, SwitchMapoutDepo, SwitchMapoutV,
    SwitchMapoutInf, SwitchMapoutSs, SwitchMapoutChvol, SwitchWritePCRnames, SwitchWriteCommaDelimited, SwitchWritePCRtimeplot,
    SwitchNoErosionOutlet, SwitchDrainage, SwitchPestout, SwitchSeparateOutput,
    SwitchInterceptionLAI, SwitchTwoLayer, SwitchSimpleSedKinWave, SwitchSOBEKoutput,
    SwitchPCRoutput, SwitchWriteHeaders, SwitchGeometric, SwitchIncludeTile, SwitchKETimebased, SwitchHouses, SwitchChannelFlood, SwitchRaindrum,
    Switchheaderpest, SwitchPesticide, SwitchRainfallFlood, SwitchFloodSedimentMethod, SwitchStoninessDET,
    SwitchFloodExplicit, SwitchFloodSWOForder1, SwitchFloodSWOForder2, SwitchMUSCL, SwitchLevees, SwitchFloodInitial, SwitchWatershed;

    int SwitchFlood1D2DCoupling;
    int SwitchKinematic2D;
    int SwitchEfficiencyDET;

    // multiple options that are set in interface or runfile, see defines above

    /// Interception storage function based on LAI
    int InterceptionLAIType;

    /// infiltration method
    int InfilMethod;

    /// erosion units in output: to/ha; kg/cell; kg/m2
    int ErosionUnits;

    /// type of kinetic energy equation;
    int KEequationType;

    /// parameters in KE equations
    double KEParamater_a1, KEParamater_b1, KEParamater_c1;
    double KEParamater_a2, KEParamater_b2;
    double KEParamater_a3, KEParamater_b3;

    /// calibration parameters
    double ksatCalibration;
    double nCalibration;
    double thetaCalibration;
    double psiCalibration;
    double ChnCalibration;
    double ChKsatCalibration;
    double SplashDelivery;
    double DepositedCohesion;
    double StripN;
    double StemflowFraction;

    double CanopyOpeness; // VJ 110209 added Aston factor as user input
    double waterRep_a;
    double waterRep_b;
    double waterRep_c;
    double waterRep_d;

    ///rainfall to flood max gradient
 //   double rainFloodingGradient;

    /// totals for mass balance checks and output
    /// Water totals for mass balance and output (in m3)
    double MB, MBeM3, Qtot,QtotT, QtotOutlet, IntercTot, WaterVolTot, floodVolTot, floodVolTotInit, floodVolTotMax, floodAreaMax, WaterVolSoilTot, InfilTot, RainTot, SnowTot, SurfStoremm, InfilKWTot,BaseFlowTot,BaseFlow;
    double floodBoundaryTot;
    //houses
    double IntercHouseTot, IntercHouseTotmm;
    double ChannelVolTot, ChannelSedTot, ChannelDepTot, ChannelDetTot, TileVolTot;
    /// Sediment totals for mass balance and output (in kg)
    double MBs, DetTot, DetSplashTot, DetFlowTot, DepTot, SoilLossTot,SoilLossTotT, SoilLossTotOutlet, SedTot, SoilLossTotSub,
           FloodDetTot, FloodDepTot, FloodSedTot;
    /// Water totals for output in file and UI (in mm), copied to 'op' structure
    double RainTotmm, SnowTotmm, IntercTotmm, WaterVolTotmm,WaterVolRunoffmm, InfilTotmm, Qtotmm, RainAvgmm, SnowAvgmm;
    double floodTotmm, floodTotmmInit;
    /// peak times (min)
    double RainstartTime, RainpeakTime, SnowpeakTime, QpeakTime, Qpeak, Rainpeak, Snowpeak;
    bool rainStarted;
    double BufferVolTot, BufferSedTot, BufferVolTotInit, BufferSedTotInit, BulkDens, BufferVolin;
    double nrCells, CatchmentArea, nrFloodedCells;

    ///pesticides
    double MBp,PestMassApplied, PestLossTotOutlet, PestFluxTotOutlet, PestRunoffSpatial, PestDisMixing, PestSorMixing, PestInfilt, PestStorage, Pestdetach, PestCinfilt,PestCfilmexit;
    double MBpex,PestRunoffSpatialex,PestDisMixingex,PestSorMixingex,PestInfiltex,PestLossTotOutletex;
    int N_SPK;
    double Maxsolubility;
    double MaxVup;

    int c_outlet;  /// copy of outlet col number
    int r_outlet;  /// copy of outlet row number

    int c_plot;  /// copy of col number of hydrograph plotted on screen
    int r_plot;  /// copy of row number of hydrograph plotted on screen

    /// time and dx parameters
    double time, BeginTime, EndTime;
    double _dt, _dx;
    long runstep, printstep, printinterval;
    int fZone;

    /// timeseries variables and output strings
    int nrRainfallseries;
    int nrSnowmeltseries;
    QVector <RAIN_LIST> RainfallSeriesM;  // rainfall vector of records
    QVector <RAIN_LIST> SnowmeltSeriesM;

    // output formatting for SOBEK flood model input
    QString SOBEKdatestring;
    int SOBEKnrlines;

    // file and directory names
    QString resultDir;
    QString inputDir;
    QString outflowFileName;
    QString totalErosionFileName;
    QString totalDepositionFileName;
    QString totalChanErosionFileName;
    QString totalChanDepositionFileName;
    QString totalSoillossFileName;
    QString totalSedFileName;
    QString totalLandunitFileName;

    QString rainfallMapFileName;
    QString interceptionMapFileName;
    QString infiltrationMapFileName;
    QString runoffMapFileName;
    QString runoffFractionMapFileName;
    QString channelDischargeMapFileName;

    QString floodLevelFileName;
    QString floodTimeFileName;
    QString floodStatsFileName;
    QString floodMaxQFileName;
    QString floodMaxChanWHFileName;
    QString floodFEWFileName;
    QString floodMaxVFileName;
    QString floodWHmaxFileName;
    QString timestamp;

    QString rainFileName;
    QString rainFileDir;
    QString snowmeltFileName;
    QString snowmeltFileDir;
    QString resultFileName;
    QString temprunname;
    QStringList outputcheck;
    /// standard names of output map series
    QString Outrunoff, Outconc, Outwh, Outrwh, Outvelo, Outinf, Outss, Outchvol,
    Outtc, Outeros, Outdepo, OutSL, OutSed,
    OutTiledrain, OutHmx, OutVf, OutQf, OutHmxWH;
    QString errorFileName;


    // list with class values of land unit map
    UNIT_LIST unitList[512]; // just a fixed number for 512 classes, who cares!
    UNIT_LIST floodList[512]; // just a fixed number for 512 classes, who cares!
    QVector <UNIT_LIST> unitListM;
    int landUnitNr;
    // data initialization, runfile reading and parsing
    NAME_LIST runnamelist[NUMNAMES]; // structure for runfile variables and names
    int nrrunnamelist;


    QList <WS_LIST> WS;
    QList <COORD> FA;


    // list of pointers for substance maps: sediment, sed classes, nutrients etc.
    // used in kin wave for routing of substances
    MapListStruct SubsMaps[32];
    int nrSubsMaps;

    // functions in lisDataInit.cpp
    void InitMapList(void);
    cTMap *NewMap(double value);
    cTMap *ReadMap(cTMap *Mask, QString name);
    void DestroyData(void);
    cTMap *InitMask(QString name);
    cTMap *InitMaskChannel(QString name);
    cTMap *InitMaskTiledrain(QString name);
    void InitTiledrains(void); //VJ 110112
    void InitBuffers(void); //VJ 110112
    void InitChannel(void); //VJ 110112
    void InitShade(void); //VJ 130301
    void InitMulticlass(void); //VJ 110511
    void GetInputData(void);      // get and make input maps
    void IntializeData(void);     // make all non-input maps
    void IntializeOptions(void);  // set all options to false etc

    // functions in lisRunfile.cpp
    QString getvaluename(QString vname);
    double getvaluedouble(QString vname);
    int getvalueint(QString vname);
    QString CheckDir(QString p, bool makeit = false);
    QString GetName(QString p);
    QString checkOutputMapName(QString p, QString S, int i);
    void ParseRunfileData(void);
    void GetRunFile(void);
    //MapListStruct qx[9];


    //SEDIMENT TRANSPORT
    int SS_Method;
    int BL_Method;
    double SigmaDiffusion;

    int OF_Method;

    //int GrainSizeDistributionType;

    bool SwitchUseMaterialDepth,SwitchUse2Layer,SwitchUseGrainSizeDistribution, SwitchEstimateGrainSizeDistribution,SwitchReadGrainSizeDistribution;
//SwithEstimated90,
    int numgrainclasses;
    QString GrainMaps;
    QList<double> graindiameters;
    QList<double> settlingvelocities;
    double distD50;
    double distD90;

    double LogNormalDist(double d50,double sigma, double d);
    double DetachMaterial(int r,int c, int d,bool channel,bool flood,bool bl, double detachment);
    void SedimentSetMaterialDistribution();//(int r,int c);
    QList<cTMap *> IW_D;

    QList<cTMap *> W_D;
    QList<cTMap *> RW_D;


    QList<cTMap *> StorageDep_D;
    QList<cTMap *> Storage_D;
    QList<cTMap *> RStorageDep_D;
    QList<cTMap *> RStorage_D;

    cTMap *unity;


    double GetDpMat(int r, int c,double p,QList<cTMap *> *M);
    double GetMpMat(int r, int c,double p,QList<cTMap *> *M, QList<double> *V);
    double GetDp(int r, int c,double p);
    double GetTotalDW(int r, int c,QList<cTMap *> *M);
    double GetSV(double d);
    void SplashDetachment(void);
    void FlowDetachment(void);
    double MaxConcentration(double watvol, double sedvol);
    void ChannelFlowDetachment(int r, int c);

    void SumSedimentClasses();

    void FindBaseFlow(); //search for channel inflow from groundwater
    bool addedbaseflow = false;

    void Pestmobilisation(void);
//    void TransPesticide(int pitRowNr, int pitColNr,cTMap *_LDD,cTMap *_Qn, cTMap *_Vup, cTMap *_Vupold,cTMap *_WHoutavg,
//                         cTMap *_WHoutavgold,cTMap *_RainNet,cTMap *_CM_N,cTMap *_C_N,cTMap *_CS_N,cTMap *_InfilVol,cTMap *_InfilVolold,
//                         cTMap *_DX,cTMap *_C,cTMap *_Cold,cTMap *_CS,cTMap *_CM,cTMap *_Kfilm,cTMap *_epsil,
//                         cTMap *_KD,cTMap *_poro,cTMap *_rhob,cTMap *_kr,cTMap *_Qin, cTMap *_Sin,cTMap *_Q,cTMap *_Alpha,cTMap *_Qpn);

    double cmx_analytique(double t, double dKfi, double dpestiinf, double depsil, double drhob, double dkr, double dKD, double dn, double CM0, double CS0,double Cr);
    double csx_analytique(double t, double dKfi,double dpestiinf,double depsil,double drhob,double dkr,double dKD,double dn, double CM0,double CS0,double Cr);
    double **Factorize(double **A, int n, int m);
    double *Solve(int n,int m, double **A_LU, double *B);
    double Implicitscheme(double Qj1i1, double Qj1i, double Qji1,double Pj1i, double Pji1, double alpha, double dt,double dx, double Kfilm, double CMi1j1);
    double ConcentrationP(double watvol, double pest);

    // 1D hydro processes
    //input timeseries
    void GetRainfallDataM(QString name, bool israinfall);   // get input timeseries
    //void GetSnowmeltData(void);   // get input timeseries
    /// convert rainfall of a timestep into a map
    void RainfallMap(void);
    /// convert snowmelt of a timestep into a map
    void SnowmeltMap(void);
    /// interception of vegetation canopy resulting in rainnet
    void Interception(void);
    /// infiltration function calling all infiltration methods
    /// add rainnet to WH and calculating new WH
    void InterceptionLitter(void);
    /// subtract water retained in Litter under forest e.g.
    void InterceptionHouses(void);
    /// subtract water retained on houses, for urban projects    
    void addRainfallWH(void);
    /// add net rainfall to WH, WHroads and WHgrass

    void InfilEffectiveKsat();
    void Infiltration(void);

    void InfilSwatre(cTMap *_WH);
//    void InfilGreenAmpt1(cTMap *_WH);  //OBSOLETE
//    void InfilSmithParlange1(cTMap *_WH); //OBSOLETE
    void InfilMorelSeytoux1(cTMap *_WH);
//    void InfilKsat(cTMap *_WH); //OBSOLETE
    double IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p, REAL8 *FFull);
    void SoilWater(void);
    void InfilMethods(cTMap *_Ksateff, cTMap *_WH, cTMap *_fpot, cTMap *_fact, cTMap *_L1, cTMap *_L2, cTMap *_FFull);
    void RainfallToFlood(void);
    void SurfaceStorage(void);

    void ToTiledrain(void);
    void TileFlow(void);
    void CalcVelDischTile(void);
    void GridCell(void);

    //flood
    //QVector<int> cellRow;
    //QVector<int> cellCol;
    int *cellRow;
    int *cellCol;
    int *floodRow;
    int *floodCol;
    long nrGridcells;
    long nrFloodcells;
    void ChannelFlood(void);
    void FloodMaxandTiming(void);
    void FloodBoundary(void);
    void FloodSpuriousValues(void);
    void ChannelFloodStatistics(void);

    double courant_factor;
    double courant_factor_diffusive;
    double mixing_coefficient, runoff_partitioning;
   // double cfl_fix;
    double minReportFloodHeight;
    double correctMassBalance(double sum1, cTMap *M, double minV);

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
    double UF_MANNINGCOEFFICIENT;

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
    static const int UF_SCHEME_CENTRALSIMPLE =1;
    static const int UF_SCHEME_BOUNDARYMUSCLE =2;

    cTMap * UF2D_Test;

    //just for display
    cTMap * UF2D_h;
    cTMap * UF2D_q;
    cTMap * UF2D_qs;
    cTMap * UF2D_fsConc;
    cTMap * UF2D_sConc;
    cTMap * UF2D_tConc;
    cTMap * UF2D_velocity;
    cTMap * UF2D_u;
    cTMap * UF2D_v;

    cTMap * UF1D_h;
    cTMap * UF1D_q;
    cTMap * UF1D_qs;
    cTMap * UF1D_fsConc;
    cTMap * UF1D_sConc;
    cTMap * UF1D_tConc;
    cTMap * UF1D_velocity;


    //internal slope functions
    cTMap * UF2D_Slope;
    cTMap * UF2D_SlopeX;
    cTMap * UF2D_SlopeY;
    cTMap * UF1D_LDDs;

    //actual calculation variables
    ////2D
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
    cTMap * UF2D_sqx1;
    cTMap * UF2D_sqy1;
    cTMap * UF2D_sqx2;
    cTMap * UF2D_sqy2;

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
    cTMap * UF1D_LDD;
    cTMap * UF1D_LDDw;
    cTMap * UF1D_LDDh;
    cTMap * UF1D_T;
    cTMap * UF1D_DT;
    cTMap * UF1D_Slope;
    cTMap * UF1D_DTStep;
    cTMap * UF1D_Courant;
    //fluid phase
    cTMap * UF1D_f;
    cTMap * UF1D_fstore;
    cTMap * UF1D_visc;
    cTMap * UF1D_fu;
    cTMap * UF1D_fa;
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

    //temporary maps for generic advection functions
    cTMap * UF_t1;
    cTMap * UF_t2;
    cTMap * UF_t3;
    cTMap * UF_t4;
    cTMap * UF_t5;

    double UF2D_foutflow = 0;
    double UF2D_fsoutflow = 0;
    double UF2D_soutflow = 0;
    double UF1D_foutflow = 0;
    double UF1D_fsoutflow = 0;
    double UF1D_soutflow = 0;


    //General Function
    void UnifiedFlow();
    void UF_Init();
    void UF_SetInput();
    void UF_SetOutput();

    bool UF_Input_first = true;

    //2D version
    double UF2D_Source(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv);

    double UF2D_Scheme(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv);

    void UF2D_Advect_Momentum(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_fu, cTMap * out_fv,  cTMap * out_su, cTMap * out_sv );
    double UF2D_Advect_mass(cTMap* dt, cTMap * _dem,cTMap * _m, cTMap * _mu, cTMap * _mv, cTMap * out_m);
    void UF2D_Advect_prop(cTMap* dt, cTMap * _dem,cTMap * _m, cTMap * _mu, cTMap * _mv,cTMap *_prop, cTMap * out_prop = 0);

    void UF2D_Advect2_Momentum(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_fu, cTMap * out_fv,  cTMap * out_su, cTMap * out_sv,cTMap * out_qfx1,cTMap *out_qfx2,cTMap * out_qfy1,cTMap *out_qfy2,cTMap * out_qsx1,cTMap *out_qsx2,cTMap * out_qsy1,cTMap *out_qsy2);
    double UF2D_Advect2_mass(cTMap* dt, cTMap * _dem,cTMap * _m, cTMap * f,cTMap * _qx1, cTMap * _qy1,cTMap * _qx2, cTMap * _qy2, cTMap * out_m);
    void UF2D_Advect2_prop(cTMap* dt, cTMap * _dem,cTMap * _m, cTMap * f,cTMap * _qx1, cTMap * _qy1,cTMap * _qx2, cTMap * _qy2,cTMap *_prop, cTMap * out_prop = 0);

    void UF2D_FluidSource(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_f);
    void UF2D_SolidSource(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv, cTMap * out_s);

    void UF2D_FluidMomentumSource(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv);
    void UF2D_SolidMomentumSource(cTMap* dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv);

    void UF2D_FluidApplyMomentum(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_fu, cTMap * out_fv);
    void UF2D_SolidApplyMomentum(cTMap * dt, cTMap * _dem,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _fv,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * _sv,cTMap * out_su, cTMap * out_sv);

    void UF2D_Diffuse_mass(cTMap* dt, cTMap * _dem,cTMap * _m,cTMap * _f, cTMap * _s, cTMap * _fu, cTMap * _fv, cTMap * _su, cTMap * _sv, cTMap * out_m);

    void UF2D_Stored_mass(cTMap * dt, cTMap* dem, cTMap *_f,cTMap * _s, cTMap * out_f, cTMap * out_s);
    //1D version
    void UF1D_Source(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su);

    void UF1D_Scheme(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su);

    void UF1D_Advect_Momentum(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_fu,  cTMap * out_su );
    double UF1D_Advect_mass(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _m, cTMap * _mu, cTMap * out_m);
    void UF1D_Advect_prop(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _m, cTMap * _mu,cTMap *_prop, cTMap * out_prop = 0);

    void UF1D_Advect2_Momentum(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_fu,  cTMap * out_su , cTMap * out_fq1, cTMap * out_fq2,cTMap * out_sq1, cTMap * out_sq2);
    double UF1D_Advect2_mass(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _m, cTMap * _f, cTMap * _q1,cTMap * _q2, cTMap * out_m);
    void UF1D_Advect2_prop(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _m, cTMap * _f, cTMap * _q1,cTMap * _q2,cTMap *_prop, cTMap * out_prop = 0);

    void UF1D_FluidSource(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_f);
    void UF1D_SolidSource(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su, cTMap * out_s);

    void UF1D_SolidMomentumSource(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu);
    void UF1D_FluidMomentumSource(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_su);

    void UF1D_SolidApplyMomentum(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_su);
    void UF1D_FluidApplyMomentum(cTMap * dt, cTMap * _ldd,cTMap * _lddw,cTMap *_lddh,cTMap * _f,cTMap * _visc,cTMap * _fu,cTMap * _s,cTMap * _d,cTMap * _ifa,cTMap * _rocksize,cTMap * _su,cTMap * out_fu);
    void UF1D_Diffuse_mass(cTMap* dt, cTMap * _ldd,cTMap * _lddw,cTMap * _m, cTMap * _f,cTMap * _fu,cTMap * _s,cTMap * _su, cTMap * out_m);

    void UF1D_Stored_mass(cTMap * dt, cTMap * _ldd,cTMap * _lddw, cTMap *_f,cTMap * _s, cTMap * out_f, cTMap * out_s);

    ////boundary conditions
    double UF_BoundaryFlux2D(double dt, double cellx, double celly, double f, double s, double fu, double fv, double su, double sv, double slopeX, double slopeY,double NN, int dr, int dc );
    double UF_BoundaryFlux1D(double dt, double width, double f, double s, double fu, double su, double slope, double NN, bool front );

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
    void UF2D1D_Connection(cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                           cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                           cTMap * _fu1D,cTMap * _s1D,
                           cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                           cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                           cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                           cTMap * _su2D,cTMap * _sv2D);

    void UF2D1D_ChannelWater(cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                                     cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                                     cTMap * _fu1D,cTMap * _s1D,
                                     cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                                     cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                                     cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                                     cTMap * _su2D,cTMap * _sv2D);

    void UF2D1D_Infiltration(cTMap * dt, cTMap * _dem,cTMap * _ldd,cTMap * _lddw,
                             cTMap * _lddh,cTMap * _f1D,cTMap * _visc1D,
                             cTMap * _fu1D,cTMap * _s1D,
                             cTMap * _d1D,cTMap * _ifa1D,cTMap * _rocksize1D,cTMap * _su1D,
                             cTMap * _f2D,cTMap * _visc2D,cTMap * _fu2D,
                             cTMap * _fv2D,cTMap * _s2D,cTMap * _d2D,cTMap * _ifa2D,cTMap * _rocksize2D,
                             cTMap * _su2D,cTMap * _sv2D);

    ////common functions
    double UF_Friction(double dt,double velx,double vely, double NN, double h);

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

    void UF2D_MUSCLE(cTMap * _mask,cTMap *dt, cTMap * in, int target, int dir = UF_DIRECTION_XY);
    void UF2D_MUSCLE_operate(cTMap * _mask,cTMap *dt,int in_map1, int in_map2,int operation, int target);
    void UF2D_MUSCLE_operate(cTMap * _dem,cTMap * dt,int in_map1, int in_map2x,int in_map2y,int operation, int target);
    void UF2D_MUSCLE_operate(cTMap * _mask,cTMap *dt, int in_map,double in_double,int operation, int target);
    void UF2D_MUSCLE_operate(cTMap * _mask,cTMap * dt,int in_map,double in_double, int operation, cTMap * outx1,cTMap * outx2,cTMap * outy1,cTMap * outy2);
    void UF2D_MUSCLE_operate(cTMap * _mask,cTMap * dt,int in_map, int in_map2, int operation, cTMap * outx1,cTMap * outx2,cTMap * outy1,cTMap * outy2);
    void UF2D_MUSCLE_operate(cTMap * _mask,cTMap * dt,int in_map, int in_map2x,int in_map2y, int operation, cTMap * outx1,cTMap * outx2,cTMap * outy1,cTMap * outy2);

    void UF1D_MUSCLE(cTMap * _ldd,cTMap * _lddw,cTMap *dt,cTMap * in, int target);
    void UF1D_MUSCLE_operate(cTMap * _ldd,cTMap * _lddw,cTMap *dt,int in_map1, int in_map2,int operation, int target);
    void UF1D_MUSCLE_operate(cTMap * _ldd,cTMap * _lddw,cTMap *dt,int in_map,double in_1, int operation, int target);
    void UF1D_MUSCLE_operate(cTMap * _ldd,cTMap * _lddw,cTMap *dt,int in_map,double in_1, int operation,cTMap * outx1,cTMap * outx2);
    void UF1D_MUSCLE_operate(cTMap * _ldd,cTMap * _lddw,cTMap *dt,int in_map, int in_map2, int operation, cTMap * outx1,cTMap * outx2);

    //slope analysis and map derivative functions
    void UF_DEMLDDAnalysis(cTMap * _dem, cTMap * _ldd,cTMap * _lddw,cTMap * _lddh,cTMap * _f1D,cTMap * _s1D,cTMap * _f2D,cTMap * _s2D);

    double UF2D_Derivative(cTMap * _dem, cTMap * _in, int r, int c, int direction);
    double UF1D_Derivative(cTMap * _ldd,cTMap * _lddw, cTMap * _in, int r, int c, bool minmod = false);
    double UF2D_Derivative2(cTMap * _dem, cTMap * _in, int r, int c, int direction);
    double UF1D_Derivative2(cTMap * _ldd,cTMap * _lddw, cTMap * _in, int r, int c);

    double UF_MinMod(double a, double b);

    bool UF_OUTORMV(cTMap * mask, int r, int c);
    bool UF_NOTIME(cTMap * mask,cTMap * dt, int r, int c);

    void upstream(cTMap *_LDD, cTMap *_M, cTMap *out);

    ////Unified Flow Soil Interactions
    //General Function
    void UnifiedFlowSediment();
    void UF_FlowDetachment(double dt);
    void UF_FlowEntrainment(double dt);
    void UF_FlowDetachment(double dt, int r, int c,int d, bool channel);
    void UF_FlowEntrainment(double dt, int r, int c,int d, bool channel);

    double UF_SoilTake(int r, int c, int d, double potential,bool channel,bool bedload);
    void UF_SoilAdd(int r, int c, int d, double mass, bool channel);

    void UF_SumGrainClasses();

    //transport capacity
    double UnifiedFlowTransportCapacity(int r, int c, int d, bool channel, bool bedload);
    //active entrainment
    double UnifiedFlowActiveEntrainment(double hf, double vf,double hs, double vs, double visc, int d);

    //connection to the digital elevation model
    void UFDEMLDD_Connection(cTMap *  dt,cTMap * RemovedMaterial1D, cTMap * RemovedMaterial2D, cTMap * out_DEM,cTMap * out_LDD);

    //////////////
    //SWATRE
    //////////////
    /// filenames for Swatre soil information
    QString SwatreTableDir;
    QString SwatreTableName;
    QString initheadName;

    double swatreDT;
    bool initSwatreStructure;

    /// SWATRE infiltration model 3D soil structure
    SOIL_MODEL *SwatreSoilModel;
    SOIL_MODEL *SwatreSoilModelCrust;
    SOIL_MODEL *SwatreSoilModelCompact;
    SOIL_MODEL *SwatreSoilModelGrass;
    PROFILE **profileList;
    int nrProfileList, sizeProfileList;
    ZONE *zone;
    double precision;
    int tnode; //VJ 110122 node nr in profile with tile drains

    SOIL_MODEL *InitSwatre(cTMap *profileMap);//, QString initHeadMaps, cTMap *tiledepthMap, double dtMin);
    void SwatreStep(SOIL_MODEL *s, cTMap *_WH, cTMap *_fpot, cTMap *_drain, cTMap *_theta, cTMap *where);
    void CloseSwatre(SOIL_MODEL *s);
    void FreeSwatreInfo(void);
    //VJ 111104 old stuff, no longer used but kept for now
    int ReadSwatreInput(QString fileName, QString tablePath);
    ZONE *ReadNodeDefinition(FILE *f);
    PROFILE *ReadProfileDefinition(FILE *f, ZONE *z, const char *tablePath);
    PROFILE *ProfileNr(int profileNr);
    // VJ 111104 constructing profile with Qt commands
    QStringList swatreProfileDef;
    QList<int> swatreProfileNr;
    int *profileNr;
    void ReadSwatreInputNew(void);
    ZONE *ReadNodeDefinitionNew(void);
    PROFILE *ReadProfileDefinitionNew(int pos, ZONE *z);

    HORIZON *ReadHorizon(const char *tablePath,	const  char *tableName);
    double *ReadSoilTable(const char *fileName, int *nrRows);
    void ReadCols(const char *fileName, double *inLut, const char *buf, int   n);
    void InitializeProfile(void);
    void HeadCalc(double *h, bool *ponded, const PROFILE *p ,const double  *thetaPrev,
                  const double  *hPrev, const double  *kavg, const double  *dimoca,
                  bool fltsat, double dt, double pond, double qtop, double qbot);
    double  NewTimeStep(double prevDt, const double *hLast, const double *h, int nrNodes,
                        double precParam, double dtMin, double dtMax);
    void ComputeForPixel(PIXEL_INFO *pixel, double *waterHeightIO, double *infil, double *drain,
                         double drainfraction, double *repel, double *Theta, SOIL_MODEL *s);

    void Totals(void);
    void MassBalance(void);
    void OutputUI(void);
    void reportAll(void);
    void ReportTimeseriesNew(void);
    //void ReportTotals(void);
    void ReportMaps(void);
    void ReportTotalsNew(void);
    void ReportLandunits(void); //VJ 110107 report erosion stats per land unit
    void CountLandunits(void); //VJ 110107 report erosion stats per land unit

    double itercount;
    // thread management variables
    bool stopRequested;
    bool waitRequested;
    bool noInterface;
    bool noOutput;
    bool batchmode;
    QMutex mutex;
    QWaitCondition condition;
    void stop();

protected:
    void run();
    QTime time_ms;
    // talk to the interface

    void setupDisplayMaps();
    void setupHydrographData();
    void ClearHydrographData();

    //combobox map selection
    void GetComboMaps();
    void ClearComboMaps();
    void AddComboMap(int listn, QString name, QString unit,cTMap * map,QList<double> ColorMap,
                     QList<QString> Colors, bool log = false,bool symcol = false, double scale = 1.0, double step = 1.0);

signals:
    void done(const QString &results);
    void debug(const QString &results);
    void show(); //use the output structure "op" declared in global.h and LisUIoutput.h

private slots:   //note, was private loop but dixygen does not recognize that
    /// the main model loop, from here all processes are called in a time loop
    void DoModel();

};



#endif
