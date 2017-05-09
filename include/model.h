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
#include "lisUnifiedFlowThread.h"
#include "lisUnifiedFlowThreadPool.h"


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

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_2DMT for(int rc = 0; rc < _nrRows; rc++)\
    {for (int cc = 0; cc < _nrCols; cc++)\
    {int r = (int) (ThreadPool->CellRDerListOrdered2d.at(thread)->data[rc][cc]);\
    int c = (int) (ThreadPool->CellCDerListOrdered2d.at(thread)->data[rc][cc]);\
    if(!INSIDE(r,c)){break;}\


/// shortcut for LDD row and col loop
#define FOR_ROW_COL_UF2DMT for(int rc = 0; rc < _nrRows; rc++)\
    {for (int cc = 0; cc < _nrCols; cc++)\
    {int r = (int) (ThreadPool->CellRMaskListOrdered2d.at(thread)->data[rc][cc]);\
    int c = (int) (ThreadPool->CellCMaskListOrdered2d.at(thread)->data[rc][cc]);\
    if(!INSIDE(r,c)){break;}\

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_UF2DMTDER for(int rc = 0; rc < _nrRows; rc++)\
    {for (int cc = 0; cc < _nrCols; cc++)\
    {int r = (int) (ThreadPool->CellRDerListOrdered2d.at(thread)->data[rc][cc]);\
    int c = (int) (ThreadPool->CellCDerListOrdered2d.at(thread)->data[rc][cc]);\
    if(!INSIDE(r,c)){break;}\

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_UF2DMT_DT for(int rc = 0; rc < _nrRows; rc++)\
    {for (int cc = 0; cc < _nrCols; cc++)\
    {int r = (int) (ThreadPool->CellRListOrdered2d.at(thread)->data[rc][cc]);\
    int c = (int) (ThreadPool->CellCListOrdered2d.at(thread)->data[rc][cc]);\
    if(!INSIDE(r,c)){break;}\


/// shortcut for LDD row and col loop
#define FOR_ROW_COL_UF1DMT  for(int rc = 0; rc < _nrRows; rc++)\
    {for (int cc = 0; cc < _nrCols; cc++)\
    {int r = (int) (ThreadPool->CellRMaskListOrdered1d.at(thread)->data[rc][cc]);\
    int c = (int) (ThreadPool->CellCMaskListOrdered1d.at(thread)->data[rc][cc]);\
    if(!INSIDE(r,c)){break;}\

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_UF1DMTDER  for(int rc = 0; rc < _nrRows; rc++)\
    {for (int cc = 0; cc < _nrCols; cc++)\
    {int r = (int) (ThreadPool->CellRDerListOrdered1d.at(thread)->data[rc][cc]);\
    int c = (int) (ThreadPool->CellCDerListOrdered1d.at(thread)->data[rc][cc]);\
    if(!INSIDE(r,c)){break;}\

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_UF1DMT_DT for(int rc = 0; rc < _nrRows; rc++)\
    {for (int cc = 0; cc < _nrCols; cc++)\
    {int r = (int) (ThreadPool->CellRListOrdered1d.at(thread)->data[rc][cc]);\
    int c = (int) (ThreadPool->CellCListOrdered1d.at(thread)->data[rc][cc]);\
    if(!INSIDE(r,c)){break;}\

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

#define NUMNAMES 3000   /// \def NUMNAMES runfile namelist max
#define NUMMAPS 1000    /// \def max nr maps
#define MIN_FLUX 1e-12 /// \def minimum flux (m3/s) in kinematic wave
#define MIN_HEIGHT 1e-6 /// \def minimum water height (m) for transport of sediment
#define MAXCONC 848.0    /// \def max concentration susp. sed. in kg/m3 0.32 * 2650 = max vol conc from experiments Govers x bulk density
#define MAXCONCBL 848.0    /// \def max concentration susp. sed. in kg/m3 0.32 * 2650 = max vol conc from experiments Govers x bulk density
#define UF_VERY_SMALL 1e-12 /// \def min timestep/flux/height in unified flow equations

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
typedef struct lCOORD {
    int _r;
    int _c;
} lCOORD;
//--------------------------------------------------------------------------
/// structure for watershed coordinates for flooding
typedef struct WS_LIST {
    bool flood;
   QList <lCOORD> cr;
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

    /// include the variables and member fuctions related to unified flow equations
#include "UFModel.h"

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
    SwitchPCRoutput, SwitchWriteHeaders, SwitchGeometric, SwitchIncludeTile, SwitchKETimebased, SwitchHouses, SwitchChannelFlood, SwitchRaindrum,SwitchLitter,
    SwitchRainfallFlood, SwitchFloodSedimentMethod, SwitchStoninessDET, SwitchLevees, SwitchFlowBarriers, SwitchBarriers, SwitchMaxVolume, SwitchChannelMaxVolume, SwitchUFInitial,SwitchUFForced;

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

    double st_scCalibration;
    double st_sifaCalibration;
    double st_sdCalibration;
    bool st_csdCalibration;
    double st_csdsfCalibration;

    double CanopyOpeness; // VJ 110209 added Aston factor as user input
    double waterRep_a;
    double waterRep_b;
    double waterRep_c;
    double waterRep_d;

    double Calibrate_EP;
    double Calibrate_TC;
    double Calibrate_SV;
    double Calibrate_YS;
    double Calibrate_DV;
    double Calibrate_DF;
    double Calibrate_SPF;
    double Calibrate_DC;
    double Calibrate_ESC;
    double Calibrate_EGC;

    bool SF_Calibrate_LF;
    bool SF_Calibrate_Initial;
    double SF_Calibrate_Margin;
    bool SF_Calibrate_First;

    ///rainfall to flood max gradient
 //   double rainFloodingGradient;
    bool firstssreport,sfset;

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
    double nrCells, CatchmentArea, nrFloodedCells;
    double LitterSmax;

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

    QString FileName_MaxDebrisFlowHeight;
    QString FileName_MaxDebrisFlowVelocity;
    QString FileName_DebrisFlowStart;
    QString FileName_Entrainment;
    QString FileName_SlopeFailure;
    QString FileName_MinimumSafetyFactor;

    QString FinalFluidPhaseFileName;
    QString FinalSolidPhaseFileName;
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
    OutTiledrain, OutHmx, OutVf, OutQf, OutHmxWH, OutSafetyFactor,OutSlopeFailure,OutDFH,OutDFV,OutFPH, OutSPH,OutEntrainment, OutTimestep;
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

    bool SwitchUseMaterialDepth,SwitchUse2Layer,SwitchUseGrainSizeDistribution, SwitchEstimateGrainSizeDistribution,SwitchReadGrainSizeDistribution, SwitchSolidPhase,SwitchEntrainment, SwitchSlopeStability, SwitchSlopeFailure;
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

    double SFMIN;
    double SFMAX;
    double SFDISPLAYMAX;
    double SFFAILMIN;

    bool OUTORMV(int r, int c);
    void SlopeStability();
    void SlopeFailure();
    void SafetyFactor();

    void CalculateSafetyFactors(cTMap * _DEM,cTMap * _SoilDepth,
                                       cTMap * _OverlandWater, cTMap * _SoilCohesion,
                                       cTMap * _InternalFrictionAngle,cTMap * _SoilWaterHeight,
                                       cTMap * _SoilWaterSuction, cTMap * _SoilDensity,
                                       cTMap * _PlantCohesion,cTMap * _PlantPressure,
                                       cTMap * _SafetyFactor,cTMap * _Threshold,cTMap * _Threshold1,cTMap * _InititationHeight,cTMap * _Initiated,cTMap * _SFCalibration);
    double CalculateSafetyFactorAt(int r, int c);

    double CalculateSafetyFactorAt(int r, int c, double slope, cTMap * _SoilDepth,
                                       cTMap * _OverlandWater, cTMap * _SoilCohesion,
                                       cTMap * _InternalFrictionAngle,cTMap * _SoilWaterHeight,
                                       cTMap * _SoilWaterSuction, cTMap * _SoilDensity,
                                       cTMap * _PlantCohesion,cTMap * _PlantPressure,
                                       cTMap * _SFCalibration);

    double CalculateSafetyFactor( double slope, double _SoilDepth,
                                       double _OverlandWater, double _SoilCohesion,
                                       double _InternalFrictionAngle,double _SoilWaterHeight,
                                       double _SoilWaterSuction, double _SoilDensity,
                                       double _PlantCohesion, double _PlantPressure,
                                       double _SFCalibration);

    double SolveStableDepthAt(int r, int c);

    double SolveStableDepthAt(int r, int c, double slope, cTMap * _SoilDepth,
                                       cTMap * _OverlandWater, cTMap * _SoilCohesion,
                                       cTMap * _InternalFrictionAngle,cTMap * _SoilWaterHeight,
                                       cTMap * _SoilWaterSuction, cTMap * _SoilDensity,
                                       cTMap * _PlantCohesion,cTMap * _PlantPressure,
                                       cTMap * _Threshold,cTMap * _Threshold1,cTMap * _SFCalibration);

    double SolveStableDepth( double slope, double _SoilDepth,
                                       double _OverlandWater, double _SoilCohesion,
                                       double _InternalFrictionAngle,double _SoilWaterHeight,
                                       double _SoilWaterSuction, double _SoilDensity,
                                       double _PlantCohesion, double _PlantPressure,
                                       double _Threshold, double _Threshold1, double _SFCalibration);

    double GetTotalSoilDepth(int r, int c);

    void InitiateDebrisFlow();

    double GetDpMat(int r, int c,double p,QList<cTMap *> *M);
    double GetMpMat(int r, int c,double p,QList<cTMap *> *M, QList<double> *V);
    double GetDp(int r, int c,double p);
    double GetTotalDW(int r, int c,QList<cTMap *> *M);
    double GetSV(double d);
    void SplashDetachment(int thread);

    double MaxConcentration(double watvol, double sedvol);

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
    void Interception(int thread);
    /// infiltration function calling all infiltration methods
    /// add rainnet to WH and calculating new WH
    void InterceptionLitter(int thread);
    /// subtract water retained in Litter under forest e.g.
    void InterceptionHouses(int thread);
    /// subtract water retained on houses, for urban projects    
    void addRainfallWH(int thread);
    /// add net rainfall to WH, WHroads and WHgrass

    void InfilEffectiveKsat();
    void Infiltration(int thread);

    void InfilSwatre(int thread, cTMap *_WH);
//    void InfilGreenAmpt1(cTMap *_WH);  //OBSOLETE
//    void InfilSmithParlange1(cTMap *_WH); //OBSOLETE
    void InfilMorelSeytoux1(cTMap *_WH);
//    void InfilKsat(cTMap *_WH); //OBSOLETE
    double IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p, REAL8 *FFull);
    void SoilWater(int thread);
    void InfilMethods(int thread,cTMap *_Ksateff, cTMap *_WH, cTMap *_fpot, cTMap *_fact, cTMap *_L1, cTMap *_L2, cTMap *_FFull);

    void SurfaceStorage(int thread);

    void GridCell(int thread);

    void CellProcessWrapper();
    void CellProcesses(int thread);
    void DynamicProcessWrapper();

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
    void ChannelFloodStatistics(int not_used);

    double courant_factor;
    double courant_factor_diffusive;
    double mixing_coefficient, runoff_partitioning;
   // double cfl_fix;
    double minReportFloodHeight;
    double correctMassBalance(double sum1, cTMap *M, double minV);


    ////FlowBarriers
    QList<int> FBid;

    QList<double> FBHeightN;
    QList<double> FBHeightS;
    QList<double> FBHeightE;
    QList<double> FBHeightW;

    QList<double> FBTimeN;
    QList<double> FBTimeS;
    QList<double> FBTimeE;
    QList<double> FBTimeW;


    void InitFlowBarriers(void);
    void SetFlowBarriers();
    void GetFlowBarrierData(QString name);
    double GetFlowBarrierHeight(int r, int c, int dr, int dc);
    //DEM based barriers
    cTMap * Barriers;

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
    void SwatreStep(int thread,SOIL_MODEL *s, cTMap *_WH, cTMap *_fpot, cTMap *_drain, cTMap *_theta, cTMap *where);
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

    QList<double> TSList_point;
    QList<double> TSList_rainav;
    QList<double> TSList_snowav;
    QList<double> TSList_q;
    QList<double> TSList_h;
    QList<double> TSList_qs;
    QList<double> TSList_c;

    void FillTimeSeriesData();

    void Totals(void);
    void MassBalance(void);
    void OutputUI(void);
    void reportAll(void);
    void reportWrapper(int not_used);

    void ReportTimeseriesNew(int not_used);
    //void ReportTotals(void);
    void ReportMaps(int not_used);
    void ReportTotalsNew(int not_used);
    void ReportLandunits(int not_used); //VJ 110107 report erosion stats per land unit
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




    ////MULTITHREADING STUFF
    //LisemThreadPool *ThreadPool;
 std::function<void(int)> fcompute;
 std::function<void(int)> flowcompute;

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
