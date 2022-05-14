/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
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


#include <QtGui>
#include <QMutex>
#include "C:\\Qt\\msys64\\mingw64\\lib\\gcc\\x86_64-w64-mingw32\\10.3.0\\include\\omp.h"
//#include <omp.h>

#include "CsfMap.h"
#include "io.h"
#include "error.h"
#include "swatre_p.h"
#include "swatre_g.h"

//#define OLDSWATRE 1

//---------------------------------------------------------------------------
//#define NULL nullptr

#define PI 3.14159265

#define HMIN 1e-6

#define DEBUG(s) emit debug(QString(s))
#define TIMEDB(s) emit timedb(QString(s))

//#define mwrite(name) writeRaster(QString(resultDir+name))
#define report(raster, name) WriteMapSeries(raster, resultDir,QString(name), printstep, mapFormat)

// defines to make life easier

/// shortcut to access data
#define Drc     data[r][c]
#define Drcr    data[rr][cr]
#define Drcd    at(d)->data[r][c]
#define Drcdr   at(d)->data[rr][cr]
#define Drci   data[r+dy[i]][c+dx[i]]

#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
    ( ldd != 0 &&  rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

/// shortcuts missing value and inside map
#define MV(r,c) pcr::isMV(LDD->data[r][c])
#define notMVIn(r,c) (!pcr::isMV(LDD->data[r][c]) && r < _nrRows && c < _nrCols && r >= 0 && c >= 0)
#define INSIDE(r, c) (r>=0 && r<_nrRows && c>=0 && c<_nrCols)
#define OUTORMV(r, c)  (INSIDE(r,c) && !pcr::isMV(LDD->data[r][c]) ? false : true)

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_MV for(int r = 0; r < _nrRows; r++)\
    for (int c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(LDD->data[r][c]))

//SLOW!
//#define FOR_ROW_COL_MV_L for(int r = 0; r < _nrRows; r++)\
//    for (int c = 0; c < _nrCols; c++)\
//    if(!pcr::isMV(LDD->data[r][c]))

// faster
//#define FOR_ROW_COL_MV_L for(QVector <LDD_COOR>::iterator crit_ = cr_.begin(); crit_ != cr_.end();  ++crit_)\
//{int r = crit_->r; int c = crit_->c;

// faster
//#define FOR_ROW_COL_MV_L for(long i_ = nrValidCells-1; i_ >= 0; i_--)\
//{long _i_ = cri_[i_]; int r = (int)(_i_/_nrCols); int c = (int)(_i_ % _nrCols);

// fastest, QVector stores all elements in the same consequtive memory!
//#define FOR_ROW_COL_MV_L for(long i_ = 0; i_< nrValidCells; i_++)\
//{int r = cr_[i_].r; int c = cr_[i_].c;

#define FOR_ROW_COL_MV_L for(long i_ = nrValidCells-1; i_ >= 0; i_--)\
 {int r = cr_[i_].r; int c = cr_[i_].c;

//#define FOR_ROW_COL_MV_L for(long i_ = 0; i_ < _nrCols*_nrRows; i_++)\
//{int r = i_/_nrCols; int c=i_%_nrCols;\
//if(!pcr::isMV(LDD->data[r][c]))

//nrValidCellsWS
#define FOR_ROW_COL_MV_LWS(nr_) for(long i_ = WScr.at(nr_).size()-1; i_ >= 0; i_--)\
{int r = WScr.at(nr_)[i_].r; int c = WScr.at(nr_)[i_].c;


#define FOR_ROW_COL_LDD5 for(long i_ = nrValidCellsLDD5-1; i_ >= 0; i_--)\
{int r = crldd5_[i_].r; int c = crldd5_[i_].c;

#define FOR_ROW_COL_LDDCH5 for(long i_ = nrValidCellsLDDCH5-1; i_ >= 0; i_--)\
{int r = crlddch5_[i_].r; int c = crlddch5_[i_].c;

#define FOR_ROW_COL_MV_CHL for(long i_ = nrValidCellsCH-1; i_ >= 0; i_--)\
{int r = crch_[i_].r; int c = crch_[i_].c;

#define FOR_ROW_COL_MV_OUTL for(int i_ = 0; i_ < crout_.size(); i_++)\
{int r = crout_[i_].r; int c = crout_[i_].c;


#define FOR_GRAIN_CLASSES for(int d  = 0 ; d < numgrainclasses;d++)

/// shortcut for channel row and col loop
#define FOR_ROW_COL_MV_CH for (int  r = 0; r < _nrRows; r++)\
    for (int  c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(LDDChannel->data[r][c]))

/// shortcut for tile network row and col loop
#define FOR_ROW_COL_MV_TILE for (int  r = 0; r < _nrRows; r++)\
    for (int  c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(LDDTile->data[r][c]))

#define NRUNITS 512  /// \def max number of landunits or depth classes in flooding
#define NUMNAMES 512   /// \def NUMNAMES runfile namelist max
#define NUMMAPS 512    /// \def max nr maps
#define MIN_FLUX 1e-6 /// \def minimum flux (m3/s) in kinematic wave
#define MIN_HEIGHT 1e-6 /// \def minimum water height (m) for transport of sediment
#define MAXCONC 848.0    /// \def max concentration susp. sed. in kg/m3 0.32 * 2650 = max vol conc from experiments Govers x bulk density
#define MAXCONCBL 848.0    /// \def max concentration susp. sed. in kg/m3 0.32 * 2650 = max vol conc from experiments Govers x bulk density
#define MIN_SLOPE 1e-6

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
#define VANALBEDA 3
#define VANLEER 2

#define FSGOVERS 0
#define FSRIJN 1
#define FSRIJNFULL 2
#define FENGELUND 3

#define FSHAIRSINEROSE 11
#define FSWUWANGJIA 10
//#define FSWUWANGJIABL 3

#define RGOVERS 0
#define RRIJN 1
#define RRIJNFULL 2
#define RWUWANGJIA 3

#define K2D_METHOD_KIN   1
#define K2D_METHOD_KINDYN  2
#define K2D_METHOD_DYN   3


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
typedef struct IDI_POINT {
    int r;
    int c;
    double V;
}  IDI_POINT;
//---------------------------------------------------------------------------list
typedef struct LDD_COOR {
    int r;
    int c;
}  LDD_COOR;
//---------------------------------------------------------------------------list
typedef struct LDD_COORi {
    int r;
    int c;
    int ldd;
}  LDD_COORi;
//---------------------------------------------------------------------------
typedef struct LDD_COORIN {
    int r;
    int c;
    int nr;
    int ldd;
    LDD_COOR *inn;
    //QVector <LDD_COOR> in;
}  LDD_COORIN;
//---------------------------------------------------------------------------
typedef struct LDD_COORloc {
    int r;
    int c;
    int loc;
    int nr;
}  LDD_COORloc;
//---------------------------------------------------------------------------
typedef struct LDD_COORout {
    int r;
    int c;
    int nr;
    QString code;
}  LDD_COORout;
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
    bool gotit;
} NAME_LIST;
//---------------------------------------------------------------------------
/// structure for output of land unit stats
typedef struct UNIT_LIST {
    long nr; // can have any class number
    double var0;
    double var1;
    double var2;
    double var3;
    double var4;
    double var5;
    double var6;
    double var7;
} UNIT_LIST;
//---------------------------------------------------------------------------
/// vec4 used for HLL
typedef struct vec4 { double v[4]; } vec4;
//---------------------------------------------------------------------------
/// Structure to store rain station values of rainfile mapnames
typedef struct RAIN_LIST {
    double time;
    QVector <double> intensity;
} RAIN_LIST;
//---------------------------------------------------------------------------
/// Structure to store rain station values of rainfile mapnames
typedef struct METEO_LIST {
    double time;
    QString name;
    double calib;
} METEO_LIST;
//---------------------------------------------------------------------------
/// Structure to store rain station values of rainfile mapnames
typedef struct Q_LIST {
    double time;
    QVector <double> Qin;
} Q_LIST;
//---------------------------------------------------------------------------
typedef struct BUFFER_LIST {
    double area;
    double h;
    int ID;
} BUFFER_LIST;

typedef struct ExtCH {
    QList <int> childRow;
    QList <int> childCol;
    int chRow;
    int chCol;
    bool isExtended;
} ExtCH;


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
    TWorld(QObject *parent = nullptr);
    ~TWorld();

    QLocale loc;

    /// copy of overall rows and columns, set in initmask
    int _nrRows;
    int _nrCols;

    long nrValidCells;
    long nrValidCellsLDD5;
    long nrValidCellsCH;
    long nrValidCellsLDDCH5;
    long nrValidCellsWS;
    int nrWatersheds;
    QVector <LDD_COOR> crldd5_;
    QVector <LDD_COOR> crlddch5_;
    QVector <LDD_COOR> cr_;
    QVector <LDD_COOR> crws_;
    QList< QVector <LDD_COOR> > WScr;
    QVector <LDD_COORi> dcr_;
    QVector <LDD_COOR> crch_;
    QVector <LDD_COORIN> crlinkedldd_;
    QVector <LDD_COORIN> crlinkedlddch_;
    QVector <LDD_COORIN> crlinkedlddbase_;
    QVector <LDD_COORout> crout_;

    /// map management structure, automatic adding and deleting of all cTMap variables
    MapListStruct maplistCTMap[NUMNAMES];
    int maplistnr;

    /// variable declaration list of all maps with comments:
#include "TMmapVariables.h"

    /// SwitchXXX are boolean options that are set in interface and runfile, mainly corrsponding to checkboxes in the UI

/* OBSOLETE
 SwitchLimitDepTC,SwatreInitialized, SwitchLimitTC, SwitchWheelPresent,SwitchChannelExtended,startbaseflowincrease,SwitchAllinChannel,
 SwitchErosionInsideLoop, SwitchSimpleDepression,SwitchRunoffPerM,SwitchWheelAsChannel, SwitchGullies, SwitchGullyEqualWD, SwitchGullyInfil,
 SwitchGullyInit,SwitchMapoutRunoff, SwitchMapoutConc, SwitchWritePCRnames,SwitchNoErosionOutlet,SwitchSimpleSedKinWave, SwitchSOBEKoutput,
 SwitchMapoutWH, SwitchMapoutWHC, SwitchMapoutTC, SwitchMapoutEros, SwitchMapoutDepo, SwitchMapoutV,SwitchMapoutInf, SwitchMapoutSs, SwitchMapoutChvol
 SwitchChannelFlood, SwitchFloodSedimentMethod, SwitchStoninessDET,SwitchRainfallFlood,SwitchLevees,SwitchWatershed,SwitchMinDTfloodMethod,SwitchNeedD90,
 int SwitchFlood1D2DCoupling;SwitchPercolation, SwitchInfilGA2, SwitchCompactPresent, SwitchNutrients, SwitchPestout, SwitchDrainage,
*/

    bool SwitchRoadsystem, SwitchHardsurface, SwitchIncludeChannel, SwitchChannelBaseflow,SwitchChannelInflow, SwitchChannelAdjustCHW,
    SwitchChannelBaseflowStationary, SwitchRainfallSatellite, SwitchIncludeET, SwitchETSatellite, SwitchSnowmelt, SwitchSnowmeltSatellite,
    SwitchRainfall, SwitchEventbased, SwitchInverseDistanceRainfall,
    SwitchDailyET, SwitchChannelInfil,  SwitchErosion, SwitchLinkedList, SwitchSedtrap, SwitchInfilCompact,
    SwitchInfilCrust, SwitchGrassStrip, SwitchImpermeable, SwitchDumphead, SwitchWaterRepellency,
    SwitchMulticlass,  SwitchOutputTimeStep, SwitchOutputTimeUser, SwitchWriteCommaDelimited, SwitchWritePCRtimeplot,
    SwitchSeparateOutput, SwitchEndRun, SwitchInterceptionLAI, SwitchTwoLayer,  SwitchChannelKinWave,
    SwitchPCRoutput, SwitchWriteHeaders, SwitchGeometric, SwitchIncludeTile, SwitchIncludeStormDrains, SwitchKETimebased,
    SwitchHouses, SwitchRaindrum, SwitchLitter, Switchheaderpest, SwitchPesticide,
    SwitchTimeavgV, SwitchCorrectDEM, Switch2DDiagonalFlow, Switch2DDiagonalFlowNew, SwitchSWOFopen, SwitchMUSCL,  SwitchFloodInitial, SwitchFlowBarriers, SwitchBuffers,
    SwitchCulverts, SwitchUserCores, SwitchVariableTimestep,  SwitchHeun,  SwitchImage, SwitchResultDatetime,SwitchOutputTimestamp,
    SwitchChannelKinwaveDt, SwitchChannelKinwaveAvg,SwitchSWOFWatersheds,SwitchGravityToChannel,
    SwitchDumpH,SwitchDumpTheta,SwitchDumpK, SwitchIncludeDiffusion, SwitchIncludeRiverDiffusion, SwitchAdvancedOptions, SwitchFixedAngle;

    int SwitchKinematic2D;
    int SwitchEfficiencyDET; // detachment efficiency
    int SwitchEfficiencyDETCH; // channel detachment efficiency
    int SwitchSplashEQ; //use lisem (1) or eurosem (1)
    int ReportDigitsOut;
    int FlowBoundaryType; // open, closed
    int userCores;
    int SwitchSV; //ettling velocity
    double splashb; // splash strength coef b limburg equtions,

    // flow bloundaries
    QList<int> FBid;
    QList<double> FBHeightN;
    QList<double> FBHeightS;
    QList<double> FBHeightE;
    QList<double> FBHeightW;
    QList<double> FBTimeN;
    QList<double> FBTimeS;
    QList<double> FBTimeE;
    QList<double> FBTimeW;

    // extended channel obsolete
    QVector <ExtCH> ExtChannel;
    ExtCH ExtendCH;

    // multiple options that are set in interface or runfile, see defines above

    /// Interception storage function based on LAI
    int InterceptionLAIType;

    /// infiltration method
    int InfilMethod;

    /// erosion units in output: to/ha; kg/cell; kg/m2
    int ErosionUnits;

    /// discharge units in output: l/s or m3/s
    int QUnits;

    /// type of kinetic energy equation;
    int KEequationType;

    /// parameters in KE equations
    double KEParamater_a1, KEParamater_b1, KEParamater_c1;
    double KEParamater_a2, KEParamater_b2;
    double KEParamater_a3, KEParamater_b3;

    double PBiasCorrection;
    double ETBiasCorrection;
    double SoilETMBcorrection;

    double GW_recharge;
    double GW_flow;
    double GW_slope;
    double GW_lag;
    double GW_bypass;
    double GW_threshold;
    double GW_initlevel;

    double rainfallETa_threshold;

    double totetafac;

    /// calibration parameters
    double gsizeCalibrationD50;
    double gsizeCalibrationD90;
    double SmaxCalibration;
    double ksatCalibration;
    double ksatCalibration2;
    double ksat2Calibration;
    double nCalibration;
    double thetaCalibration;
    double psiCalibration;
    double ChnCalibration;
    double ChnTortuosity;
    double ChKsatCalibration;
    double COHCalibration;
    double COHCHCalibration;
    double UcrCHCalibration;
    double SVCHCalibration;
    double ASCalibration;
    double SplashDelivery;
    double DepositedCohesion;
    double StripN, SedTrapN;
    double StemflowFraction;
    double DirectEfficiency;

    double CanopyOpeness; // VJ 110209 added Aston factor as user input
    double waterRep_a;
    double waterRep_b;
    double waterRep_c;
    double waterRep_d;

    /// totals for mass balance checks and output
    /// Water totals for mass balance and output (in m3)
    double MB, MBeM3, Qtot, Qtot_dt, QTiletot, IntercTot, IntercETaTot, WaterVolTot, WaterVolSoilTileTot, InfilTot, RainTot, SnowTot, theta1tot, theta2tot;
    double SurfStoremm, InfilKWTot,BaseFlowTot,BaseFlowInit, BaseFlowTotmm, Qfloodout,QfloodoutTot;
    double floodBoundaryTot, floodVolTot, floodVolTotInit, floodVolTotMax, floodAreaMax, floodArea, floodBoundarySedTot, ChannelVolTot, ChannelVolTotmm, WHinitVolTot,StormDrainVolTot;
    double IntercHouseTot, IntercHouseTotmm, IntercLitterTot, IntercLitterTotmm;
    double ChannelSedTot, ChannelDepTot, ChannelDetTot, TileVolTot;
    /// Sediment totals for mass balance and output (in kg)
    double MBs, DetTot, DetSplashTot, DetFlowTot, DepTot, SoilLossTot, SoilLossTot_dt, SedTot,
           FloodDetTot, FloodDepTot, FloodSedTot;
    /// Water totals for output in file and UI (in mm), copied to 'op' structure
    double RainTotmm, SnowTotmm, IntercTotmm, IntercETaTotmm, WaterVolTotmm, WaterVolRunoffmm, FloodBoundarymm, InfilTotmm, Qtotmm, RainAvgmm, SnowAvgmm, WaterVolRunoffmm_F;
    double StormDrainTotmm, floodVolTotmm, floodTotmmInit;
    /// peak times (min)
    double RainstartTime, RainpeakTime, SnowpeakTime, QpeakTime, Qpeak, Rainpeak, Snowpeak;
    bool rainStarted;    
    bool ETStarted;
    double ETstartTime;
    double BulkDens;
    double nrCells, CatchmentArea, nrFloodedCells;
    double LitterSmax, ETaTot, ETaTotmm, ETaTotVol, GWlevel;
    double thetai1tot, thetai2tot, thetai1cur, thetai2cur;

    double maxRainaxis;
    double latitude;
    double CulvertWidth, CulvertHeight, CulvertN, CulvertS;
    //double D50CH, D90CH;

    ///pesticides
    double MBp,PestMassApplied, PestLossTotOutlet, PestFluxTotOutlet, PestRunoffSpatial, PestDisMixing, PestSorMixing, PestInfilt, PestStorage, Pestdetach, PestCinfilt,PestCfilmexit;
    double MBpex,PestRunoffSpatialex,PestDisMixingex,PestSorMixingex,PestInfiltex,PestLossTotOutletex;
    int N_SPK;
    double Maxsolubility;
    double MaxVup;

    /// time and dx parameters
    double time, BeginTime, EndTime;
    double _dt, _dx;
    double _dt_user, _dtCHkin;
    long runstep, printstep, printinterval;

    QString mapFormat;

    /// timeseries variables and output strings
    int nrDischargeseries;
    int nrRainfallseries;
    int nrETseries;
    int nrSnowmeltseries;
    int currentRainfallrow;
    int currentETrow;
    int currentSnowmeltrow;
    int rainplace;
    int ETplace;
    int snowmeltplace;
    QVector <RAIN_LIST> RainfallSeries;  // rainfall vector of records
    QVector <RAIN_LIST> ETSeries;
    QVector <RAIN_LIST> SnowmeltSeries;
    QVector <METEO_LIST> RainfallSeriesMaps;  // rainfall vector of records
    bool calibRainfallinFile;
    QVector <METEO_LIST> ETSeriesMaps;  // rainfall vector of records
    QVector <METEO_LIST> SnowmeltSeriesMaps;  // rainfall vector of records
    QVector <Q_LIST> DischargeInSeries;
    QVector <LDD_COORloc> crQin_;
    QVector <int> locationnnrsrec;

    // output formatting for SOBEK flood model input
    QString SOBEKdatestring;
    int SOBEKnrlines;

    //flow barriers table filename
    QString FlowBarriersFileName;

    // file and directory names
    QString resultDir;
    QString inputDir;
    QString dumpDir;
    QString outflowFileName;
    QString totalSeriesFileName;
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
   // QString runoffFractionMapFileName;
    QString channelDischargeMapFileName;

    QString floodLevelFileName;
    QString floodTimeFileName;
    QString floodStatsFileName;
    QString floodMaxQFileName;
    QString floodMaxChanWHFileName;
    QString floodFEWFileName;
    QString floodMaxVFileName;
    QString floodMaxVHFileName;
    QString floodWHmaxFileName;
    QString tileWaterVolfilename;
    QString tileQmaxfilename;
    QString timestamp;

    QString rainFileDir;
    QString rainFileName;
    QString rainSatFileDir;
    QString rainSatFileName;
    QString ETFileName;
    QString ETSatFileDir;
    QString ETSatFileName;
    QString ETFileDir;
    QString snowmeltFileName;
    QString snowmeltSatFileName;
    QString snowmeltFileDir;
    QString snowmeltSatFileDir;
    QString dischargeinFileName;
    QString dischargeinFileDir;
    QString resultFileName;
    QString temprunname;
    /// standard names of output map series
    QString Outrunoff, Outconc, Outwh, Outrwh, Outvelo, Outinf, Outss, Outchvol,
    Outtc, Outeros, Outdepo, OutSL, OutSed, OutInt,OutSedSS, OutSedBL,
    OutTiledrain, OutTileVol,OutTileV, OutHmx, OutVf, OutQf, OutHmxWH, OutTheta1, OutTheta2;
    bool  SwitchOutrunoff, SwitchOutconc, SwitchOutwh, SwitchOutrwh, SwitchOutvelo, SwitchOutinf, SwitchOutss, SwitchOutchvol,
    SwitchOutConc, SwitchOutTC, SwitchOutDet, SwitchOutDep, SwitchOutSL, SwitchOutSed, SwitchOutInt, SwitchOutSedSS, SwitchOutSedBL,
    SwitchOutTiledrain, SwitchOutTileVol, SwitchOutHmx, SwitchOutVf, SwitchOutQf, SwitchOutHmxWH, SwitchOutTheta;
    QString errorFileName;
    QString errorSedFileName;
    QString satImageFileName;
    QString satImageFileDir;

    // list with class values of land unit map
    UNIT_LIST unitList[NRUNITS];
    UNIT_LIST floodList[NRUNITS];
    int landUnitNr;
    // data initialization, runfile reading and parsing
    NAME_LIST runnamelist[NUMNAMES]; // structure for runfile variables and names
    int nrrunnamelist;

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
    void InitBoundary(void); //VJ 110112
    void InitShade(void); //VJ 130301
    void InitImages(void);
    void InitErosion(void); //VJ 110511
    void GetInputData(void);      // get and make input maps
    void InitParameters(void);
    void IntializeData(void);     // make all non-input maps
    void IntializeOptions(void);  // set all options to false etc
    void InitStandardInput(void);
    void InitLULCInput(void);
    void InitSoilInput(void);
    void InitFlood(void);
    void InitScreenChanNetwork();
    void FindChannelAngles();
    void CorrectDEM(cTMap *h, cTMap * g);
    void DiagonalFlowDEM();


    // functions in lisRunfile.cpp
    QString getvaluename(QString vname);
    double getvaluedouble(QString vname);
    int getvalueint(QString vname);
    QString getvaluestring(QString vname);
    QString CheckDir(QString p, bool makeit = false);
    QString GetName(QString p);
    QString checkOutputMapName(QString p, QString S, int i);
    void ParseRunfileData(void);
    void GetRunFile(void);
    //MapListStruct qx[9];

    //FLOOD according to FULLSWOF2D
    double Flood_DTMIN;
    int F_scheme, F_fluxLimiter, F_MaxIter, F_AddGravity;
    double F_Angle, F_minWH;
   // double HLL2_f1, HLL2_f2, HLL2_f3, HLL2_cfl, HLL_tmp;
    double F_pitValue;
    bool prepareFlood, startFlood;
    int iter_n;
    int F_SWOFSolution;

    double fullSWOF2RO(cTMap *h, cTMap *u, cTMap *v, cTMap *z);
    double fullSWOF2open(cTMap *h, cTMap *vx, cTMap *vy, cTMap *z);
    double fullSWOF2openWS(int nr_, cTMap *h, cTMap *vx, cTMap *vy, cTMap *z);
    void ChannelSWOFopen();

    void KinematicSWOFopen(cTMap *_h, cTMap *_V);

    void prepareFloodZ(cTMap *z);
    void setFloodMask(cTMap * h);
    void setFloodMaskDT(cTMap * DT);
    void setFloodDT(cTMap * h);
    double checkforMinMaxV(double Ves1);

    double limiter(double a, double b);
    vec4 F_HLL3(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    vec4 F_HLL2(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    vec4 F_HLL(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    vec4 F_Rusanov(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    vec4 F_Riemann(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    void correctSpuriousVelocities(int r, int c, cTMap *hes, cTMap *ves1, cTMap *ves2);

    //runoff dynamic
    double maincalcfluxOF(cTMap *_h, double dt, double dt_max);
    void maincalcschemeOF(double dt, cTMap *he, cTMap *ve1, cTMap *ve2,cTMap *hes, cTMap *ves1, cTMap *ves2);
    void dynOutflowPoints(void);
    void OverlandFlow2Ddyn(void);
    void SolveDeepWH(void);
    void Boundary2Ddyn();
    void MUSCLOF(cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z);
    void setZeroOF(cTMap *_h, cTMap *_u, cTMap *_v);
    void simpleSchemeOF(cTMap *_h,cTMap *_u,cTMap *_v);
    void SWOFDiagonalFlow(double dt_req_min, cTMap *h, cTMap *vx, cTMap *vy);
    void SWOFDiagonalFlowNew(double dt_req_min, cTMap *h, cTMap *vx, cTMap *vy);

    void infilInWave(cTMap *_h, double dt1);

    void MUSCL(cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z);
    void maincalcscheme(double dt, cTMap *he, cTMap *ve1, cTMap *ve2,cTMap *hes, cTMap *ves1, cTMap *ves2);
    double maincalcflux(cTMap *_h,double dt, double dt_max);

    //SEDIMENT TRANSPORT
    int FS_SS_Method;
    int FS_BL_Method;
    double FS_SigmaDiffusion;

    int R_SS_Method;
    int R_BL_Method;
    double R_SigmaDiffusion;

    //int GrainSizeDistributionType;

    bool SwitchAdvancedSed, SwitchUseMaterialDepth,SwitchNoBoundarySed, SwitchUse2Phase,SwitchUseGrainSizeDistribution,
        SwitchEstimateGrainSizeDistribution,SwitchReadGrainSizeDistribution, SwitchD50CHavg;

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
    //flood sediment
    QList<cTMap *> BL_D; //bed load sediment for a certain grain size (see graindiameters)
    QList<cTMap *> SS_D; //suspended sediment for a certain grain size
    QList<cTMap *> BLC_D; //concentration
    QList<cTMap *> SSC_D; //concentration
    QList<cTMap *> BLTC_D; //transport capacity
    QList<cTMap *> SSTC_D; //transport capacity
    QList<cTMap *> BLD_D; //layer depth
    QList<cTMap *> SSD_D; //layer depth

    //river sediment
    QList<cTMap *> RBL_D;
    QList<cTMap *> RSS_D;
    QList<cTMap *> RBLC_D;
    QList<cTMap *> RSSC_D;
    QList<cTMap *> RBLTC_D;
    QList<cTMap *> RSSTC_D;
    QList<cTMap *> RBLD_D;
    QList<cTMap *> RSSD_D;

    //overland flow
    QList<cTMap *> Sed_D;
    QList<cTMap *> TC_D;
    QList<cTMap *> Conc_D;

    //used for advection in the 1d kinematic method
    QList<cTMap *> Tempa_D;
    QList<cTMap *> Tempb_D;
    QList<cTMap *> Tempc_D;
    QList<cTMap *> Tempd_D;

    //material that is available for detachment
    QList<cTMap *> StorageDep_D;
    QList<cTMap *> Storage_D;
    cTMap *Storage;
    cTMap *StorageDep;
    cTMap *SedimentMixingDepth;
    cTMap *maxDetachment;

    QList<cTMap *> RStorageDep_D;
    QList<cTMap *> RStorage_D;
    cTMap *RStorage;
    cTMap *RStorageDep;
    cTMap *RSedimentMixingDepth;

    cTMap *unity;

    //keep track of any dissolved substances that need to be advected by the kinematic wave
    //not used!!!
    QList<cTMap *> OF_Advect;
    QList<cTMap *> R_Advect;
    QList<cTMap *> F_Advect;


    //sediment for SWOF flood model
   // void SWOFSedimentBalance();
    void SWOFSedimentCheckZero(int r, int c, cTMap * h);
    void SWOFSedimentSetConcentration(int r, int c, cTMap * h);
    void SWOFSedimentDiffusion(double dt, cTMap * h,cTMap * u,cTMap * v, cTMap * _SS,cTMap * _SSC);
    void SWOFSedimentFlowInterpolation(double dt, cTMap * h,cTMap * u,cTMap * v, cTMap * _SS,cTMap * _SSC);
    void SWOFSedimentDet(cTMap *dt,int r,int c, cTMap * h,cTMap * u,cTMap * v);
    void SWOFSedimentDetNew(double dt, cTMap * h,cTMap * u,cTMap * v);
    void SWOFSediment(double dt, cTMap * h,cTMap * u,cTMap * v);
    void SWOFSedimentLayerDepth(int r , int c, double h, double velocity);//cTMap * u,cTMap * v);

    double simpleSedCalc(double Qj1i1, double Qj1i, double Sj1i, double vol, double sed);
    double complexSedCalc(double Qj1i1, double Qj1i, double Qji1, double Sj1i,double Sji1, double alpha, double dx);

    void routeSubstance(int pitRowNr, int pitColNr, cTMap *_LDD,
                                cTMap *_Q, cTMap *_Qn, cTMap *_Qs, cTMap *_Qsn,
                                cTMap *_Alpha, cTMap *_DX, cTMap*_Sed);//,cTMap*_VolStore, cTMap*_SedStore);

  //  double K2DSolvebyInterpolationSed( cTMap *M, cTMap *MC);


    double GetDpMat(int r, int c,double p,QList<cTMap *> *M);
    double GetMpMat(int r, int c,double p,QList<cTMap *> *M, QList<double> *V);
    double GetDp(int r, int c,double p);
    double GetTotalDW(int r, int c,QList<cTMap *> *M);
    double GetSV(double d);
    void SplashDetachment();
    double MaxConcentration(double watvol, double *sedvol, double *dep);
    void ChannelFlowDetachmentNew();


    void RiverSedimentDiffusion(double dt, cTMap * _SS,cTMap * _SSC);
    void RiverSedimentLayerDepth(int r , int c);
     void RiverSedimentMaxC(int r, int c);

    double calcTCSuspended(int r,int c, int _d, int method, double h, double U, int type);
    double calcTCBedload(int r,int c, int _d, int method, double h, double U, int type);


    void FindBaseFlow(); //search for channel inflow from groundwater
    bool addedbaseflow;

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
    void GetInputTimeseries();
    void GetDischargeData(QString name);
    void GetRainfallData(QString name);   // get input timeseries
    void GetSpatialMeteoData(QString name, int type);   // get input timeseries
    void GetETData(QString name);   // get input timeseries
    void GetSnowmeltData(QString name);   // get input timeseries
    double getTimefromString(QString sss);
    double getmaxRainfall();
    void IDInterpolation(QVector <IDI_POINT> *pointlist, double IDIpower);
    /// convert rainfall of a timestep into a map
    void GetRainfallMapfromStations(void);
    void GetRainfallMap(void);
    /// convert ET of a timestep into a map
    void GetETSatMap(void);
    void GetETMap(void);
    /// convert snowmelt of a timestep into a map
    void GetSnowmeltMap(void);
    void DischargeInflow(void);
    /// interception of vegetation canopy resulting in rainnet
    void Interception();
    /// infiltration function calling all infiltration methods
    /// add rainnet to WH and calculating new WH
    void InterceptionLitter();
    /// subtract water retained in Litter under forest e.g.
    void InterceptionHouses();
    /// subtract water retained on houses, for urban projects    
    void addRainfallWH();
    /// add net rainfall to WH, WHroads and WHgrass

    void cell_Interception(int r, int c);
    double cell_Percolation(int r, int c, double factor);
    double cell_Percolation1(int r, int c, double factor);
    void cell_Redistribution(int r, int c);
    void cell_SurfaceStorage(int r, int c);
    void cell_InfilMethods(int r, int c);
    void cell_InfilSwatre(int r, int c);
    void cell_depositInfil(int r, int c);
    void cell_SplashDetachment(int r, int c);
    void cell_FlowDetachment(int r, int c);
    void MoistureContent();

    void InfilEffectiveKsat();
    void Infiltration();
    void InfilSwatre();

    double IncreaseInfiltrationDepthNew(double fact_, int r, int c);

    void SoilWater();
    void InfilMethods(cTMap *_Ksateff, cTMap *_WH, cTMap *_fpot, cTMap *_fact, cTMap *_L1, cTMap *_L2, cTMap *_FFull);
    void SurfaceStorage();
    void doETa();
    void avgTheta();
    void OverlandFlow();
    void OverlandFlow2D();
    void correctWH(cTMap *_WH);

    void HydrologyProcesses();
    void ChannelandTileflow();

    void OverlandFlow1D(void);
    void ChannelFlow();
    void ChannelBaseflow();
    void ChannelRainandInfil();
    void ChannelSedimentFlow();
    void ChannelFlowandErosion();

    double getMassCH(cTMap *M);
    void correctMassBalanceCH(double sum1, cTMap *M);
    void ToChannel();//int r, int c);
    void ToFlood();
    void CalcVelDisch(); //(int r, int c);
    void CalcVelDischChannel();
    void fromChannelVoltoWH(int r, int c);
    double channelVoltoWH(double vol, int r, int c);
    void fromChannelWHtoVol(int r, int c);
    void ToTiledrain();//);
    void TileFlow(void);
    void StormDrainFlow(void);
    void CalcVelDischTile(void);
    void CalcVelDischDrain(void);
    void GridCell();

    void doExtendRow(int r, int c, int n,  double w2, double adx);
    void doExtendCol(int r, int c, int n,  double w2, double adx);
    void extendRow(int r, int c, int i, double w);
    void extendCol(int r, int c, int i, double w);

    void ExtendChannel();
    bool ExtendChannelNew();
    bool IsExtendedChannel(int r, int c, int dr, int dc);
    void DistributeOverExtendedChannel(cTMap * _In, cTMap * _Out);

    void InitFlowBarriers(void);
    double DEMFB(int r, int c, int rd, int cd, bool addwh);
    double FB(int r, int c, int rd, int cd);
    void SetFlowBarriers();
    void GetFlowBarrierData(QString name);
    double FBW(double h, int r, int c, int dr, int dc);

    long nrGridcells;
    long nrFloodcells;
    void ChannelFlood(void);
    void FloodMaxandTiming(cTMap *_h, cTMap *_UV, double threshold);
    void ChannelFloodStatistics(void);
    void ChannelOverflow(cTMap *_h, cTMap *_V);
    void ChannelOverflowNew(cTMap *_h, cTMap *_V, bool doOF);

    double courant_factor;
    double courant_factorSed;
    double mixing_coefficient, runoff_partitioning;
    double minReportFloodHeight;
    void correctMassBalance(double sum1, cTMap *M, double th);
    void correctMassBalanceWS(int nr_, double sum1, cTMap *M, double th);
    void correctMassBalanceSed(double sum1, cTMap *M, double th);
    double getMass(cTMap *M, double th);
    double getMassWS(int nr_, cTMap *M, double th);
    double getMassSed(cTMap *M, double th);
    void Kinematic(int pitRowNr, int pitColNr, cTMap *_LDD, cTMap *_Q, cTMap *_Qn, cTMap *_q, cTMap *_Alpha, cTMap *_DX, cTMap *_Vol);
    double IterateToQnew(double Qin, double Qold, double q, double alpha, double deltaT, double deltaX, double maxQ);
    void upstream(cTMap *_LDD, cTMap *_M, cTMap *out);
    void KinematicExplicit(QVector<LDD_COORIN> _crlinked, cTMap *_Q, cTMap *_Qn, cTMap *_q, cTMap *_Alpha,cTMap *_DX, cTMap *_Qmax);
    void KinematicSubstance(QVector<LDD_COORIN> _crlinked_, cTMap *_LDD, cTMap *_Q, cTMap *_Qn, cTMap *_Qs, cTMap *_Qsn, cTMap *_Alpha,cTMap *_DX, cTMap *_Sed);
    void AccufluxGW(QVector <LDD_COORIN>_crlinked_ , cTMap *_Q, cTMap *_Qn, cTMap *_CW);
    //LDD_COOR *_crlinked_
    QVector <LDD_COORIN> MakeLinkedList(cTMap *_LDD);

    // kinematic 2D
    double BoundaryQ;
    double BoundaryQs;
    double TimestepfloodMin, TimestepfloodLast;
    //SWATRE
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
    void SwatreStep(int step, int r, int c, SOIL_MODEL *s, cTMap *_WH, cTMap *_fpot, cTMap *_drain, cTMap *_theta);//, cTMap *where);
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


    double MapTotal(cTMap &M);
    void Totals(void);
    void MassBalance(void);
    void OutputUI(void);
    void reportAll(void);
    void ReportTimeseriesNew(void);
    void ReportTotalSeries(void);

    void ReportMaps(void);
    void ReportDump(void);
    void ReportMapSeries(void);
    void ReportTotalsNew(void);
    void ReportLandunits(void); //VJ 110107 report erosion stats per land unit
    void CountLandunits(void); //VJ 110107 report erosion stats per land unit
    void saveMBerror2file(bool doError, bool start);

    int nrBuffers;
    QVector <BUFFER_LIST> bufferarea;
    double itercount;
    // thread management variables
    bool stopRequested;
    bool waitRequested;
    bool noInterface;
    bool noInfo;
    bool noOutput;
    bool batchmode;
    QMutex mutex;
    QWaitCondition condition;
    void stop();

//    QList<double> TSList_point;
//    QList<double> TSList_rainav;
//    QList<double> TSList_snowav;
//    QList<double> TSList_q;
//    QList<double> TSList_h;
//    QList<double> TSList_qs;
//    QList<double> TSList_c;

protected:
    void run();
  //  QTime time_ms;
    // talk to the interface
    QElapsedTimer time_ms;
    double startTime;
    void setupDisplayMaps();
    void setupHydrographData();
    void ClearHydrographData();

    //combobox map selection
    void GetComboMaps();
    void ClearComboMaps();
    void AddComboMap(int listn, QString name, QString unit,cTMap * map,QList<double> ColorMap,
                     QList<QString> Colors, bool log = false,bool symcol = false, double scale = 1.0, double step = 1.0);
    void setLegendColors();

    QList<double> Colormap;
    QList<QString> Colors;

    QList <QStringList> Legend;
    QList <QList <double>> LegendMap;

signals:
    void done(const QString &results);
    void debug(const QString &results);
    void timedb(const QString &results);
    void show(bool showall); //use the output structure "op" declared in global.h and LisUIoutput.h

private slots:   //note, was private loop but dixygen does not recognize that
    /// the main model loop, from here all processes are called in a time loop
    void DoModel();

};



#endif
