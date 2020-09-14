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


#include <QtGui>
#include <QMutex>
#include "C:/Qt/msys64/mingw64/lib/gcc/x86_64-w64-mingw32/10.1.0/include/omp.h"
//#include "omp.h"

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

//#define mwrite(name) writeRaster(QString(resultDir+name))
#define report(raster, name) WriteMapSeries(raster, resultDir,QString(name), printstep, mapFormat)

// defines to make life easier

/// shortcut to access data
#define Drc     data[r][c]
#define Drcr    data[rr][cr]
#define Drcd    at(d)->data[r][c]
#define Drcdr   at(d)->data[rr][cr]
#define Drci    data[r+dr[i]][c+dc[i]]
#define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
    ( ldd != 0 && rFrom >= 0 && cFrom >= 0 && rFrom+dy[ldd]==rTo && cFrom+dx[ldd]==cTo )

/// shortcut missing value in map
#define MV(r,c) pcr::isMV(LDD->data[r][c])
#define MVin(r,c) (pcr::isMV(LDD->data[r][c]) && r < _nrRows && c < _nrCols && r >= 0 && c >= 0)

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_MV for(int r = 0; r < _nrRows; r++)\
    for (int c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(LDD->data[r][c]))

#define FOR_ROW_COL_MV_L for(int r = 0; r < _nrRows; r++)\
    for (int c = 0; c < _nrCols; c++)\
    if(!pcr::isMV(LDD->data[r][c]))

#define FOR_GRAIN_CLASSES for(int d  = 0 ; d < numgrainclasses;d++)

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
#define OUTORMV(r, c)  (INSIDE(r,c) && !pcr::isMV(LDD->data[r][c]) ? false : true)

#define NUMNAMES 2000   /// \def NUMNAMES runfile namelist max
#define NUMMAPS 1000    /// \def max nr maps
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
#define FSHAIRSINEROSE 3
#define FSWUWANGJIA 4
//#define FSWUWANGJIABL 3

#define RGOVERS 0
#define RRIJN 1
#define RRIJNFULL 2
#define RWUWANGJIA 3

#define K2D_METHOD_KIN   1
#define K2D_METHOD_DIFF  2
#define K2D_METHOD_DYN   3
#define K2D_METHOD_KINDYN  4


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
typedef struct LDD_COOR {
    int r;
    int c;
}  LDD_COOR;
//---------------------------------------------------------------------------
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
    long nr;
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
/// Strunture to store rain station values of rainfile mapnames
typedef struct RAIN_LIST {
    double time;
    //double *intensity;
    QVector <double> intensity;
    bool isMap;
    QString name;
} RAIN_LIST;
//---------------------------------------------------------------------------

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
    SwitchWheelPresent, SwitchCompactPresent, SwitchIncludeChannel, SwitchChannelBaseflow, SwitchChannelExtended,
    startbaseflowincrease, SwitchChannelInfil, SwitchAllinChannel, SwitchErosion, SwitchAltErosion,
    SwitchSimpleDepression,SwitchSedtrap, SwitchSnowmelt, SwitchRainfall, /*SwitchRunoffPerM,*/ SwitchInfilCompact,
    SwitchInfilCrust, SwitchGrassStrip, SwitchImpermeable, SwitchPercolation, SwitchDumphead, SwitchWaterRepellency,
    SwitchWheelAsChannel, SwitchMulticlass, SwitchNutrients, SwitchGullies, SwitchGullyEqualWD, SwitchGullyInfil,
    SwitchGullyInit, SwitchOutputTimeStep, SwitchOutputTimeUser, SwitchMapoutRunoff, SwitchMapoutConc,
    SwitchMapoutWH, SwitchMapoutWHC, SwitchMapoutTC, SwitchMapoutEros, SwitchMapoutDepo, SwitchMapoutV,
    SwitchMapoutInf, SwitchMapoutSs, SwitchMapoutChvol, SwitchWritePCRnames, SwitchWriteCommaDelimited, SwitchWritePCRtimeplot,
    SwitchNoErosionOutlet, SwitchDrainage, SwitchPestout, SwitchSeparateOutput, SwitchEndRun,
    SwitchInterceptionLAI, SwitchTwoLayer, SwitchSimpleSedKinWave, SwitchSOBEKoutput, SwitchChannelKinWave,
    SwitchPCRoutput, SwitchWriteHeaders, SwitchGeometric, SwitchIncludeTile, SwitchIncludeStormDrains, SwitchKETimebased, SwitchHouses, SwitchChannelFlood, SwitchRaindrum, SwitchLitter,
    Switchheaderpest, SwitchPesticide, SwitchRainfallFlood, SwitchFloodSedimentMethod, SwitchStoninessDET,
    SwitchTimeavgV, SwitchSWOFopen, SwitchMUSCL, SwitchLevees, SwitchFloodInitial, SwitchWatershed,SwitchFlowBarriers, SwitchBuffers,
    SwitchCulverts, SwitchUserCores, SwitchVariableTimestep, SwitchMinDTfloodMethod, SwitchHeun, SwitchNeedD90, SwitchImage, SwitchResultDatetime,
    SwitchDumpH,SwitchDumpTheta,SwitchDumpK, SwitchIncludeDiffusion, SwitchIncludeRiverDiffusion, SwitchAdvancedOptions, SwitchFixedAngle;
    int SwitchFlood1D2DCoupling;
    int SwitchKinematic2D;
    int SwitchEfficiencyDET;
    int ReportDigitsOut;
    int FlowBoundaryType;
    int userCores;
    int SwitchSV;

    QList<int> FBid;

    QList<double> FBHeightN;
    QList<double> FBHeightS;
    QList<double> FBHeightE;
    QList<double> FBHeightW;

    QList<double> FBTimeN;
    QList<double> FBTimeS;
    QList<double> FBTimeE;
    QList<double> FBTimeW;

    QList <ExtCH> ExtChannel;
    ExtCH ExtendCH;


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
double sum1;
    /// calibration parameters
    double gsizeCalibration;
    double ksatCalibration;
    double nCalibration;
    double thetaCalibration;
    double psiCalibration;
    double ChnCalibration;
    double ChKsatCalibration;
    double COHCalibration;
    double COHCHCalibration;
    double ASCalibration;
    double SplashDelivery;
    double DepositedCohesion;
    double StripN, SedTrapN;
    double StemflowFraction;

    double CanopyOpeness; // VJ 110209 added Aston factor as user input
    double waterRep_a;
    double waterRep_b;
    double waterRep_c;
    double waterRep_d;

    /// totals for mass balance checks and output
    /// Water totals for mass balance and output (in m3)
    double MB, MBeM3, Qtot,QtotT,QTiletot, IntercTot, WaterVolTot, WaterVolSoilTot, InfilTot, RainTot, SnowTot, SurfStoremm, InfilKWTot,BaseFlowTot,Qfloodout,QfloodoutTot;
    double floodBoundaryTot, floodVolTot, floodVolTotInit, floodVolTotMax, floodAreaMax, floodBoundarySedTot, ChannelVolTot, ChannelVolTotmm, WHinitVolTot,StormDrainVolTot;
    double IntercHouseTot, IntercHouseTotmm, IntercLitterTot, IntercLitterTotmm;
    double ChannelSedTot, ChannelDepTot, ChannelDetTot, TileVolTot;
    /// Sediment totals for mass balance and output (in kg)
    double MBs, DetTot, DetSplashTot, DetFlowTot, DepTot, SoilLossTot,SoilLossTotT, SedTot,
           FloodDetTot, FloodDepTot, FloodSedTot;
    /// Water totals for output in file and UI (in mm), copied to 'op' structure
    double RainTotmm, SnowTotmm, IntercTotmm, WaterVolTotmm,WaterVolRunoffmm, FloodBoundarymm, InfilTotmm, Qtotmm, RainAvgmm, SnowAvgmm, WaterVolRunoffmm_F;
    double StormDrainTotmm, floodVolTotmm, floodTotmmInit;
    /// peak times (min)
    double RainstartTime, RainpeakTime, SnowpeakTime, QpeakTime, Qpeak, Rainpeak, Snowpeak;
    bool rainStarted;
    double BulkDens;
    double nrCells, CatchmentArea, nrFloodedCells;
    double LitterSmax, ETaTot, ETaTotmm;

    ///pesticides
    double MBp,PestMassApplied, PestLossTotOutlet, PestFluxTotOutlet, PestRunoffSpatial, PestDisMixing, PestSorMixing, PestInfilt, PestStorage, Pestdetach, PestCinfilt,PestCfilmexit;
    double MBpex,PestRunoffSpatialex,PestDisMixingex,PestSorMixingex,PestInfiltex,PestLossTotOutletex;
    int N_SPK;
    double Maxsolubility;
    double MaxVup;

    /// time and dx parameters
    double time, BeginTime, EndTime;
    double _dt, _dx;
    long runstep, printstep, printinterval;

    QString mapFormat;

    QString timestampRun; // for output filenames

    /// timeseries variables and output strings
    int nrRainfallseries;
    int nrSnowmeltseries;
    QVector <RAIN_LIST> RainfallSeriesM;  // rainfall vector of records
    QVector <RAIN_LIST> SnowmeltSeriesM;

    // output formatting for SOBEK flood model input
    QString SOBEKdatestring;
    int SOBEKnrlines;

    //flow barriers table filename
    QString FlowBarriersFileName;

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

    QString rainFileName;
    QString rainFileDir;
    QString snowmeltFileName;
    QString snowmeltFileDir;
    QString resultFileName;
    QString temprunname;
    /// standard names of output map series
    QString Outrunoff, Outconc, Outwh, Outrwh, Outvelo, Outinf, Outss, Outchvol,
    Outtc, Outeros, Outdepo, OutSL, OutSed, OutInt,OutSedSS, OutSedBL,
    OutTiledrain, OutTileVol,OutTileV, OutHmx, OutVf, OutQf, OutHmxWH;
    bool  SwitchOutrunoff, SwitchOutconc, SwitchOutwh, SwitchOutrwh, SwitchOutvelo, SwitchOutinf, SwitchOutss, SwitchOutchvol,
    SwitchOutConc, SwitchOutTC, SwitchOutDet, SwitchOutDep, SwitchOutSL, SwitchOutSed, SwitchOutInt, SwitchOutSedSS, SwitchOutSedBL,
    SwitchOutTiledrain, SwitchOutTileVol, SwitchOutHmx, SwitchOutVf, SwitchOutQf, SwitchOutHmxWH;
    QString errorFileName;
    QString errorSedFileName;
    QString satImageFileName;
    QString satImageFileDir;

    // list with class values of land unit map
    UNIT_LIST unitList[512]; // just a fixed number for 512 classes, who cares!
    UNIT_LIST floodList[512]; // just a fixed number for 512 classes, who cares!
    QVector <UNIT_LIST> unitListM;
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
    void IntializeData(void);     // make all non-input maps
    void IntializeOptions(void);  // set all options to false etc
    void InitStandardInput(void); //VJ 170211
    void InitFlood(void);
    void InitChanNetwork();
    void FindChannelAngles();

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

    //FLOOD according to FULLSWOF2D
    double Flood_DTMIN;
    int F_scheme, F_fluxLimiter, F_MaxIter, F_AddGravity;
    double F_Angle;
   // double HLL2_f1, HLL2_f2, HLL2_f3, HLL2_cfl, HLL_tmp;
    double F_pitValue;
    bool prepareFlood, startFlood;
    int iter_n;
    int F_SWOFSolution;

    double fullSWOF2RO(cTMap *h, cTMap *u, cTMap *v, cTMap *z);
    double fullSWOF2open(cTMap *h, cTMap *vx, cTMap *vy, cTMap *z);
    void ChannelSWOFopen();
    void makeChannelList();

    void prepareFloodZ(cTMap *z);
    void setFloodMask(cTMap * h);
    void setFloodMaskDT(cTMap * DT);
    void setFloodDT(cTMap * h);

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
    void Boundary2Ddyn();
    void MUSCLOF(cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z);
    void setZeroOF(cTMap *_h, cTMap *_u, cTMap *_v);
    void simpleSchemeOF(cTMap *_h,cTMap *_u,cTMap *_v);

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

    bool SwitchAdvancedSed, SwitchUseMaterialDepth,SwitchNoBoundarySed, SwitchUse2Layer,SwitchUseGrainSizeDistribution, SwitchEstimateGrainSizeDistribution,SwitchReadGrainSizeDistribution;

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
    void SWOFSedimentBalance();
    void SWOFSedimentCheckZero(int r, int c, cTMap * h);
    void SWOFSedimentSetConcentration(int r, int c, cTMap * h);
    void SWOFSedimentDiffusion(cTMap* DT, cTMap * h,cTMap * u,cTMap * v, cTMap * _SS,cTMap * _SSC);
    void SWOFSedimentFlowInterpolation(cTMap* dt, cTMap * h,cTMap * u,cTMap * v, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC);
    void SWOFSedimentDet(cTMap *dt,int r,int c, cTMap * h,cTMap * u,cTMap * v);
    void SWOFSediment(cTMap* DT, cTMap * h,cTMap * u,cTMap * v);
    void SWOFSedimentLayerDepth(int r , int c, double h, double velocity);//cTMap * u,cTMap * v);

    double simpleSedCalc(double Qj1i1, double Qj1i, double Sj1i, double vol, double sed);
    double complexSedCalc(double Qj1i1, double Qj1i, double Qji1, double Sj1i,double Sji1, double alpha, double dx);

    void routeSubstance(int pitRowNr, int pitColNr, cTMap *_LDD,
                                cTMap *_Q, cTMap *_Qn, cTMap *_Qs, cTMap *_Qsn,
                                cTMap *_Alpha, cTMap *_DX, cTMap*_Vol, cTMap*_Sed);//,cTMap*_VolStore, cTMap*_SedStore);

    double K2DSolvebyInterpolationSed( cTMap *M, cTMap *MC);


    double GetDpMat(int r, int c,double p,QList<cTMap *> *M);
    double GetMpMat(int r, int c,double p,QList<cTMap *> *M, QList<double> *V);
    double GetDp(int r, int c,double p);
    double GetTotalDW(int r, int c,QList<cTMap *> *M);
    double GetSV(double d);
    void SplashDetachment();
    void FlowDetachment();
    double MaxConcentration(double watvol, double *sedvol, double *dep);
    void ChannelFlowDetachment();

    void RiverSedimentDiffusion(double dt, cTMap * _SS,cTMap * _SSC);
    void RiverSedimentLayerDepth(int r , int c);
     void RiverSedimentMaxC(int r, int c);

    double calcTCSuspended(int r,int c, int _d, int method, double h, double U, int type);
    double calcTCBedload(int r,int c, int _d, int method, double h, double U, int type);


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
    /// convert rainfall of a timestep into a map
    void RainfallMap(void);
    /// convert snowmelt of a timestep into a map
    void SnowmeltMap(void);
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

    void InfilEffectiveKsat();
    void Infiltration();
    void InfilSwatre(cTMap *_WH);

    double IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p, REAL8 *FFull);
    double IncreaseInfiltrationDepthNew(int r, int c);

    void SoilWater();
    void InfilMethods(cTMap *_Ksateff, cTMap *_WH, cTMap *_fpot, cTMap *_fact, cTMap *_L1, cTMap *_L2, cTMap *_FFull);
    void InfilMethodsNew();
    void SurfaceStorage();
    void doETa();
    void OverlandFlow();
    void OverlandFlow2D();
    void correctWH(cTMap *_WH);

    void OverlandFlow1D(void);

    void ChannelFlow();
    double ChannelIterateWH(double _h, int r, int c);
    void ChannelAddBaseandRain();
    void ChannelWaterHeightFromVolume();
    void ToChannel();//);
    void ToFlood();//);
    void CalcVelDisch();
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
    bool OUTORMVc(int r, int c);  // replaced by define

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
    void ChannelOverflow(cTMap *_h, cTMap *_V, bool doOF);
    void ChannelOverflowNew(cTMap *_h, cTMap *_V, bool doOF);

    double courant_factor;
    double courant_factor_diffusive;
    double courant_factor_sed;
    double mixing_coefficient, runoff_partitioning;
    double minReportFloodHeight;
    void correctMassBalance(double sum1, cTMap *M, double th);
    double getMass(cTMap *M, double th);
    void Kinematic(int pitRowNr, int pitColNr, cTMap *_LDD,
                   cTMap *_Q, cTMap *_Qn,
                   cTMap *_q, cTMap *_Alpha, cTMap *_DX,
                   cTMap *_Vol);
    double IterateToQnew(double Qin, double Qold, double q, double alpha, double deltaT, double deltaX, double maxQ);
    void upstream(cTMap *_LDD, cTMap *_M, cTMap *out);
    void KinematicExplicit(cTMap *_LDD, cTMap *Q, cTMap *Qn, cTMap *q, cTMap *Alpha,cTMap *_DX, cTMap *Vol);
    // kinematic 2D
    double K2DFlux(double t, double tmax);
    void K2DPreSolve();
    void K2DSolve();
    void K2DSolvebyInterpolation();
    void K2DInit();
    void K2DCalcVelDisch();
    void K2DDEMA();
    void K2DDEMAInitial();
    double K2DQOut, K2DQOutBoun;
    double K2DQSOut, K2DQSOutBoun;
    double K2DQPOut;
    double ConcentrateKin;
    double TimestepKinMin;
    double KinematicBoundaryFraction = 0.05;
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
    void SwatreStep(int step, SOIL_MODEL *s, cTMap *_WH, cTMap *_fpot, cTMap *_drain, cTMap *_theta, cTMap *where);
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

    void ReportMaps(void);
    void ReportMapSeries(void);
    void ReportTotalsNew(void);
    void ReportLandunits(void); //VJ 110107 report erosion stats per land unit
    void CountLandunits(void); //VJ 110107 report erosion stats per land unit
    void saveMBerror2file(bool doError, bool start);


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

    QList<double> TSList_point;
    QList<double> TSList_rainav;
    QList<double> TSList_snowav;
    QList<double> TSList_q;
    QList<double> TSList_h;
    QList<double> TSList_qs;
    QList<double> TSList_c;


    void CellProcesses();
    void CellProcesses2();
    void OrderedProcesses();

protected:
    void run();
  //  QTime time_ms;
    // talk to the interface
    QElapsedTimer time_ms;
    void setupDisplayMaps();
    void setupHydrographData();
    void ClearHydrographData();

    //combobox map selection
    void GetComboMaps();
    void setColor(int i);
    void ClearComboMaps();
    void AddComboMap(int listn, QString name, QString unit,cTMap * map,QList<double> ColorMap,
                     QList<QString> Colors, bool log = false,bool symcol = false, double scale = 1.0, double step = 1.0);

    QList<double> Colormap;
    QList<QString> Colors;
signals:
    void done(const QString &results);
    void debug(const QString &results);
    void show(); //use the output structure "op" declared in global.h and LisUIoutput.h

private slots:   //note, was private loop but dixygen does not recognize that
    /// the main model loop, from here all processes are called in a time loop
    void DoModel();

};



#endif
