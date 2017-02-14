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
#define report(raster, name) WriteMapSeries(raster, resultDir,QString(name), printstep, mapFormat)

// defines to make life easier

/// shortcut to access data
#define Drc     data[r][c]
#define Drcr    data[rr][cr]
#define Drcd    at(d)->data[r][c]
#define Drcdr    at(d)->data[rr][cr]
#define Drci data[r+dr[i]][c+dc[i]]

// obsolete /// shortcut to access the outlet point data
//#define DrcOutlet     data[r_outlet][c_outlet]

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

/*
/// shortcut for all cell in watershed with nr wsnr
#define FOR_WATERSHED_ROW_COL(wsnr) for (long k = 0; k < WS[wsnr].cr.count(); k++) {\
    int c = WS[wsnr].cr[k]._c;\
    int r = WS[wsnr].cr[k]._r;\
    */

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


#define NUMNAMES 2000   /// \def NUMNAMES runfile namelist max
#define NUMMAPS 1000    /// \def max nr maps
#define MIN_FLUX 1e-12 /// \def minimum flux (m3/s) in kinematic wave
#define MIN_HEIGHT 1e-12 /// \def minimum water height (m) for transport of sediment
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
    double var6;
    double var7;
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
    SwitchPCRoutput, SwitchWriteHeaders, SwitchGeometric, SwitchIncludeTile, SwitchKETimebased, SwitchHouses, SwitchChannelFlood, SwitchRaindrum, SwitchLitter,
    Switchheaderpest, SwitchPesticide, SwitchRainfallFlood, SwitchFloodSedimentMethod, SwitchStoninessDET,
    SwitchFloodSWOForder1, SwitchFloodSWOForder2, SwitchMUSCL, SwitchLevees, SwitchFloodInitial, SwitchWatershed,SwitchFlowBarriers;

    int SwitchFlood1D2DCoupling;
    int SwitchKinematic2D;
    int SwitchEfficiencyDET;
    int ReportDigitsOut;
    int FlowBoundaryType;

    QList<int> FBid;

    QList<double> FBHeightN;
    QList<double> FBHeightS;
    QList<double> FBHeightE;
    QList<double> FBHeightW;

    QList<double> FBTimeN;
    QList<double> FBTimeS;
    QList<double> FBTimeE;
    QList<double> FBTimeW;


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
    double COHCalibration;
    double ASCalibration;
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
    double MB, MBeM3, Qtot,QtotT,QTiletot, IntercTot, WaterVolTot, WaterVolSoilTot, InfilTot, RainTot, SnowTot, SurfStoremm, InfilKWTot,BaseFlowTot,BaseFlow;
   //obsolete double QtotTileOutlet, QtotOutlet, QtotChannelOutlet;
    double floodBoundaryTot, floodVolTot, floodVolTotInit, floodVolTotMax, floodAreaMax;
    //houses
    double IntercHouseTot, IntercHouseTotmm;
    double ChannelSedTot, ChannelDepTot, ChannelDetTot, TileVolTot;
    /// Sediment totals for mass balance and output (in kg)
    double MBs, DetTot, DetSplashTot, DetFlowTot, DepTot, SoilLossTot,SoilLossTotT, /*SoilLossTotOutlet,*/ SedTot, /*SoilLossTotSub,*/
           FloodDetTot, FloodDepTot, FloodSedTot;
    /// Water totals for output in file and UI (in mm), copied to 'op' structure
    double RainTotmm, SnowTotmm, IntercTotmm, WaterVolTotmm,WaterVolRunoffmm, InfilTotmm, Qtotmm, RainAvgmm, SnowAvgmm;
    double ChannelVolTotmm, floodTotmm, floodTotmmInit;
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

//    int c_outlet;  /// copy of outlet col number OBSOLETE
//    int r_outlet;  /// copy of outlet row number

//    int c_plot;  /// copy of col number of hydrograph plotted on screen
//    int r_plot;  /// copy of row number of hydrograph plotted on screen

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
    //MapListStruct SubsMaps[32];
    //int nrSubsMaps;

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
    void InitStandardInput(void); //VJ 170211

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
    void prepareFloodZ(cTMap *z);
    double fullSWOF2Do2(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct);
    double fullSWOF2Do2light(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct);
    double fullSWOF2Do1(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct);
    double limiter(double a, double b);
    void MUSCL(cTMap *ah, cTMap *au, cTMap *av, cTMap *az);
    void ENO(cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z);
    void simpleScheme(cTMap *_h, cTMap *_u, cTMap *_v);
    double maincalcflux(double dt, double dt_max);
    void maincalcscheme(double dt, cTMap *he, cTMap *ve1, cTMap *ve2,cTMap *hes, cTMap *ves1, cTMap *ves2);
    void Fr_Manning(double uold, double vold, double hnew, double q1new, double q2new, double dt, double N);
    void Fr_ManningSf(double h, double u, double v, double cf);
    void setZero(cTMap *_h, cTMap *_u, cTMap *_v);
    void F_HLL2(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    void F_HLL(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    void F_Rusanov(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    int F_scheme, F_fluxLimiter, F_reconstruction, F_replaceV, F_MaxIter;
    double F_maxVelocity;
    double F_extremeHeight;
    double F_extremeDiff;

    double F_levee;
    double HLL2_f1, HLL2_f2, HLL2_f3, HLL2_cfl, HLL_tmp;
    double q1man, q2man;
    //double dt_max, dt1;
    bool prepareFlood, startFlood;
    int verif, iter_n;
/*
    // floods watershed based
    double fullSWOF2Do2ws(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct);
    double fullSWOF2Do1ws(cTMap *h, cTMap *u, cTMap *v, cTMap *z, bool correct);
    void MUSCLws(int wsnr, cTMap *ah, cTMap *au, cTMap *av, cTMap *az);
    void ENOws(int wsnr, cTMap *_h, cTMap *_u, cTMap *_v, cTMap *_z);
    void simpleSchemews(int wsnr, cTMap *_h, cTMap *_u, cTMap *_v);
    void maincalcfluxws(int wsnr);
    void findDTws(int wsnr, bool two);
    void maincalcschemews(int wsnr, cTMap *he, cTMap *ve1, cTMap *ve2,cTMap *hes, cTMap *ves1, cTMap *ves2);
    void setZerows(int wsnr, cTMap *_h, cTMap *_u, cTMap *_v);
    void MakeWatersheds(void);
*/

    //SEDIMENT TRANSPORT
    int FS_SS_Method;
    int FS_BL_Method;
    double FS_SigmaDiffusion;

    int OF_Method;

    int R_SS_Method;
    int R_BL_Method;
    double R_SigmaDiffusion;

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
    void FS_Flux(cTMap * _sbl,cTMap * _sss,cTMap * _h1d,cTMap * _h1g,cTMap * _h2d,cTMap * _h2g,cTMap * _u1r,cTMap * _u1l,cTMap * _v1r,cTMap * _v1l,cTMap * _u2r,cTMap * _u2l,cTMap * _v2r,cTMap * _v2l);
    void FS_MUSCLE(cTMap * _sbl,cTMap * _sss);
    void FS_ENO(cTMap * _sbl,cTMap * _sss);
    void FS_Simple(cTMap * _sbl,cTMap * _sss);
    void FS_MainCalc(cTMap * _h, cTMap * _sbl,cTMap * _sbln,cTMap * _sss,cTMap * _sssn, double dt);

    void FS_FluxWS(int wsnr, cTMap * _sbl,cTMap * _sss,cTMap * _h1d,cTMap * _h1g,cTMap * _h2d,cTMap * _h2g,cTMap * _u1r,cTMap * _u1l,cTMap * _v1r,cTMap * _v1l,cTMap * _u2r,cTMap * _u2l,cTMap * _v2r,cTMap * _v2l);
    void FS_MUSCLEWS(int wsnr, cTMap * _sbl,cTMap * _sss);
    void FS_ENOWS(int wsnr, cTMap * _sbl,cTMap * _sss);
    void FS_SimpleWS(int wsnr, cTMap * _sbl,cTMap * _sss);
    void FS_MainCalcWS(int wsnr, cTMap * _h, cTMap * _sbl,cTMap * _sbln,cTMap * _sss,cTMap * _sssn, double dt);

    void FS_HLL(double h_L,double bl_L,double ss_L,double u_L,double v_L,double h_R, double bl_R,double ss_R,double u_R,double v_R);
    void FS_HLL2(double h_L,double bl_L,double ss_L,double u_L,double v_L,double h_R, double bl_R,double ss_R,double u_R,double v_R);
    void FS_Rusanov(double h_L,double bl_L,double ss_L,double u_L,double v_L,double h_R, double bl_R,double ss_R,double u_R,double v_R);

    void SWOFSedimentBalance();
    void SWOFSedimentBalanceWS(int l);

    void SWOFSedimentMaxC(int r, int c);//, cTMap * h,cTMap * u,cTMap * v);
    void SWOFSedimentCheckZero(int r, int c, cTMap * h,cTMap * u,cTMap * v);
    void SWOFSedimentSetConcentration(int r, int c, cTMap * h,cTMap * u,cTMap * v);

    void SWOFSedimentDiffusion(double dt, cTMap * h,cTMap * u,cTMap * v, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC);
    void SWOFSedimentDiffusionWS(int wsnr, double dt, cTMap * h,cTMap * u,cTMap * v, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC);


    double SWOFSedimentTCBL(int r,int c, int d, cTMap * h,cTMap * u,cTMap * v);
    double SWOFSedimentTCSS(int r,int c, int d, cTMap * h,cTMap * u,cTMap * v);

    void SWOFSedimentFlow(double dt, cTMap * h,cTMap * u,cTMap * v, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC);
    void SWOFSedimentFlowInterpolation(double dt, cTMap * h,cTMap * u,cTMap * v, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC);
    void SWOFSedimentDet(double dt,int r,int c, cTMap * h,cTMap * u,cTMap * v);
    void SWOFSedimentFlowWS(int wsnr, double dt, cTMap * h,cTMap * u,cTMap * v, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC);
    void SWOFSedimentFlowInterpolationWS(int wsnr, double dt, cTMap * h,cTMap * u,cTMap * v, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC);
    void SWOFSediment(double dt, cTMap * h,cTMap * u,cTMap * v);
    void SWOFSedimentWS(int l,double dt, cTMap * h,cTMap * u,cTMap * v);
    void SWOFSedimentLayerDepth(int r , int c, cTMap * h,cTMap * u,cTMap * v);

    double simpleSedCalc(double Qj1i1, double Qj1i, double Sj1i, double dt, double vol, double sed);
    double complexSedCalc(double Qj1i1, double Qj1i, double Qji1, double Sj1i,
                          double Sji1, double alpha, double dt, double dx);

    void routeSubstance(int pitRowNr, int pitColNr, cTMap *_LDD,
                                cTMap *_Q, cTMap *_Qn, cTMap *_Qs, cTMap *_Qsn,
                                cTMap *_Alpha, cTMap *_DX, cTMap*_Vol, cTMap*_Sed,cTMap*_VolStore, cTMap*_SedStore);


    double K2DSolvebyFluxSed(double dt, cTMap *M, cTMap *MC);
    double K2DSolvebyInterpolationSed(double dt, cTMap *M, cTMap *MC);


    double OFTC(int r, int c, int d);
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


    void RiverSedimentDiffusion(double dt, cTMap * _BL,cTMap * _BLC, cTMap * _SS,cTMap * _SSC);
    void RiverSedimentLayerDepth(int r , int c);
    double RiverSedimentTCBL(int r,int c, int d);
    double RiverSedimentTCSS(int r,int c, int d);
    void RiverSedimentMaxC(int r, int c);


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
    void InfiltrationFlood(void);
    void InfiltrationFloodNew(void);
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
    void OverlandFlow(void);
    void OverlandFlowNew(void);
    void ChannelFlow(void);
    double ChannelIterateWH(double _h, int r, int c);
    void ChannelWaterHeight(void);
    void ChannelWaterHeightFromVolume();
    void ToChannel(void);
    void ToFlood(void);
    void CalcVelDisch();
    void CalcVelDischChannel(void);
    void ToTiledrain(void);
    void TileFlow(void);
    void CalcVelDischTile(void);
    void GridCell(void);

    void ExtendChannel();
    bool IsExtendedChannel(int r, int c, int dr, int dc);
    void DistributeOverExtendedChannel(cTMap * _In, cTMap * _Out, bool do_not_div = false, bool proportional = true);
    bool OUTORMV(int r, int c);
    void InitFlowBarriers(void);
    double DEMFB(int r, int c, int rd, int cd, bool addwh);
    double FB(int r, int c, int rd, int cd);
    void SetFlowBarriers();
    void GetFlowBarrierData(QString name);
    double FBW(double h, int r, int c, int dr, int dc);

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
    void ChannelOverflow();
    void getFloodParameters(void);
    double courant_factor;
    double courant_factor_diffusive;
    double mixing_coefficient, runoff_partitioning;
   // double cfl_fix;
    double minReportFloodHeight;
    double correctMassBalance(double sum1, cTMap *M, double minV);
    void Kinematic(int pitRowNr, int pitColNr, cTMap *_LDD,
                   cTMap *_Q, cTMap *_Qn,
                   cTMap *_q, cTMap *_Alpha, cTMap *_DX,
                   cTMap *_Vol,
                   cTMap *_StorVol);
    double IterateToQnew(double Qin, double Qold, double q, double alpha, double deltaT, double deltaX);


    void upstream(cTMap *_LDD, cTMap *_M, cTMap *out);
   // void KinWave(cTMap *_LDD,cTMap *_Q, cTMap *_Qn,cTMap *_q, cTMap *_Alpha, cTMap *_DX);

    // kinematic 2D
    double K2DFlux();
    void K2DSolve(double dt);
    void K2DSolvebyFlux(double dt);
    void K2DSolvebyInterpolation(double dt);
    void K2DInit();
    void K2DCalcVelDisch();
    void K2DDEMA();
    double K2DQOut;
    double K2DQSOut;
    double K2DQPOut;
    double CourantKin;
    int ConcentrateKin;
    double TimestepKinMin;
    double KinematicBoundaryFraction = 0.05;
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
