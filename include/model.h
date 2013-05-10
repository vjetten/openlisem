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
#include <QMutex.h>

#include "csfmap.h"
#include "mmath.h"
#include "error.h"
#include "swatre_p.h"
#include "swatre_g.h"

#define OLDSWATRE 1

//---------------------------------------------------------------------------
#define PI 3.14159265

#define DEBUG(s) emit debug(QString(s))

//#define mwrite(name) WriteMap(QString(resultDir+name))
#define report(name) WriteMapSeries(resultDir,QString(name), printstep)

// defines to make life easier

/// shortcut to access data
#define Drc     Data[r][c]
#define Drci Data[r+dr[i]][c+dc[i]]

/// shortcut to access the outlet point data
#define DrcOutlet     Data[r_outlet][c_outlet]

/// shortcut to access the point data plotted on screen
#define DrcPlot     Data[r_plot][c_plot]
// VJ 110630 show hydrograph for selected output point

/// shortcut missing value in map
#define MV(r,c) IS_MV_REAL8(&LDD->Data[r][c])

/// shortcut for LDD row and col loop
#define FOR_ROW_COL_MV for (int r = 0; r < _nrRows; r++)\
    for (int c = 0; c < _nrCols; c++)\
    if(!IS_MV_REAL8(&LDD->Data[r][c]))

/// shortcut for LDD row and col loop in SWOF, rows/cols 1 to nrRows/nrCols-1
#define FOR_ROW_COL_MV_MV for (int r = 1; r < _nrRows-1; r++)\
    for (int c = 1; c < _nrCols-1; c++)\
    if(!IS_MV_REAL8(&LDD->Data[r][c]) && \
    !IS_MV_REAL8(&LDD->Data[r-1][c]) && \
    !IS_MV_REAL8(&LDD->Data[r+1][c]) && \
    !IS_MV_REAL8(&LDD->Data[r][c-1]) && \
    !IS_MV_REAL8(&LDD->Data[r][c+1]))
      //if (FloodDomain->Drc > 0)

/// shortcut for channel row and col loop
#define FOR_ROW_COL_MV_CH for (int  r = 0; r < _nrRows; r++)\
    for (int  c = 0; c < _nrCols; c++)\
    if(!IS_MV_REAL8(&LDDChannel->Data[r][c]))

/// shortcut for tile network row and col loop.
#define FOR_ROW_COL_MV_TILE for (int  r = 0; r < _nrRows; r++)\
    for (int  c = 0; c < _nrCols; c++)\
    if(!IS_MV_REAL8(&LDDTile->Data[r][c]))

/// shortcut to check if r,c is inside map boundaries, used in kinematic and flooding
#define INSIDE(r, c) (r>=0 && r<_nrRows && c>=0 && c<_nrCols)


#define NUMNAMES 360   /// \def NUMNAMES runfile namelist max
#define NUMMAPS 256    /// \def max nr maps
#define MIN_FLUX 1e-12 /// \def minimum flux (m3/s) in kinematic wave
#define MIN_HEIGHT 1e-6 /// \def minimum water height (m) for transport of sediment
#define MAXCONC 848    /// \def max concentration susp. sed. in kg/m3 0.32 * 2650 = max vol conc from experiments Govers x bulk density


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


//---------------------------------------------------------------------------
/// structure containing pointers to all maps

/** structure containing pointers to all maps for automatic destruction after runs
 so memory doesn't have to be freed for each map. The functions Newmap(double) and
ReadMap(cTMap *Mask, QString name) put a map on this list
*/
typedef struct MapListStruct {
    TMMap *m;
}  MapListStruct;
//---------------------------------------------------------------------------
/// linked list structure for network in kin wave
typedef struct LDD_LINKEDLIST {
    int rowNr;
    int colNr;
    struct LDD_LINKEDLIST *prev;
}  LDD_LINKEDLIST;
//---------------------------------------------------------------------------
/// structure used for sorting of the LDD
typedef struct LDD_POINT {
    int rowNr;
    int colNr;
}  LDD_POINT;
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
    double area;
    double totdet;
    double totdep;
    double totsl;
} UNIT_LIST;
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
  Every timestep the mass alance is calculated and output is reported to the UI and disk.
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

    /// map management structure, automatic adding and deleting of all TMMap variables
    MapListStruct maplistTMMap[NUMNAMES];
    int maplistnr;

    /// variable declaration list of all maps with comments:
#include "TMmapVariables.h"

    /// SwitchXXX are boolean options that are set in interface and runfile, mainly corrsponding to checkboxes in the UI
    bool SwitchHardsurface, SwatreInitialized, SwitchInfilGA2, SwitchLimitTC, SwitchLimitDepTC,
    SwitchWheelPresent, SwitchCompactPresent, SwitchIncludeChannel, SwitchChannelBaseflow,
    startbaseflowincrease, SwitchChannelInfil, SwitchAllinChannel, SwitchErosion, SwitchAltErosion,
    SwitchSimpleDepression, SwitchBuffers, SwitchBuffersImpermeable, SwitchSedtrap, SwitchSnowmelt, SwitchRainfall, SwitchRunoffPerM, SwitchInfilCompact,
    SwitchInfilCrust, SwitchGrassStrip, SwitchImpermeable, SwitchDumphead, SwitchWaterRepellency,
    SwitchWheelAsChannel, SwitchMulticlass, SwitchNutrients, SwitchGullies, SwitchGullyEqualWD, SwitchGullyInfil,
    SwitchGullyInit, SwitchOutputTimeStep, SwitchOutputTimeUser, SwitchMapoutRunoff, SwitchMapoutConc,
    SwitchMapoutWH, SwitchMapoutWHC, SwitchMapoutTC, SwitchMapoutEros, SwitchMapoutDepo, SwitchMapoutV,
    SwitchMapoutInf, SwitchMapoutSs, SwitchMapoutChvol, SwitchWritePCRnames, SwitchWriteCommaDelimited, SwitchWritePCRtimeplot,
    SwitchNoErosionOutlet, SwitchDrainage, SwitchPestout, SwitchSeparateOutput,
    SwitchInterceptionLAI, SwitchTwoLayer, SwitchSimpleSedKinWave, SwitchSoilwater, SwitchSOBEKoutput,
    SwitchPCRoutput, SwitchWriteHeaders, SwitchGeometric, SwitchIncludeTile, SwitchKETimebased, SwitchHouses, SwitchChannelFlood, SwitchRaindrum,
    SwitchFloodExplicit, SwitchFloodSWOForder1, SwitchFloodSWOForder1a, SwitchFloodSWOForder2, SwitchRoadsystem;

    // multiple options that are set in interface or runfile, see defines above
    /// Interception storage function based on LAI
    int InterceptionLAIType;
    /// infiltration method
    int InfilMethod;
    /// erosion units in output: to/ha; kg/cell; kg/m2
    int ErosionUnits;
    /// calibration factors
    int KEequationType;
    /// type of kinetic energy equation;
    double KEParamater_a1, KEParamater_b1, KEParamater_c1;
    double KEParamater_a2, KEParamater_b2;
    double KEParamater_a3, KEParamater_b3;
    /// parameters in KE equations

    double ksatCalibration;
    double nCalibration;
    double thetaCalibration;
    double psiCalibration;
    double ChnCalibration;
    double ChKsatCalibration;
    double SplashDelivery;
    double StripN;
    double StemflowFraction;
    double CanopyOpeness; // VJ 110209 added Aston factor as user input
    double waterRep_a;
    double waterRep_b;
    double waterRep_c;
    double waterRep_d;

    /// totals for mass balance checks and output
    /// Water totals for mass balance and output (in m3)
    double MB, Qtot, QtotOutlet, IntercTot, WaterVolTot, FloodVolTot, WaterVolSoilTot, InfilTot, RainTot, SnowTot, SurfStoremm, InfilKWTot;
    //houses
    double IntercHouseTot, IntercHouseTotmm;
    double ChannelVolTot, ChannelSedTot, ChannelDepTot, ChannelDetTot, TileVolTot;
    /// Sediment totals for mass balance and output (in kg)
    double MBs, DetTot, DetSplashTot, DetFlowTot, DepTot, SoilLossTot, SoilLossTotOutlet, SedTot, SoilLossTotSub;
    /// Water totals for output in file and UI (in mm), copied to 'op' structure
    double RainTotmm, SnowTotmm, IntercTotmm, WaterVolTotmm, InfilTotmm, Qtotmm, RainAvgmm, SnowAvgmm;
    double floodTotmm;
    /// peak times (min)
    double RainpeakTime, SnowpeakTime, QpeakTime, Qpeak, Rainpeak, Snowpeak;
    double BufferVolTot, BufferSedTot, BufferVolTotInit, BufferSedTotInit, BulkDens, BufferVolin;
    double nrCells, CatchmentArea;
    double QPlot, QtotPlot, QpeakPlot, SoilLossTotPlot;

    int c_outlet;  /// copy of outlet col number
    int r_outlet;  /// copy of outlet row number

    int c_plot;  /// copy of col number of hydrograph plotted on screen
    int r_plot;  /// copy of row number of hydrograph plotted on screen

    /// time and dx parameters
    double time, BeginTime, EndTime;
    double _dt, _dx;
    long runstep, printstep, printinterval;

    /// timeseries variables and output strings
    //double **RainfallSeries;
    int nrRainfallseries;
    int nrSnowmeltseries;
    //double **SnowmeltSeries;
    //int nrSnowmeltstations, nrSnowmeltseries;
    //RAIN_LIST *RainfallSeriesM;
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
    QString totalSoillossFileName;
    QString totalLandunitFileName;
    QString rainFileName;
    QString rainFileDir;
    QString snowmeltFileName;
    QString snowmeltFileDir;
    QString resultFileName;
    QString temprunname;
    QStringList outputcheck;
    /// standard names of output map series
    QString Outrunoff, Outconc, Outwh, Outrwh, Outtc, Outeros, Outdepo, Outvelo, Outinf, Outss, Outchvol,
    OutTiledrain, OutHmx, OutVf, OutQf;

    // list with class values of land unit map
    UNIT_LIST unitList[512]; // just a fixed number for 512 classes, who cares!
    QVector <UNIT_LIST> unitListM;
    int landUnitNr;
    // data initialization, runfile reading and parsing
    NAME_LIST runnamelist[NUMNAMES]; // structure for runfile variables and names
    int nrrunnamelist;

    // list of pointers for substance maps: sediment, sed classes, nutrients etc.
    // used in kin wave for routing of substances
    MapListStruct SubsMaps[32];
    int nrSubsMaps;

    // functions in lisDataInit.cpp
    void InitMapList(void);
    TMMap *NewMap(double value);
    TMMap *ReadMap(cTMap *Mask, QString name);
    void DestroyData(void);
    TMMap *InitMask(QString name);
    TMMap *InitMaskChannel(QString name);
    TMMap *InitMaskTiledrain(QString name);
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
    void ParseRunfileData(void);
    void GetRunFile(void);
    //MapListStruct qx[9];

    //FLOOD according to LISFLOOD
    double floodExplicit();
    //FLOOD according to FULLSWOF2D
    double fullSWOF2D(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2);
    double fullSWOF2Do1(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2);
    double fullSWOF2Do1a(TMMap *h, TMMap *u, TMMap *v, TMMap *z, TMMap *q1, TMMap *q2);
    double limiter(double a, double b);
    void MUSCL(TMMap *h,TMMap *u,TMMap *v,TMMap *z);
    void ENO(TMMap *h,TMMap *u,TMMap *v,TMMap *z);
    void simpleScheme(TMMap *_h,TMMap *_u,TMMap *_v);
    double bloc1(double dt, double dt_max);
    void bloc2(double dt, TMMap *he, TMMap *ve1, TMMap *ve2,TMMap *hes, TMMap *ves1, TMMap *ves2);
    void Fr_Manning(double uold, double vold, double hnew, double q1new, double q2new, double dt, double N);
    void Fr_ManningSf(double h, double u, double v, double cf);
    void setZero(TMMap *_h, TMMap *_u, TMMap *_v);
    void F_HLL2(double hg,double ug,double vg,double hd,double ud,double vd);
    void F_HLL(double hg,double ug,double vg,double hd,double ud,double vd);
    void F_Rusanov(double hg,double ug,double vg,double hd,double ud,double vd);
    int F_scheme;
    double F_levee;
    double HLL2_f1, HLL2_f2, HLL2_f3, HLL2_cfl;
    double dt1, dx, dy, dt_max, tx, ty;
    double q1mod, q2mod, Sf1, Sf2;
    bool prepareFlood, startFlood;
    int verif;
    int iter_n;

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
    void InterceptionHouses(void);
    /// subtract water retained on houses, for urban projects    
    void addRainfallWH(void);
    /// add net rainfall to WH, WHroads and WHgrass

    void Infiltration(void);
    void InfiltrationFlood(void);
    void InfilSwatre(void);
//    void InfilGreenAmpt1(TMMap *_WH);
//    void InfilSmithParlange1(TMMap *_WH);
    void InfilMorelSeytoux1(TMMap *_WH);
//    void InfilKsat(TMMap *_WH);
    double IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p, REAL8 *FFull);
    void SoilWater(void);
    void InfilMethods(TMMap *_Ksateff, TMMap *_WH, TMMap *_fpot, TMMap *_fact, TMMap *_L1, TMMap *_L2, TMMap *_FFull);
    void SurfaceStorage(void);
    void OverlandFlow(void);
    void ChannelFlow(void);
    void ChannelWaterHeight(void);
    void ToChannel(void);
    void ToFlood(void);
    void CalcVelDisch(void);
    void CalcVelDischChannel(void);
    void ToTiledrain(void);
    void TileFlow(void);
    void CalcVelDischTile(void);
    void GridCell(void);
    void SplashDetachment(void);
    void FlowDetachment(void);
    double MaxConcentration(double watvol, double sedvol);
    void ChannelFlowDetachment(void);
    //flood
    void ChannelFlood(void);
    void ChannelOverflow(void);
    double courant_factor;
    double cfl_fix;

    void Kinematic(int pitRowNr, int pitColNr, TMMap *_LDD, TMMap *_Q, TMMap *_Qn, TMMap *_Qs,
                   TMMap *_Qsn, TMMap *_q, TMMap *_Alpha, TMMap *_DX, TMMap *Vol, TMMap*SedVol,
                   TMMap *_StorVol, TMMap*_StorVolSed);

    double simpleSedCalc(double Qj1i1, double Qj1i, double Sj1i, double dt, double vol, double sed);
    double complexSedCalc(double Qj1i1, double Qj1i, double Qji1, double Sj1i,
                          double Sji1, double alpha, double dt, double dx);
    double IterateToQnew(double Qin, double Qold, double q, double alpha, double deltaT, double deltaX);

    void routeSubstance(int pitRowNr, int pitColNr, TMMap *_LDD,
                                TMMap *_Q, TMMap *_Qn, TMMap *_Qs, TMMap *_Qsn,
                                TMMap *_Alpha, TMMap *_DX, TMMap*_Vol, TMMap*_Sed);
    void findFlood(int pitRowNr, int pitColNr, TMMap *_LDD);

    // alternative kin wave based on a pre-sorted network
    bool useSorted;
    LDD_POINT **makeSortedNetwork(TMMap *_LDD, long *lddlistnr);
    void KinematicSorted(LDD_POINT **_lddlist, long _lddlistnr,
                         TMMap *_Q, TMMap *_Qn, TMMap *_Qs, TMMap *_Qsn,
                         TMMap *_q, TMMap *_Alpha, TMMap *_DX, TMMap *Vol, TMMap*SedVol,
                         TMMap *_StorVol, TMMap*_StorVolSed);

    //   QList <LDD_POINT *> listldd;
    //VJ 110123 sorted networks for faster kin wave
    LDD_POINT **lddlist;
    long lddlistnr;
    LDD_POINT **lddlistch;
    long lddlistchnr;
    LDD_POINT **lddlisttile;
    long lddlisttilenr;
    // QVector <TMMap> Substance;

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

    SOIL_MODEL *InitSwatre(TMMap *profileMap);//, QString initHeadMaps, TMMap *tiledepthMap, double dtMin);
    void SwatreStep(SOIL_MODEL *s, TMMap *_WH, TMMap *_fpot, TMMap *_drain, TMMap *_theta, TMMap *where);
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
signals:
    void done(const QString &results);
    void debug(const QString &results);
    void show(); //use the output structure "op" declared in global.h and LisUIoutput.h

private slots:   //note, was private loop but dixygen does not recognize that
    /// the main model loop, from here all processes are called in a time loop
    void DoModel();

};



#endif
