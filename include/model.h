/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 1992, 2003, 2016, 2024  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program. If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout, Meindert Commelin
**  Developed in: MingW/Qt/, GDAL, PCRaster
**  website, information and code: https://github.com/vjetten/openlisem
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
#include <omp.h>

#include "CsfMap.h"
#include "pcrtypes.h"
#include "io.h"
#include "lerror.h"
#include "swatre_p.h"



//---------------------------------------------------------------------------
//#define NULL nullptr

#define SHOWDEBUG (showr == 270 && showc == 318)

#define PI 3.14159265

#define HMIN 1e-6
#define DO_SEDDEP 1
#define GRAV 9.8067

#define he_ca 1e-12
#define ve_ca 1e-12

#define EPSILON 1e-10

#define GRAV_DEM 4.90335

#define BETArect 0.6
#define BETAcirc 0.6


#define Aavg(a,b)  (0.5*(a+b))
#define Savg(a,b)  sqrt(a*b)
#define Havg(a,b,w1,w2)  ((w1+w2)/(w1/a+w2/b))  //  sum (weight/variable) / sum weights
#define Mavg(a,b)  std::min(a,b)
#define SQR(a) ((a)*(a))

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
// slower
//#define FOR_ROW_COL_MV_L for(long i_ = 0; i_ < _nrCols*_nrRows; i_++)\
//{int r = i_/_nrCols; int c=i_%_nrCols;\
//if(!pcr::isMV(LDD->data[r][c]))


// #define FOR_ROW_COL_MV_L for(long i_ = nrValidCells-1; i_ >= 0; i_--)\
//  {int r = cr_[i_]->r; int c = cr_[i_]->c;
#define FOR_ROW_COL_MV_L for(long i_ = 0; i_ < nrValidCells; i_++)\
 {int r = cr_[i_].r; int c = cr_[i_].c;

#define FOR_ROW_COL_LDD5 for(long i_ = nrValidCellsLDD5-1; i_ >= 0; i_--)\
{int r = crldd5_[i_].r; int c = crldd5_[i_].c;

#define FOR_ROW_COL_LDDCH5 for(long i_ = nrValidCellsLDDCH5-1; i_ >= 0; i_--)\
{int r = crlddch5_[i_].r; int c = crlddch5_[i_].c;

// #define FOR_ROW_COL_MV_CHL for(long i_ = nrValidCellsCH-1; i_ >= 0; i_--)\
// {int r = crch_[i_]->r; int c = crch_[i_]->c;
#define FOR_ROW_COL_MV_CHL for(long i_ = 0; i_ < nrValidCellsCH; i_++)\
{int r = crch_[i_].r; int c = crch_[i_].c;

#define FOR_ROW_COL_MV_TILEL for(long i_ = nrValidCellsTile-1; i_ >= 0; i_--)\
{int r = crtile_[i_].r; int c = crtile_[i_].c;

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
#define MIN_FLUX 1e-6 /// \def minimum flux (m3/s)
#define MIN_HEIGHT 1e-6 /// \def minimum water height (m) for transport of sediment
#define MAXCONC 848.0    /// \def max concentration susp. sed. in kg/m3 0.32 * 2650 = max vol conc from experiments Govers x bulk density
#define MAXCONCBL 848.0    /// \def max concentration susp. sed. in kg/m3 0.32 * 2650 = max vol conc from experiments Govers x bulk density
#define MIN_SLOPE 1e-3

#define INFIL_NONE 0
#define INFIL_SWATRE 1
#define INFIL_SOAP 2
#define INFIL_GREENAMPT 3
#define INFIL_GREENAMPT2 4
//#define INFIL_KSAT 5
//#define INFIL_MOREL 21
#define INFIL_SMITH 22
#define INFIL_SMITH2 23

#define KE_EXPFUNCTION 0
#define KE_LOGFUNCTION 1
#define KE_POWERFUNCTION 2

#define MINMOD 1
#define VANALBEDA 3
#define VANLEER 2

#define FSGOVERS 0
#define FSHAIRSINEROSE 1

#define FSRIJN 1
#define FSRIJNFULL 2
#define FENGELUND 3
#define FSWUWANGJIA 10
//#define FSWUWANGJIABL 3

#define K2D_METHOD_KIN   1
#define K2D_METHOD_KINDYN  3
#define K2D_METHOD_DYN   2


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
    int nr;
    double V;
}  IDI_POINT;
//---------------------------------------------------------------------------list
typedef struct LDD_COOR {
    int r;
    int c;
}  LDD_COOR;
//---------------------------------------------------------------------------list
typedef struct LDD_COORldd {
    int r;
    int c;
    int ldd;
}  LDD_COORldd;
//---------------------------------------------------------------------------
typedef struct LDD_COORIN {
    int r;
    int c;
    int ldd;
    QVector <LDD_COOR> inn;
    int nr;
    //LDD_COOR *inn;
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
/// vec6 used for muscl
typedef struct vec6 { double v[6]; } vec6;
//---------------------------------------------------------------------------
/// Structure to store rain station values of rainfile mapnames
typedef struct RAIN_LIST {
    double time;
    QList <int> stationnr;
    QVector <double> intensity;
} RAIN_LIST;
//---------------------------------------------------------------------------
/// Structure to store meteo station values of rainfile mapnames
typedef struct METEO_LIST {
    double time;
    QString name;
    double calib;
} METEO_LIST;
//---------------------------------------------------------------------------
/// Structure to store discharge station values
typedef struct Q_LIST {
    double time;
    QList <int> stationnr;
    QVector <double> Qin;
} Q_LIST;
//---------------------------------------------------------------------------
/// Structure to store boundary water level
typedef struct WH_LIST {
    double time;
    double WH;
} WH_LIST;
//---------------------------------------------------------------------------
typedef struct BUFFER_LIST {
    double area;
    double h;
    int ID;
} BUFFER_LIST;

// typedef struct ExtCH {
//     QList <int> childRow;
//     QList <int> childCol;
//     int chRow;
//     int chCol;
//     bool isExtended;
// } ExtCH;

typedef struct SOIL_LIST {
    int c;
    int r;
    double dts;
    double dtsum;
    bool ponded;
    double drain;
    double Infact;
    double InfPot;
    double SD;

    QVector <double> pore;
    QVector <double> Ks;
    QVector <double> z;
    QVector <double> dz;
    QVector <double> h;
    QVector <double> hb;
    QVector <double> lambda;
    QVector <double> theta;
    QVector <double> thetar;
    QVector <double> rootz;
    QVector <double> vg_alpha;
    QVector <double> vg_n;

} SOIL_LIST;
//---------------------------------------------------------------------------
typedef struct DRAIN_PROP {
    int r;
    int c;
    int ldd;
    double Afull;
    double Qfull;
    double beta, Beta1;
    double sMax, sFull;
    double dxdt;
    double ain, aout;
    double qin, qout;
    double C1, C2;
    double a1, a2;
    double q1, q2;
    double diam;
}  DRAIN_PROP;

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
    int nNodes, nN1_, nN2_, nN3_;
    int nrSoilLayers;
    double SoilWBdtfactor;
    int KavgType;

    long nrValidCells;
    long nrValidCellsLDD5;
    long nrValidCellsCH;
    long nrValidCellsLDDCH5;
    long nrValidCellsWS;
    long nrValidCellsTile;
    QVector <LDD_COOR> cr_;
    QVector <LDD_COOR> crch_;
    QVector <LDD_COORIN> crlinkedldd_;
    QVector <LDD_COORIN> crlinkedlddch_;
    QVector <LDD_COORIN> crlinkedlddbase_;
    QVector <LDD_COORIN> crlinkedlddtile_;

    QVector <LDD_COOR> crldd5_;
    QVector <LDD_COOR> crlddch5_;
    QVector <LDD_COOR> crtile_;
    QVector <LDD_COORout> crout_;
    QVector <LDD_COORldd> dcr_;

    // vector of soil structure
    QVector <SOIL_LIST> crSoil;

  //  QVector <IDI_POINT> IDIpoints;
    QVector <IDI_POINT> IDIpointsRC;
    QList <int> stationID;
    QList <int> stationQID;
    QVector <double> IDIpointsV;

    /// map management structure, automatic adding and deleting of all cTMap variables
    //MapListStruct maplistCTMap[NUMNAMES];
    QVector <cTMap*> maplistCTMap;
    int maplistnr;

    /// variable declaration list of all maps with comments:
#include "TMmapVariables.h"

    /// SwitchXXX are boolean options that are set in interface and runfile, mainly corrsponding to checkboxes in the UI

    bool
        // input rniafall and inflow
        SwitchRainfall,
        SwitchEventbased,
        SwitchIDinterpolation,
        SwitchUseIDmap,
        SwitchDailyET,
        SwitchRainfallSatellite,
        SwitchIncludeET,
        SwitchETSatellite,
        SwitchSnowmelt,
        SwitchSnowmeltSatellite,

        SwitchInterception,
        SwitchInfiltration,
        // channel and Overland flow
        SwitchIncludeChannel,
       // SwitchChannelBaseflow,
        SwitchChannelBaseflowStationary,
        SwitchChannelAdjustCHW,
        SwitchChannelInfil,
        SwitchGWflow,
        SwitchGW2Dflow,
        SwitchGWSWOFflow,
        SwitchLDDGWflow,
        SwitchSWATGWflow,
        SwitchChannel2DflowConnect,
        SwitchChannelWFinflow,
        SwitchGWChangeSD,
        SwitchDischargeUser,
        SwitchWaveUser,
        SwitchIncludeDiffusion,
        SwitchIncludeRiverDiffusion,
        SwitchFloodInitial,
        SwitchFlowBarriers,
        SwitchBuffers,
        SwitchCulverts,
        SwitchLitter,

        // output
        //SwitchOutputTimeUser,
        SwitchWriteCommaDelimited,
        SwitchWritePCRtimeplot,
        SwitchSeparateOutput,
        SwitchWriteHeaders,
        SwitchEndRun,
        SwitchResultDatetime,
        SwitchOutputTimestamp,

        // erosion,
        SwitchErosion,
        SwitchSlopeStability,
        SwitchKETimebased,
        SwitchUse2Phase,


        // infiltration,
        SwitchInfilCompact,
        SwitchInfilCrust,
        SwitchDynamicCrusting,
        SwitchImpermeable,
        SwitchDumphead,
        SwitchGeometric,
        SwitchTwoLayer,
        SwitchThreeLayer,
        SwitchHinit4all,
        SwitchOMCorrection,
        SwitchDensCorrection,
        //SwitchWaterRepellency,
        //SwitchInterceptionLAI,
        SwitchPsiUser,
        SwitchNrLayers,
        //SwitchDumpH,
        //SwitchDumpTheta,
        //SwitchDumpK,
        SwitchVanGenuchten,
        SwitchBrooksCorey,

        // tiles and drains, buildings
        SwitchRoadsystem,
        SwitchHardsurface,
        SwitchIncludeTile,
        SwitchIncludeStormDrains,
        SwitchStormDrainCircular,
        SwitchUseSWMMflow,
        SwitchHouses,
        SwitchInfrastructure,
        SwitchRaindrum,
        SwitchAddBuildingsDEM,

        SwitchConservation,
        SwitchGridRetention,
        SwitchSedtrap,
        SwitchGrassStrip,


        // advanced
        SwitchAdvancedOptions,
        SwitchTimeavgV,
        SwitchCorrectMB_WH,
        SwitchCorrectDEM,
        Switch2DDiagonalFlow,
        SwitchSWOFopen,
        SwitchMUSCL,
        SwitchUserCores,
        SwitchVariableTimestep,
        SwitchHeun,
        SwitchImage,
        SwitchChannelKinwaveDt,
        SwitchChannelKinwaveAvg,
        SwitchLinkedList,
        SwitchPerimeterKW,
        SwitchChannelKinWave,
        SwitchChannelMaxV;

    // // TODO multi class sed
    bool SwitchAdvancedSed,
         SwitchUseMaterialDepth,
         SwitchNoBoundarySed,
         SwitchUseGrainSizeDistribution,
         SwitchEstimateGrainSizeDistribution,
         SwitchReadGrainSizeDistribution,
         SwitchD50CHavg;

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

    // extended channel not used
    // QVector <ExtCH> ExtChannel;
    // ExtCH ExtendCH;

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

    //Groundwater flow parameters
    double GW_recharge;
    double GW_flow;
    double GW_slope;
    double GW_deep;
    double GW_threshold;
    double GW_initlevel;

    double rainfallETa_threshold;
    double rainIDIfactor;

    // calibration parameters
    double gsizeCalibrationD50;
    double gsizeCalibrationD90;
    double SmaxCalibration;
    double RRCalibration;
    double ksatCalibration;
    double ksat2Calibration;
    double ksat3Calibration;
    double nCalibration;
    double thetaCalibration;
    double psiCalibration;
//    double SD1Calibration;
//    double SD2Calibration;
    double ChnCalibration;
    double WaveCalibration;
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
    double CanopyOpeness;
    double TurbulenceFactor;
    double TileDrainDistance;

    //sed transport equations
    int FS_SS_Method;
    int FS_BL_Method;
    double FS_SigmaDiffusion;
    int R_SS_Method;
    int R_BL_Method;
    double R_SigmaDiffusion;

    /// totals for mass balance checks and output
    /// Water totals for mass balance and output (in m3)
    double MB, MBeM3, Qtot, Qtot_dt, QTiletot, QTile, IntercTot, IntercETaTot, WaterVolTot, RetentionVolTot, WaterVolSoilTileTot, InfilTot, RainTot, SnowTot, theta1tot, theta2tot;
    double SurfStoremm, InfilKWTot,BaseFlowTot,BaseFlowInit, BaseFlowInitmm, BaseFlowTotmm, PeakFlowTotmm, Qfloodout, QfloodoutTot, QuserInTot;
    double QBoundaryTot, floodVolTot, floodVolTotInit, floodVolTotMax, floodAreaMax, floodArea, floodBoundarySedTot, ChannelVolTot, ChannelVolTotmm, WHinitVolTot,StormDrainVolTot;
    double IntercHouseTot, IntercHouseTotmm, IntercLitterTot, IntercLitterTotmm;
    double ChannelSedTot, ChannelDepTot, ChannelDetTot, SoilMoistTot, SoilMoistDiff, SoilMoistTotmm, QSideVolTot;
    /// Sediment totals for mass balance and output (in kg)
    double MBs, DetTot, DetSplashTot, DetFlowTot, DepTot, SoilLossTot, SoilLossTot_dt, SedTot,
           FloodDetTot, FloodDepTot, FloodSedTot;
    /// Water totals for output in file and UI (in mm), copied to 'op' structure
    double RainTotmm, SnowTotmm, IntercTotmm, IntercETaTotmm, WaterVolTotmm, WaterVolRunoffmm, Qboundtotmm, InfilTotmm, Qtotmm, RainAvgmm, SnowAvgmm, GWdeeptot;
    double StormDrainTotmm, floodVolTotmm, floodTotmmInit;
    /// peak times (min)
    double RainstartTime, RainpeakTime, SnowpeakTime, QpeakTime, Qpeak, Rainpeak, Snowpeak;
    bool rainStarted;
    bool ETStarted;
    double ETstartTime;
    double BulkDens;
    double nrCells, CatchmentArea, nrFloodedCells;
    double LitterSmax, ETaTot, ETaTotmm, ETaTotVol, GWlevel, GWleveltot;
    double thetai1tot, thetai2tot, thetai1cur, thetai2cur;
    double maxRainaxis;
    double latitude;
    double HinitValue;

    ///pesticides
    double MBp,PestMassApplied, PestLossTotOutlet, PestFluxTotOutlet, PestRunoffSpatial, PestDisMixing, PestSorMixing, PestInfilt, PestStorage, Pestdetach, PestCinfilt,PestCfilmexit;
    double MBpex,PestRunoffSpatialex,PestDisMixingex,PestSorMixingex,PestInfiltex,PestLossTotOutletex;
    int N_SPK;
    double Maxsolubility;
    double MaxVup;

    /// swatre
    double SwatrePrecision;

    /// time and dx parameters
    double time, BeginTime, EndTime;
    double _dt, _dx;
    double _dt_user, _dtCHkin, _CHMaxV;
    long runstep, printstep, printinterval;
    double _llx, _lly;

    QString mapFormat; //Gtiff or pcraster

    /// timeseries variables and output strings
    int nrDischargeseries;
    int nrRainfallseries;
    int nrETseries;
    int nrSnowmeltseries;
    int nrWHseries;
    int currentRainfallrow;
    int currentETrow;
    int currentSnowmeltrow;
    int currentDischargerow;
    int currentWHrow;
    QVector <RAIN_LIST> RainfallSeries;  // rainfall vector of records
    QVector <RAIN_LIST> ETSeries;
    QVector <RAIN_LIST> SnowmeltSeries;
    QVector <Q_LIST> DischargeSeries;
    QVector <WH_LIST> WHSeries;
    QVector <METEO_LIST> RainfallSeriesMaps;  // rainfall vector of records
    QVector <METEO_LIST> ETSeriesMaps;  // rainfall vector of records
    QVector <METEO_LIST> SnowmeltSeriesMaps;  // rainfall vector of records
    QVector <LDD_COORloc> crQin_;
    QVector <int> locationnnrsrec;
    QVector <double> raintime; // times in seconds of records for fast searching where we are
    QVector <double> ETtime;
    QVector <double> WHtime;
    QVector <double> dischargetime;
    QVector <double> snowmelttime;

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
    QString WaveinFileName;
    QString WaveinFileDir;
    QString resultFileName;
    QString temprunname;
    /// standard names of output map series
    QString Outrunoff, Outconc, Outwh, Outrwh, Outvelo, Outinf, Outss, Outchvol,
    Outtc, Outeros, Outdepo, OutSL, OutSed, OutInt,OutSedSS, OutSedBL,
    OutTiledrain, OutTileVol,OutTileV, OutHmx, OutVf, OutQf, OutHmxWH, OutTheta1, OutTheta2, OutGW;
    bool  SwitchOutrunoff, SwitchOutconc, SwitchOutwh, SwitchOutrwh, SwitchOutvelo, SwitchOutinf, SwitchOutss, SwitchOutchvol,
    SwitchOutConc, SwitchOutTC, SwitchOutDet, SwitchOutDep, SwitchOutSL, SwitchOutSed, SwitchOutInt, SwitchOutSedSS, SwitchOutSedBL,
    SwitchOutTiledrain, SwitchOutTileVol, SwitchOutHmx, SwitchOutVf, SwitchOutQf, SwitchOutHmxWH, SwitchOutTheta, SwitchOutGW;
    QString errorFileName;
    QString errorSedFileName;
    QString satImageFileName;
    QString satImageFileDir;

    // list with class values of land unit map
    UNIT_LIST unitList[NRUNITS];
    UNIT_LIST floodList[NRUNITS];
    int landUnitNr;

    // CENTRAL STRUCTURE WITH ALL MAP POINTERS
    // data initialization, runfile reading and parsing
    NAME_LIST runnamelist[NUMNAMES]; // structure for runfile variables and names
    int nrrunnamelist;

    // => functions in lisRunfile.cpp
    QString getvaluename(QString vname);
    double getvaluedouble(QString vname);
    int getvalueint(QString vname);
    QString getvaluestring(QString vname);
    QString CheckDir(QString p, bool makeit = false);
    QString GetName(QString p);
    QString checkOutputMapName(QString p, QString S, int i);
    void ParseRunfileData(void);
    void GetRunFile(void);
    // <= runfile

    // => functions in lisDataInit.cpp
    void InitMapList(void);
    cTMap *NewMap(double value);
    cTMap *ReadMap(cTMap *Mask, QString name);
    void DestroyData(void);
    cTMap *InitMask(QString name);
    cTMap *InitMaskChannel(QString name);
    cTMap *InitMaskTiledrain(QString name);
    void InitTiledrains(void);
    void InitBuffers(void);
    void InitChannel(void);
    void InitBoundary(void);
    void InitShade(void);
    void InitImages(void);
    void InitErosion(void);
    void GetInputData(void);      // get and make input maps
    void InitParameters(void);
    void IntializeData(void);     // make all non-input maps
    void IntializeOptions(void);  // set all options to false etc
    void InitStandardInput(void);
    void InitLULCInput(void);
    void InitSoilInput(void);
    void InitFlood(void);
    void InitMeteoInput(void);
    void InitScreenChanNetwork();
    void CorrectDEM(cTMap *h, cTMap * g);
    void DiagonalFlowDEM();
    void calcSoilPhysics(cTMap *Ksat, cTMap *lambda, cTMap *thfc, cTMap *thr,
                                 cTMap *psi, cTMap *psiae, double calk, double calpsi);
    // <= initiatlisation


    //int GrainSizeDistributionType;

    double LogNormalDist(double d50,double sigma, double d); // not used
    double DetachMaterial(int r,int c, int d,bool channel,bool flood,bool bl, double detachment); //not used

    //material that is available for detachment
    QList<cTMap *> StorageDep_D;
    QList<cTMap *> Storage_D;
    cTMap *Storage;
    cTMap *StorageDep;
    cTMap *SedimentMixingDepth;
    cTMap *maxDetachment;

    //QList<cTMap *> RStorageDep_D;
    //QList<cTMap *> RStorage_D;
    cTMap *RStorage;
    cTMap *RStorageDep;
    cTMap *RSedimentMixingDepth;

    //keep track of any dissolved substances that need to be advected by the kinematic wave
    //not used!!!
    //QList<cTMap *> OF_Advect;
    //QList<cTMap *> R_Advect;
    //QList<cTMap *> F_Advect;
    bool addedbaseflow;

    // TODO PEST stuff, replace with work Meindert
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
    // <= PEST

    // 1D hydro processes
    // => input timeseries
    void GetInputTimeseries();
    void GetUserDischargeData(QString name);
    void GetWHboundaryData(QString name);
    void GetRainfallStationData(QString name);   // get input timeseries
    void GetSpatialMeteoData(QString name, int type);   // get input timeseries
    void GetETStationData(QString name);   // get input timeseries
    void GetSnowmeltData(QString name);   // get input timeseries
    double getTimefromString(QString sss);
    double getmaxRainfall();
    void IDInterpolation();
    void GetRainfallMapfromStations(double currenttime);
    void GetRainfallMapfromSat(double currenttime);
    void GetETSatMap(double currenttime);
    void GetETMapfromStations(double currenttime);
    void GetSnowmeltMap(void);
    // user defined water height at boundary
    void GetWHboundaryMap(double currenttime);
    // user defined input discharge, e.g. dam spill
    void GetDischargeMapfromStations(double currenttime);
    // <= input timeseries

    // => not used, replaced by cell_[process]
    void Interception();
    void SoilWater();
    // <= not used

    double SoilWaterMass();

    // => TODO: SOAP infil model, swatre works better for now
    void cell_Soilwater(long i_); //SOAP
    double calcSinkterm(long i_,  double WH, double *S);
    void calcSinktermSWATRE(long i_,  PIXEL_INFO *pixel, double *h, double *S);
    double calculateDayLength(double latitude, int dayNumber);
    void VanGenuchten(SOIL_LIST s, double Hnew[], double K[], double C1[], bool analytical);
    void BrooksCorey(SOIL_LIST s, double Hnew[], double K[], double C1[], bool analytical);
    void getThetafromH(int j, SOIL_LIST s);
    void getHfromTheta(int j, SOIL_LIST s);
    // <= SOAP

    // => vertical hydro processes, OMP
    void GridCell();
    void HydrologyProcesses();
    void cell_Interception(int r, int c);
    void cell_SurfaceStorage(int r, int c);
    void cell_InfilMethods(int r, int c);
    void cell_SWATRECalc(long i_); // not used, too complex
    double cell_Percolation(int r, int c, double factor);
    double cell_PercolationMulti(int r, int c, double factor);
    void cell_Redistribution0(int r, int c);
    void cell_Redistribution1(int r, int c);
    void cell_Redistribution2(int r, int c);
    void cell_Tiledrain1(int r, int c);
    void cell_Tiledrain2(int r, int c);
    void cell_Channelinfow1(int r, int c);
    void cell_Channelinfow2(int r, int c);
    void cell_SplashDetachment(int r, int c);
    void cell_FlowDetachment(int r, int c);
    void cell_FlowDetachmentContinuous(int r, int c);
    void cell_ETa(int r, int c);
    double getETaFactor();
    double ETafactor;
    double ETafactorTot;
    void InfilEffectiveKsat();
    void InfilDynamicCrusting();
    void InfilSwatre();
    void InfilMethods(cTMap *_Ksateff, cTMap *_WH, cTMap *_fpot, cTMap *_fact, cTMap *_L1, cTMap *_L2, cTMap *_FFull);
    double IncreaseInfiltrationDepthNew1(double fact_, int r, int c);
    double IncreaseInfiltrationDepthNew2(double fact_, int r, int c);
    double IncreaseInfiltrationDepthNew3(double fact_, int r, int c);
    void avgTheta();
    // <= vertical processes

    // => 1D and 2D overlandflow
    void OverlandFlow();
    void CalcVelDisch(); //(int r, int c);
    void OverlandFlow1D(void);
    void OverlandFlow2D();
    void ToChannel();//int r, int c);
    void ToFlood();
    void ToTiledrainAll();
    // <= OF

    //SWMM pipe flow
    void PipeFlowSWMM();
    double getAfromS(DRAIN_PROP *dr, double s);
    int findroot_Newton(DRAIN_PROP *dr, double x1, double x2);
    int solveContinuity(DRAIN_PROP *dr);
    double solve_theta(double r, double psi_target);
    double psi_rel(double r, double theta);

    // => 1D flow on network
    void FindStationaryBaseFlow();
    void ChannelFlow();
    void IterateCulvert(DRAIN_PROP *dr);
    double findroot_NewtonCH(DRAIN_PROP *dr,double x1, double x2, double xacc);
    void ChannelBaseflow();
    void ChannelRainandInfil();
    void ChannelSedimentFlow();
    void ChannelFlowandErosion();
    void ChannelVelocityandDischarge();
    double pipeThetafroma(int r, int c, double a);
    void ChannelFlood(void);
    void ChannelOverflow(cTMap *_h, cTMap *_V);
    void ChannelOverflowIteration(cTMap *_h, cTMap *_V);

    // tiles/stormdrains
    void TileFlow(void);
    void TileFlowSWMM(void);
    void CalcVelDischRectangular(void);
    void CalcVelDischCircular(void);
    double getMassCH(cTMap *M);
    void correctMassBalanceCH(double sum1, cTMap *M);

    // <= 1D flow

    // => 2D flow according to FULLSWOF2D
    double Flood_DTMIN;
    int F_scheme, F_fluxLimiter, F_MaxIter, F_AddGravity;
    double F_minWH;
    double F_pitValue;
    bool prepareFlood, startFlood;
    int iter_n;
    double fullSWOF2openMUSCL(cTMap *h, cTMap *vx, cTMap *vy, cTMap *z);
    double doSWOFMUSCLdt(double dt, double timesum, cTMap *h, cTMap *u, cTMap *v, cTMap *z);

    void doSWOFStV(double dt, cTMap *h, cTMap *u, cTMap *v);

    void ChannelSWOFopen();  //TODO not used
    void KinematicSWOFopen(cTMap *_h, cTMap *_V);
    double limiter(double a, double b);
    vec4 F_ROE(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    vec4 F_HLL4(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    vec4 F_HLL3(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    vec4 F_HLL2(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    vec4 F_HLL(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    vec4 F_Rusanov(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
    vec4 F_Riemann(double h_L,double u_L,double v_L,double h_R,double u_R,double v_R);
   // void dynOutflowPoints(cTMap *h);
    void OverlandFlow2Ddyn(void);
    void updateWHandHmx(void);
    void Boundary2Ddyn(double dt, cTMap *h, cTMap *u, cTMap *v);
    void SWOFDiagonalFlow(double dt_req_min, cTMap *h, cTMap *vx, cTMap *vy);  //OBSOLETE
    void SWOFDiagonalFlowNew(double dt_req_min, cTMap *h, cTMap *vx, cTMap *vy);
    // <= 2D flow

    // <= groundwater
    void GroundwaterRecharge();
    void GroundwaterFlow();
    void GWFlow2D(double factor);
    void GWFlowSWAT();
    void GWFlowLDDKsat();
    double fullSWOF2GW(cTMap *h, cTMap *u, cTMap *v, cTMap *z);
    // => groundwater

    void cell_SlopeStability(int r, int c); // TODO

    // => extend channel, not used for now
    void doExtendRow(int r, int c, int n,  double w2, double adx);
    void doExtendCol(int r, int c, int n,  double w2, double adx);
    void extendRow(int r, int c, int i, double w);
    void extendCol(int r, int c, int i, double w);
    void ExtendChannel();
    bool ExtendChannelNew();
    bool IsExtendedChannel(int r, int c, int dr, int dc);
    void DistributeOverExtendedChannel(cTMap * _In, cTMap * _Out);
    // <= extend channel

    void InitFlowBarriers(void);
    double DEMFB(int r, int c, int rd, int cd, bool addwh);
    double FB(int r, int c, int rd, int cd);
    void SetFlowBarriers();
    void GetFlowBarrierData(QString name);
    double FBW(double h, int r, int c, int dr, int dc);

    double courant_factor;
    double courant_factorSed;
    double mixing_coefficient, runoff_partitioning;
    double minReportFloodHeight;
    // boundary in 2D flow
    double QBoundary;
    double QsBoundary;
    double TimestepfloodMin, TimestepfloodLast;
    QVector <double> Qout;

    // => kinematic, linked lists routing
    void Kinematic(int pitRowNr, int pitColNr, cTMap *_LDD, cTMap *_Q, cTMap *_Qn, cTMap *_Alpha, cTMap *_DX, cTMap *_Qmax, cTMap *_Amax);
    void KinematicExplicit(QVector<LDD_COORIN> _crlinked, cTMap *_Q, cTMap *_Qn, cTMap *_Alpha,cTMap *_DX, cTMap *_Qmax, cTMap *_Amax);
    void routeSubstance(int pitRowNr, int pitColNr, cTMap *_LDD,
                                cTMap *_Q, cTMap *_Qn, cTMap *_Qs, cTMap *_Qsn,
                                cTMap *_Alpha, cTMap *_DX, cTMap*_Sed);//,cTMap*_VolStore, cTMap*_SedStore);
    void KinematicSubstance(QVector<LDD_COORIN> _crlinked_, cTMap *_LDD, cTMap *_Q, cTMap *_Qn, cTMap *_Qs, cTMap *_Qsn, cTMap *_Alpha,cTMap *_DX, cTMap *_Sed);
    double IterateToQnew(double Qin, double Qold, double alpha, double deltaT, double deltaX, double Qm, double Am);
    double simpleSedCalc(double Qj1i1, double Qj1i, double Sj1i, double vol, double sed);
    double complexSedCalc(double Qj1i1, double Qj1i, double Qji1, double Sj1i,double Sji1, double alpha, double dx);
    void upstream(QVector <LDD_COORIN>_crlinked_, cTMap *_Q, cTMap *_Qn);
    void downstream(QVector <LDD_COORIN>_crlinked_, cTMap *_Q, cTMap *_Qn);
    void upstreamMax(QVector <LDD_COORIN>_crlinked_, cTMap *MaxQ, cTMap *Q, cTMap *_Qn);
    void UpstreamAvg(QVector <LDD_COORIN>_crlinked_ , cTMap *_Q, cTMap *_Qn);
    void AccufluxGW(QVector <LDD_COORIN>_crlinked_ , cTMap *_Q, cTMap *_Qn, cTMap *_CW);
    QVector <LDD_COORIN> MakeLinkedList(cTMap *_LDD);
    double itercount;
    // <= kinematic

    // => sediment stuff
    double rillfactor;
    double GetSV(double d);
    void SplashDetachment();
    double MaxConcentration(double watvol, double sedvol);
    void ChannelFlowDetachmentNew();
    void RiverSedimentDiffusion(double dt, cTMap * _SS,cTMap * _SSC);
    void RiverSedimentLayerDepth(int r , int c);
    void RiverSedimentMaxC(int r, int c);
    double calcTCSuspended(int r,int c, int _d, int method, double h, double w,  double U, int type);
    double calcTCBedload(int r,int c, int _d, int method, double h, double w, double U, int type);
    void SWOFSedimentCheckZero(int r, int c, cTMap * h);
    void SWOFSedimentSetConcentration(int r, int c, double h, double w);
    void SWOFSedimentDiffusion(double dt, cTMap * h,cTMap * u,cTMap *v, cTMap * _SS,cTMap * _SSC);
    void SWOFSedimentFlowInterpolation(double dt, cTMap * h, cTMap * u,cTMap * v, cTMap * _SS,cTMap * _SSC);
    void SWOFSedimentDetNew(double dt, cTMap * h, cTMap *w, cTMap * u,cTMap * v);
    void SWOFSediment(double dt, cTMap * h, cTMap *w, cTMap * u,cTMap * v);
    void SWOFSedimentLayerDepth(int r , int c, double h, double velocity);//cTMap * u,cTMap * v);
    void correctMassBalance(double sum1, cTMap *M);
    void correctMassBalanceSed(double sum1, cTMap *M, double th);
    double getMass(cTMap *M);
    double getMassSed(cTMap *M, double th);

    // => SWATRE
    /// filenames for Swatre soil information
    QList <int> ProfileIDList;
    QString SwatreTableDir;
    QString SwatreTableName;
    QString initheadName;
    void InitNewSoilProfile();
    double swatreDT;
    bool initSwatreStructure;
    SOIL_MODEL *SwatreSoilModel;
    SOIL_MODEL *SwatreSoilModelCrust;
    SOIL_MODEL *SwatreSoilModelCompact;
    SOIL_MODEL *SwatreSoilModelGrass;
    PROFILE **profileList = nullptr;
    HORIZON **horizonList = nullptr;
    ZONE *zone = nullptr;
    int tnode; //VJ 110122 node nr in profile with tile drains
    SOIL_MODEL *InitSwatre(cTMap *profileMap);//, QString initHeadMaps, cTMap *tiledepthMap, double dtMin);
    void CloseSwatre(SOIL_MODEL *s);
    void FreeSwatreInfo(void);
    QStringList swatreProfileDef;
    QList<int> swatreProfileNr;
    int sizeProfileList; // to free mem
    int nrHorizonList;
    int sizeHorizonList;
    void ReadSwatreInputNew(void);
    PROFILE *ReadProfileDefinitionNew(int pos, ZONE *z);
    HORIZON *ReadHorizonNew(QString tablePath, QString tableName);
    LUT *ReadSoilTableNew(QString fileName);
    void checkFileForInvalidLetters(const QString &filePath);
    double SwatreStep(long i_, int r, int c, SOIL_MODEL *s, double _WH, cTMap *_drain);
    void HeadCalc(const PROFILE *p, double *h, bool *isPonded, bool fltsat,
                  const double *thetaPrev, const double *hPrev, const double *kavg, const double *dimoca,
                  double dt, double pond, double qtop, double qbot);
    double  NewTimeStep(double prevDt, const double *hLast, const double *h, int nrNodes, double dtMin, double precParam);
    void ComputeForPixel(long i_, SOIL_MODEL *s);
    double DmcNode(double head,const  HORIZON *hor,bool on_dmch);
    double FindValue(double value,const  HORIZON *hor, int colv, int col);
    double HNode(double theta,const  HORIZON *hor); // obsolete
    double TheNode(double head,const  HORIZON *hor);// obsolete
    double HcoNode(double head,const HORIZON *hor); // obsolete
    // <= SWATRE

int showr;// for debugging
int showc;

    void Fill(cTMap &M, double value);
    void Copy(cTMap &M, cTMap &M1);
    double MapTotal(cTMap &M);
    void Average3x3(cTMap &M, cTMap &mask, bool only);
    void Average2x2(cTMap &M, cTMap &mask);

    // => mass balance and output
    void TotalsHydro(void);
    void TotalsFlow(void);
    void TotalsSediment(void);
    void MassBalance(void);
    void reportToUI(void);
    void reportToFile(void);
    void ReportTimeseriesPCR(void);
    void ReportTimeseriesCSV(void);
    void ReportTotalSeries(void);
    void ReportMaps(void);
    void ReportMapSeries(void);
    void ReportTotalsNew(void);
    void ReportErosionLandunits(void); //VJ 110107 report erosion stats per land unit
    void CountLandunits(void); //VJ 110107 report erosion stats per land unit
    void saveMBerror2file(bool start);
    void FloodMaxandTiming();
    void FloodStatistics(void);
    // => mass balance and output

    // thread management variables
    bool stopRequested;
    bool waitRequested;
    bool noInterface;
    bool showInfo;
    bool noOutput;
    bool batchmode;
    QMutex mutex;
    QWaitCondition condition;
    void stop();

protected:
    void run();

    // talk to the interface
    QElapsedTimer time_ms;
    double startTime;
    void setupDisplayMaps();
    void setupHydrographData();


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
