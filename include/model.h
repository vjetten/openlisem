/*---------------------------------------------------------------------------
project: openLISEM
name: model.h
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website SVN: http://sourceforge.net/projects/lisem

Functionality in model.h:
- TWorld class that combines all model variables and processes
- global defines Drc, MV, FOR_ROW_COL_MV etc
- global defines for lisem type; infiltration type etc.
---------------------------------------------------------------------------*/

#ifndef modelH
#define modelH

#include <math.h>
#include <stdlib.h>

#include <QtGui>
#include <QMutex.h>

#include "csfmap.h"
#include "mmath.h"
#include "error.h"

#define DEBUGv(x) {QString sss; sss.setNum(x);emit debug("debug: " + sss);}
//msleep(300)
#define DEBUG(s) emit debug(QString("debug: "+s));msleep(10)

#define mwrite(name) WriteMap(QString(resultDir+name))
#define report(name) WriteMapSeries(resultDir,QString(name), printstep+1)

#define Drc     Data[r][c]
#define MV(r,c) IS_MV_REAL4(&Mask->Data[r][c])
#define FOR_ROW_COL_MV for (int r = 0; r < nrRows; r++)\
                            for (int c = 0; c < nrCols; c++)\
                              if(!IS_MV_REAL4(&Mask->Data[r][c]))

#define FOR_ROW_COL_MV_CH for (int r = 0; r < nrRows; r++)\
                            for (int c = 0; c < nrCols; c++)\
                              if(!IS_MV_REAL4(& ChannelMask->Data[r][c]))

#define NUMNAMES 300
#define MIN_FLUX 1e-12
#define MIN_HEIGHT 1e-6
#define MAXCONC 848
//0.32 * 2650 = max vol conc * bulk density

#define LISEMBASIC 0
#define LISEMWHEELTRACKS 1
#define LISEMMULTICLASS 2
#define LISEMNUTRIENTS 3
#define LISEMGULLIES 4

#define INFIL_NONE 0
#define INFIL_SWATRE 1
#define INFIL_HOLTAN 2
#define INFIL_GREENAMPT 3
#define INFIL_GREENAMPT2 4
#define INFIL_KSAT 5
#define INFIL_MOREL 21
#define INFIL_SMITH 22


//---------------------------------------------------------------------------
typedef struct MapListStruct {
   TMMap *m;
}  MapListStruct;
//---------------------------------------------------------------------------
typedef struct Liststruct {
    int rowNr;
    int colNr;
    struct Liststruct *prev;
}  Liststruct;
//---------------------------------------------------------------------------
typedef struct _nameList{
     QString name;
     QString value;
} _nameList;
//---------------------------------------------------------------------------
class TWorld: public QThread
{
	Q_OBJECT

public:
      TWorld(QObject *parent = 0);
     ~TWorld();

     // copy of overall rows and columns, set in initmask
     long nrRows, nrCols;

     // map management structure
     MapListStruct maplist[NUMNAMES];
     int maplistnr;

     // All maps are declared here, no lacal declarations
     TMMap *tm, *Mask, *MaskChannel, *DEM, *DX, *Grad, *LDD, *Outlet, *RainZone, *N, *RR, *MDS,
      *Rain, *RainCum, *RainNet, *LeafDrain, *RainIntensity, *RainM3, *CStor, *Interc,
      *WH, *WHinf, *WHroad, *WHrunoff, *WHstore, *WaterVolrunoff, *WaterVolin, *WaterVolall, *InfilVolKinWave, *InfilVol, *fpa,
      *FlowWidth, *V, *Alpha, *Q, *Qoutflow, *Qn, *Qoutput, *Qs, *Qsn, *q, *R, *Perim, *WheelWidthDX, *StoneWidthDX,
      *SoilWidthDX, *GullyWidthDX, *RoadWidthDX, *WheelWidth, *StoneFraction, *CompactFraction, *CrustFraction,
      *PlantHeight, *Cover, *CanopyStorage, *LAI,
      *Cohesion, *RootCohesion, *CohesionSoil, *Y, *AggrStab, *D50, *DETSplash, *DETFlow, *HardSurface,
      *DEP, *TC, *Conc, *SedVol, *Qsoutflow, *CG, *DG, *SettlingVelocity, *Fcum, *FSurplus, *fact, *fpot,
      *ThetaS1, *ThetaI1, *Psi1, *Ksat1, *SoilDepth1, *L1, *Soilwater,
      *ThetaS2, *ThetaI2, *Psi2, *Ksat2, *SoilDepth2, *L2, *Soilwater2,
      *KsatCrust, *KsatCompact, *KsatGrass, *Ksateff, *L1gr, *L2gr, *factgr, *fpotgr,
      *WHGrass, *Fcumgr, *GrassFraction, *GrassWidthDX, *GrassPresent,
      *RunoffVolinToChannel, *LDDChannel, *ChannelWidth, *ChannelSide, *ChannelQ, *ChannelQn, *ChannelQs, *ChannelQsn,
      *ChannelQoutflow, *ChannelGrad, *ChannelV, *ChannelN, *ChannelWH, *ChannelWaterVol, *Channelq,
      *ChannelAlpha, *ChannelWidthUpDX, *ChannelMask, *ChannelDX, *ChannelDetFlow, *ChannelDep, *ChannelKsat,
      *ChannelSedVol, *ChannelConc, *ChannelTC, *SedToChannel, *ChannelQsoutflow, *ChannelCohesion, *ChannelY,
      *PointMap, *TotalDetMap, *TotalDepMap, *TotalSoillossMap;

     // boolean options that are set in interface and runfile, initialized in DataInit
     bool SwitchHardsurface, SwatreInitialized, SwitchInfilGA2, SwitchCrustPresent,
      SwitchWheelPresent, SwitchCompactPresent, SwitchIncludeChannel, SwitchChannelBaseflow,
      startbaseflowincrease, SwitchChannelInfil, SwitchAllinChannel, SwitchErosion, SwitchAltErosion,
      SwitchSimpleDepression, SwitchBuffers, SwitchSedtrap, SwitchSnowmelt, SwitchRunoffPerM, SwitchInfilCompact,
      SwitchInfilCrust, SwitchInfilGrass, SwitchImpermeable, SwitchDumphead, SwitchGeometricMean,
      SwitchWheelAsChannel, SwitchMulticlass, SwitchNutrients, SwitchGullies, SwitchGullyEqualWD, SwitchGullyInfil,
      SwitchGullyInit, SwitchOutputTimeStep, SwitchOutputTimeUser, SwitchMapoutRunoff, SwitchMapoutConc,
      SwitchMapoutWH, SwitchMapoutWHC, SwitchMapoutTC, SwitchMapoutEros, SwitchMapoutDepo, SwitchMapoutV,
      SwitchMapoutInf, SwitchMapoutSs, SwitchMapoutChvol, SwitchWritePCRnames, SwitchWritePCRtimeplot,
      SwitchNoErosionOutlet, SwitchDrainage, SwitchPestout, SwitchSeparateOutput, SwitchSOBEKOutput,
      SwitchInterceptionLAI, SwitchTwoLayer, SwitchSimpleSedKinWave, SwitchSoilwater, SwitchSOBEKoutput,
      SwitchPCRoutput;

     // multiple options that are set in interface or runfile, see defines above
    int InterceptionLAIType;
    int InfilMethod;

    // calibration factors
    double ksatCalibration;
    double nCalibration;
    double ChnCalibration;
    double ChKsatCalibration;
    double SplashDelivery;
    double StripN;
    double StemflowFraction;

    // totals for mass balance checks and output
    double MB, Qtot, IntercTot, WaterVolTot, InfilTot, RainTot, SurfStorTot, InfilKWTot;
    double MBs, DetTot, DetTotSplash, DetTotFlow, DepTot, SoilLossTot, SedVolTot;
    double ChannelVolTot, ChannelSedTot, ChannelDepTot, ChannelDetTot;
    double RainTotmm, IntercTotmm, WaterVolTotmm, InfilTotmm, Qtotmm, Qpeak, Rainpeak;
    double nrCells, CatchmentArea, RainpeakTime, QpeakTime;

    // time and dx parameters
    double time, BeginTime, EndTime;
    double _dt, _dx;
    long runstep, printstep;

    // timeseries variables and output strings
    double **RainfallSeries;
    int nrstations, nrrainfallseries, placerainfallseries;
    QString SOBEKdatestring;
    int SOBEKnrlines;
    char ErosionUnits;

    // file and directory names
    QString resultDir;
    QString inputDir;
    QString outflowFileName;
    QString totalErosionFileName;
    QString totalDepositionFileName;
    QString totalSoillossFileName;
    QString rainFileName;
    QString rainFileDir;
    QString snowmeltFileName;
    QString snowmeltFileDir;
    QString tableDir;
    QString resultFileName;
    QString temprunname;
    QStringList outputcheck;
    QString Outrunoff, Outconc, Outwh, Outrwh, Outtc, Outeros, Outdepo, Outvelo, Outinf, Outss, Outchvol;
    QString baseNameMap;

    // data initialization, runfile reading and parsing
    _nameList namelist[NUMNAMES]; // structire for runfile variables and names
    int nrnamelist;
    void IntializeOptions(void);  // set all options to false etc
    void IntializeData(void);     // make all non-input maps
    void GetRainfallData(void);   // get input timeseries
    void GetInputData(void);      // get and make input maps
    TMMap *InitMask(QString name);
    TMMap *InitMaskChannel(QString name);
    TMMap *ReadMapMask(QString name);
    TMMap *ReadMap(cTMap *Mask, QString name);
    TMMap *NewMap(double value);
    void InitMapList(void);
    void DestroyData(void);
    void ParseInputData();
    void GetRunFile();
    QString GetName(QString p);
    QString getvaluename(const char *vname);
    double getvaluedouble(const char *vname);
    int getvalueint(const char *vname);
    QString CheckDir(QString p, QString p1);

    // LISEM model processes
    void Rainfall(void);
    void Interception(void);
    void Infiltration(void);
    void InfilGreenAmpt1(void);
    void InfilSmithParlange1(void);
    void InfilMorelSeytoux1(void);
    void InfilKsat(void);
    void SoilWater(void);
    double IncreaseInfiltrationDepth(int r, int c, double fact, REAL4 *L1p, REAL4 *L2p);
    void SurfaceStorage(void);
    void OverlandFlow(void);
    void ChannelFlow(void);
    void ToChannel(void);
    void CalcVelDisch(void);
    void CalcVelDischChannel();
    void GridCell(void);
    void SplashDetachment(void);
    void FlowDetachment(void);
    double MaxConcentration(double watvol, double sedvol, double dep);
    void ChannelFlowDetachment(void);
    void Kinematic(int pitRowNr, int pitColNr, TMMap *_LDD, TMMap *_Q, TMMap *_Qn, TMMap *_Qs,
                   TMMap *_Qsn, TMMap *_q, TMMap *_Alpha, TMMap *_DX, TMMap *Vol, TMMap*SedVol);
    void MassBalance(void);
    void Output(void);
    //void ReportTimeseries();
    void ReportTimeseriesNew();
    void ReportTotals();
    void ReportMaps();
    void ReportTotalsNew();

    // thread management variables
    bool stopRequested;
    QMutex mutex;
    void stop();

protected:
    void run();
    QTime time_ms;

// talk to the interface with the output structure "op" declared in ifacebasic
signals:
    void done(const QString &results);
    void debug(const QString &results);
    void show();

private slots:
    void DoModel();

};

#endif
