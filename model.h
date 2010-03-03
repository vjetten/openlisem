//---------------------------------------------------------------------------

#ifndef modelH
#define modelH

#include <math.h>
#include <stdlib.h>

#include <QtGui>
#include <QMutex.h>

#include "csfmap.h"
#include "mmath.h"
#include "error.h"

#define DEBUGv(x) emit debug("debug: "+QString::setNum(x)+"<=");msleep(300)
#define DEBUG(s) emit debug(s);msleep(300)

#define mwrite(name) WriteMap(QString(resultdir+name))
#define report(name, step) WriteMapSeries(QString(resultDir+name), step)

#define    Drc     Data[r][c]
#define    MV(r,c) IS_MV_REAL4(&Mask->Data[r][c])
#define    FOR_ROW_COL_MV for (int r = 0; r < nrRows; r++)\
                            for (int c = 0; c < nrCols; c++)\
                              if(!IS_MV_REAL4(&Mask->Data[r][c]))

#define    FOR_ROW_COL_MV_CH for (int r = 0; r < nrRows; r++)\
                            for (int c = 0; c < nrCols; c++)\
                              if(!IS_MV_REAL4(& ChannelMask->Data[r][c]))

#define NUMNAMES 300

#define MIN_FLUX 1e-12
#define MIN_HEIGHT 1e-6
#define MAXCONC 850

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

     TMMap *tm, *Mask, *DEM, *DX, *Grad, *LDD, *Outlet, *RainZone, *N, *RR, *MDS,
      *Rain, *RainCum, *RainNet, *RainIntensity, *RainM3, *CStor, *Interc,
      *WH, *WHinf, *WHroad, *WHrunoff, *WHstore, *WaterVol, *WaterVolRunoff, *InfilVolKinWave, *InfilVol, *fpa,
      *FlowWidth, *V, *Alpha, *Q, *Qoutflow, *Qn, *Qs, *Qsn, *q, *R, *Perim, *WheelWidthDX, *StoneWidthDX,
      *SoilWidthDX, *GullyWidthDX, *RoadWidthDX, *WheelWidth, *StoneFraction, *CompactFraction, *CrustFraction,
      *PlantHeight, *Cover, *CanopyStorage, *LAI,
      *Cohesion, *RootCohesion, *CohesionSoil, *Y, *AggrStab, *D50, *DETSplash, *DETFlow,
      *DEP, *TC, *Conc, *SedVol, *Qsoutflow, *CG, *DG, *SettlingVelocity, *Fcum, *FSurplus, *fact, *fpot,
      *ThetaS1, *ThetaI1, *Psi1, *Ksat1, *SoilDepth1, *L1, *Soilwater,
      *ThetaS2, *ThetaI2, *Psi2, *Ksat2, *SoilDepth2, *L2, *Soilwater2,
      *KsatCrust, *KsatCompact,
      *KsatGrass, *Ksateff, *L1gr, *L2gr, *factgr, *fpotgr, *WHGrass, *Fcumgr, *GrassFraction, *GrassWidthDX,
      *RunoffVolinToChannel, *LDDChannel, *ChannelWidth, *ChannelSide, *ChannelQ, *ChannelQn, *ChannelQs, *ChannelQsn,
      *ChannelQoutflow, *ChannelGrad, *ChannelV, *ChannelN, *ChannelWH, *ChannelWaterVol, *Channelq,
      *ChannelAlpha, *ChannelWidthUpDX, *ChannelMask, *ChannelDX, *ChannelDetFlow, *ChannelDep,
      *ChannelSedVol, *ChannelConc, *ChannelTC, *SedToChannel, *ChannelQsoutflow, *ChannelCohesion, *ChannelY,
      *PointMap;

    bool SwitchHardsurface, SwatreInitialized, SwitchInfilGA2, SwitchCrustPresent, SwitchGrassPresent,
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

    QString SOBEKdatestring;
    char ErosionUnits;
    int SOBEKnrlines;
    int InterceptionLAIType;
    int InfilMethod;

    // calibration factors
    double ksatCalibration;
    double nCalibration;
    double ChnCalibration;
    double ChKsatCalibration;
    double SplashDelivery;

    // totals for mass balance checks
    double MB, Qtot, IntercTot, WaterVolTot, InfilTot, RainTot, SurfStorTot, InfilKWTot;
    double MBs, DetTot, DetTotSplash, DetTotFlow, DepTot, SoilLossTot, SedVolTot;
    double ChannelVolTot, ChannelSedTot, ChannelDepTot;

    double time, BeginTime, EndTime;
    double _dt, _dx;
    long nrRows, nrCols;
    long runstep;

    double **RainfallSeries;
    int nrstations, nrrainfallseries, placerainfallseries;


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

    // data initialization, runfile reading and parsing
    void IntializeOptions(void);
    void GetRainfallData(void);
    void IntializeData(void);
    void GetInputData(void);
    void InitMask(cTMap *M);
    TMMap *ReadMap(QString name);
    TMMap *NewMap(double value);
    void InitMapList(void);
    void DestroyData(void);
    void ParseInputData();
    void GetRunFile();
    QString getvaluename(const char *vname);
    double getvaluedouble(const char *vname);
    int getvalueint(const char *vname);

    // process functions
    void Rainfall(void);
    void Interception(void);
    void Infiltration(void);
    void InfilGreenAmpt1(void);
    void InfilSmithParlange1(void);
    void InfilMorelSeytoux1(void);
    void InfilKsat(void);
    void SoilWater(void);
    double IncreaseInfiltrationDepth(int r, int c, double fact, double L1, double L2);
    void SurfaceStorage(void);
    void OverlandFlow(void);
    void ChannelFlow(void);
    void ToChannel(void);
    void CalcVelDisch(void);
    void CalcVelDischChannel();
    void GridCell(void);
    void SplashDetachment(void);
    void FlowDetachment(void);
    double MaxConcentration(int r, int c);
    double MaxChannelConcentration(int r, int c);
    void ChannelFlowDetachment(void);
    void Kinematic(int pitRowNr, int pitColNr, TMMap *_LDD, TMMap *_Q, TMMap *_Qn, TMMap *_Qs,
                   TMMap *_Qsn, TMMap *_q, TMMap *_Alpha, TMMap *_DX, TMMap *Vol, TMMap*SedVol);
    void MassBalance(void);
    void Output(void);
    void ReportMap();

    MapListStruct maplist[NUMNAMES];
    int maplistnr;

    _nameList namelist[NUMNAMES];
    int nrnamelist;

    // thread management variables
    bool stopRequested;
    QMutex mutex;
    void stop();

protected:
    void run();
    QTime time_ms;

signals:
    void done(const QString &results);
    void debug(const QString &results);
    void show(const int i);

private slots:
    void DoModel();

};

#endif
