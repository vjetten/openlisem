/*---------------------------------------------------------------------------
project: openLISEM
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website, information and code: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/

/*
Functionality in model.h:
- TWorld class that combines all model variables and processes
- global defines Drc, MV, FOR_ROW_COL_MV etc
- global defines for lisem type; infiltration type etc.
 */
#ifndef modelH
#define modelH

#include <math.h>
#include <stdlib.h>

#include <QtGui>
#include <QMutex.h>

#include "csfmap.h"
#include "mmath.h"
#include "error.h"
#include "swatre_g.h"

#define DEBUGv(x) {QString sss; sss.setNum(x);emit debug("debug: " + sss);}
//msleep(300)
#define DEBUG(s) emit debug(QString(s))

#define mwrite(name) WriteMap(QString(resultDir+name))
#define report(name) WriteMapSeries(resultDir,QString(name), printstep)

#define Drc     Data[r][c]
#define DrcOutlet     Data[r_outlet][c_outlet]
#define MV(r,c) IS_MV_REAL4(&Mask->Data[r][c])
#define FOR_ROW_COL_MV for (int r = 0; r < nrRows; r++)\
		for (int c = 0; c < nrCols; c++)\
		if(!IS_MV_REAL4(&Mask->Data[r][c]))

#define FOR_ROW_COL_MV_CH for (int  r = 0; r < nrRows; r++)\
		for (int  c = 0; c < nrCols; c++)\
		if(!IS_MV_REAL4(& ChannelMask->Data[r][c]))

#define NUMNAMES 300
#define NUMMAPS 110
#define MIN_FLUX 1e-12
#define MIN_HEIGHT 1e-6
#define MAXCONC 848
//0.32 * 2650 = max vol conc * bulk density

#define INFIL_NONE 0
#define INFIL_SWATRE 1
#define INFIL_HOLTAN 2
#define INFIL_GREENAMPT 3
#define INFIL_GREENAMPT2 4
#define INFIL_KSAT 5
#define INFIL_MOREL 21
#define INFIL_SMITH 22
#define INFIL_SMITH2 23


//---------------------------------------------------------------------------
// structure containing pointers to all maps for automatic destruction after runs
// so memory doesn't have to be freed for each map
typedef struct MapListStruct {
	TMMap *m;
}  MapListStruct;
//---------------------------------------------------------------------------
// linked list structure for network in kin wave
typedef struct Liststruct {
	int rowNr;
	int colNr;
	struct Liststruct *prev;
}  Liststruct;
//---------------------------------------------------------------------------
// name list structure used to read run file
typedef struct _nameList{
	QString name;
	QString value;
} _nameList;
//---------------------------------------------------------------------------
// map name list structure for interaction with interface
typedef struct _mapList{
   QString name;
   QString dir;
   QString id;
   int groupnr;
   int varnr;
} _mapList;
//---------------------------------------------------------------------------
typedef struct _unitList{
   long nr;
   double area;
   double totdet;
   double totdep;
   double totsl;
} _unitList;
//---------------------------------------------------------------------------


// The world: main class containing all variables, maps, options, filenames etc etc
class TWorld: public QThread
{
	Q_OBJECT

public:
	TWorld(QObject *parent = 0);
	~TWorld();

	// copy of overall rows and columns, set in initmask
	long nrRows, nrCols;

   // map management structure, automatic adding and deleting of all TMMap variables
   MapListStruct maplistTMMap[NUMNAMES];
	int maplistnr;

	// All maps are declared here, no lacal declarations of maps are done
	TMMap
	// terrain
	*tm, *tma, *Mask, *MaskChannel, *DEM, *DX, *CellArea, *Grad, *LDD, *Outlet,*PointMap,
	// rainfall and interception
	*RainZone, *Rain, *Rainc, *RainCum, *RainNet, *LeafDrain, *RainIntensity, *RainM3, *CStor, *Interc,
	*SnowmeltZone, *Snowcover, *Snowmelt, *Snowmeltc, *SnowmeltCum, *SnowmeltNet,
	// runoff
	*WH, *WHroad, *WHrunoff, *WHrunoffCum, *WHstore, *WaterVolrunoff, *WaterVolin, *WaterVolall,
	*FlowWidth, *V, *Alpha, *Q, *Qoutflow, *Qn, *Qoutput, *Qs, *Qsn, *Qsoutput, *q, *R, *Perim,
	// soil surface, storage and cover
	*N, *RR, *MDS,*fpa,
	*SoilWidthDX, *RoadWidthDX, *StoneFraction, *CompactFraction, *CrustFraction,
	*PlantHeight, *Cover, *CanopyStorage, *LAI,
   *LandUnit, //VJ 110107 land unit map for erosion output
	// not used yet:
	*WheelWidth, *WheelWidthDX, *GullyWidthDX,
	//erosion
	*Cohesion, *RootCohesion, *CohesionSoil, *Y, *AggrStab, *D50,
	*DETSplash, *DETFlow, *HardSurface,*DEP, *TC, *Conc, *Sed, *Qsoutflow, *CG, *DG, *SettlingVelocity,
	//infiltration
   *Fcum, *FSurplus, *FFull, *fact, *fpot, *InfilVolKinWave, *InfilVol, *InfilVolCum,
	*ThetaS1, *ThetaI1, *Psi1, *Ksat1, *SoilDepth1, *L1, *Soilwater,
	*ThetaS2, *ThetaI2, *Psi2, *Ksat2, *SoilDepth2, *L2, *Soilwater2,
	*KsatCrust, *KsatCompact, *KsatGrass, *Ksateff, *L1gr, *L2gr, *factgr, *fpotgr,
	*WHGrass, *Fcumgr, *GrassFraction, *GrassWidthDX, *GrassPresent,
	*ProfileID, *ProfileIDCrust, *ProfileIDCompact, *ProfileIDGrass,
	// Channels
	*RunoffVolinToChannel, *LDDChannel, *ChannelWidth, *ChannelSide, *ChannelQ, *ChannelQn, *ChannelQs, *ChannelQsn,
	*ChannelQoutflow, *ChannelGrad, *ChannelV, *ChannelN, *ChannelWH, *ChannelWaterVol, *Channelq,
   *ChannelAlpha, *ChannelWidthUpDX, *ChannelPerimeter, *ChannelMask, *ChannelDX, *ChannelDetFlow, *ChannelDep, *ChannelKsat,
	*ChannelSed, *ChannelConc, *ChannelTC, *SedToChannel, *ChannelQsoutflow, *ChannelCohesion, *ChannelY,
	// buffers
	*BufferID, *BufferVol, *BufferSed, *ChannelBufferSed, *ChannelBufferVol,
	*BufferVolInit, *BufferSedInit, *ChannelBufferSedInit, *ChannelBufferVolInit,
	//      *BufferSedVolMax, *ChannelBufferSedVolMax,
	// totals for output
	*TotalDetMap, *TotalDepMap, *TotalSoillossMap, *TotalSed, *TotalWatervol, *TotalConc;

	// boolean options that are set in interface and runfile, initialized in DataInit
	bool SwitchHardsurface, SwatreInitialized, SwitchInfilGA2, SwitchCrustPresent,
	SwitchWheelPresent, SwitchCompactPresent, SwitchIncludeChannel, SwitchChannelBaseflow,
	startbaseflowincrease, SwitchChannelInfil, SwitchAllinChannel, SwitchErosion, SwitchAltErosion,
	SwitchSimpleDepression, SwitchBuffers, SwitchSedtrap, SwitchSnowmelt, SwitchRainfall, SwitchRunoffPerM, SwitchInfilCompact,
	SwitchInfilCrust, SwitchInfilGrass, SwitchImpermeable, SwitchDumphead,
	SwitchWheelAsChannel, SwitchMulticlass, SwitchNutrients, SwitchGullies, SwitchGullyEqualWD, SwitchGullyInfil,
	SwitchGullyInit, SwitchOutputTimeStep, SwitchOutputTimeUser, SwitchMapoutRunoff, SwitchMapoutConc,
	SwitchMapoutWH, SwitchMapoutWHC, SwitchMapoutTC, SwitchMapoutEros, SwitchMapoutDepo, SwitchMapoutV,
	SwitchMapoutInf, SwitchMapoutSs, SwitchMapoutChvol, SwitchWritePCRnames, SwitchWritePCRtimeplot,
	SwitchNoErosionOutlet, SwitchDrainage, SwitchPestout, SwitchSeparateOutput, SwitchSOBEKOutput,
	SwitchInterceptionLAI, SwitchTwoLayer, SwitchSimpleSedKinWave, SwitchSoilwater, SwitchSOBEKoutput,
	SwitchPCRoutput, SwitchWriteHeaders, SwitchGeometric;

	// multiple options that are set in interface or runfile, see defines above
	int InterceptionLAIType;
	int InfilMethod;
   int ErosionUnits;

	// calibration factors
	double ksatCalibration;
	double nCalibration;
	double ChnCalibration;
	double ChKsatCalibration;
	double SplashDelivery;
	double StripN;
	double StemflowFraction;

	// totals for mass balance checks and output
	double MB, Qtot, QtotOutlet, IntercTot, WaterVolTot, InfilTot, RainTot, SnowTot, SurfStoremm, InfilKWTot;
	double MBs, DetTot, DetSplashTot, DetFlowTot, DepTot, SoilLossTot, SoilLossTotOutlet, SedTot;
	double ChannelVolTot, ChannelSedTot, ChannelDepTot, ChannelDetTot;
	double RainTotmm, SnowTotmm, IntercTotmm, WaterVolTotmm, InfilTotmm, Qtotmm, RainAvgmm, SnowAvgmm;
	double RainpeakTime, SnowpeakTime, QpeakTime, Qpeak, Rainpeak, Snowpeak;
	double BufferVolTot, BufferSedTot, BufferVolTotInit, BufferSedTotInit, BulkDens;
	double nrCells, CatchmentArea;

	int c_outlet;  // outlet row and col
	int r_outlet;

	//SWATRE
	SOIL_MODEL *SwatreSoilModel;
	SOIL_MODEL *SwatreSoilModelCrust;
	SOIL_MODEL *SwatreSoilModelCompact;
	SOIL_MODEL *SwatreSoilModelGrass;
	double swatreDT;
   bool initSwatreStructure;

	// time and dx parameters
	double time, BeginTime, EndTime;
	double _dt, _dx;
	long runstep, printstep, printinterval;

	// timeseries variables and output strings
	double **RainfallSeries;
	int nrstations, nrrainfallseries, placerainfallseries;
	double **SnowmeltSeries;
	int nrSnowmeltstations, nrSnowmeltseries, placeSnowmeltseries;
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
	QString SwatreTableDir;
	QString SwatreTableName;
	QString initheadName;
	QString resultFileName;
	QString temprunname;
	QStringList outputcheck;
	QString Outrunoff, Outconc, Outwh, Outrwh, Outtc, Outeros, Outdepo, Outvelo, Outinf, Outss, Outchvol;

   _unitList unitList[512]; // just a fixed number for 1024 classes who cares!
   int landUnitNr;

	// data initialization, runfile reading and parsing
   _nameList runnamelist[NUMNAMES]; // structire for runfile variables and names
   int nrrunnamelist;
	void IntializeOptions(void);  // set all options to false etc
	void IntializeData(void);     // make all non-input maps
	void GetRainfallData(void);   // get input timeseries
	void GetSnowmeltData(void);   // get input timeseries
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
	//QString getvaluename(const char *vname);
	QString getvaluename(QString vname);
	double getvaluedouble(QString vname);
	int getvalueint(QString vname);
	QString CheckDir(QString p);

	// LISEM model processes
	void RainfallMap(void);
	void SnowmeltMap(void);
	void Interception(void);
	void Infiltration(void);
	void InfilSwatre(void);
	void InfilGreenAmpt1(void);
	void InfilSmithParlange1(void);
	void InfilMorelSeytoux1(void);
	void InfilKsat(void);
	void SoilWater(void);
	double IncreaseInfiltrationDepth(int r, int c, double fact, REAL8 *L1p, REAL8 *L2p);
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
			TMMap *_Qsn, TMMap *_q, TMMap *_Alpha, TMMap *_DX, TMMap *Vol, TMMap*SedVol,
			TMMap *_StorVol, TMMap*_StorVolSed);
	void Totals(void);
	void MassBalance(void);
	void OutputUI(void);
	//void ReportTimeseries();
	void ReportTimeseriesNew();
	void ReportTotals();
	void ReportMaps();
	void ReportTotalsNew();
   void ReportLandunits(); //VJ 110107 report erosion stats per land unit
   void CountLandunits(); //VJ 110107 report erosion stats per land unit

	// thread management variables
	bool stopRequested;
	bool waitRequested;
	QMutex mutex;
	QWaitCondition condition;
	void stop();

//	int ReadSwatreInput(const char *fileName,	const char *tablePath);
	int ReadSwatreInput(QString fileName, QString tablePath);

	SOIL_MODEL *InitSwatre(
			TMMap *profileMap, // r- profile id map
			QString initHeadMaps,    // init head maps path
			double dtMin,   	      		// minumum timestep, is also initial timestep
			double precis,                // precision factor to adapt timestep
			double calibration,
			bool geom,
			bool bottomClosed);
			//double curtime);

	void SwatreStep(SOIL_MODEL *s, TMMap *_WH, TMMap *_fpot, TMMap *where);
	void CloseSwatre(SOIL_MODEL *s);

protected:
	void run();
	QTime time_ms;
	void DEBUGs(QString SSS);

	// talk to the interface with the output structure "op" declared in ifacebasic
	signals:
	void done(const QString &results);
	void debug(const QString &results);
	void show();

private slots:
void DoModel();

};

#endif
