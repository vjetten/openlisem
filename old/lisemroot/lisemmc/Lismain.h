//---------------------------------------------------------------------------
#ifndef LisnewH
#define LisnewH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include "mprolog.h"   // <=== links to all header files, cps.h etc
//---------------------------------------------------------------------------

#define INFIL_NONE 0
#define INFIL_SWATRE 1
#define INFIL_HOLTAN 2
#define INFIL_GREENAMPT 3
#define INFIL_GREENAMPT2 4
#define INFIL_KSAT 5 //7
//VJ 080521
#define INFIL_MOREL 21
#define INFIL_SMITH 22

//========comma delimited output===========
#define commaoutputstr "\"t\",\"P\",\"Q\",\"Qsed\",\"Conc\"\n"\
                       "\"min\",\"mm/h\",\"l/s\",\"kg/s\",\"g/l\"\n"
#define commaoutputstrMC "\"t\",\"P\",\"Q\",\"Qsed\",\"Conc\"\
                         ,\"mu0\",\"mu1\",\"mu2\",\"mu3\",\"mu4\",\"mu5\"\n"\
                         "\"min\",\"mm/h\",\"l/s\",\"kg/s\",\"g/l\",\
                          \"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\"\n"
#define commaoutputstrNUT "\"t\",\"P\",\"Q\",\"Qsed\",\"Conc\"\
                          ,\"mu0\",\"mu1\",\"mu2\",\"mu3\",\"mu4\",\"mu5\"\
                          ,\"Psol\",\"NH4sol\",\"NO3sol\",\"Psus\",\"NH4sus\",\"NO3sus\"\n"\
                          "\"min\",\"mm/h\",\"l/s\",\"kg/s\",\"g/l\"\
                          ,\"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\"\
                          ,\"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\",\"g/l\"\n"

#define commaoutputstrBuffer "\"t\",\"P\",\"Water Vol\",\"Sed Vol\"\n" \
                             "\"min\", \"mm/h\", \"m3\", \"kg\"\n"

//VJ 050822 changed rainfall precision to 4 digits
#define commaformat4      "%.6g,%.4g,%.6g,%.6g\n"
#define commaformat       "%.6g,%.4g,%.6g,%.6g,%.6g\n"
#define commaformatMC     "%.6g,%.4g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n"
#define commaformatNUT    "%.6g,%.4g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n"

//========space delimited output===========
#define ncommaoutputstr "t   P    Q   Qsed Conc \n"     \
                        "min mm/h l/s kg/s g/l \n"
#define ncommaoutputstrMC "t   P     Q     Qsed   Conc  mu0   mu1   mu2   mu3   mu4   mu5\n"     \
                          "min mm/h  l/s   kg/s   g/l   g/l   g/l   g/l   g/l   g/l   g/l\n"
#define ncommaoutputstrNUT "t   P    Q   Qsed Conc  mu0   mu1   mu2   mu3   mu4   mu5  Psol  NH4sol  NO3sol Psus  NH4sus  NO3sus\n"\
                           "min mm/h l/s kg/s g/l   g/l   g/l   g/l   g/l   g/l   g/l  g/l   g/l     g/l    g/l   g/l     g/l\n"

#define ncommaoutputstrBuffer "t    P     Water Vol  Sed Vol\n" \
                              "min  mm/h  m3         kg\n"

//VJ 050822 changed rainfall precision to 4 digits
#define ncommaformat4      "%.6g %.4g %.6g %.6g\n"
#define ncommaformat       "%.6g %.4g %.6g %.6g %.6g\n"
#define ncommaformatMC     "%.6g %.4g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n"
#define ncommaformatNUT    "%.6g %.4g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n"

//========timeplot headers
#define timeplotoutputstr "LISEM output\n 5 \n"\
                       "t (min)\nP (mm/h)\nQ (l/s)\nQsed (kg/s)\nConc (g/l)\n"
#define timeplotoutputstrMC "LISEM output\n 11 \n"\
                       "t (min)\nP (mm/h)\nQ (l/s)\nQsed (kg/s)\nConc (g/l)\n"\
                       "mu0 (g/l)\nmu1 (g/l)\nmu2 (g/l)\nmu3 (g/l)\nmu4 (g/l)\nmu5 (g/l)\n"
#define timeplotoutputstrNUT "LISEM output\n 11 \n"\
                          "t (min)\nP (mm/h)\nQ (l/s)\nQsed (kg/s)\nConc (g/l)\n"\
                          "Psol (g/l)\nNH4sol (g/l)\nNO3sol (g/l)\nPsus (g/l)\nNH4sus (g/l)\nNO3sus (g/l)\n"
#define timeplotoutputstrBuffer "LISEM output\n 4 \n"\
                          "t (min)\nP (mm/h)\nWater Vol(m3)\nSed Vol (kg)\n"


//#define WRITE_END_RESULT       (writeEachTime || lastStep)
// obsolete

#define ReportV(S)    ReportText = S;

#define LisemWarning    LisemWarn

#define writeList(v, s) if (PERIODMAPS && PERIODCOUNT == PERIOD){\
                        strcpy(FileName, CatPath(s, RESPATH));\
                        write(v, MakeFileName(FileName,timestepindex));}

#define writepath(v, s) strcpy(FileName, CatPath(s, RESPATH));write(v, FileName)

#define writeTimeseries(v, s) if (SwitchWritePCRnames) write(v, MakePCRFileName(s,timeIndex+1));\
                              else write(v, MakeFileName(s,timestepindex))

//#define mapname(s) LisIFace->Lmapname(s)
#define mapname(s) Lmapname(s, 1)
#define mapnameonly(s) Lmapname(s, 0)

bool swatreBottomClosed = false;//(bool)SwitchImpermeable;


typedef struct {
   char name[64];
   char value[256];
   int type;
   int iii;
   float vvv;
   bool done;
} nameList;


class LisThread : public TThread
{
private:
protected:
    void __fastcall Execute();
    AnsiString WarningText;
    AnsiString ReportText;
public:

     int StepCounter;
     int errorlevel;
     bool InitDone;
     bool RunDone;
     nameList *namelist;
     int nrnamelist;
     char temprunname[512];
     int runlen;

     char rainFileName[128];
//VJ 080423 add snowmelt
     char snowmeltFileName[128];
     char outflowFileName[128];
   //  char outflowFileName2[128];
   //  char outflowFileName3[128];
   //  char outPointFileName[128];

//VJ 040823 add buffer output
     char BufferFileName[128];
//VJ 050913 add pest output
     char PestoutFileName[128];

//VJ 100115 add total runoff     
//     char totalRunoffFileName[128];
     char totalErosionFileName[128];
     char totalDepositionFileName[128];
     char totalSoillossFileName[128];
     char resultFileName[128];
//     char periodFileName[100];
     char tableDir[128];
     char PATH[128];
     char RESPATH[128];

     char FileName[128]; //used in writeList

     double BEGINTime;
     double ENDTime;
     double DTSec;
     float WRITETIME[100];
     int NrWriteTime;
     // output moments in the ineterface, so write maps at 1, 10, 24, 36 etc.
     int PERIOD;
     bool PERIODMAPS;
     int PERIODCOUNT;
     int LastPCRTimestep;

     void *soilModelcrackV;
     void *soilModelwheelV;
     void *soilModelcrustV;
     void *soilModelcompactV;
     void *soilModelgrassV;

     double mu_cl[6];

     int INFIL_METHOD; //infiltration method number, see lisparaminfil.cpp

     FILE *fpoutflow1;
     FILE *fpoutflow2;
     FILE *fpoutflow3;
     FILE *BufferFout;
     MEM_HANDLE *Mout;
      // hydrograph output files declaration


//VJ 050913 commandline param pestout
	  float pestreport;
     int pestcounter;
     FILE *fpestout;

     float timestepindex; // the main counter for the timeloop

     void __fastcall LisThread::LisemWarningV();
     void __fastcall LisThread::LisemWarn(AnsiString s, AnsiString s1 = "");
     void __fastcall LisThread::UpdateTotalsSync();
     void __fastcall LisThread::UpdateTotalsPestSync();
     void __fastcall LisThread::InitializeTotalsSync();
     void __fastcall LisThread::DumpScreenSync();
     char* __fastcall LisThread::Lmapname(char *vname, int nr);
     float __fastcall LisThread::GetFloat(char *vname);
     int __fastcall LisThread::GetInt(char *vname);
     void __fastcall LisThread::ParseInputData();
     void __fastcall LisThread::GetRunFile();
     void __fastcall LisThread::SetTimeseriesMinmax(char *filename);
     void __fastcall LisThread::DrawMapSync(void);//MEM_HANDLE *M, double timestep);
     void __fastcall LisThread::DrawMapInitSync(void);//(MEM_HANDLE *M, double timestep);

      AnsiString Convert(double f, int dig = 3);

    __fastcall LisThread(bool CreateSuspended);

    //INPUT variables
      bool SwitchCorrectMass;
      bool SwitchCorrectMassSED;
      bool SwatreInitialized;
      bool SwitchCrustPresent;
      bool SwitchGrassPresent;
      bool SwitchWheelPresent;
      bool SwitchCompactPresent;
      bool SwitchInfilGA2;
      bool SwitchIncludeChannel;
      bool SwitchChannelBaseflow;
      bool startbaseflowincrease;
      bool SwitchChannelInfil;
      bool SwitchAllinChannel;
      bool SwitchKinwaveInfil;
      bool SwitchNoErosion;
      bool SwitchAltErosion;
      bool SwitchSimpleDepression;
      bool SwitchHardsurface;
      bool SwitchBuffers;
      bool SwitchSedtrap;
      bool SwitchSnowmelt;
      bool SwitchRunoffPerM;
      bool SwitchInfilCompact;
      bool SwitchInfilCrust;
      bool SwitchInfilGrass;
      bool SwitchImpermeable;
      bool SwitchDumphead;
      bool SwitchGeometricMean;
      bool SwitchWheelAsChannel;
      bool SwitchMulticlass;
      bool SwitchNutrients;
      bool SwitchGullies;
      bool SwitchOutputTimeStep;
      bool SwitchOutputTimeUser;
      bool SwitchMapoutRunoff;
      bool SwitchMapoutConc;
      bool SwitchMapoutWH;
      bool SwitchMapoutWHC;
      bool SwitchMapoutTC;
      bool SwitchMapoutEros;
      bool SwitchMapoutDepo;
      bool SwitchMapoutV;
      bool SwitchMapoutInf;
      bool SwitchMapoutSs;
      bool SwitchMapoutChvol;
      bool SwitchMapoutMC0;
      bool SwitchMapoutMC1;
      bool SwitchMapoutMC2;
      bool SwitchMapoutMC3;
      bool SwitchMapoutMC4;
      bool SwitchMapoutMC5;
      bool SwitchMapoutMC6;
      bool SwitchMapoutPsol;
      bool SwitchMapoutPsus;
      bool SwitchMapoutPinf;
      bool SwitchMapoutNH4sol;
      bool SwitchMapoutNH4sus;
      bool SwitchMapoutNH4inf;
      bool SwitchMapoutNO3sol;
      bool SwitchMapoutNO3sus;
      bool SwitchMapoutNO3inf;
      bool SwitchMapoutPdep;
      bool SwitchMapoutNH4dep;
      bool SwitchMapoutNO3dep;
      bool SwitchMapoutPdet;
      bool SwitchMapoutNH4det;
      bool SwitchMapoutNO3det;
      bool SwitchMapoutGul0;
      bool SwitchMapoutGul1;
      bool SwitchMapoutGul2;
      bool SwitchMapoutGul3;
      bool SwitchMapoutGul4;
      bool SwitchWritePCRnames;
      bool SwitchWritePCRtimeplot;
      bool SwitchOutlet1;
      bool SwitchOutlet2;
      bool SwitchOutlet3;
      bool SwitchNoErosionOutlet;
      bool SwitchGullyEqualWD;
      bool SwitchGullyInit;
      bool SwitchGullyInfil;
      bool SwitchDrainage;
      char ErosionUnits;
      bool SwitchPestout;
      bool SwitchSeparateOutput;
      bool SwitchSOBEKOutput;
      bool SwitchMinimumdisplay;
      char SOBEKdatestring[12];
      int SOBEKnrlines;
      bool SwitchInterceptionLAI;
      int  InterceptionLAIType;
      bool runv3;

      float ksatCalibration;
      float nCalibration;
      float ChnCalibration;
      float ChKsatCalibration;
      double runtimelength;

  // these are vars that link the model to the interface
  // purpose is to make the model and the interface independent
     double StartTime;
     double EndTime;
     double CurrentTime;
     double timestep;
     double CellSize;
     double CatchmentAreaHa;
     double TOTRainMM;
     double AVGRainMM;
     double TOTInterceptionMM;
     double TOTInfilMM;
     double TOTSurfStorMM;
     double TOTRunoffMM;
     double TOTDischargeMM;
     double TOTDischargeM3;
     double P_QPercentage;
     double PeakDischargeQ;
     double BaseflowDischargeQ;
     double PeakRainfallTime;
     double PeakDischargeTime;
     double TOTBufferVolume;
     double TOTBufferSedVolume;
     double TOTFlowErosion;
     double TOTSplashErosion;
     double TOTDeposition;
     double TOTChanFlowErosion;
     double TOTChanDeposition;
     double TOTSoilLoss;
     double AVGSoilLoss;
     double OutputDischarge;
     double OutputDischarge1;
     double OutputDischarge2;
     double OutputSnow;
     double OutputSedConc;
     double OutputSediment;
     double OutputSedSusp;
     double OutputChanSedSusp;
     double MassBalance;
     double SedMassBalance;

};

#endif
