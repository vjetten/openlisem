//---------------------------------------------------------------------------
#ifndef LisnewH
#define LisnewH
//---------------------------------------------------------------------------
#include <Classes.hpp>

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



typedef struct {
   char ID[64];
   char fullname[128];
} nameList;


class LisThread : public TThread
{
private:
protected:
    void __fastcall Execute();
    AnsiString WarningText;
    AnsiString ReportText;
public:
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

   int StepCounter;
   int errorlevel;
   bool InitDone;
   bool RunDone;
   nameList *namelist;
   int nrnamelist;
   char temprunname[512];

   void __fastcall LisemWarningV();
   void __fastcall LisThread::LisemWarn(AnsiString s, AnsiString s1 = "");
   void __fastcall UpdateTotals();
   void __fastcall InitializeTotals();
   void __fastcall LisThread::ShowHelp();
   void __fastcall LisThread::DoHelp();
   void __fastcall LisThread::DumpScreenSync();
   char* __fastcall LisThread::Lmapname(char *vname);
   void __fastcall LisThread::MakeMapnameList();
   void __fastcall LisThread::MakeSwitchList();
   
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
      bool SwitchAltDepression;
      bool SwitchBuffers;
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
      bool SwitchMapoutRWH;
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

      double ksatCalibration;
      double runtimelength; 

};

#endif
