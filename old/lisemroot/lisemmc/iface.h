//---------------------------------------------------------------------------
#ifndef ifaceH
#define ifaceH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Dialogs.hpp>
#include <ExtCtrls.hpp>
#include <Menus.hpp>
#include <ComCtrls.hpp>
#include <DBCtrls.hpp>
#include <ToolWin.hpp>
#include <Buttons.hpp>
#include "cgauges.h"
#include <graphsv3.hpp>
#include <OleCtrls.hpp>
#include <Chart.hpp>
#include <TeEngine.hpp>
#include <TeeProcs.hpp>
#include <Series.hpp>
#include <jpeg.hpp>
#include <Filectrl.hpp>
//#include <systdate.h>

#include "lismain.h"
#include "CGAUGES.h"
#include <ExtDlgs.hpp>
#include "CSPIN.h"
#include <Grids.hpp>
#include <Mask.hpp>

#define ProgName "LISEM"
#define ProgVersion "version 2.58"
#define DateVersion "Compiled: 14/02/09"

#define LISEMBASIC 0
#define LISEMWHEELTRACKS 1
#define LISEMMULTICLASS 2
#define LISEMNUTRIENTS 3
#define LISEMGULLIES 4
#define wWidth 800
#define wHeight 600

//---------------------------------------------------------------------------
class TLisIFace : public TForm
{
__published:	// IDE-managed Components
   TOpenDialog *OpenDialog;
   TSaveDialog *SaveDialog;
   TToolBar *ToolBar;
   TSpeedButton *FileOpenButton;
   TSpeedButton *FileSaveButton;
   TSpeedButton *FileSaveAsButton;
   TSpeedButton *ButtonRunprog;
   TSplitter *Splitter1;
   TSpeedButton *ButtonStopprog;
        TPanel *Panel5;
        TCGauge *CGauge;
        TSpeedButton *ButtonPauseprog;
        TSpeedButton *PrintButton;
        TSpeedButton *DumpScreenButton;
	TSpeedButton *AboutButton;
        TSplitter *Splitter3;
        TSpeedButton *StopAllButton;
	TSplitter *Splitter4;
	TSplitter *Splitter5;
	TSpeedButton *HelpButton;
        TPageControl *PageControl;
        TTabSheet *TabSheetTot;
        TBevel *Bevel2;
        TBevel *Bevel3;
        TBevel *Bevel5;
        TBevel *Bevel1;
        TLabel *Label1;
        TLabel *Label2;
        TLabel *Label3;
        TLabel *Label4;
        TLabel *Label5;
        TLabel *Label6;
        TLabel *Label8;
        TLabel *Label9;
        TLabel *Label10;
        TLabel *Label11;
        TLabel *Label12;
        TLabel *Label13;
        TLabel *Label17;
        TLabel *Label18;
        TBevel *Bevel8;
        TBevel *Bevel9;
        TBevel *Bevel11;
        TBevel *Bevel20;
        TBevel *Bevel21;
        TLabel *warning;
        TLabel *Infil;
        TBevel *Bevel24;
        TLabel *Label25;
        TLabel *Label26;
        TLabel *Label28;
        TLabel *Label29;
        TLabel *Label27;
        TLabel *LO_infilmethod;
        TLabel *Label30;
        TLabel *Label32;
        TBevel *Bevel23;
        TLabel *Label23;
        TLabel *Label24;
        TLabel *Label33;
        TLabel *Label34;
        TLabel *Label35;
        TPanel *Panel2;
        TPanel *Panel4;
        TPanel *Panel6;
        TMemo *Messages;
        TLabel *LabelRunfile;
        TSavePictureDialog *SavePictureDialog;
        TPrintDialog *PrintDialog;
        TTabSheet *TabMapNames;
    TPageControl *Names;
    TTabSheet *Area;
    TLabel *Label109;
    TLabel *Label65;
        TStringGrid *MapsCatchment;
    TStringGrid *MapsLanduse;
    TTabSheet *Surface;
    TLabel *Label92;
    TStringGrid *MapsErosion;
    TTabSheet *Infiltration;
        TLabel *LabelInfilName;
    TStringGrid *MapsInfilExtra;
    TStringGrid *MapsInfilHoltan;
        TStringGrid *MapsInfilGA2;
        TStringGrid *MapsInfilGA1;
    TStringGrid *MapsInfilSwatre;
    TTabSheet *Channels;
    TLabel *Label42;
    TStringGrid *MapsChannelinfil;
    TStringGrid *MapsChannels;
    TTabSheet *TabSheetRun;
    TLabel *Label90;
    TStringGrid *MapsSurface;
    TGroupBox *GroupOptionsModel;
    TCheckBox *CheckNoErosion;
    TCheckBox *CheckIncludeChannel;
    TGroupBox *GroupDirsInput;
    TBevel *Bevel29;
    TLabel *Label51;
    TSpeedButton *LoadMapDir;
    TLabel *Label22;
    TSpeedButton *LoadWorkDir;
    TBevel *Bevel16;
    TEdit *E_Workdir;
    TEdit *E_MapDir;
    TGroupBox *GroupOptionsSimTime;
    TLabel *Label38;
    TLabel *Label37;
    TLabel *Label36;
    TBevel *Bevel25;
    TEdit *E_begintime;
    TEdit *E_Endtime;
    TEdit *E_Timestep;
    TGroupBox *GroupRunbatch;
    TBevel *Bevel28;
    TSpeedButton *EraseListButton;
    TBitBtn *RunFileButton;
    TListBox *ListRunfiles;
    TButton *ButtonRestore;
    TTabSheet *TabMapsOutput;
        TCheckBox *CheckWritePCRnames;
        TGroupBox *GroupOptionsInfil;
        TRadioGroup *E_InfilMethod;
        TGroupBox *GroupOptionsInfilExtra;
        TLabel *Label64;
        TCSpinEdit *CSpinEditKsat;
        TSpeedButton *SpeedButton1;
        TTabSheet *TabWheeltracks;
        TStringGrid *MapsWheeltrack;
        TTabSheet *TabMulticlass;
        TLabel *Label39;
        TStringGrid *MapsTexture;
        TStringGrid *TextureClass;
        TTabSheet *TabNutrients;
        TLabel *Label69;
        TLabel *Label70;
        TLabel *Label71;
        TStringGrid *MapsNutsP;
        TStringGrid *MapsNutsNO3;
        TStringGrid *MapsNutsNH4;
        TTabSheet *TabGullies;
        TGroupBox *GroupSwatreOptions;
        TBevel *Bevel31;
        TLabel *Label56;
        TBevel *Bevel26;
        TLabel *Label52;
        TSpeedButton *LoadTableDir;
        TBevel *Bevel33;
        TEdit *E_SwatreDTSEC;
        TCheckBox *CheckDumphead;
        TCheckBox *CheckGeometric;
        TEdit *E_TableDir;
        TStringGrid *MapsInfilMorel;
        TStringGrid *MapsInfilSmith;
        TGroupBox *GroupSurfacetype;
        TCheckBox *CheckInfilCompact;
        TCheckBox *CheckInfilCrust;
        TCheckBox *CheckInfilGrass;
        TCheckBox *CheckRunoffPerM;
        TLabel *Label53;
    TGroupBox *GroupOptionsOutTime;
    TMemo *E_OutputTimes;
    TRadioButton *E_OutputTimeUser;
    TRadioButton *E_OutputTimeStep;
    TCSpinEdit *E_OutputTimeSteps;
    TLabel *Label62;
    TLabel *Label68;
    TBevel *Bevel36;
    TRadioGroup *E_QWrelation;
    TStringGrid *MapsGully;
    TRadioGroup *E_Fcritical;
    TEdit *E_ThresholdGrad;
    TPanel *Panel7;
    TLabel *Label74;
    TBevel *Bevel38;
    TLabel *Label76;
    TEdit *E_QWparama;
    TEdit *E_QWparamb;
    TPanel *Panel9;
    TStringGrid *MapsNutsBD;
    TBevel *Bevel39;
    TBevel *Bevel40;
    TBevel *Bevel41;
        TPanel *Panel10;
        TChart *Qgraph;
        TLineSeries *Series1;
        TLineSeries *Series2;
        TLineSeries *Series3;
        TLineSeries *Series4;
        TLineSeries *Series6;
        TLineSeries *Series7;
        TPanel *Panel1;
        TCheckBox *ShowRainfall;
        TCheckBox *ShowOutlet1;
        TCheckBox *ShowOutlet2;
        TCheckBox *ShowOutlet3;
        TCheckBox *ShowQsed;
        TCheckBox *ShowSedConcentration;
        TRichEdit *OutletValues;
        TCheckBox *CheckChannelInfil;
        TPageControl *MapOutControl;
        TTabSheet *TabMapOutMC;
        TTabSheet *TabMapOutNut;
        TTabSheet *TabMapOutGul;
        TStringGrid *MapsOutputMC;
        TStringGrid *MapsOutputGul;
        TPanel *GroupMapsoutGul;
        TCheckBox *CheckMapoutGul1;
        TCheckBox *CheckMapoutGul3;
        TCheckBox *CheckMapoutGul2;
        TCheckBox *CheckMapoutGul0;
        TCheckBox *CheckMapoutGul4;
        TPanel *GroupMapsoutNut;
        TCheckBox *CheckMapoutNut1;
        TCheckBox *CheckMapoutNut3;
        TCheckBox *CheckMapoutNut2;
        TCheckBox *CheckMapoutNut5;
        TCheckBox *CheckMapoutNut4;
        TCheckBox *CheckMapoutNut6;
        TCheckBox *CheckMapoutNut8;
        TCheckBox *CheckMapoutNut7;
        TCheckBox *CheckMapoutNut0;
        TStringGrid *MapsOutputNut;
        TPanel *GroupMapsoutMC;
        TCheckBox *CheckMapoutMC1;
        TCheckBox *CheckMapoutMC3;
        TCheckBox *CheckMapoutMC2;
        TCheckBox *CheckMapoutMC5;
        TCheckBox *CheckMapoutMC4;
        TCheckBox *CheckMapoutMC0;
        TLabel *LO_catcharea;
        TLabel *LO_Tstart;
        TLabel *LO_dtcurr;
        TLabel *LO_totrain;
        TLabel *LO_interception;
        TLabel *LO_infil;
        TLabel *LO_surfstor;
        TLabel *LO_totdisch;
        TLabel *LO_raindisch;
        TLabel *LO_Runoff;
        TLabel *LO_gridcell;
        TLabel *LO_Tend;
        TLabel *LO_timestep;
        TLabel *LO_peaktimeP;
        TLabel *LO_peaktime;
        TLabel *LO_peakdisch;
        TLabel *LO_masserror;
        TLabel *Label49;
        TLabel *LO_chanflow;
        TLabel *LO_chandepo;
        TLabel *LO_ChSedSuspended;
        TLabel *LO_Sedbal;
    TTabSheet *TabMapOutBasic;
    TStringGrid *MapsOutputBASIC;
    TCheckBox *CheckMapout7;
    TCheckBox *CheckMapout4;
    TCheckBox *CheckMapout0;
    TCheckBox *CheckMapout1;
    TCheckBox *CheckMapout2;
    TCheckBox *CheckMapout5;
    TCheckBox *CheckMapout3;
    TCheckBox *CheckMapout6;
    TCheckBox *CheckMapout8;
    TCheckBox *CheckMapout9;
    TCheckBox *CheckWritePCRtimeplot;
    TCheckBox *CheckMapoutMC6;
    TStringGrid *MapsOutputNutErosDep;
    TLabel *Label50;
    TPanel *GroupMaspoutNutErosDepo;
    TCheckBox *CheckMapoutNut10;
    TCheckBox *CheckMapoutNut12;
    TCheckBox *CheckMapoutNut11;
    TCheckBox *CheckMapoutNut14;
    TCheckBox *CheckMapoutNut13;
    TCheckBox *CheckMapoutNut9;
//	TLMDSpinEdit *LMDSpinOutletValues;
    TGroupBox *GroupDirsOutput;
    TLabel *Label41;
    TLabel *Label43;
    TLabel *Label45;
    TLabel *Label46;
    TLabel *Label47;
    TLabel *Label48;
    TBevel *Bevel15;
    TBevel *Bevel6;
    TBevel *Bevel32;
    TEdit *E_ErosionName;
    TEdit *E_DepositionName;
    TEdit *E_OutletName;
    TEdit *E_Outlet1Name;
    TEdit *E_Outlet2Name;
    TEdit *E_TotalName;
    TCSpinEdit *SpinOutletValues;
    TLabel *Label40;
    TSpeedButton *LoadResultDir;
    TEdit *E_ResultDir;
    TBevel *Bevel19;
   TCheckBox *CheckNoErosionOutlet;
   TStringGrid *MapsGullyInit;
   TCheckBox *CheckGullyInit;
   TStringGrid *MapsBuffers;
   TLabel *Label55;
   TCheckBox *CheckImpermeable;
   TCheckBox *CheckBuffers;
   TLabel *Label58;
   TLabel *Label60;
   TLabel *Label57;
   TBevel *Bevel7;
   TEdit *E_SedBulkDensity;
   TLabel *Label7;
   TLabel *Label19;
   TLabel *Label20;
   TLabel *Label21;
   TBevel *Bevel18;
   TBevel *Bevel12;
   TLabel *LO_outflow;
   TLabel *LO_soilloss;
   TLabel *LO_dischM3;
   TPanel *Panel8;
   TLabel *LO_avgsoilloss;
   TBevel *Bevel13;
   TLabel *Label63;
   TLabel *LO_buffervol;
   TLabel *Label61;
   TLabel *LO_buffersedvol;
   TBevel *Bevel43;
   TBevel *Bevel44;
	TCheckBox *CheckAltErosion;
    TCheckBox *CheckSubsoilDrainage;
    TStringGrid *MapsInfilDrainage;
	TRadioGroup *E_OutputUnits;
   TBevel *Bevel45;
   TCheckBox *CheckGullyInfil;
   TTabSheet *Macropores;
   TStringGrid *MapsMacropore;
   TLabel *Label66;
   TCheckBox *CheckMacroporeFlow;
   TCheckBox *CheckChannelBaseflow;
   TStringGrid *MapsChannelBaseflow;
   TLabel *LO_baseflow;
   TLabel *Label72;
   TCheckBox *CheckMapout10;
   TCheckBox *CheckSnowmelt;
   TBevel *Bevel30;
   TLabel *Label54;
   TSpeedButton *RainViewButton;
   TSpeedButton *LoadRainfile;
   TBevel *Bevel34;
   TEdit *E_RainfallName;
   TLabel *Label73;
   TSpeedButton *LoadSnowmeltfile;
   TEdit *E_SnowmeltName;
   TSpeedButton *SnowmeltViewButton;
   TStringGrid *MapsInfilKsat;
   TLabel *Label67;
   TStringGrid *MapsSnowmelt;
   TLabel *Label75;
   TPanel *Panel3;
   TLabel *Label14;
   TLabel *Label15;
   TLabel *Label16;
   TBevel *Bevel10;
   TLabel *Label31;
   TLabel *LO_flow;
   TLabel *LO_depo;
   TLabel *LO_SedSuspended;
   TLabel *LO_splash;
   TLabel *Label77;
   TLabel *Label59;
   TLabel *Label44;
   TEdit *E_SplashDelivery;
   TEdit *E_ManningsNGrass;
   TBevel *Bevel27;
   TBevel *Bevel4;
   TCheckBox *CheckAltDepression;
   TCheckBox *CheckAllinChannel;
//        TNMPOP3 *NMPOP31;
    void __fastcall ButtonRunprogClick(TObject *Sender);
    void __fastcall FileOpenClick(TObject *Sender);
    void __fastcall FileSaveClick(TObject *Sender);
    void __fastcall FileSaveAsClick(TObject *Sender);
    void __fastcall FileCloseClick(TObject *Sender);
    void __fastcall ButtonStopprogClick(TObject *Sender);
    void __fastcall ButtonPauseprogClick(TObject *Sender);
    void __fastcall MakeIniButtonClick(TObject *Sender);
    void __fastcall PrintButtonClick(TObject *Sender);
    void __fastcall DumpScreenButtonClick(TObject *Sender);
    void __fastcall RunThread();
    void __fastcall ThreadDone(TObject *Sender);
    void __fastcall AboutButtonClick(TObject *Sender);
    void __fastcall StopAllButtonClick(TObject *Sender);
	void __fastcall CheckIncludeChannelClick(TObject *Sender);
	void __fastcall FormClose(TObject *Sender, TCloseAction &Action);
	void __fastcall HelpButtonClick(TObject *Sender);
	void __fastcall SizeButtonClick(TObject *Sender);
	void __fastcall FormResize(TObject *Sender);
        void __fastcall CheckNoErosionClick(TObject *Sender);
	void __fastcall E_TableDirClick(TObject *Sender);
	void __fastcall E_RainfallNameClick(TObject *Sender);
	void __fastcall RunFileButtonClick(TObject *Sender);
        void __fastcall E_OutputTimeStepClick(TObject *Sender);
        void __fastcall E_OutputTimeUserClick(TObject *Sender);
        void __fastcall RainViewButtonClick(TObject *Sender);
    void __fastcall LoadRainfileClick(TObject *Sender);
	void __fastcall LoadTableDirClick(TObject *Sender);
        void __fastcall LoadWorkDirClick(TObject *Sender);
	void __fastcall LoadMapDirClick(TObject *Sender);
        void __fastcall LoadResultDirClick(TObject *Sender);
        void __fastcall E_InfilMethodClick(TObject *Sender);
        void __fastcall EraseListButtonClick(TObject *Sender);
    void __fastcall ListRunfilesDblClick(TObject *Sender);
        void __fastcall MapsCatchmentDblClick(TObject *Sender);
        void __fastcall MapsErosionDblClick(TObject *Sender);
        void __fastcall MapsSurfaceDblClick(TObject *Sender);
        void __fastcall MapsInfilGA1DblClick(TObject *Sender);
    void __fastcall MapsWheeltrackDblClick(TObject *Sender);
    void __fastcall MapsChannelsDblClick(TObject *Sender);
    void __fastcall MapsChannelinfilDblClick(TObject *Sender);
    void __fastcall MapsOutputBASICDblClick(TObject *Sender);
    void __fastcall MapsInfilExtraDblClick(TObject *Sender);
    void __fastcall MapsInfilHoltanDblClick(TObject *Sender);
    void __fastcall MapsInfilGA2DblClick(TObject *Sender);
    void __fastcall MapsInfilSwatreDblClick(TObject *Sender);
    void __fastcall MapsTextureDblClick(TObject *Sender);
    void __fastcall MapsLanduseDblClick(TObject *Sender);
    void __fastcall CheckInfilCrustClick(TObject *Sender);
    void __fastcall CheckInfilCompactClick(TObject *Sender);
    void __fastcall TextureClassClick(TObject *Sender);
    void __fastcall TextureClassKeyDown(TObject *Sender, WORD &Key,
          TShiftState Shift);
    void __fastcall ButtonRestoreClick(TObject *Sender);
    void __fastcall SpeedButton1Click(TObject *Sender);
        void __fastcall MapsInfilSmithDblClick(TObject *Sender);
        void __fastcall MapsInfilMorelDblClick(TObject *Sender);
        void __fastcall CheckInfilGrassClick(TObject *Sender);
        void __fastcall MapsGullyDblClick(TObject *Sender);
        void __fastcall MapsNutsPDblClick(TObject *Sender);
        void __fastcall MapsNutsNO3DblClick(TObject *Sender);
        void __fastcall MapsNutsNH4DblClick(TObject *Sender);
        void __fastcall MapsNutsBDDblClick(TObject *Sender);
    void __fastcall E_MapDirDblClick(TObject *Sender);
    void __fastcall E_WorkdirDblClick(TObject *Sender);
    void __fastcall E_ResultDirDblClick(TObject *Sender);
    void __fastcall E_ResultDirKeyPress(TObject *Sender, char &Key);
    void __fastcall E_MapDirKeyPress(TObject *Sender, char &Key);
    void __fastcall E_WorkdirKeyPress(TObject *Sender, char &Key);
   void __fastcall E_ResultDirExit(TObject *Sender);
   void __fastcall E_MapDirExit(TObject *Sender);
   void __fastcall E_WorkdirExit(TObject *Sender);
   void __fastcall MapsGullyInitDblClick(TObject *Sender);
   void __fastcall MapsBuffersDblClick(TObject *Sender);
   void __fastcall CheckBuffersClick(TObject *Sender);
    void __fastcall CheckSubsoilDrainageClick(TObject *Sender);
    void __fastcall CheckImpermeableClick(TObject *Sender);
    void __fastcall MapsInfilDrainageDblClick(TObject *Sender);
   void __fastcall CheckMacroporeFlowClick(TObject *Sender);
   void __fastcall MapsMacroporeDblClick(TObject *Sender);
   void __fastcall CheckChannelBaseflowClick(TObject *Sender);
   void __fastcall CheckChannelInfilClick(TObject *Sender);
   void __fastcall CheckSnowmeltClick(TObject *Sender);
   void __fastcall MapsSnowmeltDblClick(TObject *Sender);
   void __fastcall LoadSnowmeltfileClick(TObject *Sender);
   void __fastcall SnowmeltViewButtonClick(TObject *Sender);
   void __fastcall MapsInfilKsatDblClick(TObject *Sender);
   void __fastcall MapsChannelBaseflowDblClick(TObject *Sender);
   void __fastcall CheckAllinChannelClick(TObject *Sender);

private:	// User declarations
   AnsiString __fastcall TLisIFace::CheckDir(AnsiString Comment, AnsiString dir);
   AnsiString __fastcall TLisIFace::CheckFile(AnsiString dir);
   bool blockaction;
   AnsiString OldMapDir;

public:		// User declarations
   int nrruns;
   int thisrun;
//   bool done;
   bool closeapp;
   bool take5;
   bool notrunning;
   bool multipleRuns;
   bool batchrun;
//   bool __fastcall TLisIFace::CheckRunFile(TMemo *M);
   void __fastcall TLisIFace::LoadRunFile(AnsiString RunFilename);
   void __fastcall TLisIFace::MakeIni(AnsiString name);
   void __fastcall TLisIFace::MakeNewRunfile(AnsiString name);
   bool __fastcall TLisIFace::ReadNewRunfile(AnsiString name);
   void __fastcall TLisIFace::ParseOldRunfile();
   bool __fastcall TLisIFace::ProcessIni(AnsiString name);
   void __fastcall TLisIFace::ShowHelp();
   bool __fastcall TLisIFace::CheckError(AnsiString S);
   bool __fastcall TLisIFace::DumpScreen(bool saveio);
   void __fastcall TLisIFace::GetMapname(TStrings *S);
   void __fastcall TLisIFace::InitMapNames();
   void __fastcall TLisIFace::InitOutMapNames();
   void __fastcall TLisIFace::FillMapNames(TStringGrid *N, int i,
   AnsiString Varname, AnsiString Filename, AnsiString Desc, AnsiString ID);
   void __fastcall TLisIFace::FillMapNamesO(TStringGrid *N, int i,
   AnsiString Varname, AnsiString Filename, AnsiString Desc, AnsiString ID);
   void __fastcall TLisIFace::SizeMapNames(TStringGrid *N, int w1, int w2,int w3);
   void __fastcall TLisIFace::SizeMapNamesL(TStringGrid *N, int w1, int w2,int w3, int w4);

   void __fastcall TLisIFace::GetNewRunfile();

   void __fastcall TLisIFace::SetTimeseriesMinmax(char *filename);

   void __fastcall TLisIFace::CheckMapsEnabled();
   void __fastcall TLisIFace::TextureClassCheck();
   void __fastcall TLisIFace::CheckLisemType();
   void __fastcall TLisIFace::GetDirEdit(TEdit *E);
   void __fastcall CheckMulticlassModel();
   void __fastcall CheckNutrientsModel();
   void __fastcall CheckWheelAsChannelModel();
   void __fastcall CheckGulliesModel();
   bool CheckWheelAsChannel;
   bool CheckGullies;
   bool CheckNutrients;
   bool CheckMulticlass;

   bool CheckAdjustMapDirectoryName;

   // VJ 031218 check dumpscreen to halt processes to finish dumpscreen
   bool CheckDumpScreen;

   float ClassMu[6];
   char DirIni[128];
   char DirProg[128];
   char DirLog[128];
   char workdir[128];
   char RainfallNamePath[128];
   char SnowmeltNamePath[128];
   bool ForceMapdir;
   bool CheckPestout;
   float PestoutTimeinterval;
   LisThread *LisemRun;

   AnsiString RunFilename;
   AnsiString RainfallDir;
   AnsiString SnowmeltDir;
   __fastcall TLisIFace(TComponent* Owner);
   __fastcall ~TLisIFace();

   int LisemType;

//   void __fastcall TLisIFace::MakeNamelist();
//   void __fastcall TLisIFace::KillNamelist();
   int LastPCRTimestep;

   void __fastcall TLisIFace::ResetIFace();
   void __fastcall TLisIFace::ResetInfil();
   void __fastcall TLisIFace::ResetMain();

   void __fastcall TLisIFace::DoMapDir(bool CheckName, bool CheckForce);

   TDateTime t_begin;
   TDateTime t_end;
   AnsiString runduration;
   double runlength;
};
//---------------------------------------------------------------------------
extern PACKAGE TLisIFace *LisIFace;
//---------------------------------------------------------------------------
#endif
