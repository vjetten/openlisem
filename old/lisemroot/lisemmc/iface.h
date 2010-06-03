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

#include "CGAUGES.h"
#include <ExtDlgs.hpp>
#include "CSPIN.h"
#include <Grids.hpp>
#include <Mask.hpp>
#include <ImgList.hpp>
#include "JvExMask.hpp"
#include "JvSpin.hpp"
#include "JvBaseDlg.hpp"
#include "JvBrowseFolder.hpp"
#include "JvSelectDirectory.hpp"

#include "lismain.h"
#include "lishelp.h"
#include "lisabout.h"
#include "lisrunf.h"
#include "lisrainf.h"
#include "lisstart.h"
#include "ifaceinit.h"
#include "ifacethread.h"
#include "lishint.h"
#include "lisscreenoutput.h"
#include "CsfMapDraw.h"
#include "JvButton.hpp"
#include "JvExStdCtrls.hpp"
#include "JvRecentMenuButton.hpp"
#include "JvCheckedMaskEdit.hpp"
#include "JvMaskEdit.hpp"
#include "JvToolEdit.hpp"
#include "JvExControls.hpp"
#include "JvExExtCtrls.hpp"
#include "JvOutlookBar.hpp"
#include "JvRadioGroup.hpp"
#include "JvExtComponent.hpp"
#include "JvItemsPanel.hpp"
#include "JvAVICapture.hpp"
#include "JvLookOut.hpp"
#include "JvStaticText.hpp"
#include "JvComponentBase.hpp"
#include "JvFormAutoSize.hpp"
#include "JvSpecialProgress.hpp"
#include "JvProgressBar.hpp"
#include "JvMenus.hpp"
#include "JvFormMagnet.hpp"
#include "JvRadioButton.hpp"


#define ProgName "LISEM"
#define ProgVersion "version 2.64"
#define DateVersion "Compiled: 10/01/21"

#define LISEMBASIC 0
#define LISEMWHEELTRACKS 1
#define LISEMMULTICLASS 2
#define LISEMNUTRIENTS 3
#define LISEMGULLIES 4
#define wWidth 920
#define wHeight 700


//---------------------------------------------------------------------------
typedef struct _VPoint {int x, y;} VPoint;
typedef struct _EPoint {
                 int x, y;
                 int r, c;
                } EPoint;

//---------------------------------------------------------------------------
class TLisIFace : public TForm
{
__published:	// IDE-managed Components
   TOpenDialog *OpenDialog;
   TSaveDialog *SaveDialog;
   TToolBar *ToolBar;
   TSplitter *Splitter1;
	TSpeedButton *AboutButton;
    TSplitter *Splitter3;
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
    TGroupBox *GroupDirsInput;
    TLabel *Label51;
    TSpeedButton *LoadMapDir;
    TBevel *Bevel16;
    TEdit *E_MapDir;
    TGroupBox *GroupOptionsSimTime;
    TLabel *Label38;
    TLabel *Label37;
    TLabel *Label36;
    TBevel *Bevel25;
    TEdit *E_begintime;
    TEdit *E_Endtime;
    TEdit *E_Timestep;
    TButton *ButtonRestore;
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
        TStringGrid *MapsInfilMorel;
        TStringGrid *MapsInfilSmith;
        TLabel *Label53;
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
        TCheckBox *ShowQsed;
        TCheckBox *ShowSedConcentration;
        TRichEdit *OutletValues;
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
    TCSpinEdit *SpinOutletValues;
   TStringGrid *MapsGullyInit;
   TCheckBox *CheckGullyInit;
   TStringGrid *MapsBuffers;
   TLabel *Label55;
   TLabel *Label58;
   TLabel *Label60;
   TLabel *Label7;
   TLabel *Label19;
   TLabel *Label20;
   TLabel *Label21;
   TBevel *Bevel18;
   TBevel *Bevel12;
   TLabel *LO_outflow;
   TLabel *LO_soilloss;
   TLabel *LO_dischM3;
   TLabel *LO_avgsoilloss;
   TBevel *Bevel13;
   TLabel *Label63;
   TLabel *LO_buffervol;
   TLabel *Label61;
   TLabel *LO_buffersedvol;
   TBevel *Bevel43;
   TBevel *Bevel44;
    TStringGrid *MapsInfilDrainage;
   TBevel *Bevel45;
   TCheckBox *CheckGullyInfil;
   TLabel *Label66;
   TStringGrid *MapsChannelBaseflow;
   TLabel *LO_baseflow;
   TLabel *Label72;
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
   TImageList *ImageList1;
   TToolButton *ToolButtonFileopen;
   TToolButton *ToolButtonFilesave;
   TToolButton *ToolButtonFilesaveAs;
   TToolButton *ToolButtonLisemstart;
   TToolButton *ToolButtonDumpScreen;
   TToolButton *ButtonRunProg;
   TToolButton *ButtonPauseProg;
   TToolButton *ButtonStopProg;
   TToolButton *ButtonStopAll;
   TBevel *Bevel29;
   TGroupBox *GroupDirsOutput;
   TLabel *Label41;
   TLabel *Label43;
   TLabel *Label45;
   TLabel *Label48;
   TBevel *Bevel15;
   TBevel *Bevel6;
   TBevel *Bevel32;
   TSpeedButton *LoadResultDir;
   TLabel *Label40;
   TBevel *Bevel22;
   TEdit *E_ErosionName;
   TEdit *E_DepositionName;
   TEdit *E_OutletName;
   TEdit *E_TotalName;
   TEdit *E_ResultDir;
   TGroupBox *GroupOptionsOutTime;
   TMemo *E_OutputTimes;
   TRadioButton *E_OutputTimeUser;
   TRadioButton *E_OutputTimeStep;
   TCSpinEdit *E_OutputTimeSteps;
   TGroupBox *GroupBox1;
   TCheckBox *CheckWritePCRtimeplot;
   TTabSheet *TabSheetSpare;
   TCheckBox *CheckAllinChannel;
   TEdit *E_Workdir;
   TLabel *Label46;
   TLabel *Label78;
   TLabel *Label50;
   TEdit *E_SoillossName;
   TToolButton *ToolButtonHint;
   TLabel *Label79;
   TLabel *Label82;
   TJvBrowseForFolderDialog *JvGetDirDialog;
   TLabel *Label83;
   TCheckBox *CheckSOBEKOutput;
   TCheckBox *CheckRunoffPerM;
   TCheckBox *CheckWritePCRnames;
   TGroupBox *GroupBox3;
   TPageControl *MapOutControl;
   TTabSheet *TabMapOutBasic;
   TCheckBox *CheckMapout10;
   TCheckBox *CheckMapout9;
   TCheckBox *CheckMapout8;
   TCheckBox *CheckMapout7;
   TCheckBox *CheckMapout4;
   TCheckBox *CheckMapout0;
   TCheckBox *CheckMapout1;
   TCheckBox *CheckMapout2;
   TCheckBox *CheckMapout5;
   TCheckBox *CheckMapout3;
   TCheckBox *CheckMapout6;
   TStringGrid *MapsOutputBASIC;
   TTabSheet *TabMapOutMC;
   TPanel *GroupMapsoutMC;
   TCheckBox *CheckMapoutMC1;
   TCheckBox *CheckMapoutMC3;
   TCheckBox *CheckMapoutMC2;
   TCheckBox *CheckMapoutMC5;
   TCheckBox *CheckMapoutMC4;
   TCheckBox *CheckMapoutMC0;
   TCheckBox *CheckMapoutMC6;
   TStringGrid *MapsOutputMC;
   TTabSheet *TabMapOutNut;
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
   TTabSheet *TabMapOutNut2;
   TPanel *GroupMaspoutNutErosDepo;
   TLabel *Label84;
   TLabel *Label85;
   TCheckBox *CheckMapoutNut10;
   TCheckBox *CheckMapoutNut12;
   TCheckBox *CheckMapoutNut11;
   TCheckBox *CheckMapoutNut14;
   TCheckBox *CheckMapoutNut13;
   TCheckBox *CheckMapoutNut9;
   TStringGrid *MapsOutputNutErosDep;
   TTabSheet *TabMapOutGul;
   TPanel *GroupMapsoutGul;
   TCheckBox *CheckMapoutGul1;
   TCheckBox *CheckMapoutGul3;
   TCheckBox *CheckMapoutGul2;
   TCheckBox *CheckMapoutGul0;
   TCheckBox *CheckMapoutGul4;
   TStringGrid *MapsOutputGul;
   TTabSheet *TabSheetDrawmap;
   TSplitter *Splitter2;
   TScrollBox *ScrollImage;
   TPaintBox *MapImage;
   TToolBar *ToolBar1;
   TToolButton *ToolButtonZoomin;
   TToolButton *ToolButtonZoomout;
   TToolButton *ToolButtonZoom0;
   TToolButton *ToolButton2;
   TToolButton *ToolButtonGrid;
   TToolButton *ToolButtonPalette;
   TToolButton *ToolButtonClassType;
   TPanel *Panel11;
   TPaintBox *LegendBox;
   TImageList *ImageList2;
   TToolButton *ToolButtonDisplay;
   TPanel *Panel5;
   TRadioGroup *GroupDisplayMap;
   TMaskEdit *SOBEKDateString;
   TCGauge *CGauge;
   TStatusBar *StatusBar1;
   TBevel *Bevel19;
   TSpeedButton *RunFileButton;
   TLabel *Label86;
   TComboBox *ListRunfilesa;
   TToolButton *ToolButton1;
   TBevel *Bevel4;
   TLabel *Label87;
   TEdit *E_TotalRunoffName;
   TCheckBox *CheckSeparateOutput;
   TRadioGroup *E_OutputUnits;
   TJvLookOut *JvLookOut1;
   TJvLookOutPage *LookOutPage1;
   TJvLookOutPage *LookOutPage2;
   TJvLookOutPage *LookOutPage3;
   TJvLookOutPage *LookOutPage4;
   TLabel *Label44;
   TLabel *Label57;
   TBevel *Bevel7;
   TCheckBox *CheckBuffers;
   TCheckBox *CheckSedtrap;
   TCheckBox *CheckInfilGrass;
   TEdit *E_SedBulkDensity;
   TJvSpinEdit *CalibrateGrassN;
   TCheckBox *CheckNoErosion;
   TCheckBox *CheckIncludeChannel;
   TCheckBox *CheckChannelInfil;
   TCheckBox *CheckNoErosionOutlet;
   TCheckBox *CheckAltErosion;
   TCheckBox *CheckChannelBaseflow;
   TCheckBox *CheckSnowmelt;
   TCheckBox *CheckSimpleDepression;
   TBevel *Bevel33;
   TLabel *Label56;
   TBevel *Bevel26;
   TSpeedButton *LoadTableDir;
   TLabel *Label81;
   TComboBox *E_InfilMethoda;
   TCheckBox *CheckInfilCrust;
   TCheckBox *CheckInfilCompact;
   TCheckBox *CheckSubsoilDrainage;
   TCheckBox *CheckImpermeable;
   TEdit *E_SwatreDTSEC;
   TCheckBox *CheckDumphead;
   TCheckBox *CheckGeometric;
   TEdit *E_TableDir;
   TLabel *Label64;
   TLabel *Label22;
   TLabel *Label47;
   TLabel *Label80;
   TLabel *Label59;
   TJvSpinEdit *CalibrateKsat;
   TJvSpinEdit *CalibrateN;
   TJvSpinEdit *CalibrateChN;
   TJvSpinEdit *CalibrateChKsat;
   TJvSpinEdit *CalibrateSplashDelivery;
   TJvLookOutPage *LookOutPage5;
   TLabel *Label52;
   TRadioGroup *E_InterceptionLAIType;
   TCheckBox *CheckInterceptionLAI;
   TJvSpinEdit *FractionStemflow;
   TLabel *Label88;
   TBevel *Bevel14;
   TToolButton *ToolButton3;
   TLabel *Label89;
   TToolButton *ToolButton4;
   TJvCheckedMaskEdit *E_ClassifyMax;
   TToolButton *ToolButton5;
   TLabel *Label91;
   TLabel *Label93;
   TCheckBox *CheckHardsurface;
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
    void __fastcall ButtonDumpScreenClick(TObject *Sender);
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
    void __fastcall LisemstartButtonClick(TObject *Sender);
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
   void __fastcall CheckChannelBaseflowClick(TObject *Sender);
   void __fastcall CheckChannelInfilClick(TObject *Sender);
   void __fastcall CheckSnowmeltClick(TObject *Sender);
   void __fastcall MapsSnowmeltDblClick(TObject *Sender);
   void __fastcall LoadSnowmeltfileClick(TObject *Sender);
   void __fastcall SnowmeltViewButtonClick(TObject *Sender);
   void __fastcall MapsInfilKsatDblClick(TObject *Sender);
   void __fastcall MapsChannelBaseflowDblClick(TObject *Sender);
   void __fastcall CheckAllinChannelClick(TObject *Sender);
   void __fastcall ListRunfilesaClick(TObject *Sender);
   void __fastcall E_InfilMethodaClick(TObject *Sender);
   void __fastcall CheckSedtrapClick(TObject *Sender);
   void __fastcall DisplayHint(TObject *Sender);
   void __fastcall FormCreate(TObject *Sender);
   void __fastcall SetHelpIFace();
   void __fastcall ToolButtonHintClick(TObject *Sender);
   void __fastcall CheckWritePCRtimeplotClick(TObject *Sender);
   void __fastcall CheckSOBEKOutputClick(TObject *Sender);
   void __fastcall PageControlChange(TObject *Sender);
   void __fastcall MapImagePaint(TObject *Sender);
   void __fastcall ToolButtonGridClick(TObject *Sender);
   void __fastcall ToolButtonPaletteMouseDown(TObject *Sender,
          TMouseButton Button, TShiftState Shift, int X, int Y);
   void __fastcall ToolButtonNrClassClick(TObject *Sender);
   void __fastcall ToolButtonZoominClick(TObject *Sender);
   void __fastcall ToolButtonZoomoutClick(TObject *Sender);
   void __fastcall ToolButtonZoom0Click(TObject *Sender);
   void __fastcall MapImageMouseDown(TObject *Sender, TMouseButton Button,
          TShiftState Shift, int X, int Y);
   void __fastcall MapImageMouseMove(TObject *Sender, TShiftState Shift,
          int X, int Y);
   void __fastcall MapImageMouseUp(TObject *Sender, TMouseButton Button,
          TShiftState Shift, int X, int Y);
   void __fastcall ToolButtonDisplayClick(TObject *Sender);
   void __fastcall GroupDisplayMapClick(TObject *Sender);
   void __fastcall E_OutputUnitsClick(TObject *Sender);
   void __fastcall ToolButton1Click(TObject *Sender);
   void __fastcall ToolButton3Click(TObject *Sender);
   void __fastcall E_ClassifyMaxCheckClick(TObject *Sender);
   void __fastcall E_ClassifyMaxChange(TObject *Sender);
   void __fastcall CheckInterceptionLAIClick(TObject *Sender);
/*
   void __fastcall CalibrateKsatChange(TObject *Sender);
   void __fastcall CalibrateNChange(TObject *Sender);
   void __fastcall CalibrateChKsatChange(TObject *Sender);
   void __fastcall CalibrateChNChange(TObject *Sender);
   void __fastcall CalibrateSplashDeliveryChange(TObject *Sender);
*/
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
   LisThread *LisemRun;
   bool CheckPestout; // to see if lisem is called with for pest batch running
   float PestoutTimeinterval;


   AnsiString RunFilename;
   AnsiString RainfallDir;
   AnsiString SnowmeltDir;
   __fastcall TLisIFace(TComponent* Owner);
   __fastcall ~TLisIFace();

   int LisemType;

   void __fastcall TLisIFace::ResetIFace();
   void __fastcall TLisIFace::ResetInfil();
   void __fastcall TLisIFace::ResetMain();

   void __fastcall TLisIFace::DoMapDir(bool CheckName, bool CheckForce);

   TDateTime t_begin;
   TDateTime t_end;
   AnsiString runduration;
   double runlength;
   TMouse *mouse;

   void __fastcall TLisIFace::DumpScreenMaps();
   AnsiString MapdumpFilename;
   bool InitDumpMap;

   //Drawing
   TMapDraw *CsfMap;
   void __fastcall TLisIFace::Draw();
   void __fastcall TLisIFace::EraseScreen();
   void __fastcall TLisIFace::DrawMapandLegend();
   void __fastcall TLisIFace::ShowCellInfo(int X, int Y);
   void __fastcall TLisIFace::GetSelection();


   void __fastcall TLisIFace::StartDrawMap(MEM_HANDLE *M, double timestep);
   void __fastcall TLisIFace::CleanUpAll();
   void __fastcall TLisIFace::KillDrawMapStructure();
   int MapWidth, MapHeight;
   VPoint Pzoom;
   VPoint Cursor[2];
   EPoint Pdraw[2];
   REAL4 ClassMinV, oldClassMinV;
   REAL4 ClassMaxV, oldClassMaxV;
   bool zoomin;
   bool SwitchDisplayMaps;
   AnsiString DisplayMapItems;

};
//---------------------------------------------------------------------------
extern PACKAGE TLisIFace *LisIFace;
//---------------------------------------------------------------------------
#endif
