//--------------------------------------------------------------------------
#include <vcl.h>
#include <dir.h>
#pragma hdrstop

#include "iface.h"
#include "lismain.h"
#include "lishelp.h"
#include "lisabout.h"
#include "lisrunf.h"
#include "lisrainf.h"
//#include "lisdirv.h"
#include "lisstart.h"
#include "ifaceinit.h"
#include "ifacethread.h"
#include "lishint.h"


#include "mprolog.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma link "cgauges"
#pragma link "CSPIN"

#pragma link "JvExMask"
#pragma link "JvSpin"
#pragma link "JvBaseDlg"
#pragma link "JvBrowseFolder"
#pragma link "JvSelectDirectory"
#pragma resource "*.dfm"
TLisIFace *LisIFace;
//---------------------------------------------------------------------------
__fastcall TLisIFace::TLisIFace(TComponent* Owner)
    : TForm(Owner)
{
//VJ 050907 no more lisemwin.ini, just run batch files directly from the commandline
// syntax: lisemwin [run filename] [lisemtype 0-4] [pestout] [pest time interval in min]
    LisemType = -1;
    batchrun = false;
    CheckPestout = false;
    PestoutTimeinterval = 0;
    ListRunfilesa->Items->Clear();
    //ListRunfilesa->Items->Append(" ");

    // run directly from commandline, get name and lisemtype
//VJ 080925 behaviour when double clicking on run file
    if (ParamCount() == 1)
    {
  	     RunFilename = ParamStr(1);
        ListRunfilesa->Items->Strings[0] = RunFilename;
        batchrun = true;
        LisemType = 0;
    }
//VJ 050913 cleaned pestout possibility
    if (ParamCount() >= 2){
  	     RunFilename = ParamStr(1);
  	     if (FileExists(RunFilename)){
  	        LisemType = ParamStr(2).ToIntDef(-1);
//VJ 060109 fixed bug, changed to < basic instead of <= basic
           if (LisemType < LISEMBASIC || LisemType > LISEMGULLIES) {
         	  Application->MessageBox("Invalid Lisem type (must be 0-4)", "LISEM Error", MB_OK+MB_ICONERROR);
              RunFilename = "";
        	  } else {
             // runfile exists and lisemtype is good
			    batchrun = true;
	   	    ListRunfilesa->Items->Strings[0] = RunFilename;
//             CheckLisemType();
             //check for pestout
             if (ParamCount() >= 3)
				    CheckPestout = ParamStr(3).ToIntDef(0) == 1;
             if (CheckPestout) {
	             if (ParamCount() == 4)
   	             PestoutTimeinterval = ParamStr(4).ToDouble();
                else {
      	     	    Application->MessageBox("Pest output asked but no timeinterval given, syntax: \nlisemwin [runfilename] [lisem type] [pestout 0/1] [pest time in min]", "LISEM Error", MB_OK+MB_ICONERROR);
                   CheckPestout = false;
                }
             }
           }
        } else {
      	  AnsiString S = "File \"" + RunFilename + "\" not found.";
           Application->MessageBox(S.c_str(), "LISEM Error", MB_OK+MB_ICONERROR);
	        RunFilename = "";
        }
    }//paramcount >= 2

    ResetMain();

	 CheckLisemType();

    if (batchrun) {
        LoadRunFile(RunFilename);
        nrruns = 1;
        thisrun = 0;
        multipleRuns = false;
        ButtonRunProg->Enabled = true;
        ButtonStopProg->Enabled = true;
        ButtonPauseProg->Enabled = true;
	     Messages->Clear();
//        if (ParamCount() != 1)
  //      {
       	  PageControl->ActivePage = TabSheetTot;
   	     ThreadDone(NULL);
    //    }
    }

 //   Application->HintPause = 200;
 //   Application->HintHidePause = 50000;

    HorzScrollBar->Visible = true;
    VertScrollBar->Visible = true;
    HorzScrollBar->Range = 100;
    VertScrollBar->Range = 100;

    SetHelpIFace();

}
//---------------------------------------------------------------------------
__fastcall TLisIFace::~TLisIFace()
{
//    done = true;
    if (LisemRun)
    {
       LisemRun->RunDone = true;
       LisemRun->Terminate();
    }   
    delete StartForm;

}
//---------------------------------------------------------------------------
bool __fastcall TLisIFace::CheckError(AnsiString S)
{
     Beep();
     Application->MessageBox(S.c_str(), "LISEM Error", MB_OK+MB_ICONERROR);
     return false;
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::LoadRunFile(AnsiString RunFilename)
{
	 TStringList *lrun = new TStringList;

    SaveDialog->FileName = RunFilename;
    lrun->LoadFromFile(RunFilename);
    //RunForm->RunEdit->Lines->LoadFromFile(RunFilename);
    //gaat mis bij automatic run (pest) want runform bestaat nog niet
    //RunFileButton->Caption = "Acitive: " + ExtractFileName(RunFilename);

    if(lrun->Strings[0] == "[LISEM for WINDOWS run file]")
    {
       ReadNewRunfile(RunFilename);
    }
    else
    {
       if(Application->MessageBox("Loading a DOS type run file, will be saved as new type",
       "LISEM Warning", MB_OKCANCEL+MB_ICONWARNING)==IDOK)
         ParseOldRunfile();
    }
    delete lrun;
}

//---------------------------------------------------------------------------
AnsiString FromDouble(AnsiString txt)
{
    double v = txt.ToDouble();
    AnsiString S = v;
    return (S);
}
//---------------------------------------------------------------------------
AnsiString FromInt(AnsiString txt)
{
    int v = txt.ToInt();
    AnsiString S = v;
    return (S);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::ParseOldRunfile()
{

    TStrings *tmp = RunForm->RunEdit->Lines;
    int vi;

    E_MapDir->Text = CheckDir("Map Directory",tmp->Strings[0]);
    E_TableDir->Text = CheckDir("Table Directory",tmp->Strings[1]);
    E_ResultDir->Text = tmp->Strings[2];
    E_RainfallName->Text = ExtractFileName(tmp->Strings[3]);
    RainfallDir = ExtractFileDir(tmp->Strings[3]);
    RainfallDir = CheckDir("Rainfall Directory",RainfallDir);
    E_begintime->Text = FromDouble(tmp->Strings[4]);
    E_Endtime->Text = FromDouble(tmp->Strings[5]);
    E_Timestep->Text = FromDouble(tmp->Strings[6]);
    // 7 is not used
    vi = tmp->Strings[8].ToInt();
    if (vi >= 2)
       vi--;

    double value = tmp->Strings[9].ToDouble();
    value *= 24*60*60;
    if (value == 0)
         value = E_Timestep->Text.ToDouble()/10;
    E_SwatreDTSEC->Text = value;

    E_InfilMethoda->ItemIndex = vi;

    E_ErosionName->Text = tmp->Strings[17];
    E_DepositionName->Text = tmp->Strings[18];
    E_TotalName->Text = tmp->Strings[19];
    E_OutletName->Text = tmp->Strings[20];
//    E_Outlet1Name->Text = tmp->Strings[21];
//    E_Outlet2Name->Text = tmp->Strings[22];
//        tmp->Strings[24] = tmp->Strings[24].LowerCase();
    int i = 23;
    tmp->Strings[i] = tmp->Strings[i].LowerCase();
    if (tmp->Strings[i].AnsiPos("y")>0)
    {
       i++;
       E_OutputTimeStep->Checked = true;
       E_OutputTimeUser->Checked = false;
       E_OutputTimeSteps->Value = tmp->Strings[i].ToInt();
    }
    i++;
    tmp->Strings[i] = tmp->Strings[i].LowerCase();
    if (tmp->Strings[i].AnsiPos("y")>0)
    {
       E_OutputTimeUser->Checked = true;
       E_OutputTimeStep->Checked = false;
       E_OutputTimes->Lines->Clear();
       i++;
       while(tmp->Strings[i].ToInt() > 0)
       {
          E_OutputTimes->Lines->Append(tmp->Strings[i]);
          i++;
       }
    }
    i++;
    //E_RunoffName->Text = tmp->Strings[i];

}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::FileOpenClick(TObject *Sender)
{
  PageControl->ActivePage = TabSheetRun;

  OpenDialog->FileName = "";
  OpenDialog->InitialDir = E_Workdir->Text;
  OpenDialog->Options.Clear();
  OpenDialog->Options << ofAllowMultiSelect << ofFileMustExist << ofHideReadOnly << ofEnableSizing;
  OpenDialog->Filter = "Run files (*.run)|*.run|All files (*.*)|*.*";
  OpenDialog->FilterIndex = 1; // start the dialog showing all files
  if (OpenDialog->Execute())
  {
    RunForm->RunEdit->Clear();
    E_OutputTimes->Clear();

    ListRunfilesa->Items->Assign(OpenDialog->Files);
    ListRunfilesa->ItemIndex = 0;

    RunFilename = ListRunfilesa->Items->Strings[0];
      // strings are added to the top added
    LoadRunFile(RunFilename);

    nrruns = ListRunfilesa->Items->Count;
    thisrun = 0;
    multipleRuns = nrruns > 1;

    ButtonRunProg->Enabled = true;
    ButtonStopProg->Enabled = true;
    ButtonPauseProg->Enabled = true;
    ButtonStopAll->Enabled = nrruns > 1;
  }
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::FileSaveClick(TObject *Sender)
{
  if (SaveDialog->FileName != "")
    MakeNewRunfile(SaveDialog->FileName);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::FileSaveAsClick(TObject *Sender)
{
  SaveDialog->Title = "Save As";
  if (SaveDialog->Execute())
    MakeNewRunfile(SaveDialog->FileName);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::FileCloseClick(TObject *Sender)
{
  if (RunForm->RunEdit->Modified)
  {
     int result = Application->MessageBox(
      "The current file has changed. Save changes?","", MB_YESNOCANCEL+MB_ICONEXCLAMATION);
     if (result == IDYES) FileSaveClick(0);
     if (result == IDCANCEL) return;
  }
  if (RunForm->RunEdit->Lines->Count > 0)
     RunForm->RunEdit->Clear();
}

//---------------------------------------------------------------------------

void __fastcall TLisIFace::PrintButtonClick(TObject *Sender)
{
//     PrintScale = poPrintToFit;
 //    if (PrintDialog->Execute())
   //     Print();
}
//---------------------------------------------------------------------------
bool __fastcall TLisIFace::DumpScreen(bool saveio)
{
        bool doit = false;

        AnsiString savename = ExtractFileName(RunFilename);
        if (savename.IsEmpty())
        {
           if(Application->MessageBox("Load a runfile first","LISEM Warning", MB_OK+MB_ICONWARNING)==IDOK)
           return (doit);
        }

        PageControl->ActivePage = TabSheetTot;

        try{
          TJPEGImage *J = new TJPEGImage;
          J->Assign(GetFormImage());

          savename.Delete(savename.Length() - 3, 4);
          if (!E_ResultDir->Text.IsDelimiter("\\/",E_ResultDir->Text.Length()))
            E_ResultDir->Text = E_ResultDir->Text + "\\";
          SavePictureDialog->FileName = ExtractFilePath(E_ResultDir->Text)+savename+".jpg";

          if (saveio)
          {
             if (SavePictureDialog->Execute())
                J->SaveToFile(SavePictureDialog->FileName);
          }
          else
                J->SaveToFile(SavePictureDialog->FileName);
          delete J;
        }catch(...)
        {
           CheckError("Could not make jpeg");
        }
        doit = true;

        return (doit);
        // VJ 031218 boolean check dumpscreen to halt processes to finish dumpscreen in lismain
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::ButtonDumpScreenClick(TObject *Sender)
{
    bool check = DumpScreen(true);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::MakeIni(AnsiString name)
{
    TStringList *IniOp = new TStringList;

    IniOp->Clear();
    if (E_Workdir->Text.IsEmpty())
    {
       getcwd(workdir, 127);
       IniOp->Add("WorkDir="+(AnsiString)workdir);
    }
    else
       IniOp->Add("WorkDir="+E_Workdir->Text);

    IniOp->Add("LisemType=-1");
    IniOp->Add("[Startup]");
    IniOp->Add("Runfile0=");
    IniOp->Add("Runfile1=");
    IniOp->Add("Runfile2=");
    IniOp->Add("Runfile3=");

    IniOp->SaveToFile(name);

    delete IniOp;
}
//---------------------------------------------------------------------------
bool __fastcall TLisIFace::ProcessIni(AnsiString name)
{
    TStringList *IniOp = new TStringList;
    AnsiString _workdir;
    try
    {
        if (FileExists(name))
           IniOp->LoadFromFile(name);
        else
        {
           delete IniOp;
           return false;
        }
        _workdir = IniOp->Values[(AnsiString)"WorkDir"];
        _workdir = CheckDir("Work Directory",_workdir);
        if (DirectoryExists(_workdir))
        {
           chdir(_workdir.c_str());
           strcpy(workdir, _workdir.c_str());
           E_Workdir->Text = _workdir;
        }

        LisemType = IniOp->Values[(AnsiString)"LisemType"].ToIntDef(-1);
        if (LisemType == -1)
           return true;


        AnsiString S = IniOp->Values[(AnsiString)"Runfile0"];
        if (!S.IsEmpty()) ListRunfilesa->Items->Insert(0,S);
        S = IniOp->Values[(AnsiString)"Runfile1"];
        if (!S.IsEmpty()) ListRunfilesa->Items->Insert(0,S);
        S = IniOp->Values[(AnsiString)"Runfile2"];
        if (!S.IsEmpty()) ListRunfilesa->Items->Insert(0,S);
        S = IniOp->Values[(AnsiString)"Runfile3"];
        if (!S.IsEmpty()) ListRunfilesa->Items->Insert(0,S);

//VJ 050817 fixed batch runs
        if (ListRunfilesa->Items->Count > 0)
        {
            RunFilename = ListRunfilesa->Items->Strings[0];
            SaveDialog->FileName = RunFilename;
        }
        delete IniOp;
    }
    catch(...)
    {
    	delete IniOp;
        return false;
    }
	return true;

}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::AboutButtonClick(TObject *Sender)
{
       AboutBox->Version->Caption = ProgVersion;
       AboutBox->ProductDate->Caption = DateVersion;
       AboutBox->ShowModal();
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::StopAllButtonClick(TObject *Sender)
{
        Beep();
        LisemRun->RunDone = true;
        thisrun = nrruns;
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::MakeIniButtonClick(TObject *Sender)
{
    MakeIni((AnsiString)DirIni);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::CheckIncludeChannelClick(TObject *Sender)
{
    MapsChannels->Visible = CheckIncludeChannel->Checked;
    if (!CheckIncludeChannel->Checked)
       CheckChannelInfil->Checked = false;
//VJ 080412 baseflow an infil mutually exclusive
    if (!CheckChannelInfil->Checked)
       CheckChannelBaseflow->Checked = false;
    MapsChannelinfil->Visible = CheckChannelInfil->Checked && CheckIncludeChannel->Checked;
    MapsChannelBaseflow->Visible = CheckChannelBaseflow->Checked && CheckIncludeChannel->Checked;
    Label49->Visible = !CheckIncludeChannel->Checked;
    Label49->Caption = "Activate \"Include main channels\" on the start screen";
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::FormClose(TObject *Sender, TCloseAction &Action)
{
   if (LisemRun)
   {
      LisemRun->RunDone = true;
   //   done = true;
      Action=caNone;
      closeapp = true;
      Messages->Lines->Append("PROGRAM HALT, PLEASE WAIT ...");
   }
   else
      Action=caFree;   
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::HelpButtonClick(TObject *Sender)
{
       HelpForm->ShowModal();
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::SizeButtonClick(TObject *Sender)
{
    Width = wWidth;
    Height = wHeight;

    HorzScrollBar->Range = 100;
    VertScrollBar->Range = 100;
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::FormResize(TObject *Sender)
{
//    HorzScrollBar->Visible = true;
//    VertScrollBar->Visible = true;
    if (Width < wWidth) HorzScrollBar->Range = wWidth;
       else HorzScrollBar->Range = 100;
    if (Height < wHeight) VertScrollBar->Range = wHeight;
       else VertScrollBar->Range = 100;
    Left = (Screen->Width-Width)/2;
    Top = (Screen->Height-Height)/2;
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::CheckNoErosionClick(TObject *Sender)
{
      //VJ 080217 changed behaviour of sedigraph on screen
        //Panel3->Visible=!CheckNoErosion->Checked;
        if (!CheckNoErosion->Checked)
        {
            ShowQsed->Checked = false;
        }
        ShowSedConcentration->Checked = (!CheckNoErosion->Checked);

        if (CheckNoErosion->Checked && CheckMulticlass)
        {
           CheckError("Simulation of multiple class sediment or nutrient losses only with erosion");
           CheckNutrients = false;
           CheckMulticlass = false;
        }
        MapsErosion->Visible = !CheckNoErosion->Checked;
        //Label92->Visible = !CheckNoErosion->Checked;
        Label60->Visible = CheckNoErosion->Checked;
        Label60->Caption = "Deactivate \"Runoff only\" calculation on Start screen";
//        CheckMapsEnabled();
}
//---------------------------------------------------------------------------
AnsiString __fastcall TLisIFace::CheckDir(AnsiString Comment, AnsiString dir)
{
    if (dir != "")
    {
      AnsiString OldDir = dir;
//VJ 030626 changed
      if (dir.c_str()[dir.Length()-1] != '\\')
         dir=dir+"\\";
//VJ 080613 added expandfilename         
      dir = ExtractFileDir(ExpandFileName(dir));
      if (dir.c_str()[dir.Length()-1] != '\\')
         dir=dir+"\\";
      if (!DirectoryExists(dir))
      {
          AnsiString S = Comment+": \"" + dir + "\" does not exist";
          Application->MessageBox(S.c_str(), "LISEM Error", MB_OK+MB_ICONERROR);
          dir = "DIR_NOT_EXIST";
      }
    }
    return(dir);
}
//---------------------------------------------------------------------------
AnsiString __fastcall TLisIFace::CheckFile(AnsiString dir)
{
    if (dir != "")
    {
      dir = ExpandFileName(dir);
      if (!FileExists(dir))
      {
          AnsiString S = "File \"" + dir + "\" not found";
          Application->MessageBox(S.c_str(), "LISEM Error", MB_OK+MB_ICONERROR);
          dir = "";
      }
      else
        dir = ExtractFileName(dir);
    }
    return(dir);
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::E_TableDirClick(TObject *Sender)
{
    E_TableDir->Text = CheckDir("Table Directory",E_TableDir->Text);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::E_RainfallNameClick(TObject *Sender)
{
    E_RainfallName->Text = CheckFile(RainfallDir+E_RainfallName->Text);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::RunFileButtonClick(TObject *Sender)
{
       if (ListRunfilesa->Items->Count == 0 || ListRunfilesa->Items->Strings[0].IsEmpty())
          (void) CheckError("Load RUN file first !");
       else
       {
          RunForm->RunEdit->Lines->LoadFromFile(RunFilename);
          RunForm->Caption = RunFilename;
          RunForm->ShowModal();
       }
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::E_OutputTimeStepClick(TObject *Sender)
{
     E_OutputTimeStep->Checked = !E_OutputTimeUser->Checked;
     E_OutputTimeSteps->Enabled = E_OutputTimeStep->Checked;
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::E_OutputTimeUserClick(TObject *Sender)
{
     E_OutputTimeUser->Checked = !E_OutputTimeStep->Checked;
     E_OutputTimeSteps->Enabled = E_OutputTimeStep->Checked;
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::RainViewButtonClick(TObject *Sender)
{
   if (RainViewButton->Tag == 1)
   {
     RainViewButton->Tag = 0;
     AnsiString filename = CheckDir("Rainfall Directory",RainfallDir);
     if (FileExists(filename))
     {
       RainForm->RainViewList->Lines->Clear();
       RainForm->RainViewList->Lines->LoadFromFile(filename);
       RainForm->ShowModal();
     }
     else
     {
       CheckError("File \""+filename+"\" not found");
     }
   }
   else
   {
     AnsiString filename = CheckDir("Rainfall Directory",RainfallDir);
     filename = filename+E_RainfallName->Text;
     if (FileExists(filename))
     {
       RainForm->RainViewList->Lines->Clear();
       RainForm->RainViewList->Lines->LoadFromFile(filename);
       RainForm->Caption = E_RainfallName->Text;
       RainForm->ShowModal();
     }
     else
     {
       CheckError("File \""+filename+"\" not found");
     }
   }
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::GetDirEdit(TEdit *E)
{
     ActiveControl = E;
/*
     if (E->Text.IsEmpty())
         DirView->DirectoryListBox->Directory = (AnsiString)workdir;
     else
         DirView->DirectoryListBox->Directory = E->Text;

     DirView->ShowModal();
     DirView->DirectoryListBox->OpenCurrent();
     E->Text = DirView->DirectoryListBox->Directory;
*/
     if (E->Text.IsEmpty())
         JvGetDirDialog->Directory = (AnsiString)workdir;
     else
         JvGetDirDialog->Directory = E->Text;
     JvGetDirDialog->Execute();
     E->Text = JvGetDirDialog->Directory;
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::LoadRainfileClick(TObject *Sender)
{
    OpenDialog->FileName = "";
    OpenDialog->Options.Clear();
    OpenDialog->Filter = "Rainfall files (*.txt; *.dat; *.prn)|*.txt;*.dat;*.prn|All files (*.*)|*.*";
    OpenDialog->FilterIndex = 1;
    OpenDialog->InitialDir = RainfallDir;
    if(OpenDialog->Execute())
    {
       E_RainfallName->Text = ExtractFileName(OpenDialog->FileName);
       RainfallDir = ExtractFileDir(OpenDialog->FileName)+"\\";
       strcpy(RainfallNamePath,OpenDialog->FileName.c_str());
    }
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::LoadTableDirClick(TObject *Sender)
{
     GetDirEdit(E_TableDir);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::LoadMapDirClick(TObject *Sender)
{
    GetDirEdit(E_MapDir);
    DoMapDir(true, true);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::LoadWorkDirClick(TObject *Sender)
{
     GetDirEdit(E_Workdir);
     chdir(E_Workdir->Text.c_str());
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::LoadResultDirClick(TObject *Sender)
{
//VJ 030626 Changed
   GetDirEdit(E_ResultDir);

   if (!DirectoryExists(E_ResultDir->Text))
   {
/*
       int button = Application->MessageBox("Directory does not exist, create ?",
       "LISEM Warning", MB_OKCANCEL+MB_ICONWARNING+MB_DEFBUTTON1);
       if (button == IDOK)
*/
       ForceDirectories(E_ResultDir->Text);
       return;
    }
    E_ResultDirDblClick(Sender);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::E_InfilMethodClick(TObject *Sender)
{
     ResetInfil();
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::EraseListButtonClick(TObject *Sender)
{
   if (!LisemRun)// || done)
   {
       ResetIFace();
       //RunFileButton->Caption = "Show active";
   }
   else
      Application->MessageBox("Cannot apply changes while running, stop first!",
       "LISEM Warning", MB_OK+MB_ICONWARNING);

}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::GetMapname(TStrings *S)
{
 //   PageControl->ActivePage = TabMapNames;

    OpenDialog->FileName = S->Strings[1];
//VJ 040328 corrected small bug, wrong string
    if (S->Strings[3] != "")
       OpenDialog->InitialDir = S->Strings[3];
    else
       OpenDialog->InitialDir = E_MapDir->Text;
    OpenDialog->Options.Clear();
    OpenDialog->Options << ofFileMustExist << ofHideReadOnly;
    OpenDialog->Filter = "PCRaster maps (*.map; *.csf)|*.map;*.csf|All files (*.*)|*.*";
    OpenDialog->FilterIndex = 1; // start the dialog showing all files
    if (OpenDialog->Execute())
    {
        S->Strings[1] = ExtractFileName(OpenDialog->FileName);
        S->Strings[3] = ExtractFileDir(OpenDialog->FileName);
        S->Strings[4] = OpenDialog->FileName;
    }
}


//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsCatchmentDblClick(TObject *Sender)
{
        GetMapname(MapsCatchment->Rows[MapsCatchment->Row]);
}


void __fastcall TLisIFace::MapsErosionDblClick(TObject *Sender)
{
        GetMapname(MapsErosion->Rows[MapsErosion->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsSurfaceDblClick(TObject *Sender)
{
        GetMapname(MapsSurface->Rows[MapsSurface->Row]);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::MapsInfilGA1DblClick(TObject *Sender)
{
        GetMapname(MapsInfilGA1->Rows[MapsInfilGA1->Row]);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::MapsWheeltrackDblClick(TObject *Sender)
{
        GetMapname(MapsWheeltrack->Rows[MapsWheeltrack->Row]);
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::MapsChannelsDblClick(TObject *Sender)
{
        GetMapname(MapsChannels->Rows[MapsChannels->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsChannelinfilDblClick(TObject *Sender)
{
        GetMapname(MapsChannelinfil->Rows[MapsChannelinfil->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsOutputBASICDblClick(TObject *Sender)
{
        GetMapname(MapsOutputBASIC->Rows[MapsOutputBASIC->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsInfilExtraDblClick(TObject *Sender)
{
        GetMapname(MapsInfilExtra->Rows[MapsInfilExtra->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsInfilHoltanDblClick(TObject *Sender)
{
                GetMapname(MapsInfilHoltan->Rows[MapsInfilHoltan->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsInfilGA2DblClick(TObject *Sender)
{
                GetMapname(MapsInfilGA2->Rows[MapsInfilGA2->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsInfilSwatreDblClick(TObject *Sender)
{
        GetMapname(MapsInfilSwatre->Rows[MapsInfilSwatre->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsTextureDblClick(TObject *Sender)
{
         GetMapname(MapsTexture->Rows[MapsTexture->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsLanduseDblClick(TObject *Sender)
{
        GetMapname(MapsLanduse->Rows[MapsLanduse->Row]);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::CheckMapsEnabled()
{
   for (int i = 0; i < ComponentCount; i++)
    if (Components[i]->Name.SubString(0,4) == "Maps")
    {
         TStringGrid* S =(TStringGrid *)Components[i];
         if (!S->Enabled)
         {
            S->Font->Color = clGrayText;
            S->Color = clBtnFace;
         }
         else
         {
            S->Font->Color = clBlack;
            S->Color = clWindow;
         }
    }
         if (!TextureClass->Enabled)
         {
            TextureClass->Font->Color = clGrayText;
            TextureClass->Color = clBtnFace;
         }
         else
         {
            TextureClass->Font->Color = clBlack;
            TextureClass->Color = clWindow;
         }

}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckInfilCrustClick(TObject *Sender)
{
    MapsInfilExtra->Visible = CheckInfilCrust->Checked;
//    MapsInfilExtra->Enabled = CheckInfilCrust->Checked;
//    CheckMapsEnabled();
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::CheckInfilCompactClick(TObject *Sender)
{
    MapsInfilExtra->Visible = CheckInfilCompact->Checked;
//    CheckMapsEnabled();
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::TextureClassCheck()
{
   for (int i = 0; i < 6; i++)
   {
       if (TextureClass->Cells[0][i+1].IsEmpty())
          return;
       sscanf(TextureClass->Cells[0][i+1].c_str(),"%f",&ClassMu[i]);
       if (ClassMu[i] == 0)
       {
           CheckError("Class boundary must be larger than 0 mu.");
           TextureClass->Cells[0][i+1] = "";
           return;
       }
   }
   for (int j = 1; j < 6; j++)
   for (int i = 0; i < j; i++)
   {
       if (ClassMu[j] <= ClassMu[i])
       {
           CheckError("Class boundary smaller than lower boundary");
           TextureClass->Cells[0][j+1] = "";
//           ClassMu[j] = 0;
           return;
       }
   }
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::TextureClassClick(TObject *Sender)
{
      TextureClassCheck();
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::TextureClassKeyDown(TObject *Sender, WORD &Key,
      TShiftState Shift)
{
   if (Key == VK_RETURN)
      TextureClassCheck();
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::ButtonRestoreClick(TObject *Sender)
{
    DoMapDir(false, true);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::LisemstartButtonClick(TObject *Sender)
{
    if (!LisemRun)
    {
      ResetMain();
      ResetIFace();

      LisemType = LISEMBASIC;

      CheckLisemType();
      if (!StartForm)
         StartForm = new TStartForm(static_cast<void *>(NULL));
      else
         StartForm->Show();
      LisIFace->Enabled = false;
    }
    else
      Application->MessageBox("Please stop run first","LISEM message",MB_OK+MB_ICONEXCLAMATION);

}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::CheckLisemType()
{
   CheckWheelAsChannel = false;
   CheckMulticlass = false;
   CheckNutrients = false;
   CheckGullies = false;
   MapOutControl->Visible = true;
   switch (LisemType) {
     case LISEMBASIC :  //MapOutControl->Visible = false;
                   CheckWheelAsChannelModel();
                   CheckMulticlassModel();
                   CheckNutrientsModel();
                   CheckGulliesModel();
                   break;
     case LISEMWHEELTRACKS : CheckWheelAsChannel = true; CheckWheelAsChannelModel();break;
     case LISEMMULTICLASS : CheckMulticlass = true; CheckMulticlassModel(); break;
     case LISEMNUTRIENTS : CheckNutrients = true; CheckMulticlass = true; CheckMulticlassModel();CheckNutrientsModel();break;
     case LISEMGULLIES : CheckGullies = true;CheckGulliesModel(); break;
   }
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckWheelAsChannelModel()
{
   TabWheeltracks->TabVisible=CheckWheelAsChannel;
   MapsWheeltrack->Enabled = CheckWheelAsChannel;
   CheckMapsEnabled();
   if (CheckWheelAsChannel)
      Panel4->Caption = "Wheeltrack and Chan. erosion";
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckMulticlassModel()
{

       TabMulticlass->TabVisible = CheckMulticlass;
       MapsOutputMC->Visible = false;//CheckMulticlass;
       GroupMapsoutMC->Visible = CheckMulticlass;
       TabMapOutMC->TabVisible = CheckMulticlass;
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckNutrientsModel()
{
       TabMulticlass->TabVisible = CheckNutrients;
       TabNutrients->TabVisible = CheckNutrients;
       MapsOutputNut->Visible = false;//CheckNutrients;
       GroupMapsoutNut->Visible = CheckNutrients;
       TabMapOutNut->TabVisible = CheckNutrients;
       TabMapOutNut2->TabVisible = CheckNutrients;

       MapsOutputMC->Visible = false;//CheckNutrients;//false;
       GroupMapsoutMC->Visible = CheckNutrients;//false;
       TabMapOutMC->TabVisible = CheckNutrients;//false;
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::CheckGulliesModel()
{
     TabGullies->TabVisible = CheckGullies;
     MapsOutputGul->Visible = false;//CheckGullies;
     GroupMapsoutGul->Visible = CheckGullies;
     TabMapOutGul->TabVisible = CheckGullies;
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsInfilSmithDblClick(TObject *Sender)
{
        GetMapname(MapsInfilSmith->Rows[MapsInfilSmith->Row]);
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::MapsInfilMorelDblClick(TObject *Sender)
{
        GetMapname(MapsInfilMorel->Rows[MapsInfilMorel->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckInfilGrassClick(TObject *Sender)
{
    MapsInfilExtra->Visible = CheckInfilGrass->Checked;
    E_ManningsNGrass->Enabled = CheckInfilGrass->Checked;
    Label44->Enabled = CheckInfilGrass->Checked;
}

//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsGullyDblClick(TObject *Sender)
{
        GetMapname(MapsGully->Rows[MapsGully->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsNutsPDblClick(TObject *Sender)
{
        GetMapname(MapsNutsP->Rows[MapsNutsP->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsNutsNO3DblClick(TObject *Sender)
{
        GetMapname(MapsNutsNO3->Rows[MapsNutsNO3->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsNutsNH4DblClick(TObject *Sender)
{
        GetMapname(MapsNutsNH4->Rows[MapsNutsNH4->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsNutsBDDblClick(TObject *Sender)
{
        GetMapname(MapsNutsBD->Rows[MapsNutsBD->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::E_MapDirDblClick(TObject *Sender)
{
//VJ 030626 changed behaviour
     DoMapDir(true, false);
/*
     E_MapDir->Text = CheckDir("Map Directory",E_MapDir->Text);
     if (OldMapDir != E_MapDir->Text && E_MapDir->Text != "DIR_NOT_EXIST")
     {
        CheckAdjustMapDirectoryName = true;
        ButtonRestoreClick(Sender);
     }
     if (E_MapDir->Text == "DIR_NOT_EXIST")
        E_MapDir->Text = OldMapDir;
     OldMapDir = E_MapDir->Text;
*/
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::E_WorkdirDblClick(TObject *Sender)
{
//VJ 030626 changed behaviour
   E_Workdir->Text = CheckDir("Work Directory",E_Workdir->Text);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::E_ResultDirDblClick(TObject *Sender)
{
//VJ 030626 changed behaviour
//VJ 091216 empty result dir check
   if (!E_ResultDir->Text.IsEmpty() && !DirectoryExists(E_ResultDir->Text))
   {
       AnsiString S = "Result Directory: \"" + E_ResultDir->Text + "\" does not exist, create it?";
       if (Application->MessageBox(S.c_str(), "LISEM Warning", MB_OKCANCEL+MB_ICONERROR) == IDOK)
          ForceDirectories(E_ResultDir->Text);
       else
          E_ResultDir->Text = "";
   }

   if (!E_ResultDir->Text.IsEmpty());
     InitOutMapNames();
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::E_ResultDirKeyPress(TObject *Sender, char &Key)
{
//VJ 030626 changed behaviour
    if (Key == VK_RETURN)
        E_ResultDirDblClick(Sender);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::E_MapDirKeyPress(TObject *Sender, char &Key)
{
//VJ 030626 changed behaviour
    if (Key == VK_RETURN)
        E_MapDirDblClick(Sender);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::E_WorkdirKeyPress(TObject *Sender, char &Key)
{
//VJ 030626 changed behaviour
    if (Key == VK_RETURN)
        E_WorkdirDblClick(Sender);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::DoMapDir(bool CheckName, bool CheckForce)
{
   int res = 0;
    if (CheckName)
    {
    if (OldMapDir != E_MapDir->Text)
       res = Application->MessageBox("Caution: new map Directory name added to old mapnames","LISEM message",MB_OKCANCEL+MB_ICONEXCLAMATION);
    }
    else
	    res = Application->MessageBox("Caution: default names and values are restored with map directory name","LISEM message",MB_OKCANCEL+MB_ICONEXCLAMATION);
//VJ 050907 added cancel option
    if (res == 2)
      	return;
     E_MapDir->Text = CheckDir("Map Directory",E_MapDir->Text);
     if (CheckForce || (OldMapDir != E_MapDir->Text && E_MapDir->Text != "DIR_NOT_EXIST"))
     {
        ForceMapdir = CheckForce;
        CheckAdjustMapDirectoryName = CheckName;
        InitMapNames();
     }
     if (E_MapDir->Text == "DIR_NOT_EXIST")
        E_MapDir->Text = OldMapDir;
     OldMapDir = E_MapDir->Text;
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::E_ResultDirExit(TObject *Sender)
{
        E_ResultDirDblClick(Sender);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::E_MapDirExit(TObject *Sender)
{
        E_MapDirDblClick(Sender);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::E_WorkdirExit(TObject *Sender)
{
        E_WorkdirDblClick(Sender);
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::MapsGullyInitDblClick(TObject *Sender)
{
        GetMapname(MapsGullyInit->Rows[MapsGullyInit->Row]);
//VJ 040331 Included init gully dimensions
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsBuffersDblClick(TObject *Sender)
{
        GetMapname(MapsBuffers->Rows[MapsBuffers->Row]);
// VJ 040514 include buffers
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::CheckBuffersClick(TObject *Sender)
{
//VJ 040823 include buffers
   Label55->Visible = CheckBuffers->Checked;
   MapsBuffers->Visible = CheckBuffers->Checked;
 //  CheckSedtrap->Checked = !CheckBuffers->Checked;
   //output screen
   Label61->Enabled = CheckBuffers->Checked;
   Label63->Enabled = CheckBuffers->Checked;
   Label75->Enabled = CheckBuffers->Checked;
   LO_buffervol->Enabled = CheckBuffers->Checked;
   LO_buffersedvol->Enabled = CheckBuffers->Checked;

}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::CheckSubsoilDrainageClick(TObject *Sender)
{
     //VJ 050812 included impermeable and drainage checks
     if (E_InfilMethoda->ItemIndex != INFIL_GREENAMPT && E_InfilMethoda->ItemIndex != INFIL_GREENAMPT2)
        CheckSubsoilDrainage->Checked = false;
    MapsInfilDrainage->Visible = CheckSubsoilDrainage->Checked;
    if (CheckSubsoilDrainage->Checked)
       CheckImpermeable->Checked = false;
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckImpermeableClick(TObject *Sender)
{
     //VJ 050812 included impermeable and drainage checks 
    if (CheckImpermeable->Checked)
       CheckSubsoilDrainage->Checked = false;
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsInfilDrainageDblClick(TObject *Sender)
{
     GetMapname(MapsInfilDrainage->Rows[MapsInfilDrainage->Row]);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckChannelBaseflowClick(TObject *Sender)
{
//VJ 080412 baseflow is mutually exclusive with infil
//VJ 080614 why?????
    if (!CheckIncludeChannel->Checked)
       CheckChannelBaseflow->Checked = false;
 //   if (CheckChannelBaseflow->Checked)
   //    CheckChannelInfil->Checked = false;
    MapsChannelBaseflow->Visible = CheckChannelBaseflow->Checked && CheckIncludeChannel->Checked;
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckChannelInfilClick(TObject *Sender)
{
//VJ 080412 baseflow is mutually exclusive with infil
//VJ 080616 why??????
    if (!CheckIncludeChannel->Checked)
       CheckChannelInfil->Checked = false;
 //   if (CheckChannelInfil->Checked)
   //    CheckChannelBaseflow->Checked = false;
    MapsChannelinfil->Visible = CheckChannelInfil->Checked && CheckIncludeChannel->Checked;
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckSnowmeltClick(TObject *Sender)
{
//VJ 080423 include snowmelt
   Label73->Enabled = CheckSnowmelt->Checked;
   E_SnowmeltName->Enabled = CheckSnowmelt->Checked;
   E_SnowmeltName->Visible = CheckSnowmelt->Checked;
   SnowmeltViewButton->Enabled = CheckSnowmelt->Checked;
   LoadSnowmeltfile->Enabled = CheckSnowmelt->Checked;

   Label67->Visible = CheckSnowmelt->Checked;
   MapsSnowmelt->Visible = CheckSnowmelt->Checked;
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::MapsSnowmeltDblClick(TObject *Sender)
{
        GetMapname(MapsSnowmelt->Rows[MapsSnowmelt->Row]);
//VJ 0080423 snowmelt
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::LoadSnowmeltfileClick(TObject *Sender)
{
    OpenDialog->FileName = "";
    OpenDialog->Options.Clear();
    OpenDialog->Filter = "Snowmelt files (*.txt; *.dat; *.prn)|*.txt;*.dat;*.prn|All files (*.*)|*.*";
    OpenDialog->FilterIndex = 1;
    OpenDialog->InitialDir = SnowmeltDir;
    if(OpenDialog->Execute())
    {
       E_SnowmeltName->Text = ExtractFileName(OpenDialog->FileName);
       SnowmeltDir = ExtractFileDir(OpenDialog->FileName)+"\\";
       strcpy(SnowmeltNamePath,OpenDialog->FileName.c_str());
    }
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::SnowmeltViewButtonClick(TObject *Sender)
{
   if (RainViewButton->Tag == 1)
   {
     RainViewButton->Tag = 0;
     AnsiString filename = CheckDir("Snowmelt Directory",SnowmeltDir);
     if (FileExists(filename))
     {
       RainForm->RainViewList->Lines->Clear();
       RainForm->RainViewList->Lines->LoadFromFile(filename);
       RainForm->ShowModal();
     }
     else
     {
       CheckError("File \""+filename+"\" not found");
     }
   }
   else
   {
     AnsiString filename = CheckDir("Snowmelt Directory",SnowmeltDir);
     filename = filename+E_SnowmeltName->Text;
     if (FileExists(filename))
     {
       RainForm->RainViewList->Lines->Clear();
       RainForm->RainViewList->Lines->LoadFromFile(filename);
       RainForm->Caption = E_SnowmeltName->Text;
       RainForm->ShowModal();
     }
     else
     {
       CheckError("File \""+filename+"\" not found");
     }
   }

}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::MapsInfilKsatDblClick(TObject *Sender)
{
     GetMapname(MapsInfilKsat->Rows[MapsInfilKsat->Row]);
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::MapsChannelBaseflowDblClick(TObject *Sender)
{
        GetMapname(MapsChannelBaseflow->Rows[MapsChannelBaseflow->Row]);
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::CheckAllinChannelClick(TObject *Sender)
{
    if (!CheckIncludeChannel->Checked)
       CheckAllinChannel->Checked = false;
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::GetNewRunfile()
{
   LisIFace->RunFilename = LisIFace->ListRunfilesa->Items->Strings[thisrun];
   LisIFace->SaveDialog->FileName = LisIFace->RunFilename;
   RunForm->RunEdit->Lines->LoadFromFile(LisIFace->RunFilename);
   LisIFace->ReadNewRunfile(LisIFace->RunFilename);
   LisIFace->PageControl->ActivePage = LisIFace->TabSheetTot;
   LisIFace->LabelRunfile->Caption = "Active run file: " + LisIFace->RunFilename;
}
//---------------------------------------------------------------------------



void __fastcall TLisIFace::ListRunfilesaClick(TObject *Sender)
{
         thisrun = ListRunfilesa->ItemIndex;
         RunFilename = ListRunfilesa->Items->Strings[thisrun];
         LoadRunFile(RunFilename);
}
//---------------------------------------------------------------------------



void __fastcall TLisIFace::E_InfilMethodaClick(TObject *Sender)
{
     ResetInfil();
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckSedtrapClick(TObject *Sender)
{

//VJ 090527 include sedtrap
   Label55->Visible = CheckSedtrap->Checked;
   MapsBuffers->Visible = CheckSedtrap->Checked;
   CheckBuffers->Checked = CheckSedtrap->Checked;
   //output screen
   Label61->Enabled = CheckSedtrap->Checked;
   Label63->Enabled = CheckSedtrap->Checked;
   Label75->Enabled = CheckSedtrap->Checked;
   LO_buffervol->Enabled = CheckSedtrap->Checked;
   LO_buffersedvol->Enabled = CheckSedtrap->Checked;

}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::DisplayHint(TObject *Sender)
{
    HintForm->expl->Text = GetLongHint(Application->Hint);
}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::FormCreate(TObject *Sender)
{
    Application->OnHint = DisplayHint;
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::ToolButtonHintClick(TObject *Sender)
{
   if (ToolButtonHint->Down)
   {
   HintForm->Show();
   HintForm->expl->Color = 0x00bbffff;
   }
   else
   HintForm->Close();
}
//---------------------------------------------------------------------------


void __fastcall TLisIFace::CheckWritePCRtimeplotClick(TObject *Sender)
{
   if(CheckWritePCRtimeplot->Checked)
   CheckSOBEKOutput->Checked = false;

}
//---------------------------------------------------------------------------

void __fastcall TLisIFace::CheckSOBEKOutputClick(TObject *Sender)
{
   if(CheckSOBEKOutput->Checked)
   CheckWritePCRtimeplot->Checked = false;

}
//---------------------------------------------------------------------------

