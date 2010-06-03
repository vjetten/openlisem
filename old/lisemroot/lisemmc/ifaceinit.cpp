//---------------------------------------------------------------------------

#pragma hdrstop

#include "ifaceinit.h"
#include "iface.h"
#include "lisstart.h"



//---------------------------------------------------------------------------

#pragma package(smart_init)


void __fastcall TLisIFace::ResetMain()
{
    Caption = (AnsiString)ProgName + (AnsiString)" - " + (AnsiString)ProgVersion;

    multipleRuns = false;
    take5 = false;
    closeapp = false;
    ForceMapdir = false;
    LisemRun = NULL;
    nrruns = 0;
    thisrun = 0;
    Width = wWidth;
    Height = wHeight;
    PageControl->ActivePage = TabSheetRun;
    Application->Title = "LISEMWIN "+(AnsiString)ProgVersion;
    //namelist = NULL;

    ClassMu[0] = 2;
    ClassMu[1] = 16;
    ClassMu[2] = 32;
    ClassMu[3] = 53;
    ClassMu[4] = 75;
    ClassMu[5] = 105;

    InitMapNames();

    TabWheeltracks->TabVisible=false;
    TabMulticlass->TabVisible=false;
    TabNutrients->TabVisible=false;
    TabGullies->TabVisible=false;

    MapsOutputNut->Visible = false;
    GroupMapsoutNut->Visible = false;
    TabMapOutNut->TabVisible = false;

    MapsOutputGul->Visible = false;
    GroupMapsoutGul->Visible = false;
    TabMapOutGul->TabVisible = false;

    MapsOutputMC->Visible = false;
    GroupMapsoutMC->Visible = false;
    TabMapOutMC->TabVisible = false;

    CheckWheelAsChannel = false;
    CheckGullies = false;
    CheckNutrients = false;
    CheckMulticlass = false;

    E_InfilMethodClick(NULL);
    CheckIncludeChannelClick(NULL);
    Names->ActivePage = Area;

    if (LisemType== -1)
    {
        LisIFace->Enabled = false;
        StartForm = new TStartForm(static_cast<void *>(NULL));
        //StartForm->Show();
//    Application->ProcessMessages(); // this is necessary to show the Splash now

    }
    else
    {
        LisIFace->Enabled = true;
    }

    CheckAdjustMapDirectoryName = false;

//VJ 040823 include buffers
   Label55->Visible = CheckBuffers->Checked;
   MapsBuffers->Visible = CheckBuffers->Checked;
   //output screen
   Label61->Enabled = CheckBuffers->Checked;
   Label63->Enabled = CheckBuffers->Checked;
   Label75->Enabled = CheckBuffers->Checked;
   LO_buffervol->Enabled = CheckBuffers->Checked;
   LO_buffersedvol->Enabled = CheckBuffers->Checked;

//VJ 080423 snowmelt
   Label73->Enabled = CheckSnowmelt->Checked;
   E_SnowmeltName->Enabled = CheckSnowmelt->Checked;
   E_SnowmeltName->Visible = CheckSnowmelt->Checked;
   Label67->Visible = CheckSnowmelt->Checked;
   MapsSnowmelt->Visible = CheckSnowmelt->Checked;

   SnowmeltViewButton->Enabled = CheckSnowmelt->Checked;
   LoadSnowmeltfile->Enabled = CheckSnowmelt->Checked;
}

//---------------------------------------------------------------------------

void __fastcall TLisIFace::ResetInfil()
{
      int _top = 22, _left = 8;

     //VJ 050812 included impermeable check here
     if (E_InfilMethod->ItemIndex == INFIL_NONE){
        CheckImpermeable->Checked = false;
     }

     //VJ 050812 included drainage check here
     if (E_InfilMethod->ItemIndex != INFIL_GREENAMPT && E_InfilMethod->ItemIndex != INFIL_GREENAMPT2)
        CheckSubsoilDrainage->Checked = false;

     Label58->Visible = E_InfilMethod->ItemIndex == 0;
     Label58->Caption = "Choose an infiltration method on the Start screen";

     switch (E_InfilMethod->ItemIndex) {
       case INFIL_NONE : LabelInfilName->Caption = "Infiltration maps";break;
       case INFIL_SWATRE : LabelInfilName->Caption = "SWATRE";break;
       case INFIL_HOLTAN : LabelInfilName->Caption = "Holtan";break;
       case INFIL_GREENAMPT : LabelInfilName->Caption = "Green && Ampt";break;
       case INFIL_GREENAMPT2 : LabelInfilName->Caption = "Green && Ampt";break;
       case INFIL_MOREL : LabelInfilName->Caption = "Morel-Seytoux && Verdin";break;
       case INFIL_SMITH : LabelInfilName->Caption = "Smith && Parlange";break;
       case INFIL_KSAT : LabelInfilName->Caption = "Ksat"; break;
       }
//VJ 080521 replaced values by defines to avoid morel and smith
/*
#define INFIL_NONE 0
#define INFIL_SWATRE 1
#define INFIL_HOLTAN 2
#define INFIL_GREENAMPT 3
#define INFIL_GREENAMPT2 4
#define INFIL_KSAT 5
#define INFIL_MOREL 21
#define INFIL_SMITH 22
*/
     MapsInfilSwatre->Visible = E_InfilMethod->ItemIndex == INFIL_SWATRE;
     GroupSwatreOptions->Visible = E_InfilMethod->ItemIndex == INFIL_SWATRE;
     MapsInfilSwatre->Top = _top;
     MapsInfilSwatre->Left = _left;
     GroupSwatreOptions->Top = _top+12 + MapsInfilSwatre->Height;
     GroupSwatreOptions->Left = _left;

     MapsInfilHoltan->Visible = E_InfilMethod->ItemIndex == INFIL_HOLTAN;
     MapsInfilHoltan->Top = _top;
     MapsInfilHoltan->Left = _left;

     MapsInfilGA1->Visible = E_InfilMethod->ItemIndex == INFIL_GREENAMPT ||
                             E_InfilMethod->ItemIndex == INFIL_GREENAMPT2;
     MapsInfilGA1->Top = _top;
     MapsInfilGA1->Left = _left;

     MapsInfilGA2->Visible = E_InfilMethod->ItemIndex == INFIL_GREENAMPT2;
     MapsInfilGA2->Top = _top+12 + MapsInfilGA1->Height;
     MapsInfilGA2->Left = _left;

     MapsInfilMorel->Visible = E_InfilMethod->ItemIndex == INFIL_MOREL;
     MapsInfilMorel->Top = _top;
     MapsInfilMorel->Left = _left;

     MapsInfilSmith->Visible = E_InfilMethod->ItemIndex == INFIL_SMITH;
     MapsInfilSmith->Top = _top;
     MapsInfilSmith->Left = _left;

     MapsInfilKsat->Visible = E_InfilMethod->ItemIndex == INFIL_KSAT;
     MapsInfilKsat->Top = _top;
     MapsInfilKsat->Left = _left;

     MapsInfilExtra->Top = 280;
     MapsInfilExtra->Left = _left;
     if (!CheckSubsoilDrainage->Checked)
        MapsInfilExtra->Visible = false;
//VJ 050812 drainage with GA
     MapsInfilDrainage->Top = 380;
     MapsInfilDrainage->Left = _left;
     if (E_InfilMethod->ItemIndex == 1)
        MapsInfilDrainage->Visible = false;
     MapsInfilDrainage->Visible = CheckSubsoilDrainage->Checked;

     if (E_InfilMethod->ItemIndex == INFIL_MOREL || E_InfilMethod->ItemIndex == INFIL_SMITH)
     {
        CheckError("Not implemented yet");
     }

     if (E_InfilMethod->ItemIndex == INFIL_SWATRE && E_SwatreDTSEC->Text.IsEmpty())
     {
        Application->MessageBox("Swatre options (min dt) not specified",
        "LISEM Warning", MB_OK+MB_ICONWARNING);
     }
}

//---------------------------------------------------------------------------

void __fastcall TLisIFace::ResetIFace()
{
    ListRunfiles->Items->Clear();
    LabelRunfile->Caption = "Active run file: NONE";

    CheckRunoffPerM->Checked = false;
    CheckChannelInfil->Checked = false;
    CheckChannelBaseflow->Checked = false;
    CheckIncludeChannel->Checked = false;
    CheckNoErosion->Checked = false;
    CheckSnowmelt->Checked = false;

    E_begintime->Text = "";
    E_Endtime->Text = "";
    E_Timestep->Text = "";

    E_InfilMethod->ItemIndex = 0;
    CheckImpermeable->Checked = false;
    CSpinEditKsat->Value = 100;
    E_SwatreDTSEC->Text = "";
    CheckDumphead->Checked = false;
    CheckGeometric->Checked = false;

    E_MapDir->Text = "";
    E_MapDir->Text = "";
    E_TableDir->Text = "";
    RainfallDir = "";

    E_RainfallName->Text = "";
    E_ResultDir->Text = "";
    E_TotalName->Text = "";
    E_OutletName->Text = "";
    E_Outlet1Name->Text = "";
    E_Outlet2Name->Text = "";
    E_ErosionName->Text = "";
    E_DepositionName->Text = "";
//        E_RunoffName->Text = "";

    E_OutputTimeStep->Checked = true;
    if (E_OutputTimeStep->Checked)
       E_OutputTimeSteps->Value = 1;

    E_OutputTimeUser->Checked = false;
    E_OutputTimes->Clear();

    TabWheeltracks->TabVisible=false;

    ResetInfil();

    MapsChannels->Visible = CheckIncludeChannel->Checked;
    MapsChannelinfil->Visible = CheckChannelInfil->Checked;
    Label49->Visible = !CheckIncludeChannel->Checked;
    Label49->Caption = "Activate \"Include main channels\" on the start screen";

    MapsInfilExtra->Visible = false;
    if(CheckInfilCrust->Checked)  MapsInfilExtra->Visible = true;
    if(CheckInfilCompact->Checked)  MapsInfilExtra->Visible = true;
    if(CheckInfilGrass->Checked)  MapsInfilExtra->Visible = true;
    Label55->Visible = CheckBuffers->Checked;
    Label57->Visible = CheckBuffers->Checked;
    MapsBuffers->Visible = CheckBuffers->Checked;
    E_SedBulkDensity->Visible = CheckBuffers->Checked;
    Bevel7->Visible = CheckBuffers->Checked;

}

