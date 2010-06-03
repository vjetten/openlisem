//---------------------------------------------------------------------------
#include <vcl.h>
//#include <algorithm>
#pragma hdrstop

#include "mprolog.h"
#include "lisscreenoutput.h"
#include "lisrunf.h"
#include "lismain.h"
#include "iface.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
static float maxrain;
static float MaxOutputSediment;
static float MaxOutlet2;
static float MaxOutlet3;

#define max(a,b) (a>b?a:b)
//---------------------------------------------------------------------------
__fastcall LisThread::LisThread(bool CreateSuspended)
    : TThread(CreateSuspended)
{
}
//---------------------------------------------------------------------------
void __fastcall LisThread::LisemWarningV()
{
//    Application->MessageBox(WarningText.c_str(),"Lisem warning",MB_OK);
    LisIFace->Messages->Lines->Append(WarningText);
}
//---------------------------------------------------------------------------
void __fastcall LisThread::InitializeTotals()
{
     char s[255];

     sprintf(s,"%.3f",StartTime);
     LisIFace->LO_Tstart->Caption = s;
     sprintf(s,"%.3f",EndTime);
     LisIFace->LO_Tend->Caption = s;
     sprintf(s,"%.3f",CatchmentAreaHa);
     LisIFace->LO_catcharea->Caption = s;
     sprintf(s,"%.3f",CellSize);
     LisIFace->LO_gridcell->Caption = s;
     sprintf(s,"%.3f",timestep*60);
     LisIFace->LO_timestep->Caption = s;

      sprintf(s,"%.3f",0.0);

      LisIFace->LO_dtcurr->Caption = "";
      LisIFace->LO_totrain->Caption = "";
      LisIFace->LO_totdisch->Caption = "";
      LisIFace->LO_raindisch->Caption = "";
      LisIFace->LO_interception->Caption = "";
      LisIFace->LO_infil->Caption = "";
      LisIFace->LO_surfstor->Caption = "";
      LisIFace->LO_dischM3->Caption = "";
      LisIFace->LO_peakdisch->Caption = "";
//VJ 080214 add baseflow
      LisIFace->LO_baseflow->Caption = "";
      LisIFace->LO_peaktime->Caption = "";
      LisIFace->LO_peaktimeP->Caption = "";
      LisIFace->LO_masserror->Caption = "";
      LisIFace->LO_splash->Caption = "";
      LisIFace->LO_flow->Caption = "";
      LisIFace->LO_depo->Caption = "";
      LisIFace->LO_chanflow->Caption = "";
      LisIFace->LO_chandepo->Caption = "";
//VJ 040823 buffers
      LisIFace->LO_buffervol->Caption = "";
      LisIFace->LO_buffersedvol->Caption = "";

      LisIFace->LO_soilloss->Caption = "";
      LisIFace->LO_avgsoilloss->Caption = "";
      LisIFace->LO_outflow->Caption = "";
      LisIFace->LO_Sedbal->Caption = "";

     strcpy(LisIFace->LO_dtcurr->Caption.c_str(), s);
     strcpy(LisIFace->LO_totrain->Caption.c_str(), s);
     strcpy(LisIFace->LO_totdisch->Caption.c_str(),  s);
     strcpy(LisIFace->LO_raindisch->Caption.c_str(),  s);
     strcpy(LisIFace->LO_interception->Caption.c_str(),  s);
     strcpy(LisIFace->LO_infil->Caption.c_str(),  s);
     strcpy(LisIFace->LO_surfstor->Caption.c_str(),  s);

     strcpy(LisIFace->LO_peakdisch->Caption.c_str(),  s);
     strcpy(LisIFace->LO_baseflow->Caption.c_str(),  s);
     strcpy(LisIFace->LO_dischM3->Caption.c_str(),  s);
     strcpy(LisIFace->LO_peaktime->Caption.c_str(),  s);
     strcpy(LisIFace->LO_peaktimeP->Caption.c_str(),  s);
     strcpy(LisIFace->LO_masserror->Caption.c_str(),  s);
     strcpy(LisIFace->LO_splash->Caption.c_str(),  s);
     strcpy(LisIFace->LO_flow->Caption.c_str(),  s);
     strcpy(LisIFace->LO_depo->Caption.c_str(),  s);
     strcpy(LisIFace->LO_SedSuspended->Caption.c_str(),  s);
     strcpy(LisIFace->LO_chanflow->Caption.c_str(),  s);
     strcpy(LisIFace->LO_chandepo->Caption.c_str(),  s);
//VJ 040823 buffers
     strcpy(LisIFace->LO_buffervol->Caption.c_str(),  s);
     strcpy(LisIFace->LO_buffersedvol->Caption.c_str(),  s);
     strcpy(LisIFace->LO_soilloss->Caption.c_str(),  s);
     strcpy(LisIFace->LO_avgsoilloss->Caption.c_str(),  s);
     strcpy(LisIFace->LO_outflow->Caption.c_str(),  s);
     strcpy(LisIFace->LO_Sedbal->Caption.c_str(),  s);

     LisIFace->CGauge->Progress = (CurrentTime-StartTime)/(EndTime-StartTime)*100;

     LisIFace->Qgraph->BottomAxis->Maximum = EndTime;
     LisIFace->Qgraph->BottomAxis->Minimum = StartTime;
     LisIFace->Qgraph->BottomAxis->Increment = 10*int((EndTime-StartTime)/100.0);
     LisIFace->Series1->Clear();
     LisIFace->Series2->Clear();
     LisIFace->Series3->Clear();
     LisIFace->Series4->Clear();
     LisIFace->Series6->Clear();
     LisIFace->Series7->Clear();
     LisIFace->Series1->AddXY(CurrentTime, 0,"",clBlue);
     LisIFace->Series2->AddXY(CurrentTime, 0,"",clRed);
     LisIFace->Series3->AddXY(CurrentTime, 0,"",clLime);
     LisIFace->Series4->AddXY(CurrentTime, 0,"",clBlack);
     LisIFace->Series6->AddXY(CurrentTime, 0,"",clBlack);
     LisIFace->Series7->AddXY(CurrentTime, 0,"",clBlack);
     LisIFace->Qgraph->RightAxis->Maximum = 0;
     LisIFace->Qgraph->LeftAxis->Maximum = 0;

     LisIFace->OutletValues->Lines->Clear();
     LisIFace->FocusControl(LisIFace->OutletValues);
     LisIFace->LO_infilmethod->Caption = ReportText;
     maxrain = 0;
     MaxOutputSediment = 1;
     MaxOutlet2 = 0;
     MaxOutlet3 = 0;
     InitDone = true;
}
//---------------------------------------------------------------------------
void __fastcall LisThread::UpdateTotals()
{

     char s[255];
     sprintf(s, "%.3f",CurrentTime);
     LisIFace->LO_dtcurr->Caption = s;
     sprintf(s, "%.3f",TOTRainMM);
     LisIFace->LO_totrain->Caption = s;
     sprintf(s,"%.3f",TOTInterceptionMM);
     LisIFace->LO_interception->Caption = s;
     sprintf(s,"%.3f",TOTInfilMM);
     LisIFace->LO_infil->Caption = s;
     sprintf(s,"%.3f",TOTSurfStorMM);
     LisIFace->LO_surfstor->Caption = s;
     sprintf(s,"%.3f",TOTRunoffMM);
     LisIFace->LO_Runoff->Caption = s;
     sprintf(s,"%.3f",TOTDischargeMM);
     LisIFace->LO_totdisch->Caption = s;
     sprintf(s,"%.2f",P_QPercentage);
     LisIFace->LO_raindisch->Caption = s;
     sprintf(s,"%.1f",TOTDischargeM3);
     LisIFace->LO_dischM3->Caption = s;
     sprintf(s,"%.2f",PeakDischargeQ);
     LisIFace->LO_peakdisch->Caption = s;
//VJ 080217 add baseflow
//     sprintf(s,"%.2f",BaseflowDischargeQ);
//     LisIFace->LO_baseflow->Caption = s;

     sprintf(s,"%.2f",PeakDischargeTime);
     LisIFace->LO_peaktime->Caption = s;
     sprintf(s,"%.2f",PeakRainfallTime);
     LisIFace->LO_peaktimeP->Caption = s;
     sprintf(s,"%.5f",MassBalance);
     LisIFace->LO_masserror->Caption = s;
     sprintf(s,"%.3f",TOTSplashErosion);
     LisIFace->LO_splash->Caption = s;
     sprintf(s,"%.3f",TOTFlowErosion);
     LisIFace->LO_flow->Caption = s;
     sprintf(s,"%.3f",TOTDeposition);
     LisIFace->LO_depo->Caption = s;
     sprintf(s,"%.3f",TOTChanFlowErosion);
     LisIFace->LO_chanflow->Caption = s;
     sprintf(s,"%.3f",TOTChanDeposition);
     LisIFace->LO_chandepo->Caption = s;
//VJ 040823 include buffers
     sprintf(s,"%.2f",TOTBufferVolume);
     LisIFace->LO_buffervol->Caption = s;
     sprintf(s,"%.2f",TOTBufferSedVolume);
     LisIFace->LO_buffersedvol->Caption = s;
     sprintf(s,"%.3f",TOTSoilLoss);
     LisIFace->LO_soilloss->Caption = s;
     sprintf(s,"%.1f",AVGSoilLoss);
     LisIFace->LO_avgsoilloss->Caption = s;
     sprintf(s,"%.2f",OutputDischarge);
     LisIFace->LO_outflow->Caption = s;
//     LisIFace->LO_sedoutflow->Caption = Convert(OutputSedConc,2);
     sprintf(s,"%.5f",SedMassBalance);
     LisIFace->LO_Sedbal->Caption = s;
     sprintf(s,"%.3f",OutputSedSusp);
     LisIFace->LO_SedSuspended->Caption = s;
     sprintf(s,"%.3f",OutputChanSedSusp);
     LisIFace->LO_ChSedSuspended->Caption = s;

     LisIFace->CGauge->Progress = (CurrentTime-StartTime)/(EndTime-StartTime)*100;

     maxrain = max(AVGRainMM, maxrain);

     MaxOutputSediment = max(MaxOutputSediment, OutputSediment);
     MaxOutlet2 = max(MaxOutlet2, OutputDischarge1);
     MaxOutlet3 = max(MaxOutlet3, OutputDischarge1);

     //add factor 1.05 to show top axis above top graph value
      if (LisIFace->ShowOutlet1->Checked)
           LisIFace->Qgraph->LeftAxis->Maximum = max(PeakDischargeQ,maxrain) * 1.05;
      else
      if (LisIFace->ShowOutlet2->Checked)
           LisIFace->Qgraph->LeftAxis->Maximum = max(MaxOutlet2,maxrain) * 1.05;
      else
      if (LisIFace->ShowOutlet3->Checked)
           LisIFace->Qgraph->LeftAxis->Maximum = max(MaxOutlet3,maxrain) * 1.05;

      if (LisIFace->ShowRainfall->Checked)
           LisIFace->Series4->AddXY(CurrentTime, AVGRainMM,"",clBlack);
      if (LisIFace->ShowOutlet1->Checked)
           LisIFace->Series1->AddXY(CurrentTime, OutputDischarge,"",clNavy);

     if (!LisIFace->CheckNoErosion->Checked)
     {
        if (LisIFace->ShowQsed->Checked)
        {
             LisIFace->Series2->AddXY(CurrentTime, OutputSediment,"",clMaroon);
             LisIFace->Qgraph->RightAxis->Maximum = MaxOutputSediment * 1.05;
        }
        if (LisIFace->ShowSedConcentration->Checked)
        {
             LisIFace->Series3->AddXY(CurrentTime, OutputSedConc,"",clRed);
             LisIFace->Qgraph->RightAxis->Maximum = max(LisIFace->Qgraph->RightAxis->Maximum,
                         OutputSedConc * 1.05);
        }
     }

      if (LisIFace->ShowOutlet2->Checked)
           LisIFace->Series6->AddXY(CurrentTime, OutputDischarge1,"",clBlue);
      if (LisIFace->ShowOutlet3->Checked)
           LisIFace->Series7->AddXY(CurrentTime, OutputDischarge2,"",clFuchsia);

//VJ 030624 added outlet values lijnes displayed, saves memory, -1 is display all


     if (LisIFace->SpinOutletValues->Value > -1)
	     if (LisIFace->OutletValues->Lines->Count+1 > LisIFace->SpinOutletValues->Value)
     		     LisIFace->OutletValues->Lines->Delete(0);

     if (LisIFace->SpinOutletValues->Value != 0)
     {
   	  sprintf(s,"    %5d %8.3f\t%6.2f\t%12.3f\t%8.3f",StepCounter, CurrentTime,AVGRainMM,OutputDischarge,OutputSedConc);
	     LisIFace->OutletValues->Lines->Add((AnsiString)s);
     }
}
//---------------------------------------------------------------------------
AnsiString LisThread::Convert(double f, int dig)
{
    char s[255];
    char fmt[5];
    sprintf(fmt,"%%.%df",dig);
    sprintf(s,fmt,f);
    return ((AnsiString)s);
    // beetje lomp
}
//---------------------------------------------------------------------------
void __fastcall LisThread::DumpScreenSync()
{
// VJ 031218 check dumpscreen to halt processes to finish dumpscreen
   LisIFace->CheckDumpScreen = LisIFace->DumpScreen(false);
}
//---------------------------------------------------------------------------
/*
void __fastcall LisThread::GetNewRunfile()
{
   LisIFace->RunFilename = LisIFace->ListRunfiles->Items->Strings[LisIFace->thisrun];
   LisIFace->SaveDialog->FileName = LisIFace->RunFilename;
   RunForm->RunEdit->Lines->LoadFromFile(LisIFace->RunFilename);
   LisIFace->ReadNewRunfile(LisIFace->RunFilename);
   LisIFace->PageControl->ActivePage = LisIFace->TabSheetTot;
   LisIFace->LabelRunfile->Caption = "Active run file: " + LisIFace->RunFilename;
}
//---------------------------------------------------------------------------
*/


