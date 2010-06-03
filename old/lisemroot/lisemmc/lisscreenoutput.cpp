//---------------------------------------------------------------------------
#include <vcl.h>
//#include <algorithm>
#pragma hdrstop

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
void __fastcall LisThread::InitializeTotalsSync()
{

     LisIFace->LO_Tstart->Caption = Format("%.3f",&TVarRec(StartTime), 0);
     LisIFace->LO_Tend->Caption = Format("%.3f",&TVarRec(EndTime), 0);
     LisIFace->LO_catcharea->Caption = Format("%.3f",&TVarRec(CatchmentAreaHa), 0);
     LisIFace->LO_gridcell->Caption = Format("%.3f",&TVarRec(CellSize), 0);
     LisIFace->LO_timestep->Caption = Format("%.3f",&TVarRec(timestep*60), 0);

      LisIFace->LO_dtcurr->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_totrain->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_totdisch->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_raindisch->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_interception->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_infil->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_surfstor->Caption = Format("%.3f",&TVarRec(0.0), 0);

      LisIFace->LO_peakdisch->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_baseflow->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_dischM3->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_peaktime->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_peaktimeP->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_masserror->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_splash->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_flow->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_depo->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_SedSuspended->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_chanflow->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_chandepo->Caption = Format("%.3f",&TVarRec(0.0), 0);
      //VJ 040823 buffers
      LisIFace->LO_buffervol->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_buffersedvol->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_soilloss->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_avgsoilloss->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_outflow->Caption = Format("%.3f",&TVarRec(0.0), 0);
      LisIFace->LO_Sedbal->Caption = Format("%.3f",&TVarRec(0.0), 0);

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
void __fastcall LisThread::UpdateTotalsSync()
{
     LisIFace->LO_dtcurr->Caption = Format("%.3f",&TVarRec(CurrentTime), 0);
     LisIFace->LO_totrain->Caption = Format("%.3f",&TVarRec(TOTRainMM), 0);
     LisIFace->LO_interception->Caption = Format("%.3f",&TVarRec(TOTInterceptionMM), 0);
     LisIFace->LO_infil->Caption = Format("%.3f",&TVarRec(TOTInfilMM), 0);
     LisIFace->LO_surfstor->Caption = Format("%.3f",&TVarRec(TOTSurfStorMM), 0);
     LisIFace->LO_Runoff->Caption = Format("%.3f",&TVarRec(TOTRunoffMM), 0);
     LisIFace->LO_totdisch->Caption = Format("%.3f",&TVarRec(TOTDischargeMM), 0);
     LisIFace->LO_raindisch->Caption = Format("%.2f",&TVarRec(P_QPercentage), 0);
     LisIFace->LO_dischM3->Caption = Format("%.1f",&TVarRec(TOTDischargeM3), 0);
     LisIFace->LO_peakdisch->Caption = Format("%.2f",&TVarRec(PeakDischargeQ), 0);
      //VJ 080217 add baseflow
      //
      //     LisIFace->LO_baseflow->Caption = Format("%.2f",&TVarRec(BaseflowDischargeQ), 0);
      LisIFace->LO_peaktime->Caption = Format("%.2f",&TVarRec(PeakDischargeTime), 0);
      LisIFace->LO_peaktimeP->Caption = Format("%.2f",&TVarRec(PeakRainfallTime), 0);
      LisIFace->LO_masserror->Caption = Format("%.5f",&TVarRec(MassBalance), 0);
      LisIFace->LO_splash->Caption = Format("%.3f",&TVarRec(TOTSplashErosion), 0);
      LisIFace->LO_flow->Caption = Format("%.3f",&TVarRec(TOTFlowErosion), 0);
      LisIFace->LO_depo->Caption = Format("%.3f",&TVarRec(TOTDeposition), 0);
      LisIFace->LO_chanflow->Caption = Format("%.3f",&TVarRec(TOTChanFlowErosion), 0);
      LisIFace->LO_chandepo->Caption = Format("%.3f",&TVarRec(TOTChanDeposition), 0);
      //vj 040823 include buffers
      LisIFace->LO_buffervol->Caption = Format("%.2f",&TVarRec(TOTBufferVolume), 0);
      LisIFace->LO_buffersedvol->Caption = Format("%.2f",&TVarRec(TOTBufferSedVolume), 0);
      LisIFace->LO_soilloss->Caption = Format("%.3f",&TVarRec(TOTSoilLoss), 0);
      LisIFace->LO_avgsoilloss->Caption = Format("%.1f",&TVarRec(AVGSoilLoss), 0);
      LisIFace->LO_outflow->Caption = Format("%.2f",&TVarRec(OutputDischarge), 0);
  //    LisIFace->LO_sedoutflow->Caption = Convert(OutputSedConc,2);

      LisIFace->LO_Sedbal->Caption = Format("%.5f",&TVarRec(SedMassBalance), 0);
      LisIFace->LO_SedSuspended->Caption = Format("%.3f",&TVarRec(OutputSedSusp), 0);
      LisIFace->LO_ChSedSuspended->Caption = Format("%.3f",&TVarRec(OutputChanSedSusp), 0);

      LisIFace->CGauge->Progress = (CurrentTime-StartTime)/(EndTime-StartTime)*100;

     maxrain = max(AVGRainMM, maxrain);

     MaxOutputSediment = max(MaxOutputSediment, OutputSediment);
     MaxOutlet2 = max(MaxOutlet2, OutputDischarge1);
     MaxOutlet3 = max(MaxOutlet3, OutputDischarge1);

     //add factor 1.05 to show top axis above top graph value
      if (LisIFace->ShowOutlet1->Checked)
           LisIFace->Qgraph->LeftAxis->Maximum = max(PeakDischargeQ,maxrain) * 1.05;
 /*     else
      if (LisIFace->ShowOutlet2->Checked)
           LisIFace->Qgraph->LeftAxis->Maximum = max(MaxOutlet2,maxrain) * 1.05;
      else
      if (LisIFace->ShowOutlet3->Checked)
           LisIFace->Qgraph->LeftAxis->Maximum = max(MaxOutlet3,maxrain) * 1.05;
 */
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
/*
      if (LisIFace->ShowOutlet2->Checked)
           LisIFace->Series6->AddXY(CurrentTime, OutputDischarge1,"",clBlue);
      if (LisIFace->ShowOutlet3->Checked)
           LisIFace->Series7->AddXY(CurrentTime, OutputDischarge2,"",clFuchsia);
*/
//VJ 030624 added outlet values lijnes displayed, saves memory, -1 is display all


     if (LisIFace->SpinOutletValues->Value > -1)
	     if (LisIFace->OutletValues->Lines->Count+1 > LisIFace->SpinOutletValues->Value)
     		     LisIFace->OutletValues->Lines->Delete(0);

     if (LisIFace->SpinOutletValues->Value != 0)
     {
        TVarRec args[6] = {StepCounter, CurrentTime,AVGRainMM,OutputDischarge,OutputSediment,OutputSedConc};
        AnsiString s = Format("%5d %8.3f\t%6.2f\t%10.3f\t%10.3f\t%8.3f",args, 5);
	     LisIFace->OutletValues->Lines->Add(s);
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
void __fastcall LisThread::DrawMapSync()//MEM_HANDLE *M, double timestep)
{
   double timestep = 0;
   LisIFace->StartDrawMap(Mout, timestep);
}
//---------------------------------------------------------------------------
void __fastcall LisThread::DrawMapInitSync()//MEM_HANDLE *M, double timestep)
{
   LisIFace->StartDrawMap(Mout, -1);
}
//---------------------------------------------------------------------------

