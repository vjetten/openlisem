//---------------------------------------------------------------------------
// Thread relate functions, run, stop, pause

#pragma hdrstop
#include <dir.h>
#include "ifacethread.h"
#include "iface.h"

//---------------------------------------------------------------------------

#pragma package(smart_init)

//---------------------------------------------------------------------------
void __fastcall TLisIFace::ButtonStopprogClick(TObject *Sender)
{
    Beep();
    if (LisemRun)
    {
       LisemRun->RunDone = true;
       if (LisemRun->Suspended)
          LisemRun->Resume();
    }      

//    done = true;
    LisIFace->Messages->Lines->Append("User Interrupt, please wait ...");
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::ButtonPauseprogClick(TObject *Sender)
{
    if (!LisemRun)
       return;
    if (!take5)
    {
       LisemRun->Suspend();
       take5 = true;
       LisIFace->Messages->Lines->Append("Simulation paused, press run to continue...");
    }
}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::ButtonRunprogClick(TObject *Sender)
{
    if (E_OutletName->Text.IsEmpty())
      if (Application->MessageBox("No \'Point Output\' filename, please give a name", "LISEM Warning", MB_OK+MB_ICONERROR) == IDOK)
          return;
    if (E_SoillossName->Text.IsEmpty())
      if (Application->MessageBox("No \'Soil Loss map\' filename, please give a name", "LISEM Warning", MB_OK+MB_ICONERROR) == IDOK)
          return;
    if (take5)
    {
       LisemRun->Resume();
       take5 = false;
    }
    else
    if (!LisemRun)
    {
        // clear message windows
       Messages->Clear();

       PageControl->ActivePage = TabSheetTot;
       ThreadDone(Sender);
    }
    else
      LisIFace->Messages->Lines->Append("Run active, press stop first ...");

}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::ThreadDone(TObject *Sender)
{
    ButtonRunProg->Down = false;
    ButtonStopProg->Down = true;
    take5 = false;

    KillDrawMapStructure();


    if (thisrun < nrruns)
    {
       if (multipleRuns)
          GetNewRunfile();    //load first or next runfile in list if multiple runs

       ButtonRunProg->Down = true;
       ButtonStopProg->Down = false;

       t_begin = Time();
       runduration = TimeToStr(t_begin);
       // fill a Tstringlist with in and output map names


       LisemRun = new LisThread(true);    // mke and run the thread, the old LISEM model!
       LisemRun->FreeOnTerminate = true;
       LisemRun->RunDone = false;
       LisemRun->OnTerminate = ThreadDone;
       LisIFace->LabelRunfile->Caption = "Active run file: " + LisIFace->RunFilename;
       LisIFace->Messages->Clear();
       LisIFace->Messages->Lines->Append("Preparing input data ...");
       strcpy(LisemRun->temprunname, ExtractFilePath(ParamStr(0)).c_str());
       strcat(LisemRun->temprunname,"lisemtemp.run");
       MakeNewRunfile(LisemRun->temprunname);
       LisemRun->Resume();
       // dump the current configuration to a temp run file to be read by the model itself

       thisrun++;

    }
    else
    {
       if (multipleRuns && !batchrun)
       {
          Application->MessageBox("Finished run files","LISEM message",MB_OK+MB_ICONEXCLAMATION);
          thisrun = 0;
         // reset the nr of runs
       }
       else
          thisrun = 0;

       LisemRun = NULL;

       if (closeapp || (batchrun && thisrun == 0))
          Close();
    }
}
//---------------------------------------------------------------------------
/*
void __fastcall TLisIFace::MakeNamelist()
{
// makes a tstringlist with all map names from the interface for easy reference
   TStringList *nl = new TStringList();

   for(int i = 0; i < ComponentCount; i++)
    if (Components[i]->Name.SubString(0,4) == "Maps")
     {
         TStringGrid* S =(TStringGrid *)LisIFace->Components[i];
         for (int j = 1; j < S->RowCount; j++)
         {
             if (S->Cells[3][j].IsPathDelimiter(S->Cells[3][j].Length()))
                 S->Cells[4][j] = S->Cells[3][j]+S->Cells[1][j];
              else
                 S->Cells[4][j] = S->Cells[3][j]+"\\"+S->Cells[1][j];
             nl->Add(S->Cells[5][j]+"="+S->Cells[4][j]);
         }
    }
    if (LisemRun->namelist)
       KillNamelist();
    LisemRun->nrnamelist = nl->Count;
    LisemRun->namelist = (nameList *)malloc(LisemRun->nrnamelist*sizeof(nameList));

    for (int i = 0; i < LisemRun->nrnamelist; i++)
    {
       strcpy(LisemRun->namelist[i].ID, nl->Names[i].UpperCase().c_str());
       strcpy(LisemRun->namelist[i].fullname, nl->Values[nl->Names[i]].c_str());
    }


    delete nl;


}
//---------------------------------------------------------------------------
void __fastcall TLisIFace::KillNamelist()
{
    free(LisemRun->namelist);
    LisemRun->namelist = NULL;
}
*/
