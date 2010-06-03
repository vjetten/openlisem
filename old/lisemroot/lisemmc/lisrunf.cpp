//---------------------------------------------------------------------------
#include <vcl.h>
#include <algorithm.h>
#pragma hdrstop

#include "lisrunf.h"
#include "iface.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TRunForm *RunForm;
//---------------------------------------------------------------------------
__fastcall TRunForm::TRunForm(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TRunForm::FormResize(TObject *Sender)
{
   RunForm->Height = max(RunForm->Height, 120);
   RunEdit->Height = RunForm->Height-74;
   RunEdit->Width = RunForm->Width-16;
   BitBtn1->Top = RunEdit->Height + 12;
   BitBtn1->Left = RunEdit->Left+RunEdit->Width-BitBtn1->Width-BitBtn2->Width-20;
   BitBtn2->Top = RunEdit->Height + 12;
   BitBtn2->Left = RunEdit->Left+RunEdit->Width-BitBtn2->Width-16;
}
//---------------------------------------------------------------------------

void __fastcall TRunForm::BitBtn1Click(TObject *Sender)
{

     if (RunEdit->Modified)
     {
       if(Application->MessageBox("Pressing OK will save changes and overwrite existing run file!",
        "LISEM Error", MB_OKCANCEL+MB_ICONWARNING) == IDOK)
        {
            RunForm->RunEdit->Lines->SaveToFile(LisIFace->RunFilename);
            RunForm->RunEdit->Modified = false;
//          LisIFace->MakeNewRunfile(LisIFace->RunFilename);
          LisIFace->ReadNewRunfile(LisIFace->RunFilename);
        }  
     }
}
//---------------------------------------------------------------------------

