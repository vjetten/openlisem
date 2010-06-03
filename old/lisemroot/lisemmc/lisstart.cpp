//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop
#include <process.h>

#include "lisstart.h"
#include "iface.h"


//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TStartForm *StartForm;
//---------------------------------------------------------------------------
__fastcall TStartForm::TStartForm(TComponent* Owner)
        : TForm(Owner)
{
    Caption = (AnsiString)ProgName + (AnsiString)" - " + (AnsiString)ProgVersion;
}
//---------------------------------------------------------------------------
void Visibility()
{
    LisIFace->Enabled = true;
    StartForm->Close();
    LisIFace->CheckLisemType();
}
//---------------------------------------------------------------------------


void __fastcall TStartForm::BitBtn1Click(TObject *Sender)
{
     LisIFace->LisemType = LISEMBASIC;
     Visibility();
}
//---------------------------------------------------------------------------

void __fastcall TStartForm::BitBtn2Click(TObject *Sender)
{
     LisIFace->LisemType = LISEMWHEELTRACKS;
     Visibility();
}
//---------------------------------------------------------------------------

void __fastcall TStartForm::BitBtn3Click(TObject *Sender)
{
     LisIFace->LisemType = LISEMMULTICLASS;
     Visibility();
}
//---------------------------------------------------------------------------

void __fastcall TStartForm::BitBtn4Click(TObject *Sender)
{
     LisIFace->LisemType = LISEMNUTRIENTS;
     Visibility();
}
//---------------------------------------------------------------------------

void __fastcall TStartForm::BitBtn5Click(TObject *Sender)
{
     LisIFace->LisemType = LISEMGULLIES;
     Visibility();
}
//---------------------------------------------------------------------------

void __fastcall TStartForm::Image2Click(TObject *Sender)
{
   int hoi = spawnl(P_NOWAIT, "c:\\program files\\internet explorer\\iexplore.exe",
   "c:\\program files\\internet explorer\\iexplore.exe","http://www.itc.nl/lisem",NULL);
   if (hoi == -1)
    spawnl(P_NOWAIT, "c:\\program files\\netscape\\netscape.exe",
   "c:\\program files\\netscape\\netscape.exe","http://www.itc.nl/lisem",NULL);
}
//---------------------------------------------------------------------------

void __fastcall TStartForm::BitBtn6Click(TObject *Sender)
{
   StartForm->Close();
   LisIFace->Enabled;
   LisIFace->Close();
}
//---------------------------------------------------------------------------




