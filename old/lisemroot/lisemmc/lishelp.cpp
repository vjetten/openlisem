//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

//#include "lishelp.h"
#include "iface.h"


//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
THelpForm *HelpForm;
//---------------------------------------------------------------------------
__fastcall THelpForm::THelpForm(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall THelpForm::BitBtn1Click(TObject *Sender)
{
        Close();
}
//---------------------------------------------------------------------------



void __fastcall THelpForm::FormCreate(TObject *Sender)
{
      Memo1->Lines->Strings[0] = "LISEM for Windows - "+(AnsiString)ProgVersion;

}
//---------------------------------------------------------------------------





