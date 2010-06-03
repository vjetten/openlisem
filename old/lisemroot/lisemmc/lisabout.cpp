//---------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "lisabout.h"
//--------------------------------------------------------------------- 
#pragma link "JvExControls"
#pragma link "JvLinkLabel"
#pragma resource "*.dfm"
TAboutBox *AboutBox;
//--------------------------------------------------------------------- 
__fastcall TAboutBox::TAboutBox(TComponent* AOwner)
	: TForm(AOwner)
{
}
//---------------------------------------------------------------------




void __fastcall TAboutBox::Button1Click(TObject *Sender)
{
    Close();        
}
//---------------------------------------------------------------------------






