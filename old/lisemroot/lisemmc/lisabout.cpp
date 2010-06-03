//---------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "lisabout.h"
//--------------------------------------------------------------------- 
#pragma link "JvExControls"
#pragma link "JvLinkLabel"
#pragma link "LMDCustomComponent"
#pragma link "LMDStarter"
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





void __fastcall TAboutBox::Label4Click(TObject *Sender)
{
    LMDStarter1->Command = "http://www.itc.nl/lisem";
    LMDStarter1->Execute();
}
//---------------------------------------------------------------------------

