//----------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "lisdirv.h"
//----------------------------------------------------------------------------
#pragma resource "*.dfm"
TDirView *DirView;
//----------------------------------------------------------------------------
__fastcall TDirView::TDirView(TComponent *Owner)
	: TForm(Owner)
{
}
//----------------------------------------------------------------------------
void __fastcall TDirView::Button2Click(TObject *Sender)
{
    DirView->DirectoryListBox->Directory =backupdir;
}
//---------------------------------------------------------------------------


void __fastcall TDirView::FormActivate(TObject *Sender)
{
     backupdir = DirView->DirectoryListBox->Directory;
}
//---------------------------------------------------------------------------


void __fastcall TDirView::FormCreate(TObject *Sender)
{
     DriveComboBox1->DirList = DirectoryListBox;        
}
//---------------------------------------------------------------------------

