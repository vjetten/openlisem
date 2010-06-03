//----------------------------------------------------------------------------
#ifndef lisdirvH
#define lisdirvH
//----------------------------------------------------------------------------
#include <SysUtils.hpp>
#include <Windows.hpp>
#include <Messages.hpp>
#include <Classes.hpp>
#include <Graphics.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <ExtCtrls.hpp>
#include <Forms.hpp>
#include <FileCtrl.hpp>
//----------------------------------------------------------------------------
class TDirView : public TForm
{
__published:
	TButton *Button1;
	TButton *Button2;
        TDirectoryListBox *DirectoryListBox;
        TDriveComboBox *DriveComboBox1;
	void __fastcall Button2Click(TObject *Sender);
	void __fastcall FormActivate(TObject *Sender);
        void __fastcall FormCreate(TObject *Sender);
private:
public:
	virtual __fastcall TDirView(TComponent *Owner);
        AnsiString backupdir;
};
//----------------------------------------------------------------------------
extern TDirView *DirView;
//----------------------------------------------------------------------------
#endif
