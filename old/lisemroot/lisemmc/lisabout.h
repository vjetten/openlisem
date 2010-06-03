//----------------------------------------------------------------------------
#ifndef lisaboutH
#define lisaboutH
//----------------------------------------------------------------------------
#include <vcl\System.hpp>
#include <vcl\Windows.hpp>
#include <vcl\SysUtils.hpp>
#include <vcl\Classes.hpp>
#include <vcl\Graphics.hpp>
#include <vcl\Forms.hpp>
#include <vcl\Controls.hpp>
#include <vcl\StdCtrls.hpp>
#include <vcl\Buttons.hpp>
#include <vcl\ExtCtrls.hpp>
#include "JvExControls.hpp"
#include "JvLinkLabel.hpp"
#include "LMDCustomComponent.hpp"
#include "LMDStarter.hpp"
//----------------------------------------------------------------------------
class TAboutBox : public TForm
{
__published:
	TLabel *ProductName;
	TLabel *Version;
	TLabel *ProductDate;
	TLabel *Copyright1;
	TLabel *Copyright2;
	TLabel *Authors1;
	TLabel *Contrib;
	TLabel *Authors2;
	TPanel *Panel2;
	TImage *Image1;
	TBevel *Bevel1;
        TLabel *Label3;
    TLabel *Label4;
    TBevel *Bevel3;
        TButton *Button1;
   TLabel *Label1;
   TLabel *Label2;
   TLMDStarter *LMDStarter1;
        void __fastcall Button1Click(TObject *Sender);
   void __fastcall Label4Click(TObject *Sender);
private:
public:
	virtual __fastcall TAboutBox(TComponent* AOwner);
};
//----------------------------------------------------------------------------
extern PACKAGE TAboutBox *AboutBox;
//----------------------------------------------------------------------------
#endif    
