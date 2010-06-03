//---------------------------------------------------------------------------
#ifndef lisrunfH
#define lisrunfH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
//---------------------------------------------------------------------------
class TRunForm : public TForm
{
__published:	// IDE-managed Components
	TMemo *RunEdit;
	TBitBtn *BitBtn1;
    TBitBtn *BitBtn2;
	void __fastcall FormResize(TObject *Sender);
    void __fastcall BitBtn1Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TRunForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TRunForm *RunForm;
//---------------------------------------------------------------------------
#endif
