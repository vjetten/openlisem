//---------------------------------------------------------------------------
#ifndef lisrainfH
#define lisrainfH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Buttons.hpp>
//---------------------------------------------------------------------------
class TRainForm : public TForm
{
__published:	// IDE-managed Components
        TMemo *RainViewList;
        TBitBtn *BitBtn1;
private:	// User declarations
public:		// User declarations
        __fastcall TRainForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TRainForm *RainForm;
//---------------------------------------------------------------------------
#endif
