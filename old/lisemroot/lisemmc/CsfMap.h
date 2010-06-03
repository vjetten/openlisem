//---------------------------------------------------------------------------
#ifndef CsfMapH
#define CsfMapH
//---------------------------------------------------------------------------
#include <SysUtils.hpp>
#include <Controls.hpp>
#include <Classes.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>

#include "csf.h"
#include "cps.h"

#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))
//---------------------------------------------------------------------------

class TMap
{
protected:
    char error[64];
public:
    CSF_RASTER_HEADER MRH;
    CSF_MAIN_HEADER MMH;
    
    void **CreateMatrix(CSF_CR cr);
    REAL4 **Data;
    AnsiString MapName;
    int nrRows, nrCols;

    CSF_VS GetValueScale();
    CSF_CR GetCellRepr();
    float GetCellSize();

    bool Created;
    void KillMap();
    
    void LoadFromFile(MEM_HANDLE *M, bool start);
    void ResetMinMax(void);
    __fastcall TMap();
    __fastcall ~TMap();
};

#endif
