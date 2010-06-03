//---------------------------------------------------------------------------
#include <vcl.h>
#include <stdlib.h>
#include <math.h>
//---------------------------------------------------------------------------
#pragma hdrstop

#include "CsfMap.h"
#define ErrorMsg Application->MessageBox((char *)MstrError(),"MapEdit error",MB_OK)

//---------------------------------------------------------------------------
__fastcall TMap::TMap()
{
    Data = NULL;
    Created = false;
    nrRows = 0;
    nrCols = 0;    
}
//---------------------------------------------------------------------------
__fastcall TMap::~TMap()
{
    KillMap();
}
//---------------------------------------------------------------------------
void TMap::KillMap()
{
    if (Data)
    {
        for(int r=0; r < nrRows; r++)
            delete[] Data[r];
        delete[] Data;
      Data = NULL;
    }
    Created = false;
}
//---------------------------------------------------------------------------
void TMap::LoadFromFile(MEM_HANDLE *M, bool start)
{
    CONST_REAL4_MAP MMap;
  	 MMap=(CONST_REAL4_MAP)CpsNormalHandle(M,CR_REAL4);
    nrRows = RgiveNrRows();
	 nrCols = RgiveNrCols();

    if (start)
    {
       KillMap();

       Data = new REAL4*[nrRows];
       for(int r=0; r < nrRows; r++)
          Data[r] = new REAL4[nrCols];

       if (Data == NULL)
          return;
       Created = true;
     }

    if (!Created)
       return;
/*
    MapName = Name;

    m = Mopen(MapName.c_str(), M_READ);
    RuseAs(m, CR_REAL4);
    for(int r=0; r < nrRows; r++)
       RgetSomeCells(m, (UINT4)r*nrCols, (UINT4)nrCols, Data[r]);

    Mclose(m);
*/

   MRH.valueScale = VS_SCALAR;

	for(int r = 0; r < nrRows; r++)
	   SetMemMV(Data[r],nrCols,CR_REAL4);

    for(int r=0; r < nrRows; r++)
      for(int c=0; c < nrCols; c++)
        if (!IS_MV_REAL4(&MMap[r][c]))
         Data[r][c] = MMap[r][c];


    ResetMinMax();
}
//---------------------------------------------------------------------------
void TMap::ResetMinMax(void)
{
     REAL8 minv = 1e30, maxv = -1e30;

   //  if (!Created)
     //   return;

      for(int r=0; r < nrRows; r++)
        for(int c=0; c < nrCols; c++)
        if (!IS_MV_REAL4(&Data[r][c]))
        {
           if (maxv <Data[r][c]) maxv = Data[r][c];
           if (minv >Data[r][c]) minv = Data[r][c];
        }

     MRH.minVal = minv;
     MRH.maxVal = maxv;
}


