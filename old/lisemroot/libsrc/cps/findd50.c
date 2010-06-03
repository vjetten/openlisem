#include "stddefx.h"
#include "cps.h"
#include "misc.h"

#include "chkdata.h"



static void doD50(
	      REAL4   **D50out,
    CONST_REAL4_MAP Si0,
	CONST_REAL4_MAP Si1,
	CONST_REAL4_MAP Si2,
	CONST_REAL4_MAP Si3,
	CONST_REAL4_MAP Si4,
	CONST_REAL4_MAP Si5,
    double *mu)
{
	size_t nrRows   =(int) RgiveNrRows();
	size_t nrCols= RgiveNrCols();
	size_t r, c; /* (r,c) becomes co-ordinate of outflowpoint */

        /* mark all cells as not done: */
	for(r = 0; r < nrRows; r++)
	   SetMemMV(D50out[r],nrCols,CR_REAL4);


	for(r = 0; r < nrRows; r++)
	 for(c = 0; c < nrCols; c++)
       if(!IS_MV_REAL4(&Si0[r][c]))
     {
         double sedtot=0;
         double fr[6];
         int i = 0;
         double f = 0;

         fr[0] = Si0[r][c];
         fr[1] = Si1[r][c];
         fr[2] = Si2[r][c];
         fr[3] = Si3[r][c];
         fr[4] = Si4[r][c];
         fr[5] = Si5[r][c];
         for (i = 0; i < 6; i++)
            sedtot += fr[i];

         if (sedtot > 0)
         {
          for (i = 0; i < 6; i++)
              fr[i] /= sedtot;

           for (i = 0; i < 6; i++)
           {
               f += fr[i];
               if (f > 0.5) break;
           }
           if (i==6)
              i = 5;
           f -= fr[i];
           if (i == 0)
              D50out[r][c] = fr[0];
           else
              D50out[r][c] = mu[i-1]+(0.5-f)/(fr[i])*(mu[i]-mu[i-1]);
         }
         else
           D50out[r][c] = 0;
     }
}


void FindD50(
	const char *srcFile, int srcLineNr,
	MEM_HANDLE *D50out,
   MEM_HANDLE *Si0, MEM_HANDLE *Si1, MEM_HANDLE *Si2,
   MEM_HANDLE *Si3, MEM_HANDLE *Si4, MEM_HANDLE *Si5,
   double *mu)
{
	REAL4 **D50outMap;
	CONST_REAL4_MAP Si0Map, Si1Map, Si2Map, Si3Map, Si4Map, Si5Map;

	SetCaller(srcFile,srcLineNr,"findD50");

	D50outMap=(REAL4 **)CpsNormalHandle(D50out,CR_REAL4);

	Si0Map=(CONST_REAL4_MAP)CpsNormalHandle(Si0,CR_REAL4);
	Si1Map=(CONST_REAL4_MAP)CpsNormalHandle(Si1,CR_REAL4);
	Si2Map=(CONST_REAL4_MAP)CpsNormalHandle(Si2,CR_REAL4);
	Si3Map=(CONST_REAL4_MAP)CpsNormalHandle(Si3,CR_REAL4);
	Si4Map=(CONST_REAL4_MAP)CpsNormalHandle(Si4,CR_REAL4);
	Si5Map=(CONST_REAL4_MAP)CpsNormalHandle(Si5,CR_REAL4);


    doD50(D50outMap, Si0Map, Si1Map, Si2Map, Si3Map, Si4Map, Si5Map, mu);
}



