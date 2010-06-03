#include "stddefx.h" 


/********/
/* USES */
/********/

/* libs ext. <>, our ""  */
#include <math.h> /* sqrt */
#include "cps.h"
#include "misc.h"

#include "chkdata.h"


/* global header (opt.) and test's prototypes "" */


/* headers of this app. modules called */ 

/***************/
/* EXTERNALS   */
/***************/

/**********************/ 
/* LOCAL DECLARATIONS */
/**********************/ 

/*********************/ 
/* LOCAL DEFINITIONS */
/*********************/ 

/******************/
/* IMPLEMENTATION */
/******************/
void Edge(const char *srcFile,
	  int    srcLineNr,
          MEM_HANDLE *OutMap,
          MEM_HANDLE *InMap,
          REAL4 value)
{
	REAL4 **out;
	CONST_REAL4_MAP In;
	int 	r, c;
 	int nrCols = (int)RgiveNrCols();
 	int nrRows = (int)RgiveNrRows();

	SetCaller(srcFile,srcLineNr,"edge");

        out=(REAL4 **)CpsNormalHandle(OutMap,CR_REAL4);
        In=(CONST_REAL4_MAP)CpsNormalHandle(InMap,CR_REAL4);

	/* Fill out with missing values. This is the initial value. */
	for(r=0; r < nrRows; r++)
           SetMemMV(out[r],(size_t)nrCols,CR_REAL4);

	/* for every cell check where it flows to */
	for(r = 0; r < nrRows; r++)
	{
           for(c = 0; c < nrCols; c++)
	   {
	      if (!IS_MV_REAL4(In[r]+c))
              {
                  out[r][c]=In[r][c];
                  if (c == 0 || c == nrCols-1 || r == 0 || r == nrRows-1)
                     out[r][c] = value;
                  if (c > 0 && r > 0 && c < nrCols-1 && r < nrRows -1)
                  {
                    if (
                        IS_MV_REAL4(In[r-1]+c-1)||
                        IS_MV_REAL4(In[r-1]+c)||
                        IS_MV_REAL4(In[r-1]+c+1)||
                        IS_MV_REAL4(In[r]+c-1)||
                        IS_MV_REAL4(In[r]+c+1)||
                        IS_MV_REAL4(In[r+1]+c-1)||
                        IS_MV_REAL4(In[r+1]+c)||
                        IS_MV_REAL4(In[r+1]+c+1)
                       )
                     out[r][c] = value;
                  }
              }    
           }
	}
	return;
}


