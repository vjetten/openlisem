#include "stddefx.h" 


/********/
/* USES */
/********/

/* libs ext. <>, our ""  */
#include <math.h> /* sqrt */
#include "cps.h"
#include "misc.h"

#include "lddutil.h"
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

/* these are copies taken from pcrcalc.lib
 */
/*
 #define FLOWS_TO(ldd, rFrom, cFrom, rTo, cTo) \
		( ldd[rFrom][cFrom]!=MV_UINT1 &&\
		  rFrom+LddData[ldd[rFrom][cFrom]].deltaY==rTo &&\
		  cFrom+LddData[ldd[rFrom][cFrom]].deltaX==cTo )
*/

static void DoUpstream(
           REAL4 **out,			/* write-only output map  */
     CONST_REAL4_MAP in,		/* value map */
     CONST_UINT1_MAP ldd) 		/* ldd map		*/
{
//	UINT1 	l;
   	int 	r, c;
 	int nrCols = (int)RgiveNrCols();
 	int nrRows = (int)RgiveNrRows();

	for(r = 0; r < nrRows; r++)
	{
	 for(c = 0; c < nrCols; c++)
	{
	 if(ldd[r][c] != MV_UINT1)
	   {
	    REAL8 sum = 0;
	    int d;
	   CheckReal4Cell(in,r,c, "arg.nr 2:amount-in");
	    for(d=1; d <= 9; d++)
	    {
          int rNext = ((int)r)+LddData[d].deltaY;
          int cNext = ((int)c)+LddData[d].deltaX;
          if (d == 5)
              continue;
          if (rNext>=0 && rNext<nrRows &&
		      cNext>=0 && cNext<nrCols &&
		      FLOWS_TO(ldd,rNext,cNext,r, c))
		    {
		     if (!IS_MV_REAL4(in[rNext]+cNext))
		     	sum += in[rNext][cNext];
		     else
		     {  /* mv in in-map make it's downstream MV
		         */
		         goto putMV;
		     }
		    }
	    } /* eofor all directions */
	    out[r][c] = sum;
	  }
	 else
putMV:
	  SET_MV_REAL4(out[r]+c);
	} /* eofor cols */
       } /* eofor rows */
}

void Upstream(
		const char *srcFile,
		int    srcLineNr,
               MEM_HANDLE *amountOutMap,
	       MEM_HANDLE *amountInMap,
	       MEM_HANDLE *lddMap)
{
	REAL4 **amountOut;
	CONST_REAL4_MAP amountIn;
	CONST_UINT1_MAP ldd;

	SetCaller(srcFile,srcLineNr,"upstream");
        amountOut=(REAL4 **)CpsNormalHandle(amountOutMap,CR_REAL4);
        amountIn=(CONST_REAL4_MAP)CpsNormalHandle(amountInMap,CR_REAL4);
        ldd=(CONST_UINT1_MAP)CpsNormalHandle(lddMap,CR_UINT1);

	DoUpstream(amountOut, amountIn, ldd);
}

/* Gives each cell the value of its first downstream element.
 * When the cell is a pit, it keeps its own input value.
 */
int Downstream(
		const char *srcFile,
		int    srcLineNr,
        MEM_HANDLE *amountOutMap,
	    MEM_HANDLE *amountInMap,
	    MEM_HANDLE *lddMap)
{
	REAL4 **out;
	CONST_REAL4_MAP amountIn;
	CONST_UINT1_MAP ldd;
	REAL8 	ownVal;
	int 	r, c;
 	int nrCols = (int)RgiveNrCols();
 	int nrRows = (int)RgiveNrRows();

	SetCaller(srcFile,srcLineNr,"downstream");

        out=(REAL4 **)CpsNormalHandle(amountOutMap,CR_REAL4);
        amountIn=(CONST_REAL4_MAP)CpsNormalHandle(amountInMap,CR_REAL4);
        ldd=(CONST_UINT1_MAP)CpsNormalHandle(lddMap,CR_UINT1);


	/* Fill out with missing values. This is the initial value. */
	for(r=0; r < nrRows; r++)
		SetMemMV(out[r],(size_t)nrCols,CR_REAL4);

	/* for every cell check where it flows to */
	for(r = 0; r < nrRows; r++)
	{
		for(c = 0; c < nrCols; c++)
		{
		      if( ldd[r][c] != MV_UINT1)
		      {	/* determine 1st cell downstream */
			   int lddVal = ldd[r][c];
			   int rNext = ((int)r)+LddData[lddVal].deltaY;
		       int cNext = ((int)c)+LddData[lddVal].deltaX;
			   CheckReal4Cell(amountIn,r,c, "arg.nr 2:amount-in");
		       ownVal = amountIn[r][c];

	           if ((!IS_MV_REAL4(amountIn[rNext]+cNext))
		           && ( ldd[rNext][cNext] != MV_UINT1) )
			   {
				if(lddVal != 5)
					out[r][c] = amountIn[rNext][cNext];
				else
					out[r][c] = amountIn[r][c];
			   }
		      }
			/* else-> (r,c) keeps MV as output value */
		}
	}
	return 0;
}


