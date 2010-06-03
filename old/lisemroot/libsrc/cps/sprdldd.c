#include "stddefx.h" 


/********/
/* USES */
/********/

/* libs ext. <>, our ""  */
#include <math.h> /* sqrt */
#include "cps.h"
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
#define IN_MAP(r,c,nrRows, nrCols) \
	( 0 <= (r) && (r) < (nrRows) && 0 <= (c) && (c) < (nrCols) ) 
#define NR_LDD_DIR 9
#define LDD_PIT    5
static struct {
     int deltaX;
     int deltaY;
       } LddData[10] = {{ 0,   0},
               {-1,   1},    /* 1 */
               { 0,   1},    /* 2 */
               { 1,   1},    /* 3 */
               {-1,   0},    /* 4 */
               { 0,   0},    /* LDD_PIT */
               { 1,   0},    /* 6 */
               {-1,  -1},    /* 7 */
               { 0,  -1},    /* 8 */
               { 1,  -1} };  /* 9 = NR_LDD_DIR */
/******************/
/* IMPLEMENTATION */
/******************/

/* Calculates the row number of the neighbor.
 * An index and the row number of the current cell are input.
 * According to the LddData structure, the
 * row number of the neighbor is calculated.
 * Returns the row number of the neighbor.
 */
static int RNeighbor( int rowNr,	/* rowNr from current cell */
 		 int index)	/* determines which neighbor */
{
 	int rNext;
 	PRECOND(0 <= index && index <= NR_LDD_DIR);
 	rNext = rowNr + LddData[index].deltaY;
 	return rNext;
}

/* Calculates the column number of the neighbor.
 * An index and the column number of the current cell are input.
 * According to the LddData structure, the
 * column number of the neighbor is calculated.
 * Returns the column number of the neighbor.
 */
static int CNeighbor( int colNr,	/* colNr from current cell */
 		 int index)	/* determines which neighbor */
{
 	int cNext;
 	PRECOND(0 <= index && index <= NR_LDD_DIR);
 	cNext = colNr + LddData[index].deltaX;
 	return cNext;
}			 

/* Spreads down from one point.
 * For this point its first downstream element is taken.
 * This is repeated until a downstream element is found with
 * already a lower cost output value in the valRes map.
 * Returns nothing, just modifies the output maps.
 */

static void SpreadDown(
     	REAL4 **valRes,		/* write-only output map  */
     	UINT1 **idRes,	        /* read-write output map */
     	CONST_UINT1_MAP points,	/* points	*/
     	CONST_UINT1_MAP ldd, 	/* ldd map	*/
     	int rowNr,		/* point from which is spread */
     	int colNr)		/* point from which is spread */
{
	REAL4 costVal1;
	UINT1 lddVal1, id1, pointVal1;
	int  nrCols, nrRows;
     	int r = rowNr;
     	int c = colNr;

	nrRows=(int)RgiveNrRows();
	nrCols=(int)RgiveNrCols();

     	while( IN_MAP(r,c,nrRows, nrCols) &&
     	       ((lddVal1 = ldd[r][c]) != MV_UINT1) &&
     	       ((pointVal1 = points[r][c]) != MV_UINT1) &&
	      (lddVal1 != LDD_PIT)
	     )
	{
		UINT1 lddVal2;
		UINT1 pntVal2, id2;
		int rNext = RNeighbor(r, lddVal1);  /* downstream elt.*/
		int cNext = CNeighbor(c, lddVal1);

		PRECOND(idRes[r][c] != MV_UINT1);

		/* Check on MV in input maps */
		if ( IN_MAP(rNext,cNext, nrRows, nrCols) &&
		    (lddVal2 = ldd[rNext][cNext]) != MV_UINT1 &&
                    (id2 = idRes[rNext][cNext]) != MV_UINT1 &&
                    (id1 = idRes[r][c]) != MV_UINT1 &&
		    (! IS_MV_REAL4(valRes[r]+c)) &&
                    (pntVal2 = points[rNext][cNext]) != MV_UINT1
                   )
		{
			REAL4 totalcost;

		        costVal1 = valRes[r][c]; 

			/* distance (so total cost also) depends on 
			 * neighbor being a corner neighbor or not.
			 * fricVal = 1
			 */
			totalcost = 1 * (lddVal1 % 2) ? sqrt(2.0) : 1;
			totalcost += costVal1; /* cost current cell */

			if((id2 == 0) ||
			   ((! IS_MV_REAL4(valRes[rNext]+cNext)) &&
			    (totalcost < valRes[rNext][cNext])) )
			{	
			        /* this path is first path or cheaper path 
			         */
				valRes[rNext][cNext] = totalcost;
				idRes[rNext][cNext] = id1;
				r = rNext; c = cNext;
			}	
			else return; /* old path was cheaper, stop */
		} else   return;	/* MV in input, stop */	
	}
}

/* Spreads along the ldd from each nonzero point in points.map.
 * Returns 0.
 */
static int SpreadLdd(
     UINT1 **idRes,			/* read-write output map  */
     CONST_UINT1_MAP points,		/* points	*/
     CONST_UINT1_MAP ldd, 		/* ldd map		*/
     double bufSize)
{
/*     REAL4 initCostVal, fricVal; always 0 and 1 */
       size_t r, c;
	size_t nrRows=RgiveNrRows();
	size_t nrCols=RgiveNrCols();


     /* create valRes (the spread value map) local */
     REAL4 **valRes = (REAL4 **)CpsMalloc2d(CR_REAL4);

     /* compute in pixels: */
     bufSize /= RgiveCellSizeX();


/* Fill idResBuf with 0, this is the initial value */
	for (r = 0; r < nrRows; r++)
	 for (c = 0; c < nrCols; c++)
	{
		if( ldd[r][c] != MV_UINT1 )
		{
		        CheckUint1Cell(points,(int)r,(int)c,"arg.nr.2: points");
			if(points[r][c] == 0)	
			{	
				idRes[r][c] = 0;
				SET_MV_REAL4(valRes[r]+c);
			}	
			else
			{	
				valRes[r][c] = 0; /* initial costs */
				idRes[r][c] = points[r][c];
			}	
		}	
	   	else
	   	{
			SET_MV_REAL4(valRes[r]+c);
			idRes[r][c] = MV_UINT1;
	   	}
	}	

	/* For every nonzero point in the pointmap perform the 
	 * spread function.
	 */
	for (r = 0; r < nrRows; r++)
	 for (c = 0; c < nrCols; c++)
	   if ( ! IS_MV_REAL4(valRes[r]+c ) )
     		 SpreadDown(valRes, idRes, points, ldd, (int)r, (int)c);

     /*
      *  now classify everything with a cost value larger than  
      */
	for (r = 0; r < nrRows; r++)
	 for (c = 0; c < nrCols; c++)
	   if ( idRes[r][c] != MV_UINT1 ) 
	   {
	        if ( (! IS_MV_REAL4(valRes[r]+c))
		          && (valRes[r][c] < bufSize) ) 
		{
	   	       idRes[r][c] = 1;
	   		if (valRes[r][c] == 0) /* point where spread is started */
	   	   		idRes[r][c] = 0;
	   	}   		
	   	else
	   	       idRes[r][c] = 2; 
           }

     /* remove/free valRes (the spread value map) local */
     CpsFree2d((void **)valRes);
    return 0;			/* successful terminated */ 
}

void SortLdd(
		const char *srcFile,
		int    srcLineNr,
      MEM_HANDLE *m_ldd,
      MEM_HANDLE *m_lddid)
{
/*
	UINT4 **lddid;
	CONST_UINT1_MAP ldd;
   size_t r, c;
	size_t nrRows=RgiveNrRows();
	size_t nrCols=RgiveNrCols();
   UINT1 value, value1, value2, value3, value4, value5, value6, value7, value8, value9;
   UINT1 checked = MV_UINT1;
   UINT1 nodata = MV_UINT1;
   long x = 0, i = 0;

	SetCaller(srcFile,srcLineNr,"sortldd");

   lddid=(UINT4 **)CpsNormalHandle(m_lddid,CR_UINT4);
   ldd=(CONST_UINT1_MAP)CpsNormalHandle(m_ldd,CR_UINT1);


   do {
 	for (r = 0; r < nrRows; r++)
	 for (c = 0; c < nrCols; c++)
	   if ( ldd[r][c] != MV_UINT1 )
      {
          if (lddid[r][c] > 0)
              continue;
          value = ldd[r][c];
          value7 = ldd[r-1][c-1];
          value8 = ldd[r-1][c];
          value9 = ldd[r-1][c+1];
          value4 = ldd[r][c-1];
          value6 = ldd[r][c+1];
          value1 = ldd[r+1][c-1];
          value2 = ldd[r+1][c];
          value3 = ldd[r+1][c+1];
          if (value7 != 3 &&
              value8 != 2 &&
              value9 != 1 &&
              value4 != 6 &&
              value6 != 4 &&
              value1 != 9 &&
              value2 != 8 &&
              value3 != 7 &&
              value != nodata)
              {
                  //map2array.cell.push_back(cells(r,c,value));
                  lddid[r][c] = value;
              }
      }

 */

}
