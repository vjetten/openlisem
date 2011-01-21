/*---------------------------------------------------------------------------
project: openLISEM
name: lisModel.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/ 
website SVN: http://sourceforge.net/projects/lisem
---------------------------------------------------------------------------*/
/**************************************************************************/
/* soillut.c                                                              */
/*                                                                        */
/*   Soil characteristics look-up table for Swatre component of LISEM     */
/*                                                                        */
/**************************************************************************/

#include <string.h>
#include <ctype.h>
#include <math.h>

#include "swatresoillut.h"
//#include "swatremisc.h"

#define INC_LUT 90

/* error messages: */
#define LUT_LOW    "In lookup table: lowest value of column %d is larger than "\
   " requested value (%g)"
#define LUT_HIGH    "In lookup table: highest value of column %d is smaller "\
   "than requested value (%g)"

static int keyCol;

<<<<<<< .mine
//--------------------------------------------------------------------------------
/* Comparison function for double
 * Usable for qsort(),bsearch(), lfind() type comparison arguments
 * returns
 * < 0 if e1 is less than e2
 * = 0 if e1 is equivalent to e2
 * > 0 if e1 is greater than e2
 */
int CmpDouble(
   const double *e1,  /* pointer to single double */
   const double *e2)  /* pointer to single double */
{
/* see cmpdoubl.s for GNU def */
//  register double e1_min_e2 = (*e1)-(*e2);
  double e1_min_e2 = (*e1)-(*e2);
  if (e1_min_e2 < 0)
   return(-1);
  return (e1_min_e2 > 0);
}
//--------------------------------------------------------------------------------

int intervalBsearch(
		/* arguments like std. ANSI. bsearch() */
		const void *key, const void *base, size_t num, size_t
		width, int (*cmp)(const void *e1, const void *e2))
=======
>>>>>>> .r114
/* RETURNS -1     if first element of array 'base' is bigger  then 'key'
 *          n     where n is the element of array 'base' that is smaller
 *                than 'key'
 */
int intervalBsearch(
   /* arguments like std. ANSI. bsearch() */
   const void *key, const void *base, size_t num, size_t
   width, int (*cmp)(const void *e1, const void *e2))
{
	int l,r,x, cmpVal;
	void *a;
	l=0;
	r=(int)num;
	do {
      //      lVal = *(double *)(((char *)base)+(l*width));
      //      rX = (r == (int)num) ? r-1 : r;
      //      rVal = *(double *)(((char *)base)+(rX*width));

		x = (l+r)/2;
		if (x == (int)num) /* reached the end, return last element */
		   return x-1;

		a = ((char *)base)+(x*width);
      //      xVal = *((double *)a);
		/* debug code patch here */
		cmpVal = cmp(key, a);
      //                cmpVal = 0;
      //              if (key < a) cmpVal = -1;
      //              if (key > a) cmpVal = 1;
      if (cmpVal < 0)
		   r = x-1;
		else
         l = x+1;
	} while (cmpVal != 0 && l <= r);
	if (cmpVal < 0 || x == (int)num)
		return(x-1);
	return(x);
}


LUT *CreateLutFromContents(
   const double *lutCont,  /* array with nrRows * nrCols values
                           this pointer is grabbed, space freed by FreeLut() */
	bool  gotoMinMax, /*  see struct LUT definition */
	int nrRows, 
	int nrCols)
{
	LUT *l;
	int i;

	l = (LUT *)malloc(sizeof(LUT));
	l->key = (double *)malloc(sizeof(double)*nrCols);
	l->lut = (const double **)malloc(sizeof(double *)*nrRows);
	for (i=0; i < nrRows; i++)
	{  /* make indirect 2d-array */
		l->lut[i] = lutCont;
		lutCont  += nrCols; /* next row, skip nrCol items */
	}
	l->gotoMinMax = gotoMinMax;
	l->nrRows = nrRows;
	l->nrCols = nrCols;

	return(l);
}


void FreeLut(LUT *l)
{

	free((void *)(l->lut[0])); /* this ptr is grabbed 
        * in CreateLutFromContents:lutCont
        */
	free((void *)(l->lut));
	free(l->key);
	free(l);
}

double LUT_ValueAt(
	const LUT *l, 
	int   indexCol,
	int   indexRow)
{
	return  l->lut[indexRow][indexCol];
}
//--------------------------------------------------------------------------------

/* Comparison function for double
 * Usable for qsort(),bsearch(), lfind() type comparison arguments
 * returns
 * < 0 if e1 is less than e2
 * = 0 if e1 is equivalent to e2
 * > 0 if e1 is greater than e2
 */
static int Cmp(double *e1, double *e2)
{
   double e1_min_e2 = (*e1+keyCol)-(*e2+keyCol);
   if (e1_min_e2 < 0)
      return(-1);
   return (e1_min_e2 > 0);
}

// return CmpDouble((e1+keyCol), (e2+keyCol));


//--------------------------------------------------------------------------------
double  LUT_LinIntPol(
   const LUT *l,     /* lookup in this table */
   int   wantedCol,  /* Column number requested */
   double keyVal,    /* value of index Column */
   int   indexCol    /* find index through this column */
   )
{
	int e;
	keyCol = indexCol;
	l->key[keyCol] = keyVal;
	e = intervalBsearch(l->key, l->lut[0], (size_t)l->nrRows,l->nrCols*sizeof(double),
                       (QSORT_CMP)Cmp);

	if (e == -1)
	{
		if (l->gotoMinMax)
			return  l->lut[0][ wantedCol]; /* lowest val */
      //	else
		//	Error(LUT_LOW, keyCol+1, keyVal);
	}
	if (e == (l->nrRows-1))
	{
		/* test for bigger. so equivalence is else-part */
		if (keyVal > l->lut[e][ keyCol] && !l->gotoMinMax)
         ;//     		Error(LUT_HIGH, keyCol+1, keyVal);
		else
			return  l->lut[e][wantedCol]; /* highest val */
	}

	{ /* else: interpolate between e and e+1 */
		double lowWanted = l->lut[e][wantedCol];
		double highWanted = l->lut[(e+1)][wantedCol];
		double lowKey = l->lut[e][keyCol];
		double highKey = l->lut[(e+1)][keyCol];
		double dRel;
		if (highKey == lowKey)
			return (lowWanted+highWanted)/2;
		dRel = (keyVal - lowKey)/(highKey - lowKey);
		return lowWanted+ ((highWanted-lowWanted)*dRel);
	}
}


double  LUT_LinIntPol1(
   const LUT *l,     /* lookup in this table */
   int   wantedCol,  /* Column number requested */
   double keyVal,    /* value of index Column */
   int   indexCol    /* find index through this column */
   )
{
	int e;

	keyCol = indexCol;
	l->key[keyCol] = keyVal;
	e = intervalBsearch(l->key, l->lut[0], (size_t)l->nrRows,l->nrCols*sizeof(double),
                       (QSORT_CMP)Cmp);

	if (e == -1)
	{
		if (l->gotoMinMax)
			return  l->lut[0][ wantedCol]; /* lowest val */
		else
         ;//		Error(LUT_LOW, keyCol+1, keyVal);
	}
	if (e == (l->nrRows-1))
	{
		/* test for bigger. so equivalence is else-part */
		if (keyVal > l->lut[e][ keyCol] && !l->gotoMinMax)
         ;//		Error(LUT_HIGH, keyCol+1, keyVal);
		else
			return  l->lut[e][wantedCol]; /* highest val */
	}

	{ /* else: interpolate between e and e+1 */
		double lowWanted = l->lut[e][wantedCol];
		double highWanted = l->lut[(e+1)][wantedCol];
		double lowKey = l->lut[e][keyCol];
		double highKey = l->lut[(e+1)][keyCol];
		double dRel;
		if (highKey == lowKey)
			return (lowWanted+highWanted)/2;
		dRel = (keyVal - lowKey)/(highKey - lowKey);
      /*
 *		if(dRel < 0)
 *		{
 *			int r,c;
 *			printf("keyCol %d, wantCol %d, e %d, keyVal %g, "
 *			       " high %g, low %g, dRel %g\n", keyCol, 
 *			       wantedCol, e, keyVal, highKey, lowKey, dRel);
 *			  DumpLut("bug.lut", l);
 *		}
 */
		return lowWanted+ ((highWanted-lowWanted)*dRel);
	}
}

int  LUT_Index_LE(
	const LUT *l,     
	double keyVal,  
	int    indexCol)
/*  RETURNS index in LUT for that value that is most near keyVal
 *          but is less or equal than keyVal
 *          -1 means keyVal < first entry keyVal
 */
{

	keyCol = indexCol;
	l->key[keyCol] = keyVal;
	return intervalBsearch(l->key, l->lut[0], (size_t)l->nrRows,
	                       l->nrCols*sizeof(double), (QSORT_CMP)Cmp);
}


