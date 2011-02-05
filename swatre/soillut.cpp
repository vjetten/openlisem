/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Author: Victor Jetten
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
  \file soillut.cpp
  \brief SWATRE functions to deal with SWATRE look up tables (LUT)

  functions:
  - int CmpDouble(const double *e1, const double *e2); \n
  - int intervalBsearch(const void *key, const void *base, size_t num, size_t \n
         width, int (*cmp)(const void *e1, const void *e2))\n
  - void FreeLut(LUT *l); \n
  - static int Cmp(double *e1, double *e2) \n
  - double  LUT_LinIntPol(const LUT *l, int wantedCol,double keyVal,int indexCol) \n
  - int  LUT_Index_LE(const LUT *l, double keyVal, int indexCol) \n
*/

#include "swatresoillut.h"
#include "error.h"

static int keyCol;


//--------------------------------------------------------------------------------
/*! Comparison function for double
 * Usable for qsort(),bsearch(), lfind() type comparison arguments
 * returns
 * < 0 if e1 is less than e2
 * = 0 if e1 is equivalent to e2
 * > 0 if e1 is greater than e2
 */
int CmpDouble(
   const double *e1,  // pointer to single double
   const double *e2)  // pointer to single double
{
   double e1_min_e2 = (*e1)-(*e2);

   if (e1_min_e2 < 0)
      return(-1);
   return (e1_min_e2 > 0);
}
//--------------------------------------------------------------------------------
/*! RETURNS -1    if first element of array 'base' is bigger than 'key'
 *           n    where n is the element of array 'base' that is smaller than 'key'
 */
int intervalBsearch(
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

//--------------------------------------------------------------------------------
LUT *CreateLutFromContents(
   const double *lutCont,  // array with nrRows * nrCols values this pointer is grabbed, space freed by FreeLut()
   bool  gotoMinMax, //  see struct LUT definition
	int nrRows, 
	int nrCols)
{
	LUT *l;
	int i;

	l = (LUT *)malloc(sizeof(LUT));
   l->key = (double *)malloc(sizeof(double)*nrCols);
	l->lut = (const double **)malloc(sizeof(double *)*nrRows);

	for (i=0; i < nrRows; i++)
   {
      // make indirect 2d-array
		l->lut[i] = lutCont;
      // next row, skip nrCol items
      lutCont  += nrCols;
	}
	l->gotoMinMax = gotoMinMax;
	l->nrRows = nrRows;
	l->nrCols = nrCols;

	return(l);
}
//--------------------------------------------------------------------------------
void FreeLut(LUT *l)
{
   free((void *)(l->lut[0])); /* this ptr is grabbed in CreateLutFromContents:lutCont */
	free((void *)(l->lut));
	free(l->key);
	free(l);
}
//--------------------------------------------------------------------------------
/*! Comparison function for double
  Usable for qsort(),bsearch(), lfind() type comparison arguments
  returns
  < 0 if e1 is less than e2
  = 0 if e1 is equivalent to e2
  > 0 if e1 is greater than e2
 */
static int Cmp(double *e1, double *e2)
{
   return CmpDouble((e1+keyCol), (e2+keyCol));
}
//--------------------------------------------------------------------------------
double  LUT_LinIntPol(
   const LUT *l,     // lookup in this table
   int   wantedCol,  // Column number requested
   double keyVal,    // value of index Column
   int   indexCol    // find index through this column
   )
{
	int e;

   keyCol = indexCol;
	l->key[keyCol] = keyVal;
   e = intervalBsearch(l->key, l->lut[0], (size_t)l->nrRows,
                       l->nrCols*sizeof(double), (QSORT_CMP)Cmp);

	if (e == -1)
	{
		if (l->gotoMinMax)
			return  l->lut[0][ wantedCol]; /* lowest val */
      else
      {
         ErrorString = LUT_LOW(keyCol+1, keyVal);
         throw 1;
      }
	}

	if (e == (l->nrRows-1))
	{
		/* test for bigger. so equivalence is else-part */
		if (keyVal > l->lut[e][ keyCol] && !l->gotoMinMax)
      {
         ErrorString = LUT_HIGH(keyCol+1, keyVal);
         throw 1;
      }
		else
			return  l->lut[e][wantedCol]; /* highest val */
	}

   /* else: interpolate between e and e+1 */
   double lowWanted = l->lut[e][wantedCol];
   double highWanted = l->lut[(e+1)][wantedCol];
   double lowKey = l->lut[e][keyCol];
   double highKey = l->lut[(e+1)][keyCol];

   if (highKey == lowKey)
      return (lowWanted + highWanted)/2;

   return lowWanted + (highWanted-lowWanted)*(keyVal - lowKey)/(highKey - lowKey);

}
//--------------------------------------------------------------------------------
/*!  RETURNS index in LUT for that value that is most near keyVal
                   but is less or equal than keyVal
            -1     means keyVal < first entry keyVal
 */
int  LUT_Index_LE(
	const LUT *l,     
	double keyVal,  
	int    indexCol)
{

	keyCol = indexCol;
	l->key[keyCol] = keyVal;
	return intervalBsearch(l->key, l->lut[0], (size_t)l->nrRows,
	                       l->nrCols*sizeof(double), (QSORT_CMP)Cmp);
}
//--------------------------------------------------------------------------------


