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
  \file swatresoillut.h
  \brief SWATRE local declarations to deal with input tables
*/

#ifndef  SOILLUT_H
#define  SOILLUT_H

#include <stdio.h>

#define THETA_COL	0
#define H_COL           1
#define K_COL           2
#define DMCH_COL        3
#define DMCC_COL        4
#define NR_COL          (DMCC_COL+1)

//---------------------------------------------------------------------------

#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))


extern int CmpDouble(const double *e1, const double *e2);

typedef int (*QSORT_CMP)(const void *e1, const void *e2);
// sort function for looking in lut tables

//---------------------------------------------------------------------------

/// SWATRE Land use tables, nrRows and nrCols mean rows and cols (3) in the table
typedef struct LUT {
	const double **lut;
	double *key;  /* buffer for search key */
	int    nrRows, nrCols;
	bool  gotoMinMax;
} LUT;


LUT *CreateLutFromContents(
	const double *lutCont, /* array with nrCols*nrRows items with values 
                             this pointer is grabbed, space is freed by FreeLut() */
   bool  gotoMinMax,      /*  see struct LUT definition */
	int nrRows, 
	int nrCols) ;

void FreeLut(LUT *l);
	/* frees space allocated in ReadLut() call */

double  LUT_LinIntPol(
	 const LUT *l,     /* lookup in this table */
	 int   wantedCol,  /* Column number requested */
	 double keyVal,    /* value of index Column */
	 int   indexCol);   /* find index through this column */
/* RET linear interpolatation between two values found or 
 *       Exits withs Error() call if keyVal is not in range if LUT 
 *         if LUT is created with gotoMinMax set to FALSE
 */

int  LUT_Index_LE(
	const LUT *l,     
	double keyVal,  
	int    indexCol);
/*  RETURNS index in LUT for that value that is most near keyVal
 *          but is less or equal than keyVal
 *          -1 means keyVal < first entry keyVal
 */

#define LUT_nrRows(l)	(l->nrRows)

double LUT_ValueAt(
	const LUT *l, 
	int   indexCol,
	int   indexRow);
/* returns value of index[Row.Col] */

#define  LUT_Highest( l, indexCol)\
         (l->lut[(LUT_nrRows(l)-1)][indexCol])
//	(LUT_ValueAt(l, indexCol, (LUT_nrRows(l)-1)))


#endif /* SOILLUT_H */
