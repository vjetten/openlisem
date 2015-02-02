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


#define THETA_COL	0
#define H_COL           1
#define K_COL           2
#define DMCH_COL        3
#define DMCC_COL        4
#define NR_COL          (DMCC_COL+1)

#define LUT_nrRows(l)	           (l->nrRows)
#define LUT_ValueAt(l,indexCol,indexRow)  (l->lut[indexRow][indexCol])
#define LUT_Highest(l, indexCol)   (l->lut[(LUT_nrRows(l)-1)][indexCol])

#define INC_LUT 90

/* error messages: */
#define LUT_LOW(a,b)  QString("In lookup table: lowest value of column %1 is larger than requested value (%2)").arg(a).arg(b)
#define LUT_HIGH(a,b) QString("In lookup table: highest value of column %1 is smaller than requested value (%2)").arg(a).arg(b)


//---------------------------------------------------------------------------

/// comparison function for sort
extern int CmpDouble(const double *e1, const double *e2);

/// sort function for looking in lut tables
typedef int (*QSORT_CMP)(const void *e1, const void *e2);

//---------------------------------------------------------------------------

/// SWATRE Land use tables, nrRows and nrCols mean rows and cols (3) in the table
typedef struct LUT {
	const double **lut;
   double *key;  // buffer for search key
   int   nrRows, nrCols;
	bool  gotoMinMax;
} LUT;


LUT *CreateLutFromContents(
   const double *lutCont, // array with nrCols*nrRows items, this pointer is grabbed, space is freed by FreeLut()
   bool  gotoMinMax,      //  see struct LUT definition
	int nrRows, 
   int nrCols);

/*! frees space allocated in ReadLut() call */
void FreeLut(LUT *l);

/*! linear interpolatation between two values found or exits
    withs Error() call if keyVal is not in range if LUT.
    if LUT is created with gotoMinMax set to FALSE
 */
double  LUT_LinIntPol(
   const LUT *l,     /* lookup in this table */
   int   wantedCol,  /* Column number requested */
   double keyVal,    /* value of index Column */
   int   indexCol);   /* find index through this column */

/*!  RETURNS index in LUT for that value that is most near keyVal
                   but is less or equal than keyVal
            -1     means keyVal < first entry keyVal
 */
int  LUT_Index_LE(
	const LUT *l,     
	double keyVal,  
	int    indexCol);


#endif /* SOILLUT_H */
