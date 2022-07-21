/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011,2020  Victor Jetten
**  contact:
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
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
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/

/*!
\file lutio.cpp
\brief  SWATRE read theta, h, k-table

functions:
- double* TWorld::ReadSoilTable(const char *fileName, int *nrRows); \n
- void TWorld::ReadCols(const char *fileName, double *inLut, const char *buf); \n
*/


#include "swatresoillut.h"
#include "lerror.h"
#include "model.h"


///  number of elements added to the lut when initializing (malloc) or resizing (realloc)
#define EXP_NR_COLS 3
#define INC_STEP (45*3)

/// error messages:
#define OPEN_ERRORs "SWATRE: Can't open %1"
#define READ_ERRORs "SWATRE: Read error on %1"
#define COL_ERRORs  "SWATRE: Table %1, entry nr. %2 contains %3 than 3 columns"
#define EOF_MESSs   "SWATRE: Unexpected end of file while reading lookup table"
#define NO_ENTRIESs "SWATRE: Table %1 contains no entries"
#define NR_COLSs    "SWATRE: Encountered line containing %1 while first row had %2 columns"
#define NAN_MESSs   "SWATRE: Table %1 contains a non number symbol: %2"
#define SMALLERs    "SWATRE: Table %1 column %2 on entry %3 has smaller value (%4) than previous element"

static const char *colName[3] = { "theta", "h", "k" };

//----------------------------------------------------------------------------------------
double *TWorld::ReadSoilTable(
		const char *fileName,
      int   *nrRows
      )
{
	char buf[1024];
	double  *l;
	int     sizeL, nrL;
	FILE    *f;

	sizeL = INC_STEP;
	nrL = 0;
	l = (double *)malloc(sizeof(double) * sizeL);
	f = fopen(fileName, "r");
    if (f == nullptr)
		Error(QString(OPEN_ERRORs).arg(fileName));

	do {
		int currNrCols;
        if (fgets(buf, 1024, f) == nullptr)
		{
			if (feof(f))
				break; /* OK, END OF FILE */
			Error(QString(READ_ERRORs).arg(fileName));
		}

      QStringList SL = QString(buf).split(QRegExp("\\s+"),Qt::SkipEmptyParts);
      currNrCols = SL.count();
      strcpy(buf, SL.join(" ").toLatin1());
      // trim spaces and count columns

		if (currNrCols == 0)
			continue;  /* EMPTY LINE, next one please */

		if (currNrCols != EXP_NR_COLS)
			Error(QString(COL_ERRORs).arg(fileName).arg(nrL/EXP_NR_COLS+1).arg(currNrCols < EXP_NR_COLS?"less":"more"));
		/* increase l if neccessary */
		if (sizeL <= (nrL + EXP_NR_COLS))
		{
			sizeL += INC_STEP;
			l = (double *)realloc(l, sizeL*sizeof(double));
		}

		ReadCols(fileName, l+nrL,buf,EXP_NR_COLS);
		/* test if lut is monotonous increasing */

		if (nrL > 0)
			for (int i=0; i< EXP_NR_COLS; i++)
				if ( l[nrL+i] < l[(nrL-EXP_NR_COLS)+i])
               Error(QString(SMALLERs).arg(fileName).arg(colName[i]).arg(nrL/EXP_NR_COLS).arg(l[nrL+i]));

		nrL += EXP_NR_COLS;
	} while (/* infinite: */ nrL > -1 );

	/* loop contains one break and one continue */
	if (nrL == 0)
		Error(QString(NO_ENTRIESs).arg(fileName));

	*nrRows = nrL / EXP_NR_COLS;
	fclose(f);
	return l;
}
//----------------------------------------------------------------------------------------
void TWorld::ReadCols(
        const char * /* fileName */, /* for error reporting only */
		double *inLut,   /* -w current position in lut that will be filled */
		const char *buf, /* buffer to read from */
		int   n)         /* number of items to read */
{
	const char *cp; /* current ptr */
	char  *rp;      /* result prt */
	int i;

	cp = buf;
	for (i = 0; i < n; i++)
	{
		inLut[i] = strtod(cp, &rp);
		//	if ( (!isspace(*rp)) && (*rp != '\0') )
		//	Error(NAN_MESS, fileName, cp);
		cp = rp; /* advance ptr to next number */
	}
}
//----------------------------------------------------------------------------------------
