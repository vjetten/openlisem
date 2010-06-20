/*---------------------------------------------------------------------------
project: openLISEM
name: lisModel.cpp
author: Victor Jetten
licence: GNU General Public License (GPL)
Developed in: MingW/Qt/Eclipse
website SVN: http://sourceforge.net/projects/lisem

---------------------------------------------------------------------------*/
/**************************************************************************/
/* lutio.c                                                                */
/*                                                                        */
/*   read a theta,h,k-table for Swatre component of LISEM                 */
/*                                                                        */
/**************************************************************************/
#include <ctype.h>

#include "error.h"
#include "lutio.h"
#include "soillut.h"

/**********************/
/*  number of elements added to the lut
 *  when initializing (malloc) or resizing
 *  (realloc)
 */
#define EXP_NR_COLS 3
#define INC_STEP (45*3)

/* error messages: */
#define OPEN_ERRORs "SWATRE: Can't open '%s'"
#define READ_ERRORs "SWATRE: Read error on '%s'"
#define COL_ERRORs  "SWATRE: Table '%s', entry nr. %d contains %s than 3 columns"
#define EOF_MESSs   "SWATRE: Unexpected end of file while reading lookup table"
#define NO_ENTRIESs "SWATRE: Table '%s' contains no entries"
#define NR_COLSs    "SWATRE: Encountered line containing '%s' while first row "\
"had %d columns"
#define NAN_MESSs   "SWATRE: Table: '%s' contains a non number symbol: '%s'"
#define SMALLERs    "SWATRE: Table '%s' column '%s' on entry '%d' has smaller "\
		"value ('%g') than previous element "

//----------------------------------------------------------------------------------------
		static void ReadCols(
				const char *fileName, /* for error reporting only */
				double *inLut,   /* -w current position in lut that will be filled */
				const char *buf, /* buffer to read from */
				int   n);        /* number of items to read */
//----------------------------------------------------------------------------------------
/* Trim string and replace space by single space, and count tokens.
 * Removes leading and trailing isspace() characters and
 * substitutes sequences of isspace() chararacters with one
 * space (that is ' ' not '\t').
 * A token is a string of non-isspace() characters terminated by a space
 * or '\0'
 * Returns the number of tokens.
 */
int TokenSpaceTrim(
		char *s)  /* read-write. String to be modified and counted */
{
	int i;    /* index over s */
	int d;    /* destination index */
	int t=0;  /* #tokens */

	/* remove leading spaces */
	for(i=0; isspace(s[i]); i++)
	{
		/* inc i is all we want */;
	}
	/* copy string */
	for(d=0; s[i] != '\0'; )
	{
		if (isspace(s[i]))
		{
		   s[d++] = ' ';
		   t++;
		   while (isspace(s[i]) )
				i++;
		}
		else
			s[d++] = s[i++];
	}
	/* adjust for trailing spaces */
	if (isspace(s[d-1]))
	{
		d--;
		t--;
	}
	/* put string terminator */
	s[d] = '\0';
	t = ( d == 0 ) ? 0 : t + 1;
	return(t);
} /* TokenSpaceTrim */
//----------------------------------------------------------------------------------------
double *ReadSoilTable(
		const char *fileName,
		int   *nrRows          /* nr of rows read */
      )
		/* every row is on a single line
 * EOF is end of table 
 * EXAMPLE
   AT THIS MOMENT ONLY TOTALLY SORTED LUTS ARE SUPPORTED
  */
{
	char buf[1024];
	double  *l;
	int     sizeL, nrL; 
	FILE    *f;

	sizeL = INC_STEP;
	nrL = 0;
	l = (double *)malloc(sizeof(double) * sizeL);
	f = fopen(fileName, "r");
	if (f == NULL)
		Error(OPEN_ERRORs,fileName,0);

	do {
		int currNrCols;
		if (fgets(buf, 1024, f) == NULL)
		{
			if (feof(f))
				break; /* OK, END OF FILE */
			Error(READ_ERRORs, fileName,0);
		}

		currNrCols = TokenSpaceTrim(buf);
		// trims leading trailing spaces
		if (currNrCols == 0)
			continue;  /* EMPTY LINE, next one please */
		if (currNrCols != EXP_NR_COLS)
			Error(COL_ERRORs,fileName, (nrL/EXP_NR_COLS)+1);
		//(currNrCols < EXP_NR_COLS) ? "less" : "more");
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
					Error(SMALLERs, fileName, i);//colName[i]);//,(nrL/EXP_NR_COLS), l[nrL+i] );

		nrL += EXP_NR_COLS;
	} while (/* infinite: */ nrL > -1 );
	/* loop contains one break and one continue */
	if (nrL == 0)
		Error(NO_ENTRIESs,fileName,0);

	*nrRows = nrL / EXP_NR_COLS;
	fclose(f);
	return l;
}
//----------------------------------------------------------------------------------------
static void ReadCols(
		const char *fileName, /* for error reporting only */
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
