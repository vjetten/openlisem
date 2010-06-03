#include "stddefx.h"
#include "cps.h"
#include "chkdata.h"

#include "misc.h"
#include <stdlib.h>
#include <string.h>   /* strcmp */



/* doit in floats, not doubles
 * otherwise roundoff errros
 * works as longs as we don't use INT4
 */
static BOOL ValidRange(
	double v_d,
	const char *ops,
	double low_d,
	double high_d)
{
	float v   = (float)v_d;
	float low = (float)low_d;
	float high = (float)high_d;
	switch(ops[0]) {
	 case '_': break;
	 case '[': if (v < low) return FALSE;
	 	   break;
	 case '<': if (v <=low) return FALSE;
	 	   break;
	 default : POSTCOND(FALSE);
       }
	switch(ops[1]) {
	 case '_': break;
	 case ']': if (v > high) return FALSE;
	 	   break;
	 case '>': if (v >=high) return FALSE;
	 	   break;
       }
       return TRUE;
}

static const char *RangeStr(
	const char *ops,
	double low,
	double high)
{
	static char buf[256];
	if (ops[0] != '_' && ops[1] != '_')
		sprintf(buf,"%c %g , %g %c",ops[0],low,high,ops[1]);
	else
	{
		switch(ops[0]) {
		 case '[': sprintf(buf,">= %g",low); break;
		 case '<': sprintf(buf,"> %g",low); break;
		 default : switch(ops[1]) {
		 	case ']': sprintf(buf,"<= %g",high); break;
		 	case '<': sprintf(buf,"< %g",high); break;
		 	default : POSTCOND(FALSE);;
		 	}
		 }
	}
	return buf;
}

static void ErrorSpatial(
	VARINFO *v,
	const char *ops,
	double low,
	double high,
	const char *parameterDescription)
{
	char mapName[128];
	char paramDescr[256];

	LeftRightTabTrim(strcpy(paramDescr,parameterDescription));

        if (strlen(paramDescr) >= 2 && paramDescr[0] == '@' && paramDescr[1] != '@')
        	Error(paramDescr+1); /* one @, print message as is */

        	Error("'%s'%s, some values are not in valid range: %s",
        		paramDescr+1,mapName,RangeStr(ops,low,high));


	if (EmptyString(v->mapName))
		mapName[0] = '\0';
	else
		sprintf(mapName," (map: '%s')",v->mapName);

        if (strlen(paramDescr) >= 2 && paramDescr[0] == '@' && paramDescr[1] == '@')
        {
                char lowTxt[128];
                char highTxt[128];
        	switch(ops[0]) {
        	 case '_' : lowTxt[0] = '\0'; break;
        	 case '<' : sprintf(lowTxt," or less than or equal to %g",low); break;
        	 case '[' : sprintf(lowTxt," or less than %g",low); break;
        	}
        	switch(ops[1]) {
        	 case '_' : POSTCOND(FALSE); break;
        	 case '>' : strcpy(highTxt," or equal to "); break;
        	 case ']' : highTxt[0] = '\0';
        	}
        	Error("'%s%s' is larger than%s cellwidth (= %g )%s",
        		paramDescr+2, mapName, highTxt, high, lowTxt);
        }
        else
        	Error("'%s'%s, some values are not in valid range: %s",
        		paramDescr,mapName,RangeStr(ops,low,high));
}

static void ErrorNonSpatial(
	const char *ops,
	double low,
	double high,
	const char *paramDescr)
{
        Error("'%s', value is not in valid range: %s",
        	paramDescr,RangeStr(ops,low,high));
}

static void CheckSpatUint1(
	VARINFO *v,
	const char *ops,
	double low,
	double high,
	const char *paramDescr)
{
	const UINT1 *d = (const UINT1 *)CpsCompressed(v);
	size_t i,n = GiveNrDefinedCells();

	for (i=0; i < n; i++)
	 if (d[i] != MV_UINT1)
	 	if (! ValidRange((double)d[i],ops,low,high))
	   		ErrorSpatial(v,ops,low,high,paramDescr);
}

static void CheckNonSpatUint1(
	VARINFO *v,
	const char *ops,
	double low,
	double high,
	const char *paramDescr)
{
	double d = (double)*(v->addr.uint1NonSpatial);

	if (! ValidRange(d,ops,low,high))
	   ErrorNonSpatial(ops,low,high,paramDescr);
}
static void CheckSpatReal4(
	VARINFO *v,
	const char *ops,
	double low,
	double high,
	const char *paramDescr)
{
	const REAL4 *d = (const REAL4 *)CpsCompressed(v);
	size_t i,n = GiveNrDefinedCells();

	for (i=0; i < n; i++)
	 if (! IS_MV_REAL4(d+i))
	 	if (! ValidRange(d[i],ops,low,high))
	   		ErrorSpatial(v,ops,low,high,paramDescr);
}

static void CheckNonSpatReal4(
	VARINFO *v,
	const char *ops,
	double low,
	double high,
	const char *paramDescr)
{
	REAL4 d = *(v->addr.real4NonSpatial);

	if (! ValidRange(d,ops,low,high))
	   ErrorNonSpatial(ops,low,high,paramDescr);
}

void CpsRangeTest(
	void *a,
	const char *ops,
	double low,
	double high,
	const char *paramDescr,
	const char *fileName,
	int  lineNr)
{
	VARINFO *v;
	if (strlen(ops) != 2)
		Error("At %s line %d: range code is not 2 symbols (e.g. [])",
			fileName, lineNr);
	if (strchr("<[_",ops[0]) == NULL)
		Error("At %s line %d: low range code is '%c', only _[< valid",
			fileName, lineNr);
	if (strchr(">]_",ops[1]) == NULL)
		Error("At %s line %d: high range code is '%c', only _[< valid",
			fileName, lineNr);
	if (ops[0] != '_' && ops[1] != '_' && low > high)
		Error("At %s line %d: low range value greater than high range value",
			fileName, lineNr);
	v = CpsFindByAddress(a);
	if (v->spatial) {
	  switch(v->type) {
	   case CR_REAL4 : CheckSpatReal4(v,ops,low,high,paramDescr); break;
	   case CR_UINT1 : CheckSpatUint1(v,ops,low,high,paramDescr); break;
	   default: PRECOND(FALSE);
	  }
	}
	else {
	  switch(v->type) {
	   case CR_REAL4 : CheckNonSpatReal4(v,ops,low,high,paramDescr); break;
	   case CR_UINT1 : CheckNonSpatUint1(v,ops,low,high,paramDescr); break;
	   default: PRECOND(FALSE);
	  }
	}
}


void CpsCellTest(
		MEM_HANDLE *ldd,
		MEM_HANDLE *m1,
		const char *srcFile,
		int    srcLineNr
      )
{
	CONST_UINT1_MAP lddMap;
	REAL4 **m1Map;
  	size_t nrRows = RgiveNrRows();
	size_t nrCols = RgiveNrCols();
	size_t r, c;
	VARINFO *v;

	SetCaller(srcFile,srcLineNr,"CellTest");

	m1Map=(REAL4 **)CpsNormalHandle(m1,CR_REAL4);
	lddMap=(CONST_UINT1_MAP)CpsNormalHandle(ldd,CR_UINT1);

	v = CpsFindByAddress(m1);

	for(r = 0; r < nrRows; r++)
	 for(c = 0; c < nrCols; c++)
  	   if ( lddMap[r][c]  != MV_UINT1) // ldd map has cell but other map may not have
      {
         switch(v->type) {
           case CR_REAL4 : CheckReal4Cell(m1Map,r,c, v->mapName); break;
           case CR_UINT1 : CheckUint1Cell(m1Map,r,c, v->mapName); break;
           default: PRECOND(FALSE);
         }
      }
      /*
      else // other map has cell but ldd has MV
      {
         switch(v->type) {
           case CR_REAL4 : if (!(IS_MV_REAL4(m1Map[r]+c)))
            Error("\n %s \nhas a value on (r,c):(%d,%d), LDD is not defined here!",v->mapName,r,c); break;
           case CR_UINT1 : if (m1Map[r][c] != MV_UINT1)
            Error("\n %s \nhas a value on (r,c):(%d,%d), LDD is not defined here!",v->mapName,r,c); break;
           default: PRECOND(FALSE);
         }
      }
      */
}


