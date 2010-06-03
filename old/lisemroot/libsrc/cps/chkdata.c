#include "stddefx.h"
#include "cps.h"
#include "misc.h"

#include "chkdata.h"

static const char *currentSrcFile;
static int  currentSrcLineNr;
static const char *currentInFunction;

void SetCaller(
		const char *srcFile,
		int    srcLineNr,
		const char *inFunction)
{
	currentSrcFile = srcFile;
	currentSrcLineNr = srcLineNr;
	currentInFunction =inFunction;
   sprintf(ErrorMessage,"%s - %ld : %s",srcFile,srcLineNr,inFunction);
}

void CheckReal4Cell(
	CONST_REAL4_MAP m,
	int r,
	int c,
	const char *argStr)
{
	if (IS_MV_REAL4(m[r]+c))
		Error("\n %s \nhas a MV on (r,c):(%d,%d), LDD has a value here!",argStr,r,c);
//VJ 090208 simplifying error message
//		Error("%s: (%s,line %d)\n %s has a MV on (r,c):(%d,%d)\nLDD has a value here!",
//		 currentInFunction,currentSrcFile,currentSrcLineNr, argStr,r,c);

}

void CheckUint1Cell(
	CONST_UINT1_MAP m,
	int r,
	int c,
	const char *argStr)
{
	if (m[r][c] == MV_UINT1)
		Error("\n %s \nhas a MV on (r,c):(%d,%d), LDD has a value here!",argStr,r,c);
//VJ 090208 simplifying error message
//		Error("%s: (%s,line %d)\n %s has a MV on (r,c):(%d,%d)\nLDD has a value here!",
//		 currentInFunction,currentSrcFile,currentSrcLineNr, argStr,r,c);
}
