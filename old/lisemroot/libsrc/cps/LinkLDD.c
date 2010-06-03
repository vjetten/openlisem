#include "stddefx.h"
#include "cps.h"
#include "misc.h"

#include "lddutil.h"
#include "chkdata.h"


static void CalculateLink(
      int pitRowNr,
      int pitColNr,
      ListStruct *list,
      RECMEM_HEAP *heap,
      CONST_UINT1_MAP ldd,
      UNIT1 ** check)
{
 int nrRows, nrCols;
 ListStruct *temp;
 heap = NewRecMemHeap(sizeof(ListStruct), 200, NULL, NULL);
 list = (ListStruct *)NewRecord(heap);
 list->prev = NULL;
 list->rowNr = pitRowNr;
 list->colNr = pitColNr;

 nrRows = (int)RgiveNrRows();
 nrCols = (int)RgiveNrCols();

 while (list != NULL)
 {
	int i;
	int rowNr = list->rowNr;
	int colNr = list->colNr;

	for (i=1; i<=9; i++)
	{
		int r,c;

		if (i==5)  /* this is the current cell*/
			continue;
		r = ((int)rowNr)+LddData[i].deltaY;
		c = ((int)colNr)+LddData[i].deltaX;
		if (r>=0 && r<nrRows &&
		    c>=0 && c<nrCols &&
		    FLOWS_TO(ldd,r,c,rowNr, colNr) &&
                    check[r][c] == 0)
		{
		    temp = NewRecord(heap);
		    temp->prev = list;
		    list = temp;
		    list->rowNr = r;
		    list->colNr = c;
                    check[r][c] = 1;
		}
	}

 }

// FreeAllRecords(heap);

}


void LinkLDD(
		const char *srcFile,
		int    srcLineNr,
		ListStruct *LDDlist,
                RECMEM_HEAP *heap,
		MEM_HANDLE *ldd)
{
	size_t nrRows   =(int) RgiveNrRows();
	size_t nrCols= RgiveNrCols();
	size_t r, c; /* (r,c) becomes co-ordinate of outflowpoint */
	CONST_UINT1_MAP lddMap;
        UINT1 **check = (UINT1 **)CpsMalloc2d(CR_UINT1);

	SetCaller(srcFile,srcLineNr,"linkLDD");

	lddMap=(CONST_UINT1_MAP)CpsNormalHandle(ldd,CR_UINT1);

	for(r = 0; r < nrRows; r++)
	 for(c = 0; c < nrCols; c++)
            check[r][c] = 0;

	for(r = 0; r < nrRows; r++)
	 for(c = 0; c < nrCols; c++)
	  if ( ldd[r][c] == 5 ) {
	   CalculateLink(r,c, LDDList,(const UINT1 **) ldd, check);

}


