#include "stddefx.h"
#include "cps.h"
#include "misc.h"

#include "lddutil.h"
#include "chkdata.h"


static void CalculateAccuflux(
      int pitRowNr,    /* r- Y-coordinate of OutflowPoint */
      int pitColNr,    /* r- X-coordinate of OutflowPoint */
      REAL4   **m_Q1, /* -w debiet at time j+1 */
CONST_REAL4_MAP m_Q0, /* debiet at time j   */
CONST_UINT1_MAP ldd)
{
/***************************************************************************/
/* this funtion creates a list with a LDD path of all cells draining to the*/
/* cell colNr,rowNr, starting with the outflow point. For each cell in that*/
/* path the fluxes Q1 of the 8 neighbour cells is summed in Qin and a new  */
/* Q1 is estimated. The last cell is the outflow point.                    */
/***************************************************************************/
 int nrRows, nrCols;
 Liststruct *list, *temp;
 RECMEM_HEAP *heap = NewRecMemHeap(sizeof(Liststruct), 200, NULL, NULL);
 list = (Liststruct *)NewRecord(heap);
 list->prev = NULL;
 list->rowNr = pitRowNr;
 list->colNr = pitColNr;

 nrRows = (int)RgiveNrRows();
 nrCols = (int)RgiveNrCols();

 while (list != NULL)
 {
	int i;
	BOOL  subCachDone = TRUE; /* are sub-catchment cells done ? */
	int rowNr = list->rowNr;
	int colNr = list->colNr;
       /* put all points that have to be calculated to
        * calculate the current point in the list,
	* before the current point
	*/
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
		    IS_MV_REAL4(&m_Q1[r][c]) ) /* cell not computed */
		{
		    temp = NewRecord(heap);
		    temp->prev = list;
		    list = temp;
		    list->rowNr = r;
		    list->colNr = c;
		    subCachDone = FALSE;
		}
	}


	/* no other points are found that drain to point [rowNr,colNr] */

	if (subCachDone)
	{
		REAL8 Qin=0.0;

        	/* find for each point in list the neighbour
         	 * LDD points and sum Q1, starting with last point of prev
         	 */

		for (i=1;i<=9;i++) /* for all incoming cells */
		{
			int r = ((int)rowNr)+LddData[i].deltaY;
			int c = ((int)colNr)+LddData[i].deltaX;

			if (i==5)  /* Skip current cell */
				continue;

			if (r>=0 && r < nrRows &&
			    c>=0 && c < nrCols &&
			    FLOWS_TO(ldd,r,c,rowNr, colNr) &&
			    !IS_MV_REAL4(&m_Q0[r][c]))
			 {
				Qin += m_Q1[r][c];
             }
		} /* eof all incoming cells */

		temp=list;
   		list=list->prev;
		FreeRecord(temp,heap);

         // for one cell in the list iterate Qj + sum Qupstreamj to Qtimej+1
		m_Q1[rowNr][colNr] = (REAL4)Qin+m_Q0[rowNr][colNr];

		/* cell rowNr, colNr is now done */
	}/* eof subcatchment done */
 } /* eowhile list != NULL */

 FreeAllRecords(heap);
}


void Accuflux(
		const char *srcFile,
		int    srcLineNr,
		MEM_HANDLE *Q1,
		MEM_HANDLE *Q0,
		MEM_HANDLE *ldd)
{
	REAL4 **Q1Map;
	CONST_REAL4_MAP Q0Map;
	CONST_UINT1_MAP lddMap;
	size_t nrRows   =(int) RgiveNrRows();
	size_t nrCols= RgiveNrCols();
	size_t r, c; /* (r,c) becomes co-ordinate of outflowpoint */

	SetCaller(srcFile,srcLineNr,"accuflux");

	Q1Map=(REAL4 **)CpsNormalHandle(Q1,CR_REAL4);
	Q0Map=(CONST_REAL4_MAP)CpsNormalHandle(Q0,CR_REAL4);
	lddMap=(CONST_UINT1_MAP)CpsNormalHandle(ldd,CR_UINT1);

        /* mark all cells as not done: */
	for(r = 0; r < nrRows; r++)
	   SetMemMV(Q1Map[r],nrCols,CR_REAL4);

	for(r = 0; r < nrRows; r++)
	 for(c = 0; c < nrCols; c++)
	  if ( lddMap[r][c] == 5 )
      {
	   CalculateAccuflux(r,c, Q1Map, Q0Map, (const UINT1 **) lddMap);
      }

}



