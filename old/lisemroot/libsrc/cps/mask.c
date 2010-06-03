#include "stddefx.h"
#include <stdio.h>
#include <stdlib.h>
#include "misc.h"
#include "cps.h"


/***********************************************************************/
/***********************************************************************/

#define MASKMAPERROR "Error while reading mask-map\n"
#define MASK_INIT_TEST if (!maskInitialized) Error("InitMask() not yet called\n")
#define REALLOCFAILED "\nreallocation of maskblock failed\n"
#define MASKOUTOFMEMORY "\nallocating memory for mask failed\n"
#define NOPOINTS " not a single point on mask map\nDon't be so silly\nI refuse to work with this mask\nExit\n"

/************************MASK INFO*********************/
static UINT1 *mask;  
static size_t nrMask;
static size_t nrDefinedCells=0;
static BOOL maskInitialized=FALSE;
static MAP  *maskM;
/******************************************************/

UINT1 *GiveMask()
/*function returns the mask*/
{
	MASK_INIT_TEST;
	return(mask);
}

size_t GiveNrDefinedCells()
/*function returns number of defined cells on 2dmap*/
{
	MASK_INIT_TEST;
	return(nrDefinedCells);
}

size_t GiveNrMask()
{
	MASK_INIT_TEST;
	return(nrMask);
}

void FreeMask()
/*function free's all memory allocated to hold the mask-info
  called at exit*/
{
	if (maskInitialized)
	{
		mask = GiveMask();
		Free(mask);
//VJ??
//        maskInitialized = FALSE;
	}
	Mclose(maskM);
}


void InitMask(char *name)
/*P name r- is the name of the file which will be read*/
{
      size_t i,nrCells;
      size_t number,offset ;
      BOOL definedSearch;
      void **maskMap;
      REAL4 *linBuf;
      IFDEBUG(size_t totalFound = 0);

      PRECOND(!maskInitialized); /* CW RUNTIME ERROR */
      	/*there can only be one mask at this time*/
      maskInitialized=TRUE;

      
      maskMap = CpsOpenAndReadAs2d(name,CR_REAL4);
      if ((maskM = Mopen(name,M_READ)) == NULL)
       MperrorExit(name,1);	
	
      linBuf= ((REAL4 **)maskMap)[0];
      nrCells = RgiveNrCells();

      mask=(UINT1 *) ChkMalloc(nrCells*sizeof(UINT1 *));

      nrDefinedCells=0;
      /* global static variable which holds 
       * the number of defined cells
       */
      definedSearch=TRUE; /* flip between searches 
                           * start with defined
                           */
      number=0;
      offset=0;

      for (i=0;i < nrCells; i++)
      {
	    nrDefinedCells += IS_MV_REAL4(linBuf+i) ? 0 : 1;
	    if (IS_MV_REAL4(linBuf+i) && definedSearch)
	    {
		mask[offset]=number;
		IFDEBUG(totalFound += number);
		number=0;
		offset++;
		definedSearch=FALSE;
	    }
	    if (!(IS_MV_REAL4(linBuf+i)) && (!definedSearch))
	    {
		mask[offset]=number;
		IFDEBUG(totalFound += number);
		number=0;
		offset++;
		definedSearch=TRUE;
	    }
	    number++;
	    if (number==255) /*maximimum of variable number*/
	    {
		mask[offset]=number;
		IFDEBUG(totalFound += number);
		number=0;
		offset++;
		definedSearch=!definedSearch;
	    }
	}
	if (number != 0) /* stuff pending */
	{
		mask[offset]=number;
		IFDEBUG(totalFound += number);
		number=0;
		offset++;
	}
	POSTCOND(totalFound == nrCells);

      CpsFree2d(maskMap);
      nrMask = offset;
      mask = (UINT1 *)ChkRealloc(mask,sizeof(UINT1)*nrMask);
      POSTCOND(mask != NULL);
}

