#include "stddefx.h"

#include "cps.h"
#include "misc.h"

#ifdef DEBUG
 int nrNormalMapsAllocated=0;
 int nrCompressedMapsAllocated=0;
#endif

#ifdef DEVELOP_DEBUG
# define  MAX_PTRS 500
static void *compPtrs[MAX_PTRS];
static BOOL firstTime = TRUE;
static size_t maxPtrs = 0;

static  void RegPtr(void *p)
{
	size_t i;
	if (firstTime)
	{
		printf("HEAVY DEBUG ACTIVE\n");
		for(i=0; i < MAX_PTRS; i++)
		 compPtrs[i] = NULL;
		firstTime = FALSE;
	}
	for(i=0; i < MAX_PTRS; i++)
		if (compPtrs[i] == NULL)
		{ 
		 compPtrs[i] = p; 
		 maxPtrs = MAX(maxPtrs, i);
		 return;
		}
	POSTCOND(FALSE);
}

static void UnRegPtr(void *p)
{
	size_t i;
	PRECOND(!firstTime);
	for(i=0; i < MAX_PTRS; i++)
		if (compPtrs[i] == p)
		{ compPtrs[i] = NULL; return;}
	POSTCOND(FALSE);
}

BOOL ChkCompPtr(const void *p)
{
	size_t  i;
	PRECOND(!firstTime);
	for(i=0; i < MAX_PTRS; i++)
		if (compPtrs[i] == p)
			return TRUE;
	return FALSE;
}

#else
# define UnRegPtr(x)	
# define RegPtr(x)	
#endif

void *CpsMallocCompressed(size_t elementSize)
/*P elementSize r- is the number of bytes required for each element*/
/*function returns LOCKED memory to space for a compressed map*/
{
	size_t size = GiveNrDefinedCells();
	void *p = ChkMalloc(size*elementSize);
	IFDEBUG(nrCompressedMapsAllocated++);
	RegPtr(p);
	return p;
}

void CpsFreeCompressed(void *hhMap)
{
	UnRegPtr(hhMap);
	Free(hhMap);
	IFDEBUG(nrCompressedMapsAllocated--);
}


void **CpsMalloc2d(CSF_CR cr)
{
	IFDEBUG(nrNormalMapsAllocated++);
	return Malloc2d(RgiveNrRows(),RgiveNrCols(), CELLSIZE(cr));
}

void CpsFree2d(void **data)
{
	IFDEBUG(nrNormalMapsAllocated--);
	Free2d(data,RgiveNrRows());
}

/* free data contents (a map) of an variable
 */
void CpsFreeContents(
	VARINFO *v) 
{
	if (v->spatial)
	{
	  switch(v->state) {
	   case CPS_NORMAL:
		CpsFree2d(v->addr.memH->voidSpatialNorm);
		break;
	   case CPS_COMPRESSED:
		CpsFreeCompressed( v->addr.memH->voidSpatialComp);
		break;
	   case CPS_EMPTY:
	      /*
	       *	Warning("static spatial '%s' never used",v->name);
	       */
	       break;
	  }
	}
	v->state = CPS_EMPTY;
}
