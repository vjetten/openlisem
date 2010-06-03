#include "stddefx.h"

#include <string.h>

/**************************************************************************/
/*  compress.c                                                            */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/**************************************************************************/

/********/
/* USES */
/********/
#include "cps.h"

/***************/
/* EXTERNALS   */
/***************/


/*********************/
/* LOCAL DEFINITIONS */
/*********************/

/**********************/
/* LOCAL DECLARATIONS */
/**********************/

/******************/
/* IMPLEMENTATION */
/******************/

static void Compressed2Normal(VARINFO *vi)
{
	char *com = (char *)(vi->addr.memH->voidSpatialComp);/* compressed data */
	void **nor = CpsMalloc2d(vi->type); /* new normal data */
	char *p=nor[0]; /* moves over nor */
	char *l = com; /* moves over com */
	size_t i,cs = CELLSIZE(vi->type); /* cell size */
	UINT1 *mask = GiveMask();
	BOOL definedSearch = TRUE;
#       ifdef DEBUG
	size_t cellsFilled = 0;
#       endif

	for(i=0; i < GiveNrMask(); i++)
	{
		if (definedSearch)
		{
			memcpy(p,l,cs*mask[i]);
			l += cs * mask[i];
		}
		else
			SetMemMV(p,mask[i],vi->type);
		p += cs * mask[i];
		definedSearch = ! definedSearch;
#       	ifdef DEBUG
		cellsFilled += mask[i];
#       	endif
	}
	POSTCOND(cellsFilled == RgiveNrCells());
	CpsFreeCompressed(com);
	vi->addr.memH->voidSpatialNorm = nor;
}

void *CpsNorm2Comp(void **nor, CSF_CR type)
{
	char *com = (char *)CpsMallocCompressed(CELLSIZE(type)); 
	            /* new compressed data */
	char *c=com; /* moves over com */
	char *n=(char *)nor[0]; /* moves over nor */
	size_t i,cs = CELLSIZE(type); /* cell size */
	UINT1 *mask = GiveMask();
	BOOL definedSearch = TRUE;

	for(i=0; i < GiveNrMask(); i++)
	{
		if (definedSearch)
		{
			memcpy(c,n,cs*mask[i]);
			c += cs * mask[i];
		}
		n += cs * mask[i];
		definedSearch = ! definedSearch;
	}
	POSTCOND(c == (com + GiveNrDefinedCells()*cs));
	POSTCOND(n == (((char *)(nor[0])) + RgiveNrCells()*cs));
	CpsFree2d(nor);
	return com;
}

static void Normal2Compressed(VARINFO *vi)
{
	char *com = CpsNorm2Comp( vi->addr.memH->voidSpatialNorm, vi->type);
	vi->addr.memH->voidSpatialComp = com;
}

void **CpsNormal(VARINFO *vi)
{

	PRECOND(vi->spatial);

	switch(vi->state) {
	 case CPS_COMPRESSED:
		Compressed2Normal(vi); break;
	 case CPS_EMPTY:
		vi->addr.memH->voidSpatialNorm= CpsMalloc2d(vi->type);
		break;
	 default: 
	 	POSTCOND(vi->state == CPS_NORMAL);
	}
	vi->state=CPS_NORMAL;
	return vi->addr.memH->voidSpatialNorm;
}

void *CpsCompressed(VARINFO *vi)
{
	PRECOND(vi->spatial);

	switch(vi->state) {
	 case CPS_COMPRESSED:
		break;
	 case CPS_EMPTY:
		vi->addr.memH->voidSpatialComp = 
			CpsMallocCompressed(CELLSIZE(vi->type));
		break;
	 default: 
		Normal2Compressed(vi); break;
	 	POSTCOND(vi->state == CPS_NORMAL);
	}
	vi->state=CPS_COMPRESSED;
	return vi->addr.memH->voidSpatialComp;
}

void CpsInitWithDataCompressed(VARINFO *vi, MEM_HANDLE data)
{
	PRECOND(vi->spatial);

	if (vi->addr.memH->voidSpatialComp == data.voidSpatialComp)
	{ /* See Save() in calc, can happen with A += 1
	   */
		POSTCOND(vi->state == CPS_COMPRESSED);
		return;
	}
	if(vi->state != CPS_EMPTY)
		CpsFreeContents(vi);
	vi->state=CPS_COMPRESSED;
	vi->addr.memH->voidSpatialComp = data.voidSpatialComp;
}

void *CpsCompressedHandle(MEM_HANDLE *m, CSF_CR cr)
{
	VARINFO *v = CpsFindByMemoryHandle(m);
	PRECOND(v->type == cr);
	return CpsCompressed(v);
}

void *CpsNormalHandle(MEM_HANDLE *m, CSF_CR cr)
{
	VARINFO *v = CpsFindByMemoryHandle(m);
	PRECOND(v->type == cr);
	return CpsNormal(v);
}
