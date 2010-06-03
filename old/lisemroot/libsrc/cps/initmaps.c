#include "stddefx.h"
#include "csf.h"
/**************************************************************************/
/*      init.c                                                            */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/**************************************************************************/

/********/
/* USES */
/********/
#include "cps.h"
#include "misc.h"

/***************/
/* EXTERNALS   */
/***************/


/*********************/
/* LOCAL DEFINITIONS */
/*********************/
#define OK 1
#define NOTOK 0
#define MAPHANDLENULL "Error at initMaps\nNo valid mapHandle\nprogram halted"
#define RE_INITIALIZATION "Reinitialization of map specifications\nShould only be called once\n"
#define TOMUCHCOLS "\nThe number of columns is to big. In this case not all columns of a row\ncan be in one block. This is a problem of DOS machines which have \nsegments of 64K. If you have a 386SX-machine or higher, you can use\nthe DOS-extended 32-bits version. This version has not this problem.\nSorry.\n"
#define IS_INIT_TEST  if (! isInitialized)\
		      {\
			printf("Programmers Fault\nMap specifications not initialized\nUse InitMaps\n");\
			exit(1);\
		      };


/**********************/
/* LOCAL DECLARATIONS */
/**********************/


/******************/
/* IMPLEMENTATION */
/******************/

static CSF_VS valueScale;
static REAL8 x0;
static REAL8 y0;
static size_t nrRows;
static size_t nrCols;
static REAL8 cellSizeX;
static REAL8 cellSizeY;
static CSF_PT projection;

static UINT1 isInitialized=FALSE;

void InitMaps(
 MAP *idMapHandle) /* to map that will be used for initialization*/
{
	if (idMapHandle==NULL)
		Error("%s\n",MAPHANDLENULL);

	if (isInitialized)
		return;

	valueScale=RgetValueScale(idMapHandle);
	x0=RgetX0(idMapHandle);
	y0=RgetY0(idMapHandle);
	nrRows=RgetNrRows(idMapHandle);
	nrCols=RgetNrCols(idMapHandle);
	cellSizeX=RgetCellSize(idMapHandle);
	cellSizeY=RgetCellSize(idMapHandle);

	projection=MgetProjection(idMapHandle);

	if (Merrno!=NOERROR)
		Error("%s",MAPHANDLENULL);
	isInitialized=TRUE;
}

void ResetInitMaps(void)
{
	isInitialized=FALSE;
}

BOOL MapSpecsInit(void)
/*test if map specifications are init*/
{
	return(isInitialized);
}

CSF_VS RgiveValueScale(void)
{
	IS_INIT_TEST;
	return(valueScale);
}

REAL8 RgiveX0(void)
{
	IS_INIT_TEST;
	return(x0);
}

REAL8 RgiveY0(void)
{
	IS_INIT_TEST;
	return(y0);
}

size_t RgiveNrRows(void)
{
	IS_INIT_TEST;
	return(nrRows);
}
size_t RgiveNrCells(void)
{
	IS_INIT_TEST;
	return(nrRows *nrCols);
}

size_t RgiveNrCols(void)
{
	IS_INIT_TEST;
	return(nrCols);
}

REAL8 RgiveCellSizeX(void)
{
	IS_INIT_TEST;
	return(cellSizeX);
}

REAL8 RgiveCellSizeY(void)
{
	IS_INIT_TEST;
	return(cellSizeY);
}

CSF_PT MgiveProjection(void)
{
	IS_INIT_TEST;
	return(projection);
}

BOOL RtestSpecs(mapHandle)
MAP *mapHandle;
/*P mapHandle r- handle to file which has to be verified*/
/*function returns TRUE if map specifications are equal to
  the map specifications which are already initialized*/
{
	IS_INIT_TEST;
	if ( (RgetX0(mapHandle) != RgiveX0()) ||
	     (RgetY0(mapHandle) != RgiveY0()) ||
	     (RgetNrRows(mapHandle) != RgiveNrRows()) ||
	     (RgetNrCols(mapHandle) != RgiveNrCols()) ||
	     (RgetCellSize(mapHandle) != RgiveCellSizeX()) ||
	     (RgetCellSize(mapHandle) != RgiveCellSizeY()) ||
	     (MgetProjection(mapHandle) != MgiveProjection())
	   )   return(FALSE);
	else return(TRUE);
}
