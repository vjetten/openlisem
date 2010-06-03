#include "stddefx.h"
#include "csf.h"


/**************************************************************************/
/*  timeout.c                                                             */
/*                                                                        */
/*                                                                        */
/*                                                                        */
/**************************************************************************/

/********/
/* USES */
/********/
#include <stdlib.h>
#include <string.h>
#include "cps.h"
#include "misc.h"


/***************/
/* EXTERNALS   */
/***************/


/*********************/ 
/* LOCAL DEFINITIONS */
/*********************/ 


/**********************/
/* LOCAL DECLARATIONS */
/**********************/
#define NOTCSFFILE "Given seriefile is not a RUU CSF TIMESERIE file\nSorry.\n"

/******************/
/* IMPLEMENTATION */
/******************/


static TIMESERIE *Init_Stations(
	UINT1 **idMap)

/*P idMap  r- the pointers to the map in memory */
/* Init_Stations will return a pointer to a timeserie with in it the tabel
	with stations*/
{
size_t nrRows,nrCols,colNr,rowNr;
TIMESERIE *tSerie;

	nrRows=(size_t) RgiveNrRows();
	nrCols=(size_t) RgiveNrCols();

	tSerie=(TIMESERIE *) ChkMalloc (sizeof(TIMESERIE));

	tSerie->numberStations=0;    /* init all members of structure */
	tSerie->station=NULL;
	tSerie->cumulatief=TSER_UNDEF;
	tSerie->type=TSER_UNDEF;
	tSerie->lines=0;
	tSerie->time=NULL;
	tSerie->timeSerie=NULL;


	for (rowNr=0;rowNr<nrRows;rowNr++)
	  for (colNr=0;colNr<nrCols;colNr++)
	  {
		if (idMap[rowNr][colNr] != MV_UINT1)
      if (idMap[rowNr][colNr] > 0)
		{
		     tSerie->station=ChkRealloc(tSerie->station,((tSerie->numberStations+1) * sizeof(COORDINATES *)));
		     tSerie->station[tSerie->numberStations]=(COORDINATES *) ChkMalloc(sizeof(COORDINATES));
		     tSerie->station[tSerie->numberStations]->rowNr=rowNr;
		     tSerie->station[tSerie->numberStations]->colNr=colNr;
           tSerie->station[tSerie->numberStations]->id = (UINT4) idMap[rowNr][colNr];
		     tSerie->numberStations++;
		}
	  }

	return(tSerie);
}




static void AddInterval(TIMESERIE *tSerie,REAL4 **dataMap,double time)

/*P tSerie  rw serie with station_coordinates and timeserie*/
/*P dataMap r- map with data to be stored */
/*P time    r- time to be stored*/
{
   size_t stationNr;


	tSerie->timeSerie=ChkRealloc(tSerie->timeSerie,((tSerie->lines+1)*sizeof(REAL4 *)));
	tSerie->time=ChkRealloc(tSerie->time,((tSerie->lines+1)*sizeof(REAL4)));

	tSerie->timeSerie[tSerie->lines]=(REAL4 *) ChkMalloc (tSerie->numberStations * sizeof(REAL4));

	tSerie->time[tSerie->lines]=time;


	for (stationNr=0;stationNr<tSerie->numberStations;stationNr++)
		tSerie->timeSerie[tSerie->lines][stationNr]=dataMap[tSerie->station[stationNr]->rowNr][tSerie->station[stationNr]->colNr];
	tSerie->lines++;
}




static BOOL IsCSFTimeSerie(FILE *fp)

/*P fp r- is the filepointer to the file which has to be checked to be a CSF timeserie*/
{
char word[18];

	if (fgets(word,18,fp)!=NULL)     /*fgets reads 18-1 characters plus the end of string*/
	{
		StrToUpper(word);
		if(strcmp("RUU CSF TIMESERIE",word)==0) return(TRUE);
	}

	/*in all other cases*/
	return(FALSE);

}



static void WriteHeader(FILE *fp,TIMESERIE *tSerie)

/*P fp r- the filepointer to the timeserie file*/
/*P tSerie r- the timeSerie which has to be used to get the number of stations*/
{
	size_t stationNr;

	fprintf(fp,"RUU CSF TIMESERIE NOTYPE NOTYPE \n%d\n",tSerie->numberStations);
	for (stationNr=0;stationNr<tSerie->numberStations;stationNr++)
    	fprintf(fp,"location id %d\n",tSerie->station[stationNr]->id);

/*
	   fprintf(fp,"%d %d\n",
	     (int)(tSerie->station[stationNr]->rowNr),
	     (int)(tSerie->station[stationNr]->colNr));
   	fprintf(fp,"\n");
*/

}


static void TimeSerieToFile(const char *fileName,TIMESERIE *tSerie)
/*P fileName r- name of TIMESERIE-file*/
/*P tSerie r- TIMESERIE which has to be written*/
{
FILE *fp;
char dummychar;
size_t stationNr;

	fp=fopen(fileName,"a+");

	if (fp==NULL)
		Error("can't open file '%s'\n",fileName);

	/*check if file exists. if not than write header else compare to be CSF*/
	if ((dummychar=getc(fp))==EOF)
		WriteHeader(fp,tSerie);
	else
	{
		ungetc(dummychar,fp);
		if (!IsCSFTimeSerie(fp)) Error(NOTCSFFILE);
	}

	fseek(fp,0l,SEEK_END);

	fprintf(fp,"%g",tSerie->time[0]);
	for (stationNr=0;stationNr<tSerie->numberStations;stationNr++)
	   fprintf(fp," %g",tSerie->timeSerie[0][stationNr]);
	fprintf(fp,"\n");
	fclose(fp);
}

static void FreeTimeSerie(TIMESERIE *serie)
/*P serie r- TIMESERIE of which memory has to be freed*/
{
size_t index;

	/*for (index=0;index<serie->lines;index++)
		Free(serie->timeSerie[index]);*/
	Free(serie->timeSerie[0]); /*this is one block, now reallocated*/
	Free(serie->timeSerie);

	for (index=0;index<serie->numberStations;index++)
		Free(serie->station[index]);
	Free(serie->station);

	Free(serie->time);

	Free(serie);
}


void CpsTimeOutput(const char *fileName,MEM_HANDLE *vIdMap,MEM_HANDLE *vDataMap,double time)
/*P fileName r- timeserie file*/
/*P vIdMap r- idmap where stations are on*/
/*P vDataMap r- map with data to be stored*/
/*P time r- time of writing*/
{
	TIMESERIE *serie;
	REAL4 **dataMap;
	UINT1 **idMap;

	/* CW could be compressed if COORDINATE stuff is changed
	 */
	idMap=(UINT1 **) CpsNormalHandle(vIdMap,CR_UINT1);
	dataMap=(REAL4 **)CpsNormalHandle(vDataMap,CR_REAL4);

	serie=Init_Stations(idMap);
	AddInterval(serie,dataMap,time);

	TimeSerieToFile(fileName,serie);

	FreeTimeSerie(serie);
}
