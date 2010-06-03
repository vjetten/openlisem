#include "stddefx.h"

#include "cps.h"
#include "misc.h"
#include "string.h"


#define HEADERFAILURE "timeserie header must be RUU CSF TIMESERIE INTENSITY NORMAL <nr stations>\n"
#define NOTSUPPORTED "Timein does not yet support timeserie's of: "
#define TIME 0
#define NOTIMEINTERVAL "No time interval. Can't generate map.\n"
#define NODATAINTERVAL "No (more) input data for this time interval, stopping run.\n"
#define NOTIMPLEMENTED "Datatype of timeserie is not implemented.\nCan't generate map.\n"
#define WRONGIDMAP "Wrong ID map.\nUnknown station number: "
#define OUTOFMEMORY "Not enough memory to handle make map.\n"


TIMESERIE *timeSerie[3];


static TIMESERIE *InstallTimeSerie(
	const char *fileName)
{
	FILE *fp;
	//char word[12];
	char dummyBuf[256] ="\0";
   char *p, *p1;
	TIMESERIE *serie;
	size_t index;
   int res;

	serie=(TIMESERIE *)ChkMalloc(sizeof(TIMESERIE));

	fp=fopen(fileName,"r");
	if (fp==NULL)
		Error("can't open file: %s\nfile doesn't exist\n",fileName);

//VJ 080628 more relaxed header reading of timeseries and better checking
/*
	res = fscanf(fp,"%s",word);
	StrToUpper(word);

	if (strcmp(word,"RUU")!=0)
	  Error("%s",HEADERFAILURE);

	res = fscanf(fp,"%s",word);
	StrToUpper(word);
	if (strcmp(word,"CSF")!=0)
	  Error("%s",HEADERFAILURE);

	res = fscanf(fp,"%s",word);
	StrToUpper(word);
	if (strcmp(word,"TIMESERIE")!=0)
	  Error("%s %s",HEADERFAILURE,word);

	res = fscanf(fp,"%s",word);
	StrToUpper(word);
	if (strcmp(word,"INTENSITY")!=0)
	  Error("%s%s\n",NOTSUPPORTED,word);
	else serie->type=TSER_INTENSITY;

	res = fscanf(fp,"%s",word);
	StrToUpper(word);
	if (strcmp(word,"NORMAL")!=0)
	  Error("%s%s\n",NOTSUPPORTED,word);
	else serie->cumulatief=TSER_NORMAL;
*/
   serie->type=TSER_INTENSITY;
   serie->cumulatief=TSER_NORMAL;
   serie->numberStations = 0;

   //VJ scan first line
	fscanf(fp,"%[^\n]\n",dummyBuf);
   p = strtok(dummyBuf, " \t");
   while (p)
   {
      p1 = strdup(p);
      p = strtok(NULL, " \t");
   }
   serie->numberStations = atoi(p1);

	if (serie->numberStations <= 0)
      Error("Incorrect timeseries header format in file [%s].\nMust be\
       \"<some header text><space><nr stations>\"\nNr stations read=%d\n",
       fileName,(int)serie->numberStations);

	//res = fscanf(fp,"\n",dummyBuf);

	/*scan names of stations, do nothing with them */
	for (index=1;index<=serie->numberStations;index++)
		fscanf(fp,"%[^\n]\n",dummyBuf);

	serie->lines=0;

	if (feof(fp)==0)
	{
	  serie->timeSerie=(REAL4 **) ChkMalloc (sizeof(REAL4 *));
	  serie->time=(REAL4 *) ChkMalloc (sizeof(REAL4));
	  serie->timeSerie[0]=(REAL4 *) ChkMalloc ((serie->numberStations) * sizeof(REAL4));

	  res = fscanf(fp,"%f",&(serie->time[0]));

	  for (index=0;index<serie->numberStations;index++)
			fscanf(fp,"%f",&(serie->timeSerie[0][index]));

	  res = fscanf(fp,"\n",dummyBuf);
	}


	while (feof(fp)==0)
	{
	  (serie->lines)++;
	  serie->timeSerie=ChkRealloc(serie->timeSerie,((serie->lines+1)*sizeof(REAL4 *)));
	  serie->time=ChkRealloc(serie->time,((serie->lines+1)*sizeof(REAL4)));

	  serie->timeSerie[0]=(REAL4 *)ChkRealloc(serie->timeSerie[0],((serie->numberStations) * sizeof(REAL4) * (serie->lines+1)));

	  for (index=0;index<(serie->lines+1);index++)
		serie->timeSerie[index]=(REAL4 *) ((char *) serie->timeSerie[0]+(serie->numberStations * sizeof(REAL4) * index));


	  res = fscanf(fp,"%f",&(serie->time[serie->lines]));

	  for (index=0;index < (serie->numberStations);index++)
		fscanf(fp,"%f",&(serie->timeSerie[(serie->lines)][index]));

	  res = fscanf(fp,"\n",dummyBuf);
	}

	fclose(fp);
 /* for debugging
   fp = fopen("try.txt","w");
   for (index = 0; index < serie->lines; index++)
    fprintf(fp,"%f %f\n",serie->time[index],serie->timeSerie[index][0]);
   fclose(fp);
 */  

	(serie->lines)++; //???? why?
	return(serie);
}


static void Serie2map(REAL4 *intervalMap,TIMESERIE *serie,UINT1 *idMap,double intervalEnd, double intervalLenght)
/*P serie  r- pointer to timeserie structure */
/*P idMap  r- map with for each point whre it belongs to */
/*P intervalEnd   r- end of interval */
/*P intervalLenght r- lenght of interal */
{
	REAL4 *valueArea;
	size_t index,count;
	double startInterval;
	size_t i,nrCells = GiveNrDefinedCells();

	if (serie->lines==0)
		Error("%s",NOTIMEINTERVAL);

        startInterval=intervalEnd-intervalLenght;

   //VJ 050809 waar slaat dit op?
   //     si = (startInterval*1000000)+1;
   //	    startInterval=si;
   //     startInterval/=1000000;

	if (startInterval<serie->time[0])
		Error("%s",NODATAINTERVAL);


    //added 1e-6 to avoid rounding of errors
	index=0;
	while ((index < serie->lines) && (startInterval > serie->time[index] + 1e-6))
		index++;

	if (index >= serie->lines)
		Error("%s",NODATAINTERVAL);




	valueArea=(REAL4 *) ChkMalloc (serie->numberStations * sizeof(REAL4));
	/* array with values for stations */

	if (serie->type==TSER_INTENSITY && serie->cumulatief==TSER_NORMAL)
	{

	  if (serie->time[index] > intervalEnd)
	  {
		for (count=0;count<serie->numberStations;count++)
		 valueArea[count]=serie->timeSerie[index][count]*intervalLenght;
	  }
	  else
	  {
		for (count=0;count<serie->numberStations;count++)
		 valueArea[count]=serie->timeSerie[index][count]*(serie->time[index]-startInterval);

		index++;

		while ((index < serie->lines) && (intervalEnd > serie->time[index]) )
		{
		 for (count=0;count < serie->numberStations;count++)
		  valueArea[count]+=serie->timeSerie[index][count]* (serie->time[index] - serie->time[index-1]);
		  index++;
		}

		if (index < serie->lines)
		{
		 for (count=0;count<serie->numberStations;count++)
		  valueArea[count]+=serie->timeSerie[index][count] * (intervalEnd-serie->time[index-1]);
		}
		else
		{
		  Free(valueArea);
		  Error("%s",NODATAINTERVAL);
		}
	  }
	  for (count=0;count<serie->numberStations;count++)
		valueArea[count]=valueArea[count]/intervalLenght;

		/*to calculate the value for the interval*/
	}
	else
	{
		/* if not INTENSITY and NORMAL */
		printf("intensity %d    normal %d\n",serie->type,serie->cumulatief);
		Error("%s",NOTIMPLEMENTED);
	}




	/**********generate map for interval************/


	for (i=0;i<nrCells;i++)
	{
		if (idMap[i]==MV_UINT1)
		   SET_MV_REAL4(intervalMap+i);
		else
		{
		  if (idMap[i] > ((UINT1) serie->numberStations))
		  {
		    /* for (count=0;count<nrRows;count++) CW WHY HERE ?
		     *  Free(intervalMap[count]);
		     * Free(intervalMap);
		     */
		     Error("%s cell (rownr,colnr):(%d,%d) %d\n",
		        WRONGIDMAP,(int)i,(int)i, /* CW convert to row,col */
			(int)idMap[i]);
		  }
		  else
		     intervalMap[i]=valueArea[idMap[i]-1];
		}
	}
	Free(valueArea);
}

void CpsTimeIn(MEM_HANDLE *vResultMap,MEM_HANDLE *vIdMap,
     const char *serieFileName,double intervalEnd, double intervalLength, BOOL FirstTime, int TSnumber)
/*P vResultMap -w is the map which will contain the new data*/
/*P vIdMap r- is the idmap .. has to be UINT1 !!!!!!!*/
/*P serieFileName r- is the used timeserie filename*/
/*P intervalEnd r- end of wanted interval*/
/*P intervalLength r- interval Length*/
{
	UINT1 *idMap;
	REAL4 *resultMap;
//	static BOOL FIRSTTIME=TRUE;


	if (FirstTime)
	{
		 timeSerie[TSnumber]=InstallTimeSerie(serieFileName);
	}

	idMap=(UINT1 *) CpsCompressedHandle(vIdMap,CR_UINT1);
	resultMap=(REAL4 *)CpsCompressedHandle(vResultMap,CR_REAL4);

	Serie2map(resultMap, timeSerie[TSnumber], idMap, intervalEnd,intervalLength);
}

void FreeTimeSeries(TIMESERIE *serie)
/*P serie r- TIMESERIE of which memory has to be freed*/
{
    size_t index;
    if (serie)
    {
//        for (index=0;index<serie->lines;index++)
//                Free(serie->timeSerie[index]);
        if(serie->timeSerie)
        {
          Free(serie->timeSerie[0]); /*this is one block, now reallocated*/
          Free(serie->timeSerie);
        }

        if(serie->time)
          Free(serie->time);

        Free(serie);
        serie = NULL;
    }
}

