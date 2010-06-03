#ifndef __TIMEOUTP
#define __TIMEOUTP


TIMESERIE *Init_Stations(UINT1 **idMap);

/*P idMap  r- the pointers to the map in memory */
/* Init_Stations will return a pointer to a timeserie with in it the tabel
	with stations*/

BOOL IsCSFTimeSerie(FILE *fp);
/*P fp r- is the filepointer to the file which has to be checked to be a CSF timeserie*/

void AddInterval(TIMESERIE *tSerie,REAL4 **dataMap,double time);

/*P tSerie  rw serie with station_coordinates and timeserie*/
/*P dataMap r- map with data to be stored */
/*P time    r- the time to be stored*/

void WriteHeader(FILE *fp,TIMESERIE *tSerie);

/*P fp r- the filepointer to the timeserie file*/
/*P tSerie r- the timeSerie which has to be used to get the number of stations*/

void TimeSerieToFile(char *fileName,TIMESERIE *tSerie);
/*P fileName r- name of TIMESERIE-file*/
/*P tSerie r- TIMESERIE which has to be written*/


void FreeTimeSerie(TIMESERIE *serie);
/*P serie r- TIMESERIE of which memory has to be freed*/

#endif /*__TIMEOUTP*/
