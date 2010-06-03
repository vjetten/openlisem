#ifndef __TOUTMAIN
#define __TOUTMAIN


TIMESERIE *Init_Stations_Smart(char *fileName);

/*P fileName r- is the name of the seriefile where the station-coordinates are in*/



BOOL IsCSFTimeSerie(FILE *fp);

/*P fp r- the filepointer of the file which has to be checked on CSF format */
/*        filepointer has to be at beginning of file */

void WriteHeader(FILE *fp,TIMESERIE *tSerie);

/* fp r- filepointer of the file where the header should be written to */
/*       fp has to be at right position*/

#endif /*__TOUTMAIN*/
