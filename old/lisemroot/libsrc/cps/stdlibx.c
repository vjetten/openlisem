#include "stddefx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>     /*strlen*/
#include <ctype.h>      /*islower,toupper*/
#include <math.h>	/*floor*/

#include "cps.h"
#include "misc.h"



void StrToUpper(char *str)
{
       size_t i,n;
       n = strlen(str);
	for (i=0;i<n;i++)
		str[i]= (char)toupper((str[i]));
}


char *CatPath(const char *fileName,const char *path)
/*P fileName rw fileName which has to be concatenated to the path*/
/*P path r- name of path which has to be in front of filename*/
/*  function adds a path to a filename in a buffer*/
/*  returned filename exists until!! next call CatPath*/
{
static char buffer[128];
       int i = 0;
	strcpy(buffer,path);
        if (buffer[strlen(path)-1] != DIR_PATH_DELIM_CHAR)
        {
	   buffer[strlen(path)] = DIR_PATH_DELIM_CHAR;
           i = 1;
        }
	strcpy(buffer+strlen(path)+i,fileName);

	return(buffer);
}

char *MakeFileNameSeq(char *fileName,int nr)
/*P fileName r- is the base filename (first 5 characters of filename) */
/*P minutes r- is the rest of the filename (extention is number of seconds*/
/* function adds number to filename*/
/* returned filename exists until!! next call to MakeFileName*/
{
    static char buffer[128];
    char number[5];
	sprintf(number,".%03u",nr);
	strcpy(buffer,fileName);

    strcat(buffer,number);

	return(buffer);
}


char *MakeFileName(char *fileName,float minutes)
/*P fileName r- is the base filename (first 5 characters of filename) */
/*P minutes r- is the rest of the filename (extention is number of seconds*/
/* function adds number to filename*/
/* returned filename exists until!! next call to MakeFileName*/
{
static char buffer[128];
char number[9];
unsigned int min,frac;
float minf;

	min=(unsigned int) floor(minutes);
        minf = floor((minutes-min)*1000 + 0.5);
	frac=(unsigned int) minf;
	sprintf(number,"%04u.%03u",min,frac);

	strcpy(buffer,fileName);
        strcat(buffer,number);

	return(buffer);
}

char *__MakePCRFileName(char *iname,int  atTime)
{
      /*
       *	----------1---
       *	012345678901
       *	xxxxxxxx.xxx
       */
	char filePartBuf[128];
	char dirPartBuf[128];
	static char result[128], nrs[14];
	char *dirPart,*filePart,*lastPoint;
        SplitFilePathName(iname, &dirPart,&filePart);
        strcpy(filePartBuf,filePart);
        strcpy(dirPartBuf,dirPart);
	lastPoint = strrchr(filePartBuf,'.');
	if (strlen(filePartBuf) > 11 || (lastPoint != NULL && strlen(lastPoint) >= 4))
	{ /* this is a hack, feature if a calc script
	   * has a dynamic report of 8+3 length then at
	   * each timestep the same map is written,
	   * result: only keep the last step
           */
		return MakeFilePathName(result,dirPartBuf,filePartBuf);
	}
	(void)sprintf(nrs,"%011d",atTime);
	/* move overlapping last 3 digits + '\0'
	 * to insert '.'
	 */
	(void)memmove(nrs+9,nrs+8,4);
	nrs[8] = '.';
	/* overprint nrs with the name part
	 * excluding '\0'
	 */
	(void)memcpy(nrs,filePartBuf,strlen(filePartBuf));
	return MakeFilePathName(result,dirPartBuf,nrs);
}

char *MakePCRFileName(char *iname,int  atTime)
{
      /*
       *	----------1---
       *	012345678901
       *	xxxxxxxx.xxx
       */
	char filePartBuf[128];
	char dirPartBuf[128];
	static char result[128], nrs[14];
	char *dirPart,*filePart,*lastPoint;

        SplitFilePathName(iname, &dirPart,&filePart);
        strcpy(filePartBuf,filePart);
        strcpy(dirPartBuf,dirPart);
	lastPoint = strchr(filePartBuf,'.');
        if (lastPoint != NULL)
           lastPoint='\0';
        if (strlen(filePartBuf) > 8)
           filePartBuf[8]='\0';
    //VJ: cutof filename at 8 chars for new display

	(void)sprintf(nrs,"%011d",atTime);
	/* move overlapping last 3 digits + '\0'
	 * to insert '.'
	 */
	(void)memmove(nrs+9,nrs+8,4);
	nrs[8] = '.';
	/* overprint nrs with the name part
	 * excluding '\0'
	 */
	(void)memcpy(nrs,filePartBuf,strlen(filePartBuf));
	return MakeFilePathName(result,dirPartBuf,nrs);
}
