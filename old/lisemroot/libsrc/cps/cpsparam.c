#include "stddefx.h"
#include "cps.h"

#include "misc.h"

#include <stdlib.h>
#include <string.h>   /* strcmp */

static FILE *f = NULL;
static size_t lineNr = 1;
static char fileName[256];
static char buf[256];
BOOL silentMode = FALSE;

void OpenParamFile(
	BOOL *writeEachTime,
	int argc,
	char *argv[],
	const char *modelName,
	const char *version)
{
   int i;
   int argsNoOptions=argc;

   exitOnError=1;
/* a realy silly way to find options, but we do NOT have
 * getopt on WIN32, so be it
*/
   for (i = 1 ; i < argc; i++) 
   	if (argv[i][0] == '-') {
   	  int j;
   	  for (j=1; argv[i][j] != '\0'; j++) 
   	       switch (argv[i][j]) {
   		case 'w' : *writeEachTime = TRUE;
   		           break;
   		case 's' : silentMode = TRUE;
   		           break;
 		default:
			Error("Unknown option '-%c'\n",argv[i][j]);
		}
          argsNoOptions--;
	} else {
	   f = fopen(argv[i],"r");
	   if (f == NULL)
		Error("Can't open '%s'", argv[i]);
	   strcpy(fileName, argv[i]);
	   lineNr = 1;
	}

   if (argsNoOptions != 2)
   {
   	fprintf(stderr,"Usage: %s [options] pre-file \n"
   	 " options:\n"
   	 "  -w : write end result each time step\n"
   	 "  -s : silent mode, no report of each timestep\n"
   	 "  Note: no '<' between %s and pre-file (as in older version of lisem)\n"
   	 , modelName, modelName);
   	exit(1);
   }
   fprintf(stderr,"%s (version %s)\n",modelName,version);
}

static void ReadString(void)
{
	size_t i=0;
	int c;
	while( (c = fgetc(f)) != EOF)
	{
		buf[i++] = c;
		if (c == '\n')
		{
			buf[i-1] = '\0';
			break;
		}
	}
	buf[i] = '\0';
	TokenSpaceTrim(buf);
	if (strchr(buf,' ') != NULL)
	    Error("file '%s' line '%u': contains two parameters ('%s')",
	    	fileName, lineNr, buf);
	if (buf[0] == '\0')
	    Error("file '%s' line '%u' is an empty line",
	    	fileName, lineNr);
	lineNr++;
}

void ReadStringParamFile(
	char *str)	/* string to read */
{
	ReadString();
	strcpy(str,buf);
}

void ReadNumberParamFile(
	REAL4 *val)	/* string to read */
{
	double v;
	ReadString();
	if (! CnvrtDouble(&v,buf))
	    Error("file '%s' line '%u': '%s' is not a valid number",
	    	fileName, lineNr, buf);
	*val = (float)v;
}
