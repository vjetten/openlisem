#include "stddefx.h"
#include "cps.h"

#include "misc.h"

#include <stdlib.h>
#include <string.h>   /* strcmp */

static BOOL warningPrinted = FALSE;


void LisemError(const char *fmt, ... )
{
    va_list marker;
    va_start(marker,VA_START_ARG(fmt));
    vfError(fmt, marker);
    va_end(marker);
}

void LisemWarning(
	const char *fmt,  /* Format control, see printf() for a description */
	... )             /* Optional arguments */
{
    va_list marker;

    warningPrinted = TRUE;
    (void)fprintf(stderr,"WARNING: ");

    /* Write text to a string and output the string. */

    va_start( marker,VA_START_ARG(fmt));
    (void)vfprintf(stderr, fmt, marker );
    va_end( marker );

    if (fmt[strlen(fmt)-1] != '\n')
		(void)fprintf(stderr, "\n");
    StartTimer();
    while (ReadTimer() < 4)
		;
}

void Report(
	const char *fmt,  /* Format control, see printf() for a description */
	... )             /* Optional arguments */
{
    va_list marker;

    if (silentMode)
    	return;

    /* Write text to a string and output the string. */

    va_start( marker,VA_START_ARG(fmt));
    (void)vprintf(fmt, marker );
    va_end( marker );
}
