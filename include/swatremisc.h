#ifndef MISC__H
#define MISC__H

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#ifdef __cplusplus
 extern "C" {
#endif

#define DIR_PATH_DELIM_CHAR  '\\'
#ifndef _MAX_
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))
#endif

/* qsortcmp.c */
 typedef int (*QSORT_CMP)(const void *e1, const void *e2);

extern int CmpUchar(const unsigned char *e1, const unsigned char *e2);
extern int CmpInt(const int *e1, const int *e2);
extern int CmpFloat(const float *e1, const float *e2);
extern int CmpDouble(const double *e1, const double *e2);


#ifdef __cplusplus
 }
#endif

#endif /* MISC__H */
