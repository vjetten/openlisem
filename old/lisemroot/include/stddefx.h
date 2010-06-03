#ifndef __STDDEFX__
#define __STDDEFX__

#ifdef __cplusplus
 extern "C" {
#endif

/************************************************************************/
/*                                                                      */
/* Things that should be in every module                                */
/*   Only one time include                                              */
/*                                                                      */
/*                                                                      */
/************************************************************************/


#ifdef DJGPP
#define PLATFORM_TXT "dos32/dpmi"
# ifndef DOS_FS
#  define DOS_FS
# endif
# define INTEL32
# ifndef CPU_LITTLE_ENDIAN
#  define CPU_LITTLE_ENDIAN
# endif
# define  VA_START_ARG(x)	x
# include <unistd.h>
#else
# define  VA_START_ARG(x)	x
#endif

#ifdef __BORLANDC__
# ifndef BORLANDC
#  define BORLANDC
# endif
#endif

#ifdef BORLANDC
#define PLATFORM_TXT "win32"
# ifndef DOS_FS
#  define DOS_FS
# endif
# ifndef CPU_LITTLE_ENDIAN
# define CPU_LITTLE_ENDIAN
# endif
# ifndef INTEL32
#  define INTEL32
# endif
#endif


#ifdef __linux__
#define PLATFORM_TXT "linux"
# ifndef INTEL32
#  define INTEL32
# endif
# ifndef UNIX_FS
#  define UNIX_FS
# endif
# ifndef CPU_LITTLE_ENDIAN
#  define CPU_LITTLE_ENDIAN
# endif
#endif

#ifdef BORLANDC
/* including standard headers complains in C type linking */
# ifdef __cplusplus
  }
# endif
#endif

/* ANSI */
#include <stddef.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>

#ifdef BORLANDC
# ifdef __cplusplus
 extern "C" {
# endif
#endif

/* SELF */
#include "environ.h"
#include "typedef.h"
#include "debug.h"

/***********************************************/
/* Extended library functions                  */
/***********************************************/

#ifdef MIN
# if   MIN(4,2) != 2
#  error Expected that MIN macro is the minimum of two numbers
# endif
#else
# define MIN(a,b)	(((a) < (b)) ? (a) : (b))
#endif

#ifdef MAX
# if   MAX(4,2) != 4
#  error Expected that MAX macro is the maximum of two numbers
#endif
#else
# define MAX(a,b)	(((a) > (b)) ? (a) : (b))
#endif

#define ABS(a) 		(((a) > 0) ? (a) : -(a))
#define XOR		^

/* typedef for fourth argument qsort(),bsearch(), etc.. */
typedef int (*QSORT_CMP)(const void *e1, const void *e2);

/* OBSOLOTE?: #define VOID void */
/* retype VOID to char in a module
 if calculation is required */
typedef char *VPTR;
#define ARRAY_SIZE(arrayName) (sizeof(arrayName)/sizeof(arrayName[0]))

/* STUFF to shut up compiler warnings
 */
#define USED_UNINIT_ZERO   0
#define USED_UNINIT_NULL   NULL

#ifdef __cplusplus
 }
#endif

#endif /* __STDDEFX__ */
