#ifndef __ENVIRON
#define __ENVIRON
/***************************************************************************/
/*                                                                         */
/*  PURPOSE:  To define environment and compiler specific things           */
/***************************************************************************/

/***************************************************************************/
/*                                                                         */
/*  FILE SYSTEM STUFF                                                      */
/***************************************************************************/
/* Cygnus defines: */
#ifdef _WIN32
# ifndef DOS_FS
#  define DOS_FS
# endif
# ifndef FN_LIMIT_NONE
#  define FN_LIMIT_NONE
# endif
#endif
/* Cygnus defines: */
#ifdef i386
# ifndef CPU_LITTLE_ENDIAN
#  define CPU_LITTLE_ENDIAN
# endif
#endif

#ifdef DOS_FS  
# define DIR_PATH_DELIM_CHAR  '\\'
# ifdef UNIX_FS
#  error  TWO FILESYSTEMS SPECIFIED (DOS_FS, UNIX_FS)
# endif
#else
# ifdef UNIX_FS
#  define DIR_PATH_DELIM_CHAR  '/'
#  ifdef DOS_FS
#   error  TWO FILESYSTEMS SPECIFIED (DOS_FS, UNIX_FS)
#  endif
# else
#  error  NO FILESYSTEM SPECIFIED (DOS_FS, UNIX_FS)
# endif
#endif

/***************************************************************************/
/*                                                                         */
/*  CPU STUFF                                                              */
/***************************************************************************/
#ifdef INTEL32
# ifndef CPU_LITTLE_ENDIAN
#  define CPU_LITTLE_ENDIAN
# endif
#endif

#ifdef INTEL16
# ifndef CPU_LITTLE_ENDIAN
#  define CPU_LITTLE_ENDIAN
# endif
#endif

#ifdef MOTOROLA
# ifndef CPU_BIG_ENDIAN
#  define CPU_BIG_ENDIAN
# endif
#endif

#ifndef CPU_BIG_ENDIAN
# ifndef CPU_LITTLE_ENDIAN
# error NO CPU_ENDIAN SPECIFIED (CPU_BIG_ENDIAN, CPU_LITTLE_ENDIAN)
# endif
#endif

#ifdef CPU_BIG_ENDIAN
# ifdef CPU_LITTLE_ENDIAN
#  error TWO CPU_ENDIANS SPECIFIED (CPU_BIG_ENDIAN, CPU_LITTLE_ENDIAN)
# endif
#endif

/***************************************************************************/
/*                                                                         */
/*  OS STUFF, NOT USED AS FAR AS I KNOW (CW Oct 95)                        */
/***************************************************************************/

#ifdef THINK_C       /* predefined in Think C 4.0 & 5.0*/
# define MAC 0
#endif

/***************************************************************************/
/*                                                                         */
/*  Types and constants who are equal for all environments                 */
/*                                                                         */
/***************************************************************************/

/***************************************************************************/
/* library functions that are not part of the ANSI standard but are        */
/* available on most platforms if the platform does NOT                    */
/* support these functions then use their implementation from our libraries*/
/* prototypes are given here but the misc/mathx-library must be linked in  */
/* the *_AVAILABLE functions are used to control their inclusion in the    */
/* misc or mathx library                                                   */
/* functions included and library place of our implementation are:         */
/* lfind misc                                                              */
/* hypot mathx                                                             */
/* rint/Rint  mathx re-implemented always see doc of Rint                  */
/***************************************************************************/

/* hypot must be undeffed if not available
 * on virtual platform available
 */
#	define HYPOT_AVAILABLE

#ifdef _MSC_VER
#	define lfind	_lfind
#	define LFIND_AVAILABLE
#	define hypot   _hypot
#	define stricmp   _stricmp
#	define STRICMP_AVAILABLE    /* NOT USED YET */
#endif

#ifdef BORLANDC
#	define LFIND_AVAILABLE
#endif

#ifdef DJGPP
#	define RINT_AVAILABLE 
#endif

#ifdef __linux__
#	define RINT_AVAILABLE 
#endif

#ifdef _HPUX_SOURCE
#	define RINT_AVAILABLE 
#	define LFIND_AVAILABLE 
#endif

#endif /* __ENVIRON */
