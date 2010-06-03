
#ifndef __TYPEDEF
#define __TYPEDEF

/***************************************************************/
/*                                                             */
/*  SOME GENERAL STUFF                                         */
/*                                                             */
/***************************************************************/

#define global    /* prefix of variable declaration who are   */
                  /* referenced as extern in other modules    */
/*
 * UINT_T used for proper cast to %u printf argument
 */
typedef unsigned int  UINT_T;

typedef unsigned long ULONG_T;
typedef int    BOOL;

#define FALSE   0  /* FAILURE, NO, YOU POOR LOOSER! etc. */
#define TRUE    1  /* SUCCES , YES, drinking beer, etc */

#endif /* __TYPEDEF */
