#ifndef ERRORF__H
#define ERRORF__H

//#ifdef __cplusplus
// extern "C" {
//#endif

#include "csf.h"

#ifdef COMPASDOS
  #define stderr_v stderr
#else
  extern FILE *errorfile;
  #define stderr_v errorfile
#endif

//#ifdef __cplusplus
// }
//#endif

#endif /* ERRORF__H */
