#ifndef PHSOLVER_PHIO_H
#define PHSOLVER_PHIO_H

#include <FCMangle.h>

#define phio_readheader FortranCInterface_GLOBAL_(phio_readheader, PHIO_READHEADER)
#define phio_readdatablock FortranCInterface_GLOBAL_(phio_readdatablock, PHIO_READDATABLOCK)

#ifdef __cplusplus
extern "C" {
#endif
  void phio_readheader(
      int* fileDescriptor,
      const  char keyphrase[],
      void* valueArray,
      int*  nItems,
      const char  datatype[],
      const char  iotype[] );
  void phio_readdatablock(
      int* fileDescriptor,
      const  char keyphrase[],
      void* valueArray,
      int*  nItems,
      const char  datatype[],
      const char  iotype[] );
#ifdef __cplusplus
}
#endif

#endif

