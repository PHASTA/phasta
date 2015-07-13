#ifndef PHSOLVER_PHIO_H
#define PHSOLVER_PHIO_H

#include <FCMangle.h>

#define phio_readheader FortranCInterface_GLOBAL_(phio_readheader, PHIO_READHEADER)
#define phio_readdatablock FortranCInterface_GLOBAL_(phio_readdatablock, PHIO_READDATABLOCK)
#define phio_openfile_read FortranCInterface_GLOBAL_(phio_openfile_read, PHIO_OPENFILE_READ)
#define phio_openfile_write FortranCInterface_GLOBAL_(phio_openfile_write, PHIO_OPENFILE_WRITE)
#define phio_restartname FortranCInterface_GLOBAL_(phio_restartname, PHIO_RESTARTNAME)

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
  void phio_openfile_read(
      const char filename[],
      int* numFiles,
      int* fileDescriptor);
  void phio_openfile_write(
      const char filename[],
      int* numFiles,
      int* numFields,
      int* numPPF,
      int* fileDescriptor);
  void phio_restartname(int* step, char* filename);
#ifdef __cplusplus
}
#endif

#endif

