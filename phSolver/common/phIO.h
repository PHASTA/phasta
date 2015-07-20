#ifndef PHSOLVER_PHIO_H
#define PHSOLVER_PHIO_H

#include <FCMangle.h>

#define phio_restartname \
  FortranCInterface_GLOBAL_(phio_restartname, PHIO_RESTARTNAME)

#ifdef __cplusplus
extern "C" {
#endif
  typedef struct phio_file* phio_fp;
  void phio_readheader(
      phio_fp fileDescriptor,
      const  char keyphrase[],
      void* valueArray,
      int*  nItems,
      const char  datatype[],
      const char  iotype[] );
  void phio_writeheader( 
      phio_fp fileDescriptor,
      const char keyphrase[],
      const void* valueArray,
      const int* nItems,
      const int* ndataItems,
      const char datatype[],
      const char iotype[] );
  void phio_readdatablock(
      phio_fp fileDescriptor,
      const  char keyphrase[],
      void* valueArray,
      int*  nItems,
      const char  datatype[],
      const char  iotype[] );
  void phio_writedatablock(
      phio_fp fileDescriptor,
      const char keyphrase[],
      const void* valueArray,
      const int* nItems,
      const char datatype[],
      const char iotype[]);
  void phio_openfile_read(
      const char filename[],
      int* numFiles,
      phio_fp* fileDescriptor);
  void phio_openfile_write(
      const char filename[],
      int* numFiles,
      int* numFields,
      int* numPPF,
      phio_fp* fileDescriptor);
  void phio_restartname(int* step, char* filename);
  void phio_closefile_read(phio_fp fileDescriptor);
  void phio_closefile_write(phio_fp fileDescriptor);
#ifdef __cplusplus
}
#endif

#endif

