#ifndef PHSOLVER_PHIO_H
#define PHSOLVER_PHIO_H

#ifdef __cplusplus
extern "C" {
#endif
  typedef struct phio_file* phio_fp;
  void phio_openfile(
      const char filename[],
      phio_fp fileDescriptor);
  void phio_closefile(phio_fp fileDescriptor);
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
  void phio_constructName(
      phio_fp fileDescriptor,
      const char* inName,
      char* outName);
#ifdef __cplusplus
}
#endif

#endif

