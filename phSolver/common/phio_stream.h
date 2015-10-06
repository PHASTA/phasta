#ifndef PHSOLVER_PHIO_STREAM_H
#define PHSOLVER_PHIO_STREAM_H
#include "phio_base.h"
#include "phstream.h"
struct streamio_file : phio_file {
  RStream* rs;
  GRStream* grs;
};
typedef struct streamio_file* stream_fp;
void stream_openfile(
    const char filename[],
    phio_fp fileDescriptor);
void stream_readheader(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] );
void stream_writeheader( 
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const int* ndataItems,
    const char datatype[],
    const char iotype[] );
void stream_readdatablock(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] );
void stream_writedatablock(
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[]);
void stream_closefile(
    phio_fp fileDescriptor);
void stream_constructname(
    const char* in,
    char* out);
#endif
