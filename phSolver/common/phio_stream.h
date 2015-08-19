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
void stream_closefile(
    phio_fp fileDescriptor);
void stream_constructname(
    const char* in,
    char* out);
#endif
