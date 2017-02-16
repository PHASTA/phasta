#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <sstream>
#include "phIO.h"
#include "phComm.h"
#include "phio_base.h"
#include <mpi.h>

#ifndef PHASTAIO_TIMERS_ON
#define PHASTAIO_TIMERS_ON 0
#endif

namespace {
  inline double getTime() {
    return MPI_Wtime();
  }
  inline bool isRankZero() {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return !rank;
  }
  inline void printTime(const char* key, double t) {
#if PHASTAIO_TIMERS_ON==1
    if( isRankZero() )
      fprintf(stderr, "%s %f\n", key, t);
#endif
  }
}

#define PHIO_TRACING 0
namespace {
  void trace(const char* key, const char* aux="", void* obj=NULL) {
    if(PHIO_TRACING)
      fprintf(stderr, "PHIO_TRACE entering %s %s %p\n", key, aux, obj);
  }
}

#ifdef __cplusplus
extern "C" {
#endif

void phio_readheader(
    phio_fp f,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  const double t0 = getTime();
  f->ops->readheader(f->file, keyphrase, valueArray,
      nItems, datatype, iotype);
  printTime(__func__, getTime()-t0);
}
void phio_writeheader(
    phio_fp f,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const int* ndataItems,
    const char datatype[],
    const char iotype[] ) {
  const double t0 = getTime();
  f->ops->writeheader(f->file, keyphrase, valueArray,
      nItems, ndataItems, datatype, iotype);
  printTime(__func__, getTime()-t0);
}
void phio_readdatablock(
    phio_fp f,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  const double t0 = getTime();
  f->ops->readdatablock(f->file, keyphrase, valueArray,
      nItems, datatype, iotype);
  printTime(__func__, getTime()-t0);
}
void phio_writedatablock(
    phio_fp f,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[]) {
  const double t0 = getTime();
  f->ops->writedatablock(f->file, keyphrase, valueArray,
      nItems, datatype, iotype);
  printTime(__func__, getTime()-t0);
}

void phio_constructName(
    phio_fp f,
    const char inName[],
    char* outName) {
  const double t0 = getTime();
  f->ops->constructname(inName, outName);
  printTime(__func__, getTime()-t0);
}

void phio_openfile(
    const char filename[],
    phio_fp f) {
  const double t0 = getTime();
  trace("openfile",filename,f);
  f->ops->openfile(filename, f);
  printTime(__func__, getTime()-t0);
}

void phio_closefile(phio_fp f) {
  const double t0 = getTime();
  trace("closefile","unknown",f);
  f->ops->closefile(f);
  printTime(__func__, getTime()-t0);
}

void phio_appendInt(char* dest, int v) {
  std::stringstream ss;
  ss << dest << v << '.';
  std::string s = ss.str();
  strcpy(dest, s.c_str());
}

#ifdef __cplusplus
}
#endif
