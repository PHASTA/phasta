#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <cstring>
#include <string>
#include <sstream>
#include "phIO.h"
#include "phComm.h"
#include "phio_base.h"
#include "ph_mpi_help.h"


#ifndef PHASTAIO_TIMERS_ON
#define PHASTAIO_TIMERS_ON 0
#endif
struct phastaio_stats {
  double readTime;
  double writeTime;
  double openTime;
  double closeTime;
  size_t readBytes;
  size_t writeBytes;
};
phastaio_stats phio_global_stats;

namespace {
  inline double getTime() {
    return MPI_Wtime();
  }
  inline size_t getSize(const char* t) {
    std::string type(t);
    if(type == "integer")
      return sizeof(int);
    else if(type == "double")
      return sizeof(double);
    else {
      assert(0);
      exit(EXIT_FAILURE);
    }
  }
}

#define PHIO_TRACING 0
namespace {
  void trace(const char* key, const char* aux="", void* obj=NULL) {
    if(PHIO_TRACING)
      fprintf(stderr, "PHIO_TRACE entering %s %s %p\n", key, aux, obj);
  }
  void printMinMaxAvg(const char* key, size_t v) {
    int val = static_cast<int>(v);
    int min = ph_min_int(val);
    int max = ph_max_int(val);
    long tot = ph_add_long(static_cast<long>(val));
    double avg = tot/static_cast<double>(ph_peers());
    if(!ph_self())
      fprintf(stderr, "phio_%s min max avg %d %d %f\n",
          key, min, max, avg);
  }

  void printMinMaxAvg(const char* key, double v) {
    double min = ph_min_double(v);
    double max = ph_max_double(v);
    double tot = ph_add_double(v);
    double avg = tot/ph_peers();
    if(!ph_self())
      fprintf(stderr, "phio_%s min max avg %f %f %f\n",
          key, min, max, avg);
  }
}

#ifdef __cplusplus
extern "C" {
#endif

void phio_printStats() {
  const int mebi=1024*1024;
  printMinMaxAvg("readTime (s)",phio_getReadTime());
  printMinMaxAvg("writeTime (s)", phio_getWriteTime());
  printMinMaxAvg("openTime (s)", phio_getOpenTime());
  printMinMaxAvg("closeTime (s)", phio_getCloseTime());
  printMinMaxAvg("readBytes (B)", phio_getReadBytes());
  printMinMaxAvg("writeBytes (B)", phio_getWriteBytes());
  printMinMaxAvg("readBandwidth (MiB/s)",
      (phio_getReadBytes()/phio_getReadTime())/mebi);
  printMinMaxAvg("writeBandwidth (MiB/s)",
      (phio_getWriteBytes()/phio_getWriteTime())/mebi);
}

void phio_initStats() {
  phio_global_stats.readTime = 0;
  phio_global_stats.writeTime = 0;
  phio_global_stats.openTime = 0;
  phio_global_stats.closeTime = 0;
  phio_global_stats.readBytes = 0;
  phio_global_stats.writeBytes = 0;
}

double phio_getReadTime() {
  return phio_global_stats.readTime;
}

double phio_getWriteTime() {
  return phio_global_stats.writeTime;
}

double phio_getOpenTime() {
  return phio_global_stats.openTime;
}

double phio_getCloseTime() {
  return phio_global_stats.closeTime;
}

size_t phio_getReadBytes() {
  return phio_global_stats.readBytes;
}

size_t phio_getWriteBytes() {
  return phio_global_stats.writeBytes;
}

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
  phio_global_stats.readBytes += (*nItems)*getSize("integer");
  phio_global_stats.readTime += getTime()-t0;
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
  phio_global_stats.writeBytes += (*nItems)*getSize("integer");
  phio_global_stats.writeTime += getTime()-t0;
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
  phio_global_stats.readBytes += (*nItems)*getSize(datatype);
  phio_global_stats.readTime += getTime()-t0;
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
  phio_global_stats.writeBytes += (*nItems)*getSize(datatype);
  phio_global_stats.writeTime += getTime()-t0;
}

void phio_constructName(
    phio_fp f,
    const char inName[],
    char* outName) {
  f->ops->constructname(inName, outName);
}

void phio_openfile(
    const char filename[],
    phio_fp f) {
  const double t0 = getTime();
  trace("openfile",filename,f);
  f->ops->openfile(filename, f);
  phio_global_stats.openTime += getTime()-t0;
}

void phio_closefile(phio_fp f) {
  const double t0 = getTime();
  trace("closefile","unknown",f);
  f->ops->closefile(f);
  phio_global_stats.closeTime += getTime()-t0;
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
