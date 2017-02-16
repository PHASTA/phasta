#include<stdio.h>
#include"phiostats.h"
#include"phiompi.h"

struct phastaio_stats {
  double readTime;
  double writeTime;
  double openTime;
  double closeTime;
  size_t readBytes;
  size_t writeBytes;
  size_t reads;
  size_t writes;
  size_t opens;
  size_t closes;
};
struct phastaio_stats phastaio_global_stats;

void phastaio_printSzt(const char* key, size_t v) {
  int val = (int)v;
  int min = phio_min_int(val);
  int max = phio_max_int(val);
  long tot = phio_add_long((long)val);
  double avg = tot/(double)phio_peers();
  if(!phio_self())
    fprintf(stderr, "phastaio_%s min max avg %d %d %f\n",
        key, min, max, avg);
}

void phastaio_printDbl(const char* key, double v) {
  double min = phio_min_double(v);
  double max = phio_max_double(v);
  double tot = phio_add_double(v);
  double avg = tot/phio_peers();
  if(!phio_self())
    fprintf(stderr, "phastaio_%s min max avg %f %f %f\n",
        key, min, max, avg);
}

double phastaio_getOpenTime() {
  return phastaio_global_stats.openTime;
}
double phastaio_getCloseTime() {
  return phastaio_global_stats.closeTime;
}
double phastaio_getReadTime() {
  return phastaio_global_stats.readTime;
}
double phastaio_getWriteTime() {
  return phastaio_global_stats.writeTime;
}
size_t phastaio_getReadBytes() {
  return phastaio_global_stats.readBytes;
}
size_t phastaio_getWriteBytes() {
  return phastaio_global_stats.writeBytes;
}
void phastaio_addReadBytes(size_t bytes) {
  phastaio_global_stats.readBytes+=bytes;
}
void phastaio_addWriteBytes(size_t bytes) {
  phastaio_global_stats.writeBytes+=bytes;
}
void phastaio_addReadTime(double time) {
  phastaio_global_stats.reads++;
  phastaio_global_stats.readTime+=time;
}
void phastaio_addWriteTime(double time) {
  phastaio_global_stats.writes++;
  phastaio_global_stats.writeTime+=time;
}
void phastaio_addOpenTime(double time) {
  phastaio_global_stats.opens++;
  phastaio_global_stats.openTime+=time;
}
void phastaio_addCloseTime(double time) {
  phastaio_global_stats.closes++;
  phastaio_global_stats.closeTime+=time;
}
size_t phastaio_getReads() {
  return phastaio_global_stats.reads;
}
size_t phastaio_getWrites() {
  return phastaio_global_stats.writes;
}
size_t phastaio_getOpens() {
  return phastaio_global_stats.opens;
}
size_t phastaio_getCloses() {
  return phastaio_global_stats.closes;
}

void phastaio_printStats() {
  const int mebi=1024*1024;
  const int reads = phio_max_int(phastaio_getReads());
  const int writes = phio_max_int(phastaio_getWrites());
  const int opens = phio_max_int(phastaio_getOpens());
  const int closes = phio_max_int(phastaio_getCloses());
  if(opens) {
    phastaio_printSzt("opens", phastaio_getOpens());
    phastaio_printDbl("openTime (s)",phastaio_getOpenTime());
  }
  if(closes) {
    phastaio_printSzt("closes", phastaio_getCloses());
    phastaio_printDbl("closeTime (s)",phastaio_getCloseTime());
  }
  if(reads) {
    phastaio_printSzt("reads", phastaio_getReads());
    phastaio_printDbl("readTime (s)",phastaio_getReadTime());
    phastaio_printSzt("readBytes (B)", phastaio_getReadBytes());
    phastaio_printDbl("readBandwidth (MiB/s)",
        (phastaio_getReadBytes()/phastaio_getReadTime())/mebi);
  }
  if(writes) {
    phastaio_printSzt("writes", phastaio_getWrites());
    phastaio_printDbl("writeTime (s)", phastaio_getWriteTime());
    phastaio_printSzt("writeBytes (B)", phastaio_getWriteBytes());
    phastaio_printDbl("writeBandwidth (MiB/s)",
        (phastaio_getWriteBytes()/phastaio_getWriteTime())/mebi);
  }
}

void phastaio_initStats() {
  phastaio_global_stats.readTime = 0;
  phastaio_global_stats.writeTime = 0;
  phastaio_global_stats.readBytes = 0;
  phastaio_global_stats.writeBytes = 0;
  phastaio_global_stats.reads = 0;
  phastaio_global_stats.writes = 0;
  phastaio_global_stats.opens = 0;
  phastaio_global_stats.closes = 0;
}
