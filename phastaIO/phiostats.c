#include<stdio.h>
#include"phiostats.h"
#define __STDC_FORMAT_MACROS /* c++ requires this for the print macros */
#include <inttypes.h> /* PRIu64 */
#include <unistd.h> /* usleep */
#include"phiompi.h"
#include"phiotmrc.h"

struct phastaio_stats {
  size_t readTime;
  size_t writeTime;
  size_t openTime;
  size_t closeTime;
  size_t readBytes;
  size_t writeBytes;
  size_t reads;
  size_t writes;
  size_t opens;
  size_t closes;
};
struct phastaio_stats phastaio_global_stats;

void phastaio_printSzt(const char* key, size_t v) {
  size_t min = phio_min_sizet(v);
  size_t max = phio_max_sizet(v);
  size_t tot = phio_add_sizet(v);
  double avg = ((double)tot)/phio_peers();
  if(!phio_self())
    fprintf(stderr, "phastaio_%s min max avg %" PRIu64 " %" PRIu64 " %f\n", key, min, max, avg);
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

size_t phastaio_getOpenTime() {
  return phastaio_global_stats.openTime;
}
size_t phastaio_getCloseTime() {
  return phastaio_global_stats.closeTime;
}
size_t phastaio_getReadTime() {
  return phastaio_global_stats.readTime;
}
size_t phastaio_getWriteTime() {
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
void phastaio_addReadTime(size_t time) {
  phastaio_global_stats.reads++;
  phastaio_global_stats.readTime+=time;
}
void phastaio_addWriteTime(size_t time) {
  phastaio_global_stats.writes++;
  phastaio_global_stats.writeTime+=time;
}
void phastaio_addOpenTime(size_t time) {
  phastaio_global_stats.opens++;
  phastaio_global_stats.openTime+=time;
}
void phastaio_addCloseTime(size_t time) {
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
  if(!phio_self()) {
    const size_t us = 1000;
    phioTime t0,t1;
    size_t elapsed;
    phastaio_time(&t0);
    usleep(us);
    phastaio_time(&t1);
    elapsed = phastaio_time_diff(&t0,&t1);
    fprintf(stderr, "%" PRIu64 " us measured as %" PRIu64 " us\n", us, elapsed);
  }
  const size_t reads = phio_max_sizet(phastaio_getReads());
  const size_t writes = phio_max_sizet(phastaio_getWrites());
  const size_t opens = phio_max_sizet(phastaio_getOpens());
  const size_t closes = phio_max_sizet(phastaio_getCloses());
  if(opens) {
    phastaio_printSzt("opens", phastaio_getOpens());
    phastaio_printSzt("openTime (us)",phastaio_getOpenTime());
  }
  if(closes) {
    phastaio_printSzt("closes", phastaio_getCloses());
    phastaio_printSzt("closeTime (us)",phastaio_getCloseTime());
  }
  if(reads) {
    phastaio_printSzt("reads", phastaio_getReads());
    phastaio_printSzt("readTime (us)",phastaio_getReadTime());
    phastaio_printSzt("readBytes (B)", phastaio_getReadBytes());
    /* B  * 10^6us *  1MB   = MB
     * -    ------   -----    --
     * us     1s     10^6B    s
     */
    const double bw = ((double)phastaio_getReadBytes())/phastaio_getReadTime();
    phastaio_printDbl("readBandwidth (MB/s)", bw);
    if( phastaio_getReadTime() == 0 ) {
      fprintf(stderr, "%d ZERO read time reads %" PRIu64 " readTime %" PRIu64 " readBytes %" PRIu64 " bw %f\n",
          phio_self(), phastaio_getReads(), phastaio_getReadTime(), phastaio_getReadBytes(), bw);
    }

  }
  if(writes) {
    phastaio_printSzt("writes", phastaio_getWrites());
    phastaio_printSzt("writeTime (us)", phastaio_getWriteTime());
    phastaio_printSzt("writeBytes (B)", phastaio_getWriteBytes());
    phastaio_printDbl("writeBandwidth (MB/s)",
        ((double)phastaio_getWriteBytes())/phastaio_getWriteTime());
  }
}

void phastaio_initStats() {
#ifdef __INTEL_COMPILER
  phastaio_setCyclesPerMicroSec();
#endif
  phastaio_global_stats.readTime = 0;
  phastaio_global_stats.writeTime = 0;
  phastaio_global_stats.readBytes = 0;
  phastaio_global_stats.writeBytes = 0;
  phastaio_global_stats.reads = 0;
  phastaio_global_stats.writes = 0;
  phastaio_global_stats.opens = 0;
  phastaio_global_stats.closes = 0;
}
