#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <cassert>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#define __STDC_FORMAT_MACROS /* c++ requires this for the print macros */
#include <inttypes.h> /* PRIu64 */

#ifdef __bgq__
#include "hwi/include/bqc/A2_inlines.h"
#endif

#include "phiotmrc.h"
#include "phiompi.h"

#define BILLION 1000L*1000L*1000L
#define MILLION 1000L*1000L

#ifdef __INTEL_COMPILER
size_t phastaio_global_cpus;
size_t phastaio_time_diff(phioTime* start, phioTime* end);
/* return the cycle count */
void phastaio_time(phioTime* t) {
  *t = _rdtsc(); //intel intrinsic
}
/* determine the reference clock frequency */
void phastaio_setCyclesPerMicroSec() {
  const size_t usec = 5*MILLION;
  size_t cpus, cycles;
  phioTime t0, t1;
  phastaio_time(&t0);
  /* Testing on Theta indicates that 5s is long enough
   * to get a stable value for the reference frequency.
   */
  usleep(usec);
  phastaio_time(&t1);
  cycles = t1 - t0;
  cpus = ((double)cycles)/(usec);
  if(!phio_self())
    fprintf(stderr, "cycles %" PRIu64 " us %" PRIu64 " cycles per micro second %" PRIu64"\n", cycles, usec, cpus);
  phastaio_global_cpus = cpus;
}
/*return elapsed time in micro seconds*/
size_t phastaio_time_diff(phioTime* start, phioTime* end) {
  size_t cycles = *end - *start;
  size_t us = ((double)cycles)/phastaio_global_cpus;
  return us;
}
#else
void phastaio_time(phioTime* t) {
  int err;
  err = clock_gettime(CLOCK_MONOTONIC,t);
  assert(!err);
}
/*return elapsed time in micro seconds*/
size_t phastaio_time_diff(phioTime* start, phioTime* end) {
  assert(sizeof(size_t)==8);
  size_t elapsed = 0;
  phioTime diff;
  if ((end->tv_nsec-start->tv_nsec)<0) {
    diff.tv_sec = end->tv_sec-start->tv_sec-1;
    diff.tv_nsec = BILLION+end->tv_nsec-start->tv_nsec;
  } else {
    diff.tv_sec = end->tv_sec-start->tv_sec;
    diff.tv_nsec = end->tv_nsec-start->tv_nsec;
  }
  elapsed = (diff.tv_sec)*MILLION + (diff.tv_nsec)/1000L;
  return elapsed;
}
#endif

double phiotmrc (void)
{

#ifdef __bgq__

   // use the GetTimeBase function available on BGQ 
   uint64_t TB  = GetTimeBase();
   double t1 = 6.25e-10*TB; // = 1/1.6e9

#else

  // use the gettimeofday function available on any Linux plateform

  int rc;
  struct timeval tv;

  rc = gettimeofday (&tv, NULL);
  if (rc == -1) {
    fprintf(stderr,"tmrc: gettimeofday\n");
    return 0.;
  }
  double t1 =  ((double) tv.tv_sec) + 1.e-6 * ((double) tv.tv_usec);

#endif

  return t1;

}

