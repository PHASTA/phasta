#include <stdio.h>
#include <sys/types.h>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#ifdef __bgq__
#include "hwi/include/bqc/A2_inlines.h"
#endif

#include "phiotmrc.h"

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

