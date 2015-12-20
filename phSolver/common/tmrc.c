#include "tmrc.h"
#include <stdio.h>
#include <sys/types.h>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#ifdef __bgq__
#include "hwi/include/bqc/A2_inlines.h"
#endif

double TMRC (void)
{

#ifdef __bgq__

   /* use the GetTimeBase function available on BGQ */
   uint64_t TB  = GetTimeBase();
   double t1 = 6.25e-10*TB; /* = 1/1.6e9 */

#else

  /* use the gettimeofday function available on any Linux platform */

  int rc;
  struct timeval tv;
  double t1 = 0;

  rc = gettimeofday (&tv, NULL);
  if (rc == -1) {
    fprintf(stderr,"tmrc: gettimeofday\n");
    return 0.;
  }
  t1 =  ((double) tv.tv_sec) + 1.e-6 * ((double) tv.tv_usec);

#endif

  return t1;
}
