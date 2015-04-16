#include <FCMangle.h>
#define TMRC FortranCInterface_GLOBAL_(tmrc, TMRC)

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


  // Old stuff
  /*
  struct rusage cputimePoints;
  getrusage(RUSAGE_SELF,&cputimePoints);
  double t1=((cputimePoints.ru_utime).tv_sec+(cputimePoints.ru_stime).tv_sec);
  t1+=((cputimePoints.ru_utime).tv_usec+(cputimePoints.ru_stime).tv_usec)/1000000.0;
  */

   /* double t1=rts_get_timebase()/( 2500000000.0); */

   /* return time(NULL); */
   /* return clock()/CLOCKS_PER_SEC; */
}




