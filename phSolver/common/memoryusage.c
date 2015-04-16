/* #include <stdio.h> */
#include <FCMangle.h> 
#include "memoryusage.h"
#include "common_c.h"
#include "mpi.h"

#ifdef __bgp__  // For Intrepid/Challenger
#include <spi/kernel_interface.h>
#elif __bgq__  // For Mira/Cetus
#include <spi/include/kernel/memory.h>
#endif


void initmemstat( void) {

  memstats.rheap = 0.0;
  memstats.rheapavail = 0.0;
  memstats.rstack = 0.0;
  memstats.rstackavail = 0.0;
  memstats.rshared = 0.0;
  memstats.rpersist = 0.0;
  memstats.rguard = 0.0;
  memstats.rmmap = 0.0;

} 

void printmeminfo( const char *msg) {


#ifdef __bgp__ 
    unsigned int heapnew, heapavailnew, stacknew, stackavailnew;
//    unsigned int sharednew, persistnew, guardnew, mmapnew;
#elif __bgq__
    uint64_t heapnew, heapavailnew, stacknew, stackavailnew;
//    uint64_t sharednew, persistnew, guardnew, mmapnew;
#endif



#ifdef __bgp__ 
// copy the same lines
#elif __bgq__

  // Print info for the following ranks
  if (workfc.myrank == 0 || workfc.myrank == (workfc.numpe/2) || workfc.myrank == (workfc.numpe-1) ) {
  //for (int irank=0 ; irank<workfc.numpe ; irank++){
  // if(irank == workfc.myrank) {

    Kernel_GetMemorySize( KERNEL_MEMSIZE_HEAP, &heapnew);    
    Kernel_GetMemorySize( KERNEL_MEMSIZE_HEAPAVAIL, &heapavailnew);    
    Kernel_GetMemorySize( KERNEL_MEMSIZE_STACK, &stacknew);    
    Kernel_GetMemorySize( KERNEL_MEMSIZE_STACKAVAIL, &stackavailnew);    

    double rheapnew, rheapavailnew, rstacknew, rstackavailnew;
    rheapnew = (double) heapnew*inv1024sq;
    rheapavailnew = (double) heapavailnew*inv1024sq;
    rstacknew = (double) stacknew*inv1024sq;
    rstackavailnew = (double) stackavailnew*inv1024sq;

    double rdiffheap, rdiffheapavail, rdiffstack, rdiffstackavail;
    rdiffheap = rheapnew - memstats.rheap;
    rdiffheapavail = rheapavailnew - memstats.rheapavail;
    rdiffstack = rstacknew - memstats.rstack;
    rdiffstackavail = rstackavailnew - memstats.rstackavail;

    printf("MEM %s (rank %d ): allocated heap: %6.2f MB ( %6.2f MB), avail. heap: %6.2f MB (%6.2f MB), allocated stack: %6.2f MB ( %6.2f MB), avail. stack: %6.2f MB ( %6.2f MB),\n", 
                msg, workfc.myrank, rheapnew, rdiffheap, rheapavailnew, rdiffheapavail, rstacknew, rdiffstack, rstackavailnew, rdiffstackavail );
  
    memstats.rheap = rheapnew;
    memstats.rheapavail = rheapavailnew;
    memstats.rstack = rstacknew;
    memstats.rstackavail = rstackavailnew ;

/*
    Kernel_GetMemorySize( KERNEL_MEMSIZE_SHARED, &sharednew);
    Kernel_GetMemorySize( KERNEL_MEMSIZE_PERSIST, &persistnew);    
    Kernel_GetMemorySize( KERNEL_MEMSIZE_GUARD, &guardnew);    
    Kernel_GetMemorySize( KERNEL_MEMSIZE_MMAP, &mmapnew);    

    double rsharednew, rpersistnew, rguardnew, rmmapnew;
    rsharednew = (double) sharednew*inv1024sq;
    rpersistnew = (double) persistnew*inv1024sq;
    rguardnew = (double) guardnew*inv1024sq;
    rmmapnew = (double) mmapnew*inv1024sq;

    double  rdiffshared, rdiffpersist, rdiffguard, rdiffmmap;
    rdiffshared = rsharednew - memstats.rshared;
    rdiffpersist = rpersistnew - memstats.rpersist;
    rdiffguard = rguardnew - memstats.rguard;
    rdiffmmap = rmmapnew - memstats.rmmap;

    printf("MEM %s (rank %d ): allocated heap: %6.2f MB ( %6.2f MB), avail. heap: %6.2f MB (%6.2f MB)\n",msg, workfc.myrank, rheapnew, rdiffheap, rheapavailnew, rdiffheapavail );
    printf("MEM %s (rank %d ): shared: %6.2f MB, persist.: %6.2f MB ( %6.2f MB), guard: %6.2f MB ( %6.2f MB), mmap: %6.2f MB ( %6.2f MB) \n", msg, workfc.myrank, rsharednew, rdiffshared, 
          rpersistnew, rdiffpersist, rguardnew, rdiffguard, rmmapnew, rdiffmmap);

    memstats.rshared = rsharednew;
    memstats.rpersist = rpersistnew;
    memstats.rguard = rguardnew;
    memstats.rmmap = rmmapnew;
*/

//   } // if irank == myrank
//   MPI_Barrier(MPI_COMM_WORLD);
//  } // for irank=1 -> numpe

  } // if irank == chosen ranks

#endif

}

