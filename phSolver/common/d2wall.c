#include <stdlib.h>
#include <FCMangle.h>
#include <new_interface.h>
#include <stdio.h>
#include <string.h> //memset
#include <assert.h>
#include "common_c.h"
#include "phastaIO.h"
#include "phIO.h"
#include "setsyncioparam.h"

void
read_d2wall(  int* pid,
              int* numnp,
              double* array1,
              int* foundd2wall ) {
    int isize, nitems;
    int iarray[10];
    int j;
    for ( j = 0; j < 10; j++) { 
       //Initialize iarray to 0 so that we can assess the result of readheader
       iarray[j] = 0;
    }

    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

    //  Retrieve and compute the parameters required for SyncIO
    nfiles = outpar.nsynciofiles;
    numparts = workfc.numpe;
    irank = *pid; // workfc.myrank;
    nprocs = workfc.numpe;

    // Calculate number of parts each proc deal with and where it start and end ...
    int nppp = numparts/nprocs;// nppp : Number of parts per proc ...
    assert(nppp==1);
    int startpart = irank * nppp +1;// Part id from which I (myrank) start ...
    int endpart = startpart + nppp - 1;// Part id to which I (myrank) end ...

    phio_fp handle;
    char filename[255],path[255];
    memset((void*)filename,0,255);
    *foundd2wall = 0;
    ////////////////////////////////////////////////////
    // First we try to read dwal from the restart files.
    ////////////////////////////////////////////////////

    if ( nfiles == 0 ) {
      sprintf(filename,"restart.%d.", timdat.lstep);
    }
    else {
      sprintf(filename,"restart-dat.%d.", timdat.lstep);
    }
    phio_openfile_read(filename, &nfiles, &handle);

    int i;
    for ( i = 0; i < nppp; i++) { //This loop is useful only if several parts per processor
      nitems = 2;
      phio_readheader(handle, "dwal", (void*)iarray, &nitems, "double", phasta_iotype);

      if (iarray[0] == (*numnp)) {
        if (irank==0) {
          printf("d2wall field found in %s\n",filename);
        }
        *foundd2wall = 1;
        isize = (*numnp);
        phio_readdatablock(handle, "dwal", (void*)(array1), &isize, "double", phasta_iotype );
      }
      else { //d2wall fields was not found in the restart file
        *foundd2wall = 0;
        if (irank==0) {
          printf("d2wall field not found in %s - trying d2wall files now\n",filename);
        }
      }
    }
    phio_closefile_read(handle);

    if (irank==0) {
      printf("\n");
    }
}
