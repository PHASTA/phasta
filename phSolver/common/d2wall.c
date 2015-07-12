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

/*
void
read_d2wall(  int* pid,
               int* numnp,
               double* array1 ) {

//    time_t timenow = time ( &timenow);
    int isize, nitems;
    int iarray[10];

    //MPI_Barrier(MPI_COMM_WORLD);

    /////////////////////////////// Start of writing using new-lib ////////////////////////////

    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

    //  First, count the number of fields to write and store the result in
    //countfieldstowriterestart();

    //  Retrieve and compute the parameters required for SyncIO
    nfiles = outpar.nsynciofiles;
    //    nfields = 1; //outpar.nsynciofieldswriterestart;  // Only the distance to the walls in d2wall
    numparts = workfc.numpe;
    irank = *pid; // workfc.myrank;
    nprocs = workfc.numpe;
    int nppf = numparts/nfiles;
    int GPID;

    // Calculate number of parts each proc deal with and where it start and end ...
    int nppp = numparts/nprocs;// nppp : Number of parts per proc ...
    int startpart = irank * nppp +1;// Part id from which I (myrank) start ...
    int endpart = startpart + nppp - 1;// Part id to which I (myrank) end ...

    int descriptor;
    char filename[255],path[255],fieldtag_s[255];
    bzero((void*)filename,255);
    bzero((void*)fieldtag_s,255);

    sprintf(filename,"d2wall.%d",((int)(irank/(nprocs/nfiles))+1));

    queryphmpiio(filename, &nfields, &nppf);

    initphmpiio(&nfields, &nppf, &nfiles, &f_descriptor, "read");

    if (irank==0) {
      printf("Filename is %s \n",filename);
    }
    openfile(filename, "read", &f_descriptor);

    field_flag=0;

    int i;
    for ( i = 0; i < nppp; i++) { //This loop is useful only if several parts per processor
    // GPID : global part id, corresponds to rank ...
      // e.g : (in this example)
      // proc 0 : 1--4
      // proc 1 : 5--8 ...
      GPID = startpart + i;

      // Write solution field ...
      sprintf(fieldtag_s,"d2wall@%d",GPID);

      nitems = 2;
      readheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, "double", phasta_iotype);
      //iarray[ 0 ] = (*numnp); What we should get from readheader
      //iarray[ 1 ] = 1;

      if (iarray[0] != (*numnp)) {
        printf("ERROR - numnp in d2wall.dat not coherent %d %d/n", iarray[0], *numnp);
      }

      isize = (*numnp);
      readdatablock( &f_descriptor, fieldtag_s, (void*)(array1), &isize, "double", phasta_iotype );

    }
    field_flag++;

    if (field_flag==1){

      closefile(&f_descriptor, "read");

      finalizephmpiio(&f_descriptor);

      if (irank==0) {
        printf("\n");
      }
    }
}
*/


void
read_d2wall(  int* pid,
              int* numnp,
              double* array1,
              int* foundd2wall ) {

//    time_t timenow = time ( &timenow);

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
    int nppf = numparts/nfiles;

    // Calculate number of parts each proc deal with and where it start and end ...
    int nppp = numparts/nprocs;// nppp : Number of parts per proc ...
    assert(nppp==1);
    int startpart = irank * nppp +1;// Part id from which I (myrank) start ...
    int endpart = startpart + nppp - 1;// Part id to which I (myrank) end ...

    int descriptor;
    char filename[255],path[255];
    memset((void*)filename,0,255);
    *foundd2wall = 0;
    ////////////////////////////////////////////////////
    // First we try to read dwal from the restart files.
    ////////////////////////////////////////////////////

    phio_restartname(&(timdat.lstep), filename);
    phio_openfile(filename, "read", &nfiles, &f_descriptor);

    int i;
    for ( i = 0; i < nppp; i++) { //This loop is useful only if several parts per processor
      nitems = 2;
      phio_readheader( &f_descriptor, "dwal", (void*)iarray, &nitems, "double", phasta_iotype);
      //iarray[ 0 ] = (*numnp); What we should get from readheader
      //iarray[ 1 ] = 1;

      if (iarray[0] == (*numnp)) {
        if (irank==0) {
          printf("d2wall field found in %s\n",filename);
        }
        *foundd2wall = 1;
        isize = (*numnp);
        phio_readdatablock( &f_descriptor, "dwal", (void*)(array1), &isize, "double", phasta_iotype );
      }
      else { //d2wall fields was not found in the restart file
        *foundd2wall = 0;
        if (irank==0) {
          printf("d2wall field not found in %s - trying d2wall files now\n",filename);
        }
      }
    }
    closefile(&f_descriptor, "read");
    finalizephmpiio(&f_descriptor);

    ////////////////////////////////////////////////////
    // We try to read dwal from the d2wall files if not found in the restart files
    ////////////////////////////////////////////////////

    int numd2wallfiles;
    if (*foundd2wall == 0) {

      detectd2wallfiles(&numd2wallfiles);

      if (numd2wallfiles == outpar.nsynciofiles ) {
        // Read the d2wall field from the d2wall files
        phio_openfile("d2wall.", "read", &nfiles, &f_descriptor);

        int i;
        for ( i = 0; i < nppp; i++) { //This loop is useful only if several parts per processor
          nitems = 2;
          phio_readheader( &f_descriptor, "d2wall", (void*)iarray, &nitems, "double", phasta_iotype);
          //iarray[ 0 ] = (*numnp); What we should get from readheader
          //iarray[ 1 ] = 1;

          if (iarray[0] == (*numnp)) {
            if (irank==0) {
              printf("d2wall field found in %s\n",filename);
            }
            *foundd2wall = 1;
            isize = (*numnp);
            phio_readdatablock( &f_descriptor, "d2wall", (void*)(array1), &isize, "double", phasta_iotype );
          }
          else {
            *foundd2wall = 0;
              printf("WARNING - numnp not coherent in d2wall files: %d - %d\n",iarray[0],*numnp);
              printf("WARNING - Recomputing the d2wall field for safety\n");
            }
          }

          closefile(&f_descriptor, "read");
          finalizephmpiio(&f_descriptor);
      }
      else if (numd2wallfiles != 0) {
        // The number of d2wall file should be either 0 or outpar.nsynciofiles
        if (irank==0) {
          printf("WARNING - Number of d2wall files not coherent: %d - %d\n",numd2wallfiles,outpar.nsynciofiles);
          printf("WARNING - Recomputing the d2wall field for safety\n");
          *foundd2wall = 0;
        }
      }
    } // end of tentative reading from d2wall files

    if (irank==0) {
      printf("\n");
    }
}
