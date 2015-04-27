/* This file provides interface functions for 'partial ' random
   access into the PHASTA input files

   Anil Karanam March 2001 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "phastaIO.h"
#include "rdtsc.h"
#include <FCMangle.h>
#include "new_interface.h"

//MR CHANGE
#include "common_c.h"
//MR CHANGE END

#ifdef intel
#include <winsock2.h>
#else
#include <unistd.h>
#include <strings.h>
#endif

//extern double cpu_speed = 2600000000.0; //for Jaguar XT5
//extern double cpu_speed = 2100000000.0;  //for Jaguar xt4
//extern double cpu_speed =  850000000.0;  //for Intrepid

void igetMinMaxAvg(int *ivalue, double *stats, int *statRanks) {
  int isThisRank;

  double *value = (double*)malloc(sizeof(double));
  *value = 1.0*(*ivalue);

  rgetMinMaxAvg(value,stats,statRanks);

  /* MPI_Allreduce(value,&stats[0],1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  isThisRank=workfc.numpe+1;
  if(*value==stats[0])
    isThisRank=workfc.myrank; 
  MPI_Allreduce(&isThisRank,&statRanks[0],1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

  MPI_Allreduce(value,&stats[1],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  isThisRank=workfc.numpe+1;
  if(*value==stats[1])
    isThisRank=workfc.myrank; 
  MPI_Allreduce(&isThisRank,&statRanks[1],1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

  MPI_Allreduce(value,&stats[2],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  stats[2] /= workfc.numpe; */

  free(value);
}

void rgetMinMaxAvg(double *value, double *stats, int *statRanks) {
  int isThisRank;

  MPI_Allreduce(value,&stats[0],1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  isThisRank=workfc.numpe+1;
  if(*value==stats[0])
    isThisRank=workfc.myrank; 
  MPI_Allreduce(&isThisRank,&statRanks[0],1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

  MPI_Allreduce(value,&stats[1],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  isThisRank=workfc.numpe+1;
  if(*value==stats[1])
    isThisRank=workfc.myrank; 
  MPI_Allreduce(&isThisRank,&statRanks[1],1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

  MPI_Allreduce(value,&stats[2],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  stats[2] /= workfc.numpe;

  double sqValue = (*value)*(*value), sqValueAvg = 0.;
  MPI_Allreduce(&sqValue,&sqValueAvg,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  sqValueAvg /= workfc.numpe;
  // stats[3] = sqValueAvg;

  stats[3] = sqrt(sqValueAvg-stats[2]*stats[2]);
}

void print_mesh_stats(void) {
  int statRanks[2];
  double iStats[4], rStats[4];

  igetMinMaxAvg(&conpar.nshg,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("nshg    : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&conpar.numel,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("numel   : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&conpar.numelb,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("numelb  : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&conpar.nnz_tot,iStats,statRanks);
  if(workfc.myrank==workfc.master) {
    printf("nnz_tot : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
    printf("\n");
  }
}

void print_mpi_stats(void) {
  int statRanks[2];
  double iStats[4], rStats[4];

// NS equations
  igetMinMaxAvg(&mpistats.iISend,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iISend : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&mpistats.iIRecv,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iIRecv : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&mpistats.iWaitAll,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iWtAll : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&mpistats.iAllR,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iAllR  : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);

  rgetMinMaxAvg(&mpistats.rISend,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rISend : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rIRecv,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rIRecv : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rWaitAll,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rWtAll : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rCommu,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rCommu : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rAllR,rStats,statRanks);
  if(workfc.myrank==workfc.master) {
    printf("rAllR  : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
    printf("\n");
  }
// Scalars
  igetMinMaxAvg(&mpistats.iISendScal,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iISendScal : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&mpistats.iIRecvScal,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iIRecvScal : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&mpistats.iWaitAllScal,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iWtAllScal : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);
  igetMinMaxAvg(&mpistats.iAllRScal,iStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("iAllRScal : min [%d,%d], max[%d,%d] and avg[.,%d] (rms=%d)\n",statRanks[0],(int)iStats[0],statRanks[1],(int)iStats[1],(int)iStats[2],(int)iStats[3]);

  rgetMinMaxAvg(&mpistats.rISendScal,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rISendScal : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rIRecvScal,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rIRecvScal : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rWaitAllScal,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rWtAllScal : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rCommuScal,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rCommuScal : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
  rgetMinMaxAvg(&mpistats.rAllRScal,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("rAllRScal  : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);


}

//void print_system_stats(double tcorecp[2]) {
void print_system_stats(double *tcorecp, double *tcorecpscal) {
  int statRanks[2];
  double iStats[4], rStats[4];
  double syst_assembly, syst_solve;

// NS equations
  syst_assembly = tcorecp[0];
  syst_solve = tcorecp[1];

  rgetMinMaxAvg(&syst_assembly,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("Elm. form. : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);

  rgetMinMaxAvg(&syst_solve,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("Lin. alg. sol : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);

// Scalars
  syst_assembly = tcorecpscal[0];
  syst_solve = tcorecpscal[1];

  rgetMinMaxAvg(&syst_assembly,rStats,statRanks);
  if(workfc.myrank==workfc.master)
    printf("Elm. form. Scal. : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);

  rgetMinMaxAvg(&syst_solve,rStats,statRanks);
  if(workfc.myrank==workfc.master) {
    printf("Lin. alg. sol Scal. : min [%d,%2.5f], max[%d,%2.5f] and avg[.,%2.5f] (rms=%2.5f)\n",statRanks[0],rStats[0],statRanks[1],rStats[1],rStats[2],rStats[3]);
    printf("\n");
  }
  //printf("rank %d - syst_assembly %f - syst_solve %f\n",workfc.myrank,syst_assembly,syst_solve);
}



void countfieldstowriterestart()
{
  int nfields;

//     printf("TEST: %d %d %d %d %d\n",timdat.istep,timdat.itseq,inpdat.nstep[0],inpdat.nstep[1],timdat.lstep);

  nfields = 2; //solution, time derivatives

  if(outpar.ivort == 1){
    nfields++; //vorticity
  }

  if(abs(turbvar.itwmod) != 1 && outpar.iowflux == 1) { 
    nfields++; //instantaneous wss in bflux.f
  }

//   if(ideformwall.eq.1) not handled yet

  if(timdat.istep == inpdat.nstep[timdat.itseq-1]){ //Last time step of the computation

    //projection vectors and pressure projection vectors (call saveLesRestart in itrdrv)
    if(matdat.matflg[0][0]==-1) nfields = nfields +2;

    //if Print Error Indicators = true (call write_error in itrdrv)
    if(turbvar.ierrcalc == 1){
      nfields++;
    }

    //if Print ybar = True (call write_field(myrank,'a','ybar',4,... in itrdrv)
    if(outpar.ioybar == 1){
      nfields++;  //ybar

      //phase average fields
      if(outpar.nphasesincycle >0) {
        nfields = nfields + outpar.nphasesincycle;
      }

      if(abs(turbvar.itwmod) != 1 && outpar.iowflux == 1) { 
        nfields++; //wssbar
      }

    }

    if(turbvari.irans < 0) {
      nfields++; //dwal
    }

  }

  outpar.nsynciofieldswriterestart = nfields;

  if(workfc.myrank == 0) {
    printf("Number of fields to write in restart files: %d\n", nfields);
  }
}


void
Write_Restart(  int* pid,
                int* stepno,
                int* nshg,
                int* numVars,
                double* array1,
                double* array2 ) {

    char fname[255];
    char rfile[60];
    char existingfile[30], linkfile[30];
    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
//    time_t timenow = time ( &timenow);
    double version=0.0;
    int isize, nitems;
    int iarray[10];

    /*sprintf(rfile,"restart.%d.%d",*stepno,*pid+1);
    openfile_(rfile,"write", &irstou);

    // writing the top ascii header for the restart file

    writestring_( &irstou,"# PHASTA Input File Version 2.0\n");
    writestring_( &irstou,
                  "# format \"keyphrase : sizeofnextblock usual headers\"\n");

    bzero( (void*)fname, 255 );
    sprintf(fname,"# Output generated by phasta version (NOT YET CURRENT): %lf \n", version);
    writestring_( &irstou, fname );

    bzero( (void*)fname, 255 );
    gethostname(fname,255);
    writestring_( &irstou,"# This result was produced on: ");
    writestring_( &irstou, fname );
    writestring_( &irstou,"\n");

    bzero( (void*)fname, 255 );
    sprintf(fname,"# %s\n", ctime( &timenow ));
    writestring_( &irstou, fname );

    isize = 1;
    nitems = 1;
    iarray[ 0 ] = 1;
    writeheader_( &irstou, "byteorder magic number ",
                  (void*)iarray, &nitems, &isize, "integer", phasta_iotype );

    nitems = 1;
    writedatablock_( &irstou, "byteorder magic number ",
                     (void*)mptr, &nitems, "integer", phasta_iotype );


    bzero( (void*)fname, 255 );
    sprintf(fname,"number of modes : < 0 > %d\n", *nshg);
    writestring_( &irstou, fname );

    bzero( (void*)fname, 255 );
    sprintf(fname,"number of variables : < 0 > %d\n", *numVars);
    writestring_( &irstou, fname );


    isize = (*nshg)*(*numVars);
    nitems = 3;
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numVars);
    iarray[ 2 ] = (*stepno);
    writeheader_( &irstou, "solution ",
                  (void*)iarray, &nitems, &isize, "double", phasta_iotype );


    nitems = (*nshg)*(*numVars);
    writedatablock_( &irstou, "solution ",
                     (void*)(array1), &nitems, "double", phasta_iotype );



    nitems = 3;
    writeheader_( &irstou, "time derivative of solution ",
                  (void*)iarray, &nitems, &isize, "double", phasta_iotype );


    nitems = (*nshg)*(*numVars);
    writedatablock_( &irstou, "time derivative of solution ",
                     (void*)(array2), &nitems, "double", phasta_iotype );


    closefile_( &irstou, "write" );
    */
    //MPI_Barrier(MPI_COMM_WORLD);

    /////////////////////////////// Start of writing using new-lib ////////////////////////////

//MR CHANGE
    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

    //  First, count the number of fields to write and store the result in
    countfieldstowriterestart();

    //  Retrieve and compute the parameters required for SyncIO
    nfiles = outpar.nsynciofiles;
    nfields = outpar.nsynciofieldswriterestart;
    numparts = workfc.numpe;
    irank = *pid; //workfc.myrank;
    nprocs = workfc.numpe;
//MR CHANGE END
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

    sprintf(filename,"restart-dat.%d.%d",*stepno,((int)(irank/(nprocs/nfiles))+1));

//    unsigned long long timer_start;
//    unsigned long long timer_end;
//    double time_span;

//MR CHANGE
//  Measure the time - Start the timer
//    MPI_Barrier(MPI_COMM_WORLD);
//    timer_start = rdtsc();
//MR CHANGE END

    initphmpiio(&nfields, &nppf, &nfiles, &f_descriptor, "write");

//MR CHANGE
//  Measure the time - End of timer
//    timer_end = rdtsc();
//    time_span=(double)((timer_end-timer_start)/cpu_speed);
    if (*pid==0) {
//      printf("\n*****************************\n");
//      printf("Time: 'initphmpiio' of %s with %d fields and %d files is:    %f s\n",filename,nfields,nfiles,time_span);
      printf("Filename is %s \n",filename);
    }
//MR CHANGE END


//MR CHANGE
//  Measure the time - Start the timer
//    MPI_Barrier(MPI_COMM_WORLD);
//    timer_start = rdtsc();
//MR CHANGE END

    openfile(filename, "write", &f_descriptor);

//MR CHANGE
//  Measure the time - End of timer
//    MPI_Barrier(MPI_COMM_WORLD);
//    timer_end = rdtsc();
//    time_span=(double)((timer_end-timer_start)/cpu_speed);
//    if (*pid==0) {
//      printf("Time: 'openfile' of %s with %d fields and %d files is:    %f s\n",filename,nfields,nfiles,time_span);
//      printf("*****************************\n");
//    }
//MR CHANGE END

    field_flag=0;

     int i;
     for ( i = 0; i < nppp; i++) { //This loop is useful only if several parts per processor
     // GPID : global part id, corresponds to rank ...
        // e.g : (in this example)
        // proc 0 : 1--4
        // proc 1 : 5--8 ...
        GPID = startpart + i;

        // Write solution field ...
        sprintf(fieldtag_s,"solution@%d",GPID);

        isize = (*nshg)*(*numVars);
        nitems = 3;
        iarray[ 0 ] = (*nshg);
        iarray[ 1 ] = (*numVars);
        iarray[ 2 ] = (*stepno);

//MR CHANGE
//  Measure the time - Start the timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_start = rdtsc();
//MR CHANGE END

        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);

//MR CHANGE
//  Measure the time - End of timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_end = rdtsc();
//        time_span=(double)((timer_end-timer_start)/cpu_speed);
//        if (*pid==0) {
//          printf("\n*****************************\n");
//          printf("Time: header for 'Solution':    %f s\n",time_span);
//        }
//MR CHANGE END

        nitems = (*nshg)*(*numVars);

//MR CHANGE
//  Measure the time - Start the timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_start = rdtsc();
//MR CHANGE END

        writedatablock( &f_descriptor, fieldtag_s, (void*)(array1), &isize, "double", phasta_iotype );

//MR CHANGE
//  Measure the time - End of timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_end = rdtsc();
//        time_span=(double)((timer_end-timer_start)/cpu_speed);

//        int isizemin,isizemax,isizetot;
//        double sizemin,sizemax,sizeavg,sizetot,rate;

//        MPI_Allreduce(&isize,&isizemin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
//        MPI_Allreduce(&isize,&isizemax,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
//        MPI_Allreduce(&isize,&isizetot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

//        sizemin=(double)(8.0*isizemin/1024.0/1024.0);
//        sizemax=(double)(8.0*isizemax/1024.0/1024.0);
//        sizetot=(double)(8.0*isizetot/1024.0/1024.0);
//        sizeavg=(double)(1.0*sizetot/workfc.numpe);
//        rate=(double)(1.0*sizetot/time_span);

//        if (*pid==0) {
//          printf("Time: block for 'Solution':    %f s\n",time_span);
//          printf("Time: block:   Min= %f MB; Max= %f MB; Avg= %f MB; Tot= %f MB; Rate= %f MB/s; \n",sizemin,sizemax,sizeavg,sizetot,rate);

//        }
//MR CHANGE END

    }
    field_flag++;

    for ( i = 0; i < nppp; i++) {

        // GPID : global part id, corresponds to rank ...
        // e.g : (in this example)
        // proc 0 : 1--4
        // proc 1 : 5--8 ...
        GPID = startpart + i;

        // Write solution field ...
        sprintf(fieldtag_s,"time derivative of solution@%d",GPID);

        isize = (*nshg)*(*numVars);
        nitems = 3;
        iarray[ 0 ] = (*nshg);
        iarray[ 1 ] = (*numVars);
        iarray[ 2 ] = (*stepno);

//MR CHANGE
//  Measure the time - Start the timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_start = rdtsc();
//MR CHANGE END

        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);

//MR CHANGE
//  Measure the time - End of timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_end = rdtsc();
//        time_span=(double)((timer_end-timer_start)/cpu_speed);
//        if (*pid==0) {
//          printf("Time: header for 'Time Derivative of solution':    %f s\n",time_span);
//        }
//MR CHANGE END

        nitems = (*nshg)*(*numVars);

//MR CHANGE
//  Measure the time - Start the timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_start = rdtsc();
//MR CHANGE END

        writedatablock( &f_descriptor, fieldtag_s, (void*)(array2), &isize, "double", phasta_iotype );

//MR CHANGE
//  Measure the time - End of timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_end = rdtsc();
//        time_span=(double)((timer_end-timer_start)/cpu_speed);

//        int isizemin,isizemax,isizetot;
//        double sizemin,sizemax,sizeavg,sizetot,rate;

//        MPI_Allreduce(&isize,&isizemin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
//        MPI_Allreduce(&isize,&isizemax,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
//        MPI_Allreduce(&isize,&isizetot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

//        sizemin=(double)(8.0*isizemin/1024.0/1024.0);
//        sizemax=(double)(8.0*isizemax/1024.0/1024.0);
//        sizetot=(double)(8.0*isizetot/1024.0/1024.0);
//        sizeavg=sizetot/workfc.numpe;
//        rate=sizetot/time_span;

//        if (*pid==0) {
//          printf("Time: block for 'Time Derivative of Solution':    %f s\n",time_span);
//          printf("Time: block:   Min= %f MB; Max= %f MB; Avg= %f MB; Tot= %f MB; Rate= %f MB/s; \n",sizemin,sizemax,sizeavg,sizetot,rate);
//          printf("*****************************\n");

//        }
//MR CHANGE END

    }
    field_flag++;

    if (field_flag==nfields){

//MR CHANGE
//  Measure the time - Start the timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_start = rdtsc();
//MR CHANGE END

      closefile(&f_descriptor, "write");

//MR CHANGE
//    Measure the time - End of timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_end = rdtsc();
//      time_span=(double)((timer_end-timer_start)/cpu_speed);
//      if (*pid==0) {
//        printf("\n*****************************\n");
//        printf("Time: 'closefile' is:    %f s\n",time_span);
//      }
//MR CHANGE END

//MR CHANGE
//  Measure the time - Start the timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_start = rdtsc();
//MR CHANGE END

      finalizephmpiio(&f_descriptor);

//MR CHANGE
//    Measure the time - End of timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_end = rdtsc();
//      time_span=(double)((timer_end-timer_start)/cpu_speed);
      if (*pid==0) {
//        printf("Time: 'finalizephmpiio' is:    %f s\n",time_span);
//        printf("Last field %d '%s' finished! \n",nfields, fieldtag_s);
        printf("\n");
//        printf("*****************************\n");
      }
    }
//MR CHANGE END



    ///////////////////////////////////////////////////////////////////////////////////////////

    /* create a soft link of the restart we just wrote to restart.latest
     this is the file the next run will always try to start from */

/*    sprintf( linkfile, "restart.latest.%d", *pid+1 );
    unlink( linkfile );
    sprintf( existingfile, "restart.%d.%d", *stepno, *pid+1 );
    link( existingfile, linkfile );
*/
}

void
Write_Error(  int* pid,
              int* stepno,
              int* nshg,
              int* numVars,
              double* array1 ) {


    char fname[255];
    char rfile[60];
    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    //printf("Time is commented\n");
    //time_t timenow = time ( &timenow);
    //printf("Yes\n");
    double version=0.0;
    int isize, nitems;
    int iarray[10];

    /*sprintf(rfile,"restart.%d.%d",*stepno,*pid+1);
    openfile_(rfile,"append", &irstou);

    isize = (*nshg)*(*numVars);
    nitems = 3;
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numVars);
    iarray[ 2 ] = (*stepno);
    writeheader_( &irstou, "errors", (void*)iarray, &nitems, &isize, "double", phasta_iotype );


    nitems = (*nshg)*(*numVars);
    writedatablock_( &irstou, "errors ", (void*)(array1), &nitems, "double", phasta_iotype );

    closefile_( &irstou, "append" );*/

    /////////////////////////////// Start of writing using new-lib ////////////////////////////

    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

//    unsigned long long timer_start;
//    unsigned long long timer_end;
//    double time_span;

    nfiles = outpar.nsynciofiles;
    nfields = outpar.nsynciofieldswriterestart;
    numparts = workfc.numpe;
    irank = *pid; //workfc.myrank;
    nprocs = workfc.numpe;

    int nppf = numparts/nfiles;
    int GPID;

    // Calculate number of parts each  proc deal with and where it start and end ...
    int nppp = numparts/nprocs;// nppp : Number of parts per proc ...
    int startpart = irank * nppp +1;// Part id from which I (myrank) start ...
    int endpart = startpart + nppp - 1;// Part id to which I (myrank) end ...

    field_flag++;

    char fieldtag[255];

    int i;
    for ( i = 0; i < nppp; i++  ) {
        GPID = startpart + i;
        sprintf(fieldtag,"errors@%d",GPID);

        if(*pid==0) {
//          printf("\n*****************************\n");
          printf("\n");
          printf("The %d/%d th field to be written is '%s'\n",field_flag,nfields,fieldtag);
        }

        isize = (*nshg)*(*numVars);
        nitems = 3;
        iarray[ 0 ] = (*nshg);
        iarray[ 1 ] = (*numVars);
        iarray[ 2 ] = (*stepno);

//MR CHANGE
//  Measure the time - Start the timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_start = rdtsc();
//MR CHANGE END

        writeheader( &f_descriptor, fieldtag, (void*)iarray, &nitems, &isize, "double", phasta_iotype);

//MR CHANGE
//  Measure the time - End of timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_end = rdtsc();
//        time_span=(double)((timer_end-timer_start)/cpu_speed);
//        if (*pid==0) {
//          printf("Time: header for 'error':    %f s\n",time_span);
//        }
//MR CHANGE END

//MR CHANGE
//  Measure the time - Start the timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_start = rdtsc();
//MR CHANGE END

        writedatablock( &f_descriptor, fieldtag, (void*)array1, &isize, "double", phasta_iotype );

//MR CHANGE
//  Measure the time - End of timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_end = rdtsc();
//        time_span=(double)((timer_end-timer_start)/cpu_speed);

//        int isizemin,isizemax,isizetot;
//        double sizemin,sizemax,sizeavg,sizetot,rate;

//        MPI_Allreduce(&isize,&isizemin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
//        MPI_Allreduce(&isize,&isizemax,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
//        MPI_Allreduce(&isize,&isizetot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

//        sizemin=(double)(8.0*isizemin/1024.0/1024.0);
//        sizemax=(double)(8.0*isizemax/1024.0/1024.0);
//        sizetot=(double)(8.0*isizetot/1024.0/1024.0);
//        sizeavg=sizetot/workfc.numpe;
//        rate=sizetot/time_span;

//        if (*pid==0) {
//          printf("Time: block for 'error':    %f s\n",time_span);
//          printf("Time: block:   Min= %f MB; Max= %f MB; Avg= %f MB; Tot= %f MB; Rate= %f MB/s; \n",sizemin,sizemax,sizeavg,sizetot,rate);
//          printf("*****************************\n");
//        }
//MR CHANGE END

    }

//     MPI_Barrier(MPI_COMM_WORLD);
//     timer_end = rdtsc();
//     time_span=(double)(timer_end-timer_start)/cpu_speed;

//     if (*pid==0) {
//         printf("Field 'error' written in:     %f s\n",time_span);
//         printf("Write field '%s' finished! \n",fieldtag);
//     }

    if (field_flag==nfields){

//MR CHANGE
//  Measure the time - Start the timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_start = rdtsc();
//MR CHANGE END

      closefile(&f_descriptor, "write");

//MR CHANGE
//    Measure the time - End of timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_end = rdtsc();
//      time_span=(double)((timer_end-timer_start)/cpu_speed);
//      if (*pid==0) {
//        printf("\n*****************************\n");
//        printf("Time: 'closefile' is:    %f s\n",time_span);
//      }
//MR CHANGE END

//MR CHANGE
//  Measure the time - Start the timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_start = rdtsc();
//MR CHANGE END

      finalizephmpiio(&f_descriptor);

//MR CHANGE
//    Measure the time - End of timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_end = rdtsc();
//      time_span=(double)((timer_end-timer_start)/cpu_speed);
      if (*pid==0) {
//        printf("Time: 'finalizephmpiio' is:    %f s\n",time_span);
        printf("Last field %d '%s' finished! \n",nfields, fieldtag);
        printf("\n");
//        printf("*****************************\n");
      }
    }
//MR CHANGE END

    ///////////////////////////////////////////////////////////////////////////////////////////


}


void
Write_Displ(  int* pid,
              int* stepno,
              int* nshg,
              int* numVars,
              double* array1 ) {


    char fname[255];
    char rfile[60];
    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    time_t timenow = time ( &timenow);
    double version=0.0;
    int isize, nitems;
    int iarray[10];

    sprintf(rfile,"restart.%d.%d",*stepno,*pid+1);
    openfile(rfile,"append", &irstou);

    isize = (*nshg)*(*numVars);
    nitems = 3;
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numVars);
    iarray[ 2 ] = (*stepno);
    writeheader( &irstou, "displacement", (void*)iarray, &nitems, &isize, "double", phasta_iotype );

    nitems = (*nshg)*(*numVars);
    writedatablock( &irstou, "displacement", (void*)(array1), &nitems, "double", phasta_iotype );

    closefile( &irstou, "append" );
}

void
Write_Debug(	int *pid,			//ramk of mpi process
				char* filetag,		//Name of the file to write. This will be appended with %d.%d to track step and processor. WARNING: this must be null terminated. 
				char* fieldtag,		//name of field to write data into. WARNING: this must be null terminated. When calling from Fortran, use //char(0) 
				void* array,		//data array to write to file
				char* arraytype,	//arraytype[0] == 'd': double, arraytype[0] == 'i': int
				int* nshg,			//total number of shape functions
				int* numvars,		//number of variables per state (i.e. size of the matrix in dimension 2)
				int* stepno) {		//step number
//Wites the field in array to a set of files of the form [filetag].[time step].[file number]. 
//
//WARNING:  Both filetag and fieldtag must be null terminated. Fortran does not normally do this when storing strings. 
//          To get around this, the string MUST be concatenated with char(0) before or in the function call, e.g. 
//          Write_FieldDebug(myrank, 'foo'//char(0) ...)
//
  
//	int irstou;
//	int magic_number = 362436;
//	int* mptr = &magic_number;
//	double version=0.0;
//	char fmode[10];
//	strcpy(fmode, "write");
//    memset((void*)filename,   0, 255);	//unnecessary with if arrays are zeroed at allocation using 
//    memset((void*)fieldtag_s, 0, 255);	// name[size] = {0}

    int nfields = 1; // outpar.nsynciofieldswriterestart;
    int numparts = workfc.numpe;
    int irank = *pid;
    int nprocs = workfc.numpe;
#if(WRITE_DEBUG_OUTPUT_TYPE == POSIX)		//set in new_interface.h
	int nfiles = nprocs;
#elif(WRITE_DEBUG_OUTPUT_TYPE == SYNCIO)
    int nfiles = outpar.nsynciofiles;
#endif

    // Calculate number of parts each  proc deal with and where it start and end ...
    int nppp = numparts/nprocs;// nppp : Number of parts per proc ...
    int nppf = numparts/nfiles;
    int GPID;
    int startpart = irank * nppp +1;// Part id from which I (myrank) start ...
//	int endpart = startpart + nppp - 1;// Part id to which I (myrank) end ...

	int magic_number = 362436;
	int* mptr = &magic_number;
	int isize, nitems;
	int iarray[10];

    char filename[255] = {0}, fieldtag_s[255] = {0};
	sprintf(filename, "%s.%d.%d",filetag, *stepno, (int)(irank/(nprocs/nfiles))+1 );

    char datatype[10];
    if(arraytype[0] == 'i')		
		strcpy(datatype,"int");
    else 								// default is double
		strcpy(datatype,"double");

//    initphmpiio(&nfields, &nppf, &nfiles, &f_descriptor, "write");
	if (irank == 0) 
		printf("Filename is %s \n",filename);
    openfile(filename, "write", &f_descriptor);

    /////////////////////////////// Start of writing using new-lib ////////////////////////////

    if(irank == 0)
		printf("\nStarting to write the field '%s'\n", fieldtag);

	int i;
    for (i = 0; i < nppp; i++  ) {			//loop over all prarts per processor
        GPID = startpart + i;

        // Write solution field ...
#if(WRITE_DEBUG_OUTPUT_TYPE == POSIX)
		//write the byte order magic number before the main field
		isize = 1;
		nitems = 1;
		iarray[ 0 ] = 1;
		writeheader(    &f_descriptor, "byteorder magic number ", (void*)iarray, &nitems, &isize, "integer", phasta_iotype );
		writedatablock( &f_descriptor, "byteorder magic number ", (void*)mptr, &nitems, "integer", phasta_iotype );

		sprintf(fieldtag_s, "%s", fieldtag);
#elif(WRITE_DEBUG_OUTPUT_TYPE == SYNCIO)
        sprintf(fieldtag_s, "%s@%d",fieldtag, GPID);
#endif

        isize = (*nshg)*(*numvars);
        nitems = 3;
        iarray[ 0 ] = (*nshg);
        iarray[ 1 ] = (*numvars);
        iarray[ 2 ] = (*stepno);

        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, datatype, phasta_iotype);
//		nitems = (*nshg)*(*numvars);
        writedatablock( &f_descriptor, fieldtag_s, array, &isize, datatype, phasta_iotype );
    }

    closefile(&f_descriptor, "write");
//	finalizephmpiio(&f_descriptor);

	if (irank == 0) 
		printf("Finish writting the field '%s'! \n\n", fieldtag);
}
//N Mati change end

void
Write_Field(  int *pid,
              char* filemode,
              char* fieldtag,
              int* tagsize,
              void* array,
              char* arraytype,
              int* nshg,
              int* numvars,
              int* stepno) {

    //printf("Rank is %d, field is %s, tagsize is %d, nshg is %d, numvars is %d\n",*pid,fieldtag,*tagsize,*nshg,*numvars);

//     char rfile[32];
    // assuming restart.sn.(pid+1)
//     sprintf(rfile,"restart.%d.%d",*stepno,*pid+1);

    char *fieldlabel = (char *)malloc((*tagsize+1)*sizeof(char));
    strncpy(fieldlabel, fieldtag, *tagsize);
    fieldlabel[*tagsize] = '\0';

    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    double version=0.0;
    int isize, nitems;
    int iarray[10];

    char fmode[10];
    if(!strncmp(filemode,"w",1))
      strcpy(fmode,"write");
    else // default is append
      strcpy(fmode,"append");

    char datatype[10];
    if(!strncmp(arraytype,"i",1))
      strcpy(datatype,"int");
    else // default is double
      strcpy(datatype,"double");

/*     openfile_(rfile, fmode, &irstou);

     nitems = 3; // assuming field will write 3 items in iarray
     iarray[ 0 ] = (*nshg);
     iarray[ 1 ] = (*numvars);
     iarray[ 2 ] = (*stepno);

     isize = (*nshg)*(*numvars);
     writeheader_( &irstou, fieldlabel, (void*)iarray, &nitems, &isize, datatype, phasta_iotype );

     nitems = (*nshg)*(*numvars);
     writedatablock_( &irstou, fieldlabel, array, &nitems, datatype, phasta_iotype );
     closefile_( &irstou, fmode);
*/
    /////////////////////////////// Start of writing using new-lib ////////////////////////////

    int nfiles = outpar.nsynciofiles;
    int nfields = outpar.nsynciofieldswriterestart;
    int numparts = workfc.numpe;
    int irank = *pid;
    int nprocs = workfc.numpe;

//    unsigned long long timer_start;
//    unsigned long long timer_end;
//    double time_span;


    int nppf = numparts/nfiles;
    int GPID;

    // Calculate number of parts each  proc deal with and where it start and end ...
    int nppp = numparts/nprocs;// nppp : Number of parts per proc ...
    int startpart = irank * nppp +1;// Part id from which I (myrank) start ...
    int endpart = startpart + nppp - 1;// Part id to which I (myrank) end ...

    char filename[255],path[255],fieldtag_s[255];
    bzero((void*)filename,255);
    bzero((void*)fieldtag_s,255);

    strncpy(fieldlabel, fieldtag, *tagsize);

    field_flag++;
    if(*pid==0) {
//      printf("\n*****************************\n");
      printf("\n");
      printf("The %d/%d the field to be written is '%s'\n",field_flag,nfields,fieldtag);
    }

    sprintf(filename,"restart-dat.%d.%d",*stepno,((int)(irank/(nprocs/nfiles))+1));

//     MPI_Barrier(MPI_COMM_WORLD);
//     timer_start = rdtsc();

    int i;
    for ( i = 0; i < nppp; i++  ) {
        GPID = startpart + i;

        // Write solution field ...
        sprintf(fieldtag_s,"%s@%d",fieldlabel,GPID);

        isize = (*nshg)*(*numvars);
        nitems = 3;
        iarray[ 0 ] = (*nshg);
        iarray[ 1 ] = (*numvars);
        iarray[ 2 ] = (*stepno);

//MR CHANGE
//  Measure the time - Start the timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_start = rdtsc();
//MR CHANGE END

        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, datatype, phasta_iotype);

//MR CHANGE
//  Measure the time - End of timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_end = rdtsc();
//        time_span=(double)((timer_end-timer_start)/cpu_speed);
//        if (*pid==0) {
//          printf("Time: header for '%s':    %f s\n",fieldtag_s,time_span);
//        }
//MR CHANGE END

        nitems = (*nshg)*(*numvars);

//MR CHANGE
//  Measure the time - Start the timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_start = rdtsc();
//MR CHANGE END

        writedatablock( &f_descriptor, fieldtag_s, array, &isize, datatype, phasta_iotype );

//MR CHANGE
//  Measure the time - End of timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_end = rdtsc();
//        time_span=(double)((timer_end-timer_start)/cpu_speed);

//        int isizemin,isizemax,isizetot;
//        double sizemin,sizemax,sizeavg,sizetot,rate;

//        MPI_Allreduce(&isize,&isizemin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
//        MPI_Allreduce(&isize,&isizemax,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
//        MPI_Allreduce(&isize,&isizetot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

//        sizemin=(double)(8.0*isizemin/1024.0/1024.0);
//        sizemax=(double)(8.0*isizemax/1024.0/1024.0);
//        sizetot=(double)(8.0*isizetot/1024.0/1024.0);
//        sizeavg=sizetot/workfc.numpe;
//        rate=sizetot/time_span;

//        if (*pid==0) {
//          printf("Time: block for '%s':    %f s\n",fieldtag_s,time_span);
//          printf("Time: block:   Min= %f MB; Max= %f MB; Avg= %f MB; Tot= %f MB; Rate= %f MB/s; \n",sizemin,sizemax,sizeavg,sizetot,rate);
//          printf("*****************************\n");
//        }
//MR CHANGE END

    }

//     MPI_Barrier(MPI_COMM_WORLD);
//     timer_end = rdtsc();
//     time_span=(double)(timer_end-timer_start)/cpu_speed;

//     if (*pid==0) {
//         printf("Field '%s' written in:     %f s\n",fieldtag,time_span);
//         printf("Write field '%s' finished! \n",fieldtag_s);
//     }

//     if (field_flag==nfields){
//       closefile_(&f_descriptor, "write");
//       finalizephmpiio_(&f_descriptor);
//       if(*pid==0) {
//         printf("Last field %d '%s' finished! \n",nfields, fieldtag_s);
//         printf("\n*****************************\n");
//       }
//     }

    if (field_flag==nfields){

//MR CHANGE
//  Measure the time - Start the timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_start = rdtsc();
//MR CHANGE END

      closefile(&f_descriptor, "write");

//MR CHANGE
//    Measure the time - End of timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_end = rdtsc();
//      time_span=(double)((timer_end-timer_start)/cpu_speed);
//      if (*pid==0) {
//        printf("\n*****************************\n");
//        printf("Time: 'closefile' is:    %f s\n",time_span);
//      }
//MR CHANGE END

//MR CHANGE
//  Measure the time - Start the timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_start = rdtsc();
//MR CHANGE END

      finalizephmpiio(&f_descriptor);

//MR CHANGE
//    Measure the time - End of timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_end = rdtsc();
//      time_span=(double)((timer_end-timer_start)/cpu_speed);
      if (*pid==0) {
//        printf("Time: 'finalizephmpiio' is:    %f s\n",time_span);
        printf("Last field %d '%s' finished! \n",nfields, fieldtag);
        printf("\n");
//        printf("*****************************\n");
      }
    }
//MR CHANGE END

    ///////////////////////////////////////////////////////////////////////////////////////////

    free(fieldlabel);
}

//MR CHANGE
void
Write_PhAvg(  int* pid,
              char* filemode,
              char* fieldtag,
              int* tagsize,
              int* iphase,
              void* array,
              char* arraytype,
              int* nshg,
              int* numvars,
              int* stepno) {

    char rfile[32];
    // assuming restart_phase_avg_<sn>.<iphase>.<pid+1>
    sprintf(rfile,"restart_phase_avg_%d.%d.%d",*stepno,*iphase,*pid+1);

    char *fieldlabel = (char *)malloc((*tagsize+1)*sizeof(char));
    strncpy(fieldlabel, fieldtag, *tagsize);
    fieldlabel[*tagsize] = '\0';

    int irstou;
    int isize, nitems;
    int iarray[10];

    char fmode[10];
    if(!strncmp(filemode,"w",1))
      strcpy(fmode,"write");
    else // default is append
      strcpy(fmode,"append");

    char datatype[10];
    if(!strncmp(arraytype,"i",1))
      strcpy(datatype,"int");
    else // default is double
      strcpy(datatype,"double");

    openfile(rfile, fmode, &irstou);

    if(!strcmp(fmode,"write")) {
      // may be create a routine for 'top' portion under write mode
      int magic_number = 362436;
      int* mptr = &magic_number;
      time_t timenow = time ( &timenow);
      double version=0.0;

      /* writing the top ascii header for the restart file */

      writestring( &irstou,"# PHASTA Input File Version 2.0\n");
      writestring( &irstou,
                    "# format \"keyphrase : sizeofnextblock usual headers\"\n");

      char fname[255];
      bzero( (void*)fname, 255 );
      sprintf(fname,"# Output generated by phasta version (NOT YET CURRENT): %lf \n", version);
      writestring( &irstou, fname );

      bzero( (void*)fname, 255 );
      gethostname(fname,255);
      writestring( &irstou,"# This result was produced on: ");
      writestring( &irstou, fname );
      writestring( &irstou,"\n");

      bzero( (void*)fname, 255 );
      sprintf(fname,"# %s\n", ctime( &timenow ));
      writestring( &irstou, fname );

      isize = 1;
      nitems = 1;
      iarray[ 0 ] = 1;
      writeheader( &irstou, "byteorder magic number ",
                    (void*)iarray, &nitems, &isize, "integer", phasta_iotype );
      nitems = 1;
      writedatablock( &irstou, "byteorder magic number ",
                       (void*)mptr, &nitems, "integer", phasta_iotype );
    }

    nitems = 3; // assuming field will write 3 items in iarray
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numvars);
    iarray[ 2 ] = (*stepno);

    isize = (*nshg)*(*numvars);
    writeheader( &irstou, fieldlabel, (void*)iarray, &nitems, &isize, datatype, phasta_iotype );

    nitems = (*nshg)*(*numvars);
    writedatablock( &irstou, fieldlabel, array, &nitems, datatype, phasta_iotype );

    closefile( &irstou, fmode);

    free(fieldlabel);
}

//MR CHANGE
void
Write_PhAvg2( int* pid,
              char* filemode,
              char* fieldtag,
              int* tagsize,
              int* iphase,
              int* nphasesincycle,
              void* array,
              char* arraytype,
              int* nshg,
              int* numvars,
              int* stepno) {

//     char rfile[32];
    // assuming restart.sn.(pid+1)
//     sprintf(rfile,"restart.%d.%d",*stepno,*pid+1);

    int addtagsize; // phase number is added to the name of the field
    if(*iphase<10)
      addtagsize=1;
    else if(*iphase<100)
      addtagsize=2;
    else if(*iphase<1000)
      addtagsize=3;

    int tagsize2;
    tagsize2=*tagsize+addtagsize;

//     char *fieldlabel = (char *)malloc((*tagsize+1)*sizeof(char));
//     strncpy(fieldlabel, fieldtag, *tagsize);
//     fieldlabel[*tagsize] = '\0';

    char *fieldlabel = (char *)malloc((tagsize2+1)*sizeof(char));
    strncpy(fieldlabel, fieldtag, *tagsize);
    fieldlabel[tagsize2] = '\0';

    char straddtagsize[10];
    sprintf(straddtagsize,"%d",*iphase);

    if(*iphase<10) {
      fieldlabel[tagsize2-1]=straddtagsize[0];
    }
    else if(*iphase<100) {
      fieldlabel[tagsize2-2]=straddtagsize[0];
      fieldlabel[tagsize2-1]=straddtagsize[1];
    }
    else if(*iphase<1000) {
      fieldlabel[tagsize2-3]=straddtagsize[0];
      fieldlabel[tagsize2-2]=straddtagsize[1];
      fieldlabel[tagsize2-1]=straddtagsize[2];
    }

    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    double version=0.0;
    int isize, nitems;
    int iarray[10];

    char fmode[10];
    if(!strncmp(filemode,"w",1))
      strcpy(fmode,"write");
    else // default is append
      strcpy(fmode,"append");

    char datatype[10];
    if(!strncmp(arraytype,"i",1))
      strcpy(datatype,"int");
    else // default is double
      strcpy(datatype,"double");

//
// //     if(*iphase==1) //open the file but then keep it open for the remaining cycles
//     openfile_(rfile, fmode, &irstou);
//
// //     printf("iphase: %d - pid: %d - irstou %d\n",*iphase,*pid,irstou);
//
//
//     nitems = 3; // assuming field will write 3 items in iarray
//     iarray[ 0 ] = (*nshg);
//     iarray[ 1 ] = (*numvars);
//     iarray[ 2 ] = (*stepno);
//
//     isize = (*nshg)*(*numvars);
//     writeheader_( &irstou, fieldlabel, (void*)iarray, &nitems, &isize, datatype, phasta_iotype );
//
//     nitems = (*nshg)*(*numvars);
//     writedatablock_( &irstou, fieldlabel, array, &nitems, datatype, phasta_iotype );
//
// //     if(*iphase==*nphasesincycle) //close the file after nphasesincycle
//       closefile_( &irstou, fmode);
//

    /////////////////////////////// Start of writing using new-lib ////////////////////////////

    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;
//    unsigned long long timer_start;
//    unsigned long long timer_end;
//    double time_span;

    nfiles = outpar.nsynciofiles;
    nfields = outpar.nsynciofieldswriterestart;
    numparts = workfc.numpe;
    irank = *pid; //workfc.myrank;
    nprocs = workfc.numpe;

    int nppf = numparts/nfiles;
    int GPID;

    // Calculate number of parts each  proc deal with and where it start and end ...
    int nppp = numparts/nprocs;// nppp : Number of parts per proc ...
    int startpart = irank * nppp +1;// Part id from which I (myrank) start ...
    int endpart = startpart + nppp - 1;// Part id to which I (myrank) end ...

    //int descriptor;
    char filename[255],path[255],fieldtag_s[255];
    bzero((void*)filename,255);
    bzero((void*)fieldtag_s,255);

//     char * namer;
//     namer = strtok(fieldlabel," ");
//     strncpy(fieldlabel, fieldtag, *tagsize);

    field_flag++;
    if(*pid==0) {
//      printf("\n*****************************\n");
      printf("\n");
      printf("The %d/%d the field to be written is '%s'\n",field_flag,nfields,fieldlabel);
    }

    sprintf(filename,"restart-dat.%d.%d",*stepno,((int)(irank/(nprocs/nfiles))+1));

    int i;
    for ( i = 0; i < nppp; i++  ) {
        GPID = startpart + i;

        // Write solution field ...
        sprintf(fieldtag_s,"%s@%d",fieldlabel,GPID);

        //printf("This is %d and fieldtag_s is %s \n",myrank,fieldtag_s);

        isize = (*nshg)*(*numvars);
        nitems = 3;
        iarray[ 0 ] = (*nshg);
        iarray[ 1 ] = (*numvars);
        iarray[ 2 ] = (*stepno);

//MR CHANGE
//  Measure the time - Start the timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_start = rdtsc();
//MR CHANGE END

        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);

//MR CHANGE
//  Measure the time - End of timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_end = rdtsc();
//        time_span=(double)((timer_end-timer_start)/cpu_speed);
//        if (*pid==0) {
//          printf("Time: header for '%s':    %f s\n",fieldtag_s,time_span);
//        }
//MR CHANGE END

        nitems = (*nshg)*(*numvars);

//MR CHANGE
//  Measure the time - Start the timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_start = rdtsc();
//MR CHANGE END

        writedatablock( &f_descriptor, fieldtag_s, array, &isize, "double", phasta_iotype );

//MR CHANGE
//  Measure the time - End of timer
//        MPI_Barrier(MPI_COMM_WORLD);
//        timer_end = rdtsc();
//        time_span=(double)((timer_end-timer_start)/cpu_speed);

//        int isizemin,isizemax,isizetot;
//        double sizemin,sizemax,sizeavg,sizetot,rate;

//        MPI_Allreduce(&isize,&isizemin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
//        MPI_Allreduce(&isize,&isizemax,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
//        MPI_Allreduce(&isize,&isizetot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

//        sizemin=(double)(8.0*isizemin/1024.0/1024.0);
//        sizemax=(double)(8.0*isizemax/1024.0/1024.0);
//        sizetot=(double)(8.0*isizetot/1024.0/1024.0);
//        sizeavg=sizetot/workfc.numpe;
//        rate=sizetot/time_span;

//        if (*pid==0) {
//          printf("Time: block for '%s':    %f s\n",fieldtag_s,time_span);
//          printf("Time: block:   Min= %f MB; Max= %f MB; Avg= %f MB; Tot= %f MB; Rate= %f MB/s; \n",sizemin,sizemax,sizeavg,sizetot,rate);
//          printf("*****************************\n");
//        }
//MR CHANGE END

    }

//     if (*pid==0) {
//         printf("Field '%s' written in:     %f s\n",fieldtag,time_span);
//         printf("Write field '%s' finished! \n",fieldtag_s);
//     }

//
//     if (field_flag==nfields){
//       closefile_(&f_descriptor, "write");
//       finalizephmpiio_(&f_descriptor);
//       if(*pid==0) {
//         printf("Last field %d '%s' finished! \n",nfields, fieldtag_s);
//         printf("\n*****************************\n");
//       }
//     }

    if (field_flag==nfields){

//MR CHANGE
//  Measure the time - Start the timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_start = rdtsc();
//MR CHANGE END

      closefile(&f_descriptor, "write");

//MR CHANGE
//    Measure the time - End of timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_end = rdtsc();
//      time_span=(double)((timer_end-timer_start)/cpu_speed);
//     if (*pid==0) {
//        printf("\n*****************************\n");
//        printf("Time: 'closefile' is:    %f s\n",time_span);
//      }
//MR CHANGE END

//MR CHANGE
//  Measure the time - Start the timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_start = rdtsc();
//MR CHANGE END

      finalizephmpiio(&f_descriptor);

//MR CHANGE
//    Measure the time - End of timer
//      MPI_Barrier(MPI_COMM_WORLD);
//      timer_end = rdtsc();
//      time_span=(double)((timer_end-timer_start)/cpu_speed);
      if (*pid==0) {
//        printf("Time: 'finalizephmpiio' is:    %f s\n",time_span);
//        printf("Last field %d '%s' finished! \n",nfields, fieldtag);
        printf("\n");
//        printf("*****************************\n");
      }
    }
//MR CHANGE END

    ///////////////////////////////////////////////////////////////////////////////////////////

    free(fieldlabel);
}


void
Write_d2wall(   int* pid,
                int* numnp,
                double* array1 ) {

//    time_t timenow = time ( &timenow);
    int isize, nitems;
    int iarray[10];

//    MPI_Barrier(MPI_COMM_WORLD);

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
    nfields = 1; //outpar.nsynciofieldswriterestart;  // Only the distance to the walls in d2wall
    numparts = workfc.numpe;
    irank = *pid; //workfc.myrank;
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

    if (irank==0) {
      printf("Filename is %s \n",filename);
    }

    initphmpiio(&nfields, &nppf, &nfiles, &f_descriptor, "write");

    openfile(filename, "write", &f_descriptor);

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

        isize = (*numnp);
        nitems = 2;
        iarray[ 0 ] = (*numnp);
        iarray[ 1 ] = 1; //numVars = 1

        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);

        //nitems = (*nshg)*(*numVars);
        //nitems = (*numnp);

        writedatablock( &f_descriptor, fieldtag_s, (void*)(array1), &isize, "double", phasta_iotype );


    }
    field_flag++;

    if (field_flag==nfields){

      closefile(&f_descriptor, "write");

      finalizephmpiio(&f_descriptor);

      if (irank==0) {
        printf("\n");
      }
    }
}

