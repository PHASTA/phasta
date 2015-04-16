/* This file provides interface functions for 'partial ' random
   access into the PHASTA input files

   Anil Karanam March 2001 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "phastaIO.h"
#include <FCMangle.h>
#include "new_interfaceAcuStat.h"

//MR CHANGE
#include "commonAcuStat_c.h"
//MR CHANGE END

#ifdef intel
#include <winsock2.h>
#else
#include <unistd.h>
#include <strings.h>
#endif


void
Write_AcuStat(  int* pid,
                int* stepstart,
                int* stepend,
                int* steptot,
                int* nshg,
                int* numVars,
                double* array1 ) {

    char fname[255];
    char rfile[60];
    char existingfile[30], linkfile[30];
    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    double version=0.0;
    int isize, nitems;
    int iarray[10];
    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;


    irank = *pid;

    //  Retrieve and compute the parameters required for SyncIO
    nfiles = outpar.nsynciofiles; //We use the same number of files as before
    nfields = outpar.nsynciofieldswriterestart; // >=4
    numparts = workfc.numpe;
    nprocs= workfc.numpe;

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

    //sprintf(filename,"restartAcuStat-dat.%d.%d",*stepend,((int)(irank/(nprocs/nfiles))+1));
    sprintf(filename,"restartAcuStat.%d.%d.%d.%d",*stepstart,*stepend,*steptot,((int)(irank/(nprocs/nfiles))+1));
    if (*pid==0) {
	    printf("Filename is %s \n",filename);
    }

    initphmpiio(&nfields, &nppf, &nfiles, &f_descriptor, "write");
    openfile(filename, "write", &f_descriptor);

    field_flag=0;
    // solution 
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
	    iarray[ 2 ] = (*stepend);
	    writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
	    writedatablock( &f_descriptor, fieldtag_s, (void*)(array1), &isize, "double", phasta_iotype );
    }
    field_flag++;

    if (field_flag==nfields){
      closefile(&f_descriptor, "write");
      finalizephmpiio(&f_descriptor);
      printf("file closed\n");
      if (*pid==0) {
        printf("\n");
      }
    }

    return;
}



void
Write_AcuStat2(  int* pid,
                int* stepstart,
                int* stepend,
                int* steptot,
                int* nshg,
                int* numVars,
                int* ndoferrors,
                int* ndofybar,
                double* array1,
                double* array2,
                double* array3,
                double* array4 ) {

    char fname[255];
    char rfile[60];
    char existingfile[30], linkfile[30];
    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    double version=0.0;
    int isize, nitems;
    int iarray[10];
    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;


    irank = *pid;

    //  Retrieve and compute the parameters required for SyncIO
    nfiles = outpar.nsynciofiles; //We use the same number of files as before
    nfields = outpar.nsynciofieldswriterestart; // >=4
    numparts = workfc.numpe;
    nprocs= workfc.numpe;

    int nppf = numparts/nfiles;
    int GPID;

    // Calculate number of parts each proc deal with and where it start and end ...
    int nppp = numparts/nprocs;// nppp : Number of parts per proc ...
    int startpart = irank * nppp +1;// Part id from which I (myrank) start ...
    int endpart = startpart + nppp - 1;// Part id to which I (myrank) end ...

    int descriptor;
    char filename[255],fieldtag_s[255];
    bzero((void*)filename,255);
    bzero((void*)fieldtag_s,255);

    sprintf(filename,"restartAcuStat.%d.%d.%d.%d",*stepstart,*stepend,*steptot,((int)(irank/(nprocs/nfiles))+1));
    if (*pid==0) {
	    printf("Filename is %s \n",filename);
    }

    initphmpiio(&nfields, &nppf, &nfiles, &f_descriptor, "write");
    openfile(filename, "write", &f_descriptor);

    field_flag=0;
    // solution 
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
	    iarray[ 2 ] = (*stepend);
	    writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
	    writedatablock( &f_descriptor, fieldtag_s, (void*)(array1), &isize, "double", phasta_iotype );
    }
    field_flag++;

    // dwal
    for ( i = 0; i < nppp; i++) { //This loop is useful only if several parts per processor
	    // GPID : global part id, corresponds to rank ...
	    // e.g : (in this example)
	    // proc 0 : 1--4
	    // proc 1 : 5--8 ...
	    GPID = startpart + i;

	    // Write solution field ...
	    sprintf(fieldtag_s,"dwal@%d",GPID);

	    isize = (*nshg)*1;
	    nitems = 3;
	    iarray[ 0 ] = (*nshg);
	    iarray[ 1 ] = 1;
	    iarray[ 2 ] = (*stepend);
	    writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
	    writedatablock( &f_descriptor, fieldtag_s, (void*)(array2), &isize, "double", phasta_iotype );
    }
    field_flag++;

    // errors
    for ( i = 0; i < nppp; i++) { //This loop is useful only if several parts per processor
	    // GPID : global part id, corresponds to rank ...
	    // e.g : (in this example)
	    // proc 0 : 1--4
	    // proc 1 : 5--8 ...
	    GPID = startpart + i;

	    // Write solution field ...
	    sprintf(fieldtag_s,"errors@%d",GPID);

	    isize = (*nshg)*(*ndoferrors);
	    nitems = 3;
	    iarray[ 0 ] = (*nshg);
	    iarray[ 1 ] = (*ndoferrors);
	    iarray[ 2 ] = (*stepend);
	    writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
	    writedatablock( &f_descriptor, fieldtag_s, (void*)(array3), &isize, "double", phasta_iotype );
    }
    field_flag++;

    // ybar
    for ( i = 0; i < nppp; i++) { //This loop is useful only if several parts per processor
	    // GPID : global part id, corresponds to rank ...
	    // e.g : (in this example)
	    // proc 0 : 1--4
	    // proc 1 : 5--8 ...
	    GPID = startpart + i;

	    // Write solution field ...
	    sprintf(fieldtag_s,"ybar@%d",GPID);

	    isize = (*nshg)*(*ndofybar);
	    nitems = 3;
	    iarray[ 0 ] = (*nshg);
	    iarray[ 1 ] = (*ndofybar);
	    iarray[ 2 ] = (*stepend);
	    writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
	    writedatablock( &f_descriptor, fieldtag_s, (void*)(array4), &isize, "double", phasta_iotype );
    }
    field_flag++;


    if (field_flag==nfields){
      closefile(&f_descriptor, "write");
      finalizephmpiio(&f_descriptor);
      printf("file closed\n");
      if (*pid==0) {
        printf("\n");
      }
    }


    return;
}

void
Write_Ybar(     int* pid,
                int* stepstart,
                int* stepend,
                int* steptot,
                int* nshg,
                int* numVars,
                double* array1) {

    char fname[255];
    char rfile[60];
    char existingfile[30], linkfile[30];
    int irstou;
    int magic_number = 362436;
    int* mptr = &magic_number;
    double version=0.0;
    int isize, nitems;
    int iarray[10];

    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

    //  First, count the number of fields to write and store the result in
    //countfieldstowriterestart();

    //  Retrieve and compute the parameters required for SyncIO
    nfiles = outpar.nsynciofiles;
    //nfields = outpar.nsynciofieldswriterestart;
    nfields=1;
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
    char filename[255],fieldtag_s[255];
    bzero((void*)filename,255);
    bzero((void*)fieldtag_s,255);

//    sprintf(filename,"ybar.%d.%d",*stepno,((int)(irank/(nprocs/nfiles))+1));
    sprintf(filename,"restartAcuStat.%d.%d.%d.%d",*stepstart,*stepend,*steptot,((int)(irank/(nprocs/nfiles))+1));

    initphmpiio(&nfields, &nppf, &nfiles, &f_descriptor, "write");

    if (*pid==0) {
      printf("Filename is %s \n",filename);
    }

    openfile(filename, "write", &f_descriptor);

    field_flag=0;

     int i;
     for ( i = 0; i < nppp; i++) { //This loop is useful only if several parts per processor
     // GPID : global part id, corresponds to rank ...
        // e.g : (in this example)
        // proc 0 : 1--4
        // proc 1 : 5--8 ...
        GPID = startpart + i;

        // Write ybar field ...
        sprintf(fieldtag_s,"ybar@%d",GPID);

        isize = (*nshg)*(*numVars);
        nitems = 3;
        iarray[ 0 ] = (*nshg);
        iarray[ 1 ] = (*numVars);
        iarray[ 2 ] = (*stepend);
        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
        writedatablock( &f_descriptor, fieldtag_s, (void*)(array1), &isize, "double", phasta_iotype );
    }
    field_flag++;

    if (field_flag==nfields){
      closefile(&f_descriptor, "write");
      finalizephmpiio(&f_descriptor);
      if (*pid==0) {
        printf("\n");
      }
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
    double version=0.0;
    int isize, nitems;
    int iarray[10];
    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

    //  First, count the number of fields to write and store the result in
    //countfieldstowriterestart();

    //  Retrieve and compute the parameters required for SyncIO
    nfiles = outpar.nsynciofiles;
    nfields = outpar.nsynciofieldswriterestart;
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
    char filename[255],fieldtag_s[255];
    bzero((void*)filename,255);
    bzero((void*)fieldtag_s,255);

    sprintf(filename,"restart-dat.%d.%d",*stepno,((int)(irank/(nprocs/nfiles))+1));

    initphmpiio(&nfields, &nppf, &nfiles, &f_descriptor, "write");

    if (*pid==0) {
      printf("Filename is %s \n",filename);
    }

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
        sprintf(fieldtag_s,"solution@%d",GPID);

        isize = (*nshg)*(*numVars);
        nitems = 3;
        iarray[ 0 ] = (*nshg);
        iarray[ 1 ] = (*numVars);
        iarray[ 2 ] = (*stepno);
        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
        writedatablock( &f_descriptor, fieldtag_s, (void*)(array1), &isize, "double", phasta_iotype );
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
        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
        writedatablock( &f_descriptor, fieldtag_s, (void*)(array2), &isize, "double", phasta_iotype );
     }
     field_flag++;

     if (field_flag==nfields){
       closefile(&f_descriptor, "write");
       finalizephmpiio(&f_descriptor);
       if (*pid==0) {
         printf("\n");
       }
     }
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
    double version=0.0;
    int isize, nitems;
    int iarray[10];

    /////////////////////////////// Start of writing using new-lib ////////////////////////////

    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

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
          printf("\n");
          printf("The %d/%d th field to be written is '%s'\n",field_flag,nfields,fieldtag);
        }

        isize = (*nshg)*(*numVars);
        nitems = 3;
        iarray[ 0 ] = (*nshg);
        iarray[ 1 ] = (*numVars);
        iarray[ 2 ] = (*stepno);
        writeheader( &f_descriptor, fieldtag, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
        writedatablock( &f_descriptor, fieldtag, (void*)array1, &isize, "double", phasta_iotype );
    }

    if (field_flag==nfields){
      closefile(&f_descriptor, "write");
      finalizephmpiio(&f_descriptor);
      if (*pid==0) {
        printf("Last field %d '%s' finished! \n",nfields, fieldtag);
        printf("\n");
      }
    }

}

void
Write_Displ(  int* pid,
              int* stepno,
              int* nshg,
              int* numVars,
              double* array1 ) { //TO BE UPDATED FOR SYNCIO


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
    writedatablock( &irstou, "displacement", (void*)(array1), &isize, "double", phasta_iotype );

    closefile( &irstou, "append" );
}

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

//   Old posix format
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

    char fieldtag_s[255];
    bzero((void*)fieldtag_s,255);

    strncpy(fieldlabel, fieldtag, *tagsize);

    field_flag++;
    if(*pid==0) {
      printf("\n");
      printf("The %d/%d th field to be written is '%s'\n",field_flag,nfields,fieldlabel);
    }

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
        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, datatype, phasta_iotype);
        writedatablock( &f_descriptor, fieldtag_s, array, &isize, datatype, phasta_iotype );
    }

    if (field_flag==nfields){
      closefile(&f_descriptor, "write");
      finalizephmpiio(&f_descriptor);
      if (*pid==0) {
        printf("Last field %d '%s' finished! \n",nfields, fieldtag);
        printf("\n");
      }
    }

    free(fieldlabel);
}

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
      writedatablock( &irstou, "byteorder magic number ",
                       (void*)mptr, &isize, "integer", phasta_iotype );
    }

    isize = (*nshg)*(*numvars);
    nitems = 3; // assuming field will write 3 items in iarray
    iarray[ 0 ] = (*nshg);
    iarray[ 1 ] = (*numvars);
    iarray[ 2 ] = (*stepno);
    writeheader( &irstou, fieldlabel, (void*)iarray, &nitems, &isize, datatype, phasta_iotype );
    writedatablock( &irstou, fieldlabel, array, &isize, datatype, phasta_iotype );

    closefile( &irstou, fmode);

    free(fieldlabel);
}

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

    int addtagsize; // phase number is added to the name of the field
    if(*iphase<10)
      addtagsize=1;
    else if(*iphase<100)
      addtagsize=2;
    else if(*iphase<1000)
      addtagsize=3;

    int tagsize2;
    tagsize2=*tagsize+addtagsize;

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

    /////////////////////////////// Start of writing using new-lib ////////////////////////////

    int nfiles;
    int nfields;
    int numparts;
    int irank;
    int nprocs;

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

    char fieldtag_s[255];
    bzero((void*)fieldtag_s,255);

    field_flag++;
    if(*pid==0) {
      printf("\n");
      printf("The %d/%d th field to be written is '%s'\n",field_flag,nfields,fieldlabel);
    }

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
        writeheader( &f_descriptor, fieldtag_s, (void*)iarray, &nitems, &isize, "double", phasta_iotype);
        writedatablock( &f_descriptor, fieldtag_s, array, &isize, "double", phasta_iotype );
    }

    if (field_flag==nfields){
      closefile(&f_descriptor, "write");
      finalizephmpiio(&f_descriptor);
      if (*pid==0) {
        printf("\n");
      }
    }

    free(fieldlabel);
}


void
Write_d2wall(   int* pid,
                int* numnp,
                double* array1 ) {

    int isize, nitems;
    int iarray[10];

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
    char filename[255],fieldtag_s[255];
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

