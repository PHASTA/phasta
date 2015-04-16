#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <cstring>
#include <sys/types.h>
#include <dirent.h> // For opendir(),readdir(),...
#include <mpi.h>
#include "phastaIO.h"

//#include "../include/Input.h"
#include "commonM2N_c.h"
#include "setsyncioparamM2N.h"


using namespace std; //useful for ifstream. Other solution is std::ifstream

// extern void queryphmpiio_(const char[], int*, int*);
// void queryphmpiio_(const char filename[],int *nfields, int *nppf)

void setIOparam()
{
  int count;
  int countRed;
  int nfields;
  int nppf;
  int stepno;
  int part;
  char fname[255];
  DIR *d;
  struct dirent *filename;


//Number of geombc files
  if(workfc.myrank == 0){
    if( (d = opendir(".")) != NULL) {
      count=0;
      countRed=0;
      while(filename=readdir(d)) {
        //printf("%s\n", filename->d_name);
        if(strncmp(filename->d_name,"geombc-dat",10)==0) {
          count=count+1;
        }
        if(strncmp(filename->d_name,"geombcRed-dat",13)==0) {
          countRed=countRed+1;
        }
      }
      closedir(d);

      // Sanity check
      if (count == 0) {
        printf("ERROR: Could not find any geombc-dat file in the directory\n");
      }
      if (countRed == 0) {
        printf("WARNING: Could not find any geombcRed-dat file in the directory\n");
        printf("Either you are reducing down to 1 part (geombcRed not needed) or links are missing\n");
        printf("Number of output SyncIO files will be set to either 1 or to the number of input SyncIO files read by M2N accordingly\n");
      }

    }
    else {
      printf("ERROR when counting geombc-dat\n");
    }
  }

  MPI_Bcast( &count, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &countRed, 1, MPI_INT, 0, MPI_COMM_WORLD );
  //printf("Here we gogo: %d %d\n",workfc.myrank,count);

  outpar.nsynciofiles = count;
  outpar.nsynciofilesred = countRed;
  if(workfc.myrank == 0) {
    printf("Number of geombc-dat and restart-dat files to read: %d\n", count);
    printf("Number of geombcRed-dat and restartRed-dat files to read: %d\n", countRed);
  }

}

void detectd2wallfiles(int* numd2wallfiles)
{
  int count;
  int nfields;
  int nppf;
  int stepno;
  int part;
  char fname[255];
  DIR *d;
  struct dirent *filename;

//Number of d2wall files
  if(workfc.myrank == 0){
    if( (d = opendir(".")) != NULL) {
      count=0;
      while(filename=readdir(d)) {
        //printf("%s\n", filename->d_name);
        if(strncmp(filename->d_name,"d2wall",6)==0) {
          count=count+1;
        }
      }
      closedir(d);
    }
  }

  MPI_Bcast( &count, 1, MPI_INT, 0, MPI_COMM_WORLD );
  //printf("Here we gogo: %d %d\n",workfc.myrank,count);

  *numd2wallfiles = count;
  if(workfc.myrank == 0) {
    printf("Number of d2wall files present in the proc_case directory: %d\n", count);
  }

}


// void countfieldstowriterestart() //See new_interface.c
