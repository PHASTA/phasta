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

#include "commonM2NFixBnd_c.h"
#include "setsyncioparamM2NFixBnd.h"


using namespace std; //useful for ifstream. Other solution is std::ifstream

// extern void queryphmpiio_(const char[], int*, int*);
// void queryphmpiio_(const char filename[],int *nfields, int *nppf)

void setIOparam()
{
  int count;
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
      while(filename=readdir(d)) {
        //printf("%s\n", filename->d_name);
        if(strncmp(filename->d_name,"geombcRed-dat",13)==0) {
          count=count+1;
        }
      }
      closedir(d);
    }
    else {
      printf("ERROR when counting geombcRed-dat\n");
    }
  }

  MPI_Bcast( &count, 1, MPI_INT, 0, MPI_COMM_WORLD );
  //printf("Here we gogo: %d %d\n",workfc.myrank,count);

  outpar.nsynciofiles = count;
  if(workfc.myrank == 0) {
    printf("Number of geombcRed-dat and restartRed-dat files to read: %d\n", count);
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
