#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "phIO.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if( argc != 2 ) {
    fprintf(stderr, "Usage: %s <numSyncFiles>\n",argv[0]);
    MPI_Finalize();
    return 1;
  }
  const char* phrase = "number of fishes";
  const char* type = "integer";
  const char* iotype = "binary";
  int fish = 2;
  int numFish[2] = {0,0};
  int nfiles[2] = {atoi(argv[1]), 1};
  const char* filename[2] = {"water-dat.", "water.dat."};
  phio_fp file;
  int one = 1;
  int ppf = size/nfiles[0];
  for(int i=0; i<2; i++) {
    phio_openfile_write(filename[i], &(nfiles[i]), &one, &ppf, &file);
    phio_writeheader(file, phrase, &fish, &one, &one, type, iotype);
    phio_readheader(file, phrase, &(numFish[i]), &one, type, iotype);
    phio_closefile_write(file);
  }
  int match = (numFish[0] == numFish[1]);
  if(!rank && match)
    fprintf(stderr, "number of fish match!\n");
  if(!rank && !match)
    fprintf(stderr, "number of fish don't match... :(\n");
  MPI_Finalize();
  return !match;
}
