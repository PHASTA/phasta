#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "phIO.h"
#include "syncio.h"
#include "posixio.h"

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
  const char* type = "double";
  const char* iotype = "binary";
  int one = 1;
  int fish = 2;
  int numFish[2] = {0,0};
  double fishWeight[2] = {1.23,1.23};
  int nfiles = atoi(argv[1]);
  int ppf = size/nfiles;
  const char* filename[2] = {"water-dat.", "water.dat."};
  phio_fp file[2];
  syncio_setup_write(nfiles, one, ppf, &(file[0]));
  posixio_setup(&(file[1]), 'w');
  for(int i=0; i<2; i++) {
    phio_openfile(filename[i], file[i]);
    phio_writeheader(file[i], phrase, &fish, &one, &one, type, iotype);
    phio_writedatablock(file[i], phrase, &(fishWeight[i]), &one, type, iotype);
    phio_closefile(file[i]);
  }
  syncio_setup_read(nfiles, &(file[0]));
  posixio_setup(&(file[1]), 'r');
  for(int i=0; i<2; i++) {
    phio_openfile(filename[i], file[i]);
    phio_readheader(file[i], phrase, &(numFish[i]), &one, type, iotype);
    phio_readdatablock(file[i], phrase, &(fishWeight[i]), &one, type, iotype);
    phio_closefile(file[i]);
  }
  int match = (numFish[0] == numFish[1]) && (fishWeight[0] == fishWeight[1]);
  if(!rank && match)
    fprintf(stderr, "number of fish match!\n");
  if(!rank && !match)
    fprintf(stderr, "number of fish don't match... :(\n");
  MPI_Finalize();
  return !match;
}
