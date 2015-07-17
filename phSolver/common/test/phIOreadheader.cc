#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "phIO.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if( argc != 2 ) {
    fprintf(stderr, "Usage: %s <numSyncFiles>\n",argv[0]);
    MPI_Finalize();
    return 1;
  }
  const char* iotype = "binary";
  int numberOfNodes[2] = {0,0};
  int nfiles[2] = {atoi(argv[1]), 1};
  const char* dir[2] = {"4-procs_case-SyncIO-2", "4-procs_case-Posix"};
  const char* filename[2] = {"geombc-dat.", "geombc.dat."};
  phio_fp file;
  int one = 1;
  for(int i=0; i<2; i++) {
    chdir(dir[i]);
    MPI_Barrier(MPI_COMM_WORLD);
    phio_openfile_read(filename[i], &(nfiles[i]), &file);
    phio_readheader(file, "number of nodes", &(numberOfNodes[i]),
        &one, "integer", iotype);
    phio_closefile_read(file);
    chdir("..");
    MPI_Barrier(MPI_COMM_WORLD);
  }
  int match = (numberOfNodes[0] == numberOfNodes[1]);
  if(!rank && match)
    fprintf(stderr, "number of nodes match!\n");
  if(!rank && !match)
    fprintf(stderr, "number of nodes don't match... :(\n");
  MPI_Finalize();
  return !match;
}
