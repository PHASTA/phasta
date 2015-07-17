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
  int nfiles[2] = {atoi(argv[1]), 1};
  const char* phrase = "co-ordinates";
  const char* type = "double";
  const char* iotype = "binary";
  const char* dir[2] = {"4-procs_case-SyncIO-2", "4-procs_case-Posix"};
  const char* filename[2] = {"geombc-dat.", "geombc.dat."};
  double* coords[2] = {NULL, NULL};
  int len[2] = {0, 0};
  int numpts[2];
  phio_fp file;
  int two = 2;
  for(int i=0; i<2; i++) {
    chdir(dir[i]);
    MPI_Barrier(MPI_COMM_WORLD);
    phio_openfile_read(filename[i], &(nfiles[i]), &file);
    phio_readheader(file, phrase, numpts, &two, type, iotype);
    len[i] = numpts[0]*3; //numPts * 3 dimensions
    fprintf(stderr, "%d %s len %d\n", rank, __func__, len[i]);
    coords[i] = (double*) malloc(len[i]*sizeof(double));
    phio_readdatablock(file, phrase, coords[i], &(len[i]), type, iotype);
    phio_closefile_read(file);
    chdir("..");
    MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) 
      fprintf(stderr, "-------------------\n");
  }
  int match = (len[0] == len[1]);
  if(!rank && match)
    fprintf(stderr, "number of points match!\n");
  if(!rank && !match) {
    fprintf(stderr, "number of points don't match... :(\n");
    return 1;
  }

  match = true;
  for(int i=0; i<len[0]; i++)
    match = match && ( coords[0][i] == coords[1][i] );

  if(!rank && match)
    fprintf(stderr, "points match!\n");
  if(!rank && !match) {
    fprintf(stderr, "points don't match... :(\n");
    return 1;
  }

  for(int i=0; i<2; i++)
    free(coords[i]);

  MPI_Finalize();
  return !match;
}
