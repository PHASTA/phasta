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
  const char* iotype = "binary";
  int seven = 7;
  int headerData[7] = {0,0,0,0,0,0,0};
  int blocksRead = 0;
  phio_fp file;
  posixio_setup(&file, 'r');
  phio_openfile("geombc.dat.", file);
  MPI_Barrier(MPI_COMM_WORLD);
  do {
    phio_readheader(file, "connectivity interior",
        headerData, &seven, "integer", iotype);
    fprintf(stderr, "rank %d data[0] %d\n", rank, headerData[0]);
    blocksRead += (headerData[0] > 0);
  } while( headerData[0] > 0 );
  phio_closefile(file);
  MPI_Barrier(MPI_COMM_WORLD);
  fprintf(stderr, "rank %d number of blocks read %d\n", rank, blocksRead);
  MPI_Finalize();
  return (blocksRead==2);
}
