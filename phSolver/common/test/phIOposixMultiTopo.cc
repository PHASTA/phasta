#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "phIO.h"
#include "syncio.h"
#include "posixio.h"
#include <assert.h>

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
  for(int i=0;i<2;i++) {
    phio_readheader(file, "connectivity interior",
        headerData, &seven, "integer", iotype);
    assert(headerData[0] > 0 && headerData[3] > 0);
    int size = headerData[0]*headerData[3]; /* neltp*nshl */
    int* vals = (int*) calloc(size,sizeof(int));
    phio_readdatablock(file,"connectivity interior",vals,&size,"integer",iotype);
    free(vals);
    blocksRead += (headerData[0] > 0);
  }
  phio_closefile(file);
  MPI_Finalize();
  return !(blocksRead==2);
}
