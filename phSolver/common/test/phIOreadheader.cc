#include <mpi.h>
#include <assert.h>
#include "phIO.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  int numberOfNodes[2];
  const char* iotype = "binary";
  int nfiles[2] = {1,2};
  int file = 0;
  int one = 1;
  for(int i=0; i<2; i++) {
    phio_openfile_read("geombc", nfiles, &file);
    phio_readheader(&file, "number of nodes", &(numberOfNodes[i]),
        &one, "integer", iotype);
    phio_closefile_read(&file);
  }
  int match = (numberOfNodes[0] == numberOfNodes[1]);
  MPI_Finalize();
  return match;
}
