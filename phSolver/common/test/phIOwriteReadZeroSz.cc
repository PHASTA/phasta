#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cassert>
#include "phIO.h"
#include "phstream.h" //for makeRStream and makeGRStream
#include "syncio.h"
#include "posixio.h"
#include "streamio.h"

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
  int zero = 0;
  int one = 1;
  int numFish = 0;
  double fishWeight = 1.23;
  int nfiles = atoi(argv[1]);
  int ppf = size/nfiles;
  const char* filename[2] = {"water-dat.", "water.dat."};
  rstream rs = makeRStream();
  grstream grs = makeGRStream();
  phio_fp file[3];
  syncio_setup_write(nfiles, one, ppf, &(file[0]));
  posixio_setup(&(file[1]), 'w');
  streamio_setup_write(&(file[2]), rs);
  for(int i=0; i<3; i++) {
    phio_openfile(filename[i], file[i]);
    phio_writeheader(file[i], phrase, &zero, &one, &zero, type, iotype);
    phio_writedatablock(file[i], phrase, &fishWeight, &zero, type, iotype);
    phio_closefile(file[i]);
  }
  syncio_setup_read(nfiles, &(file[0]));
  posixio_setup(&(file[1]), 'r');
  streamio_read_r(&(file[2]), rs);
  for(int i=0; i<3; i++) {
    phio_openfile(filename[i], file[i]);
    phio_readheader(file[i], phrase, &numFish, &one, type, iotype);
    assert(!numFish);
    phio_readdatablock(file[i], phrase, &fishWeight, &numFish, type, iotype);
    assert(fishWeight == 1.23);
    phio_closefile(file[i]);
  }
  MPI_Finalize();
  return 0;
}
