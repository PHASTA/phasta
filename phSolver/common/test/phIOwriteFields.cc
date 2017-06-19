#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cassert>
#include "phiostats.h"
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
  double blockArray[4] = {2.0,4.0,8.0,16.0};
  int blockEntries = 4;
  int one = 1;
  int nfields = 2;
  int numFish = 0;
  double fishWeight = 1.23;
  int nfiles = 1;
  int ppf = size/nfiles;
  const char* filename[3] = {"water.dat.", "water.","water-dat."};
  rstream rs = makeRStream();
  phio_fp file[3];
  const char* modes[3]={"syncio","posixio","streamio"};
  syncio_setup_write(nfiles, nfields, ppf, &(file[0]));
  posixio_setup(&(file[1]), 'w');
  streamio_setup_r(&(file[2]), rs, 'w');
  for(int i=0; i<3; i++) {
    if(!rank) fprintf(stderr, "%s\n", modes[i]);
    phastaio_initStats();
    phio_openfile(filename[i], file[i]);
    char str [50];
    int n;
    for (int j = 0; j < 2 ; j++) {
      n = sprintf(str, " Number of times %d ", j);
      assert(n);
      phio_writeheader(file[i], str, &blockEntries, &one, &blockEntries, type, iotype);
      phio_writedatablock(file[i], str, blockArray, &blockEntries, type, iotype);
    }
    phio_closefile(file[i]);
    phastaio_printStats();
  }
  syncio_setup_read(nfiles, &(file[0]));
  posixio_setup(&(file[1]), 'r');
  streamio_setup_r(&(file[2]), rs, 'r');
  for(int i=0; i<3; i++) {
    if(!rank) fprintf(stderr, "%s\n", modes[i]);
    phastaio_initStats();
    phio_openfile(filename[i], file[i]);
    phio_readheader(file[i], phrase, &numFish, &one, type, iotype);
    assert(!numFish);
    phio_readdatablock(file[i], phrase, &fishWeight, &numFish, type, iotype);
    assert(fishWeight == 1.23);
    phio_closefile(file[i]);
    phastaio_printStats();
  }
  MPI_Finalize();
  return 0;
}
