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
  int zero = 0;
  int one = 1;
  int numFish = 0;
  double fishWeight = 1.23;
  int nfiles = 1;
  int ppf = size/nfiles;
  const char* filename[2] = {"water-dat.", "water.dat."};
  rstream rs = makeRStream();
     phio_fp file[3];
  const char* modes[3]={"syncio", "posixio", "streamio"};
  syncio_setup_write(nfiles, one, ppf, &(file[0]));
  posixio_setup(&(file[1]), 'w');
  streamio_setup_r(&(file[2]), rs, 'w');
  fprintf(stderr, "%s\n" ,"Outside loop 1.0"); 
  for(int i=0; i<3; i++) {
    fprintf(stderr, "%s\n" ,"Within the i loop");
    if(!rank) fprintf(stderr, "%s\n", modes[i]);
    fprintf(stderr, "%s\n" ,"Before phastaio");
    phastaio_initStats();
    fprintf(stderr, "%s\n" ,"Opening files with ", filename[i], file[i] );
    phio_openfile(filename[i], file[i]);
    fprintf(stderr, "%s\n" ,"Entering for loop for ", atoi(argv[1]) );
   // const char* str = "Number of times "+ nfiles ;
    for (int j = 0; j < nfiles ; j++) {
      fprintf(stderr,"%s\n", "Inside loop");
      fprintf(stderr,"%d\n",atoi(argv[1]));
      const char* str = "Number of times " +j;
      fprintf(stderr,"%s\n", "Writing the header time - " );
      fprintf(stderr,"%d\n",zero);
      fprintf(stderr,"%d\n",one);
      fprintf(stderr,"%s\n",file[i] );
      fprintf(stderr,"%s\n","Should have printed the file" );
      fprintf(stderr,"%s\n",type);
      fprintf(stderr,"%s\n",iotype);
      fprintf(stderr,"%s\n","Opening the file" );
      phio_writeheader(file[i], str, &zero, &one, &zero, type, iotype);
      fprintf(stderr,"%s\n",str );
      fprintf(stderr,"%d\n",j );
      fprintf(stderr,"%s\n", "Writing the data block time - ");
      fprintf(stderr,"%d\n",j );
      phio_writedatablock(file[i], str, &fishWeight, &zero, type, iotype);
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
    //Str was added
   // const char* str = "Number of times "+ nfiles ;
   // for (int j = 0; j < nfiles ; j++) {
    phio_readheader(file[i], phrase, &numFish, &one, type, iotype);
    //  phio_readheader(file[i], str, &numFish, &one, type, iotype);
     //Changing argument from file[i] to file[j]
    assert(!numFish);
    phio_readdatablock(file[i], phrase, &fishWeight, &numFish, type, iotype);
    //  phio_readdatablock(file[i], str, &fishWeight, &numFish, type, iotype);
    assert(fishWeight == 1.23);
    phio_closefile(file[i]);
    phastaio_printStats();
  }
  MPI_Finalize();
  return 0;
}
