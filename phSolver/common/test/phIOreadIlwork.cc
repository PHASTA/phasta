#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "phIO.h"
#include "posixio.h"
#include "phio_posix.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  if( argc != 2 ) {
    fprintf(stderr, "Usage: %s <geombc posix file>\n",argv[0]);
    MPI_Finalize();
    return 1;
  }
  const char* filename = argv[1];
  const char* phrase = "ilwork";
  const char* type = "integer";
  const char* iotype = "binary";
  int* ilwork = NULL;
  int len = 0;
  int one = 2;

  phio_fp file;
  posixio_setup(&(file), 'r');
  posix_openfile_single(filename, file);
  phio_readheader(file, phrase, &len, &one, type, iotype);
  fprintf(stderr, "len %d\n", len);
  ilwork = (int*) malloc(len*sizeof(int));
  phio_readdatablock(file, phrase, ilwork, &len, type, iotype);
  phio_closefile(file);
  free(ilwork);

  MPI_Finalize();
  return 0;
}
