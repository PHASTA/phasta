#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <set>
#include "phIO.h"
#include "posixio.h"
#include "phio_posix.h"

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  if( argc != 3 ) {
    fprintf(stderr, "Usage: %s <geombc posix file> <write ilwork=0|1>\n",argv[0]);
    MPI_Finalize();
    return 1;
  }
  const char* filename = argv[1];
  const int writeilwork = atoi(argv[2]);
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

  // Initialization 
  int itkbeg = 0;
  int m = 0;
  int idl = 0;
  std::set<int> neighbors;

  // Compute summary info first such as number of communication tasks, neighboring parts, onwer parts, etc
  int numtask = ilwork[0];
  int numowner = 0;
  for(int itask=0;itask<numtask;itask++) {
    int iacc   = ilwork [itkbeg + 2];
    int iother = ilwork [itkbeg + 3];
    int numseg = ilwork [itkbeg + 4];
    if(iacc == 1) numowner++;
    neighbors.insert(iother);
    itkbeg = itkbeg + 4 + 2*numseg;
  }
  printf("Number of neighboring parts: %d\n",neighbors.size());
  printf("Number of communication tasks: %d\n",numtask);
  printf("Number of owner communications: %d\n",numowner);
  printf("Number of non-owner communications: %d\n",numtask-numowner);

  // Print now communication info
  printf("\n");
  printf("Communication info for this part:\n");
  itkbeg = 0;
  for(int itask=0;itask<numtask;itask++) {
    int itag   = ilwork [itkbeg + 1];
    int iacc   = ilwork [itkbeg + 2];
    int iother = ilwork [itkbeg + 3];
    int numseg = ilwork [itkbeg + 4];
    int isgbeg = ilwork [itkbeg + 5];

    // Toal number of nodes involved in this communication (lfront), including all segments
    int lfront = 0;
    for(int is=1;is<=numseg;is++) {
      int lenseg = ilwork [itkbeg + 4 + 2*is];
      lfront = lfront + lenseg;
    }
    printf("Communication %d:\ttag %d\townership %d\trank %d\tnumseg %d\tnumvtx %d\n",itask,itag,iacc,iother,numseg,lfront);
    itkbeg = itkbeg + 4 + 2*numseg;
  }

  // Print now the raw ilwork array
  printf("\n");
  if( writeilwork ) {
    printf("ilwork array:\n");
    for(int i=0;i<len;i++) {
      printf("%d\n",ilwork[i]);
    }
  }

  free(ilwork);

  MPI_Finalize();
  return 0;
}
