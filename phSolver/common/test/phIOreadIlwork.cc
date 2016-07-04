#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <set>
#include <string>
#include <sstream>
#include <cassert>
#include "phIO.h"
#include "posixio.h"
#include "phio_posix.h"

int openfile(const char* geomfilename, phio_fp& file) {
  int commrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
  int err = posix_openfile_single(geomfilename, file);
  int globalerr = 0;
  MPI_Allreduce(&err, &globalerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(err > 0 && !commrank)
    fprintf(stderr, "failed to open file %s\n", geomfilename);
  return globalerr;
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  int commrank,commsize;
  MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
  MPI_Comm_size(MPI_COMM_WORLD,&commsize);
  if( argc != 5 ) {
    if( !commrank )
      fprintf(stderr, "Usage: %s <geombc posix file> <verbosity=0|1|2> <rankoffset> <outfile>\n",argv[0]);
      fprintf(stderr, "verbosity=0 will only write to the specified \'outfile\'\n");
      fprintf(stderr, "verbosity>0 will write to the specified \'outfile\' and to stdout\n");
    MPI_Finalize();
    return 1;
  }
  const char* filename = argv[1];
  const int verbosity = atoi(argv[2]);
  const int rankoffset = atoi(argv[3]);
  char* outfilename = argv[4];
  const char* phrase = "ilwork";
  const char* type = "integer";
  const char* iotype = "binary";
  int* ilwork = NULL;
  int len = 0;
  int one = 2;

  const int phastarank = commrank+1+rankoffset;
  const int group = (commrank+rankoffset)/2048;

  char geomfilename[4096];
  sprintf(geomfilename, "%s/%d/geombc.dat.%d",filename,group,phastarank);

  phio_fp file;
  posixio_setup(&(file), 'r');
  int err = openfile(geomfilename, file);
  if(err > 0) {
    if(!commrank) fprintf(stderr, "trying again without the sub-directory...\n");
    sprintf(geomfilename, "%s/geombc.dat.%d",filename,phastarank);
    err = openfile(geomfilename, file);
    assert(!err);
    if(!commrank)
      fprintf(stderr, "geombc files opened successfully\n");
  }
  phio_readheader(file, phrase, &len, &one, type, iotype);
  ilwork = (int*) malloc(len*sizeof(int));
  phio_readdatablock(file, phrase, ilwork, &len, type, iotype);
  phio_closefile(file);

  // Initialization
  int itkbeg = 0;
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
  assert(neighbors.size() != 0);
  MPI_Status status;
  MPI_File outfile;
  MPI_File_open(MPI_COMM_WORLD,outfilename,
      MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&outfile);
  std::string header("rank,peers,tasks,owned,notowned\n");
  if( !commrank ) //write header
    MPI_File_write_at(outfile,0,(void*)header.c_str(),header.size(),MPI_BYTE,&status);
  std::stringstream ss;
  ss << phastarank << "," 
     << neighbors.size() << ","
     << numtask << ","
     << numowner << ","
     << numtask-numowner << "\n";
  std::string s = ss.str();
  int size = s.size();
  int offset = 0;
  MPI_Exscan(&size,&offset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  offset += header.size();
  int ret = MPI_File_write_at(outfile,offset,(void*)s.c_str(),s.size(),MPI_BYTE,&status);
  assert(ret == MPI_SUCCESS);
  if( verbosity > 0 ) {
    // Print now communication info
    printf("\n");
    printf("Communication info for this part:\n");
    itkbeg = 0;
    for(int itask=0;itask<numtask;itask++) {
      int itag   = ilwork [itkbeg + 1];
      int iacc   = ilwork [itkbeg + 2];
      int iother = ilwork [itkbeg + 3];
      int numseg = ilwork [itkbeg + 4];

      // Toal number of nodes involved in this communication (lfront), including all segments
      int lfront = 0;
      for(int is=1;is<=numseg;is++) {
        int lenseg = ilwork [itkbeg + 4 + 2*is];
        lfront = lfront + lenseg;
      }
      printf("Communication %d:\ttag %d\townership %d\trank %d\tnumseg %d\tnumvtx %d\n",itask,itag,iacc,iother,numseg,lfront);
      itkbeg = itkbeg + 4 + 2*numseg;
    }
    printf("\n");
  }

  // Print now the raw ilwork array
  if( verbosity > 1 ) {
    printf("ilwork array:\n");
    for(int i=0;i<len;i++) {
      printf("%d\n",ilwork[i]);
    }
  }

  free(ilwork);
  MPI_File_close(&outfile);
  MPI_Finalize();
  return 0;
}
