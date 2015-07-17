#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "phIO.h"
#include "phComm.h"
#include "phio_base.h"
#include "phio_sync.h"
#include "phio_posix.h"

void phio_readheader(
    phio_fp f,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  f->ops->readheader(f->file, keyphrase, valueArray,
      nItems, datatype, iotype);
}
void phio_writeheader( 
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const int* ndataItems,
    const char datatype[],
    const char iotype[] ) {
}
void phio_readdatablock(
    phio_fp f,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  f->ops->readdatablock(f->file, keyphrase, valueArray,
      nItems, datatype, iotype);
}
void phio_writedatablock(
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[]) {
}
void phio_openfile_read(
    const char filename[],
    int* numFiles,
    phio_fp* fileDescriptor) {
  std::string fn(filename);
  std::string syncSuffix("-dat");
  std::string posixSuffix(".dat");
  if( fn.find(syncSuffix) != std::string::npos ) 
    sync_openfile_read(filename, numFiles, fileDescriptor);
  else if( fn.find(posixSuffix) != std::string::npos ) 
    posix_openfile_read(filename, numFiles, fileDescriptor);
  else {
    fprintf(stderr,
        "type of file %s is unknown... exiting\n", filename);
    exit(1);
  }
}
void phio_openfile_write(
    const char filename[],
    int* numFiles,
    int* numFields,
    int* numPPF,
    int* fileDescriptor) {
}
void phio_restartname(int* step, char* filename) {
}
void phio_closefile_read(phio_fp f) {
  f->ops->closefile_read(f);
}
void phio_closefile_write(int* fileDescriptor) {
}

