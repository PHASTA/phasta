#include "phIO.h"
#include "phio_base.h"
#include "phio_posix.h"
#include "phComm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <phastaIO.h>
#include <sstream>
#include <string>

namespace {
  std::string appendPosix(const char* phrase) {
    std::stringstream ss;
    ss << phrase << "?";
    return ss.str();
  }
  std::string appendRank(const char* phrase) {
    std::stringstream ss;
    ss << phrase << phcomm_rank()+1;
    return ss.str();
  }
}

static phio_ops posix_ops = {
  posix_readheader,
  posix_writeheader,
  posix_readdatablock,
  posix_writedatablock,
  posix_restartname,
  posix_closefile_read,
  posix_closefile_write
};

void posix_readheader(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  std::string posixPhrase = appendPosix(keyphrase);
  readheader(fileDescriptor, posixPhrase.c_str(),
      valueArray, nItems, datatype, iotype);
}

void posix_writeheader(
      const int* fileDescriptor,
      const char keyphrase[],
      const void* valueArray,
      const int* nItems,
      const int* ndataItems,
      const char datatype[],
      const char iotype[] ) {
  std::string posixPhrase = appendPosix(keyphrase);
  writeheader(fileDescriptor, posixPhrase.c_str(),
      valueArray, nItems, ndataItems, datatype, iotype);
}


void posix_readdatablock(
    int*  fileDescriptor,
    const char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  std::string posixPhrase = appendPosix(keyphrase);
  readdatablock(fileDescriptor, posixPhrase.c_str(),
      valueArray, nItems, datatype, iotype);
}

void posix_writedatablock(
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[]) {
  std::string posixPhrase = appendPosix(keyphrase);
  writedatablock(fileDescriptor, posixPhrase.c_str(),
      valueArray, nItems, datatype, iotype);
}

void posix_openfile_read(
    const char filename[],
    int* numFiles,
    phio_fp* fileDescriptor) {
  *fileDescriptor =
    (struct phio_file*) malloc(sizeof(struct phio_file));
  (*fileDescriptor)->ops = &posix_ops; 
  (*fileDescriptor)->file = (int*) malloc(sizeof(int*));
  const char* mode = "read";
  std::string posixName = appendRank(filename);
  openfile(posixName.c_str(), mode, (*fileDescriptor)->file);
}

void posix_openfile_write(
    const char filename[],
    int* numFiles,
    int* numFields,
    int* numPPF,
    int* fileDescriptor) {
  //TODO - define a good upper bound
  assert(*numFields > 0 && *numFields < 1024);
  assert(*numPPF > 0 && *numPPF < 1024);
  const char* mode = "write";
  std::string posixName = appendRank(filename);
  openfile(posixName.c_str(), mode, fileDescriptor);
}

void posix_restartname(int* step, char* filename) {
  std::stringstream ss;
  ss << "restart.dat." << *step << '.';
  std::string s = ss.str();
  strcpy(filename, s.c_str()); 
}

void posix_closefile_read(phio_fp f) {
  const char* mode = "read";
  closefile(f->file, mode);
  free(f->file);
  free(f);
}

void posix_closefile_write(int* fileDescriptor) {
  const char* mode = "write";
  closefile(fileDescriptor, mode);
}
