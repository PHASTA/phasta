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
  std::string appendRank(const char* phrase) {
    std::stringstream ss;
    ss << phrase << phcomm_rank()+1;
    return ss.str();
  }
  void close(phio_fp f, const char* mode) {
    closefile(f->file, mode);
    free(f->file);
    free(f);
  }
}

void posix_openfile(const char filename[], phio_fp f) {
  std::string posixName = appendRank(filename);
  int err = posix_openfile_single(posixName.c_str(),f);
  assert(!err);
}

int posix_openfile_single(const char filename[], phio_fp f) {
  assert(f->mode == 'r' || f->mode == 'w');
  std::string posixName(filename);
  if(f->mode == 'r')
    openfile(posixName.c_str(), "read", f->file);
  else if(f->mode == 'w')
    openfile(posixName.c_str(), "write", f->file);
  if( ! *(f->file) )
    return 1;
  else
    return 0;
}

void posix_closefile(phio_fp f) {
  assert(f->mode == 'r' || f->mode == 'w');
  if(f->mode == 'r')
    close(f, "read");
  else if(f->mode == 'w')
    close(f, "write");
}

void posix_readheader(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  readheader(fileDescriptor, keyphrase,
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
  writeheader(fileDescriptor, keyphrase, valueArray,
      nItems, ndataItems, datatype, iotype);
}

void posix_readdatablock(
    int*  fileDescriptor,
    const char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  readdatablock(fileDescriptor, keyphrase,
      valueArray, nItems, datatype, iotype);
}

void posix_writedatablock(
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[]) {
  writedatablock(fileDescriptor, keyphrase, valueArray,
      nItems, datatype, iotype);
}

void posix_constructname(
    const char* in,
    char* out) {
  std::string fullname(in);
  std::string gname("geombc");
  if( fullname.find(gname) != std::string::npos )
    fullname.append(".dat");
  fullname.append(".");
  sprintf(out, "%s", fullname.c_str());
}
