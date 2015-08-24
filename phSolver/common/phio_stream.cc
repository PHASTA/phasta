#include <string>
#include <sstream>
#include <phastaIO.h>
#include "phio_stream.h"

#define PHIO_STREAM_TRACING 0
namespace {
  std::string appendPosix(const char* phrase) {
    std::stringstream ss;
    ss << phrase << "?";
    return ss.str();
  }
  void traceEnter(const char* key, const char* aux="") {
    if(PHIO_STREAM_TRACING)
      fprintf(stderr, "CAKE entering %s %s\n", key, aux);
  }
  void traceExit(const char* key) {
    if(PHIO_STREAM_TRACING)
      fprintf(stderr, "CAKE exiting %s\n", key);
  }
}

void stream_openfile(
    const char filename[],
    phio_fp f) {
  traceEnter(__func__, filename);
  stream_fp sf = (stream_fp) f;
  if(sf->mode == 'w' && sf->rs != NULL)
    sf->file = (int*) openRStreamWrite(sf->rs);
  else if(sf->mode == 'r' && sf->grs != NULL)
    sf->file = (int*) openGRStreamRead(sf->grs, filename);
  else
    fprintf(stderr,
        "ERROR %s type of stream %s is unknown... exiting\n",
        __func__, filename);
  traceExit(__func__);
}

void stream_readheader(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  traceEnter(__func__, keyphrase);
  std::string posixPhrase = appendPosix(keyphrase);
  readHeader((FILE*)fileDescriptor, posixPhrase.c_str(),
      (int*)valueArray, *nItems, iotype);
  traceExit(__func__);
}

void stream_writeheader(
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const int* ndataItems,
    const char datatype[],
    const char iotype[] ) {
  traceEnter(__func__, keyphrase);
  std::string posixPhrase = appendPosix(keyphrase);
  writeHeader((FILE*)fileDescriptor, posixPhrase.c_str(),
      (int*)valueArray, *nItems, *ndataItems, datatype, iotype);
  traceExit(__func__);
}

void stream_readdatablock(
    int* fileDescriptor,
    const  char*,
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  traceEnter(__func__);
  readDataBlock((FILE*)fileDescriptor, valueArray, *nItems,
      datatype, iotype);
  traceExit(__func__);
}

void stream_writedatablock(
    const int* fileDescriptor,
    const char*,
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[]) {
  traceEnter(__func__);
  writeDataBlock((FILE*)fileDescriptor, valueArray,
      *nItems, datatype, iotype);
  traceExit(__func__);
}

void stream_closefile(phio_fp f) {
  traceEnter(__func__);
  stream_fp sf = (stream_fp) f;
  fclose((FILE*)sf->file);
  traceExit(__func__);
}

void stream_constructname(const char* in, char* out) {
  traceEnter(__func__);
  sprintf(out, "%s", in); 
  traceExit(__func__);
}
