#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <sstream>
#include "phIO.h"
#include "phComm.h"
#include "phio_base.h"

#define PHIO_TRACING 0
namespace {
  void trace(const char* key, const char* aux="", void* obj=NULL) {
    if(PHIO_TRACING)
      fprintf(stderr, "PHIO_TRACE entering %s %s %p\n", key, aux, obj);
  }
}

#ifdef __cplusplus
extern "C" {
#endif

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
    phio_fp f,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const int* ndataItems,
    const char datatype[],
    const char iotype[] ) {
  f->ops->writeheader(f->file, keyphrase, valueArray,
      nItems, ndataItems, datatype, iotype);
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
    phio_fp f,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[]) {
  f->ops->writedatablock(f->file, keyphrase, valueArray,
      nItems, datatype, iotype);
}

void phio_constructName(
    phio_fp f,
    const char inName[],
    char* outName) {
  f->ops->constructname(inName, outName);
}

void phio_openfile(
    const char filename[],
    phio_fp f) {
  trace("openfile",filename,f);
  f->ops->openfile(filename, f);
}

void phio_closefile(phio_fp f) {
  trace("closefile","unknown",f);
  f->ops->closefile(f);
}

void phio_appendInt(char* dest, int v) {
  std::stringstream ss;
  ss << dest << v << '.';
  std::string s = ss.str();
  strcpy(dest, s.c_str());
}

#ifdef __cplusplus
}
#endif
