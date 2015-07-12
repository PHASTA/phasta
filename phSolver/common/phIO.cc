#include "phIO.h"
#include "phComm.h"
#include <stdio.h>
#include <string.h>
#include <phastaIO.h>
#include <sstream>
#include <string>

namespace {
  std::string appendSync(const char* phrase) {
    std::stringstream ss;
    ss << phrase << "@" << phcomm_rank()+1 << "?";
    std::string s = ss.str();
    return s;
  }
  std::string appendColor(const char* phrase, int numFiles) {
    const int color = computeColor(phcomm_rank(), phcomm_size(), numFiles);
    std::stringstream ss;
    ss << phrase << color+1;
    std::string s = ss.str();
    return s;
  }
}

void phio_readheader(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  std::string syncPhrase = appendSync(keyphrase);
  readheader(fileDescriptor, syncPhrase.c_str(),
      valueArray, nItems, datatype, iotype);
}

void phio_readdatablock(
    int*  fileDescriptor,
    const char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  std::string syncPhrase = appendSync(keyphrase);
  readdatablock(fileDescriptor, syncPhrase.c_str(),
      valueArray, nItems, datatype, iotype);
}

void phio_openfile(
    const char filename[],
    const char mode[],
    int* numFiles,
    int* fileDescriptor) {
  std::string syncName = appendColor(filename, *numFiles);
  fprintf(stderr, "filename %s syncName %s\n", filename, syncName.c_str());
  int nfields, nppf;
  queryphmpiio(syncName.c_str(), &nfields, &nppf);
  fprintf(stderr, "nfields %d nppf %d\n", nfields, nppf);
  initphmpiio(&nfields, &nppf, numFiles, fileDescriptor, mode); 
  openfile(syncName.c_str(), mode, fileDescriptor);
}

void phio_restartname(int* step, char* filename) {
  std::stringstream ss;
  ss << "restart-dat." << *step << '.';
  std::string s = ss.str();
  strcpy(filename, s.c_str()); 
}
