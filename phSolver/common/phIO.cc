#include "phIO.h"
#include "phComm.h"
#include <stdio.h>
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
}

void phio_readheader( int* fileDescriptor,
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


