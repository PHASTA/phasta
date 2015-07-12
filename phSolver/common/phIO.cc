#include "phIO.h"
#include "phComm.h"
#include <phastaIO.h>
#include <sstream>
#include <string>

void phio_readheader( int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  std::stringstream ss;
  ss << keyphrase << "@" << phcomm_rank()+1 << "?";
  std::string s = ss.str();
  readheader(fileDescriptor, s.c_str(),
      valueArray, nItems, datatype, iotype);
}

