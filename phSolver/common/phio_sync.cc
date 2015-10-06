#include "phIO.h"
#include "phio_sync.h"
#include "phComm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <phastaIO.h>
#include <sstream>
#include <string>

namespace {
  void appendRank(std::stringstream& ss, const char* phrase) {
    ss << phrase << "@" << phcomm_rank()+1;
  }
  std::string appendSync(const char* phrase) {
    std::stringstream ss;
    appendRank(ss,phrase);
    ss << "?";
    return ss.str();
  }
  std::string appendSyncWrite(const char* phrase) {
    std::stringstream ss;
    appendRank(ss,phrase);
    return ss.str();
  }
  std::string appendColor(const char* phrase, int numFiles) {
    const int color = computeColor(phcomm_rank(), phcomm_size(), numFiles);
    std::stringstream ss;
    ss << phrase << color+1;
    return ss.str();
  }
  void close(sync_fp f, const char* mode) {
    int* file = f->file;
    closefile(file, mode);
    finalizephmpiio(file);
    free(file);
    free(f);
  }
}

void sync_openfile_read(
    const char filename[],
    phio_fp f) {
  sync_fp sf = (sync_fp) f;
  std::string syncName = appendColor(filename, sf->nfiles);
  int nfields=0;
  int nppf=0;
  queryphmpiio(syncName.c_str(), &nfields, &nppf);
  const char* mode = "read";
  int* file = sf->file;
  initphmpiio(&nfields, &nppf, &(sf->nfiles), file, mode); 
  openfile(syncName.c_str(), mode, file);
}

void sync_openfile_write(
    const char filename[],
    phio_fp f) {
  sync_fp sf = (sync_fp) f;
  std::string syncName = appendColor(filename, sf->nfiles);
  const char* mode = "write";
  int* file = sf->file;
  initphmpiio(&(sf->nfields), &(sf->nppf),
      &(sf->nfiles), file, mode); 
  openfile(syncName.c_str(), mode, file);
}

void sync_closefile(phio_fp f) {
  sync_fp sf = (sync_fp) f;
  const char m = sf->mode;
  if(m == 'r')
    close(sf, "read");
  else if(m == 'w')
    close(sf, "write");
  else {
    fprintf(stderr, "ERROR unsupported file mode in %s on line %d"
        "... exiting", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

void sync_readheader(
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

void sync_writeheader(
      const int* fileDescriptor,
      const char keyphrase[],
      const void* valueArray,
      const int* nItems,
      const int* ndataItems,
      const char datatype[],
      const char iotype[] ) {
  std::string syncPhrase = appendSyncWrite(keyphrase);
  writeheader(fileDescriptor, syncPhrase.c_str(),
      valueArray, nItems, ndataItems, datatype, iotype);
}

void sync_readdatablock(
    int* fileDescriptor,
    const char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] ) {
  std::string syncPhrase = appendSync(keyphrase);
  readdatablock(fileDescriptor, syncPhrase.c_str(),
      valueArray, nItems, datatype, iotype);
}

void sync_writedatablock(
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[]) {
  std::string syncPhrase = appendSyncWrite(keyphrase);
  writedatablock(fileDescriptor, syncPhrase.c_str(),
      valueArray, nItems, datatype, iotype);
}

void sync_constructname(
    const char* in,
    char* out) {
  std::string fullname(in);
  fullname.append("-dat.");
  sprintf(out, "%s", fullname.c_str());
}
