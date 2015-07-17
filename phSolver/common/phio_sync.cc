#include "phIO.h"
#include "phio_base.h"
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
  void close(phio_fp f, const char* mode) {
    closefile(f->file, mode);
    finalizephmpiio(f->file);
    free(f->file);
    free(f);
  }
}

static struct phio_ops sync_ops = {
  sync_readheader,
  sync_writeheader,
  sync_readdatablock,
  sync_writedatablock,
  sync_restartname,
  sync_closefile_read,
  sync_closefile_write
};

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

void sync_openfile_read(
    const char filename[],
    int* numFiles,
    phio_fp* fileDescriptor) {
  *fileDescriptor =
    (struct phio_file*) malloc(sizeof(struct phio_file));
  (*fileDescriptor)->ops = &sync_ops; 
  (*fileDescriptor)->file = (int*) malloc(sizeof(int*));
  std::string syncName = appendColor(filename, *numFiles);
  int nfields=0;
  int nppf=0;
  queryphmpiio(syncName.c_str(), &nfields, &nppf);
  const char* mode = "read";
  initphmpiio(&nfields, &nppf, numFiles, (*fileDescriptor)->file, mode); 
  openfile(syncName.c_str(), mode, (*fileDescriptor)->file);
}

void sync_openfile_write(
    const char filename[],
    int* numFiles,
    int* numFields,
    int* numPPF,
    phio_fp* fileDescriptor) {
  *fileDescriptor =
    (struct phio_file*) malloc(sizeof(struct phio_file));
  (*fileDescriptor)->ops = &sync_ops; 
  (*fileDescriptor)->file = (int*) malloc(sizeof(int*));
  std::string syncName = appendColor(filename, *numFiles);
  //TODO - define a good upper bound
  assert(*numFields > 0 && *numFields < 1024);
  assert(*numPPF > 0 && *numPPF < 1024);
  const char* mode = "write";
  initphmpiio(numFields, numPPF, numFiles, (*fileDescriptor)->file, mode); 
  openfile(syncName.c_str(), mode, (*fileDescriptor)->file);
}

void sync_restartname(int* step, char* filename) {
  std::stringstream ss;
  ss << "restart-dat." << *step << '.';
  std::string s = ss.str();
  strcpy(filename, s.c_str()); 
}

void sync_closefile_read(phio_fp f) {
  close(f, "read");
}

void sync_closefile_write(phio_fp f) {
  close(f, "write");
}
