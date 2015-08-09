#ifndef PHSOLVER_PHIO_SYNC_H
#define PHSOLVER_PHIO_SYNC_H
typedef struct phio_file* phio_fp;
struct syncio_file {
  phio_fp base;
  int nfiles;
  int nfields;
  int nppf;
};
typedef struct syncio_file* sync_fp;
void sync_openfile_read(
    const char filename[],
    phio_fp fileDescriptor);
void sync_openfile_write(
    const char filename[],
    phio_fp fileDescriptor);
void sync_closefile(
    phio_fp fileDescriptor);
void sync_readheader(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] );
void sync_writeheader( 
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const int* ndataItems,
    const char datatype[],
    const char iotype[] );
void sync_readdatablock(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] );
void sync_writedatablock(
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[]);
#endif
