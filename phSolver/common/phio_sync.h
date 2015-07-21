#ifndef PHSOLVER_PHIO_SYNC_H
#define PHSOLVER_PHIO_SYNC_H
typedef struct phio_file* phio_fp;
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
void sync_openfile_read(
    const char filename[],
    int* numFiles,
    phio_fp* fileDescriptor);
void sync_openfile_write(
    const char filename[],
    int* numFiles,
    int* numFields,
    int* numPPF,
    phio_fp* fileDescriptor);
void sync_closefile_read(phio_fp fileDescriptor);
void sync_closefile_write(phio_fp fileDescriptor);
#endif
