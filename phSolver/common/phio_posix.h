#ifndef PHSOLVER_PHIO_POSIX_H
#define PHSOLVER_PHIO_POSIX_H
typedef struct phio_file* phio_fp;
void posix_openfile(
    const char filename[],
    phio_fp fileDescriptor);
int posix_openfile_single(
    const char filename[],
    phio_fp fileDescriptor);
void posix_closefile(
    phio_fp fileDescriptor);
void posix_readheader(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] );
void posix_writeheader( 
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const int* ndataItems,
    const char datatype[],
    const char iotype[] );
void posix_readdatablock(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] );
void posix_writedatablock(
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[]);
void posix_constructname(
    const char* in,
    char* out);
#endif
