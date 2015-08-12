#ifndef PHSOLVER_PHIO_BASE_H
#define PHSOLVER_PHIO_BASE_H
struct phio_file {
  struct phio_ops const* ops;
  int* file;
  char mode;
};
typedef struct phio_file* phio_fp;
struct phio_ops {
  void (*openfile)(
      const  char keyphrase[],
      phio_fp fileDescriptor);
  void (*closefile)(
      phio_fp fileDescriptor);
  void (*readheader)(
      int* fileDescriptor,
      const  char keyphrase[],
      void* valueArray,
      int*  nItems,
      const char  datatype[],
      const char  iotype[] );
  void (*writeheader)( 
      const int* fileDescriptor,
      const char keyphrase[],
      const void* valueArray,
      const int* nItems,
      const int* ndataItems,
      const char datatype[],
      const char iotype[] );
  void (*readdatablock)(
      int* fileDescriptor,
      const  char keyphrase[],
      void* valueArray,
      int*  nItems,
      const char  datatype[],
      const char  iotype[] );
  void (*writedatablock)(
      const int* fileDescriptor,
      const char keyphrase[],
      const void* valueArray,
      const int* nItems,
      const char datatype[],
      const char iotype[]);
  void (*constructname)(
      const char inName[],
      char* outName);
};

#endif
