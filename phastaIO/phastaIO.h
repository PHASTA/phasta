/*  Primary interface for the Phasta Binary read and write routines these*/
/*  functions are 'C' callable.( All arguments have been kept as pointers to*/
/*  facilitate calling from Fortran )*/
/*  Anil Kumar Karanam  Spring 2003*/
#ifndef _PHASTAIO_H_
#define _PHASTAIO_H_

#include <FCMangle.h>
#include <mpi.h>

#ifdef intel
#define ios_base ios
#endif

#define queryphmpiio FortranCInterface_GLOBAL_(queryphmpiio,QUERYPHMPIIO)
#define initphmpiio FortranCInterface_GLOBAL_(initphmpiio,INITPHMPIIO)
#define initphmpiiosub FortranCInterface_GLOBAL_(initphmpiiosub,INITPHMPIIOSUB)
#define finalizephmpiio FortranCInterface_GLOBAL_(finalizephmpiio,FINALIZEPHMPIIO)

#define openfile FortranCInterface_GLOBAL_(openfile, OPENFILE)
#define closefile FortranCInterface_GLOBAL_(closefile, CLOSEFILE)
#define readheader FortranCInterface_GLOBAL_(readheader, READHEADER)
#define readdatablock FortranCInterface_GLOBAL_(readdatablock, READDATABLOCK)
#define writeheader FortranCInterface_GLOBAL_(writeheader, WRITEHEADER)
#define writedatablock FortranCInterface_GLOBAL_(writedatablock, WRITEDATABLOCK)
#define writestring FortranCInterface_GLOBAL_(writestring, WRITESTRING)
#define togglestrictmode FortranCInterface_GLOBAL_(togglestrictmode, TOGGLESTRICTMODE)
#define SwapArrayByteOrder FortranCInterface_GLOBAL_(swaparraybyteorder, SWAPARRAYBYTEORDER)
#define isLittleEndian FortranCInterface_GLOBAL_(islittleendian, ISLITTLEENDIAN)


#if defined (__cplusplus)
extern "C" {
#endif

  void mem_alloc( void* p, char* type, int size );

  void
  queryphmpiio( const char filename[],
		 int *nfields,
		 int *nppf );

  int
  initphmpiio( int *nfields,
		int *nppf,
		int *nfiles,
		int *filehandle,
		const char mode[] );
  int
  initphmpiiosub( int *nfields,
		  int *nppf,
		  int *nfiles,
		  int *filehandle,
		  const char mode[],
                  MPI_Comm my_local_comm );

  void
  finalizephmpiio( int *fileDescriptor );

    void
    SwapArrayByteOrder( void* array,
                         int   nbytes,
                         int   nItems ) ;
    void

    openfile( const char filename[],
               const char mode[],
               int* fileDescriptor );

    void
    closefile( int* fileDescriptor,
                const char mode[] );

    void
    readheader( int*   fileDescriptor,
                 const char  keyphrase[],
                 void*  valueArray,
                 int*   nItems,
                 const char   datatype[],
                 const char   iotype[] );

    void
    writeheader( const int*  fileDescriptor,
                  const char keyphrase[],
                  const void* valueArray,
                  const int*  nItems,
                  const int*  ndataItems,
                  const char  datatype[],
                  const char  iotype[] );

    void
    readdatablock( int*  fileDescriptor,
                    const char keyphrase[],
                    void* valueArray,
                    int*  nItems,
                    const char  datatype[],
                    const char  iotype[] );


    void
    writedatablock( const int*   fileDescriptor,
                     const char  keyphrase[],
                     const void*  valueArray,
                     const int*   nItems,
                     const char   datatype[],
                     const char   iotype[]  );

    void
    writestring( int* fileDescriptor,
                  const char inString[] );

    void
    togglestrictmode( );

  int
  isLittleEndian( ) ;

  int computeColor( int myrank, int numprocs, int nfiles);

#ifdef __cplusplus
} // end of extern "C".

#include <string>

namespace PHASTA {
template<class T>
struct PhastaIO_traits;

template<>
struct PhastaIO_traits<int> {
  static const char* const type_string;
};


template<>
struct PhastaIO_traits<double> {
  static const char* const type_string;
};


template<class T>
void
write_data_block( const int   fileDescriptor,
		  const std::string keyphrase,
		  const T* const  valueArray,
		  const int   nItems,
		  const bool inBinary = false) {
  writedatablock(&fileDescriptor,
		  keyphrase.c_str(),
		  reinterpret_cast<const void*>(valueArray),
		  &nItems,
		  PhastaIO_traits<T>::type_string,
		  inBinary ? "binary" : "text");

}
template<class T>
void
write_header( const int  fileDescriptor,
		 const std::string& keyphrase,
		 const T* valueArray,
		 const int  nItems,
		 const int  nDataItems,
		 const bool inBinary = false) {
      writeheader(&fileDescriptor,
		   keyphrase.c_str(),
		   reinterpret_cast<const void*>(valueArray),
		   &nItems,
		   &nDataItems,
		   PhastaIO_traits<T>::type_string,
		   inBinary ? "binary" : "text");
}


} // namespace PHASTA

#endif // __cplusplus

#endif // _PHASTAIO_H_
