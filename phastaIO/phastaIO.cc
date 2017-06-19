/*
 *
 * This is the SyncIO library that uses MPI-IO collective functions to
 * implement a flexible I/O checkpoint solution for a large number of
 * processors.
 *
 * Previous developer: Ning Liu         (liun2@cs.rpi.edu)
 *                     Jing Fu          (fuj@cs.rpi.edu)
 * Current developers: Michel Rasquin   (Michel.Rasquin@colorado.edu),
 *                     Ben Matthews     (benjamin.a.matthews@colorado.edu)
 *
 */

#include <map>
#include <vector>
#include <string>
#include <cstdarg>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include "phastaIO.h"
#include "phiotmrc.h"
#include "phiompi.h"
#include "mpi.h"
#include "phiostats.h"
#include "phiotimer.h"
#include <assert.h>

/* OS-specific things try to stay here */
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>


#define VERSION_INFO_HEADER_SIZE 8192
#define DB_HEADER_SIZE 1024
#define ONE_MEGABYTE 1048576
#define TWO_MEGABYTE 2097152
#define ENDIAN_TEST_NUMBER 12180 // Troy's Zip Code!!
#define MAX_PHASTA_FILES 64
#define MAX_PHASTA_FILE_NAME_LENGTH 1024
#define MAX_FIELDS_NUMBER ((VERSION_INFO_HEADER_SIZE/MAX_FIELDS_NAME_LENGTH)-4) // The meta data include - MPI_IO_Tag, nFields, nFields*names of the fields, nppf
// -3 for MPI_IO_Tag, nFields and nppf, -4 for extra security (former nFiles)
#define MAX_FIELDS_NAME_LENGTH 128
#define DefaultMHSize (4*ONE_MEGABYTE)
//#define DefaultMHSize (8350) //For test
#define LATEST_WRITE_VERSION 1
#define inv1024sq 953.674316406e-9 // = 1/1024/1024
int MasterHeaderSize = -1;

bool PRINT_PERF = false; // default print no perf results
int irank = -1; // global rank, should never be manually manipulated
int mysize = -1;

// Static variables are bad but used here to store the subcommunicator and associated variables
// Prevent MPI_Comm_split to be called more than once, especially on BGQ with the V1R2M1 driver (leak detected in MPI_Comm_split - IBM working on it)
static int s_assign_local_comm = 0;
static MPI_Comm s_local_comm;
static int s_local_size = -1;
static int s_local_rank = -1;

// the following defines are for debug printf
#define PHASTAIO_DEBUG 0 //default to not print any debugging info

void phprintf(const char* fmt, ...) {
  (void)fmt;
#if PHASTAIO_DEBUG
  char format[1024];
  snprintf(format, sizeof(format), "phastaIO debug: %s", fmt);
  va_list ap;
  va_start(ap,fmt);
  vprintf(format,ap);
  va_end(ap);
#endif
}

void phprintf_0(const char* fmt, ...) {
  (void)fmt;
#if PHASTAIO_DEBUG
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0){
    char format[1024];
    snprintf(format, sizeof(format), "phastaIO debug: irank=0 %s", fmt);
    va_list ap;
    va_start(ap,fmt);
    vprintf(format, ap);
    va_end(ap);
  }
#endif
}

enum PhastaIO_Errors
{
  MAX_PHASTA_FILES_EXCEEDED = -1,
  UNABLE_TO_OPEN_FILE = -2,
  NOT_A_MPI_FILE = -3,
  GPID_EXCEEDED = -4,
  DATA_TYPE_ILLEGAL = -5
};

using namespace std;

namespace{

  map<int, std::string> LastHeaderKey;
  vector< FILE* > fileArray;
  vector< bool > byte_order;
  vector< int > header_type;
  int DataSize=0;
  bool LastHeaderNotFound = false;
  bool Wrong_Endian = false ;
  bool Strict_Error = false ;
  bool binary_format = true;

  /***********************************************************************/
  /***************** NEW PHASTA IO CODE STARTS HERE **********************/
  /***********************************************************************/

  typedef struct
  {
    char filename[MAX_PHASTA_FILE_NAME_LENGTH];   /* defafults to 1024 */
    unsigned long my_offset;
    unsigned long next_start_address;
    unsigned long **my_offset_table;
    unsigned long **my_read_table;

    double * double_chunk;
    double * read_double_chunk;

    int field_count;
    int part_count;
    int read_field_count;
    int read_part_count;
    int GPid;
    int start_id;

    int mhsize;

    int myrank;
    int numprocs;
    int local_myrank;
    int local_numprocs;

    int nppp;
    int nPPF;
    int nFiles;
    int nFields;

    int * int_chunk;
    int * read_int_chunk;

    int Wrong_Endian; /* default to false */
    char * master_header;
    MPI_File file_handle;
    MPI_Comm local_comm;
  } phastaio_file_t;

  typedef struct
  {
    int nppf, nfields;
    char * masterHeader;
  }serial_file;

  serial_file *SerialFile;
  phastaio_file_t *PhastaIOActiveFiles[MAX_PHASTA_FILES];
  int PhastaIONextActiveIndex = 0; /* indicates next index to allocate */

  // the caller has the responsibility to delete the returned string
  // TODO: StringStipper("nbc value? ") returns NULL?
  char* StringStripper( const char  istring[] ) {
    int length = strlen( istring );
    char* dest = (char *)malloc( length + 1 );
    strcpy( dest, istring );
    dest[ length ] = '\0';
    if ( char* p = strpbrk( dest, " ") )
      *p = '\0';
    return dest;
  }

  inline int cscompare( const char teststring[], const char targetstring[] ) {
    char* s1 = const_cast<char*>(teststring);
    char* s2 = const_cast<char*>(targetstring);

    while( *s1 == ' ') s1++;
    while( *s2 == ' ') s2++;
    while( ( *s1 )
        && ( *s2 )
        && ( *s2 != '?')
        && ( tolower( *s1 )==tolower( *s2 ) ) ) {
      s1++;
      s2++;
      while( *s1 == ' ') s1++;
      while( *s2 == ' ') s2++;
    }
    if ( !( *s1 ) || ( *s1 == '?') ) return 1;
    else return 0;
  }

  inline void isBinary( const char iotype[] ) {
    char* fname = StringStripper( iotype );
    if ( cscompare( fname, "binary" ) ) binary_format = true;
    else binary_format = false;
    free (fname);

  }

  inline size_t typeSize( const char typestring[] ) {
    char* ts1 = StringStripper( typestring );
    if ( cscompare( "integer", ts1 ) ) {
      free (ts1);
      return sizeof(int);
    } else if ( cscompare( "double", ts1 ) ) {
      free (ts1);
      return sizeof( double );
    } else {
      free (ts1);
      fprintf(stderr,"unknown type : %s\n",ts1);
      return 0;
    }
  }

  int readHeader(
      FILE*       fileObject,
        const char  phrase[],
        int*        params,
        int         expect ) {
    char* text_header;
    char* token;
    char Line[1024] = "\0";
    char junk;
    bool FOUND = false ;
    int real_length;
    int skip_size, integer_value;
    int rewind_count=0;

    if( !fgets( Line, 1024, fileObject ) && feof( fileObject ) ) {
      rewind( fileObject );
      clearerr( fileObject );
      rewind_count++;
      fgets( Line, 1024, fileObject );
    }

    while( !FOUND  && ( rewind_count < 2 ) )  {
      if ( ( Line[0] != '\n' ) && ( real_length = strcspn( Line, "#" )) ) {
        text_header = (char *)malloc( real_length + 1 );
        strncpy( text_header, Line, real_length );
        text_header[ real_length ] =static_cast<char>(NULL);
        token = strtok ( text_header, ":" );
        assert(token);
        if( cscompare( phrase , token ) ) {
          FOUND = true ;
          token = strtok( NULL, " ,;<>" );
          assert(token);
          skip_size = atoi( token );
          int i;
          for( i=0; i < expect && ( token = strtok( NULL," ,;<>") ); i++) {
            assert(token);
            params[i] = atoi( token );
          }
          if ( i < expect ) {
            fprintf(stderr,
                "Aloha Expected # of ints not found for: %s\n",phrase );
          }
        } else if ( cscompare(token,"byteorder magic number") ) {
          if ( binary_format ) {
            fread((void*)&integer_value,sizeof(int),1,fileObject);
            fread( &junk, sizeof(char), 1 , fileObject );
            if ( 362436 != integer_value ) Wrong_Endian = true;
          } else{
            fscanf(fileObject, "%d\n", &integer_value );
          }
        } else {
          /* some other header, so just skip over */
          token = strtok( NULL, " ,;<>" );
          assert(token);
          skip_size = atoi( token );
          if ( binary_format)
            fseek( fileObject, skip_size, SEEK_CUR );
          else
            for( int gama=0; gama < skip_size; gama++ )
              fgets( Line, 1024, fileObject );
        }
        free (text_header);
      } // end of if before while loop

      if ( !FOUND )
        if( !fgets( Line, 1024, fileObject ) && feof( fileObject ) ) {
          rewind( fileObject );
          clearerr( fileObject );
          rewind_count++;
          fgets( Line, 1024, fileObject );
        }
    }

    if ( !FOUND ) {
      //fprintf(stderr, "Error: Could not find: %s\n", phrase);
      if(irank == 0) printf("WARNING: Could not find: %s\n", phrase);
      return 1;
    }
    return 0;
  }

} // end unnamed namespace


// begin of publicly visible functions

/**
 * This function takes a long long pointer and assign (start) phiotmrc value to it
 */
void startTimer(double* start) {
  if( !PRINT_PERF ) return;
  MPI_Barrier(MPI_COMM_WORLD);
  *start =  phiotmrc();
}

/**
 * This function takes a long long pointer and assign (end) phiotmrc value to it
 */
void endTimer(double* end) {
  if( !PRINT_PERF ) return;
  *end = phiotmrc();
  MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * choose to print some performance results (or not) according to
 * the PRINT_PERF macro
 */
void printPerf(
    const char* func_name,
    double start,
    double end,
    unsigned long datasize,
    int printdatainfo,
    const char* extra_msg) {
  if( !PRINT_PERF ) return;
  unsigned long data_size = datasize;
  double time = end - start;
  unsigned long isizemin,isizemax,isizetot;
  double sizemin,sizemax,sizeavg,sizetot,rate;
  double tmin, tmax, tavg, ttot;

  MPI_Allreduce(&time, &tmin,1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&time, &tmax,1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&time, &ttot,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  tavg = ttot/mysize;

  if(irank == 0) {
    if ( PhastaIONextActiveIndex == 0 ) printf("** 1PFPP ");
    else  printf("** syncIO ");
    printf("%s(): Tmax = %f sec, Tmin = %f sec, Tavg = %f sec", func_name, tmax, tmin, tavg);
  }

  if(printdatainfo == 1) { // if printdatainfo ==1, compute I/O rate and block size
    MPI_Allreduce(&data_size,&isizemin,1,MPI_LONG_LONG_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&data_size,&isizemax,1,MPI_LONG_LONG_INT,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&data_size,&isizetot,1,MPI_LONG_LONG_INT,MPI_SUM,MPI_COMM_WORLD);

    sizemin=(double)(isizemin*inv1024sq);
    sizemax=(double)(isizemax*inv1024sq);
    sizetot=(double)(isizetot*inv1024sq);
    sizeavg=(double)(1.0*sizetot/mysize);
    rate=(double)(1.0*sizetot/tmax);

    if( irank == 0) {
      printf(", Rate = %f MB/s [%s] \n \t\t\t"
             " block size: Min= %f MB; Max= %f MB; Avg= %f MB; Tot= %f MB\n",
             rate, extra_msg, sizemin,sizemax,sizeavg,sizetot);
    }
  }
  else {
    if(irank == 0) {
      printf(" \n");
      //printf(" (%s) \n", extra_msg);
    }
  }
}

/**
 * This function is normally called at the beginning of a read operation, before
 * init function.
 * This function (uses rank 0) reads out nfields, nppf, master header size,
 * endianess and allocates for masterHeader string.
 * These values are essential for following read operations. Rank 0 will bcast
 * these values to other ranks in the commm world
 *
 * If the file set is of old POSIX format, it would throw error and exit
 */
void queryphmpiio(const char filename[],int *nfields, int *nppf)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);
  MPI_Comm_size(MPI_COMM_WORLD, &mysize);

  if(irank == 0) {
    FILE * fileHandle;
    char* fname = StringStripper( filename );

    PHASTAIO_OPENTIME(fileHandle = fopen (fname,"rb");)
    if (fileHandle == NULL ) {
      printf("\nError: File %s doesn't exist! Please check!\n",fname);
    }
    else {
      SerialFile =(serial_file *)calloc( 1,  sizeof( serial_file) );
      int meta_size_limit = VERSION_INFO_HEADER_SIZE;
      SerialFile->masterHeader = (char *)malloc( meta_size_limit );
      fread(SerialFile->masterHeader, 1, meta_size_limit, fileHandle);

      char read_out_tag[MAX_FIELDS_NAME_LENGTH];
      char version[MAX_FIELDS_NAME_LENGTH/4];
      int mhsize;
      char * token;
      int magic_number;

      memcpy( read_out_tag,
          SerialFile->masterHeader,
          MAX_FIELDS_NAME_LENGTH-1 );

      if ( cscompare ("MPI_IO_Tag",read_out_tag) ) {
        // Test endianess ...
        memcpy (&magic_number,
            SerialFile->masterHeader + sizeof("MPI_IO_Tag : ")-1, //-1 sizeof returns the size of the string+1 for "\0"
            sizeof(int) );                                        // masterheader should look like "MPI_IO_Tag : 12180 " with 12180 in binary format

        if ( magic_number != ENDIAN_TEST_NUMBER ) {
          printf("Endian is different!\n");
          // Will do swap later
        }

        // test version, old version, default masterheader size is 4M
        // newer version masterheader size is read from first line
        memcpy(version,
            SerialFile->masterHeader + MAX_FIELDS_NAME_LENGTH/2,
            MAX_FIELDS_NAME_LENGTH/4 - 1); //TODO: why -1?

        if( cscompare ("version",version) ) {
          // if there is "version" tag in the file, then it is newer format
          // read master header size from here, otherwise use default
          // Note: if version is "1", we know mhsize is at 3/4 place...

          token = strtok(version, ":");
          token = strtok(NULL, " ,;<>" );
          int iversion = atoi(token);

          if( iversion == 1) {
            memcpy( &mhsize,
                SerialFile->masterHeader + MAX_FIELDS_NAME_LENGTH/4*3 + sizeof("mhsize : ")-1,
                sizeof(int));
            if ( magic_number != ENDIAN_TEST_NUMBER )
              SwapArrayByteOrder(&mhsize, sizeof(int), 1);

            if( mhsize > DefaultMHSize ) {
              //if actual headersize is larger than default, let's re-read
              free(SerialFile->masterHeader);
              SerialFile->masterHeader = (char *)malloc(mhsize);
              fseek(fileHandle, 0, SEEK_SET); // reset the file stream position
              fread(SerialFile->masterHeader,1,mhsize,fileHandle);
            }
          }
          //TODO: check if this is a valid int??
          MasterHeaderSize = mhsize;
        }
        else { // else it's version 0's format w/o version tag, implicating MHSize=4M
          MasterHeaderSize = DefaultMHSize;
        }

        memcpy( read_out_tag,
            SerialFile->masterHeader+MAX_FIELDS_NAME_LENGTH+1,
            MAX_FIELDS_NAME_LENGTH ); //TODO: why +1

        // Read in # fields ...
        token = strtok( read_out_tag, ":" );
        token = strtok( NULL," ,;<>" );
        *nfields = atoi( token );
        if ( *nfields > MAX_FIELDS_NUMBER) {
          printf("Error queryphmpiio: nfields is larger than MAX_FIELDS_NUMBER!\n");
        }
        SerialFile->nfields=*nfields; //TODO: sanity check of this int?

        memcpy( read_out_tag,
            SerialFile->masterHeader + MAX_FIELDS_NAME_LENGTH * 2
            + *nfields * MAX_FIELDS_NAME_LENGTH,
            MAX_FIELDS_NAME_LENGTH);

        token = strtok( read_out_tag, ":" );
        token = strtok( NULL," ,;<>" );
        *nppf = atoi( token );
        SerialFile->nppf=*nppf; //TODO: sanity check of int
      } // end of if("MPI_IO_TAG")
      else {
        printf("Error queryphmpiio: The file you opened is not of syncIO new format, please check! read_out_tag = %s\n",read_out_tag);
        exit(1);
      }
      PHASTAIO_CLOSETIME(fclose(fileHandle);)
      free(SerialFile->masterHeader);
      free(SerialFile);
    } //end of else
    free(fname);
  }

  // Bcast value to every one
  MPI_Bcast( nfields, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( nppf, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &MasterHeaderSize, 1, MPI_INT, 0, MPI_COMM_WORLD );
  phprintf("Info queryphmpiio: myrank = %d, MasterHeaderSize = %d\n", irank, MasterHeaderSize);
}

/**
 * This function computes the right master header size (round to size of 2^n).
 * This is only needed for file format version 1 in "write" mode.
 */
int computeMHSize(int nfields, int nppf, int version) {
  int mhsize=0;
  if(version == 1) {
    //int meta_info_size = (2+nfields+1) * MAX_FIELDS_NAME_LENGTH; // 2 is MPI_IO_TAG and nFields, the others 1 is nppf
    int meta_info_size = VERSION_INFO_HEADER_SIZE;
    int actual_size =  nfields * nppf * sizeof(unsigned long) + meta_info_size;
    //printf("actual_size = %d, offset table size = %d\n", actual_size,  nfields * nppf * sizeof(long long));
    if (actual_size > DefaultMHSize) {
      mhsize = (int) ceil( (double) actual_size/DefaultMHSize); // it's rounded to ceiling of this value
      mhsize *= DefaultMHSize;
    }
    else {
      mhsize = DefaultMHSize;
    }
  } else {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(!rank) {
      fprintf(stderr,
          "ERROR invalid version passed to %s... exiting\n", __func__);
      exit(EXIT_FAILURE);
    }
  }
  return mhsize;
}

/**
 * Computes correct color of a rank according to number of files.
 */
extern "C" int computeColor( int myrank, int numprocs, int nfiles) {
  int color =
    (int)(myrank / (numprocs / nfiles));
  return color;
}


/**
 * Check the file descriptor.
 */
void checkFileDescriptor(const char fctname[], int*  fileDescriptor ) {
  if ( *fileDescriptor < 0 ) {
    printf("Error: File descriptor = %d in %s\n",*fileDescriptor,fctname);
    exit(1);
  }
}

/**
 * Initialize the file struct members and allocate space for file struct
 * buffers.
 *
 * Note: this function is only called when we are using new format. Old POSIX
 * format should skip this routine and call openfile() directly instead.
 */
int initphmpiio( int *nfields, int *nppf, int *nfiles, int *filehandle, const char mode[])
{
  // we init irank again in case query not called (e.g. syncIO write case)
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);
  MPI_Comm_size(MPI_COMM_WORLD, &mysize);

  phprintf("Info initphmpiio: entering function, myrank = %d, MasterHeaderSize = %d, nfields %d, nppf %d, nfiles %d\n", irank, MasterHeaderSize, *nfields, *nppf, *nfiles);

  double timer_start, timer_end;
  startTimer(&timer_start);

  char* imode = StringStripper( mode );

  // Note: if it's read, we presume query was called prior to init and
  // MasterHeaderSize is already set to correct value from parsing header
  // otherwise it's write then it needs some computation to be set
  if ( cscompare( "read", imode ) ) {
    // do nothing
  }
  else if( cscompare( "write", imode ) ) {
    MasterHeaderSize =  computeMHSize(*nfields, *nppf, LATEST_WRITE_VERSION);
  }
  else {
    printf("Error initphmpiio: can't recognize the mode %s", imode);
    exit(1);
  }
  free ( imode );

  phprintf("Info initphmpiio: myrank = %d, MasterHeaderSize = %d\n", irank, MasterHeaderSize);

  int i, j;

  if( PhastaIONextActiveIndex == MAX_PHASTA_FILES ) {
    printf("Error initphmpiio: PhastaIONextActiveIndex = MAX_PHASTA_FILES");
    endTimer(&timer_end);
    printPerf("initphmpiio", timer_start, timer_end, 0, 0, "");
    return MAX_PHASTA_FILES_EXCEEDED;
  }
  //		else if( PhastaIONextActiveIndex == 0 )  //Hang in debug mode on Intrepid
  //		{
  //			for( i = 0; i < MAX_PHASTA_FILES; i++ );
  //			{
  //				PhastaIOActiveFiles[i] = NULL;
  //			}
  //		}


  PhastaIOActiveFiles[PhastaIONextActiveIndex] = (phastaio_file_t *)calloc( 1,  sizeof( phastaio_file_t) );

  i = PhastaIONextActiveIndex;
  PhastaIONextActiveIndex++;

  //PhastaIOActiveFiles[i]->next_start_address = 2*TWO_MEGABYTE;

  PhastaIOActiveFiles[i]->next_start_address = MasterHeaderSize;  // what does this mean??? TODO

  PhastaIOActiveFiles[i]->Wrong_Endian = false;

  PhastaIOActiveFiles[i]->nFields = *nfields;
  PhastaIOActiveFiles[i]->nPPF = *nppf;
  PhastaIOActiveFiles[i]->nFiles = *nfiles;
  MPI_Comm_rank(MPI_COMM_WORLD, &(PhastaIOActiveFiles[i]->myrank));
  MPI_Comm_size(MPI_COMM_WORLD, &(PhastaIOActiveFiles[i]->numprocs));


  if( *nfiles > 1 ) { // split the ranks according to each mpiio file

    if ( s_assign_local_comm == 0) { // call mpi_comm_split for the first (and only) time

      if (PhastaIOActiveFiles[i]->myrank == 0) printf("Building subcommunicator\n");

      int color = computeColor(PhastaIOActiveFiles[i]->myrank, PhastaIOActiveFiles[i]->numprocs, PhastaIOActiveFiles[i]->nFiles);
      MPI_Comm_split(MPI_COMM_WORLD,
          color,
          PhastaIOActiveFiles[i]->myrank,
          &(PhastaIOActiveFiles[i]->local_comm));
      MPI_Comm_size(PhastaIOActiveFiles[i]->local_comm,
          &(PhastaIOActiveFiles[i]->local_numprocs));
      MPI_Comm_rank(PhastaIOActiveFiles[i]->local_comm,
          &(PhastaIOActiveFiles[i]->local_myrank));

      // back up now these variables so that we do not need to call comm_split again
      s_local_comm = PhastaIOActiveFiles[i]->local_comm;
      s_local_size = PhastaIOActiveFiles[i]->local_numprocs;
      s_local_rank = PhastaIOActiveFiles[i]->local_myrank;
      s_assign_local_comm = 1;
    }
    else { // recycle the subcommunicator
      if (PhastaIOActiveFiles[i]->myrank == 0) printf("Recycling subcommunicator\n");
      PhastaIOActiveFiles[i]->local_comm = s_local_comm;
      PhastaIOActiveFiles[i]->local_numprocs = s_local_size;
      PhastaIOActiveFiles[i]->local_myrank = s_local_rank;
    }
  }
  else { // *nfiles == 1 here - no need to call mpi_comm_split here

    if (PhastaIOActiveFiles[i]->myrank == 0) printf("Bypassing subcommunicator\n");
    PhastaIOActiveFiles[i]->local_comm = MPI_COMM_WORLD;
    PhastaIOActiveFiles[i]->local_numprocs = PhastaIOActiveFiles[i]->numprocs;
    PhastaIOActiveFiles[i]->local_myrank = PhastaIOActiveFiles[i]->myrank;

  }

  PhastaIOActiveFiles[i]->nppp =
    PhastaIOActiveFiles[i]->nPPF/PhastaIOActiveFiles[i]->local_numprocs;

  PhastaIOActiveFiles[i]->start_id = PhastaIOActiveFiles[i]->nPPF *
    (int)(PhastaIOActiveFiles[i]->myrank/PhastaIOActiveFiles[i]->local_numprocs) +
    (PhastaIOActiveFiles[i]->local_myrank * PhastaIOActiveFiles[i]->nppp);

  PhastaIOActiveFiles[i]->my_offset_table =
    ( unsigned long ** ) calloc( MAX_FIELDS_NUMBER , sizeof( unsigned long *) );

  PhastaIOActiveFiles[i]->my_read_table =
    ( unsigned long ** ) calloc( MAX_FIELDS_NUMBER , sizeof( unsigned long *) );

  for (j=0; j<*nfields; j++)
  {
    PhastaIOActiveFiles[i]->my_offset_table[j] =
      ( unsigned long * ) calloc( PhastaIOActiveFiles[i]->nppp , sizeof( unsigned long) );

    PhastaIOActiveFiles[i]->my_read_table[j] =
      ( unsigned long * ) calloc( PhastaIOActiveFiles[i]->nppp , sizeof( unsigned long) );
  }
  *filehandle = i;

  PhastaIOActiveFiles[i]->master_header = (char *)calloc(MasterHeaderSize,sizeof( char ));
  PhastaIOActiveFiles[i]->double_chunk = (double *)calloc(1,sizeof( double ));
  PhastaIOActiveFiles[i]->int_chunk = (int *)calloc(1,sizeof( int ));
  PhastaIOActiveFiles[i]->read_double_chunk = (double *)calloc(1,sizeof( double ));
  PhastaIOActiveFiles[i]->read_int_chunk = (int *)calloc(1,sizeof( int ));

  // Time monitoring
  endTimer(&timer_end);
  printPerf("initphmpiio", timer_start, timer_end, 0, 0, "");

  phprintf_0("Info initphmpiio: quiting function");

  return i;
}

/**
 * Destruct the file struct and free buffers allocated in init function.
 */
void finalizephmpiio( int *fileDescriptor )
{
  double timer_start, timer_end;
  startTimer(&timer_start);

  int i, j;
  i = *fileDescriptor;
  //PhastaIONextActiveIndex--;

  /* //free the offset table for this phasta file */
  //for(j=0; j<MAX_FIELDS_NUMBER; j++) //Danger: undefined behavior for my_*_table.[j] not allocated or not initialized to NULL
  for(j=0; j<PhastaIOActiveFiles[i]->nFields; j++)
  {
    free( PhastaIOActiveFiles[i]->my_offset_table[j]);
    free( PhastaIOActiveFiles[i]->my_read_table[j]);
  }
  free ( PhastaIOActiveFiles[i]->my_offset_table );
  free ( PhastaIOActiveFiles[i]->my_read_table );
  free ( PhastaIOActiveFiles[i]->master_header );
  free ( PhastaIOActiveFiles[i]->double_chunk );
  free ( PhastaIOActiveFiles[i]->int_chunk );
  free ( PhastaIOActiveFiles[i]->read_double_chunk );
  free ( PhastaIOActiveFiles[i]->read_int_chunk );

  if( PhastaIOActiveFiles[i]->nFiles > 1 && s_assign_local_comm ) { // the comm was split
    if (PhastaIOActiveFiles[i]->myrank == 0) printf("Freeing subcommunicator\n");
    s_assign_local_comm = 0;
    MPI_Comm_free(&(PhastaIOActiveFiles[i]->local_comm));
  }

  free( PhastaIOActiveFiles[i]);

  endTimer(&timer_end);
  printPerf("finalizempiio", timer_start, timer_end, 0, 0, "");

  PhastaIONextActiveIndex--;
}


/**
 * Special init for M2N in order to create a subcommunicator for the reduced solution (requires PRINT_PERF to be false for now)
 * Initialize the file struct members and allocate space for file struct buffers.
 *
 * Note: this function is only called when we are using new format. Old POSIX
 * format should skip this routine and call openfile() directly instead.
 */
int initphmpiiosub( int *nfields, int *nppf, int *nfiles, int *filehandle, const char mode[],MPI_Comm my_local_comm)
{
  // we init irank again in case query not called (e.g. syncIO write case)

  MPI_Comm_rank(my_local_comm, &irank);
  MPI_Comm_size(my_local_comm, &mysize);

  phprintf("Info initphmpiio: entering function, myrank = %d, MasterHeaderSize = %d\n", irank, MasterHeaderSize);

  double timer_start, timer_end;
  startTimer(&timer_start);

  char* imode = StringStripper( mode );

  // Note: if it's read, we presume query was called prior to init and
  // MasterHeaderSize is already set to correct value from parsing header
  // otherwise it's write then it needs some computation to be set
  if ( cscompare( "read", imode ) ) {
    // do nothing
  }
  else if( cscompare( "write", imode ) ) {
    MasterHeaderSize =  computeMHSize(*nfields, *nppf, LATEST_WRITE_VERSION);
  }
  else {
    printf("Error initphmpiio: can't recognize the mode %s", imode);
    exit(1);
  }
  free ( imode );

  phprintf("Info initphmpiio: myrank = %d, MasterHeaderSize = %d\n", irank, MasterHeaderSize);

  int i, j;

  if( PhastaIONextActiveIndex == MAX_PHASTA_FILES ) {
    printf("Error initphmpiio: PhastaIONextActiveIndex = MAX_PHASTA_FILES");
    endTimer(&timer_end);
    printPerf("initphmpiio", timer_start, timer_end, 0, 0, "");
    return MAX_PHASTA_FILES_EXCEEDED;
  }
  //		else if( PhastaIONextActiveIndex == 0 )  //Hang in debug mode on Intrepid
  //		{
  //			for( i = 0; i < MAX_PHASTA_FILES; i++ );
  //			{
  //				PhastaIOActiveFiles[i] = NULL;
  //			}
  //		}


  PhastaIOActiveFiles[PhastaIONextActiveIndex] = (phastaio_file_t *)calloc( 1,  sizeof( phastaio_file_t) );

  i = PhastaIONextActiveIndex;
  PhastaIONextActiveIndex++;

  //PhastaIOActiveFiles[i]->next_start_address = 2*TWO_MEGABYTE;

  PhastaIOActiveFiles[i]->next_start_address = MasterHeaderSize;  // what does this mean??? TODO

  PhastaIOActiveFiles[i]->Wrong_Endian = false;

  PhastaIOActiveFiles[i]->nFields = *nfields;
  PhastaIOActiveFiles[i]->nPPF = *nppf;
  PhastaIOActiveFiles[i]->nFiles = *nfiles;
  MPI_Comm_rank(my_local_comm, &(PhastaIOActiveFiles[i]->myrank));
  MPI_Comm_size(my_local_comm, &(PhastaIOActiveFiles[i]->numprocs));

  int color = computeColor(PhastaIOActiveFiles[i]->myrank, PhastaIOActiveFiles[i]->numprocs, PhastaIOActiveFiles[i]->nFiles);
  MPI_Comm_split(my_local_comm,
      color,
      PhastaIOActiveFiles[i]->myrank,
      &(PhastaIOActiveFiles[i]->local_comm));
  MPI_Comm_size(PhastaIOActiveFiles[i]->local_comm,
      &(PhastaIOActiveFiles[i]->local_numprocs));
  MPI_Comm_rank(PhastaIOActiveFiles[i]->local_comm,
      &(PhastaIOActiveFiles[i]->local_myrank));
  PhastaIOActiveFiles[i]->nppp =
    PhastaIOActiveFiles[i]->nPPF/PhastaIOActiveFiles[i]->local_numprocs;

  PhastaIOActiveFiles[i]->start_id = PhastaIOActiveFiles[i]->nPPF *
    (int)(PhastaIOActiveFiles[i]->myrank/PhastaIOActiveFiles[i]->local_numprocs) +
    (PhastaIOActiveFiles[i]->local_myrank * PhastaIOActiveFiles[i]->nppp);

  PhastaIOActiveFiles[i]->my_offset_table =
    ( unsigned long ** ) calloc( MAX_FIELDS_NUMBER , sizeof( unsigned long *) );

  PhastaIOActiveFiles[i]->my_read_table =
    ( unsigned long ** ) calloc( MAX_FIELDS_NUMBER , sizeof( unsigned long *) );

  for (j=0; j<*nfields; j++)
  {
    PhastaIOActiveFiles[i]->my_offset_table[j] =
      ( unsigned long * ) calloc( PhastaIOActiveFiles[i]->nppp , sizeof( unsigned long) );

    PhastaIOActiveFiles[i]->my_read_table[j] =
      ( unsigned long * ) calloc( PhastaIOActiveFiles[i]->nppp , sizeof( unsigned long) );
  }
  *filehandle = i;

  PhastaIOActiveFiles[i]->master_header = (char *)calloc(MasterHeaderSize,sizeof( char ));
  PhastaIOActiveFiles[i]->double_chunk = (double *)calloc(1,sizeof( double ));
  PhastaIOActiveFiles[i]->int_chunk = (int *)calloc(1,sizeof( int ));
  PhastaIOActiveFiles[i]->read_double_chunk = (double *)calloc(1,sizeof( double ));
  PhastaIOActiveFiles[i]->read_int_chunk = (int *)calloc(1,sizeof( int ));

  // Time monitoring
  endTimer(&timer_end);
  printPerf("initphmpiiosub", timer_start, timer_end, 0, 0, "");

  phprintf_0("Info initphmpiiosub: quiting function");

  return i;
}

namespace {

  enum {
    DIR_MODE = S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH
  };

  bool my_mkdir(std::string name) {
    if(name.empty())
      return true;
    errno = 0;
    int err = mkdir(name.c_str(), DIR_MODE);
    if ((err == -1) && (errno == EEXIST)) {
      errno = 0;
      err = 0;
      return false;
    }
    assert(!err);
    return true;
  }

  enum {
    DIR_FANOUT = 2048
  };

  std::string getSubDirPrefix() {
    if (phio_peers() <= DIR_FANOUT)
      return string("");
    int self = phio_self();
    int subSelf = self % DIR_FANOUT;
    int subGroup = self / DIR_FANOUT;
    std::stringstream ss;
    ss << subGroup << '/';
    return ss.str();
  }
}

/** open file for both POSIX and MPI-IO syncIO format.
 *
 * If it's old POSIX format, simply call posix fopen().
 *
 * If it's MPI-IO foramt:
 * in "read" mode, it builds the header table that points to the offset of
 * fields for parts;
 * in "write" mode, it opens the file with MPI-IO open routine.
 */
void openfile(const char filename[], const char mode[], int*  fileDescriptor )
{
  phprintf_0("Info: entering openfile");

  double timer_start, timer_end;
  startTimer(&timer_start);

  if ( PhastaIONextActiveIndex == 0 )
  {
    FILE* file=NULL ;
    *fileDescriptor = 0;
    char* fname = StringStripper( filename );
    char* imode = StringStripper( mode );

    std::string posixname = getSubDirPrefix();
    if (!phio_self())
      my_mkdir(posixname);
    phio_barrier();
    posixname += string(fname);
    if ( cscompare( "read", imode ) )
      PHASTAIO_OPENTIME(file = fopen(posixname.c_str(), "rb" );)
    else if( cscompare( "write", imode ) )
      PHASTAIO_OPENTIME(file = fopen(posixname.c_str(), "wb" );)
    else if( cscompare( "append", imode ) )
      PHASTAIO_OPENTIME(file = fopen(posixname.c_str(), "ab" );)

    if ( !file ){
      fprintf(stderr,"Error openfile: unable to open file %s\n",fname);
    } else {
      fileArray.push_back( file );
      byte_order.push_back( false );
      header_type.push_back( sizeof(int) );
      *fileDescriptor = fileArray.size();
    }
    free (fname);
    free (imode);
  }
  else // else it would be parallel I/O, opposed to posix io
  {
    char* fname = StringStripper( filename );
    char* imode = StringStripper( mode );
    int rc;
    int i = *fileDescriptor;
    checkFileDescriptor("openfile",&i);
    char* token;

    if ( cscompare( "read", imode ) )
    {
      //	      if (PhastaIOActiveFiles[i]->myrank == 0)
      //                printf("\n **********\nRead open ... ... regular version\n");

      PHASTAIO_OPENTIME(
        rc = MPI_File_open( PhastaIOActiveFiles[i]->local_comm,
            fname,
            MPI_MODE_RDONLY,
            MPI_INFO_NULL,
            &(PhastaIOActiveFiles[i]->file_handle) );
      )

      if(rc)
      {
        *fileDescriptor = UNABLE_TO_OPEN_FILE;
        int error_string_length;
        char error_string[4096];
        MPI_Error_string(rc, error_string, &error_string_length);
        fprintf(stderr, "Error openfile: Unable to open file %s! MPI reports \"%s\"\n",fname,error_string);
        endTimer(&timer_end);
        printPerf("openfile", timer_start, timer_end, 0, 0, "");
        return;
      }

      MPI_Status read_tag_status;
      char read_out_tag[MAX_FIELDS_NAME_LENGTH];
      int j;
      int magic_number;

      if ( PhastaIOActiveFiles[i]->local_myrank == 0 ) {
        MPI_File_read_at( PhastaIOActiveFiles[i]->file_handle,
            0,
            PhastaIOActiveFiles[i]->master_header,
            MasterHeaderSize,
            MPI_CHAR,
            &read_tag_status );
      }

      MPI_Bcast( PhastaIOActiveFiles[i]->master_header,
          MasterHeaderSize,
          MPI_CHAR,
          0,
          PhastaIOActiveFiles[i]->local_comm );

      memcpy( read_out_tag,
          PhastaIOActiveFiles[i]->master_header,
          MAX_FIELDS_NAME_LENGTH-1 );

      if ( cscompare ("MPI_IO_Tag",read_out_tag) )
      {
        // Test endianess ...
        memcpy ( &magic_number,
            PhastaIOActiveFiles[i]->master_header+sizeof("MPI_IO_Tag : ")-1, //-1 sizeof returns the size of the string+1 for "\0"
            sizeof(int) );                                                   // masterheader should look like "MPI_IO_Tag : 12180 " with 12180 in binary format

        if ( magic_number != ENDIAN_TEST_NUMBER )
        {
          PhastaIOActiveFiles[i]->Wrong_Endian = true;
        }

        memcpy( read_out_tag,
            PhastaIOActiveFiles[i]->master_header+MAX_FIELDS_NAME_LENGTH+1, // TODO: WHY +1???
            MAX_FIELDS_NAME_LENGTH );

        // Read in # fields ...
        token = strtok ( read_out_tag, ":" );
        token = strtok( NULL," ,;<>" );
        PhastaIOActiveFiles[i]->nFields = atoi( token );

        unsigned long **header_table;
        header_table = ( unsigned long ** )calloc(PhastaIOActiveFiles[i]->nFields, sizeof(unsigned long *));

        for ( j = 0; j < PhastaIOActiveFiles[i]->nFields; j++ )
        {
          header_table[j]=( unsigned long * ) calloc( PhastaIOActiveFiles[i]->nPPF , sizeof( unsigned long));
        }

        // Read in the offset table ...
        for ( j = 0; j < PhastaIOActiveFiles[i]->nFields; j++ )
        {
          if ( PhastaIOActiveFiles[i]->local_myrank == 0 ) {
            memcpy( header_table[j],
                PhastaIOActiveFiles[i]->master_header +
                VERSION_INFO_HEADER_SIZE +
                j * PhastaIOActiveFiles[i]->nPPF * sizeof(unsigned long),
                PhastaIOActiveFiles[i]->nPPF * sizeof(unsigned long) );
          }

          MPI_Scatter( header_table[j],
              PhastaIOActiveFiles[i]->nppp,
              MPI_LONG_LONG_INT,
              PhastaIOActiveFiles[i]->my_read_table[j],
              PhastaIOActiveFiles[i]->nppp,
              MPI_LONG_LONG_INT,
              0,
              PhastaIOActiveFiles[i]->local_comm );

          // Swap byte order if endianess is different ...
          if ( PhastaIOActiveFiles[i]->Wrong_Endian ) {
            SwapArrayByteOrder( PhastaIOActiveFiles[i]->my_read_table[j],
                sizeof(unsigned long),
                PhastaIOActiveFiles[i]->nppp );
          }
        }

        for ( j = 0; j < PhastaIOActiveFiles[i]->nFields; j++ ) {
          free ( header_table[j] );
        }
        free (header_table);

      } // end of if MPI_IO_TAG
      else //else not valid MPI file
      {
        *fileDescriptor = NOT_A_MPI_FILE;
        printf("Error openfile: The file %s you opened is not in syncIO new format, please check again! File descriptor = %d, MasterHeaderSize = %d, read_out_tag = %s\n",fname,*fileDescriptor,MasterHeaderSize,read_out_tag);
        //Printing MasterHeaderSize is useful to test a compiler bug on Intrepid BGP
        endTimer(&timer_end);
        printPerf("openfile", timer_start, timer_end, 0, 0, "");
        return;
      }
    } // end of if "read"
    else if( cscompare( "write", imode ) )
    {
      PHASTAIO_OPENTIME(
      rc = MPI_File_open( PhastaIOActiveFiles[i]->local_comm,
          fname,
          MPI_MODE_WRONLY | MPI_MODE_CREATE,
          MPI_INFO_NULL,
          &(PhastaIOActiveFiles[i]->file_handle) );
      )
      if(rc != MPI_SUCCESS)
      {
        *fileDescriptor = UNABLE_TO_OPEN_FILE;
        int error_string_length;
        char error_string[4096];
        MPI_Error_string(rc, error_string, &error_string_length);
        fprintf(stderr, "Error openfile: Unable to open file %s! MPI reports \"%s\"\n",fname,error_string);
        return;
      }
    } // end of if "write"
    free (fname);
    free (imode);
  } // end of if FileIndex != 0

  endTimer(&timer_end);
  printPerf("openfile", timer_start, timer_end, 0, 0, "");
}

/** close file for both POSIX and MPI-IO syncIO format.
 *
 * If it's old POSIX format, simply call posix fclose().
 *
 * If it's MPI-IO foramt:
 * in "read" mode, it simply close file with MPI-IO close routine.
 * in "write" mode, rank 0 in each group will re-assemble the master header and
 * offset table and write to the beginning of file, then close the file.
 */
void closefile( int* fileDescriptor, const char mode[] )
{
  double timer_start, timer_end;
  startTimer(&timer_start);

  int i = *fileDescriptor;
  checkFileDescriptor("closefile",&i);

  if ( PhastaIONextActiveIndex == 0 ) {
    char* imode = StringStripper( mode );

    if( cscompare( "write", imode )
        || cscompare( "append", imode ) ) {
      fflush( fileArray[ *fileDescriptor - 1 ] );
    }

    PHASTAIO_CLOSETIME(fclose( fileArray[ *fileDescriptor - 1 ] );)
    free (imode);
  }
  else {
    char* imode = StringStripper( mode );

    //write master header here:
    if ( cscompare( "write", imode ) ) {
      //	      if ( PhastaIOActiveFiles[i]->nPPF * PhastaIOActiveFiles[i]->nFields < 2*ONE_MEGABYTE/8 )  //SHOULD BE CHECKED
      //		MasterHeaderSize = 4*ONE_MEGABYTE;
      //	      else
      //		MasterHeaderSize = 4*ONE_MEGABYTE + PhastaIOActiveFiles[i]->nPPF * PhastaIOActiveFiles[i]->nFields * 8 - 2*ONE_MEGABYTE;

      MasterHeaderSize = computeMHSize( PhastaIOActiveFiles[i]->nFields, PhastaIOActiveFiles[i]->nPPF, LATEST_WRITE_VERSION);
      phprintf_0("Info closefile: myrank = %d, MasterHeaderSize = %d\n", PhastaIOActiveFiles[i]->myrank, MasterHeaderSize);

      MPI_Status write_header_status;
      char mpi_tag[MAX_FIELDS_NAME_LENGTH];
      char version[MAX_FIELDS_NAME_LENGTH/4];
      char mhsize[MAX_FIELDS_NAME_LENGTH/4];
      int magic_number = ENDIAN_TEST_NUMBER;

      if ( PhastaIOActiveFiles[i]->local_myrank == 0 )
      {
        bzero((void*)mpi_tag,MAX_FIELDS_NAME_LENGTH);
        sprintf(mpi_tag, "MPI_IO_Tag : ");
        memcpy(PhastaIOActiveFiles[i]->master_header,
            mpi_tag,
            MAX_FIELDS_NAME_LENGTH);

        bzero((void*)version,MAX_FIELDS_NAME_LENGTH/4);
        // this version is "1", print version in ASCII
        sprintf(version, "version : %d",1);
        memcpy(PhastaIOActiveFiles[i]->master_header + MAX_FIELDS_NAME_LENGTH/2,
            version,
            MAX_FIELDS_NAME_LENGTH/4);

        // master header size is computed using the formula above
        bzero((void*)mhsize,MAX_FIELDS_NAME_LENGTH/4);
        sprintf(mhsize, "mhsize : ");
        memcpy(PhastaIOActiveFiles[i]->master_header + MAX_FIELDS_NAME_LENGTH/4*3,
            mhsize,
            MAX_FIELDS_NAME_LENGTH/4);

        bzero((void*)mpi_tag,MAX_FIELDS_NAME_LENGTH);
        sprintf(mpi_tag,
            "\nnFields : %d\n",
            PhastaIOActiveFiles[i]->nFields);
        memcpy(PhastaIOActiveFiles[i]->master_header+MAX_FIELDS_NAME_LENGTH,
            mpi_tag,
            MAX_FIELDS_NAME_LENGTH);

        bzero((void*)mpi_tag,MAX_FIELDS_NAME_LENGTH);
        sprintf(mpi_tag, "\nnPPF : %d\n", PhastaIOActiveFiles[i]->nPPF);
        memcpy( PhastaIOActiveFiles[i]->master_header+
            PhastaIOActiveFiles[i]->nFields *
            MAX_FIELDS_NAME_LENGTH +
            MAX_FIELDS_NAME_LENGTH * 2,
            mpi_tag,
            MAX_FIELDS_NAME_LENGTH);

        memcpy( PhastaIOActiveFiles[i]->master_header+sizeof("MPI_IO_Tag : ")-1, //-1 sizeof returns the size of the string+1 for "\0"
            &magic_number,
            sizeof(int));

        memcpy( PhastaIOActiveFiles[i]->master_header+sizeof("mhsize : ") -1 + MAX_FIELDS_NAME_LENGTH/4*3,
            &MasterHeaderSize,
            sizeof(int));
      }

      int j = 0;
      unsigned long **header_table;
      header_table = ( unsigned long ** )calloc(PhastaIOActiveFiles[i]->nFields, sizeof(unsigned long *));

      for ( j = 0; j < PhastaIOActiveFiles[i]->nFields; j++ ) {
        header_table[j]=( unsigned long * ) calloc( PhastaIOActiveFiles[i]->nPPF , sizeof( unsigned long));
      }

      //if( irank == 0 ) printf("gonna mpi_gather, myrank = %d\n", irank);
      for ( j = 0; j < PhastaIOActiveFiles[i]->nFields; j++ ) {
        MPI_Gather( PhastaIOActiveFiles[i]->my_offset_table[j],
            PhastaIOActiveFiles[i]->nppp,
            MPI_LONG_LONG_INT,
            header_table[j],
            PhastaIOActiveFiles[i]->nppp,
            MPI_LONG_LONG_INT,
            0,
            PhastaIOActiveFiles[i]->local_comm );
      }

      if ( PhastaIOActiveFiles[i]->local_myrank == 0 ) {

        //if( irank == 0 ) printf("gonna memcpy for every procs, myrank = %d\n", irank);
        for ( j = 0; j < PhastaIOActiveFiles[i]->nFields; j++ ) {
          memcpy ( PhastaIOActiveFiles[i]->master_header +
              VERSION_INFO_HEADER_SIZE +
              j * PhastaIOActiveFiles[i]->nPPF * sizeof(unsigned long),
              header_table[j],
              PhastaIOActiveFiles[i]->nPPF * sizeof(unsigned long) );
        }

        //if( irank == 0 ) printf("gonna file_write_at(), myrank = %d\n", irank);
        MPI_File_write_at( PhastaIOActiveFiles[i]->file_handle,
            0,
            PhastaIOActiveFiles[i]->master_header,
            MasterHeaderSize,
            MPI_CHAR,
            &write_header_status );
      }

      ////free(PhastaIOActiveFiles[i]->master_header);

      for ( j = 0; j < PhastaIOActiveFiles[i]->nFields; j++ ) {
        free ( header_table[j] );
      }
      free (header_table);
    }

    //if( irank == 0 ) printf("gonna file_close(), myrank = %d\n", irank);
    PHASTAIO_CLOSETIME(
      MPI_File_close( &( PhastaIOActiveFiles[i]->file_handle ) );
    )
    free ( imode );
  }

  endTimer(&timer_end);
  printPerf("closefile_", timer_start, timer_end, 0, 0, "");
}

int readHeader( FILE* f, const char phrase[],
    int* params, int numParams, const char* iotype) {
  isBinary(iotype);
  return readHeader(f,phrase,params,numParams);
}

void readheader(
    int* fileDescriptor,
    const  char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] )
{
  double timer_start, timer_end;

  startTimer(&timer_start);

  int i = *fileDescriptor;
  checkFileDescriptor("readheader",&i);

  if ( PhastaIONextActiveIndex == 0 ) {
    int filePtr = *fileDescriptor - 1;
    FILE* fileObject;
    int* valueListInt;

    if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
      fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
      fprintf(stderr,"openfile_ function has to be called before \n") ;
      fprintf(stderr,"acessing the file\n ") ;
      fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
      endTimer(&timer_end);
      printPerf("readheader", timer_start, timer_end, 0, 0, "");
      return;
    }

    LastHeaderKey[filePtr] = std::string(keyphrase);
    LastHeaderNotFound = false;

    fileObject = fileArray[ filePtr ] ;
    Wrong_Endian = byte_order[ filePtr ];

    isBinary( iotype );
    typeSize( datatype );   //redundant call, just avoid a compiler warning.

    // right now we are making the assumption that we will only write integers
    // on the header line.

    valueListInt = static_cast< int* >( valueArray );
    int ierr = readHeader( fileObject ,
        keyphrase,
        valueListInt,
        *nItems ) ;

    byte_order[ filePtr ] = Wrong_Endian ;

    if ( ierr ) LastHeaderNotFound = true;

    //return ; // don't return, go to the end to print perf
  }
  else {
    int* valueListInt;
    valueListInt = static_cast <int*>(valueArray);
    char* token = NULL;
    bool FOUND = false ;
    isBinary( iotype );

    MPI_Status read_offset_status;
    char read_out_tag[MAX_FIELDS_NAME_LENGTH];
    memset(read_out_tag, '\0', MAX_FIELDS_NAME_LENGTH);
    char readouttag[MAX_FIELDS_NUMBER][MAX_FIELDS_NAME_LENGTH];
    int j;

    int string_length = strlen( keyphrase );
    char* buffer = (char*) malloc ( string_length+1 );

    strcpy ( buffer, keyphrase );
    buffer[ string_length ] = '\0';

    char* st2 = strtok ( buffer, "@" );
    st2 = strtok (NULL, "@");
    PhastaIOActiveFiles[i]->GPid = atoi(st2);
    if ( char* p = strpbrk(buffer, "@") )
      *p = '\0';

    // Check if the user has input the right GPid
    if ( ( PhastaIOActiveFiles[i]->GPid <=
          PhastaIOActiveFiles[i]->myrank *
          PhastaIOActiveFiles[i]->nppp )||
        ( PhastaIOActiveFiles[i]->GPid >
          ( PhastaIOActiveFiles[i]->myrank + 1 ) *
          PhastaIOActiveFiles[i]->nppp ) )
    {
      *fileDescriptor = NOT_A_MPI_FILE;
      printf("Error readheader: The file is not in syncIO new format, please check! myrank = %d, GPid = %d, nppp = %d, keyphrase = %s\n", PhastaIOActiveFiles[i]->myrank, PhastaIOActiveFiles[i]->GPid, PhastaIOActiveFiles[i]->nppp, keyphrase);
      // It is possible atoi could not generate a clear integer from st2 because of additional garbage character in keyphrase
      endTimer(&timer_end);
      printPerf("readheader", timer_start, timer_end, 0, 0, "");
      return;
    }

    // Find the field we want ...
    //for ( j = 0; j<MAX_FIELDS_NUMBER; j++ )
    for ( j = 0; j<PhastaIOActiveFiles[i]->nFields; j++ )
    {
      memcpy( readouttag[j],
          PhastaIOActiveFiles[i]->master_header + j*MAX_FIELDS_NAME_LENGTH+MAX_FIELDS_NAME_LENGTH*2+1,
          MAX_FIELDS_NAME_LENGTH-1 );
    }

    for ( j = 0; j<PhastaIOActiveFiles[i]->nFields; j++ )
    {
      token = strtok ( readouttag[j], ":" );

      //if ( cscompare( buffer, token ) )
      if ( cscompare( token , buffer ) && cscompare( buffer, token ) )
        // This double comparison is required for the field "number of nodes" and all the other fields that start with "number of nodes" (i.g. number of nodes in the mesh").
        // Would be safer to rename "number of nodes" by "number of nodes in the part" so that the name are completely unique. But much more work to do that (Nspre, phParAdapt, etc).
        // Since the field name are unique in SyncIO (as it includes part ID), this should be safe and there should be no issue with the "?" trailing character.
      {
        PhastaIOActiveFiles[i]->read_field_count = j;
        FOUND = true;
        //printf("buffer: %s | token: %s | j: %d\n",buffer,token,j);
        break;
      }
    }
    free(buffer);

    if (!FOUND)
    {
      //if(irank==0) printf("Warning readheader: Not found %s \n",keyphrase); //PhastaIOActiveFiles[i]->myrank is certainly initialized here.
      if(PhastaIOActiveFiles[i]->myrank == 0) printf("WARNING readheader: Not found %s\n",keyphrase);
      endTimer(&timer_end);
      printPerf("readheader", timer_start, timer_end, 0, 0, "");
      return;
    }

    // Find the part we want ...
    PhastaIOActiveFiles[i]->read_part_count = PhastaIOActiveFiles[i]->GPid -
      PhastaIOActiveFiles[i]->myrank * PhastaIOActiveFiles[i]->nppp - 1;

    PhastaIOActiveFiles[i]->my_offset =
      PhastaIOActiveFiles[i]->my_read_table[PhastaIOActiveFiles[i]->read_field_count][PhastaIOActiveFiles[i]->read_part_count];

    // 	  printf("****Rank %d offset is %d\n",PhastaIOActiveFiles[i]->myrank,PhastaIOActiveFiles[i]->my_offset);

    // Read each datablock header here ...

    MPI_File_read_at_all( PhastaIOActiveFiles[i]->file_handle,
        PhastaIOActiveFiles[i]->my_offset+1,
        read_out_tag,
        MAX_FIELDS_NAME_LENGTH-1,
        MPI_CHAR,
        &read_offset_status );
    token = strtok ( read_out_tag, ":" );

    // 	  printf("&&&&Rank %d read_out_tag is %s\n",PhastaIOActiveFiles[i]->myrank,read_out_tag);

    if( cscompare( keyphrase , token ) ) //No need to compare also token with keyphrase like above. We should already have the right one. Otherwise there is a problem.
    {
      FOUND = true ;
      token = strtok( NULL, " ,;<>" );
      for( j=0; j < *nItems && ( token = strtok( NULL," ,;<>") ); j++ )
        valueListInt[j] = atoi( token );

      if ( j < *nItems )
      {
        fprintf( stderr, "Expected # of ints not found for: %s\n", keyphrase );
      }
    }
    else {
      //if(irank==0)
      if(PhastaIOActiveFiles[i]->myrank == 0)
        // If we enter this if, there is a problem with the name of some fields
      {
        printf("Error readheader: Unexpected mismatch between keyphrase = %s and token = %s\n",keyphrase,token);
      }
    }
  }

  endTimer(&timer_end);
  printPerf("readheader", timer_start, timer_end, 0, 0, "");

}

void readDataBlock(
    FILE* fileObject,
    void* valueArray,
    int nItems,
    const char  datatype[],
    const char  iotype[] )
{
  isBinary(iotype);
  size_t type_size = typeSize( datatype );
  phastaioTime t0,t1;
  phastaio_time(&t0);
  if ( binary_format ) {
    char junk = '\0';
    fread( valueArray, type_size, nItems, fileObject );
    fread( &junk, sizeof(char), 1 , fileObject );
    if ( Wrong_Endian ) SwapArrayByteOrder( valueArray, type_size, nItems );
  } else {
    char* ts1 = StringStripper( datatype );
    if ( cscompare( "integer", ts1 ) ) {
      for( int n=0; n < nItems ; n++ )
        fscanf(fileObject, "%d\n",(int*)((int*)valueArray+n) );
    } else if ( cscompare( "double", ts1 ) ) {
      for( int n=0; n < nItems ; n++ )
        fscanf(fileObject, "%lf\n",(double*)((double*)valueArray+n) );
    }
    free (ts1);
  }
  phastaio_time(&t1);
  const size_t elapsed = phastaio_time_diff(&t0,&t1);
  phastaio_addReadTime(elapsed);
  phastaio_addReadBytes(nItems*type_size);
}

void readdatablock(
    int*  fileDescriptor,
    const char keyphrase[],
    void* valueArray,
    int*  nItems,
    const char  datatype[],
    const char  iotype[] )
{
  //if(irank == 0) printf("entering readdatablock()\n");
  unsigned long data_size = 0;
  double timer_start, timer_end;
  startTimer(&timer_start);

  int i = *fileDescriptor;
  checkFileDescriptor("readdatablock",&i);

  if ( PhastaIONextActiveIndex == 0 ) {
    int filePtr = *fileDescriptor - 1;
    FILE* fileObject;

    if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
      fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
      fprintf(stderr,"openfile_ function has to be called before\n") ;
      fprintf(stderr,"acessing the file\n ") ;
      fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
      endTimer(&timer_end);
      printPerf("readdatablock", timer_start, timer_end, 0, 0, "");
      return;
    }

    // error check..
    // since we require that a consistant header always preceed the data block
    // let us check to see that it is actually the case.

    if ( ! cscompare( LastHeaderKey[ filePtr ].c_str(), keyphrase ) ) {
      fprintf(stderr, "Header not consistant with data block\n");
      fprintf(stderr, "Header: %s\n", LastHeaderKey[ filePtr ].c_str() );
      fprintf(stderr, "DataBlock: %s\n ", keyphrase );
      fprintf(stderr, "Please recheck read sequence \n");
      if( Strict_Error ) {
        fprintf(stderr, "fatal error: cannot continue, returning out of call\n");
        endTimer(&timer_end);
        printPerf("readdatablock", timer_start, timer_end, 0, 0, "");
        return;
      }
    }

    if ( LastHeaderNotFound ) {
      endTimer(&timer_end);
      printPerf("readdatablock", timer_start, timer_end, 0, 0, "");
      return;
    }
    fileObject = fileArray[ filePtr ];
    Wrong_Endian = byte_order[ filePtr ];
    LastHeaderKey.erase(filePtr);
    readDataBlock(fileObject,valueArray,*nItems,datatype,iotype);

    //return;
  }
  else {
    // 	  printf("read data block\n");
    MPI_Status read_data_status;
    size_t type_size = typeSize( datatype );
    int nUnits = *nItems;
    isBinary( iotype );

    // read datablock then
    //MR CHANGE
    //             if ( cscompare ( datatype, "double"))
    char* ts2 = StringStripper( datatype );
    if ( cscompare ( "double" , ts2))
      //MR CHANGE END
    {

      phastaioTime t0,t1;
      phastaio_time(&t0);
      MPI_File_read_at_all_begin( PhastaIOActiveFiles[i]->file_handle,
          PhastaIOActiveFiles[i]->my_offset + DB_HEADER_SIZE,
          valueArray,
          nUnits,
          MPI_DOUBLE );
      MPI_File_read_at_all_end( PhastaIOActiveFiles[i]->file_handle,
          valueArray,
          &read_data_status );
      data_size=8*nUnits;
      phastaio_time(&t1);
      const size_t elapsed = phastaio_time_diff(&t0,&t1);
      phastaio_addReadTime(elapsed);
      phastaio_addReadBytes(nUnits*sizeof(double));
    }
    //MR CHANGE
    //             else if ( cscompare ( datatype, "integer"))
    else if ( cscompare ( "integer" , ts2))
      //MR CHANGE END
    {
      phastaioTime t0,t1;
      phastaio_time(&t0);
      MPI_File_read_at_all_begin(PhastaIOActiveFiles[i]->file_handle,
          PhastaIOActiveFiles[i]->my_offset + DB_HEADER_SIZE,
          valueArray,
          nUnits,
          MPI_INT );
      MPI_File_read_at_all_end( PhastaIOActiveFiles[i]->file_handle,
          valueArray,
          &read_data_status );
      data_size=4*nUnits;
      phastaio_time(&t1);
      const size_t elapsed = phastaio_time_diff(&t0,&t1);
      phastaio_addReadTime(elapsed);
      phastaio_addReadBytes(nUnits*sizeof(int));
    }
    else
    {
      *fileDescriptor = DATA_TYPE_ILLEGAL;
      printf("readdatablock - DATA_TYPE_ILLEGAL - %s\n",datatype);
      endTimer(&timer_end);
      printPerf("readdatablock", timer_start, timer_end, 0, 0, "");
      return;
    }
    free(ts2);


    // 	  printf("%d Read finishe\n",PhastaIOActiveFiles[i]->myrank);

    // Swap data byte order if endianess is different ...
    if ( PhastaIOActiveFiles[i]->Wrong_Endian )
    {
      SwapArrayByteOrder( valueArray, type_size, nUnits );
    }
  }

  endTimer(&timer_end);
  char extra_msg[1024];
  memset(extra_msg, '\0', 1024);
  char* key = StringStripper(keyphrase);
  sprintf(extra_msg, " field is %s ", key);
  printPerf("readdatablock", timer_start, timer_end, data_size, 1, extra_msg);
  free(key);

}

void writeHeader(
    FILE* f,
    const char keyphrase[],
    const void* valueArray,
    const int nItems,
    const int ndataItems,
    const char datatype[],
    const char iotype[])
{
  isBinary( iotype );

  const int _newline =
    ( ndataItems > 0 ) ? sizeof( char ) : 0;
  int size_of_nextblock =
    ( binary_format ) ? typeSize(datatype) * ndataItems + _newline : ndataItems;

  fprintf( f, "%s : < %d > ", keyphrase, size_of_nextblock );
  for( int i = 0; i < nItems; i++ )
    fprintf( f, "%d ", *((int*)((int*)valueArray+i)));
  fprintf( f, "\n");
}

void writeheader(
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const int* ndataItems,
    const char datatype[],
    const char iotype[])
{

  //if(irank == 0) printf("entering writeheader()\n");

  double timer_start, timer_end;
  startTimer(&timer_start);

  int i = *fileDescriptor;
  checkFileDescriptor("writeheader",&i);

  if ( PhastaIONextActiveIndex == 0 ) {
    int filePtr = *fileDescriptor - 1;
    if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
      fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
      fprintf(stderr,"openfile_ function has to be called before \n") ;
      fprintf(stderr,"acessing the file\n ") ;
      fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
      endTimer(&timer_end);
      printPerf("writeheader", timer_start, timer_end, 0, 0, "");
      return;
    }

    LastHeaderKey[filePtr] = std::string(keyphrase);
    DataSize = *ndataItems;
    FILE* fileObject = fileArray[ filePtr ] ;
    header_type[ filePtr ] = typeSize( datatype );
    writeHeader(fileObject,keyphrase,valueArray,*nItems,
        *ndataItems,datatype,iotype);
  }
  else { // else it's parallel I/O
    DataSize = *ndataItems;
    size_t type_size = typeSize( datatype );
    isBinary( iotype );
    char mpi_tag[MAX_FIELDS_NAME_LENGTH];

    int string_length = strlen( keyphrase );
    char* buffer = (char*) malloc ( string_length+1 );

    strcpy ( buffer, keyphrase);
    buffer[ string_length ] = '\0';

    char* st2 = strtok ( buffer, "@" );
    st2 = strtok (NULL, "@");
    PhastaIOActiveFiles[i]->GPid = atoi(st2);

    if ( char* p = strpbrk(buffer, "@") )
      *p = '\0';
    assert(PhastaIOActiveFiles[i]->field_count < MAX_FIELDS_NUMBER);
    assert(PhastaIOActiveFiles[i]->part_count < PhastaIOActiveFiles[i]->nppp);
    bzero((void*)mpi_tag,MAX_FIELDS_NAME_LENGTH);
    sprintf(mpi_tag, "\n%s : %d\n", buffer, PhastaIOActiveFiles[i]->field_count);
    unsigned long offset_value;

    int temp = *ndataItems;
    unsigned long number_of_items = (unsigned long)temp;
    MPI_Barrier(PhastaIOActiveFiles[i]->local_comm);

    MPI_Scan( &number_of_items,
        &offset_value,
        1,
        MPI_LONG_LONG_INT,
        MPI_SUM,
        PhastaIOActiveFiles[i]->local_comm );

    offset_value = (offset_value - number_of_items) * type_size;

    offset_value += PhastaIOActiveFiles[i]->local_myrank *
      DB_HEADER_SIZE +
      PhastaIOActiveFiles[i]->next_start_address;
    // This offset is the starting address of each datablock header...
    PhastaIOActiveFiles[i]->my_offset = offset_value;

    // Write in my offset table ...
    PhastaIOActiveFiles[i]->my_offset_table[PhastaIOActiveFiles[i]->field_count][PhastaIOActiveFiles[i]->part_count] = PhastaIOActiveFiles[i]->my_offset;

    // Update the next-start-address ...
    PhastaIOActiveFiles[i]->next_start_address = offset_value +
      number_of_items * type_size +
      DB_HEADER_SIZE;
    MPI_Bcast( &(PhastaIOActiveFiles[i]->next_start_address),
        1,
        MPI_LONG_LONG_INT,
        PhastaIOActiveFiles[i]->local_numprocs-1,
        PhastaIOActiveFiles[i]->local_comm );

    // Prepare datablock header ...
    int _newline = (*ndataItems>0)?sizeof(char):0;
    unsigned int size_of_nextblock = type_size * (*ndataItems) + _newline;

    //char datablock_header[255];
    //bzero((void*)datablock_header,255);
    char datablock_header[DB_HEADER_SIZE];
    bzero((void*)datablock_header,DB_HEADER_SIZE);

    PhastaIOActiveFiles[i]->GPid = PhastaIOActiveFiles[i]->nppp*PhastaIOActiveFiles[i]->myrank+PhastaIOActiveFiles[i]->part_count;
    sprintf( datablock_header,
        "\n%s : < %u >",
        keyphrase,
        size_of_nextblock );

    for ( int j = 0; j < *nItems; j++ )
    {
      sprintf( datablock_header,
          "%s %d ",
          datablock_header,
          *((int*)((int*)valueArray+j)));
    }
    sprintf( datablock_header,
        "%s\n ",
        datablock_header );

    // Write datablock header ...
    //MR CHANGE
    // 	if ( cscompare(datatype,"double") )
    char* ts1 = StringStripper( datatype );
    if ( cscompare("double",ts1) )
      //MR CHANGE END
    {
      free ( PhastaIOActiveFiles[i]->double_chunk );
      PhastaIOActiveFiles[i]->double_chunk = ( double * )malloc( (sizeof( double )*number_of_items+ DB_HEADER_SIZE));

      double * aa = ( double * )datablock_header;
      memcpy(PhastaIOActiveFiles[i]->double_chunk, aa, DB_HEADER_SIZE);
    }
    //MR CHANGE
    // 	if  ( cscompare(datatype,"integer") )
    else if ( cscompare("integer",ts1) )
      //MR CHANGE END
    {
      free ( PhastaIOActiveFiles[i]->int_chunk );
      PhastaIOActiveFiles[i]->int_chunk = ( int * )malloc( (sizeof( int )*number_of_items+ DB_HEADER_SIZE));

      int * aa = ( int * )datablock_header;
      memcpy(PhastaIOActiveFiles[i]->int_chunk, aa, DB_HEADER_SIZE);
    }
    else {
      //             *fileDescriptor = DATA_TYPE_ILLEGAL;
      printf("writeheader - DATA_TYPE_ILLEGAL - %s\n",datatype);
      endTimer(&timer_end);
      printPerf("writeheader", timer_start, timer_end, 0, 0, "");
      return;
    }
    free(ts1);

    PhastaIOActiveFiles[i]->part_count++;
    if ( PhastaIOActiveFiles[i]->part_count == PhastaIOActiveFiles[i]->nppp ) {
      //A new field will be written
      if ( PhastaIOActiveFiles[i]->local_myrank == 0 ) {
        memcpy( PhastaIOActiveFiles[i]->master_header +
            PhastaIOActiveFiles[i]->field_count *
            MAX_FIELDS_NAME_LENGTH +
            MAX_FIELDS_NAME_LENGTH * 2,
            mpi_tag,
            MAX_FIELDS_NAME_LENGTH-1);
      }
      PhastaIOActiveFiles[i]->field_count++;
      PhastaIOActiveFiles[i]->part_count=0;
    }
    free(buffer);
  }

  endTimer(&timer_end);
  printPerf("writeheader", timer_start, timer_end, 0, 0, "");
}

void writeDataBlock(
    FILE* f,
    const void* valueArray,
    const int   nItems,
    const char  datatype[],
    const char  iotype[]  )
{
  isBinary( iotype );
  size_t type_size = typeSize( datatype );
  phastaioTime t0,t1;
  phastaio_time(&t0);
  if ( binary_format ) {
    fwrite( valueArray, type_size, nItems, f );
    fprintf( f,"\n");
  } else {
    char* ts1 = StringStripper( datatype );
    if ( cscompare( "integer", ts1 ) ) {
      const int* vals = (int*) valueArray;
      for( int n=0; n < nItems ; n++ )
        fprintf(f,"%d\n",vals[n]);
    } else if ( cscompare( "double", ts1 ) ) {
      const double* vals = (double*) valueArray;
      for( int n=0; n < nItems ; n++ )
        fprintf(f,"%f\n",vals[n]);
    }
    free (ts1);
  }
  phastaio_time(&t1);
  const size_t elapsed = phastaio_time_diff(&t0,&t1);
  phastaio_addWriteTime(elapsed);
  phastaio_addWriteBytes(nItems*type_size);
}

void writedatablock(
    const int* fileDescriptor,
    const char keyphrase[],
    const void* valueArray,
    const int* nItems,
    const char datatype[],
    const char iotype[] )
{
  //if(irank == 0) printf("entering writedatablock()\n");

  unsigned long data_size = 0;
  double timer_start, timer_end;
  startTimer(&timer_start);

  int i = *fileDescriptor;
  checkFileDescriptor("writedatablock",&i);

  if ( PhastaIONextActiveIndex == 0 ) {
    int filePtr = *fileDescriptor - 1;

    if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
      fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
      fprintf(stderr,"openfile_ function has to be called before \n") ;
      fprintf(stderr,"acessing the file\n ") ;
      fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
      endTimer(&timer_end);
      printPerf("writedatablock", timer_start, timer_end, 0, 0, "");
      return;
    }
    // since we require that a consistant header always preceed the data block
    // let us check to see that it is actually the case.

    if ( ! cscompare( LastHeaderKey[ filePtr ].c_str(), keyphrase ) ) {
      fprintf(stderr, "Header not consistant with data block\n");
      fprintf(stderr, "Header: %s\n", LastHeaderKey[ filePtr ].c_str() );
      fprintf(stderr, "DataBlock: %s\n ", keyphrase );
      fprintf(stderr, "Please recheck write sequence \n");
      if( Strict_Error ) {
        fprintf(stderr, "fatal error: cannot continue, returning out of call\n");
        endTimer(&timer_end);
        printPerf("writedatablock", timer_start, timer_end, 0, 0, "");
        return;
      }
    }

    FILE* fileObject =  fileArray[ filePtr ] ;
    size_t type_size=typeSize( datatype );
    isBinary( iotype );

    LastHeaderKey.erase(filePtr);

    if ( header_type[filePtr] != (int)type_size ) {
      fprintf(stderr,"header and datablock differ on typeof data in the block for\n");
      fprintf(stderr,"keyphrase : %s\n", keyphrase);
      if( Strict_Error ) {
        fprintf(stderr,"fatal error: cannot continue, returning out of call\n" );
        endTimer(&timer_end);
        printPerf("writedatablock", timer_start, timer_end, 0, 0, "");
        return;
      }
    }

    int nUnits = *nItems;

    if ( nUnits != DataSize ) {
      fprintf(stderr,"header and datablock differ on number of data items for\n");
      fprintf(stderr,"keyphrase : %s\n", keyphrase);
      if( Strict_Error ) {
        fprintf(stderr,"fatal error: cannot continue, returning out of call\n" );
        endTimer(&timer_end);
        printPerf("writedatablock", timer_start, timer_end, 0, 0, "");
        return;
      }
    }
    writeDataBlock(fileObject,valueArray,*nItems,datatype,iotype);
  }
  else {  // syncIO case
    MPI_Status write_data_status;
    isBinary( iotype );
    int nUnits = *nItems;

    //MR CHANGE
    // 	if ( cscompare(datatype,"double") )
    char* ts1 = StringStripper( datatype );
    if ( cscompare("double",ts1) )
      //MR CHANGE END
    {
      memcpy((PhastaIOActiveFiles[i]->double_chunk+DB_HEADER_SIZE/sizeof(double)), valueArray, nUnits*sizeof(double));
      phastaioTime t0,t1;
      phastaio_time(&t0);
      MPI_File_write_at_all_begin( PhastaIOActiveFiles[i]->file_handle,
          PhastaIOActiveFiles[i]->my_offset,
          PhastaIOActiveFiles[i]->double_chunk,
          //BLOCK_SIZE/sizeof(double),
          nUnits+DB_HEADER_SIZE/sizeof(double),
          MPI_DOUBLE );
      MPI_File_write_at_all_end( PhastaIOActiveFiles[i]->file_handle,
          PhastaIOActiveFiles[i]->double_chunk,
          &write_data_status );
      data_size=8*nUnits;
      phastaio_time(&t1);
      const size_t elapsed = phastaio_time_diff(&t0,&t1);
      phastaio_addWriteTime(elapsed);
      phastaio_addWriteBytes((nUnits*sizeof(double))+DB_HEADER_SIZE);
    }
    //MR CHANGE
    // 	else if ( cscompare ( datatype, "integer"))
    else if ( cscompare("integer",ts1) )
      //MR CHANGE END
    {
      memcpy((PhastaIOActiveFiles[i]->int_chunk+DB_HEADER_SIZE/sizeof(int)), valueArray, nUnits*sizeof(int));
      phastaioTime t0,t1;
      phastaio_time(&t0);
      MPI_File_write_at_all_begin( PhastaIOActiveFiles[i]->file_handle,
          PhastaIOActiveFiles[i]->my_offset,
          PhastaIOActiveFiles[i]->int_chunk,
          nUnits+DB_HEADER_SIZE/sizeof(int),
          MPI_INT );
      MPI_File_write_at_all_end( PhastaIOActiveFiles[i]->file_handle,
          PhastaIOActiveFiles[i]->int_chunk,
          &write_data_status );
      data_size=4*nUnits;
      phastaio_time(&t1);
      const size_t elapsed = phastaio_time_diff(&t0,&t1);
      phastaio_addWriteTime(elapsed);
      phastaio_addWriteBytes((nUnits*sizeof(int))+DB_HEADER_SIZE);
    }
    else {
      printf("Error: writedatablock - DATA_TYPE_ILLEGAL - %s\n",datatype);
      endTimer(&timer_end);
      printPerf("writedatablock", timer_start, timer_end, 0, 0, "");
      return;
    }
    free(ts1);
  }

  endTimer(&timer_end);
  char extra_msg[1024];
  memset(extra_msg, '\0', 1024);
  char* key = StringStripper(keyphrase);
  sprintf(extra_msg, " field is %s ", key);
  printPerf("writedatablock", timer_start, timer_end, data_size, 1, extra_msg);
  free(key);

}

void SwapArrayByteOrder( void* array, int nbytes, int nItems )
{
  /* This swaps the byte order for the array of nItems each
     of size nbytes , This will be called only locally  */
  int i,j;
  unsigned char* ucDst = (unsigned char*)array;
  for(i=0; i < nItems; i++) {
    for(j=0; j < (nbytes/2); j++)
      std::swap( ucDst[j] , ucDst[(nbytes - 1) - j] );
    ucDst += nbytes;
  }
}

void writestring( int* fileDescriptor, const char inString[] )
{
  int filePtr = *fileDescriptor - 1;
  FILE* fileObject = fileArray[filePtr] ;
  fprintf(fileObject,"%s",inString );
  return;
}

void Gather_Headers( int* fileDescriptor, vector< string >& headers )
{
  FILE* fileObject;
  char Line[1024];

  fileObject = fileArray[ (*fileDescriptor)-1 ];

  while( !feof(fileObject) ) {
    fgets( Line, 1024, fileObject);
    if ( Line[0] == '#' ) {
      headers.push_back( Line );
    } else {
      break;
    }
  }
  rewind( fileObject );
  clearerr( fileObject );
}

void isWrong( void ) {
  (Wrong_Endian) ? fprintf(stdout,"YES\n") : fprintf(stdout,"NO\n");
}

void togglestrictmode( void ) { Strict_Error = !Strict_Error; }

int isLittleEndian( void )
{
  // this function returns a 1 if the current running architecture is
  // LittleEndian Byte Ordered, else it returns a zero
  union  {
    long a;
    char c[sizeof( long )];
  } endianUnion;
  endianUnion.a = 1 ;
  if ( endianUnion.c[sizeof(long)-1] != 1 ) return 1 ;
  else return 0;
}

namespace PHASTA {
  const char* const PhastaIO_traits<int>::type_string = "integer";
  const char* const PhastaIO_traits<double>::type_string = "double";
}
