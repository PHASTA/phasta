//MR CHANGE
#define OMPI_SKIP_MPICXX 1
//MR CHANGE END
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "commonM2N_c.h"

#if !(defined IOSTREAMH)
#include <iostream>
#include <strstream>
using namespace std;
#endif

#include <FCMangle.h>
#define input FortranCInterface_GLOBAL_(input,INPUT)
#define timer FortranCInterface_GLOBAL_(timer,TIMER)

#ifdef intel
#include <direct.h>
#define chdir _chdir
#else
#include <unistd.h>
#endif

extern "C" char phasta_iotype[80];
char phasta_iotype[80];

extern "C" void input();
//MR CHANGE
extern void setIOparam(); // For SyncIO
//MR CHANGE END

int myrank; /* made file global for ease in debugging */

void
catchDebugger() {
    while (1) { 
      int debuggerPresent=0;
      int fakeSTOP = 1; // please stop HERE and assign as next line 
      // assign or set debuggerPresent=1
      if(debuggerPresent) {
        break;
      }
    }
}

// some useful debugging functions

void
pdarray( void* darray , int start, int end ) {
    for( int i=start; i < end; i++ ){
        cout << ((double*)darray)[i] << endl;
    }
}

void
piarray( void* iarray , int start, int end ) {
    for( int i=start; i < end; i++ ){
        cout << ((int*)iarray)[i] << endl;
    }
}

extern "C" int 
phasta( int argc,   
        char *argv[] ) {
  
    int size,ierr;
    char inpfilename[100];
    char* pauseDebugger = getenv("catchDebugger");
    //cout << "pauseDebugger" << pauseDebugger << endl;

    MPI_Init(&argc,&argv);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

    workfc.numpe = size;
    workfc.myrank = myrank;

#if (defined WIN32)
    if(argc > 2 ){
      catchDebugger();
    }
#endif
#if (1) // ALWAYS ( defined LAUNCH_GDB ) && !( defined WIN32 )

    if ( pauseDebugger ) {

        int parent_pid = getpid();
        int gdb_child = fork();
        cout << "gdb_child" << gdb_child << endl;

        if( gdb_child == 0 ) {
     
            cout << "Debugger Process initiating" << endl;
            strstream exec_string;
         
#if ( defined decalp )
            exec_string <<"xterm -e idb " 
                        << " -pid "<< parent_pid <<" "<< argv[0] << endl;
#endif
#if ( defined LINUX )
            exec_string <<"xterm -e gdb"
                        << " -pid "<< parent_pid <<" "<< argv[0] << endl;
#endif
#if ( defined SUN4 )
            exec_string <<"xterm -e dbx " 
                        << " - "<< parent_pid <<" "<< argv[0] << endl;
#endif
#if ( defined IRIX )
            exec_string <<"xterm -e dbx " 
                        << " -p "<< parent_pid <<" "<< argv[0] << endl;
#endif
            system( exec_string.str() );
            exit(0);
        }
        catchDebugger();
    }

#endif

//MR CHANGE
        setIOparam();
//MR CHANGE END

        input();
//MR CHANGE

    MPI_Barrier(MPI_COMM_WORLD);
    if ( myrank == 0 ) {
      printf("phasta.cc - last call before finalize!\n");
    }
//MR CHANGE
 
    MPI_Finalize();
    return 0;
}
