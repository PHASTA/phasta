#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cassert>
#ifdef HAVE_PETSC
#include <petscsys.h>
#include <petscviewer.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>

#include "common_c.h"
#include "Input.h"
#include "phstream.h"
#include "streamio.h"

#if !(defined IOSTREAMH)
#include <iostream>
#include <sstream>
using namespace std;
#endif

#include <FCMangle.h>
#define input FortranCInterface_GLOBAL_(input,INPUT)
#define proces FortranCInterface_GLOBAL_(proces,PROCES)
#define timer FortranCInterface_GLOBAL_(timer,TIMER)

#ifdef intel
#include <direct.h>
#define chdir _chdir
#else
#include <unistd.h>
#endif

extern "C" char phasta_iotype[80];
char phasta_iotype[80];

extern int SONFATH;
extern "C" void proces();
extern "C" void input();
extern int input_fform(phSolver::Input&);
extern void setIOparam(); // For SyncIO
extern "C" void initPhastaCommonVars();

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

namespace {
  int cdToParent() {
    if( chdir("..") ) {
      fprintf(stderr,"could not change to the parent directory\n");
      return 1;
    } else {
      return 0;
    }
  }
  int run(phSolver::Input& ctrl) {
    int size,ierr;
    char inpfilename[100];
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

    workfc.numpe = size;
    workfc.myrank = myrank;

    initPhastaCommonVars();
    /* Input data  */
    ierr = input_fform(ctrl);
    if(!ierr){
      sprintf(inpfilename,"%d-procs_case/",size);
      if( chdir( inpfilename ) ) {
        cerr << "could not change to the problem directory "
          << inpfilename << endl;
        return -1;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      input();
      /* now we can start the solver */
      proces();
    }
    else{
      printf("error during reading ascii input \n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if ( myrank == 0 ) {
      printf("phasta.cc - last call before finalize!\n");
    }
    if( cdToParent() )
      return -1;
    return timdat.lstep;
  }
}

int phasta(phSolver::Input& ctrl) {
  outpar.input_mode = 0; //FIXME magic value for posix
  outpar.output_mode = 0; //FIXME magic value for posix
  return run(ctrl);
}

int phasta(phSolver::Input& ctrl, grstream grs) {
  assert(grs);
  outpar.input_mode = -1; //FIXME magic value for streams
  outpar.output_mode = 1; //FIXME magic value for syncio
  streamio_set_gr(grs);
  return run(ctrl);
}

int phasta(phSolver::Input& ctrl, RStream* rs) {
  fprintf(stderr, "HEY! if you see this email Cameron and tell him "
      "to implement %s(...) on line %d of %s "
      "... returning an error\n", __func__, __LINE__, __FILE__);
  return -1;
}

int phasta(phSolver::Input& ctrl, GRStream* grs, RStream* rs) {
  outpar.input_mode = -1; //FIXME magic value for streams
  outpar.output_mode = -1; //FIXME magic value for streams
  assert(grs);
  assert(rs);
  streamio_set_gr(grs);
  streamio_set_r(rs);
  return run(ctrl);
}

int phasta( int argc, char *argv[] ) {
    int size,ierr;
    char inpfilename[100];
    char* pauseDebugger = getenv("catchDebugger");
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

#ifdef HAVE_PETSC
    PETSC_COMM_WORLD=MPI_COMM_WORLD;
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    PetscInitializeFortran();
    PetscPopSignalHandler(); //Let us segfault in peace ;-)
    PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD);
// ok with Master    PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD);
// ok with 3.6x    PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD);
    if(sizeof(PetscInt) != sizeof(long long int))
    {
      //PetscInt and gcorp_t (gen_ncorp.c)
      //must be the same size. hard-coded for now
      //FIXME
	    if(myrank == 0)
	    {
		    printf("WARNING: PETSc Index Size Mismatch\n");
		    printf("WARNING: Proceed at your own risk\n");
	    }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == 0)
    {
	    printf("PETSc Initialized\n");
	    fflush(stdout);
    }
#endif
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
            stringstream exec_string;
         
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
            string s = exec_string.str();
            system( s.c_str() );
            exit(0);
        }
        catchDebugger();
    }

#endif

    /* Input data  */
    if(argc > 1 ){
        strcpy(inpfilename,argv[1]);
    } else {
        strcpy(inpfilename,"solver.inp");
    }
    string defaultConf = ".";
    const char* path_to_config = getenv("PHASTA_CONFIG");
    if(path_to_config) 
      defaultConf = path_to_config;
    defaultConf.append("/input.config");
    string userConf(inpfilename);
    phSolver::Input ctrl(userConf, defaultConf);
    initPhastaCommonVars();
    ierr = input_fform(ctrl);
    if(!ierr){
      sprintf(inpfilename,"%d-procs_case/",size);
      if( chdir( inpfilename ) ) {
        cerr << "could not change to the problem directory "
          << inpfilename << endl;
        return -1;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      setIOparam();
      outpar.input_mode = outpar.nsynciofiles; //FIXME this is awful
      outpar.output_mode = outpar.nsynciofiles; //FIXME this is awful
      input();
      /* now we can start the solver */
      proces();
    }
    else{
        printf("error during reading ascii input \n");
    }
#ifdef HAVE_PETSC
    PetscFinalize();
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    if ( myrank == 0 ) {
      printf("phasta.cc - last call before finalize!\n");
    }
    if( cdToParent() )
      return -1;
    return timdat.lstep;
}
