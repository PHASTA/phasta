#include<phasta.h>
#include<mpi.h>

int 
main( int argc,   
      char* argv[] ) {
  MPI_Init(&argc,&argv);
  phasta ( argc, argv);  
  MPI_Finalize();
  return 0;
}
