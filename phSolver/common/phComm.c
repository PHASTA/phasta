#include "phComm.h"
#include <mpi.h>

int phcomm_rank() {
  int r;
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  return r;
}
