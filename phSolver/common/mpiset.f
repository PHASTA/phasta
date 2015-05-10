      subroutine mpiset

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      logical     reorder

      call MPI_COMM_RANK (MPI_COMM_WORLD, myrank)  ! we ditched the ierr that fortran 
	                                             ! normally have here to pacify Digital

      return
      end

