      subroutine rerun_check(stopjob)

      include "common.h"
      include "mpif.h"

      integer stopjob

      stopjob = -1

      if(myrank.eq.master) then
        open(unit=772,file='rerun-check',status='old',iostat=ierr)
        if(ierr.eq.0) read(772,*)stopjob
        close(772)
      endif

      call MPI_BCAST(stopjob,1,MPI_INTEGER,master,
     .               MPI_COMM_WORLD,ierr)

      return
      end
