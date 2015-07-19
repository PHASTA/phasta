      program readwrite
      use iso_c_binding
      include "common.h"
      include "mpif.h"

      integer :: rank, ierror, ione, nfiles
      type(c_ptr), TARGET :: igeom
      type(c_ptr) :: p

      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

      ione = 1
      nfiles = 2
      p = c_loc(igeom)
      write (*,*) 'rank numfiles', rank, numfiles
      call phio_openfile_read('geombc-dat.' // char(0), nfiles, p)
      write (*,*) 'rank numfiles', rank, numfiles
      !      call phio_readheader(igeom,'number of nodes' // char(0),numnp,ione,
      !     & 'integer' // char(0), iotype)
      write (*,*) rank, ' calling closefile_read'

      call phio_closefile_read(igeom)


      call MPI_Finalize(ierror)
      end
