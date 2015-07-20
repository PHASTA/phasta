      program readwrite
      use iso_c_binding
      use phio
      include "mpif.h"

      integer :: rank, ierror, ione, nfiles, numnp
      type(c_ptr) :: handle
      character(len=30) :: dataInt, iotype
      dataInt = c_char_"integer"//c_null_char
      iotype =  c_char_"binary"//c_null_char

      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

      ione = 1
      nfiles = 2
      call phio_openfile_read(c_char_"geombc-dat."//c_null_char, nfiles, handle)
      call phio_readheader(handle,
     &  c_char_"number of nodes"//char(0),
     &  c_loc(numnp), ione, dataInt, iotype)
      write (*,*) rank, ' calling closefile_read'

      call phio_closefile_read(handle)

      call MPI_Finalize(ierror)
      end
