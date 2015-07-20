      program readheaderFtn
      use iso_c_binding
      use phio
      use chdir_mod
      include "mpif.h"

      integer :: rank, ierror, one
      type(c_ptr) :: handle
      character(len=30) :: dataInt, iotype
      character(len=256), dimension(2) :: dir, fname
      integer, dimension(2) :: nfiles, numnp

      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

      dataInt = c_char_"integer"//c_null_char
      iotype =  c_char_"binary"//c_null_char
      one = 1

      dir(1) = c_char_"4-procs_case-SyncIO-2"//c_null_char
      dir(2) = c_char_"4-procs_case-Posix"//c_null_char
      fname(1) = c_char_"geombc-dat."//c_null_char
      fname(2) = c_char_"geombc.dat."//c_null_char
      nfiles(1) = 2
      nfiles(2) = 1
      do i=1,2
        call chdir(dir(i))
        call MPI_Barrier(MPI_COMM_WORLD, ierror)
        nfiles = 2
        call phio_openfile_read(fname(i), nfiles(i), handle)
        call phio_readheader(handle,
     &      c_char_"number of nodes"//char(0),
     &      c_loc(numnp(i)), one, dataInt, iotype)
        call phio_closefile_read(handle)
        call chdir(c_char_'..'//c_null_char)
      end do
      if( numnp(1) .ne. numnp(2) ) then
        write (*,*) 'rank numnp', rank, numnp
        stop 1  
      endif
      call MPI_Finalize(ierror)
      end
