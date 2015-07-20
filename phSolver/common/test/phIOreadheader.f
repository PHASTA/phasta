      program readheaderFtn
      use iso_c_binding
      use phio
      use chdir_mod
      include "mpif.h"

      integer :: rank, ierror, one, nfiles, numnp
      type(c_ptr) :: handle
      character(len=30) :: dataInt, iotype
      character(len=256) :: dir, fname

      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

      dataInt = c_char_"integer"//c_null_char
      iotype =  c_char_"binary"//c_null_char
      one = 1

      dir = c_char_"4-procs_case-SyncIO-2"//c_null_char
      call chdir(dir)
      call MPI_Barrier(MPI_COMM_WORLD, ierror)
      if( rank .eq. 0 ) then
        write (*,*) 'sync'
      endif
      nfiles = 2
      fname = c_char_"geombc-dat."//c_null_char
      call phio_openfile_read(fname, nfiles, handle)
      call phio_readheader(handle,
     &  c_char_"number of nodes"//char(0),
     &  c_loc(numnp), one, dataInt, iotype)
      call phio_closefile_read(handle)
      call chdir(c_char_'..'//c_null_char)

      dir = c_char_"4-procs_case-Posix"//c_null_char
      call chdir(dir)
      call MPI_Barrier(MPI_COMM_WORLD, ierror)
      if( rank .eq. 0 ) then
        write (*,*) 'posix'
      endif
      nfiles = 1
      fname = c_char_"geombc.dat."//c_null_char
      call phio_openfile_read(fname, nfiles, handle)
      call phio_readheader(handle,
     &  c_char_"number of nodes"//char(0),
     &  c_loc(numnp), one, dataInt, iotype)
      call phio_closefile_read(handle)
      call chdir(c_char_'..'//c_null_char)

      call MPI_Finalize(ierror)
      end
