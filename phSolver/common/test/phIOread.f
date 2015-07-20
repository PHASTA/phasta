      program readheaderFtn
      use iso_c_binding
      use phio
      use chdir_mod
      include "mpif.h"

      type :: ptrarr
        real(c_double), pointer :: ptr(:,:)
      end type ptrarr

      integer :: rank, ierror, two
      type(c_ptr) :: handle
      character(len=30) :: dataDbl, iotype
      character(len=256) :: phrase
      character(len=256), dimension(2) :: dir, fname
      integer, dimension(2) :: nfiles, numpts, ncoords
      real(c_double), allocatable, target :: syncCoords(:,:), posixCoords(:,:)
      type(ptrarr), dimension(2) :: coords

      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

      coords(1)%ptr => syncCoords
      coords(2)%ptr => posixCoords

      phrase = c_char_"co-ordinates"//c_null_char
      dataDbl = c_char_"double"//c_null_char
      iotype =  c_char_"binary"//c_null_char
      two = 2 

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
        call phio_readheader(handle, phrase, c_loc(numpts), 
     &      two, dataDbl, iotype)
        ncoords(i) = numpts(1)*numpts(2)
        allocate( coords(i)%ptr(numpts(1),numpts(2)) )
        call phio_readdatablock(handle, phrase, 
     &      c_loc(coords(i)%ptr), ncoords(i), dataDbl, iotype)
        call phio_closefile_read(handle)
        call chdir(c_char_'..'//c_null_char)
      end do
      if( ncoords(1) .ne. ncoords(2) ) then
        write (*,*) 'rank ncoords', rank, ncoords
        stop 1  
      endif
      do i=1,numpts(1)
        do j=1,numpts(2)
          if( coords(1)%ptr(i,j) .ne. coords(2)%ptr(i,j) ) then
            write (*,*) 'rank coordinate mismatch i,j', rank, i, j
            stop 1
          end if
        end do
      end do
      deallocate(coords(1)%ptr)
      deallocate(coords(2)%ptr)
      call MPI_Finalize(ierror)
      end
