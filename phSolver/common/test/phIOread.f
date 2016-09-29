      program readheaderFtn
      use iso_c_binding
      use phio
      use syncio
      use posixio
      use chdir_mod
      include "mpif.h"

      type :: ptrarr
        real(c_double), pointer :: ptr(:,:)
      end type ptrarr

      integer :: rank, ierror, two, nfiles
      type(c_ptr), dimension(2) :: handle
      character(len=30) :: dataDbl, iotype
      character(len=256) :: phrase
      character(len=256), dimension(2) :: dir, fname
      integer, target, dimension(2) :: numpts, ncoords
      real(c_double), allocatable, target :: syncCoords(:,:), posixCoords(:,:)
      type(ptrarr), target, dimension(2) :: coords

      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

      coords(1)%ptr => syncCoords
      coords(2)%ptr => posixCoords

      phrase = c_char_"co-ordinates"//c_null_char
      dataDbl = c_char_"double"//c_null_char
      iotype =  c_char_"binary"//c_null_char
      two = 2 
      nfiles = 2

      dir(1) = c_char_"4-procs_case-SyncIO-2_ref"//c_null_char
      dir(2) = c_char_"4-procs_case-Posix_ref"//c_null_char
      fname(1) = c_char_"geombc-dat."//c_null_char
      fname(2) = c_char_"geombc.dat."//c_null_char
      call syncio_setup_read(nfiles, handle(1))
      call posixio_setup(handle(2), c_char_"r"//c_null_char)
      do i=1,2
        call chdir(dir(i))
        call MPI_Barrier(MPI_COMM_WORLD, ierror)
        call phio_openfile(fname(i), handle(i))
        call phio_readheader(handle(i), phrase, c_loc(numpts), 
     &      two, dataDbl, iotype)
        ncoords(i) = numpts(1)*numpts(2)
        allocate( coords(i)%ptr(numpts(1),numpts(2)) )
        call phio_readdatablock(handle(i), phrase, 
     &      c_loc(coords(i)%ptr(1,1)), ncoords(i), dataDbl, iotype)
        call phio_closefile(handle(i))
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
