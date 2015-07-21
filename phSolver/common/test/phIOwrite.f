      program readheaderFtn
      use iso_c_binding
      use phio
      include "mpif.h"

      integer, target :: rank, ierror, one, ppf, peers, fish
      type(c_ptr) :: handle
      character(len=30) :: dataDbl, iotype
      character(len=256) :: phrase
      character(len=256), dimension(2) :: fname
      integer, target, dimension(2) :: nfiles, fishweight, numFish

      call MPI_Init(ierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
      call MPI_Comm_size(MPI_COMM_WORLD, peers, ierror)

      phrase = c_char_"number of fishes"//c_null_char
      dataDbl = c_char_"double"//c_null_char
      iotype =  c_char_"binary"//c_null_char
      one = 1 
      fish = 2

      fishweight(1) = 1.23
      fishweight(2) = 1.23

      fname(1) = c_char_"fortranWater-dat."//c_null_char
      fname(2) = c_char_"fortranWater.dat."//c_null_char
      nfiles(1) = 2
      nfiles(2) = 1
      ppf = peers/nfiles(1)
      do i=1,2
        call phio_openfile_write(fname(i), nfiles(i), one, ppf, handle)
        call phio_writeheader(handle, phrase, c_loc(fish), one, one,
     &      dataDbl, iotype)
        call phio_writedatablock(handle, phrase, c_loc(fishweight(i)),
     &      one, dataDbl, iotype)
        call phio_closefile_write(handle)
      end do
      do i=1,2
        call phio_openfile_read(fname(i), nfiles(i), handle)
        call phio_readheader(handle, phrase, c_loc(numFish(i)),
     &      one, dataDbl, iotype)
        call phio_readdatablock(handle, phrase, c_loc(fishweight(i)),
     &      one, dataDbl, iotype)
        call phio_closefile_read(handle)
      end do
      if( numFish(1) .ne. numFish(2) .or.
     &    fishweight(1) .ne. fishweight(2) ) then
        write (*,*) "fish don\'t match"
        stop 1
      end if
      call MPI_Finalize(ierror)
      end
