      program readheaderFtn
      use iso_c_binding
      use phio
      use syncio
      use posixio
      include "mpif.h"

      integer, target :: rank, ierror, one, ppf, peers, fish, nfiles
      type(c_ptr), dimension(2) :: handle
      character(len=30) :: dataDbl, iotype
      character(len=256) :: phrase
      character(len=256), dimension(2) :: fname
      integer, target, dimension(2) :: fishweight, numFish

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
      nfiles = 2
      ppf = peers/nfiles
      ! handle(1) is the file for syncio writing
      call syncio_setup_write(nfiles, one, ppf, handle(1))
      ! handle(2) is the file for posix writing
      call posixio_setup(handle(2), c_char_"w"//c_null_char)
      ! if there were a handle(3) for streams we would do the following
      ! call streamio_setup_write(handle(3), <stream obj>)
      ! after the handles are setup the function calls are the same
      ! write the same garbage to posix and syncio files
      do i=1,2
        call phio_openfile(fname(i), handle(i))
        call phio_writeheader(handle(i), phrase, c_loc(fish), one, one,
     &      dataDbl, iotype)
        call phio_writedatablock(handle(i), phrase, c_loc(fishweight(i)),
     &      one, dataDbl, iotype)
        call phio_closefile(handle(i))
      end do
      ! this is the read side... less interesting for us
      call syncio_setup_read(nfiles, handle(1))
      call posixio_setup(handle(2), c_char_"r"//c_null_char)
      do i=1,2
        call phio_openfile(fname(i), handle(i))
        call phio_readheader(handle(i), phrase, c_loc(numFish(i)),
     &      one, dataDbl, iotype)
        call phio_readdatablock(handle(i), phrase, c_loc(fishweight(i)),
     &      one, dataDbl, iotype)
        call phio_closefile(handle(i))
      end do
      if( numFish(1) .ne. numFish(2) .or.
     &    fishweight(1) .ne. fishweight(2) ) then
        write (*,*) "fish don\'t match"
        stop 1
      end if
      call MPI_Finalize(ierror)
      end
