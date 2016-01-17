        subroutine genblk (IBKSZ)
c
c----------------------------------------------------------------------
c
c  This routine reads the interior elements and generates the
c  appropriate blocks.
c
c Zdenek Johan, Fall 1991.
c----------------------------------------------------------------------
c
        use pointer_data
        use phio
        use iso_c_binding
        include "common.h"
        include "mpif.h" !Required to determine the max for itpblk

        integer, target, allocatable :: ientp(:,:)
        integer mater(ibksz)
        integer, target :: intfromfile(50) ! integers read from headers
        character*255 fname1
        integer :: descriptor, descriptorG, GPID, color
        integer ::  numparts, writeLock
        integer :: ierr_io, numprocs
        integer, target :: itpblktot,ierr,iseven
        character*255 fname2
        character(len=30) :: dataInt
        dataInt = c_char_'integer'//c_null_char
        numparts = numpe !This is the common settings. Beware if you try to compute several parts per process
        ione=1
        itwo=2
        iseven=7
        ieleven=11
        iel=1
        itpblk=nelblk

        ! Get the total number of different interior topologies in the whole domain. 
        ! Try to read from a field. If the field does not exist, scan the geombc file.
          itpblktot=1  ! hardwired to montopology for now
        call phio_readheader(fhandle,
     &   c_char_'total number of interior tpblocks' // char(0),
     &   c_loc(itpblktot), ione, dataInt, iotype) 

        if (itpblktot == -1) then 
          ! The field 'total number of different interior tpblocks' was not found in the geombc file.
          ! Scan all the geombc file for the 'connectivity interior' fields to get this information.
          iblk=0
          neltp=0
          do while(neltp .ne. -1) 

            ! intfromfile is reinitialized to -1 every time.
            ! If connectivity interior@xxx is not found, then 
            ! readheader will return intfromfile unchanged

            intfromfile(:)=-1
            iblk = iblk+1
            if(input_mode.ge.1) then
              write (fname2,"('connectivity interior',i1)") iblk
            else
              write (fname2,"('connectivity interior linear tetrahedron')") 
!              write (fname2,"('connectivity interior?')") 
            endif

            !write(*,*) 'rank, fname2',myrank, trim(adjustl(fname2))
            call phio_readheader(fhandle, fname2 // char(0),
     &       c_loc(intfromfile), iseven, dataInt, iotype)
            neltp = intfromfile(1) ! -1 if fname2 was not found, >=0 otherwise
          end do
          itpblktot = iblk-1   
        end if

        if (myrank == 0) then
          write(*,*) 'Number of interior topologies: ',itpblktot
        endif

        nelblk=0
        mattyp = 0
        ndofl = ndof
        nsymdl = nsymdf

        do iblk = 1, itpblktot
           writeLock=0;
            if(input_mode.ge.1) then
              write (fname2,"('connectivity interior',i1)") iblk
            else
              write (fname2,"('connectivity interior linear tetrahedron')") 
!              write (fname2,"('connectivity interior?')") 
            endif

           ! Synchronization for performance monitoring, as some parts do not include some topologies
           call MPI_Barrier(MPI_COMM_WORLD,ierr) 
           call phio_readheader(fhandle, fname2 // char(0),
     &      c_loc(intfromfile), iseven, dataInt, iotype)
           neltp  =intfromfile(1)
           nenl   =intfromfile(2)
           ipordl =intfromfile(3)
           nshl   =intfromfile(4)
           ijunk  =intfromfile(5)
           ijunk  =intfromfile(6)
           lcsyst =intfromfile(7)
           allocate (ientp(neltp,nshl))
           iientpsiz=neltp*nshl

           if (neltp==0) then
              writeLock=1;
           endif

           call phio_readdatablock(fhandle,fname2 // char(0),
     &      c_loc(ientp), iientpsiz, dataInt, iotype)

           if(writeLock==0) then
             do n=1,neltp,ibksz 
                nelblk=nelblk+1
                npro= min(IBKSZ, neltp - n + 1)
                lcblk(1,nelblk)  = iel
                lcblk(3,nelblk)  = lcsyst
                lcblk(4,nelblk)  = ipordl
                lcblk(5,nelblk)  = nenl
                lcblk(6,nelblk)  = nfacel
                lcblk(7,nelblk)  = mattyp
                lcblk(8,nelblk)  = ndofl
                lcblk(9,nelblk)  = nsymdl 
                lcblk(10,nelblk) = nshl ! # of shape functions per elt
c
c.... allocate memory for stack arrays
c
                allocate (mmat(nelblk)%p(npro))
c
                allocate (mien(nelblk)%p(npro,nshl))
                allocate (mxmudmi(nelblk)%p(npro,maxsh))
                if(usingpetsc.eq.0) then
                    allocate (mienG(nelblk)%p(1,1))
                else
                    allocate (mienG(nelblk)%p(npro,nshl))
                endif
                ! note mienG will be passed to gensav but nothing filled if not 
                ! using PETSc so this is safe
c
c.... save the element block
c
                n1=n
                n2=n+npro-1
                mater=1   ! all one material for now
                call gensav (ientp(n1:n2,1:nshl),
     &                       mater,           mien(nelblk)%p,
     &                       mienG(nelblk)%p,
     &                       mmat(nelblk)%p)
                iel=iel+npro
             enddo
           endif
           deallocate(ientp)
        enddo

        lcblk(1,nelblk+1) = iel
        return
1000    format(a80,//,
     &  ' N o d a l   C o n n e c t i v i t y',//,
     &  '   Elem  ',/,
     &  '  Number  ',7x,27('Node',i2,:,2x))
1100    format(2x,i5,6x,27i8)
        end
