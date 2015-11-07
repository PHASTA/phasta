        subroutine genbkb (ibksz)
c
c----------------------------------------------------------------------
c
c  This routine reads the boundary elements, reorders them and
c  generates traces for the gather/scatter operations.
c
c Zdenek Johan, Fall 1991.
c----------------------------------------------------------------------
c
        use dtnmod
        use pointer_data
        use phio
        use iso_c_binding
        include "common.h"
        include "mpif.h" !Required to determine the max for itpblk

        integer, target, allocatable :: ientp(:,:),iBCBtp(:,:)
        real*8, target, allocatable :: BCBtp(:,:)
        integer materb(ibksz)
        integer, target :: intfromfile(50) ! integers read from headers
        character*255 fname1
        integer :: descriptor, descriptorG, GPID, color, nfields
        integer :: numparts, nppp, nprocs, writeLock
        integer :: ierr_io, numprocs, itmp, itmp2 
        integer, target :: itpblktot,ierr
        character*255 fname2
        character(len=30) :: dataInt, dataDbl
        dataInt = c_char_'integer'//c_null_char
        dataDbl = c_char_'double'//c_null_char

        nfields = nsynciofieldsreadgeombc
        numparts = numpe !This is the common settings. Beware if you try to compute several parts per process
        nppp = numparts/numpe
        ione=1
        itwo=2
        ieight=8
        ieleven=11
        itmp = int(log10(float(myrank+1)))+1
        iel=1
        itpblk=nelblb

        ! Get the total number of different interior topologies in the whole domain. 
        ! Try to read from a field. If the field does not exist, scan the geombc file.
        itpblktot=1  ! hardwired to monotopology for now
        call phio_readheader(fhandle,
     &   c_char_'total number of boundary tpblocks' // char(0),
     &   c_loc(itpblktot), ione, dataInt, iotype)

        if (itpblktot == -1) then 
          ! The field 'total number of different boundary tpblocks' was not found in the geombc file.
          ! Scan all the geombc file for the 'connectivity interior' fields to get this information.
          iblk=0
          neltp=0
          do while(neltp .ne. -1) 

            ! intfromfile is reinitialized to -1 every time.
            ! If connectivity boundary@xxx is not found, then 
            ! readheader will return intfromfile unchanged

            intfromfile(:)=-1
            iblk = iblk+1
            if(input_mode.ge.1)then
               write (fname2,"('connectivity boundary',i1)") iblk
            else
               write (fname2,"('connectivity boundary linear tetrahedron')")
            endif
 
            call phio_readheader(fhandle, fname2 // char(0),
     &       c_loc(intfromfile), ieight, dataInt, iotype)
            neltp = intfromfile(1) ! -1 if fname2 was not found, >=0 otherwise
          end do
          itpblktot = iblk-1   
        end if

        if (myrank == 0) then
          write(*,*) 'Number of boundary topologies: ',itpblktot
        endif

        nelblb=0
        mattyp=0
        ndofl = ndof

        do iblk = 1, itpblktot
           writeLock=0;
            if(input_mode.ge.1)then
               write (fname2,"('connectivity boundary',i1)") iblk
            else
               write (fname2,"('connectivity boundary linear tetrahedron')")
            endif
           
           ! Synchronization for performance monitoring, as some parts do not include some topologies
           call MPI_Barrier(MPI_COMM_WORLD,ierr) 
           call phio_readheader(fhandle, fname2 // char(0),
     &      c_loc(intfromfile), ieight, dataInt, iotype)
           neltp =intfromfile(1)
           nenl  =intfromfile(2)
           ipordl=intfromfile(3)
           nshl  =intfromfile(4)
           nshlb =intfromfile(5)
           nenbl =intfromfile(6)
           lcsyst=intfromfile(7)
           numnbc=intfromfile(8)

           allocate (ientp(neltp,nshl))
           allocate (iBCBtp(neltp,ndiBCB))
           allocate (BCBtp(neltp,ndBCB))
           iientpsiz=neltp*nshl

           if (neltp==0) then
              writeLock=1;
           endif

           call phio_readdatablock(fhandle, fname2 // char(0),
     &      c_loc(ientp),iientpsiz,dataInt,iotype)
c     
c.... Read the boundary flux codes
c
           call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            if(input_mode.ge.1)then
               write (fname2,"('nbc codes',i1)") iblk
            else
               write (fname2,"('nbc codes linear tetrahedron')")
            endif

           call phio_readheader(fhandle, fname2 // char(0),
     &      c_loc(intfromfile), ieight, dataInt, iotype)
           iiBCBtpsiz=neltp*ndiBCB
           call phio_readdatablock(fhandle, fname2 // char(0),
     &      c_loc(iBCBtp),iiBCBtpsiz,dataInt,iotype)
c     
c.... read the boundary condition data
c     
           call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            if(input_mode.ge.1)then
               write (fname2,"('nbc values',i1)") iblk
            else
               write (fname2,"('nbc values linear tetrahedron')")
            endif

           call phio_readheader(fhandle, fname2 // char(0),
     &      c_loc(intfromfile), ieight, dataInt, iotype)
           BCBtp    = zero
           iBCBtpsiz=neltp*ndBCB
           call phio_readdatablock(fhandle, fname2 // char(0),
     &      c_loc(BCBtp),iBCBtpsiz,dataDbl,iotype)
c
c This is a temporary fix until NSpre properly zeros this array where it
c is not set.  DEC has indigestion with these arrays though the
c result is never used (never effects solution).
c
           if(writeLock==0) then
              where(.not.btest(iBCBtp(:,1),0)) BCBtp(:,1)=zero
              where(.not.btest(iBCBtp(:,1),1)) BCBtp(:,2)=zero
              where(.not.btest(iBCBtp(:,1),3)) BCBtp(:,6)=zero
              if(ndBCB.gt.6) then
                 do i=6,ndof
                    where(.not.btest(iBCBtp(:,1),i-1)) BCBtp(:,i+1)=zero
                 enddo
              endif
              where(.not.btest(iBCBtp(:,1),2)) 
                 BCBtp(:,3)=zero
                 BCBtp(:,4)=zero
                 BCBtp(:,5)=zero
              endwhere
              
              do n=1,neltp,ibksz 
                 nelblb=nelblb+1
                 npro= min(IBKSZ, neltp - n + 1)
                 lcblkb(1,nelblb)  = iel
                 lcblkb(3,nelblb)  = lcsyst
                 lcblkb(4,nelblb)  = ipordl
                 lcblkb(5,nelblb)  = nenl
                 lcblkb(6,nelblb)  = nenbl
                 lcblkb(7,nelblb)  = mattyp
                 lcblkb(8,nelblb)  = ndofl
                 lcblkb(9,nelblb)  = nshl 
                 lcblkb(10,nelblb) = nshlb ! # of shape functions per elt
c     
c.... save the element block
c     
                 n1=n
                 n2=n+npro-1
                 materb=1       ! all one material for now
c     
c.... allocate memory for stack arrays
c
                 allocate (mienb(nelblb)%p(npro,nshl))
                 allocate (miBCB(nelblb)%p(npro,ndiBCB))
                 allocate (mBCB(nelblb)%p(npro,nshlb,ndBCB))
                 allocate (mmatb(nelblb)%p(npro))
c 
c.... save the boundary element block
c 
                 call gensvb (ientp(n1:n2,1:nshl),
     &                iBCBtp(n1:n2,:),      BCBtp(n1:n2,:),
     &                materb,        mienb(nelblb)%p,
     &                miBCB(nelblb)%p,        mBCB(nelblb)%p,
     &                mmatb(nelblb)%p)
                 iel=iel+npro
              enddo
           endif
           deallocate(ientp)
           deallocate(iBCBtp)
           deallocate(BCBtp)

        enddo
        lcblkb(1,nelblb+1) = iel
        return
c
c.... end of file error handling
c
 911    call error ('genbcb  ','end file',igeomBAK)
1000    format(a80,//,
     &  ' B o u n d a r y   E l e m e n t   C o n n e c t i v i t y',//,
     &  '   Elem   BC codes',/,
     &  '  Number  C P V H ',5x,27('Node',i1,:,2x))
1100    format(2x,i5,2x,4i2,3x,27i7)
2100    format(2x,i5,1p,1x,6e12.4)
        end
