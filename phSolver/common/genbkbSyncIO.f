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
c
        include "common.h"
        include "mpif.h" !Required to determine the max for itpblk
c

        integer, allocatable :: ientp(:,:),iBCBtp(:,:)
        real*8, allocatable :: BCBtp(:,:)
        integer materb(ibksz)
        integer intfromfile(50) ! integers read from headers
        character*255 fname1

cccccccccccccc New Phasta IO starts here cccccccccccccccccccccccccccccc

        integer :: descriptor, descriptorG, GPID, color, nfiles, nfields
        integer :: numparts, nppf, nppp, nprocs, writeLock
        integer :: ierr_io, numprocs, itmp, itmp2 
        integer :: itpblktot,ierr

        character*255 fnamer, fname2, temp2
        character*64 temp1, temp3

        nfiles = nsynciofiles
        nfields = nsynciofieldsreadgeombc
        numparts = numpe !This is the common settings. Beware if you try to compute several parts per process

        nppp = numparts/numpe
        nppf = numparts/nfiles

        color = int(myrank/(numparts/nfiles))
        itmp2 = int(log10(float(color+1)))+1
        write (temp2,"('(''geombc-dat.'',i',i1,')')") itmp2
        temp2=trim(temp2)
        write (fnamer,temp2) (color+1)
        fnamer=trim(fnamer)

        ione=1
        itwo=2
        ieight=8
        ieleven=11
        itmp = int(log10(float(myrank+1)))+1


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        iel=1
        itpblk=nelblb


!MR CHANGE
        ! Get the total number of different interior topologies in the whole domain. 
        ! Try to read from a field. If the field does not exist, scan the geombc file.
        itpblktot=-1
        write(temp1,
     &   "('(''total number of boundary tpblocks@'',i',i1,',A1)')") itmp
        write (fname2,temp1) (myrank+1),'?'
        call readheader(igeom,fname2 // char(0) ,itpblktot,ione,
     &  'integer' // char(0),iotype) 

!        write (*,*) 'Rank: ',myrank,' boundary itpblktot intermediate:',
!     &               itpblktot

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
            write (temp1,"('connectivity boundary',i1)") iblk
            temp1 = trim(temp1)
            write (temp3,"('(''@'',i',i1,',A1)')") itmp
            write (fname2, temp3) (myrank+1), '?'
            fname2 = trim(temp1)//trim(fname2)
            !write(*,*) 'rank, fname2',myrank, trim(adjustl(fname2))
            call readheader(igeom,fname2 // char(0),intfromfile,
     &      ieight,'integer' // char(0),iotype)
            neltp = intfromfile(1) ! -1 if fname2 was not found, >=0 otherwise
          end do
          itpblktot = iblk-1   
        end if

        if (myrank == 0) then
          write(*,*) 'Number of boundary topologies: ',itpblktot
        endif
!        write (*,*) 'Rank: ',myrank,' boundary itpblktot final:',
!     &               itpblktot

        nelblb=0
        mattyp=0
        ndofl = ndof
        do iblk = 1, itpblktot
           writeLock=0;

           fname1='connectivity boundary?'

!           print *, "Loop ",iblk, myrank, itpblk, trim(fnamer)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           write (temp1,"('connectivity boundary',i1)") iblk
           temp1 = trim(temp1)
           write (temp3,"('(''@'',i',i1,',A1)')") itmp
           write (fname2, temp3) (myrank+1), '?'
           fname2 = trim(temp1)//trim(fname2)
           fname2 = trim(fname2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           ! Synchronization for performance monitoring, as some parts do not include some topologies
           call MPI_Barrier(MPI_COMM_WORLD,ierr) 
           call readheader(igeom,fname2 // char(0),intfromfile,ieight,
     &                     'integer' // char(0),iotype)
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

!           print *, "neltp is ", neltp

           call readdatablock(igeom,fname2 // char(0),ientp,iientpsiz,
     &                     'integer' // char(0),iotype)


c     
c.... Read the boundary flux codes
c
     


           fname1='nbc codes?'

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           write (temp1,"('nbc codes',i1)") iblk
           temp1=trim(temp1)
           write (temp3,"('(''@'',i',i1,',A1)')") itmp
           write (fname2, temp3) (myrank+1), '?'
           fname2 = trim(temp1)//trim(fname2)
           call MPI_BARRIER(MPI_COMM_WORLD, ierr)
           
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           call readheader(igeom,fname2 // char(0) ,intfromfile,
     &      ieight,'integer' // char(0),iotype)
           iiBCBtpsiz=neltp*ndiBCB
           call readdatablock(igeom,fname2 // char(0) ,iBCBtp,
     &      iiBCBtpsiz,'integer' // char(0),iotype)

c     
c.... read the boundary condition data
c     
           fname1='nbc values?'
           
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           write (temp1,"('nbc values',i1)") iblk
           temp1=trim(temp1)
           write (temp3,"('(''@'',i',i1,',A1)')") itmp
           write (fname2, temp3) (myrank+1), '?'
           fname2 = trim(temp1)//trim(fname2)
           call MPI_BARRIER(MPI_COMM_WORLD, ierr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           call readheader(igeom,fname2 // char(0) ,intfromfile,
     &      ieight,'integer' // char(0) ,iotype)
           BCBtp    = zero
           iBCBtpsiz=neltp*ndBCB
           call readdatablock(igeom,fname2 // char(0),BCBtp,iBCBtpsiz,
     &                     'double' // char(0) ,iotype)


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
c     
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
c
                 allocate (miBCB(nelblb)%p(npro,ndiBCB))
c 
                 allocate (mBCB(nelblb)%p(npro,nshlb,ndBCB))
c 
                 allocate (mmatb(nelblb)%p(npro))
c 
c.... save the boundary element block
c 
                 call gensvb (ientp(n1:n2,1:nshl),
     &                iBCBtp(n1:n2,:),      BCBtp(n1:n2,:),
     &                materb,        mienb(nelblb)%p,
     &                miBCB(nelblb)%p,        mBCB(nelblb)%p,
     &                mmatb(nelblb)%p)
c 
                 iel=iel+npro
              enddo

           endif
           deallocate(ientp)
           deallocate(iBCBtp)
           deallocate(BCBtp)

        enddo
        lcblkb(1,nelblb+1) = iel

c
c.... return
c
        return
c
c.... end of file error handling
c
 911    call error ('genbcb  ','end file',igeomBAK)
c
1000    format(a80,//,
     &  ' B o u n d a r y   E l e m e n t   C o n n e c t i v i t y',//,
     &  '   Elem   BC codes',/,
     &  '  Number  C P V H ',5x,27('Node',i1,:,2x))
1100    format(2x,i5,2x,4i2,3x,27i7)
c$$$2000    format(a80,//,
c$$$     &  ' B o u n d a r y   E l e m e n t   B C   D a t a ',//,
c$$$     &  '   Node   ',3x,'mass',/,
c$$$     &  '  Number  ',3x,'flux',6x,'Pressure',6x,'Heat',6x,
c$$$     &  3('Viscous',i1,:,4x))
2100    format(2x,i5,1p,1x,6e12.4)
c
        end

