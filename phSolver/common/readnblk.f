c  readnblk.f (pronounce "Reed and Block Dot Eff") contains:
c
c    module readarrays ("Red Arrays") -- contains the arrays that
c     are read in from binary files but not immediately blocked 
c     through pointers.
c
c    subroutine readnblk ("Reed and Block") -- allocates space for
c     and reads data to be contained in module readarrays.  Reads
c     all remaining data and blocks them with pointers.
c


      module readarrays
      
      real*8, allocatable :: point2x(:,:)
      real*8, allocatable :: qold(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: acold(:,:)
      integer, allocatable :: iBCtmp(:)
      real*8, allocatable :: BCinp(:,:)

      integer, allocatable :: point2ilwork(:)
      integer, allocatable :: nBC(:)
      integer, allocatable :: point2iper(:)
      integer, target, allocatable :: point2ifath(:)
      integer, target, allocatable :: point2nsons(:)
      
      end module



      subroutine readnblk
c
      use iso_c_binding 
      use readarrays
      use phio
      include "common.h"
c
      real*8, target, allocatable :: xread(:,:), qread(:,:), acread(:,:)
      real*8, target, allocatable :: uread(:,:)
      real*8, target, allocatable :: BCinpread(:,:)
      integer, target, allocatable :: iperread(:), iBCtmpread(:)
      integer, target, allocatable :: ilworkread(:), nBCread(:)
      character*10 cname2
      character*8 mach2
!MR CHANGE
!      character*20 fmt1
      character*30 fmt1
!MR CHANGE END
      character*255 fname1,fnamer,fnamelr
      character*255 warning
      integer igeomBAK, ibndc, irstin, ierr
      integer, target :: intfromfile(50) ! integers read from headers

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc New PhastaIO Definition Part ccccccccccccccccccccccccccccccccccccccccc

      integer :: descriptor, descriptorG, GPID, color, nfiles, nfields
      integer ::  numparts, nppf
      integer :: ierr_io, numprocs, itmp, itmp2
      integer :: ignored
      character*255 fname2, temp2
      character*64 temp1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      type(c_ptr) :: handle
      character(len=1024) :: dataInt, dataDbl
      dataInt = c_char_'integer'//c_null_char
      dataDbl = c_char_'double'//c_null_char

c
c.... determine the step number to start with
c
      open(unit=72,file='numstart.dat',status='old')
      read(72,*) irstart
      close(72)
      lstep=irstart ! in case restart files have no fields

c
      fname1='geombc.dat'
      fname1= trim(fname1)  // cname2(myrank+1)
c      fnamelr='restart.latest'

      itmp=1
      if (irstart .gt. 0) itmp = int(log10(float(irstart)))+1
      write (fmt1,"('(''restart.'',i',i1,',1x)')") itmp
      write (fnamer,fmt1) irstart
      fnamer = trim(fnamer) // cname2(myrank+1)
c      fnamelr = trim(fnamelr) // cname2(myrank+1)

c
c.... open input files
c


c      call openfile(  fname1,  'read?', igeomBAK );


c
c.... try opening restart.latest.proc before trying restart.stepno.proc
c
c      call openfile(  fnamelr,  'read?', irstin );
c      if ( irstin .eq. 0 ) 

!MR CHANGE
!       call openfile( fnamer, 'read?', irstin );
!MR CHANGE END

! either one will work
c
c.... input the geometry parameters
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!MR CHANGE

!      nfiles = 2
!      nfields = 31
!      numparts = 8
!      nppp = numparts/numpe
!      nppf = numparts/nfiles

      nfiles = nsynciofiles
!      nfields = nsynciofieldsreadgeombc
      numparts = numpe !This is the common settings. Beware if you try to compute several parts per process

!      nppp = numparts/numpe
!      nppf = numparts/nfiles
!MR CHANGE END

      itwo=2
      ione=1
      ieleven=11
      itmp = int(log10(float(myrank+1)))+1

      call phio_openfile_read(c_char_'geombc-dat.' // char(0), nfiles, fhandle);

      call phio_readheader(fhandle,c_char_'number of nodes' // char(0),
     & c_loc(numnp),ione, dataInt, iotype)

      call phio_readheader(fhandle,c_char_'number of modes' // char(0),
     & c_loc(nshg),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of interior elements' // char(0),
     &  c_loc(numel),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of boundary elements' // char(0),
     &  c_loc(numelb),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'maximum number of element nodes' // char(0),
     &  c_loc(nen),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of interior tpblocks' // char(0),
     &  c_loc(nelblk),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of boundary tpblocks' // char(0),
     & c_loc(nelblb),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of nodes with Dirichlet BCs' // char(0),
     & c_loc(numpbc),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of shape functions' // char(0),
     & c_loc(ntopsh),ione, dataInt, iotype)

c      call closefile( igeom, "read" )
c      call finalizephmpiio( igeom )

!       if(myrank==0) then
!          print *, "Reading Finished, ********* : ", trim(fname2)
!       endif


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      ieleven=11
c      ione=1
c      fname1='number of nodes?'
c      call readheader(igeomBAK,fname1,numnp,ione,'integer', iotype)
c      fname1='number of modes?'
c      call readheader(igeomBAK,fname1,nshg,ione,'integer', iotype)
cc      fname1='number of global modes?'
cc      call readheader(igeomBAK,fname1,nshgt,ione,'integer', iotype)
c      fname1='number of interior elements?'
c      call readheader(igeomBAK,fname1,numel,ione,'integer', iotype)
c      fname1='number of boundary elements?'
c      call readheader(igeomBAK,fname1,numelb,ione,'integer', iotype)
c      fname1='maximum number of element nodes?'
c      call readheader(igeomBAK,fname1,nen,ione,'integer', iotype)
c      fname1='number of interior tpblocks?'
c      call readheader(igeomBAK,fname1,nelblk,ione,'integer', iotype)
c      fname1='number of boundary tpblocks?'
c      call readheader(igeomBAK,fname1,nelblb,ione,'integer', iotype)
c      fname1='number of nodes with Dirichlet BCs?'
c      call readheader(igeomBAK,fname1,numpbc,ione,'integer', iotype)
c      fname1='number of shape functions?'
c      call readheader(igeomBAK,fname1,ntopsh,ione,'integer', iotype)

c
c.... calculate the maximum number of boundary element nodes
c     
      nenb = 0
      do i = 1, melCat
         if (nen .eq. nenCat(i,nsd)) nenb = max(nenCat(i,nsd-1), nenb)
      enddo
c     
      if (myrank == master) then
         if (nenb .eq. 0) call error ('input   ','nen     ',nen)
      endif
c
c.... setup some useful constants
c
      I3nsd  = nsd / 3          ! nsd=3 integer flag
      E3nsd  = float(I3nsd)     ! nsd=3 real    flag
c    
      if(matflg(1,1).lt.0) then
         nflow = nsd + 1
      else
         nflow = nsd + 2
      endif 
      ndof   = nsd + 2
      nsclr=impl(1)/100
      ndof=ndof+nsclr           ! number of sclr transport equations to solve
      
      ndofBC = ndof + I3nsd     ! dimension of BC array
      ndiBCB = 2                ! dimension of iBCB array
      ndBCB  = ndof + 1         ! dimension of BCB array
c     
      nsymdf = (ndof*(ndof + 1)) / 2 ! symm. d.o.f.'s
c
c.... ----------------------> Communication tasks <--------------------
c

cc still read in new

      if(numpe > 1) then

cc MR CHANGE
         call phio_readheader(fhandle,
     &    c_char_'size of ilwork array' // char(0),
     &    c_loc(nlwork),ione, dataInt, iotype)

         call phio_readheader(fhandle,
     &    c_char_'ilwork' //char(0),
     &    c_loc(nlwork),ione, dataInt, iotype)

         allocate( point2ilwork(nlwork) )
         allocate( ilworkread(nlwork) )
         call phio_readdatablock(fhandle, c_char_'ilwork' // char(0),
     &      c_loc(ilworkread), nlwork, dataInt , iotype)

c      call closefile( igeom, "read" )
c      call finalizephmpiio( igeom )

c         fname1='size of ilwork array?'
c         call readheader(igeomBAK,fname1,nlwork,ione,'integer', iotype)

c         ione=1
c         fname1='ilwork?'
c         call readheader(igeomBAK,fname1,nlwork,ione,'integer', iotype)

c         allocate( point2ilwork(nlwork) )
c         allocate( ilworkread(nlwork) )
c         call readdatablock(igeomBAK,fname1,ilworkread,
c     &                      nlwork,'integer', iotype)
cc MR CHANGE


         point2ilwork = ilworkread
         call ctypes (point2ilwork)
      else
           nlwork=1
           allocate( point2ilwork(1))
           nshg0 = nshg
      endif

cccccccccccccccccccccccccccccccccccccccccc

      itwo=2

c      print *, "fname2 is", fname2

cc MR CHANGE
c      call initphmpiio(nfields,nppf,nfiles,igeom,'read')
c      call openfile( fnamer, 'read', igeom )
CC MR CHANGE

      call phio_readheader(fhandle,
     & c_char_'co-ordinates' // char(0),
     & c_loc(intfromfile),itwo, dataDbl, iotype)
      numnp=intfromfile(1)
c      print *, "read out @@@@@@ is ", numnp
      allocate( point2x(numnp,nsd) )
      allocate( xread(numnp,nsd) )
      ixsiz=numnp*nsd
      call phio_readdatablock(fhandle,
     & c_char_'co-ordinates' // char(0),
     & c_loc(xread),ixsiz, dataDbl, iotype)
      point2x = xread

!      call closefile( igeom, "read" )
!      call finalizephmpiio( igeom )

!       print *, "Finished finalize"

c      deallocate (point2x)
c      deallocate (xread)

cccccccccccccccccccccccccccccccccccccccccc

c     
c.... read the node coordinates
c

c      itwo=2
c      fname1='co-ordinates?'
c      call readheader(igeomBAK,fname1,intfromfile,itwo, 'double', iotype)
c      numnp=intfromfile(1)
cc      nsd=intfromfile(2)
c      allocate( point2x(numnp,nsd) )
c      allocate( xread(numnp,nsd) )
c      ixsiz=numnp*nsd
c      call readdatablock(igeomBAK,fname1,xread,ixsiz, 'double',iotype)
c      point2x = xread

c
c.... read in and block out the connectivity
c

! !MR CHANGE
!     This is not necessary but this avoids to have the geombc files opend two times.
!     A better way consists in pasisng the file handler to genblk or make it global or use igeomBAK instead of igeom
!      call closefile( igeom, "read" )
!      call finalizephmpiio( igeom )
! !MR CHANGE END

      call genblk (IBKSIZ)

! !MR CHANGE
!      call initphmpiio( nfields, nppf, nfiles, igeom, 'read')
!      call openfile( fnamer, 'read', igeom )
! !MR CHANGE END

c
c.... read the boundary condition mapping array
c

cc MR CHANGE
!      call initphmpiio(nfields,nppf,nfiles,igeom, 'read')
!      call openfile( fnamer, 'read', igeom )
cc MR CHANGE

      ione=1
      call phio_readheader(fhandle,
     & c_char_'bc mapping array' // char(0),
     & c_loc(nshg),ione, dataInt, iotype)

c      fname1='bc mapping array?'
c      call readheader(igeomBAK,fname1,nshg,
c     &     ione,'integer', iotype)

      allocate( nBC(nshg) )

      allocate( nBCread(nshg) )

c      call readdatablock(igeomBAK,fname1,nBCread,nshg,'integer',iotype)
      call phio_readdatablock(fhandle,
     & c_char_'bc mapping array' // char(0),
     & c_loc(nBCread), nshg, dataInt, iotype)

      nBC=nBCread

c
c.... read the temporary iBC array
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'bc codes array' // char(0),
     & c_loc(numpbc),ione, dataInt, iotype)

c      ione = 1
c      fname1='bc codes array?'
c      call readheader(igeomBAK,fname1,numpbc,
c     &     ione, 'integer', iotype)

!MR CHANGE
!       if ( numpbc > 0 ) then
!          allocate( iBCtmp(numpbc) )
!          allocate( iBCtmpread(numpbc) )
! c         call readdatablock(igeomBAK,fname1,iBCtmpread,numpbc,
! c     &                      'integer',iotype)
!         call readdatablock(igeom,fname2,iBCtmpread,numpbc,
!      &                      'integer',iotype)
!          iBCtmp=iBCtmpread
!       else  ! sometimes a partition has no BC's
!          allocate( iBCtmp(1) )
!          iBCtmp=0
!       endif

      if ( numpbc > 0 ) then
        allocate( iBCtmp(numpbc) )
        allocate( iBCtmpread(numpbc) )
      else
        allocate( iBCtmp(1) )
        allocate( iBCtmpread(1) )
      endif
c         call readdatablock(igeomBAK,fname1,iBCtmpread,numpbc,
c     &                      'integer',iotype)
      call phio_readdatablock(fhandle,
     & c_char_'bc codes array' // char(0),
     & c_loc(iBCtmpread), numpbc, dataInt, iotype)

      if ( numpbc > 0 ) then
         iBCtmp=iBCtmpread
      else  ! sometimes a partition has no BC's
         deallocate( iBCtmpread)
         iBCtmp=0
      endif
!MR CHANGE END

c
c.... read boundary condition data
c

      ione=1

c      ione=1
c      fname1='boundary condition array?'
c      call readheader(igeomBAK,fname1,intfromfile,
c     &     ione, 'double', iotype)
      call phio_readheader(fhandle,
     & c_char_'boundary condition array' // char(0),
     & c_loc(intfromfile),ione, dataDbl, iotype)
c here intfromfile(1) contains (ndof+7)*numpbc
!MR CHANGE
!       if ( numpbc > 0 ) then
!          if(intfromfile(1).ne.(ndof+7)*numpbc) then
!            warning='WARNING more data in BCinp than needed: keeping 1st'
!            write(*,*) warning, ndof+7
!          endif
!          allocate( BCinp(numpbc,ndof+7) )
!          nsecondrank=intfromfile(1)/numpbc
!          allocate( BCinpread(numpbc,nsecondrank) )
!          iBCinpsiz=intfromfile(1)
! c         call readdatablock(igeomBAK,fname1,BCinpread,iBCinpsiz,
! c     &                      'double',iotype)
!          call readdatablock(igeom,fname2,BCinpread,iBCinpsiz,
!      &                      'double',iotype)
!          BCinp(:,1:(ndof+7))=BCinpread(:,1:(ndof+7))
!       else  ! sometimes a partition has no BC's
!          allocate( BCinp(1,ndof+7) )
!          BCinp=0
!       endif

      if ( numpbc > 0 ) then
!         if(intfromfile(1).ne.(ndof+7)*numpbc) then
!           warning='WARNING more data in BCinp than needed: keeping 1st'
!           write(*,*) warning, ndof+7
!         endif
         allocate( BCinp(numpbc,ndof+7) )
         nsecondrank=intfromfile(1)/numpbc
         allocate( BCinpread(numpbc,nsecondrank) )
         iBCinpsiz=intfromfile(1)
      else
         allocate( BCinp(1,ndof+7) )
         allocate( BCinpread(0,0) ) !dummy
         iBCinpsiz=intfromfile(1)
      endif
c         call readdatablock(igeomBAK,fname1,BCinpread,iBCinpsiz,
c     &                      'double',iotype)

      call phio_readdatablock(fhandle,
     & c_char_'boundary condition array' // char(0),
     & c_loc(BCinpread), iBCinpsiz, dataDbl, iotype)

      if ( numpbc > 0 ) then
         BCinp(:,1:(ndof+7))=BCinpread(:,1:(ndof+7))
      else  ! sometimes a partition has no BC's
         deallocate(BCinpread)
         BCinp=0
      endif
!MR CHANGE END

c
c.... read periodic boundary conditions
c

      ione=1
c      fname1='periodic masters array?'
c      call readheader(igeomBAK,fname1,nshg,
c     &     ione, 'integer', iotype)
      call phio_readheader(fhandle,
     & c_char_'periodic masters array' // char(0),
     & c_loc(nshg), ione, dataInt, iotype)
      allocate( point2iper(nshg) )
      allocate( iperread(nshg) )
c      call readdatablock(igeomBAK,fname1,iperread,nshg,
c     &                      'integer',iotype)
      call phio_readdatablock(fhandle,
     & c_char_'periodic masters array' // char(0),
     & c_loc(iperread), nshg, dataInt, iotype)
      point2iper=iperread


! !MR CHANGE
!      call closefile( igeom, "read" )
!      call finalizephmpiio( igeom )
! !MR CHANGE END

c
c.... generate the boundary element blocks
c
      call genbkb (ibksiz)


! !MR CHANGE
!       write(*,*) 'HELLO 12 from ', myrank
! !MR CHANGE END

c
c  Read in the nsons and ifath arrays if needed
c
c  There is a fundamental shift in the meaning of ifath based on whether
c  there exist homogenous directions in the flow.  
c
c  HOMOGENOUS DIRECTIONS EXIST:  Here nfath is the number of inhomogenous
c  points in the TOTAL mesh.  That is to say that each partition keeps a 
c  link to  ALL inhomogenous points.  This link is furthermore not to the
c  sms numbering but to the original structured grid numbering.  These 
c  inhomogenous points are thought of as fathers, with their sons being all
c  the points in the homogenous directions that have this father's 
c  inhomogeneity.  The array ifath takes as an arguement the sms numbering
c  and returns as a result the father.
c
c  In this case nsons is the number of sons that each father has and ifath
c  is an array which tells the 
c
c  NO HOMOGENOUS DIRECTIONS.  In this case the mesh would grow to rapidly
c  if we followed the above strategy since every partition would index its
c  points to the ENTIRE mesh.  Furthermore, there would never be a need
c  to average to a node off processor since there is no spatial averaging.
c  Therefore, to properly account for this case we must recognize it and
c  inerrupt certain actions (i.e. assembly of the average across partitions).
c  This case is easily identified by noting that maxval(nsons) =1 (i.e. no
c  father has any sons).  Reiterating to be clear, in this case ifath does
c  not point to a global numbering but instead just points to itself.
c
      nfath=1  ! some architectures choke on a zero or undeclared
                 ! dimension variable.  This sets it to a safe, small value.
      if(((iLES .lt. 20) .and. (iLES.gt.0))
     &                   .or. (itwmod.gt.0)  ) then ! don't forget same
                                                    ! conditional in proces.f

c           read (igeomBAK) nfath  ! nfath already read in input.f,
                                     ! needed for alloc
         ione=1
c         call creadlist(igeomBAK,ione,nfath)
c         fname1='keyword sonfath?'
         if(nohomog.gt.0) then


!             fname1='number of father-nodes?'
!             call readheader(igeomBAK,fname1,nfath,ione,'integer', iotype)

            call phio_readheader(fhandle,
     &       c_char_'number of father-nodes' // char(0),
     &       c_loc(nfath), ione, dataInt, iotype)

c
c     fname1='keyword nsons?'
!             fname1='number of son-nodes for each father?'
!             call readheader(igeomBAK,fname1,nfath,ione,'integer', iotype)

            call phio_readheader(fhandle,
     &       c_char_'number of son-nodes for each father' // char(0),
     &       c_loc(nfath), ione, dataInt, iotype)

            allocate (point2nsons(nfath))

!             call readdatablock(igeomBAK,fname1,point2nsons,nfath,
!      &                      'integer',iotype)
            call phio_readdatablock(fhandle,
     &       c_char_'number of son-nodes for each father' // char(0),
     &       c_loc(point2nsons),nfath, dataInt, iotype)

c
!             fname1='keyword ifath?'
!             call readheader(igeomBAK,fname1,nshg,ione,'integer', iotype)

            call phio_readheader(fhandle,
     &       c_char_'keyword ifath' // char(0),
     &       c_loc(nshg), ione, dataInt, iotype);

            allocate (point2ifath(nshg))

!             call readdatablock(igeomBAK,fname1,point2ifath,nshg,
!      &                      'integer',iotype)
            call phio_readdatablock(fhandle,
     &       c_char_'keyword ifath' // char(0),
     &       c_loc(point2ifath), nshg, dataInt, iotype)
     
c     
            nsonmax=maxval(point2nsons)
c
         else  ! this is the case where there is no homogeneity
               ! therefore ever node is a father (too itself).  sonfath
               ! (a routine in NSpre) will set this up but this gives
               ! you an option to avoid that.
            nfath=nshg
            allocate (point2nsons(nfath))
            point2nsons=1
            allocate (point2ifath(nshg))
            do i=1,nshg
               point2ifath(i)=i
            enddo
            nsonmax=1
c
         endif
      else
         allocate (point2nsons(1))
         allocate (point2ifath(1))
      endif

      call phio_closefile_read(fhandle);

! !MR CHANGE
!       write(*,*) 'HELLO 13 from ', myrank
! !MR CHANGE END

c
c  renumber the master partition for SPEBC
c
c      if((myrank.eq.master).and.(irscale.ge.0)) then
c         call setSPEBC(numnp, nfath, nsonmax)
c         call renum(point2x,point2ifath,point2nsons)
c      endif
c
c.... Read restart files

c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      nfields = 11
c      numparts = 512
c      nppp = numparts/numpe
c      startpart = myrank * nppp +1
c      int endpart = startpart + nppp - 1
c      nppf = numparts/nfiles
cc      fnamer = "/users/nliu/PIG4/4-procs_case/restart-file"
cc      fnamer="./4-procs_case/restart-file"

      fnamer = c_char_"restart-dat."//c_null_char
      call phio_appendStep(fnamer, irstart)
      call phio_openfile_read(fnamer, nfiles, fhandle)

      ithree=3
c      call creadlist(irstin,ithree,nshg2,ndof2,lstep)

      itmp = int(log10(float(myrank+1)))+1

c      print *, "Solution is : ", fname1

      intfromfile=0
      call phio_readheader(fhandle,
     & c_char_'solution' // char(0), 
     & c_loc(intfromfile), ithree, dataInt, iotype)
c
c.... read the values of primitive variables into q
c

c      print *, "intfromfile(1) is ", intfromfile(1)
c      print *, "intfromfile(2) is ", intfromfile(2)
c      print *, "intfromfile(3) is ", intfromfile(3)

      allocate( qold(nshg,ndof) )
      if(intfromfile(1).ne.0) then
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)
         if(ndof2.ne.ndof) then

         endif
c
        if (nshg2 .ne. nshg)
     &        call error ('restar  ', 'nshg   ', nshg)
         allocate( qread(nshg,ndof2) )
         iqsiz=nshg*ndof2
         call phio_readdatablock(fhandle,
     &    c_char_'solution' // char(0),
     &    c_loc(qread),iqsiz, dataDbl,iotype)
         qold(:,1:ndof)=qread(:,1:ndof)
         deallocate(qread)
      else
         if (myrank.eq.master) then
            if (matflg(1,1).eq.0) then ! compressible
               warning='Solution is set to zero (with p and T to one)'
            else
               warning='Solution is set to zero'
            endif
            write(*,*) warning
         endif
         qold=zero
         if (matflg(1,1).eq.0) then ! compressible
            qold(:,1)=one ! avoid zero pressure
            qold(:,nflow)=one ! avoid zero temperature
         endif
      endif


! !MR CHANGE
!       write(*,*) 'HELLO 16-8 from ', myrank
! !MR CHANGE END

c      itmp=1
c      if (myrank .gt. 0) itmp = int(log10(float(myrank)))+1
      intfromfile=0
      call phio_readheader(fhandle,
     & c_char_'time derivative of solution' // char(0),
     & c_loc(intfromfile), ithree, dataInt, iotype)
      allocate( acold(nshg,ndof) )
      if(intfromfile(1).ne.0) then
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)

c      print *, "intfromfile(1) is ", intfromfile(1)
c      print *, "intfromfile(2) is ", intfromfile(2)
c      print *, "intfromfile(3) is ", intfromfile(3)

         if (nshg2 .ne. nshg)
     &        call error ('restar  ', 'nshg   ', nshg)
c
         allocate( acread(nshg,ndof2) )
         acread=zero
         iacsiz=nshg*ndof2
         call phio_readdatablock(fhandle,
     &    c_char_'time derivative of solution' // char(0),
     &    c_loc(acread), iacsiz, dataDbl,iotype)
         acold(:,1:ndof)=acread(:,1:ndof)
         deallocate(acread)
      else
         if (myrank.eq.master) then
            warning='Time derivative of solution is set to zero (SAFE)'
            write(*,*) warning
         endif
         acold=zero
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
c
cc
cc.... read the header and check it against the run data
cc


c      ithree=3
ccc      call creadlist(irstin,ithree,nshg2,ndof2,lstep)
c      fname1='solution?'
c      intfromfile=0
c      call readheader(irstin,fname1,intfromfile,
c     &     ithree,'integer', iotype)
cc
cc.... read the values of primitive variables into q
cc
c      allocate( qold(nshg,ndof) )
c      if(intfromfile(1).ne.0) then 
c         nshg2=intfromfile(1)
c         ndof2=intfromfile(2)
c         lstep=intfromfile(3)
c         if(ndof2.ne.ndof) then
c           warning='WARNING more data in restart than needed: keeping 1st '
c           write(*,*) warning , ndof
c         endif
cc
c         if (nshg2 .ne. nshg) 
c     &        call error ('restar  ', 'nshg   ', nshg)
c         allocate( qread(nshg,ndof2) )

c         iqsiz=nshg*ndof2
c         call readdatablock(irstin,fname1,qread,iqsiz,
c     &                         'double',iotype)
c         qold(:,1:ndof)=qread(:,1:ndof)
c         deallocate(qread)
c      else
c         if (myrank.eq.master) then
c            if (matflg(1,1).eq.0) then ! compressible
c               warning='Solution is set to zero (with p and T to one)'
c            else
c               warning='Solution is set to zero'
c            endif
c            write(*,*) warning
c         endif
c         qold=zero
c         if (matflg(1,1).eq.0) then ! compressible
c            qold(:,1)=one ! avoid zero pressure
c            qold(:,nflow)=one ! avoid zero temperature
c         endif
c      endif
cc 
c      fname1='time derivative of solution?'
c      intfromfile=0
c      call readheader(irstin,fname1,intfromfile,
c     &     ithree,'integer', iotype)
c      allocate( acold(nshg,ndof) )
c      if(intfromfile(1).ne.0) then 
c         nshg2=intfromfile(1)
c         ndof2=intfromfile(2)
c         lstep=intfromfile(3)
c         
c         if (nshg2 .ne. nshg) 
c     &        call error ('restar  ', 'nshg   ', nshg)
cc     
c         allocate( acread(nshg,ndof2) )
c         acread=zero
c
c         iacsiz=nshg*ndof2
c         call readdatablock(irstin,fname1,acread,iacsiz,
c     &                   'double',iotype)
c         acold(:,1:ndof)=acread(:,1:ndof)
c         deallocate(acread)
c      else
c         if (myrank.eq.master) then
c            warning='Time derivative of solution is set to zero (SAFE)'
c            write(*,*) warning
c         endif
c         acold=zero
c      endif

c      call creadlist(irstin,ithree,nshg2,ndisp,lstep)

      if (ideformwall.eq.1) then
!          fname1='displacement?'
!          call readheader(irstin,fname1,intfromfile,
!      &        ithree,'integer', iotype)

          intfromfile=0
          call phio_readheader(fhandle,
     &     c_char_'displacement' // char(0),
     &     c_loc(intfromfile), ithree, dataInt, iotype)

         nshg2=intfromfile(1)
         ndisp=intfromfile(2)
         lstep=intfromfile(3)
         if(ndisp.ne.nsd) then
            warning='WARNING ndisp not equal nsd'
            write(*,*) warning , ndisp
         endif
c
         if (nshg2 .ne. nshg) 
     &        call error ('restar  ', 'nshg   ', nshg)
c
c.... read the values of primitive variables into uold
c

         allocate( uold(nshg,nsd) )
         allocate( uread(nshg,nsd) )
         
         iusiz=nshg*nsd

!          call readdatablock(irstin,fname1,uread,iusiz,
!      &        'double',iotype)
         call phio_readdatablock(fhandle,
     &    c_char_'displacement' // char(0),
     &    c_loc(uread), iusiz, dataDbl, iotype)

         uold(:,1:nsd)=uread(:,1:nsd)
       else
         allocate( uold(nshg,nsd) )
         uold(:,1:nsd) = zero
       endif

c
c.... close c-binary files
c
!MR CHANGE
!      call closefile( irstin, "read" )

      call phio_closefile_read(fhandle)

!MR CHANGE
!      call closefile( igeomBAK,  "read" )
c
      deallocate(xread)
      if ( numpbc > 0 )  then
         deallocate(bcinpread)
         deallocate(ibctmpread)
      endif
      deallocate(iperread)
      if(numpe.gt.1)
     &     deallocate(ilworkread)
      deallocate(nbcread)

      return
c
 994  call error ('input   ','opening ', igeomBAK)
 995  call error ('input   ','opening ', igeomBAK)
 997  call error ('input   ','end file', igeomBAK)
 998  call error ('input   ','end file', igeomBAK)
c
      end

c
c No longer called but kept around in case....
c
      subroutine genpzero(iBC)

      use pointer_data
c
      include "common.h"
      integer iBC(nshg)
c
c....  check to see if any of the nodes have a dirichlet pressure
c
      pzero=1
      if (any(btest(iBC,2))) pzero=0  
c
      do iblk = 1, nelblb
         npro = lcblkb(1,iblk+1)-lcblkb(1,iblk)
         do i=1, npro
            iBCB1=miBCB(iblk)%p(i,1)
c     
c.... check to see if any of the nodes have a Neumann pressure 
c     but not periodic (note that 
c     
            if(btest(iBCB1,1)) pzero=0
         enddo
c     
c.... share results with other processors
c     
         pzl=pzero
         if (numpe .gt. 1)
     &        call MPI_ALLREDUCE (pzl, pzero, 1,
     &        MPI_DOUBLE_PRECISION,MPI_MIN, MPI_COMM_WORLD,ierr)
           
      enddo
c
c.... return
c
      return
c
      end

