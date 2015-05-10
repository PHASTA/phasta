c  readnblk.f (pronounce "Reed and Block Dot Eff") contains:
c
c    module readarrays ("Read Arrays") -- contains the arrays that
c     are read in from binary files but not immediately blocked 
c     through pointers.
c
c    subroutine readnblk ("Read and Block") -- allocates space for
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
      integer, allocatable :: point2ifath(:)
      integer, allocatable :: point2nsons(:)
      
      end module



      subroutine readnblk
c
      use readarrays
      include "common.h"
c
      real*8, allocatable :: xread(:,:), qread(:,:), acread(:,:)
      real*8, allocatable :: uread(:,:)
      real*8, allocatable :: BCinpread(:,:)
      integer, allocatable :: iperread(:), iBCtmpread(:)
      integer, allocatable :: ilworkread(:), nBCread(:)
      character*10 cname2
      character*8 mach2
      character*20 fmt1
      character*255 fname1,fnamer,fnamelr
      character*255 warning
      integer igeom, ibndc, irstin, ierr
      integer intfromfile(50) ! integers read from headers
c
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
      call openfile(  fname1,  'read?', igeom );
c
c.... try opening restart.latest.proc before trying restart.stepno.proc
c.... no restart.latest.proc !!
c      call openfile(  fnamelr,  'read?', irstin );
c      if ( irstin .eq. 0 ) then
          call openfile( fnamer, 'read?', irstin );
c      endif

! either one will work
c
c.... input the geometry parameters
c

      ieleven=11
      ione=1
      fname1='number of nodes?'
      call readheader(igeom,fname1,numnp,ione,'integer', iotype)
      fname1='number of modes?'
      call readheader(igeom,fname1,nshg,ione,'integer', iotype)
c      fname1='number of global modes?'
c      call readheader(igeom,fname1,nshgt,ione,'integer', iotype)
      fname1='number of interior elements?'
      call readheader(igeom,fname1,numel,ione,'integer', iotype)
      fname1='number of boundary elements?'
      call readheader(igeom,fname1,numelb,ione,'integer', iotype)
      fname1='maximum number of element nodes?'
      call readheader(igeom,fname1,nen,ione,'integer', iotype)
      fname1='number of interior tpblocks?'
      call readheader(igeom,fname1,nelblk,ione,'integer', iotype)
      fname1='number of boundary tpblocks?'
      call readheader(igeom,fname1,nelblb,ione,'integer', iotype)
      fname1='number of nodes with Dirichlet BCs?'
      call readheader(igeom,fname1,numpbc,ione,'integer', iotype)
      fname1='number of shape functions?'
      call readheader(igeom,fname1,ntopsh,ione,'integer', iotype)
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
      if(numpe > 1) then

         fname1='size of ilwork array?'
         call readheader(igeom,fname1,nlwork,ione,'integer', iotype)

         ione=1
         fname1='ilwork?'
         call readheader(igeom,fname1,nlwork,ione,'integer', iotype)

         allocate( point2ilwork(nlwork) )
         allocate( ilworkread(nlwork) )
         call readdatablock(igeom,fname1,ilworkread,
     &                      nlwork,'integer', iotype)
         point2ilwork = ilworkread
         call ctypes (point2ilwork)
      else
           nlwork=1
           allocate( point2ilwork(1))
           nshg0 = nshg
      endif
c     
c.... read the node coordinates
c
      itwo=2
      fname1='co-ordinates?'
      call readheader(igeom,fname1,intfromfile,itwo, 'double', iotype)
      numnp=intfromfile(1)
c      nsd=intfromfile(2)
      allocate( point2x(numnp,nsd) )
      allocate( xread(numnp,nsd) )
      ixsiz=numnp*nsd
      call readdatablock(igeom,fname1,xread,ixsiz, 'double',iotype)
      point2x = xread

c..............................for NGC
      if(istretchOutlet.eq.1)then
         
c...geometry6
        if(iNGCgeometryType .eq. 6) then
          xmaxn = 1.276
          xmaxo = 0.848
          xmin  = 0.42
c...geometry8
        elseif(iNGCgeometryType .eq. 8)then
          xmaxn=1.6*4.5*0.0254+0.85*1.5
          xmaxo=1.6*4.5*0.0254+0.85*1.0
          xmin =1.6*4.5*0.0254+0.85*0.5
        endif
c...
        alpha=(xmaxn-xmaxo)/(xmaxo-xmin)**2
        where (point2x(:,1) .ge. xmin)
c..... N=# of current elements from .42 to exit(~40)
c..... (x_mx-x_mn)/N=.025
c..... alpha=3    3*.025=.075
           point2x(:,1)=point2x(:,1)+
     &     alpha*(point2x(:,1)-xmin)**2
c..... ftn to stretch x at exit
        endwhere
      endif

        

c.....................................



c
c.... read in and block out the connectivity
c
      call genblk (IBKSIZ)
c
c.... read the boundary condition mapping array
c
      ione=1
      fname1='bc mapping array?'
      call readheader(igeom,fname1,nshg,
     &     ione,'integer', iotype)
      allocate( nBC(nshg) )

      allocate( nBCread(nshg) )
      call readdatablock(igeom,fname1,nBCread,nshg,'integer',iotype)
      nBC=nBCread
c
c.... read the temporary iBC array
c
      ione = 1
      fname1='bc codes array?'
      call readheader(igeom,fname1,numpbc,
     &     ione, 'integer', iotype)
      if ( numpbc > 0 ) then
         allocate( iBCtmp(numpbc) )
         allocate( iBCtmpread(numpbc) )
         call readdatablock(igeom,fname1,iBCtmpread,numpbc,
     &                      'integer',iotype)
         iBCtmp=iBCtmpread
      else  ! sometimes a partition has no BC's
         allocate( iBCtmp(1) )
         iBCtmp=0
      endif
c
c.... read boundary condition data
c
      ione=1
      fname1='boundary condition array?'
      call readheader(igeom,fname1,intfromfile,
     &     ione, 'double', iotype)
c here intfromfile(1) contains (ndof+7)*numpbc
      if ( numpbc > 0 ) then
         if(intfromfile(1).ne.(ndof+7)*numpbc) then
           warning='WARNING more data in BCinp than needed: keeping 1st'
           write(*,*) warning, ndof+7
         endif
         allocate( BCinp(numpbc,ndof+7) )
         nsecondrank=intfromfile(1)/numpbc
         allocate( BCinpread(numpbc,nsecondrank) )
         iBCinpsiz=intfromfile(1)
         call readdatablock(igeom,fname1,BCinpread,iBCinpsiz,
     &                      'double',iotype)
         BCinp(:,1:(ndof+7))=BCinpread(:,1:(ndof+7))
      else  ! sometimes a partition has no BC's
         allocate( BCinp(1,ndof+7) )
         BCinp=0
      endif
c
c.... read periodic boundary conditions
c
      fname1='periodic masters array?'
      call readheader(igeom,fname1,nshg,
     &     ione, 'integer', iotype)
      allocate( point2iper(nshg) )
      allocate( iperread(nshg) )
      call readdatablock(igeom,fname1,iperread,nshg,
     &                      'integer',iotype)
      point2iper=iperread
c
c.... generate the boundary element blocks
c
      call genbkb (ibksiz)

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

c           read (igeom) nfath  ! nfath already read in input.f,
                                     ! needed for alloc
         ione=1
c         call creadlist(igeom,ione,nfath)
c         fname1='keyword sonfath?'
         if(nohomog.gt.0) then
            fname1='number of father-nodes?'
            call readheader(igeom,fname1,nfath,ione,'integer', iotype)
c
c     fname1='keyword nsons?'
            fname1='number of son-nodes for each father?'
            call readheader(igeom,fname1,nfath,ione,'integer', iotype)
            allocate (point2nsons(nfath))
            call readdatablock(igeom,fname1,point2nsons,nfath,
     &                      'integer',iotype)
c
            fname1='keyword ifath?'
            call readheader(igeom,fname1,nshg,ione,'integer', iotype)
            allocate (point2ifath(nshg))
            call readdatablock(igeom,fname1,point2ifath,nshg,
     &                      'integer',iotype)
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
c
c  renumber the master partition for SPEBC
c
c      if((myrank.eq.master).and.(irscale.ge.0)) then
c         call setSPEBC(numnp, nfath, nsonmax)
c         call renum(point2x,point2ifath,point2nsons)
c      endif
c
c.... Read restart files
c
c.... read the header and check it against the run data
c

      ithree=3
c      call creadlist(irstin,ithree,nshg2,ndof2,lstep)
      fname1='solution?'
      intfromfile=0
      call readheader(irstin,fname1,intfromfile,
     &     ithree,'integer', iotype)
c
c.... read the values of primitive variables into q
c
      allocate( qold(nshg,ndof) )
      if(intfromfile(1).ne.0) then 
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)
         if(ndof2.ne.ndof) then
           warning='WARNING more data in restart than needed: 
     $              keeping 1st'
           write(*,*) warning , ndof
         endif
c
         if (nshg2 .ne. nshg) 
     &        call error ('restar  ', 'nshg   ', nshg)
         allocate( qread(nshg,ndof2) )

         iqsiz=nshg*ndof2
         call readdatablock(irstin,fname1,qread,iqsiz,
     &                         'double',iotype)
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
c 
      fname1='time derivative of solution?'
      intfromfile=0
      call readheader(irstin,fname1,intfromfile,
     &     ithree,'integer', iotype)
      allocate( acold(nshg,ndof) )
      if(intfromfile(1).ne.0) then 
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)
         
         if (nshg2 .ne. nshg) 
     &        call error ('restar  ', 'nshg   ', nshg)
c     
         allocate( acread(nshg,ndof2) )
         acread=zero

         iacsiz=nshg*ndof2
         call readdatablock(irstin,fname1,acread,iacsiz,
     &                   'double',iotype)
         acold(:,1:ndof)=acread(:,1:ndof)
         deallocate(acread)
      else
         if (myrank.eq.master) then
            warning='Time derivative of solution is set to zero (SAFE)'
            write(*,*) warning
         endif
         acold=zero
      endif

c      call creadlist(irstin,ithree,nshg2,ndisp,lstep)
      if (ideformwall.eq.1) then
         fname1='displacement?'
         call readheader(irstin,fname1,intfromfile,
     &        ithree,'integer', iotype)
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
         call readdatablock(irstin,fname1,uread,iusiz,
     &        'double',iotype)
         uold(:,1:nsd)=uread(:,1:nsd)
       else
         allocate( uold(nshg,nsd) )
         uold(:,1:nsd) = zero
       endif

c 
c
c.... close c-binary files
c
      call closefile( irstin, "read" )
      call closefile( igeom,  "read" )
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
 994  call error ('input   ','opening ', igeom)
 995  call error ('input   ','opening ', igeom)
 997  call error ('input   ','end file', igeom)
 998  call error ('input   ','end file', igeom)
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
