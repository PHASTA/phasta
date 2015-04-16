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
      real*8, allocatable :: dwal(:)
      real*8, allocatable :: errors(:,:)
      real*8, allocatable :: ybar(:,:)
      real*8, allocatable :: yphbar(:,:,:)
      real*8, allocatable :: vort(:,:)
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
      include "commonM2NFixBnd.h"
      include "mpif.h"
c
      real*8, allocatable :: xread(:,:), qread(:,:), qread1(:)
      real*8, allocatable :: uread(:,:), acread(:,:)
      real*8, allocatable :: BCinpread(:,:)
      real*8 globmax,globmin
      integer, allocatable :: iperread(:), iBCtmpread(:)
      integer, allocatable :: ilworkread(:), nBCread(:)
      character*10 cname2
      character*30 fmt1
      character*255 fname1,fnamer,fnamelr
      character*255 warning
      
      integer :: descriptor, color, nfiles, nfields
      integer ::  numparts, nppf
      character*255 fname2, temp2
      character*64 temp1
      integer :: igeom, ibndc, irstin, ierr
      integer :: ndof, ndoferrors, ndofybar
      integer :: itmp, itmp2
      integer :: irstart, irstartmap, iybar
      integer :: ierror, numphavg, ivort, idwal, idebug, iphavg

      integer intfromfile(50) ! integers read from headers
      logical exinput
c
c
c.... determine the step number to start with
c
!      open(unit=72,file='numstart.dat',status='old')
!      read(72,*) irstart
!      close(72)

      if(myrank == 0) then 
        fnamer='M2N_input.dat'
        fnamer = trim(fnamer) // char(0)
        inquire(file=fnamer,exist=exinput)
        if(exinput) then 
          open(unit=72,file=fnamer,status='old')
          read(72,*) irstart
          read(72,*) irstartmap
          read(72,*) iybar
          read(72,*) ierror
          read(72,*) numphavg
          read(72,*) ivort
          read(72,*) idwal
          read(72,*) idebug
          close(72)
        else
           write(*,*) 'ERROR: Input file ',
     &                 trim(fnamer),' does not exist!'
        endif
      endif
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(exinput,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
      if(.not. exinput) then ! M2NFixBnd_input.dat does not exist. Quit
!        call mpi_abort(mpi_comm_world,911,ierr)
!        call mpi_finalize(ierr)
        return
      else ! broadcast the information read by rank 0
        call mpi_bcast(irstart,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(iybar,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(ierror,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(numphavg,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(ivort,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(idwal,1,MPI_INTEGER,0,mpi_comm_world,ierr)
      endif

      if(myrank == 0 ) then
      ! Print some info
        write(*,*) 'The solution field is reduced by default'

        if(iybar .gt. 0) then
          write(*,*) 'The ybar field will also be reduced'
          iybar = 1 ! security
        else
          write(*,*) 'The ybar field will NOT be reduced'
        endif

        if(ierror .gt. 0) then
          write(*,*) 'The error field will also be reduced'
          ierror = 1 ! security
        else
          write(*,*) 'The error field will NOT be reduced'
        endif

        if(numphavg .gt. 0) then
          write(*,*) 'The phase average fields (',numphavg,
     &                ') will also be reduced'
        else
          write(*,*) 'The phase average fields will NOT be reduced'
        endif

        if(ivort .gt. 0) then
          write(*,*) 'The vorticity field will also be reduced'
          ivort = 1 ! security
        else
          write(*,*) 'The vorticity field will NOT be reduced'
        endif

        if(idwal .gt. 0) then
          write(*,*) 'The dwal field will also be reduced'
          idwal = 1 ! security
        else
          write(*,*) 'The dwal field will NOT be reduced'
        endif
        write(*,*) ''
      endif


      call mpi_barrier(mpi_comm_world, ierr)
      if(myrank == 0) then
        write(*,*) 'Reading the geombcRed-dat files'
      endif

      lstep=irstart ! in case restart files have no fields

      nfiles = nsynciofiles
      numparts = numpe !This is the common settings. Beware if you try to compute several parts per process
c
c.... input the geometry parameters
c

      color = int(myrank/(numparts/nfiles)) !Should call the color routine in SyncIO here
      itmp2 = int(log10(float(color+1)))+1
      write (temp2,"('(''geombcRed-dat.'',i',i1,')')") itmp2
      write (fnamer,temp2) (color+1)
      fnamer = trim(fnamer)//char(0)

      ieleven=11
      ione=1

      itmp = int(log10(float(myrank+1)))+1

      call queryphmpiio(fnamer, nfields, nppf);
      if (myrank == 0) then
        write(*,*) 'Number of fields in geombcRed-dat: ',nfields
        write(*,*) 'Number of parts per file geombcRed-dat: ',nppf
      endif
      call initphmpiio( nfields, nppf, nfiles, igeom,
     & 'read'//char(0))
      call openfile( fnamer, 'read'//char(0), igeom )

      write (temp1,"('(''number of nodes@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2//char(0),numnp,ione,
     & 'integer'//char(0), iotype)

      write (temp1,"('(''number of modes@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2//char(0),nshg,ione,
     & 'integer'//char(0), iotype)

      write (temp1,"('(''number of interior elements@'',i',i1,',A1)')")
     &       itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2//char(0),numel,ione,
     & 'integer'//char(0), iotype)

      write (temp1,"('(''number of boundary elements@'',i',i1,',A1)')")
     &       itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2//char(0),numelb,ione,
     & 'integer'//char(0),iotype)

      write (temp1,
     & "('(''maximum number of element nodes@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2//char(0),nen,ione,
     &'integer'//char(0),iotype)

      write (temp1,"('(''number of interior tpblocks@'',i',i1,',A1)')")
     &       itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2//char(0),nelblk,ione,
     & 'integer'//char(0) ,iotype)

      write (temp1,"('(''number of boundary tpblocks@'',i',i1,',A1)')")
     &       itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2//char(0),nelblb,ione,
     & 'integer'//char(0), iotype)

      write (temp1,
     & "('(''number of nodes with Dirichlet BCs@'',i',i1,',A1)')") itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2//char(0),numpbc,ione,
     & 'integer'//char(0),iotype)

      write (temp1,"('(''number of shape functions@'',i',i1,',A1)')")
     &       itmp
      write (fname2,temp1) (myrank+1),'?'
      call readheader(igeom,fname2//char(0),ntopsh,ione,
     & 'integer'//char(0),iotype)

      call closefile( igeom, "read"//char(0) )
      call finalizephmpiio( igeom )

c
c.... calculate the maximum number of boundary element nodes
c     
      nenb = 3 !was initialized to 0 but
      do i = 1, melCat !melCat is 0 here
         if (nen .eq. nenCat(i,nsd)) nenb = max(nenCat(i,nsd-1), nenb)
      enddo
c     
      if (myrank == master) then
         if (nenb .eq. 0) call error ('input   '//char(0),
     &                                'nen     '//char(0),nen)
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
c.... Read restart files
c
c.... read the header and check it against the run data
c
      call mpi_barrier(mpi_comm_world, ierr)
      if(myrank == 0) then
        write(*,*) 'Reading the RestartRedTmp-dat files'
      endif

!     Beware in what follows! nshg2 read from the header of the solution,
!     dwal, errors and ybar can be different from nshg read from the geombc files
!     This is due to the fact that nshg2 represents the largest vID referenced 
!     by the mapping and not the largest vID present in the mesh.
!     qold, dwal, errors, ybar must be allocated to nshg and initialized
!     to a large negative value. The unmapped vertices will be updated
!     through communication with ilwork.

      itmp=1
      if (irstart .gt. 0) itmp = int(log10(float(irstart+1)))+1

      write (fmt1,"('(''restartRedTmp-dat.'',i',i1,',1x)')") itmp
      write (fnamer,fmt1) irstart
      fnamer = trim(fnamer) // cname2(color+1)
      
      call queryphmpiio(fnamer//char(0), nfields, nppf);
      if (myrank == 0) then
        write(*,*) 'Number of fields in restartRedTmp-dat: ',nfields
        write(*,*) 'Number of parts per file restartRedTmp-dat: ',nppf
      endif

      call initphmpiio(nfields,nppf,nfiles,descriptor,
     & 'read'//char(0))
      call openfile( fnamer//char(0) , 
     & 'read'//char(0), descriptor )

c
c     Read the solution 
c
      ithree=3
      itmp = int(log10(float(myrank+1)))+1
      write (temp1,"('(''solution@'',i',i1,',A1)')") itmp
      write (fname1,temp1) (myrank+1),'?'
      fname1 = trim(fname1)

      intfromfile=0
      call readheader(descriptor,fname1//char(0) ,intfromfile,
     & ithree,'integer'//char(0), iotype)

c
c.... read the values of primitive variables into q
c
      if(intfromfile(1).ne.0) then 
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         ndof=ndof2 !This must be the same anyway
         allocate( qold(nshg,ndof) )
         qold(:,:) = -9.87654321e32
         lstep=intfromfile(3)
         allocate( qread(nshg2,ndof2) )

         if (nshg2 .ne. nshg) then 
           write(*,*) 'nshg from geombc and nshg2 from restart differ'
     &                //' on rank', myrank, ' :',nshg,nshg2,
     &               ' - Probably mixing phasta files'
         endif
         call mpi_barrier(mpi_comm_world, ierr)

         iqsiz=nshg2*ndof2

         call readdatablock(descriptor,fname1//char(0),qread,iqsiz,
     &                         'double'//char(0),iotype)
         qold(1:nshg2,1:ndof2)=qread(1:nshg2,1:ndof2)
         deallocate(qread)
      else
         if (myrank.eq.master) then
            if (matflg(1,1).eq.0) then ! compressible
               warning='Solution is set to zero (with p and T to one)'
            else
               warning='Solution is set to zero'
            endif
            write(*,*) warning// char(0)
         endif
         qold=zero
         if (matflg(1,1).eq.0) then ! compressible
            qold(:,1)=one ! avoid zero pressure
            qold(:,nflow)=one ! avoid zero temperature
         endif
      endif

c
c    Read the ybar
c
      ndofybar = 0
      if(iybar==1) then
        itmp = int(log10(float(myrank+1)))+1
        write (temp1,"('(''ybar@'',i',i1,',a1)')") itmp
        write (fname1,temp1) (myrank+1),'?'
        fname1 = trim(fname1)
        intfromfile=0
        call readheader(descriptor,fname1//char(0),intfromfile,
     &     ithree,'integer'//char(0),iotype)
        nshg2=intfromfile(1)
        ndofybar=intfromfile(2)
        !lstep=intfromfile(3)
        if(ndofybar .ne. 0) then 
          allocate( ybar(nshg,ndofybar) )
          ybar(:,:) = -9.87654321e32
          allocate( qread(nshg2,ndofybar) )

          iqsiz=nshg2*ndofybar
          call readdatablock(descriptor,fname1//char(0) ,qread,iqsiz,
     &                   'double'//char(0),iotype)
          ybar(1:nshg2,1:ndofybar)=qread(1:nshg2,1:ndofybar)
          deallocate(qread)
        else
          write(*,*) 'WARNING: ybar is missing in the restart files'
          iybar = 0
        endif
      endif

c
c    Read the errors
c
      ndoferrors=0
      if(ierror == 1) then
        write (temp1,"('(''errors@'',i',i1,',a1)')")
     &       itmp
        write (fname1,temp1) (myrank+1),'?'
        fname1 = trim(fname1)
        intfromfile=0
        call readheader(descriptor,fname1//char(0),intfromfile,
     &     ithree,'integer'//char(0),iotype)
        nshg2=intfromfile(1)
        ndoferrors=intfromfile(2)
        !lstep=intfromfile(3)
        if(ndoferrors .ne. 0) then 
          allocate( errors(nshg,ndoferrors) )
          errors(:,:) = -9.87654321e32
          allocate( qread(nshg2,ndoferrors) )
          iqsiz=nshg2*ndoferrors
          call readdatablock(descriptor,fname1//char(0),qread,iqsiz,
     &                   'double'//char(0),iotype)
          errors(1:nshg2,1:ndoferrors)=qread(1:nshg2,1:ndoferrors)
          deallocate(qread)
        else
          if(myrank==0) then
            write(*,*) 'WARNING: errors is missing in the restart files'
          endif
          ierror = 0
        endif
      endif


!
!     Read the phase_average fields
!
      ndofyphbar=0
      if(numphavg .gt. 0) then
        do iphavg = 1,numphavg
          itmp = int(log10(float(myrank+1)))+1
          itmp2 = int(log10(float(iphavg)))+1
          write (temp1,
     &            "('(''phase_average'',i',i1,',''@'',i',i1,',A1)')")
     &             itmp2, itmp
          write (fname1,temp1) iphavg,(myrank+1),'?'
          fname1 = trim(fname1)
          intfromfile=0
          call readheader(descriptor,fname1//char(0),intfromfile,
     &                  ithree,'integer'//char(0),iotype)

          nshg2=intfromfile(1)
          ndofyphbar=intfromfile(2)
          !lstep=intfromfile(3)

          if(ndofyphbar.ne.0) then
            ! Allocate some memory for the first ts only
            if(iphavg==1) then
              allocate( yphbar(nshg,ndofyphbar,numphavg) )
              yphbar(:,:,:) = -9.87654321e32
            endif

            allocate( qread(nshg2,ndofyphbar) )
            iqsiz = nshg2*ndofyphbar
            call readdatablock(descriptor,fname1//char(0),qread,iqsiz,
     &                     'double'//char(0),iotype)
            yphbar(1:nshg2,1:ndofyphbar,iphavg) = 
     &                         qread(1:nshg2,1:ndofyphbar)
            deallocate(qread)
          else
            if(myrank==0) then
              write(*,*)'WARNING: phase_average is missing '//
     &            'in the restart files'
            endif
            numphavg = 0
            if(iphavg > 1) then
              deallocate(yphbar)
            endif
            exit
          endif
        enddo
      endif

!
! follow the usual convention for loading the vorticity field
!
      ndofvort=0
      if(ivort == 1) then
        itmp = int(log10(float(myrank+1)))+1
        write (temp1,"('(''vorticity@'',i',i1,',a1)')") itmp
        write (fname1,temp1) (myrank+1),'?'
        fname1 = trim(fname1)
        intfromfile=0
        call readheader(descriptor,fname1//char(0),intfromfile,
     &     ithree,'integer'//char(0),iotype)

        nshg2=intfromfile(1)
        ndofvort=intfromfile(2)
        !lstep=intfromfile(3)

        if(ndofvort .ne. 0) then 
          allocate(vort(nshg,ndofvort))
          vort(:,:) = -9.87654321e32
          allocate(qread(nshg2,ndofvort))
          iqsiz = nshg2*ndofvort
          call readdatablock(descriptor,fname1//char(0),qread,iqsiz,
     &                   'double'//char(0),iotype)
          vort(1:nshg2,1:ndofvort)=qread(1:nshg2,1:ndofvort)
          deallocate(qread)
        else
          if(myrank==0) then
            write(*,*) 'WARNING: vorticity is missing '//
     &                 'in the restart files'
          endif
          ivort = 0
        endif
      endif

c
c     Read the dwal
c
      ndofdwal=0
      if(idwal==1) then
        write (temp1,"('(''dwal@'',i',i1,',a1)')") itmp
        write (fname1,temp1) (myrank+1),'?'
        fname1 = trim(fname1)
        intfromfile=0
        call readheader(descriptor,fname1//char(0),intfromfile,
     &     ithree,'integer'//char(0),iotype)

        nshg2=intfromfile(1)
        ndofdwal=intfromfile(2)
        if(ndofdwal .ne. 0) then 
          if(ndofdwal.ne.1) then
            warning='WARNING: ndofdwal not equal 1'
            write(*,*) warning, ndofdwal
          endif
          allocate( dwal(nshg) )
          dwal(:) = -9.87654321e32
          allocate( qread1(nshg2) )
          iqsiz=nshg2*1
          call readdatablock(descriptor,fname1//char(0),qread1,iqsiz,
     &                   'double'//char(0),iotype)
          dwal(1:nshg2)=qread1(1:nshg2)
          deallocate(qread1)
        else
          if(myrank==0) then
            write(*,*) 'WARNING: dwal is missing in the restart files'
          endif
          idwal = 0
        endif
      endif

!
!.... close c-binary files
!
      call closefile( descriptor, "read"//char(0) )
      call finalizephmpiio( descriptor )
c.... ----------------------> Communication tasks <--------------------
c

      if(numpe > 1) then

        call mpi_barrier(mpi_comm_world, ierr)
        if(myrank == 0) then
          write(*,*) 'Reading the geombc-dat files again for ilwork'
        endif

        color = int(myrank/(numparts/nfiles)) !Should call the color routine in SyncIO here
        itmp2 = int(log10(float(color+1)))+1
        write (temp2,"('(''geombcRed-dat.'',i',i1,')')") itmp2
        write (fnamer,temp2) (color+1)
        fnamer = trim(fnamer)//char(0)

        ieleven=11
        ione=1

        itmp = int(log10(float(myrank+1)))+1

        call queryphmpiio(fnamer, nfields, nppf);
        if (myrank == 0) then
          write(*,*) 'Number of fields in geombcRed-dat: ',nfields
          write(*,*) 'Number of parts per file geombcRed-dat: ',nppf
        endif
        call initphmpiio( nfields, nppf, nfiles, igeom,
     &       'read'//char(0))
        call openfile( fnamer, 'read'//char(0), igeom )

        write (temp1,"('(''size of ilwork array@'',i',i1,',A1)')") itmp
        write (fname2,temp1) (myrank+1),'?'
        call readheader(igeom,fname2//char(0),nlwork,ione,
     &   'integer' // char(0) ,iotype)

        write (temp1,"('(''ilwork@'',i',i1,',A1)')") itmp
        write (fname2,temp1) (myrank+1),'?'
        call readheader(igeom,fname2//char(0) ,nlwork,ione,
     &   'integer'//char(0) , iotype)

        allocate( point2ilwork(nlwork) )
        allocate( ilworkread(nlwork) )
        call readdatablock(igeom,fname2//char(0),ilworkread,
     &                      nlwork,'integer'//char(0) , iotype)

        point2ilwork = ilworkread
        deallocate(ilworkread)

        call closefile( igeom, "read"//char(0) )
        call finalizephmpiio( igeom )

        call ctypes (point2ilwork)
         
      else
        nlwork=1
        allocate( point2ilwork(1))
        nshg0 = nshg
      endif

c
c.... -------------------->   communications <-------------------------
c
      call mpi_barrier(mpi_comm_world, ierr)  ! make sure every rank is synced here
      if(myrank == 0) then
        write(*,*) 'Updating the vertices on the part boundaries'
      endif

      if (numpe > 1) then
         ! solution
          call commuMax (qold, point2ilwork, ndof, 'in '//char(0))
          call commuMax (qold, point2ilwork, ndof, 'out'//char(0))
          call mpi_barrier(mpi_comm_world, ierr)  ! make sure everybody is done with ilwork

          ! ybar
          if(iybar == 1) then
            call commuMax (ybar, point2ilwork, ndofybar, 'in '//char(0))
            call commuMax (ybar, point2ilwork, ndofybar, 'out'//char(0))
            call mpi_barrier(mpi_comm_world, ierr)  ! make sure everybody is done with ilwork
          endif

         ! errors
          if(ierror == 1) then
            call commuMax (errors, point2ilwork, ndoferrors, 
     &                                             'in '//char(0))
            call commuMax (errors, point2ilwork, ndoferrors, 
     &                                             'out'//char(0))
            call mpi_barrier(mpi_comm_world, ierr)  ! make sure everybody is done with ilwork
          endif

          ! phase_average
          if(numphavg .gt. 0) then
            do iphavg = 1,numphavg
              call commuMax (yphbar(:,:,iphavg), point2ilwork,
     &                                  ndofyphbar, 'in '//char(0))
              call commuMax (yphbar(:,:,iphavg), point2ilwork,
     &                                  ndofyphbar, 'out'//char(0))
              call mpi_barrier(mpi_comm_world, ierr)  ! make sure everybody is done with ilwork
            enddo
          endif

          ! vorticity
          if(ivort == 1) then
            call commuMax (vort, point2ilwork, ndofvort, 'in '//char(0))
            call commuMax (vort, point2ilwork, ndofvort, 'out'//char(0))
            call mpi_barrier(mpi_comm_world, ierr)  ! make sure everybody is done with ilwork
          endif

         ! dwal
          if(idwal == 1) then
            call commuMax (dwal, point2ilwork, 1, 'in '//char(0))
            call commuMax (dwal, point2ilwork, 1, 'out'//char(0))
            call mpi_barrier(mpi_comm_world, ierr)  ! make sure everybody is done with ilwork
          endif
      endif

c
c.... -------------------->  Print Min and Max of each field  <-------------------------
c
      call mpi_barrier(mpi_comm_world, ierr)
      if (myrank == 0) then
        write(*,*) 'Printing min and max of each field component'
        write(*,*) ''
      endif

      ! qold
      do idof=1,ndof
        call mpi_allreduce(maxval(qold(:,idof)),globmax,1,MPI_DOUBLE, 
     &                                  MPI_MAX, mpi_comm_world, ierr ) 
        call mpi_allreduce(minval(qold(:,idof)),globmin,1,MPI_DOUBLE, 
     &                                  MPI_MIN, mpi_comm_world, ierr ) 
        if (myrank == 0) then
          write(*,925) 'max/min qold(',idof,'): ',globmax,globmin
        endif
      enddo
      if (myrank == 0) then
        write(*,*) ''
      endif

      ! ybar
      if(iybar == 1) then
        do idof=1,ndofybar
          call mpi_allreduce(maxval(ybar(:,idof)),globmax,1,MPI_DOUBLE, 
     &                                  MPI_MAX, mpi_comm_world, ierr ) 
          call mpi_allreduce(minval(ybar(:,idof)),globmin,1,MPI_DOUBLE, 
     &                                  MPI_MIN, mpi_comm_world, ierr ) 
          if (myrank == 0) then
            write(*,925) 'max/min ybar(',idof,'): ',globmax,globmin
          endif
        enddo
        if (myrank == 0) then
          write(*,*) ''
        endif
      endif

      ! errors
      if(ierror == 1) then
        do idof=1,ndoferrors
          call mpi_allreduce(maxval(errors(:,idof)),globmax,1,
     &                      MPI_DOUBLE,MPI_MAX, mpi_comm_world, ierr ) 
          call mpi_allreduce(minval(errors(:,idof)),globmin,1,
     &                      MPI_DOUBLE,MPI_MIN, mpi_comm_world, ierr ) 
          if (myrank == 0) then
            write(*,925) 'max/min errors(',idof,'): ',globmax,globmin
          endif
        enddo
        if (myrank == 0) then
          write(*,*) ''
        endif
      endif

      ! phase_average
      if(numphavg .gt. 0) then
        do iphavg = 1,numphavg
          do idof=1,ndofyphbar
            call mpi_allreduce(maxval(yphbar(:,idof,iphavg)), globmax,1,
     &                     MPI_DOUBLE, MPI_MAX, mpi_comm_world, ierr ) 
            call mpi_allreduce(minval(yphbar(:,idof,iphavg)), globmin,1,
     &                     MPI_DOUBLE, MPI_MIN, mpi_comm_world, ierr ) 
            if (myrank == 0) then
              write(*,926) 'max/min yphbar(',idof,iphavg,'):',
     &                                globmax,globmin
            endif
          enddo
          if (myrank == 0) then
            write(*,*) ''
          endif
        enddo
      endif

      ! vorticity
      if(ivort == 1) then
        do idof=1,ndofvort
          call mpi_allreduce(maxval(vort(:,idof)),globmax,1,
     &                      MPI_DOUBLE,MPI_MAX, mpi_comm_world, ierr ) 
          call mpi_allreduce(minval(vort(:,idof)),globmin,1,
     &                      MPI_DOUBLE,MPI_MIN, mpi_comm_world, ierr ) 
          if (myrank == 0) then
            write(*,925) 'max/min vort(',idof,'): ',globmax,globmin
          endif
        enddo
        if (myrank == 0) then
          write(*,*) ''
        endif
      endif

      ! dwal
      if(idwal == 1) then
        call mpi_allreduce(maxval(dwal(:)), globmax, 1, MPI_DOUBLE, 
     &                                 MPI_MAX, mpi_comm_world, ierr ) 
        call mpi_allreduce(minval(dwal(:)), globmin, 1, MPI_DOUBLE, 
     &                                 MPI_MIN, mpi_comm_world, ierr ) 
        if (myrank == 0) then
          write(*,925) 'max/min dwal(',1,'):',globmax,globmin
          write(*,*) ''
        endif
      endif

925   format(A,i2,A,2(e24.17,1x))
926   format(A,i2,i2,A,2(e24.17,1x))

c
c.... e------------------->  Write data to disks  <-------------------------
c
      call mpi_barrier(mpi_comm_world, ierr)
      if(myrank == 0) then
        write(*,*) 'Writing the reduced restartRed-dat files'
      endif

!      call Write_M2NFixBnd(myrank, lstep, nshg, 
!     &          ndof, ndofybar, ndoferrors,
!     &          qold, ybar, errors, dwal)

      nsynciofieldswriterestart = 1+iybar+ierror+numphavg+ivort+idwal
      call Write_M2NFixBnd_SolOnly(myrank, lstep, nshg, 
     &          ndof, qold)
      deallocate(qold)

      ! ybar
      if(iybar == 1) then
        call Write_Field(myrank,'a','ybar',4,ybar,'d',
     &                    nshg,ndofybar,lstep) 
        deallocate(ybar)
      endif

      ! errors
      if(ierror == 1) then
        call Write_Field(myrank,'a','errors',6,errors,'d',
     &                    nshg,ndoferrors,lstep) 
        deallocate(errors)
      endif

      ! phase_average
      if(numphavg .gt. 0) then
        do iphavg=1,numphavg
          call write_phavg2(myrank,'a','phase_average',13,
     &             iphavg,numphavg,yphbar(:,:,iphavg),'d',
     &             nshg,ndofyphbar,lstep)
        enddo
        deallocate(yphbar)
      endif

      ! vorticity 
      if(ivort == 1) then
        call Write_Field(myrank,'a','vorticity',9,vort,'d',
     &                    nshg,ndofvort,lstep) 
        deallocate(vort)
      endif

      ! dwal
      if(idwal == 1) then
        call Write_Field(myrank,'a','dwal',4,dwal,'d',
     &                    nshg,1,lstep) 
        deallocate(dwal)
      endif


!
!     Deallocate some remaining memory
!
      if(numpe.gt.1) then
        deallocate(point2ilwork)
      endif

      return
c
 994  call error ('input   ','opening ', igeom)
 995  call error ('input   ','opening ', igeom)
 997  call error ('input   ','end file', igeom)
 998  call error ('input   ','end file', igeom)
c
      end
