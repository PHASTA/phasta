!  readnblk.f (pronounce "Reed and Block Dot Eff") contains:
!
!    module readarrays ("Red Arrays") -- contains the arrays that
!     are read in from binary files but not immediately blocked 
!     through pointers.
!
!    subroutine readnblk ("Reed and Block") -- allocates space for
!     and reads data to be contained in module readarrays.  Reads
!     all remaining data and blocks them with pointers.
!


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
!
      use readarrays
      include "commonM2N.h"
      include "mpif.h"

!
      real*8, allocatable :: xread(:,:), qread(:,:), acread(:,:)
      real*8, allocatable ::  qreadN(:), qoldN(:,:)
      real*8, allocatable :: pID(:), vID(:), dwalN(:), ybarN(:,:)
      real*8, allocatable :: errorsN(:,:), yphbarN(:,:,:), vortN(:,:)
      real*8, allocatable :: qTota2as(:,:), qTota2ar(:,:)
      real*8, allocatable :: qrecv(:), qsend(:)
      real*8 :: qreadN1
      integer, allocatable :: iperread(:), iBCtmpread(:)
      integer, allocatable :: ilworkread(:), nBCread(:)
      integer, allocatable :: indexpart(:) 
      integer, allocatable :: getinfo(:)
      character*10 cname2
      character*8 mach2
      character*30 fmt1
      character*255 fname1,fnamer, fnamelr
      character*255 warning
      integer :: igeomBAK, ibndc, irstin, ierr
      integer ::  icountN(numpe) ! integers read from headers
      integer :: intfromfile(50) ! integers read from headers
      integer :: maxvID(numpe),maxvIDout(numpe)
      integer :: nmap, ndim
      integer :: maxicountN, maxicountNglob, numvar, sumvar
      integer :: ncount, ifill, tag
      integer :: stat (MPI_STATUS_SIZE)
      integer :: ndofdwal,ndofybar,ndoferrors,ndofyphbar,ndofvort
      logical :: exmap, exinput
      integer ::  my_local_comm, my_local_size, my_local_rank
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccc New PhastaIO Definition Part ccccccccccccccccccccccccccccccccccccccccc

      integer :: descriptor, descriptorG, GPID, color, nfiles, nfields
      integer :: numparts, nppf, islesseqN
      integer :: ierr_io, numprocs, itmp, itmp2
      integer :: i, nonzero 
      integer :: ireducemethod, isbinary, irstart, irstartmap, iybar
      integer :: ierror, numphavg, ivort, idwal, idebug, iphavg, irankN
      character*255 fname2, temp2
      character*64 temp1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!
!  Set some parameters for the reduce operation
!
!     0 for communication based approach with mpi_send and mpi_recv (recommended), 1 for file based approach (not maintained any more)
      ireducemethod = 0

!     1 to write the restartMap* files in binary format, 0 in ascii (in conjunction with ireducemethod = 1)
      isbinary = 1
!
!     Read the input parameter for the the M2N reduction
!
!      open(unit=72,file='numstart.dat',status='old')
!      read(72,*) irstart
!      close(72)

      if(myrank == 0) then
        write(*,*) 'Reading M2N_input.dat'
        fnamer='M2N_input.dat'
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
      if(.not. exinput) then ! M2N_input.dat does not exist. Quit
        return
      else ! broadcast the information read by rank 0
        call mpi_bcast(irstart,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(irstartmap,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(iybar,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(ierror,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(numphavg,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(ivort,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(idwal,1,MPI_INTEGER,0,mpi_comm_world,ierr)
        call mpi_bcast(idebug,1,MPI_INTEGER,0,mpi_comm_world,ierr)
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
        write(*,*) 'Reading the geombc-dat files'
      endif

      lstep=irstart ! in case restart files have no fields

      nfiles = nsynciofiles
      numparts = numpe !This is the common settings. Beware if you try to compute several parts per process

      color = int(myrank/(numparts/nfiles)) !Should call the color routine in SyncIO here
      itmp2 = int(log10(float(color+1)))+1
      write (temp2,"('(''geombc-dat.'',i',i1,')')") itmp2
      write (fnamer,temp2) (color+1)
      fnamer = trim(fnamer)//char(0)

      itwo=2
      ione=1
      ieleven=11
      itmp = int(log10(float(myrank+1)))+1

      call queryphmpiio(fnamer, nfields, nppf);
      if (myrank == 0) then
        write(*,*) 'Number of fields in geombc-dat: ',nfields
        write(*,*) 'Number of parts per file geombc-dat: ',nppf
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

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!.... calculate the maximum number of boundary element nodes
!     
      nenb = 0
      do i = 1, melCat
         if (nen .eq. nenCat(i,nsd)) nenb = max(nenCat(i,nsd-1), nenb)
      enddo
!     
!
!.... setup some useful constants
!
      I3nsd  = nsd / 3          ! nsd=3 integer flag
      E3nsd  = float(I3nsd)     ! nsd=3 real    flag
!    
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
!     
      nsymdf = (ndof*(ndof + 1)) / 2 ! symm. d.o.f.'s
!
!.... ----------------------> Communication tasks <--------------------
!
      call closefile( igeom, "read"//char(0) )
      call finalizephmpiio( igeom )

!
!.... Read restart files for the solution, dwal, error and ybar
!
      call mpi_barrier(mpi_comm_world, ierr)
      if(myrank == 0) then
        write(*,*)'Reading the restart-dat files for the s-d-e-y fields'
      endif

      itmp=1
      if (irstart .gt. 0) itmp = int(log10(float(irstart+1)))+1

      write (fmt1,"('(''restart-dat.'',i',i1,',1x)')") itmp

      write (fnamer,fmt1) irstart
      fnamer = trim(fnamer) // cname2(color+1)

      call queryphmpiio(fnamer//char(0), nfields, nppf);
      if (myrank == 0) then
        write(*,*) 'Number of fields in restart-dat: ',nfields
        write(*,*) 'Number of parts per file restart-dat: ',nppf
      endif
      call initphmpiio(nfields,nppf,nfiles,descriptor,
     & 'read'//char(0))
      call openfile( fnamer//char(0) , 
     & 'read'//char(0), descriptor )

      ithree=3
!      call creadlist(irstin,ithree,nshg2,ndof2,lstep)

      itmp = int(log10(float(myrank+1)))+1
      write (temp1,"('(''solution@'',i',i1,',A1)')") itmp
      write (fname1,temp1) (myrank+1),'?'
      fname1 = trim(fname1)

!      print *, "Solution is : ", fname1

      intfromfile=0
      call readheader(descriptor,fname1//char(0) ,intfromfile,
     & ithree,'integer'//char(0), iotype)
!
!.... read the values of primitive variables into q
!

!      print *, "intfromfile(1) is ", intfromfile(1)
!      print *, "intfromfile(2) is ", intfromfile(2)
!      print *, "intfromfile(3) is ", intfromfile(3)

      if(intfromfile(1).ne.0) then
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)
         if(nshg2.ne.nshg) then
           write(*,*) 'ERROR: mistmatch between nshg '//
     &      'from the geombc and restart files',nshg,nshg2
         endif

         ndof=ndof2
         allocate( qold(nshg,ndof2) )
         allocate( qread(nshg,ndof2) )
         iqsiz=nshg*ndof2
         call readdatablock(descriptor,fname1//char(0),qread,iqsiz,
     &                         'double'//char(0),iotype)
         qold(:,1:ndof)=qread(:,1:ndof)
         deallocate(qread)
      else
         if (myrank.eq.master) then
            if (matflg(1,1).eq.0) then ! compressible
               warning='WARNING: Solution is set to zero '//
     &                   '(with p and T to one)'
            else
               warning='WARNING: Solution is set to zero'
            endif
            write(*,*) warning
         endif
         qold=zero
         if (matflg(1,1).eq.0) then ! compressible
            qold(:,1)=one ! avoid zero pressure
            qold(:,nflow)=one ! avoid zero temperature
         endif
      endif

!
!  Above this line is the usual loading of fields in readnblk.  Quite a bit of stuff has been removed.
!
 

!
! follow the usual convention for loading the ybar
!
      ndofybar=0
      if(iybar == 1) then
        itmp = int(log10(float(myrank+1)))+1
        write (temp1,"('(''ybar@'',i',i1,',a1)')") itmp
        write (fname1,temp1) (myrank+1),'?'
        fname1 = trim(fname1)
        intfromfile=0
        call readheader(descriptor,fname1//char(0),intfromfile,
     &     ithree,'integer'//char(0),iotype)

        !nshg2=intfromfile(1)
        ndofybar=intfromfile(2)
        !lstep=intfromfile(3)

        if(ndofybar .ne. 0) then 
          allocate( ybar(nshg,ndofybar) )
          iqsiz=nshg*ndofybar
          call readdatablock(descriptor,fname1//char(0) ,ybar,iqsiz,
     &                   'double'//char(0),iotype)
        else
          if(myrank==0) then
            write(*,*) 'WARNING: ybar is missing in the restart files'
          endif
          iybar = 0
        endif
      endif

!
! follow the usual convention for loading the error
!
      ndoferrors=0
      if(ierror == 1) then
        itmp = int(log10(float(myrank+1)))+1
        write (temp1,"('(''errors@'',i',i1,',a1)')") itmp
        write (fname1,temp1) (myrank+1),'?'
        fname1 = trim(fname1)
        intfromfile=0
        call readheader(descriptor,fname1//char(0),intfromfile,
     &     ithree,'integer'//char(0),iotype)

        !nshg2=intfromfile(1)
        ndoferrors=intfromfile(2)
        !lstep=intfromfile(3)

        if(ndoferrors .ne. 0) then 
          allocate( errors(nshg,ndoferrors) )
          iqsiz=nshg*ndoferrors
          call readdatablock(descriptor,fname1//char(0),errors,iqsiz,
     &                   'double'//char(0),iotype)
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

          !nshg=intfromfile(1) !Do not use nshg and ndof from common.h here!
          ndofyphbar=intfromfile(2)
          !lstep=intfromfile(3)

          if(ndofyphbar.ne.0) then
            ! Allocate some memory for the first ts only
            if(iphavg==1) then
              allocate( yphbar(nshg,ndofyphbar,numphavg) )
            endif

            allocate( qread(nshg,ndofyphbar)  ) ; qread = 0.d0
            iqsiz = nshg*ndofyphbar
            call readdatablock(descriptor,fname1//char(0),qread,iqsiz,
     &                     'double'//char(0),iotype)
            yphbar(:,:,iphavg) = qread(:,:)
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

        !nshg2=intfromfile(1)
        ndofvort=intfromfile(2)
        !lstep=intfromfile(3)

        if(ndofvort .ne. 0) then 
          allocate( vort(nshg,ndofvort) )
          iqsiz=nshg*ndofvort
          call readdatablock(descriptor,fname1//char(0) ,vort ,iqsiz,
     &                   'double'//char(0),iotype)
        else
          if(myrank==0) then
            write(*,*) 'WARNING: vorticity is missing '//
     &                 'in the restart files'
          endif
          ivort = 0
        endif
      endif


!
! follow the usual convention for loading the d2wal
!
      ndofdwal=0
      if(idwal == 1 ) then
        write (temp1,"('(''dwal@'',i',i1,',a1)')")
     &       itmp
        write (fname1,temp1) (myrank+1),'?'
        fname1 = trim(fname1)
        intfromfile=0
        call readheader(descriptor,fname1//char(0),intfromfile,
     &     ithree,'integer'//char(0),iotype)

         !nshg2=intfromfile(1)
         ndofdwal=intfromfile(2)
         !lstep=intfromfile(3)

        if(ndofdwal .ne. 0) then 
          if(ndofdwal.ne.1) then
            warning='WARNING: ndofdwal not equal 1'
            write(*,*) warning, ndofdwal
          endif
          allocate( dwal(nshg) )
          iqsiz=nshg
          call readdatablock(descriptor,fname1//char(0), dwal, iqsiz,
     &                   'double'//char(0),iotype)
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

!
!.... Read restart files for the vtx mapping between M and N parts
!
      call mpi_barrier(mpi_comm_world, ierr)
      if(myrank == 0) then
        write(*,*) 'Reading the restart-dat files for vtx mapping'
      endif

      itmp=1
      if (irstartmap .gt. 0) itmp = int(log10(float(irstartmap+1)))+1

      write (fmt1,"('(''restart-dat.'',i',i1,',1x)')") itmp

      write (fnamer,fmt1) irstartmap
      fnamer = trim(fnamer) // cname2(color+1)

      call queryphmpiio(fnamer//char(0), nfields, nppf);
      if (myrank == 0) then
        write(*,*) 'Number of fields in restart-dat: ',nfields
        write(*,*) 'Number of parts per file restart-dat: ',nppf
      endif
      call initphmpiio(nfields,nppf,nfiles,descriptor,
     & 'read'//char(0))
      call openfile( fnamer//char(0) , 
     & 'read'//char(0), descriptor )

!
! follow the usual convention for loading the partid
!
      itmp = int(log10(float(myrank+1)))+1
       write (temp1,"('(''mapping_partid@'',i',i1,',A1)')")
     &          itmp
       write (fname1,temp1) (myrank+1),'?'
       fname1 = trim(fname1)
       intfromfile=0
       call readheader(descriptor,fname1//char(0),intfromfile,
     &        ithree,'integer'//char(0),iotype)

       nshg2=intfromfile(1)
       nmap=intfromfile(2)
       !lstep=intfromfile(3) !lstep comes from another restart file
       if(nmap.ne.1) then
          warning='WARNING mapping_partid not equal 1'
          write(*,*) warning , nmap, intfromfile(1), intfromfile(2)
       endif

       allocate( pID(nshg) )
       iqsiz=nshg
       call readdatablock(descriptor,fname1//char(0) ,pID,iqsiz,
     &                   'double'//char(0),iotype)

!
! follow the usual convention for loading the vtxid
!
        write (temp1,"('(''mapping_vtxid@'',i',i1,',a1)')")
     &           itmp
        write (fname1,temp1) (myrank+1),'?'
        fname1 = trim(fname1)
        intfromfile=0
        call readheader(descriptor,fname1//char(0),intfromfile,
     &     ithree,'integer'//char(0),iotype)

       nshg2=intfromfile(1)
       nmap=intfromfile(2)
       !lstep=intfromfile(3) !lstep comes from another restart file
       if(nmap.ne.1) then
          warning='warning mapping_partid not equal 1'
          write(*,*) warning , nmap
       endif

       allocate( vid(nshg) )
       iqsiz=nshg
       call readdatablock(descriptor,fname1//char(0) ,vid,iqsiz,
     &                   'double'//char(0),iotype)

!
!.... close c-binary files
!
      call closefile( descriptor, "read"//char(0) )
      call finalizephmpiio( descriptor )

!
!  Debugging information
!
!      minvid=minval(vid)
!      write(*,*) minvid,myrank
!      do i=1,nshg
!      write(800+myrank,*) pID(i),vID(i),qold(i,1), qold(i,2)
!      ivID=vID(i)+0.1
!      ipID=pID(i)+0.1
!      write(900+myrank,*) ipID,ivID,qold(i,1), qold(i,2)
!      enddo
     
!
! here starts the mapping code
! 

!  first, on each of the M parts, count the number of vertices on each of the 
!  N parts and put the result in icountN for each of the M parts
!
      call mpi_barrier(mpi_comm_world, ierr) 
      if(myrank == 0) then
        write(*,*) 'Computing the size of the mapped partition'
      endif

      maxPID=0
      icountN=0
      do i=1,nshg
          ifileN=pID(i)+1.1
          icountN(ifileN)=icountN(ifileN)+1
          maxPID=max(maxPID,ifileN)
      enddo
!
!     ALLREDUCE so that every processor knows N the size of the mapped partition
!
      call MPI_ALLREDUCE (maxPID, maxPIDout, 1,
     &          MPI_INT, MPI_MAX, mpi_comm_world, ierr)

      irankN=maxPIDout

      call mpi_barrier(mpi_comm_world, ierr) 
      if(myrank.eq.0) then
        write (*,*) 'M2N reduction with M =',
     &                        numpe, 'and N =',irankN
      endif
!
!     DEBUG echo to the screen the number of number of vtx that will be mapped to 
!     each of the N for each of the M parts  
!
      if (idebug == 1) then
        call mpi_barrier(mpi_comm_world, ierr) 
        do i=1,numpe
          if(myrank == i) then
            write(*,926) myrank, (icountN(j),j=1,irankN)
          endif
          call mpi_barrier(mpi_comm_world, ierr) 
        enddo
      endif

!
!     Compute the max vertex ID for every N part
!
      call mpi_barrier(mpi_comm_world, ierr) 
      if(myrank == 0) then
        write(*,*) 'Computing the max vertex ID for every N part'
      endif

      maxvID(:)=0 !only maxvID(1:irankN) should have non-zero entries at the end of the program
      do i=1,nshg
        ifileN=pID(i)+1.1    !shift pID to start from 1
        intvID=(vID(i)+1.1)  !shift vID to start from 1
        maxvID(ifileN)=max(maxvID(ifileN),intvID)
      enddo

!
!     Sanity check for maxvID
!
      do i=irankN+1,numpe
        if (maxvID(i) .ne. 0) then
           write(*,*) 'ERROR: maxvID(',i,') not equal to 0 on rank',
     &                 myrank
        endif 
      enddo
!
!     ALLREDUCE so that every processor knows the number of vtx on each of the N 
!     parts in the reduced partition. No need to reduce maxvID(irankN+1:numpe), as it should be 0
!
      call MPI_ALLREDUCE (maxvID, maxvIDout, irankN, 
     &          MPI_INT, MPI_MAX, mpi_comm_world, ierr)

!
!    echo debug information to the screen
!
      if (idebug == 1) then
        if(myrank.eq.0) write(*,*) 'Writing maxvID on each of N ranks'
        if(myrank.eq.0) write(*,926) (maxvIDout(j),j=1,irankN)
      endif

!     Beware! maxvIDout does not represent the largest vID present in the parts
!     but the largest vID referenced in the parts through the mapping.
!     On part boundaries, it is possible to have unreferenced vID that will need
!     to be fixed through P2P communication.

!
!     Allocate the memory for qoldN on N parts handled by ranks 1..N.
!     Memory available per comupte node could be improved by chosing others ranks
!
      call mpi_barrier(mpi_comm_world, ierr) 
      if(myrank == 0) then
        write(*,*) 'Allocating memory for the reduced s-d-e-y vectors'
      endif

      if(myrank.lt.irankN) then
        nshgN=maxvIDout(myrank+1) ! number of vtx on my part computed above

        allocate(qoldN(nshgN,ndof))! size of my final mapped field array
        qoldN=-9.87654321e32       ! set to absurd, easily recognized value to find holes

        if(iybar == 1) then
          allocate(ybarN(nshgN,ndofybar))   ! size of my final mapped field array
          ybarN=-9.87654321e32       ! set to absurd, easily recognized value to find holes
        endif

        if(ierror == 1) then
          allocate(errorsN(nshgN,ndoferrors))   ! size of my final mapped field array
          errorsN=-9.87654321e32       ! set to absurd, easily recognized value to find holes
        endif

        if(numphavg .gt. 0) then
          allocate(yphbarN(nshgN,ndofyphbar,numphavg))   ! size of my final mapped field array
          yphbarN=-9.87654321e32       ! set to absurd, easily recognized value to find holes
        endif

        if(ivort == 1) then
          allocate( vortN(nshgN,ndofvort))
          vortN=-9.87654321e32
        endif

        if(idwal == 1) then
          allocate(dwalN(nshgN))   ! size of my final mapped field array
          dwalN=-9.87654321e32       ! set to absurd, easily recognized value to find holes
        endif
 

      else ! Not used but sane, as passed to write_M2N
        allocate(qoldN(1,1))
        if(iybar == 1) then
          allocate(ybarN(1,1))
        endif
        if(ierror == 1) then
          allocate(errorsN(1,1))
        endif
        if(numphavg .gt. 0) then
          allocate(yphbarN(1,1,1))   
        endif
        if(ivort == 1) then
          allocate(vortN(1,1))
        endif
        if(idwal == 1) then
          allocate(dwalN(1))
        endif
      endif

!
!     Most important routine: build the qoldN array for the reduction of the solution from M to N
!
      call mpi_barrier(mpi_comm_world, ierr) 
      if(myrank == 0) then
        write(*,*) 'Building the reduced field vectors'
      endif

      if (ireducemethod  == 0) then
!
!       Communication approach through mpi_isend and mpi_irecv to limit memory consumption
!       and avoid disk access
!
        numvar = 1+ndof+ndofybar+ndoferrors
     &           +(numphavg*ndofyphbar)+ndofvort+ndofdwal
!        write(*,*) numvar,ndof,ndofybar,ndoferrors,numphavg,ndofyphbar,
!     &             ndofvort
        call mpi_barrier(mpi_comm_world, ierr) 
        !1 for vtxid, ndof for solution, ndofdwal for dwal, ndoferrors for error, ndofybar for ybar

        do irN = 0,irankN-1  ! every irN rank 1..N receives information from


          ! Get which of the numpe-1 ranks have some information to communicate to rank irN
          if (myrank == irN) then
            allocate(getinfo(numpe))
            getinfo(:) = 0
          endif

          ifill = icountN(irN+1)
          call mpi_gather(ifill, 1, MPI_INT, getinfo, 1, MPI_INT, 
     &                   irN, mpi_comm_world, ierr)

!          if (myrank == irN) then
!            write(*,*) getinfo(:)
!          endif
          
          do irM = 0,numpe-1 ! all the irM ranks 1..M

! DEBUG
!            if (myrank == 0) then
!              write(*,*) 'NEW data exchange between irN:',
!     &                  irN,'and irM:',irM
!            endif
!            call mpi_barrier(mpi_comm_world, ierr) 
!DEBUG
      
            if (myrank == irN .and. myrank == irM) then ! no need to send the data. Already there so get them
!             write(*,*) '  myrank=irN=irM=',irN,' so treat data locally'
              do i=1,nshg
                ifileN = pID(i)+1.1 !shift to start from 1
                if(ifileN == irN+1) then ! these data need to be packed because rank irN is expecting them
                  intvID = vID(i)+1.1 ! shift to start from 1
                  ! Solution qoldN
                  do idof=1,ndof
                    qoldN(intvID,idof) = qold(i,idof)
                  enddo
                  ! ybarN
                  if(iybar==1) then
                    do idof=1,ndofybar
                      ybarN(intvID,idof) = ybar(i,idof)
                    enddo
                  endif
                  ! errorsN
                  if(ierror==1) then
                    do idof=1,ndoferrors
                      errorsN(intvID,idof) = errors(i,idof)
                    enddo
                  endif
                  ! phase average
                  if(numphavg .gt. 0) then
                    do iphavg = 1,numphavg
                      do idof=1,ndofyphbar
                        yphbarN(intvID,idof,iphavg) = 
     &                        yphbar(i,idof,iphavg)
                      enddo
                    enddo
                  endif
                  ! vortN
                  if(ivort==1) then
                    do idof=1,ndofvort
                      vortN(intvID,idof) = vort(i,idof)
                    enddo
                  endif
                 ! Distance to the wall dwalN
                  if(idwal==1) then
                    dwalN(intvID) = dwal(i)
                  endif
               endif 
              enddo

            else !myrank ne to irN and irM
              tag=irN
              if(myrank == irN) then 
                ! rank irN receives data from rank j and update qoldN. 
!                call mpi_recv(ifill, 1, MPI_INTEGER, irM, tag,  
!     &                                    mpi_comm_world, stat, ierr)
!                write(*,*) '  ', irN,'receives',ifill,getinfo(irM+1),
!     &                     'data from',irM
                ifill = getinfo(irM+1)
                if(ifill .gt. 0) then
                  ! there are really data to receive
!                  write(*,*) '  ', irN,'receives',ifill,'data from',irM
                  allocate(qrecv(ifill*numvar))
                  qrecv = -9.87654321e32 
                  call mpi_recv(qrecv, ifill*numvar, MPI_DOUBLE, irM,
     &                                tag, mpi_comm_world, stat, ierr)
                  do i=1,ifill
                    intvID = qrecv( (i-1)*numvar+1 )+1.1 ! shift to start from 1
                    ! Solution qold
                    do idof=1,ndof
                      sumvar = 1+idof
                      qoldN(intvID,idof) = qrecv( (i-1)*numvar+sumvar )
                    enddo
                    ! ybarN
                    if(iybar==1) then
                      do idof=1,ndofybar
                        sumvar = 1+ndof+idof
                        ybarN(intvID,idof) = qrecv( (i-1)*numvar+sumvar)
                      enddo
                    endif
                    ! errorN
                    if(ierror==1) then
                      do idof=1,ndoferrors
                        sumvar = 1+ndof+ndofybar+idof
                        errorsN(intvID,idof) = 
     &                            qrecv( (i-1)*numvar+sumvar)
                      enddo
                    endif
                    ! yphbarN
                    if(numphavg .gt. 0) then
                      do iphavg = 1,numphavg
                        do idof=1,ndofyphbar
                          sumvar = 1+ndof+ndofybar+ndoferrors+
     &                                    (iphavg-1)*ndofyphbar+idof
                          yphbarN(intvID,idof,iphavg) =  
     &                                     qrecv( (i-1)*numvar+sumvar)
                        enddo
                      enddo
                    endif
                    ! vortN
                    if(ivort==1) then
                      do idof=1,ndofvort
                        sumvar = 1+ndof+ndofybar+ndoferrors+
     &                                    numphavg*ndofyphbar+idof
                        vortN(intvID,idof) = 
     &                            qrecv( (i-1)*numvar+sumvar)
                      enddo
                    endif
                    ! Distance to the wall dwalN
                    if(idwal==1) then
                      sumvar = 1+ndof+ndofybar+ndoferrors+
     &                         numphavg*ndofyphbar+ndofvort+1
                      dwalN(intvID) = qrecv ( (i-1)*numvar+sumvar)
                    endif
                  enddo
                  deallocate(qrecv)
                endif

              elseif(myrank == irM) then ! rank irM sends data to rank irN
                ! First build the data to send to rank i for qoldN
                ifill = 0
                if (icountN(irN+1) .gt. 0) then !there are data to transfer
                  allocate(qsend(icountN(irN+1)*numvar))
                  qsend(:) = -9.87654321e32
                  do i=1,nshg
                    ifileN = pID(i)+1.1 !shift to start from 1
                    if(ifileN == irN+1) then 
                      ! these data need to be packed because rank irN is expecting them
                      ifill = ifill+1
                      ! vtx ID
                      qsend( (ifill-1)*numvar+1) = vID(i) ! no shift to start from 1 yet
                      ! Solution qold
                      do idof=1,ndof
                        sumvar = 1+idof
                        qsend( (ifill-1)*numvar+sumvar) = qold(i,idof)
                      enddo
                      ! ybar
                      if(iybar==1) then
                        do idof=1,ndofybar
                          sumvar = 1+ndof+idof
                          qsend( (ifill-1)*numvar+sumvar) = ybar(i,idof)
                        enddo
                      endif
                     ! errors
                      if(ierror==1) then
                        do idof=1,ndoferrors
                          sumvar = 1+ndof+ndofybar+idof
                          qsend( (ifill-1)*numvar+sumvar) = 
     &                                                 errors(i,idof)
                        enddo
                      endif
                      ! yphbarN
                      if(numphavg .gt. 0) then
                        do iphavg = 1,numphavg
                          do idof=1,ndofyphbar
                            sumvar = 1+ndof+ndofybar+ndoferrors+
     &                                    (iphavg-1)*ndofyphbar+idof
                            qsend( (ifill-1)*numvar+sumvar) = 
     &                                            yphbar(i,idof,iphavg)
                          enddo
                        enddo
                      endif
                      ! vortN
                      if(ivort==1) then
                        do idof=1,ndofvort
                          sumvar = 1+ndof+ndofybar+ndoferrors+
     &                                    numphavg*ndofyphbar+idof
                          qsend( (ifill-1)*numvar+sumvar) = 
     &                                          vort(i,idof)
                        enddo
                      endif
                      ! Distance to the wall dwal
                      if(idwal==1) then
                        sumvar = 1+ndof+ndofybar+ndoferrors+
     &                           numphavg*ndofyphbar+ndofvort+1
                        qsend( (ifill-1)*numvar+sumvar) = dwal(i)
                      endif
                    endif
                  enddo
                endif  
                if(ifill .ne. icountN(irN+1)) then
                  write(*,*) 'ERROR with data set:', irM, irN, 
     &                                   ifill, icountN(irN+1)
                endif
!                write(*,*) '  ', irM,'sends', ifill,'data to',irN
!                call mpi_send(ifill, 1 ,MPI_INTEGER, irN, tag,   
!     &                                          mpi_comm_world, ierr) ! Send the size of the data
                if(ifill .gt. 0) then !if the size of the data is >0, send it for good
!                  write(*,*) '  ', irM,'sends', ifill,'data to',irN
                  call mpi_send(qsend(1), ifill*numvar, MPI_DOUBLE, irN,
     &                                      tag, mpi_comm_world, ierr)
                  deallocate(qsend)
                endif

              endif !if(myrank==irN)
            endif !if (myrank == irN = irM)

          enddo !loop over irM

          call mpi_barrier(mpi_comm_world,ierr) 
          if (myrank == irN) then
            nonzero = 0
            do i = 1,numpe
              if (getinfo(i) .gt. 0) then
                nonzero = nonzero +1
              endif
            enddo
            write(*,*)"Part ",irN,"out of N = ",irankN," updated from",
     &                nonzero, "parts"
            deallocate(getinfo)
          endif          
        enddo !loop over irN

!
!       Do it without files through mpialltoall
!
!      maxicountN=maxval(icountN(:))
!      call MPI_ALLREDUCE (maxicountN, maxicountNglob, 1,
!     &          MPI_INT, MPI_MAX, mpi_comm_world, ierr)
!
!      numvar = 1+ndof ! vID + ndof from solution
!      write(*,*) 'maxicountN:',myrank, maxicountN, maxicountNglob, nshg
!      call mpi_barrier(mpi_comm_world, ierr)
!      allocate(qTota2as(numvar*maxicountNglob,numpe))
!      qTota2as(:,:) = -9.87654321e32

!      allocate(indexpart(numpe)); indexpart(:)=0
!      do i=1,nshg
!        ifileN = pID(i)+1.1                              !shift pID to start from 1
!        indexpart(ifileN) = indexpart(ifileN)+1
!        qTota2as( (indexpart(ifileN)-1)*numvar+1 , ifileN ) = vID(i)   ! vID not yet shifted to start from 1 here
!        do j=1,ndof
!          
!         if((indexpart(ifileN)-1)*numvar+1+j>numvar*maxicountNglob) then
!            write(*,*) 'ERROR DIMENSION1:',myrank,
!     &       ifileN, indexpart(ifileN),
!     &       (indexpart(ifileN)-1)*numvar+1+j, numvar*maxicountNglob   
!          endif
!          if (ifileN > numpe ) then
!              write(*,*) 'ERROR DIMENSION2:',myrank, ifileN, numpe
!          endif
!          qTota2as( (indexpart(ifileN)-1)*numvar+1+j, ifileN )=qold(i,j)
!        enddo
!      enddo
      
!      allocate(qTota2ar(numvar*maxicountNglob,numpe))
!      qTota2ar(:,:) = -9.87654321e32

!      ncount = numvar*maxicountNglob
!      call mpi_alltoall(qTota2as(1,1), ncount, MPI_DOUBLE, 
!     &                  qTota2ar(1,1), ncount, MPI_DOUBLE, 
!     &                  mpi_comm_world, ierr)
!      if (ierr .ne. 0) then 
!        write(*,*) 'Error with mpi_alltoall:', myrank, ierr
!      endif
!      deallocate(qTota2as)


      elseif (ireducemethod == 1) then ! file based approach
!
!       Now for each process on M parts push the fields into an i,j pair named
!       ascii file  restartMap.i.j   1<=i<=N  1<=j<=M.  This is embarassingly 
!       parallel with no data collision but writes a lot of files
!
        do l=1,irankN
          if(icountN(l).gt.0) then  !  I have data to write for rank l
!
!           generate the filename since this is a file we have to write
!           note this avoids opening files we don't need but each rank does open
!           all the files it will need to write to....lots of files=weakness
! 
!
!           write the number of vtx being mapped from part j to part i
!           note ndof will need to change to the total size of vector field being mapped
!           and this also should be computed in some routine above for allocation 
!           purposes.  currently it is at risk of being confused with ndof
!
            itmp=1
            if (l .gt. 0) itmp = int(log10(float(l)))+1
            write (fmt1,"('(''restartMap.'',i',i1,',1x)')") itmp
            write (fnamer,fmt1) l
            fnamer = trim(fnamer) // cname2(myrank+1)
            if (isbinary == 1) then      
              open(unit=700+l,file=fnamer,status='unknown', 
     &                              form='unformatted')
            else
              open(unit=700+l,file=fnamer,status='unknown')
            endif
          endif
        enddo
!
!       so that we can read directly into the mapped array later, it is helpful to 
!       find the number of verts on each part in the N partition which is easily 
!       obtained  by looking at the max of each vtx id going to each N rank
!       and then doing an allReduce to find global maximum
!
!       Write solution
!
        do l=1,irankN
          if(icountN(l).gt.0) then  !  I have data to write for rank l
            if (isbinary == 1) then      
               write(700+l) icountN(l), ndof
            else
               write(700+l,*) icountN(l), ndof
            endif
          endif
        enddo
 
        do i=1,nshg
          ifileN=pID(i)+1.1    !shift pID to start from 1
          intvID=(vID(i)+1.1)  !shift vID to start from 1
          if (isbinary == 1) then      
            write(700+ifileN) intvID,(qold(i,j),j=1,ndof)
          else
            write(700+ifileN,925) intvID,(qold(i,j),j=1,ndof)
          endif
!         the above "dealing" of records to files based on pID is safe because each 
!         process 1..M has a unique set of files to write to
        enddo

!
!       Write dwal
!
        ione = 1
        do l=1,irankN
          if(icountN(l).gt.0) then  !  I have data to write for rank l
            if (isbinary == 1) then      
               write(700+l) icountN(l), ione
            else
               write(700+l,*) icountN(l), ione
            endif
          endif
        enddo
 
        do i=1,nshg
          ifileN=pID(i)+1.1    !shift pID to start from 1
          intvID=(vID(i)+1.1)  !shift vID to start from 1
          if (isbinary == 1) then      
            write(700+ifileN) intvID, dwal(i)
          else
            write(700+ifileN,925) intvID, dwal(i)
          endif
!         the above "dealing" of records to files based on pID is safe because each 
!         process 1..M has a unique set of files to write to
        enddo
!
!       close all the files on each of the M ranks and deallocate as we are done
!       with the ascii file creation step
!
        do l=1,irankN
          if(icountN(l).gt.0) close(700+l)
        enddo

      endif ! ireducemethod 
          
!
!  Deallocate some memory
!
      call mpi_barrier(mpi_comm_world, ierr) 
      if(myrank == 0) then
        write(*,*) 'Deallocating some memory'
      endif
      deallocate(vID)
      deallocate(pID)
      deallocate(qold)
      if(iybar == 1) then
        deallocate(ybar)
      endif
      if(ierror == 1) then
        deallocate(errors)
      endif
      if(numphavg .gt. 0) then
        deallocate(yphbar)
      endif
      if(ivort == 1) then
        deallocate(vort)
      endif
      if(idwal == 1) then
        deallocate(dwal)
      endif
!
!  Make sure all the ranks are done with their restartMap files before
!  trying to read them for ireducemethod = 1
!
      call mpi_barrier(mpi_comm_world, ierr)

!
!  Create a subcommunicator for the first N ranks
!
 
!      islesseqN = 0
!      if(myrank.lt.irankN) then
!        islesseqN = 1
!      endif
!      call mpi_comm_split(mpi_comm_world, islesseqN, myrank, 
!     &                    my_local_comm, ierr)

!      my_local_size = -1
!      call mpi_comm_size(my_local_comm, my_local_size, ierr)
!      my_local_rank = -1
!      call mpi_comm_rank(my_local_comm, my_local_rank, ierr)
!      write(*,*) "RANKS:", myrank, my_local_rank, my_local_size, irankN

!
!     Now for all N ranks, load the partial files written above and merge them
!
!      if(myrank.lt.irankN) then   ! only parts 1..N do work in this phase

        if(ireducemethod == 1) then ! Read the files back for the file based approach. 
                                    ! Beware, this method is not maintained any more

          l=myrank+1                ! help me find my mapped files 
          do j=1,numpe              ! potentially all M have a map for me
!
!           create the filename for all M of my potential Maps
!
            itmp=1
            if (l .gt. 0) itmp = int(log10(float(l)))+1
            write (fmt1,"('(''restartMap.'',i',i1,',1x)')") itmp
            write (fnamer,fmt1) l
            fnamer = trim(fnamer) // cname2(j)
!
!           most won't exist so inquire for existance
!
            inquire(file=fnamer,exist=exmap)
            if(exmap) then  ! so this map does exist so we open the Map files

              if(isbinary == 1) then
                open(unit=700+l,file=fnamer,status='unknown', 
     &                                    form='unformatted')
              else
                open(unit=700+l,file=fnamer,status='unknown')
              endif
!
!             Read solution first
!
              if(isbinary == 1) then
                read(700+l) icountRead, ndofRead  ! read it's header
              else
                read(700+l,*) icountRead, ndofRead  ! read it's header
              endif
              allocate(qreadN(ndofRead))
              do i=1, icountRead                  ! loop over its contents
                if(isbinary == 1) then
                  read(700+l) intvID,(qreadN(k),k=1,ndofRead) ! read a line in binary format
                else
                  read(700+l,925) intvID,(qreadN(k),k=1,ndofRead) ! read a line in ascii format
                endif
! DEBUG
!              if((i.eq.1).and.(myrank.eq.0)) 
!     &          write(*,*) vidRead, ndofRead, (qreadN(k),k=1,ndofRead)
! DEBUG
                qoldN(intvID,1:ndof)=qreadN(1:ndof) ! fill directly into qoldN
! 
!               the above avoids sorting and sifting as duplicated vtx values just overwrite 
!               and you get the last one as a final value.  Assuming the field was commu-ed
!               all duplicat values are the same so this is safe.  For error fields this 
!               might not be true but it should be o.k.
!
! DEBUG
!              if((i.eq.1).and.(myrank.eq.0)) 
!     &        write(*,*) intvID,(qoldN(intvID,k),k=1,ndofRead)
! DEBUG
              enddo ! end of this map file so close the file and look for more
              deallocate(qreadN)
!
!             Read dwalN second
!
              if(isbinary == 1) then
                read(700+l) icountRead, ndofRead  ! read it's header
              else
                read(700+l,*) icountRead, ndofRead  ! read it's header
              endif
              do i=1, icountRead
                if(isbinary == 1) then
                  read(700+l) intvID,qreadN1 
                else
                  read(700+l,925) intvID,qreadN1
                endif
                dwalN(intvID)=qreadN1 ! fill directly into dwalN
              enddo ! end of this map file so close the file and look for more

              close(700+l)
            endif ! map existed
          enddo ! loop over all M ranks

        endif ! ireducemethod = 1
!
!       At this point I should have a complete field for my part so I write the 
!       mapped file
!

! DEBUG
!        itmp=1
!        if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
!        write (fmt1,"('(''restartMapped.'',i',i1,',1x)')") itmp
!        write (fnamer,fmt1) lstep
!        fnamer = trim(fnamer) // cname2(myrank+1)
!        open(unit=700,file=fnamer,status='unknown') !,
!        write(700,*) maxvIDout(myrank+1),ndof
!        do i=1,maxVIDout(myrank+1)+1
!          write(700,927) i, (qoldN(i,j),j=1,ndof)
!          enddo
!        close(700)
! DEBUG

!
!      Write the N posix files. Note that the nodes on the part boundaries still need to be updated
!
!      Posix format
!        call write_restartonly(myrank, lstep, nshgN, ndof, qoldN)
!        call write_field (myrank,'a','dwal',4,dwalN,'d',nshgN,1,lstep)
!        call write_field (myrank,'a','errors',6,errorsN,'d',
!     &                                        nshgN,ndoferrors,lstep)
!        call write_field (myrank,'a','ybar',4,ybarN,'d',
!     &                                           nshgN,ndofybar,lstep)

!        SyncIO format
          

!
!      Test the result from mpi_alltoall
!
!      qoldN=-9.87654321e32
!      l=myrank+1
!      do j=1,numpe ! should really go to irankN
!        do i=1,maxicountNglob
!          intvID = qTota2ar( (i-1)*numvar+1 , j) + 1.1
!          if (intvID > 0) then
!            do idof=1,ndof 
!              qoldN(intvID,idof) = qTota2ar( (i-1)*numvar+1+idof, j)
!            enddo
!          endif
!        enddo
!      enddo

!      itmp=1
!      if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
!      write (fmt1,"('(''restartMapped2.'',i',i1,',1x)')") itmp
!      write (fnamer,fmt1) lstep
!      fnamer = trim(fnamer) // cname2(myrank+1)
!      open(unit=700,file=fnamer,status='unknown') !,
!      write(700,*) maxvIDout(myrank+1),ndof
!      do i=1,maxVIDout(myrank+1)+1
!        write(700,927) i, (qoldN(i,j),j=1,ndof)
!      enddo
!      close(700)
     
!      endif ! my rank is within range 1..N

      call mpi_barrier(mpi_comm_world, ierr) 
      if(myrank == 0) then
        write(*,*) 'Writing the reduced restartRedTmp-dat files'
      endif

!      call Write_M2N(myrank, irankN, lstep, nshgN, 
!     &          ndof, ndofybar, ndoferrors,
!     &          qoldN, ybarN, errorsN, dwalN)


      nsynciofieldswriterestart = 1+iybar+ierror+numphavg+ivort+idwal

!     We check here how many geombcRed files have been red.
!     If 0, then two possibilities:
!       - Reduction to 1 part so geombcRed are not required. Set nsynciofilesred to 1
!       - User forgot to create the link to geombcRed. 
!         Then, assumes the number of geombc and geombcRed files ia the same  

      if(nsynciofilesred == 0) then
        if(irankN == 1) then ! 1 part reduction -> set to 1
          nsynciofilesred = 1
        else ! Links are missing
          nsynciofilesred = nsynciofiles
        endif
      endif

      call Write_M2N_SolOnly(myrank,irankN,lstep,nshgN,ndof,qoldN)
      deallocate(qoldN)

      if(iybar == 1) then
        call Write_M2N_Field(myrank,irankN,'a','ybar',4,ybarN,'d',
     &                    nshgN,ndofybar,lstep) 
        deallocate(ybarN)
      endif

      if(ierror == 1) then
        call Write_M2N_Field(myrank,irankN,'a','errors',6,errorsN,'d',
     &                    nshgN,ndoferrors,lstep) 
        deallocate(errorsN)
      endif

      if(numphavg .gt. 0) then
        do iphavg=1,numphavg
          call write_M2N_phavg2(myrank,irankN,'a','phase_average',13,
     &             iphavg,numphavg,yphbarN(:,:,iphavg),'d',
     &             nshgN,ndofyphbar,lstep)
        enddo
        deallocate(yphbarN)
      endif

      if(ivort == 1) then
        call Write_M2N_Field(myrank,irankN,'a','vorticity',9,
     &                    vortN,'d',nshgN,ndofvort,lstep) 
        deallocate(vortN)
      endif

      if(idwal == 1) then
        call Write_M2N_Field(myrank,irankN,'a','dwal',4,dwalN,'d',
     &                    nshgN,1,lstep) 
        deallocate(dwalN)
      endif


!
!     Below is debug stuff trying to figure out why the code crashes on exit
!

!       call mpi_barrier(mpi_comm_world, ierr)
!       write(*,*) myrank,'leaving readnblk'

!925   format(i5,2x,6(e14.7,2x)) 
925   format(i9,1x,6(e24.17,1x)) ! used for file based approach

926   format(9(i5,2x))  !used for debugging
927   format(i5,2x,6(e14.7,2x))  !used for debugging

      return
!
!
      end

