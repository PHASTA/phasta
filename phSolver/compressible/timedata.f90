      !--------------------------------------------------------
      ! Initialize the Probe Point Arrays and write the Header
      ! -------------------------------------------------------
      subroutine initProbePoints()
        !Tests if the probe point file xyzts.dat exists, loads probe
        !point locations, initializes a number of arrays used by
        !timedataC, and writes the initial header for the output file. 
        !
        ! Rewritten by:             Nicholas Mati  2014-04-18
        ! Revision history:
        !  2014-04-18   Code moved from itrdrv to here

        use timedataC
        include "common.h"
        include "mpif.h"

        logical :: exVarts 
        
        !Test if xyzts.dat exists and broadcast the result. 
        if(myrank.eq.master) inquire(file='xyzts.dat',exist=exts)
        if(numpe .gt. 1) then
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call MPI_BCAST(exts,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
        endif
        
        if(.not. exts) return
        call readProbePoints
          
        allocate (statptts(ntspts,2))
        allocate (parptts( ntspts,nsd))
        allocate (varts(   ntspts,ndof))

        statptts(:,:) = 0
        parptts(:,:) = zero
        varts(:,:) = zero     
        ivartsbuff = 0
        vartsResetbuffer = .false.
        
        allocate (ivarts(    ntspts*ndof))
        allocate (ivartsg(   ntspts*ndof))
        allocate (vartssoln( ntspts*ndof))
        allocate (vartssolng(ntspts*ndof))
        allocate (vartsbuff( ntspts,ndof,nbuff))
        allocate (vartsbuffstep(nbuff))

        !test if the varts folder exists. If it doesn't create it. 
        if(myrank .eq. master) then
          inquire(file="./varts/.", exist=exVarts)
          if(.not. exVarts) then
            call system("mkdir varts")    !Only works on *nix, but we
                                          !never really run on Windows
                                          !anymore so...
          endif
        endif

!       initProbePoints = exts
!     end function
      end subroutine


      !------------------------
      ! Read Probe Point Input
      !------------------------
      subroutine readProbePoints
        ! Original Code written by: ??             ????-??-??
        ! Rewritten by:             Nicholas Mati  2014-04-18
        ! Revision history:
        !  2014-04-18   Rewritten code moved from itrdrv to here. 
        !
        !Reads the file xyzts.dat for probe point locations, write
        !frequency, tolerance, ... The file is expected to have the
        !form:
        ! ntspts freq tolpt iterat nbuff
        ! x1 y1 z1
        ! x2 y2 z2
        ! ...
        ! xN yN zN
        !
        ! ... where ntspts is the number of probe points and freq is the
        ! number of steps to take before flushing data. If nbuff is
        ! zero, the number of time steps between restarts, ntout, is
        ! used. 

        use timedataC
        include "common.h"
        include "mpif.h"
      
        if(myrank.eq.master) then
          open(unit=626,file='xyzts.dat',status='old')
          read(626,*) ntspts, freq, tolpt, iterat, nbuff
        endif
        
        !Broadcase out the header of xyzts.dat. These should probably
        !be combined into two calls, but this is quick and dirty. 
        if(numpe .gt. 1) then
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call MPI_Bcast(ntspts, 1, MPI_INTEGER,          master,
     &                              MPI_COMM_WORLD,       ierr)
          call MPI_Bcast(freq,   1, MPI_INTEGER,          master, 
     &                              MPI_COMM_WORLD,       ierr)
          call MPI_Bcast(tolpt,  1, MPI_DOUBLE_PRECISION, master, 
     &                              MPI_COMM_WORLD,       ierr)
          call MPI_Bcast(iterat, 1, MPI_INTEGER,          master, 
     &                              MPI_COMM_WORLD,       ierr)
          call MPI_Bcast(nbuff,  1, MPI_INTEGER,          master, 
     &                              MPI_COMM_WORLD,       ierr)
        endif

        allocate (ptts(    ntspts,nsd))

        !Read probe point coordinates and broadcast to the rest of the
        !processors
        if(myrank .eq. master) then
          do jj=1,ntspts        ! read coordinate data where solution desired
             read(626,*) ptts(jj,1), ptts(jj,2), ptts(jj,3)
           enddo
         close(626)
        endif
        
        if(numpe .gt. 1) then
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call MPI_Bcast(ptts, ntspts*nsd, MPI_DOUBLE_PRECISION,
     &                          master,     MPI_COMM_WORLD,       ierr)
        endif
      
        if (nbuff .eq. 0) 
     &    nbuff=ntout
      end subroutine

      
      !-----------------------------
      ! Write the Header varts file
      !-----------------------------
      subroutine TD_writeHeader(fvarts)
        !Creates the file fvarts and writes the data header.  
        !fvarts:    Name The file to create
        
        use timedataC
        include "common.h"

        character(len=*) fvarts

        !Open the output varts file and write the header 
         if (myrank .eq. master) then
           
           !fvarts='varts/varts'
           !fvarts=trim(fvarts)//trim(cname2(lstep))
           !fvarts=trim(fvarts)//'.dat'
           open(unit=1001, file=fvarts, status='unknown')
           
           !Write the time step        
           write(1001, *) "Time Step: ", Delt(1) 
           write(1001, *)

           !Write the probe locations to varts.ts.dat so that post
           !processing tools actually know what point goes where.
           !From experience, it's difficult to keep this straight.
           write(1001, *) 
     &                 "Probe ID   x              y              z"
           do jj = 1, ntspts 
             write(1001, "(I5, T12, 3(F16.12))") jj, ptts(jj,1:3)
           enddo
           write(1001, *)      
                
           !write the output format string. This can't be hard
           !coded because ntspts is not known in advance. 
           write(vartsIOFrmtStr, '("(I8, ", I4, "(E15.7e2))")')
     &           ndof*ntspts
           
           !Header to delinieate the probe locations with the data.
           write(1001, *) "Probe Data:"
           close(unit=1001)
         endif  ! if(myrank .eq. master)
      end subroutine
       
      

      !------------------------
      ! Accumulate Probe Data
      !------------------------
      subroutine TD_bufferData()

        use timedataC
        include "common.h"
        include "mpif.h"

        integer :: icheck, istp, nstp
                  
        if (mod(lstep,freq).eq.0) then
          if(vartsResetBuffer) then
            ivartsBuff = 0
            vartsResetBuffer = .false.
          endif

          !------------------------
          !Merge Data across parts
          !------------------------
          if (numpe > 1) then
            !load the contents of varts into vartssoln
            do jj = 1, ntspts
               vartssoln((jj-1)*ndof+1:jj*ndof)=varts(jj,:)
               ivarts=zero
            enddo

            !mark which points have been computed on this processor
            do k=1,ndof*ntspts
               if(vartssoln(k).ne.zero) ivarts(k)=1
            enddo

            !merge the solution
            call MPI_REDUCE(vartssoln, vartssolng, ndof*ntspts,
     &           MPI_DOUBLE_PRECISION, MPI_SUM, master,
     &           MPI_COMM_WORLD, ierr)
              
            call MPI_REDUCE(ivarts, ivartsg, ndof*ntspts,
     &           MPI_INTEGER, MPI_SUM, master,
     &           MPI_COMM_WORLD, ierr)
                 
             !if the probe point happened to span multiple partitions,
             !divide the sum by the number of contributing partitions. 
             if (myrank.eq.master) then
               do jj = 1, ntspts
                 indxvarts = (jj-1)*ndof
                 do k=1,ndof
                   if(ivartsg(indxvarts+k).ne.0) then ! none of the vartssoln(parts) were non zero
                      varts(jj,k) = 
     &                    vartssolng(indxvarts+k) / ivartsg(indxvarts+k)
                   endif
                 enddo !loop over states / DoF
               enddo !loop over probe points
             endif !only on master
          endif !only if numpe > 1
          
          ivartsBuff = ivartsBuff + 1
          if (myrank.eq.master) then
            if(ivartsBuff .gt. nbuff) then
              write(*,*) "WARNING: vartsbuff has overflowed"
              ivartsBuff = nbuff
            endif

            vartsBuffStep(ivartsBuff) = lstep
            do jj = 1, ntspts
              vartsbuff(jj,1:ndof, ivartsBuff) = varts(jj,1:ndof)
            enddo
          endif
        endif

        varts(:,:) = zero

      end subroutine 
      

      !------------
      ! Write Data 
      !------------
      subroutine TD_writeData(fvarts, forceFlush)
        !writes the probe point data accumulated durring calls to
        !TD_bufferData(). Note that actual file IO only occurs when the
        !buffer is full or when DT_writeData is called with forceFlush
        !set to true. Also note that TD_writeHeader must be called prior
        !to calling DT_writeData. 
        use timedataC
        include "common.h"

        character(len=*) :: fvarts
        logical :: forceFlush
!        logical, optional :: forceflush    
        logical :: flush
        integer :: k, ibuf

        if (myrank.eq.master) then

          !if provided, use the default value passed in to determine
          !wheather to flush the buffer
!          if(present(forceFlush)) then   !optional version breaks the
            flush = forceFlush            !compiler on Bluegene? 
!          else
!            flush = .false.    !set the default value 
!          endif

          !make sure incomplete buffers get purged at the end of a run
          !regardless of the default.
!         if(ivartsBuff .eq. nbuff) flush = .true.
          if(mod(lstep, nbuff) .eq. 0) flush = .true.
          if(vartsResetBuffer) flush = .false.  !Prevent repeated calls without updating 
                                                          !the buffer from writting multiple times. 
                  
          if(flush) then   !flush the buffer to disc
            open(unit=1001, file      = fvarts,   status = "old",
     &                  position = "append", action = "write")
            do ibuf = 1,ivartsBuff    
              write(1001, vartsIOFrmtStr)
     &                    vartsBuffStep(ibuf),                  !write the time step in the first column. 
     &                   ((vartsbuff(jj,k,ibuf),  k=1, ndof)   !loop over the variables that you actually want to output. 
     &                                         , jj=1, ntspts) !loop over probe points
            enddo
                      
            close(1001)

            vartsResetBuffer = .true.
!            ivartsBuff = 0      !need to reset ivartsBuff because
!                                !writeDate can be called consecutively
          endif !only dump when buffer full
        endif !only on master
                      
!              call flush(1001)
!              call fsync(1001)
                     
               !Code for writting one file per probe point
!              do jj = 1, ntspts !loop through probe points
!                ifile = 1000+jj
!                do ibuf=1,nbuff 
!                  write(ifile,555) lstep-1 -nbuff+ibuf, 
!     &               (vartsbuff(jj,k,ibuf) , k=1, ndof)
!!     &              vartsbuff(jj,:,ibuf)
!                 
!                enddo ! buff empty
!                  
!                call flush(ifile) 
!              enddo ! jj ntspts
              
            
!         varts(:,:) = zero ! reset the array for next step   !MOVED FOR Mach Control                 
! 555     format(i6,6(2x,E12.5e2))

      end subroutine


      subroutine TD_finalize()
       use timedataC
         
        deallocate(ivarts)
        deallocate(ivartsg)
        deallocate(vartssoln)
        deallocate(vartssolng)
        deallocate(vartsbuff)
        deallocate(vartsbuffstep)

        deallocate(ptts)
        deallocate(varts)
      end subroutine


      !---------------------
      ! allocate the arrays
      !---------------------
      subroutine sTD 
        !Allocates the arrays statptts, ptts, parptts, and varts for use
        !in itrdrv and ??
        !Subroutine is Depricated. 
       
       use timedataC
        include "common.h"
  
        allocate (statptts(ntspts,2))
        allocate (ptts(ntspts,nsd))
        allocate (parptts(ntspts,nsd))
        allocate (varts(ntspts,ndof))
  
        return
      end

      !-------------------
      ! delete the arrays
      !-------------------
      subroutine dTD 
        !Deallocates ptts and varts
       use timedataC
  
        deallocate (ptts)
        deallocate (varts)
  
        return
      end
