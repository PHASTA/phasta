      subroutine forces (y ,  ilwork)
c
c----------------------------------------------------------------------
c
c This subroutine calculates and updates the aerodynamic forces
c
c----------------------------------------------------------------------
c
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
c     
      integer ilwork(nlwork)

      real*8  y(nshg,ndof)
      real*8  Forin(5), Forout(5),spmasss
      real*8  ftots(3,0:MAXSURF),ftot(3)

!MR CHANGE
        integer :: iv_rankpernode, iv_totnodes, iv_totcores
        integer :: iv_node, iv_core, iv_thread
!MR CHANGE
 
c
c  START OF DISTURBANCE BLOCK
c
      dimension diste(1), distesum(1)
      distcalc=0  ! we would want to activate this for T-S instability studies
      if(distcalc.eq.1) then
        if(iter.eq.nitr) then  ! we have completed the last iteration
          if (numpe > 1) then
          call MPI_REDUCE (dke,dkesum,1,MPI_DOUBLE_PRECISION,
     &                     MPI_SUM, master, MPI_COMM_WORLD,ierr)
          endif
          if (numpe.eq.1)
     &      write (76,1000) lstep+1, dke,dke

          if ((myrank .eq. master).and.(numpe > 1)) then 
             write (76,1000) lstep+1, dkesum,dkesum
c
             call flush(76)
c     
          endif
          dke=0.0               ! we must zero it back out for next step

       endif
       endif

C
C END OF DISTURBANCE BLOCK
C

c
c.... -------------------->  Aerodynamic Forces  <----------------------
c
c.... output the forces and the heat flux
c
      if (numpe > 1) then
         Forin = (/ Force(1), Force(2), Force(3), HFlux, 
     &              entrop /)
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE (Forin(1), Forout(1),5,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         Force = Forout(1:3)
         HFlux = Forout(4)
         entrop= Forout(5)
      endif

      if (numpe > 1) then
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE (flxID(1,isrfIM), Atots,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE (flxID(2,isrfIM), spmasss,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE (flxID(3,:), Ftots(1,:),MAXSURF+1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE (flxID(4,:), Ftots(2,:),MAXSURF+1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE (flxID(5,:), Ftots(3,:),MAXSURF+1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
      else
         Ftots=flxID(3:5,:)
         Atots=flxID(1,isrfIM)
         spmasss=flxID(2,isrfIM)
      endif
      call sforce()

      ftot(1)=sum(Ftots(1,0:MAXSURF))
      ftot(2)=sum(Ftots(2,0:MAXSURF))
      ftot(3)=sum(Ftots(3,0:MAXSURF))

      if (myrank .eq. master) then
         write (iforce,1000) lstep, (Force(i), i=1,nsd), 
     &        HFlux,spmasss,(ftot(i),i=1,nsd)
         call flush(iforce)
      endif
   
      if((spmasss.ne.0) .and. (matflg(5,1).ne.0))then ! we are doing
                                                      ! force adjustment 
         tmp=Dtgl*0.025         ! .025= 1/2/tmp1 when tmp1=20
         
         select case ( matflg(5,1) )
         case ( 1 )             ! standard linear body force
            fwall=Force(1)
            vlngth=xlngth
         case ( 2 )       
         case ( 3 )     
            fwall=Force(3)
            vlngth=zlngth
         end select

         write(222,*) spmasss
         fnew=(vel-spmasss/(ro*Atots))*tmp + 
     &        fwall/(vlngth*Atots)
         if(myrank.eq.0) then
            write(880,*) datmat(1,5,1),fnew
            call flush(880)
         endif
         datmat(1,5,1)=fnew     !* (one - tmp2) + datmat(1,5,1) * tmp2
      endif

      if (pzero .eq. 1) then
c

c  also adjust the pressure here so that the mean stays at zero (prevent drift)
c
c  This is only valid if there are NO Dirichlet Boundary conditions on p!!!
c
c  (method below counts partition boundary nodes twice. When we get ntopsh
c  into the method it will be easy to fix this.  Currently it would require
c  some work and it is probably not worth it 
c
c  switched to avg of the on-processor averages to get around hierarchic
c  modes problem
c
c              pave=sum(y(1:numnp,1)) 
c              call MPI_ALLREDUCE (pave, pavet, 1,
c     &             MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
c     pavet=pavet/nshgt
         pave=sum(y(1:numnp,4)) 
         if (numpe .gt. 1) then
            xnpts=numnp
           call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE (pave, pavet, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE (xnpts, xnptst, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            pavet=pavet/xnptst 
         else
            pavet = pave/numnp
         endif

         y(1:numnp,4)=y(1:numnp,4)-pavet
         pave=sum(y(1:numnp,4))
c         if(myrank.eq.0) then
c            write(90+myrank,*) pavet,pave
c            call flush(900)
c         endif
c
      endif

      return
c
! 1000 format(1p,i6,8e13.5)
 1000 format(1p,i6,8e20.12)
c
      end

      subroutine sforce()

      include "common.h"
      include "mpif.h"

      real*8  SFlux(10),SFluxg(10)
      character*20 fname1
      character*5  cname
      character*10  cname2

      integer icalled
      data icalled /0/
      save icalled

      if(icalled.eq.0) then
        icalled = 1

!MR CHANGE
        iv_rankpernode = iv_rankpercore*iv_corepernode
        iv_totnodes = numpe/iv_rankpernode
        iv_totcores = iv_corepernode*iv_totnodes
        irankfilesforce(:) = -1
!MR CHANGE

        do isrf = 0,MAXSURF

!MR CHANGE
          ! Compute the adequate rank which will take care of isrf
          iv_node = (iv_totnodes-1)-mod(isrf,iv_totnodes)
          iv_core = (iv_corepernode-1) - mod((isrf - 
     &              mod(isrf,iv_totnodes))/iv_totnodes,iv_corepernode)
          iv_thread = (iv_rankpercore-1) - mod((isrf- 
     &              (mod(isrf,iv_totcores)))/iv_totcores,iv_rankpercore)
          irankfilesforce(isrf) = iv_node*iv_rankpernode 
     &                     + iv_core*iv_rankpercore
     &                     + iv_thread
                 
          if(myrank == 0 .and. isrf.le.iv_rankpernode) then
              write(*,*) '  sforce ', isrf, 'if any handled by rank',
     &                      irankfilesforce(isrf), ' on node', iv_node
          endif

          ! Verification just in case
          if(irankfilesforce(isrf) .lt.0 .or. 
     &                 irankfilesforce(isrf) .ge. numpe) then 
            write(*,*) 'WARNING: irankfilesforce(',isrf,') is ',
     &                                  irankfilesforce(isrf),
     &                      ' and reset to numpe-1'
            irankfilesforce(isrf) = numpe-1
          endif
!MR CHANGE

!          if ( nsrflist(isrf).ne.0 .and. myrank.eq.master) then
          if ( nsrflist(isrf).ne.0 .and.
     &                     myrank.eq.irankfilesforce(isrf)) then
            iunit=60+isrf
            fname1 = 'forces_s'//trim(cname(isrf))// cname2(lskeep+1)
            open(unit=iunit,file=trim(fname1),status='unknown')
          endif
        enddo
      endif

      do isrf = 0,MAXSURF
        if ( nsrflist(isrf).ne.0 ) then
          SFlux(:) = flxID(:,isrf) ! Area, Mass Flux, Force 1, Force 2, Force 3, Scalars
          if ( numpe > 1 ) then
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_ALLREDUCE (SFlux, SFluxg,10,
     .        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
            SFlux = SFluxg
          endif
!          if(myrank.eq.master) then
          if(myrank.eq.irankfilesforce(isrf)) then
            iunit=60+isrf
!            write(iunit,"(i7,1p10e14.5)") lstep, SFlux
            write(iunit,"(i7,1p10e20.12)") lstep, SFlux
            call flush(iunit)
! done in itrdrv.f
!            if(istep.eq.nstep(1)) then !end step then close file
!              close(iunit)
!            endif
          endif
        endif
      enddo

      return
      end
