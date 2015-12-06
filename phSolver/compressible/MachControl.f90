      module MachControl
        logical exMC            !Mach Control enable / Does the MC file
                                !exist

        integer :: MC_iPP       !Probe Point index
        integer :: MC_surfID    !Surface ID to identify where to apply
                                ! the boundary condition 

        real*8 :: MC_dt         !simulaiton time step
        real*8 :: MC_MtTarget   !Target throat Mach number
        real*8 :: MC_dMtCO      !M_throat Cut Off (frequency) (for
                                ! filtering the state)
        real*8 :: MC_ulim(2)    !u_inlet velocity limit
                                !u_inlet(1)     lower limit
                                !u_inlet(2)     upper limit
        !States
        real*8 :: MC_Mt         !actual Mach number
        real*8 :: MC_dMtdt      !time derivative of MC_dMt
        real*8 :: MC_dMt        !filtered state deviation from target
        real*8 :: MC_IdMt       !integreal of MC_dMt

        !Gains
        real*8 :: MC_Kd         !derivative gain
        real*8 :: MC_Kp         !proportional gain
        real*8 :: MC_KI         !integral gain

        !Variables to apply the boundary condition
        integer, allocatable, dimension(:) :: MC_imapInlet(:) 
        integer :: MC_icountInlet   !number of inlet boundary points 
        real*8, allocatable, dimension(:) :: MC_vscale(:)   !scaling factor for BL
        real*8 :: MC_BLdelta        !boundary layer thickness
        real*8 :: MC_uInlet     !x-velocity at the inlet 
        
      end module

      subroutine MC_init(dt, tstep, BC)
        use MachControl
        use turbSA
        include "common.h"
        include "mpif.h"

        dimension :: BC(nshg, ndofBC)
        logical :: MC_existDir
        real*8 :: dt
        real*8 :: vBC, vBCg
        integer :: tstep
        vBCg = 0
        
        if(myrank == master) then 
          inquire(file="./MachControl/.", exist=MC_existDir)
          if(.not. MC_existDir) then
            call system("mkdir ./MachControl")  !Doesn't seem to work on BGQ.
          endif
        endif

        call MC_readState(tstep)
        
        if(exMC) then
          !find the node associated with isetBlowerID_Duct and write the
          !scale vector to linear ramp velocity near the walls. 
          allocate(MC_imapInlet(nshg))
          allocate(MC_vscale(nshg))
          call sfID2np(MC_surfID, MC_icountInlet, MC_imapInlet)       
             
          if(MC_BLdelta > 0) then
            do i = 1, MC_icountInlet
              imapped = MC_imapInlet(i)
              
              MC_vscale(i) = min(one, d2wall(imapped)/MC_BLdelta)
              !MC_vscale(i) = d2wall(imapped)/MC_BLdelta
              !if(MC_vscale(i) > 1) MC_vscale(i) = 1
  
              !Also look for the inlet velocity. This will get 
              !overwritten in the next time step, but it's nice to know
              !what it started with.
              vBCg = max(vBCg, BC(imapped, 3))
            enddo 
          endif
       
          !Look for the maximum global inlet velocity among all parts.  
          if(numpe > 1) then
            call MPI_ALLREDUCE(vBCg, MC_uInlet, 1, 
     &           MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
          endif
  
          MC_dt = dt
        endif
      end subroutine


      subroutine MC_updateState()
        !Updates the derivative, proportional and integral values used
        !when calculating MC_uInlet. 
        !Note: This function uses the most recent values in the TimeData
        !buffer vartsbuff. TD_bufferData should be called prior to this
        !subroutine. 
        
        use MachControl
        use timedataC
        include "common.h"
        include "mpif.h"
        
        real*8 :: u2, T
        real*8 :: Mt, dMt, MC_alpha
        logical :: bangbang
        
        if(myrank == master) then !MC_* are not set on other processors
          T  = vartsbuff(MC_iPP, 5, ivartsBuff) 
          u2 = vartsbuff(MC_iPP, 2, ivartsBuff)**2 + 
     &         vartsbuff(MC_iPP, 3, ivartsBuff)**2 +
     &         vartsbuff(MC_iPP, 4, ivartsBuff)**2  

          MC_Mt = sqrt(u2/(1.4*Rgas*T))
          dMt = MC_Mt - MC_MtTarget
         
          !Note: MC_uInlet should always equal the previous time step;
          !on the first time step of a run, MC_uInlet is computed based
          !off of the maximum inlet velocity in the solution. 
          bangbang = MC_uInlet <= MC_uLim(1) .or.
     &               MC_uInlet >= MC_uLim(2)

          !Low pass filter dMt
          MC_alpha = MC_dt/(MC_dMtCO + MC_dt)
          dMt = MC_alpha*dMt + (1 - MC_alpha)*MC_dMt
                  
          !Calculate derivatives and integrals
          if(.not. bangbang) then  
            !Only integrate when not hitting the limits. This introduces
            !a nonlinearity in the control scheme which minimizes the
            !overshoot when starting from zero initial conditions
            MC_IdMt  = (dMt + MC_dMt)*MC_dt/2 + MC_IdMt !trapezoidal rule
          endif
          MC_dMtdt = (dMt - MC_dMt)/MC_dt !first order BE FD derivative
          MC_dMt = dMt !update the actual state
        endif
      end subroutine


      subroutine MC_applyBC(BC)
        use MachControl
        include "common.h"
        include "mpif.h"

        dimension :: BC(nshg, ndofBC)
        real*8 :: vInlet
        integer :: imapped
        real*8 :: tmp

        if(myrank .eq. master) then
          MC_uInlet = MC_Kd*MC_dMtdt + 
     &                MC_Kp*MC_dMt + 
     &                MC_KI*MC_IdMt
               
          !clip to [MC_uLim(1),MC_uLim(2)] 
          if(MC_uLim(1) < MC_ulim(2)) then 
            MC_uInlet = min(max(MC_uInlet, MC_uLim(1)), MC_uLim(2))
          endif
        endif
        if(numpe .gt. 1) 
     &    call MPI_Bcast(MC_uInlet, 1, MPI_DOUBLE_PRECISION, 
     &                         master, MPI_COMM_WORLD, ierr)
      
        if(MC_BLdelta .gt. 0) then !Apply a velocity profile with a BL
          do i = 1, MC_icountInlet
            imapped = MC_imapInlet(i)
            BC(imapped,3) = MC_uInlet*MC_vscale(i)
          enddo
        else   !Apply a uniform inlet velocity. 
          do i = 1, MC_icountInlet
            BC(MC_imapInlet(i), 3) = MC_uInlet
          enddo
        endif
        
      end subroutine
         
      
      subroutine MC_printState()
      !Prints the Mach Control states into the output
        
        use MachControl
        include "common.h"
        
        if(myrank == master) then
          write(*, "('M_throat: ',F8.6,'  dM_throat: ',F8.6)")
     &             MC_Mt, MC_Mt - MC_MtTarget
          write(*,"('dM/dt = ',E14.6,'  M  = ',E14.6,'  I(M) =',E14.6)")
     &             MC_dMtdt, MC_dMt, MC_IdMt
          write(*,"('Kd    = ',E14.6,'  Kp = ',E14.6,'  KI   =',E14.6)")
     &            MC_Kd, MC_Kp, MC_KI
          write(*, "('u_inlet: ',F9.6)") MC_uInlet
        endif
      end subroutine
 
      subroutine MC_readState(tstep)
        !Reads the intput / state file for MachControl. 
        !Written by: Nicholas Mati      2014-04-19
        !Revision History:
        ! - 2014-04-19: created
        !
        !The file is expected to be of the form:
        ! MC_iPP   MC_surfID  MC_MtTarget
        ! MC_dMtdt MC_dMt MC_IdMt
        ! MC_Kd    MC_Kp  MC_KI 
        ! MC_dMtCO MC_BLdelta

        use MachControl
        include "common.h"
        include "mpif.h"
        
        integer :: tstep 
        character(len=60) :: fname
         
        if(myrank == master) then
          write(fname, "('./MachControl/MachControl.',I0,'.dat')") tstep
          inquire(file=fname, exist=exMC)
  
          if(exMC) then 
            open(unit=1002, file=fname, status='old') 
                     
            read(1002, *)  MC_iPP,   MC_surfID
            read(1002, *)  MC_Mt,    MC_MtTarget
            read(1002, *)  MC_dMtdt, MC_dMt,    MC_IdMt
            read(1002, *)  MC_Kd,    MC_Kp,     MC_KI
            read(1002, *)  MC_dMtCO, MC_BLdelta
            read(1002, *)  MC_uLim(1), MC_uLim(2)
            close(1002)
                  
            !If d2wall is not available, then set the inlet BL thickness 
            !to zero to disable wall scalling. 
            if(iRANS >= 0) MC_BLdelta = 0
          endif 
        endif
           
        !Only a few values are needed across other processors. 
         if(numpe > 1) then
          call MPI_BCAST(MC_BLdelta, 1, MPI_DOUBLE_PRECISION, 
     &                          master, MPI_COMM_WORLD, ierr)
          call MPI_BCAST(MC_surfID,  1, MPI_INTEGER, 
     &                          master, MPI_COMM_WORLD, ierr)
          call MPI_BCAST(exMC,       1, MPI_INTEGER,  
     &                          master, MPI_COMM_WORLD, ierr)
        endif
              
      end subroutine

      subroutine MC_writeState(tstep)
        !Writes a file with the name MachControl.[ts].dat that can later
        !be read in and restarted from. 
             
        use MachControl
        include "common.h"
            
        integer tstep 
        character(len=60) :: fname
       
        if(myrank == master) then 
          write(fname, "('./MachControl/MachControl.',I0,'.dat')") tstep
          open(unit=1002, file=fname, status='unknown') 
            
          write(1002, "(2(i9))")         MC_iPP,   MC_surfID 
          write(1002, "(2(F18.14))")     MC_Mt,    MC_MtTarget
          write(1002, "(3(F18.14))")     MC_dMtdt, MC_dMt,    MC_IdMt
          write(1002, "(3(F18.10))")     MC_Kd,    MC_Kp,     MC_KI
          write(1002, "(2(F18.14))")     MC_dMtCO, MC_BLdelta
          write(1002, "(2(F18.14))")     MC_uLim(1), MC_uLim(2)
          close(1002)
        endif
      end subroutine
      
      subroutine MC_finalize()
     
        use MachControl 

        deallocate(MC_imapInlet)
        deallocate(MC_vscale)
      end subroutine


