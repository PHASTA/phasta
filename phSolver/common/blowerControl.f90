      module blowerControl
        logical :: BC_enable
        integer :: nBlower
        real*8  :: BC_dt
        real*8  :: BC_istep0
        real*8  :: BC_t

        type :: blowerData
          integer :: surfID,                !surface ID to use for blower inlet
     &               nmap,                  !number of inlet nodes
     &               mode                   !constant = 0, trapezoid = 1, sinusoid = 2
          integer, allocatable :: map(:)    !mapping to the inlet nodes
          real*8,  allocatable :: vscale(:) !scalling map used to create an inlet boundary layer
          logical :: updateBC, !true if the BC dynamically updated
     &               enable    !true to enable blower on with the particular surface ID
          
          !The wave form is expected to look something like:
          !      |--------- (a) ---------|
          !       ________                _______    _ vmax
          !      /        \              /          |
          !     /          \            /           |
          !____/            \__________/            |_ vmin
          !    |-|--------|-|
          !    (b)   (c)  (d)

          !Temporal parameters
          real*8 :: t_cycle                 !period of the cycle (a)
          real*8 :: t_riseTime, t_fallTime  !time of rising edge (b), falling edge (d)
          real*8 :: t_fullOn                !duration pulse is on full (c)
          real*8 :: t0                      !t = 0 reference for wave form
            
          !State values
          real*8 :: vmax, vmin              !maximum, minimum blower velocity
          real*8 :: T                       !uniform temperature at blower inlet
          real*8 :: nu                      !uniform eddy viscosity at blower inlet
          real*8 :: deltaBL                 !Boundary layer thickness
          real*8 :: deltaBLsclr             !Boundary layer thickness for scalar
          
        end type blowerData

        type(blowerData), allocatable :: blower(:)

      contains
        subroutine BC_setVars(nBlowers,     blowerMode,
     &                        surfID,       enable, 
     &                        t_cycle,      t_fullOn,
     &                        t_riseTime,   t_fallTime,
     &                        vmax,         vmin, 
     &                        T,            nu,
     &                        deltaBL,      deltaBLsclr )

          integer :: nBlowers
          integer :: blowerMode(nBlowers),  
     &               surfID(    nBlowers),  enable(      nBlowers)
          real*8  :: t_cycle(   nBlowers),  t_fullOn(    nBlowers),
     &               t_riseTime(nBlowers),  t_fallTime(  nBlowers),
     &               vmax(      nBlowers),  vmin(        nBlowers),
     &               T(         nBlowers),  nu(          nblowers),
     &               deltaBL(   nBlowers),  deltaBLsclr(nBlowers)
        
          nBlower = nBlowers
          
          if(.not. allocated(blower)) then
            allocate(blower(nBlower))
          endif

          do i = 1, nBlower
            blower(i)%mode        = blowerMode(i)
            blower(i)%surfID      = surfID(i)
            blower(i)%enable      = (enable(i) == 1)
            blower(i)%t_cycle     = t_cycle(i)
            blower(i)%t_fullOn    = t_fullOn(i)
            blower(i)%t_riseTime  = t_riseTime(i)
            blower(i)%t_fallTime  = t_fallTime(i)
            blower(i)%vmax        = vmax(i)
            blower(i)%vmin        = vmin(i)
            blower(i)%T           = T(i)
            blower(i)%nu          = nu(i)
            blower(i)%deltaBL     = deltaBL(i)
            blower(i)%deltaBLsclr = deltaBLsclr(i)
          enddo

        end subroutine
         
  
        subroutine BC_iter(BC)
        !Updates BC with the correct velocities for the blower. 
        !Note: itrBC still needs to be called after blowerIter to update
        !      values in yold. 
  
!         use blowerControl
          use wallData
          include "common.h"
          
          integer :: i, iBlower, nNodeMap
          real*8  :: vmagref, vmag
          real*8  :: t
          real*8  :: BC(nshg, ndofBC)
  
          !local variables to make the code shorter
          real*8  :: t0, t1, t2, t3, t4, vmax, vmin
          
          if(.not. BC_enable) return
  
          BC_t = BC_dt*(istep - BC_istep0)  !Update time
  
          do iBlower = 1, nBlower
            if(.not. blower(iBlower)%enable .or. 
     &         .not. blower(iBlower)%updateBC ) cycle

            !------------------------            
            !Setup usefull variables
            !------------------------
            vmax = blower(iBlower)%vmax 
            vmin = blower(iBlower)%vmin

            if(blower(iBlower)%mode == 1) then !Trapezoid
              !The wave form is expected to look something like:
              !       ________                _______    _ vmax
              !      /        \              /          |
              !     /          \            /           |
              !____/            \__________/            |_ vmin
              !    | |        | |          |
              !   t0 t1      t2 t3         t4
    
              t0 = 0.0d00
              t1 = t0 + blower(iBlower)%t_riseTime
              t2 = t1 + blower(iBlower)%t_fullOn
              t3 = t2 + blower(iBlower)%t_fallTime
              t4 = blower(iBlower)%t_cycle
    
              !Calculate the velocity magnitude
              !t = mod(BC_dt*istep + blower(iBlower)%t0, t4)
              t = mod(BC_t - blower(iBlower)%t0, t4)
    
              if(t0 <= t .and. t < t1) then     !rising edge
                vmagref = (vmax - vmin)*(t - t0)/(t1 - t0) + vmin
              elseif(t1 <= t .and. t < t2) then !high
                vmagref = vmax
              elseif(t2 <= t .and. t < t3) then !falling edge
                vmagref = (vmax - vmin)*(t - t3)/(t2 - t3) + vmin
              elseif(t3 <= t .and. t < t4) then !low
                vmagref = vmin
              endif
            elseif(blower(iBlower)%mode == 2) then !sinusoid
              !         _.--._                _.--._   |- vmax
              !      _-'      '-_          _-'      '  | 
              !__..-'            '-..__..-'            |- vmin
              !                                        |
              !      |     |------ t_cycle ------|
              !      t0                  
              
              !t is number of waves since t0
              t = (BC_t - blower(iBlower)%t0)/blower(iBlower)%t_cycle
              vmagref = 0.5*(vmax - vmin)*sin(2*pi*t) + 0.5*(vmax + vmin)
            elseif(blower(iBlower)%mode == 0) then !Constant 
              vmagref = max(vmax, vmin) 
            endif 
            
            !---------- 
            !update BC
            !----------
            do i = 1,blower(iBlower)%nmap
              imapped =      blower(iBlower)%map(i)
              vmag = vmagref*blower(iBlower)%vscale(i) !scale the reference 
                                                       !to create a BL
              BC(imapped, 3) = -vmag*wnorm(imapped, 1)  !x velocity
              BC(imapped, 4) = -vmag*wnorm(imapped, 2)  !y velocity
              BC(imapped, 5) = -vmag*wnorm(imapped, 3)  !z velocity
            enddo  !end loop over blower surface nodes
          enddo
  
        end subroutine
  
  
        subroutine BC_Finalize()
        !Deallocates allocatable variables in blower and then deallocates blower. 
          
          integer :: iBlower
  
          if(BC_enable) then 
            do iBlower = 1,nBlower
              deallocate(blower(iBlower)%map)
              deallocate(blower(iBlower)%vscale)
            enddo
          endif
          
          deallocate(blower) !note that blower gets allocated regardless
                             !of whether BC_enable is true
        end subroutine
  
  
        subroutine BC_ReadPhase(tstep)
          include "common.h"
          include "mpif.h"
           
          integer :: tstep
          integer :: i, nBlowerIn, surfIDin
          character(len=128) :: fname
          logical :: existFname
          real*8  :: phi
          
          if(.not. BC_enable) return  !nothing needs to be written
  
          BC_t = BC_dt*(istep - istep0)  !Update time
          if(myrank == master) then
            !Even if blowerPhase.%i.dat exists, it may not have the same 
            !number of blowers; zero all t0 in blower before proceeding.
            do i = 1,nBlower
              blower(i)%t0 = 0.0d00
            enddo
                     
            !Assemble the file name and test if it exists. 
            write(fname, "('blowerPhase.', I0, '.dat')") tstep
            inquire(file=fname, exist=existFname)
            
            if(existFname) then
              open(unit=1003, file=fname, status='old')
    
              read(1003, *) nBlowerIn
  
              !In the future, the values should probably be matched by
              !surface ID, but I don't have time for that now. 
              do i = 1,nBlowerIn  
                read(1003, *) surfIDin, phi
                blower(i)%t0 = -phi*blower(i)%t_cycle 
                !Note the negative sign; t0 is the reference starting
                !time. For positive phase, the cycle started sometime in
                !the past. 
              enddo
              
              close(1003)
            endif !existFname
          endif !myrank == master
           
          !At this point, t0 on rank 0 should be correct; 
          !broadcast the contents of the file
          if(numpe > 1) then 
            do i = 1,nBlower
              call MPI_BCAST(blower(i)%t0, 1, MPI_DOUBLE_PRECISION, 
     &                              master, MPI_COMM_WORLD, ierr)
            enddo
          endif
        end subroutine
  
  
        subroutine BC_writePhase(tstep)
          include "common.h"
           
          integer :: tstep    !time step to use in the file name
               
          character(len=128) :: fname
          integer :: i
          real*8  :: phi
            
          if(.not. BC_enable) return  !nothing needs to be written
  
          BC_t = BC_dt*(istep - BC_istep0)  !Update time
  
          if(myrank == master) then
            write(fname, "('blowerPhase.', I0, '.dat')") tstep
            open(unit=1003, file=fname, status='unknown')
  
            write(1003, "(i9)") nBlower
            do i = 1,nBlower
              !Calculate the present possition in the blowing cycle
              if(blower(i)%updateBC) then
                phi = mod((BC_t - blower(i)%t0)/blower(i)%t_cycle,1.0d0)
              else
                phi = 0.0d00
              endif
  
              write(1003, "(i9, F18.14)") blower(i)%surfID, phi
            enddo
            
            close(1003)
          endif
        end subroutine
  
!        subroutine BC_setNBlower(nBlowers)
!          integer :: nBlowers
!          nBlower = nBlowers
!        end subroutine
! 
!        subroutine BC_setSurfID(surfID)
!          integer :: surfID(nBlower)
!          do i = 1, nBlower
!            blower(i)%surfID = surfID(i)
!          enddo
!        end subroutine
! 
!        subroutine BC_setEnable(enable)
!          integer :: enable(nBlower)
!          do i = 1, nBlower
!            blower(i)%enable = enable(i)
!          enddo
!        end subroutine
! 
!        subroutine BC_setSurfID(surfID)
!          real*8 :: surfID(nBlower)
!          do i = 1, nBlower
!            blower(i)%surfID = surfID(i)
!          enddo
!        end subroutine
! 
!        subroutine BC_setSurfID(surfID)
!          integer :: surfID(nBlower)
!          do i = 1, nBlower
!            blower(i)%surfID = surfID(i)
!          enddo
!        end subroutine
 
        
      end module 
      
       subroutine BC_init(dt, tstep, BC)
          use blowerControl
          use wallData !wnorm, 
          use turbSA !d2wall
          include "common.h"    
   
          real*8 :: BC(nshg, ndofBC)
  
          integer :: i, iBlower, imapped, tstep 
          real*8  :: deltaBLinv, vmag, dt
          
          !local variables to shorten the code a bit
          integer :: nNodeMap
          real*8  :: vmax !T, nu
          
          !variables that should be set in solver.inp:
          !nBlower + surfID, t_cycle, t_riseTime, t_fallTime, t_fullOn, vmax, vmin, T, nu, deltaBL, enable
          
          !Test if any of the blowers are enabled. If not, then set 
          !BC_enable to false and exit. 
          BC_enable = .false.
          do iBlower = 1, nBlower
            if(blower(iBlower)%enable) then
              BC_enable = .true.
              exit
            endif
          enddo
          if(.not. BC_enable) return
          
          BC_istep0 = istep
          BC_dt = dt
                     
          !check to see if a blowerControl.dat has been written with phase data
          call BC_ReadPhase(tstep)  !t0 is updated in the call. 
          
          do iBlower = 1, nBlower
            !------------------------------------
            !Nodal Mapping and Referrence Arrays
            !------------------------------------
            
            !create the nodal mapping
            allocate(blower(iBlower)%map(nshg))
            call sfID2np(blower(iBlower)%surfID, !analog to sfID2np, except 
     &                   blower(iBlower)%nmap,   !that map is allocated in 
     &                   blower(iBlower)%map)    !the subroutine call. 
!            call asfID2np(blower(iBlower)%surfID, !analog to sfID2np, except 
!     &                   blower(iBlower)%nmap,   !that map is allocated in 
!     &                   blower(iBlower)%map)    !the subroutine call. 
            nNodeMap = blower(iBlower)%nmap                
  
            !Calculate the boundary layer scalling      
            allocate(blower(iBlower)%vscale(nNodeMap))
            
            if(blower(iBlower)%deltaBL > 0.0) then   !use the provided BL
              deltaBLinv = 1/blower(iBlower)%deltaBL 
            else                     !use a default of essentially zero thickness 
              deltaBLinv = 1e16      
            endif
                
            do i = 1, nNodeMap 
              !loop over each surface node and calculate vscale with a 
              !linear ramp. Add more options in the future. 
              imapped = blower(iBlower)%map(i)
              blower(iBlower)%vscale(i) = 
     &                          min(1.0d0, d2wall(imapped)*deltaBLinv)
            enddo
            
            !-------------------------------------
            !Apply persistent boundary conditions
            !-------------------------------------
                 
            !temperature          
            if(blower(iBlower)%T > 0.0) then 
              do i = 1, nNodeMap
                imapped = blower(iBlower)%map(i)
                BC(imapped, 2) = blower(iBlower)%T
              enddo
            endif
          
            !eddy viscosity
            if(blower(iBlower)%nu >= 0.0) then
              if(blower(iBlower)%deltaBLsclr > 0.0) then   !use the provided BL
                deltaBLinv = 1/blower(iBlower)%deltaBLsclr
              else                     !use a default of essentially zero thickness 
                deltaBLinv = 1e16
              endif

              do i = 1, nNodeMap
                imapped = blower(iBlower)%map(i)
                BC(imapped, 7) = blower(iBlower)%nu * 
     &                           min(1.0d0, d2wall(imapped)*deltaBLinv)
              enddo
            endif
                
            !velocity
            blower(iBlower)%updateBC = .true.  !update velocity later with a call to BC_iter. 
          enddo !loop over blower surfaces

          call BC_iter(BC)

          do iBlower = 1, nBlower
            if(blower(iBlower)%t_cycle >= 1.0 .or.  !Unreasonable cycle 
     &         blower(iBlower)%t_cycle <= 0.0 .or.  !periods are provided, 
     &         blower(iBlower)%mode == 0) then      !so assume const. 
              blower(iBlower)%updateBC = .false.
            else
              blower(iBlower)%updateBC = .true.
            endif
          enddo
          
        end subroutine

      
