              subroutine itrdrv (y,         ac,   uold, x,         
     &                   iBC,       BC,         
     &                   iper,      ilwork,     shp,       
     &                   shgl,      shpb,       shglb,
     &                   ifath,     velbar,     nsons ) 
c
c----------------------------------------------------------------------
c
c This iterative driver is the semi-discrete, predictor multi-corrector 
c algorithm. It contains the Hulbert Generalized Alpha method which
c is 2nd order accurate for Rho_inf from 0 to 1.  The method can be
c made  first-order accurate by setting Rho_inf=-1. It uses a
c GMRES iterative solver.
c
c working arrays:
c  y      (nshg,ndof)           : Y variables
c  x      (nshg,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodicity table
c
c shape functions:
c  shp    (nshape,ngauss)        : interior element shape functions
c  shgl   (nsd,nshape,ngauss)    : local shape function gradients
c  shpb   (nshapeb,ngaussb)      : boundary element shape functions
c  shglb  (nsd,nshapeb,ngaussb)  : bdry. elt. shape gradients
c
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use pvsQbi        !gives us splag (the spmass at the end of this run 
      use specialBC     !gives us itvn
      use timedataC      !allows collection of time series
      use MachControl   !PID to control the inlet velocity. 
      use blowerControl !gives us BC_enable 
      use turbSA
      use wallData
      use fncorpmod

        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
      
c
        dimension y(nshg,ndof),            ac(nshg,ndof),  
     &            yold(nshg,ndof),         acold(nshg,ndof),           
     &            x(numnp,nsd),            iBC(nshg),
     &            BC(nshg,ndofBC),         ilwork(nlwork),
     &            iper(nshg),              uold(nshg,nsd)
c
        dimension res(nshg,nflow),         
     &            rest(nshg),              solinc(nshg,ndof)
c     
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
        real*8   almit, alfit, gamit

        dimension ifath(numnp),    velbar(nfath,ndof),  nsons(nfath)
        real*8 rerr(nshg,10),ybar(nshg,ndof+8) ! 8 is for avg. of square as uu, vv, ww, pp, TT, uv, uw, and vw
        real*8, allocatable, dimension(:,:) :: vortG
        real*8, allocatable, dimension(:,:,:) :: BDiag

!       integer, allocatable, dimension(:) :: iv_rank  !used with MRasquin's version of probe points
!       integer :: iv_rankpernode, iv_totnodes, iv_totcores
!       integer :: iv_node,        iv_core,     iv_thread

! assuming three profiles to control (inlet, bottom FC and top FC)
! store velocity profile set via BC
        real*8 vbc_prof(nshg,3)
        character(len=60) fvarts
        integer ifuncs(6), iarray(10)
        integer BCdtKW, tsBase

        real*8 elDw(numel) ! element average of DES d variable

        real*8, allocatable, dimension(:,:) :: HBrg
        real*8, allocatable, dimension(:) :: eBrg
        real*8, allocatable, dimension(:) :: yBrg
        real*8, allocatable, dimension(:) :: Rcos, Rsin
c
c  Here are the data structures for sparse matrix GMRES
c
        integer, allocatable, dimension(:,:) :: rowp
        integer, allocatable, dimension(:) :: colm
        real*8, allocatable, dimension(:,:) :: lhsK
        real*8, allocatable, dimension(:,:) :: EGmass
        real*8, allocatable, dimension(:,:) :: EGmasst
 
        integer iTurbWall(nshg) 
        real*8 yInlet(3), yInletg(3)
        integer imapped, imapInlet(nshg)  !for now, used for setting Blower conditions
!        real*8 M_th, M_tc, M_tt
!        logical  exMc
!        real*8 vBC, vBCg
        real*8 vortmax, vortmaxg

       iprec=0 !PETSc - Disable PHASTA's BDiag. TODO: Preprocssor Switch

       call findTurbWall(iTurbWall)

!-------
! SETUP
!-------

       !HACK for debugging suction
!       call Write_Debug(myrank, 'wallNormal'//char(0), 
!     &                          'wnorm'//char(0), wnorm, 
!     &                          'd', nshg, 3, lstep)

       !Probe Point Setup
       call initProbePoints()
       if(exts) then  !exts is set in initProbePoints 
         write(fvarts, "('./varts/varts.', I0, '.dat')") lstep    
         fvarts = trim(fvarts)  

         if(myrank .eq. master) then 
           call TD_writeHeader(fvarts)
         endif
       endif
       
       !Mach Control Setup
       call MC_init(Delt, lstep, BC)
       exMC = exMC .and. exts   !If probe points aren't available, turn
                                !the Mach Control off
       if(exMC) then 
         call MC_applyBC(BC)
         call MC_printState()
       endif


c
c.... open history and aerodynamic forces files
c
        if (myrank .eq. master) then
          open (unit=ihist,  file=fhist,  status='unknown')
          open (unit=iforce, file=fforce, status='unknown')
        endif
c
c
c.... initialize
c
        ifuncs  = 0                      ! func. evaluation counter
        istep  = 0
        ntotGM = 0                      ! number of GMRES iterations
        time   = 0
        yold   = y
        acold  = ac

!Blower Setup
       call BC_init(Delt, lstep, BC)  !Note: sets BC_enable
! fix the yold values to the reset BC
      if(BC_enable) call itrBC (yold,  ac,  iBC,  BC,  iper, ilwork)
! without the above, second solve of first steps is fouled up
!
        allocate(HBrg(Kspace+1,Kspace))
        allocate(eBrg(Kspace+1))
        allocate(yBrg(Kspace))
        allocate(Rcos(Kspace))
        allocate(Rsin(Kspace))

        if (mod(impl(1),100)/10 .eq. 1) then
c
c     generate the sparse data fill vectors
c
           allocate  (rowp(nshg,nnz))
           allocate  (colm(nshg+1))
           call genadj(colm, rowp, icnt ) ! preprocess the adjacency list

           nnz_tot=icnt         ! this is exactly the number of non-zero 
                                ! blocks on this proc
           if(usingpetsc.eq.1) then
             allocate (BDiag(1,1,1))
           else
             allocate (BDiag(nshg,nflow,nflow))
             allocate (lhsK(nflow*nflow,nnz_tot))
           endif
        endif
        if (mod(impl(1),100)/10 .eq. 3) then
c
c     generate the ebe data fill vectors
c
           nedof=nflow*nshape
           allocate  (EGmass(numel,nedof*nedof))
        endif

c..........................................
        rerr = zero
        ybar(:,1:ndof) = y(:,1:ndof)
        do idof=1,5
           ybar(:,ndof+idof) = y(:,idof)*y(:,idof)
        enddo
        ybar(:,ndof+6) = y(:,1)*y(:,2)
        ybar(:,ndof+7) = y(:,1)*y(:,3)
        ybar(:,ndof+8) = y(:,2)*y(:,3)
c.........................................

!  change the freestream and inflow eddy viscosity to reflect different
!  levels of freestream turbulence
!
! First override boundary condition values
!  USES imapInlet from Mach Control so if that gets conditionaled away
!  it has to know if it is needed here
!
      if(isetEV_IC_BC.eq.1) then
        allocate(vortG(nshg, 4))
        call vortGLB(yold, x, shp, shgl, ilwork, vortG)
        vortmax=maxval(abs(vortG(:,4)))  !column 4 is the magnitude of the shear tensor - should actually probably be calld shearmax instead of vortmax
        write(*,*) "vortmax = ", vortmax
        
        !Find the maximum shear in the simulation
        if(numpe.gt.1) then
           call MPI_ALLREDUCE(vortmax, vortmaxg, 1, 
     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
           vortmax = vortmaxg
        endif

        !Apply eddy viscosity at the inlet
        do i=1,icountInlet ! for now only coded to catch primary inlet, not blower
           BC(imapInlet(i),7)=evis_IC_BC
        enddo
    
        !Apply eddy viscosity through the quasi-inviscid portion of the domain
        do i = 1,nshg
          if(abs(vortG(i,4)).ge.vortmax*0.01) yold(i,6)=evis_IC_BC
        enddo
        isclr=1 ! fix scalar
        call itrBCsclr(yold,ac,iBC,BC,iper,ilwork)
        
        deallocate(vortG)
      endif
c
c.... loop through the time sequences
c
        do 3000 itsq = 1, ntseq
        itseq = itsq
c
c.... set up the current parameters
c
        nstp   = nstep(itseq)
        nitr   = niter(itseq)
        LCtime = loctim(itseq)
c
        call itrSetup ( y,  acold)
        isclr=0

        niter(itseq)=0          ! count number of flow solves in a step
                                !  (# of iterations)
        do i=1,seqsize
           if(stepseq(i).eq.0) niter(itseq)=niter(itseq)+1
        enddo
        nitr = niter(itseq)
c
c.... determine how many scalar equations we are going to need to solve
c
        nsclrsol=nsclr          ! total number of scalars solved. At
                                ! some point we probably want to create
                                ! a map, considering stepseq(), to find
                                ! what is actually solved and only
                                ! dimension EGmasst to the appropriate
                                ! size.

        if(nsclrsol.gt.0)allocate  (EGmasst(numel*nshape*nshape
     &                              ,nsclrsol))
 
c
c.... loop through the time steps
c
        ttim(1) = REAL(secs(0.0)) / 100.
        ttim(2) = secs(0.0)

c        tcorecp1 = REAL(secs(0.0)) / 100.
c        tcorewc1 = secs(0.0)
        if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if(myrank.eq.master)  then
           tcorecp1 = TMRC()
        endif

        rmub=datmat(1,2,1)
        if(rmutarget.gt.0) then
           rmue=rmutarget
           xmulfact=(rmue/rmub)**(1.0/nstp)
           if(myrank.eq.master) then
              write(*,*) 'viscosity will by multiplied by ', xmulfact
              write(*,*) 'to bring it from ', rmub,' down to ', rmue
           endif
           datmat(1,2,1)=datmat(1,2,1)/xmulfact ! make first step right
        else
           rmue=datmat(1,2,1)   ! keep constant
           xmulfact=one
        endif
        if(iramp.eq.1) then
                call initBCprofileScale(vbc_prof,BC,yold,x)
! fix the yold values to the reset BC
                call itrBC (yold,  ac,  iBC,  BC,  iper, ilwork)
                isclr=1 ! fix scalar
                call itrBCsclr(yold,ac,iBC,BC,iper,ilwork)
        endif   
c  Time Varying BCs------------------------------------(Kyle W 6-6-13)
c        BCdtKW=0
        if(BCdtKW.gt.0) then
           call BCprofileInitKW(PresBase,VelBase,BC)
        endif
c  Time Varying BCs------------------------------------(Kyle W 6-6-13)
                                                        
867     continue


c============ Start the loop of time steps============================c
!        edamp2=.99
!        edamp3=0.05
        deltaInlInv=one/(0.125*0.0254)
        do 2000 istp = 1, nstp
           rerr=zero !extreme limit of 1 step window or error stats....later a variable

!           if (myrank.eq.master) write(*,*) 'Time step of current run', 
!     &                                    istp

c  Time Varying BCs------------------------------------(Kyle W 6-6-13)
           if(BCdtKW.gt.0) then
              call BCprofileScaleKW(PresBase,VelBase,BC,yold)
           endif
c  Time Varying BCs------------------------------------(Kyle W 6-6-13)

           if(iramp.eq.1) 
     &        call BCprofileScale(vbc_prof,BC,yold)

c           call rerun_check(stopjob)
cc          if(stopjob.ne.0) goto 2001
c
c Decay of scalars
c
           if(nsclr.gt.0 .and. tdecay.ne.1) then
              yold(:,6:ndof)=y(:,6:ndof)*tdecay
              BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*tdecay
           endif

           if(nosource.eq.1) BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*0.8


c           xi=istp*one/nstp
c           datmat(1,2,1)=rmub*(1.0-xi)+xi*rmue
           datmat(1,2,1)=xmulfact*datmat(1,2,1)

c.... if we have time varying boundary conditions update the values of BC.
c     these will be for time step n+1 so use lstep+1
c
           if(itvn.gt.0) call BCint((lstep+1)*Delt(1), shp, shgl,
     &                              shpb, shglb, x, BC, iBC)

            if(iLES.gt.0) then
c
c.... get dynamic model coefficient
c
            ilesmod=iLES/10  
c
c digit bit set filter rule, 10 bit set model
c
            if (ilesmod.eq.0) then ! 0 < iLES < 10 => dyn. model calculated
                                   ! at nodes based on discrete filtering
               call getdmc (yold,       shgl,      shp, 
     &                      iper,       ilwork,    nsons,
     &                      ifath,      x)
            endif
            if (ilesmod .eq. 1) then ! 10 < iLES < 20 => dynamic-mixed
                                     ! at nodes based on discrete filtering
               call bardmc (yold,       shgl,      shp, 
     &                      iper,       ilwork,    
     &                      nsons,      ifath,     x) 
            endif
            if (ilesmod .eq. 2) then ! 20 < iLES < 30 => dynamic at quad
                                     ! pts based on lumped projection filt. 
               call projdmc (yold,       shgl,      shp, 
     &                       iper,       ilwork,    x) 
            endif
c
           endif ! endif of iLES


c
c.... set traction BCs for modeled walls
c
            if (itwmod.ne.0) then   !wallfn check
               call asbwmod(yold,   acold,   x,      BC,     iBC,
     &                      iper,   ilwork,  ifath,  velbar)
            endif
c
c.... -----------------------> predictor phase <-----------------------
c
            call itrPredict(   yold,    acold,    y,   ac )
            call itrBC (y,  ac,  iBC,  BC,  iper, ilwork)
            isclr = zero
            if (nsclr.gt.zero) then
            do isclr=1,nsclr
               call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
            enddo
            endif
c
c.... --------------------> multi-corrector phase <--------------------
c
           iter=0
            ilss=0  ! this is a switch thrown on first solve of LS redistance
c
cHACK to make it keep solving RANS until tolerance is reached
c
       istop=0
! blocking extra RANS steps for now       iMoreRANS=0 
       iMoreRANS=6 
c 
c find the the RANS portion of the sequence
c
       do istepc=1,seqsize
          if(stepseq(istepc).eq.10) iLastRANS=istepc
       enddo

       iseqStart=1
9876   continue
c
            do istepc=iseqStart,seqsize
               icode=stepseq(istepc)
               if(mod(icode,10).eq.0) then ! this is a solve
                  isolve=icode/10
                  if(isolve.eq.0) then   ! flow solve (encoded as 0)
c
                     etol=epstol(1)
                     iter   = iter+1
                     ifuncs(1)  = ifuncs(1) + 1
c     
c.... reset the aerodynamic forces
c     
                     Force(1) = zero
                     Force(2) = zero
                     Force(3) = zero
                     HFlux    = zero
c     
c.... form the element data and solve the matrix problem
c     
c.... explicit solver
c     
                     if (impl(itseq) .eq. 0) then
                        if (myrank .eq. master)
     &                       call error('itrdrv  ','impl ',impl(itseq))
                     endif
                     if (mod(impl(1),100)/10 .eq. 1) then  ! sparse solve
c     
c.... preconditioned sparse matrix GMRES solver
c     
                        lhs = 1 - min(1,mod(ifuncs(1)-1,LHSupd(1))) 
                        iprec=lhs
                        nedof = nflow*nshape
c                        write(*,*) 'lhs=',lhs
                    if(usingpetsc.eq.1) then
#if (HAVE_PETSC)
               call SolGMRp (y,             ac,            yold,
     &                       x,
     &                       iBC,           BC,
     &                       colm,          rowp,          lhsk,
     &                       res,
     &                       BDiag,         
     &                       iper,          ilwork,
     &                       shp,           shgl,
     &                       shpb,          shglb,         solinc,
     &                       rerr,          fncorp )
#else
                     if(myrank.eq.0) write(*,*) 'exiting because run time input asked for PETSc, not linked in exec'
                     call error('itrdrv  ','noPETSc',usingpetsc)
#endif
                     else
                      call SolGMRs (y,             ac,            yold,
     &                       acold,         x,
     &                       iBC,           BC,
     &                       colm,          rowp,          lhsk,
     &                       res,
     &                       BDiag,         hBrg,          eBrg,
     &                       yBrg,          Rcos,          Rsin,
     &                       iper,          ilwork,
     &                       shp,           shgl,
     &                       shpb,          shglb,         solinc,
     &                       rerr)
                    endif
                      else if (mod(impl(1),100)/10 .eq. 2) then ! mfg solve
c     
c.... preconditioned matrix-free GMRES solver
c     
                        lhs=0
                        iprec = 1 - min(1,mod(ifuncs(1)-1,LHSupd(1))) 
                        nedof = 0
                        call SolMFG (y,             ac,            yold,
     &                       acold,         x,
     &                       iBC,           BC,
     &                       res,           
     &                       BDiag,         HBrg,          eBrg,
     &                       yBrg,          Rcos,          Rsin,
     &                       iper,          ilwork,
     &                       shp,           shgl,
     &                       shpb,          shglb,         solinc, 
     &                       rerr)
c     
                     else if (mod(impl(1),100)/10 .eq. 3) then ! ebe solve
c.... preconditioned ebe matrix GMRES solver
c     

                        lhs = 1 - min(1,mod(ifuncs(1)-1,LHSupd(1))) 
                        iprec = lhs
                        nedof = nflow*nshape
c                        write(*,*) 'lhs=',lhs
                      call SolGMRe (y,             ac,            yold,
     &                       acold,         x,
     &                       iBC,           BC,
     &                       EGmass,        res,
     &                       BDiag,         HBrg,          eBrg,
     &                       yBrg,          Rcos,          Rsin,
     &                       iper,          ilwork,
     &                       shp,           shgl,
     &                       shpb,          shglb,         solinc,
     &                       rerr)
                     endif
c     
                else          ! solve a scalar  (encoded at isclr*10)
                     isclr=isolve
                     etol=epstol(isclr+2) ! note that for both epstol and LHSupd 1 is flow 2 temp isclr+2 for scalars
                     ifuncs(isclr+2)  = ifuncs(isclr+2) + 1
                     if((iLSet.eq.2).and.(ilss.eq.0)
     &                    .and.(isclr.eq.2)) then 
                        ilss=1  ! throw switch (once per step)
c     
c... copy the first scalar at t=n+1 into the second scalar of the 
c... level sets
c     
                     y(:,6)    = yold(:,6)  + (y(:,6)-yold(:,6))/alfi
                     y(:,7)    =  y(:,6)
                     yold(:,7) = y(:,7)
                     ac(:,7)   = zero
c     
                     call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
c     
c....store the flow alpha, gamma parameter values and assigm them the 
c....Backward Euler parameters to solve the second levelset scalar
c     
                        alfit=alfi
                        gamit=gami
                        almit=almi
                        alfi = 1
                        gami = 1
                        almi = 1
                     endif
c     
                     lhs = 1 - min(1,mod(ifuncs(isclr+2)-1,
     &                                       LHSupd(isclr+2)))
                     iprec = lhs
                     istop=0
                 if(usingPETSc.eq.1) then
#if (HAVE_PETSC)
                     call SolGMRpSclr(y,             ac,  
     &                    x,             elDw,
     &                    iBC,           BC,          
     &                    colm,           rowp, 
     &                    iper,          ilwork,
     &                    shp,           shgl,
     &                    shpb,          shglb,     rest,
     &                    solinc(1,isclr+5),fncorp)
#else
                     write(*,*) 'exiting because run time input asked for PETSc, not linked in exec'
                     call error('itrdrv  ','noPETSc',usingpetsc)
#endif
                 else
                     call SolGMRSclr(y,             ac,         yold,
     &                    acold,         EGmasst(1,isclr),
     &                    x,             elDw,
     &                    iBC,           BC,          
     &                    rest,           
     &                    HBrg,          eBrg,
     &                    yBrg,          Rcos,      Rsin,
     &                    iper,          ilwork,
     &                    shp,           shgl,
     &                    shpb,          shglb, solinc(1,isclr+5))
                  endif
c     
                  endif         ! end of scalar type solve
c     
c     
c.... end of the multi-corrector loop
c     
 1000             continue      !check this

               else             ! this is an update  (mod did not equal zero)
                  iupdate=icode/10 ! what to update
                  if(iupdate.eq.0) then !update flow  
                     call itrCorrect ( y, ac, yold, acold, solinc)
                     call itrBC (y,  ac,  iBC,  BC, iper, ilwork)
                     call tnanq(y, 5, 'y_updbc')
c Elaine-SPEBC
                     if((irscale.ge.0).and.(myrank.eq.master)) then
                        call genscale(y, x, iBC)
c                       call itrBC (y,  ac,  iBC,  BC, iper, ilwork)
                     endif
                  else          ! update scalar
                     isclr=iupdate !unless
                     if(iupdate.eq.nsclr+1) isclr=0
                     call itrCorrectSclr ( y, ac, yold, acold,
     &                                     solinc(1,isclr+5))
                     if (ilset.eq.2 .and. isclr.eq.2)  then
                        fct2=one/almi
                        fct3=one/alfi
                        acold(:,7) = acold(:,7) 
     &                             + (ac(:,7)-acold(:,7))*fct2
                        yold(:,7)  = yold(:,7)  
     &                             + (y(:,7)-yold(:,7))*fct3  
                        call itrBCSclr (  yold,  acold,  iBC,  BC, 
     &                                    iper,  ilwork)
                        ac(:,7) = acold(:,7)*(one-almi/gami)
                        y(:,7)  = yold(:,7)
                        ac(:,7) = zero 
                        if (ivconstraint .eq. 1) then
     &                       
c ... applying the volume constraint
c
                           call solvecon (y,    x,      iBC,  BC, 
     &                                    iper, ilwork, shp,  shgl)
c
                        endif   ! end of volume constraint calculations
                     endif
                     call itrBCSclr (  y,  ac,  iBC,  BC, iper, ilwork)
                  endif
               endif            !end of switch between solve or update
            enddo               ! loop over sequence in step
        if((istop.lt.0).and.(iMoreRANS.lt.5)) then
            iMoreRANS=iMoreRANS+1
            if(myrank.eq.master) write(*,*) 'istop =', istop
       iseqStart=iLastRANS
       goto 9876
       endif
c     
c     Find the solution at the end of the timestep and move it to old
c  
c.... First to reassign the parameters for the original time integrator scheme
c
            if((iLSet.eq.2).and.(ilss.eq.1)) then 
               alfi =alfit
               gami =gamit
               almi =almit  
            endif          
            call itrUpdate( yold,  acold,   y,    ac)
            call itrBC (yold, acold,  iBC,  BC, iper,ilwork)  
c Elaine-SPEBC      
            if((irscale.ge.0).and.(myrank.eq.master)) then
                call genscale(yold, x, iBC)
c               call itrBC (y,  ac,  iBC,  BC, iper, ilwork)
            endif           
            do isclr=1,nsclr
               call itrBCSclr (yold, acold,  iBC, BC, iper, ilwork)
            enddo
c     
            istep = istep + 1
            lstep = lstep + 1
            ntoutv=max(ntout,100) 
            !boundary flux output moved after the error calculation so
            !everything can be written out in a single chunk of code -
            !Nicholas Mati
            
            !dump TIME SERIES
            if (exts) then
              !Write the probe data to disc. Note that IO to disc only
              !occurs when mod(lstep, nbuff) == 0. However, this
              !function also does data buffering and must be called
              !every time step. 
              call TD_bufferData()
              call TD_writeData(fvarts, .false.)
            endif
            
            !Update the Mach Control
            if(exts .and. exMC) then
              !Note: the function MC_updateState must be called after
              !the function TD_bufferData due to dependencies on
              !vartsbuff(:,:,:). 
              call MC_updateState()
              call MC_applyBC(BC)
              call MC_printState()
                 
              !Write the state if a restart is also being written. 
              if(mod(lstep,ntout).eq.0 ) call MC_writeState(lstep)
            endif

            !update blower control
            if(BC_enable) then
              !Update the blower boundary conditions for the next 
              !iteration. 
              call BC_iter(BC)

              !Also write the current phases of the blowers if a 
              !restart is also being written. 
              if(mod(lstep, ntout) == 0) call BC_writePhase(lstep)
            endif
            
            !.... Yi Chen Duct geometry8
            if(isetBlowing_Duct.gt.0)then
              if(ifixBlowingVel_Duct.eq.0)then
                if(nstp.gt.nBlowingStepsDuct)then
                  nBlowingStepsDuct = nstp-2
                endif
                call setBlowing_Duct2(x,BC,yold,iTurbWall,istp)
              endif
            endif
          !... Yi Chen Duct geometry8

c
c.... -------------------> error calculation  <-----------------
            if(ierrcalc.eq.1.or.ioybar.eq.1) then
               tfact=one/istep
               do idof=1,ndof
                 ybar(:,idof) =tfact*yold(:,idof) +
     &                         (one-tfact)*ybar(:,idof)
               enddo
c....compute average
c...  ybar(:,ndof+1:ndof+8) is for avg. of square as uu, vv, ww, pp, TT, uv, uw, and vw
               do idof=1,5 ! avg. of square for 5 flow variables
                   ybar(:,ndof+idof) = tfact*yold(:,idof)**2 +
     &                             (one-tfact)*ybar(:,ndof+idof)
               enddo
               ybar(:,ndof+6) = tfact*yold(:,1)*yold(:,2) + !uv
     &                          (one-tfact)*ybar(:,ndof+6)
               ybar(:,ndof+7) = tfact*yold(:,1)*yold(:,3) + !uw
     &                          (one-tfact)*ybar(:,ndof+7)
               ybar(:,ndof+8) = tfact*yold(:,2)*yold(:,3) + !vw
     &                          (one-tfact)*ybar(:,ndof+8)
c... compute err
c hack ShockError
c  
               errmax=maxval(rerr(:,6))
               errswitch=0.1*errmax  
!
! note this scalefactor will govern the thickness of the refinement band around the shock.  
! Higher values make it thinner (less refinement), lower -> thicker
! what is below is specific to SAM adapt's expectation to adapt when the 6th field is > 1.0e-6
! note also that this field was altered in e3.f and e3ls.f to no longer be the traditional error 
! indicator, rather it is based on element jump of Temperature to identify shocks
!
               where(rerr(:,6).gt.errswitch)
                    rerr(:,6)=1.0
               elsewhere
                    rerr(:,6)=1e-10
               endwhere
               rerr(:, 7)=rerr(:, 7)+(yold(:,1)-ybar(:,1))**2
               rerr(:, 8)=rerr(:, 8)+(yold(:,2)-ybar(:,2))**2
               rerr(:, 9)=rerr(:, 9)+(yold(:,3)-ybar(:,3))**2
               rerr(:,10)=rerr(:,10)+(yold(:,4)-ybar(:,4))**2
            endif
           
c.. writing ybar field if requested in each restart file		
            
            !here is where we save our averaged field.  In some cases we want to
            !write it less frequently		
            if( (irs >= 1) .and. (
     &        mod(lstep, ntout) == 0 .or. !Checkpoint
     &        istep == nstp) )then        !End of simulation
              if(output_mode .eq. -1 ) then ! this is an in-memory adapt case
                if(istep == nstp) then ! go ahead and take care of it
                  call checkpoint (nstp,yold, acold, ybar, rerr,  velbar, 
     &                       x, iper, ilwork, shp, shgl, iBC )
                endif
                if(ntout.le.lstep) then ! user also wants file output
                  output_mode=0
                  call checkpoint (nstp,yold, acold, ybar, rerr,  velbar, 
     &                       x, iper, ilwork, shp, shgl, iBC )
                  output_mode=-1 ! reset to stream 
                endif
              else
                call checkpoint (nstp,yold, acold, ybar, rerr,  velbar, 
     &                       x, iper, ilwork, shp, shgl, iBC )
              endif   
             endif   

 2000    continue  !end of NSTEP loop
 2001    continue  

         ttim(1) = REAL(secs(0.0)) / 100. - ttim(1)
         ttim(2) = secs(0.0)              - ttim(2)

c         tcorecp2 = REAL(secs(0.0)) / 100.
c         tcorewc2 = secs(0.0)
         if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         if(myrank.eq.master)  then
            tcorecp2 = TMRC()
            write(6,*) 'T(core) cpu = ',tcorecp2-tcorecp1
         endif
        
c     call wtime

      call destroyWallData
      call destroyfncorp

 3000 continue !end of NTSEQ loop
c     
c.... ---------------------->  Post Processing  <----------------------
c     
c.... print out the last step
c     
!      if ( (irs .ge. 1) .and. ((mod(lstep, ntout) .ne. 0) .or.
!     &    (nstp .eq. 0))) then
!          if( (mod(lstep, ntoutv) .eq. 0) .and.
!     &        ((irscale.ge.0).or.(itwmod.gt.0) .or. 
!     &        ((nsonmax.eq.1).and.(iLES.gt.0))))
!     &        call rwvelb  ('out ',  velbar  ,ifail)
!
!          call Bflux  (yold,  acold,     x,
!     &         shp,           shgl,      shpb,
!     &         shglb,         nodflx,    ilwork)
!      endif



c      if(ioybar.eq.1) then
c         call write_field(myrank,'a'//char(0),'ybar'//char(0),4,
c     &                      ybar,'d'//char(0),nshg,ndof+8,lstep)
c      endif

c     if(iRANS.lt.0 .and. idistcalc.eq.1) then
c        call write_field(myrank,'a'//char(0),'DESd'//char(0),4,
c     &                      elDw,'d'//char(0),numel,1,lstep)
c
c         call write_field(myrank,'a'//char(0),'dwal'//char(0),4,
c     &                    d2wall,'d'//char(0),nshg,1,lstep)
c     endif 

c
c.... close history and aerodynamic forces files
c     
      if (myrank .eq. master) then
         close (ihist)
         close (iforce)
             
         if(exMC) then 
           call MC_writeState(lstep) 
           call MC_finalize
         endif
             
         if(exts) then 
           call TD_writeData(fvarts, .true.)    !force the flush of the buffer. 
           call TD_finalize
         endif
      endif

      if(BC_enable) then  !blower is allocated on all processes. 
        if(mod(lstep, ntout) /= 0) then !May have already written file.
           call BC_writePhase(lstep)
        endif

        call BC_finalize
      endif

      close (iecho)
      if(iabc==1) deallocate(acs)
c
c.... end
c
      return
      end

      subroutine checkpoint (nstp,yold, acold, ybar, rerr,  velbar, 
     &                       x, iper, ilwork, shp, shgl, iBC )
c
      use turbSA
      include "common.h"
      dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            iper(nshg),              iBC(nshg),
     &            x(nshg,nsd),         ilwork(nlwork)
      real*8  velbar(nfath,ndof),
     &        yold(nshg,ndof),      acold(nshg,ndof),           
     &        rerr(nshg,10),        ybar(nshg,ndof+8) 
! 8 is for avg. of square as uu, vv, ww, pp, TT, uv, uw, and vw

      if( (mod(lstep, ntout) .eq. 0) .and.
     &              ((irscale.ge.0).or.(itwmod.gt.0) .or. 
     &              ((nsonmax.eq.1).and.(iLES.gt.0))))
     &              call rwvelb  ('out ',  velbar  ,ifail)

!BUG: need to update new_interface to work with SyncIO.
      !Bflux is presently completely crippled. Note that restar
      !has also been moved here for readability. 
!              call Bflux  (yold,          acold,     x,  compute boundary fluxes and print out
!    &              shp,           shgl,      shpb,
!    &              shglb,         nodflx,    ilwork)
                  
      call timer ('Output  ')      !set up the timer

      !write the solution and time derivative 
      call restar ('out ',  yold, acold)  

      !Write the distance to wall field in each restart
      if((istep==nstp) .and. (irans < 0 )) then !d2wall is allocated
                 call write_field(myrank,'a'//char(0),'dwal'//char(0),4,
     &                            d2wall,'d'//char(0), nshg, 1, lstep)
      endif 
           
      !Write the time average in each restart. 
      if(ioybar.eq.1)then
                 call write_field(myrank,'a'//char(0),'ybar'//char(0),4,
     &                              ybar,'d'//char(0),nshg,ndof+8,lstep)
      endif
                 
      !Write the error feild at the end of each step sequence
      if(ierrcalc.eq.1 .and. istep == nstp) then 
        !smooth the error indicators
      
        do i=1,ierrsmooth
         call errsmooth( rerr, x, iper, ilwork, shp, shgl, iBC )
        enddo
                   
        call write_field( myrank, 'a'//char(0), 'errors'//char(0), 6, 
     &                        rerr, 'd'//char(0), nshg, 10, lstep)
      endif

c the following is a future way to have the number of steps in the header...done for posix but not yet for syncio
c
c              call write_field2(myrank,'a'//char(0),'ybar'//char(0),
c     &                          4,ybar,'d'//char(0),nshg,ndof+8,
c     &                         lstep,istep)
c
c.... end
c
      return
      end

