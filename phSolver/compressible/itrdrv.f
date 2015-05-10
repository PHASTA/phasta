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
      use pvsQbi     !gives us splag (the spmass at the end of this run 
      use specialBC  !gives us itvn
      use timedata   !allows collection of time series
      use turbSA

        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
      
c
        dimension y(nshg,ndof),            ac(nshg,ndof),           
     &           yold(nshg,ndof),         acold(nshg,ndof),           
     &            x(numnp,nsd),            iBC(nshg),
     &            BC(nshg,ndofBC),         ilwork(nlwork),
     &            iper(nshg),              uold(nshg,nsd)
c
        dimension res(nshg,nflow),         BDiag(nshg,nflow,nflow),
     &            rest(nshg),              solinc(nshg,ndof)
c     
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
        real*8   almit, alfit, gamit
        dimension ifath(numnp),    velbar(nfath,ndof),  nsons(nfath)
        real*8 rerr(nshg,10),ybar(nshg,ndof+8) ! 8 is for avg. of square as uu, vv, ww, pp, TT, uv, uw, and vw
        integer, allocatable, dimension(:) :: ivarts
        integer, allocatable, dimension(:) :: ivartsg
        real*8, allocatable, dimension(:) :: vartssoln
        real*8, allocatable, dimension(:) :: vartssolng
        real*8, allocatable, dimension(:,:,:) :: vartsbuff
! assuming three profiles to control (inlet, bottom FC and top FC)
! store velocity profile set via BC
        real*8 vbc_prof(nshg,3)
        real*8 PresBase, VelBase
        character*20 fname1, fmt1, fname2, fmt2
        character*4 fname4c ! 4 characters
        character*60 fvarts
        character*5  cname
        character*10  cname2
        integer ifuncs(6), iarray(10)
        integer BCdtKW, tsBase

        real*8 elDw(numel) ! element average of DES d variable

c
c  Here are the data structures for sparse matrix GMRES
c
       integer, allocatable, dimension(:,:) :: rowp
       integer, allocatable, dimension(:) :: colm
       real*8, allocatable, dimension(:,:) :: lhsK
       real*8, allocatable, dimension(:,:) :: EGmass
       real*8, allocatable, dimension(:,:) :: EGmasst

       integer iTurbWall(nshg)
 
       call findTurbWall(iTurbWall)

       inquire(file='xyzts.dat',exist=exts)
       lskeep=lstep 
       if(exts) then
         
          open(unit=626,file='xyzts.dat',status='old')
          read(626,*) ntspts, freq, tolpt, iterat, varcod
          call sTD              ! sets data structures
          do jj=1,ntspts        ! read coordinate data where solution desired
             read(626,*) ptts(jj,1),ptts(jj,2),ptts(jj,3)
          enddo
          close(626)

           statptts(:,:) = 0
           parptts(:,:) = zero
           varts(:,:) = zero           

           allocate (ivarts(ntspts*ndof))
           allocate (ivartsg(ntspts*ndof))
           allocate (vartssoln(ntspts*ndof))
           allocate (vartssolng(ntspts*ndof))

           nbuff=ntout
           allocate (vartsbuff(ntspts,ndof,nbuff))

           if (myrank .eq. master) then
              do jj=1,ntspts
                 fvarts='varts/varts'
                 fvarts=trim(fvarts)//trim(cname2(jj))
                 fvarts=trim(fvarts)//trim(cname2(lstep))
                 fvarts=trim(fvarts)//'.dat'
                 fvarts=trim(fvarts)
                 open(unit=1000+jj, file=fvarts, status='unknown')
              enddo
           endif 
          
       endif
 
c
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
        if (mod(impl(1),100)/10 .eq. 1) then
c
c     generate the sparse data fill vectors
c
           allocate  (rowp(nshg,nnz))
           allocate  (colm(nshg+1))
           call genadj(colm, rowp, icnt ) ! preprocess the adjacency list

           nnz_tot=icnt         ! this is exactly the number of non-zero 
                                ! blocks on this proc
           allocate (lhsK(nflow*nflow,nnz_tot))
        endif
        if (mod(impl(1),100)/10 .eq. 3) then
c
c     generate the ebe data fill vectors
c
           nedof=nflow*nshape
           allocate  (EGmass(numel,nedof*nedof))
        endif
cc

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


        vbc_prof(:,1:3) = BC(:,3:5)

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
        if(myrank.eq.0)  then
           tcorecp1 = TMRC()
        endif

        rmub=datmat(1,2,1)
        if(rmutarget.gt.0) then
           rmue=rmutarget
           xmulfact=(rmue/rmub)**(1.0/nstp)
           if(myrank.eq.0) then
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
        do 2000 istp = 1, nstp

           if (myrank.eq.master) write(*,*) 'Time step of current run', 
     &                                    istp

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
           iMoreRANS=0 
c 
c find the the RANS portion of the sequence
c
           do istepc=1,seqsize
              if(stepseq(istepc).eq.10) iLastRANS=istepc
           enddo

           iseqStart= 1
9876       continue
c

c         open(unit=72,file='ien_solgmrs.dat',
c     &            status='unknown')

c=========== Start the loop of the non-linear iteration s=========c
           do istepc=iseqStart,seqsize

c                      if(myrank.eq.master)then
c                          write(*,*) 'icode, start sequence',icode
c                          write(*,*) mien(409)%p(:,1)
c                      endif

c                     do iblk = 407,409
c                      if( myrank.eq.master)then
c                        write(*,*) 'generate ien.dat'
c                        write(72,*)'step construction',istepc
c                        write(72,*)'ien',iblk
c                        write(72,*) mien(iblk)%p(:,1)
c                      endif
c                     enddo


               icode=stepseq(istepc)
c               if(myrank.eq.master)write(*,*) 'Stagger:',istepc,icode

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


c------------- Duct debug, output initial condition confined with boundary conditions --------
c         lstep_orig=lstep
c         lstep = 300+istepc
c         call restar('out ',y,ac)
c         lstep=lstep_orig
c         if (myrank.eq.master) then
c           open(unit=72,file='numstart.dat',status='old')
c           write(72,*) lstep
c           close(72)
c         endif
c         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c.................
c                     if(myrank.eq.master)then
c                       write(*,*)'icode,start solve', icode
c                       write(*,*)mien(409)%p(:,1)
c                     endif


                      call SolGMRs (y,             ac,            yold,
     &                       acold,         x,
     &                       iBC,           BC,
     &                       colm,          rowp,          lhsk,
     &                       res,
     &                       BDiag,         a(mHBrg),      a(meBrg),
     &                       a(myBrg),      a(mRcos),      a(mRsin),
     &                       iper,          ilwork,
     &                       shp,           shgl,
     &                       shpb,          shglb,         solinc,
     &                       rerr)

c                      if(myrank.eq.master)then
c                          write(*,*) 'icode, end solve',icode
c                          write(*,*) mien(409)%p(:,1)
c                      endif
 
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
     &                       BDiag,         a(mHBrg),      a(meBrg),
     &                       a(myBrg),      a(mRcos),      a(mRsin),
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
     &                       BDiag,         a(mHBrg),      a(meBrg),
     &                       a(myBrg),      a(mRcos),      a(mRsin),
     &                       iper,          ilwork,
     &                       shp,           shgl,
     &                       shpb,          shglb,         solinc,
     &                       rerr)
                     endif
c     
                  else          ! solve a scalar  (encoded at isclr*10)
                     ifuncs(isclr+2)  = ifuncs(isclr+2) + 1
                     etol=epstol(isclr+1)
                     isclr=isolve
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
                     call SolGMRSclr(y,             ac,         yold,
     &                    acold,         EGmasst(1,isclr),
     &                    x,             elDw,
     &                    iBC,           BC,          
     &                    rest,           
     &                    a(mHBrg),      a(meBrg),
     &                    a(myBrg),      a(mRcos),    a(mRsin),
     &                    iper,          ilwork,
     &                    shp,           shgl,
     &                    shpb,          shglb, solinc(1,isclr+5))
c     
                  endif         ! end of scalar type solve
c     
c     
c.... end of the multi-corrector loop
c     
 1000             continue      !check this

               else             ! this is an correct  (mod did not equal zero)
                  iupdate=icode/10 ! what to correct
                  if(iupdate.eq.0) then !correct flow

                     call itrCorrect ( y, ac, yold, acold, solinc)
c------------- NASA debug, output before itrBC --------
c                     if(myrank.eq.master)then
c                       write(*,*)'Creating debug restart file...'
c                     endif
c                     lstep_orig=lstep
c                     lstep = 400+lstep
c                     call restar('out ',y,ac)
c                     lstep=lstep_orig
c                     if (myrank.eq.master) then
c                       open(unit=72,file='numstart.dat',status='old')
c                       write(72,*) lstep
c                       close(72)
c                     endif
c                     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c.................

                     call itrBC (y,  ac,  iBC,  BC, iper, ilwork)
                     call tnanq(y, 5, 'y_updbc')
c Elaine-SPEBC
                     if((irscale.ge.0).and.(myrank.eq.master)) then
                        call genscale(y, x, iBC)
c                       call itrBC (y,  ac,  iBC,  BC, iper, ilwork)
                     endif

c                     if(myrank.eq.master)then
c                       write(*,*)'icode,end correct', icode
c                       write(*,*)mien(409)%p(:,1)
c                     endif

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

c            if(myrank.eq.master)then 
c              write(*,*)'icode, end loop over sequence', icode
c              write(*,*)mien(49)%p(:,1)
c            endif
 
            enddo               ! loop over sequence in step

            if((istop.lt.0).and.(iMoreRANS.lt.5)) then
                iMoreRANS=iMoreRANS+1
            if(myrank.eq.0) write(*,*) 'istop =', istop
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
c     
c.... compute boundary fluxes and print out
c     
            lstep = lstep + 1

            ntoutv=max(ntout,100) 
            if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
c
c here is where we save our averaged field.  In some cases we want to
c write it less frequently
               if( (mod(lstep, ntoutv) .eq. 0) .and.
     &              ((irscale.ge.0).or.(itwmod.gt.0) .or. 
     &              ((nsonmax.eq.1).and.(iLES.gt.0))))
     &              call rwvelb  ('out ',  velbar  ,ifail)

               call Bflux  (yold,          acold,     x,
     &              shp,           shgl,      shpb,
     &              shglb,         nodflx,    ilwork)

               call vortGLB(yold, x, shp, shgl, ilwork)
            endif
c...  dump TIME SERIES
            
            if (exts) then
               if (mod(lstep-1,freq).eq.0) then
                  
                  if (numpe > 1) then
                     do jj = 1, ntspts
                        vartssoln((jj-1)*ndof+1:jj*ndof)=varts(jj,:)
                        ivarts=zero
                     enddo
                     do k=1,ndof*ntspts
                        if(vartssoln(k).ne.zero) ivarts(k)=1
                     enddo
                     call MPI_REDUCE(vartssoln, vartssolng, ndof*ntspts,
     &                    MPI_DOUBLE_PRECISION, MPI_SUM, master,
     &                    MPI_COMM_WORLD, ierr)

                     call MPI_REDUCE(ivarts, ivartsg, ndof*ntspts,
     &                    MPI_INTEGER, MPI_SUM, master,
     &                    MPI_COMM_WORLD, ierr)

                     if (myrank.eq.zero) then
                        do jj = 1, ntspts

                           indxvarts = (jj-1)*ndof
                           do k=1,ndof
                              if(ivartsg(indxvarts+k).ne.0) then ! none of the vartssoln(parts) were non zero
                                 varts(jj,k)=vartssolng(indxvarts+k)/
     &                                ivartsg(indxvarts+k)
                              endif
                           enddo
                        enddo
                     endif !only on master
                  endif !only if numpe > 1
        
                 if( istp.eq. nstp) then !make sure incomplete buffers get purged.
                    icheck=mod(nstp,nbuff)
                    if(icheck.ne.0) nbuff=icheck
                 endif

                  if (myrank.eq.zero) then
                     k=mod(lstep,nbuff)
                     if(k.eq.0) k=nbuff
                     do jj = 1, ntspts
                       vartsbuff(jj,1:5,k)=varts(jj,1:5)
                     enddo
                     if(k.eq. nbuff) then
                       do jj = 1, ntspts
                        ifile = 1000+jj
                        do ibuf=1,nbuff
                        write(ifile,555) lstep-1 -nbuff+ibuf, 
     &                  (vartsbuff(jj,k,ibuf), k=1,5) ! assuming ndof to be 5
                        enddo ! buff empty
        
c                        call flush(ifile)
                       enddo ! jj ntspts
                        endif !only dump when buffer full
                  endif !only on master

cc.... Yi Chen Duct geometry8
         if(isetBlowing_Duct.gt.0)then
           if(ifixBlowingVel_Duct.eq.0)then
             if(nstp.gt.nBlowingStepsDuct)then
               nBlowingStepsDuct = nstp-2
             endif
             call setBlowing_Duct2(x,BC,yold,iTurbWall,istp)
           endif
         endif
         if(isetSuction_Duct.gt.0)then
            call setSuction_Duct2(x,BC,yold)
         endif
cc... Yi Chen Duct geometry8


                  varts(:,:) = zero ! reset the array for next step

 555              format(i6,5(2x,E12.5e2))

               endif
            endif
            
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
            endif
           
c.. writing ybar field if requested in each restart file
            if(ioybar.eq.1 .and. mod(lstep,ntout).eq.0 )then
              call write_field2(myrank,'a','ybar',4,
     &             ybar,'d',nshg,ndof+8,lstep,istep)
            endif   

c.... end of the NSTEP and NTSEQ loops...................................
 2000    continue
 2001    continue

         ttim(1) = REAL(secs(0.0)) / 100. - ttim(1)
         ttim(2) = secs(0.0)                 - ttim(2)

c         tcorecp2 = REAL(secs(0.0)) / 100.
c         tcorewc2 = secs(0.0)
         if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
         if(myrank.eq.0)  then
            tcorecp2 = TMRC()
            write(6,*) 'T(core) cpu = ',tcorecp2-tcorecp1
         endif
        
c     call wtime

 3000 continue
c     
c.... ---------------------->  Post Processing  <----------------------
c     
c.... print out the last step
c     
      if ((irs .ge. 1) .and. ((mod(lstep, ntout) .ne. 0) .or.
     &     (nstp .eq. 0))) then
         if( (mod(lstep, ntoutv) .eq. 0) .and.
     &        ((irscale.ge.0).or.(itwmod.gt.0) .or. 
     &        ((nsonmax.eq.1).and.(iLES.gt.0))))
     &        call rwvelb  ('out ',  velbar  ,ifail)

         call Bflux  (yold,  acold,     x,
     &        shp,           shgl,      shpb,
     &        shglb,         nodflx,    ilwork)

      endif
c     
      if(ierrcalc.eq.1) then
c
c.....smooth the error indicators
c
c ! errsmooth is currently available only for incompressible code
c        do i=1,ierrsmooth
c            call errsmooth( rerr, x, iper, ilwork, shp, shgl, iBC )
c        end do

         call write_error(myrank, lstep, nshg, 10, rerr )
      endif


      if(iRANS.lt.0 .and. idistcalc.eq.1) then
         call write_field(myrank,'a','DESd',4,
     &                    elDw,'d',numel,1,lstep)

         call write_field(myrank,'a','dwal',4,
     &                    d2wall,'d',nshg,1,lstep)
      endif 

c
c.... close history and aerodynamic forces files
c     
      if (myrank .eq. master) then
         close (ihist)
         close (iforce)
         if(exts) then
            deallocate(ivarts)
            deallocate(ivartsg)
            deallocate(vartssoln)
            deallocate(vartssolng)
            do jj=1,ntspts
               close(1000+jj)
            enddo
         endif
      endif
      close (iecho)
      if(iabc==1) deallocate(acs)
c
c.... end
c
        return
        end
