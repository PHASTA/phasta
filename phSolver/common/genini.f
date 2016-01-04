        subroutine genini (iBC, BC, y, ac, iper, ilwork,
     &               ifath, velbar,  nsons,x,
     &               shp,     shgl,    shpb,    shglb ) 
c
c----------------------------------------------------------------------
c This routine reads the initial values in primitive form (density,
c velocity and temperature), satisfies the boundary conditions and 
c converts them to Y-variables.
c
c input:
c  iBC    (nshg)               : boundary condition code
c  BC     (nshg,ndofBC)        : boundary condition constrain data
c   x     (numnp,nsd)	       : locations of nodes, numnp-> # of node
c				  nsd-> space dimension, 1=x, 2=y, 3=z 
C			
c
c output:
c  y      (nshg,ndof)          : initial values of Y variables
c
c
c Farzin Shakib, Winter 1986.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use specialBC   ! gets itvn from here
        use convolImpFlow ! brings in ntimeptpT and other variables
        use convolRCRFlow ! brings RCR variables
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension iBC(nshg),                iper(nshg),
     &            BC(nshg,ndofBC),          y(nshg,ndof),
     &            ac(nshg,ndof),            x(numnp,nsd),
     &            shp(MAXTOP,maxsh,MAXQPT),
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)

c
        dimension ilwork(nlwork)
c
c  stuff for dynamic model s.w.avg and wall model
c
        dimension ifath(numnp),    velbar(nfath,nflow),
     &           nsons(nfath) 

        character*20 fname1
        character*10 cname2
        character*5 cname
c
c.... -------------------------->  Restart  <---------------------------
c
c.... read q from [RESTAR.INP], reset LSTEP
c
        call restar ('in  ',  y,  ac)
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(matflg(1,1).eq.0)then ! compressible code
          call INIprofile(BC,y,x)
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if((itwmod.gt.0) 
     &     .or. (nsonmax.eq.1 .and. iLES.gt.0) ) then 
           call rwvelb('in  ',velbar,ifail)
c
c if the read failed calculate velbar
c
           if(ifail.eq.1) then
              call getvel (y,     ilwork, iBC,
     &                     nsons, ifath, velbar)
           endif
 
        endif   ! for the itwmod or irscale
c
c.... time varying boundary conditions as set from file bct.dat and impt.dat 
c     (see function for format in file bctint.f)
c
        if (itvn .gt. 0 ) then !for inlet velocities
           call initBCt( x, iBC, BC)
           call BCint(lstep*Delt(1),shp,shgl,shpb,shglb,x, BC, iBC)
        endif
        if (impfile .gt. 0 ) then !for impedance BC
           if(myrank.eq.master) then 
              write(*,*) 'reading Qhistor.dat'
           endif
           open(unit=816, file='Qhistor.dat',status='old')
           read (816,*) ntimeptpT
           allocate (QHistImp(ntimeptpT+1,numImpSrfs)) 
           allocate (QHistTry(ntimeptpT,numImpSrfs)) 
           allocate (QHistTryF(ntimeptpT,numImpSrfs)) 
           do j=1,ntimeptpT+1
              read(816,*) (QHistImp(j,n),n=1,numImpSrfs) !read flow history
           enddo
           close(816)              
           call initImpt() !read impedance data and initialize begin/end values
           do i=2,ntimeptpT
              call Impint((ntimeptpT-i+1)*Delt(1),i) !return Imp values in reverse order ZN->Z0
           enddo

           allocate (poldImp(0:MAXSURF)) !for pressure part that depends on the history only
        endif
        if (ircrfile .gt. 0 ) then !for RCR BC
           call initRCRt()      !read RCR data
           dtRCR(:) = Delt(1)/(ValueListRCR(2,:)*ValueListRCR(3,:))
           allocate (RCRConvCoef(nstep(1)+lstep+2,numRCRSrfs)) !for convolution coeff
           !last one needed only to have array of same size as imp BC
           allocate (QHistRCR(nstep(1)+lstep+1,numRCRSrfs)) !for flow history
           QHistRCR = zero
           if(lstep.ne.0) then
              fname1 = 'Qhistor'
              fname1 = trim(fname1)//trim(cname2(lstep))//'.dat'
              fname1 = trim(fname1)
              if(myrank.eq.master) then
                 write(*,*) 'reading ', fname1
              endif
              open(unit=816, file=fname1, status='old')
              read(816,*) k
              do j=1,lstep+1
                 read(816,*) (QHistRCR(j,n),n=1,numRCRSrfs) !read flow history
              enddo
              close(816) 
           endif
           allocate (poldRCR(0:MAXSURF)) !for pressure part that depends on the history only
           allocate (HopRCR(0:MAXSURF)) !for H operator contribution
        endif


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  this subroutine initBCprofileScale is called to generate the correct
c  BC array, including the siuation of Take BC from IC for inlet.
c  so this subroutine must be called before BCs are applied to ICs
c  as those following lines do
c        call INIBCprofile(BC,y,x) 
c        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c.... satisfy the boundary conditions
c
        call itrBC (y, ac,  iBC, BC, iper, ilwork)

        itempr=mod(impl(1),2)  ! tempr solve if impl odd
        if(itempr.eq.1) then
           isclr=0
           call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
        endif
        do isclr=1,nsclr
           call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
        enddo

        if((irscale.ge.0) .and. (myrank.eq.master)) then
           call genscale(y, x, iBC)
        endif
c
c.... --------------------------->  Echo  <----------------------------
c
c.... echo the initial data
c
        if ((necho .lt. 0).and.(myrank.eq.master)) then
          do n = 1, nshg
            if (mod(n,50) .eq. 1) write(iecho,1000) ititle,(i,i=1,ndof)
            write (iecho,1100) n, (y(n,i),i=1,ndof)
          enddo
        endif
c
c.... return
c
        return
c
1000    format(a80,//,
     &  ' I n i t i a l   V a l u e s                        ',//,
     &  '    Node ',/,
     &  '   Number ',6x,6('dof',i1,:,10x))
1100    format(1p,2x,i5,5x,5(e12.5,2x))
c
        end
