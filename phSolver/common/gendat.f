        subroutine gendat (y,       ac,       x,      iBC,     BC,
     &                     iper,    ilwork,
     &                     shp,     shgl,    shpb,    shglb,
     &                     ifath,   velbar,   nsons ) 
c
c----------------------------------------------------------------------
c
c This routine inputs the geometry and the boundary conditions.
c
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      
        use readarrays          ! used to acess nBC
        use dtnmod
        use pointer_data
        use wallData            ! give access to wnorm, findWallNOrm
        include "common.h"
        include "mpif.h"

c
c arrays in the following line are now dimensioned in readnblk
c        dimension nBC(nshg)
c
        dimension y(nshg,ndof),      ac(nshg,ndof),
     &            x(numnp,nsd),      iBC(nshg),
     &            BC(nshg,ndofBC),
     &            nodflx(numflx),    ilwork(nlwork),
     &            iper(nshg)
c
c.... shape function declarations
c     
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
c  stuff for dynamic model s.w.avg and wall model
c
        dimension ifath(numnp),    velbar(nfath,nflow), nsons(nfath) 

! Hack to get suction right on part boundaries
!       dimension BC3(numnp,5)

c... Duct
        integer iTurbWall(nshg)
        integer nTopNodes
        integer TopSurface(nshg)

c...

c.... start the timer
c
        
CAD        call timer ('PrProces')

c
c  put a barrier here so that all of the files are done reading
c  This SHOULD shield any mpi_profile information from io bottlenecks
c
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

c... updated, truly useful things for Duct.........
         if(isetInitial_Duct.gt.0)then
            call setInitial_Duct(x)    ! in INIprofile.f
         endif

c
c.... ---------------------------->  Nodes  <--------------------------
c
c.... compute length scales
c
        call xyzbound(x)
c
c.... echo the coordinates
c
        if ((necho .lt. 2).and.(myrank.eq.master)) then
          do n = 1, numnp
            if (mod(n,50) .eq. 1) write (iecho,1000) ititle,(i,i=1,nsd)
            write (iecho,1100) n, (x(n,i),i=1,nsd)
          enddo
        endif
c
c.... prepare periodic boundary conditions
c
        do i = 1,nshg
          if (iper(i) .ne. 0) then
            nshg0 = nshg0 - 1
          else
            iper(i) = i
          endif
        enddo
c
c.... ---------------------->  Interior Elements  <--------------------
c
        ibound = 0
c
c.... generate the interior nodal mapping
c
        call genshp ( shp, shgl, nshape, nelblk)
c
c.... --------------------->  Boundary Conditions  <-------------------
c
c.... read and generate the boundary condition codes (iBC array)
c
        call geniBC (iBC)
c
c.... read and generate the essential boundary conditions (BC array)
c
        call genBC  (iBC,   BC,   point2x,
     &               point2ilwork, point2iper)
        deallocate(nBC)

c=========================================================================================
c Yi Chen 
c Duct geometry8
c... updated, truly useful things.........


c===== specify outlet pressure
         if(isetOutletID .gt. 0)then
           call SetUniOutPres(BC)     ! in BCprofile2.f
         endif

c==== specify inlet boundary conditions
c         if(isetInlet_Duct.gt.0)then
c           call setInlet_Duct(x,BC,iTurbWall) ! in BCprofile2.f
c         endif

c==== specify blowing conditions 
         if(isetBlowing_Duct.gt.0)then
            if (ifixBlowingVel_Duct.eq.0)then
              call setBlowing_Duct(BC,iTurbWall)
            else
              call setBlowing_Duct3(x,BC,iTurbWall) ! in setBlowing_Duct3.f, fixed jet inlet velocity
            endif
         endif

c====== specify wall conditions 
         call findTurbWall(iTurbWall)

c==== apply suction patch on sides
c suction is applied at the end so it will overwrite the velocity at any nodes shared by the no-slip walls
         call findWallNorm(x,iBC,ilwork,iper)
         if(isetSuctionID_Duct.gt.0)then
            call setSuction_Duct3(x, BC, y, ilwork)
         endif

c==== 
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c.... Yi Chen
c==========================================================================================
c
c.... ---------------------->  Boundary Elements  <--------------------
c
        ibound = 1
        call gtnods
c
c  We now take care of Direchlet to Neumann BC's.  It had to move here
c  so that the IBC array was of size nshg and ready to be marked.
c

        if(nsclr.gt.0) then 
           call initDtN         ! Dirichlet to Neumann module: 
                                     ! initialize this once only
           do iblk = 1, nelblb  ! number of blocks
              iel    = lcblkb(1,iblk)
              npro   = lcblkb(1,iblk+1) - iel
c
c  for the DtN BC we need to mark all of the nodes that are involved.
c
              do i=1,npro
c
c if this element has the BCB AND it has not been found yet then mark it
c
                 if(miBCB(iblk)%p(i,2).lt.0) then  
                    idtn = 1    !set the flag for dtn bc's
                    do j=1,nshapeb
                       do isclr=1,nsclr
                          ignd=mienb(iblk)%p(i,j)
                             ifeature(ignd) = abs(miBCB(iblk)%p(i,2))       
                             iBC(ignd)=ior(iBC(ignd),2**13)
                                ! must mark this as a Neumann BC now
                             miBCB(iblk)%p(i,1)=
     &                       ior(miBCB(iblk)%p(i,1),2**(4+isclr))
                       end do
                    end do
                 endif
              end do
           end do
        endif
c
c.... generate the boundary element shape functions
c
        call genshpb ( shpb, shglb, nshapeb, nelblb)
c.... Evaluate the shape funcs. and their gradients at the desired quadrature
c.... for filtering. Save these evaluations using a module
c
c KEJ moved them to this point because cdelsq now passed with module
c     and it is read in with velb.<stepnum>.<proc#> now
c
        if (iLES .gt. 0) then

           call setfilt         ! For setting quad. rule to use for integrating
           call filtprep        ! the hat filter.
           if(iLES/10 .eq. 2) then
              call setave       ! For averaging cdelsq computed at quad pts
              call aveprep(shp,x)
           endif
        endif
c
c User sets request pzero in solver.inp now
c
c        call genpzero(iBC,iper)
c
      if((myrank.eq.master).and.(irscale.ge.0)) then
         call setSPEBC(numnp,nsd) 
         call eqn_plane(point2x, iBC)
      endif
c
c.... --------------------->  Initial Conditions  <--------------------
c
c.... generate the initial conditions and initialize time varying BC
c
        call genini (iBC,      BC,         y, 
     &               ac,       iper, 
     &               ilwork,   ifath,      velbar,  
     &               nsons,    x,
     &               shp,     shgl,    shpb,    shglb) 
c
c.... close the geometry, boundary condition and material files
c
!MR CHANGE
!        close (igeom)
!MR CHANGE END
        close (ibndc)
        if (mexist) close (imat)
c
c.... return
c
CAD        call timer ('Back    ')
        return
c
c.... end of file error handling
c
999     call error ('gendat  ','end file',igeom)
c
1000    format(a80,//,
     &  ' N o d a l   C o o r d i n a t e s                  ',//,
     &  '    Node     ',12x,3('x',i1,:,17x))
1100    format(1p,2x,i5,13x,3(1e12.5,7x))
2000    format(a80,//,
     &  ' B o u n d a r y   F l u x   N o d e s              '//,
     &  '   index          Node          ')
2100    format(1x,i5,5x,i10)
c
        end


        subroutine xyzbound(x)

        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        dimension x(numnp,3)

        real*8   Forout(3), Forin(3)

        xlngth=maxval(x(:,1))
        ylngth=maxval(x(:,2))
        zlngth=maxval(x(:,3))
        if(numpe. gt. 1) then
           Forin=(/xlngth,ylngth,zlngth/)
           call MPI_ALLREDUCE (Forin, Forout, 3,
     &       MPI_DOUBLE_PRECISION,MPI_MAX, MPI_COMM_WORLD,ierr)
           xmax = Forout(1)
           ymax = Forout(2)
           zmax = Forout(3)
        else
           xmax = xlngth
           ymax = ylngth
           zmax = zlngth
        endif
        xlngth=minval(x(:,1))
        ylngth=minval(x(:,2))
        zlngth=minval(x(:,3))
        if(numpe .gt. 1) then
           Forin=(/xlngth,ylngth,zlngth/)
           call MPI_ALLREDUCE (Forin, Forout, 3,
     &       MPI_DOUBLE_PRECISION,MPI_MIN, MPI_COMM_WORLD,ierr)
        else
           Forout(1) = xlngth
           Forout(2) = ylngth
           Forout(3) = zlngth
        endif

        xlngth = xmax-Forout(1)
        ylngth = ymax-Forout(2)
        zlngth = zmax-Forout(3)

        if(myrank.eq.master) then
           print 108,  xlngth,ylngth,zlngth
        endif
 108    format(' Domain size (x,y,z):',2x,3f15.10)
        return
        end
