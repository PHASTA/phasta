c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc       SPARSE   
c_______________________________________________________________

        subroutine ElmGMRPETSc (y,         ac,        x,         
     &                     shp,       shgl,      iBC,
     &                     BC,        shpb,      shglb,
     &                     res,       rmes,      
     &                     iper,      ilwork,   
     &                     rerr, lhsP) 
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c----------------------------------------------------------------------
c
        use pointer_data
        use timedataC
c
        include "common.h"
        include "mpif.h"
c
!        integer col(nshg+1), row(nnz*nshg)
!        real*8 lhsK(nflow*nflow,nnz_tot)
        real*8 BDiag(1,1,1)
        
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),               
     &            iBC(nshg),
     &            BC(nshg,ndofBC),      
     &            res(nshg,nflow),
     &            rmes(nshg,nflow),      
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg, idflx),     rmass(nshg)
c
        dimension ilwork(nlwork)

        real*8 Bdiagvec(nshg,nflow), rerr(nshg,10)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)
        real*8, allocatable :: EGmass(:,:,:)
        integer gnode(numnp)
	ttim(80) = ttim(80) - secs(0.0)
        iprec=0
c
c.... set up the timer
c
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!        if(myrank.eq.0) write (*,*) 'top of elmgmrpetsc '

         call timer ('Elm_Form')
c
c.... -------------------->   interior elements   <--------------------
c
c.... set up parameters
c
        ires   = 1
c
        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction
                                                       ! of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
        qres = zero
        rmass = zero
        
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          ngauss = nint(lcsyst)
c
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIq (y,                x,                       
     &               tmpshp,              
     &               tmpshgl,
     &               mien(iblk)%p,     mxmudmi(iblk)%p,
     &               qres,                   
     &               rmass)

          deallocate ( tmpshp )
          deallocate ( tmpshgl ) 
       enddo
       
c
c.... take care of periodic boundary conditions
c

       call qpbc( rmass, qres, iBC, iper, ilwork )       
c
      endif                     ! computation of global diffusive flux
c
c.... loop over element blocks to compute element residuals
c
c
c.... initialize the arrays
c
        res    = zero
        rmes   = zero ! to avoid trap_uninitialized
        !if (lhs. eq. 1)   lhsK = zero
        flxID = zero
c
c.... loop over the element-blocks
c
!        call commu (res  , ilwork, nflow  , 'in ') !FOR TEST
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!        if(myrank.eq.0) write (*,*) 'after res zeroed '
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          iblkts = iblk          ! used in timeseries
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
c

          if(lhs.eq.1) then
             allocate (EGmass(npro,nedof,nedof))
             EGmass = zero
          else
             allocate (EGmass(1,1,1))
          endif

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIGMR (y,                   ac,
     &                 x,                   mxmudmi(iblk)%p,
     &                 tmpshp,
     &                 tmpshgl,             mien(iblk)%p,
     &                 mmat(iblk)%p,        res,
     &                 rmes,                BDiag,
     &                 qres,                EGmass,
     &                 rerr )
          if(lhs.eq.1) then
c
c.... satisfy the BCs on the implicit LHS
c     
             call bc3LHS (iBC,                  BC,  mien(iblk)%p, 
     &                    EGmass  ) 

c
c.... Fill-up the global sparse LHS mass matrix
c
!        if(myrank.eq.0) write (*,*) 'before fillsparsepetscc ',iblk
             call cycle_count_start()
             call fillsparsecpetscc( mieng(iblk)%p, EGmass, lhsP)
             call cycle_count_stop()
          endif
c
          deallocate ( EGmass )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c.... end of interior element loop
c
       enddo
       if(lhs.eq.1)   call cycle_count_print()
!        call commu (res  , ilwork, nflow  , 'in ') !FOR TEST
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!        if(myrank.eq.0) write (*,*) 'after fillsparsepetscc '
c
c.... -------------------->   boundary elements   <--------------------
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          if(lcsyst.eq.3) lcsyst=nenbl
          ngaussb = nintb(lcsyst)
          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c

          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBMFG (y,                       x,
     &                 tmpshpb,                 tmpshglb, 
     &                 mienb(iblk)%p,           mmatb(iblk)%p,
     &                 miBCB(iblk)%p,           mBCB(iblk)%p,
     &                 res,                     rmes)

          deallocate (tmpshpb)
          deallocate (tmpshglb)
c
c.... end of boundary element loop
c
        enddo
c
      ttim(80) = ttim(80) + secs(0.0)
c
c before the commu we need to rotate the residual vector for axisymmetric
c boundary conditions (so that off processor periodicity is a dof add instead
c of a dof combination).  Take care of all nodes now so periodicity, like
c commu is a simple dof add.
c
      if(iabc==1) then               !are there any axisym bc's
          call rotabc(res(1,2), iBC,  'in ')
       endif

c.... -------------------->   communications <-------------------------
c

!      if (numpe > 1) then
      if (numpe < 1) then
        call commu (res  , ilwork, nflow  , 'in ')

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3Res (y,  iBC,  BC,  res,  iper, ilwork)
c
c
c.... return
c
!      if (numpe > 1) then
      if (numpe < 1) then
        call commu (res  , ilwork, nflow  , 'out')
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif
      call timer ('Back    ')
      return
      end
c
c

c

        subroutine ElmGMRPETScSclr(y,      ac,
     &                        x,      elDw,         
     &                        shp,    shgl,   iBC,
     &                        BC,     shpb,   shglb,
     &                        rest,   rmest,  
     &                        iper,   ilwork, lhsPs)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c----------------------------------------------------------------------
c
        use pointer_data
c
        include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),              
     &            iBC(nshg),
     &            BC(nshg,ndofBC),      
     &            rest(nshg),           Diag(1),
     &            rmest(nshg),          
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qrest(nshg),          rmasst(nshg)
        real*8 elDw(numel)
c
        dimension ilwork(nlwork)
c
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)
        real*8, allocatable :: EGMasst(:,:,:)
        real*8, allocatable :: ElDwl(:)
c
	ttim(80) = ttim(80) - tmr()
        iprec=0
c
c.... set up the timer
c

        call timer ('Elm_Form')
c
c.... -------------------->   interior elements   <--------------------
c
c.... set up parameters
c
        intrul = intg  (1,itseq)
        intind = intpt (intrul)
c
        ires   = 1

c
c.... initialize the arrays
c
        rest    = zero
        rmest   = zero ! to avoid trap_uninitialized
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... compute and assemble the residual and tangent matrix
c
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
          if(lhs.eq.1) then
             allocate (EGMasst(npro,nshl,nshl))
             EGMasst=zero
          endif

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
          
          allocate (elDwl(npro))
c
          call AsIGMRSclr(y,                   
     &                    ac,
     &                    x,               elDwl,   
     &                    tmpshp,          tmpshgl,
     &                    mien(iblk)%p,
     &                    mmat(iblk)%p,    rest,
     &                    rmest,               
     &                    qrest,           EGmasst,
     &                    Diag )
c
           elDw(iel:inum) = elDwl(1:npro)
          deallocate ( elDwl )

           if(lhs.eq.1) then 
c.... satisfy the BCs on the implicit LHS
c     
          call bc3LHSSclr (iBC, mien(iblk)%p, EGmasst )
c
c
c.... Fill-up the global sparse LHS mass matrix
c
             call fillsparsecpetscs( mieng(iblk)%p, EGmasst, lhsPs)
          endif
          if(lhs.eq.1) deallocate ( EGmasst )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c.... end of interior element loop
c
       enddo
c
c.... -------------------->   boundary elements   <--------------------
c
c.... set up parameters
c
        intrul = intg   (2,itseq)
        intind = intptb (intrul)
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          if(lcsyst.eq.3) lcsyst=nenbl
          ngaussb = nintb(lcsyst)          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c

          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)
c
          call AsBMFGSclr (y,                  x,
     &                     tmpshpb,
     &                     tmpshglb, 
     &                     mienb(iblk)%p,      mmatb(iblk)%p,
     &                     miBCB(iblk)%p,      mBCB(iblk)%p,
     &                     rest,               rmest)
c
          deallocate ( tmpshpb )
          deallocate ( tmpshglb )

c.... end of boundary element loop
c
        enddo


      ttim(80) = ttim(80) + tmr()
c
c.... -------------------->   communications <-------------------------
c

      if (numpe < 1) then
        call commu (rest  , ilwork, 1  , 'in ')

        call MPI_BARRIER (MPI_COMM_WORLD,ierr)

      endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResSclr (y,  iBC,  BC,  rest,  iper, ilwork)
      if (numpe < 1) then
        call commu (rest  , ilwork, 1  , 'out')
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif
c
c.... return
c
      call timer ('Back    ')
      return
      end
