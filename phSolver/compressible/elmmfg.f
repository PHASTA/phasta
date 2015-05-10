        subroutine ElmMFG (y,         ac,        x,         
     &                     shp,       shgl,      iBC,
     &                     BC,        shpb,      shglb,
     &                     res,       rmes,      BDiag,
     &                     iper,      ilwork,    rerr)
c
c----------------------------------------------------------------------
c
c This routine calculates the preconditioning matrix and the
c R.H.S. residual vector for the Matrix-free Iterative solver.
c
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
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
     &            res(nshg,nflow),
     &            rmes(nshg,nflow),      BDiag(nshg,nflow,nflow),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)  
c   
        dimension qres(nshg,idflx),     rmass(nshg),
     &            tmp(nshg)
c
        dimension ilwork(nlwork)

        real*8  rerr(nshg,10)
	real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
	real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

	ttim(80) = ttim(80) - secs(0.0)
c
c.... set up the timer
c

        call timer ('Elm_Form')
c
c.... -------------------->   interior elements   <--------------------
c
c.... set up parameters
c
        ires   = 3

        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction
                                                       ! of qdiff
c     
c     loop over element blocks for the global reconstruction
c     of the diffusive flux vector, q, and lumped mass matrix, rmass
c     
           qres = zero
           rmass = zero
           
           do iblk = 1, nelblk
c     
c.... set up the parameters
c     
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
     &             tmpshp,
     &             tmpshgl,
     &             mien(iblk)%p,      mxmudmi(iblk)%p,
     &             qres,                   
     &             rmass)

	      deallocate (tmpshp)
	      deallocate (tmpshgl)
              
           enddo
c
c.... take care of periodic boundary conditions
c

       call qpbc( rmass, qres, iBC,  iper, ilwork ) 

c     
        endif                   ! computation of global diffusive flux
c     
c.... loop over element blocks to compute element residuals
c     
c     
c.... initialize the arrays
c     
        res  = zero
        rmes = zero
c
        if (iprec .ne. 0) BDiag = zero
c     
c.... loop over the element-blocks
c     
        do iblk = 1, nelblk
c     
c.... set up the parameters
c
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
c.... compute and assemble the residuals and the preconditioner
c     
	   allocate (tmpshp(nshl,MAXQPT))
	   allocate (tmpshgl(nsd,nshl,MAXQPT))

	   tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
	   tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

           call AsIMFG (y,                       ac,
     &                  x,                       mxmudmi(iblk)%p,
     &                  tmpshp,                  tmpshgl,
     &                  mien(iblk)%p,
     &                  mmat(iblk)%p,            res,
     &                  rmes,                    BDiag,
     &                  qres,                    rerr )

	   deallocate (tmpshp)
	   deallocate (tmpshgl)

c     
c.... end of interior element loop
c     
        enddo
c     
c.... -------------------->   boundary elements   <--------------------
c     
c.... loop over the elements
c     
        do iblk = 1, nelblb
c     
c.... set up the parameters
c     
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
     &          tmpshpb,                 tmpshglb, 
     &          mienb(iblk)%p,           mmatb(iblk)%p,
     &          miBCB(iblk)%p,           mBCB(iblk)%p,
     &          res,                     rmes)
  
	   deallocate (tmpshpb)
	   deallocate (tmpshglb)
c     
c.... end of boundary element loop
c     
        enddo
c     
        ttim(80) = ttim(80) + secs(0.0)
c     
c.... -------------------->   communications <-------------------------
c     
        if((iabc==1)) then      !are there any axisym bc's
           call rotabc(res(1,2), iBC, 'in ')     
           call rotabc(rmes(1,2), iBC, 'in ')
        endif

        if (numpe > 1) then
c     

            call commu (res  , ilwork, nflow  , 'in ')

            call MPI_BARRIER (MPI_COMM_WORLD,ierr)
            
            call commu (rmes , ilwork, nflow  , 'in ')
            if(iprec.ne.0) call commu(BDiag,ilwork, nflow*nflow, 'in ')
        endif

c     
c.... ---------------------->   post processing  <----------------------
c     
c.... satisfy the BCs on the residual and the modified residual
c     
        call bc3Res (y,  iBC,  BC,  res,  iper, ilwork)
        call bc3Res (y,  iBC,  BC,  rmes, iper, ilwork)

c     
c.... satisfy the BCs on the block-diagonal preconditioner
c     
        if (iprec .ne. 0) then
           call bc3BDg (y,  iBC,  BC,  BDiag, iper, ilwork)
        endif


c     
c.... return
c     
        call timer ('Back    ')
        return
        end


