        subroutine ElmGMR (u,         y,         ac,        x,     
     &                     shp,       shgl,      iBC,
     &                     BC,        shpb,      shglb,
     &                     res,       iper,      ilwork,
     &                     rowp,      colm,     lhsK,      
     &                     lhsP,      rerr,     GradV)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c Alberto Figueroa, Winter 2004.  CMM-FSI
c Irene Vignon, Spring 2004.
c----------------------------------------------------------------------
c
        use pvsQbi  ! brings in NABI
        use stats   !  
        use pointer_data  ! brings in the pointers for the blocked arrays
        use local_mass
        use timedata
c
        include "common.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            u(nshg,nsd),
     &            x(numnp,nsd),               
     &            iBC(nshg),           
     &            BC(nshg,ndofBC),  
     &            res(nshg,nflow),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg,idflx),     rmass(nshg)
        dimension GradV(nshg,nsdsq)
c
        dimension ilwork(nlwork)

        integer rowp(nshg*nnz),         colm(nshg+1)

        real*8 lhsK(9,nnz_tot), lhsP(4,nnz_tot)

        real*8, allocatable, dimension(:,:,:,:) :: xKebe, xGoC

        real*8  rerr(nshg,10)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

        real*8 spmasstot(20),  ebres(nshg)
c
c.... set up the timer
c

CAD        call timer ('Elm_Form')
c
c.... -------------------->   diffusive flux   <--------------------
c
c.... set up parameters
c
        ires   = 1

        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction
                                                       ! of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
           qres = zero
           if(iter == nitr .and. icomputevort == 1 ) then
             !write(*,*) 'iter:',iter,' - nitr:',nitr
             !write(*,*) 'Setting GradV to zero'
             GradV = zero
           endif
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)
c     
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

            if(iter == nitr .and. icomputevort == 1 ) then
               !write(*,*) 'Calling AsIqGradV'
               call AsIqGradV (y,                x,
     &                   shp(lcsyst,1:nshl,:), 
     &                   shgl(lcsyst,:,1:nshl,:),
     &                   mien(iblk)%p,
     &                   GradV)
             endif
              call AsIq (y,                x,                       
     &                   shp(lcsyst,1:nshl,:), 
     &                   shgl(lcsyst,:,1:nshl,:),
     &                   mien(iblk)%p,     mxmudmi(iblk)%p,  
     &                   qres,             rmass )
           enddo
       
c
c.... form the diffusive flux approximation
c
           call qpbc( rmass, qres, iBC, iper, ilwork )       
           if(iter == nitr .and. icomputevort == 1 ) then
             !write(*,*) 'Calling solveGradV'
             call solveGradV( rmass, GradV, iBC, iper, ilwork )
           endif
c
        endif 
c
c.... -------------------->   interior elements   <--------------------
c
        res    = zero
        if (stsResFlg .ne. 1) then
           flxID = zero
        endif

        if (lhs .eq. 1) then
           lhsp   = zero
           lhsk   = zero
        endif
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
          iblock = iblk         ! used in local mass inverse (p>2)
          iblkts = iblk         ! used in timeseries
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
c.... allocate the element matrices
c
          allocate ( xKebe(npro,9,nshl,nshl) )
          allocate ( xGoC (npro,4,nshl,nshl) )
c
c.... compute and assemble the residual and tangent matrix
c
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIGMR (y,                   ac,
     &                 x,                   mxmudmi(iblk)%p,      
     &                 tmpshp, 
     &                 tmpshgl,
     &                 mien(iblk)%p,
     &                 res,
     &                 qres,                xKebe,
     &                 xGoC,                rerr)
c
c.... satisfy the BC's on the implicit LHS
c     
          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             if(ipord.eq.1) 
     &         call bc3lhs (iBC, BC,mien(iblk)%p, xKebe)  
             call fillsparseI (mien(iblk)%p, 
     &                 xKebe,            lhsK,
     &                 xGoC,             lhsP,
     &                 rowp,                      colm)
          endif

          deallocate ( xKebe )
          deallocate ( xGoC  )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c.... end of interior element loop
c
       enddo
c$$$       if(ibksiz.eq.20 .and. iwrote.ne.789) then
c$$$          do i=1,nshg
c$$$             write(789,*) 'eqn block ',i 
c$$$             do j=colm(i),colm(i+1)-1
c$$$                write(789,*) 'var block',rowp(j)
c$$$
c$$$                do ii=1,3
c$$$                   write(789,111) (lhsK((ii-1)*3+jj,j),jj=1,3)
c$$$                enddo
c$$$             enddo
c$$$          enddo
c$$$          close(789)
c$$$          iwrote=789
c$$$       endif
c$$$ 111   format(3(e14.7,2x))
c$$$c
c.... add in lumped mass contributions if needed
c
       if((flmpr.ne.0).or.(flmpl.ne.0)) then
          call lmassadd(ac,res,rowp,colm,lhsK,gmass)
       endif

       have_local_mass = 1
c
c.... time average statistics
c       
       if ( stsResFlg .eq. 1 ) then

          if (numpe > 1) then
             call commu (stsVec, ilwork, nResDims  , 'in ')
          endif
          do j = 1,nshg
             if (btest(iBC(j),10)) then
                i = iper(j)
                stsVec(i,:) = stsVec(i,:) + stsVec(j,:)
             endif
          enddo
c     
          do i = 1,nshg
             stsVec(i,:) = stsVec(iper(i),:)
          enddo

          if (numpe > 1) then
             call commu (stsVec, ilwork, nResDims  , 'out')
          endif
          return
          
       endif
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
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel 


          if(lcsyst.eq.3) lcsyst=nenbl
c
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
c
c.... allocate the element matrices
c
          allocate ( xKebe(npro,9,nshl,nshl) )
          allocate ( xGoC (npro,4,nshl,nshl) )
          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBMFG (u,                       y,
     &                 ac,                      x,
     &                 tmpshpb,
     &                 tmpshglb,
     &                 mienb(iblk)%p,           mmatb(iblk)%p,
     &                 miBCB(iblk)%p,           mBCB(iblk)%p,
     &                 res,                     xKebe)

c
c.... satisfy (again, for the vessel wall contributions) the BC's on the implicit LHS
c
c.... first, we need to make xGoC zero, since it doesn't have contributions from the 
c.... vessel wall elements

          xGoC = zero

          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             if(ipord.eq.1)
     &         call bc3lhs (iBC, BC,mienb(iblk)%p, xKebe)
             call fillsparseI (mienb(iblk)%p,
     &                 xKebe,           lhsK,
     &                 xGoC,             lhsP,
     &                 rowp,                      colm)
          endif

          deallocate ( xKebe )
          deallocate ( xGoC )
          deallocate (tmpshpb)
          deallocate (tmpshglb)
c
c.... end of boundary element loop
c
       enddo

       if(ipvsq.ge.1) then
c
c....  pressure vs. resistance boundary condition sets pressure at
c      outflow to linearly increase as flow through that face increases
c      (routine is at bottom of this file)
c
          call ElmpvsQ (res,y,-1.0d0)     
       endif
           
c
c before the commu we need to rotate the residual vector for axisymmetric
c boundary conditions (so that off processor periodicity is a dof add instead
c of a dof combination).  Take care of all nodes now so periodicity, like
c commu is a simple dof add.
c
       if(iabc==1)              !are there any axisym bc's
     &       call rotabc(res, iBC,  'in ')
c
c
c.... -------------------->   communications <-------------------------
c

       if (numpe > 1) then
          call commu (res  , ilwork, nflow  , 'in ')
       endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3Res (iBC,  BC,  res,  iper, ilwork)
c
c.... return
c
c      call timer ('Back    ')
      return
      end


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!********************************************************************
!--------------------------------------------------------------------

      subroutine ElmGMRSclr (y,         ac,        x,     
     &                       shp,       shgl,      iBC,
     &                       BC,        shpb,      shglb,
     &                       res,       iper,      ilwork,
     &                       rowp,      colm,      lhsS    )
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c----------------------------------------------------------------------
c
        use pointer_data
        use local_mass
c
        include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),         iBC(nshg),           
     &            BC(nshg,ndofBC),      res(nshg),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg,nsd),     rmass(nshg)
c
        integer ilwork(nlwork), rowp(nshg*nnz),   colm(nshg+1)

        real*8 lhsS(nnz_tot)

        real*8, allocatable, dimension(:,:,:) :: xSebe
c
c.... set up the timer
c

CAD        call timer ('Elm_Form')
c
c.... -------------------->   diffusive flux   <--------------------
c
        ires   = 1

        if (idiff==1 .or. idiff==3) then ! global reconstruction of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
           qres = zero
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lcsyst = lcblk(3,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              
              ngauss = nint(lcsyst)
c     
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

              call AsIqSclr (y,                   x,                       
     &                       shp(lcsyst,1:nshl,:), 
     &                       shgl(lcsyst,:,1:nshl,:),
     &                       mien(iblk)%p,     qres,                   
     &                       rmass )
       
           enddo
       
c
c.... form the diffusive flux approximation
c
           call qpbcSclr ( rmass, qres, iBC, iper, ilwork )       
c
        endif 
c
c.... -------------------->   interior elements   <--------------------
c
        res    = zero
        spmass = zero

        if (lhs .eq. 1) then
           lhsS   = zero
        endif

        if ((impl(1)/10) .eq. 0) then   ! no flow solve so flxID was not zeroed
           flxID = zero
        endif
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
          iblock = iblk         ! used in local mass inverse (p>2)
          iel    = lcblk(1,iblk)
          lcsyst = lcblk(3,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          ndofl  = lcblk(8,iblk)
          npro   = lcblk(1,iblk+1) - iel

          ngauss = nint(lcsyst)
c
c.... allocate the element matrices
c
          allocate ( xSebe(npro,nshl,nshl) )
c
c.... compute and assemble the residual and tangent matrix
c
          call AsIGMRSclr(y,                   ac,
     &                 x,
     &                 shp(lcsyst,1:nshl,:), 
     &                 shgl(lcsyst,:,1:nshl,:),
     &                 mien(iblk)%p,        res,
     &                 qres,                xSebe, mxmudmi(iblk)%p )
c
c.... satisfy the BC's on the implicit LHS
c     
          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             call fillsparseSclr (mien(iblk)%p, 
     &                 xSebe,             lhsS,
     &                 rowp,              colm)
          endif

          deallocate ( xSebe )
c
c.... end of interior element loop
c
       enddo

c
c.... add in lumped mass contributions if needed
c
       if((flmpr.ne.0).or.(flmpl.ne.0)) then
          call lmassaddSclr(ac(:,isclr), res,rowp,colm,lhsS,gmass)
       endif

       have_local_mass = 1
c
c
c  call DtN routine which updates the flux to be consistent with the
c  current solution values.  We will put the result in the last slot of
c  BC (we added a space in input.f).  That way we can localize this
c  value to the boundary elements.  This is important to keep from calling
c  the DtN evaluator more than once per node (it can be very expensive).
c
         if(idtn.eq.1)  call DtN(iBC,BC,y)
c
c.... -------------------->   boundary elements   <--------------------
c
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lcsyst = lcblkb(3,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel

          if(lcsyst.eq.3) lcsyst=nenbl
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
c
c localize the dtn boundary condition
c

          if(idtn.eq.1)   call dtnl(   iBC, BC, mienb(iblk)%p,
     &              miBCB(iblk)%p,  mBCB(iblk)%p)

c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          call AsBSclr (y,                       x,
     &                  shpb(lcsyst,1:nshl,:),
     &                  shglb(lcsyst,:,1:nshl,:),
     &                  mienb(iblk)%p,           mmatb(iblk)%p,
     &                  miBCB(iblk)%p,           mBCB(iblk)%p,
     &                  res)
c
c.... end of boundary element loop
c
        enddo
c
c
c.... -------------------->   communications <-------------------------
c

      if (numpe > 1) then
        call commu (res  , ilwork, 1  , 'in ')
      endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResSclr (iBC,  res,  iper, ilwork)
c
c.... return
c
CAD      call timer ('Back    ')
      return
      end

        
c
c....routine to compute and return the flow rates for coupled surfaces of a given type
c        
      subroutine GetFlowQ (qsurf,y,srfIdList,numSrfs)
        
      use pvsQbi  ! brings in NABI
c
      include "common.h"
      include "mpif.h"

      real*8  y(nshg,3)
      real*8  qsurf(0:MAXSURF),qsurfProc(0:MAXSURF)
      integer numSrfs, irankCoupled, srfIdList(0:MAXSURF)

c note we only need the first three entries (u) from y

      qsurfProc=zero
      
      do i = 1,nshg
      
        if(numSrfs.gt.zero) then
          do k = 1,numSrfs
            irankCoupled = 0
            if (srfIdList(k).eq.ndsurf(i)) then
              irankCoupled=k
              do j = 1,3              
                 qsurfProc(irankCoupled) = qsurfProc(irankCoupled)
     &                            + NABI(i,j)*y(i,j)
              enddo
              exit
            endif      
          enddo       
        endif
      
      enddo
c      
c     at this point, each qsurf has its "nodes" contributions to Q
c     accumulated into qsurf. Note, because NABI is on processor this
c     will NOT be Q for the surface yet
c
c.... reduce integrated Q for each surface, push on qsurf
c
       npars=MAXSURF+1
       if(impistat.eq.1) then
         iAllR = iAllR+1
       elseif(impistat.eq.2) then
          iAllR = iAllR+1
       endif
       if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
       if(impistat.gt.0) rmpitmr = TMRC()
       call MPI_ALLREDUCE (qsurfProc, qsurf(:), npars,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
       if(impistat.eq.1) then 
         rAllR = rAllR+TMRC()-rmpitmr
       elseif(impistat.eq.2) then
         rAllRScal = rAllRScal+TMRC()-rmpitmr
       endif
  
c
c.... return
c
      return
      end


        
c
c... routine to couple pressure with flow rate for each coupled surface
c
      subroutine ElmpvsQ (res,y,sign)     

      use pvsQbi  ! brings in NABI
      use convolImpFlow !brings in the current part of convol coef for imp BC
      use convolRCRFlow !brings in the current part of convol coef for RCR BC

      include "common.h"
      include "mpif.h"

      real*8 res(nshg,ndof),y(nshg,3)
      real*8 p(0:MAXSURF)
      integer irankCoupled

c
c... get p for the resistance BC
c           
      if(numResistSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistResist,numResistSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        p(:)=sign*p(:)*ValueListResist(:) ! p=QR  now we have the true pressure on each
                                        ! outflow surface.  Note sign is -1
                                        ! for RHS, +1 for LHS
c
c....  multiply it by integral NA n_i
c     
       do i = 1,nshg
          do k = 1,numResistSrfs
              irankCoupled = 0
              if (nsrflistResist(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)     
                  exit 
              endif
          enddo   
       enddo
       
      endif !end of coupling for Resistance BC

      
c
c... get p for the impedance BC
c     
      if(numImpSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistImp,numImpSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        do j = 1,numImpSrfs
            if(sign.lt.zero) then ! RHS so -1
               p(j)= sign*(poldImp(j) + p(j)*ImpConvCoef(ntimeptpT+2,j))  !pressure p=pold+ Qbeta
            elseif(sign.gt.zero) then ! LHS so sign is positive
                p(j)= sign*p(j)*ImpConvCoef(ntimeptpT+2,j)
            endif
        enddo
             
c
c....  multiply it by integral NA n_i
c     
       do i = 1,nshg
          do k = 1,numImpSrfs
              irankCoupled = 0
              if (nsrflistImp(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)      
                  exit
              endif
          enddo   
       enddo
       
      endif !end of coupling for Impedance BC
c
c... get p for the RCR BC
c     
      if(numRCRSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistRCR,numRCRSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        do j = 1,numRCRSrfs
            if(sign.lt.zero) then ! RHS so -1
                p(j)= sign*(poldRCR(j) + p(j)*RCRConvCoef(lstep+2,j)) !pressure p=pold+ Qbet
                p(j)= p(j) - HopRCR(j) ! H operator contribution 
            elseif(sign.gt.zero) then ! LHS so sign is positive
                p(j)= sign*p(j)*RCRConvCoef(lstep+2,j)
            endif
        enddo
             
c
c....  multiply it by integral NA n_i
c     
       do i = 1,nshg
          do k = 1,numRCRSrfs
              irankCoupled = 0
              if (nsrflistRCR(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)
                  exit
              endif
          enddo   
       enddo
       
      endif !end of coupling for RCR BC

      return
      end







