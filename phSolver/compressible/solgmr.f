        subroutine SolGMRe (y,         ac,        yold,      acold,
     &			   x,         iBC,       BC,        EGmass,    
     &                     res,       BDiag,     HBrg,      eBrg,
     &                     yBrg,      Rcos,      Rsin,      iper,
     &                     ilwork,    shp,       shgl,      shpb,
     &                     shglb,     Dy, rerr)
c
c----------------------------------------------------------------------
c
c  This is the preconditioned GMRES driver routine.
c
c input:
c  y      (nshg,ndof)           : Y-variables at n+alpha_v
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  yold   (nshg,ndof)           : Y-variables at beginning of step
c  acold  (nshg,ndof)           : Primvar. accel. variable at begng step
c  x      (numnp,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  HBrg   (Kspace+1,Kspace)      : Hessenberg matrix (LHS matrix)
c  eBrg   (Kspace+1)             : RHS      of Hessenberg minim. problem
c  yBrg   (Kspace)               : solution of Hessenberg minim. problem
c  Rcos   (Kspace)               : Rotational cosine of QR algorithm
c  Rsin   (Kspace)               : Rotational sine   of QR algorithm
c  shp(b) (nen,maxsh,melCat)     : element shape functions (boundary)
c  shgl(b)(nsd,nen,maxsh,melCat) : local gradients of shape functions
c
c output:
c  res    (nshg,nflow)           : preconditioned residual
c  BDiag  (nshg,nflow,nflow)      : block-diagonal preconditioner
c
c  
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use pointer_data
        
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension y(nshg,ndof),             ac(nshg,ndof),
     &            yold(nshg,ndof),          acold(nshg,ndof),
     &            x(numnp,nsd),
     &            iBC(nshg),                BC(nshg,ndofBC),
     &            res(nshg,nflow),
     &            BDiag(nshg,nflow,nflow),
     &            HBrg(Kspace+1,*),         eBrg(*),
     &            yBrg(*),                  Rcos(*),
     &            Rsin(*),                  ilwork(nlwork),
     &            iper(nshg),               EGmass(numel,nedof,nedof)!,
ctoomuchmem     &            Binv(numel,nedof,nedof)
c
        dimension Dy(nshg,nflow),            rmes(nshg,nflow),
     &            temp(nshg,nflow),
     &            uBrg(nshg,nflow,Kspace+1)
c        
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
      real*8    rerr(nshg,10)
c     
c.... *******************>> Element Data Formation <<******************
c
c
c.... set the parameters for flux and surface tension calculations
c
c
        idflx = 0 
        if(idiff >= 1)  idflx= idflx + (nflow-1) * nsd
        if (isurf == 1) idflx=idflx + nsd
c
c.... form the LHS matrices, the residual vector, and the block
c     diagonal preconditioner
c
        call ElmGMRe(y,             ac,            x,
     &               shp,           shgl,          iBC,
     &               BC,            shpb,
     &               shglb,         res,
     &               rmes,          BDiag,         iper,      
     &               ilwork,        EGmass,        rerr )
      rmes=res  ! saving the b vector (residual)
c
c.... **********************>>    EBE - GMRES    <<********************
c
        call timer ('Solver  ')
c
c.... ------------------------> Initialization <-----------------------
c
c
c.... LU decompose the block diagonals
c
        if (iprec .ne. 0)
     &  call i3LU (BDiag, res,  'LU_Fact ')
        
c
c.... block diagonal precondition residual
c
        call i3LU (BDiag, res,  'forward ')
c
c.... initialize Dy
c
        Dy = zero
c
c.... Pre-precondition the LHS mass matrix and set up the element
c     by element preconditioners
c
ctoomuchmemory note that Binv is demoted from huge array to just one
c        real*8 in i3pre because it takes too much memory

        call i3pre (BDiag,    Binv,   EGmass,  ilwork)
c     
c.... left EBE precondition the residual
c     
ctoomuchmem        call i3Pcond (Binv,  res,  ilwork,   'L_Pcond ')
c     
c.... copy res in uBrg(1)
c     
        uBrg(:,:,1) = res
c     
c.... calculate norm of residual
c
        temp  = res**2

        call sumgat (temp, nflow, summed, ilwork)
        unorm = sqrt(summed)
c
c.... check if GMRES iterations are required
c
        iKs    = 0
        lGMRES = 0
c
c.... if we are down to machine precision, don't bother solving
c
        if (unorm .lt. 100.*epsM**2) goto 3000 
c
c.... set up tolerance of the Hessenberg's problem
c
        epsnrm = etol * unorm
c
c.... ------------------------>  GMRES Loop  <-------------------------
c
c.... loop through GMRES cycles
c
        do 2000 mGMRES = 1, nGMRES
        lGMRES = mGMRES - 1
c
        if (lGMRES .gt. 0) then
c
c.... if GMRES restarts are necessary, calculate  R - A x
c
c
c.... right precondition Dy
c
           temp = Dy
           
ctoomuchmem           call i3Pcond (Binv,  temp,  ilwork,  'R_Pcond ')
c
c.... perform the A x product
c
           call Au1GMR (EGmass,  temp,  ilwork, iBC,iper)
c
c.... periodic nodes have to assemble results to their partners
c
           call bc3per (iBC,  temp,  iper, ilwork, nflow)
c
c.... left preconditioning
c
ctoomuchmem           call i3Pcond (Binv,  temp,  ilwork,  'L_Pcond ')
c
c.... subtract A x from residual and calculate the norm
c           
           temp = res - temp
           uBrg(:,:,1) = temp
c
c.... calculate the norm
c
           temp  = temp**2
           call sumgat (temp, nflow, summed, ilwork)
           unorm = sqrt(summed)
c     
c.... flop count
c     
      !      flops = flops + nflow*nshg+nshg
c     
        endif
c
c.... set up RHS of the Hessenberg's problem
c
        call clear (eBrg, Kspace+1)
        eBrg(1) = unorm
c
c.... normalize the first Krylov vector
c
        uBrg(:,:,1) = uBrg(:,:,1) / unorm
c
c.... loop through GMRES iterations
c
        do 1000 iK = 1, Kspace
           iKs = iK

           uBrg(:,:,iKs+1) = uBrg(:,:,iKs)
c
c.... right EBE precondition the LHS ( u_{i+1} <-- inverse(U) u_i )
c
ctoomuchmem           call i3Pcond (Binv,  uBrg(:,:,iKs+1), ilwork,  'R_Pcond ')
c
c.... Au product  ( u_{i+1} <-- EGmass u_{i+1} )
c
           call Au1GMR ( EGmass, uBrg(:,:,iKs+1),  ilwork, iBC,iper)
c
c.... periodic nodes have to assemble results to their partners
c
           call bc3per (iBC,  uBrg(:,:,iKs+1),  iper, ilwork, nflow)

c
c.... left EBE precondition the LHS ( u_{i+1} <-- inverse(L) u_{i+1} )
c
ctoomuchmem           call i3Pcond (Binv,  uBrg(:,:,iKs+1), ilwork, 'L_Pcond ')
c
c.... orthogonalize and get the norm
c
          do jK = 1, iKs+1  
c
            if (jK .eq. 1) then
c
              temp = uBrg(:,:,iKs+1) * uBrg(:,:,1)  ! {u_{i+1}*u_1} vector 
              call sumgat (temp, nflow, beta, ilwork) ! sum vector=(u_{i+1},u_1)
c
            else
c
c project off jK-1 vector
c
              uBrg(:,:,iKs+1) = uBrg(:,:,iKs+1) - beta * uBrg(:,:,jK-1)
c
              temp = uBrg(:,:,iKs+1) * uBrg(:,:,jK) !{u_{i+1}*u_j} vector
              call sumgat (temp, nflow, beta, ilwork) ! sum vector=(u_{i+1},u_j)
c
            endif
c
            HBrg(jK,iKs) = beta   ! put this in the Hessenberg Matrix
c
        enddo
c
   !      flops = flops + (3*iKs+1)*nflow*numnp+(iKs+1)*numnp
c
c  the last inner product was with what was left of the vector (after
c  projecting off all of the previous vectors
c
        unorm           = sqrt(beta)
        HBrg(iKs+1,iKs) = unorm   ! this fills the 1 sub diagonal band
c
c.... normalize the Krylov vector
c
        uBrg(:,:,iKs+1) = uBrg(:,:,iKs+1) / unorm  ! normalize the next Krylov
c vector
c
c.... construct and reduce the Hessenberg Matrix
c  since there is only one subdiagonal we can use a Givens rotation to 
c  rotate off each subdiagonal AS IT IS FORMED.   We do this because it
c  allows us to check progress of solution and quit when satisfied.  Note
c  that all future K vects will put a subdiagonal in the next column so
c  there is no penalty to work ahead as  the rotation for the next vector
c  will be unaffected by this rotation.
        
c     
c     H Y = E ========>   R_i H Y = R_i E
c     
           do jK = 1, iKs-1
              tmp            =  Rcos(jK) * HBrg(jK,  iKs) +
     &                          Rsin(jK) * HBrg(jK+1,iKs)
              HBrg(jK+1,iKs) = -Rsin(jK) * HBrg(jK,  iKs) +
     &                          Rcos(jK) * HBrg(jK+1,iKs)
              HBrg(jK,  iKs) =  tmp
           enddo
c     
           tmp             = sqrt(HBrg(iKs,iKs)**2 + HBrg(iKs+1,iKs)**2)
           Rcos(iKs)       = HBrg(iKs,  iKs) / tmp
           Rsin(iKs)       = HBrg(iKs+1,iKs) / tmp
           HBrg(iKs,  iKs) = tmp
           HBrg(iKs+1,iKs) = zero
c     
c.... rotate eBrg    R_i E
c     
           tmp         = Rcos(iKs) * eBrg(iKs) + Rsin(iKs) * eBrg(iKs+1)
           eBrg(iKs+1) =-Rsin(iKs) * eBrg(iKs) + Rcos(iKs) * eBrg(iKs+1)
           eBrg(iKs)   = tmp
c     
c.... check for convergence
c     
           ntotGM = ntotGM + 1
           echeck=abs(eBrg(iKs+1))
           if (echeck .le. epsnrm) exit
c     
c.... end of GMRES iteration loop
c     
 1000   continue
c
c.... ------------------------->   Solution   <------------------------
c
c.... if converged or end of Krylov space
c
c.... solve for yBrg
c
        do jK = iKs, 1, -1
           yBrg(jK) = eBrg(jK) / HBrg(jK,jK)
           do lK = 1, jK-1
              eBrg(lK) = eBrg(lK) - yBrg(jK) * HBrg(lK,jK)
           enddo
        enddo
c     
c.... update Dy
c
        do jK = 1, iKs
           Dy = Dy + yBrg(jK) * uBrg(:,:,jK)
        enddo
c     
c.... flop count
c
   !      flops = flops + (3*iKs+1)*nflow*nshg
c
c.... check for convergence
c     

        echeck=abs(eBrg(iKs+1))
        if (echeck .le. epsnrm) exit
        if(myrank.eq.master) write(*,*)'solver tolerance %satisfaction',
     &  (one-echeck/unorm)/(one-etol)*100
c     
c.... end of mGMRES loop
c
 2000 continue
c     
c.... ------------------------>   Converged   <------------------------
c     
c.... if converged
c     
 3000 continue
c
c.... back EBE precondition the results
c
ctoomuchmem      call i3Pcond (Binv,  Dy, ilwork, 'R_Pcond ')
c     
c.... back block-diagonal precondition the results 
c
      call i3LU (BDiag, Dy, 'backward')
c     
c
c.... output the statistics
c
              call rstat (res, ilwork,rmes) 
c    
c.... stop the timer
c     
 3002 continue                  ! no solve just res.
      call timer ('Back    ')
c     
c.... end
c     
      return
      end





      subroutine SolGMRs(y,         ac,        yold,      acold,
     &			 x,         iBC,       BC,  
     &                   col,       row,       lhsk,         
     &                   res,       BDiag,     HBrg,      eBrg,
     &                   yBrg,      Rcos,      Rsin,      iper,
     &                   ilwork,    shp,       shgl,      shpb,
     &                   shglb,     Dy,        rerr)
c
c----------------------------------------------------------------------
c
c  This is the preconditioned GMRES driver routine.
c
c input:
c  y      (nshg,ndof)           : Y-variables at n+alpha_v
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  yold   (nshg,ndof)           : Y-variables at beginning of step
c  acold  (nshg,ndof)           : Primvar. accel. variable at begng step
c  x      (numnp,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  HBrg   (Kspace+1,Kspace)      : Hessenberg matrix (LHS matrix)
c  eBrg   (Kspace+1)             : RHS      of Hessenberg minim. problem
c  yBrg   (Kspace)               : solution of Hessenberg minim. problem
c  Rcos   (Kspace)               : Rotational cosine of QR algorithm
c  Rsin   (Kspace)               : Rotational sine   of QR algorithm
c  shp(b) (nen,maxsh,melCat)     : element shape functions (boundary)
c  shgl(b)(nsd,nen,maxsh,melCat) : local gradients of shape functions
c
c output:
c  res    (nshg,nflow)           : preconditioned residual
c  BDiag  (nshg,nflow,nflow)      : block-diagonal preconditioner
c
c  
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use pointer_data
        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
c
      integer col(nshg+1), row(nnz*nshg)
      real*8 lhsK(nflow*nflow,nnz_tot)


      dimension y(nshg,ndof),             ac(nshg,ndof),
     &          yold(nshg,ndof),          acold(nshg,ndof),
     &          x(numnp,nsd),
     &          iBC(nshg),                BC(nshg,ndofBC),
     &          res(nshg,nflow),
     &          BDiag(nshg,nflow,nflow),
     &          HBrg(Kspace+1,Kspace),    eBrg(Kspace+1),
     &          yBrg(Kspace),             Rcos(Kspace),
     &          Rsin(Kspace),             ilwork(nlwork),
     &          iper(nshg)
c
      dimension Dy(nshg,nflow),            rmes(nshg,nflow),
     &          temp(nshg,nflow),
     &          uBrg(nshg,nflow,Kspace+1)
c        
      dimension shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
      real*8    rerr(nshg,10)
c      
c     
c.... *******************>> Element Data Formation <<******************
c
c
c.... set the parameters for flux and surface tension calculations
c
c
      idflx = 0 
      if(idiff >= 1)  idflx= idflx + (nflow-1) * nsd
      if (isurf == 1) idflx=idflx + nsd
c
c.... form the LHS matrices, the residual vector, and the block
c     diagonal preconditioner
c
      call ElmGMRs(y,             ac,            x,
     &             shp,           shgl,          iBC,
     &             BC,            shpb,
     &             shglb,         res,
     &             rmes,          BDiag,         iper,      
     &             ilwork,        lhsK,          col, 
     &             row,           rerr )
      rmes=res  ! saving the b vector (residual) 
c    

	call tnanq(res,5, 'res_egmr')
	call tnanq(BDiag,25, 'bdg_egmr')
c
c.... **********************>>    EBE - GMRES    <<********************
c
      call timer ('Solver  ')
c
c.... ------------------------> Initialization <-----------------------
c
c
c.... LU decompose the block diagonals
c
      if (iprec .ne. 0) then
         call i3LU (BDiag, res,  'LU_Fact ')
         if (numpe > 1) then
            call commu (BDiag  , ilwork, nflow*nflow  , 'out')
         endif
      endif
c
c.... block diagonal precondition residual
c
      call i3LU (BDiag, res,  'forward ')
!  from this point forward b is btilde (Preconditioned residual)
c
c Check the residual for divering trend
c
	call rstatCheck(res,ilwork,y,ac)
c
c.... initialize Dy
c
      Dy = zero
c
c.... Pre-precondition the LHS mass matrix and set up the sparse 
c     preconditioners
c

      if(lhs.eq.1) call Spsi3pre (BDiag,    lhsK,  col, row)
c     
c.... copy res in uBrg(1)
c     
      uBrg(:,:,1) = res
c     
c.... calculate norm of residual
c
      temp  = res**2

      call sumgat (temp, nflow, summed, ilwork)
      unorm = sqrt(summed)
c
c.... check if GMRES iterations are required
c
      iKs    = 0
      lGMRESs = 0
c
c.... if we are down to machine precision, don't bother solving
c
      if (unorm .lt. 100.*epsM**2) goto 3000 
c
c.... set up tolerance of the Hessenberg's problem
c
      epsnrm = etol * unorm
c
c.... ------------------------>  GMRES Loop  <-------------------------
c
c.... loop through GMRES cycles
c
      do 2000 mGMRES = 1, nGMRES
         lGMRESs = mGMRES - 1
c
         if (lGMRES .gt. 0) then
c
c.... if GMRES restarts are necessary, calculate  R - A x
c
c
c.... right precondition Dy
c
            temp = Dy
           
c
c.... perform the A x product
c
            call SparseAp (iper,ilwork,iBC, col, row, lhsK,  temp)
c           call tnanq(temp,5, 'q_spAPrs')

c     
c.... periodic nodes have to assemble results to their partners
c
            call bc3per (iBC,  temp,  iper, ilwork, nflow) 
c           call tnanq(temp,5, 'q_BCprs')
c
c.... subtract A x from residual and calculate the norm
c           
            temp = res - temp
            uBrg(:,:,1) = temp
c
c.... calculate the norm
c
            temp  = temp**2
            call sumgat (temp, nflow, summed, ilwork)
            unorm = sqrt(summed)
c     
c.... flop count
c     
       !      flops = flops + nflow*nshg+nshg
c     
         endif
c
c.... set up RHS of the Hessenberg's problem
c
         call clear (eBrg, Kspace+1)
         eBrg(1) = unorm
c
c.... normalize the first Krylov vector
c
         uBrg(:,:,1) = uBrg(:,:,1) / unorm
c
c.... loop through GMRES iterations
c
         do 1000 iK = 1, Kspace
            iKs = iK

            uBrg(:,:,iKs+1) = uBrg(:,:,iKs)
c
c.... Au product  ( u_{i+1} <-- EGmass u_{i+1} )
c
            call SparseAp (iper, ilwork, iBC,
     &                     col,  row,    lhsK,
     &                     uBrg(:,:,iKs+1) )
c           call tnanq(uBrg(:,:,iKS+1),5, 'q_spAP')

c     
c.... periodic nodes have to assemble results to their partners
c
            call bc3per (iBC,  uBrg(:,:,iKs+1),  iper, ilwork, nflow)
c           call tnanq(uBrg(:,:,iKS+1),5, 'q_bc')

c
c.... orthogonalize and get the norm
c
            do jK = 1, iKs+1  
c
               if (jK .eq. 1) then
c
                  temp = uBrg(:,:,iKs+1) * uBrg(:,:,1) ! {u_{i+1}*u_1} vector 
                  call sumgat (temp, nflow, beta, ilwork) ! sum vector=(u_{i+1},u_1)
c
               else
c
c project off jK-1 vector
c
                  uBrg(:,:,iKs+1)=uBrg(:,:,iKs+1)-beta * uBrg(:,:,jK-1)
c
                  temp = uBrg(:,:,iKs+1) * uBrg(:,:,jK) !{u_{i+1}*u_j} vector
                  call sumgat (temp, nflow, beta, ilwork) ! sum vector=(u_{i+1},u_j)
c
               endif
c
               HBrg(jK,iKs) = beta ! put this in the Hessenberg Matrix
c
            enddo
c
       !      flops = flops + (3*iKs+1)*nflow*numnp+(iKs+1)*numnp
c
c  the last inner product was with what was left of the vector (after
c  projecting off all of the previous vectors
c
        if(beta.le.0) write(*,*) 'beta in solgmr non-positive'
            unorm           = sqrt(beta)
            HBrg(iKs+1,iKs) = unorm ! this fills the 1 sub diagonal band
c
c.... normalize the Krylov vector
c
            uBrg(:,:,iKs+1) = uBrg(:,:,iKs+1) / unorm ! normalize the next Krylov
c vector
c
c.... construct and reduce the Hessenberg Matrix
c  since there is only one subdiagonal we can use a Givens rotation to 
c  rotate off each subdiagonal AS IT IS FORMED.   We do this because it
c  allows us to check progress of solution and quit when satisfied.  Note
c  that all future K vects will put a subdiagonal in the next column so
c  there is no penalty to work ahead as  the rotation for the next vector
c  will be unaffected by this rotation.
        
c     
c     H Y = E ========>   R_i H Y = R_i E
c     
            do jK = 1, iKs-1
               tmp            =  Rcos(jK) * HBrg(jK,  iKs) +
     &                           Rsin(jK) * HBrg(jK+1,iKs)
               HBrg(jK+1,iKs) = -Rsin(jK) * HBrg(jK,  iKs) +
     &                           Rcos(jK) * HBrg(jK+1,iKs)
               HBrg(jK,  iKs) =  tmp
            enddo
c     
            tmp            = sqrt(HBrg(iKs,iKs)**2 + HBrg(iKs+1,iKs)**2)
            Rcos(iKs)      = HBrg(iKs,  iKs) / tmp
            Rsin(iKs)      = HBrg(iKs+1,iKs) / tmp
            HBrg(iKs,  iKs)= tmp
            HBrg(iKs+1,iKs)= zero
c     
c.... rotate eBrg    R_i E
c     
            tmp        = Rcos(iKs) * eBrg(iKs) + Rsin(iKs) * eBrg(iKs+1)
            eBrg(iKs+1)=-Rsin(iKs) * eBrg(iKs) + Rcos(iKs) * eBrg(iKs+1)
            eBrg(iKs)  = tmp
c     
c.... check for convergence
c     
            ntotGM = ntotGM + 1
            echeck=abs(eBrg(iKs+1))
            if (echeck .le. epsnrm.and. iKs .ge. minIters) exit
c     
c.... end of GMRES iteration loop
c     
 1000    continue
c
c.... ------------------------->   Solution   <------------------------
c
c.... if converged or end of Krylov space
c
c.... solve for yBrg
c
         do jK = iKs, 1, -1
            yBrg(jK) = eBrg(jK) / HBrg(jK,jK)
            do lK = 1, jK-1
               eBrg(lK) = eBrg(lK) - yBrg(jK) * HBrg(lK,jK)
            enddo
         enddo
c     
c.... update Dy
c
         do jK = 1, iKs
            Dy = Dy + yBrg(jK) * uBrg(:,:,jK)
         enddo
c     
c.... flop count
c
    !      flops = flops + (3*iKs+1)*nflow*nshg
c
c.... check for convergence
c     
        echeck=abs(eBrg(iKs+1))
        if (echeck .le. epsnrm) exit
!        if(myrank.eq.master) write(*,*)'solver tolerance %satisfaction',
!     &  (one-echeck*etol/epsnrm)/(one-etol)*100

c     
c.... end of mGMRES loop
c
 2000 continue
c     
c.... ------------------------>   Converged   <------------------------
c     
c.... if converged
c     
 3000 continue

c
c     
c.... back block-diagonal precondition the results 
c
      call i3LU (BDiag, Dy, 'backward')
c     
c
c.... output the statistics
c
      call rstat (res, ilwork,rmes) 
      
      if(myrank.eq.master) then
        if (echeck .le. epsnrm) then
            write(*,*)
        else
            write(*,*)'solver tolerance %satisfaction',
     &  (one-echeck*etol/epsnrm)/(one-etol)*100
        endif
      endif
c    
c.... stop the timer
c     
 3002 continue                  ! no solve just res.
      call timer ('Back    ')
c     
c.... end
c     
      return
      end

        subroutine SolGMRSclr(y,       ac,      yold,
     &                        acold,   EGmasst,
     &                        x,       elDw, 
     &                        iBC,      BC,           
     &                        rest,     HBrg,     eBrg,
     &                        yBrg,     Rcos,     Rsin,    iper,
     &                        ilwork,
     &                        shp,      shgl,
     &                        shpb,     shglb,  Dyt)
c
c----------------------------------------------------------------------
c
c  This is the preconditioned GMRES driver routine.
c
c input:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  yold   (nshg,ndof)           : Y-variables at beginning of step
c  acold  (nshg,ndof)           : Primvar. accel. variable at begng step
c  x      (numnp,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodic nodal information
c
c output:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  HBrg   (Kspace+1,Kspace)      : Hessenberg matrix (LHS matrix)
c  eBrg   (Kspace+1)             : RHS      of Hessenberg minim. problem
c  yBrg   (Kspace)               : solution of Hessenberg minim. problem
c  Rcos   (Kspace)               : Rotational cosine of QR algorithm
c  Rsin   (Kspace)               : Rotational sine   of QR algorithm
c output:
c  rest   (numnp)           : preconditioned residual
c
c  
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use pointer_data
        
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension y(nshg,ndof),      ac(nshg,ndof),
     &            yold(nshg,ndof),   acold(nshg,ndof),
     &            x(numnp,nsd),
     &            iBC(nshg),         BC(nshg,ndofBC),
     &            rest(nshg),
     &            Diag(nshg),
     &            HBrg(Kspace+1,*),  eBrg(*),
     &            yBrg(*),           Rcos(*),
     &            Rsin(*),           ilwork(nlwork),
     &            iper(nshg),        EGmasst(numel,nshape,nshape)
c
        dimension Dyt(nshg),         rmest(nshg),
     &            tempt(nshg),       Dinv(nshg),
     &            uBrgt(nshg,Kspace+1)
c        
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)
        real*8    elDw(numel) 
c     
c.... *******************>> Element Data Formation <<******************
c
c
c.... form the LHS matrices, the residual vector, and the block
c     diagonal preconditioner
c
        call ElmGMRSclr(y,             ac, 
     &                  x,             elDw,      shp,        shgl, 
     &                  iBC,           BC,
     &                  shpb,          shglb,
     &                  rest,
     &                  rmest,         Diag,       iper,      
     &                  ilwork,        EGmasst)
c
c.... **********************>>    EBE - GMRES    <<********************
c
        call timer ('Solver  ')
c
c.... ------------------------> Initialization <-----------------------
c
c
      id = isclr+5
c.... initialize Dy
c
        Dyt = zero
c
c.... Right preconditioner
c
        call i3preSclr(Diag, Dinv, EGmassT, ilwork)
c
c Check the residual for divering trend
c

        call rstatCheckSclr(rest,ilwork,y,ac)

c     
c.... copy rest in uBrgt(1)
c     
        uBrgt(:,1) = rest
c     
c.... calculate norm of residual
c
        tempt  = rest**2

        call sumgat (tempt, 1, summed, ilwork)
        unorm = sqrt(summed)
c
c.... check if GMRES iterations are required
c
        iKss    = 0
        lGMRESt = 0
c
c.... if we are down to machine precision, don't bother solving
c
        if (unorm .lt. 100.*epsM**2) goto 3000 
c
c.... set up tolerance of the Hessenberg's problem
c
        epsnrm = etol * unorm 
c
c.... ------------------------>  GMRES Loop  <-------------------------
c
c.... loop through GMRES cycles
c
        do 2000 mGMRES = 1, nGMRES
        lGMRESt = mGMRES - 1
c
        if (lGMRESt .gt. 0) then
c
c.... if GMRES restarts are necessary, calculate  R - A x
c
c

c.... perform the A x product
c
           call Au1GMRSclr (EGmasst,  tempt,  ilwork, iper)
c
c.... periodic nodes have to assemble results to their partners
c
c          call bc3perSclr (iBC,  tempt,  iper)
c
c.... subtract A x from residual and calculate the norm
c           
           tempt = rest - tempt
           uBrgt(:,1) = tempt
c
c.... calculate the norm
c
           tempt  = tempt**2
           call sumgat (tempt, 1, summed, ilwork)
           unorm = sqrt(summed)
c     
c.... flop count
c     
      !      flops = flops + ndof*numnp+numnp
c     
        endif
c
c.... set up RHS of the Hessenberg's problem
c
        call clear (eBrg, Kspace+1)
        call clear (HBrg, Kspace+1)
        eBrg(1) = unorm
c
c.... normalize the first Krylov vector
c
        uBrgt(:,1) = uBrgt(:,1) / unorm
c
c.... loop through GMRES iterations
c
        do 1000 iK = 1, Kspace
           iKss = iK

           uBrgt(:,iKss+1) = uBrgt(:,iKss)

c.... Au product  ( u_{i+1} <-- EGmass u_{i+1} )
c
           call Au1GMRSclr ( EGmasst, uBrgt(:,iKss+1),  ilwork, iper )

c
c.... periodic nodes have to assemble results to their partners
c
           call bc3perSclr (iBC,  uBrgt(:,iKss+1),  iper)
c
c.... orthogonalize and get the norm
c
          do jK = 1, iKss+1  
c
            if (jK .eq. 1) then
c
              tempt = uBrgt(:,iKss+1) * uBrgt(:,1)  ! {u_{i+1}*u_1} vector 
              call sumgat (tempt, 1, beta, ilwork) ! sum vector=(u_{i+1},u_1)
c
            else
c
c project off jK-1 vector
c
          uBrgt(:,iKss+1) = uBrgt(:,iKss+1) - beta * uBrgt(:,jK-1)
c
              tempt = uBrgt(:,iKss+1) * uBrgt(:,jK) !{u_{i+1}*u_j} vector
              call sumgat (tempt, 1, beta, ilwork) ! sum vector=(u_{i+1},u_j)
c
            endif
c
            HBrg(jK,iKss) = beta   ! put this in the Hessenberg Matrix
c
        enddo
c
   !      flops = flops + (3*iKss+1)*ndof*numnp+(iKss+1)*numnp
c
c  the last inner product was with what was left of the vector (after
c  projecting off all of the previous vectors
c
        unorm           = sqrt(beta)
        HBrg(iKss+1,iKss) = unorm   ! this fills the 1 sub diagonal band
c
c.... normalize the Krylov vector
c
        uBrgt(:,iKss+1) = uBrgt(:,iKss+1) / unorm  ! normalize the next Krylov
c vector
c
c.... construct and reduce the Hessenberg Matrix
c  since there is only one subdiagonal we can use a Givens rotation to 
c  rotate off each subdiagonal AS IT IS FORMED.   We do this because it
c  allows us to check progress of solution and quit when satisfied.  Note
c  that all future K vects will put a subdiagonal in the next column so
c  there is no penalty to work ahead as  the rotation for the next vector
c  will be unaffected by this rotation.
        
c     
c     H Y = E ========>   R_i H Y = R_i E
c     
           do jK = 1, iKss-1
              tmp            =  Rcos(jK) * HBrg(jK,  iKss) +
     &                          Rsin(jK) * HBrg(jK+1,iKss)
              HBrg(jK+1,iKss) = -Rsin(jK) * HBrg(jK,  iKss) +
     &                          Rcos(jK) * HBrg(jK+1,iKss)
              HBrg(jK,  iKss) =  tmp
           enddo
c     
           tmp        = sqrt(HBrg(iKss,iKss)**2 + HBrg(iKss+1,iKss)**2)
           Rcos(iKss) = HBrg(iKss,  iKss) / tmp
           Rsin(iKss) = HBrg(iKss+1,iKss) / tmp
           HBrg(iKss,  iKss) = tmp
           HBrg(iKss+1,iKss) = zero
c     
c.... rotate eBrg    R_i E
c     
           tmp         = Rcos(iKss)*eBrg(iKss) + Rsin(iKss)*eBrg(iKss+1)
           eBrg(iKss+1)=-Rsin(iKss)*eBrg(iKss) + Rcos(iKss)*eBrg(iKss+1)
           eBrg(iKss)  = tmp
c     
c.... check for convergence
c     
           ercheck=eBrg(iKss+1)
           ntotGMs = ntotGMs + 1
           if (abs(eBrg(iKss+1)) .le. epsnrm) exit
c     
c.... end of GMRES iteration loop
c     
 1000   continue
c
c.... ------------------------->   Solution   <------------------------
c
c.... if converged or end of Krylov space
c
c.... solve for yBrg
c
        do jK = iKss, 1, -1
           yBrg(jK) = eBrg(jK) / HBrg(jK,jK)
           do lK = 1, jK-1
              eBrg(lK) = eBrg(lK) - yBrg(jK) * HBrg(lK,jK)
           enddo
        enddo
c     
c.... update Dy
c
        do jK = 1, iKss
           Dyt = Dyt + yBrg(jK) * uBrgt(:,jK)
        enddo
c     
c.... flop count
c
   !      flops = flops + (3*iKss+1)*ndof*numnp
c
c.... check for convergence
c     
        if (abs(eBrg(iKss+1)) .le. epsnrm) exit
c     
c.... end of mGMRES loop
c
 2000 continue
c     
c.... ------------------------>   Converged   <------------------------
c     
c.... if converged
c     
 3000 continue
c
c.... back precondition the result
c
      Dyt(:) = Dyt(:) * Dinv(:)
c     
c.... output the statistics
c
      call rstatSclr(rest, ilwork)
c.... stop the timer
c     
 3002 continue                  ! no solve just res.
      call timer ('Back    ')
c     
c.... end
c     
      return
      end


