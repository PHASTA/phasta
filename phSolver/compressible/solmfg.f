        subroutine SolMFG (y,         ac,        yold,      acold,
     &			   x,         
     &                     iBC,       BC,        res,
     &                     BDiag,     HBrg,      eBrg,
     &                     yBrg,      Rcos,      Rsin,      iper,
     &                     ilwork,    shp,       shgl,     
     &                     shpb,      shglb,     Dy, rerr)
c
c----------------------------------------------------------------------
c
c  This is the preconditioned matrix-free GMRES driver routine.
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
c  res    (nshg,ndof)           : preconditioned residual
c  BDiag  (nshg,ndof,ndof)      : block-diagonal preconditioner
c
c  
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
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
     &            HBrg(Kspace+1,Kspace),    eBrg(Kspace+1),
     &            yBrg(Kspace),             Rcos(Kspace),
     &            Rsin(Kspace),             ilwork(nlwork),
     &            iper(nshg)
c
        dimension Dy(nshg,nflow),            rmes(nshg,nflow),
     &            ypre(nshg,nflow),          temp(nshg,nflow),
     &            uBrg(nshg,nflow,Kspace+1),  ytmp(nshg,nflow)

c        
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)       
c
c.... *******************>> Element Data Formation <<******************
c
c
c.... set the parameters for flux and surface tension calculations
c
        idflx = zero
        if(idiff >= 1)  idflx= idflx + (nflow-1) * nsd
        if (isurf == 1) idflx=idflx + nsd
c
        call ElmMFG (y,             ac,            x,
     &               shp,           shgl,
     &               iBC,           BC,
     &               shpb,          shglb,
     &               res,           rmes,
     &               BDiag,         iper,          ilwork,
     &               rerr)
c
c.... **********************>> Matrix-Free GMRES <<********************
c
c
c.... start the timer
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
c  This is a feature that allows one to take an extra pass just to
c  find the residual at the end of the last solve.
c
c$$$        if(iter.ge.(press*nitr) ) then
c$$$          iKs=0
c$$$          lGMRES=0
c$$$          goto 3002
c$$$        endif

c
c.... block diagonal precondition modified residual
c
        call i3LU (BDiag, rmes, 'forward ')

c
c.... copy res in uBrg(1)
c
        uBrg(:,:,1) = res
c
c.... initialize Dy
c
        Dy = zero
c
c.... block diagonal precondition y-variables
c
        ypre(:,:) = y(:,1:nflow) ! ypre is the pre-conditioned,
                                 ! unperturbed, base vector
c
        call yshuffle(ypre,'new2old ')
c
        call i3LU (BDiag, ypre,   'product ')
c
c  since we will never use ypre in the "new" form again, leave it
c  shuffled
c
c        call yshuffle(ypre, 'old2new ')
c
c.... calculate norm of residual
c
        temp  = res**2

c        call tnanq(temp,5,"res**2  ")
        
        call sumgat (temp, nflow, summed, ilwork)
        unorm = sqrt(summed)
c
c.... flop count
c
!        flops = flops + ndof*nshg+nshg
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
c.... compute the finite difference interval
c
        if ((iter .eq. 1) .and. (mod(istep,20) .eq. 0)) then
           call itrFDI (ypre,            y,            ac,
     &                  x,               rmes,
     &                  res,             BDiag,         iBC,
     &                  BC,              iper,
     &                  ilwork,          shp,           shgl,
     &                  shpb,            shglb)
        endif
        ires=2
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
c.... calculate  R - A u
c
          call Au2MFG (ypre,        y,        ac,
     &                 x,           rmes,
     &                 res,         Dy,          temp,
     &                 BDiag,       iBC,         BC,
     &                 iper,        ilwork,
     &                 shp,         shgl,        
     &                 shpb,        shglb)

c
          uBrg(:,:,1) = temp
c
c.... calculate the norm
c
          temp  = temp**2
          call sumgat (temp, ndof, summed, ilwork)
          unorm = sqrt(summed)

c
c.... flop count
c
!          flops = flops + ndof*nshg+nshg
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
c
c.... Au product
c
        temp = uBrg(:,:,iKs)
c
        call Au1MFG (ypre,        y,           ac,
     &               x,           rmes,
     &               res,         temp,        BDiag,
     &               iBC,         BC,          
     &               iper,        ilwork,
     &               shp,         shgl,         
     &               shpb,        shglb)
c
        uBrg(:,:,iKs+1) = temp   ! u_{i+1}= J u_i  In Johan Thesis p 15c
c$$$c
c$$$c.... debug
c$$$c
c$$$           do i=1,nshg
c$$$              write(78,'(5(f14.7))')(uBrg(i,j,iKs+1),j=1,5)
c$$$           enddo
c$$$           stop
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
!        flops = flops + (3*iKs+1)*nflow*nshg+(iKs+1)*nshg
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
c   will be unaffected by this rotation.

c
c   H Y = E ========>   R_i H Y = R_i E
c
        do jK = 1, iKs-1
          tmp            =  Rcos(jK) * HBrg(jK,  iKs) +
     &                      Rsin(jK) * HBrg(jK+1,iKs)
          HBrg(jK+1,iKs) = -Rsin(jK) * HBrg(jK,  iKs) +
     &                      Rcos(jK) * HBrg(jK+1,iKs)
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
        tmp         =  Rcos(iKs) * eBrg(iKs) + Rsin(iKs) * eBrg(iKs+1)
        eBrg(iKs+1) = -Rsin(iKs) * eBrg(iKs) + Rcos(iKs) * eBrg(iKs+1)
        eBrg(iKs)   =  tmp
c
c.... check for convergence
c
        ercheck=eBrg(iKs+1)
        ntotGM = ntotGM + 1
        if (abs(eBrg(iKs+1)) .le. epsnrm) exit
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
!        flops = flops + (3*iKs+1)*nflow*nshg
c
c.... check for convergence
c
        if (abs(eBrg(iKs+1)) .le. epsnrm) exit
c
c.... end of mGMRES loop
c
2000    continue
c
c.... ------------------------>   Converged   <------------------------
c
c.... if converged
c
3000    continue

c
c.... back precondition the results 
c
        call i3LU (BDiag, Dy, 'backward')
c
c.... output the statistics
c
              call rstat (res, ilwork) 
c 
c ... reset ires to 3 again (asires changed ires to 2)
c    
              ires = 3
c
c.... stop the timer
c
3002    continue  ! no solve just res.
        call timer ('Back    ')
c
c.... end
c
        return
        end
