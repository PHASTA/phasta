        subroutine e3LS   (A1,        A2,          A3,    
     &                     rho,       rmu,         cp,
     &                     cv,        con,         T,
     &                     u1,        u2,          u3, 
     &                     rLyi,      dxidx,       tau,   
     &                     ri,        rmi,         rk,
     &                     dui,       aci,         A0,
     &                     divqi,     shape,       shg,
     &                     EGmass,    stiff,       WdetJ,
     &                     giju,      rTLS,        raLS,
     &                     A0inv,     dVdY,        rerrl,
     &                     compK,     pres,        PTau)
c                                                                
c----------------------------------------------------------------------
c
c This routine calculates the contribution of the least-squares 
c operator to the RHS vector and LHS tangent matrix. The temporary 
c results are put in ri.
c
c input:
c  A1    (npro,nflow,nflow)     : A_1
c  A2    (npro,nflow,nflow)     : A_2
c  A3    (npro,nflow,nflow)     : A_3
c  rho   (npro)               : density
c  T     (npro)               : temperature
c  cp    (npro)               : specific heat at constant pressure
c  u1    (npro)               : x1-velocity component
c  u2    (npro)               : x2-velocity component
c  u3    (npro)               : x3-velocity component
c  rLyi  (npro,nflow)          : least-squares residual vector
c  dxidx (npro,nsd,nsd)       : inverse of deformation gradient
c  tau   (npro,3)             : stability parameter
c  PTau  (npro,5,5)           : matrix of stability parameters
c  rLyi  (npro,nflow)          : convective portion of least-squares
c                               residual vector 
c  divqi (npro,nflow-1)        : divergence of diffusive flux
c  shape (npro,nshl)        : element shape functions
c  shg   (npro,nshl,nsd)    : global element shape function grads
c  WdetJ (npro)               : weighted jacobian determinant
c  stiff (npro,nsd*nflow,nsd*nflow) : stiffness matrix
c  EGmass(npro,nedof,nedof)   : partial mass matrix
c  compK (npro,10)             : K_ij in (eq.134) A new ... III 
c
c output:
c  ri     (npro,nflow*(nsd+1)) : partial residual
c  rmi    (npro,nflow*(nsd+1)) : partial modified residual
c  EGmass (npro,nedof,nedof)  : partial mass matrix
c
c
c Zdenek Johan, Summer 1990. (Modified from e2ls.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997. Prim. Var. with Diag Tau
c----------------------------------------------------------------------
c
        include "common.h"

c
c  passed arrays
c
        dimension A1(npro,nflow,nflow),    A2(npro,nflow,nflow),
     &            A3(npro,nflow,nflow),    cv(npro),
     &            A0(npro,nflow,nflow),    rho(npro),
     &            con(npro),               cp(npro),
     &            dui(npro,nflow),         aci(npro,nflow),
     &            u1(npro),                u2(npro),
     &            u3(npro),                rk(npro),
     &            rLyi(npro,nflow),        dxidx(npro,nsd,nsd),
     &            tau(npro,5),             giju(npro,6),
     &            rTLS(npro),              raLS(npro),
     &            ri(npro,nflow*(nsd+1)),  rmi(npro,nflow*(nsd+1)),
     &            divqi(npro,nflow-1),     stiff(npro,3*nflow,3*nflow),
     &            EGmass(npro,nedof,nedof),shape(npro,nshl),
     &            shg(npro,nshl,nsd),      WdetJ(npro),
     &            PTau(npro,5,5),          T(npro),
     &            pres(npro)
c
c local arrays
c
        dimension rLymi(npro,nflow),         Atau(npro,nflow,nflow),
     &            A1tauA0(npro,nflow,nflow), A2tauA0(npro,nflow,nflow),
     &            A3tauA0(npro,nflow,nflow), fact(npro),
     &            A0inv(npro,15),            dVdY(npro,15),
     &            compK(npro,10),            ac1(npro),
     &            ac2(npro),                 ac3(npro)
c
        real*8    rerrl(npro,nshl,6), tmp(npro), tmp1(npro)
        ttim(24) = ttim(24) - secs(0.0)
c
c
c last step to the least squares is adding the time term.  So that we
c only have to localize one vector for each Krylov vector the modified
c residual is quite different from the total residual.
c
c
c the modified residual
c
       fct1=almi/gami/alfi*dtgl
c
c  add possibility of not including time term
c
c       if(idiff.ne.-1)
        
       if(ires.ne.1) rLymi = rLyi + fct1*dui
c       
       if(ires.eq.1 .or. ires .eq. 3) then
c       rLymi = rLyi

        rLyi(:,1) = rLyi(:,1) 
     &            + A0(:,1,1)*aci(:,1)
c    &            + A0(:,1,2)*aci(:,2)
c    &            + A0(:,1,3)*aci(:,3)
c    &            + A0(:,1,4)*aci(:,4)
     &            + A0(:,1,5)*aci(:,5)
c
        rLyi(:,2) = rLyi(:,2) 
     &            + A0(:,2,1)*aci(:,1)
     &            + A0(:,2,2)*aci(:,2)
c    &            + A0(:,2,3)*aci(:,3)
c    &            + A0(:,2,4)*aci(:,4)
     &            + A0(:,2,5)*aci(:,5)
c
        rLyi(:,3) = rLyi(:,3) 
     &            + A0(:,3,1)*aci(:,1)
c    &            + A0(:,3,2)*aci(:,2)
     &            + A0(:,3,3)*aci(:,3)
c    &            + A0(:,3,4)*aci(:,4)
     &            + A0(:,3,5)*aci(:,5)
c
        rLyi(:,4) = rLyi(:,4) 
     &            + A0(:,4,1)*aci(:,1)
c    &            + A0(:,4,2)*aci(:,2)
c    &            + A0(:,4,3)*aci(:,3)
     &            + A0(:,4,4)*aci(:,4)
     &            + A0(:,4,5)*aci(:,5)
c
        rLyi(:,5) = rLyi(:,5) 
     &            + A0(:,5,1)*aci(:,1)
     &            + A0(:,5,2)*aci(:,2)
     &            + A0(:,5,3)*aci(:,3)
     &            + A0(:,5,4)*aci(:,4)
     &            + A0(:,5,5)*aci(:,5)
c
      endif
c
c.... subtract div(q) from the least squares term
c
      if ((idiff >= 1).and.(ires==3 .or. ires==1)) then
c
      if (isurf.eq.zero) then
         rLyi(:,2) = rLyi(:,2) - divqi(:,1)
         rLyi(:,3) = rLyi(:,3) - divqi(:,2)
         rLyi(:,4) = rLyi(:,4) - divqi(:,3)       
         rLyi(:,5) = rLyi(:,5) - divqi(:,4)
      endif
      endif
c
c.... -------------------> error calculation  <-----------------
c     
       if((ierrcalc.eq.1).and.(nitr.eq.iter)) then
          do ia=1,nshl
             tmp=shape(:,ia)*WdetJ(:)
             tmp1=shape(:,ia)*Qwt(lcsyst,intp)
             rerrl(:,ia,1) = rerrl(:,ia,1) +
     &                       tmp1(:)*rLyi(:,1)*rLyi(:,1)
             rerrl(:,ia,2) = rerrl(:,ia,2) +
     &                       tmp1(:)*rLyi(:,2)*rLyi(:,2)
             rerrl(:,ia,3) = rerrl(:,ia,3) +
     &                       tmp1(:)*rLyi(:,3)*rLyi(:,3)

             rerrl(:,ia,4) = rerrl(:,ia,4) +
     &                       tmp(:)*divqi(:,1)*divqi(:,1)
             rerrl(:,ia,5) = rerrl(:,ia,5) +
     &                       tmp(:)*divqi(:,2)*divqi(:,2)
             rerrl(:,ia,6) = rerrl(:,ia,6) +
     &                       tmp(:)*divqi(:,3)*divqi(:,3)
          enddo
       endif
c      
c
c.... --------------------------->  Tau  <-----------------------------
c
c.... calculate the tau matrix
c
c.... in the first incarnation we will ONLY have a diagonal tau here

       if (itau .lt. 10) then    ! diagonal tau

          call e3tau  (rho,             cp,		rmu,
     &         u1,              u2,             u3,
     &         con,             dxidx,          rLyi,  
     &         rLymi,           tau,            rk,
     &         giju,            rTLS,           raLS,
     &         A0inv,           dVdY,           cv)	
          
       else

c.... looks like need a non-diagonal, discribed in "A new ... III"
c.... by Hughes and Mallet. Original work - developed diffusive (and as
c.... well advective) correction factors for 1-D a-d equation w/ hier. b. 

          
          ac1(:) = aci(:,2)
          ac2(:) = aci(:,3)
          ac3(:) = aci(:,4)

          call e3tau_nd  (rho,       cp,  rmu,   T,
     &         u1,              u2,             u3,
     &         ac1,             ac2,             ac3,
     &         con,             dxidx,          rLyi,  
     &         rLymi,           PTau,           rk,
     &         giju,            rTLS,           raLS,
     &         cv,              compK,          pres,
     &         A0inv,           dVdY)

       endif
       

        ttim(25) = ttim(25) + secs(0.0)
c
c Warning:  to save space I have put the tau times the modified residual 
c           in rLymi and the tau times the total residual back in rLyi
c
c
c.... ---------------------------->  RHS  <----------------------------
c
c.... calculate (A_i^T tau Ly)
c

       if(ires.ne.1) then
c
c  A1 * Tau L(Y):  to be hit on left with Na,x in e3wmlt
c
        rmi(:,1) =  
     &               A1(:,1,1) * rLymi(:,1) 
     &             + A1(:,1,2) * rLymi(:,2)
c    &             + A1(:,1,3) * rLymi(:,3) 
c    &             + A1(:,1,4) * rLymi(:,4)
     &             + A1(:,1,5) * rLymi(:,5)
     &             + rmi(:,1)
        rmi(:,2) =
     &               A1(:,2,1) * rLymi(:,1) 
     &             + A1(:,2,2) * rLymi(:,2)
c    &             + A1(:,2,3) * rLymi(:,3) 
c    &             + A1(:,2,4) * rLymi(:,4)
     &             + A1(:,2,5) * rLymi(:,5)
     &             + rmi(:,2)
        rmi(:,3) =
     &               A1(:,3,1) * rLymi(:,1) 
     &             + A1(:,3,2) * rLymi(:,2)
     &             + A1(:,3,3) * rLymi(:,3) 
c    &             + A1(:,3,4) * rLymi(:,4)
     &             + A1(:,3,5) * rLymi(:,5)
     &             + rmi(:,3)
        rmi(:,4) =
     &               A1(:,4,1) * rLymi(:,1) 
     &             + A1(:,4,2) * rLymi(:,2)
c    &             + A1(:,4,3) * rLymi(:,3) 
     &             + A1(:,4,4) * rLymi(:,4)
     &             + A1(:,4,5) * rLymi(:,5)
     &             + rmi(:,4)
        rmi(:,5) =
     &               A1(:,5,1) * rLymi(:,1) 
     &             + A1(:,5,2) * rLymi(:,2)
     &             + A1(:,5,3) * rLymi(:,3) 
     &             + A1(:,5,4) * rLymi(:,4)
     &             + A1(:,5,5) * rLymi(:,5)
     &             + rmi(:,5)
c
c  A2 * Tau L(Y),  to be hit on left with Na,y 
c
        rmi(:,6) = 
     &               A2(:,1,1) * rLymi(:,1) 
c    &             + A2(:,1,2) * rLymi(:,2)
     &             + A2(:,1,3) * rLymi(:,3) 
c    &             + A2(:,1,4) * rLymi(:,4)
     &             + A2(:,1,5) * rLymi(:,5)
     &             + rmi(:,6)
        rmi(:,7) =
     &               A2(:,2,1) * rLymi(:,1) 
     &             + A2(:,2,2) * rLymi(:,2)
     &             + A2(:,2,3) * rLymi(:,3) 
c    &             + A2(:,2,4) * rLymi(:,4)
     &             + A2(:,2,5) * rLymi(:,5)
     &             + rmi(:,7)
        rmi(:,8) =
     &               A2(:,3,1) * rLymi(:,1) 
c    &             + A2(:,3,2) * rLymi(:,2)
     &             + A2(:,3,3) * rLymi(:,3) 
c    &             + A2(:,3,4) * rLymi(:,4)
     &             + A2(:,3,5) * rLymi(:,5)
     &             + rmi(:,8)
        rmi(:,9) =
     &               A2(:,4,1) * rLymi(:,1) 
c    &             + A2(:,4,2) * rLymi(:,2)
     &             + A2(:,4,3) * rLymi(:,3) 
     &             + A2(:,4,4) * rLymi(:,4)
     &             + A2(:,4,5) * rLymi(:,5)
     &             + rmi(:,9)
        rmi(:,10) =
     &               A2(:,5,1) * rLymi(:,1) 
     &             + A2(:,5,2) * rLymi(:,2)
     &             + A2(:,5,3) * rLymi(:,3) 
     &             + A2(:,5,4) * rLymi(:,4)
     &             + A2(:,5,5) * rLymi(:,5)
     &             + rmi(:,10)
c
c  A3 * Tau L(Y)  to be hit on left with Na,z
c
        rmi(:,11) = 
     &               A3(:,1,1) * rLymi(:,1) 
c    &             + A3(:,1,2) * rLymi(:,2)
c    &             + A3(:,1,3) * rLymi(:,3) 
     &             + A3(:,1,4) * rLymi(:,4)
     &             + A3(:,1,5) * rLymi(:,5)
     &             + rmi(:,11)
        rmi(:,12) =
     &               A3(:,2,1) * rLymi(:,1) 
     &             + A3(:,2,2) * rLymi(:,2)
c    &             + A3(:,2,3) * rLymi(:,3) 
     &             + A3(:,2,4) * rLymi(:,4)
     &             + A3(:,2,5) * rLymi(:,5)
     &             + rmi(:,12)
        rmi(:,13) =
     &               A3(:,3,1) * rLymi(:,1) 
c    &             + A3(:,3,2) * rLymi(:,2)
     &             + A3(:,3,3) * rLymi(:,3) 
     &             + A3(:,3,4) * rLymi(:,4)
     &             + A3(:,3,5) * rLymi(:,5)
     &             + rmi(:,13)
        rmi(:,14) =
     &               A3(:,4,1) * rLymi(:,1) 
c    &             + A3(:,4,2) * rLymi(:,2)
c    &             + A3(:,4,3) * rLymi(:,3) 
     &             + A3(:,4,4) * rLymi(:,4)
     &             + A3(:,4,5) * rLymi(:,5)
     &             + rmi(:,14)
        rmi(:,15) =
     &               A3(:,5,1) * rLymi(:,1) 
     &             + A3(:,5,2) * rLymi(:,2)
     &             + A3(:,5,3) * rLymi(:,3) 
     &             + A3(:,5,4) * rLymi(:,4)
     &             + A3(:,5,5) * rLymi(:,5)
     &             + rmi(:,15)
      endif  !ires.ne.1

c
c  same thing for the real residual
c
       if(ires.eq.3 .or. ires .eq. 1) then  ! we need the total residual
        ri(:,1) = 
     &               A1(:,1,1) * rLyi(:,1) 
     &             + A1(:,1,2) * rLyi(:,2)
c    &             + A1(:,1,3) * rLyi(:,3) 
c    &             + A1(:,1,4) * rLyi(:,4)
     &             + A1(:,1,5) * rLyi(:,5)
     &             + ri(:,1)
        ri(:,2) =
     &               A1(:,2,1) * rLyi(:,1) 
     &             + A1(:,2,2) * rLyi(:,2)
c    &             + A1(:,2,3) * rLyi(:,3) 
c    &             + A1(:,2,4) * rLyi(:,4)
     &             + A1(:,2,5) * rLyi(:,5)
     &             + ri(:,2)
        ri(:,3) =
     &               A1(:,3,1) * rLyi(:,1) 
     &             + A1(:,3,2) * rLyi(:,2)
     &             + A1(:,3,3) * rLyi(:,3) 
c    &             + A1(:,3,4) * rLyi(:,4)
     &             + A1(:,3,5) * rLyi(:,5)
     &             + ri(:,3)
        ri(:,4) =
     &               A1(:,4,1) * rLyi(:,1) 
     &             + A1(:,4,2) * rLyi(:,2)
c    &             + A1(:,4,3) * rLyi(:,3) 
     &             + A1(:,4,4) * rLyi(:,4)
     &             + A1(:,4,5) * rLyi(:,5)
     &             + ri(:,4)
        ri(:,5) =
     &               A1(:,5,1) * rLyi(:,1) 
     &             + A1(:,5,2) * rLyi(:,2)
     &             + A1(:,5,3) * rLyi(:,3) 
     &             + A1(:,5,4) * rLyi(:,4)
     &             + A1(:,5,5) * rLyi(:,5)
     &             + ri(:,5)
c
        ri(:,6) = 
     &               A2(:,1,1) * rLyi(:,1) 
c    &             + A2(:,1,2) * rLyi(:,2)
     &             + A2(:,1,3) * rLyi(:,3) 
c    &             + A2(:,1,4) * rLyi(:,4)
     &             + A2(:,1,5) * rLyi(:,5)
     &             + ri(:,6)
        ri(:,7) =
     &               A2(:,2,1) * rLyi(:,1) 
     &             + A2(:,2,2) * rLyi(:,2)
     &             + A2(:,2,3) * rLyi(:,3) 
c    &             + A2(:,2,4) * rLyi(:,4)
     &             + A2(:,2,5) * rLyi(:,5)
     &             + ri(:,7)
        ri(:,8) =
     &               A2(:,3,1) * rLyi(:,1) 
c    &             + A2(:,3,2) * rLyi(:,2)
     &             + A2(:,3,3) * rLyi(:,3) 
c    &             + A2(:,3,4) * rLyi(:,4)
     &             + A2(:,3,5) * rLyi(:,5)
     &             + ri(:,8)
        ri(:,9) =
     &               A2(:,4,1) * rLyi(:,1) 
c    &             + A2(:,4,2) * rLyi(:,2)
     &             + A2(:,4,3) * rLyi(:,3) 
     &             + A2(:,4,4) * rLyi(:,4)
     &             + A2(:,4,5) * rLyi(:,5)
     &             + ri(:,9)
        ri(:,10) =
     &               A2(:,5,1) * rLyi(:,1) 
     &             + A2(:,5,2) * rLyi(:,2)
     &             + A2(:,5,3) * rLyi(:,3) 
     &             + A2(:,5,4) * rLyi(:,4)
     &             + A2(:,5,5) * rLyi(:,5)
     &             + ri(:,10)
        ri(:,11) = 
     &               A3(:,1,1) * rLyi(:,1) 
c    &             + A3(:,1,2) * rLyi(:,2)
c    &             + A3(:,1,3) * rLyi(:,3) 
     &             + A3(:,1,4) * rLyi(:,4)
     &             + A3(:,1,5) * rLyi(:,5)
     &             + ri(:,11)
        ri(:,12) =
     &               A3(:,2,1) * rLyi(:,1) 
     &             + A3(:,2,2) * rLyi(:,2)
c    &             + A3(:,2,3) * rLyi(:,3) 
     &             + A3(:,2,4) * rLyi(:,4)
     &             + A3(:,2,5) * rLyi(:,5)
     &             + ri(:,12)
        ri(:,13) =
     &               A3(:,3,1) * rLyi(:,1) 
c    &             + A3(:,3,2) * rLyi(:,2)
     &             + A3(:,3,3) * rLyi(:,3) 
     &             + A3(:,3,4) * rLyi(:,4)
     &             + A3(:,3,5) * rLyi(:,5)
     &             + ri(:,13)
        ri(:,14) =
     &               A3(:,4,1) * rLyi(:,1) 
c    &             + A3(:,4,2) * rLyi(:,2)
c    &             + A3(:,4,3) * rLyi(:,3) 
     &             + A3(:,4,4) * rLyi(:,4)
     &             + A3(:,4,5) * rLyi(:,5)
     &             + ri(:,14)
        ri(:,15) =
     &               A3(:,5,1) * rLyi(:,1) 
     &             + A3(:,5,2) * rLyi(:,2)
     &             + A3(:,5,3) * rLyi(:,3) 
     &             + A3(:,5,4) * rLyi(:,4)
     &             + A3(:,5,5) * rLyi(:,5)
     &             + ri(:,15)
c
       endif ! for ires=3 or 1

c     
c.... ---------------------------->  LHS  <----------------------------
c
       if (lhs .eq. 1) then
c
c.... calculate (Atau) <-- (A_1 tau) (Recall that we are using a 
c                                     diagonal tau here)          
c
          
          if (itau.lt.10) then

             do i = 1, nflow
                Atau(:,i,1) = A1(:,i,1)*tau(:,1)
                Atau(:,i,2) = A1(:,i,2)*tau(:,2)
                Atau(:,i,3) = A1(:,i,3)*tau(:,2)
                Atau(:,i,4) = A1(:,i,4)*tau(:,2)
                Atau(:,i,5) = A1(:,i,5)*tau(:,3)
             enddo
             
          else

             Atau = zero
             do i = 1, nflow
                do j = 1, nflow
                   do k = 1, nflow
                      Atau(:,i,j) =Atau(:,i,j) + A1(:,i,k)*PTau(:,k,j)
                   enddo
                enddo
             enddo
             
          endif
c     
c.... calculate (A_1 tau A_0) (for L.S. time term of EGmass)
c
       do j = 1, nflow
          do i = 1, nflow
             A1tauA0(:,i,j) = 
     &            Atau(:,i,1)*A0(:,1,j) +
     &            Atau(:,i,2)*A0(:,2,j) +
     &            Atau(:,i,3)*A0(:,3,j) +
     &            Atau(:,i,4)*A0(:,4,j) +
     &            Atau(:,i,5)*A0(:,5,j)
          enddo
       enddo
c
c.... add (A_1 tau A_1) to stiff [1,1]
c
       do j = 1, nflow
          do i = 1, nflow
             stiff(:,i,j) = stiff(:,i,j) + (
     &              Atau(:,i,1)*A1(:,1,j)
     &            + Atau(:,i,2)*A1(:,2,j)
     &            + Atau(:,i,3)*A1(:,3,j)
     &            + Atau(:,i,4)*A1(:,4,j)
     &            + Atau(:,i,5)*A1(:,5,j)
     &            )
          enddo
       enddo
c
c.... add (A_1 tau A_2) to stiff [1,2]
c
       do j = 1, nflow
          do i = 1, nflow
             stiff(:,i,j+5) = stiff(:,i,j+5) + (
     &              Atau(:,i,1)*A2(:,1,j)
     &            + Atau(:,i,2)*A2(:,2,j)
     &            + Atau(:,i,3)*A2(:,3,j)
     &            + Atau(:,i,4)*A2(:,4,j)
     &            + Atau(:,i,5)*A2(:,5,j)
     &            )
          enddo
       enddo
c
c.... add (A_1 tau A_3) to stiff [1,3]
c
       do j = 1, nflow
          do i = 1, nflow
             stiff(:,i,j+10) = stiff(:,i,j+10) + (
     &              Atau(:,i,1)*A3(:,1,j)
     &            + Atau(:,i,2)*A3(:,2,j)
     &            + Atau(:,i,3)*A3(:,3,j)
     &            + Atau(:,i,4)*A3(:,4,j)
     &            + Atau(:,i,5)*A3(:,5,j)
     &            )
          enddo
       enddo
c
c.... calculate (Atau) <-- (A_2 tau) (Recall that we are using a 
c                                     diagonal tau here)          
c     
       if (itau.lt.10) then

          do i = 1, nflow
             Atau(:,i,1) = A2(:,i,1)*tau(:,1)
             Atau(:,i,2) = A2(:,i,2)*tau(:,2)
             Atau(:,i,3) = A2(:,i,3)*tau(:,2)
             Atau(:,i,4) = A2(:,i,4)*tau(:,2)
             Atau(:,i,5) = A2(:,i,5)*tau(:,3)
          enddo
          
       else
          Atau = zero
          do i = 1, nflow
             do j = 1, nflow
                do k = 1, nflow
                   Atau(:,i,j) = Atau(:,i,j) + A2(:,i,k)*PTau(:,k,j)
                enddo
             enddo
          enddo
          
       endif
c     
c.... calculate (A_2 tau A_0) (for L.S. time term of EGmass)
c
       do j = 1, nflow
          do i = 1, nflow
             A2tauA0(:,i,j) = 
     &            Atau(:,i,1)*A0(:,1,j) +
     &            Atau(:,i,2)*A0(:,2,j) +
     &            Atau(:,i,3)*A0(:,3,j) +
     &            Atau(:,i,4)*A0(:,4,j) +
     &            Atau(:,i,5)*A0(:,5,j)
          enddo
       enddo
c
c.... add (A_2 tau A_1) to stiff [2,1]
c
       do j = 1, nflow
          do i = 1, nflow
             stiff(:,i+5,j) = stiff(:,i+5,j) + (
     &              Atau(:,i,1)*A1(:,1,j)
     &            + Atau(:,i,2)*A1(:,2,j)
     &            + Atau(:,i,3)*A1(:,3,j)
     &            + Atau(:,i,4)*A1(:,4,j)
     &            + Atau(:,i,5)*A1(:,5,j)
     &            )
          enddo
       enddo
c
c.... add (A_2 tau A_2) to stiff [2,2]
c
       do j = 1, nflow
          do i = 1, nflow
             stiff(:,i+5,j+5) = stiff(:,i+5,j+5) + (
     &              Atau(:,i,1)*A2(:,1,j)
     &            + Atau(:,i,2)*A2(:,2,j)
     &            + Atau(:,i,3)*A2(:,3,j)
     &            + Atau(:,i,4)*A2(:,4,j)
     &            + Atau(:,i,5)*A2(:,5,j)
     &            )
          enddo
       enddo
c
c.... add (A_2 tau A_3) to stiff [2,3]
c
       do j = 1, nflow
          do i = 1, nflow
             stiff(:,i+5,j+10) = stiff(:,i+5,j+10) + (
     &              Atau(:,i,1)*A3(:,1,j)
     &            + Atau(:,i,2)*A3(:,2,j)
     &            + Atau(:,i,3)*A3(:,3,j)
     &            + Atau(:,i,4)*A3(:,4,j)
     &            + Atau(:,i,5)*A3(:,5,j)
     &            )
          enddo
       enddo
c
c.... calculate (Atau) <-- (A_3 tau) (Recall that we are using a 
c                                     diagonal tau here)          
c     
       if (itau.lt.10) then
          
          do i = 1, nflow
             Atau(:,i,1) = A3(:,i,1)*tau(:,1)
             Atau(:,i,2) = A3(:,i,2)*tau(:,2)
             Atau(:,i,3) = A3(:,i,3)*tau(:,2)
             Atau(:,i,4) = A3(:,i,4)*tau(:,2)
             Atau(:,i,5) = A3(:,i,5)*tau(:,3)
          enddo
       
       else
          Atau = zero
          do i = 1, nflow
             do j = 1, nflow
                do k = 1, nflow
                   Atau(:,i,j) = Atau(:,i,j) + A3(:,i,k)*PTau(:,k,j)
                enddo
             enddo
          enddo
          
       endif
c
c.... calculate (A_3 tau A_0) (for L.S. time term of EGmass)
c
       do j = 1, nflow
          do i = 1, nflow
             A3tauA0(:,i,j) = 
     &            Atau(:,i,1)*A0(:,1,j) +
     &            Atau(:,i,2)*A0(:,2,j) +
     &            Atau(:,i,3)*A0(:,3,j) +
     &            Atau(:,i,4)*A0(:,4,j) +
     &            Atau(:,i,5)*A0(:,5,j)
          enddo
       enddo
c
c.... add (A_3 tau A_1) to stiff [3,1]
c
       do j = 1, nflow
          do i = 1, nflow
             stiff(:,i+10,j) = stiff(:,i+10,j) + (
     &              Atau(:,i,1)*A1(:,1,j)
     &            + Atau(:,i,2)*A1(:,2,j)
     &            + Atau(:,i,3)*A1(:,3,j)
     &            + Atau(:,i,4)*A1(:,4,j)
     &            + Atau(:,i,5)*A1(:,5,j)
     &            )
          enddo
       enddo
c
c.... add (A_3 tau A_2) to stiff [3,2]
c
       do j = 1, nflow
          do i = 1, nflow
             stiff(:,i+10,j+5) = stiff(:,i+10,j+5) + (
     &              Atau(:,i,1)*A2(:,1,j)
     &            + Atau(:,i,2)*A2(:,2,j)
     &            + Atau(:,i,3)*A2(:,3,j)
     &            + Atau(:,i,4)*A2(:,4,j)
     &            + Atau(:,i,5)*A2(:,5,j)
     &            )
          enddo
       enddo
c
c.... add (A_3 tau A_3) to stiff [3,3]
c
       do j = 1, nflow
          do i = 1, nflow
             stiff(:,i+10,j+10) = stiff(:,i+10,j+10) + (
     &              Atau(:,i,1)*A3(:,1,j)
     &            + Atau(:,i,2)*A3(:,2,j)
     &            + Atau(:,i,3)*A3(:,3,j)
     &            + Atau(:,i,4)*A3(:,4,j)
     &            + Atau(:,i,5)*A3(:,5,j)
     &            )
          enddo
       enddo
c
c.... add least squares time term to the LHS tangent mass matrix
c
c
c.... loop through rows (nodes i)
c
       do i = 1, nshl
          i0 = nflow * (i - 1)
c
c.... first calculate (Atau) <-- (N_a,i A_i tau A_0)
c     ( use Atau to conserve space )
c
          do idof = 1, nflow
             do jdof = 1, nflow
                Atau(:,idof,jdof) =
     &               shg(:,i,1) * A1tauA0(:,idof,jdof) +
     &               shg(:,i,2) * A2tauA0(:,idof,jdof) +
     &               shg(:,i,3) * A3tauA0(:,idof,jdof)
             enddo
          enddo
c
c.... loop through column nodes, add (N_a,i A_i tau N_b) to EGmass
c
          do j = 1, nshl
             j0 = nflow * (j - 1)
c
c.... compute the factors
c
            fact = shape(:,j) * WdetJ * almi/gami/alfi*dtgl
c
c.... loop through d.o.f.'s
c
            do idof = 1, nflow
               il = i0 + idof
               
               EGmass(:,il,j0+1) = EGmass(:,il,j0+1) +
     &                             fact * Atau(:,idof,1)
               EGmass(:,il,j0+2) = EGmass(:,il,j0+2) +
     &                             fact * Atau(:,idof,2)
               EGmass(:,il,j0+3) = EGmass(:,il,j0+3) +
     &                             fact * Atau(:,idof,3)
               EGmass(:,il,j0+4) = EGmass(:,il,j0+4) +
     &                             fact * Atau(:,idof,4)
               EGmass(:,il,j0+5) = EGmass(:,il,j0+5) +
     &                             fact * Atau(:,idof,5)
            enddo
c
c.... end loop on column nodes
c
         enddo
c
c.... end loop on row nodes
c
       enddo
c
c.... end LHS computation
c      
       endif
       
       ttim(24) = ttim(24) + secs(0.0)
c
c.... return
c
        return
        end
c
c
c
        subroutine e3LSSclr  (A1t,     A2t,     A3t,    
     &                        rho,     rmu,     rTLSt,
     &                        u1,      u2,      u3, 
     &                        rLyti,   dxidx,   raLSt,
     &                        rti,     rk,      giju,
     &                        acti,    A0t,
     &                        shape,   shg,
     &                        EGmasst, stifft,  WdetJ,
     &                        srcp)
c                                                                
c----------------------------------------------------------------------
c
c This routine calculates the contribution of the least-squares 
c operator to the RHS vector and LHS tangent matrix. The temporary 
c results are put in ri.
c
c input:
c  A0t    (npro)              : A_0
c  A1t    (npro)              : A_1
c  A2t    (npro)              : A_2
c  A3t    (npro)              : A_3
c  acti   (npro)              : time-deriv. of Sclr
c  rho    (npro)              : density
c  rmu    (npro)              : molecular viscosity
c  rk     (npro)              : kinetic energy
c  u1     (npro)              : x1-velocity component
c  u2     (npro)              : x2-velocity component
c  u3     (npro)              : x3-velocity component
c  rLyti  (npro)              : least-squares residual vector
c  dxidx  (npro,nsd,nsd)      : inverse of deformation gradient
c  taut   (npro)              : stability parameter
c  rLyti  (npro)              : convective portion of least-squares
c                               residual vector 
c  divqti (npro,1)            : divergence of diffusive flux
c  shape  (npro,nshl)         : element shape functions
c  shg    (npro,nshl,nsd)     : global element shape function grads
c  WdetJ  (npro)              : weighted jacobian determinant
c  stifft (npro,nsd,nsd)      : stiffness matrix
c  EGmasst(npro,nshape,nshape): partial mass matrix
c
c output:
c  rti    (npro,nsd+1)        : partial residual
c  EGmasst(npro,nshape,nshape): partial mass matrix
c
c
c Zdenek Johan, Summer 1990. (Modified from e2ls.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997. Prim. Var. with Diag Tau
c----------------------------------------------------------------------
c
        include "common.h"

c
c  passed arrays
c
        dimension A1t(npro),                 A2t(npro),
     &            A3t(npro),                   
     &            A0t(npro),                 rho(npro),
     &            acti(npro),                rmu(npro),
     &            u1(npro),                  u2(npro),
     &            u3(npro),                  rk(npro),
     &            rLyti(npro),               dxidx(npro,nsd,nsd),
     &            taut(npro),                raLSt(npro),
     &            rti(npro,nsd+1),           rTLSt(npro),
     &            stifft(npro,3,3),          giju(npro,6),
     &            EGmasst(npro,nshape,nshape),  
     &            shape(npro,nshl),
     &            shg(npro,nshl,nsd),        WdetJ(npro),
     &            srcp(npro)
c
c local arrays
c
        dimension rLymti(npro),          Ataut(npro),
     &            A1tautA0(npro),        A2tautA0(npro),
     &            A3tautA0(npro),        fact(npro)

        ttim(24) = ttim(24) - tmr()
c
       if(ivart.lt.2) return
c
c last step to the least squares is adding the time term.  So that we
c only have to localize one vector for each Krylov vector the modified
c residual is quite different from the total residual.
c
c
c the modified residual
c
       fct1=almi/gami/alfi*dtgl
c
c  add possibility of not including time term
c
c       if(idiff.ne.-1) 
c       rLymti = rLyti + fct1*duti

       if((ires.eq.1 .or. ires .eq. 3).and. idiff.ne.-1) then

        rLyti(:) = rLyti(:) + A0t(:)*acti(:)

      endif
c
c.... subtract div(q) from the least squares term
c
c      if ((idiff >= 1).and.(ires==3 .or. ires==1)) then
c         rLyi(:) = rLyi(:) - divqti(:)
c      endif
c
c.... --------------------------->  Tau  <-----------------------------
c
c.... calculate the tau matrix
c

c
c.... we will use the same tau as used for momentum equations here
c 
        ttim(25) = ttim(25) - tmr()

          call e3tauSclr(rho,         rmu,        A0t, 
     &                   u1,          u2,         u3,
     &                   dxidx,       rLyti,      rLymti,
     &                   taut,        rk,         raLSt,
     &                   rTLSt,       giju)	

        ttim(25) = ttim(25) + tmr()
c
c Warning:  to save space I have put the tau times the modified residual 
c           in rLymi and the tau times the total residual back in rLyi
c
c
c.... ---------------------------->  RHS  <----------------------------
c
c.... calculate (A_i^T tau Ly)
c
        
c      if(ires.ne.1) then
c
c  A1 * Tau L(Y):  to be hit on left with Na,x in e3wmlt
c
c        rmti(:,1) =   A1t(:) * rLymti(:) 
c
c
c  A2 * Tau L(Y),  to be hit on left with Na,y 
c
c        rmti(:,2) = A2t(:) * rLymti(:) 
c
c
c  A3 * Tau L(Y)  to be hit on left with Na,z
c
c        rmti(:,3) = A3t(:) * rLymti(:) 
c
c      endif  !ires.ne.1

c
c  same thing for the real residual
c
       if(ires.eq.3 .or. ires .eq. 1) then  ! we need the total residual
        rti(:,1) = rti(:,1) + A1t(:) * rLyti(:) 

        rti(:,2) = rti(:,2) + A2t(:) * rLyti(:) 

        rti(:,3) = rti(:,3) + A3t(:) * rLyti(:) 

       endif ! for ires=3 or 1

c     
c.... ---------------------------->  LHS  <----------------------------
c
       if (lhs .eq. 1) then
c       
c
c.... calculate (Atau) <-- (A_1 tau)                                    
c     

          Ataut(:) = A1t(:)*taut(:)

c
c.... calculate (A_1 tau (A_0-srcp)) (for L.S. time term of EGmass)
c

             A1tautA0(:) = Ataut(:)*(a0t(:)*fct1-srcp(:))

c
c.... add (A_1 tau A_1) to stiff [1,1]
c

             stifft(:,1,1) = stifft(:,1,1) + Ataut(:)*A1t(:)
c            stifft(:,1,1) = Ataut(:)*A1t(:)
c
c.... add (A_1 tau A_2) to stiff [1,2]
c

             stifft(:,1,2) = stifft(:,1,2) + Ataut(:)*A2t(:)
c            stifft(:,1,2) =  Ataut(:)*A2t(:)
c
c.... add (A_1 tau A_3) to stiff [1,3]
c

             stifft(:,1,3) = stifft(:,1,3) + Ataut(:)*A3t(:)
c            stifft(:,1,3) =  Ataut(:)*A3t(:)
c
c.... calculate (Atau) <-- (A_2 tau)
c     

          Ataut(:) = A2t(:)*taut(:)

c
c.... calculate (A_2 tau (A_0-srcp)) (for L.S. time term of EGmass)
c

             A2tautA0(:) = Ataut(:)*(a0t(:)*fct1-srcp(:))

c
c.... add (A_2 tau A_1) to stiff [2,1]
c

             stifft(:,2,1) = stifft(:,1,2)
c
c.... add (A_2 tau A_2) to stiff [2,2]
c

             stifft(:,2,2) = stifft(:,2,2) + Ataut(:)*A2t(:)

c
c.... add (A_2 tau A_3) to stiff [2,3]
c

             stifft(:,2,3) = stifft(:,2,3) + Ataut(:)*A3t(:)

c
c.... calculate (Atau) <-- (A_3 tau)
c     

          Ataut(:) = A3t(:)*taut(:)

c
c.... calculate (A_3 tau (A_0-srcp)) (for L.S. time term of EGmass)
c

             A3tautA0(:) = Ataut(:)*(a0t(:)*fct1-srcp(:))

c
c.... add (A_3 tau A_1) to stiff [3,1]
c

             stifft(:,3,1) = stifft(:,1,3)

c
c.... add (A_3 tau A_2) to stiff [3,2]
c

             stifft(:,3,2) = stifft(:,2,3)

c
c.... add (A_3 tau A_3) to stiff [3,3]
c

             stifft(:,3,3) = stifft(:,3,3) + Ataut(:)*A3t(:)

c
c.... add least squares time term to the LHS tangent mass matrix
c
c
c.... loop through rows (nodes i)
c
       do ia = 1, nshl
c
c.... first calculate (Atau) <-- (N_a,i A_i tau A_0)
c     ( use Atau to conserve space )
c

                Ataut(:) =
     &               shg(:,ia,1) * A1tautA0(:) +
     &               shg(:,ia,2) * A2tautA0(:) +
     &               shg(:,ia,3) * A3tautA0(:)

c
c.... loop through column nodes, add (N_a,i A_i tau N_b) to EGmass
c
          do jb = 1, nshl

            fact = shape(:,jb) * WdetJ

            EGmasst(:,ia,jb) = EGmasst(:,ia,jb) + fact * Ataut(:)

c
c.... end loop on column nodes
c

          enddo
c
c.... end loop on row nodes
c
       enddo
c
c.... end LHS computation
c      
       endif
       
       ttim(24) = ttim(24) + tmr()
c
c.... return
c
        return
        end

