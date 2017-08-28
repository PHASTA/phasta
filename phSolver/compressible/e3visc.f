        subroutine e3visc (g1yi,    g2yi,    g3yi,
     &                     dxidx,
     &                     rmu,     rlm,     rlm2mu,
     &                     u1,      u2,      u3,      
     &                     ri,      rmi,     stiff,
     &                     con,     rlsli,   compK, T)
c
c----------------------------------------------------------------------
c
c This routine calculates the contribution of viscous and heat fluxes
c to both RHS and LHS.
c
c input:
c  g1yi   (npro,nflow)         : grad-y in direction 1
c  g2yi   (npro,nflow)         : grad-y in direction 2
c  g3yi   (npro,nflow)         : grad-y in direction 3
c  u1     (npro)              : x1-velocity component
c  u2     (npro)              : x2-velocity component
c  u3     (npro)              : x3-velocity component
c  rmu    (npro)              : viscosity
c  rlm    (npro)              : Lambda
c  rlm2mu (npro)              : Lambda + 2 Mu
c  con    (npro)              : Conductivity
c  cp     (npro)              : specific heat at constant pressure
c  rlsli  (npro,6)              : Resolved Reynold's stresses
c output:
c  ri     (npro,nflow*(nsd+1)) : partial residual
c  rmi    (npro,nflow*(nsd+1)) : partial modified residual
c  stiff  (npro,nsd*nflow,nsd*nflow) : stiffness matrix
c  compK (npro,10)             : K_ij in (eq.134) A new ... III 
c
c
c Zdenek Johan, Summer 1990. (Modified from e2visc.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997 Primitive Variables
c----------------------------------------------------------------------
c
      include "common.h"
c
c     passed arrays
c
      dimension g1yi(npro,nflow),           g2yi(npro,nflow),
     &     g3yi(npro,nflow),   
     &     Sclr(npro),                dxidx(npro,nsd,nsd),
     &     u1(npro),                  u2(npro),
     &     u3(npro),                  rho(npro),
     &     ri(npro,nflow*(nsd+1)),     rmi(npro,nflow*(nsd+1))
c
c
      dimension rmu(npro),                 rlm(npro),
     &     rlm2mu(npro),                   con(npro),
     &     stiff(npro,3*nflow,3*nflow),
     &     rlsli(npro,6),                  compK(npro,10),
     &     f1(npro), f2(npro), f3(npro), f4(npro), 
     &     f5(npro), f6(npro),  T(npro), rk(npro) 

      ttim(23) = ttim(23) - secs(0.0)
c
c... dynamic model is now being accounted for in getdiff.f
c
c
c.... --------------------->  Diffusivity Matrix  <-------------------
c

      if (lhs .eq. 1) then
c
c.... K11
c
         stiff(:, 2, 2) = rlm2mu
         stiff(:, 3, 3) = rmu
         stiff(:, 4, 4) = rmu
         stiff(:, 5, 2) = rlm2mu * u1
         stiff(:, 5, 3) = rmu    * u2
         stiff(:, 5, 4) = rmu    * u3
         stiff(:, 5, 5) = con
c     
c.... K12
c     
         stiff(:, 2, 8) = rlm
         stiff(:, 3, 7) = rmu
         stiff(:, 5, 7) = rmu    * u2
         stiff(:, 5, 8) = rlm    * u1
c     
c.... K13
c     
         stiff(:, 2,14) = rlm
         stiff(:, 4,12) = rmu
         stiff(:, 5,12) = rmu    * u3
         stiff(:, 5,14) = rlm    * u1
           
c     
c.... K21
c     
         stiff(:, 7, 3) = rmu
         stiff(:, 8, 2) = rlm
         stiff(:,10, 2) = rlm    * u2
         stiff(:,10, 3) = rmu    * u1
           
c     
c.... K22
c     
         stiff(:, 7, 7) = rmu
         stiff(:, 8, 8) = rlm2mu
         stiff(:, 9, 9) = rmu
         stiff(:,10, 7) = rmu    * u1
         stiff(:,10, 8) = rlm2mu * u2
         stiff(:,10, 9) = rmu    * u3
         stiff(:,10,10) = con
c     
c.... K23
c     
         stiff(:, 8,14) = rlm
         stiff(:, 9,13) = rmu
         stiff(:,10,13) = rmu    * u3
         stiff(:,10,14) = rlm    * u2
c     
c.... K31
c     
         stiff(:,12, 4) = rmu
         stiff(:,14, 2) = rlm
         stiff(:,15, 2) = rlm    * u3
         stiff(:,15, 4) = rmu    * u1
c     
c.... K32
c     
         stiff(:,13, 9) = rmu
         stiff(:,14, 8) = rlm
         stiff(:,15, 8) = rlm    * u3
         stiff(:,15, 9) = rmu    * u2
c     
c.... K33
c     
         stiff(:,12,12) = rmu
         stiff(:,13,13) = rmu
         stiff(:,14,14) = rlm2mu
         stiff(:,15,12) = rmu    * u1
         stiff(:,15,13) = rmu    * u2
         stiff(:,15,14) = rlm2mu * u3
         stiff(:,15,15) = con

      endif

      if (itau .ge. 10) then     ! non-diag tau, buld compK
      
c     
c.... calculate element factors
c     
         f1 = dxidx(:,1,1) * dxidx(:,1,1) +
     &           dxidx(:,2,1) * dxidx(:,2,1) +
     &           dxidx(:,3,1) * dxidx(:,3,1)
         f2 = dxidx(:,1,1) * dxidx(:,1,2) +
     &           dxidx(:,2,1) * dxidx(:,2,2) +
     &           dxidx(:,3,1) * dxidx(:,3,2)
         f3 = dxidx(:,1,2) * dxidx(:,1,2) +
     &           dxidx(:,2,2) * dxidx(:,2,2) +
     &           dxidx(:,3,2) * dxidx(:,3,2)
         f4 = dxidx(:,1,1) * dxidx(:,1,3) +
     &           dxidx(:,2,1) * dxidx(:,2,3) +
     &           dxidx(:,3,1) * dxidx(:,3,3)
         f5 = dxidx(:,1,2) * dxidx(:,1,3) +
     &           dxidx(:,2,2) * dxidx(:,2,3) +
     &           dxidx(:,3,2) * dxidx(:,3,3)
         f6 = dxidx(:,1,3) * dxidx(:,1,3) +
     &           dxidx(:,2,3) * dxidx(:,2,3) +
     &           dxidx(:,3,3) * dxidx(:,3,3)
c     
c.... correction for tetrahedra (invariance w.r.t triangular coord)
c     
         if (lcsyst .eq. 1) then
            f1 = ( f1 + (dxidx(:,1,1) +
     &           dxidx(:,2,1) + dxidx(:,3,1))**2 ) * pt39
            f2 = ( f2 + (dxidx(:,1,1) +
     &              dxidx(:,2,1) + dxidx(:,3,1)) *
     &              (dxidx(:,1,2) +
     &              dxidx(:,2,2) + dxidx(:,3,2)) ) * pt39
            f3 = ( f3 + (dxidx(:,1,2) +
     &              dxidx(:,2,2) + dxidx(:,3,2))**2 ) * pt39
            f4 = ( f4 + (dxidx(:,1,1) +
     &              dxidx(:,2,1) + dxidx(:,3,1)) *
     &              (dxidx(:,1,3) +
     &              dxidx(:,2,3) + dxidx(:,3,3)) ) * pt39
            f5 = ( f5 + (dxidx(:,1,2) +
     &              dxidx(:,2,2) + dxidx(:,3,2)) *
     &              (dxidx(:,1,3) +
     &              dxidx(:,2,3) + dxidx(:,3,3)) ) * pt39
            f6 = ( f6 + (dxidx(:,1,3) +
     &              dxidx(:,2,3) + dxidx(:,3,3))**2 ) * pt39
c     
       !      flops = flops + 36*npro
         endif
c     
c.... correction for wedges
c     
         if ((lcsyst .eq. 3) .or. (lcsyst .eq. 4)) then
            f1 = ( dxidx(:,1,1) * dxidx(:,1,1) +
     &           dxidx(:,2,1) * dxidx(:,2,1) +
     &           ( dxidx(:,1,1) + dxidx(:,2,1) )**2 ) * pt57
     &           + dxidx(:,3,1) * dxidx(:,3,1)
            f2 = ( dxidx(:,1,1) * dxidx(:,1,2) +
     &           dxidx(:,2,1) * dxidx(:,2,2) +
     &           ( dxidx(:,1,1) + dxidx(:,2,1) ) *
     &           ( dxidx(:,1,2) + dxidx(:,2,2) ) ) * pt57
     &           + dxidx(:,3,1) * dxidx(:,3,2)
            f3 = ( dxidx(:,1,2) * dxidx(:,1,2) +
     &           dxidx(:,2,2) * dxidx(:,2,2) +
     &           ( dxidx(:,1,2) + dxidx(:,2,2) )**2 ) * pt57
     &           + dxidx(:,3,2) * dxidx(:,3,2)
            f4 = ( dxidx(:,1,1) * dxidx(:,1,3) +
     &           dxidx(:,2,1) * dxidx(:,2,3) +
     &           ( dxidx(:,1,1) + dxidx(:,2,1) ) *
     &           ( dxidx(:,1,3) + dxidx(:,2,3) ) ) * pt57
     &           + dxidx(:,3,1) * dxidx(:,3,3)
            f5 = ( dxidx(:,1,2) * dxidx(:,1,3) +
     &           dxidx(:,2,2) * dxidx(:,2,3) +
     &           ( dxidx(:,1,2) + dxidx(:,2,2) ) *
     &           ( dxidx(:,1,3) + dxidx(:,2,3) ) ) * pt57
     &           + dxidx(:,3,2) * dxidx(:,3,3)
            f6 = ( dxidx(:,1,3) * dxidx(:,1,3) +
     &           dxidx(:,2,3) * dxidx(:,2,3) +
     &           ( dxidx(:,1,3) + dxidx(:,2,3) )**2 ) * pt57
     &           + dxidx(:,3,3) * dxidx(:,3,3)
c     
       !      flops = flops + 51*npro
         endif
c     
c.... calculate compact K matrix in local parent coordinates
c.... equation 134 in "A new ... III" only w/ K^tilde_jj. Might need
c.... complete Kij.
      
         compK(:, 1) = f1 * T * rlm2mu + f3 * T * rmu
     &        + f6 * T * rmu
c     
         compK(:, 2) = f2 * T * (rlm + rmu)
         compK(:, 3) = f1 * T * rmu + f3 * T * rlm2mu
     &        + f6 * T * rmu
c     
         compK(:, 4) = f4 * T * (rlm + rmu)
         compK(:, 5) = f5 * T * (rlm + rmu)
         compK(:, 6) = f1 * T * rmu + f3 * T * rmu
     &        + f6 * T * rlm2mu
c     
         compK(:, 7) = f1 * T * rlm2mu  * u1 + f2 * T * (rlm + rmu) * u2
     &        + f3 * T * rmu * u1 + f4 * T * (rlm + rmu) * u3
     &        + f6 * T * rmu * u1
         compK(:, 8) = f1 * T * rmu * u2 + f2 * T * (rlm + rmu) * u1
     &        + f3 * T * rlm2mu  * u2 + f5 * T * (rlm + rmu) * u3
     &        + f6 * T * rmu * u2
         compK(:, 9) = f1 * T * rmu * u3 + f3 * T * rmu * u3
     &        + f4 * T * (rlm + rmu) * u1 + f5 * T * (rlm + rmu) * u2
     &        + f6 * T * rlm2mu  * u3

         rk=pt5*(u1**2+u2**2+u3**2)
         
         compK(:,10) = f1 * T * (con    * T + two * rmu * rk + (rlm +
     &        rmu) * u1**2) + f2 * T * (rlm + rmu) * two * u1 * u2 
     &        + f3 * T * (con    * T + two * rmu * rk + (rlm + rmu) *
     &        u2**2) + f4 * T * (rlm + rmu) * two * u1 * u3 
     &        + f5 * T * (rlm + rmu) * two * u2 * u3 + f6 * T * (con
     &        * T + two * rmu * rk + (rlm + rmu) * u3**2)  
c     
c.... flop count
c
    !      flops = flops + 86*npro
c     
c.... end of GLS
c     
      
      endif
c     
c.... --------------------------->  RHS  <-----------------------------
c     
c.... compute diffusive fluxes and add them to ri and rmi

c
c.... diffusive flux in x1-direction
c
c       rmi(:,1) = zero ! already initialized
        rmi(:,2) =  rlm2mu      * g1yi(:,2) 
     &               +      rlm * g2yi(:,3) 
     &               +      rlm * g3yi(:,4)
     &               -      rlsli(:,1)
        rmi(:,3) =  rmu         * g1yi(:,3) 
     &               +      rmu * g2yi(:,2) 
     &               -      rlsli(:,4)
        rmi(:,4) =  rmu         * g1yi(:,4)
     &               +      rmu * g3yi(:,2)
     &               -      rlsli(:,5)
        rmi(:,5) =  rlm2mu * u1 * g1yi(:,2) + rmu * u2 * g1yi(:,3)
     &                                   +    rmu * u3 * g1yi(:,4)
     &               + rmu * u2 * g2yi(:,2) + rlm * u1 * g2yi(:,3)
     &               + rmu * u3 * g3yi(:,2) + rlm * u1 * g3yi(:,4)
     &               + con      * g1yi(:,5)

c
      ri (:,2:5) = ri (:,2:5) + rmi(:,2:5)
c     rmi(:,2:5) = rmi(:,2:5) + qdi(:,2:5)
c
c!      flops = flops + 74*npro
c
c.... diffusive flux in x2-direction
c
c       rmi(:, 6) = zero
        rmi(:, 7) =       rmu * g1yi(:,3) 
     &             +      rmu * g2yi(:,2)
     &             -      rlsli(:,4)
        rmi(:, 8) =       rlm * g1yi(:,2)
     &             +   rlm2mu * g2yi(:,3)
     &             +      rlm * g3yi(:,4)
     &             -      rlsli(:,2)
        rmi(:, 9) =       rmu * g2yi(:,4)
     &             +      rmu * g3yi(:,3)
     &             -      rlsli(:,6)
        rmi(:,10) =  rlm * u2 * g1yi(:,2) +    rmu * u1 * g1yi(:,3)
     &             + rmu * u1 * g2yi(:,2) + rlm2mu * u2 * g2yi(:,3)
     &             + rmu * u3 * g2yi(:,4)
     &             + rmu * u3 * g3yi(:,3) +    rlm * u2 * g3yi(:,4)
     &             +    con   * g2yi(:,5)
c
      ri (:,7:10) = ri (:,7:10) + rmi(:,7:10)
c     rmi(:,7:10) = rmi(:,7:10) + qdi(:,2:5)
c
c!      flops = flops + 74*npro
c
c.... diffusive flux in x3-direction
c
c       rmi(:,11) = zero
        rmi(:,12) =       rmu * g1yi(:,4)
     &             +      rmu * g3yi(:,2)
     &             -      rlsli(:,5)
        rmi(:,13) =       rmu * g2yi(:,4)
     &             +      rmu * g3yi(:,3)
     &             -      rlsli(:,6)
        rmi(:,14) =       rlm * g1yi(:,2)
     &             +      rlm * g2yi(:,3)
     &             +   rlm2mu * g3yi(:,4)
     &             -   rlsli(:,3)
        rmi(:,15) =     rlm * u3 * g1yi(:,2) + rmu * u1 * g1yi(:,4)
     &             +    rlm * u3 * g2yi(:,3) + rmu * u2 * g2yi(:,4)
     &             +    rmu * u1 * g3yi(:,2) + rmu * u2 * g3yi(:,3)
     &             + rlm2mu * u3 * g3yi(:,4)
     &             +    con      * g3yi(:,5) 
c
      ri (:,12:15) = ri (:,12:15) + rmi(:,12:15)
c
c!      flops = flops + 74*npro
c
c  stiff for preconditioner has been eliminated
c  all preconditioner stuff is in e3bdg.f
c

      ttim(23) = ttim(23) + secs(0.0)
      
c
c.... return
c
        return
        end
c     
c     
c     
      subroutine e3viscSclr (g1yti,    g2yti,    g3yti,
     &     rmu,      Sclr,     rho,
     &     rti,      rmti,     stifft)
c     
c----------------------------------------------------------------------
c     
c     This routine calculates the contribution of viscous and heat fluxes
c     to both RHS and LHS.
c     
c     input:
c     g1yti  (npro)              : grad-y in direction 1
c     g2yti  (npro)              : grad-y in direction 2
c     g3yti  (npro)              : grad-y in direction 3
c     rmu    (npro)              : viscosity
c     Sclr   (npro)              : scalar variable
c     
c     output:
c     rti     (npro,nsd+1)       : partial residual
c     rmti    (npro,nsd+1)       : partial modified residual
c     stifft  (npro,nsd,nsd)     : stiffness matrix
c     
c     
c     
c     Zdenek Johan, Summer 1990. (Modified from e2visc.f)
c     Zdenek Johan, Winter 1991. (Fortran 90)
c     Kenneth Jansen, Winter 1997 Primitive Variables
c----------------------------------------------------------------------
c     
      use turbSA                ! for saSigma
      include "common.h"
c     
c     passed arrays
c     
      dimension g1yti(npro),           g2yti(npro),
     &     g3yti(npro),           rmu(npro),
     &     Sclr(npro),            rho(npro),
     &     rti(npro,nsd+1),       rmti(npro,nsd+1),
     &     stifft(npro,3,3)

      ttim(23) = ttim(23) - tmr()

      if ((ilset.ne.zero) .and. (isclr.lt.3)) return 
c     
c.... --------------------------->  RHS  <-----------------------------
c     
c.... --------------------->  Diffusivity Matrix  <-------------------
c     
c      if (lhs .eq. 1) then

         stifft = zero
c     
c.... K11
c     
         stifft(:,1,1)=rmu
         if (iRANS .lt. 0) then
            stifft(:,1,1)=saSigmaInv*((rmu/rho)+max(zero,Sclr))
            if(iconvsclr.eq.1) stifft(:,1,1)=rho*stifft(:,1,1)
         endif
c     
c.... K22
c     
         stifft(:,2,2)=stifft(:,1,1)
c     
c.... K33
c     
         stifft(:,3,3)=stifft(:,1,1)

c      endif
c     
c.... --------------------------->  RHS  <-----------------------------
c     
c.... compute diffusive fluxes and add them to ri and rmi
c     
c.... diffusive fluxes
c     
      rmti(:,1) = stifft(:,1,1) * g1yti(:) 
      rmti(:,2) = stifft(:,2,2) * g2yti(:) 
      rmti(:,3) = stifft(:,3,3) * g3yti(:) 

      rti (:,:) = rti (:,:) + rmti(:,:)
c     rmi(:,2:5) = rmi(:,2:5) + qdi(:,2:5)
c     
c!      flops = flops + 74*npro
c     
      ttim(23) = ttim(23) + tmr()

c     
c.... return
c     
      return
      end
