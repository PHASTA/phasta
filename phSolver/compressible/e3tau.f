      subroutine e3tau  (rho,    cp,     rmu, 	
     &                     u1,     u2,     u3,
     &                     con,    dxidx,  rLyi,
     &                     rLymi,  tau,    rk, 
     &                     giju,   rTLS,   raLS,
     &                     A0inv,  dVdY,   cv)
c
c----------------------------------------------------------------------
c
c This routine computes the diagonal Tau for least-squares operator.  
c We receive the regular residual L operator and a
c modified residual L operator, calculate tau, and return values for
c tau and tau times these operators to maintain the format of past
c ENSA codes
c
c input:
c  rho    (npro)           : density
c  T      (npro)           : temperature
c  cp     (npro)           : specific heat at constant pressure
c  u1     (npro)           : x1-velocity component
c  u2     (npro)           : x2-velocity component
c  u3     (npro)           : x3-velocity component
c  dxidx  (npro,nsd,nsd)   : inverse of deformation gradient
c  rLyi   (npro,nflow)      : least-squares residual vector
c  rLymi   (npro,nflow)     : modified least-squares residual vector
c
c output:
c  rLyi   (npro,nflow)      : least-squares residual vector
c  rLymi   (npro,nflow)     : modified least-squares residual vector
c  tau    (npro,3)         : 3 taus
c
c
c Zdenek Johan, Summer 1990.  (Modified from e2tau.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      include "common.h"
c
      dimension rho(npro),                 con(npro), 
     &            cp(npro),                  u1(npro),
     &            u2(npro),                  u3(npro),
     &            dxidx(npro,nsd,nsd),       rk(npro),
     &            tau(npro,5),               rLyi(npro,nflow),
     &            rLymi(npro,nflow),         dVdY(npro,15), 
     &            rTLS(npro),                raLS(npro),
     &            rLyitemp(npro,nflow),      detgijI(npro)
c
      dimension   rmu(npro),	 cv(npro),
     &		  gijd(npro,6),  uh1(npro),
     &		  fact(npro),	 h2o2u(npro),   giju(npro,6),
     &            A0inv(npro,15),gijdu(npro,6)
c
      call e3gijd( dxidx, gijd )
c
c  next form the diffusive length scale |u| h_1 = 2 ( ui (gijd)-1 uj)^{1/2}
c
c   dividing factor for the inverse of gijd
c
         fact = gijd(:,1) * gijd(:,3) * gijd(:,6)
     &        - gijd(:,1) * gijd(:,5) * gijd(:,5)
     &        - gijd(:,3) * gijd(:,4) * gijd(:,4)
     &        - gijd(:,6) * gijd(:,2) * gijd(:,2)
     &        + gijd(:,2) * gijd(:,4) * gijd(:,5) * two
c
         uh1=    u1*u1*(gijd(:,3)*gijd(:,6)-gijd(:,5)*gijd(:,5))
     &         + u2*u2*(gijd(:,1)*gijd(:,6)-gijd(:,4)*gijd(:,4))
     &         + u3*u3*(gijd(:,1)*gijd(:,3)-gijd(:,2)*gijd(:,2))
     &   + two *(u1*u2*(gijd(:,4)*gijd(:,5)-gijd(:,2)*gijd(:,6))
     &         + u1*u3*(gijd(:,2)*gijd(:,5)-gijd(:,4)*gijd(:,3))
     &         + u2*u3*(gijd(:,4)*gijd(:,2)-gijd(:,1)*gijd(:,5)))
c
c   at this point we have (u h1)^2 * inverse coefficient^2 / 4
c
c                                    ^ fact
c
      uh1= two * sqrt(uh1/fact)
c
c  next form the advective length scale |u|/h_2 = 2 ( ui (gijd) uj)^{1/2}
c
      h2o2u =   u1*u1*gijd(:,1)
     &     + u2*u2*gijd(:,3)
     &     + u3*u3*gijd(:,6)
     &     +(u1*u2*gijd(:,2)
     &     + u1*u3*gijd(:,4)
     &     + u2*u3*gijd(:,5))*two  + 1.0e-15 !FIX FOR INVALID MESHES
c
c  at this point we have (2 u / h_2)^2
c

c       call tnanqe(h2o2u,1,"riaconv ")

      h2o2u = one / sqrt(h2o2u) ! this flips it over leaves it h_2/(2u)
c  
c.... diffusive corrections

      if(itau.eq.1) then        ! tau proposed by  for nearly incompressible
c                                 flows by Guillermo Hauke
c     
c.... cell Reynold number
c     
         fact=pt5*rho*uh1/rmu

c     
c.... continuity tau
c     
         tau(:,1)=pt5*uh1*min(one,fact)*taucfct
c     
c...  momentum tau
c     
         dts=one/(Dtgl*dtsfct)
         tau(:,2)=min(dts,h2o2u)
         tau(:,2)=tau(:,2)/rho
c     
c.... energy tau   cv=cp/gamma  assumed
c   
c         tau(:,3)=gamma*tau(:,2)/cp
          tau(:,3)=tau(:,2)/cv
c     
c.... diffusive corrections
c     
         if (ipord == 1) then
            celt = pt66
         else if (ipord == 2) then
            celt = pt33
c            celt = pt33*0.04762
         else if (ipord == 3) then
            celt = pt33         !.02  just a guess...
         else if (ipord >= 4) then
            celt = .008         ! yet another guess...
         endif
c     
c          fact=h2o2u*h2o2u*rk*pt66/rmu
         fact=h2o2u*h2o2u*rk*celt/rmu
c
         tau(:,2)=min(tau(:,2),fact)
         tau(:,3)=min(tau(:,3),fact*rmu/con)*temper
c     
      else if(itau.eq.0)  then  ! tau proposed by Farzin and Shakib
c     
c...  momentum tau
c     

c     
c...  higher order element diffusive correction
c
         if (ipord == 1) then
            fff = 36.0d0
         else if (ipord == 2) then
            fff = 60.0d0
c     fff = 36.0d0
         else if (ipord == 3) then
            fff = 128.0d0
c     fff = 144.0d0
      endif           
         dts = dtsfct*Dtgl
         tau(:,2)=rho**2*((two*dts)**2
     &        + u1*(u1*gijd(:,1) + two*(u2*gijd(:,2)+u3*gijd(:,4)))
     &        + u2*(u2*gijd(:,3) + two*u3*gijd(:,5))
     &        + u3*u3*gijd(:,6))
     &        +fff*rmu**2*(gijd(:,1)**2 + gijd(:,3)**2 + gijd(:,6)**2 +
     &        two*(gijd(:,2)**2 + gijd(:,4)**2 + gijd(:,5)**2))
         fact=sqrt(tau(:,2))
         tau(:,1)=pt125*fact/(gijd(:,1)+gijd(:,3)+gijd(:,6))*taucfct
         tau(:,2)=one/fact
c     
c.... energy tau   cv=cp/gamma  assumed
c     
         tau(:,3)=tau(:,2)/cv*temper
c
      endif
c     
c...  finally multiply this tau matrix times the
c     two residual vectors
c
c ... calculate (tau Ly) --> (rLyi)
c ... storing rLyi for the DC term
        if(iDC .ne. 0) rLyitemp=rLyi

      if(ires.eq.3 .or. ires .eq. 1) then
         rLyi(:,1) = tau(:,1) * rLyi(:,1) 
         rLyi(:,2) = tau(:,2) * rLyi(:,2)
         rLyi(:,3) = tau(:,2) * rLyi(:,3)
         rLyi(:,4) = tau(:,2) * rLyi(:,4)
         rLyi(:,5) = tau(:,3) * rLyi(:,5)
      endif
c
      if(iDC .ne. 0) then
c..... calculation of rTLS & raLS for DC term
c
c..... calculation of (Ly - S).tau^tilda*(Ly - S) 
c 
         rTLS = rLYItemp(:,1)       * (rLyi(:,1)*dVdY(:,1)
     &        + dVdY(:,2)*rLyi(:,2) + dVdY(:,4)*rLyi(:,3)
     &        + rLyi(:,4)*dVdY(:,7) + dVdY(:,11)*rLyi(:,5))
     &        + rLyitemp(:,2)       * (rLyi(:,2)*dVdY(:,3)
     &        + rLyi(:,3)*dVdY(:,5) + dVdY(:,8)*rLyi(:,4)
     &        + rLyi(:,5)*dVdY(:,12))
     &        + rLyitemp(:,3)       * (rLyi(:,3)*dVdY(:,6)
     &        + dVdY(:,9)*rLyi(:,4) + dVdY(:,13)*rLyi(:,5))
     &        + rLyitemp(:,4)       * (rLyi(:,4)*dVdY(:,10)
     &        + dVdY(:,14)*rLyi(:,5))
     &        + rLyitemp(:,5)       * (dVdY(:,15)*rLyi(:,5))

c
c...... calculation of (Ly - S).A0inv*(Ly - S)
c
           raLS = two*rLyitemp(:,4)*rLyitemp(:,5)*A0inv(:,15)
     &          + two*rLyitemp(:,3)*rLyitemp(:,5)*A0inv(:,14)
     &          + two*rLyitemp(:,1)*rLyitemp(:,2)*A0inv( :,6)
     &          + two*rLyitemp(:,2)*rLyitemp(:,3)*A0inv(:,10)
     &          + two*rLyitemp(:,2)*rLyitemp(:,4)*A0inv(:,11)
     &          + two*rLyitemp(:,1)*rLyitemp(:,3)*A0inv( :,7)
     &          + two*rLyitemp(:,3)*rLyitemp(:,4)*A0inv(:,13)
     &          + two*rLyitemp(:,2)*rLyitemp(:,5)*A0inv(:,12)
     &          + two*rLyitemp(:,1)*rLyitemp(:,4)*A0inv( :,8)
     &          + two*rLyitemp(:,1)*rLyitemp(:,5)*A0inv( :,9)
     &          + rLyitemp(:,1)**2*A0inv(:,1)
     &          + rLyitemp(:,2)**2*A0inv(:,2)
     &          + rLyitemp(:,3)**2*A0inv(:,3)
     &          + rLyitemp(:,4)**2*A0inv(:,4)
     &          + rLyitemp(:,5)**2*A0inv(:,5)
c
c.....****************calculation of giju for DC term***************
c     
c.... for the notation of different numbering
c     
           gijdu(:,1)=gijd(:,1)
           gijdu(:,2)=gijd(:,3)
           gijdu(:,3)=gijd(:,6)
           gijdu(:,4)=gijd(:,2)
           gijdu(:,5)=gijd(:,4)
           gijdu(:,6)=gijd(:,5)
c     
c     
           detgijI = one/(
     &          gijdu(:,1) * gijdu(:,2) * gijdu(:,3)
     &          - gijdu(:,1) * gijdu(:,6) * gijdu(:,6)
     &          - gijdu(:,4) * gijdu(:,4) * gijdu(:,3)
     &          + gijdu(:,4) * gijdu(:,5) * gijdu(:,6) * two
     &          - gijdu(:,5) * gijdu(:,5) * gijdu(:,2) 
     &          )
           giju(:,1) = detgijI * (gijdu(:,2)*gijdu(:,3) 
     &               - gijdu(:,6)**2)
           giju(:,2) = detgijI * (gijdu(:,1)*gijdu(:,3) 
     &               - gijdu(:,5)**2)
           giju(:,3) = detgijI * (gijdu(:,1)*gijdu(:,2)
     &               - gijdu(:,4)**2)
           giju(:,4) = detgijI * (gijdu(:,5)*gijdu(:,6)
     &               - gijdu(:,4)*gijdu(:,3) )
           giju(:,5) = detgijI * (gijdu(:,4)*gijdu(:,6)
     &               - gijdu(:,5)*gijdu(:,2) )
           giju(:,6) = detgijI * (gijdu(:,4)*gijdu(:,5)
     &               - gijdu(:,1)*gijdu(:,6) )
c
        endif                   ! end of iDC.ne.0
c
c.... calculate (tau Lym) --> (rLymi)
c
      if(ires.ne.1 ) then
         rLymi(:,1) = tau(:,1) * rLymi(:,1) 
         rLymi(:,2) = tau(:,2) * rLymi(:,2)
         rLymi(:,3) = tau(:,2) * rLymi(:,3)
         rLymi(:,4) = tau(:,2) * rLymi(:,4)
         rLymi(:,5) = tau(:,3) * rLymi(:,5)
      endif
c
c  INCORRECT NOW    !      flops = flops + 255*npro
c
c
c.... return
c
        return
        end
c
c


      subroutine e3tau_nd  (rho,    cp,	    rmu,   T,
     &     u1,              u2,             u3,
     &     a1,              a2,             a3,
     &     con,             dxidx,          rLyi,  
     &     rLymi,           Tau,            rk,
     &     giju,            rTLS,           raLS,
     &     cv,              compK,          pres,
     &     A0inv,           dVdY)
c
c----------------------------------------------------------------------
c
c This routine computes the diagonal Tau for least-squares operator.  
c We receive the regular residual L operator and a
c modified residual L operator, calculate tau, and return values for
c tau and tau times these operators to maintain the format of past
c ENSA codes
c
c input:
c  rho    (npro)           : density
c  T      (npro)           : temperature
c  cp     (npro)           : specific heat at constant pressure
c  u1     (npro)           : x1-velocity component
c  u2     (npro)           : x2-velocity component
c  u3     (npro)           : x3-velocity component
c  dxidx  (npro,nsd,nsd)   : inverse of deformation gradient
c  rLyi   (npro,nflow)      : least-squares residual vector
c  rLymi   (npro,nflow)     : modified least-squares residual vector
c
c output:
c  rLyi   (npro,nflow)      : least-squares residual vector
c  rLymi   (npro,nflow)     : modified least-squares residual vector
c  Mtau    (npro,5,5)       : Matrix of Stabilized Parameters
c
c
c Zdenek Johan, Summer 1990.  (Modified from e2tau.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      include "common.h"
c
      dimension rho(npro),                 con(npro), 
     &            cp(npro),                  u1(npro),
     &            u2(npro),                  u3(npro),
     &            dxidx(npro,nsd,nsd),       rk(npro),
     &            rLyi(npro,nflow),
     &            rLymi(npro,nflow),         dVdY(npro,15), 
     &            rTLS(npro),                raLS(npro),
     &            rLyitemp(npro,nflow),      detgijI(npro)
c
      dimension   rmu(npro),	 cv(npro),
     &		  gijd(npro,6),  
     &		  fact(npro),	 giju(npro,6),
     &            A0inv(npro,15),gijdu(npro,6), compK(npro,10),
     &            a1(npro),    a2(npro),      a3(npro),
     &            T(npro),      pres(npro),     alphap(npro),
     &            betaT(npro),  gamb(npro),     c(npro),
     &            u(npro),      eb1(npro),      Q(npro,5,5),
     &            rlam(npro,5), eigmax(npro,5),   Pe(npro),
     &            Ak(npro),    B(npro),    D(npro),   E(npro),
     &            STau(npro,15),  Tau(npro,nflow,nflow)


c... build some necessary quantities at integration point:

c      alfap = one / T
c      betaT = one / pres     
      gamb  = gamma1
      c  = sqrt( (gamma * Rgas) * T )
      u = sqrt(u1**2+u2**2+u3**2)
      eb1 = cp*T - pt5*(u1**2+u2**2+u3**2) 

c... get eigenvectors for diagonalizing Tau_adv (e.v) before final rotations
c... Note: the ordering of eigenvectors corresponds to the following
c... ordering of the eigenvalues in the 1-st Euler Jacobian rotated to
c... the streamline coordinates

c     |u+c      | 
c     |   u     | 
c     |    u    |
c     |     u   |  ! for origins of this see Beam, Warming, Hyett
c     |      u-c|  ! Mathematics of Computation vol. 29 No. 132 p. 1037

c... different ordering assumed in Shakib/Johan entropy vars. code


      
      call e3eig1 (rho,             T,               cp,
     &               gamb,            c,               u1,
     &               u2,              u3,              a1,
     &               a2,              a3,              eb1,
     &               dxidx,           u,               Q)

c... Perform a Jacobi rotation on the Lambda matrix as well as adjust 
c... the eigenvectors. Tau still remains in entropy variables.



      call e3eig2 (u,               c,               a1,
     &             a2,              a3,              rlam,
     &             Q,               eigmax)

c
c.... invert the eigenvalues
c


      where (rlam .gt. ((epsM**2) * eigmax))
         rlam = one / sqrt(rlam)
      elsewhere
         rlam = zero
      endwhere

c
c.... flop count
c
   !      flops = flops + 65*npro

c.... reduce the diffusion contribution
c
 
        if (Navier .eq. 1) then
c
c.... calculate sigma
c
           
           do i = 1, nflow       ! diff. corr for every mode of Tau

              c(:) = pt33 * (
     &    Q(:,2,i) * ( compK(:, 1) * Q(:,2,i) + compK(:, 2) * Q(:,3,i)
     &               + compK(:, 4) * Q(:,4,i) + compK(:, 7) * Q(:,5,i) )
     &  + Q(:,3,i) * ( compK(:, 2) * Q(:,2,i) + compK(:, 3) * Q(:,3,i)
     &               + compK(:, 5) * Q(:,4,i) + compK(:, 8) * Q(:,5,i) )
     &  + Q(:,4,i) * ( compK(:, 4) * Q(:,2,i) + compK(:, 5) * Q(:,3,i)
     &               + compK(:, 6) * Q(:,4,i) + compK(:, 9) * Q(:,5,i) )
     &  + Q(:,5,i) * ( compK(:, 7) * Q(:,2,i) + compK(:, 8) * Q(:,3,i)
     &               + compK(:, 9) * Q(:,4,i) + compK(:,10) * Q(:,5,i) )
     &                  )

c... get Peclet Number
              

              Pe(:) = one / (c(:)*rlam(:,i)) 


c
c...  branch out into different polynomial orders
c


              if (ipord == 1) then ! linears - l = l^a*(coth(Pe)-1/Pe)
                                   ! approx. l = l^a*min(1/3*Pe,1)
                 rlam(:,i) = rlam(:,i)*min(pt33*Pe(:),one)
                 
              endif
              
              if (ipord == 2) then ! quads - Codina, CMAME' 92
                                ! approx. l = l^a*min(1/6*Pe,1/2)
                 rlam(:,i) = rlam(:,i)*min(pt33*pt5*Pe(:),pt5)
                 
              endif
              
              if (ipord == 3) then ! cubics - Recent Work
                                ! l = l^a*1/Pe*b
                                ! b is a real root of the
                                ! following polynomial:
             !  b^3+(-Pe*coth(Pe)b^2+(6/15*Pe^2-Pe*coth(Pe)+1)b+
             !  +(-1/15*Pe^3*coth(Pe)+6/15*Pe^2-Pe*coth(Pe)+1) = 0
             !  proceed setting up newton iteration w/ initial
             !  guess coming from the quad estimate.
                 
                 Ak(:)=1.0      ! get parameters for iteration

                 where(Pe.lt.5.0) ! approx. to hyp. cothangent
                    alphap = exp(Pe)
                    betaT = exp(-Pe)
                    c = (alphap+betaT)/(alphap-betaT)
                 elsewhere
                    c = one
                 endwhere

                 B= -Pe*c + Ak
                 D= 0.4 * (Pe**2) + B
                 E=-0.066666667 * (Pe**3) * c +D
                 
                                ! initial guess, use betaT
                 betaT(:) = Pe(:)*min(pt33*pt5*Pe(:),pt5)
                 
                                ! iteration - 3 iterations sufficient
                 do j = 1, 3

                    betaT=betaT-(Ak*betaT**3+B*betaT**2+D*betaT+E)/(3*Ak
     &                   *betaT**2+2*B*betaT+D)
                 enddo
                 
                 rlam(:,i) = rlam(:,i)*(1/Pe(:))*betaT(:)
                 
              endif             ! done w/ different polynomial orders
              
           enddo                ! done with loop over dof's
           
        endif                   ! done with diffusive correction
        

c
c.... form Tau (symmetric - entropy variables - then convert)
c
        STau(:, 1) = rlam(:,1) * Q(:,1,1) * Q(:,1,1) + 
     &                rlam(:,2) * Q(:,1,2) * Q(:,1,2) +
     &                rlam(:,3) * Q(:,1,3) * Q(:,1,3) +
     &                rlam(:,4) * Q(:,1,4) * Q(:,1,4) +
     &                rlam(:,5) * Q(:,1,5) * Q(:,1,5)
c
        STau(:, 2) = rlam(:,1) * Q(:,1,1) * Q(:,2,1) + 
     &                rlam(:,2) * Q(:,1,2) * Q(:,2,2) +
     &                rlam(:,3) * Q(:,1,3) * Q(:,2,3) +
     &                rlam(:,4) * Q(:,1,4) * Q(:,2,4) +
     &                rlam(:,5) * Q(:,1,5) * Q(:,2,5)
c
        STau(:, 3) = rlam(:,1) * Q(:,2,1) * Q(:,2,1) + 
     &                rlam(:,2) * Q(:,2,2) * Q(:,2,2) +
     &                rlam(:,3) * Q(:,2,3) * Q(:,2,3) +
     &                rlam(:,4) * Q(:,2,4) * Q(:,2,4) +
     &                rlam(:,5) * Q(:,2,5) * Q(:,2,5)
c
        STau(:, 4) = rlam(:,1) * Q(:,1,1) * Q(:,3,1) + 
     &                rlam(:,2) * Q(:,1,2) * Q(:,3,2) +
     &                rlam(:,3) * Q(:,1,3) * Q(:,3,3) +
     &                rlam(:,4) * Q(:,1,4) * Q(:,3,4) +
     &                rlam(:,5) * Q(:,1,5) * Q(:,3,5)
c
        STau(:, 5) = rlam(:,1) * Q(:,2,1) * Q(:,3,1) + 
     &                rlam(:,2) * Q(:,2,2) * Q(:,3,2) +
     &                rlam(:,3) * Q(:,2,3) * Q(:,3,3) +
     &                rlam(:,4) * Q(:,2,4) * Q(:,3,4) +
     &                rlam(:,5) * Q(:,2,5) * Q(:,3,5)
c
        STau(:, 6) = rlam(:,1) * Q(:,3,1) * Q(:,3,1) + 
     &                rlam(:,2) * Q(:,3,2) * Q(:,3,2) +
     &                rlam(:,3) * Q(:,3,3) * Q(:,3,3) +
     &                rlam(:,4) * Q(:,3,4) * Q(:,3,4) +
     &                rlam(:,5) * Q(:,3,5) * Q(:,3,5)
c
        STau(:, 7) = rlam(:,1) * Q(:,1,1) * Q(:,4,1) + 
     &                rlam(:,2) * Q(:,1,2) * Q(:,4,2) +
     &                rlam(:,3) * Q(:,1,3) * Q(:,4,3) +
     &                rlam(:,4) * Q(:,1,4) * Q(:,4,4) +
     &                rlam(:,5) * Q(:,1,5) * Q(:,4,5)
c
        STau(:, 8) = rlam(:,1) * Q(:,2,1) * Q(:,4,1) + 
     &                rlam(:,2) * Q(:,2,2) * Q(:,4,2) +
     &                rlam(:,3) * Q(:,2,3) * Q(:,4,3) +
     &                rlam(:,4) * Q(:,2,4) * Q(:,4,4) +
     &                rlam(:,5) * Q(:,2,5) * Q(:,4,5)
c
        STau(:, 9) = rlam(:,1) * Q(:,3,1) * Q(:,4,1) + 
     &                rlam(:,2) * Q(:,3,2) * Q(:,4,2) +
     &                rlam(:,3) * Q(:,3,3) * Q(:,4,3) +
     &                rlam(:,4) * Q(:,3,4) * Q(:,4,4) +
     &                rlam(:,5) * Q(:,3,5) * Q(:,4,5)
c
        STau(:,10) = rlam(:,1) * Q(:,4,1) * Q(:,4,1) + 
     &                rlam(:,2) * Q(:,4,2) * Q(:,4,2) +
     &                rlam(:,3) * Q(:,4,3) * Q(:,4,3) +
     &                rlam(:,4) * Q(:,4,4) * Q(:,4,4) +
     &                rlam(:,5) * Q(:,4,5) * Q(:,4,5)
c
        STau(:,11) = rlam(:,1) * Q(:,1,1) * Q(:,5,1) + 
     &                rlam(:,2) * Q(:,1,2) * Q(:,5,2) +
     &                rlam(:,3) * Q(:,1,3) * Q(:,5,3) +
     &                rlam(:,4) * Q(:,1,4) * Q(:,5,4) +
     &                rlam(:,5) * Q(:,1,5) * Q(:,5,5)
c
        STau(:,12) = rlam(:,1) * Q(:,2,1) * Q(:,5,1) + 
     &                rlam(:,2) * Q(:,2,2) * Q(:,5,2) +
     &                rlam(:,3) * Q(:,2,3) * Q(:,5,3) +
     &                rlam(:,4) * Q(:,2,4) * Q(:,5,4) +
     &                rlam(:,5) * Q(:,2,5) * Q(:,5,5)
c
        STau(:,13) = rlam(:,1) * Q(:,3,1) * Q(:,5,1) + 
     &                rlam(:,2) * Q(:,3,2) * Q(:,5,2) +
     &                rlam(:,3) * Q(:,3,3) * Q(:,5,3) +
     &                rlam(:,4) * Q(:,3,4) * Q(:,5,4) +
     &                rlam(:,5) * Q(:,3,5) * Q(:,5,5)
c
        STau(:,14) = rlam(:,1) * Q(:,4,1) * Q(:,5,1) + 
     &                rlam(:,2) * Q(:,4,2) * Q(:,5,2) +
     &                rlam(:,3) * Q(:,4,3) * Q(:,5,3) +
     &                rlam(:,4) * Q(:,4,4) * Q(:,5,4) +
     &                rlam(:,5) * Q(:,4,5) * Q(:,5,5)
c
        STau(:,15) = rlam(:,1) * Q(:,5,1) * Q(:,5,1) + 
     &                rlam(:,2) * Q(:,5,2) * Q(:,5,2) +
     &                rlam(:,3) * Q(:,5,3) * Q(:,5,3) +
     &                rlam(:,4) * Q(:,5,4) * Q(:,5,4) +
     &                rlam(:,5) * Q(:,5,5) * Q(:,5,5)


c
c... Form Primitive Variable Tau as [dY/dV]*tilde{Tau},
c... see G.Hauke's thesis appendix for [dY/dV] matrix
c
        betaT = cp*T + pt5*(u1**2+u2**2+u3**2) !reuse betaT
          
        Tau(:,1,1) = rho*T*(STau(:,1)+u1*STau(:,2)+
     &         u2*STau(:,4)+u3*STau(:,7)+betaT*STau(:,11))
        Tau(:,1,2) = rho*T*(STau(:,2)+u1*STau(:,3)+
     &         u2*STau(:,5)+u3*STau(:,8)+betaT*STau(:,12))
        Tau(:,1,3) = rho*T*(STau(:,4)+u1*STau(:,5)+
     &         u2*STau(:,6)+u3*STau(:,9)+betaT*STau(:,13))
        Tau(:,1,4) = rho*T*(STau(:,7)+u1*STau(:,8)+
     &         u2*STau(:,9)+u3*STau(:,10)+betaT*STau(:,14))
        Tau(:,1,5) = rho*T*(STau(:,11)+u1*STau(:,12)+
     &         u2*STau(:,13)+u3*STau(:,14)+betaT*STau(:,15))

          
        Tau(:,2,1) = T(:)*(STau(:,2) + u1(:)*STau(:,11))
        Tau(:,2,2) = T(:)*(STau(:,3) + u1(:)*STau(:,12))
        Tau(:,2,3) = T(:)*(STau(:,5) + u1(:)*STau(:,13))
        Tau(:,2,4) = T(:)*(STau(:,8) + u1(:)*STau(:,14))
        Tau(:,2,5) = T(:)*(STau(:,12) + u1(:)*STau(:,15))
        
        Tau(:,3,1) = T(:)*(STau(:,4) + u2(:)*STau(:,11))
        Tau(:,3,2) = T(:)*(STau(:,5) + u2(:)*STau(:,12))
        Tau(:,3,3) = T(:)*(STau(:,6) + u2(:)*STau(:,13))
        Tau(:,3,4) = T(:)*(STau(:,9) + u2(:)*STau(:,14))
        Tau(:,3,5) = T(:)*(STau(:,13) + u2(:)*STau(:,15))
        
        Tau(:,4,1) = T(:)*(STau(:,7) + u3(:)*STau(:,11))
        Tau(:,4,2) = T(:)*(STau(:,8) + u3(:)*STau(:,12))
        Tau(:,4,3) = T(:)*(STau(:,9) + u3(:)*STau(:,13))
        Tau(:,4,4) = T(:)*(STau(:,10) + u3(:)*STau(:,14))
        Tau(:,4,5) = T(:)*(STau(:,14) + u3(:)*STau(:,15))
        
        betaT = T**2
        
        Tau(:,5,1) = betaT(:)*STau(:,11)
        Tau(:,5,2) = betaT(:)*STau(:,12)
        Tau(:,5,3) = betaT(:)*STau(:,13)
        Tau(:,5,4) = betaT(:)*STau(:,14)
        Tau(:,5,5) = betaT(:)*STau(:,15)


c     
c...  done with conversion to pressure primitive variables
c...  now need to interface with the rest of the computations
c     
        
c...  finally multiply this tau matrix times the
c     two residual vectors
c
c ... calculate (tau Ly) --> (rLyi)
c ... storing rLyi for the DC term
          

        if(iDC .ne. 0) rLyitemp=rLyi

        if(ires.eq.3 .or. ires .eq. 1) then
           eigmax = rLyi  !reuse
           rLyi = zero
           do i = 1, nflow
              do  j = 1, nflow
                 rLyi(:,i) = rLyi(:,i) + Tau(:,i,j) * eigmax(:,j) 
              enddo
           enddo
        endif


        if(iDC .ne. 0) then
c.....calculation of rTLS & raLS for DC term
c
c.....calculation of (Ly - S).tau^tilda*(Ly - S) 
c 
           rTLS = rLYItemp(:,1)       * (rLyi(:,1)*dVdY(:,1)
     &        + dVdY(:,2)*rLyi(:,2) + dVdY(:,4)*rLyi(:,3)
     &        + rLyi(:,4)*dVdY(:,7) + dVdY(:,11)*rLyi(:,5))
     &        + rLyitemp(:,2)       * (rLyi(:,2)*dVdY(:,3)
     &        + rLyi(:,3)*dVdY(:,5) + dVdY(:,8)*rLyi(:,4)
     &        + rLyi(:,5)*dVdY(:,12))
     &        + rLyitemp(:,3)       * (rLyi(:,3)*dVdY(:,6)
     &        + dVdY(:,9)*rLyi(:,4) + dVdY(:,13)*rLyi(:,5))
     &        + rLyitemp(:,4)       * (rLyi(:,4)*dVdY(:,10)
     &        + dVdY(:,14)*rLyi(:,5))
     &        + rLyitemp(:,5)       * (dVdY(:,15)*rLyi(:,5))

c
c...... calculation of (Ly - S).A0inv*(Ly - S)
c
           raLS = two*rLyitemp(:,4)*rLyitemp(:,5)*A0inv(:,15)
     &          + two*rLyitemp(:,3)*rLyitemp(:,5)*A0inv(:,14)
     &          + two*rLyitemp(:,1)*rLyitemp(:,2)*A0inv( :,6)
     &          + two*rLyitemp(:,2)*rLyitemp(:,3)*A0inv(:,10)
     &          + two*rLyitemp(:,2)*rLyitemp(:,4)*A0inv(:,11)
     &          + two*rLyitemp(:,1)*rLyitemp(:,3)*A0inv( :,7)
     &          + two*rLyitemp(:,3)*rLyitemp(:,4)*A0inv(:,13)
     &          + two*rLyitemp(:,2)*rLyitemp(:,5)*A0inv(:,12)
     &          + two*rLyitemp(:,1)*rLyitemp(:,4)*A0inv( :,8)
     &          + two*rLyitemp(:,1)*rLyitemp(:,5)*A0inv( :,9)
     &          + rLyitemp(:,1)**2*A0inv(:,1)
     &          + rLyitemp(:,2)**2*A0inv(:,2)
     &          + rLyitemp(:,3)**2*A0inv(:,3)
     &          + rLyitemp(:,4)**2*A0inv(:,4)
     &          + rLyitemp(:,5)**2*A0inv(:,5)
c
c.....****************calculation of giju for DC term***************
c     
c.... for the notation of different numbering
c     
           call e3gijd( dxidx, gijd )

           gijdu(:,1)=gijd(:,1)
           gijdu(:,2)=gijd(:,3)
           gijdu(:,3)=gijd(:,6)
           gijdu(:,4)=gijd(:,2)
           gijdu(:,5)=gijd(:,4)
           gijdu(:,6)=gijd(:,5)
c     
c     
           detgijI = one/(
     &          gijdu(:,1) * gijdu(:,2) * gijdu(:,3)
     &          - gijdu(:,1) * gijdu(:,6) * gijdu(:,6)
     &          - gijdu(:,4) * gijdu(:,4) * gijdu(:,3)
     &          + gijdu(:,4) * gijdu(:,5) * gijdu(:,6) * two
     &          - gijdu(:,5) * gijdu(:,5) * gijdu(:,2) 
     &          )
           giju(:,1) = detgijI * (gijdu(:,2)*gijdu(:,3) 
     &               - gijdu(:,6)**2)
           giju(:,2) = detgijI * (gijdu(:,1)*gijdu(:,3) 
     &               - gijdu(:,5)**2)
           giju(:,3) = detgijI * (gijdu(:,1)*gijdu(:,2)
     &               - gijdu(:,4)**2)
           giju(:,4) = detgijI * (gijdu(:,5)*gijdu(:,6)
     &               - gijdu(:,4)*gijdu(:,3) )
           giju(:,5) = detgijI * (gijdu(:,4)*gijdu(:,6)
     &               - gijdu(:,5)*gijdu(:,2) )
           giju(:,6) = detgijI * (gijdu(:,4)*gijdu(:,5)
     &               - gijdu(:,1)*gijdu(:,6) )
c
        endif                   ! end of iDC.ne.0
c
c.... calculate (tau Lym) --> (rLymi)
c
        if(ires.ne.1 ) then
           eigmax = rLymi
           rLymi = zero
           do i = 1, nflow
              do j = 1, nflow
                 rLymi(:,i) = rLymi(:,i) + Tau(:,i,j) * eigmax(:,j) 
              enddo
           enddo
        endif
c
c  INCORRECT NOW    !      flops = flops + 255*npro
c
c
c.... return
c
      return
      end
c


      subroutine e3tau_nd_modal  (rho,    cp,	    rmu,   T,
     &     u1,              u2,             u3,
     &     a1,              a2,             a3,
     &     con,             dxidx,          rLyi,  
     &     rLymi,           Tau,            rk,
     &     giju,            rTLS,           raLS,
     &     cv,              compK,          pres,
     &     A0inv,           dVdY)
c     
c----------------------------------------------------------------------
c     
c     This routine computes the diagonal Tau for least-squares operator.
c     
c We receive the regular residual L operator and a
c modified residual L operator, calculate tau, and return values for
c tau and tau times these operators to maintain the format of past
c ENSA codes
c
c input:
c  rho    (npro)           : density
c  T      (npro)           : temperature
c  cp     (npro)           : specific heat at constant pressure
c  u1     (npro)           : x1-velocity component
c  u2     (npro)           : x2-velocity component
c  u3     (npro)           : x3-velocity component
c  dxidx  (npro,nsd,nsd)   : inverse of deformation gradient
c  rLyi   (npro,nflow)      : least-squares residual vector
c  rLymi   (npro,nflow)     : modified least-squares residual vector
c
c output:
c  rLyi   (npro,nflow)      : least-squares residual vector
c  rLymi   (npro,nflow)     : modified least-squares residual vector
c  Mtau    (npro,5,5)       : Matrix of Stabilized Parameters
c
c
c Zdenek Johan, Summer 1990.  (Modified from e2tau.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      include "common.h"
c
      dimension rho(npro),                 con(npro), 
     &            cp(npro),                  u1(npro),
     &            u2(npro),                  u3(npro),
     &            dxidx(npro,nsd,nsd),       rk(npro),
     &            rLyi(npro,nflow,ipord),
     &            rLymi(npro,nflow,ipord),   dVdY(npro,15), 
     &            rTLS(npro),                raLS(npro),
     &            rLyitemp(npro,nflow),      detgijI(npro)
c
      dimension   rmu(npro),	 cv(npro),
     &		  gijd(npro,6),  
     &		  fact(npro),	 giju(npro,6),
     &            A0inv(npro,15),gijdu(npro,6), compK(npro,10),
     &            a1(npro),    a2(npro),      a3(npro),
     &            T(npro),      pres(npro),     alphap(npro),
     &            betaT(npro),  gamb(npro),     c(npro),
     &            u(npro),      eb1(npro),      Q(npro,5,5),
     &            rlam(npro,5), eigmax(npro,5),   Pe(npro),
     &            STau(npro,15, ipord),  Tau(npro,nflow,nflow,ipord),
     &            rlamtmp(npro,5,ipord)


c... build some necessary quantities at integration point:

c      alfap = one / T
c      betaT = one / pres     
      gamb  = gamma1
      c  = sqrt( (gamma * Rgas) * T )
      u = sqrt(u1**2+u2**2+u3**2)
      eb1 = cp*T - pt5*(u1**2+u2**2+u3**2) 

c... get eigenvectors for diagonalizing Tau_adv (e.v) before final rotations
c... Note: the ordering of eigenvectors corresponds to the following
c... ordering of the eigenvalues in the 1-st Euler Jacobian rotated to
c... the streamline coordinates

c     |u+c      | 
c     |   u     | 
c     |    u    |
c     |     u   |  ! for origins of this see Beam, Warming, Hyett
c     |      u-c|  ! Mathematics of Computation vol. 29 No. 132 p. 1037

c... different ordering assumed in Shakib/Johan entropy vars. code


      
      call e3eig1 (rho,             T,               cp,
     &               gamb,            c,               u1,
     &               u2,              u3,              a1,
     &               a2,              a3,              eb1,
     &               dxidx,           u,               Q)

c... Perform a Jacobi rotation on the Lambda matrix as well as adjust 
c... the eigenvectors. Tau still remains in entropy variables.



      call e3eig2 (u,               c,               a1,
     &             a2,              a3,              rlam,
     &             Q,               eigmax)

c
c.... invert the eigenvalues
c


      where (rlam .gt. ((epsM**2) * eigmax))
         rlam = one / sqrt(rlam)
      elsewhere
         rlam = zero
      endwhere

      do i = 1, ipord
         rlamtmp(:,:,i) = rlam(:,:)
      enddo
      
c
c.... flop count
c
   !      flops = flops + 65*npro

c.... reduce the diffusion contribution
c
 
        if (Navier .eq. 1) then
c
c.... calculate sigma
c
           
           do i = 1, nflow       ! diff. corr for every mode of Tau

              c(:) = pt33 * (
     &    Q(:,2,i) * ( compK(:, 1) * Q(:,2,i) + compK(:, 2) * Q(:,3,i)
     &               + compK(:, 4) * Q(:,4,i) + compK(:, 7) * Q(:,5,i) )
     &  + Q(:,3,i) * ( compK(:, 2) * Q(:,2,i) + compK(:, 3) * Q(:,3,i)
     &               + compK(:, 5) * Q(:,4,i) + compK(:, 8) * Q(:,5,i) )
     &  + Q(:,4,i) * ( compK(:, 4) * Q(:,2,i) + compK(:, 5) * Q(:,3,i)
     &               + compK(:, 6) * Q(:,4,i) + compK(:, 9) * Q(:,5,i) )
     &  + Q(:,5,i) * ( compK(:, 7) * Q(:,2,i) + compK(:, 8) * Q(:,3,i)
     &               + compK(:, 9) * Q(:,4,i) + compK(:,10) * Q(:,5,i) )
     &                  )

c... get Peclet Number
              

              Pe(:) = one / (c(:)*rlam(:,i)) 


c
c...  branch out into different polynomial orders
c


              if (ipord == 1) then ! linears - l = l^a*(coth(Pe)-1/Pe)
                                   ! approx. l = l^a*min(1/3*Pe,1)
              rlamtmp(:,i,1) = rlamtmp(:,i,1)*min(pt33*Pe(:),one)
                 
              endif
              
              if (ipord == 2) then 

              rlamtmp(:,i,1) = rlamtmp(:,i,1)*min((1.0/15.0)*Pe(:),pt33)
              rlamtmp(:,i,2) = rlamtmp(:,i,2)*min((1.0/12.0)*Pe(:),pt5)
                 
              endif
              
              if (ipord == 3) then ! cubics - Recent Work
                 
                 do ii = 1, npro
                    if (Pe(ii).lt.3.0) then
                       rlamtmp(ii,i,1) = rlamtmp(ii,i,1)*
     &                      0.00054*Pe(ii)**2
                    endif
                    
                    if ((Pe(ii).ge.3).and.(Pe(ii).lt.17.20)) then
                       rlamtmp(ii,i,1) = rlamtmp(ii,i,1)*(0.0114*Pe(ii)
     &                      -0.0294)
                    endif
                    
                    if (Pe(ii).ge.17.20) then
                       rlamtmp(ii,i,1) = rlamtmp(ii,i,1)*(1.0/6.0)
                    endif
                    
                 enddo
                 
                 rlamtmp(:,i,2) = rlamtmp(:,i,2)*min((1.0/45.0)*Pe(:)
     &                ,0.2d0)
                 rlamtmp(:,i,3) = rlamtmp(:,i,3)*min((1.0/25.0)*Pe(:)
     &                ,pt33)   
                
                 
              endif             ! done w/ different polynomial orders
              
           enddo                ! done with loop over dof's
           
        endif                   ! done with diffusive correction
        

c
c.... form Tau (symmetric - entropy variables - then convert)
c
        do i = 1, ipord

        STau(:, 1, i) = rlamtmp(:,1,i) * Q(:,1,1) * Q(:,1,1) + 
     &                rlamtmp(:,2,i) * Q(:,1,2) * Q(:,1,2) +
     &                rlamtmp(:,3,i) * Q(:,1,3) * Q(:,1,3) +
     &                rlamtmp(:,4,i) * Q(:,1,4) * Q(:,1,4) +
     &                rlamtmp(:,5,i) * Q(:,1,5) * Q(:,1,5)
c
        STau(:, 2, i) = rlamtmp(:,1,i) * Q(:,1,1) * Q(:,2,1) + 
     &                rlamtmp(:,2,i) * Q(:,1,2) * Q(:,2,2) +
     &                rlamtmp(:,3,i) * Q(:,1,3) * Q(:,2,3) +
     &                rlamtmp(:,4,i) * Q(:,1,4) * Q(:,2,4) +
     &                rlamtmp(:,5,i) * Q(:,1,5) * Q(:,2,5)
c
        STau(:, 3, i) = rlamtmp(:,1,i) * Q(:,2,1) * Q(:,2,1) + 
     &                rlamtmp(:,2,i) * Q(:,2,2) * Q(:,2,2) +
     &                rlamtmp(:,3,i) * Q(:,2,3) * Q(:,2,3) +
     &                rlamtmp(:,4,i) * Q(:,2,4) * Q(:,2,4) +
     &                rlamtmp(:,5,i) * Q(:,2,5) * Q(:,2,5)
c
        STau(:, 4, i) = rlamtmp(:,1,i) * Q(:,1,1) * Q(:,3,1) + 
     &                rlamtmp(:,2,i) * Q(:,1,2) * Q(:,3,2) +
     &                rlamtmp(:,3,i) * Q(:,1,3) * Q(:,3,3) +
     &                rlamtmp(:,4,i) * Q(:,1,4) * Q(:,3,4) +
     &                rlamtmp(:,5,i) * Q(:,1,5) * Q(:,3,5)
c
        STau(:, 5, i) = rlamtmp(:,1,i) * Q(:,2,1) * Q(:,3,1) + 
     &                rlamtmp(:,2,i) * Q(:,2,2) * Q(:,3,2) +
     &                rlamtmp(:,3,i) * Q(:,2,3) * Q(:,3,3) +
     &                rlamtmp(:,4,i) * Q(:,2,4) * Q(:,3,4) +
     &                rlamtmp(:,5,i) * Q(:,2,5) * Q(:,3,5)
c
        STau(:, 6, i) = rlamtmp(:,1,i) * Q(:,3,1) * Q(:,3,1) + 
     &                rlamtmp(:,2,i) * Q(:,3,2) * Q(:,3,2) +
     &                rlamtmp(:,3,i) * Q(:,3,3) * Q(:,3,3) +
     &                rlamtmp(:,4,i) * Q(:,3,4) * Q(:,3,4) +
     &                rlamtmp(:,5,i) * Q(:,3,5) * Q(:,3,5)
c
        STau(:, 7, i) = rlamtmp(:,1,i) * Q(:,1,1) * Q(:,4,1) + 
     &                rlamtmp(:,2,i) * Q(:,1,2) * Q(:,4,2) +
     &                rlamtmp(:,3,i) * Q(:,1,3) * Q(:,4,3) +
     &                rlamtmp(:,4,i) * Q(:,1,4) * Q(:,4,4) +
     &                rlamtmp(:,5,i) * Q(:,1,5) * Q(:,4,5)
c
        STau(:, 8, i) = rlamtmp(:,1,i) * Q(:,2,1) * Q(:,4,1) + 
     &                rlamtmp(:,2,i) * Q(:,2,2) * Q(:,4,2) +
     &                rlamtmp(:,3,i) * Q(:,2,3) * Q(:,4,3) +
     &                rlamtmp(:,4,i) * Q(:,2,4) * Q(:,4,4) +
     &                rlamtmp(:,5,i) * Q(:,2,5) * Q(:,4,5)
c
        STau(:, 9, i) = rlamtmp(:,1,i) * Q(:,3,1) * Q(:,4,1) + 
     &                rlamtmp(:,2,i) * Q(:,3,2) * Q(:,4,2) +
     &                rlamtmp(:,3,i) * Q(:,3,3) * Q(:,4,3) +
     &                rlamtmp(:,4,i) * Q(:,3,4) * Q(:,4,4) +
     &                rlamtmp(:,5,i) * Q(:,3,5) * Q(:,4,5)
c
        STau(:,10, i) = rlamtmp(:,1,i) * Q(:,4,1) * Q(:,4,1) + 
     &                rlamtmp(:,2,i) * Q(:,4,2) * Q(:,4,2) +
     &                rlamtmp(:,3,i) * Q(:,4,3) * Q(:,4,3) +
     &                rlamtmp(:,4,i) * Q(:,4,4) * Q(:,4,4) +
     &                rlamtmp(:,5,i) * Q(:,4,5) * Q(:,4,5)
c
        STau(:,11, i) = rlamtmp(:,1,i) * Q(:,1,1) * Q(:,5,1) + 
     &                rlamtmp(:,2,i) * Q(:,1,2) * Q(:,5,2) +
     &                rlamtmp(:,3,i) * Q(:,1,3) * Q(:,5,3) +
     &                rlamtmp(:,4,i) * Q(:,1,4) * Q(:,5,4) +
     &                rlamtmp(:,5,i) * Q(:,1,5) * Q(:,5,5)
c
        STau(:,12, i) = rlamtmp(:,1,i) * Q(:,2,1) * Q(:,5,1) + 
     &                rlamtmp(:,2,i) * Q(:,2,2) * Q(:,5,2) +
     &                rlamtmp(:,3,i) * Q(:,2,3) * Q(:,5,3) +
     &                rlamtmp(:,4,i) * Q(:,2,4) * Q(:,5,4) +
     &                rlamtmp(:,5,i) * Q(:,2,5) * Q(:,5,5)
c
        STau(:,13, i) = rlamtmp(:,1,i) * Q(:,3,1) * Q(:,5,1) + 
     &                rlamtmp(:,2,i) * Q(:,3,2) * Q(:,5,2) +
     &                rlamtmp(:,3,i) * Q(:,3,3) * Q(:,5,3) +
     &                rlamtmp(:,4,i) * Q(:,3,4) * Q(:,5,4) +
     &                rlamtmp(:,5,i) * Q(:,3,5) * Q(:,5,5)
c
        STau(:,14, i) = rlamtmp(:,1,i) * Q(:,4,1) * Q(:,5,1) + 
     &                rlamtmp(:,2,i) * Q(:,4,2) * Q(:,5,2) +
     &                rlamtmp(:,3,i) * Q(:,4,3) * Q(:,5,3) +
     &                rlamtmp(:,4,i) * Q(:,4,4) * Q(:,5,4) +
     &                rlamtmp(:,5,i) * Q(:,4,5) * Q(:,5,5)
c
        STau(:,15, i) = rlamtmp(:,1,i) * Q(:,5,1) * Q(:,5,1) + 
     &                rlamtmp(:,2,i) * Q(:,5,2) * Q(:,5,2) +
     &                rlamtmp(:,3,i) * Q(:,5,3) * Q(:,5,3) +
     &                rlamtmp(:,4,i) * Q(:,5,4) * Q(:,5,4) +
     &                rlamtmp(:,5,i) * Q(:,5,5) * Q(:,5,5)

      enddo

c
c... Form Primitive Variable Tau as [dY/dV]*tilde{Tau},
c... see G.Hauke's thesis appendix for [dY/dV] matrix
c
      do k = 1, ipord

         betaT = cp*T + pt5*(u1**2+u2**2+u3**2) !reuse betaT
         
         Tau(:,1,1,k) = rho*T*(STau(:,1,k)+u1*STau(:,2,k)+
     &        u2*STau(:,4,k)+u3*STau(:,7,k)+betaT*STau(:,11,k))
         Tau(:,1,2,k) = rho*T*(STau(:,2,k)+u1*STau(:,3,k)+
     &        u2*STau(:,5,k)+u3*STau(:,8,k)+betaT*STau(:,12,k))
         Tau(:,1,3,k) = rho*T*(STau(:,4,k)+u1*STau(:,5,k)+
     &        u2*STau(:,6,k)+u3*STau(:,9,k)+betaT*STau(:,13,k))
         Tau(:,1,4,k) = rho*T*(STau(:,7,k)+u1*STau(:,8,k)+
     &        u2*STau(:,9,k)+u3*STau(:,10,k)+betaT*STau(:,14,k))
         Tau(:,1,5,k) = rho*T*(STau(:,11,k)+u1*STau(:,12,k)+
     &        u2*STau(:,13,k)+u3*STau(:,14,k)+betaT*STau(:,15,k))
         
         
         Tau(:,2,1,k) = T(:)*(STau(:,2,k) + u1(:)*STau(:,11,k))
         Tau(:,2,2,k) = T(:)*(STau(:,3,k) + u1(:)*STau(:,12,k))
         Tau(:,2,3,k) = T(:)*(STau(:,5,k) + u1(:)*STau(:,13,k))
         Tau(:,2,4,k) = T(:)*(STau(:,8,k) + u1(:)*STau(:,14,k))
         Tau(:,2,5,k) = T(:)*(STau(:,12,k) + u1(:)*STau(:,15,k))
         
         Tau(:,3,1,k) = T(:)*(STau(:,4,k) + u2(:)*STau(:,11,k))
         Tau(:,3,2,k) = T(:)*(STau(:,5,k) + u2(:)*STau(:,12,k))
         Tau(:,3,3,k) = T(:)*(STau(:,6,k) + u2(:)*STau(:,13,k))
         Tau(:,3,4,k) = T(:)*(STau(:,9,k) + u2(:)*STau(:,14,k))
         Tau(:,3,5,k) = T(:)*(STau(:,13,k) + u2(:)*STau(:,15,k))
         
         Tau(:,4,1,k) = T(:)*(STau(:,7,k) + u3(:)*STau(:,11,k))
         Tau(:,4,2,k) = T(:)*(STau(:,8,k) + u3(:)*STau(:,12,k))
         Tau(:,4,3,k) = T(:)*(STau(:,9,k) + u3(:)*STau(:,13,k))
         Tau(:,4,4,k) = T(:)*(STau(:,10,k) + u3(:)*STau(:,14,k))
         Tau(:,4,5,k) = T(:)*(STau(:,14,k) + u3(:)*STau(:,15,k))
         
         betaT = T**2
         
         Tau(:,5,1,k) = betaT(:)*STau(:,11,k)
         Tau(:,5,2,k) = betaT(:)*STau(:,12,k)
         Tau(:,5,3,k) = betaT(:)*STau(:,13,k)
         Tau(:,5,4,k) = betaT(:)*STau(:,14,k)
         Tau(:,5,5,k) = betaT(:)*STau(:,15,k)
         
      enddo
      
c     
c...  done with conversion to pressure primitive variables
c...  now need to interface with the rest of the computations
c     
        
c...  finally multiply this tau matrix times the
c     two residual vectors
c
c ... calculate (tau Ly) --> (rLyi)
c ... storing rLyi for the DC term
          

        if(iDC .ne. 0) rLyitemp(:,:)=rLyi(:,:,1)

        if(ires.eq.3 .or. ires .eq. 1) then
           eigmax(:,:) = rLyi(:,:,1) !reuse
           rLyi = zero
           do k = 1, ipord
              do i = 1, nflow
                 do  j = 1, nflow
                    rLyi(:,i,k) = rLyi(:,i,k)+Tau(:,i,j,k)*eigmax(:,j) 
                 enddo
              enddo
           enddo
        endif
        
        
        if(iDC .ne. 0) then
c.....calculation of rTLS & raLS for DC term
c
c.....calculation of (Ly - S).tau^tilda*(Ly - S) 
c 
           rTLS = rLYItemp(:,1)     * (rLyi(:,1,1)*dVdY(:,1)
     &        + dVdY(:,2)*rLyi(:,2,1) + dVdY(:,4)*rLyi(:,3,1)
     &        + rLyi(:,4,1)*dVdY(:,7) + dVdY(:,11)*rLyi(:,5,1))
     &        + rLyitemp(:,2)       * (rLyi(:,2,1)*dVdY(:,3)
     &        + rLyi(:,3,1)*dVdY(:,5) + dVdY(:,8)*rLyi(:,4,1)
     &        + rLyi(:,5,1)*dVdY(:,12))
     &        + rLyitemp(:,3)       * (rLyi(:,3,1)*dVdY(:,6)
     &        + dVdY(:,9)*rLyi(:,4,1) + dVdY(:,13)*rLyi(:,5,1))
     &        + rLyitemp(:,4)       * (rLyi(:,4,1)*dVdY(:,10)
     &        + dVdY(:,14)*rLyi(:,5,1))
     &        + rLyitemp(:,5)       * (dVdY(:,15)*rLyi(:,5,1))

c
c...... calculation of (Ly - S).A0inv*(Ly - S)
c
           raLS = two*rLyitemp(:,4)*rLyitemp(:,5)*A0inv(:,15)
     &          + two*rLyitemp(:,3)*rLyitemp(:,5)*A0inv(:,14)
     &          + two*rLyitemp(:,1)*rLyitemp(:,2)*A0inv( :,6)
     &          + two*rLyitemp(:,2)*rLyitemp(:,3)*A0inv(:,10)
     &          + two*rLyitemp(:,2)*rLyitemp(:,4)*A0inv(:,11)
     &          + two*rLyitemp(:,1)*rLyitemp(:,3)*A0inv( :,7)
     &          + two*rLyitemp(:,3)*rLyitemp(:,4)*A0inv(:,13)
     &          + two*rLyitemp(:,2)*rLyitemp(:,5)*A0inv(:,12)
     &          + two*rLyitemp(:,1)*rLyitemp(:,4)*A0inv( :,8)
     &          + two*rLyitemp(:,1)*rLyitemp(:,5)*A0inv( :,9)
     &          + rLyitemp(:,1)**2*A0inv(:,1)
     &          + rLyitemp(:,2)**2*A0inv(:,2)
     &          + rLyitemp(:,3)**2*A0inv(:,3)
     &          + rLyitemp(:,4)**2*A0inv(:,4)
     &          + rLyitemp(:,5)**2*A0inv(:,5)
c
c.....****************calculation of giju for DC term***************
c     
c.... for the notation of different numbering
c     
           gijdu(:,1)=gijd(:,1)
           gijdu(:,2)=gijd(:,3)
           gijdu(:,3)=gijd(:,6)
           gijdu(:,4)=gijd(:,2)
           gijdu(:,5)=gijd(:,4)
           gijdu(:,6)=gijd(:,5)
c     
c     
           detgijI = one/(
     &          gijdu(:,1) * gijdu(:,2) * gijdu(:,3)
     &          - gijdu(:,1) * gijdu(:,6) * gijdu(:,6)
     &          - gijdu(:,4) * gijdu(:,4) * gijdu(:,3)
     &          + gijdu(:,4) * gijdu(:,5) * gijdu(:,6) * two
     &          - gijdu(:,5) * gijdu(:,5) * gijdu(:,2) 
     &          )
           giju(:,1) = detgijI * (gijdu(:,2)*gijdu(:,3) 
     &               - gijdu(:,6)**2)
           giju(:,2) = detgijI * (gijdu(:,1)*gijdu(:,3) 
     &               - gijdu(:,5)**2)
           giju(:,3) = detgijI * (gijdu(:,1)*gijdu(:,2)
     &               - gijdu(:,4)**2)
           giju(:,4) = detgijI * (gijdu(:,5)*gijdu(:,6)
     &               - gijdu(:,4)*gijdu(:,3) )
           giju(:,5) = detgijI * (gijdu(:,4)*gijdu(:,6)
     &               - gijdu(:,5)*gijdu(:,2) )
           giju(:,6) = detgijI * (gijdu(:,4)*gijdu(:,5)
     &               - gijdu(:,1)*gijdu(:,6) )
c
        endif                   ! end of iDC.ne.0
c
c.... calculate (tau Lym) --> (rLymi)
c
        if(ires.ne.1 ) then
           eigmax(:,:) = rLymi(:,:,1)
           rLymi = zero
           do k = 1, ipord
              do i = 1, nflow
                 do j = 1, nflow
         rLymi(:,i,k) = rLymi(:,i,k)+Tau(:,i,j,k)*eigmax(:,j)                 
                 enddo
              enddo
           enddo
        endif
c
c  INCORRECT NOW    !      flops = flops + 255*npro
c
c
c.... return
c
      return
      end
c



        subroutine e3tauSclr(rho,    rmu,    A0t, 
     &                       u1,     u2,     u3,
     &                       dxidx,  rLyti,  rLymti,
     &                       taut,   rk,     raLSt,
     &                       rTLSt,  giju)
c
c----------------------------------------------------------------------
c
c This routine computes the diagonal Tau for least-squares operator.  
c This Tau is the one proposed for nearly incompressible flows by
c Franca and Frey.  We receive the regular residual L operator and a
c modified residual L operator, calculate tau, and return values for
c tau and tau times these operators to maintain the format of past
c ENSA codes
c
c input:
c  rho    (npro)           : density
c  T      (npro)           : temperature
c  cp     (npro)           : specific heat at constant pressure
c  u1     (npro)           : x1-velocity component
c  u2     (npro)           : x2-velocity component
c  u3     (npro)           : x3-velocity component
c  dxidx  (npro,nsd,nsd)   : inverse of deformation gradient
c  rLyti   (npro)          : least-squares residual vector
c  rLymti   (npro)         : modified least-squares residual vector
c
c output:
c  rLyti   (npro,nflow)     : least-squares residual vector
c  rLymti   (npro,nflow)    : modified least-squares residual vector
c  tau    (npro,3)         : 3 taus
c
c
c Zdenek Johan, Summer 1990.  (Modified from e2tau.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use turbSA
      include "common.h"
c
      dimension rho(npro),                 T(npro),
     &            cp(npro),                  u1(npro),
     &            u2(npro),                  u3(npro),
     &            dxidx(npro,nsd,nsd),       rk(npro),
     &            taut(npro),                rLyti(npro),
     &            rLymti(npro)
c
        dimension rmu(npro),                 A0t(npro),
     &		  gijd(npro,6),              uh1(npro),
     &		  fact(npro),	             h2o2u(npro),
     &            rlytitemp(npro),           raLSt(npro),
     &            rTLSt(npro),               gijdu(npro,6),
     &            giju(npro,6),              detgijI(npro)
c
c      
      call e3gijd( dxidx, gijd )

c
c  next form the diffusive length scale |u| h_1 = 2 ( ui (gijd)-1 uj)^{1/2}
c
c   dividing factor for the inverse of gijd
c
      fact = gijd(:,1) * gijd(:,3) * gijd(:,6)
     &     - gijd(:,1) * gijd(:,5) * gijd(:,5)
     &     - gijd(:,3) * gijd(:,4) * gijd(:,4)
     &     - gijd(:,6) * gijd(:,2) * gijd(:,2)
     &     + gijd(:,2) * gijd(:,4) * gijd(:,5) * two
c
      uh1=    u1*u1*(gijd(:,3)*gijd(:,6)-gijd(:,5)*gijd(:,5))
     &     + u2*u2*(gijd(:,1)*gijd(:,6)-gijd(:,4)*gijd(:,4))
     &     + u3*u3*(gijd(:,1)*gijd(:,3)-gijd(:,2)*gijd(:,2))
     &     + two *(u1*u2*(gijd(:,4)*gijd(:,5)-gijd(:,2)*gijd(:,6))
     &     + u1*u3*(gijd(:,2)*gijd(:,5)-gijd(:,4)*gijd(:,3))
     &     + u2*u3*(gijd(:,4)*gijd(:,2)-gijd(:,1)*gijd(:,5)))
c
c   at this point we have (u h1)^2 * inverse coefficient^2 / 4
c
c                                    ^ fact
c

      uh1= two * sqrt(uh1/fact)

c
c  next form the advective length scale |u|/h_2 = 2 ( ui (gijd) uj)^{1/2}
c
      h2o2u =   u1*u1*gijd(:,1)
     &     + u2*u2*gijd(:,3)
     &     + u3*u3*gijd(:,6)
     &     +(u1*u2*gijd(:,2)
     &     + u1*u3*gijd(:,4)
     &     + u2*u3*gijd(:,5))*two  + 1.0e-15 !FIX FOR INVALID MESHES
c
c  at this point we have (2 u / h_2)^2
c

c       call tnanqe(h2o2u,1,"riaconv ")

      h2o2u = one / sqrt(h2o2u) ! this flips it over leaves it h_2/(2u)
c
c...  momentum tau
c 
c
c... rmu will now hold the total (molecular plus eddy) viscosity...
      dts=one/(Dtgl*dtsfct)
      if(iremoveStabTimeTerm.gt.0) dts = dts*100000  ! remove time term from scalar
! Duct code had this       dts=1.0e16
      taut(:)=min(dts,h2o2u)
      taut(:)=taut(:)/rho
      taut(:)=min(taut(:),h2o2u*h2o2u*rk*pt66*saSigma/rmu)
c
c...  finally multiply this tau matrix times the
c     two residual vectors
c
c.... calculate (tau Lyt) --> (rLyti)
c
c.... storing rLyi for the DC term
          rLytitemp=rLyti

	if(ires.eq.3 .or. ires .eq. 1) then
          rLyti(:) = taut(:) * rLyti(:) 

        endif
        if(iDCSclr(1) .ne. 0) then
c..... calculation of rTLS & raLS for DC term
c..... calculation of (Ly - S).tau^tilda*(Ly - S) 
c
         rTLSt = rLYtItemp(:)*rLyti(:)
c
c...... calculation of (Ly - S).A0inv*(Ly - S)
c
         raLSt = rLYtItemp(:) * rLYtItemp(:)
c.....*****************calculation of giju for DC term******************
c
c.... for the notation of different numbering
c
           gijdu(:,1)=gijd(:,1)
           gijdu(:,2)=gijd(:,3)
           gijdu(:,3)=gijd(:,6)
           gijdu(:,4)=gijd(:,2)
           gijdu(:,5)=gijd(:,4)
           gijdu(:,6)=gijd(:,5)
c
c  we are going to need this in the dc factor later so we calculate it.
c
         detgijI = one/(
     &             gijdu(:,1) * gijdu(:,2) * gijdu(:,3)
     &           - gijdu(:,1) * gijdu(:,6) * gijdu(:,6)
     &           - gijdu(:,4) * gijdu(:,4) * gijdu(:,3)
     &           + gijdu(:,4) * gijdu(:,5) * gijdu(:,6) * two
     &           - gijdu(:,5) * gijdu(:,5) * gijdu(:,2) 
     &                 )
          giju(:,1) = detgijI * (gijdu(:,2)*gijdu(:,3) 
     &              - gijdu(:,6)**2)
          giju(:,2) = detgijI * (gijdu(:,1)*gijdu(:,3) 
     &              - gijdu(:,5)**2)
          giju(:,3) = detgijI * (gijdu(:,1)*gijdu(:,2)
     &              - gijdu(:,4)**2)
          giju(:,4) = detgijI * (gijdu(:,5)*gijdu(:,6)
     &              - gijdu(:,4)*gijdu(:,3) )
          giju(:,5) = detgijI * (gijdu(:,4)*gijdu(:,6)
     &              - gijdu(:,5)*gijdu(:,2) )
          giju(:,6) = detgijI * (gijdu(:,4)*gijdu(:,5)
     &              - gijdu(:,1)*gijdu(:,6) )
c
         endif    ! end of iDCSclr(1).ne.0
c
c.... calculate (tau Lym) --> (rLymi)
c
c        if(ires.ne.1 ) then
c          rLymi(:,1) = tau(:,1) * rLymi(:,1) 
c          rLymi(:,2) = tau(:,2) * rLymi(:,2)
c          rLymi(:,3) = tau(:,2) * rLymi(:,3)
c          rLymi(:,4) = tau(:,2) * rLymi(:,4)
c          rLymi(:,5) = tau(:,3) * rLymi(:,5)
c        endif
c
c  INCORRECT NOW    !      flops = flops + 255*npro
c
c
c.... return
c
      return
      end

c-----------------------------------------------------------------------
c get the metric tensor g_{ij}=xi_{k,i} xi_{k,j}.  
c-----------------------------------------------------------------------
      subroutine e3gijd( dxidx,  gijd )
      
      include "common.h"
      
      real*8  dxidx(npro,nsd,nsd),  gijd(npro,6),
     &        tmp1(npro),           tmp2(npro),
     &        tmp3(npro)
c  form metric tensor g_{ij}=xi_{k,i} xi_{k,j}.  It is a symmetric
c  tensor so we only form 6 components and use symmetric matrix numbering.
c  (d for down since giju=[gijd]^{-1})
c  (Note FARZIN and others use numbering of 1,2,3 being diagonal 456 off)
      if (lcsyst .ge. 2) then

         gijd(:,1) = dxidx(:,1,1) * dxidx(:,1,1)
     &            + dxidx(:,2,1) * dxidx(:,2,1)
     &            + dxidx(:,3,1) * dxidx(:,3,1)
c
         gijd(:,2) = dxidx(:,1,1) * dxidx(:,1,2)
     &            + dxidx(:,2,1) * dxidx(:,2,2)
     &            + dxidx(:,3,1) * dxidx(:,3,2)
c
         gijd(:,3) = dxidx(:,1,2) * dxidx(:,1,2)
     &            + dxidx(:,2,2) * dxidx(:,2,2)
     &            + dxidx(:,3,2) * dxidx(:,3,2)
c
         gijd(:,4) = dxidx(:,1,1) * dxidx(:,1,3)
     &            + dxidx(:,2,1) * dxidx(:,2,3)
     &            + dxidx(:,3,1) * dxidx(:,3,3)
c
         gijd(:,5) = dxidx(:,1,2) * dxidx(:,1,3)
     &            + dxidx(:,2,2) * dxidx(:,2,3)
     &            + dxidx(:,3,2) * dxidx(:,3,3)
c
         gijd(:,6) = dxidx(:,1,3) * dxidx(:,1,3)
     &            + dxidx(:,2,3) * dxidx(:,2,3)
     &        + dxidx(:,3,3) * dxidx(:,3,3)
c
      else   if (lcsyst .eq. 1) then   
c
c  There is an invariance problem with tets 
c  It is fixed by the following modifications to gijd 
c

         c1 = 1.259921049894873D+00
         c2 = 6.299605249474365D-01
c
         tmp1(:) = c1 * dxidx(:,1,1) + c2 *(dxidx(:,2,1)+dxidx(:,3,1))
         tmp2(:) = c1 * dxidx(:,2,1) + c2 *(dxidx(:,1,1)+dxidx(:,3,1))
         tmp3(:) = c1 * dxidx(:,3,1) + c2 *(dxidx(:,1,1)+dxidx(:,2,1))
         gijd(:,1) = dxidx(:,1,1) * tmp1
     1             + dxidx(:,2,1) * tmp2
     2             + dxidx(:,3,1) * tmp3
c
         tmp1(:) = c1 * dxidx(:,1,2) + c2 *(dxidx(:,2,2)+dxidx(:,3,2))
         tmp2(:) = c1 * dxidx(:,2,2) + c2 *(dxidx(:,1,2)+dxidx(:,3,2))
         tmp3(:) = c1 * dxidx(:,3,2) + c2 *(dxidx(:,1,2)+dxidx(:,2,2))
         gijd(:,2) = dxidx(:,1,1) * tmp1
     1             + dxidx(:,2,1) * tmp2
     2             + dxidx(:,3,1) * tmp3
c
         gijd(:,3) = dxidx(:,1,2) * tmp1
     1             + dxidx(:,2,2) * tmp2
     2             + dxidx(:,3,2) * tmp3
c
         tmp1(:) = c1 * dxidx(:,1,3) + c2 *(dxidx(:,2,3)+dxidx(:,3,3))
         tmp2(:) = c1 * dxidx(:,2,3) + c2 *(dxidx(:,1,3)+dxidx(:,3,3))
         tmp3(:) = c1 * dxidx(:,3,3) + c2 *(dxidx(:,1,3)+dxidx(:,2,3))
         gijd(:,4) = dxidx(:,1,1) * tmp1
     1             + dxidx(:,2,1) * tmp2
     2             + dxidx(:,3,1) * tmp3
c
         gijd(:,5) = dxidx(:,1,2) * tmp1
     1             + dxidx(:,2,2) * tmp2
     2             + dxidx(:,3,2) * tmp3
c
         gijd(:,6) = dxidx(:,1,3) * tmp1
     1             + dxidx(:,2,3) * tmp2
     2             + dxidx(:,3,3) * tmp3
c
      else  
c  This is just the hex copied from above.  I have
c  to find my notes on invariance factors for wedges
c         

         gijd(:,1) = dxidx(:,1,1) * dxidx(:,1,1)
     &            + dxidx(:,2,1) * dxidx(:,2,1)
     &            + dxidx(:,3,1) * dxidx(:,3,1)
c
         gijd(:,2) = dxidx(:,1,1) * dxidx(:,1,2)
     &            + dxidx(:,2,1) * dxidx(:,2,2)
     &            + dxidx(:,3,1) * dxidx(:,3,2)
c
         gijd(:,3) = dxidx(:,1,2) * dxidx(:,1,2)
     &            + dxidx(:,2,2) * dxidx(:,2,2)
     &            + dxidx(:,3,2) * dxidx(:,3,2)
c
         gijd(:,4) = dxidx(:,1,1) * dxidx(:,1,3)
     &            + dxidx(:,2,1) * dxidx(:,2,3)
     &            + dxidx(:,3,1) * dxidx(:,3,3)
c
         gijd(:,5) = dxidx(:,1,2) * dxidx(:,1,3)
     &            + dxidx(:,2,2) * dxidx(:,2,3)
     &            + dxidx(:,3,2) * dxidx(:,3,3)
c
         gijd(:,6) = dxidx(:,1,3) * dxidx(:,1,3)
     &            + dxidx(:,2,3) * dxidx(:,2,3)
     &            + dxidx(:,3,3) * dxidx(:,3,3)
      endif
c
      return
      end
        
