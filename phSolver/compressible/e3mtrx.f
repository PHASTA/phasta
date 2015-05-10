        subroutine e3mtrx (rho,    pres,   T,
     &                     ei,     h,	   alfap,
     &                     betaT,  cp,     rk,
     &                     u1,     u2,   u3,
     &                     A0,     A1,    
     &                     A2,     A3,
     &                     e2p,    e3p,    e4p,
     &                     drdp,   drdT,   A0DC, 
     &                     A0inv,  dVdY)
c
c----------------------------------------------------------------------
c
c This routine sets up the necessary matrices at the integration point.
c 
c input:
c  rho   (npro)         : density
c  pres  (npro)         : pressure
c  T     (npro)         : temperature
c  ei    (npro)         : internal energy
c  h     (npro)         : enthalpy
c  alfap (npro)         : expansivity
c  betaT (npro)         : isothermal compressibility
c  cp    (npro)         : specific heat at constant pressure
c  c     (npro)         : speed of sound
c  rk    (npro)         : kinetic energy
c  u1    (npro)         : x1-velocity component
c  u2    (npro)         : x2-velocity component
c  u3    (npro)         : x3-velocity component
c
c output:
c  A0    (npro,nflow,nflow)  : A0 matrix   
c  A1   (npro,nflow,nflow)  : A_1 matrix   
c  A2   (npro,nflow,nflow)  : A_2 matrix     
c  A3   (npro,nflow,nflow)  : A_3 matrix    
c
c Note: the definition of the matrices can be found in 
c       thesis by Hauke. 
c
c Zdenek Johan, Summer 1990.  (Modified from e2mtrx.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c Kenneth Jansen, Winter 1997 Primitive Variables
c----------------------------------------------------------------------
c
        include "common.h"
c
c
c  passed arrays
c
        dimension rho(npro),                 pres(npro),
     &            T(npro),                   ei(npro),
     &            h(npro),                   alfap(npro),
     &            betaT(npro),               
     &            cp(npro),                  
     &            rk(npro),
     &            u1(npro),                 u2(npro),
     &            u3(npro),                 fact1(npro),
     &            A0(npro,nflow,nflow),     dVdY(npro,15),   
     &            A1(npro,nflow,nflow),     A2(npro,nflow,nflow),
     &            A3(npro,nflow,nflow),     A0DC(npro,4),
     &            A0inv(npro,15),           d(npro),
     &            fact2(npro),               s1(npro),
     &            e1bar(npro),              e2bar(npro),
     &            e3bar(npro),              e4bar(npro),
     &            e5bar(npro),              c1bar(npro),
     &            c2bar(npro),              cv(npro),
     &            c3bar(npro),              u12(npro),
     &            u31(npro),                u23(npro)
c
c  local work arrays that are passed shared space
c
        dimension e2p(npro),                  
     &            e3p(npro),                 e4p(npro),
     &            drdp(npro),                drdT(npro)

	ttim(21) = ttim(21) - secs(0.0)
c
c.... initialize
c
        A0 = zero
        A1 = zero
        A2 = zero
        A3 = zero
c
c.... set up the constants
c
c
        drdp = rho * betaT
        drdT = -rho * alfap
        A0(:,5,1) = drdp * (h + rk)  - alfap * T    ! e1p
c        A0(:,5,1) = drdp * (ei + rk) + betaT * pres - alfap * T    ! e1p
          e2p  = A0(:,5,1) + one
          e3p  = rho * ( h + rk)
          e4p  = drdT * (h + rk) + rho * cp
c
c
c.... Calculate A0
c
        A0(:,1,1) = drdp 
c       A0(:,1,2) = zero 
c       A0(:,1,3) = zero 
c       A0(:,1,4) = zero 
        A0(:,1,5) = drdT 
c
        A0(:,2,1) = drdp * u1 
        A0(:,2,2) = rho 
c       A0(:,2,3) = zero 
c       A0(:,2,4) = zero 
        A0(:,2,5) = drdT * u1 
c
        A0(:,3,1) = drdp * u2 
c       A0(:,3,2) = zero 
        A0(:,3,3) = rho 
c       A0(:,3,4) = zero 
        A0(:,3,5) = drdT * u2 
c
        A0(:,4,1) = drdp * u3 
c       A0(:,4,2) = zero 
c       A0(:,4,3) = zero 
        A0(:,4,4) = rho 
        A0(:,4,5) = drdT * u3 
c
covered above       A0(:,5,1) = drdp * u1 
        A0(:,5,2) = rho * u1 
        A0(:,5,3) = rho * u2 
        A0(:,5,4) = rho * u3 
        A0(:,5,5) = e4p
c
        flops = flops + 67*npro
c
c.... Calculate A-tilde-1, A-tilde-2 and A-tilde-3
c
        A1(:,1,1) = drdp * u1
        A1(:,1,2) = rho
c       A1(:,1,3) = zero
c       A1(:,1,4) = zero
        A1(:,1,5) = drdT * u1
c
        A1(:,2,1) = drdp * u1 * u1 +1
        A1(:,2,2) = two *rho  * u1
c       A1(:,2,3) = zero
c       A1(:,2,4) = zero
        A1(:,2,5) = drdT * u1 * u1
c
        A1(:,3,1) = drdp * u1 * u2 
        A1(:,3,2) = rho  * u2
        A1(:,3,3) = rho  * u1
c       A1(:,3,4) = zero
        A1(:,3,5) = drdT * u1 * u2
c
        A1(:,4,1) = drdp * u1 * u3 
        A1(:,4,2) = rho  * u3
c       A1(:,4,3) = zero
        A1(:,4,4) = rho  * u1
        A1(:,4,5) = drdT * u1 * u3
c
        A1(:,5,1) = u1 * e2p
        A1(:,5,2) = e3p + rho * u1 * u1
        A1(:,5,3) = rho * u1 * u2
        A1(:,5,4) = rho * u1 * u3
        A1(:,5,5) = u1 * e4p
c
        flops = flops + 35*npro
c
        A2(:,1,1) = drdp * u2
c       A2(:,1,2) = zero
        A2(:,1,3) = rho
c       A2(:,1,4) = zero
        A2(:,1,5) = drdT * u2
c
        A2(:,2,1) = drdp * u1 * u2 
        A2(:,2,2) = rho  * u2
        A2(:,2,3) = rho  * u1
c       A2(:,2,4) = zero
        A2(:,2,5) = drdT * u1 * u2
c
        A2(:,3,1) = drdp * u2 * u2 +1
c       A2(:,3,2) = zero
        A2(:,3,3) = two * rho  * u2
c       A2(:,3,4) = zero
        A2(:,3,5) = drdT * u2 * u2
c
        A2(:,4,1) = drdp * u2 * u3 
c       A2(:,4,2) = zero
        A2(:,4,3) = rho  * u3
        A2(:,4,4) = rho  * u2
        A2(:,4,5) = drdT * u2 * u3
c
        A2(:,5,1) = u2 * e2p
        A2(:,5,2) = rho * u1 * u2
        A2(:,5,3) = e3p + rho * u2 * u2
        A2(:,5,4) = rho * u2 * u3
        A2(:,5,5) = u2 * e4p
c
        flops = flops + 35*npro
c
        A3(:,1,1) = drdp * u3
c       A3(:,1,2) = zero
c       A3(:,1,3) = zero
        A3(:,1,4) = rho
        A3(:,1,5) = drdT * u3
c
        A3(:,2,1) = drdp * u1 * u3 
        A3(:,2,2) = rho  * u3
c       A3(:,2,3) = zero
        A3(:,2,4) = rho  * u1
        A3(:,2,5) = drdT * u1 * u3
c
        A3(:,3,1) = drdp * u3 * u2 
c       A3(:,3,2) = zero
        A3(:,3,3) = rho  * u3
        A3(:,3,4) = rho  * u2
        A3(:,3,5) = drdT * u3 * u2
c
        A3(:,4,1) = drdp * u3 * u3 +1
c       A3(:,4,2) = zero
c       A3(:,4,3) = zero
        A3(:,4,4) = two *rho  * u3
        A3(:,4,5) = drdT * u3 * u3
c
        A3(:,5,1) = u3 * e2p
        A3(:,5,2) = rho * u1 * u3
        A3(:,5,3) = rho * u2 * u3
        A3(:,5,4) = e3p + rho * u3 * u3
        A3(:,5,5) = u3 * e4p
c
        flops = flops + 35*npro
	ttim(21) = ttim(21) + secs(0.0)

c
c.... return
c
      if (idc .ne. 0) then
c.... for Discountinuity Capturing Term
c
c.... calculation of A0^DC matrix
c
c.... Ref P-163 of the handout
c
       s1 = one/(rho**2 * betaT * T)
       cv = cp - (alfap**2 * T/rho/betaT)
       A0DC(:,1) = (rho * betaT)**2*s1
       A0DC(:,2) = -rho*alfap*rho*betaT*s1
       A0DC(:,3) = rho/T
       A0DC(:,4) = (-rho*alfap)**2 * s1 + (rho*cv/T**2)
c
c.... calculation of A0^tilda^inverse matrix
c.... ref P-169 of the hand out 
c
       fact1 = one/(rho*cv*T**2)
       d = alfap*T/rho/betaT
       e1bar = h - rk
       e2bar = e1bar - d
       e3bar = e2bar - cv * T
       e4bar = e2bar - 2* cv * T
       e5bar = e1bar**2 - 2*e1bar*d + 2*rk*cv*T + cp*T/rho/betaT
       c1bar = u1**2 + cv * T      
       c2bar = u2**2 + cv * T
       c3bar = u3**2 + cv * T      
       u12 = u1 * u2
       u31 = u3 * u1
       u23 = u2 * u3
       A0inv( :,1) = e5bar*fact1
       A0inv( :,2) = c1bar*fact1
       A0inv( :,3) = c2bar*fact1
       A0inv( :,4) = c3bar*fact1
       A0inv( :,5) = 1*fact1
       A0inv( :,6) = u1*e3bar*fact1
       A0inv( :,7) = u2*e3bar*fact1
       A0inv( :,8) = u3*e3bar*fact1
       A0inv( :,9) = -e2bar*fact1
       A0inv(:,10) = u12*fact1
       A0inv(:,11) = u31*fact1
       A0inv(:,12) = -u1*fact1
       A0inv(:,13) = u23*fact1
       A0inv(:,14) = -u2*fact1
       A0inv(:,15) = -u3*fact1
c  
c.....calculation of dV/dY (derivative of entropy variables w.r.to primitive
c
      fact1 = 1/T
      fact2 = fact1/T
      dVdY( :,1) = fact1/rho                
      dVdY( :,2) = -fact1*u1                
      dVdY( :,3) =  fact1               
      dVdY( :,4) =  -fact1*u2                
      dVdY( :,5) =  zero                
      dVdY( :,6) =  fact1               
      dVdY( :,7) =  -fact1*u3                
      dVdY( :,8) =  zero               
      dVdY( :,9) =  zero               
      dVdY(:,10) =  fact1               
      dVdY(:,11) =  -(h-rk)*fact2               
      dVdY(:,12) =  -fact2*u1               
      dVdY(:,13) =  -fact2*u2                
      dVdY(:,14) =  -fact2*u3                
      dVdY(:,15) =   fact2 

      endif  !end of idc.ne.0               

        return
        end
c
c
        subroutine e3mtrxSclr (rho,   
     &                         u1,     u2,   u3,
     &                         A0t,    A1t,    
     &                         A2t,    A3t  )
c
c----------------------------------------------------------------------
c
c This routine sets up the necessary matrices at the integration point.
c 
c input:
c  rho   (npro)         : density
c  u1    (npro)         : x1-velocity component
c  u2    (npro)         : x2-velocity component
c  u3    (npro)         : x3-velocity component
c
c output:
c  A0t    (npro) : A_0 "matrix"   
c  A1t   (npro)  : A_1 "matrix"   
c  A2t   (npro)  : A_2 "matrix"     
c  A3t   (npro)  : A_3 "matrix"    
c
c Note: the definition of the matrices can be found in 
c       thesis by Hauke. 
c
c Zdenek Johan, Summer 1990.  (Modified from e2mtrx.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c Kenneth Jansen, Winter 1997 Primitive Variables
c----------------------------------------------------------------------
c
        include "common.h"
c
c
c  passed arrays
c
        dimension rho(npro),
     &            u1(npro),        u2(npro),
     &            u3(npro),                  
     &            A0t(npro),        
     &            A1t(npro),       A2t(npro),
     &            A3t(npro)
c
        if (iconvsclr.eq.2) then  !convective form
           A0t(:) = one
           A1t(:) = u1(:)
           A2t(:) = u2(:)
           A3t(:) = u3(:)
        else                    !conservative form
           A0t(:) = rho(:)
           A1t(:) = rho(:)*u1(:)
           A2t(:) = rho(:)*u2(:)
           A3t(:) = rho(:)*u3(:)
        endif

c
c.... return
c
        return
        end


