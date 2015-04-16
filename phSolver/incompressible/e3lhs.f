      subroutine e3LHS ( u1,        u2,         u3,
     &                   uBar,      WdetJ,      rho,
     &                   rLui,      rmu,       
     &                   tauC,      tauM,       tauBar,
     &                   shpfun,    shg,        xKebe,
     &                   xGoC )
c------------------------------------------------------------------------
c 
c  This routine computes the left hand side tangent matrix at an 
c  integration point.
c
c  input:
c     u1(npro)                  : x1-velocity
c     u2(npro)                  : x2-velocity
c     u3(npro)                  : x3-velocity
c     uBar(npro,3)              : u - tauM * Li
c     WdetJ(npro)               : weighted jacobian determinant
c     rLui(npro,3)              : total residual of NS equations
c     rmu(npro)                 : fluid viscosity
c     rho(npro)				  : fluid density
c     tauC(npro)                : continuity tau
c     tauM(npro)                : momentum tau
c     tauBar(npro)              : additional tau
c     shpfun(npro,nshl)         : element shpfun functions
c     shg(npro,nshl,3)          : global grad of element shape functions
c
c  output:
c     xKebe(npro,9,nshl,nshl) : left hand side
c     xGoC(npro,4,nshl,nshl)    : left hand side
c
c
c------------------------------------------------------------------------
      include "common.h"

      dimension u1(npro),         u2(npro),       u3(npro),
     &          uBar(npro,3),     WdetJ(npro),    rho(npro),
     &          rLui(npro,3),     rmu(npro),   
     &          tauC(npro),       tauM(npro),     tauBar(npro),
     &          shpfun(npro,nshl),shg(npro,nshl,3)
      
      dimension xKebe(npro,9,nshl,nshl), xGoC(npro,4,nshl,nshl)
c
c.... local declarations
c
      dimension t1(npro,3),       t2(npro,3),      t3(npro,3),
     &          tmp1(npro),       tmp2(npro),    
     &          tmp(npro),        tlW(npro)

      integer   aa, b
      
      real*8    lhmFct, lhsFct,           tsFct(npro)
      
      lhsFct = alfi * gami * Delt(itseq)
      lhmFct = almi * (one - flmpl) 
c
c.... scale variables for efficiency
c
      tlW      = lhsFct * WdetJ     
      tmp1      = tlW * rho
      tauM      = tlW * tauM 
      tauC      = tlW * tauC 
      rmu       = tlW * rmu 
      tsFct     = lhmFct * WdetJ * rho
      if(iconvflow.eq.2) then  ! 2 is ubar form 3 is cons form but ubar tang. 
         tauBar    = lhsFct * WdetJ * tauBar 
         uBar(:,1) = tmp1 * uBar(:,1)
         uBar(:,2) = tmp1 * uBar(:,2)
         uBar(:,3) = tmp1 * uBar(:,3)
      else
         tauBar = zero  !lazy tangent...if effective code it
         uBar(:,1) = tmp1 * u1(:)
         uBar(:,2) = tmp1 * u2(:)
         uBar(:,3) = tmp1 * u3(:)
      endif

c
c.... compute mass and convection terms
c
      do b = 1, nshl
         t1(:,1) = uBar(:,1) * shg(:,b,1)
     &           + uBar(:,2) * shg(:,b,2)
     &           + uBar(:,3) * shg(:,b,3)
c
c t1=ubar_k N^b,k*rho*alpha_f*gamma*deltat*WdetJ  
c
         do aa = 1, nshl
            tmp1 = tsFct * shpfun(:,aa) * shpfun(:,b)
            tmp2 = tmp1 + t1(:,1) * shpfun(:,aa)
c
c tmp1=alpha_m*(1-lmp)*WdetJ*N^aN^b*rho   the time term CORRECT
c tmp2=tmp1+N^a*ubar_k N^b,k*rho*alpha_f*gamma*deltat*WdetJ   the 
c    second term is convective term CORRECT
c            
            xKebe(:,1,aa,b) = xKebe(:,1,aa,b) + tmp2
            xKebe(:,5,aa,b) = xKebe(:,5,aa,b) + tmp2
            xKebe(:,9,aa,b) = xKebe(:,9,aa,b) + tmp2
         enddo
      enddo
c
c.... compute the rest of K (symmetric terms)
c      
      do b = 1, nshl
         
         t1(:,1) = tauC * shg(:,b,1)
         t1(:,2) = tauC * shg(:,b,2)
         t1(:,3) = tauC * shg(:,b,3)

c t1 is tauC*N^b_i,j*alpha_f*gamma*deltat*WdetJ
         
         t2(:,1) = rmu  * shg(:,b,1)
         t2(:,2) = rmu  * shg(:,b,2)
         t2(:,3) = rmu  * shg(:,b,3)
c t2 is mu*N^b_j,k*alpha_f*gamma*deltat*WdetJ
      
         tmp1 = tauM   * ( u1 * shg(:,b,1)  
     &                   + u2 * shg(:,b,2) 
     &                   + u3 * shg(:,b,3) )*rho
c tmp1 is tauM*(rho u_m N^b_j,m)*alpha_f*gamma*deltat*WdetJ

         tmp2 = tauBar * ( rLui(:,1) * shg(:,b,1)
     &                   + rLui(:,2) * shg(:,b,2)
     &                   + rLui(:,3) * shg(:,b,3) )
c tmp2 is taubar*(L_m N^b_j,m)*alpha_f*gamma*deltat*WdetJ
         t3(:,1) = t2(:,1) + tmp1 * u1 + tmp2 * rLui(:,1)
         t3(:,2) = t2(:,2) + tmp1 * u2 + tmp2 * rLui(:,2)
         t3(:,3) = t2(:,3) + tmp1 * u3 + tmp2 * rLui(:,3)

c t3 is   (mu*N^b_j,k + u_k tauM*(rho u_m N^b_j,m)+ L_k*taubar*(L_mN^b_j,m ) 
c   *alpha_f*gamma*deltat*WdetJ     which isline 2 page 40 of whiting
c   ALMOST (waiting to get hit with N^a_{i,k}
c mu correct NOW (wrong before) and rho weight on tauM term
c
c.... first do the (nodal) diagonal blocks         
c
         aa  = b
         
         tmp = t3(:,1) * shg(:,aa,1)
     &       + t3(:,2) * shg(:,aa,2)
     &       + t3(:,3) * shg(:,aa,3)
c previous command is the N^a_{i,k} dot product with t3 defined above

         xKebe(:,1,aa,b) = xKebe(:,1,aa,b) + tmp
     &                      + t1(:,1) * shg(:,aa,1)
     &                      + t2(:,1) * shg(:,aa,1)
         xKebe(:,5,aa,b) = xKebe(:,5,aa,b) + tmp
     &                      + t1(:,2) * shg(:,aa,2)
     &                      + t2(:,2) * shg(:,aa,2)
         xKebe(:,9,aa,b) = xKebe(:,9,aa,b) + tmp
     &                      + t1(:,3) * shg(:,aa,3)
     &                      + t2(:,3) * shg(:,aa,3)
c
         tmp1               = t1(:,1) * shg(:,aa,2)
     &                      + t2(:,2) * shg(:,aa,1)
         xKebe(:,2,aa,b) = xKebe(:,2,aa,b) + tmp1 
         xKebe(:,4,b,aa) = xKebe(:,4,b,aa) + tmp1 
c
         tmp1               = t1(:,1) * shg(:,aa,3)
     &                      + t2(:,3) * shg(:,aa,1)
         xKebe(:,3,aa,b) = xKebe(:,3,aa,b) + tmp1 
         xKebe(:,7,b,aa) = xKebe(:,7,b,aa) + tmp1 
c
         tmp1               = t1(:,2) * shg(:,aa,3)
     &                      + t2(:,3) * shg(:,aa,2)
         xKebe(:,6,aa,b) = xKebe(:,6,aa,b) + tmp1 
         xKebe(:,8,b,aa) = xKebe(:,8,b,aa) + tmp1 
c
c.... now the off-diagonal (nodal) blocks
c
         do aa = b+1, nshl
            tmp             = t3(:,1) * shg(:,aa,1)
     &                      + t3(:,2) * shg(:,aa,2)
     &                      + t3(:,3) * shg(:,aa,3)
c
            tmp1            = tmp
     &                      + t1(:,1) * shg(:,aa,1)
     &                      + t2(:,1) * shg(:,aa,1)
            xKebe(:,1,aa,b) = xKebe(:,1,aa,b) + tmp1
            xKebe(:,1,b,aa) = xKebe(:,1,b,aa) + tmp1
c
            tmp1            = tmp
     &                      + t1(:,2) * shg(:,aa,2)
     &                      + t2(:,2) * shg(:,aa,2)
            xKebe(:,5,aa,b) = xKebe(:,5,aa,b) + tmp1
            xKebe(:,5,b,aa) = xKebe(:,5,b,aa) + tmp1
c
            tmp1            = tmp
     &                      + t1(:,3) * shg(:,aa,3)
     &                      + t2(:,3) * shg(:,aa,3)
            xKebe(:,9,aa,b) = xKebe(:,9,aa,b) + tmp1
            xKebe(:,9,b,aa) = xKebe(:,9,b,aa) + tmp1
c
c.... ( i != j )
c
            tmp1               = t1(:,1) * shg(:,aa,2)
     &                         + t2(:,2) * shg(:,aa,1)
            xKebe(:,2,aa,b) = xKebe(:,2,aa,b) + tmp1
            xKebe(:,4,b,aa) = xKebe(:,4,b,aa) + tmp1
c
            tmp1               = t1(:,1) * shg(:,aa,3)
     &                         + t2(:,3) * shg(:,aa,1)
            xKebe(:,3,aa,b) = xKebe(:,3,aa,b) + tmp1
            xKebe(:,7,b,aa) = xKebe(:,7,b,aa) + tmp1
c
            tmp1               = t1(:,2) * shg(:,aa,1)
     &                         + t2(:,1) * shg(:,aa,2)
            xKebe(:,4,aa,b) = xKebe(:,4,aa,b) + tmp1
            xKebe(:,2,b,aa) = xKebe(:,2,b,aa) + tmp1
c
            tmp1               = t1(:,2) * shg(:,aa,3)
     &                         + t2(:,3) * shg(:,aa,2)
            xKebe(:,6,aa,b) = xKebe(:,6,aa,b) + tmp1
            xKebe(:,8,b,aa) = xKebe(:,8,b,aa) + tmp1
c
            tmp1               = t1(:,3) * shg(:,aa,1)
     &                         + t2(:,1) * shg(:,aa,3)
            xKebe(:,7,aa,b) = xKebe(:,7,aa,b) + tmp1
            xKebe(:,3,b,aa) = xKebe(:,3,b,aa) + tmp1
c
            tmp1               = t1(:,3) * shg(:,aa,2)
     &                         + t2(:,2) * shg(:,aa,3)
            xKebe(:,8,aa,b) = xKebe(:,8,aa,b) + tmp1
            xKebe(:,6,b,aa) = xKebe(:,6,b,aa) + tmp1
c
         enddo
      enddo
c
c.... compute G   Nai Nbp,j
c
      
      do b = 1, nshl
         t1(:,1) = tlW * shg(:,b,1)
         t1(:,2) = tlW * shg(:,b,2)
         t1(:,3) = tlW * shg(:,b,3)
         do aa = 1, nshl
            xGoC(:,1,aa,b) = xGoC(:,1,aa,b) + t1(:,1) * shpfun(:,aa)  
            xGoC(:,2,aa,b) = xGoC(:,2,aa,b) + t1(:,2) * shpfun(:,aa)  
            xGoC(:,3,aa,b) = xGoC(:,3,aa,b) + t1(:,3) * shpfun(:,aa)  
         enddo
      enddo
c
c.... compute C
c we divide by rho because the L on the weight space is density divided
c      form
c
      tauM=tauM/rho
      do b = 1, nshl
         t1(:,1) = tauM * shg(:,b,1)
         t1(:,2) = tauM * shg(:,b,2)
         t1(:,3) = tauM * shg(:,b,3)
         do aa = b, nshl
            xGoC(:,4,aa,b) = xGoC(:,4,aa,b) 
     &                      + t1(:,1) * shg(:,aa,1)
     &                      + t1(:,2) * shg(:,aa,2)
     &                      + t1(:,3) * shg(:,aa,3)
         enddo
      enddo
      
c
c.... return
c
      return
      end


c------------------------------------------------------------------------
c
c     calculate the tangent matrix for the advection-diffusion equation
c
c------------------------------------------------------------------------
      subroutine e3LHSSclr ( uMod,      giju,       dcFct,
     &                       Sclr,      Sdot,       gradS,  
     &                       WdetJ,     rLS,        tauS,
     &                       shpfun,    shg,        srcL,
     &                       diffus,
     &                       xSebe )

c
      include "common.h"

      real*8    uMod(npro,nsd),
     &          Sclr(npro),       Sdot(npro),   gradS(npro,nsd),
     &          WdetJ(npro),      rLS(npro),        rho(npro), 
     &          tauS(npro),       shpfun(npro,nshl),  
     &          srcL(npro),        shg(npro,nshl,3),
     &			xSebe(npro,nshl,nshl)
      
      real*8    diffus(npro),  cp,  kptmp(npro),tauSo(npro)

c
c.... local declarations
c
      dimension t1(npro,3),       tmp1(npro),       tmp2(npro),
     &          tmp(npro),        dcFct(npro),      giju(npro,6)

      integer   aa, b
      
      real*8    lhsFct,           tsFct(npro)
      
      lhsFct = alfi * gami * Delt(itseq)
c
c.... scale variables for efficiency
c     
      tauSo     = tauS
      tauS      = lhsFct * WdetJ * tauS 
      kptmp     = lhsFct * WdetJ * diffus
      tsFct     = almi   * WdetJ * (one - flmpl)
      srcL       = srcL    * WdetJ * lhsFct
c
c.... compute mass and convection terms
c
      do b = 1, nshl
         t1(:,1) = WdetJ * ( uMod(:,1) * shg(:,b,1)
     &                     + uMod(:,2) * shg(:,b,2)
     &                     + uMod(:,3) * shg(:,b,3) )
         t1(:,2) = t1(:,1) * tauSo
         do aa = 1, nshl
            tmp1 = shpfun(:,aa) * shpfun(:,b)
            tmp2 = shpfun(:,aa) * lhsFct
            xSebe(:,aa,b) = xSebe(:,aa,b) + tmp1 * (tsFct + srcL)
     &                                    + tmp2 * t1(:,1)
c
c.... compute mass term for stab u_j N_{a,j} tau N_b (note that a and b
c            flipped on both sides below)
c
            xSebe(:,b,aa) = xSebe(:,b,aa) + t1(:,2)*shpfun(:,aa)
         enddo
      enddo
c
c.... compute the rest of S (symmetric terms)
c      
      do b = 1, nshl
         tmp     = tauS(:) 
     &             * ( uMod(:,1) * shg(:,b,1)
     &               + uMod(:,2) * shg(:,b,2)
     &               + uMod(:,3) * shg(:,b,3) )

         t1(:,1) = kptmp * shg(:,b,1) + uMod(:,1) * tmp
         t1(:,2) = kptmp * shg(:,b,2) + uMod(:,2) * tmp
         t1(:,3) = kptmp * shg(:,b,3) + uMod(:,3) * tmp
         if (idcsclr(1) .ne. 0) then
            if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or. 
     &           (idcsclr(2).eq.2 .and. isclr.eq.2)) then ! scalar with dc
c
               tmp = WdetJ * dcFct * lhsFct
c
               giju(:,1)	= tmp * giju(:,1)
               giju(:,2)	= tmp * giju(:,2)
               giju(:,3)	= tmp * giju(:,3)
               giju(:,4)	= tmp * giju(:,4)
               giju(:,5)	= tmp * giju(:,5)
               giju(:,6)	= tmp * giju(:,6)
c       
               t1(:,1) = t1(:,1) + giju(:,1) * shg(:,b,1) 
     2                           + giju(:,4) * shg(:,b,2) 
     3			         + giju(:,6) * shg(:,b,3)
               t1(:,2) = t1(:,2) + giju(:,4) * shg(:,b,1) 
     2                           + giju(:,2) * shg(:,b,2) 
     3			         + giju(:,5) * shg(:,b,3)
               t1(:,3) = t1(:,3) + giju(:,6) * shg(:,b,1) 
     2                           + giju(:,5) * shg(:,b,2) 
     3			         + giju(:,3) * shg(:,b,3)
            endif
         endif                  !end of idcsclr
c
c.... first do the (nodal) diagonal blocks         
c
         aa  = b
         
         xSebe(:,aa,b) = xSebe(:,aa,b) + t1(:,1) * shg(:,aa,1)
     &                                 + t1(:,2) * shg(:,aa,2)
     &                                 + t1(:,3) * shg(:,aa,3)

c
c.... now the off-diagonal (nodal) blocks
c
         do aa = b+1, nshl
            tmp             = t1(:,1) * shg(:,aa,1)
     &                      + t1(:,2) * shg(:,aa,2)
     &                      + t1(:,3) * shg(:,aa,3)
c
            xSebe(:,aa,b) = xSebe(:,aa,b) + tmp
            xSebe(:,b,aa) = xSebe(:,b,aa) + tmp
c
         enddo
      enddo
      
c
c.... return
c
      return
      end

