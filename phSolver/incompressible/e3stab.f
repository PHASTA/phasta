      subroutine e3stab (rho,          u1,       u2,
     &                   u3,           dxidx,    rLui,   
     &                   rmu,          tauC,     tauM,   
     &                   tauBar,       uBar )  
c
c----------------------------------------------------------------------
c
c This routine computes the diagonal Tau for least-squares operator.  
c Diagonal tau proposed by Shakib.
c
c input:
c  u1     (npro)           : x1-velocity component
c  u2     (npro)           : x2-velocity component
c  u3     (npro)           : x3-velocity component
c  dxidx  (npro,nsd,nsd)   : inverse of deformation gradient
c  rLui   (npro,nsd)      : least-squares residual vector
c
c output:
c  tauC    (npro)          : continuity tau
c  tauM    (npro)          : momentum tau
c  tauBar  (npro)          : additional tau
c  uBar    (npro,nsd)      : modified velocity
c
c Zdenek Johan, Summer 1990.  (Modified from e2tau.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension rho(npro),                 u1(npro),
     &            u2(npro),                  u3(npro),
     &            dxidx(npro,nsd,nsd), 
     &            rLui(npro,nsd),
     &            tauC(npro),    tauM(npro), tauBar(npro),
     &            rmu(npro),     uBar(npro,3), unorm(npro)

c
        dimension gijd(npro,6),       fact(npro), rnu(npro),
     &       rhoinv(npro)
c
c
c.... get the metric tensor
c      
      call e3gijd( dxidx, gijd )
c
c... higher order element diffusive correction
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

      omegasq=zero
      if(matflg(6,1).eq.1) omegasq = datmat(1,6,1)**2
     .                              +datmat(2,6,1)**2
     .                              +datmat(3,6,1)**2
      rhoinv=one/rho
      rnu=rmu*rhoinv

      if(itau.eq.0)  then  ! original tau
c
c...  momentum tau
c

!MR CHANGE 
         dts=  Dtgl*dtsfct ! Dtgl = (time step)^-1
!         dts=  min(Dtgl,28800.0d0)*dtsfct      ! Dtgl = (time step)^-1 !28800 = 1600*180 / 10
!         dts=  min(Dtgl,2880.0d0)*dtsfct      ! Dtgl = (time step)^-1 !2880 = 1600*180 / 100
!         dts=  min(Dtgl,288.0d0)*dtsfct      ! Dtgl = (time step)^-1 !288 = 1600*180 / 1000
!MR CHANGE

         tauM = ( (two*dts)**2
     3        + ( u1 * ( gijd(:,1) * u1
     4                       + gijd(:,4) * u2
     5                       + gijd(:,6) * u3 )
     6          + u2 * ( gijd(:,4) * u1
     7                       + gijd(:,2) * u2
     8			             + gijd(:,5) * u3 )
     9		        + u3 * ( gijd(:,6) * u1
     a			             + gijd(:,5) * u2
     1			             + gijd(:,3) * u3 ) ) )
     2		    + fff * rnu** 2
     3		    * ( gijd(:,1) ** 2
     4		      + gijd(:,2) ** 2
     5		      + gijd(:,3) ** 2
     6		      + 2.
     7		      * ( gijd(:,4) ** 2
     8		        + gijd(:,5) ** 2
     9		        + gijd(:,6) ** 2 ) 
     b              +omegasq)
        
         fact = sqrt(tauM)
         dtsi=one/dts
         ff=taucfct/dtsfct
         tauC =rho* pt125*fact/(gijd(:,1)+gijd(:,2)+gijd(:,3))*ff
         tauM = one/fact
      else if(itau.eq.1)  then  ! new tau

c
c  determinant of gijd
c
         fact = gijd(:,1) * gijd(:,2) * gijd(:,3)
     &        - gijd(:,2) * gijd(:,6) * gijd(:,6)
     &        - gijd(:,1) * gijd(:,5) * gijd(:,5)
     &        - gijd(:,3) * gijd(:,4) * gijd(:,4)
     &        + gijd(:,6) * gijd(:,4) * gijd(:,5) * two
        
c
c put 1/2u*h1 = sqrt(u_i g^{ij} u_j) into tau_M  note inverse is calculated
c on the fly here from cofactors over the determinent dotted from left and 
c right with u
c
         
         tauM = 
     1       u1 * ( (gijd(:,2)*gijd(:,3)-gijd(:,5)*gijd(:,5))  * u1
     2     +  two * (gijd(:,5)*gijd(:,6)-gijd(:,4)*gijd(:,3))  * u2 
     3     +  two * (gijd(:,4)*gijd(:,5)-gijd(:,6)*gijd(:,2))  * u3)
     1     + u2 * ( (gijd(:,1)*gijd(:,3)-gijd(:,6)*gijd(:,6))  * u2
     3     +  two * (gijd(:,4)*gijd(:,6)-gijd(:,1)*gijd(:,5))  * u3)
     1     + u3 * ( (gijd(:,1)*gijd(:,2)-gijd(:,4)*gijd(:,4))  * u3)
         tauM=fact/taum  ! here we have (u_i g^{ij} u^j)^{-1} approx 4/u^2h^2
c
c  we can calculate tauC more efficiently now
c
         tauC=tauM*(one+tauM*rmu*rmu)
         tauC=one/tauC
         tauC=taucfct*sqrt(tauC)
c
c
c...  momentum tau
c
c
c     this tau needs a u/h instead of a u*h so we contract with g_{ij} as
c     follows  (i.e. u_i g_{ij} u_j approx u^2/(h^2)/4) 
c
         fact = 
     3          u1 * ( gijd(:,1) * u1
     4               + gijd(:,4) * u2
     5               + gijd(:,6) * u3 )
     6        + u2 * ( gijd(:,4) * u1
     7               + gijd(:,2) * u2
     8               + gijd(:,5) * u3 )
     9        + u3 * ( gijd(:,6) * u1
     a               + gijd(:,5) * u2
     1               + gijd(:,3) * u3 ) 
c 
c first limit dt effect on tau from causing trouble if user drops CFL below
c .05 (this could cause loss of spatial stability)
c
         velsq=vel*vel
         unorm = (u1*u1+u2*u2+u3*u3)/velsq
         dtsfsq=dtsfct*dtsfct
         dt=one/Dtgl
         taubar=  dtsfsq/( dt*dt + .01*unorm/fact)  ! never gets above (C_1 20*u_inf/h)^2
c
c  this means tau will never get below h/(20*C_1*u) no matter what time step 
c  you choose.  The 0.01 constant comes from minCFL=.05=> .05*.05*4 (where the 
c  4 comes from the bi-unit mapping). If you want to limit sooner the formula
c  would be  ".01-factor"=minCFL^2*4
c

         tauM = rho ** 2
     1		    * ( four*taubar + fact
     2		    + fff * rmu** 2
     3		    * ( gijd(:,1) ** 2
     4		      + gijd(:,2) ** 2
     5		      + gijd(:,3) ** 2
     6		      + 2.
     7		      * ( gijd(:,4) ** 2
     8		        + gijd(:,5) ** 2
     9		        + gijd(:,6) ** 2 ) ) 
     b              +omegasq)
         fact=sqrt(tauM)
cdebugcheck         tauBar = pt125*fact/(gijd(:,1)+gijd(:,2)+gijd(:,3)) !*dtsi
      
        tauM=one/fact           ! turn it right side up.
      else if(itau.eq.2)  then  ! new tau different continuity h

         unorm = (u1*u1+u2*u2+u3*u3)
         
         tauM=(gijd(:,1)+gijd(:,2)+gijd(:,3))/unorm ! here we have  4/u^2h^2
c
c  we can calculate tauC more efficiently now
c
         tauC=tauM*(one+tauM*rmu*rmu)
         tauC=one/tauC
         tauC=sqrt(tauC)*taucfct
c
c
c...  momentum tau
c
c
c     this tau needs a u/h instead of a u*h so we contract with g_{ij} as
c     follows  (i.e. u_i g_{ij} u_j approx u^2/(h^2)/4) 
c
         fact = 
     3          u1 * ( gijd(:,1) * u1
     4               + gijd(:,4) * u2
     5               + gijd(:,6) * u3 )
     6        + u2 * ( gijd(:,4) * u1
     7               + gijd(:,2) * u2
     8               + gijd(:,5) * u3 )
     9        + u3 * ( gijd(:,6) * u1
     a               + gijd(:,5) * u2
     1               + gijd(:,3) * u3 ) 
c 
c first limit dt effect on tau from causing trouble if user drops CFL below
c .05 (this could cause loss of spatial stability)
c
         velsq=vel*vel
         dtsfsq=dtsfct*dtsfct
         dt=one/Dtgl
         unorm=unorm/velsq
         taubar=  dtsfsq/( dt*dt + .01*unorm/fact)  ! never gets above (C_1 20*u_inf/h)^2
c
c  this means tau will never get below h/(20*C_1*u) no matter what time step 
c  you choose.  The 0.01 constant comes from minCFL=.05=> .05*.05*4 (where the 
c  4 comes from the bi-unit mapping). If you want to limit sooner the formula
c  would be  ".01-factor"=minCFL^2*4
c

         tauM = rho ** 2
     1		    * ( four*taubar + fact
     2		    + fff * rmu** 2
     3		    * ( gijd(:,1) ** 2
     4		      + gijd(:,2) ** 2
     5		      + gijd(:,3) ** 2
     6		      + 2.
     7		      * ( gijd(:,4) ** 2
     8		        + gijd(:,5) ** 2
     9		        + gijd(:,6) ** 2 ) ) 
     b              +omegasq)
         fact=sqrt(tauM)
c         tauBar = pt125*fact/(gijd(:,1)+gijd(:,2)+gijd(:,3)) !*dtsi
      
        tauM=one/fact           ! turn it right side up.
      else if(itau.eq.3)  then  ! compressible tau

c
c  determinant of gijd
c
         fact = gijd(:,1) * gijd(:,2) * gijd(:,3)
     &        - gijd(:,2) * gijd(:,6) * gijd(:,6)
     &        - gijd(:,1) * gijd(:,5) * gijd(:,5)
     &        - gijd(:,3) * gijd(:,4) * gijd(:,4)
     &        + gijd(:,6) * gijd(:,4) * gijd(:,5) * two
        
c
c put 1/2u*h1 = sqrt(u_i g^{ij} u_j) into tau_M  note inverse is calculated
c on the fly here from cofactors over the determinent dotted from left and 
c right with u
c
         
         tauM = 
     1       u1 * ( (gijd(:,2)*gijd(:,3)-gijd(:,5)*gijd(:,5))  * u1
     2     +  two * (gijd(:,5)*gijd(:,6)-gijd(:,4)*gijd(:,3))  * u2 
     3     +  two * (gijd(:,4)*gijd(:,5)-gijd(:,6)*gijd(:,2))  * u3)
     1     + u2 * ( (gijd(:,1)*gijd(:,3)-gijd(:,6)*gijd(:,6))  * u2
     3     +  two * (gijd(:,4)*gijd(:,6)-gijd(:,1)*gijd(:,5))  * u3)
     1     + u3 * ( (gijd(:,1)*gijd(:,2)-gijd(:,4)*gijd(:,4))  * u3)
c
c  we can calculate tauC more efficiently now
c
         tauM=sqrt(tauM/fact)*two
         tauC=pt5*tauM*min(one,pt5*tauM/rmu)*taucfct
c
c
c...  momentum tau
c
c
c     this tau needs a u/h instead of a u*h so we contract with g_{ij} as
c     follows  (i.e. u_i g_{ij} u_j approx u^2/(h^2)/4) 
c
         fact = 
     3          u1 * ( gijd(:,1) * u1
     4               + gijd(:,4) * u2
     5               + gijd(:,6) * u3 )
     6        + u2 * ( gijd(:,4) * u1
     7               + gijd(:,2) * u2
     8               + gijd(:,5) * u3 )
     9        + u3 * ( gijd(:,6) * u1
     a               + gijd(:,5) * u2
     1               + gijd(:,3) * u3 ) 
         fact=one/sqrt(fact)

         unorm = (u1*u1+u2*u2+u3*u3)

         dts= one/( Dtgl*dtsfct)
         tauM =min(dts,min(fact,fact*fact*unorm*pt33/rmu))
      endif
c
c.... calculate tauBar
c
      tauBar = rLui(:,1) * ( gijd(:,1) * rLui(:,1)
     &                       + gijd(:,4) * rLui(:,2)
     &                       + gijd(:,6) * rLui(:,3) )
     &         + rLui(:,2) * ( gijd(:,4) * rLui(:,1)
     &                       + gijd(:,2) * rLui(:,2)
     &                       + gijd(:,5) * rLui(:,3) ) 
     &         + rLui(:,3) * ( gijd(:,6) * rLui(:,1)
     &                       + gijd(:,5) * rLui(:,2)
     &                       + gijd(:,3) * rLui(:,3) )
      where ( tauBar .ne. 0.0 ) 
         tauBar = tauM / sqrt(tauBar)
      endwhere

c
c.... compute the modified velocity, uBar
c
        uBar(:,1) = u1 - tauM * rLui(:,1)*rhoinv
        uBar(:,2) = u2 - tauM * rLui(:,2)*rhoinv
        uBar(:,3) = u3 - tauM * rLui(:,3)*rhoinv
c     
c.... return
c
        return
        end

c-----------------------------------------------------------------------
c
c  Momentum tau
c
c-----------------------------------------------------------------------
      subroutine e3uBar (rho,          ui,         dxidx,     
     &                   rLui,         rmu,        uBar )         

      include "common.h"

      real*8     rho(npro),            ui(npro,nsd),
     &           dxidx(npro,nsd,nsd),  rLui(npro,nsd),
     &           rmu(npro),            uBar(npro,nsd)

      real*8     gijd(npro,6),         tauM(npro)

c
c.... get the metric tensor
c      
      call e3gijd( dxidx, gijd )
c
c.... higher order element diffusive correction
c
      if (ipord == 1) then
         fff = 36.0d0
      else if (ipord == 2) then
         fff = 60.0d0
      else if (ipord == 3) then
         fff = 128.0d0
      endif

!MR CHANGE 
      dts  =  (Dtgl*dtsfct)
!      dts=  min(Dtgl,28800.0d0)*dtsfct      !28800 = 1600*180 / 10
!      dts=  min(Dtgl,2880.0d0)*dtsfct      !2880 = 1600*180 / 100
!      dts=  min(Dtgl,288.0d0)*dtsfct        !288 = 1600*180 / 1000
!MR CHANGE

      tauM = rho ** 2
     1              * ( (two*dts)**2
     3                + ( ui(:,1) * ( gijd(:,1) * ui(:,1)
     4                              + gijd(:,4) * ui(:,2)
     5                              + gijd(:,6) * ui(:,3) )
     6                  + ui(:,2) * ( gijd(:,4) * ui(:,1)
     7			            + gijd(:,2) * ui(:,2)
     8			            + gijd(:,5) * ui(:,3) )
     9		        + ui(:,3) * ( gijd(:,6) * ui(:,1)
     a			            + gijd(:,5) * ui(:,2)
     1			            + gijd(:,3) * ui(:,3) ) ) )
     2		    + fff * rmu** 2
     3		    * ( gijd(:,1) ** 2
     4		      + gijd(:,2) ** 2
     5		      + gijd(:,3) ** 2
     6		      + 2.
     7		      * ( gijd(:,4) ** 2
     8		        + gijd(:,5) ** 2
     9		        + gijd(:,6) ** 2 ) )
        
      tauM = one/sqrt(tauM)
c
c.... compute the modified velocity, uBar
c
      uBar(:,1) = ui(:,1) - tauM * rLui(:,1)
      uBar(:,2) = ui(:,2) - tauM * rLui(:,2)
      uBar(:,3) = ui(:,3) - tauM * rLui(:,3)

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
c
c  form metric tensor g_{ij}=xi_{k,i} xi_{k,j}.  It is a symmetric
c  tensor so we only form 6 components and use symmetric matrix numbering.
c
      if (lcsyst .ge. 2) then  ! note this makes wedges like hexs..should
c                                be corrected later

         gijd(:,1) = dxidx(:,1,1) * dxidx(:,1,1)
     &             + dxidx(:,2,1) * dxidx(:,2,1)
     &             + dxidx(:,3,1) * dxidx(:,3,1)
c
         gijd(:,4) = dxidx(:,1,1) * dxidx(:,1,2)
     &             + dxidx(:,2,1) * dxidx(:,2,2)
     &             + dxidx(:,3,1) * dxidx(:,3,2)
c
         gijd(:,2) = dxidx(:,1,2) * dxidx(:,1,2)
     &             + dxidx(:,2,2) * dxidx(:,2,2)
     &             + dxidx(:,3,2) * dxidx(:,3,2)
c
         gijd(:,5) = dxidx(:,1,2) * dxidx(:,1,3)
     &             + dxidx(:,2,2) * dxidx(:,2,3)
     &             + dxidx(:,3,2) * dxidx(:,3,3)
c
         gijd(:,6) = dxidx(:,1,1) * dxidx(:,1,3)
     &             + dxidx(:,2,1) * dxidx(:,2,3)
     &             + dxidx(:,3,1) * dxidx(:,3,3)
c
         gijd(:,3) = dxidx(:,1,3) * dxidx(:,1,3)
     &             + dxidx(:,2,3) * dxidx(:,2,3)
     &             + dxidx(:,3,3) * dxidx(:,3,3)
c
      else   if (lcsyst .eq. 1) then
c
c  There is an invariance problem with tets 
c  It is fixed by the following modifications to gijd 
c

         c1 = 1.259921049894873D+00
         c2 = 6.299605249474365D-01
c
         tmp1(:) = c1 * dxidx(:,1,1)+c2 *(dxidx(:,2,1)+dxidx(:,3,1))
         tmp2(:) = c1 * dxidx(:,2,1)+c2 *(dxidx(:,1,1)+dxidx(:,3,1))
         tmp3(:) = c1 * dxidx(:,3,1)+c2 *(dxidx(:,1,1)+dxidx(:,2,1))
         gijd(:,1) = dxidx(:,1,1) * tmp1
     1              + dxidx(:,2,1) * tmp2
     2              + dxidx(:,3,1) * tmp3
c
         tmp1(:) = c1 * dxidx(:,1,2)+c2 *(dxidx(:,2,2)+dxidx(:,3,2))
         tmp2(:) = c1 * dxidx(:,2,2)+c2 *(dxidx(:,1,2)+dxidx(:,3,2))
         tmp3(:) = c1 * dxidx(:,3,2)+c2 *(dxidx(:,1,2)+dxidx(:,2,2))
         gijd(:,2) = dxidx(:,1,2) * tmp1
     1             + dxidx(:,2,2) * tmp2
     2             + dxidx(:,3,2) * tmp3
c
         gijd(:,4) = dxidx(:,1,1) * tmp1
     1             + dxidx(:,2,1) * tmp2
     2             + dxidx(:,3,1) * tmp3
c
         tmp1(:) = c1 * dxidx(:,1,3)+c2 *(dxidx(:,2,3)+dxidx(:,3,3))
         tmp2(:) = c1 * dxidx(:,2,3)+c2 *(dxidx(:,1,3)+dxidx(:,3,3))
         tmp3(:) = c1 * dxidx(:,3,3)+c2 *(dxidx(:,1,3)+dxidx(:,2,3))
         gijd(:,3) = dxidx(:,1,3) * tmp1
     1             + dxidx(:,2,3) * tmp2
     2             + dxidx(:,3,3) * tmp3
c
         gijd(:,5) = dxidx(:,1,2) * tmp1
     1             + dxidx(:,2,2) * tmp2
     2             + dxidx(:,3,2) * tmp3
c
         gijd(:,6) = dxidx(:,1,1) * tmp1
     1             + dxidx(:,2,1) * tmp2
     2             + dxidx(:,3,1) * tmp3
c
      else
         write(*,*) 'lcsyst eq',lcsyst,'not supported'
         stop
      endif

      return
      end

c------------------------------------------------------------------------
c
c     calculate the stabilization for the advection-diffusion equation
c
c------------------------------------------------------------------------
      subroutine e3StabSclr (uMod,  dxidx,  tauT, diffus, srcR, giju,
     &                       srcRat )
c
c
        include "common.h"
c
        real*8    rho(npro),                 uMod(npro,nsd),
     &            dxidx(npro,nsd,nsd),       diffus(npro),
     &            tauT(npro),                srcR(npro)

c
        real*8    gijd(npro,6),       giju(npro,6),   
     &            tmp1(npro),         tmp2(npro),
     &            tmp3(npro),         fact(npro),
     &            srcRat(npro)

        real*8     fff
        if(ivart.eq.1) then
           tauT=zero
           return
        endif
c
c.... get the metric tensor
c      
      call e3gijd( dxidx, gijd )
c
c...  momentum tau
c 
c
c... higher order element diffusive correction
c
        if (ipord == 1) then
           fff = 9.0d0
        else if (ipord == 2) then
           fff = 36.0d0
        else if (ipord == 3) then
           fff = 64.0d0
        endif

!MR CHANGE 
      dts  =  (Dtgl*dtsfct)
!      dts=  min(Dtgl,28800.0d0)*dtsfct      !28800 = 1600*180 / 10
!      dts=  min(Dtgl,2880.0d0)*dtsfct      !2880 = 1600*180 / 100
!      dts=  min(Dtgl,288.0d0)*dtsfct        !288 = 1600*180 / 1000
!MR CHANGE

c        if(iRANS.ne.-2) srcRat=srcR
        tauT = 
     1         (two*dts)**2 
     2       + srcRat ** 2
     3       + uMod(:,1) * ( gijd(:,1) * uMod(:,1)
     4                     + gijd(:,4) * uMod(:,2)
     5	                   + gijd(:,6) * uMod(:,3) )
     6	     + uMod(:,2) * ( gijd(:,4) * uMod(:,1)
     7	                   + gijd(:,2) * uMod(:,2)
     8	                   + gijd(:,5) * uMod(:,3) )
     9	     + uMod(:,3) * ( gijd(:,6) * uMod(:,1)
     a	                   + gijd(:,5) * uMod(:,2)
     1	                   + gijd(:,3) * uMod(:,3) )
     2	     + fff * diffus(:)** 2
     3	           * ( gijd(:,1) ** 2
     4		     + gijd(:,2) ** 2
     5		     + gijd(:,3) ** 2
     6		     + 2.
     7		      * ( gijd(:,4) ** 2
     8		        + gijd(:,5) ** 2
     9		        + gijd(:,6) ** 2 ) )
        
        tauT = one/sqrt(tauT)
c
        if(idcsclr(1) .ne. 0) then 
           if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or. 
     &          (idcsclr(2).eq.2 .and. isclr.eq.2)) then ! scalar with dc
c     
c     determinant of gijd
c     
              fact = one/(gijd(:,1) * gijd(:,2) * gijd(:,3)
     &             - gijd(:,2) * gijd(:,6) * gijd(:,6)
     &             - gijd(:,1) * gijd(:,5) * gijd(:,5)
     &             - gijd(:,3) * gijd(:,4) * gijd(:,4)
     &             + gijd(:,6) * gijd(:,4) * gijd(:,5) * two)
c
c ... note between compressible and incompressible 5 and 6 of giju 
c     are switched        
c
              giju(:,1) = fact * (gijd(:,2)*gijd(:,3) 
     &                  - gijd(:,5)**2)
              giju(:,2) = fact * (gijd(:,1)*gijd(:,3) 
     &                  - gijd(:,6)**2)
              giju(:,3) = fact * (gijd(:,1)*gijd(:,2)
     &                  - gijd(:,4)**2)
              giju(:,4) = fact * (gijd(:,5)*gijd(:,6)
     &                  - gijd(:,4)*gijd(:,3) )
              giju(:,5) = fact * (gijd(:,4)*gijd(:,6)
     &                  - gijd(:,1)*gijd(:,5) )
              giju(:,6) = fact * (gijd(:,4)*gijd(:,5)
     &                  - gijd(:,6)*gijd(:,2) )

c
           endif
        endif                   ! end of idcsclr.ne.0
c     
c.... return
c
        return
        end

