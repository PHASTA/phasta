        subroutine e3ivar (yl,      ycl,     acl,
     &                     Sclr,    shape,   shgl,    
     &                     xl,      dui,     aci,
     &                     g1yi,    g2yi,    g3yi,
     &                     shg,     dxidx,   WdetJ,   
     &                     rho,     pres,    T,
     &                     ei,      h,       alfap,   
     &                     betaT,   cp,      rk,      
     &                     u1,      u2,      u3,
     &                     ql,      divqi,   sgn, tmp,
     &                     rmu,     rlm,     rlm2mu,
     &                     con,     rlsl,    rlsli, 
     &                     xmudmi,  sforce,  cv) 
c
c----------------------------------------------------------------------
c
c  This routine computes the variables at integration point.
c
c input:
c  yl     (npro,nshl,nflow)     : primitive variables (NO SCALARS)
c  ycl    (npro,nshl,ndof)      : primitive variables at current step
c  acl    (npro,nshl,ndof)      : prim.var. accel. at the current step
c  shape  (npro,nshl)           : element shape-functions
c  shgl   (nsd,nen)             : element local-grad-shape-functions
c  xl     (npro,nenl,nsd)       : nodal coordinates at current step
c  ql     (npro,nshl,(nflow-1)*nsd) : diffusive flux vector
c  rlsl   (npro,nshl,6)       : resolved Leonard stresses
c  sgn    (npro,nshl)         : signs of shape functions
c
c output:
c  dui    (npro,nflow)           : delta U variables at current step
c  aci    (npro,nflow)           : primvar accel. variables at current step
c  g1yi   (npro,nflow)           : grad-y in direction 1
c  g2yi   (npro,nflow)           : grad-y in direction 2
c  g3yi   (npro,nflow)           : grad-y in direction 3
c  shg    (npro,nshl,nsd)       : element global grad-shape-functions
c  dxidx  (npro,nsd,nsd)        : inverse of deformation gradient
c  WdetJ  (npro)                : weighted Jacobian
c  rho    (npro)                : density
c  pres   (npro)                : pressure
c  T      (npro)                : temperature
c  ei     (npro)                : internal energy
c  h      (npro)                : enthalpy
c  alfap  (npro)                : expansivity
c  betaT  (npro)                : isothermal compressibility
c  cp     (npro)                : specific heat at constant pressure
c  rk     (npro)                : kinetic energy
c  u1     (npro)                : x1-velocity component
c  u2     (npro)                : x2-velocity component
c  u3     (npro)                : x3-velocity component
c  divqi  (npro,nflow-1)        : divergence of diffusive flux
c  rlsli  (npro,6)              : resolved Leonard stresses at quad pt
c
c Zdenek Johan, Summer 1990. (Modified from e2ivar.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997. Primitive Variables
c----------------------------------------------------------------------
c
        include "common.h"
c
c  passed arrays
c
        dimension yl(npro,nshl,nflow),        ycl(npro,nshl,ndof),
     &            acl(npro,nshl,ndof),      
     &            shape(npro,nshl),  
     &            shgl(npro,nsd,nshl),      xl(npro,nenl,nsd),
     &            dui(npro,nflow),
     &            aci(npro,nflow),            g1yi(npro,nflow),
     &            g2yi(npro,nflow),           g3yi(npro,nflow),
     &            shg(npro,nshl,nsd),        dxidx(npro,nsd,nsd),
     &            WdetJ(npro),
     &            rho(npro),                 pres(npro),
     &            T(npro),                   ei(npro),
     &            h(npro),                   alfap(npro),
     &            betaT(npro),               cp(npro),                  
     &            rk(npro),                  cv(npro),
     &            u1(npro),                  u2(npro),
     &            u3(npro),                  divqi(npro,nflow), 
     &            ql(npro,nshl,idflx),
     &            rmu(npro),                 rlm(npro),
     &            rlm2mu(npro),              con(npro), 
     &            Sclr(npro)
c
        dimension tmp(npro),  dxdxi(npro,nsd,nsd),  sgn(npro,nshl)
        dimension rlsl(npro,nshl,6),         rlsli(npro,6),
     &            xmudmi(npro,ngauss)
        dimension gyti(npro,nsd),            gradh(npro,nsd),
     &            sforce(npro,3),            weber(npro) 
        ttim(20) = ttim(20) - secs(0.0)

c
        ttim(10) = ttim(10) - secs(0.0)

        dui = zero
c
        do n = 1, nshl
           dui(:,1) = dui(:,1) + shape(:,n) * yl(:,n,1) ! p
           dui(:,2) = dui(:,2) + shape(:,n) * yl(:,n,2) ! u1
           dui(:,3) = dui(:,3) + shape(:,n) * yl(:,n,3) ! u2
           dui(:,4) = dui(:,4) + shape(:,n) * yl(:,n,4) ! u3
           dui(:,5) = dui(:,5) + shape(:,n) * yl(:,n,5) ! T
        enddo
c     
        flops = flops + 10*nshl*npro
c     
c     
c.... compute conservative variables
c
        rk = pt5 * (dui(:,2)**2 + dui(:,3)**2 + dui(:,4)**2)
c     
        if (iLSet .ne. 0)then
           Sclr = zero
           isc=abs(iRANS)+6
           do n = 1, nshl
              Sclr = Sclr + shape(:,n) * ycl(:,n,isc)
           enddo
        endif
        
        ithm = 6
        call getthm (dui(:,1),   dui(:,5),     Sclr,
     &               rk,         rho,          ei,
     &               tmp,        tmp,          tmp,
     &               tmp,        tmp,          tmp,
     &               tmp,        tmp)
c     
c
        dui(:,1) = rho
        dui(:,2) = rho * dui(:,2)
        dui(:,3) = rho * dui(:,3)
        dui(:,4) = rho * dui(:,4)
        dui(:,5) = rho * (ei + rk)
        
          
       ttim(10) = ttim(10) + secs(0.0)
c
c.... ------------->  Primitive variables at int. point  <--------------
c
c.... compute primitive variables
c
       ttim(11) = ttim(11) - secs(0.0)

       pres = zero
       u1   = zero
       u2   = zero
       u3   = zero
       T    = zero
       do n = 1, nshl 
c
c  y(int)=SUM_{a=1}^nshl (N_a(int) Ya)
c
          pres = pres + shape(:,n) * ycl(:,n,1)
          u1   = u1   + shape(:,n) * ycl(:,n,2)
          u2   = u2   + shape(:,n) * ycl(:,n,3)
          u3   = u3   + shape(:,n) * ycl(:,n,4)
          T    = T    + shape(:,n) * ycl(:,n,5)          
       enddo

       if( (iLES.gt.10).and.(iLES.lt.20))  then  ! bardina
       rlsli = zero
       do n = 1, nshl 

          rlsli(:,1) = rlsli(:,1) + shape(:,n) * rlsl(:,n,1)
          rlsli(:,2) = rlsli(:,2) + shape(:,n) * rlsl(:,n,2)
          rlsli(:,3) = rlsli(:,3) + shape(:,n) * rlsl(:,n,3)
          rlsli(:,4) = rlsli(:,4) + shape(:,n) * rlsl(:,n,4)
          rlsli(:,5) = rlsli(:,5) + shape(:,n) * rlsl(:,n,5)
          rlsli(:,6) = rlsli(:,6) + shape(:,n) * rlsl(:,n,6)

       enddo
       else
          rlsli = zero
       endif
c
       
       ttim(11) = ttim(11) + secs(0.0)

c
c.... ----------------------->  accel. at int. point  <----------------------
c
       
c       if (ires .ne. 2)  then
          ttim(12) = ttim(12) - secs(0.0)
c
c.... compute primitive variables
c
          aci = zero
c
          do n = 1, nshl
             aci(:,1) = aci(:,1) + shape(:,n) * acl(:,n,1)
             aci(:,2) = aci(:,2) + shape(:,n) * acl(:,n,2)
             aci(:,3) = aci(:,3) + shape(:,n) * acl(:,n,3)
             aci(:,4) = aci(:,4) + shape(:,n) * acl(:,n,4)
             aci(:,5) = aci(:,5) + shape(:,n) * acl(:,n,5)
          enddo
c
          flops = flops + 10*nshl*npro
          ttim(12) = ttim(12) + secs(0.0)
c       endif                    !ires .ne. 2
c
c.... ----------------->  Thermodynamic Properties  <-------------------
c
c.... compute the kinetic energy
c
        ttim(27) = ttim(27) - secs(0.0)
c
        rk = pt5 * ( u1**2 + u2**2 + u3**2 )
c
c.... get the thermodynamic properties
c
        if (iLSet .ne. 0)then
           Sclr = zero
           isc=abs(iRANS)+6
           do n = 1, nshl
              Sclr = Sclr + shape(:,n) * ycl(:,n,isc)
           enddo
        endif

        ithm = 7
        call getthm (pres,            T,               Sclr,
     &               rk,              rho,             ei,
     &               h,               tmp,             cv,
     &               cp,              alfap,           betaT,
     &               tmp,             tmp)
c
        ttim(27) = ttim(27) + secs(0.0)
c
c ........Get the material properties 
c
        call getDiff (T,        cp,       rho,        ycl, 
     &                rmu,      rlm,      rlm2mu,     con,  shape,
     &                xmudmi,   xl)

c.... --------------------->  Element Metrics  <-----------------------
c
      ttim(26) = ttim(26) - secs(0.0)
c
        call e3metric( xl,         shgl,        dxidx,  
     &                 shg,        WdetJ)
c
c       
        ttim(26) = ttim(26) + secs(0.0)
c
c.... --------------------->  Global Gradients  <-----------------------
c
        ttim(13) = ttim(13) - secs(0.0)
        
        g1yi = zero
        g2yi = zero
        g3yi = zero
c
        ttim(13) = ttim(13) + secs(0.0)
        ttim(7) = ttim(7) - secs(0.0)
c
c.... compute the global gradient of Y-variables
c
c
c  Y_{,x_i}=SUM_{a=1}^nshl (N_{a,x_i}(int) Ya)
c
        if(nshl.eq.4) then
          g1yi(:,1) = g1yi(:,1) + shg(:,1,1) * yl(:,1,1)
     &                          + shg(:,2,1) * yl(:,2,1)
     &                          + shg(:,3,1) * yl(:,3,1)
     &                          + shg(:,4,1) * yl(:,4,1)
c
          g1yi(:,2) = g1yi(:,2) + shg(:,1,1) * yl(:,1,2)
     &                          + shg(:,2,1) * yl(:,2,2)
     &                          + shg(:,3,1) * yl(:,3,2)
     &                          + shg(:,4,1) * yl(:,4,2)
c
          g1yi(:,3) = g1yi(:,3) + shg(:,1,1) * yl(:,1,3)
     &                          + shg(:,2,1) * yl(:,2,3)
     &                          + shg(:,3,1) * yl(:,3,3)
     &                          + shg(:,4,1) * yl(:,4,3)
c
          g1yi(:,4) = g1yi(:,4) + shg(:,1,1) * yl(:,1,4)
     &                          + shg(:,2,1) * yl(:,2,4)
     &                          + shg(:,3,1) * yl(:,3,4)
     &                          + shg(:,4,1) * yl(:,4,4)
c
          g1yi(:,5) = g1yi(:,5) + shg(:,1,1) * yl(:,1,5)
     &                          + shg(:,2,1) * yl(:,2,5)
     &                          + shg(:,3,1) * yl(:,3,5)
     &                          + shg(:,4,1) * yl(:,4,5)
c
          g2yi(:,1) = g2yi(:,1) + shg(:,1,2) * yl(:,1,1)
     &                          + shg(:,2,2) * yl(:,2,1)
     &                          + shg(:,3,2) * yl(:,3,1)
     &                          + shg(:,4,2) * yl(:,4,1)
c
          g2yi(:,2) = g2yi(:,2) + shg(:,1,2) * yl(:,1,2)
     &                          + shg(:,2,2) * yl(:,2,2)
     &                          + shg(:,3,2) * yl(:,3,2)
     &                          + shg(:,4,2) * yl(:,4,2)
c
          g2yi(:,3) = g2yi(:,3) + shg(:,1,2) * yl(:,1,3)
     &                          + shg(:,2,2) * yl(:,2,3)
     &                          + shg(:,3,2) * yl(:,3,3)
     &                          + shg(:,4,2) * yl(:,4,3)
c
          g2yi(:,4) = g2yi(:,4) + shg(:,1,2) * yl(:,1,4)
     &                          + shg(:,2,2) * yl(:,2,4)
     &                          + shg(:,3,2) * yl(:,3,4)
     &                          + shg(:,4,2) * yl(:,4,4)
c
          g2yi(:,5) = g2yi(:,5) + shg(:,1,2) * yl(:,1,5)
     &                          + shg(:,2,2) * yl(:,2,5)
     &                          + shg(:,3,2) * yl(:,3,5)
     &                          + shg(:,4,2) * yl(:,4,5)
c
          g3yi(:,1) = g3yi(:,1) + shg(:,1,3) * yl(:,1,1)
     &                          + shg(:,2,3) * yl(:,2,1)
     &                          + shg(:,3,3) * yl(:,3,1)
     &                          + shg(:,4,3) * yl(:,4,1)
c
          g3yi(:,2) = g3yi(:,2) + shg(:,1,3) * yl(:,1,2)
     &                          + shg(:,2,3) * yl(:,2,2)
     &                          + shg(:,3,3) * yl(:,3,2)
     &                          + shg(:,4,3) * yl(:,4,2)
c
          g3yi(:,3) = g3yi(:,3) + shg(:,1,3) * yl(:,1,3)
     &                          + shg(:,2,3) * yl(:,2,3)
     &                          + shg(:,3,3) * yl(:,3,3)
     &                          + shg(:,4,3) * yl(:,4,3)
c
          g3yi(:,4) = g3yi(:,4) + shg(:,1,3) * yl(:,1,4)
     &                          + shg(:,2,3) * yl(:,2,4)
     &                          + shg(:,3,3) * yl(:,3,4)
     &                          + shg(:,4,3) * yl(:,4,4)
c
          g3yi(:,5) = g3yi(:,5) + shg(:,1,3) * yl(:,1,5)
     &                          + shg(:,2,3) * yl(:,2,5)
     &                          + shg(:,3,3) * yl(:,3,5)
     &                          + shg(:,4,3) * yl(:,4,5)
c
        else
        do n = 1, nshl
          g1yi(:,1) = g1yi(:,1) + shg(:,n,1) * yl(:,n,1)
          g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * yl(:,n,2)
          g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * yl(:,n,3)
          g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * yl(:,n,4)
          g1yi(:,5) = g1yi(:,5) + shg(:,n,1) * yl(:,n,5)
c
          g2yi(:,1) = g2yi(:,1) + shg(:,n,2) * yl(:,n,1)
          g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * yl(:,n,2)
          g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * yl(:,n,3)
          g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * yl(:,n,4)
          g2yi(:,5) = g2yi(:,5) + shg(:,n,2) * yl(:,n,5)
c
          g3yi(:,1) = g3yi(:,1) + shg(:,n,3) * yl(:,n,1)
          g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * yl(:,n,2)
          g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * yl(:,n,3)
          g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * yl(:,n,4)
          g3yi(:,5) = g3yi(:,5) + shg(:,n,3) * yl(:,n,5)
c
        enddo
      endif


c
c     
         divqi = zero
         idflow = 0
      if (((idiff >= 1) .or. isurf==1).and.
     &     (ires == 3 .or. ires==1)) then  
         ttim(32) = ttim(32) - tmr()
c     
c     CLASS please ignore all terms related to qi until after you
c     understand EVERYTHING ELSE IN THE CODE.  You may run with
c     idiff=0 for the whole class and everything will be ok never
c     having understood this part.  (leave it for later).
c     
c.... compute divergence of diffusive flux vector, qi,i
c     
         if(idiff >= 1) then
            idflow = idflow + 4
            do n=1, nshl

               divqi(:,1) = divqi(:,1) + shg(:,n,1)*ql(:,n,1 ) 
     &              + shg(:,n,2)*ql(:,n,5 )
     &              + shg(:,n,3)*ql(:,n,9 )

               divqi(:,2) = divqi(:,2) + shg(:,n,1)*ql(:,n,2 ) 
     &              + shg(:,n,2)*ql(:,n,6 )
     &              + shg(:,n,3)*ql(:,n,10)

               divqi(:,3) = divqi(:,3) + shg(:,n,1)*ql(:,n,3 ) 
     &              + shg(:,n,2)*ql(:,n,7 )
     &              + shg(:,n,3)*ql(:,n,11)

               divqi(:,4) = divqi(:,4) + shg(:,n,1)*ql(:,n,4 ) 
     &              + shg(:,n,2)*ql(:,n,8 )
     &              + shg(:,n,3)*ql(:,n,12)

            enddo
         endif                  !end of idiff
c     
         if (isurf .eq. 1) then   
c .... divergence of normal calculation
          do n=1, nshl
             divqi(:,idflow+1) = divqi(:,idflow+1) 
     &            + shg(:,n,1)*ql(:,n,idflx-2)
     &            + shg(:,n,2)*ql(:,n,idflx-1)
     &            + shg(:,n,3)*ql(:,n,idflx)
          enddo 
c .... initialization of some variables
          Sclr = zero
          gradh= zero
          gyti = zero
          sforce=zero
          do i = 1, npro
             do n = 1, nshl      
                Sclr(i) = Sclr(i) + shape(i,n) * ycl(i,n,6) !scalar
c     
c .... compute the global gradient of Scalar variable
c     
                gyti(i,1) = gyti(i,1) + shg(i,n,1) * ycl(i,n,6) 
                gyti(i,2) = gyti(i,2) + shg(i,n,2) * ycl(i,n,6)
                gyti(i,3) = gyti(i,3) + shg(i,n,3) * ycl(i,n,6)
c     
             enddo

             if (abs (sclr(i)) .le. epsilon_ls) then
                gradh(i,1) = 0.5/epsilon_ls * (1 
     &               + cos(pi*Sclr(i)/epsilon_ls)) * gyti(i,1)
                gradh(i,2) = 0.5/epsilon_ls * (1 
     &               + cos(pi*Sclr(i)/epsilon_ls)) * gyti(i,2) 
                gradh(i,3) = 0.5/epsilon_ls * (1 
     &               + cos(pi*Sclr(i)/epsilon_ls)) * gyti(i,3)
             endif
          enddo                 !end of the loop over npro
c
c .... surface tension force calculation
c .. divide by density now as it gets multiplied in e3res.f, as surface
c    tension force is already in the form of force per unit volume
c
 
          weber(:)= Bo  
          sforce(:,1) = -(1.0/weber(:)) * divqi(:,idflow+1) !x-direction
     &         *gradh(:,1)/rho(:)
          sforce(:,2) = -(1.0/weber(:)) * divqi(:,idflow+1) !y-direction
     &         *gradh(:,2)/rho(:)
          sforce(:,3) = -(1.0/weber(:)) * divqi(:,idflow+1) !z-direction
     &         *gradh(:,3)/rho(:)          

       endif            ! end of the surface tension force calculation

         ttim(32) = ttim(32) + secs(0.0)

      endif                     ! diffusive flux computation
c
c.... return
c
       ttim(20) = ttim(20) + secs(0.0)
c
c.... ----------------------->  dist. kin energy at int. point  <--------------
c
c       
c       if (ires .ne. 2 .and. iter.eq.1)  then  !only do at beginning of step
c
c calc exact velocity for a channel at quadrature points.
c
c       dkei=0.0
c
c first guess would be this but it is wrong since it measures the
c error outside of FEM space as well.  This would be correct if we
c wanted to measure error but for disturbance we would like to get
c zero if the solution was nodally exact (i.e no disturbance added to
c the base flow which is nodally exact).
c
c      do n = 1, nshl 
c         dkei = dkei + shape(:,n) * xl(:,n,2)  ! this is just the y coord 
c      enddo
c         dkei = (1.0-dkei*dkei)  ! u_exact with u_cl=1
c
c
c  correct way
c
c       do n = 1, nshl 
c          dkei = dkei + shape(:,n) * (1.0-xl(:,n,2)**2) !u_ex^~  (in FEM space)
c       enddo
c          dkei = (u1-dkei)**2 +u2**2  ! u'^2+v'^2
c          dkei = dkei*WdetJ  ! mult function*W*det of jacobian to
c                              get this quadrature point contribution
c          dke  = dke+sum(dkei) ! we move the sum over elements inside of the
c                              sum over quadrature to save memory (we want
c                              a scalar only)
c       endif
       return
       end
        subroutine e3ivarSclr (ycl,       acl,      acti,     
     &                         shape,    shgl,     xl,      
     &                         T,        cp,
     &                         dxidx,    Sclr,               
     &                         WdetJ,    vort,  
     &                         g1yti,    g2yti,    g3yti,
     &                         rho,      rmu,      con,
     &                         rk,       u1,       u2,
     &                         u3,       shg,      dwl,
     &                         dist2w)
c
c----------------------------------------------------------------------
c
c  This routine computes the variables at integration point.
c
c input:
c  ycl     (npro,nshl,ndof)      : primitive variables
c  actl   (npro,nshl)           : time-deriv of ytl
c  dwl    (npro,nshl)           : distances to wall
c  shape  (npro,nshl)           : element shape-functions
c  shgl   (npro,nsd,nshl)       : element local-grad-shape-functions
c  xl     (npro,nenl,nsd)       : nodal coordinates at current step
c                      
c output:
c  acti   (npro)                : time-deriv of Scalar variable
c  Sclr   (npro)                : Scalar variable at intpt
c  dist2w (npro)                : distance to nearest wall at intpt
c  WdetJ  (npro)                : weighted Jacobian at intpt
c  vort   (npro)                : vorticity at intpt
c  g1yti  (npro)                : grad-Sclr in direction 1 at intpt
c  g2yti  (npro)                : grad-Sclr in direction 2 at intpt
c  g3yti  (npro)                : grad-Sclr in direction 3 at intpt
c  rho    (npro)                : density at intpt
c  rmu    (npro)                : molecular viscosity
c  con    (npro)                : conductivity
c  rk     (npro)                : kinetic energy
c  u1     (npro)                : x1-velocity component at intpt
c  u2     (npro)                : x2-velocity component at intpt
c  u3     (npro)                : x3-velocity component at intpt
c  shg    (npro,nshl,nsd)       : element global grad-shape-functions at intpt
c
c also used:
c  pres   (npro)                : pressure at intpt
c  T      (npro)                : temperature at intpt
c  ei     (npro)                : internal energy at intpt
c  h      (npro)                : enthalpy at intpt
c  alfap  (npro)                : expansivity at intpt
c  betaT  (npro)                : isothermal compressibility at intpt
c  cp     (npro)                : specific heat at constant pressure at intpt
c  dxidx  (npro,nsd,nsd)        : inverse of deformation gradient at intpt
c  divqi  (npro,nflow-1)         : divergence of diffusive flux
c
c
c Zdenek Johan, Summer 1990. (Modified from e2ivar.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997. Primitive Variables
c----------------------------------------------------------------------
c
        include "common.h"
c
c  passed arrays
c
        dimension ycl(npro,nshl,ndof),
     &            acl(npro,nshl,ndof),       acti(npro),
     &            shape(npro,nshl),  
     &            shgl(npro,nsd,nshl),       xl(npro,nenl,nsd),
     &            dui(npro,nflow),
     &            g1yi(npro,nflow),
     &            g2yi(npro,nflow),           g3yi(npro,nflow),
     &            shg(npro,nshl,nsd),        dxidx(npro,nsd,nsd),
     &            WdetJ(npro),
     &            rho(npro),                 pres(npro),
     &            T(npro),                   ei(npro),
     &            h(npro),                   alfap(npro),
     &            betaT(npro),               cp(npro),                  
     &            rk(npro),
     &            u1(npro),                  u2(npro),
     &            u3(npro),                  divqi(npro,nflow-1),
     &            ql(npro,nshl,(nflow-1)*nsd),Sclr(npro),
     &            dwl(npro,nenl),            
     &            dist2w(npro),              
     &            vort(npro),
     &            rmu(npro),                 con(npro),
     &            g1yti(npro),
     &            g2yti(npro),               g3yti(npro)
c

        dimension tmp(npro),                  dxdxi(npro,nsd,nsd)
c
        ttim(20) = ttim(20) - tmr()
c
c.... ------------->  Primitive variables at int. point  <--------------
c
c.... compute primitive variables
c
       ttim(11) = ttim(11) - tmr()

       pres   = zero
       u1     = zero
       u2     = zero
       u3     = zero
       T      = zero
       Sclr   = zero
       dist2w = zero
c
       id = isclr + 5
       do n = 1, nshl 
c
c  y(intp)=SUM_{a=1}^nshl (N_a(intp) Ya)
c
          pres   = pres   + shape(:,n) * ycl( :,n,1)
          u1     = u1     + shape(:,n) * ycl( :,n,2)
          u2     = u2     + shape(:,n) * ycl( :,n,3)
          u3     = u3     + shape(:,n) * ycl( :,n,4)
          T      = T      + shape(:,n) * ycl( :,n,5)
          Sclr   = Sclr   + shape(:,n) * ycl(:,n,id)
          if (iRANS.lt.0) then
             dist2w = dist2w + shape(:,n) * dwl(:,n)
          endif
        enddo
c
       ttim(11) = ttim(11) + tmr()
c
c.... ----------------------->  accel. at intp. point  <----------------------
c
       
       if (ires .ne. 2)  then
          ttim(12) = ttim(12) - tmr()
c
c.... compute time-derivative of Scalar variable
c
          acti = zero
          do n = 1, nshl
             acti(:)  = acti(:)  + shape(:,n) * acl(:,n,id)
          enddo
c
          flops = flops + 10*nshl*npro
          ttim(12) = ttim(12) + tmr()
       endif                    !ires .ne. 2
c
c.... ----------------->  Thermodynamic Properties  <-------------------
c
c.... compute the kinetic energy
c
        ttim(27) = ttim(27) - tmr()
c
        rk = pt5 * ( u1**2 + u2**2 + u3**2 )
c
c.... get the thermodynamic and material properties
c
        ithm = 7
        call getthm (pres,            T,               Sclr, 
     &               rk,              rho,             tmp,
     &               tmp,             tmp,             tmp,
     &               cp,              tmp,             tmp,
     &               tmp,             tmp)
c
        if (iconvsclr.eq.2) rho=one
c
        call getDiffSclr(T,            cp,          rmu,
     &                   tmp,          tmp,         con, rho, Sclr)

        ttim(27) = ttim(27) + tmr()


c
c.... --------------------->  Element Metrics  <-----------------------
c
      call e3metric( xl,         shgl,        dxidx,  
     &               shg,        WdetJ)



c.... --------------------->  Global Gradients  <-----------------------
c
        ttim(13) = ttim(13) - tmr()
        
        g1yi = zero
        g2yi = zero
        g3yi = zero
c
c.... compute the global gradient of Y-variables
c
c
c  Y_{,x_i}=SUM_{a=1}^nshl (N_{a,x_i}(intp) Ya)
c
        do n = 1, nshl
c         g1yi(:,1) = g1yi(:,1) + shg(:,n,1) * ycl(:,n,1)
c         g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * ycl(:,n,2)
          g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * ycl(:,n,3)
          g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * ycl(:,n,4)
c         g1yi(:,5) = g1yi(:,5) + shg(:,n,1) * ycl(:,n,5)
c
c         g2yi(:,1) = g2yi(:,1) + shg(:,n,2) * ycl(:,n,1)
          g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * ycl(:,n,2)
c         g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * ycl(:,n,3)
          g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * ycl(:,n,4)
c         g2yi(:,5) = g2yi(:,5) + shg(:,n,2) * ycl(:,n,5)
c
c         g3yi(:,1) = g3yi(:,1) + shg(:,n,3) * ycl(:,n,1)
          g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * ycl(:,n,2)
          g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * ycl(:,n,3)
c         g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * ycl(:,n,4)
c         g3yi(:,5) = g3yi(:,5) + shg(:,n,3) * ycl(:,n,5)
c
       enddo
       


        g1yti = zero
        g2yti = zero
        g3yti = zero
c
c.... compute the global gradient of Scalar variable
c
c
c  Y_{,x_i}=SUM_{a=1}^nshl (N_{a,x_i}(intp) Ya)
c
        do n = 1, nshl
           g1yti(:) = g1yti(:) + shg(:,n,1) * ycl(:,n,id)
           g2yti(:) = g2yti(:) + shg(:,n,2) * ycl(:,n,id)
           g3yti(:) = g3yti(:) + shg(:,n,3) * ycl(:,n,id)
c     
        enddo
c **********************************************************
c
c.... compute vorticity at this integration point
c
       vort = sqrt(
     &              (g2yi(:,4)-g3yi(:,3))**2
     &             +(g3yi(:,2)-g1yi(:,4))**2
     &             +(g1yi(:,3)-g2yi(:,2))**2
     &            ) 
c ***********************************************************

       ttim(7) = ttim(7) + tmr()
      
c
c.... return
c
       ttim(20) = ttim(20) + tmr()
       return
       end


