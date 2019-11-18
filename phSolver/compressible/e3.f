        subroutine e3 (yl,      ycl,     acl,     shp,
     &                 shgl,    xl,      rl,      rml,    xmudmi,
     &                 BDiagl,  ql,      sgn,     rlsl,   EGmass,
     &                 rerrl,   ytargetl)
c                                                                      
c----------------------------------------------------------------------
c
c This routine is the 3D element routine for the N-S equations. 
c This routine calculates the RHS residual and if requested the
c modified residual or the LHS consistent mass matrix or the block-
c diagonal preconditioner.
c
c input: 
c  yl     (npro,nshl,nflow)     : Y variables  (DOES NOT CONTAIN SCALARS)
c  ycl    (npro,nshl,ndof)      : Y variables at current step
c  acl    (npro,nshl,ndof)      : Y acceleration (Y_{,t})
c  shp    (nshl,ngauss)       : element shape-functions  N_a
c  shgl   (nsd,nshl,ngauss)   : element local-shape-functions N_{a,xi}
c  xl     (npro,nenl,nsd)       : nodal coordinates at current step (x^e_a)
c  ql     (npro,nshl,(nflow-1)*nsd) : diffusive flux vector
c  rlsl   (npro,nshl,6)       : resolved Leonard stresses
c  sgn    (npro,nshl)         : shape function sign matrix      
c
c output:
c  rl     (npro,nshl,nflow)      : element RHS residual    (G^e_a)
c  rml    (npro,nshl,nflow)      : element modified residual  (G^e_a tilde)
c  EGmass (npro,nedof,nedof)    : element LHS tangent mass matrix (dG^e_a
c                                                                  dY_b  )
c  BDiagl (npro,nshl,nflow,nflow) : block-diagonal preconditioner
c
c
c Note: This routine will calculate the element matrices for the
c        Hulbert's generalized alpha method integrator
c
c Note: nedof = nflow*nshape is the total number of degrees of freedom
c       at each element node.
c
c Mathematics done by:  Michel Mallet, Farzin Shakib (version 1)
c                       Farzin Shakib                (version 2)
c
c
c Zdenek Johan, Summer 1990.   (Modified from e2.f)
c Zdenek Johan, Winter 1991.   (Fortran 90)
c Kenneth Jansen, Winter 1997. (Primitive Variables)
c Chris Whiting, Winter 1998.  (LHS matrix formation)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension yl(npro,nshl,nflow),     ycl(npro,nshl,ndof),
     &            acl(npro,nshl,ndof),     rmu(npro),    
     &            shp(nshl,ngauss),        rlm2mu(npro),
     &            shgl(nsd,nshl,ngauss),   con(npro),
     &            xl(npro,nenl,nsd),       rlm(npro),
     &            rl(npro,nshl,nflow),     ql(npro,nshl,idflx),
     &            rml(npro,nshl,nflow),    xmudmi(npro,ngauss),
     &            BDiagl(npro,nshl,nflow,nflow),
     &            EGmass(npro,nedof,nedof),cv(npro),
     &            ytargetl(npro,nshl,nflow)
c
        dimension dui(npro,ndof),            aci(npro,ndof)
c
        dimension g1yi(npro,nflow),          g2yi(npro,nflow),
     &            g3yi(npro,nflow),          shg(npro,nshl,nsd),
     &            divqi(npro,nflow),       tau(npro,5)
c
        dimension dxidx(npro,nsd,nsd),       WdetJ(npro)
c
        dimension rho(npro),                 pres(npro),
     &            T(npro),                   ei(npro),
     &            h(npro),                   alfap(npro),
     &            betaT(npro),               DC(npro,ngauss),
     &            cp(npro),                  rk(npro),
     &            u1(npro),                  u2(npro),
     &            u3(npro),                  A0DC(npro,4),
     &            Sclr(npro),                dVdY(npro,15),
     &            giju(npro,6),              rTLS(npro),
     &            raLS(npro),                A0inv(npro,15)
c      
        dimension A0(npro,nflow,nflow),      A1(npro,nflow,nflow),
     &            A2(npro,nflow,nflow),      A3(npro,nflow,nflow)
c
        dimension rLyi(npro,nflow),          sgn(npro,nshl)
c      
        dimension ri(npro,nflow*(nsd+1)),    rmi(npro,nflow*(nsd+1)),
     &            shape(npro,nshl),          shdrv(npro,nsd,nshl),
     &            stiff(npro,nsd*nflow,nsd*nflow), 
     &            PTau(npro,5,5),
     &            sforce(npro,3),      compK(npro,10)
c
        dimension x(npro,3),              bcool(npro)

        dimension rlsl(npro,nshl,6),      rlsli(npro,6)
        
        real*8    rerrl(npro,nshl,6)
        ttim(6) = ttim(6) - secs(0.0)
c
c.... local reconstruction of diffusive flux vector
c (note: not currently included in mfg)
        if (idiff==2 .and. (ires==3 .or. ires==1)) then
           call e3ql (ycl, shp, shgl, xl, ql, xmudmi, sgn)
        endif
c
c.... loop through the integration points
c
        do intp = 1, ngauss
c
c.... if Det. .eq. 0, do not include this point
c
        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
c
c.... create a matrix of shape functions (and derivatives) for each
c     element at this quadrature point. These arrays will contain 
c     the correct signs for the hierarchic basis
c     
        call getshp(shp,          shgl,      sgn, 
     &              shape,        shdrv)
c     
c.... initialize
c
        ri  = zero
        rmi = zero
        if (lhs .eq. 1) stiff = zero
c
c
c.... calculate the integration variables
c
        ttim(8) = ttim(8) - secs(0.0)
        call e3ivar (yl,              ycl,             acl,
     &               Sclr,            shape,           shdrv,           
     &               xl,              dui,             aci,
     &               g1yi,            g2yi,            g3yi,
     &               shg,             dxidx,           WdetJ,
     &               rho,             pres,            T,
     &               ei,              h,               alfap,
     &               betaT,           cp,              rk,
     &               u1,              u2,              u3,              
     &               ql,              divqi,           sgn,
     &               rLyi,  !passed as a work array
     &               rmu,             rlm,             rlm2mu,
     &               con,             rlsl,            rlsli,
     &               xmudmi,          sforce,          cv)
        ttim(8) = ttim(8) + secs(0.0)
        
c
c.... calculate the relevant matrices
c
        ttim(9) = ttim(9) - secs(0.0)
        call e3mtrx (rho,             pres,           T,
     &               ei,              h,              alfap,
     &               betaT,           cp,             rk,
     &               u1,              u2,             u3,
     &               A0,              A1,
     &               A2,              A3,             
     &               rLyi(:,1),       rLyi(:,2),      rLyi(:,3),  ! work arrays
     &               rLyi(:,4),       rLyi(:,5),      A0DC,
     &               A0inv,           dVdY)
        ttim(9) = ttim(9) + tmr()
c
c.... calculate the convective contribution (Galerkin)
c
        ttim(14) = ttim(14) - secs(0.0)
        call e3conv (g1yi,            g2yi,            g3yi,
     &               A1,              A2,              A3,
     &               rho,             pres,            T,
     &               ei,              rk,              u1,
     &               u2,              u3,              rLyi,
     &               ri,              rmi,             EGmass,
     &               shg,             shape,           WdetJ)
        ttim(14) = ttim(14) + secs(0.0)
c
c.... calculate the diffusion contribution
c
        ttim(15) = ttim(15) - secs(0.0)
        compK = zero
        if (Navier .eq. 1) then
        call e3visc (g1yi,            g2yi,            g3yi,  
     &               dxidx,
     &               rmu,             rlm,             rlm2mu,
     &               u1,              u2,              u3,
     &               ri,              rmi,             stiff,
     &               con, rlsli,     compK, T)
        endif
        ttim(15) = ttim(15) + secs(0.0)
c
c.... calculate the body force contribution
c
        if(isurf .ne. 1 .and. matflg(5,1).gt.0) then
        call e3source (ri,            rmi,           rlyi,
     &                 rho,           u1,            u2,
     &                 u3,            pres,          sforce,
     &                 dui,           dxidx,         ytargetl,
     &                 xl,            shape,         bcool)
        else
           bcool=zero
        endif
c
c.... calculate the least-squares contribution 
c
        ttim(16) = ttim(16) - secs(0.0)
        call e3LS   (A1,              A2,            A3,
     &               rho,             rmu,           cp,
     &               cv,              con,           T,  
     &               u1,              u2,            u3,              
     &               rLyi,            dxidx,         tau,  
     &               ri,              rmi,           rk, 
     &               dui,             aci,           A0,
     &               divqi,           shape,         shg,
     &               EGmass,          stiff,         WdetJ,
     &               giju,            rTLS,          raLS,
     &               A0inv,           dVdY,          rerrl,
     &               compK,           pres,          PTau)
        ttim(16) = ttim(16) + secs(0.0)
c        
c....  Discontinuity capturing
c
        if(iDC.ne.0) then
          call e3dc  (g1yi,          g2yi,          g3yi,
     &                A0,            raLS,          rTLS,
     &                giju,          DC,            
     &                ri,            rmi,           stiff, A0DC)
        endif
! SAM wants a threshold here so we are going to take over this little used 
! error indictor for that purpose.  To revert note you will want to uncomment the original 
! form of this error indicator in e3LS.f
        if((intp.eq.1).and.(ierrcalc.eq.1).and.(nitr.eq.iter))  then
          do i=1,npro
             Tmax=maxval(yl(i,:,5))
             Tmin=minval(yl(i,:,5))
             rerrl(i,:,6)=(Tmax-Tmin)/T(i)
          enddo
! the below was somewhat suprisingly ineffective compared to above for identifying shocks.  
! it  refined on each side of the shock but left the actual shock quite coarse whereas the above
! centered well on the shock
!          do j=1,nshl
!            rerrl(:,j,6)=rerrl(:,j,6)+DC(:,intp)  !super hack to get error indicator for shock capturing
!          enddo
        endif
c
c
c.... calculate the time derivative (mass) contribution to RHS 
c
        if (ngauss .eq. 1 .and. nshl .eq. 4) then    ! trilinear tets
           ttim(17) = ttim(17) - secs(0.0)
           call e3juel (yl, acl, Sclr, A0, WdetJ, rl, rml) 
           ttim(17) = ttim(17) + secs(0.0)
        else
           call e3massr (aci, dui, ri,  rmi, A0)
        endif
        
c     
c.... calculate the time (mass) contribution to the LHS
c
        if (lhs .eq. 1) then
           call e3massl (bcool,shape,  WdetJ,   A0,  EGmass)
        endif
c
c....  calculate the preconditioner all at once now instead of in separate
c      routines
c
       if(iprec.eq.1 .and. lhs.ne.1)then
          ttim(18) = ttim(18) - secs(0.0)
          
          if (itau.lt.10) then

             call e3bdg(shape,       shg,             WdetJ,
     &		  A1,	       A2,	        A3,	
     &		  A0,          bcool,           tau,
     &            u1,          u2,              u3,    
     &            BDiagl,
     &            rmu,         rlm2mu,          con)

          else

             call e3bdg_nd(shape,       shg,             WdetJ,
     &		  A1,	       A2,	        A3,	
     &		  A0,          bcool,           PTau,
     &            u1,          u2,              u3,    
     &            BDiagl,
     &            rmu,         rlm2mu,          con)

          endif

       ttim(18) = ttim(18) + secs(0.0)
       endif
c
c
c.... multiply flux terms by shape functions and derivatives (from weight space for RHS and
c     by both the weight space and solution space for the LHS stiffness term)
c
       ttim(19) = ttim(19) - secs(0.0)
       call e3wmlt (shape,         shg,       WdetJ,
     &              ri,            rmi,       rl,        
     &              rml,           stiff,     EGmass)
       ttim(19) = ttim(19) + secs(0.0)
c
c.... end of integration loop
c
      enddo

      ttim(6) = ttim(6) + secs(0.0)
c
c.... return
c
      return
      end
c
c
c
c
       subroutine e3Sclr (ycl,          acl,  
     &                    dwl,         elDwl,            
     &                    shp,         sgn,
     &                    shgl,        xl, 
     &                    rtl,         rmtl,
     &                    qtl,         EGmasst)
c                                                                      
c----------------------------------------------------------------------
c
c This routine is the 3D element routine for the N-S equations. 
c This routine calculates the RHS residual and if requested the
c modified residual or the LHS consistent mass matrix or the block-
c diagonal preconditioner.
c
c input:    e    a   1..5   when we think of U^e_a  and U is 5 variables
c  ycl      (npro,nshl,ndof)     : Y variables
c  actl    (npro,nshl)          : scalar variable time derivative
c  dwl     (npro,nenl)          : distances to wall
c  shp     (nen,ngauss)          : element shape-functions  N_a
c  shgl    (nsd,nen,ngauss)      : element local-shape-functions N_{a,xi}
c  xl      (npro,nenl,nsd)      : nodal coordinates at current step (x^e_a)
c  qtl     (npro,nshl)          : diffusive flux vector (don't worry)
c
c output:
c  rtl     (npro,nshl)          : element RHS residual    (G^e_a)
c  rmtl    (npro,nshl)          : element modified residual  (G^e_a tilde)
c  EGmasst (npro,nshape,nshape) : element LHS tangent mass matrix (dG^e_a
c                                                                  dY_b  )
c
c
c Note: This routine will calculate the element matrices for the
c        Hulbert's generalized alpha method integrator
c
c Note: nedof = nflow*nshl is the total number of degrees of freedom
c       at each element node.
c
c Mathematics done by:  Michel Mallet, Farzin Shakib (version 1)
c                       Farzin Shakib                (version 2)
c
c
c Zdenek Johan, Summer 1990.   (Modified from e2.f)
c Zdenek Johan, Winter 1991.   (Fortran 90)
c Kenneth Jansen, Winter 1997. (Primitive Variables)
c Chris Whiting, Winter 1998.  (LHS matrix formation)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension ycl(npro,nshl,ndof),
     &            acl(npro,nshl,ndof),
     &            dwl(npro,nenl),         
     &            shp(nshl,ngauss),         shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),
     &            rtl(npro,nshl),           qtl(npro,nshl),
     &            rmtl(npro,nshl),          Diagl(npro,nshl),
     &            EGmasst(npro,nshape,nshape),
     &            dist2w(npro),             sgn(npro,nshl),   
     &            vort(npro),               gVnrm(npro),               
     &            rmu(npro),                con(npro),
     &            T(npro),                  cp(npro),
     &            g1yti(npro),              acti(npro),
     &            g2yti(npro),              g3yti(npro),
     &            Sclr(npro),               srcp (npro)

c
        dimension shg(npro,nshl,nsd) 
c
        dimension dxidx(npro,nsd,nsd),     WdetJ(npro)
c
        dimension rho(npro),               rk(npro),
     &            u1(npro),                u2(npro),
     &            u3(npro)   
c      
        dimension A0t(npro),               A1t(npro),
     &            A2t(npro),               A3t(npro)
c
        dimension rLyti(npro),             raLSt(npro),
     &            rTLSt(npro),             giju(npro,6),
     &            DCt(npro, ngauss)
c      
        dimension rti(npro,nsd+1),         rmti(npro,nsd+1),
     &            stifft(npro,nsd,nsd),
     &            shape(npro,nshl),        shdrv(npro,nsd,nshl)
        real*8    elDwl(npro)
c
        ttim(6) = ttim(6) - tmr()
c
c.... loop through the integration points
c
        elDwl(:)=zero
        do intp = 1, ngauss
c
c.... if Det. .eq. 0, do not include this point
c
        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
c
c
c.... create a matrix of shape functions (and derivatives) for each
c     element at this quadrature point. These arrays will contain 
c     the correct signs for the hierarchic basis
c     
        call getshp(shp,          shgl,      sgn, 
     &              shape,        shdrv)
c     
c.... initialize
c
        rlyti = zero
        rti   = zero
        rmti  = zero
        srcp  = zero
        stifft = zero
c        if (lhs .eq. 1) stifft = zero
c
c
c.... calculate the integration variables
c
        ttim(8) = ttim(8) - tmr()
c
        call e3ivarSclr(ycl,              acl,        acti,
     &                  shape,           shdrv,      xl,
     &                  T,               cp,
     &                  dxidx,           Sclr,        
     &                  Wdetj,           vort,       gVnrm,
     &                  g1yti,           g2yti,      g3yti,
     &                  rho,             rmu,        con,
     &                  rk,              u1,         u2,
     &                  u3,              shg,        dwl,
     &                  dist2w)
        ttim(8) = ttim(8) + tmr()

c
c.... calculate the source term contribution
c
        if(nosource.ne.1) 
     &  call e3sourceSclr (Sclr,    rho,    rmu,
     &                     dist2w,  vort,   gVnrm,   con,
     &                     g1yti,  g2yti,   g3yti,
     &                     rti,    rLyti,   srcp,
     &                     ycl,     shape,   u1,
     &                     u2,     u3,	     xl, 
     &                     elDwl)
c
         if (ilset.eq.2 .and. isclr.eq.2) then
          rk = pt5 * ( u1**2 + u2**2 + u3**2 )
         endif
c.... calculate the relevant matrices
c
        ttim(9) = ttim(9) - tmr()
        call e3mtrxSclr (rho,
     &                   u1,            u2,         u3,
     &                   A0t,           A1t,
     &                   A2t,           A3t)
        ttim(9) = ttim(9) + tmr()
c
c.... calculate the convective contribution (Galerkin)
c
        ttim(14) = ttim(14) - tmr()
        call e3convSclr (g1yti,        g2yti,       g3yti,
     &                   A1t,          A2t,         A3t,
     &                   rho,          u1,          Sclr,
     &                   u2,           u3,          rLyti,
     &                   rti,          rmti,        EGmasst,
     &                   shg,          shape,       WdetJ)
        ttim(14) = ttim(14) + tmr()
c
c.... calculate the diffusion contribution
c
        ttim(15) = ttim(15) - tmr()
        if (Navier .eq. 1) then
// Match PHASTA         if (Navier .eq. 101) then
        call e3viscSclr (g1yti,        g2yti,         g3yti,
     &                   rmu,          Sclr,          rho,
     &                   rti,          rmti,          stifft )
        endif
         ttim(15) = ttim(15) + tmr()
c
         if (ilset.eq.2)  srcp = zero
 
c
c.... calculate the least-squares contribution 
c
        ttim(16) = ttim(16) - tmr()
        call e3LSSclr(A1t,             A2t,             A3t,
     &                rho,             rmu,             rtLSt,
     &                u1,              u2,              u3,              
     &                rLyti,           dxidx,           raLSt,
     &                rti,             rk,              giju,
     &                acti,            A0t,
     &                shape,           shg,
     &                EGmasst,         stifft,          WdetJ,
     &                srcp)
        ttim(16) = ttim(16) + tmr()
c
c******************************DC TERMS*****************************
        if (idcsclr(1) .ne. 0) then
           if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or. 
     &         (idcsclr(2).eq.2 .and. isclr.eq.2)) then   ! scalar with dc
              call e3dcSclr(g1yti, g2yti,          g3yti,
     &             A0t,            raLSt,          rTLSt,
     &             DCt,            giju,         
     &             rti,            rmti,           stifft)
           endif
        endif
c     
c******************************************************************
c.... calculate the time derivative (mass) contribution to RHS
c

           call e3massrSclr (acti, rti, A0t)
c
c.... calculate the time (mass) contribution to the LHS
c
         if (lhs .eq. 1) then
            call e3masslSclr (shape,  WdetJ,   A0t,  EGmasst,srcp)
         endif
c

c.... multiply flux terms by shape functions and derivatives (from weight space for RHS and
c     by both the weight space and solution space for the LHS stiffness term)
c
       ttim(19) = ttim(19) - tmr()
       call e3wmltSclr (shape,           shg,             WdetJ,
     &                  rti,             
     &                  rtl,             stifft,          EGmasst)
       ttim(19) = ttim(19) + tmr()
c
c.... end of the loop
c
      enddo
	qptinv=one/ngauss
	elDwl(:)=elDwl(:)*qptinv


      ttim(6) = ttim(6) + tmr()
c
c.... return
c
      return
      end

