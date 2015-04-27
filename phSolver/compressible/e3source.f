        subroutine e3source  (ri,        rmi,          rlyi,    
     &                        rho,       u1,           u2,
     &                        u3,        pres,         sforce,
     &                        dui,       dxidx,        ytargetl,
     &                        xl,        shpfun,       bcool)
c                                                                
c----------------------------------------------------------------------
c
c This routine calculates the contribution of the bodyforce and surface
c tension force operator to the RHS vector and LHS tangent matrix. The 
c temporary results are put in ri.
c
c  u1    (npro)               : x1-velocity component
c  u2    (npro)               : x2-velocity component
c  u3    (npro)               : x3-velocity component
c  ri     (npro,nflow*(nsd+1)) : partial residual
c  rmi    (npro,nflow*(nsd+1)) : partial modified residual
c  rLyi  (npro,nflow)          : least-squares residual vector
c  shape (npro,nshl)          : element shape functions
c  g1yti  (npro)                : grad-Sclr in direction 1 at intpt
c  g2yti  (npro)                : grad-Sclr in direction 2 at intpt
c  g3yti  (npro)                : grad-Sclr in direction 3 at intpt
c
      use turbSA
      use specialBC
      include "common.h"
c
      dimension ri(npro,nflow*(nsd+1)),     rmi(npro,nflow*(nsd+1)),
     &            u1(npro),                  u2(npro),
     &            u3(npro),                  rho(npro),
     &            pres(npro),
     &            rLyi(npro,nflow),          sforce(npro,3),
     &            shpfun(npro,nshl),        
     &            xl(npro,nenl,3),           xx(npro,3)

      real*8 ytargeti(npro,nflow), ytargetl(npro,nshl,nflow)

      real*8 src(npro,nflow), bcool(npro),
     &       dui(npro,nflow), duitarg(npro,nflow), 
     &       dxidx( npro, nsd, nsd), xfind( npro ), delta(npro), rat
c
c......contribution of body force
c
      bcool=zero
      src=zero
c


      if(matflg(5,1).eq.1) then ! usual case
         src(:,1) = zero
         src(:,2) = rho(:) * datmat(1,5,1)
         src(:,3) = rho(:) * datmat(2,5,1)
         src(:,4) = rho(:) * datmat(3,5,1)
         src(:,5) = u1*src(:,2) + u2*src(:,3) + u3*src(:,4)
      else if(matflg(5,1).eq.3) then ! user supplied white noise

         xsor = 18
c            ampl = spamp(lstep+1)
c            rat = Delt(1)/0.1
         ampl = 0.002*exp(-(0.1248222*(lstep)-2.9957323)**2)
c            if((myrank.eq.zero).and.(intp.eq.ngauss)) write(*,*) ampl
         delta(:) = 0.5*sqrt(dxidx(:,1,1)*dxidx(:,1,1) ! 1/dx
     .            +dxidx(:,2,1)*dxidx(:,2,1)
     .            +dxidx(:,3,1)*dxidx(:,3,1))
         do i=1,npro
            xfind(i) = (xsor-minval(xl(i,:,1)))
     &               *(maxval(xl(i,:,1))-xsor)
         enddo
	
         where ( xfind .ge. 0. )
            src(:,2) = rho(:) * ampl * delta
c             scaling by element size is removed not to mess up
c             refinement  studies
c              src(:,2) = rho(:) * ampl 
            src(:,5) = u1*src(:,2)
         endwhere

      else if(matflg(5,1).ge.4) then ! cool case (sponge outside of a
                                     ! revolved box defined from zinSponge to
                                     ! zoutSponge in axial extent and 0
                                     ! to radSponge in radial extent for
                                     ! all theta)

c           determine coordinates of quadrature pt
            xx=zero
            do n  = 1,nenl
               xx(:,1) = xx(:,1)  + shpfun(:,n) * xl(:,n,1)
               xx(:,2) = xx(:,2)  + shpfun(:,n) * xl(:,n,2)
               xx(:,3) = xx(:,3)  + shpfun(:,n) * xl(:,n,3)
            enddo
            ytargeti=zero
            do j=1,nflow
               do n=1,nshl
                  ytargeti(:,j) = ytargeti(:,j)
     &                          + shpfun(:,n)*ytargetl(:,n,j)
               enddo
            enddo

            
c            we=3.0*29./682.
            rsteep=3.0
            src=zero
            radsts=radSponge*radSponge
            CoefRatioI2O = grthISponge/grthOSponge
            do id=1,npro
               radsqr=xx(id,2)**2+xx(id,1)**2
               if(xx(id,3).lt. zinSponge) then  ! map this into big outflow
                                                ! sponge to keep logic
                                                ! below simple
 
                  xx(id,3)=(zinSponge-xx(id,3))*CoefRatioI2O
     &                    + zoutSponge
 !
 !    CoefRatioI2O is the ratio of the inlet quadratic coefficient to the
 !    outlet quadratic coeficient (basically how much faster sponge
 !    coefficient grows in inlet region relative to outlet region)
 !                  
               endif
               if((xx(id,3).gt.zoutSponge).or.(radsqr.gt.radsts))  then
                  rad=sqrt(radsqr)
                  radc=max(rad,radSponge)
                  zval=max(xx(id,3),zoutSponge)
                  bcool(id)=grthOSponge*((zval-zoutSponge)**2
     &                                   +(radc-radSponge)**2)
                  bcool(id)=min(bcool(id),betamax)
c     Determine the resulting density and energies
               den   = ytargeti(id,1) / (Rgas * ytargeti(id,5))
               ei    = ytargeti(id,5) * ( Rgas / gamma1 )
               rk    = pt5 * ( ytargeti(id,2)**2+ytargeti(id,3)**2
     &                                         +ytargeti(id,4)**2 )
c     Determine the resulting conservation variables
               duitarg(id,1) = den
               duitarg(id,2) = den * ytargeti(id,2)
               duitarg(id,3) = den * ytargeti(id,3)
               duitarg(id,4) = den * ytargeti(id,4)
               duitarg(id,5) = den * (ei + rk)
c     Apply the sponge
               if(spongeContinuity.eq.1)
     &           src(id,1) = -bcool(id)*(dui(id,1) - duitarg(id,1))
               if(spongeMomentum1.eq.1)
     &           src(id,2) = -bcool(id)*(dui(id,2) - duitarg(id,2))
               if(spongeMomentum2.eq.1)
     &           src(id,3) = -bcool(id)*(dui(id,3) - duitarg(id,3))
               if(spongeMomentum3.eq.1)
     &           src(id,4) = -bcool(id)*(dui(id,4) - duitarg(id,4))
               if(spongeEnergy.eq.1)
     &           src(id,5) = -bcool(id)*(dui(id,5) - duitarg(id,5))
            endif
         enddo
      else
         if(isurf .ne. 1) then
            write(*,*) 'only vector (1) and cooling (4) implemented'
            stop
         endif
      endif
      
      if (isurf .eq. 1) then    ! add the surface tension force 
         src(:,2) = src(:,2) +  rho(:)*sforce(:,1)
         src(:,3) = src(:,3) +  rho(:)*sforce(:,2)
         src(:,4) = src(:,4) +  rho(:)*sforce(:,3)
         src(:,5) = src(:,5) + (u1*sforce(:,1)+u2*sforce(:,2)
     &                       + u3*sforce(:,3))*rho(:)
      endif

c
c==========================>>  IRES = 1 or 3  <<=======================
c
      if (ivart.gt.1) then
         rLyi(:,1) = rLyi(:,1) - src(:,1)
         rLyi(:,2) = rLyi(:,2) - src(:,2)
         rLyi(:,3) = rLyi(:,3) - src(:,3)
         rLyi(:,4) = rLyi(:,4) - src(:,4)
         rLyi(:,5) = rLyi(:,5) - src(:,5)
      endif
      
      if ((ires .eq. 1) .or. (ires .eq. 3)) then ! we need ri built
         ri (:,16) = ri (:,16) -  src(:,1)
         ri (:,17) = ri (:,17) -  src(:,2)
         ri (:,18) = ri (:,18) -  src(:,3) 
         ri (:,19) = ri (:,19) -  src(:,4)
         ri (:,20) = ri (:,20) -  src(:,5)
         
      endif
      
      if ((ires.eq.2) .or. (ires.eq.3)) then ! we need rmi built
         rmi (:,16) = rmi (:,16) -  src(:,1)
         rmi (:,17) = rmi (:,17) -  src(:,2)
         rmi (:,18) = rmi (:,18) -  src(:,3) 
         rmi (:,19) = rmi (:,19) -  src(:,4)
         rmi (:,20) = rmi (:,20) -  src(:,5)
      endif
c
      return
      end
c
c
c
      subroutine e3sourceSclr(Sclr,    rho,    rmu,
     &                          dist2w,  vort,   gVnrm, con,
     &                          g1yti,   g2yti,  g3yti,
     &                          rti,     rLyti,  srcp,
     &                          ycl,      shape,  u1,
     &                          u2,      u3,      xl, elDwl)
c
c---------------------------------------------------------------------
c
c  This routine calculates the source term indicated in the Spalart-
c  Allmaras eddy viscosity model.  After term is stored in rti(:,4), 
c  for later use by e3wmltSclr, and in rLyti(:) for later use by e3lsSclr.
c
c input:
c  Sclr   (npro)              : working turbulence variable
c  rho    (npro)              : density at intpt
c  rmu    (npro)              : molecular viscosity
c  dist2w (npro)              : distance from intpt to the nearest wall
c  vort   (npro)              : magnitude of the vorticity
c  gVnrm  (npro)              : magnitude of the velocity gradient
c  con    (npro)              : conductivity
c  g1yti  (npro)              : grad-Sclr in direction 1
c  g2yti  (npro)              : grad-Sclr in direction 2
c  g3yti  (npro)              : grad-Sclr in direction 3
c   
c output:
c  rti    (npro,4)            : components of residual at intpt
c  rLyti  (npro)              : GLS stabilization
c
c---------------------------------------------------------------------
c
      use turbSA
      include "common.h"
c
      dimension Sclr   (npro),        ycl(npro,nshl,ndof),
     &          dist2w (npro),          shape(npro,nshl),          
     &          vort   (npro), gVnrm(npro), rho   (npro),
     &          rmu    (npro),             con    (npro),
     &          g1yti  (npro),             g2yti  (npro),
     &          g3yti  (npro),             u1     (npro),
     &          u2     (npro),             u3     (npro)
c
      dimension rti    (npro,4),           rLyti  (npro)
c
      dimension ft1    (npro),       
c unfix later -- pieces used in acusim:
     &          srcrat (npro),             vdgn   (npro),
c    &          term1  (npro),             term2  (npro),
c    &          term3  (npro),
c
     &          chi    (npro),             fv1    (npro),
     &          fv2    (npro),             Stilde (npro),
     &          r      (npro),             g      (npro),
     &          fw     (npro),             ft2    (npro),
     &          fv1p   (npro),             fv2p   (npro),
     &          stp    (npro),             rp     (npro),
     &          gp     (npro),             fwp    (npro),
     &          bf     (npro),             srcp   (npro),
     &          gp6    (npro),             tmp    (npro),
     &          tmp1   (npro),             fwog   (npro)  
      real*8 elDwl(npro) ! local quadrature point DES dvar
      real*8 sclrm(npro) ! modified for non-negativity
      real*8 saCb1Scale(npro)  !Hack to change the production term and BL thickness
      real*8 xl_xbar(npro)     !Hack to store mean x location of element. 
c... for levelset 
      real*8  sign_levelset(npro), sclr_ls(npro), mytmp(npro),
     &        xl(npro,nenl,nsd)
    
c
      if(iRANS.lt.0) then    ! spalart almaras model
        sclrm=max(rmu/100.0,Sclr)
        if(iles.lt.0) then
          do i=1,npro
            dx=maxval(xl(i,:,1))-minval(xl(i,:,1))
            dy=maxval(xl(i,:,2))-minval(xl(i,:,2))
            dz=maxval(xl(i,:,3))-minval(xl(i,:,3))
            dmax=max(dx,max(dy,dz))
            dmax=0.65d0*dmax
            if( iles.eq.-1) then !original DES97
               dist2w(i)=min(dmax,dist2w(i))
            elseif(iles.eq.-2) then ! DDES
               rd=sclrm(i)*saKappaP2Inv/(dist2w(i)**2*gVnrm(i)+1.0d-12)
               fd=one-tanh((8.0000000000000000d0*rd)**3)
               dist2w(i)=dist2w(i)-fd*max(zero,dist2w(i)-dmax)
            endif
          enddo
        endif

        elDwl(:)=elDwl(:)+dist2w(:)
c
c  determine chi
        chi = rho*sclrm/rmu
c  determine f_v1
        fv1 = chi**3/(chi**3+saCv1**3)
c  determine f_v2
        fv2 = one - chi/(one+chi*fv1)
c  determine Stilde
        Stilde = vort + sclrm*fv2/(saKappa*dist2w)**2!unfix
c  determine r
        where(Stilde(:).ne.zero)
           r(:) = sclrm(:)/Stilde(:)/(saKappa*dist2w(:))**2
        elsewhere
           r(:) = 1.0d32
        endwhere
c  determine g
        saCw3l=saCw3
        g = r + saCw2*(r**6-r)
        sixth = 1.0/6.0
c            gp      = rp * (tmp + 5 * saCw2 * rP5)
c
c            gP6     = (g * g * g) ** 2
c            tmp     = 1 / (gP6 + (saCw3*saCw3*saCw3)**2)
c            tmp1    = ( (1 + (saCw3*saCw3*saCw3)**2) * tmp ) ** sixth
c            fw      = g * tmp1
c            fwp     = gp * tmp1 * (saCw3*saCw3*saCw3)**2 * tmp
c  determine f_w and f_w/g
        fwog = ((one+saCw3**6)/(g**6+saCw3**6))**sixth
        fw   = g*fwog
c  determine f_t2
c        ft2 = ct3*exp(-ct4*chi**2)
        ft2 = zero

c        srcrat=saCb1*(one-ft2)*Stilde*sclrm
c     &      -(saCw1*fw - saCb1*ft2/saKappa**2)*(sclrm/dist2w)**2
c        srcrat=srcrat/sclrm

!----------------------------------------------------------------------------
!HACK: lower the EV production rate within a region to decrease BL thickness. 
! Appear NM was not finished yet        if(scrScaleEnable) then  
        if(one.eq.zero) then  
          do i = 1,nenl !average the x-locations
            xl_xbar(:) = xl_xbar(:) + xl(:,i,1)
          enddo
          xl_xbar = xl_xbar/nenl
          
          saCb1Scale = one
          where(xl_xbar < saCb1alterXmin .and. xl_xbar > saCb1alterXmax)
            saCb1Scale(:) = seCb1alter
          endwhere
          
          srcrat = saCb1Scale*saCb1*(one-ft2)*Stilde
     &         -(saCw1*fw - saCb1*ft2/saKappa**2)*(sclrm/dist2w/dist2w)
        else
          srcrat=saCb1*(one-ft2)*Stilde
     &         -(saCw1*fw - saCb1*ft2/saKappa**2)*(sclrm/dist2w/dist2w)
        endif 

!Original:
!        srcrat=saCb1*(one-ft2)*Stilde
!     &       -(saCw1*fw - saCb1*ft2/saKappa**2)*(sclrm/dist2w/dist2w)
!End Hack
!----------------------------------------------------------------------------

c
c        term1=saCb1*(one-ft2)*Stilde*sclrm
c        term2=saCb2*saSigmaInv*(g1yti**2+g2yti**2+g3yti**2)
c        term3=-(saCw1*fw - saCb1*ft2/saKappa**2)*(sclrm/dist2w)**2
c determine d()/d(sclrm)
        fv1p = 3*(saCv1**3)*(chi**2)*rho
          fv1p = fv1p/(rmu*(chi**3+saCv1**3)**2)
        fv2p = (chi**2)*fv1p-(one/rmu)
          fv2p = fv2p/(one+chi*fv1)**2
        stp = fv2 + sclrm*fv2p
          stp = stp/(saKappa*dist2w)**2
        where(Stilde(:).ne.zero)
             rp(:) = Stilde(:) - sclrm(:)*stp(:)
             rp(:) = rp(:)/(saKappa*dist2w(:)*Stilde(:))**2
        elsewhere
             rp(:) = 1.0d32
        endwhere
        gp = one+saCw2*(6*r**5 - one)
          gp = gp*rp
        fwp = (saCw3**6)*fwog
          fwp = fwp*gp/(g**6+saCw3**6)
c  determine source term
        bf = saCb2*saSigmaInv*(g1yti**2+g2yti**2+g3yti**2)
     &      +saCb1*(one-ft2)*Stilde*sclrm
     &      -(saCw1*fw - saCb1*ft2/saKappa**2)*(sclrm/dist2w)**2
        bf = bf * rho
c determine d(source)/d(sclrm)
        srcp = rho*saCb1*(sclrm*stp+Stilde)
     &        -rho*saCw1*(fwp*sclrm**2 + 2*sclrm*fw)/dist2w**2
        do i=1, npro
          if(srcp(i).le.zero .and. srcp(i).le.srcrat(i)) then
            srcp(i)=srcp(i)
          else if(srcrat(i).lt.zero) then
            srcp(i)=srcrat(i)
          else 
            srcp(i)=zero
          endif
        enddo
c
c==========================>>  IRES = 1 or 3  <<=======================
c
c          if ((ires .eq. 1) .or. (ires .eq. 3)) then
             rti (:,4) = rti (:,4) -  bf(:) 
c          endif                 !ires

c          rmti (:,4) = rmti (:,4) - bf(:)
          rLyti(:) = rLyti(:) - bf(:)
c
       elseif (iLSet.ne.0) then
          if (isclr.eq.1)  then
             srcp = zero

          elseif (isclr.eq.2) then !we are redistancing level-sets

             sclr_ls = zero     !zero out temp variable

             do ii=1,npro

                do jj = 1, nshl ! first find the value of levelset at point ii
                   
                   sclr_ls(ii) =  sclr_ls(ii) 
     &                  + shape(ii,jj) * ycl(ii,jj,6)

                enddo

                if (sclr_ls(ii) .lt. - epsilon_ls)then
                   
                   sign_levelset(ii) = - one
                   
                elseif  (abs(sclr_ls(ii)) .le. epsilon_ls)then
c     sign_levelset(ii) = zero
c     
                   sign_levelset(ii) =sclr_ls(ii)/epsilon_ls 
     &                  + sin(pi*sclr_ls(ii)/epsilon_ls)/pi
                   
                   
                elseif (sclr_ls(ii) .gt. epsilon_ls) then
                   
                   sign_levelset(ii) = one
                   
                endif               
                srcp(ii) = sign_levelset(ii)
                
             enddo  
c     
c     ad   The redistancing equation can be written in the following form
c     ad
c     ad   d_{,t} + sign(phi)*( d_{,i}/|d_{,i}| )* d_{,i} = sign(phi)
c     ad
c     ad   This is rewritten in the form
c     ad
c     ad   d_{,t} + u * d_{,i} = sign(phi)
c     ad

c$$$  CAD   For the redistancing equation the "pseudo velocity" term is
c$$$  CAD   calculated as follows



             mytmp = srcp(:) / sqrt( g1yti(:) * g1yti(:) 
     &            + g2yti(:) * g2yti(:)
     &            + g3yti(:) * g3yti(:) ) 

             u1 = mytmp(:) * g1yti(:) 
             u2 = mytmp(:) * g2yti(:)
             u3 = mytmp(:) * g3yti(:)
c     
c==========================>>  IRES = 1 or 3  <<=======================
c     
c     if ((ires .eq. 1) .or. (ires .eq. 3)) then
             rti (:,4) = rti (:,4) - srcp(:) 
c     endif                 !ires

c     rmti (:,4) = rmti (:,4) -  srcp(:)
             rLyti(:) = rLyti(:) - srcp(:)
c     
          endif                 ! close of scalar 2 of level set

       else    ! NOT turbulence and NOT level set so this is a simple
               ! scalar. If your scalar equation has a source term
               ! then add your own like the above but leave an unforced case
               ! as an option like you see here

          srcp = zero
       endif


c
c.... Return and end
c
        return
        end

