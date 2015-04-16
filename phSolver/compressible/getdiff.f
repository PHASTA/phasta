        subroutine getDiff (T,      cp,     rho,    ycl,
     &                      rmu,    rlm,    rlm2mu, con, shp,
     &                      xmudmi, xl)

c----------------------------------------------------------------------
c
c This routine calculates the fluid material properties.
c
c input:
c  T      (npro)          : temperature
c  cp     (npro)          : specific heat at constant pressure
c **************************************************************
c  rho    (npro)          : density
c  ycl    (npro,nshl,ndof): Y variables 
c  shp    (npro,nshl)     : element shape-functions
c *************************************************************
c output:
c  rmu    (npro)        : Mu
c  rlm    (npro)        : Lambda
c  rlm2mu (npro)        : Lambda + 2 Mu
c  con    (npro)        : Conductivity
c
c Note: material type flags
c         matflg(2):
c          eq. 0, constant viscosity
c          eq. 1, generalized Sutherland viscosity
c         matflg(3):
c          eq. 0, Stokes approximation
c          eq. 1, shear proportional bulk viscosity
c
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use turbSA
      use pointer_data
      include "common.h"
c
      dimension T(npro),                   cp(npro),
     &     rho(npro),                 Sclr(npro),
     &     rmu(npro),                 rlm(npro),
     &     rlm2mu(npro),              con(npro),
     &     ycl(npro,nshl,ndof),        shp(npro,nshl),
     &     xmudmi(npro,ngauss),        xl(npro,nenl,nsd),
     &     xx(npro) 
c     
      dimension xmut(npro)      
      real*8 prop_blend(npro),test_it(npro)

      integer n, e
      integer wallmask(nshl)
      real*8  xki, xki3, fv1, evisc
c
c
c.... constant viscosity
c
      if (matflg(2,1) .eq. 0) then
c     
         if (iLSet .ne. 0)then  !two fluid properties used in this model
            Sclr = zero
            isc=abs(iRANS)+6
            do n = 1, nshl
               Sclr = Sclr + shp(:,n) * ycl(:,n,isc)
            enddo
            test_it = 0.5*(one + Sclr/epsilon_ls +
     &           (sin(pi*Sclr/epsilon_ls))/pi )
            
            prop_blend = max( min(test_it(:), one ), zero  )
            rmu = datmat(1,2,2) + (datmat(1,2,1)-datmat(1,2,2))
     &           *prop_blend
         elseif(irampViscOutlet.eq.1)then ! increase viscosity near outlet
c.............ramp rmu near outlet (for a NGC geometry)
            xx=zero
            do n=1,nenl
               xx(:)=xx(:) + shp(:,n) * xl(:,n,1)
            enddo
            fmax=10.0
            where(xx(:).le. 0.42) !healfway btwn AIP and exit
               rmu(:)=datmat(1,2,1)
            elsewhere(xx(:).ge. 0.75) !2/3 of the way to the exit
               rmu(:)=fmax*datmat(1,2,1)
            elsewhere
               rmu(:)= datmat(1,2,1)*(
     &          (55.65294821-55.65294821*fmax)*xx(:)*xx(:)*xx(:)
     &          +(-97.67092412+97.67092412*fmax)*xx(:)*xx(:)
     &          +(52.59203606-52.59203606*fmax)*xx(:)
     &          -7.982719760+8.982719760*fmax)
            endwhere
         else ! constant viscosity
            rmu = datmat(1,2,1)
         endif
c     
      else
c     
c.... generalized Sutherland viscosity
c     
         rmu = datmat(1,2,1) * (T/datmat(2,2,1))*sqrt(T/datmat(2,2,1))
     &        * ( datmat(2,2,1) + datmat(3,2,1) ) / (T + datmat(3,2,1))
c     
      endif
c     
c.... calculate the second viscosity coefficient
c     
      if (matflg(3,1) .eq. 0) then
         rlm = -pt66 * rmu
      else
         rlm = (datmat(1,3,1) - pt66) * rmu
      endif
c     
c.... calculate the remaining quantities
c
      con    = rmu * cp / pr
c
c-------------Eddy Viscosity Calculation-----------------
c
c.... dynamic model
c      
      if (iLES .gt. 0. and. iRANS.eq.0) then  ! simple LES
         xmut = xmudmi(:,intp)
      else if (iRANS .eq. 0 .and. iLES.eq.0 ) then !DNS
         xmut = zero
      else if (iRANS .lt. 0) then ! calculate RANS viscosity
c
c.... RANS
c
         do e = 1, npro
            wallmask = 0
            if(itwmod.eq.-2) then ! effective viscosity
c mark the wall nodes for this element, if there are any
               do n = 1, nshl
c
c  note that we are using ycl here so that means that these
c  terms are not perturbed for MFG difference and therefore
c  NOT in the LHS.  As they only give the evisc near the wall 
c  I doubt this is a problem.
c
                  u1=ycl(e,n,2)
                  u2=ycl(e,n,3)
                  u3=ycl(e,n,4)
                  if((u1.eq.zero).and.(u2.eq.zero).and.(u3.eq.zero))
     &                 then
                     wallmask(n)=1
                  endif
               enddo
            endif
c     
            if( any(wallmask.eq.1) ) then
c if there are wall nodes for this elt in an effective-viscosity wall
c modeled case,then eddy viscosity has been stored at the wall nodes 
c in place of the spalart-allmaras variable; the eddy viscosity for 
c the whole element is taken to be the avg of wall values
               evisc = zero
               nwnode=0
               do n = 1, nshl
                  if(wallmask(n).eq.1) then
                     evisc = evisc + ycl(e,n,6)
                     nwnode = nwnode + 1
                  endif
               enddo
               evisc = evisc/nwnode
               xmut(e)= abs(evisc)
c this is what we would use instead of the above if we were allowing
c the eddy viscosity to vary through the element based on non-wall nodes
c$$$               evisc = zero
c$$$               Turb = zero
c$$$               do n = 1, nshl
c$$$                  if(wallmask(n).eq.1) then
c$$$                     evisc = evisc + shape(e,n) * ycl(e,n,6)
c$$$                  else
c$$$                     Turb = Turb + shape(e,n) * ycl(e,n,6)
c$$$                  endif
c$$$               enddo
c$$$               xki    = abs(Turb)/rmu(e)
c$$$               xki3   = xki * xki * xki
c$$$               fv1    = xki3 / (xki3 + saCv1P3)
c$$$               rmu(e) = rmu(e) + fv1*abs(Turb)               
c$$$               rmu(e) = rmu(e) + abs(evisc)
            else
c else one of the following is the case:
c   using effective-viscosity, but no wall nodes on this elt
c   using slip-velocity
c   using no model; walls are resolved 
c in all of these cases, eddy viscosity is calculated normally
               savar = zero
               do n = 1, nshl
                  savar = savar + shp(e,n) * ycl(e,n,6)
               enddo
               xki    = abs(savar)/rmu(e)
               xki3   = xki * xki * xki
               fv1    = xki3 / (xki3 + saCv1P3)
               xmut(e) = fv1*abs(savar)
            endif
         enddo                  ! end loop over elts

         if (iLES.gt.0) then    ! this is DES so we have to blend in
                                ! xmudmi based on max edge length of
                                ! element
            call EviscDES (xl,xmut,xmudmi)
         endif
      endif                     ! check for LES or RANS
      
      rlm    = rlm - pt66*xmuT
      rmu    = rmu + xmuT  
      rlm2mu = rlm + two * rmu
      con    = con + xmuT*cp/pr
c
c.... return
c
      return
      end
c
c
c
      subroutine getDiffSclr (T,      cp,   rmu,   rlm,   
     &     rlm2mu, con, rho, Sclr)
c
c----------------------------------------------------------------------
c
c This routine calculates the fluid material properties.
c
c input:
c  T      (npro)        : temperature
c  cp     (npro)        : specific heat at constant pressure
c
c output:
c  rmu    (npro)        : Mu
c  rlm    (npro)        : Lambda
c  rlm2mu (npro)        : Lambda + 2 Mu
c  con    (npro)        : Conductivity
c
c Note: material type flags
c         matflg(2):
c          eq. 0, constant viscosity
c          eq. 1, generalized Sutherland viscosity
c         matflg(3):
c          eq. 0, Stokes approximation
c          eq. 1, shear proportional bulk viscosity
c
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension T(npro),                   cp(npro),
     &            rmu(npro),                 rlm(npro),
     &            rlm2mu(npro),              con(npro),
     &            rho(npro),                 Sclr(npro)

    

c
c
c.... constant viscosity
c
        if (matflg(2,1) .eq. 0) then
c
          rmu = datmat(1,2,1)
c
        else
c
c.... generalized Sutherland viscosity
c
          rmu = datmat(1,2,1) * (T/datmat(2,2,1))*sqrt(T/datmat(2,2,1))
     &        * ( datmat(2,2,1) + datmat(3,2,1) ) / (T + datmat(3,2,1))
c
        endif
c
*************************check****************************
c        if (iRANS(1).lt.zero) then
c           rmu = saSigmaInv*rho*((rmu/rho)+Sclr)
c        endif
c This augmentation of viscosity is performed in e3viscsclr
c The Spalart -Allmaras model will need molecular viscosity 
c  in subsequent calculations.
c.... calculate the second viscosity coefficient
c
        if (matflg(3,1) .eq. 0) then
c
          rlm = -pt66 * rmu
c
        else
c
          rlm = (datmat(1,3,1) - pt66) * rmu
c
        endif
c
c.... calculate the remaining quantities
c


        
        rlm2mu = rlm + two * rmu
        con    = rmu * cp / pr



c
c.... return
c
        return
        end
      
      subroutine EviscDES(xl,xmut,xmudmi)
     
      include "common.h"
      real*8 xmut(npro),xl(npro,nenl,nsd),xmudmi(npro,ngauss)


      do i=1,npro
         dx=maxval(xl(i,:,1))-minval(xl(i,:,1))
         dy=maxval(xl(i,:,2))-minval(xl(i,:,2))
         dz=maxval(xl(i,:,3))-minval(xl(i,:,3))
         emax=max(dx,max(dy,dz))
         if(emax.lt.eles) then  ! pure les
            xmut(i)=xmudmi(i,intp)
         else if(emax.lt.two*eles) then ! blend
            xi=(emax-eles)/(eles)
            xmut(i)=xi*xmut(i)+(one-xi)*xmudmi(1,intp)
         endif                  ! leave at RANS value as edge is twice pure les
      enddo
      
      return
      end
