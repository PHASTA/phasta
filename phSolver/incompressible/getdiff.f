      subroutine getDiff( dwl,yl, shape, xmudmi, xl, rmu,  rho)
c-----------------------------------------------------------------------
c  compute and add the contribution of the turbulent
c  eddy viscosity to the molecular viscosity.
c-----------------------------------------------------------------------
      use     turbSA
      include "common.h"

      real*8  yl(npro,nshl,ndof), rmu(npro), xmudmi(npro,ngauss),
     &        shape(npro,nshl),   rho(npro),
     &        dwl(npro,nshl),     sclr(npro),
     &        xl(npro,nenl,nsd)
      integer n, e

      real*8  epsilon_ls, kay, epsilon
     &        h_param, prop_blend(npro),test_it(npro)
c
c    
c.... get the material properties (2 fluid models will need to determine
c     the "interpolated in phase space" properties....constant for now.
c     two options exist in the interpolation 1) smooth (recommended) 
c     interpolation of nodal data, 2) discontinuous "sampling" phase at 
c     quadrature points.
c
CAD
CAD    prop_blend is a smoothing function to avoid possible large density 
CAD   gradients, e.g., water and air flows where density ratios can approach
CAD   1000.
CAD
CAD    epsilon_ls is an adjustment to control the width of the band over which
CAD    the properties are blended. 



      if (iLSet .eq. 0)then

         rho  = datmat(1,1,1)	! single fluid model, i.e., only 1 density
         rmu = datmat(1,2,1)

      else     !  two fluid properties used in this model

!        Smooth the tranistion of properties for a "distance" of epsilon_ls
!        around the interface.  Here "distance" is define as the value of the 
!        levelset function.  If the levelset function is properly defined, 
!        this is the true distance normal from the front.  Of course, the 
!        distance is in a driection normal to the front.

         Sclr = zero
         isc=abs(iRANS)+6
         do n = 1, nshl
            Sclr = Sclr + shape(:,n) * yl(:,n,isc)
         enddo
         do i= 1, npro
            if (sclr(i) .lt. - epsilon_ls)then
               prop_blend(i) = zero
            elseif  (abs(sclr(i)) .le. epsilon_ls)then
               prop_blend(i) = 0.5*(one + Sclr(i)/epsilon_ls +
     &              (sin(pi*Sclr(i)/epsilon_ls))/pi )
            elseif (sclr(i) .gt. epsilon_ls) then
               prop_blend(i) = one
            endif
         enddo
c
        rho = datmat(1,1,2) + (datmat(1,1,1)-datmat(1,1,2))*prop_blend
        rmu = datmat(1,2,2) + (datmat(1,2,1)-datmat(1,2,2))*prop_blend

      endif

CAD	At this point we have a rho that is bounded by the two values for
CAD 	density 1, datmat(1,1,1), the fluid,  and density 2, datmat(1,1,2)
CAD     the gas

c
c  The above approach evaluates all intermediate quantities at the 
c  quadrature point, then combines them to form the needed quantities there.
c  1 alternative is calculating all quanties (only rho here) at the nodes and 
c  then interpolating the result to the quadrature points.  If this is done,
c  do not forget to do the same thing for rou in e3b!!!
c  ^^^^^^^^^^
c  ||||||||||
c  WARNING
c
c.... dynamic model
c      
      if (iLES .gt. 0 .and. iRANS.eq.0) then   ! simple LES
         rmu = rmu + xmudmi(:,intp)
      else if (iRANS.lt.0) then 
         if (iRANS .eq. -1) then ! RANS (Spalart-Allmaras)
            call AddEddyViscSA(yl, shape, rmu)
         else if(iRANS.eq.-2) then ! RANS-KE
            sigmaInv=1.0        ! use full eddy viscosity for flow equations
            call AddEddyViscKE(yl, dwl, shape, rho, sigmaInv, rmu)
         endif
         if (iLES.gt.0) then    ! this is DES so we have to blend in
                                ! xmudmi based on max edge length of
                                ! element
            call EviscDESIC (xl,rmu,xmudmi)
         endif
      endif                     ! check for LES or RANS
c
      return
      end

      subroutine EviscDESIC(xl,xmut,xmudmi)
     
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
            xmut(i)=xi*xmut(i)+(one-xi)*(xmudmi(1,intp)+datmat(1,2,2))
         endif                  ! leave at RANS value as edge is twice pure les
      enddo
 !this was made messy by the addEddyVisc routines  Must clean up later.


      
      return
      end
       
      subroutine getdiffsclr(shape, dwl, yl, diffus)

      use turbSA
      use turbKE ! access to KE model constants
      include "common.h"
      
      real*8   diffus(npro), rho(npro)
      real*8   yl(npro,nshl,ndof), dwl(npro,nshl), shape(npro,nshl)
      integer n, e
      rho(:)  = datmat(1,1,1)	! single fluid model, i.e., only 1 density
      if(isclr.eq.0) then  ! solving the temperature equation
         diffus(:) = datmat(1,4,1)
      else if(iRANS.eq.-1) then ! solving SA model
         diffus(:) = datmat(1,2,1)
         call AddSAVar(yl, shape, diffus)
      else if(iRANS.eq.-2)then ! solving KE model
         diffus(:) = datmat(1,2,1)
         if(isclr.eq.2) then
            sigmaInv=1.0/ke_sigma ! different eddy viscosity for epsilon
         else
            sigmaInv=1.0 ! full eddy viscosity for solving kay equation
         endif
         call AddEddyViscKE(yl, dwl, shape, rho, sigmaInv, diffus)
      else                      ! solving scalar advection diffusion equations
         diffus = scdiff(isclr)
      endif
c
      return
      end

      function ev2sa(xmut,rm,cv1)
      implicit none
      real*8 err,ev2sa,rm,cv1,f,dfds,rat,efac
      real*8 pt5,kappa,B,xmut,chi3,denom,cv1_3
      integer iter
      pt5=0.5
      err=1.0d-6
      ev2sa=rm*cv1*1.2599       ! inflection point chi=cv1*cuberoot(2)
      kappa=0.4
c$$$        B=5.5
      efac=0.1108               ! exp(-kappa*B)
      do iter=1,50
         chi3=ev2sa/rm
         chi3=chi3*chi3*chi3
         cv1_3=cv1**3
         denom=chi3+cv1_3
         
         f=ev2sa*chi3/denom - xmut
         dfds=chi3*(chi3+4.0*cv1_3)/(denom**2)
         rat=-f/dfds
         ev2sa=ev2sa+rat
         if(abs(rat).le.err) goto 20
      enddo
      write(*,*)'ev2sa failed to converge'
      write(*,*) 'dfds,        rat,        ev2sa,        mu'
      write(*,*) dfds,rat,ev2sa,rm
 20   continue
      return
      end  
c     
      

      subroutine AddEddyViscSA(yl,shape,rmu)
      use turbSA
      include "common.h"
c     INPUTS
      double precision, intent(in), dimension(npro,nshl,ndof) ::
     &     yl
      double precision, intent(in), dimension(npro,nshl) ::
     &     shape
c     INPUT-OUTPUTS
      double precision, intent(inout), dimension(npro) ::
     &     rmu
c     LOCALS
      logical, dimension(nshl) ::
     &     wallnode
      integer e, n
      double precision xki, xki3, fv1, evisc
c
c     Loop over elements in this block
      do e = 1, npro
c        assume no wall nodes on this element
         wallnode(:) = .false.
         if(itwmod.eq.-2) then  ! effective viscosity
c           mark the wall nodes for this element, if there are any
            do n = 1, nshl
               u1=yl(e,n,2)
               u2=yl(e,n,3)
               u3=yl(e,n,4)
               if((u1.eq.zero).and.(u2.eq.zero).and.(u3.eq.zero))
     &              then
                  wallnode(n)=.true.
               endif
            enddo
         endif
c
         if( any(wallnode(:)) ) then
c if there are wall nodes for this elt, then we are using effective-
c viscosity near-wall modeling, and eddy viscosity has been stored
c at the wall nodes in place of the spalart-allmaras variable; the
c eddy viscosity for the whole element is taken to be the avg of the
c wall values
            evisc = zero
            nwnode=0
            do n = 1, nshl
               if(wallnode(n)) then
                  evisc = evisc + yl(e,n,6)
                  nwnode = nwnode + 1
               endif
            enddo
            evisc = evisc/nwnode
            rmu(e) = rmu(e) + abs(evisc)
c this is what we would use instead of the above if we were allowing
c the eddy viscosity to vary through the element based on non-wall nodes
c$$$               evisc = zero
c$$$               Turb = zero
c$$$               do n = 1, nshl
c$$$                  if(wallmask(n).eq.1) then
c$$$                     evisc = evisc + shape(e,n) * yl(e,n,6)
c$$$                  else
c$$$                     Turb = Turb + shape(e,n) * yl(e,n,6)
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
            Turb = zero
            do n = 1, nshl
               Turb = Turb + shape(e,n) * yl(e,n,6)
            enddo
            xki    = abs(Turb)/rmu(e)
            xki3   = xki * xki * xki
            fv1    = xki3 / (xki3 + saCv1P3)
            rmu(e) = rmu(e) + fv1*abs(Turb)
         endif
      enddo                     ! end loop over elts
      return
      end subroutine AddEddyViscSA



      subroutine AddSAVar(yl,shape,rmu)
      use turbSA
      include "common.h"
c     INPUTS
      double precision, intent(in), dimension(npro,nshl,ndof) ::
     &     yl
      double precision, intent(in), dimension(npro,nshl) ::
     &     shape
c     INPUT-OUTPUTS
      double precision, intent(inout), dimension(npro) ::
     &     rmu
c     LOCALS
      logical, dimension(nshl) ::
     &     wallnode
      integer e, n
      double precision savar, savarw
c     Loop over elements in this block
      do e = 1, npro
c        assume no wall nodes on this element
         wallnode(:) = .false.
         if(itwmod.eq.-2) then  ! effective viscosity
c           mark the wall nodes for this element, if there are any
            do n = 1, nshl
               u1=yl(e,n,2)
               u2=yl(e,n,3)
               u3=yl(e,n,4)
               if((u1.eq.zero).and.(u2.eq.zero).and.(u3.eq.zero))
     &              then
                  wallnode(n)=.true.
               endif
            enddo
         endif
c
         savar=zero
         do n = 1, nshl
            if( wallnode(n) ) then
c if wallmask was set, we're using effective-viscosity wall-model and
c this node is on a wall.  Eddy viscosity has been stored at the wall 
c nodes in place of the S-A variable, so we must convert it
               savarw = ev2sa(yl(e,n,6),datmat(1,2,1),saCv1)
               savar  = savar + shape(e,n) * savarw
            else
c if wallmask wasn't set, then one of the following is the case:
c   using effective-viscosity, but this isn't a wall node
c   using slip-velocity
c   using no wall model; wall is resolved
c in all these cases, the S-A variable is calculated normally
               savar  = savar + shape(e,n) * yl(e,n,6)
            endif   
         enddo
         rmu(e)=datmat(1,2,1)
         rmu(e) = (rmu(e) + abs(savar)) * saSigmaInv
      enddo                     ! end loop over elts
      return
      end subroutine AddSAVar



      subroutine AddEddyViscKE(yl, dwl, shape, rho, sigmaInv, rmu)
      use turbKE ! access to KE model constants
      include "common.h"
c     INPUTS
      double precision, intent(in), dimension(npro,nshl,ndof) ::
     &     yl
      double precision, intent(in), dimension(npro,nshl) ::
     &     shape, dwl
      double precision, intent(in), dimension(npro) ::
     &     rho
      double precision sigmaInv
c     INPUT-OUTPUTS
      double precision, intent(inout), dimension(npro) ::
     &     rmu
c     LOCALS
      double precision eviscKE, kay, epsilon, dw, CmuKE
      double precision epsInv, Rey, Ret, RetInv, tmp1, fmuKE
      integer e,n
c
      do e = 1, npro
         kay = 0.0
         epsilon = 0.0
         dw = 0.0
         do n = 1, nshl
            kay = kay + shape(e,n)*yl(e,n,6)
            epsilon = epsilon + shape(e,n)*yl(e,n,7)
            dw = dw + shape(e,n)*dwl(e,n)
         enddo
         kay = abs(kay)
         if(kay.lt.1.0e-32) kay=0.0
         epsInv	    = 0
         if ( abs(epsilon) .gt.1.e-32) then
            epsInv        = 1. / abs(epsilon)
         endif
         
         Rey                 = sqrt(kay) *    dw * rho(e) / rmu(e)       
         Ret                 = kay*kay   * epsInv * rho(e) / rmu(e)
         RetInv              = 0
         if(Ret.lt.1.d100.AND.Ret.gt.zero) RetInv=1./Ret
         tmp1     = exp(-0.0165*Rey)          
         fmuKE    = (1. -tmp1) ** 2 * (1.+20.5*RetInv) ! fmu of Lam-Bremhorst

         eviscKE=rho(e)*ke_C_mu*fmuKE*kay*kay*epsInv
         
         rmu(e) = rmu(e) + eviscKE*sigmaInv
      enddo
      return
      end subroutine AddEddyViscKE


