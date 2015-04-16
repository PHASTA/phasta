      subroutine e3bvar (yl,      acl,     ul,
     &                   shpb,    shglb,
     &                   xlb,     lnode,  
     &                   WdetJb,  bnorm,   pres,    
     &                   u1,      u2,      u3,      rmu,  
     &                   unm,     tau1n,   tau2n,   tau3n,
     &                   vdot,    rlKwall,         
     &                   xKebe,   rKwall_glob)
c
c----------------------------------------------------------------------
c
c   This routine computes the variables at integration points for 
c the boundary element routine.
c
c input:
c  yl     (npro,nshl,ndof)      : primitive variables (local)
c          ndof: 5[p,v1,v2,v3,T]+number of scalars solved 
c  acl    (npro,nshl,ndof)      : acceleration (local)
c  ul     (npro,nshlb,nsd)       : displacement (local)
c  shpb   (nen)                 : boundary element shape-functions
c  shglb  (nsd,nen)             : boundary element grad-shape-functions
c  xlb    (npro,nenl,nsd)       : nodal coordinates at current step
c  lnode  (nenb)                : local nodes on the boundary
c
c output:
c  g1yi   (npro,ndof)           : grad-v in direction 1
c  g2yi   (npro,ndof)           : grad-v in direction 2
c  g3yi   (npro,ndof)           : grad-v in direction 3
c  WdetJb (npro)                : weighted Jacobian
c  bnorm  (npro,nsd)            : outward normal
c  pres   (npro)                : pressure
c  u1     (npro)                : x1-velocity component
c  u2     (npro)                : x2-velocity component
c  u3     (npro)                : x3-velocity component
c  unm    (npro)                : BC u dot n
c  p      (npro)                : BC pressure
c  tau1n  (npro)                : BC viscous flux 1
c  tau2n  (npro)                : BC viscous flux 2
c  tau3n  (npro)                : BC viscous flux 3
c  vdot   (npro,nsd)            : acceleration at quadrature points
c  rlKwall(npro,nshlb,nsd)      : wall stiffness contribution to the local residual
c
c Zdenek Johan, Summer 1990.  (Modified from e2bvar.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c Alberto Figueroa, Winter 2004.  CMM-FSI
c----------------------------------------------------------------------
c
      use        turbsa
      include "common.h"
c
      dimension yl(npro,nshl,ndof),        rmu(npro),
     &            shpb(npro,nshl),           shglb(npro,nsd,nshl),
     &            xlb(npro,nenl,nsd),        
     &            lnode(27),                 g1yi(npro,ndof),
     &            g2yi(npro,ndof),           g3yi(npro,ndof),
     &            WdetJb(npro),              bnorm(npro,nsd),
     &            pres(npro),                
     &            u1(npro),                  u2(npro),
     &            u3(npro),
     &            unm(npro),                 
     &            tau1n(npro),               tau2n(npro),
     &            tau3n(npro),
     &            acl(npro,nshl,ndof),       ul(npro,nshl,nsd),
     &            vdot(npro,nsd),            rlKwall(npro,nshlb,nsd)
c
      dimension gl1yi(npro,ndof),          gl2yi(npro,ndof),
     &            gl3yi(npro,ndof),          dxdxib(npro,nsd,nsd),
     &            dxidxb(npro,nsd,nsd),      temp(npro),
     &            temp1(npro),               temp2(npro),
     &            temp3(npro),
     &            v1(npro,nsd),              v2(npro,nsd),
     &            v3(npro,nsd),              
     &            rotnodallocal(npro,nsd,nsd),
     &            x1rot(npro,nsd),           x2rot(npro,nsd),
     &            x3rot(npro,nsd),           detJacrot(npro),
     &            B1(npro,5,3),              B2(npro,5,3),
     &            B3(npro,5,3),              Dmatrix(npro,5,5),
     &            DtimesB1(npro,5,3),        DtimesB2(npro,5,3),
     &            DtimesB3(npro,5,3),
     &            rKwall_local11(npro,nsd,nsd),
     &            rKwall_local12(npro,nsd,nsd),
     &            rKwall_local13(npro,nsd,nsd),
     &            rKwall_local21(npro,nsd,nsd),
     &            rKwall_local22(npro,nsd,nsd),
     &            rKwall_local23(npro,nsd,nsd),
     &            rKwall_local31(npro,nsd,nsd),
     &            rKwall_local32(npro,nsd,nsd),
     &            rKwall_local33(npro,nsd,nsd),
     &            rKwall_glob11(npro,nsd,nsd),
     &            rKwall_glob12(npro,nsd,nsd),
     &            rKwall_glob13(npro,nsd,nsd),
     &            rKwall_glob21(npro,nsd,nsd),
     &            rKwall_glob22(npro,nsd,nsd),
     &            rKwall_glob23(npro,nsd,nsd),
     &            rKwall_glob31(npro,nsd,nsd),
     &            rKwall_glob32(npro,nsd,nsd),
     &            rKwall_glob33(npro,nsd,nsd)
c     
      dimension   rKwall_glob(npro,9,nshl,nshl),
     &            xKebe(npro,9,nshl,nshl)
c     
      real*8      lhmFctvw, tsFctvw(npro)

      dimension   tmp1(npro)       
c     
      real*8    Turb(npro),                xki,
     &            xki3,                      fv1
c        
      integer   e, i, j
c      
      integer   aa, b


c
c.... ------------------->  integration variables  <--------------------
c
c.... compute the primitive variables at the integration point
c
      pres = zero
      u1   = zero
      u2   = zero
      u3   = zero
c     
        
      do n = 1, nshlb
         nodlcl = lnode(n)
c     
         pres = pres + shpb(:,nodlcl) * yl(:,nodlcl,1)
         u1   = u1   + shpb(:,nodlcl) * yl(:,nodlcl,2)
         u2   = u2   + shpb(:,nodlcl) * yl(:,nodlcl,3)
         u3   = u3   + shpb(:,nodlcl) * yl(:,nodlcl,4)

      enddo
c
c.... ---------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
      dxdxib = zero
c
      do n = 1, nenl
         dxdxib(:,1,1) = dxdxib(:,1,1) + xlb(:,n,1) * shglb(:,1,n)
         dxdxib(:,1,2) = dxdxib(:,1,2) + xlb(:,n,1) * shglb(:,2,n)
         dxdxib(:,1,3) = dxdxib(:,1,3) + xlb(:,n,1) * shglb(:,3,n)
         dxdxib(:,2,1) = dxdxib(:,2,1) + xlb(:,n,2) * shglb(:,1,n)
         dxdxib(:,2,2) = dxdxib(:,2,2) + xlb(:,n,2) * shglb(:,2,n)
         dxdxib(:,2,3) = dxdxib(:,2,3) + xlb(:,n,2) * shglb(:,3,n)
         dxdxib(:,3,1) = dxdxib(:,3,1) + xlb(:,n,3) * shglb(:,1,n)
         dxdxib(:,3,2) = dxdxib(:,3,2) + xlb(:,n,3) * shglb(:,2,n)
         dxdxib(:,3,3) = dxdxib(:,3,3) + xlb(:,n,3) * shglb(:,3,n)
      enddo
c
c.... compute the normal to the boundary
c
c$$$      if (lcsyst .eq. 4) then   ! wedge-quad
c$$$         temp1 =  dxdxib(:,2,1) * dxdxib(:,3,3) -
c$$$     &            dxdxib(:,2,3) * dxdxib(:,3,1)
c$$$         temp2 =  dxdxib(:,3,1) * dxdxib(:,1,3) -
c$$$     &            dxdxib(:,3,3) * dxdxib(:,1,1)
c$$$         temp3 =  dxdxib(:,1,1) * dxdxib(:,2,3) -
c$$$     &            dxdxib(:,1,3) * dxdxib(:,2,1)
c$$$      elseif( lcyst .eq. 6) then  ! pyr-tri face
c$$$         temp1 =  dxdxib(:,2,1) * dxdxib(:,3,3) -
c$$$     &            dxdxib(:,2,3) * dxdxib(:,3,1)
c$$$         temp2 =  dxdxib(:,3,1) * dxdxib(:,1,3) -
c$$$     &            dxdxib(:,3,3) * dxdxib(:,1,1)
c$$$         temp3 =  dxdxib(:,1,1) * dxdxib(:,2,3) -
c$$$     &            dxdxib(:,1,3) * dxdxib(:,2,1)
c$$$      elseif( lcyst .eq. 1) then !usual wrong way tets
c$$$         temp1 = -dxdxib(:,2,2) * dxdxib(:,3,1) +
c$$$     &            dxdxib(:,2,1) * dxdxib(:,3,2)
c$$$         temp2 = -dxdxib(:,3,2) * dxdxib(:,1,1) +
c$$$     &            dxdxib(:,3,1) * dxdxib(:,1,2)
c$$$         temp3 = -dxdxib(:,1,2) * dxdxib(:,2,1) +
c$$$     &            dxdxib(:,1,1) * dxdxib(:,2,2)
c$$$      else
c$$$         temp1 =  dxdxib(:,2,2) * dxdxib(:,3,1) -
c$$$     &            dxdxib(:,2,1) * dxdxib(:,3,2)
c$$$         temp2 =  dxdxib(:,3,2) * dxdxib(:,1,1) -
c$$$     &            dxdxib(:,3,1) * dxdxib(:,1,2)
c$$$         temp3 =  dxdxib(:,1,2) * dxdxib(:,2,1) -
c$$$     &            dxdxib(:,1,1) * dxdxib(:,2,2)
c$$$      endif
c$$$c
c$$$      temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
c$$$      bnorm(:,1) = temp1 * temp
c$$$      bnorm(:,2) = temp2 * temp
c$$$      bnorm(:,3) = temp3 * temp
c$$$c     
c$$$      WdetJb     = Qwtb(lcsyst,intp) / temp
c$$$      if(lcsyst .eq. 3) WdetJb=WdetJb*two
c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
c
      if(lcsyst.eq.1) then      ! set to curl into element all others out
         ipt2=2
         ipt3=3
      elseif(lcsyst.eq.2) then
         ipt2=4
         ipt3=2
      elseif(lcsyst.eq.3) then
         ipt2=3
         ipt3=2
      elseif(lcsyst.eq.4) then
         ipt2=2
         ipt3=4
      elseif(lcsyst.eq.5) then
         ipt2=4
         ipt3=2
      elseif(lcsyst.eq.6) then
         ipt2=2
         ipt3=5
      endif
      v1 = xlb(:,ipt2,:) - xlb(:,1,:)
      v2 = xlb(:,ipt3,:) - xlb(:,1,:)
c
c compute cross product
c
      temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
      temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
      temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
c     
c mag is area for quads, twice area for tris
c 
      temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
      bnorm(:,1) = temp1 * temp
      bnorm(:,2) = temp2 * temp
      bnorm(:,3) = temp3 * temp
c
        
      if (lcsyst .eq. 1) then
         WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
      elseif (lcsyst .eq. 2) then
         WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
      elseif (lcsyst .eq. 3) then
         WdetJb     = Qwtb(lcsyst,intp) / (two*temp)
      elseif (lcsyst .eq. 4) then
         WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
      elseif (lcsyst .eq. 5) then
         WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
      elseif (lcsyst .eq. 6) then
         WdetJb     = Qwtb(lcsyst,intp) / (two*temp)
      endif
c
c.... -------------------------->  Grad-V  <----------------------------
c
c.... compute grad-v for Navier-Stokes terms
c
      if (Navier .eq. 1) then
c
c.... compute the inverse of deformation gradient
c
         dxidxb(:,1,1) =   dxdxib(:,2,2) * dxdxib(:,3,3) 
     &        - dxdxib(:,3,2) * dxdxib(:,2,3)
         dxidxb(:,1,2) =   dxdxib(:,3,2) * dxdxib(:,1,3) 
     &        - dxdxib(:,1,2) * dxdxib(:,3,3)
         dxidxb(:,1,3) =   dxdxib(:,1,2) * dxdxib(:,2,3) 
     &        - dxdxib(:,1,3) * dxdxib(:,2,2)
         temp          = one / ( dxidxb(:,1,1) * dxdxib(:,1,1) 
     &        + dxidxb(:,1,2) * dxdxib(:,2,1)  
     &        + dxidxb(:,1,3) * dxdxib(:,3,1) )
         dxidxb(:,1,1) =  dxidxb(:,1,1) * temp
         dxidxb(:,1,2) =  dxidxb(:,1,2) * temp
         dxidxb(:,1,3) =  dxidxb(:,1,3) * temp
         dxidxb(:,2,1) = (dxdxib(:,2,3) * dxdxib(:,3,1) 
     &        - dxdxib(:,2,1) * dxdxib(:,3,3)) * temp
         dxidxb(:,2,2) = (dxdxib(:,1,1) * dxdxib(:,3,3) 
     &        - dxdxib(:,3,1) * dxdxib(:,1,3)) * temp
         dxidxb(:,2,3) = (dxdxib(:,2,1) * dxdxib(:,1,3) 
     &        - dxdxib(:,1,1) * dxdxib(:,2,3)) * temp
         dxidxb(:,3,1) = (dxdxib(:,2,1) * dxdxib(:,3,2) 
     &        - dxdxib(:,2,2) * dxdxib(:,3,1)) * temp
         dxidxb(:,3,2) = (dxdxib(:,3,1) * dxdxib(:,1,2) 
     &        - dxdxib(:,1,1) * dxdxib(:,3,2)) * temp
         dxidxb(:,3,3) = (dxdxib(:,1,1) * dxdxib(:,2,2) 
     &        - dxdxib(:,1,2) * dxdxib(:,2,1)) * temp
c
c.... compute local-grad-Y
c
         gl1yi = zero
         gl2yi = zero
         gl3yi = zero
c     
         do n = 1, nshl
            gl1yi(:,1) = gl1yi(:,1) + shglb(:,1,n) * yl(:,n,1)
            gl1yi(:,2) = gl1yi(:,2) + shglb(:,1,n) * yl(:,n,2)
            gl1yi(:,3) = gl1yi(:,3) + shglb(:,1,n) * yl(:,n,3)
            gl1yi(:,4) = gl1yi(:,4) + shglb(:,1,n) * yl(:,n,4)
c     
            gl2yi(:,1) = gl2yi(:,1) + shglb(:,2,n) * yl(:,n,1)
            gl2yi(:,2) = gl2yi(:,2) + shglb(:,2,n) * yl(:,n,2)
            gl2yi(:,3) = gl2yi(:,3) + shglb(:,2,n) * yl(:,n,3)
            gl2yi(:,4) = gl2yi(:,4) + shglb(:,2,n) * yl(:,n,4)
c     
            gl3yi(:,1) = gl3yi(:,1) + shglb(:,3,n) * yl(:,n,1)
            gl3yi(:,2) = gl3yi(:,2) + shglb(:,3,n) * yl(:,n,2)
            gl3yi(:,3) = gl3yi(:,3) + shglb(:,3,n) * yl(:,n,3)
            gl3yi(:,4) = gl3yi(:,4) + shglb(:,3,n) * yl(:,n,4)
         enddo
c     
c.... convert local-grads to global-grads
c     
         g1yi(:,2) = dxidxb(:,1,1) * gl1yi(:,2) + 
     &        dxidxb(:,2,1) * gl2yi(:,2) +
     &        dxidxb(:,3,1) * gl3yi(:,2)
         g2yi(:,2) = dxidxb(:,1,2) * gl1yi(:,2) + 
     &        dxidxb(:,2,2) * gl2yi(:,2) +
     &        dxidxb(:,3,2) * gl3yi(:,2)
         g3yi(:,2) = dxidxb(:,1,3) * gl1yi(:,2) + 
     &        dxidxb(:,2,3) * gl2yi(:,2) +
     &        dxidxb(:,3,3) * gl3yi(:,2)
c     
         g1yi(:,3) = dxidxb(:,1,1) * gl1yi(:,3) + 
     &        dxidxb(:,2,1) * gl2yi(:,3) +
     &        dxidxb(:,3,1) * gl3yi(:,3)
         g2yi(:,3) = dxidxb(:,1,2) * gl1yi(:,3) + 
     &        dxidxb(:,2,2) * gl2yi(:,3) +
     &        dxidxb(:,3,2) * gl3yi(:,3)
         g3yi(:,3) = dxidxb(:,1,3) * gl1yi(:,3) + 
     &        dxidxb(:,2,3) * gl2yi(:,3) +
     &        dxidxb(:,3,3) * gl3yi(:,3)
c     
         g1yi(:,4) = dxidxb(:,1,1) * gl1yi(:,4) + 
     &        dxidxb(:,2,1) * gl2yi(:,4) +
     &        dxidxb(:,3,1) * gl3yi(:,4)
         g2yi(:,4) = dxidxb(:,1,2) * gl1yi(:,4) + 
     &        dxidxb(:,2,2) * gl2yi(:,4) +
     &        dxidxb(:,3,2) * gl3yi(:,4)
         g3yi(:,4) = dxidxb(:,1,3) * gl1yi(:,4) + 
     &        dxidxb(:,2,3) * gl2yi(:,4) +
     &        dxidxb(:,3,3) * gl3yi(:,4)
c     
c.... end grad-v
c     
      endif

c     
c.... mass flux
c     
      unm = bnorm(:,1) * u1 +bnorm(:,2) * u2  +bnorm(:,3) * u3
! no rho in continuity eq.


c
c.... viscous flux
c
      tau1n = bnorm(:,1) * two * rmu *  g1yi(:,2)  
     &     + bnorm(:,2) *      (rmu * (g2yi(:,2) + g1yi(:,3)))
     &     + bnorm(:,3) *      (rmu * (g3yi(:,2) + g1yi(:,4)))
      tau2n = bnorm(:,1) *      (rmu * (g2yi(:,2) + g1yi(:,3)))
     &     + bnorm(:,2) * two * rmu *  g2yi(:,3) 
     &     + bnorm(:,3) *      (rmu * (g3yi(:,3) + g2yi(:,4)))
      tau3n = bnorm(:,1) *      (rmu * (g3yi(:,2) + g1yi(:,4)))
     &     + bnorm(:,2) *      (rmu * (g3yi(:,3) + g2yi(:,4)))
     &     + bnorm(:,3) * two * rmu *  g3yi(:,4) 
c     
      temp1 = bnorm(:,1) * tau1n
     &     + bnorm(:,2) * tau2n
     &     + bnorm(:,3) * tau3n
      
      pres  = pres - temp1
      
      tau1n = tau1n - bnorm(:,1) * temp1
      tau2n = tau2n - bnorm(:,2) * temp1
      tau3n = tau3n - bnorm(:,3) * temp1

c
c.... viscous flux control
c
c     if iviscflux = 1, we consider the viscous flux on the RHS
c     otherwise, if iviscflux = 0, we eliminate this term: stability in 
c     pressure-flow coupled boundaries for cardiovascular applications in 
c     situations of flow reversal
      tau1n = tau1n * iviscflux
      tau2n = tau2n * iviscflux
      tau3n = tau3n * iviscflux

      vdot = zero
      rlKwall = zero
      if (intp.eq.ngaussb)   then    ! do this only for the last gauss point
        rKwall_glob = zero
      endif

      if(ideformwall.eq.1) then
      do n = 1, nshlb
         nodlcl = lnode(n)
c     
         vdot(:,1) = vdot(:,1) + shpb(:,nodlcl) * acl(:,nodlcl,2)
         vdot(:,2) = vdot(:,2) + shpb(:,nodlcl) * acl(:,nodlcl,3)
         vdot(:,3) = vdot(:,3) + shpb(:,nodlcl) * acl(:,nodlcl,4)

      enddo
      vdot = vdot * thicknessvw * rhovw
c     
c.... --------------------->  Stiffness matrix & residual  <-----------------
c     
c.... B^t * D * B formulation for plane stress enhanced membrane
c
c
c.... rotation matrix
c     
      v1 = xlb(:,ipt2,:) - xlb(:,1,:)
      temp       = one / sqrt ( v1(:,1)**2 + v1(:,2)**2 + v1(:,3)**2 )
      v1(:,1) = v1(:,1) * temp
      v1(:,2) = v1(:,2) * temp
      v1(:,3) = v1(:,3) * temp
      
      v2 = xlb(:,ipt3,:) - xlb(:,1,:)
      
c     compute cross product
      temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
      temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
      temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
      
      temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
      v3(:,1) = temp1 * temp
      v3(:,2) = temp2 * temp
      v3(:,3) = temp3 * temp
      
c     cross product again for v2
      temp1 = v3(:,2) * v1(:,3) - v1(:,2) * v3(:,3)
      temp2 = v1(:,1) * v3(:,3) - v3(:,1) * v1(:,3)
      temp3 = v3(:,1) * v1(:,2) - v1(:,1) * v3(:,2)
      
      temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
      v2(:,1) = temp1 * temp
      v2(:,2) = temp2 * temp
      v2(:,3) = temp3 * temp
      
      do j = 1, nsd
         rotnodallocal(:,1,j) = v1(:,j)
         rotnodallocal(:,2,j) = v2(:,j)
         rotnodallocal(:,3,j) = v3(:,j)
      enddo
      
c     
c.... rotated coordinates
c
      x1rot = zero
      x2rot = zero
      x3rot = zero
      
      do i = 1, nsd
         do j = 1, nsd 
            x1rot(:,i) = x1rot(:,i)+rotnodallocal(:,i,j)*xlb(:,1,j)
            x2rot(:,i) = x2rot(:,i)+rotnodallocal(:,i,j)*xlb(:,ipt2,j)
            x3rot(:,i) = x3rot(:,i)+rotnodallocal(:,i,j)*xlb(:,ipt3,j)
         enddo
      enddo
      
c     
c.... B matrices
c     
      B1 = zero
      B2 = zero
      B3 = zero
      detJacrot = (x2rot(:,1)-x1rot(:,1)) * (x3rot(:,2)-x1rot(:,2)) - 
     &     (x3rot(:,1)-x1rot(:,1)) * (x2rot(:,2)-x1rot(:,2))
      
      B1(:,1,1) = (x2rot(:,2)-x3rot(:,2))/detJacrot(:)
      B1(:,2,2) = (x3rot(:,1)-x2rot(:,1))/detJacrot(:)
      B1(:,3,1) = (x3rot(:,1)-x2rot(:,1))/detJacrot(:)
      B1(:,3,2) = (x2rot(:,2)-x3rot(:,2))/detJacrot(:)
      B1(:,4,3) = (x2rot(:,2)-x3rot(:,2))/detJacrot(:)
      B1(:,5,3) = (x3rot(:,1)-x2rot(:,1))/detJacrot(:)
      
      B2(:,1,1) = (x3rot(:,2)-x1rot(:,2))/detJacrot(:)
      B2(:,2,2) = (x1rot(:,1)-x3rot(:,1))/detJacrot(:)
      B2(:,3,1) = (x1rot(:,1)-x3rot(:,1))/detJacrot(:)
      B2(:,3,2) = (x3rot(:,2)-x1rot(:,2))/detJacrot(:)
      B2(:,4,3) = (x3rot(:,2)-x1rot(:,2))/detJacrot(:)
      B2(:,5,3) = (x1rot(:,1)-x3rot(:,1))/detJacrot(:)
      
      B3(:,1,1) = (x1rot(:,2)-x2rot(:,2))/detJacrot(:)
      B3(:,2,2) = (x2rot(:,1)-x1rot(:,1))/detJacrot(:)
      B3(:,3,1) = (x2rot(:,1)-x1rot(:,1))/detJacrot(:)
      B3(:,3,2) = (x1rot(:,2)-x2rot(:,2))/detJacrot(:)
      B3(:,4,3) = (x1rot(:,2)-x2rot(:,2))/detJacrot(:)
      B3(:,5,3) = (x2rot(:,1)-x1rot(:,1))/detJacrot(:)
      
C      B1 = B1 / detJacrot
C      B2 = B2 / detJacrot
C      B3 = B3 / detJacrot
      
c     
c.... D matrix
c     
      Dmatrix = zero
      temp1 = evw / (1.0d0 - rnuvw*rnuvw)
      temp2 = rnuvw * temp1
      temp3 = pt5 * (1.0d0 - rnuvw) * temp1
      Dmatrix(:,1,1) = temp1
      Dmatrix(:,1,2) = temp2
      Dmatrix(:,2,1) = temp2
      Dmatrix(:,2,2) = temp1
      Dmatrix(:,3,3) = temp3
      Dmatrix(:,4,4) = temp3*rshearconstantvw
      Dmatrix(:,5,5) = temp3*rshearconstantvw
c     
c.... D * [B1|B2|B3]
c     
      DtimesB1 = zero
      DtimesB2 = zero
      DtimesB3 = zero
      do i = 1, 5
         do j = 1, 3
            do k = 1, 5
               DtimesB1(:,i,j) = DtimesB1(:,i,j) 
     &              + Dmatrix(:,i,k) * B1(:,k,j)
               DtimesB2(:,i,j) = DtimesB2(:,i,j) 
     &              + Dmatrix(:,i,k) * B2(:,k,j)
               DtimesB3(:,i,j) = DtimesB3(:,i,j) 
     &              + Dmatrix(:,i,k) * B3(:,k,j)
            enddo
         enddo
      enddo
c     
c.... [B1|B2|B3]^T * D * [B1|B2|B3]
c     
      rKwall_local11 = zero
      rKwall_local12 = zero
      rKwall_local13 = zero
      rKwall_local21 = zero
      rKwall_local22 = zero
      rKwall_local23 = zero
      rKwall_local31 = zero
      rKwall_local32 = zero
      rKwall_local33 = zero
      
      do i = 1, 3               ! i is a node index: i=1, nenbl=3
         do j = 1, 3            ! same is true for j
            do k = 1, 5
               rKwall_local11(:,i,j)  = rKwall_local11(:,i,j)
     &              + B1(:,k,i) * DtimesB1(:,k,j)
               rKwall_local12(:,i,j)  = rKwall_local12(:,i,j)
     &              + B1(:,k,i) * DtimesB2(:,k,j)
               rKwall_local13(:,i,j)  = rKwall_local13(:,i,j)
     &              + B1(:,k,i) * DtimesB3(:,k,j)
               rKwall_local21(:,i,j)  = rKwall_local21(:,i,j)
     &              + B2(:,k,i) * DtimesB1(:,k,j)
               rKwall_local22(:,i,j)  = rKwall_local22(:,i,j)
     &              + B2(:,k,i) * DtimesB2(:,k,j)
               rKwall_local23(:,i,j)  = rKwall_local23(:,i,j)
     &              + B2(:,k,i) * DtimesB3(:,k,j)
               rKwall_local31(:,i,j)  = rKwall_local31(:,i,j)
     &              + B3(:,k,i) * DtimesB1(:,k,j)
               rKwall_local32(:,i,j)  = rKwall_local32(:,i,j)
     &              + B3(:,k,i) * DtimesB2(:,k,j)
               rKwall_local33(:,i,j)  = rKwall_local33(:,i,j)
     &              + B3(:,k,i) * DtimesB3(:,k,j)
            enddo
         enddo
      enddo
      
c     
c.... Now we need to rotate each of these submatrices to the global frame
c     
      call rotatestiff(rKwall_local11, rotnodallocal, rKwall_glob11)
      call rotatestiff(rKwall_local12, rotnodallocal, rKwall_glob12)
      call rotatestiff(rKwall_local13, rotnodallocal, rKwall_glob13)
      call rotatestiff(rKwall_local21, rotnodallocal, rKwall_glob21)
      call rotatestiff(rKwall_local22, rotnodallocal, rKwall_glob22)
      call rotatestiff(rKwall_local23, rotnodallocal, rKwall_glob23)
      call rotatestiff(rKwall_local31, rotnodallocal, rKwall_glob31)
      call rotatestiff(rKwall_local32, rotnodallocal, rKwall_glob32)
      call rotatestiff(rKwall_local33, rotnodallocal, rKwall_glob33)

c     multiply the nodal matrices by the area and the thickness
      do i =1, nsd
         do j = 1, nsd
            rKwall_glob11(:,i,j) = rKwall_glob11(:,i,j) * detJacrot(:) 
     &                           * pt5 * thicknessvw
            rKwall_glob12(:,i,j) = rKwall_glob12(:,i,j) * detJacrot(:) 
     &                           * pt5 * thicknessvw
            rKwall_glob13(:,i,j) = rKwall_glob13(:,i,j) * detJacrot(:) 
     &                           * pt5 * thicknessvw
            rKwall_glob21(:,i,j) = rKwall_glob21(:,i,j) * detJacrot(:) 
     &                           * pt5 * thicknessvw
            rKwall_glob22(:,i,j) = rKwall_glob22(:,i,j) * detJacrot(:) 
     &                           * pt5 * thicknessvw
            rKwall_glob23(:,i,j) = rKwall_glob23(:,i,j) * detJacrot(:) 
     &                           * pt5 * thicknessvw
            rKwall_glob31(:,i,j) = rKwall_glob31(:,i,j) * detJacrot(:) 
     &                           * pt5 * thicknessvw
            rKwall_glob32(:,i,j) = rKwall_glob32(:,i,j) * detJacrot(:) 
     &                           * pt5 * thicknessvw
            rKwall_glob33(:,i,j) = rKwall_glob33(:,i,j) * detJacrot(:) 
     &                           * pt5 * thicknessvw
         enddo
      enddo

c     
c.... Final K * u product (in global coordinates) to get the residual
c
      do i = 1, 3               ! now i is a spatial index: i=1, nsd=3
         rlKwall(:,1,1) = rlKwall(:,1,1) 
     &                  + rKwall_glob11(:,1,i) * ul(:,1,i) 
     &                  + rKwall_glob12(:,1,i) * ul(:,2,i) 
     &                  + rKwall_glob13(:,1,i) * ul(:,3,i) 
         rlKwall(:,1,2) = rlKwall(:,1,2)
     &                  + rKwall_glob11(:,2,i) * ul(:,1,i) 
     &                  + rKwall_glob12(:,2,i) * ul(:,2,i) 
     &                  + rKwall_glob13(:,2,i) * ul(:,3,i) 
         rlKwall(:,1,3) = rlKwall(:,1,3) 
     &                  + rKwall_glob11(:,3,i) * ul(:,1,i) 
     &                  + rKwall_glob12(:,3,i) * ul(:,2,i) 
     &                  + rKwall_glob13(:,3,i) * ul(:,3,i) 
         rlKwall(:,2,1) = rlKwall(:,2,1) 
     &                  + rKwall_glob21(:,1,i) * ul(:,1,i) 
     &                  + rKwall_glob22(:,1,i) * ul(:,2,i) 
     &                  + rKwall_glob23(:,1,i) * ul(:,3,i) 
         rlKwall(:,2,2) = rlKwall(:,2,2)
     &                  + rKwall_glob21(:,2,i) * ul(:,1,i) 
     &                  + rKwall_glob22(:,2,i) * ul(:,2,i) 
     &                  + rKwall_glob23(:,2,i) * ul(:,3,i) 
         rlKwall(:,2,3) = rlKwall(:,2,3) 
     &                  + rKwall_glob21(:,3,i) * ul(:,1,i) 
     &                  + rKwall_glob22(:,3,i) * ul(:,2,i) 
     &                  + rKwall_glob23(:,3,i) * ul(:,3,i)
         rlKwall(:,3,1) = rlKwall(:,3,1) 
     &                  + rKwall_glob31(:,1,i) * ul(:,1,i) 
     &                  + rKwall_glob32(:,1,i) * ul(:,2,i) 
     &                  + rKwall_glob33(:,1,i) * ul(:,3,i) 
         rlKwall(:,3,2) = rlKwall(:,3,2)
     &                  + rKwall_glob31(:,2,i) * ul(:,1,i) 
     &                  + rKwall_glob32(:,2,i) * ul(:,2,i) 
     &                  + rKwall_glob33(:,2,i) * ul(:,3,i) 
         rlKwall(:,3,3) = rlKwall(:,3,3) 
     &                  + rKwall_glob31(:,3,i) * ul(:,1,i) 
     &                  + rKwall_glob32(:,3,i) * ul(:,2,i) 
     &                  + rKwall_glob33(:,3,i) * ul(:,3,i)
      enddo
c     
c.... --------------> End of Stiffness matrix & residual  <-----------------
c     

c     
c.... -----> Wall Stiffness and Mass matrices for implicit LHS  <-----------
c     

c....  Here we just add the mass matrix contribution.  The stiffness contribution 
c....  is added in e3b

c      lhmFct = almi * (one - flmpl)      Maybe we have to define flmplW: lumped
                                        ! mass parameter for the wall
      lhmFctvw = almi * (one - flmpl)                                  
c
c.... scale variables for efficiency
c
      tsFctvw     = lhmFctvw * WdetJb * rhovw * thicknessvw     
c
c.... compute mass and convection terms
c
c.... NOTE:  the wall mass contributions should only have 3 nodal components 
c.... since the fourth node is an interior node... therefore, the loops should
c.... be done from 1 to nshlb=3...

      do b = 1, nshlb
         do aa = 1, nshlb
            tmp1 = tsFctvw * shpb(:,aa) * shpb(:,b)
c
c           tmp1=alpha_m*(1-lmp)*WdetJ*N^aN^b*rho*thickness   the time term 
c            
            xKebe(:,1,aa,b) = xKebe(:,1,aa,b) + tmp1
            xKebe(:,5,aa,b) = xKebe(:,5,aa,b) + tmp1
            xKebe(:,9,aa,b) = xKebe(:,9,aa,b) + tmp1
         enddo
      enddo

c
c.... assemble the nodal stiffness into the element stiffness matrix rKwall_glob
c
c.... We have passed the integer intp to make this operation only once: we are 
c.... not using the gauss points structure to compute the stiffness of the wall 
c.... elements, so we don't want to be redundant and calculate ngaussb times the 
c.... stiffness matrix which is constant for linear triangles...

c.... This is ugly, but I will fix it later...

      if (intp.eq.ngaussb)   then    ! do this only for the last gauss point
        rKwall_glob(:,1,1,1) = rKwall_glob11(:,1,1)
        rKwall_glob(:,2,1,1) = rKwall_glob11(:,1,2)
        rKwall_glob(:,3,1,1) = rKwall_glob11(:,1,3)
        rKwall_glob(:,4,1,1) = rKwall_glob11(:,2,1)
        rKwall_glob(:,5,1,1) = rKwall_glob11(:,2,2)
        rKwall_glob(:,6,1,1) = rKwall_glob11(:,2,3)
        rKwall_glob(:,7,1,1) = rKwall_glob11(:,3,1)
        rKwall_glob(:,8,1,1) = rKwall_glob11(:,3,2)
        rKwall_glob(:,9,1,1) = rKwall_glob11(:,3,3)
      
        rKwall_glob(:,1,1,2) = rKwall_glob12(:,1,1)
        rKwall_glob(:,2,1,2) = rKwall_glob12(:,1,2)
        rKwall_glob(:,3,1,2) = rKwall_glob12(:,1,3)
        rKwall_glob(:,4,1,2) = rKwall_glob12(:,2,1)
        rKwall_glob(:,5,1,2) = rKwall_glob12(:,2,2)
        rKwall_glob(:,6,1,2) = rKwall_glob12(:,2,3)
        rKwall_glob(:,7,1,2) = rKwall_glob12(:,3,1)
        rKwall_glob(:,8,1,2) = rKwall_glob12(:,3,2)
        rKwall_glob(:,9,1,2) = rKwall_glob12(:,3,3)
      
        rKwall_glob(:,1,1,3) = rKwall_glob13(:,1,1)
        rKwall_glob(:,2,1,3) = rKwall_glob13(:,1,2)
        rKwall_glob(:,3,1,3) = rKwall_glob13(:,1,3)
        rKwall_glob(:,4,1,3) = rKwall_glob13(:,2,1)
        rKwall_glob(:,5,1,3) = rKwall_glob13(:,2,2)
        rKwall_glob(:,6,1,3) = rKwall_glob13(:,2,3)
        rKwall_glob(:,7,1,3) = rKwall_glob13(:,3,1)
        rKwall_glob(:,8,1,3) = rKwall_glob13(:,3,2)
        rKwall_glob(:,9,1,3) = rKwall_glob13(:,3,3)
      
        rKwall_glob(:,1,2,1) = rKwall_glob21(:,1,1)
        rKwall_glob(:,2,2,1) = rKwall_glob21(:,1,2)
        rKwall_glob(:,3,2,1) = rKwall_glob21(:,1,3)
        rKwall_glob(:,4,2,1) = rKwall_glob21(:,2,1)
        rKwall_glob(:,5,2,1) = rKwall_glob21(:,2,2)
        rKwall_glob(:,6,2,1) = rKwall_glob21(:,2,3)
        rKwall_glob(:,7,2,1) = rKwall_glob21(:,3,1)
        rKwall_glob(:,8,2,1) = rKwall_glob21(:,3,2)
        rKwall_glob(:,9,2,1) = rKwall_glob21(:,3,3)
      
        rKwall_glob(:,1,2,2) = rKwall_glob22(:,1,1)
        rKwall_glob(:,2,2,2) = rKwall_glob22(:,1,2)
        rKwall_glob(:,3,2,2) = rKwall_glob22(:,1,3)
        rKwall_glob(:,4,2,2) = rKwall_glob22(:,2,1)
        rKwall_glob(:,5,2,2) = rKwall_glob22(:,2,2)
        rKwall_glob(:,6,2,2) = rKwall_glob22(:,2,3)
        rKwall_glob(:,7,2,2) = rKwall_glob22(:,3,1)
        rKwall_glob(:,8,2,2) = rKwall_glob22(:,3,2)
        rKwall_glob(:,9,2,2) = rKwall_glob22(:,3,3)      

        rKwall_glob(:,1,2,3) = rKwall_glob23(:,1,1)
        rKwall_glob(:,2,2,3) = rKwall_glob23(:,1,2)
        rKwall_glob(:,3,2,3) = rKwall_glob23(:,1,3)
        rKwall_glob(:,4,2,3) = rKwall_glob23(:,2,1)
        rKwall_glob(:,5,2,3) = rKwall_glob23(:,2,2)
        rKwall_glob(:,6,2,3) = rKwall_glob23(:,2,3)
        rKwall_glob(:,7,2,3) = rKwall_glob23(:,3,1)
        rKwall_glob(:,8,2,3) = rKwall_glob23(:,3,2)
        rKwall_glob(:,9,2,3) = rKwall_glob23(:,3,3)
      
        rKwall_glob(:,1,3,1) = rKwall_glob31(:,1,1)
        rKwall_glob(:,2,3,1) = rKwall_glob31(:,1,2)
        rKwall_glob(:,3,3,1) = rKwall_glob31(:,1,3)
        rKwall_glob(:,4,3,1) = rKwall_glob31(:,2,1)
        rKwall_glob(:,5,3,1) = rKwall_glob31(:,2,2)
        rKwall_glob(:,6,3,1) = rKwall_glob31(:,2,3)
        rKwall_glob(:,7,3,1) = rKwall_glob31(:,3,1)
        rKwall_glob(:,8,3,1) = rKwall_glob31(:,3,2)
        rKwall_glob(:,9,3,1) = rKwall_glob31(:,3,3)

        rKwall_glob(:,1,3,2) = rKwall_glob32(:,1,1)
        rKwall_glob(:,2,3,2) = rKwall_glob32(:,1,2)
        rKwall_glob(:,3,3,2) = rKwall_glob32(:,1,3)
        rKwall_glob(:,4,3,2) = rKwall_glob32(:,2,1)
        rKwall_glob(:,5,3,2) = rKwall_glob32(:,2,2)
        rKwall_glob(:,6,3,2) = rKwall_glob32(:,2,3)
        rKwall_glob(:,7,3,2) = rKwall_glob32(:,3,1)
        rKwall_glob(:,8,3,2) = rKwall_glob32(:,3,2)
        rKwall_glob(:,9,3,2) = rKwall_glob32(:,3,3)
      
        rKwall_glob(:,1,3,3) = rKwall_glob33(:,1,1)
        rKwall_glob(:,2,3,3) = rKwall_glob33(:,1,2)
        rKwall_glob(:,3,3,3) = rKwall_glob33(:,1,3)
        rKwall_glob(:,4,3,3) = rKwall_glob33(:,2,1)
        rKwall_glob(:,5,3,3) = rKwall_glob33(:,2,2)
        rKwall_glob(:,6,3,3) = rKwall_glob33(:,2,3)
        rKwall_glob(:,7,3,3) = rKwall_glob33(:,3,1)
        rKwall_glob(:,8,3,3) = rKwall_glob33(:,3,2)
        rKwall_glob(:,9,3,3) = rKwall_glob33(:,3,3)

        rKwall_glob = rKwall_glob*betai*Delt(itseq)*Delt(itseq)*alfi
      
      else
c....   nothing happens
        goto 123
      endif      

123   continue

      endif
c     
c.... return
c     
      return
      end
      
c---------------------------------------------------------------------
c
c     variables for boundary elements
c
c---------------------------------------------------------------------
        subroutine e3bvarSclr (yl,        shdrv,    xlb,
     &                         shape,     WdetJb,   bnorm,
     &                         flux,      dwl )

        include "common.h"
c
        dimension yl(npro,nshl,ndof),        shdrv(npro,nsd,nshl),
     &            xlb(npro,nenl,nsd),        shape(npro,nshl),
     &            WdetJb(npro),              bnorm(npro,nsd),
     &            flux(npro)
c
        dimension dxdxib(npro,nsd,nsd),
     &            dxidxb(npro,nsd,nsd),      temp(npro),
     &            temp1(npro),               temp2(npro),
     &            temp3(npro),
     &            v1(npro,nsd),              v2(npro,nsd),
     &            gradSl(npro,nsd),          gradS(npro,nsd)

        real*8    diffus(npro),              dwl(npro,nshl)
        
        call getdiffsclr(shape,dwl,yl,diffus)
c
c.... ---------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
        dxdxib = zero
c
        do n = 1, nenl
           dxdxib(:,1,1) = dxdxib(:,1,1) + xlb(:,n,1) * shdrv(:,1,n)
           dxdxib(:,1,2) = dxdxib(:,1,2) + xlb(:,n,1) * shdrv(:,2,n)
           dxdxib(:,1,3) = dxdxib(:,1,3) + xlb(:,n,1) * shdrv(:,3,n)
           dxdxib(:,2,1) = dxdxib(:,2,1) + xlb(:,n,2) * shdrv(:,1,n)
           dxdxib(:,2,2) = dxdxib(:,2,2) + xlb(:,n,2) * shdrv(:,2,n)
           dxdxib(:,2,3) = dxdxib(:,2,3) + xlb(:,n,2) * shdrv(:,3,n)
           dxdxib(:,3,1) = dxdxib(:,3,1) + xlb(:,n,3) * shdrv(:,1,n)
           dxdxib(:,3,2) = dxdxib(:,3,2) + xlb(:,n,3) * shdrv(:,2,n)
           dxdxib(:,3,3) = dxdxib(:,3,3) + xlb(:,n,3) * shdrv(:,3,n)
        enddo
c     
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
c
        v1 = xlb(:,2,:) - xlb(:,1,:)
        v2 = xlb(:,3,:) - xlb(:,1,:)
        
c     
c.....The following are done in order to correct temp1..3  
c     based on the results from compressible code.  This is done only 
c     for wedges, depending on the bounary face.(tri or quad)  
c     
        if (lcsyst .eq. 4) then
           temp1 = dxdxib(:,2,1) * dxdxib(:,3,3) -
     &             dxdxib(:,2,3) * dxdxib(:,3,1)
           temp2 = dxdxib(:,3,1) * dxdxib(:,1,3) -
     &             dxdxib(:,3,3) * dxdxib(:,1,1)
           temp3 = dxdxib(:,1,1) * dxdxib(:,2,3) -
     &             dxdxib(:,1,3) * dxdxib(:,2,1)
             
        elseif (lcsyst .eq. 1) then
           temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
           temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
           temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
        else 
           temp1 = - v1(:,2) * v2(:,3) + v2(:,2) * v1(:,3)
           temp2 = - v2(:,1) * v1(:,3) + v1(:,1) * v2(:,3)
           temp3 = - v1(:,1) * v2(:,2) + v2(:,1) * v1(:,2)
        endif
c
        temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
        bnorm(:,1) = temp1 * temp
        bnorm(:,2) = temp2 * temp
        bnorm(:,3) = temp3 * temp
c
     
        if (lcsyst .eq. 3) then
           WdetJb     = (1 - Qwtb(lcsyst,intp)) / (four*temp)
        elseif (lcsyst .eq. 4) then
           WdetJb     = Qwtb(lcsyst,intp) / temp
        else
           WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
        endif      
c
c.... -------------------------->  Grad-V  <----------------------------
c
c.... compute grad-v for Navier-Stokes terms
c
        if (Navier .eq. 1) then
c
c.... compute the inverse of deformation gradient
c
          dxidxb(:,1,1) =   dxdxib(:,2,2) * dxdxib(:,3,3) 
     &                    - dxdxib(:,3,2) * dxdxib(:,2,3)
          dxidxb(:,1,2) =   dxdxib(:,3,2) * dxdxib(:,1,3) 
     &                    - dxdxib(:,1,2) * dxdxib(:,3,3)
          dxidxb(:,1,3) =   dxdxib(:,1,2) * dxdxib(:,2,3) 
     &                    - dxdxib(:,1,3) * dxdxib(:,2,2)
          temp          = one / ( dxidxb(:,1,1) * dxdxib(:,1,1) 
     &                          + dxidxb(:,1,2) * dxdxib(:,2,1)  
     &                          + dxidxb(:,1,3) * dxdxib(:,3,1) )
          dxidxb(:,1,1) =  dxidxb(:,1,1) * temp
          dxidxb(:,1,2) =  dxidxb(:,1,2) * temp
          dxidxb(:,1,3) =  dxidxb(:,1,3) * temp
          dxidxb(:,2,1) = (dxdxib(:,2,3) * dxdxib(:,3,1) 
     &                   - dxdxib(:,2,1) * dxdxib(:,3,3)) * temp
          dxidxb(:,2,2) = (dxdxib(:,1,1) * dxdxib(:,3,3) 
     &                   - dxdxib(:,3,1) * dxdxib(:,1,3)) * temp
          dxidxb(:,2,3) = (dxdxib(:,2,1) * dxdxib(:,1,3) 
     &                   - dxdxib(:,1,1) * dxdxib(:,2,3)) * temp
          dxidxb(:,3,1) = (dxdxib(:,2,1) * dxdxib(:,3,2) 
     &                   - dxdxib(:,2,2) * dxdxib(:,3,1)) * temp
          dxidxb(:,3,2) = (dxdxib(:,3,1) * dxdxib(:,1,2) 
     &                   - dxdxib(:,1,1) * dxdxib(:,3,2)) * temp
          dxidxb(:,3,3) = (dxdxib(:,1,1) * dxdxib(:,2,2) 
     &                   - dxdxib(:,1,2) * dxdxib(:,2,1)) * temp
c
c.... compute local-grad-Y
c
c
          gradSl = zero
          isc=5+isclr
          do n = 1, nshl
            gradSl(:,1) = gradSl(:,1) + shdrv(:,1,n) * yl(:,n,isc)
            gradSl(:,2) = gradSl(:,2) + shdrv(:,2,n) * yl(:,n,isc)
            gradSl(:,3) = gradSl(:,3) + shdrv(:,3,n) * yl(:,n,isc)
          enddo
c
c.... convert local-grads to global-grads
c
          gradS(:,1) = dxidxb(:,1,1) * gradSl(:,1) + 
     &                 dxidxb(:,2,1) * gradSl(:,2) + 
     &                 dxidxb(:,3,1) * gradSl(:,3)  

c
          gradS(:,2) = dxidxb(:,1,2) * gradSl(:,1) +
     &                 dxidxb(:,2,2) * gradSl(:,2) +
     &                 dxidxb(:,3,2) * gradSl(:,3) 

          gradS(:,3) = dxidxb(:,1,3) * gradSl(:,1) +
     &                 dxidxb(:,2,3) * gradSl(:,2) +
     &                 dxidxb(:,3,3) * gradSl(:,3) 
c
c.... end grad-T
c
        endif

        flux = diffus * ( gradS(:,1) * bnorm(:,1)
     &                  + gradS(:,2) * bnorm(:,2)
     &                  + gradS(:,3) * bnorm(:,3) )
c
c.... return
c
        return
        end


c---------------------------------------------------------------------
c
c     rotates the local nodal stiffnesses to the goblal frame
c
c---------------------------------------------------------------------
      subroutine rotatestiff(rKlocal, rotation, 
     &                       rKglobal)

      include "common.h"

      dimension rKlocal(npro,nsd,nsd), rotation(npro,nsd,nsd),
     &          rKglobal(npro,nsd,nsd)

      dimension tempm(npro,nsd,nsd)

      tempm = zero
      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               tempm(:,i,j) = tempm(:,i,j) 
     &              + rKlocal(:,i,k) * rotation(:,k,j)
            enddo
         enddo
      enddo
      
      rKglobal = zero
      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               rKglobal(:,i,j) = rKglobal(:,i,j) 
     &              + rotation(:,k,i) * tempm(:,k,j)
            enddo
         enddo
      enddo

      return
      end
