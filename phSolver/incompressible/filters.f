      subroutine hfilterC (y, x, ien, hfres, shgl, shp, Qwtf)

c...  The filter operator used here uses the generalized box 
c...  kernel


      include "common.h"

      dimension y(nshg,5),             hfres(nshg,16)
      dimension x(numnp,3),            xl(npro,nenl,3)
      dimension ien(npro,nshl),        yl(npro,nshl,5),
     &          fresl(npro,16),        WdetJ(npro),
     &          u1(npro),              u2(npro),
     &          u3(npro),              
     &          S11(npro),             S22(npro),
     &          S33(npro),             S12(npro),
     &          S13(npro),             S23(npro),
     &          dxdxi(npro,nsd,nsd),   dxidx(npro,nsd,nsd),
     &          shgl(nsd,nshl,ngauss),  shg(npro,nshl,nsd),
     &          shp(nshl,ngauss),       
     &          fresli(npro,16),       Qwtf(ngaussf)

      dimension tmp(npro)

      call local (y,      yl,     ien,    5,  'gather  ')
      call localx (x,      xl,     ien,    3,  'gather  ')
c

      fresl = zero

c
      do intp = 1, ngaussf   ! Loop over quadrature points

c  calculate the metrics
c
c
c.... --------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
         dxdxi = zero
c
         do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(1,n,intp)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(2,n,intp)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(3,n,intp)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(1,n,intp)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(2,n,intp)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(3,n,intp)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(1,n,intp)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(2,n,intp)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(3,n,intp)
         enddo
c     
c.... compute the inverse of deformation gradient
c
         dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)
     &        - dxdxi(:,3,2) * dxdxi(:,2,3)
         dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)
     &        - dxdxi(:,1,2) * dxdxi(:,3,3)
         dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)
     &        - dxdxi(:,1,3) * dxdxi(:,2,2)
         tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)
     &        + dxidx(:,1,2) * dxdxi(:,2,1)
     &        + dxidx(:,1,3) * dxdxi(:,3,1) )
         dxidx(:,1,1) = dxidx(:,1,1) * tmp
         dxidx(:,1,2) = dxidx(:,1,2) * tmp
         dxidx(:,1,3) = dxidx(:,1,3) * tmp
         dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)
     &        - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
         dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)
     &        - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
         dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)
     &        - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
         dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)
     &        - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
         dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)
     &        - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
         dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)
     &        - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c     
c        wght=Qwt(lcsyst,intp)  ! may be different now
         wght=Qwtf(intp)         
         WdetJ(:) = wght / tmp(:)
         

c... compute the gradient of shape functions at the quad. point.


      do n = 1,nshl 
        shg(:,n,1) = (shgl(1,n,intp) * dxidx(:,1,1)
     &              + shgl(2,n,intp) * dxidx(:,2,1)
     &              + shgl(3,n,intp) * dxidx(:,3,1))
        shg(:,n,2) = (shgl(1,n,intp) * dxidx(:,1,2)
     &              + shgl(2,n,intp) * dxidx(:,2,2)
     &              + shgl(3,n,intp) * dxidx(:,3,2))
        shg(:,n,3) = (shgl(1,n,intp) * dxidx(:,1,3)
     &              + shgl(2,n,intp) * dxidx(:,2,3)
     &              + shgl(3,n,intp) * dxidx(:,3,3))
      enddo


c... compute the velocities and the strain rate tensor at the quad. point


         u1  = zero
         u2  = zero
         u3  = zero
         S11 = zero
         S22 = zero
         S33 = zero
         S12 = zero
         S13 = zero
         S23 = zero
         do i=1,nshl  
            u1 = u1 + shp(i,intp)*yl(:,i,2)
            u2 = u2 + shp(i,intp)*yl(:,i,3)
            u3 = u3 + shp(i,intp)*yl(:,i,4)
 
            S11 = S11 + shg(:,i,1)*yl(:,i,2)
            S22 = S22 + shg(:,i,2)*yl(:,i,3)
            S33 = S33 + shg(:,i,3)*yl(:,i,4)
            
            S12 = S12 + shg(:,i,2)*yl(:,i,2)
     &                       +shg(:,i,1)*yl(:,i,3)
            S13 = S13 + shg(:,i,3)*yl(:,i,2)
     &                       +shg(:,i,1)*yl(:,i,4)
            S23 = S23 + shg(:,i,3)*yl(:,i,3)
     &                       +shg(:,i,2)*yl(:,i,4) 
         enddo
         S12 = pt5 * S12
         S13 = pt5 * S13
         S23 = pt5 * S23


         fresli(:,1) = WdetJ * u1 !G * u1 * WdetJ
         fresli(:,2) = WdetJ * u2 !G * u2 * WdetJ
         fresli(:,3) = WdetJ * u3 !G * u3 * WdetJ

         fresli(:,4) = WdetJ * u1 * u1 ! G*u1*u1*WdetJ
         fresli(:,5) = WdetJ * u2 * u2 ! G*u2*u2*WdetJ
         fresli(:,6) = WdetJ * u3 * u3 ! G*u3*u3*WdetJ
         fresli(:,7) = WdetJ * u1 * u2 ! G*u1*u2*WdetJ
         fresli(:,8) = WdetJ * u1 * u3 ! G*u1*u3*WdetJ
         fresli(:,9) = WdetJ * u2 * u3 ! G*u2*u3*WdetJ
         
         fresli(:,10) = S11 * WdetJ ! G*S_{1,1}*WdetJ
         fresli(:,11) = S22 * WdetJ ! G*S_{2,2}*WdetJ
         fresli(:,12) = S33 * WdetJ ! G*S_{3,3}*WdetJ
         fresli(:,13) = S12 * WdetJ ! G*S_{1,1}*WdetJ
         fresli(:,14) = S13 * WdetJ ! G*S_{1,3}*WdetJ
         fresli(:,15) = S23 * WdetJ ! G*S_{2,3}*WdetJ
            
         fresli(:,16) = WdetJ !Integral of filter kernel, G,
c                                               over the element            

         do i = 1, 16           ! Add contribution of each quad. point
            fresl(:,i) = fresl(:,i) + fresli(:,i)
         enddo

      enddo                     !end of loop over integration points
c

      do j = 1,nshl
      do nel = 1,npro
        hfres(ien(nel,j),:) = hfres(ien(nel,j),:) + fresl(nel,:) 
      enddo
      enddo


      return
      end

c... Here, the filter operation (denoted w/ a tilde) uses the generalized 
c... box kernel.

      subroutine twohfilterB (y, x, strnrm, ien, fres, 
     &     hfres, shgl, shp, Qwtf)

      include "common.h"

      dimension y(nshg,ndof),            fres(nshg,33)
      dimension x(numnp,nsd),            xl(npro,nenl,nsd)
      dimension ien(npro,nshl),        yl(npro,nshl,ndof),
     &          fresl(npro,33),        WdetJ(npro),
     &          u1(npro),              u2(npro),
     &          u3(npro),              dxdxi(npro,nsd,nsd),
     &          strnrm(npro,ngauss),    dxidx(npro,nsd,nsd),
     &          shgl(nsd,nshl,ngauss),       shg(npro,nshl,nsd),
     &          shp(nshl,ngauss),       
     &          fresli(npro,33),       Qwtf(ngaussf),
     &          hfres(nshg,16),        hfresl(npro,nshl,16),
     &          S(npro,nshl),          SS11(npro,nshl),
     &          SS22(npro,nshl),       SS33(npro,nshl),
     &          SS12(npro,nshl),       SS13(npro,nshl),
     &          SS23(npro,nshl),       
     &          S11p(npro),            S22p(npro),
     &          S33p(npro),            S12p(npro),
     &          S13p(npro),            S23p(npro)

      dimension tmp(npro)

      call local (y,      yl,     ien,    5,  'gather  ')
      call localx (x,      xl,     ien,    3,  'gather  ')
      call local (hfres,  hfresl, ien,   16,  'gather  ')

      S(:,:) = sqrt(
     &     two*(hfresl(:,:,10)**2 + hfresl(:,:,11)**2 +
     &     hfresl(:,:,12)**2) + four*(
     &     hfresl(:,:,13)**2 + hfresl(:,:,14)**2 +
     &     hfresl(:,:,15)**2) )      

      SS11(:,:) = S(:,:)*hfresl(:,:,10)
      SS22(:,:) = S(:,:)*hfresl(:,:,11)
      SS33(:,:) = S(:,:)*hfresl(:,:,12)
      SS12(:,:) = S(:,:)*hfresl(:,:,13)
      SS13(:,:) = S(:,:)*hfresl(:,:,14)
      SS23(:,:) = S(:,:)*hfresl(:,:,15)

      fresl = zero

      do intp = 1, ngauss


c  calculate the metrics
c
c
c.... --------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
        dxdxi = zero
c
          do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(1,n,intp)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(2,n,intp)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(3,n,intp)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(1,n,intp)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(2,n,intp)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(3,n,intp)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(1,n,intp)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(2,n,intp)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(3,n,intp)
          enddo
c
c.... compute the inverse of deformation gradient
c
        dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)
     &                 - dxdxi(:,3,2) * dxdxi(:,2,3)
        dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)
     &                 - dxdxi(:,1,2) * dxdxi(:,3,3)
        dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)
     &                 - dxdxi(:,1,3) * dxdxi(:,2,2)
        tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)
     &                       + dxidx(:,1,2) * dxdxi(:,2,1)
     &                       + dxidx(:,1,3) * dxdxi(:,3,1) )
        dxidx(:,1,1) = dxidx(:,1,1) * tmp
        dxidx(:,1,2) = dxidx(:,1,2) * tmp
        dxidx(:,1,3) = dxidx(:,1,3) * tmp
        dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)
     &                - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
        dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)
     &                - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
        dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)
     &                - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
        dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)
     &                - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
        dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)
     &                - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
        dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)
     &                - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
        wght=Qwt(lcsyst,intp)  ! may be different now
c        wght=Qwtf(intp)
        WdetJ = wght / tmp
c

      do n = 1,nshl
        shg(:,n,1) = (shgl(1,n,intp) * dxidx(:,1,1)
     &              + shgl(2,n,intp) * dxidx(:,2,1)
     &              + shgl(3,n,intp) * dxidx(:,3,1))
        shg(:,n,2) = (shgl(1,n,intp) * dxidx(:,1,2)
     &              + shgl(2,n,intp) * dxidx(:,2,2)
     &              + shgl(3,n,intp) * dxidx(:,3,2))
        shg(:,n,3) = (shgl(1,n,intp) * dxidx(:,1,3)
     &              + shgl(2,n,intp) * dxidx(:,2,3)
     &              + shgl(3,n,intp) * dxidx(:,3,3))
      enddo

c... compute filtered quantities at the hat level evaluated at the quad. pt.
c... This should really
c... be the bar-hat level, but the bar filter is implicit so we don't
c... bother to mention it.

      fresli = zero
      S11p   = zero
      S22p   = zero
      S33p   = zero
      S12p   = zero
      S13p   = zero
      S23p   = zero

      do i = 1, nenl
         fresli(:,1) = fresli(:,1) + shp(i,intp)*hfresl(:,i,1)  ! hat{u1}
         fresli(:,2) = fresli(:,2) + shp(i,intp)*hfresl(:,i,2)  ! hat{u2}
         fresli(:,3) = fresli(:,3) + shp(i,intp)*hfresl(:,i,3)  ! hat{u3}

         fresli(:,4) = fresli(:,4) + shp(i,intp)*hfresl(:,i,1)*
     &        hfresl(:,i,1)     ! hat{u1}*hat{u1}
         fresli(:,5) = fresli(:,5) + shp(i,intp)*hfresl(:,i,2)*
     &        hfresl(:,i,2)     ! hat{u2}*hat{u2}
         fresli(:,6) = fresli(:,6) + shp(i,intp)*hfresl(:,i,3)*
     &        hfresl(:,i,3)     ! hat{u3}*hat{u3}
         fresli(:,7) = fresli(:,7) + shp(i,intp)*hfresl(:,i,1)*
     &        hfresl(:,i,2)     ! hat{u1}*hat{u2}
         fresli(:,8) = fresli(:,8) + shp(i,intp)*hfresl(:,i,1)*
     &        hfresl(:,i,3)     ! hat{u1}*hat{u3}
         fresli(:,9) = fresli(:,9) + shp(i,intp)*hfresl(:,i,2)*
     &        hfresl(:,i,3)     ! hat{u2}*hat{u3}

         fresli(:,10) = fresli(:,10) + shp(i,intp)*hfresl(:,i,10)  ! hat{S11}
         fresli(:,11) = fresli(:,11) + shp(i,intp)*hfresl(:,i,11)  ! hat{S22}
         fresli(:,12) = fresli(:,12) + shp(i,intp)*hfresl(:,i,12)  ! hat{S33}
         fresli(:,13) = fresli(:,13) + shp(i,intp)*hfresl(:,i,13)  ! hat{S12}
         fresli(:,14) = fresli(:,14) + shp(i,intp)*hfresl(:,i,14)  ! hat{S13}
         fresli(:,15) = fresli(:,15) + shp(i,intp)*hfresl(:,i,15)  ! hat{S23}

         S11p = S11p + shg(:,i,1)*hfresl(:,i,1)
         S22p = S22p + shg(:,i,2)*hfresl(:,i,2)
         S33p = S33p + shg(:,i,3)*hfresl(:,i,3)

         S12p = S12p + shg(:,i,2)*hfresl(:,i,1) +
     &        shg(:,i,1)*hfresl(:,i,2)
         S13p = S13p + shg(:,i,3)*hfresl(:,i,1) +
     &        shg(:,i,1)*hfresl(:,i,3)
         S23p = S23p + shg(:,i,3)*hfresl(:,i,2) +
     &        shg(:,i,2)*hfresl(:,i,3)

      enddo         

c... get the strain rate tensor norm at the quad. pt.

c      strnrm(:,intp) = sqrt(
c     &   two * (fresli(:,10)**2 + fresli(:,11)**2 + fresli(:,12)**2)
c     &  + four * ( fresli(:,13)**2 + fresli(:,14)**2 + 
c     &    fresli(:,15)**2 ) )

c      strnrm2(:,intp) = zero
c      do i = 1, nenl
c         strnrm2(:,intp) = strnrm2(:,intp) + S(:,i)*shp(i,intp)
c      enddo

c      strnrm3(:,intp) = sqrt(
c     &     two * (S11p(:)**2 + S22p(:)**2 + S33p(:)**2)
c     &    + four * ( S12p(:)**2 + S13p(:)**2 + S23p(:)**2) ) 

c... compute |hat{S}| * hat{Sij} 

      do i = 1, nenl
         fresli(:,16) =fresli(:,16)+shp(i,intp)*SS11(:,i)! |hat{S}|*hat{S11} 
         fresli(:,17) =fresli(:,17)+shp(i,intp)*SS22(:,i)! |hat{S}|*hat{S22}
         fresli(:,18) =fresli(:,18)+shp(i,intp)*SS33(:,i)! |hat{S}|*hat{S33}
         fresli(:,19) =fresli(:,19)+shp(i,intp)*SS12(:,i)! |hat{S}|*hat{S12}
         fresli(:,20) =fresli(:,20)+shp(i,intp)*SS13(:,i)! |hat{S}|*hat{S13}
         fresli(:,21) =fresli(:,21)+shp(i,intp)*SS23(:,i)! |hat{S}|*hat{S23}
      enddo

c... multiply fresli by WdetJ so that when we finish looping over 
c... quad. pts. and add the contributions from all the quad. points
c... we get filtered quantities at the hat-tilde level, secretly the
c... bar-hat-tilde level.

      do j = 1, 21
         fresli(:,j) = fresli(:,j)*WdetJ(:)
      enddo

c... compute volume of box kernel

      fresli(:,22) = WdetJ

c... add contributions from each quad pt to current element 

      do i = 1, 22
         fresl(:,i) = fresl(:,i) + fresli(:,i)
      enddo

      enddo ! end of loop over integration points


c... scatter locally filtered quantities to the global nodes

      do j = 1,nshl
      do nel = 1,npro
        fres(ien(nel,j),:) = fres(ien(nel,j),:) + fresl(nel,:) 
      enddo
      enddo      

      return
      end

c...  The filter operator used here uses the generalized box 
c...  kernel

      subroutine hfilterCC (y, x, ien, hfres, shgl, shp, Qwtf)

      include "common.h"

      dimension y(nshg,5),             hfres(nshg,22)
      dimension x(numnp,3),            xl(npro,nenl,3)
      dimension ien(npro,nshl),        yl(npro,nshl,5),
     &          fresl(npro,22),        WdetJ(npro),
     &          u1(npro),              u2(npro),
     &          u3(npro),              
     &          S11(npro),             S22(npro),
     &          S33(npro),             S12(npro),
     &          S13(npro),             S23(npro),
     &          dxdxi(npro,nsd,nsd),   dxidx(npro,nsd,nsd),
     &          shgl(nsd,nshl,ngauss),  shg(npro,nshl,nsd),
     &          shp(nshl,ngauss),       
     &          fresli(npro,22),       Qwtf(ngaussf),
     &          strnrmi(npro) 

      dimension tmp(npro)

      call local (y,      yl,     ien,    5,  'gather  ')
      call localx (x,      xl,     ien,    3,  'gather  ')
c

      fresl = zero

c
      do intp = 1, ngaussf   ! Loop over quadrature points

c  calculate the metrics
c
c
c.... --------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
         dxdxi = zero
c
         do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(1,n,intp)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(2,n,intp)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(3,n,intp)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(1,n,intp)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(2,n,intp)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(3,n,intp)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(1,n,intp)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(2,n,intp)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(3,n,intp)
         enddo
c     
c.... compute the inverse of deformation gradient
c
         dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)
     &        - dxdxi(:,3,2) * dxdxi(:,2,3)
         dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)
     &        - dxdxi(:,1,2) * dxdxi(:,3,3)
         dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)
     &        - dxdxi(:,1,3) * dxdxi(:,2,2)
         tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)
     &        + dxidx(:,1,2) * dxdxi(:,2,1)
     &        + dxidx(:,1,3) * dxdxi(:,3,1) )
         dxidx(:,1,1) = dxidx(:,1,1) * tmp
         dxidx(:,1,2) = dxidx(:,1,2) * tmp
         dxidx(:,1,3) = dxidx(:,1,3) * tmp
         dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)
     &        - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
         dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)
     &        - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
         dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)
     &        - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
         dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)
     &        - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
         dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)
     &        - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
         dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)
     &        - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c     
c         wght=Qwt(lcsyst,intp)  ! may be different now
         wght=Qwtf(intp)         
         WdetJ(:) = wght / tmp(:)
         

c... compute the gradient of shape functions at the quad. point.


      do n = 1,nshl 
        shg(:,n,1) = (shgl(1,n,intp) * dxidx(:,1,1)
     &              + shgl(2,n,intp) * dxidx(:,2,1)
     &              + shgl(3,n,intp) * dxidx(:,3,1))
        shg(:,n,2) = (shgl(1,n,intp) * dxidx(:,1,2)
     &              + shgl(2,n,intp) * dxidx(:,2,2)
     &              + shgl(3,n,intp) * dxidx(:,3,2))
        shg(:,n,3) = (shgl(1,n,intp) * dxidx(:,1,3)
     &              + shgl(2,n,intp) * dxidx(:,2,3)
     &              + shgl(3,n,intp) * dxidx(:,3,3))
      enddo


c... compute the velocities and the strain rate tensor at the quad. point


         u1  = zero
         u2  = zero
         u3  = zero
         S11 = zero
         S22 = zero
         S33 = zero
         S12 = zero
         S13 = zero
         S23 = zero
         do i=1,nshl  
            u1 = u1 + shp(i,intp)*yl(:,i,2)
            u2 = u2 + shp(i,intp)*yl(:,i,3)
            u3 = u3 + shp(i,intp)*yl(:,i,4)
 
            S11 = S11 + shg(:,i,1)*yl(:,i,2)
            S22 = S22 + shg(:,i,2)*yl(:,i,3)
            S33 = S33 + shg(:,i,3)*yl(:,i,4)
            
            S12 = S12 + shg(:,i,2)*yl(:,i,2)
     &                       +shg(:,i,1)*yl(:,i,3)
            S13 = S13 + shg(:,i,3)*yl(:,i,2)
     &                       +shg(:,i,1)*yl(:,i,4)
            S23 = S23 + shg(:,i,3)*yl(:,i,3)
     &                       +shg(:,i,2)*yl(:,i,4) 
         enddo
         S12 = pt5 * S12
         S13 = pt5 * S13
         S23 = pt5 * S23

c... Get the strain rate norm at the quad pts

         strnrmi = sqrt( two*(S11**2 + S22**2 + S33**2) 
     &        + four*(S12**2 + S13**2 + S23**2) )


c... Multiply quantities to be filtered by the box kernel

         fresli(:,1) = WdetJ * u1 !G * u1 * WdetJ
         fresli(:,2) = WdetJ * u2 !G * u2 * WdetJ
         fresli(:,3) = WdetJ * u3 !G * u3 * WdetJ

         fresli(:,4) = WdetJ * u1 * u1 ! G*u1*u1*WdetJ
         fresli(:,5) = WdetJ * u2 * u2 ! G*u2*u2*WdetJ
         fresli(:,6) = WdetJ * u3 * u3 ! G*u3*u3*WdetJ
         fresli(:,7) = WdetJ * u1 * u2 ! G*u1*u2*WdetJ
         fresli(:,8) = WdetJ * u1 * u3 ! G*u1*u3*WdetJ
         fresli(:,9) = WdetJ * u2 * u3 ! G*u2*u3*WdetJ
         
         fresli(:,10) = S11 * WdetJ ! G*S_{1,1}*WdetJ
         fresli(:,11) = S22 * WdetJ ! G*S_{2,2}*WdetJ
         fresli(:,12) = S33 * WdetJ ! G*S_{3,3}*WdetJ
         fresli(:,13) = S12 * WdetJ ! G*S_{1,1}*WdetJ
         fresli(:,14) = S13 * WdetJ ! G*S_{1,3}*WdetJ
         fresli(:,15) = S23 * WdetJ ! G*S_{2,3}*WdetJ
            
         fresli(:,16) = WdetJ !Integral of filter kernel, G,
c                                               over the element            

         fresli(:,17) = S11 * strnrmi * WdetJ
         fresli(:,18) = S22 * strnrmi * WdetJ
         fresli(:,19) = S33 * strnrmi * WdetJ
         fresli(:,20) = S12 * strnrmi * WdetJ         
         fresli(:,21) = S13 * strnrmi * WdetJ
         fresli(:,22) = S23 * strnrmi * WdetJ


         do i = 1, 22           ! Add contribution of each quad. point
            fresl(:,i) = fresl(:,i) + fresli(:,i)
         enddo

      enddo                     !end of loop over integration points
c

      do j = 1,nshl
      do nel = 1,npro
        hfres(ien(nel,j),:) = hfres(ien(nel,j),:) + fresl(nel,:) 
      enddo
      enddo


      return
      end


c... Here, the filter operation (denoted w/ a tilde) uses the generalized 
c... box kernel.

      subroutine twohfilterBB (y, x, strnrm, ien, fres, 
     &     hfres, shgl, shp, Qwtf)

      include "common.h"

      dimension y(nshg,ndof),            fres(nshg,33)
      dimension x(numnp,nsd),            xl(npro,nenl,nsd)
      dimension ien(npro,nshl),        yl(npro,nshl,ndof),
     &          fresl(npro,33),        WdetJ(npro),
     &          u1(npro),              u2(npro),
     &          u3(npro),              dxdxi(npro,nsd,nsd),
     &          strnrm(npro,ngauss),    dxidx(npro,nsd,nsd),
     &          shgl(nsd,nshl,ngauss),       shg(npro,nshl,nsd),
     &          shp(nshl,ngauss),       
     &          fresli(npro,33),       Qwtf(ngaussf),
     &          hfres(nshg,22),        hfresl(npro,nshl,22),
     &          S(npro,nshl),          SS11(npro,nshl),
     &          SS22(npro,nshl),       SS33(npro,nshl),
     &          SS12(npro,nshl),       SS13(npro,nshl),
     &          SS23(npro,nshl)       

      dimension tmp(npro)

      call local (y,      yl,     ien,    5,  'gather  ')
      call localx (x,      xl,     ien,    3,  'gather  ')
      call local (hfres,  hfresl, ien,   22,  'gather  ')

      S(:,:) = sqrt(
     &     two*(hfresl(:,:,10)**2 + hfresl(:,:,11)**2 +
     &     hfresl(:,:,12)**2) + four*(
     &     hfresl(:,:,13)**2 + hfresl(:,:,14)**2 +
     &     hfresl(:,:,15)**2) )      

      SS11(:,:) = S(:,:)*hfresl(:,:,10)
      SS22(:,:) = S(:,:)*hfresl(:,:,11)
      SS33(:,:) = S(:,:)*hfresl(:,:,12)
      SS12(:,:) = S(:,:)*hfresl(:,:,13)
      SS13(:,:) = S(:,:)*hfresl(:,:,14)
      SS23(:,:) = S(:,:)*hfresl(:,:,15)

      fresl = zero

      do intp = 1, ngaussf


c  calculate the metrics
c
c
c.... --------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
        dxdxi = zero
c
          do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(1,n,intp)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(2,n,intp)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(3,n,intp)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(1,n,intp)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(2,n,intp)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(3,n,intp)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(1,n,intp)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(2,n,intp)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(3,n,intp)
          enddo
c
c.... compute the inverse of deformation gradient
c
        dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)
     &                 - dxdxi(:,3,2) * dxdxi(:,2,3)
        dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)
     &                 - dxdxi(:,1,2) * dxdxi(:,3,3)
        dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)
     &                 - dxdxi(:,1,3) * dxdxi(:,2,2)
        tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)
     &                       + dxidx(:,1,2) * dxdxi(:,2,1)
     &                       + dxidx(:,1,3) * dxdxi(:,3,1) )
        dxidx(:,1,1) = dxidx(:,1,1) * tmp
        dxidx(:,1,2) = dxidx(:,1,2) * tmp
        dxidx(:,1,3) = dxidx(:,1,3) * tmp
        dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)
     &                - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
        dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)
     &                - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
        dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)
     &                - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
        dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)
     &                - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
        dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)
     &                - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
        dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)
     &                - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
c        wght=Qwt(lcsyst,intp)  ! may be different now
        wght=Qwtf(intp)
        WdetJ = wght / tmp
c

      do n = 1,nshl
        shg(:,n,1) = (shgl(1,n,intp) * dxidx(:,1,1)
     &              + shgl(2,n,intp) * dxidx(:,2,1)
     &              + shgl(3,n,intp) * dxidx(:,3,1))
        shg(:,n,2) = (shgl(1,n,intp) * dxidx(:,1,2)
     &              + shgl(2,n,intp) * dxidx(:,2,2)
     &              + shgl(3,n,intp) * dxidx(:,3,2))
        shg(:,n,3) = (shgl(1,n,intp) * dxidx(:,1,3)
     &              + shgl(2,n,intp) * dxidx(:,2,3)
     &              + shgl(3,n,intp) * dxidx(:,3,3))
      enddo

c... compute filtered quantities at the hat level evaluated at the quad. pt.
c... This should really
c... be the bar-hat level, but the bar filter is implicit so we don't
c... bother to mention it.

      fresli = zero
      S11p   = zero
      S22p   = zero
      S33p   = zero
      S12p   = zero
      S13p   = zero
      S23p   = zero

      do i = 1, nenl
         fresli(:,1) = fresli(:,1) + shp(i,intp)*hfresl(:,i,1)  ! hat{u1}
         fresli(:,2) = fresli(:,2) + shp(i,intp)*hfresl(:,i,2)  ! hat{u2}
         fresli(:,3) = fresli(:,3) + shp(i,intp)*hfresl(:,i,3)  ! hat{u3}

         fresli(:,4) = fresli(:,4) + shp(i,intp)*hfresl(:,i,4) ! hat{u1*u1}
         fresli(:,5) = fresli(:,5) + shp(i,intp)*hfresl(:,i,5) ! hat{u2*u2}
         fresli(:,6) = fresli(:,6) + shp(i,intp)*hfresl(:,i,6) ! hat{u3*u3}
         fresli(:,7) = fresli(:,7) + shp(i,intp)*hfresl(:,i,7) ! hat{u1*u2}
         fresli(:,8) = fresli(:,8) + shp(i,intp)*hfresl(:,i,8) ! hat{u1*u3}
         fresli(:,9) = fresli(:,9) + shp(i,intp)*hfresl(:,i,9) ! hat{u2*u3}

         fresli(:,10) = fresli(:,10) + shp(i,intp)*hfresl(:,i,10)  ! hat{S11}
         fresli(:,11) = fresli(:,11) + shp(i,intp)*hfresl(:,i,11)  ! hat{S22}
         fresli(:,12) = fresli(:,12) + shp(i,intp)*hfresl(:,i,12)  ! hat{S33}
         fresli(:,13) = fresli(:,13) + shp(i,intp)*hfresl(:,i,13)  ! hat{S12}
         fresli(:,14) = fresli(:,14) + shp(i,intp)*hfresl(:,i,14)  ! hat{S13}
         fresli(:,15) = fresli(:,15) + shp(i,intp)*hfresl(:,i,15)  ! hat{S23}

         fresli(:,16) = fresli(:,16) +shp(i,intp)*hfresl(:,i,17)! hat{S11*|S|}
         fresli(:,17) = fresli(:,17) +shp(i,intp)*hfresl(:,i,18)! hat{S22*|S|}
         fresli(:,18) = fresli(:,18) +shp(i,intp)*hfresl(:,i,19)! hat{S33*|S|}
         fresli(:,19) = fresli(:,19) +shp(i,intp)*hfresl(:,i,20)! hat{S12*|S|}
         fresli(:,20) = fresli(:,20) +shp(i,intp)*hfresl(:,i,21)! hat{S13*|S|}
         fresli(:,21) = fresli(:,21) +shp(i,intp)*hfresl(:,i,22)! hat{S23*|S|}


      enddo

c... multiply fresli by WdetJ so that when we finish looping over 
c... quad. pts. and add the contributions from all the quad. points
c... we get filtered quantities at the hat-tilde level, secretly the
c... bar-hat-tilde level.

      do j = 1, 21
         fresli(:,j) = fresli(:,j)*WdetJ(:)
      enddo

c... compute volume of box kernel

      fresli(:,22) = WdetJ

c... add contributions from each quad pt to current element 

      do i = 1, 22
         fresl(:,i) = fresl(:,i) + fresli(:,i)
      enddo

      enddo ! end of loop over integration points


c... scatter locally filtered quantities to the global nodes

      do j = 1,nshl
      do nel = 1,npro
        fres(ien(nel,j),:) = fres(ien(nel,j),:) + fresl(nel,:) 
      enddo
      enddo      

      return
      end


c...  The filter operator used here uses the generalized hat (witch hat) 
c...  kernel

      subroutine hfilterBB (y, x, ien, hfres, shgl, shp, Qwtf)

      include "common.h"

      dimension y(nshg,5),             hfres(nshg,24)
      dimension x(numnp,3),            xl(npro,nenl,3)
      dimension ien(npro,nshl),        yl(npro,nshl,5),
     &          fresl(npro,nshl,24),        WdetJ(npro),
     &          u1(npro),              u2(npro),
     &          u3(npro),              rho(npro),
     &          S11(npro),             S22(npro),
     &          S33(npro),             S12(npro),
     &          S13(npro),             S23(npro),
     &          dxdxi(npro,nsd,nsd),   dxidx(npro,nsd,nsd),
     &          shgl(nsd,nshl,ngauss),       shg(npro,nshl,nsd),
     &          shp(nshl,ngauss),          
     &          fresli(npro,nshl,24),   Qwtf(ngaussf),
     &          strnrmi(npro) 
    

      dimension tmp(npro)

      call local (y,      yl,     ien,    5,  'gather  ')
      call localx (x,      xl,     ien,    3,  'gather  ')
c

      fresl = zero

c
      do intp = 1, ngaussf   ! Loop over quadrature points

c  calculate the metrics
c
c
c.... --------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
         dxdxi = zero
c
         do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(1,n,intp)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(2,n,intp)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(3,n,intp)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(1,n,intp)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(2,n,intp)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(3,n,intp)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(1,n,intp)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(2,n,intp)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(3,n,intp)
         enddo
c     
c.... compute the inverse of deformation gradient
c
         dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)
     &        - dxdxi(:,3,2) * dxdxi(:,2,3)
         dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)
     &        - dxdxi(:,1,2) * dxdxi(:,3,3)
         dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)
     &        - dxdxi(:,1,3) * dxdxi(:,2,2)
         tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)
     &        + dxidx(:,1,2) * dxdxi(:,2,1)
     &        + dxidx(:,1,3) * dxdxi(:,3,1) )
         dxidx(:,1,1) = dxidx(:,1,1) * tmp
         dxidx(:,1,2) = dxidx(:,1,2) * tmp
         dxidx(:,1,3) = dxidx(:,1,3) * tmp
         dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)
     &        - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
         dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)
     &        - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
         dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)
     &        - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
         dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)
     &        - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
         dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)
     &        - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
         dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)
     &        - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c     
c        wght=Qwt(lcsyst,intp)  ! may be different now
         wght=Qwtf(intp)    
         WdetJ(:) = wght / tmp(:)
         

c... compute the gradient of shape functions at the quad. point.


      do n = 1,nshl 
        shg(:,n,1) = (shgl(1,n,intp) * dxidx(:,1,1)
     &              + shgl(2,n,intp) * dxidx(:,2,1)
     &              + shgl(3,n,intp) * dxidx(:,3,1))
        shg(:,n,2) = (shgl(1,n,intp) * dxidx(:,1,2)
     &              + shgl(2,n,intp) * dxidx(:,2,2)
     &              + shgl(3,n,intp) * dxidx(:,3,2))
        shg(:,n,3) = (shgl(1,n,intp) * dxidx(:,1,3)
     &              + shgl(2,n,intp) * dxidx(:,2,3)
     &              + shgl(3,n,intp) * dxidx(:,3,3))
      enddo


c... compute the velocities and the strain rate tensor at the quad. point


         u1  = zero
         u2  = zero
         u3  = zero
         S11 = zero
         S22 = zero
         S33 = zero
         S12 = zero
         S13 = zero
         S23 = zero
         do i=1,nshl  
            u1 = u1 + shp(i,intp)*yl(:,i,2)
            u2 = u2 + shp(i,intp)*yl(:,i,3)
            u3 = u3 + shp(i,intp)*yl(:,i,4)
 
            S11 = S11 + shg(:,i,1)*yl(:,i,2)
            S22 = S22 + shg(:,i,2)*yl(:,i,3)
            S33 = S33 + shg(:,i,3)*yl(:,i,4)
            
            S12 = S12 + shg(:,i,2)*yl(:,i,2)
     &                       +shg(:,i,1)*yl(:,i,3)
            S13 = S13 + shg(:,i,3)*yl(:,i,2)
     &                       +shg(:,i,1)*yl(:,i,4)
            S23 = S23 + shg(:,i,3)*yl(:,i,3)
     &                       +shg(:,i,2)*yl(:,i,4) 
         enddo
         S12 = pt5 * S12
         S13 = pt5 * S13
         S23 = pt5 * S23

c... Get the strain rate norm at the quad pts

         strnrmi = sqrt( two*(S11**2 + S22**2 + S33**2) 
     &        + four*(S12**2 + S13**2 + S23**2) )


c... Loop over element nodes and multiply u_{i} and S_{i,j} by the
c... hat kernel and the Jacobian over the current element evaluated at 
c... the current quad. point.

         do i = 1, nenl   ! Loop over element nodes

            fresli(:,i,1) = WdetJ * u1 * shp(i,intp) !G * u1 * WdetJ
            fresli(:,i,2) = WdetJ * u2 * shp(i,intp) !G * u2 * WdetJ
            fresli(:,i,3) = WdetJ * u3 * shp(i,intp) !G * u3 * WdetJ

            fresli(:,i,4) = WdetJ * u1 * u1 * shp(i,intp) ! G*u1*u1*WdetJ
            fresli(:,i,5) = WdetJ * u2 * u2 * shp(i,intp) ! G*u2*u2*WdetJ
            fresli(:,i,6) = WdetJ * u3 * u3 * shp(i,intp) ! G*u3*u3*WdetJ
            fresli(:,i,7) = WdetJ * u1 * u2 * shp(i,intp) ! G*u1*u2*WdetJ
            fresli(:,i,8) = WdetJ * u1 * u3 * shp(i,intp) ! G*u1*u3*WdetJ
            fresli(:,i,9) = WdetJ * u2 * u3 * shp(i,intp) ! G*u2*u3*WdetJ

            fresli(:,i,10) = S11 * shp(i,intp) * WdetJ ! G*S_{1,1}*WdetJ
            fresli(:,i,11) = S22 * shp(i,intp) * WdetJ ! G*S_{2,2}*WdetJ
            fresli(:,i,12) = S33 * shp(i,intp) * WdetJ ! G*S_{3,3}*WdetJ
            fresli(:,i,13) = S12 * shp(i,intp) * WdetJ ! G*S_{1,1}*WdetJ
            fresli(:,i,14) = S13 * shp(i,intp) * WdetJ ! G*S_{1,3}*WdetJ
            fresli(:,i,15) = S23 * shp(i,intp) * WdetJ ! G*S_{2,3}*WdetJ
            
            fresli(:,i,22) = WdetJ*shp(i,intp) !Integral of filter kernel, G,
c                                               over the element 

c...   Get G*|S|*S_{i,j}*WdetJ
           
            fresli(:,i,16) = S11 * strnrmi * shp(i,intp) * WdetJ 
            fresli(:,i,17) = S22 * strnrmi * shp(i,intp) * WdetJ            
            fresli(:,i,18) = S33 * strnrmi * shp(i,intp) * WdetJ
            fresli(:,i,19) = S12 * strnrmi * shp(i,intp) * WdetJ  
            fresli(:,i,20) = S13 * strnrmi * shp(i,intp) * WdetJ
            fresli(:,i,21) = S23 * strnrmi * shp(i,intp) * WdetJ  

            do j = 1, 22 ! Add contribution of each quad. point for each
c                          element node
               fresl(:,i,j) = fresl(:,i,j) + fresli(:,i,j)
            enddo

         enddo                  !end loop over element nodes

      enddo                     !end of loop over integration points


      call local (hfres, fresl, ien, 24, 'scatter ')

      return
      end

      subroutine ediss (y,           ac,         shgl,      
     &                  shp,         iper,       ilwork,    
     &                  nsons,       ifath,      x,
     &                  iBC,    BC,  xavegt)


      use pvsQbi           ! brings in NABI
      use stats            !  
      use pointer_data     ! brings in the pointers for the blocked arrays
      use local_mass
      use rlssave  ! Use the resolved Leonard stresses at the nodes.      
      use quadfilt ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwt contains the weights of the 
c                    quad. points.  



      include "common.h"
      include "mpif.h"
      include "auxmpi.h"


      dimension y(nshg,ndof),                  ac(nshg,ndof), 
     &          ifath(nshg),                   nsons(nshg),
     &          iper(nshg),                    ilwork(nlwork),        
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),
     &          x(numnp,3),           
     &          qres(nshg,nsd*nsd),             rmass(nshg),  
     &          iBC(nshg),                      BC(nshg,ndofBC),
     &          cdelsq(nshg),                   vol(nshg),
     &          stress(nshg,9),                 diss(nshg,3),
     &          xave(nshg,12),                  xaveg(nfath,12),
     &          xavegr(nfath,12),           xavegt(nfath,12),
     &          rnum(nfath),                rden(nfath)                
 

c.... First let us obtain cdelsq at each node in the domain.
c.... We use numNden which lives in the quadfilt module.

      rnum(ifath(:)) = numNden(:,1)
      rden(ifath(:)) = numNden(:,2)
      
c      if (myrank .eq. master) then
c         write(*,*) 'rnum25=', rnum(25), rden(25)
c         write(*,*) 'rnum26=', rnum(26), rden(26)
c      endif

      cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)       
c      cdelsq(:) = zero ! Debugging

      if (istep .eq. 0) then
         xavegt = zero  ! For averaging dissipations and SUPG stresses
      endif

        if (idiff==1 .or. idiff==3) then ! global reconstruction of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
           qres = zero
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)
c     
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

              call AsIq (y,                x,                       
     &                   shp(lcsyst,1:nshl,:), 
     &                   shgl(lcsyst,:,1:nshl,:),
     &                   mien(iblk)%p,     mxmudmi(iblk)%p,  
     &                   qres,             rmass )
           enddo
       
c
c.... form the diffusive flux approximation
c
           call qpbc( rmass, qres, iBC, BC, iper, ilwork )       
c
        endif 


c.... form the SUPG stresses well as dissipation due to eddy viscosity,
c...  and SUPG stabilization.


        stress = zero
        vol    = zero
        diss   = zero

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        call SUPGstress (y, ac, x, qres, mien(iblk)%p, mxmudmi(iblk)%p, 
     &                   cdelsq, shglf(lcsyst,:,1:nshl,:),
     &                   shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf),
     &                   shgl(lcsyst,:,1:nshl,:), shp(lcsyst,1:nshl,:),
     &                   stress, diss, vol)

      enddo

      if(numpe>1) call commu (stress, ilwork, 9, 'in ')      
      if(numpe>1) call commu (diss, ilwork, 3, 'in ')       
      if(numpe>1) call commu (vol, ilwork, 1, 'in ')      

c 
c account for periodicity 
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           stress(i,:) = stress(i,:) + stress(j,:)
           diss(i,:)   = diss(i,:)   + diss(j,:)
           vol(i)      = vol(i)      + vol(j)
        endif
      enddo

      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           stress(j,:) = stress(i,:) 
           diss(j,:)   = diss(i,:)   
           vol(j)      = vol(i)      
        endif
      enddo      

      if(numpe>1) call commu (stress, ilwork, 9, 'out ')      
      if(numpe>1) call commu (diss, ilwork, 3, 'out ')       
      if(numpe>1) call commu (vol, ilwork, 1, 'out ')      

      vol = one / vol
      do i = 1, 9
         stress(:,i) = stress(:,i)*vol(:)
      enddo
      do i = 1, 3
         diss(:,i) = diss(:,i)*vol(:)
      enddo

c---------- > Begin averaging dissipations and SUPG stress <--------------

      do i = 1, 9
         xave(:,i) = stress(:,i)
      enddo
      xave(:,10) = diss(:,1)
      xave(:,11) = diss(:,2)
      xave(:,12) = diss(:,3)

c  zero on processor periodic nodes so that they will not be added twice
        do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xave(j,:) = zero
          endif
        enddo

      if (numpe.gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xave(isgbeg:isgend,:) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
c

      xaveg = zero
      do i = 1,nshg      
         xaveg(ifath(i),:) = xaveg(ifath(i),:) + xave(i,:)
      enddo

      if(numpe .gt. 1)then
         call drvAllreduce(xaveg, xavegr,12*nfath)         

         do m = 1, 12
            xavegr(:,m) = xavegr(:,m)/nsons(:)
         enddo

         if (myrank .eq. master) then
            write(*,*)'diss=', xavegt(14,11), xavegr(14,11)
         endif

         do m = 1, 12
            xavegt(:,m) = xavegt(:,m) + xavegr(:,m)
         enddo

      else

         do m = 1, 12
            xaveg(:,m) = xaveg(:,m)/nsons(:)
         enddo         
         
         do m = 1, 12
            xavegt(:,m) = xavegt(:,m) + xaveg(:,m)
         enddo         
         
      endif

      if (myrank .eq. master) then
         write(*,*)'diss=', xavegt(14,11), xavegr(14,11)
      endif

      if ( istep .eq. (nstep(1)-1) ) then
         if ( myrank .eq. master) then

            do i = 1, nfath
               write(376,*)xavegt(i,1),xavegt(i,2),xavegt(i,3)
               write(377,*)xavegt(i,4),xavegt(i,5),xavegt(i,6)
               write(378,*)xavegt(i,7),xavegt(i,8),xavegt(i,9)
               write(379,*)xavegt(i,10),xavegt(i,11),xavegt(i,12)
            enddo

            call flush(376)
            call flush(377)
            call flush(378)
            call flush(379)
            
         endif
      endif      


      return

      end


      subroutine resSij (y, x, ien, hfres, shgl, shp, Qwtf)

      include "common.h"

      dimension y(nshg,5),             hfres(nshg,24)
      dimension x(numnp,3),            xl(npro,nenl,3)
      dimension ien(npro,nshl),        yl(npro,nshl,5),
     &          fresl(npro,nshl,24),        WdetJ(npro),
     &          u1(npro),              u2(npro),
     &          u3(npro),              rho(npro),
     &          S11(npro),             S22(npro),
     &          S33(npro),             S12(npro),
     &          S13(npro),             S23(npro),
     &          dxdxi(npro,nsd,nsd),   dxidx(npro,nsd,nsd),
     &          shgl(nsd,nshl,ngauss),       shg(npro,nshl,nsd),
     &          shp(nshl,ngauss),          
     &          fresli(npro,nshl,24),   Qwtf(ngaussf),
     &          strnrmi(npro)
    

      dimension tmp(npro)

      call local (y,      yl,     ien,    5,  'gather  ')
      call localx (x,      xl,     ien,    3,  'gather  ')
c

      fresl = zero

c
      do intp = 1, ngaussf   ! Loop over quadrature points

c  calculate the metrics
c
c
c.... --------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
         dxdxi = zero
c
         do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(1,n,intp)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(2,n,intp)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(3,n,intp)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(1,n,intp)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(2,n,intp)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(3,n,intp)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(1,n,intp)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(2,n,intp)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(3,n,intp)
         enddo
c     
c.... compute the inverse of deformation gradient
c
         dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)
     &        - dxdxi(:,3,2) * dxdxi(:,2,3)
         dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)
     &        - dxdxi(:,1,2) * dxdxi(:,3,3)
         dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)
     &        - dxdxi(:,1,3) * dxdxi(:,2,2)
         tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)
     &        + dxidx(:,1,2) * dxdxi(:,2,1)
     &        + dxidx(:,1,3) * dxdxi(:,3,1) )
         dxidx(:,1,1) = dxidx(:,1,1) * tmp
         dxidx(:,1,2) = dxidx(:,1,2) * tmp
         dxidx(:,1,3) = dxidx(:,1,3) * tmp
         dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)
     &        - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
         dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)
     &        - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
         dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)
     &        - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
         dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)
     &        - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
         dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)
     &        - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
         dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)
     &        - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c     
c        wght=Qwt(lcsyst,intp)  ! may be different now
         wght=Qwtf(intp)    
         WdetJ(:) = wght / tmp(:)
         

c... compute the gradient of shape functions at the quad. point.


      do n = 1,nshl 
        shg(:,n,1) = (shgl(1,n,intp) * dxidx(:,1,1)
     &              + shgl(2,n,intp) * dxidx(:,2,1)
     &              + shgl(3,n,intp) * dxidx(:,3,1))
        shg(:,n,2) = (shgl(1,n,intp) * dxidx(:,1,2)
     &              + shgl(2,n,intp) * dxidx(:,2,2)
     &              + shgl(3,n,intp) * dxidx(:,3,2))
        shg(:,n,3) = (shgl(1,n,intp) * dxidx(:,1,3)
     &              + shgl(2,n,intp) * dxidx(:,2,3)
     &              + shgl(3,n,intp) * dxidx(:,3,3))
      enddo


c... compute the velocities and the strain rate tensor at the quad. point


         u1  = zero
         u2  = zero
         u3  = zero
         S11 = zero
         S22 = zero
         S33 = zero
         S12 = zero
         S13 = zero
         S23 = zero
         do i=1,nshl  
            u1 = u1 + shp(i,intp)*yl(:,i,2)
            u2 = u2 + shp(i,intp)*yl(:,i,3)
            u3 = u3 + shp(i,intp)*yl(:,i,4)
 
            S11 = S11 + shg(:,i,1)*yl(:,i,2)
            S22 = S22 + shg(:,i,2)*yl(:,i,3)
            S33 = S33 + shg(:,i,3)*yl(:,i,4)
            
            S12 = S12 + shg(:,i,2)*yl(:,i,2)
     &                       +shg(:,i,1)*yl(:,i,3)
            S13 = S13 + shg(:,i,3)*yl(:,i,2)
     &                       +shg(:,i,1)*yl(:,i,4)
            S23 = S23 + shg(:,i,3)*yl(:,i,3)
     &                       +shg(:,i,2)*yl(:,i,4) 
         enddo
         S12 = pt5 * S12
         S13 = pt5 * S13
         S23 = pt5 * S23

c... Get the strain rate norm at the quad pts

         strnrmi = sqrt( two*(S11**2 + S22**2 + S33**2) 
     &        + four*(S12**2 + S13**2 + S23**2) )


c... Loop over element nodes and multiply u_{i} and S_{i,j} by the
c... hat kernel and the Jacobian over the current element evaluated at 
c... the current quad. point.

         do i = 1, nenl   ! Loop over element nodes

            fresli(:,i,1) = WdetJ * u1 * shp(i,intp) !G * u1 * WdetJ
            fresli(:,i,2) = WdetJ * u2 * shp(i,intp) !G * u2 * WdetJ
            fresli(:,i,3) = WdetJ * u3 * shp(i,intp) !G * u3 * WdetJ

            fresli(:,i,4) = WdetJ * u1 * u1 * shp(i,intp) ! G*u1*u1*WdetJ
            fresli(:,i,5) = WdetJ * u2 * u2 * shp(i,intp) ! G*u2*u2*WdetJ
            fresli(:,i,6) = WdetJ * u3 * u3 * shp(i,intp) ! G*u3*u3*WdetJ
            fresli(:,i,7) = WdetJ * u1 * u2 * shp(i,intp) ! G*u1*u2*WdetJ
            fresli(:,i,8) = WdetJ * u1 * u3 * shp(i,intp) ! G*u1*u3*WdetJ
            fresli(:,i,9) = WdetJ * u2 * u3 * shp(i,intp) ! G*u2*u3*WdetJ

            fresli(:,i,10) = S11 * shp(i,intp) * WdetJ ! G*S_{1,1}*WdetJ
            fresli(:,i,11) = S22 * shp(i,intp) * WdetJ ! G*S_{2,2}*WdetJ
            fresli(:,i,12) = S33 * shp(i,intp) * WdetJ ! G*S_{3,3}*WdetJ
            fresli(:,i,13) = S12 * shp(i,intp) * WdetJ ! G*S_{1,1}*WdetJ
            fresli(:,i,14) = S13 * shp(i,intp) * WdetJ ! G*S_{1,3}*WdetJ
            fresli(:,i,15) = S23 * shp(i,intp) * WdetJ ! G*S_{2,3}*WdetJ


            fresli(:,i,16) = strnrmi*strnrmi*strnrmi*shp(i,intp)*WdetJ

            fresli(:,i,22) = WdetJ*shp(i,intp) !Integral of filter kernel, G,
c                                               over the element 

c...   Get G*|S|*S_{i,j}*WdetJ
           
c            fresli(:,i,16) = S11 * strnrmi * shp(i,intp) * WdetJ 
            fresli(:,i,17) = S22 * strnrmi * shp(i,intp) * WdetJ            
            fresli(:,i,18) = S33 * strnrmi * shp(i,intp) * WdetJ
            fresli(:,i,19) = S12 * strnrmi * shp(i,intp) * WdetJ  
            fresli(:,i,20) = S13 * strnrmi * shp(i,intp) * WdetJ
            fresli(:,i,21) = S23 * strnrmi * shp(i,intp) * WdetJ  

            do j = 1, 22 ! Add contribution of each quad. point for each
c                          element node
               fresl(:,i,j) = fresl(:,i,j) + fresli(:,i,j)
            enddo

         enddo                  !end loop over element nodes

      enddo                     !end of loop over integration points


      call local (hfres, fresl, ien, 24, 'scatter ')

      return
      end



	subroutine sparseCG (rhsorig, trx, lhsM, row, col, iper,
     &     ilwork, iBC, BC)
c
c------------------------------------------------------------------------
c
c  This subroutine uses Conjugate Gradient,
c to solve the system of equations.
c
c Farzin Shakib,  Summer 1987.
c------------------------------------------------------------------------
c
	include "common.h"
c
	dimension rhsorig(nshg), trx(nshg)
c
	dimension d(nshg),     p(nshg),
     &            q(nshg),     ilwork(nlwork),
     &		  Dinv(nshg),  rhs(nshg),
     &            pp(nshg),    
     &            iBC(nshg),
     &            BC(nshg)
  

        integer   row(nshg*nnz),         col(nshg+1)
	integer   iBCdumb(nshg),         iper(nshg)

	real*8 BCdumb(nshg,ndofBC)

	real*8	lhsM(nnz*nshg)
c
	data nitercf / 100 /
c
c
	BCdumb  = one
	iBCdumb = 1
c
	rhs(:)=rhsorig(:)
c
c.... Calculate the inverse of the diagonal of the mass matrix
c
c       call CFDinv (Dinv, emass)
c
c.... Left precondition the mass matrix and the RHS
c
c	call CFPre (Dinv, emass, rhs)
c
c Initialize. We have a system Ax=b We have made A as
c well conditionedand as we're willing to go.  Since
c we used left preconditioning (on the old A and old b),
c we don't need to do any back-preconditioning later.
c
	rr = 0
	do n = 1, nshg
	   rr  = rr + rhs(n)*rhs(n)
	enddo
c
c  if parallel the above is only this processors part of the dot product.
c  get the total result
c
	   dotTot=zero
	   call drvAllreduce(rr,dotTot,1)
	   rr=dotTot
	rr0 = rr
c
	trx(:) = zero ! x_{0}=0
c                   ! r_{0}=b
c
c                   ! beta_{1}=0
	p(:) = rhs(:) ! p_{1}=r_{0}=b
c
c.... Iterate
c        
	do iter = 1, nitercf      ! 15  ! nitercf
c
c.... calculate alpha
c
	   pp=p   ! make a vector that we can copy masters onto slaves
		  ! and thereby make ready for an accurate Ap product

#ifdef HAVE_LESLIB
	   call commOut(pp, ilwork, 1,iper,iBCdumb,BCdumb)  !slaves= master

	   call fLesSparseApSclr(	col,	row,	lhsM,	
     &					pp,	q,	nshg,
     &                                  nnz)

	   call commIn(q, ilwork, 1,iper,iBC,BC) ! masters=masters+slaves
							 ! slaves zeroed
#else
           if(myrank.eq.0) write(*,*) 'need alt solver in filter.f'
           call error('filter   ','noleslib',1)
#endif
c	   call CFAp (p,  q) ! put Ap product in q

c	   if(nump>1) call commu (q, ilwork, 1, 'in ') 

c	   do j = 1,nshg
c	      i = iper(j)
c	      if (i .ne. j) then
c		 q(i) = q(i) + q(j)
c	      endif
c	   enddo

c	   do j = 1,nshg
c	      i = iper(j)
c	      if (i .ne. j) then
c		 q(j) = zero
c	      endif
c	   enddo

c     Need to zero off-processor slaves as well.

c      if (numpe.gt.1 .and. nsons(1).gt.1) then

c         numtask = ilwork(1)
c         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  

c         do itask = 1, numtask
            
c            iacc   = ilwork (itkbeg + 2)
c            numseg = ilwork (itkbeg + 4)
	    
c            if (iacc .eq. 0) then
c               do is = 1,numseg
c                  isgbeg = ilwork (itkbeg + 3 + 2*is)
c                  lenseg = ilwork (itkbeg + 4 + 2*is)
c                  isgend = isgbeg + lenseg - 1
c                  q(isgbeg:isgend) = zero
c               enddo
c           endif
            
c            itkbeg = itkbeg + 4 + 2*numseg
            
c         enddo

c	endif         
	 


	   pap = 0
	   do  n = 1, nshg
	      pap = pap + p(n) * q(n)
	   enddo
c
c  if parallel the above is only this processors part of the dot product.
c  get the total result
c
           dotTot=zero
	   call drvAllreduce(pap,dotTot,1)
	   pap=dotTot
	   alpha = rr / pap 
c
c.... calculate the next guess
c
	   trx(:) = trx(:) + alpha * p(:)
c
c.... calculate the new residual
c
c
	   rhs(:) = rhs(:) - alpha * q(:)
	   tmp = rr
	   rr = 0
	   do n = 1, nshg
	      rr = rr + rhs(n)*rhs(n)
	   enddo
c
c  if parallel the above is only this processors part of the dot product.
c  get the total result
c
           dotTot=zero
	   call drvAllreduce(rr,dotTot,1)
	   rr=dotTot
c
c.... check for convergence
c
	   if(rr.lt.100.*epsM**2) goto 6000
c
c.... calculate a new search direction
c
	   beta = rr / tmp
	   p(:) = rhs(:) + beta * p(:)
c
c.... end of iteration
c
	enddo
c
c.... if converged
c
6000	continue

c need a commu(out) on solution (TRX) to get slaves the correct solution AND
c on processor slaves = on processor masters

	if(numpe>1) call commu (trx, ilwork, 1, 'out ')
	trx(:) = trx(iper(:))

c
	write(*,9000) iter, rr / rr0
c	write(16,9000) iter, rr / rr0
c
c.... return
c
	return
9000	format(20x,'  number of iterations:', i10,/,
     &	       20x,'    residual reduction:', 2x,e10.2)
	end

      subroutine stdfdmc (y,      shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x)

      use pointer_data

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwt contains the weights of the 
c                    quad. points.  



      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension fres(nshg,24),         fwr(nshg),
     &          strnrm(nshg),         cdelsq(nshg),
     &          cdelsq2(nshg),
     &          xnum(nshg),           xden(nshg),
     &          xmij(nshg,6),         xlij(nshg,6),
     &          xnude(nfath,2),        xnuder(nfath,2),
     &          ynude(nfath,6),        ynuder(nfath,6)
      dimension ui(nfath,3),           snorm(nfath),
     &          uir(nfath,3),          snormr(nfath),
     &          xm(nfath,6),           xl(nfath,6),
     &          xl1(nfath,6),          xl2(nfath,6),
     &          xl1r(nfath,6),         xl2r(nfath,6),
     &          xmr(nfath,6),          xlr(nfath,6),
     &          nsons(nshg),
     &          strl(numel,ngauss)           
      dimension y(nshg,5),            yold(nshg,5),
     &          ifath(nshg),          iper(nshg),
     &          ilwork(nlwork),!        xmudmi(numel,ngauss),
     &          x(numnp,3),
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),    
     &          xnutf(nfath),         xfac(nshg,5)

      character*10 cname
      character*30 fname1, fname2, fname3, fname4, fname5, fname6,
     &             fname0
c
c
c   setup the weights for time averaging of cdelsq (now in quadfilt module)
c
      denom=max(1.0d0*(lstep),one)
      if(dtavei.lt.0) then
         wcur=one/denom
      else
         wcur=dtavei
      endif  
      whist=1.0-wcur

      if (istep .eq. 0) then
         xnd      = zero
         xmodcomp = zero
         xmcomp  = zero
         xlcomp  = zero
         xl1comp  = zero
         xl2comp  = zero
         ucomp    = zero
         scomp    = zero
      endif

     
      fres = zero
      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)
c

c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      yold(:,5) = 1.0  ! Debugging
c      yold(:,2) = 2.0*x(:,1) - 3.0*x(:,2) 
c      yold(:,2) = 2.0
c      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      yold(:,3) = 3.0
c      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
c      yold(:,4) = 4.0
c      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
c                               suitable for the



      intrul=intg(1,itseq)
      intind=intpt(intrul)

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst) 
        
        call asithf (yold, x, strl(iel:inum,:), mien(iblk)%p, fres, 
     &               shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo
c
 
c      if (ngaussf .ne. ngauss) then
      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)

        if (ngaussf .ne. ngauss) then

        call getstrl (yold, x,      mien(iblk)%p,  
     &               strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:))

        endif

      enddo
c      endif
c
c
C must fix for abc and dynamic model
c      if(iabc==1)   !are there any axisym bc's
c     &      call rotabc(res, iBC, BC,nflow, 'in ')
c
      if(numpe>1) call commu (fres, ilwork, 24, 'in ')
c 
c account for periodicity in filtered variables
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:)
        endif
      enddo
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = fres(i,:)
        endif
      enddo

      if(numpe>1)   call commu (fres, ilwork, 24, 'out')

      fres(:,23) = one / fres(:,23)
      do j = 1,22
        fres(:,j) = fres(:,j) * fres(:,23)
      enddo
c     fres(:,24) = fres(:,24) * fres(:,23)
c
c.....at this point fres is really all of our filtered quantities
c     at the nodes
c

      strnrm = sqrt( 
     &  two * (fres(:,10)**2 + fres(:,11)**2 + fres(:,12)**2)
     &  + four * ( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )

      fwr = fwr1 * fres(:,22) * strnrm

      xmij(:,1) = -fwr
     &             * fres(:,10) + fres(:,16)
      xmij(:,2) = -fwr
     &             * fres(:,11) + fres(:,17) 
      xmij(:,3) = -fwr
     &             * fres(:,12) + fres(:,18) 

      xmij(:,4) = -fwr * fres(:,13) + fres(:,19)
      xmij(:,5) = -fwr * fres(:,14) + fres(:,20)
      xmij(:,6) = -fwr * fres(:,15) + fres(:,21)

      fres(:,22) = one / fres(:,22)

      xlij(:,1) = fres(:,4) - fres(:,1) * fres(:,1) * fres(:,22)
      xlij(:,2) = fres(:,5) - fres(:,2) * fres(:,2) * fres(:,22)
      xlij(:,3) = fres(:,6) - fres(:,3) * fres(:,3) * fres(:,22)
      xlij(:,4) = fres(:,7) - fres(:,1) * fres(:,2) * fres(:,22)
      xlij(:,5) = fres(:,8) - fres(:,1) * fres(:,3) * fres(:,22)
      xlij(:,6) = fres(:,9) - fres(:,2) * fres(:,3) * fres(:,22)

      xnum =        xlij(:,1) * xmij(:,1) + xlij(:,2) * xmij(:,2) 
     &                                    + xlij(:,3) * xmij(:,3)
     &     + two * (xlij(:,4) * xmij(:,4) + xlij(:,5) * xmij(:,5)
     &                                    + xlij(:,6) * xmij(:,6))
      xden =        xmij(:,1) * xmij(:,1) + xmij(:,2) * xmij(:,2) 
     &                                    + xmij(:,3) * xmij(:,3)
     &     + two * (xmij(:,4) * xmij(:,4) + xmij(:,5) * xmij(:,5)
     &                                    + xmij(:,6) * xmij(:,6))
      xden = two * xden

c... For collectection of statistics on dyn. model components 

      xfac(:,1) = strnrm*strnrm*( fres(:,10)**2 + fres(:,11)**2 + 
     &     fres(:,12)**2
     &     + two*( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )
      
      xfac(:,2) = strnrm*( xlij(:,1)*fres(:,10) + xlij(:,2)*fres(:,11) 
     &     + xlij(:,3)*fres(:,12) + 
     &     two*(xlij(:,4)*fres(:,13) + xlij(:,5)*fres(:,14) +
     &     xlij(:,6)*fres(:,15)) )

      xfac(:,3) = strnrm*( fres(:,10)*fres(:,16) + fres(:,11)*fres(:,17)
     &     + fres(:,12)*fres(:,18) + 
     &     two*(fres(:,13)*fres(:,19) + fres(:,14)*fres(:,20) +
     &     fres(:,15)*fres(:,21)) )

      xfac(:,4) = xlij(:,1)*fres(:,16) + xlij(:,2)*fres(:,17)
     &     + xlij(:,3)*fres(:,18) + 
     &     two*(xlij(:,4)*fres(:,19) + xlij(:,5)*fres(:,20) +
     &     xlij(:,6)*fres(:,21))

      xfac(:,5) = fres(:,16)*fres(:,16) + fres(:,17)*fres(:,17)
     &     + fres(:,18)*fres(:,18) + 
     &     two*(fres(:,19)*fres(:,19) + fres(:,20)*fres(:,20) +
     &     fres(:,21)*fres(:,21))

c  zero on processor periodic nodes so that they will not be added twice
        do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xnum(j) = zero
            xden(j) = zero
            xfac(j,:) = zero
            xmij(j,:) = zero
            xlij(j,:) = zero
            fres(j,:) = zero
            strnrm(j) = zero
          endif
        enddo

      if (numpe.gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xnum(isgbeg:isgend) = zero
                  xden(isgbeg:isgend) = zero
                  strnrm(isgbeg:isgend) = zero
                  xfac(isgbeg:isgend,:) = zero
                  xmij(isgbeg:isgend,:) = zero
                  xlij(isgbeg:isgend,:) = zero
                  fres(isgbeg:isgend,:) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
c
c Description of arrays.   Each processor has an array of length equal
c to the total number of fathers times 2 xnude(nfathers,2). One to collect 
c the numerator and one to collect the denominator.  There is also an array
c of length nshg on each processor which tells the father number of each
c on processor node, ifath(nnshg).  Finally, there is an arry of length
c nfathers to tell the total (on all processors combined) number of sons
c for each father. 
c
c  Now loop over nodes and accumlate the numerator and the denominator
c  to the father nodes.  Only on processor addition at this point.
c  Note that serrogate fathers are collect some for the case where some
c  sons are on another processor
c
      xnude = zero
      ynude = zero
      xm    = zero
      xl    = zero
      xl1   = zero
      xl2   = zero
      ui    = zero
      snorm = zero

      do i = 1,nshg
         xnude(ifath(i),1) = xnude(ifath(i),1) + xnum(i)
         xnude(ifath(i),2) = xnude(ifath(i),2) + xden(i)

         ynude(ifath(i),1) = ynude(ifath(i),1) + xfac(i,1)
         ynude(ifath(i),2) = ynude(ifath(i),2) + xfac(i,2)
         ynude(ifath(i),3) = ynude(ifath(i),3) + xfac(i,3)
         ynude(ifath(i),4) = ynude(ifath(i),4) + xfac(i,4)
         ynude(ifath(i),5) = ynude(ifath(i),5) + xfac(i,5)

         xm(ifath(i),1) = xm(ifath(i),1) + xmij(i,1)
         xm(ifath(i),2) = xm(ifath(i),2) + xmij(i,2)
         xm(ifath(i),3) = xm(ifath(i),3) + xmij(i,3)
         xm(ifath(i),4) = xm(ifath(i),4) + xmij(i,4)
         xm(ifath(i),5) = xm(ifath(i),5) + xmij(i,5)
         xm(ifath(i),6) = xm(ifath(i),6) + xmij(i,6)

         xl(ifath(i),1) = xl(ifath(i),1) + xlij(i,1)
         xl(ifath(i),2) = xl(ifath(i),2) + xlij(i,2)
         xl(ifath(i),3) = xl(ifath(i),3) + xlij(i,3)
         xl(ifath(i),4) = xl(ifath(i),4) + xlij(i,4)
         xl(ifath(i),5) = xl(ifath(i),5) + xlij(i,5)
         xl(ifath(i),6) = xl(ifath(i),6) + xlij(i,6)         

         xl1(ifath(i),1) = xl1(ifath(i),1) + fres(i,4)
         xl1(ifath(i),2) = xl1(ifath(i),2) + fres(i,5)
         xl1(ifath(i),3) = xl1(ifath(i),3) + fres(i,6)
         xl1(ifath(i),4) = xl1(ifath(i),4) + fres(i,7)
         xl1(ifath(i),5) = xl1(ifath(i),5) + fres(i,8)
         xl1(ifath(i),6) = xl1(ifath(i),6) + fres(i,9)         

         xl2(ifath(i),1) = xl2(ifath(i),1) + fres(i,1)*fres(i,1) 
         xl2(ifath(i),2) = xl2(ifath(i),2) + fres(i,2)*fres(i,2)
         xl2(ifath(i),3) = xl2(ifath(i),3) + fres(i,3)*fres(i,3) 
         xl2(ifath(i),4) = xl2(ifath(i),4) + fres(i,1)*fres(i,2)
         xl2(ifath(i),5) = xl2(ifath(i),5) + fres(i,1)*fres(i,3)
         xl2(ifath(i),6) = xl2(ifath(i),6) + fres(i,2)*fres(i,3)

         ui(ifath(i),1) = ui(ifath(i),1) + fres(i,1)
         ui(ifath(i),2) = ui(ifath(i),2) + fres(i,2)
         ui(ifath(i),3) = ui(ifath(i),3) + fres(i,3)

         snorm(ifath(i)) = snorm(ifath(i)) + strnrm(i)

      enddo

c
c Now  the true fathers and serrogates combine results and update
c each other.
c       
      if(numpe .gt. 1)then
         call drvAllreduce(xnude, xnuder,2*nfath)
         call drvAllreduce(ynude, ynuder,6*nfath)
         call drvAllreduce(xm, xmr,6*nfath)
         call drvAllreduce(xl, xlr,6*nfath)
         call drvAllreduce(xl1, xl1r,6*nfath)
         call drvAllreduce(xl2, xl2r,6*nfath)
         call drvAllreduce(ui, uir,3*nfath)
         call drvAllreduce(snorm, snormr,nfath)

         do i = 1, nfath
            ynuder(i,6) = ( ynuder(i,4) - fwr1*ynuder(i,2) ) /
     &           ( two*ynuder(i,5) - four*fwr1*ynuder(i,3)
     &           + two*fwr1*fwr1*ynuder(i,1) )
         enddo

         cdelsq2(:) = ynuder(ifath(:),6)  ! For comparison w/ cdelsq
c
c  xnude is the sum of the sons for each father on this processor
c
c  xnuder is the sum of the sons for each father on all processor combined
c  (the same as if we had not partitioned the mesh for each processor)
c
c   For each father we have precomputed the number of sons (including
c   the sons off processor). 
c
c   Now divide by number of sons to get the average (not really necessary
c   for dynamic model since ratio will cancel nsons at each father)
c
         xnuder(:,1) = xnuder(:,1) / nsons(:)
         xnuder(:,2) = xnuder(:,2) / nsons(:)

         do m = 1, 5
         ynuder(:,m) = ynuder(:,m)/nsons(:)
         enddo
         do m = 1,6
         xmr(:,m) = xmr(:,m)/nsons(:)
         xlr(:,m) = xlr(:,m)/nsons(:)
         xl1r(:,m) = xl1r(:,m)/nsons(:)
         xl2r(:,m) = xl2r(:,m)/nsons(:)
         enddo

         uir(:,1) = uir(:,1)/nsons(:)
         uir(:,2) = uir(:,2)/nsons(:)
         uir(:,3) = uir(:,3)/nsons(:)

         snormr(:) = snormr(:)/nsons(:)
c
cc  the next line is c \Delta^2
cc
cc         xnuder(:,1) = xnuder(:,1) / (xnuder(:,2) + 1.d-09)
cc         do i = 1,nshg
cc            cdelsq(i) = xnuder(ifath(i),1)
cc         enddo

            numNden(:,1) = whist*numNden(:,1)+wcur*xnuder(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnuder(ifath(:),2)
            cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)
            
c            cdelsq(:) = xnuder(ifath(:),1)/(xnuder(ifath(:),2)+1.d-09)

            xnd(:,1) = xnd(:,1) + xnuder(:,1)
            xnd(:,2) = xnd(:,2) + xnuder(:,2)

            xmodcomp(:,1) = xmodcomp(:,1)+ynuder(:,1)
            xmodcomp(:,2) = xmodcomp(:,2)+ynuder(:,2)            
            xmodcomp(:,3) = xmodcomp(:,3)+ynuder(:,3)
            xmodcomp(:,4) = xmodcomp(:,4)+ynuder(:,4)
            xmodcomp(:,5) = xmodcomp(:,5)+ynuder(:,5)

            xmcomp(:,:) = xmcomp(:,:)+xmr(:,:)
            xlcomp(:,:) = xlcomp(:,:)+xlr(:,:)

            xl1comp(:,:) = xl1comp(:,:)+xl1r(:,:)
            xl2comp(:,:) = xl2comp(:,:)+xl2r(:,:)

            ucomp(:,:) = ucomp(:,:)+uir(:,:)
            u1 = uir(32,1)
            scomp(:)   = scomp(:)+snormr(:)

      else

         xnude(:,1) = xnude(:,1)/nsons(:)
         xnude(:,2) = xnude(:,2)/nsons(:)

         do m = 1, 5
         ynude(:,m) = ynude(:,m)/nsons(:)
         enddo
         do m = 1,6
         xm(:,m) = xm(:,m)/nsons(:)
         xl(:,m) = xl(:,m)/nsons(:)
         xl1(:,m) = xl1(:,m)/nsons(:)
         xl2(:,m) = xl2(:,m)/nsons(:)
         enddo

         ui(:,1) = ui(:,1)/nsons(:)
         ui(:,2) = ui(:,2)/nsons(:)
         ui(:,3) = ui(:,3)/nsons(:)

         snorm(:) = snorm(:)/nsons(:)

c     
c     the next line is c \Delta^2, not nu_T but we want to save the
c     memory
c     

cc         xnude(:,1) = xnude(:,1) / (xnude(:,2) + 1.d-09)
cc        do i = 1,nshg
cc            cdelsq(i) = xnude(ifath(i),1)
cc         enddo
cc      endif

         do i = 1, nfath
            ynude(i,6) = ( ynude(i,4) - fwr1*ynude(i,2) ) /
     &           ( two*ynude(i,5) - four*fwr1*ynude(i,3)
     &           + fwr1*fwr1*ynude(i,1) )
         enddo

            numNden(:,1) = whist*numNden(:,1)+wcur*xnude(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnude(ifath(:),2)

            xnd(:,1) = xnd(:,1)+xnude(:,1)
            xnd(:,2) = xnd(:,2)+xnude(:,2)

            cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)

c            cdelsq(:) = xnude(ifath(:),1)/(xnude(ifath(:),2))!+1.d-09)
            

          cdelsq2(:) = ynude(ifath(:),6)  ! For comparison w/ cdelsq

            xmodcomp(:,1) = xmodcomp(:,1)+ynude(:,1)
            xmodcomp(:,2) = xmodcomp(:,2)+ynude(:,2)            
            xmodcomp(:,3) = xmodcomp(:,3)+ynude(:,3)
            xmodcomp(:,4) = xmodcomp(:,4)+ynude(:,4)
            xmodcomp(:,5) = xmodcomp(:,5)+ynude(:,5)

            xmcomp(:,:) = xmcomp(:,:)+xm(:,:)
            xlcomp(:,:) = xlcomp(:,:)+xl(:,:)

            xl1comp(:,:) = xl1comp(:,:)+xl1(:,:)
            xl2comp(:,:) = xl2comp(:,:)+xl2(:,:)

            ucomp(:,:) = ucomp(:,:)+ui(:,:)
            scomp(:)   = scomp(:)+snorm(:)

         endif

c         do i = 1, nfath
c            xmodcomp(i,:) = xmodcomp(i,:)/nsons(i)
c            xmcomp(i,:) = xmcomp(i,:)/nsons(i)         
c            xlcomp(i,:) = xlcomp(i,:)/nsons(i)
c            xl2comp(i,:) = xl2comp(i,:)/nsons(i)         
c            xl1comp(i,:) = xl1comp(i,:)/nsons(i)
c            xnd(i,:) = xnd(i,:)/nsons(i)
c            scomp(i) = scomp(i)/nsons(i)
c            ucomp(i,:) = ucomp(i,:)/nsons(i)
c         enddo

         if (myrank .eq. master) then
            write(*,*)'istep, nstep=', istep, nstep(1)
         endif

         if ( istep .eq. (nstep(1)-1) ) then
         if ( myrank .eq. master) then

            do i = 1, nfath
            write(365,*)xmodcomp(i,1),xmodcomp(i,2),xmodcomp(i,3),
     &              xmodcomp(i,4),xmodcomp(i,5)

            write(366,*)xmcomp(i,1),xmcomp(i,2),xmcomp(i,3)
            write(367,*)xmcomp(i,4),xmcomp(i,5),xmcomp(i,6)            

            write(368,*)xlcomp(i,1),xlcomp(i,2),xlcomp(i,3)
            write(369,*)xlcomp(i,4),xlcomp(i,5),xlcomp(i,6)

            write(370,*)xl1comp(i,1),xl1comp(i,2),xl1comp(i,3)
            write(371,*)xl1comp(i,4),xl1comp(i,5),xl1comp(i,6) 

            write(372,*)xl2comp(i,1),xl2comp(i,2),xl2comp(i,3)
            write(373,*)xl2comp(i,4),xl2comp(i,5),xl2comp(i,6)

            write(374,*)xnd(i,1),xnd(i,2),scomp(i)
            write(375,*)ucomp(i,1),ucomp(i,2),ucomp(i,3) 
            enddo

            call flush(365)
            call flush(366)
            call flush(367)
            call flush(368)
            call flush(369)
            call flush(370)
            call flush(371)
            call flush(372)
            call flush(373)
            call flush(374)
            call flush(375)

c            close(852)
c            close(853)
c            close(854)

         endif
         endif

            if (myrank .eq. master) then
               write(*,*)'uit uic=', ucomp(32,1),u1
            endif

 555     format(e14.7,4(2x,e14.7))
 556     format(e14.7,5(2x,e14.7))

         


c $$$$$$$$$$$$$$$$$$$$$$$$$$$
      tmp1 =  MINVAL(cdelsq)
      tmp2 =  MAXVAL(cdelsq)
      if(numpe>1) then
         call MPI_REDUCE (tmp1, tmp3, 1,MPI_DOUBLE_PRECISION,
     &        MPI_MIN, master, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
     &        MPI_MAX, master, MPI_COMM_WORLD, ierr)
         tmp1=tmp3
         tmp2=tmp4
      endif
      if (myrank .EQ. master) then !print CDelta^2 range
         write(34,*)lstep,tmp1,tmp2
         call flush(34)
      endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$
      
      if (myrank .eq. master) then
         write(*,*) 'cdelsq=', cdelsq(1),cdelsq(2)
         write(*,*) 'cdelsq=', cdelsq2(1),cdelsq2(2)
         write(22,*) lstep, cdelsq(1)
         call flush(22)
      endif

      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         
         ngauss = nint(lcsyst)

         call scatnu (mien(iblk)%p, strl(iel:inum,:), 
     &        mxmudmi(iblk)%p,cdelsq,shp(lcsyst,1:nshl,:))
      enddo
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$  tmp1 =  MINVAL(xmudmi)
c$$$  tmp2 =  MAXVAL(xmudmi)
c$$$  if(numpe>1) then
c$$$  call MPI_REDUCE (tmp1, tmp3, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MIN, master, MPI_COMM_WORLD, ierr)
c$$$  call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MAX, master, MPI_COMM_WORLD, ierr)
c$$$      tmp1=tmp3
c$$$  tmp2=tmp4
c$$$  endif
c$$$  if (myrank .EQ. master) then
c$$$  write(35,*) lstep,tmp1,tmp2
c$$$  call flush(35)
c$$$  endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$

c
c  if flag set, write a restart file with info (reuse xmij's memory)
c
      if(irs.eq.11) then
         lstep=999
         xmij(:,1)=xnum(:)
         xmij(:,2)=xden(:)
         xmij(:,3)=cdelsq(:)
         xmij(:,5)=xlij(:,4)    !leave M_{12} in 4 and put L_{12} here
         call restar('out ',xmij,xlij) !also dump all of L_{ij} in ac
         stop
      endif
c
c  local clipping moved to scatnu with the creation of mxmudmi pointers
c
c$$$      rmu=datmat(1,2,1)
c$$$      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
c$$$      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0
c      stop !uncomment to test dmod
c


c  write out the nodal values of xnut (estimate since we don't calc strain
c  there and must use the filtered strain).
c

      if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
c
c  collect the average strain into xnude(2)
c
         xnude(:,2) = zero
         do i = 1,numnp
            xnude(ifath(i),2) = xnude(ifath(i),2) + strnrm(i)
         enddo

         if(numpe .gt. 1) then
             call drvAllreduce(xnude(:,2), xnuder(:,2),nfath)
          else
             xnuder=xnude
          endif
c     
c          nut= cdelsq    * |S|
c 
         xnutf=xnuder(:,1)*xnuder(:,2)/nsons(:)
c
c  collect the x and y coords into xnude
c
         xnude = zero
         do i = 1,numnp
            xnude(ifath(i),1) = xnude(ifath(i),1) + x(i,1)
            xnude(ifath(i),2) = xnude(ifath(i),2) + x(i,2)
         enddo

         if(numpe .gt. 1) 
     &        call drvAllreduce(xnude, xnuder,2*nfath)
         xnuder(:,1)=xnuder(:,1)/nsons(:)
         xnuder(:,2)=xnuder(:,2)/nsons(:)
c
c  xnude is the sum of the sons for each father on this processor
c
         if((myrank.eq.master)) then
            do i=1,nfath      ! cdelsq   * |S|
               write(444,*) xnuder(i,1),xnuder(i,2),xnutf(i)
            enddo
            call flush(444)
         endif
      endif

      return
      end
      subroutine widefdmc (y,      shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x)

      use pointer_data

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwtf contains the weights of the 
c                    quad. points.  

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension fres(nshg,33),         fwr(nshg),
     &          strnrm(nshg),         cdelsq(nshg),
     &          cdelsq2(nshg),
     &          xnum(nshg),           xden(nshg),
     &          xmij(nshg,6),         xlij(nshg,6),
     &          xnude(nfath,2),        xnuder(nfath,2),
     &          ynude(nfath,6),        ynuder(nfath,6),
     &          ui(nfath,3),           snorm(nfath),
     &          uir(nfath,3),          snormr(nfath)
      dimension xm(nfath,6),           xl(nfath,6),
     &          xl1(nfath,6),          xl2(nfath,6),
     &          xl1r(nfath,6),         xl2r(nfath,6),
     &          xmr(nfath,6),          xlr(nfath,6),
     &          nsons(nshg),
     &          strl(numel,ngauss),           
     &          y(nshg,5),            yold(nshg,5),
     &          ifath(nshg),          iper(nshg),
     &          ilwork(nlwork),
     &          x(numnp,3)
      dimension shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),    
     &          xnutf(nfath),
     &          hfres(nshg,22),
     &          xfac(nshg,5)

      character*10 cname
      character*30 fname1, fname2, fname3, fname4, fname5, fname6
c
     
c
c
c   setup the weights for time averaging of cdelsq (now in quadfilt module)
c

      denom=max(1.0d0*(lstep),one)
      if(dtavei.lt.0) then
         wcur=one/denom
      else
         wcur=dtavei
      endif  
      whist=1.0-wcur

      if (myrank .eq. master) then
         write(*,*)'istep=', istep
      endif

      if (istep .eq. 0) then
         xnd      = zero
         xmodcomp = zero
         xmcomp  = zero
         xlcomp  = zero
         xl1comp  = zero
         xl2comp  = zero
         ucomp    = zero
         scomp    = zero
      endif


      fres = zero
      hfres = zero

      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)

c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      yold(:,5) = 1.0  ! Debugging
c      yold(:,2) = 2.0*x(:,1) - 3.0*x(:,2) 
c      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
c      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
c                               suitable for the

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
c        call hfilterBB (yold, x, mien(iblk)%p, hfres, 
c     &               shglf(lcsyst,:,1:nshl,:),
c     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

        call hfilterCC (yold, x, mien(iblk)%p, hfres, 
     &               shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo

      if(numpe>1) call commu (hfres, ilwork, 22, 'in ')
c 
c... account for periodicity in filtered variables
c
      do j = 1,nshg  !    Add on-processor slave contribution to masters
        i = iper(j)
        if (i .ne. j) then
           hfres(i,:) = hfres(i,:) + hfres(j,:)
        endif
      enddo
      do j = 1,nshg ! Set on-processor slaves to be the same as masters
        i = iper(j)
        if (i .ne. j) then
           hfres(j,:) = hfres(i,:)
        endif
      enddo

c... Set off-processor slaves to be the same as their masters

      if(numpe>1)   call commu (hfres, ilwork, 22, 'out')


      hfres(:,16) = one / hfres(:,16) ! one/(volume filter kernel)

      do j = 1, 15
	hfres(:,j) = hfres(:,j) * hfres(:,16)
      enddo	    		
      do j = 17, 22
	hfres(:,j) = hfres(:,j) * hfres(:,16)
      enddo	

c... For debugging

c      hfres(:,1) = 2.0*x(:,1) - 3.0*x(:,2) 
c      hfres(:,2) = 3.0*x(:,1) + 4.0*x(:,2)
c      hfres(:,3) = 4.0*x(:,1) + x(:,2) + x(:,3)

c... Done w/ h-filtering. Begin 2h-filtering.

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        call twohfilterBB (yold, x, strl(iel:inum,:), mien(iblk)%p, 
     &               fres, hfres, shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo
c
 

      if(numpe>1) call commu (fres, ilwork, 33, 'in ')
c 
c account for periodicity in filtered variables
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:)
        endif
      enddo

      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = fres(i,:)
        endif
      enddo

      if(numpe>1)then
         call commu (fres, ilwork, 33, 'out')
      endif

      fres(:,22) = one / fres(:,22)
      do j = 1,21
        fres(:,j) = fres(:,j) * fres(:,22)
      enddo
      do j = 23,33
        fres(:,j) = fres(:,j) * fres(:,22)
      enddo

      
      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)
 
        call getstrl (yold, x,      mien(iblk)%p,  
     &               strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:))

      enddo

c
c... Obtain the hat-tilde strain rate norm at the nodes 
c

      strnrm = sqrt( 
     &  two * (fres(:,10)**2 + fres(:,11)**2 + fres(:,12)**2)
     &  + four * ( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )

      fwr = fwr1 * strnrm

      xmij(:,1) = -fwr
     &             * fres(:,10) + fres(:,16)
      xmij(:,2) = -fwr
     &             * fres(:,11) + fres(:,17) 
      xmij(:,3) = -fwr
     &             * fres(:,12) + fres(:,18) 

      xmij(:,4) = -fwr * fres(:,13) + fres(:,19)
      xmij(:,5) = -fwr * fres(:,14) + fres(:,20)
      xmij(:,6) = -fwr * fres(:,15) + fres(:,21)


      xlij(:,1) = fres(:,4) - fres(:,1) * fres(:,1) 
      xlij(:,2) = fres(:,5) - fres(:,2) * fres(:,2) 
      xlij(:,3) = fres(:,6) - fres(:,3) * fres(:,3) 
      xlij(:,4) = fres(:,7) - fres(:,1) * fres(:,2) 
      xlij(:,5) = fres(:,8) - fres(:,1) * fres(:,3) 
      xlij(:,6) = fres(:,9) - fres(:,2) * fres(:,3) 

      xnum =        xlij(:,1) * xmij(:,1) + xlij(:,2) * xmij(:,2) 
     &                                    + xlij(:,3) * xmij(:,3)
     &     + two * (xlij(:,4) * xmij(:,4) + xlij(:,5) * xmij(:,5)
     &                                    + xlij(:,6) * xmij(:,6))
      xden =        xmij(:,1) * xmij(:,1) + xmij(:,2) * xmij(:,2) 
     &                                    + xmij(:,3) * xmij(:,3)
     &     + two * (xmij(:,4) * xmij(:,4) + xmij(:,5) * xmij(:,5)
     &                                    + xmij(:,6) * xmij(:,6))
      xden = two * xden

c... For collectection of statistics on dyn. model components 

      xfac(:,1) = strnrm*strnrm*( fres(:,10)**2 + fres(:,11)**2 + 
     &     fres(:,12)**2
     &     + two*( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )
      
      xfac(:,2) = strnrm*( xlij(:,1)*fres(:,10) + xlij(:,2)*fres(:,11) 
     &     + xlij(:,3)*fres(:,12) + 
     &     two*(xlij(:,4)*fres(:,13) + xlij(:,5)*fres(:,14) +
     &     xlij(:,6)*fres(:,15)) )

      xfac(:,3) = strnrm*( fres(:,10)*fres(:,16) + fres(:,11)*fres(:,17)
     &     + fres(:,12)*fres(:,18) + 
     &     two*(fres(:,13)*fres(:,19) + fres(:,14)*fres(:,20) +
     &     fres(:,15)*fres(:,21)) )

      xfac(:,4) = xlij(:,1)*fres(:,16) + xlij(:,2)*fres(:,17)
     &     + xlij(:,3)*fres(:,18) + 
     &     two*(xlij(:,4)*fres(:,19) + xlij(:,5)*fres(:,20) +
     &     xlij(:,6)*fres(:,21))

      xfac(:,5) = fres(:,16)*fres(:,16) + fres(:,17)*fres(:,17)
     &     + fres(:,18)*fres(:,18) + 
     &     two*(fres(:,19)*fres(:,19) + fres(:,20)*fres(:,20) +
     &     fres(:,21)*fres(:,21))

c  zero on processor periodic nodes so that they will not be added twice
        do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xnum(j) = zero
            xden(j) = zero
            xfac(j,:) = zero
            xmij(j,:) = zero
            xlij(j,:) = zero
            fres(j,:) = zero
            strnrm(j) = zero
          endif
        enddo

      if (numpe.gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xnum(isgbeg:isgend) = zero
                  xden(isgbeg:isgend) = zero
                  strnrm(isgbeg:isgend) = zero
                  xfac(isgbeg:isgend,:) = zero
                  xmij(isgbeg:isgend,:) = zero
                  xlij(isgbeg:isgend,:) = zero
                  fres(isgbeg:isgend,:) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
c
c Description of arrays.   Each processor has an array of length equal
c to the total number of fathers times 2 xnude(nfathers,2). One to collect 
c the numerator and one to collect the denominator.  There is also an array
c of length nshg on each processor which tells the father number of each
c on processor node, ifath(nnshg).  Finally, there is an arry of length
c nfathers to tell the total (on all processors combined) number of sons
c for each father. 
c
c  Now loop over nodes and accumlate the numerator and the denominator
c  to the father nodes.  Only on processor addition at this point.
c  Note that serrogate fathers are collect some for the case where some
c  sons are on another processor
c
      xnude = zero
      ynude = zero
      xm    = zero
      xl    = zero
      xl1   = zero
      xl2   = zero
      ui    = zero
      snorm = zero

      do i = 1,nshg
         xnude(ifath(i),1) = xnude(ifath(i),1) + xnum(i)
         xnude(ifath(i),2) = xnude(ifath(i),2) + xden(i)

         ynude(ifath(i),1) = ynude(ifath(i),1) + xfac(i,1)
         ynude(ifath(i),2) = ynude(ifath(i),2) + xfac(i,2)
         ynude(ifath(i),3) = ynude(ifath(i),3) + xfac(i,3)
         ynude(ifath(i),4) = ynude(ifath(i),4) + xfac(i,4)
         ynude(ifath(i),5) = ynude(ifath(i),5) + xfac(i,5)

         xm(ifath(i),1) = xm(ifath(i),1) + xmij(i,1)
         xm(ifath(i),2) = xm(ifath(i),2) + xmij(i,2)
         xm(ifath(i),3) = xm(ifath(i),3) + xmij(i,3)
         xm(ifath(i),4) = xm(ifath(i),4) + xmij(i,4)
         xm(ifath(i),5) = xm(ifath(i),5) + xmij(i,5)
         xm(ifath(i),6) = xm(ifath(i),6) + xmij(i,6)

         xl(ifath(i),1) = xl(ifath(i),1) + xlij(i,1)
         xl(ifath(i),2) = xl(ifath(i),2) + xlij(i,2)
         xl(ifath(i),3) = xl(ifath(i),3) + xlij(i,3)
         xl(ifath(i),4) = xl(ifath(i),4) + xlij(i,4)
         xl(ifath(i),5) = xl(ifath(i),5) + xlij(i,5)
         xl(ifath(i),6) = xl(ifath(i),6) + xlij(i,6)         

         xl1(ifath(i),1) = xl1(ifath(i),1) + fres(i,4)
         xl1(ifath(i),2) = xl1(ifath(i),2) + fres(i,5)
         xl1(ifath(i),3) = xl1(ifath(i),3) + fres(i,6)
         xl1(ifath(i),4) = xl1(ifath(i),4) + fres(i,7)
         xl1(ifath(i),5) = xl1(ifath(i),5) + fres(i,8)
         xl1(ifath(i),6) = xl1(ifath(i),6) + fres(i,9)         

         xl2(ifath(i),1) = xl2(ifath(i),1) + fres(i,1)*fres(i,1) 
         xl2(ifath(i),2) = xl2(ifath(i),2) + fres(i,2)*fres(i,2)
         xl2(ifath(i),3) = xl2(ifath(i),3) + fres(i,3)*fres(i,3) 
         xl2(ifath(i),4) = xl2(ifath(i),4) + fres(i,1)*fres(i,2)
         xl2(ifath(i),5) = xl2(ifath(i),5) + fres(i,1)*fres(i,3)
         xl2(ifath(i),6) = xl2(ifath(i),6) + fres(i,2)*fres(i,3)

         ui(ifath(i),1) = ui(ifath(i),1) + fres(i,1)
         ui(ifath(i),2) = ui(ifath(i),2) + fres(i,2)
         ui(ifath(i),3) = ui(ifath(i),3) + fres(i,3)

         snorm(ifath(i)) = snorm(ifath(i)) + strnrm(i)

      enddo

c
c Now  the true fathers and serrogates combine results and update
c each other.
c       
      if(numpe .gt. 1)then
         call drvAllreduce(xnude, xnuder,2*nfath)
         call drvAllreduce(ynude, ynuder,6*nfath)
         call drvAllreduce(xm, xmr,6*nfath)
         call drvAllreduce(xl, xlr,6*nfath)
         call drvAllreduce(xl1, xl1r,6*nfath)
         call drvAllreduce(xl2, xl2r,6*nfath)
         call drvAllreduce(ui, uir,3*nfath)
         call drvAllreduce(snorm, snormr,nfath)

         do i = 1, nfath
            ynuder(i,6) = ( ynuder(i,4) - fwr1*ynuder(i,2) ) /
     &           ( two*ynuder(i,5) - four*fwr1*ynuder(i,3)
     &           + two*fwr1*fwr1*ynuder(i,1) )
         enddo

         cdelsq2(:) = ynuder(ifath(:),6)  ! For comparison w/ cdelsq
c
c  xnude is the sum of the sons for each father on this processor
c
c  xnuder is the sum of the sons for each father on all processor combined
c  (the same as if we had not partitioned the mesh for each processor)
c
c   For each father we have precomputed the number of sons (including
c   the sons off processor). 
c
c   Now divide by number of sons to get the average (not really necessary
c   for dynamic model since ratio will cancel nsons at each father)
c
         xnuder(:,1) = xnuder(:,1) / nsons(:)
         xnuder(:,2) = xnuder(:,2) / nsons(:)

         do m = 1, 5
         ynuder(:,m) = ynuder(:,m)/nsons(:)
         enddo
         do m = 1,6
         xmr(:,m) = xmr(:,m)/nsons(:)
         xlr(:,m) = xlr(:,m)/nsons(:)
         xl1r(:,m) = xl1r(:,m)/nsons(:)
         xl2r(:,m) = xl2r(:,m)/nsons(:)
         enddo

         uir(:,1) = uir(:,1)/nsons(:)
         uir(:,2) = uir(:,2)/nsons(:)
         uir(:,3) = uir(:,3)/nsons(:)

         snormr(:) = snormr(:)/nsons(:)

c
cc  the next line is c \Delta^2
cc
cc         xnuder(:,1) = xnuder(:,1) / (xnuder(:,2) + 1.d-09)
cc         do i = 1,nshg
cc            cdelsq(i) = xnuder(ifath(i),1)
cc         enddo

            numNden(:,1) = whist*numNden(:,1)+wcur*xnuder(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnuder(ifath(:),2)
            cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)
            
c            cdelsq(:) = xnuder(ifath(:),1)/(xnuder(ifath(:),2)+1.d-09)

            xnd(:,1) = xnd(:,1) + xnuder(:,1)
            xnd(:,2) = xnd(:,2) + xnuder(:,2)

            xmodcomp(:,1) = xmodcomp(:,1)+ynuder(:,1)
            xmodcomp(:,2) = xmodcomp(:,2)+ynuder(:,2)            
            xmodcomp(:,3) = xmodcomp(:,3)+ynuder(:,3)
            xmodcomp(:,4) = xmodcomp(:,4)+ynuder(:,4)
            xmodcomp(:,5) = xmodcomp(:,5)+ynuder(:,5)

            xmcomp(:,:) = xmcomp(:,:)+xmr(:,:)
            xlcomp(:,:) = xlcomp(:,:)+xlr(:,:)

            xl1comp(:,:) = xl1comp(:,:)+xl1r(:,:)
            xl2comp(:,:) = xl2comp(:,:)+xl2r(:,:)

            ucomp(:,:) = ucomp(:,:)+uir(:,:)
            u1 = uir(32,1)
            scomp(:)   = scomp(:)+snormr(:)

      else

         xnude(:,1) = xnude(:,1)/nsons(:)
         xnude(:,2) = xnude(:,2)/nsons(:)

         do m = 1, 5
         ynude(:,m) = ynude(:,m)/nsons(:)
         enddo
         do m = 1,6
         xm(:,m) = xm(:,m)/nsons(:)
         xl(:,m) = xl(:,m)/nsons(:)
         xl1(:,m) = xl1(:,m)/nsons(:)
         xl2(:,m) = xl2(:,m)/nsons(:)
         enddo

         ui(:,1) = ui(:,1)/nsons(:)
         ui(:,2) = ui(:,2)/nsons(:)
         ui(:,3) = ui(:,3)/nsons(:)

         snorm(:) = snorm(:)/nsons(:)
c     
c     the next line is c \Delta^2, not nu_T but we want to save the
c     memory
c     

cc         xnude(:,1) = xnude(:,1) / (xnude(:,2) + 1.d-09)
cc        do i = 1,nshg
cc            cdelsq(i) = xnude(ifath(i),1)
cc         enddo
cc      endif

         do i = 1, nfath
            ynude(i,6) = ( ynude(i,4) - fwr1*ynude(i,2) ) /
     &           ( two*ynude(i,5) - four*fwr1*ynude(i,3)
     &           + fwr1*fwr1*ynude(i,1) )
         enddo

            numNden(:,1) = whist*numNden(:,1)+wcur*xnude(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnude(ifath(:),2)

            xnd(:,1) = xnd(:,1)+xnude(:,1)
            xnd(:,2) = xnd(:,2)+xnude(:,2)

            cdelsq(:) = numNden(:,1) / (numNden(:,2)) ! + 1.d-09)

c            cdelsq(:) = xnude(ifath(:),1)/(xnude(ifath(:),2))!+1.d-09)
            

          cdelsq2(:) = ynude(ifath(:),6)  ! For comparison w/ cdelsq

            xmodcomp(:,1) = xmodcomp(:,1)+ynude(:,1)
            xmodcomp(:,2) = xmodcomp(:,2)+ynude(:,2)            
            xmodcomp(:,3) = xmodcomp(:,3)+ynude(:,3)
            xmodcomp(:,4) = xmodcomp(:,4)+ynude(:,4)
            xmodcomp(:,5) = xmodcomp(:,5)+ynude(:,5)

            xmcomp(:,:) = xmcomp(:,:)+xm(:,:)
            xlcomp(:,:) = xlcomp(:,:)+xl(:,:)

            xl1comp(:,:) = xl1comp(:,:)+xl1(:,:)
            xl2comp(:,:) = xl2comp(:,:)+xl2(:,:)

            ucomp(:,:) = ucomp(:,:)+ui(:,:)
            u1 = ui(32,1)
            scomp(:)   = scomp(:)+snorm(:)

         endif


c         do i = 1, nfath
c            xmodcomp(i,:) = xmodcomp(i,:)/nsons(i)
c            xmcomp(i,:) = xmcomp(i,:)/nsons(i)         
c            xlcomp(i,:) = xlcomp(i,:)/nsons(i)
c            xl2comp(i,:) = xl2comp(i,:)/nsons(i)         
c            xl1comp(i,:) = xl1comp(i,:)/nsons(i)
c            xnd(i,:) = xnd(i,:)/nsons(i)
c            scomp(i) = scomp(i)/nsons(i)
c            ucomp(i,:) = ucomp(i,:)/nsons(i)
c         enddo

         if ( istep .eq. (nstep(1)-1) ) then
         if ( myrank .eq. master) then

            do i = 1, nfath
            write(365,*)xmodcomp(i,1),xmodcomp(i,2),xmodcomp(i,3),
     &              xmodcomp(i,4),xmodcomp(i,5)

            write(366,*)xmcomp(i,1),xmcomp(i,2),xmcomp(i,3)
            write(367,*)xmcomp(i,4),xmcomp(i,5),xmcomp(i,6)            

            write(368,*)xlcomp(i,1),xlcomp(i,2),xlcomp(i,3)
            write(369,*)xlcomp(i,4),xlcomp(i,5),xlcomp(i,6)

            write(370,*)xl1comp(i,1),xl1comp(i,2),xl1comp(i,3)
            write(371,*)xl1comp(i,4),xl1comp(i,5),xl1comp(i,6) 

            write(372,*)xl2comp(i,1),xl2comp(i,2),xl2comp(i,3)
            write(373,*)xl2comp(i,4),xl2comp(i,5),xl2comp(i,6)

            write(374,*)xnd(i,1),xnd(i,2),scomp(i)
            write(375,*)ucomp(i,1),ucomp(i,2),ucomp(i,3) 

c            write(*,*)'uit uic=', ucomp(32,1),u1
            enddo


            call flush(365)
            call flush(366)
            call flush(367)
            call flush(368)
            call flush(369)
            call flush(370)
            call flush(371)
            call flush(372)
            call flush(373)
            call flush(374)
            call flush(375)

c            if (myrank .eq. master) then
c               write(*,*)'uit uic=', ucomp(32,1),u1
c            endif


c            close(852)
c            close(853)
c            close(854)

         endif
         endif

            if (myrank .eq. master) then
               write(*,*)'uit uic=', ucomp(32,1),u1
            endif


 555     format(e14.7,4(2x,e14.7))
 556     format(e14.7,5(2x,e14.7))

c         close(849)
c         close(850)         
c         close(851)
c         close(852)  
c         close(853)
c         close(854)          

c $$$$$$$$$$$$$$$$$$$$$$$$$$$
      tmp1 =  MINVAL(cdelsq)
      tmp2 =  MAXVAL(cdelsq)
      if(numpe>1) then
         call MPI_REDUCE (tmp1, tmp3, 1,MPI_DOUBLE_PRECISION,
     &        MPI_MIN, master, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
     &        MPI_MAX, master, MPI_COMM_WORLD, ierr)
         tmp1=tmp3
         tmp2=tmp4
      endif
      if (myrank .EQ. master) then !print CDelta^2 range
         write(34,*)lstep,tmp1,tmp2
         call flush(34)
      endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$
      
      if (myrank .eq. master) then
         write(*,*) 'cdelsq=', cdelsq(1),cdelsq(2)
         write(*,*) 'cdelsq=', cdelsq2(1),cdelsq2(2)
         write(22,*) lstep, cdelsq(1)
         call flush(22)
      endif

      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         
         ngauss = nint(lcsyst)

         call scatnu (mien(iblk)%p, strl(iel:inum,:), 
     &        mxmudmi(iblk)%p,cdelsq,shp(lcsyst,1:nshl,:))
      enddo
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$  tmp1 =  MINVAL(xmudmi)
c$$$  tmp2 =  MAXVAL(xmudmi)
c$$$  if(numpe>1) then
c$$$  call MPI_REDUCE (tmp1, tmp3, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MIN, master, MPI_COMM_WORLD, ierr)
c$$$  call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MAX, master, MPI_COMM_WORLD, ierr)
c$$$      tmp1=tmp3
c$$$  tmp2=tmp4
c$$$  endif
c$$$  if (myrank .EQ. master) then
c$$$  write(35,*) lstep,tmp1,tmp2
c$$$  call flush(35)
c$$$  endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$

c
c  if flag set, write a restart file with info (reuse xmij's memory)
c
      if(irs.eq.11) then
         lstep=999
         xmij(:,1)=xnum(:)
         xmij(:,2)=xden(:)
         xmij(:,3)=cdelsq(:)
         xmij(:,5)=xlij(:,4)    !leave M_{12} in 4 and put L_{12} here
         call restar('out ',xmij,xlij) !also dump all of L_{ij} in ac
         stop
      endif
c
c  local clipping moved to scatnu with the creation of mxmudmi pointers
c
c$$$      rmu=datmat(1,2,1)
c$$$      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
c$$$      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0
c      stop !uncomment to test dmod
c


c  write out the nodal values of xnut (estimate since we don't calc strain
c  there and must use the filtered strain).
c

      if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
c
c  collect the average strain into xnude(2)
c
         xnude(:,2) = zero
         do i = 1,numnp
            xnude(ifath(i),2) = xnude(ifath(i),2) + strnrm(i)
         enddo

         if(numpe .gt. 1) then
             call drvAllreduce(xnude(:,2), xnuder(:,2),nfath)
          else
             xnuder=xnude
          endif
c     
c          nut= cdelsq    * |S|
c 
         xnutf=xnuder(:,1)*xnuder(:,2)/nsons(:)
c
c  collect the x and y coords into xnude
c
         xnude = zero
         do i = 1,numnp
            xnude(ifath(i),1) = xnude(ifath(i),1) + x(i,1)
            xnude(ifath(i),2) = xnude(ifath(i),2) + x(i,2)
         enddo

         if(numpe .gt. 1) 
     &        call drvAllreduce(xnude, xnuder,2*nfath)
         xnuder(:,1)=xnuder(:,1)/nsons(:)
         xnuder(:,2)=xnuder(:,2)/nsons(:)
c
c  xnude is the sum of the sons for each father on this processor
c
         if((myrank.eq.master)) then
            do i=1,nfath      ! cdelsq   * |S|
               write(444,*) xnuder(i,1),xnuder(i,2),xnutf(i)
            enddo
            call flush(444)
         endif
      endif

      return
      end
