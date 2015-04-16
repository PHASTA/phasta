      subroutine getstrl( y, x, ien, strnrm, shgl, shp )

      include "common.h"

      dimension y(nshg,5)
      dimension x(numnp,3),            xl(npro,nenl,3)
      dimension ien(npro,nshl),        yl(npro,nshl,5),
     &          u1(npro),              u2(npro),
     &          u3(npro),              dxdxi(npro,nsd,nsd),
     &          strnrm(npro,maxsh),    dxidx(npro,nsd,nsd),
     &          shgl(nsd,nshl,maxsh),       shg(npro,nshl,nsd),
     &          shp(nshl,maxsh)         
      dimension tmp(npro),             fresli(npro,24)

      call localy (y,      yl,     ien,    5,  'gather  ')
      call localx (x,      xl,     ien,    3,  'gather  ')
c

      if(matflg(1,1).eq.0) then ! compressible
      yl (:,:,1) = yl(:,:,1) / (Rgas * yl(:,:,5)) 
      else
      yl (:,:,1) = one
      endif

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

      fresli=zero
      do i=1,nshl
        fresli(:,22) = fresli(:,22)+shp(i,intp)*yl(:,i,1)  ! density at qpt
c       fresli(:,24) = fresli(:,24)+shp(i,intp)*yl(:,i,5)  !temperature at qpt
      enddo
c
c
c     fresli(:,22)=fresli(:,22)*wght
c     fresli(:,24)=fresli(:,24)*wght


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

      do j=10,12  ! normal strainrate u_{i,i} no sum on i
       ig=j-9
       iv=j-8
       do i=1,nshl
        fresli(:,j) = fresli(:,j)+shg(:,i,ig)*yl(:,i,iv)
       enddo
      enddo

c shear stresses  NOTE  there may be faster ways to do this
c                  check agains CM5 code for speed WTP
       
       do i=1,nshl
        fresli(:,13) = fresli(:,13)+shg(:,i,2)*yl(:,i,2)
     &                             +shg(:,i,1)*yl(:,i,3)
        fresli(:,14) = fresli(:,14)+shg(:,i,3)*yl(:,i,2)
     &                             +shg(:,i,1)*yl(:,i,4)
        fresli(:,15) = fresli(:,15)+shg(:,i,3)*yl(:,i,3)
     &                             +shg(:,i,2)*yl(:,i,4)
       enddo


      fresli(:,13) = pt5 * fresli(:,13)
      fresli(:,14) = pt5 * fresli(:,14)
      fresli(:,15) = pt5 * fresli(:,15)

      strnrm(:,intp) = fresli(:,22) * sqrt(
     &   two * (fresli(:,10)**2 + fresli(:,11)**2 + fresli(:,12)**2)
     &  + four * ( fresli(:,13)**2 + fresli(:,14)**2 + 
     &    fresli(:,15)**2 ) )

      
      enddo !end of loop over integration points

      
      return
      end
