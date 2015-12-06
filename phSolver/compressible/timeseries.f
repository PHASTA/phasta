
      subroutine timeseries(ycl, xl, ien, sgn)

      use timedataC
      include "common.h"

      dimension shape(nshl), ycl(npro,nshl,ndofl),
     &     ien(npro,nshl), xl(npro,nenl,nsd),         
     &     sgn(npro,nshl)
      real*8 al(npro,nenl,nsd), 
     &     zi0(npro,nsd), detaij(npro), dzi0(npro,nsd),
     &     m11(npro), m12(npro), m13(npro), m21(npro), m22(npro),
     &     m23(npro), m31(npro), m32(npro), m33(npro), 
     &     r1(npro), r2(npro), r3(npro), shgradl(nshl,nsd)
     
      real*8 xts1, xts2, xts3
      real*8 soln(ndof)
      integer e, founde

      do jj = 1, ntspts
         founde = 0
         if(statptts(jj,1).gt.0) then
            if(statptts(jj,1).eq.iblkts) then
               if(lcsyst.eq.2) then ! hex
                  call shphex (ipord, parptts(jj,:),shape(:),
     &                 shgradl(:,:))
               elseif(lcsyst.eq.1) then
                  call shptet (ipord, parptts(jj,:),shape(:),
     &                 shgradl(:,:))
               endif
               founde=statptts(jj,2)
            endif
         else
            xts1 = ptts(jj,1)
            xts2 = ptts(jj,2)
            xts3 = ptts(jj,3)

            if(lcsyst.eq.2) then ! hex

               call get_a_not_hex(xl,al) ! get mapping poly. coeff.
         
c...  get initial guess for Newton Iteration procedure
         
               detaij(:) = -al(:,2,1)*al(:,3,2)*al(:,4,3) + 
     &              al(:,2,1)*al(:,4,2)*al(:,3,3) + al(:,2,2)*
     &              al(:,3,1)*al(:,4,3) - al(:,2,2)*al(:,4,1)*
     &              al(:,3,3) - al(:,2,3)*al(:,3,1)*al(:,4,2)+
     &              al(:,2,3)*al(:,4,1)*al(:,3,2)
            
               detaij = 1./detaij
            
               zi0(:,1) = detaij(:)*((al(:,4,2)*al(:,3,3)
     &              - al(:,3,2)*al(:,4,3))*(xts1-al(:,1,1)) +
     &              (al(:,3,1)*al(:,4,3)
     &              - al(:,4,1)*al(:,3,3))*(xts2-al(:,1,2)) +
     &              (al(:,4,1)*al(:,3,2)
     &              - al(:,3,1)*al(:,4,2))*(xts3-al(:,1,3)))
            
            
               zi0(:,2) = detaij(:)*((al(:,2,2)*al(:,4,3)
     &              - al(:,4,2)*al(:,2,3))*(xts1-al(:,1,1)) +
     &              (al(:,4,1)*al(:,2,3)
     &              - al(:,2,1)*al(:,4,3))*(xts2-al(:,1,2)) +
     &              (al(:,2,1)*al(:,4,2)
     &              - al(:,4,1)*al(:,2,2))*(xts3-al(:,1,3)))
            
               zi0(:,3) = detaij(:)*((al(:,3,2)*al(:,2,3)
     &              - al(:,2,2)*al(:,3,3))*(xts1-al(:,1,1)) +
     &              (al(:,2,1)*al(:,3,3)
     &              - al(:,3,1)*al(:,2,3))*(xts2-al(:,1,2)) +
     &              (al(:,3,1)*al(:,2,2)
     &              - al(:,2,1)*al(:,3,2))*(xts3-al(:,1,3)))
            
            
c...  iterate to convergence
            
               do it = 1, iterat
               
c...  build matrix
               
                  m11(:)=al(:,2,1)+al(:,5,1)*zi0(:,2)+al(:,7,1)*zi0(:,3)
     &                 +al(:,8,1)*zi0(:,2)*zi0(:,3)
                  m12(:)=al(:,3,1)+al(:,5,1)*zi0(:,1)+al(:,6,1)*zi0(:,3)
     &                 +al(:,8,1)*zi0(:,1)*zi0(:,3)
                  m13(:)=al(:,4,1)+al(:,6,1)*zi0(:,2)+al(:,7,1)*zi0(:,1)
     &                 +al(:,8,1)*zi0(:,1)*zi0(:,2)

                  m21(:)=al(:,2,2)+al(:,5,2)*zi0(:,2)+al(:,7,2)*zi0(:,3)
     &                 +al(:,8,2)*zi0(:,2)*zi0(:,3)
                  m22(:)=al(:,3,2)+al(:,5,2)*zi0(:,1)+al(:,6,2)*zi0(:,3)
     &                 +al(:,8,2)*zi0(:,1)*zi0(:,3)
                  m23(:)=al(:,4,2)+al(:,6,2)*zi0(:,2)+al(:,7,2)*zi0(:,1)
     &                 +al(:,8,2)*zi0(:,1)*zi0(:,2)

                  m31(:)=al(:,2,3)+al(:,5,3)*zi0(:,2)+al(:,7,3)*zi0(:,3)
     &                 +al(:,8,3)*zi0(:,2)*zi0(:,3)
                  m32(:)=al(:,3,3)+al(:,5,3)*zi0(:,1)+al(:,6,3)*zi0(:,3)
     &                 +al(:,8,3)*zi0(:,1)*zi0(:,3)
                  m33(:)=al(:,4,3)+al(:,6,3)*zi0(:,2)+al(:,7,3)*zi0(:,1)
     &                 +al(:,8,3)*zi0(:,1)*zi0(:,2)
               
               
c...  build rhs
               
                  r1(:)=al(:,1,1)+al(:,2,1)*zi0(:,1)+al(:,3,1)*zi0(:,2)+
     &                 al(:,4,1)*zi0(:,3)+al(:,5,1)*zi0(:,1)*zi0(:,2)+
     &                 al(:,6,1)*zi0(:,2)*zi0(:,3)+al(:,7,1)*
     &                 zi0(:,1)*zi0(:,3)+al(:,8,1)*zi0(:,1)*
     &                 zi0(:,2)*zi0(:,3) - xts1

                  r2(:)=al(:,1,2)+al(:,2,2)*zi0(:,1)+al(:,3,2)*zi0(:,2)+
     &                 al(:,4,2)*zi0(:,3)+al(:,5,2)*zi0(:,1)*zi0(:,2)+
     &                 al(:,6,2)*zi0(:,2)*zi0(:,3)+al(:,7,2)*
     &                 zi0(:,1)*zi0(:,3)+al(:,8,2)*zi0(:,1)*
     &                 zi0(:,2)*zi0(:,3) - xts2
               
                  r3(:)=al(:,1,3)+al(:,2,3)*zi0(:,1)+al(:,3,3)*zi0(:,2)+
     &                 al(:,4,3)*zi0(:,3)+al(:,5,3)*zi0(:,1)*zi0(:,2)+
     &                 al(:,6,3)*zi0(:,2)*zi0(:,3)+al(:,7,3)*
     &                 zi0(:,1)*zi0(:,3)+al(:,8,3)*zi0(:,1)*
     &                 zi0(:,2)*zi0(:,3) - xts3

c...  get solution

                  detaij = m11*m22*m33-m11*m23*m32-m21*m12*m33+
     &                 m21*m13*m32+m31*m12*m23-m31*m13*m22  

                  detaij = 1./detaij

                  dzi0(:,1) = -detaij(:)*((m22(:)*m33(:)-m23(:)*m32(:))*
     &                 r1(:) + (m13(:)*m32(:)-m12(:)*m33(:))*r2(:) + 
     &                 (m12(:)*m23(:)-m13(:)*m22(:))*r3(:))
                  dzi0(:,2) = -detaij(:)*((m23(:)*m31(:)-m21(:)*m32(:))*
     &                 r1(:) + (m11(:)*m33(:)-m13(:)*m31(:))*r2(:) + 
     &                 (m13(:)*m21(:)-m11(:)*m23(:))*r3(:))
                  dzi0(:,3) = -detaij(:)*((m21(:)*m32(:)-m22(:)*m31(:))*
     &                 r1(:) + (m12(:)*m31(:)-m11(:)*m32(:))*r2(:) + 
     &                 (m11(:)*m22(:)-m12(:)*m21(:))*r3(:))

                  zi0(:,:) = zi0(:,:) + dzi0(:,:)
               
               enddo
            
               do e = 1, npro
                  if ((abs(zi0(e,1)).lt.(one+tolpt)).and.
     &                 (abs(zi0(e,2)).lt.(one+tolpt)).and.
     &                 (abs(zi0(e,3)).lt.(one+tolpt))) then ! got the element

                     call shphex (ipord, zi0(e,:),shape(:),
     &                    shgradl(:,:))
                  
                     founde=e
                     exit
                  endif
            
               enddo
            elseif (lcsyst.eq.1) then !tet

               call get_a_not_tet(xl,al)
                       
c
c solve for r, s, t  for each elements
c 
               do e = 1, npro
                  detaij(e) = al(e,2,1)*(-al(e,3,2)*al(e,4,3) + 
     &                 al(e,4,2)*al(e,3,3)) + al(e,2,2)*
     &                 (al(e,3,1)*al(e,4,3) - al(e,4,1)*
     &                 al(e,3,3)) + al(e,2,3)*(-al(e,3,1)*al(e,4,2)+
     &                 al(e,4,1)*al(e,3,2))
            
                  detaij(e) = 1./detaij(e)

                  zi0(e,1) = detaij(e)*((al(e,4,2)*al(e,3,3)
     &                 - al(e,3,2)*al(e,4,3))*(xts1-al(e,1,1)) +
     &                 (al(e,3,1)*al(e,4,3)
     &                 - al(e,4,1)*al(e,3,3))*(xts2-al(e,1,2)) +
     &                 (al(e,4,1)*al(e,3,2)
     &                 - al(e,3,1)*al(e,4,2))*(xts3-al(e,1,3)))

                  zi0(e,2) = detaij(e)*((al(e,2,2)*al(e,4,3)
     &                 - al(e,4,2)*al(e,2,3))*(xts1-al(e,1,1)) +
     &                 (al(e,4,1)*al(e,2,3)
     &                 - al(e,2,1)*al(e,4,3))*(xts2-al(e,1,2)) +
     &                 (al(e,2,1)*al(e,4,2)
     &                 - al(e,4,1)*al(e,2,2))*(xts3-al(e,1,3)))

                  zi0(e,3) = detaij(e)*((al(e,3,2)*al(e,2,3)
     &                 - al(e,2,2)*al(e,3,3))*(xts1-al(e,1,1)) +
     &                 (al(e,2,1)*al(e,3,3)
     &                 - al(e,3,1)*al(e,2,3))*(xts2-al(e,1,2)) +
     &                 (al(e,3,1)*al(e,2,2)
     &                 - al(e,2,1)*al(e,3,2))*(xts3-al(e,1,3)))
               
                  if ((zi0(e,1)+zi0(e,2)+zi0(e,3)).lt.(one+tolpt).and.
     &                 zi0(e,1).lt.(one +tolpt).and.    !should not be necessary; the limit in the tet is already defined by the previous line
     &                 zi0(e,1).gt.(zero-tolpt).and.
     &                 zi0(e,2).lt.(one +tolpt).and.    !should not be necessary 
     &                 zi0(e,2).gt.(zero-tolpt).and.
     &                 zi0(e,3).lt.(one +tolpt).and.    !should not be necessary 
     &                 zi0(e,3).gt.(zero-tolpt)) then
                     
                     call shptet (ipord, zi0(e,:), shape(:),
     &                    shgradl(:,:))

                     founde=e
                     exit
                  endif
               enddo
            endif

            if(founde.ne.0) then !store values for next time steps
               statptts(jj,1)=iblkts
               statptts(jj,2)=founde
               parptts(jj,:)=zi0(founde,:)
            else
              if(iblkts.eq.nelblk) then !not found in last elm-blk of this part
                 statptts(jj,1)=nelblk+1 ! not searched for following steps
              endif
            endif
         endif

         if(founde.ne.0) then
            soln(1:ndofl) = zero

            do i = 1,nenl
               soln(1:ndofl) = soln(1:ndofl)
     &              +ycl(founde,i,1:ndofl)*shape(i)
            enddo
            do i = 1+nenl,nshl
               soln(1:ndofl) = soln(1:ndofl)
     &              +ycl(founde,i,1:ndofl)*shape(i)*sgn(founde,i)
            enddo

            varts(jj,:) = soln(:)
         endif
         
      enddo

      return 
      end
      
