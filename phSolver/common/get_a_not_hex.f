

      subroutine get_a_not_hex(xc,anot)

      include "common.h"

      dimension xc(npro,nenl,nsd), anot(npro,nenl,nsd)


      do i = 1, nsd

         anot(:,1,i) = pt125*(xc(:,1,i)+xc(:,2,i)+xc(:,3,i)+xc(:,4,i)
     &        +xc(:,5,i)+xc(:,6,i)+xc(:,7,i)+xc(:,8,i))
         
         anot(:,2,i) = pt125*(-xc(:,1,i)+xc(:,2,i)+xc(:,3,i)-xc(:,4,i)
     &        -xc(:,5,i)+xc(:,6,i)+xc(:,7,i)-xc(:,8,i))

         anot(:,3,i) = pt125*(-xc(:,1,i)-xc(:,2,i)+xc(:,3,i)+xc(:,4,i)
     &        -xc(:,5,i)-xc(:,6,i)+xc(:,7,i)+xc(:,8,i))
         
         anot(:,4,i) = pt125*(-xc(:,1,i)-xc(:,2,i)-xc(:,3,i)-xc(:,4,i)
     &        +xc(:,5,i)+xc(:,6,i)+xc(:,7,i)+xc(:,8,i))

         anot(:,5,i) = pt125*(xc(:,1,i)-xc(:,2,i)+xc(:,3,i)-xc(:,4,i)
     &        +xc(:,5,i)-xc(:,6,i)+xc(:,7,i)-xc(:,8,i))

         anot(:,6,i) = pt125*(xc(:,1,i)+xc(:,2,i)-xc(:,3,i)-xc(:,4,i)
     &        -xc(:,5,i)-xc(:,6,i)+xc(:,7,i)+xc(:,8,i))

         anot(:,7,i) = pt125*(xc(:,1,i)-xc(:,2,i)-xc(:,3,i)+xc(:,4,i)
     &        -xc(:,5,i)+xc(:,6,i)+xc(:,7,i)-xc(:,8,i))

         anot(:,8,i) = pt125*(-xc(:,1,i)+xc(:,2,i)-xc(:,3,i)+xc(:,4,i)
     &        +xc(:,5,i)-xc(:,6,i)+xc(:,7,i)-xc(:,8,i))

      enddo

      return
      end


      subroutine get_a_not_tet(xc,anot)

      include "common.h"

      dimension xc(npro,nenl,nsd), anot(npro,nenl,nsd)


      do i = 1, nsd

         anot(:,1,i) = xc(:,4,i)
         anot(:,2,i) = xc(:,1,i)-xc(:,4,i)
         anot(:,3,i) = xc(:,2,i)-xc(:,4,i)
         anot(:,4,i) = xc(:,3,i)-xc(:,4,i)
         
      enddo
      
      return
      end
      
