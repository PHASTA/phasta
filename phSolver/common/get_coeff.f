

c$$$      subroutine get_a_not_hex(xc,anot)
c$$$
c$$$      include "common.h"
c$$$
c$$$      dimension xc(npro,nenl,nsd), anot(npro,nenl,nsd)
c$$$
c$$$
c$$$      do i = 1, nsd
c$$$
c$$$         anot(:,1,i) = pt125*(xc(:,1,i)+xc(:,2,i)+xc(:,3,i)+xc(:,4,i)
c$$$     &        +xc(:,5,i)+xc(:,6,i)+xc(:,7,i)+xc(:,8,i))
c$$$         
c$$$         anot(:,2,i) = pt125*(-xc(:,1,i)+xc(:,2,i)+xc(:,3,i)-xc(:,4,i)
c$$$     &        -xc(:,5,i)+xc(:,6,i)+xc(:,7,i)-xc(:,8,i))
c$$$
c$$$         anot(:,3,i) = pt125*(-xc(:,1,i)-xc(:,2,i)+xc(:,3,i)+xc(:,4,i)
c$$$     &        -xc(:,5,i)-xc(:,6,i)+xc(:,7,i)+xc(:,8,i))
c$$$         
c$$$         anot(:,4,i) = pt125*(-xc(:,1,i)-xc(:,2,i)-xc(:,3,i)-xc(:,4,i)
c$$$     &        +xc(:,5,i)+xc(:,6,i)+xc(:,7,i)+xc(:,8,i))
c$$$
c$$$         anot(:,5,i) = pt125*(xc(:,1,i)-xc(:,2,i)+xc(:,3,i)-xc(:,4,i)
c$$$     &        +xc(:,5,i)-xc(:,6,i)+xc(:,7,i)-xc(:,8,i))
c$$$
c$$$         anot(:,6,i) = pt125*(xc(:,1,i)+xc(:,2,i)-xc(:,3,i)-xc(:,4,i)
c$$$     &        -xc(:,5,i)-xc(:,6,i)+xc(:,7,i)+xc(:,8,i))
c$$$
c$$$         anot(:,7,i) = pt125*(xc(:,1,i)-xc(:,2,i)-xc(:,3,i)+xc(:,4,i)
c$$$     &        -xc(:,5,i)+xc(:,6,i)+xc(:,7,i)-xc(:,8,i))
c$$$
c$$$         anot(:,8,i) = pt125*(-xc(:,1,i)+xc(:,2,i)-xc(:,3,i)+xc(:,4,i)
c$$$     &        +xc(:,5,i)-xc(:,6,i)+xc(:,7,i)-xc(:,8,i))
c$$$
c$$$      enddo
c$$$
c$$$      return
c$$$      end


      subroutine get_coeff_tet(xc,anot)

      use spebc
      include "common.h"
      

      dimension xc(nelint,nenl,nsd), anot(nelint,nenl,nsd)


      do i = 1, nsd

         anot(:,1,i) = xc(:,4,i)
         anot(:,2,i) = xc(:,1,i)-xc(:,4,i)
         anot(:,3,i) = xc(:,2,i)-xc(:,4,i)
         anot(:,4,i) = xc(:,3,i)-xc(:,4,i)
         
      enddo
      
      return
      end
      
