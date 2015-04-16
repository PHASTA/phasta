        subroutine bc3per (iBC,  res, iper, ilwork,nQs)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the periodic nodes after Ap product
c
c input:
c  iBC   (nshg)        : Boundary Condition Code
c  iper  (nshg)        : partners of periodic nodes
c  res   (nshg,nflow)        : residual before BC is applied
c
c output:
c  res   (nshg)        : residual after satisfaction of BC
c  
c
c Kenneth Jansen,  Winter 1998. 
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension iBC(nshg),
     &            res(nshg,nQs),           ilwork(nlwork),
     &            iper(nshg)
c
c
c.... local periodic boundary conditions (no communications)
c
        do j = 1,nshg
          if (btest(iBC(j),10)) then
            i = iper(j)
            res(i,:) = res(i,:) + res(j,:)
            res(j,:) = zero
          endif
        enddo
c     
c.... periodic slaves get the residual values of the masters
c
c      do i = 1,nshg
c         if (btest(iBC(i),10)) then
c            res(i,:) = res(iper(i),:)
c         endif
c      enddo       
c
c
c.... return
c
        return
        end
c
c
c
        subroutine bc3perSclr (iBC,  res, iper)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the periodic nodes after Ap product
c
c input:
c  iBC   (nshg)        : Boundary Condition Code
c  iper  (nshg)        : partners of periodic nodes
c  res   (nshg)   : residual before BC is applied
c
c output:
c  res   (nshg)   : residual after satisfaction of BC
c  
c
c Kenneth Jansen,  Winter 1998. 
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension iBC(nshg),
     &            res(nshg),
     &            iper(nshg)
c
c
c.... local periodic boundary conditions (no communications)
c
        do j = 1,nshg
          if (btest(iBC(j),10)) then
            i = iper(j)
            res(i) = res(i) + res(j)
            res(j) = zero  !changed
          endif
        enddo
c     
c.... periodic slaves get the residual values of the masters
c
c$$$      do i = 1,nshg
c$$$         if (btest(iBC(i),10)) then
c$$$            res(i) = res(iper(i))
c$$$         endif
c$$$      enddo       
c
c
c.... return
c
        return
        end

