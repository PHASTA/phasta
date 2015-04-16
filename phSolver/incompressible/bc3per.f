        subroutine bc3per (iBC,  res, iper, ilwork,nQs)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the periodic nodes after Ap product
c
c input:
c  iBC   (nshg)        : Boundary Condition Code
c  iper  (nshg)        : partners of periodic nodes
c  res   (nshg,nQs)   : residual before BC is applied
c
c output:
c  res   (nshg,nQs)   : residual after satisfaction of BC
c  
c
c Kenneth Jansen,  Winter 1998. 
c----------------------------------------------------------------------
c
        use periodicity  ! this gives you rcount(1:nshg) (real*8)
        include "common.h"
c
        dimension iBC(nshg),
     &            res(nshg,nQs),           ilwork(nlwork),
     &            iper(nshg)
c
c.... local periodic (and axisymmetric) boundary conditions (no communications)
c
           do j = 1,nshg
              if ((btest(iBC(j),10)) .or. (btest(iBC(j),12))) then
                 i = iper(j)
                 res(i,:) = res(i,:) + res(j,:)
                 res(j,:) = zero
              endif
           enddo


        if(numpe.gt.1) then
c
c.... nodes treated on another processor are eliminated
c     
           numtask = ilwork(1)
           itkbeg = 1
           
           do itask = 1, numtask
              
              iacc   = ilwork (itkbeg + 2)
              numseg = ilwork (itkbeg + 4)
              
              if (iacc .eq. 0) then
                 do is = 1,numseg
                    isgbeg = ilwork (itkbeg + 3 + 2*is)
                    lenseg = ilwork (itkbeg + 4 + 2*is)
                    isgend = isgbeg + lenseg - 1
                    res(isgbeg:isgend,:) = zero
                 enddo
              endif
              
              itkbeg = itkbeg + 4 + 2*numseg
              
           enddo
        endif
c
c.... return
c
        return
        end





