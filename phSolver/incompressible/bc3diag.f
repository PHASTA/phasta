        subroutine bc3Diag (iBC,  BC,  res)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the diagonal vector for 3D elements.
c
c input:
c  iBC   (nshg)        : Boundary Condition Code
c  BC    (nshg,ndofBC) : the boundary condition constraint parameters
c  res   (nshg,nflow)   : residual before BC is applied
c
c output:
c  res   (nshg,nflow)   : residual after satisfaction of BC
c  
c
c Thuc Bui,      Winter 1989.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension iBC(nshg),                BC(nshg,ndofBC),   
     &            res(nshg,4)
c 
c.... pressure 
c
        where (btest(iBC,2))
           res(:,4) = zero
        endwhere
c
c.... velocities
c
c ibits(n1,n2,n3) extracts bits n2+1 through n2+n3 (extending to the left
c as is traditional in binary) of the integer n1
c and returns the base 10 integer. In examples below x y z a b can 
c be 1 or zero without any effect.
c
c.... x1-velocity
c
c if iBC=4   bits of ibc =00000100 => ibits(4,3,3)=0
c if iBC=40  bits of ibc =00101000 => ibits(4,3,3)=5
c if iBC=40  bits of ibc =00101000 => ibits(4,3,2)=1
c
        where (ibits(iBC,3,3) .eq. 1)   ! bits of iBC= xy001zab 
c
c     notice that the extracted 3 bits form the number 1.  below
c     you will see the combinations which make up 2-7, all of the
c     possible velocity combinations
c
          res(:,2) = res(:,2) - BC(:,4) * res(:,1)
          res(:,3) = res(:,3) - BC(:,5) * res(:,1)
          res(:,1) = zero
        endwhere
c
c.... x2-velocity
c
        where (ibits(iBC,3,3) .eq. 2)   ! bits of iBC= xy010zab 
          res(:,1) = res(:,1) - BC(:,4) * res(:,2)
          res(:,3) = res(:,3) - BC(:,5) * res(:,2)
          res(:,2) = zero
        endwhere
c
c.... x1-velocity and x2-velocity
c
        where (ibits(iBC,3,3) .eq. 3)  ! bits of iBC= xy011zab 
          res(:,3) = res(:,3) - BC(:,4) * res(:,1) - BC(:,6) * res(:,2)
          res(:,1) = zero
          res(:,2) = zero
        endwhere
c
c.... x3-velocity
c
        where (ibits(iBC,3,3) .eq. 4)  ! bits of iBC= xy100zab 
          res(:,1) = res(:,1) - BC(:,4) * res(:,3)
          res(:,2) = res(:,2) - BC(:,5) * res(:,3)
          res(:,3) = zero
        endwhere
c
c.... x1-velocity and x3-velocity
c
        where (ibits(iBC,3,3) .eq. 5)  ! bits of iBC= xy101zab 
          res(:,2) = res(:,2) - BC(:,4) * res(:,1) - BC(:,6) * res(:,3)
          res(:,1) = zero
          res(:,3) = zero
        endwhere
c
c.... x2-velocity and x3-velocity
c
        where (ibits(iBC,3,3) .eq. 6)  ! bits of iBC= xy110zab 
          res(:,1) = res(:,1) - BC(:,4) * res(:,2) - BC(:,6) * res(:,3)
          res(:,2) = zero
          res(:,3) = zero
        endwhere
c
c.... x1-velocity, x2-velocity and x3-velocity
c
        where (ibits(iBC,3,3) .eq. 7)  ! bits of iBC= xy111zab 
          res(:,1) = zero
          res(:,2) = zero
          res(:,3) = zero
        endwhere
c
c.... scaled plane extraction boundary condition
c
        if(intpres.eq.1) then  ! interpolating pressure so zero continuity res 
           where (btest(iBC,11))
              res(:,1) = zero
              res(:,2) = zero
              res(:,3) = zero
              res(:,4) = zero
           endwhere
        else  ! leave residual in continuity equation
           where (btest(iBC,11))
              res(:,1) = zero
              res(:,2) = zero
              res(:,3) = zero
           endwhere
        endif

c.... return
c
        return
        end


        subroutine bc3SclrDiag (iBC,  res)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the diagonal vector for 3D elements.
c
c input:
c  iBC   (nshg)        : Boundary Condition Code
c  BC    (nshg,ndofBC) : the boundary condition constraint parameters
c  res   (nshg)   : residual before BC is applied
c
c output:
c  res   (nshg)   : residual after satisfaction of BC
c  
c
c Thuc Bui,      Winter 1989.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension iBC(nshg),                BC(nshg,ndofBC),   
     &            res(nshg)

       if(isclr.eq.0) then
c 
c.... temperature
c
          where (btest(iBC,1))
             res = zero
          endwhere
       else
c 
c.... scalar
c
          ib=5+isclr
          where (btest(iBC,ib))
             res = zero
          endwhere
       endif
c
c.... return
c
        return
        end
