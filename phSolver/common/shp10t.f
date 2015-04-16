        subroutine shp10T (QuadPt, nQpt,   shp,    shgl,   wght)
c
c----------------------------------------------------------------------
c
c  This subroutine generates shape functions for the 10-node
c   tetrahedra.
c
c input:
c  QuadPt (4,nQpt)              : quadrature points' local coord's
c                                   QuadPt(1,*) : r
c                                   QuadPt(2,*) : s
c                                   QuadPt(3,*) : t
c                                   QuadPt(4,*) : wght
c  nQpt                         : number of quadrature points
c
c output:
c  shp    (nen,nQpt)            : shape functions
c  shgl   (nsd,nen,nQpt)        : local-gradient of shape function 
c  wght   (nQpt)                : quadrature weights
c
c
c shape-functions:
c  N1 = 1 - r - s - t
c  N2 = r
c  N3 = s
c  N4 = t
c
c Note: To be compatible with design of Tau and DC, the local 
c       gradients are divided by 2.  This is equivalent to having
c       r=[-1,1], s=[-1,1] and t=[-1,1], without really changing
c       r, s and t points (i.e., Qpt is for r=[0,1], s=[0,1] and
c       t=[0,1] range)
c
c Zdenek Johan, Summer 1990.
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension QuadPt(4,*),                   shp(nen,*),
     &            shgl(nsd,nen,*),               wght(*)
c
        call error ('shape   ','not-impl',10)
c
c.... loop through quadrature points
c
        do m = 1, nQpt
c
c.... generate the local-shape-functions
c
          shp(1,m) = one - QuadPt(1,m) - QuadPt(2,m) - QuadPt(3,m)
          shp(2,m) = QuadPt(1,m)
          shp(3,m) = QuadPt(2,m)
          shp(4,m) = QuadPt(3,m)
c
c.... generate the grad-local-shape-functions
c
          shgl(1,1,m) = -pt5
          shgl(2,1,m) = -pt5
          shgl(3,1,m) = -pt5
          shgl(1,2,m) =  pt5
          shgl(2,2,m) =  zero
          shgl(3,2,m) =  zero
          shgl(1,3,m) =  zero
          shgl(2,3,m) =  pt5
          shgl(3,3,m) =  zero
          shgl(1,4,m) =  zero
          shgl(2,4,m) =  zero
          shgl(3,4,m) =  pt5
c
c.... copy the weight
c
          wght(m) = QuadPt(4,m)
c
c.... end of shape-function loop
c
        enddo
c
c.... return
c
        return
        end
