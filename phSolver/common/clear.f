        subroutine clear (clr, n)
c
c----------------------------------------------------------------------
c
c  This routine clears a floating point array.
c
c input:
c  n            : number of floating points to be zeroed
c
c output:
c  clr (n)      : the array to be zeroed
c
c Farzin Shakib, Summer 1985.
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension clr(n)
c
        do i = 1, n
          clr(i) = zero
        enddo
c
        return
        end
