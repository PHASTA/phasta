        subroutine i3LU (Diag, r, code)
c
c----------------------------------------------------------------------
c
c This routine performs a LU factorization/solve of a set of matrices,
c used for block diagonal preconditioning in the iterative driver.
 
c input:
c  Diag (nshg,nflow,nflow)  : block diagonal (symmetric storage)
c  r    (nshg,nflow)    : residual
c  code                 : operation code
c                           .eq. 'LU_Fact ', Cholesky Factor
c                           .eq. 'forward ', forward reduction
c                           .eq. 'backward', backward substitution
c                           .eq. 'product ', product Diag.r
c
c output:
c  Diag (nshg,nflow,nflow)  : LU decomp. of block diagonal
c  r    (nshg,nflow)    : reduced residual
c
c
c Note: the formulation used here is taken from Golub's "Matrix
c       Computations" Book (1984), pages 82-85 algorithm P5.1-3.
c
c       L and U overwrite Diag
c
c       The diagonal terms (i,i) of U are stored in inverted form.
c
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension Diag(nshg,nflow,nflow),        r(nshg,nflow)
c
        character*8 code
        
c
c.... perform LU decomposition with the Diagonal terms inverted
c
        if (code .eq. 'LU_Fact ') then
c
          Diag(:,1,1) = one   / Diag(:,1,1)
c
          Diag(:,2,1) = Diag(:,1,1) *  Diag(:,2,1)
          Diag(:,3,1) = Diag(:,1,1) *  Diag(:,3,1)
          Diag(:,4,1) = Diag(:,1,1) *  Diag(:,4,1)
          Diag(:,5,1) = Diag(:,1,1) *  Diag(:,5,1)
c
          Diag(:,2,2) = Diag(:,2,2) - Diag(:,2,1) *  Diag(:,1,2)
          Diag(:,2,3) = Diag(:,2,3) - Diag(:,2,1) *  Diag(:,1,3)
          Diag(:,2,4) = Diag(:,2,4) - Diag(:,2,1) *  Diag(:,1,4)
          Diag(:,2,5) = Diag(:,2,5) - Diag(:,2,1) *  Diag(:,1,5)
          Diag(:,2,2) = one   / Diag(:,2,2)
c
          Diag(:,3,2) = Diag(:,2,2)*(Diag(:,3,2)
     &                             - Diag(:,3,1) * Diag(:,1,2))
          Diag(:,4,2) = Diag(:,2,2)*(Diag(:,4,2)
     &                             - Diag(:,4,1) * Diag(:,1,2))
          Diag(:,5,2) = Diag(:,2,2)*(Diag(:,5,2)
     &                             - Diag(:,5,1) * Diag(:,1,2))
c

          Diag(:,3,3) = Diag(:,3,3) - Diag(:,3,1) * Diag(:,1,3)
     &                              - Diag(:,3,2) * Diag(:,2,3)
          Diag(:,3,4) = Diag(:,3,4) - Diag(:,3,1) * Diag(:,1,4)
     &                              - Diag(:,3,2) * Diag(:,2,4)
          Diag(:,3,5) = Diag(:,3,5) - Diag(:,3,1) * Diag(:,1,5)
     &                              - Diag(:,3,2) * Diag(:,2,5)
          Diag(:,3,3) = one   / Diag(:,3,3)
c
          Diag(:,4,3) = Diag(:,3,3) *(Diag(:,4,3)
     &                              - Diag(:,4,1) * Diag(:,1,3)
     &                              - Diag(:,4,2) * Diag(:,2,3))
          Diag(:,5,3) = Diag(:,3,3) *(Diag(:,5,3)
     &                              - Diag(:,5,1) * Diag(:,1,3)
     &                              - Diag(:,5,2) * Diag(:,2,3))
c
          Diag(:,4,4) = Diag(:,4,4) - Diag(:,4,1) * Diag(:,1,4)
     &                              - Diag(:,4,2) * Diag(:,2,4)
     &                              - Diag(:,4,3) * Diag(:,3,4)
          Diag(:,4,4) = one / Diag(:,4,4)
c
          Diag(:,5,4) = Diag(:,4,4) *(Diag(:,5,4)
     &                              - Diag(:,5,1) * Diag(:,1,4)
     &                              - Diag(:,5,2) * Diag(:,2,4)
     &                              - Diag(:,5,3) * Diag(:,3,4))
c
          Diag(:,5,5) = Diag(:,5,5) - Diag(:,5,1) * Diag(:,1,5)
     &                              - Diag(:,5,2) * Diag(:,2,5)
     &                              - Diag(:,5,3) * Diag(:,3,5)
     &                              - Diag(:,5,4) * Diag(:,4,5)
          Diag(:,5,5) = one / Diag(:,5,5)
c
c  INACCURATE NOW          flops = flops + 110*nshg
c
          return
        endif
c
c.... perform forward reduction
c
        if (code .eq. 'forward ') then
c
c         r(:,1) =  r(:,1) !no-op
          r(:,2) =  r(:,2) - Diag(:,2,1) *  r(:,1) 
          r(:,3) =  r(:,3) - Diag(:,3,1) *  r(:,1) 
     &                     - Diag(:,3,2) *  r(:,2) 
          r(:,4) =  r(:,4) - Diag(:,4,1) *  r(:,1) 
     &                     - Diag(:,4,2) *  r(:,2) 
     &                     - Diag(:,4,3) *  r(:,3) 
          r(:,5) =  r(:,5) - Diag(:,5,1) *  r(:,1) 
     &                     - Diag(:,5,2) *  r(:,2) 
     &                     - Diag(:,5,3) *  r(:,3) 
     &                     - Diag(:,5,4) *  r(:,4) 
c
c.... flop count
c
c  INACCURATE           flops = flops + 25*nshg
c
          return
        endif
c
c.... perform backward substitution
c
        if (code .eq. 'backward') then
c
          r(:,5) = Diag(:,5,5) *   r(:,5)
          r(:,4) = Diag(:,4,4) * ( r(:,4)
     &                         -   r(:,5) * Diag(:,4,5) )
          r(:,3) = Diag(:,3,3) * ( r(:,3) 
     &                         -   r(:,5) * Diag(:,3,5)
     &                         -   r(:,4) * Diag(:,3,4) )
          r(:,2) = Diag(:,2,2) * ( r(:,2) 
     &                         -   r(:,5) * Diag(:,2,5)
     &                         -   r(:,4) * Diag(:,2,4)
     &                         -   r(:,3) * Diag(:,2,3) )
          r(:,1) = Diag(:,1,1) * ( r(:,1) 
     &                         -   r(:,5) * Diag(:,1,5)
     &                         -   r(:,4) * Diag(:,1,4)
     &                         -   r(:,3) * Diag(:,1,3) 
     &                         -   r(:,2) * Diag(:,1,2) )
c
c.... flop count
c
          flops = flops + 25*nshg

          return
        endif
c
c.... perform product U.r
c
        if (code .eq. 'product ') then
c
          r(:,1) = r(:,1) / Diag(:,1,1) + r(:,2) * Diag(:,1,2) +
     &             r(:,3) * Diag(:,1,3) + r(:,4) * Diag(:,1,4) +
     &             r(:,5) * Diag(:,1,5)
          r(:,2) = r(:,2) / Diag(:,2,2) + 
     &             r(:,3) * Diag(:,2,3) + r(:,4) * Diag(:,2,4) +
     &             r(:,5) * Diag(:,2,5)
          r(:,3) = r(:,3) / Diag(:,3,3) + 
     &             r(:,4) * Diag(:,3,4) +
     &             r(:,5) * Diag(:,3,5)
          r(:,4) = r(:,4) / Diag(:,4,4) + 
     &             r(:,5) * Diag(:,4,5)
          r(:,5) = r(:,5) / Diag(:,5,5)
c
c.... flop count
c
          flops = flops + 40*nshg

          return
        endif
c
        call error ('i3LU    ', code, 0)
c
c.... return
c
c
111     format(5(e14.7,2x))
        return
        end
