        subroutine i3LDU (Diag, r, code)
c
c----------------------------------------------------------------------
c
c This routine preforms a Cholesky factorization/solve of a set of 
c symmetric matrices for 3-D computations, used for block diagonal
c preconditioning in the iterative driver.
c 
c input:
c  Diag (numnp,nsymdf)  : block diagonal (symmetric storage)
c  r    (numnp,nflow)    : residual
c  code                 : operation code
c                           .eq. 'LDU_Fact', Cholesky Factor
c                           .eq. 'forward ', forward reduction
c                           .eq. 'backward', backward substitution
c                           .eq. 'product ', product Diag.r
c
c output:
c  Diag (numnp,nsymdf)  : Cholesky decomp. of block diagonal
c  r    (numnp,nflow)    : reduced residual
c
c
c Note: the formulation used here to reduce the diagonal block to 
c       symmetric Cholesky triangle is taken from Golub's "Matrix
c       Computations" Book, pages 89 algorithm 5.2-1.  Followed by
c       standard solve.
c
c
c              Diag(1)   Diag(2)     Diag(4)     Diag(7)      Diag(11)
c   T             0      Diag(3)     Diag(5)     Diag(8)      Diag(12)
c  L  =  U  =     0         0        Diag(6)     Diag(9)      Diag(13)
c                 0         0           0        Diag(10)     Diag(14)
c                 0         0           0           0         Diag(15)
c
c  The diagonal terms 1, 3, 6, 10 and 15 are stored in inverted form.
c
c Farzin Shakib, Spring 1987.
c Zdenek Johan,  Fall   1989.  (Modified for option 'product')
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension Diag(nshg,nsymdf),        r(nshg,nflow)
c
        character*8 code
c
c.... perform Cholesky decomposition with the Diagonal terms inverted
c
        if (code .eq. 'LDU_Fact') then
c
          Diag(:, 1) = one   / sqrt (Diag(:, 1))
c
          Diag(:, 2) = Diag(:, 1) *  Diag(:, 2)
          Diag(:, 3) = Diag(:, 3) -  Diag(:, 2) * Diag(:, 2)
          Diag(:, 3) = one   / sqrt (Diag(:, 3))
c
          Diag(:, 4) = Diag(:, 1) *  Diag(:, 4)
          Diag(:, 5) = Diag(:, 3) * (Diag(:, 5)
     &                            -  Diag(:, 4) * Diag(:, 2))
          Diag(:, 6) = Diag(:, 6) -  Diag(:, 4) * Diag(:, 4)
     &                            -  Diag(:, 5) * Diag(:, 5)
          Diag(:, 6) = one   / sqrt (Diag(:, 6))
c
          Diag(:, 7) = Diag(:, 1) *  Diag(:, 7)
          Diag(:, 8) = Diag(:, 3) * (Diag(:, 8)
     &                            -  Diag(:, 7) * Diag(:, 2))
          Diag(:, 9) = Diag(:, 6) * (Diag(:, 9) 
     &                            -  Diag(:, 7) * Diag(:, 4) 
     &                            -  Diag(:, 8) * Diag(:, 5))
c
          Diag(:,10) = Diag(:,10) -  Diag(:, 7) * Diag(:, 7)
     &                            -  Diag(:, 8) * Diag(:, 8)
     &                            -  Diag(:, 9) * Diag(:, 9)
          Diag(:,10) = one   / sqrt (Diag(:,10))
c
          Diag(:,11) = Diag(:, 1) *  Diag(:,11)
          Diag(:,12) = Diag(:, 3) * (Diag(:,12)
     &                            -  Diag(:,11) * Diag(:, 2))
          Diag(:,13) = Diag(:, 6) * (Diag(:,13)
     &                            -  Diag(:,11) * Diag(:, 4)
     &                            -  Diag(:,12) * Diag(:, 5))
          Diag(:,14) = Diag(:,10) * (Diag(:,14)
     &                            -  Diag(:,11) * Diag(:, 7)
     &                            -  Diag(:,12) * Diag(:, 8)
     &                            -  Diag(:,13) * Diag(:, 9))
c
          Diag(:,15) = Diag(:,15) -  Diag(:,11) * Diag(:,11)
     &                            -  Diag(:,12) * Diag(:,12)
     &                            -  Diag(:,13) * Diag(:,13)
     &                            -  Diag(:,14) * Diag(:,14)
          Diag(:,15) = one   / sqrt (Diag(:,15))
c
c.... flop count
c
!      flops = flops + 110*nshg
c
          return
        endif
c
c.... perform forward reduction
c
        if (code .eq. 'forward ') then
c
          r(:,1) = Diag(:, 1) *   r(:,1)
          r(:,2) = Diag(:, 3) * ( r(:,2) 
     &                        -   r(:,1) * Diag(:, 2) )
          r(:,3) = Diag(:, 6) * ( r(:,3)
     &                        -   r(:,1) * Diag(:, 4)
     &                        -   r(:,2) * Diag(:, 5) )
          r(:,4) = Diag(:,10) * ( r(:,4)
     &                        -   r(:,1) * Diag(:, 7)
     &                        -   r(:,2) * Diag(:, 8)
     &                        -   r(:,3) * Diag(:, 9) )
          r(:,5) = Diag(:,15) * ( r(:,5)
     &                        -   r(:,1) * Diag(:,11)
     &                        -   r(:,2) * Diag(:,12)
     &                        -   r(:,3) * Diag(:,13)
     &                        -   r(:,4) * Diag(:,14) )
c
c.... flop count
c
!      flops = flops + 25*nshg
c
          return
        endif
c
c.... perform backward substitution
c
        if (code .eq. 'backward') then
c
          r(:,5) = Diag(:,15) *   r(:,5)
          r(:,4) = Diag(:,10) * ( r(:,4)
     &                        -   r(:,5) * Diag(:,14) )
          r(:,3) = Diag(:, 6) * ( r(:,3) 
     &                        -   r(:,5) * Diag(:,13)
     &                        -   r(:,4) * Diag(:, 9) )
          r(:,2) = Diag(:, 3) * ( r(:,2)
     &                        -   r(:,5) * Diag(:,12)
     &                        -   r(:,4) * Diag(:, 8)
     &                        -   r(:,3) * Diag(:, 5) )
          r(:,1) = Diag(:, 1) * ( r(:,1)
     &                        -   r(:,5) * Diag(:,11)
     &                        -   r(:,4) * Diag(:, 7)
     &                        -   r(:,3) * Diag(:, 4)
     &                        -   r(:,2) * Diag(:, 2) )
c
c.... flop count
c
!      flops = flops + 25*nshg
c
          return
        endif
c
c.... perform product U.r
c
        if (code .eq. 'product ') then
c
          r(:,1) = r(:,1) / Diag(:, 1) + r(:,2) * Diag(:, 2) +
     &             r(:,3) * Diag(:, 4) + r(:,4) * Diag(:, 7) +
     &             r(:,5) * Diag(:,11)
          r(:,2) = r(:,2) / Diag(:, 3) + r(:,3) * Diag(:, 5) +
     &             r(:,4) * Diag(:, 8) + r(:,5) * Diag(:,12)
          r(:,3) = r(:,3) / Diag(:, 6) + r(:,4) * Diag(:, 9) +
     &             r(:,5) * Diag(:,13)
          r(:,4) = r(:,4) / Diag(:,10) + r(:,5) * Diag(:,14)
          r(:,5) = r(:,5) / Diag(:,15)
c
c.... flop count
c
!      flops = flops + 40*nshg
c
          return
        endif
c
        call error ('i3LDU   ', code, 0)
c
c.... return
c
        return
        end
