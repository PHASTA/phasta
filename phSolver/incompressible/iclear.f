        subroutine iclear (iclr, n)
c
c----------------------------------------------------------------------
c
c This routine clears an integer array.
c
c input:
c  n            : number of integers to be zeroed
c
c output:
c  iclr (n)     : the array to be zeroed
c
c
c Farzin Shakib, Summer 1985.
c----------------------------------------------------------------------
c
        dimension iclr(n)
c
        do i = 1, n
          iclr(i) = 0
        enddo
c
        return
        end
