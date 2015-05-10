      subroutine tnanqe (u, n, arrname)

      include "common.h"

      dimension   u(npro,n),rnan(2)
      character*8 arrname

      nnanq = 0
      nlarge=0
      DO j = 1,n
	DO i = 1,npro
          if(abs(u(i,j)).gt.1.0e10)  nlarge=nlarge+1
	  IF (u(i,j) .ne. u(i,j)) then
	    write(*,*) myrank, i,j
             nnanq = nnanq + 1
             u(i,j)=9.876543e21
	endif
	ENDDO
      ENDDO
        rnan(1)=nnanq
        rnan(2)=nlarge
        call sumgatN(rnan,2,summed,1)
      if (summed.ge.1) then
	do i=1,npro
	write(8+myrank,245) (u(i,j), j=1,n)
	enddo 
	call error('tnanqe  ',arrname,nnanq)
      endif
245   format(10(e14.7,2x))
      return
      end
