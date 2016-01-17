      subroutine tnanq (u, n, arrname)

      include "common.h"

      dimension   u(nshg,n),rnan(2)
      character(len=*) arrname

      nnanq = 0
      nlarge = 0
      DO j = 1,n
        DO i = 1,nshg
          if (abs(u(i,j)).gt.1.0e10)  nlarge=nlarge+1
          if (u(i,j) .ne. u(i,j)) then
            nnanq = nnanq + 1
            u(i,j)=9.876543e21
          endif 
        ENDDO
      ENDDO

      rnan(1)=nnanq
      rnan(2)=nlarge
      call sumgatN(rnan,2,summed,1)
 
      if (summed.ge.1) then
        close(1001)        !Hack to close the varts files and flush the buffers. 
        call write_restart(myrank,9876543,nshg,n,u,u)
        call error('tnanq   ',arrname,nnanq)
      endif
      return

      end
