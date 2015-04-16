      subroutine ludcmp(aa,n,np,indx,d)
c---------------------------------------------------------------------
c
c  find the LU decomposition of a matrix. 
c
c  Numerical Recipes: p. 38
c
c---------------------------------------------------------------------
      include "common.h"
      
      dimension aa(np,np), indx(n), vv(500)
      
      d=1.0
      do i=1,n
         aamax=0.0
         do j=1,n
            if (abs(aa(i,j)).gt.aamax) aamax=abs(aa(i,j))
         enddo
         vv(i)=1.0/aamax
      enddo
      
      do j=1,n
         do i=1,j-1
            sum=aa(i,j)
            do k=1,i-1
               sum=sum-aa(i,k)*aa(k,j)
            enddo
            aa(i,j)=sum
         enddo
         aamax = 0.0
         do i=j,n
            sum=aa(i,j)
            do k=1,j-1
               sum=sum-aa(i,k)*aa(k,j)
            enddo
            aa(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         enddo
         if (j.ne.imax) then
            do k=1,n
               dum=aa(imax,k)
               aa(imax,k)=aa(j,k)
               aa(j,k)=dum
            enddo
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(j.ne.n)then
            dum=1.0/aa(j,j)
            do i=j+1,n
               aa(i,j)=aa(i,j)*dum
            enddo
         endif
      enddo

      return
      end
      
      
         
      
      
