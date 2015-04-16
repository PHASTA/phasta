      subroutine elem_search(xl,xts1,xts2,xts3,xsic,elmt, idfile)
c
c----------------------------------------------------------------------
c This subroutine finds for a particular point in which element it lays
c and what are its local coordinates
c
c input:
c  xl      (npro, nshl, nsd)    : local node coordinates
c  xts1, xts2, xts3		: point coordinates	
c
c output:
c  xsic   (nsd)           	: point's local coordinates 
c  elmt  			: element number
c
c
c  Elaine Bohr
c  April 2002
c----------------------------------------------------------------------
c
      use spebc
      include "common.h"

      dimension shape(nelint,nshl), xsic(nsd),
     &     xl(nelint,nenl,nsd)         
     
      real*8 al(nelint,nenl,nsd), 
     &       zi0(nelint,nsd), detaij(nelint), dzi0(nelint,nsd),
     &       neg(nelint), distance(nelint)
      integer testing(4), found(nelint), negativ(nelint)
     
      real*8 xts1, xts2, xts3, min
      integer e, elmt, result, count, counting, find

      tolpt = 0.000001
      call get_coeff_tet(xl,al)
            
      detaij(:) = -al(:,2,1)*al(:,3,2)*al(:,4,3) + 
     &           al(:,2,1)*al(:,4,2)*al(:,3,3) + al(:,2,2)*
     &           al(:,3,1)*al(:,4,3) - al(:,2,2)*al(:,4,1)*
     &           al(:,3,3) - al(:,2,3)*al(:,3,1)*al(:,4,2)+
     &           al(:,2,3)*al(:,4,1)*al(:,3,2)
            
      detaij = 1./detaij
            
      zi0(:,1) = detaij(:)*((al(:,4,2)*al(:,3,3)
     &           - al(:,3,2)*al(:,4,3))*(xts1-al(:,1,1)) +
     &           (al(:,3,1)*al(:,4,3)
     &           - al(:,4,1)*al(:,3,3))*(xts2-al(:,1,2)) +
     &           (al(:,4,1)*al(:,3,2)
     &           - al(:,3,1)*al(:,4,2))*(xts3-al(:,1,3)))
            
            
      zi0(:,2) = detaij(:)*((al(:,2,2)*al(:,4,3)
     &           - al(:,4,2)*al(:,2,3))*(xts1-al(:,1,1)) +
     &           (al(:,4,1)*al(:,2,3)
     &           - al(:,2,1)*al(:,4,3))*(xts2-al(:,1,2)) +
     &           (al(:,2,1)*al(:,4,2)
     &           - al(:,4,1)*al(:,2,2))*(xts3-al(:,1,3)))
            
      zi0(:,3) = detaij(:)*((al(:,3,2)*al(:,2,3)
     &           - al(:,2,2)*al(:,3,3))*(xts1-al(:,1,1)) +
     &           (al(:,2,1)*al(:,3,3)
     &           - al(:,3,1)*al(:,2,3))*(xts2-al(:,1,2)) +
     &           (al(:,3,1)*al(:,2,2)
     &           - al(:,2,1)*al(:,3,2))*(xts3-al(:,1,3)))
      
c      zi0(:,4) = 1 - zi0(:,1) - zi0(:,2) - zi0(:,3)
            
      elmt = 0
      counting = 0
      neg(:) = 0
      negativ(:) = 0
      found(:) = 0
c      distance(:) = 0
      do e = 1, nelint
         
	 count = 0  
	 testing(:) = 0    
         if (zi0(e,1).lt.(one+tolpt).and. 
     &          zi0(e,1).gt.(zero-tolpt)) then    
             testing(1) = 1
	     count = count + 1
	 endif
         if (zi0(e,2).lt.(one+tolpt).and. 
     &          zi0(e,2).gt.(zero-tolpt)) then    
             testing(2) = 1
	     count = count + 1
	 endif
         if (zi0(e,3).lt.(one+tolpt).and. 
     &          zi0(e,3).gt.(zero-tolpt)) then    
             testing(3) = 1
	     count = count + 1
	 endif
         if ((1-zi0(e,1)-zi0(e,2)-zi0(e,3)).lt.(one+tolpt).and. 
     &       (1-zi0(e,1)-zi0(e,2)-zi0(e,3)).gt.(zero-tolpt)) then    
             testing(4) = 1
	     count = count + 1
	 endif
     
         result = 1
	 do i = 1, 4
	    result = result*testing(i)
	 enddo
	 
         if (result .eq. 1) then      
             xsic(:) = zi0(e,:)
	     elmt = e
		  
	     return   
         
	 elseif (count .eq. 3) then
	    counting = counting + 1
	    do i = 1, 3
	       if (testing(i) .eq. 0 .and. 
     &		       zi0(e,i) .lt. 0.0) then
                  found(counting) = e
		  neg(counting) = zi0(e,i)
		  negativ(counting) = i
c		  distance(counting) = sqrt(zi0(e,1)*zi0(e,1)+
c    &		          zi0(e,2)*zi0(e,2)+zi0(e,3)*zi0(e,3))
	       endif
	    enddo
	    if (testing(4) .eq. 0) then
	       found(counting) = e
	       neg(counting) = 1-zi0(e,1)-zi0(e,2)-zi0(e,3)
	       negativ(counting) = 4
c	       distance(counting) = sqrt(zi0(e,1)*zi0(e,1)+
c     &		       zi0(e,2)*zi0(e,2)+zi0(e,3)*zi0(e,3))
	    endif
	 endif
	       
            
      enddo
           
      min = neg(1)
      elmt = found(1)
      find = 1
      do i = 2, counting
         if (min .lt. neg(i)) then 
	    min = neg(i)
	    elmt = found(i)
	    find = i
	 endif
      enddo
      if (negativ(find) .eq. 4) then
         prop = zi0(elmt,1)+zi0(elmt,2)+zi0(elmt,3)
	 xsic(:) = zi0(elmt,:) / prop
      else
         xsic(:) = zi0(elmt,:)
         xsic(negativ(find)) = 0.0
      endif
c      call error ('elem_search ', 'outrange', idfile)
         
      return

      
      end
      
