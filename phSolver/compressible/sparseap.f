      subroutine SparseAp(iper, ilwork, iBC, col, row,	lhsK,	p)   

C============================================================================
C
C "SparseAp": This routine performs the product of a sparsley stored matrix
C 		lhsK(nflow*nflow, nnz_tot) and a vector p(nshg, nflow). The
C		results of the product is returned in p(nshg, nflow).
C
C Nahid Razmara, Spring 2000. (Sparse Matrix)
C============================================================================

c
      include "common.h"
c
c
c.... Data declaration
c
	integer	col(nshg+1),	row(nnz*nshg), 
     &          iper(nshg),     ilwork(nlwork)
	real*8	lhsK(nflow*nflow,nnz_tot)
        real*8 	p(nshg,nflow),	q(nshg,nflow)
c
	real*8	tmp1,	tmp2,	tmp3,	tmp4,	tmp5
c     
c.... communicate:: copy the master's portion of uBrg to each slave
c
      if (numpe > 1) then
         call commu (p, ilwork, nflow  , 'out')
      endif
c
c.... local periodic boundary conditions (no communications)
c
        do j=1,nflow
           p(:,j)=p(iper(:),j)
        enddo
c
c       slave has masters value, for abc we need to rotate it
c        (if this is a vector only no SCALARS)
        if((iabc==1)) !are there any axisym bc's
     &     call rotabc(p(1,2), iBC,  'out')
c
c.... clear the vector
c
	q=zero
c
c.... Perform Matrix-Vector (AP) product
c
 	do i = 1, nshg 
c
	    tmp1 = 0.
	    tmp2 = 0.
	    tmp3 = 0.
	    tmp4 = 0.
	    tmp5 = 0.
c
	    do k = col(i), col(i+1)-1
		j = row(k) 
c
               tmp1 = tmp1 + lhsK(1 ,k)*p(j,1)
     1                     + lhsK(6 ,k)*p(j,2)
     2                     + lhsK(11,k)*p(j,3)
     3                     + lhsK(16,k)*p(j,4)
     4                     + lhsK(21,k)*p(j,5)
               tmp2 = tmp2 + lhsK(2 ,k)*p(j,1)
     1                     + lhsK(7 ,k)*p(j,2)
     2                     + lhsK(12,k)*p(j,3)
     3                     + lhsK(17,k)*p(j,4)
     4                     + lhsK(22,k)*p(j,5)
               tmp3 = tmp3 + lhsK(3 ,k)*p(j,1)
     1                     + lhsK(8 ,k)*p(j,2)
     2                     + lhsK(13,k)*p(j,3)
     3                     + lhsK(18,k)*p(j,4)
     4                     + lhsK(23,k)*p(j,5)
               tmp4 = tmp4 + lhsK(4 ,k)*p(j,1)
     1                     + lhsK(9 ,k)*p(j,2)
     2                     + lhsK(14,k)*p(j,3)
     3                     + lhsK(19,k)*p(j,4)
     4                     + lhsK(24,k)*p(j,5)
               tmp5 = tmp5 + lhsK(5 ,k)*p(j,1)
     1                     + lhsK(10,k)*p(j,2)
     2                     + lhsK(15,k)*p(j,3)
     3                     + lhsK(20,k)*p(j,4)
     4                     + lhsK(25,k)*p(j,5)
c  
	    enddo

	    q(i,1) = q(i,1) + tmp1
	    q(i,2) = q(i,2) + tmp2
	    q(i,3) = q(i,3) + tmp3
	    q(i,4) = q(i,4) + tmp4
	    q(i,5) = q(i,5) + tmp5
	enddo

        

c
	p =  q
c
c
c.... -------------------->   communications <-------------------------
c
c
        if((iabc==1))           !are there any axisym bc's
     &       call rotabc(p(1,2), iBC, 'in ')
c
        if (numpe > 1) then
c
c.... send slave's copy of uBrg to the master
c
           call commu (p  , ilwork, nflow  , 'in ')
c     
c.... nodes treated on another processor are eliminated
c     
           numtask = ilwork(1)
           itkbeg = 1

           do itask = 1, numtask

              iacc   = ilwork (itkbeg + 2)
              numseg = ilwork (itkbeg + 4)

              if (iacc .eq. 0) then
                 do is = 1,numseg
                    isgbeg = ilwork (itkbeg + 3 + 2*is)
                    lenseg = ilwork (itkbeg + 4 + 2*is)
                    isgend = isgbeg + lenseg - 1
                    p(isgbeg:isgend,:) = zero
                 enddo
              endif
            
              itkbeg = itkbeg + 4 + 2*numseg

           enddo
        endif
c
	return
	end

