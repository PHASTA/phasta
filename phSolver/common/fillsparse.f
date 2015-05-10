	subroutine fillsparseI(	iens, xKebe,	lhsK,
     &                               xGoC,      lhsP,
     1			             row,	col)
c
c
c
	include "common.h"
	real*8	xKebe(npro,9,nshl,nshl), xGoC(npro,4,nshl,nshl)
	integer	ien(npro,nshl),	col(nshg+1), row(nshg*nnz)
	real*8	lhsK(9,nnz_tot),	lhsP(4,nnz_tot)
c
	integer	aa,	b,	c,	e,	i,	k,	n
c
	integer sparseloc

	integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
	ien=abs(iens)
c       
c.... Accumulate the lhs
c
	do e = 1, npro ! loop over the elements
	    do aa = 1, nshl ! loop over the local equation numbers
		i = ien(e,aa) ! finds the global equation number or
			      ! block-row of our matrix
		c = col(i)    ! starting point to look for the matching column
		n = col(i+1) - c  !length of the list of entries in rowp
		do b = 1, nshl ! local variable number tangent respect
			       ! to
c function that searches row until it finds the match that gives the
c		   global equation number

		    k = sparseloc( row(c), n, ien(e,b) ) + c-1
c
c                                             *         *
c                   dimension egmass(npro,ndof,nenl,ndof,nenl)
c
c compressible      lhsT(1:5,1:5,k)=lhsT(1:5,1:5,k)+egmass(e,1:5,aa,1:5,b)
c
		    lhsK(1,k) = lhsK(1,k) + xKebe(e,1,aa,b)
		    lhsK(2,k) = lhsK(2,k) + xKebe(e,2,aa,b)
		    lhsK(3,k) = lhsK(3,k) + xKebe(e,3,aa,b)
		    lhsK(4,k) = lhsK(4,k) + xKebe(e,4,aa,b)
		    lhsK(5,k) = lhsK(5,k) + xKebe(e,5,aa,b)
		    lhsK(6,k) = lhsK(6,k) + xKebe(e,6,aa,b)
		    lhsK(7,k) = lhsK(7,k) + xKebe(e,7,aa,b)
		    lhsK(8,k) = lhsK(8,k) + xKebe(e,8,aa,b)
		    lhsK(9,k) = lhsK(9,k) + xKebe(e,9,aa,b)
c
		    lhsP(1,k) = lhsP(1,k) + xGoC(e,1,aa,b)
		    lhsP(2,k) = lhsP(2,k) + xGoC(e,2,aa,b)
		    lhsP(3,k) = lhsP(3,k) + xGoC(e,3,aa,b)
		    lhsP(4,k) = lhsP(4,k) + xGoC(e,4,aa,b)
		enddo
	    enddo
	enddo
c
c.... end
c
	return
	end


	subroutine fillsparseC(	iens, EGmass,	lhsK,
     1			             row,	col)
c
c-----------------------------------------------------------
c This routine fills up the spasely stored LHS mass matrix
c
c Nahid Razmra, Spring 2000. (Sparse Matrix)
c-----------------------------------------------------------
c
c

	include "common.h"

 	real*8	EGmass(npro,nedof,nedof)
        integer	ien(npro,nshl),	col(nshg+1), row(nnz*nshg)
        real*8 lhsK(nflow*nflow,nnz_tot)

c
        integer aa, b, c, e, i, k, n, n1
        integer f, g, r, s, t
c
	integer sparseloc

	integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
	ien=abs(iens)
c
c.... Accumulate the lhsK
c
	do e = 1, npro
	    do aa = 1, nshl
		i = ien(e,aa)
		c = col(i)
		n = col(i+1) - c
		r = (aa-1)*nflow
		do b = 1, nshl
		    s = (b-1)*nflow
                    k = sparseloc( row(c), n, ien(e,b) ) + c-1
c
		    do g = 1, nflow
		       t = (g-1)*nflow
		       do f = 1, nflow

			  lhsK(t+f,k) = lhsK(t+f,k) + EGmass(e,r+f,s+g)
c

                       enddo
                   enddo
		enddo
	    enddo
	enddo
c
c.... end
c
	return
	end

	subroutine fillsparseSclr(   iens,      xSebe,	lhsS,
     1			             row,	col)
c
c
c
	include "common.h"
	real*8	xSebe(npro,nshl,nshl)
	integer	ien(npro,nshl),	col(nshg+1), row(nshg*nnz)
	real*8	lhsS(nnz_tot)	
c
	integer	aa,	b,	c,	e,	i,	k,	n
c
	integer sparseloc

	integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
	ien=abs(iens)
c
c.... Accumulate the lhs
c
	do e = 1, npro
	    do aa = 1, nshl
		i = ien(e,aa)
		c = col(i)
		n = col(i+1) - c
		do b = 1, nshl
		    k = sparseloc( row(c), n, ien(e,b) ) + c-1
c
		    lhsS(k) = lhsS(k) + xSebe(e,aa,b)
		enddo
	    enddo
	enddo
c
c.... end
c
	return
	end

	integer	function sparseloc( list, n, target )

c-----------------------------------------------------------
c This function finds the location of the non-zero elements
c of the LHS matrix in the sparsely stored matrix 
c lhsK(nflow*nflow,nnz*numnp)
c
c Nahid Razmara, Spring 2000. 	(Sparse Matrix)
c-----------------------------------------------------------

	integer	list(n),	n,	target
	integer	rowvl,	rowvh,	rowv

c
c.... Initialize
c
	rowvl = 1
	rowvh = n + 1
c
c.... do a binary search
c
100	if ( rowvh-rowvl .gt. 1 ) then
	    rowv = ( rowvh + rowvl ) / 2
	    if ( list(rowv) .gt. target ) then
		rowvh = rowv
	    else
		rowvl = rowv
	    endif
	    goto 100
	endif
c
c.... return
c
	sparseloc = rowvl
c
	return
	end

