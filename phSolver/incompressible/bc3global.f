	subroutine bc3global (globMas, iBC)  
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of LHS mass matrix for a single 
c element.
c
c input:
c  iBC   (nshg) 	: boundary condition code
c  BC    (nshg,11)     : Dirichlet BC constraint parameters
c  ien   (npro,nshl)	: ien array for this element
c  EGmass(npro,nedof,nedof) : element consistent mass matrix before BC
c
c output:
c  EGmass(npro,nedof,nedof): LHS mass matrix after BC is satisfied
c
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Spring 1990. (Modified for general divariant gas)
c----------------------------------------------------------------------
c
        include "common.h"
c
	dimension iBC(nshg),
     &            globMas(4*nshg,4*nshg)    


	do in=1,nshg
	   i0 = (in-1)*4
c
c.... pressure
c
	   if ( btest(iBC(in),2) ) then
	      globMas(i0+1,:) = zero
	      globMas(:,i0+1) = zero
	      globMas(i0+1,i0+1) = one
	   endif
c       
c....   velocities
c       
c       
c....   x1-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 1 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	   endif
c       
c....   x2-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 2 ) then
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one

	   endif
c       
c....   x1-velocity and x2-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 3 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one
	   endif
c       
c....   x3-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 4 ) then
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one

	   endif
c       
c....   x1-velocity and x3-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 5 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one

	   endif
c       
c....   x2-velocity and x3-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 6 ) then
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one

	   endif
c       
c....   x1-velocity, x2-velocity, and x3-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 7 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one
	   endif
	enddo
	

c       
c....   return
c       
	return
	end
