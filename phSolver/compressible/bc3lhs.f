      subroutine bc3LHS (iBC,  BC,  iens,  EGmass )
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of LHS mass matrix for a single 
c element.
c
c input:
c  iBC   (nshg) 	: boundary condition code
c  BC    (nshg,11)      : Dirichlet BC constraint parameters
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
	dimension iBC(nshg),      ien(npro,nshl),
     &		  BC(nshg,11),    EGmass(npro,nedof,nedof) 
	integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
	ien=abs(iens)
c
c.... loop over elements
c
        do iel = 1, npro
c
c.... loop over number of shape functions for this element
c
	do inod = 1, nshl
c
c.... set up parameters
c
        in  = ien(iel,inod)
        if (iBC(in) .eq. 0) goto 5000
        ioff = (inod - 1) * nflow
        i1 = ioff + 1
        i2 = ioff + 2
        i3 = ioff + 3
        i4 = ioff + 4
        i5 = ioff + 5
c
c.... pressure
c
	if ( btest(iBC(in),2) ) then

           do i = 1, nedof
              EGmass(iel,i,i1) = zero
              EGmass(iel,i1,i) = zero
           enddo
           EGmass(iel,i1,i1) = one
	endif
c
c.... density  ( not activated yet for LHS matrix )
c

c
c.... velocities
c
c
c.... x1-velocity
c
        if ( ibits(iBC(in),3,3) .eq. 1) then
           do i = 1, nedof
              EGmass(iel,i3,i) = EGmass(iel,i3,i) 
     &                         - BC(in,4) * EGmass(iel,i2,i) 
              EGmass(iel,i4,i) = EGmass(iel,i4,i) 
     &                         - BC(in,5) * EGmass(iel,i2,i)
           enddo
           do i = 1, nedof
              EGmass(iel,i,i3) = EGmass(iel,i,i3) 
     &                         - BC(in,4) * EGmass(iel,i,i2) 
              EGmass(iel,i,i4) = EGmass(iel,i,i4) 
     &                         - BC(in,5) * EGmass(iel,i,i2)
           enddo
           do i = 1, nedof
              EGmass(iel,i,i2) = zero
              EGmass(iel,i2,i) = zero
           enddo
           EGmass(iel,i2,i2) = one
        endif
c
c.... x2-velocity
c
        if ( ibits(iBC(in),3,3) .eq. 2 ) then
           do i = 1, nedof
              EGmass(iel,i2,i) = EGmass(iel,i2,i) 
     &                         - BC(in,4) * EGmass(iel,i3,i) 
              EGmass(iel,i4,i) = EGmass(iel,i4,i) 
     &                         - BC(in,5) * EGmass(iel,i3,i)
           enddo
           do i = 1, nedof
              EGmass(iel,i,i2) = EGmass(iel,i,i2) 
     &                         - BC(in,4) * EGmass(iel,i,i3) 
              EGmass(iel,i,i4) = EGmass(iel,i,i4) 
     &                         - BC(in,5) * EGmass(iel,i,i3)
           enddo
           
           do i = 1, nedof
              EGmass(iel,i,i3) = zero
              EGmass(iel,i3,i) = zero
           enddo
           EGmass(iel,i3,i3) = one
        endif
c     
c.... x1-velocity and x2-velocity
c
        if ( ibits(iBC(in),3,3) .eq. 3 ) then
           do i = 1, nedof
              EGmass(iel,i4,i) = EGmass(iel,i4,i) 
     &                     - BC(in,4) * EGmass(iel,i2,i) 
     &                     - BC(in,6) * EGmass(iel,i3,i)
           enddo
           do i = 1, nedof
              EGmass(iel,i,i4) = EGmass(iel,i,i4) 
     &                     - BC(in,4) * EGmass(iel,i,i2) 
     &                     - BC(in,6) * EGmass(iel,i,i3)
           enddo

           do i = 1, nedof
              EGmass(iel,i,i2) = zero
              EGmass(iel,i2,i) = zero
              EGmass(iel,i,i3) = zero
              EGmass(iel,i3,i) = zero
           enddo
           EGmass(iel,i2,i2) = one
           EGmass(iel,i3,i3) = one
        endif
c
c.... x3-velocity
c
        if ( ibits(iBC(in),3,3) .eq. 4 ) then
           do i = 1, nedof
              EGmass(iel,i2,i) = EGmass(iel,i2,i) 
     &                         - BC(in,4) * EGmass(iel,i4,i) 
              EGmass(iel,i3,i) = EGmass(iel,i3,i) 
     &                         - BC(in,5) * EGmass(iel,i4,i)
           enddo
           do i = 1, nedof
              EGmass(iel,i,i2) = EGmass(iel,i,i2) 
     &                         - BC(in,4) * EGmass(iel,i,i4) 
              EGmass(iel,i,i3) = EGmass(iel,i,i3) 
     &                         - BC(in,5) * EGmass(iel,i,i4)
           enddo

           do i = 1, nedof
              EGmass(iel,i,i4) = zero
              EGmass(iel,i4,i) = zero
           enddo
           EGmass(iel,i4,i4) = one
        endif
c     
c.... x1-velocity and x3-velocity
c
        if ( ibits(iBC(in),3,3) .eq. 5 ) then
           do i = 1, nedof
              EGmass(iel,i3,i) = EGmass(iel,i3,i) 
     &                     - BC(in,4) * EGmass(iel,i2,i) 
     &                     - BC(in,6) * EGmass(iel,i4,i)
           enddo
           do i = 1, nedof
              EGmass(iel,i,i3) = EGmass(iel,i,i3) 
     &                     - BC(in,4) * EGmass(iel,i,i2) 
     &                     - BC(in,6) * EGmass(iel,i,i4)
           enddo

           do i = 1, nedof
              EGmass(iel,i ,i2) = zero
              EGmass(iel,i2,i ) = zero
              EGmass(iel,i ,i4) = zero
              EGmass(iel,i4,i ) = zero
           enddo
           EGmass(iel,i2,i2) = one
           EGmass(iel,i4,i4) = one
        endif
c     
c.... x2-velocity and x3-velocity
c
        if ( ibits(iBC(in),3,3) .eq. 6 ) then
           do i = 1, nedof
              EGmass(iel,i2,i) = EGmass(iel,i2,i) 
     &                     - BC(in,4) * EGmass(iel,i3,i) 
     &                     - BC(in,6) * EGmass(iel,i4,i)
           enddo
           do i = 1, nedof
              EGmass(iel,i,i2) = EGmass(iel,i,i2) 
     &                     - BC(in,4) * EGmass(iel,i,i3) 
     &                     - BC(in,6) * EGmass(iel,i,i4)
           enddo

           do i = 1, nedof
              EGmass(iel,i ,i3) = zero
              EGmass(iel,i3,i ) = zero
              EGmass(iel,i ,i4) = zero
              EGmass(iel,i4,i ) = zero
           enddo
           EGmass(iel,i3,i3) = one
           EGmass(iel,i4,i4) = one
        endif
c
c.... x1-velocity, x2-velocity, and x3-velocity
c
        if ( ibits(iBC(in),3,3) .eq. 7 ) then
           do i = 1, nedof
              EGmass(iel,i ,i2) = zero
              EGmass(iel,i2,i ) = zero
              EGmass(iel,i ,i3) = zero
              EGmass(iel,i3,i ) = zero
              EGmass(iel,i ,i4) = zero
              EGmass(iel,i4,i ) = zero
           enddo
           EGmass(iel,i2,i2) = one
           EGmass(iel,i3,i3) = one
           EGmass(iel,i4,i4) = one
        endif
c
c.... temperature
c        
	if ( btest(iBC(in),1) ) then
           do i = 1, nedof
              EGmass(iel,i,i5) = zero
              EGmass(iel,i5,i) = zero
           enddo
           EGmass(iel,i5,i5) = one
	endif
c
c	Elaine
c.... scaled plane extraction boundary conditions
c
        if(intpres.eq.1) then
	  if ( btest(iBC(in),11) ) then
            do i = 1, nedof 
              EGmass(iel,i ,i1) = zero
              EGmass(iel,i1,i ) = zero
              EGmass(iel,i ,i2) = zero
              EGmass(iel,i2,i ) = zero
              EGmass(iel,i ,i3) = zero
              EGmass(iel,i3,i ) = zero
              EGmass(iel,i ,i4) = zero
              EGmass(iel,i4,i ) = zero
              EGmass(iel,i,i5) = zero
              EGmass(iel,i5,i) = zero
            enddo
            EGmass(iel,i1,i1) = one
            EGmass(iel,i2,i2) = one
            EGmass(iel,i3,i3) = one
            EGmass(iel,i4,i4) = one
            EGmass(iel,i5,i5) = one
          endif
	else
	  if ( btest(iBC(in),11) ) then
            do i = 1, nedof 
              EGmass(iel,i ,i2) = zero
              EGmass(iel,i2,i ) = zero
              EGmass(iel,i ,i3) = zero
              EGmass(iel,i3,i ) = zero
              EGmass(iel,i ,i4) = zero
              EGmass(iel,i4,i ) = zero
              EGmass(iel,i,i5) = zero
              EGmass(iel,i5,i) = zero
            enddo
            EGmass(iel,i2,i2) = one
            EGmass(iel,i3,i3) = one
            EGmass(iel,i4,i4) = one
            EGmass(iel,i5,i5) = one
          endif
	endif
      
 5000 continue
        
c        
c.... end loop over shape functions (nodes)
c        
      enddo
c
c.... end loop over elements
c
      enddo
c     
c.... return
c
      return
      end
c
c
c
      subroutine bc3LHSSclr (iBC,  iens,  EGmasst )
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of LHS mass matrix for a single 
c element.
c
c input:
c  iBC   (nshg) 	: boundary condition code
c  BC    (nshg,11)      : Dirichlet BC constraint parameters
c  ien   (npro,nshl)	: ien array for this element
c  EGmasst(npro,nshape,nshape) : element consistent mass matrix before BC
c
c output:
c  EGmasst(npro,nshape,nshape): LHS mass matrix after BC is satisfied
c
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Spring 1990. (Modified for general divariant gas)
c----------------------------------------------------------------------
c
        include "common.h"
c
	dimension iBC(nshg),      ien(npro,nshl),
     &		  EGmasst(npro,nshape,nshape) 
c
	integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
	ien=abs(iens)
c
        id = isclr+5
c.... loop over elements
c
        do iel = 1, npro
c
c.... loop over number of shape functions for this element
c
	do inod = 1, nshl
c
c.... set up parameters
c
        in  = ien(iel,inod)
        if (iBC(in) .eq. 0) goto 5000

c
c.... scalar
c
	if ( btest(iBC(in),id) ) then

           do i = 1, nshl
              EGmasst(iel,i,inod) = zero
              EGmasst(iel,inod,i) = zero
           enddo
           EGmasst(iel,inod,inod) = one
	endif

c
      
 5000 continue
        
c        
c.... end loop over shape functions (nodes)
c        
      enddo
c
c.... end loop over elements
c
      enddo
c     
c.... return
c
      return
      end

