      subroutine bc3LHS (iBC,  BC,  iens,  xKebe )
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of LHS mass matrix for all  
c elements in this block.
c
c input:
c  iBC   (nshg)             : boundary condition code
c  BC    (nshg,ndofBC)      : Dirichlet BC constraint parameters
c  ien   (npro,nshape)      : ien array for this element
c  xKebe (npro,9,nshl,nshl) : element consistent mass matrix before BC
c
c output:
c  xKebe (npro,9,nshl,nshl) : LHS mass matrix after BC is satisfied
c
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Spring 1990. (Modified for general divariant gas)
c Ken Jansen, Summer 2000. Incompressible (only needed on xKebe)
c----------------------------------------------------------------------
c
        include "common.h"
c
      dimension iBC(nshg),      ien(npro,nshl),
     & BC(nshg,ndofBC), xKebe(npro,9,nshl,nshl)
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
c        return
        do iel = 1, npro
c
c.... loop over number of shape functions for this element
c
           do inod = 1, nshl
c
c.... set up parameters
c
              in  = abs(ien(iel,inod))
              if (ibits(iBC(in),3,3) .eq. 0) goto 5000 ! NO velocity BC's
              if (ibits(iBC(in),3,3) .eq. 7) goto 5000 ! 3 components ok

c.... 1 or 2 component velocities
c
c
c.... x1-velocity
c
              if ( ibits(iBC(in),3,3) .eq. 1) then
c
! we want to project out the x1 component of the velocity from the tangent  
! matix which is, mathematically, M^e = S^T M^e S. We will do the M^e S
! product first. It has the effect of
! subtracting the column of the block-9 matrix from each column of the block-9
! matrix  that is going to survive (weighted by the coefficient in the
! BC array associated with that row) FOR EACH column  of the
! nshl by nshl matrix FOR EACH element.  THEN the transpose of the
! operation is carried out (replace the word "column" by row
! EVERYWHERE). The following code has been set up so that we only have to 
! give the starting position in each case since we know the block-9 matrix is
! ordered like this   
!  1 2 3
!  4 5 6
!  7 8 9

c
c  adjusting the second column for the eventual removal of the first
c  column of the block-9 submatrix
c
                 irem1=1
                 irem2=irem1+3
                 irem3=irem2+3

                 iadj1=2
                 iadj2=iadj1+3
                 iadj3=iadj2+3
                 do i = 1, nshl
                    xKebe(iel,iadj1,i,inod) = xKebe(iel,iadj1,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem1,i,inod) 
                    xKebe(iel,iadj2,i,inod) = xKebe(iel,iadj2,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem2,i,inod) 
                    xKebe(iel,iadj3,i,inod) = xKebe(iel,iadj3,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem3,i,inod) 

                 enddo
! block status ' denotes colunn 1 projected off.
!  1 2' 3
!  4 5' 6
!  7 8' 9
c
c  adjusting the third column for the eventual removal of the first
c  column of the block-9 submatrix
c
                 iadj1=3
                 iadj2=iadj1+3
                 iadj3=iadj2+3
                 do i = 1, nshl
                    xKebe(iel,iadj1,i,inod) = xKebe(iel,iadj1,i,inod) 
     &                           - BC(in,5) * xKebe(iel,irem1,i,inod) 
                    xKebe(iel,iadj2,i,inod) = xKebe(iel,iadj2,i,inod) 
     &                           - BC(in,5) * xKebe(iel,irem2,i,inod) 
                    xKebe(iel,iadj3,i,inod) = xKebe(iel,iadj3,i,inod) 
     &                           - BC(in,5) * xKebe(iel,irem3,i,inod) 
! block status
!  1 2' 3'
!  4 5' 6'
!  7 8' 9'
                 enddo
                 do i=1,nshl
c
c done with the first  columns_block-9 for columns AND rows of nshl
c
                    xKebe(iel,irem1,i,inod) = zero 
                    xKebe(iel,irem2,i,inod) = zero 
                    xKebe(iel,irem3,i,inod) = zero 


! block status
!  0 2' 3'
!  0 5' 6'
!  0 8' 9'

                 enddo
c
c  now adjust the second row_block-9 for EACH row nshl for EACH element 
c

                 iadj1=4
                 iadj2=iadj1+1
                 iadj3=iadj2+1
                 irem1=1
                 irem2=irem1+1
                 irem3=irem2+1
                 do i = 1, nshl
                    xKebe(iel,iadj1,inod,i) = xKebe(iel,iadj1,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem1,inod,i) 
                    xKebe(iel,iadj2,inod,i) = xKebe(iel,iadj2,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem2,inod,i) 
                    xKebe(iel,iadj3,inod,i) = xKebe(iel,iadj3,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem3,inod,i) 

                 enddo
! block status
!  0 2' 3'
!  0 5'' 6''
!  0 8' 9'


                 iadj1=7
                 iadj2=iadj1+1
                 iadj3=iadj2+1
                 do i = 1, nshl
                    xKebe(iel,iadj1,inod,i) = xKebe(iel,iadj1,inod,i) 
     &                           - BC(in,5) * xKebe(iel,irem1,inod,i) 
                    xKebe(iel,iadj2,inod,i) = xKebe(iel,iadj2,inod,i) 
     &                           - BC(in,5) * xKebe(iel,irem2,inod,i) 
                    xKebe(iel,iadj3,inod,i) = xKebe(iel,iadj3,inod,i) 
     &                           - BC(in,5) * xKebe(iel,irem3,inod,i) 

! block status
!  0 2' 3'
!  0 5'' 6''
!  0 8'' 9''
                 enddo
                 do i=1,nshl

c
c eliminate the first row of block-9 for all rows
c 
                    xKebe(iel,irem1,inod,i) = zero 
                    xKebe(iel,irem2,inod,i) = zero 
                    xKebe(iel,irem3,inod,i) = zero 

                 enddo

! block status
!  0 0   0
!  0 5'' 6''
!  0 8'' 9''

! Be aware that this simple status of the block does not reflect that when 
! we eliminated columns we did it for columns in nshl as well for the given
! inod. Conversely when we eliminated rows in the block we did so for ALL
!  rows in nshl as can be seen by the transpose of i and inod.

                 xKebe(iel,1,inod,inod)=one
! block status
!  1 0   0
!  0 5'' 6''
!  0 8'' 9''
              endif
c
c.... x2-velocity
c
              if ( ibits(iBC(in),3,3) .eq. 2) then
c
! See comment above. Now we are eliminating the 2nd column then row of
! the block-9 matrix
!  1 2 3
!  4 5 6
!  7 8 9
c
c  adjusting the first column for the eventual removal of the second
c  column of the block-9 submatrix
c
                 irem1=2
                 irem2=irem1+3
                 irem3=irem2+3

                 iadj1=1
                 iadj2=iadj1+3
                 iadj3=iadj2+3
                 do i = 1, nshl
                    xKebe(iel,iadj1,i,inod) = xKebe(iel,iadj1,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem1,i,inod) 
                    xKebe(iel,iadj2,i,inod) = xKebe(iel,iadj2,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem2,i,inod) 
                    xKebe(iel,iadj3,i,inod) = xKebe(iel,iadj3,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem3,i,inod) 

                 enddo
c
c  adjusting the third column for the eventual removal of the second
c  column of the block-9 submatrix
c
                 iadj1=3
                 iadj2=iadj1+3
                 iadj3=iadj2+3
                 do i = 1, nshl
                    xKebe(iel,iadj1,i,inod) = xKebe(iel,iadj1,i,inod) 
     &                           - BC(in,5) * xKebe(iel,irem1,i,inod) 
                    xKebe(iel,iadj2,i,inod) = xKebe(iel,iadj2,i,inod) 
     &                           - BC(in,5) * xKebe(iel,irem2,i,inod) 
                    xKebe(iel,iadj3,i,inod) = xKebe(iel,iadj3,i,inod) 
     &                           - BC(in,5) * xKebe(iel,irem3,i,inod) 

                 enddo
                 do i=1,nshl
c
c done with the second  columns_block-9 for columns
c

                    xKebe(iel,irem1,i,inod) = zero 
                    xKebe(iel,irem2,i,inod) = zero 
                    xKebe(iel,irem3,i,inod) = zero 

                 enddo
c
c  now adjust the 1st row_block-9 for EACH row nshl for EACH element 
c

                 iadj1=1
                 iadj2=iadj1+1
                 iadj3=iadj2+1
                 irem1=4
                 irem2=irem1+1
                 irem3=irem2+1
                 do i = 1, nshl
                    xKebe(iel,iadj1,inod,i) = xKebe(iel,iadj1,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem1,inod,i) 
                    xKebe(iel,iadj2,inod,i) = xKebe(iel,iadj2,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem2,inod,i) 
                    xKebe(iel,iadj3,inod,i) = xKebe(iel,iadj3,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem3,inod,i) 

                 enddo
                 iadj1=7
                 iadj2=iadj1+1
                 iadj3=iadj2+1
                 do i = 1, nshl
                    xKebe(iel,iadj1,inod,i) = xKebe(iel,iadj1,inod,i) 
     &                           - BC(in,5) * xKebe(iel,irem1,inod,i) 
                    xKebe(iel,iadj2,inod,i) = xKebe(iel,iadj2,inod,i) 
     &                           - BC(in,5) * xKebe(iel,irem2,inod,i) 
                    xKebe(iel,iadj3,inod,i) = xKebe(iel,iadj3,inod,i) 
     &                           - BC(in,5) * xKebe(iel,irem3,inod,i) 
                 enddo
                 do i=1,nshl

c
c eliminate the second row of block-9 for all rows 
c 
                    xKebe(iel,irem1,inod,i) = zero 
                    xKebe(iel,irem2,inod,i) = zero 
                    xKebe(iel,irem3,inod,i) = zero 
                 enddo
                 xKebe(iel,5,inod,inod)=one
              endif
c
c.... x3-velocity
c
              if ( ibits(iBC(in),3,3) .eq. 4) then
c
! See comment above. Now we are eliminating the 3rd column then row of
! the block-9 matrix
!  1 2 3
!  4 5 6
!  7 8 9
c
c  adjusting the 1st column for the eventual removal of the 3rd
c  column of the block-9 submatrix
c
                 irem1=3
                 irem2=irem1+3
                 irem3=irem2+3

                 iadj1=1
                 iadj2=iadj1+3
                 iadj3=iadj2+3
                 do i = 1, nshl
                    xKebe(iel,iadj1,i,inod) = xKebe(iel,iadj1,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem1,i,inod) 
                    xKebe(iel,iadj2,i,inod) = xKebe(iel,iadj2,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem2,i,inod) 
                    xKebe(iel,iadj3,i,inod) = xKebe(iel,iadj3,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem3,i,inod) 

                 enddo
c
c  adjusting the second column for the eventual removal of the 3rd
c  column of the block-9 submatrix
c
                 iadj1=2
                 iadj2=iadj1+3
                 iadj3=iadj2+3
                 do i = 1, nshl
                    xKebe(iel,iadj1,i,inod) = xKebe(iel,iadj1,i,inod) 
     &                           - BC(in,5) * xKebe(iel,irem1,i,inod) 
                    xKebe(iel,iadj2,i,inod) = xKebe(iel,iadj2,i,inod) 
     &                           - BC(in,5) * xKebe(iel,irem2,i,inod) 
                    xKebe(iel,iadj3,i,inod) = xKebe(iel,iadj3,i,inod) 
     &                           - BC(in,5) * xKebe(iel,irem3,i,inod) 
                 enddo
                 do i=1,nshl

c
c done with the 3rd columns_block-9 for columns 
c

                    xKebe(iel,irem1,i,inod) = zero 
                    xKebe(iel,irem2,i,inod) = zero 
                    xKebe(iel,irem3,i,inod) = zero 

                 enddo
c
c  now adjust the 1st row_block-9 for EACH row nshl for EACH element 
c

                 iadj1=1
                 iadj2=iadj1+1
                 iadj3=iadj2+1
                 irem1=7
                 irem2=irem1+1
                 irem3=irem2+1
                 do i = 1, nshl
                    xKebe(iel,iadj1,inod,i) = xKebe(iel,iadj1,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem1,inod,i) 
                    xKebe(iel,iadj2,inod,i) = xKebe(iel,iadj2,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem2,inod,i) 
                    xKebe(iel,iadj3,inod,i) = xKebe(iel,iadj3,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem3,inod,i) 

                 enddo
                 iadj1=4
                 iadj2=iadj1+1
                 iadj3=iadj2+1
                 do i = 1, nshl
                    xKebe(iel,iadj1,inod,i) = xKebe(iel,iadj1,inod,i) 
     &                           - BC(in,5) * xKebe(iel,irem1,inod,i) 
                    xKebe(iel,iadj2,inod,i) = xKebe(iel,iadj2,inod,i) 
     &                           - BC(in,5) * xKebe(iel,irem2,inod,i) 
                    xKebe(iel,iadj3,inod,i) = xKebe(iel,iadj3,inod,i) 
     &                           - BC(in,5) * xKebe(iel,irem3,inod,i) 

                 enddo
                 do i=1,nshl
                    xKebe(iel,irem1,inod,i) = zero 
                    xKebe(iel,irem2,inod,i) = zero 
                    xKebe(iel,irem3,inod,i) = zero 

                 enddo
                 xKebe(iel,9,inod,inod)=one
              endif
c     
c.... x1-velocity and x2-velocity
c
              if ( ibits(iBC(in),3,3) .eq. 3 ) then
c
! See comment above. Now we are eliminating the 2nd and 1st column then
! same rows of
! the block-9 matrix
!  1 2 3
!  4 5 6
!  7 8 9
c
c  adjusting the 3rd column for the eventual removal of the first and second
c  column of the block-9 submatrix
c
                 irem1=1
                 irem2=irem1+3
                 irem3=irem2+3

                 ire21=2
                 ire22=ire21+3
                 ire23=ire22+3

                 iadj1=3
                 iadj2=iadj1+3
                 iadj3=iadj2+3
                 do i = 1, nshl
                    xKebe(iel,iadj1,i,inod) = xKebe(iel,iadj1,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem1,i,inod) 
     &                           - BC(in,6) * xKebe(iel,ire21,i,inod) 
                    xKebe(iel,iadj2,i,inod) = xKebe(iel,iadj2,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem2,i,inod) 
     &                           - BC(in,6) * xKebe(iel,ire22,i,inod) 
                    xKebe(iel,iadj3,i,inod) = xKebe(iel,iadj3,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem3,i,inod) 
     &                           - BC(in,6) * xKebe(iel,ire23,i,inod) 

! Status of the block-9 matrix
!  1 2 3'
!  4 5 6'
!  7 8 9'
                 enddo
                 do i=1,nshl
c
c done with the first and second columns_block-9 for columns AND rows of nshl
c

                    xKebe(iel,irem1,i,inod) = zero 
                    xKebe(iel,irem2,i,inod) = zero 
                    xKebe(iel,irem3,i,inod) = zero 

                    xKebe(iel,ire21,i,inod) = zero 
                    xKebe(iel,ire22,i,inod) = zero 
                    xKebe(iel,ire23,i,inod) = zero 

! Status of the block-9 matrix
!  0 0 3'
!  0 0 6'
!  0 0 9'

                 enddo
c
c  now adjust the 3rd row_block-9 for EACH row nshl for EACH element 
c

                 iadj1=7
                 iadj2=iadj1+1
                 iadj3=iadj2+1
                 irem1=1
                 irem2=irem1+1
                 irem3=irem2+1
                 ire21=4
                 ire22=ire21+1
                 ire23=ire22+1
                 do i = 1, nshl
                    xKebe(iel,iadj1,inod,i) = xKebe(iel,iadj1,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem1,inod,i) 
     &                           - BC(in,6) * xKebe(iel,ire21,inod,i) 
                    xKebe(iel,iadj2,inod,i) = xKebe(iel,iadj2,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem2,inod,i) 
     &                           - BC(in,6) * xKebe(iel,ire22,inod,i) 
                    xKebe(iel,iadj3,inod,i) = xKebe(iel,iadj3,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem3,inod,i) 
     &                           - BC(in,6) * xKebe(iel,ire23,inod,i) 


! Status of the block-9 matrix
!  0 0 3'
!  0 0 6'
!  0 0 9''
                 enddo
                 do i=1,nshl
                    xKebe(iel,irem1,inod,i) = zero 
                    xKebe(iel,irem2,inod,i) = zero 
                    xKebe(iel,irem3,inod,i) = zero 

                    xKebe(iel,ire21,inod,i) = zero 
                    xKebe(iel,ire22,inod,i) = zero 
                    xKebe(iel,ire23,inod,i) = zero 

! Status of the block-9 matrix
!  0 0 0
!  0 0 0
!  0 0 9''

                 enddo
                 xKebe(iel,1,inod,inod)=one
                 xKebe(iel,5,inod,inod)=one
              endif
c     
c.... x1-velocity and x3-velocity
c
              if ( ibits(iBC(in),3,3) .eq. 5 ) then
c
! See comment above. Now we are eliminating the 1 and 3 column then
! same rows of
! the block-9 matrix
!  1 2 3
!  4 5 6
!  7 8 9
c
c  adjusting the 3rd column for the eventual removal of the first and second
c  column of the block-9 submatrix
c
                 irem1=1
                 irem2=irem1+3
                 irem3=irem2+3

                 ire21=3
                 ire22=ire21+3
                 ire23=ire22+3

                 iadj1=2
                 iadj2=iadj1+3
                 iadj3=iadj2+3
                 do i = 1, nshl
                    xKebe(iel,iadj1,i,inod) = xKebe(iel,iadj1,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem1,i,inod) 
     &                           - BC(in,6) * xKebe(iel,ire21,i,inod) 
                    xKebe(iel,iadj2,i,inod) = xKebe(iel,iadj2,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem2,i,inod) 
     &                           - BC(in,6) * xKebe(iel,ire22,i,inod) 
                    xKebe(iel,iadj3,i,inod) = xKebe(iel,iadj3,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem3,i,inod) 
     &                           - BC(in,6) * xKebe(iel,ire23,i,inod) 

                 enddo
                 do i=1,nshl
c
c done with the first and third columns_block-9 for columns AND rows of nshl
c
                    xKebe(iel,irem1,i,inod) = zero 
                    xKebe(iel,irem2,i,inod) = zero 
                    xKebe(iel,irem3,i,inod) = zero 

                    xKebe(iel,ire21,i,inod) = zero 
                    xKebe(iel,ire22,i,inod) = zero 
                    xKebe(iel,ire23,i,inod) = zero 
                 enddo
c
c  now adjust the 2nd row_block-9 for EACH row nshl for EACH element 
c

                 iadj1=4
                 iadj2=iadj1+1
                 iadj3=iadj2+1
                 irem1=1
                 irem2=irem1+1
                 irem3=irem2+1
                 ire21=7
                 ire22=ire21+1
                 ire23=ire22+1
                 do i = 1, nshl
                    xKebe(iel,iadj1,inod,i) = xKebe(iel,iadj1,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem1,inod,i) 
     &                           - BC(in,6) * xKebe(iel,ire21,inod,i) 
                    xKebe(iel,iadj2,inod,i) = xKebe(iel,iadj2,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem2,inod,i) 
     &                           - BC(in,6) * xKebe(iel,ire22,inod,i) 
                    xKebe(iel,iadj3,inod,i) = xKebe(iel,iadj3,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem3,inod,i) 
     &                           - BC(in,6) * xKebe(iel,ire23,inod,i) 

                 enddo
                 do i=1,nshl
                    xKebe(iel,irem1,inod,i) = zero 
                    xKebe(iel,irem2,inod,i) = zero 
                    xKebe(iel,irem3,inod,i) = zero 

                    xKebe(iel,ire21,inod,i) = zero 
                    xKebe(iel,ire22,inod,i) = zero 
                    xKebe(iel,ire23,inod,i) = zero 

                 enddo
                 xKebe(iel,1,inod,inod)=one
                 xKebe(iel,9,inod,inod)=one
              endif
c     
c.... x2-velocity and x3-velocity
c
              if ( ibits(iBC(in),3,3) .eq. 6 ) then
c
! See comment above. Now we are eliminating the 2nd and 3rd column then
! same rows of
! the block-9 matrix
!  1 2 3
!  4 5 6
!  7 8 9
c
c  adjusting the 3rd column for the eventual removal of the first and second
c  column of the block-9 submatrix
c
                 irem1=2
                 irem2=irem1+3
                 irem3=irem2+3

                 ire21=3
                 ire22=ire21+3
                 ire23=ire22+3

                 iadj1=1
                 iadj2=iadj1+3
                 iadj3=iadj2+3
                 do i = 1, nshl
                    xKebe(iel,iadj1,i,inod) = xKebe(iel,iadj1,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem1,i,inod) 
     &                           - BC(in,6) * xKebe(iel,ire21,i,inod) 
                    xKebe(iel,iadj2,i,inod) = xKebe(iel,iadj2,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem2,i,inod) 
     &                           - BC(in,6) * xKebe(iel,ire22,i,inod) 
                    xKebe(iel,iadj3,i,inod) = xKebe(iel,iadj3,i,inod) 
     &                           - BC(in,4) * xKebe(iel,irem3,i,inod) 
     &                           - BC(in,6) * xKebe(iel,ire23,i,inod) 
                 enddo
                 do i=1,nshl

c
c done with the first and second columns_block-9 for columns AND rows of nshl
c
                    xKebe(iel,irem1,i,inod) = zero 
                    xKebe(iel,irem2,i,inod) = zero 
                    xKebe(iel,irem3,i,inod) = zero 

                    xKebe(iel,ire21,i,inod) = zero 
                    xKebe(iel,ire22,i,inod) = zero 
                    xKebe(iel,ire23,i,inod) = zero 

                 enddo
c
c  now adjust the 3rd row_block-9 for EACH row nshl for EACH element 
c

                 iadj1=1
                 iadj2=iadj1+1
                 iadj3=iadj2+1
                 irem1=4
                 irem2=irem1+1
                 irem3=irem2+1
                 ire21=7
                 ire22=ire21+1
                 ire23=ire22+1
                 do i = 1, nshl
                    xKebe(iel,iadj1,inod,i) = xKebe(iel,iadj1,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem1,inod,i) 
                    xKebe(iel,iadj2,inod,i) = xKebe(iel,iadj2,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem2,inod,i) 
                    xKebe(iel,iadj3,inod,i) = xKebe(iel,iadj3,inod,i) 
     &                           - BC(in,4) * xKebe(iel,irem3,inod,i) 
     &                           - BC(in,6) * xKebe(iel,ire23,inod,i) 

                 enddo
                 do i=1,nshl
                    xKebe(iel,irem1,inod,i) = zero 
                    xKebe(iel,irem2,inod,i) = zero 
                    xKebe(iel,irem3,inod,i) = zero 

c 
                    xKebe(iel,ire21,inod,i) = zero 
                    xKebe(iel,ire22,inod,i) = zero 
                    xKebe(iel,ire23,inod,i) = zero 

                 enddo
                 xKebe(iel,5,inod,inod)=one
                 xKebe(iel,9,inod,inod)=one
              endif
      
 5000         continue
        
c        
c.... end loop over shape functions (nodes)
c        
           enddo
c
c.... end loop over elements
c     
        enddo
c
c These elements should assemble to a matrix with the rows and columns 
c associated with the Dirichlet nodes zeroed out.  Note that BC3 Diag
c
c     
c.... return
c
        return
        end
