        subroutine geniBC (iBC)
c
c----------------------------------------------------------------------
c This routine reads the boundary condition codes.
c
c output: 
c  iBC   (nshg)        : Boundary Condition code
c
c         = 1 * iBC_1 + 2 * iBC_2 + 4 * iBC_3
c              density   temperature   pressure
c
c    if nsd = 3:
c
c        +  8 * iBC_4 +  16 * iBC_5 +  32 * iBC_6
c           x1-velocity   x2-velocity   x3-velocity
c
c        + 64 * iBC_7 + 128 * iBC_8 + 256 * iBC_9 + 512 * iBC_10
c          sclr1         sclr2        sclr3         sclr4
c
c        + 1024 * iBC_11  + 2048* iBC_12 + 4096* iBC_13 + 8192* iBC_14
c          perioidicity     spebc          axisym         deformwall
c
c  nBC   (nshg)        : Boundary Condition mapping array
c
c
c Farzin Shakib, Winter 1986.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
c
        use readarrays          ! used to access iBCtmp
        use pointer_data
        include "common.h"
c
c Arrays in the following 1 line are now dimensioned in readnblk
c        dimension iBCtmp(numpbc)
c
        dimension iBC(nshg)
        dimension itemp(6)
        integer, allocatable :: iBCpart(:)
c
c.... set the iBC array
c
        iBC = 0
c
        if(numpbc.eq.0) goto 9999  ! sometimes there are no BC's on a partition
        where (nBC(:) .ne. 0) iBC(:) = iBCtmp(nBC(:))
c
c.... echo the input iBC array only if other than zero
c
        if (necho .lt. 3) then
          nn = 0
          do n = 1, nshg
            if (nBC(n) .ne. 0) then
              nb = nBC(n)
              nn = nn + 1
              if (mod(nn,50).eq.1) write(iecho,1000)ititle,(j,j=1,ndof)
              itemp(   1) = mod(iBCtmp(nb)   ,2) - mod(iBCtmp(nb)/ 4,2)
              itemp(   2) = mod(iBCtmp(nb)/ 8,2)
              itemp(   3) = mod(iBCtmp(nb)/16,2)
              itemp(   4) = mod(iBCtmp(nb)/32,2)
              itemp(ndof) = mod(iBCtmp(nb)/ 2,2)
              write(iecho,1100) n,(itemp(i),i=1,ndof)
            endif
          enddo
        endif

c
c.... for deformable wall case update iBC from iBCB information
c

9999   if(ideformwall.eq.1) then
          allocate (iBCpart(nshg))
          iBCpart = 0
          do iblk = 1, nelblb
             iel    = lcblkb(1,iblk)
             iorder = lcblkb(4,iblk)
             nenl   = lcblkb(5,iblk) ! no. of vertices per element
             nenbl  = lcblkb(6,iblk) ! no. of vertices per bdry. face
             nshl   = lcblkb(9,iblk)
             nshlb  = lcblkb(10,iblk)
             npro   = lcblkb(1,iblk+1) - iel 
             call iBCupdate(iBCpart,  mienb(iblk)%p,   miBCB(iblk)%p)
          enddo
          iBC = iBC + iBCpart
          deallocate(iBCpart)
       endif

        deallocate(iBCtmp)

       
          
c
c.... return
c
        return
c
c.... end of file error handling
c
999     call error ('geniBC  ','end file',ibndc)
c
1000    format(a80,//,
     &  ' N o d a l   B o u n d a r y   C o n d i t i o n   C o d e',//,
     &  '    Node   ',13x,6('dof',i1,:,6x))
1100    format(2x,i5,10x,5i10)
c
        end
