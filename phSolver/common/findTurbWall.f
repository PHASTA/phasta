      subroutine findTurbWall(iTurbWall)

ccc
c
ccc
      use pointer_data
      include "common.h"
      include "mpif.h"
      integer, allocatable :: ienb(:)
      integer iTurbWall(nshg)

c--- Exclude all the nodes on turbulence wall----c
      iTurbWall(:)=0
      do iblk=1,nelblb ! loop element block
         npro=lcblkb(1,iblk+1)-lcblkb(1,iblk)
         nenbl=lcblkb(6,iblk)
         nshl=lcblkb(9,iblk)
         allocate(ienb(nshl))
         do i=1,npro   ! loop element
            iBCB1=miBCB(iblk)%p(i,1)
            if(.not.btest(iBCB1,4))cycle
            ienb(1:nshl)=mienb(iblk)%p(i,1:nshl)
            do j=1,nenbl  ! loop elemental vertex
               nn=ienb(j)
               iTurbWall(nn)=1   
            enddo  ! end loop elemental vertex
         enddo ! end loop element
         deallocate(ienb)
      enddo ! end loop element block 

      return
      end 
