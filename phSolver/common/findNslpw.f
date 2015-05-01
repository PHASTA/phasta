c==== This subroutine is used to find no-slip wall defined in Simulation Maker,
c     with SurfID ranging from 601 to 609 
c     About velocity BC priority
c     no-slip wall > slip wall > inlet/outlet
c     Array iNslpw(nshg) is an array defined at all the nodal points. 
c     It is a lable to identify it current point is on a no-slip wall.
c     Similar arrays are iTurbWall, islpw
c=============================================================================

      subroutine findnslpw(iNslpw)

      use pointer_data ! for miBCB and mienb
      include "common.h"
      include "mpif.h"
      integer, allocatable :: ienb(:)
      integer iNslpw(nshg)

c--- Exclude all the nodes on no-slip wall----
      iNslpw(:)=0
      do iblk=1,nelblb ! loop over boundary element blocks
         npro=lcblkb(1,iblk+1)-lcblkb(1,iblk)
         nenbl=lcblkb(6,iblk)
         nshl=lcblkb(9,iblk) 
         allocate(ienb(nshl))
         do i=1,npro   ! loop over boundary elements
            isfID=miBCB(iblk)%p(i,2) ! SurfID is stored at miBCB(2)
            if(isfID.gt.609.or.isfID.lt.601)cycle ! restart the loop over boundary elements

            ienb(1:nshl)=mienb(iblk)%p(i,1:nshl)

c           Once find the qualified boundary cell, then all the nodes belonging to that cell  
c           are assigned the label iNslpw=1

            do j=1,nenbl  ! loop elemental vertex
               nn=ienb(j)  ! ienb maps the boundary elemental nodal index back to the global nodal index
               iNslpw(nn)=1   
            enddo  ! end loop over elemental vertex
         enddo ! end loop over elements

         deallocate(ienb)
      enddo ! end loop element block 

      return
      end 
