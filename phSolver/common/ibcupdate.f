      subroutine iBCupdate(iBCpart,   ienb,    iBCB)
c
c---------------------------------------------------------------------
c  This is the subroutine update iBC for deformable wall case
c  where iBCB(:,1) = 16: Turbulence wall
c  make iBC(:) = iBC(:) + 8192
c---------------------------------------------------------------------
c
      include "common.h"
      
      dimension iBCpart(nshg),        iBCB(npro,ndiBCB),
     &          ienb(npro,nshl)

      do iel = 1, npro
         if(btest(iBCB(iel,1),4)) then ! turbulence wall
            do inode = 1, nenbl
               iglobal = ienb(iel,inode)
               iBCpart(iglobal) = 8192
            enddo
         endif
      enddo

      return
      end
