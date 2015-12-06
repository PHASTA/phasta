        subroutine gensav (ientmp, mattmp, ien,    ienG, mater)
c
c----------------------------------------------------------------------
c
c  This routine saves the element block data.
c
c input:
c  ientmp (npro,nshl)   : nodal connectivity 
c  mattmp (npro)        : material type flag
c
c output:
c  ien    (npro,nshl)   : nodal connectivity
c  mater  (npro)        : material type flag
c
c
c Zdenek Johan, Winter 1992.
c----------------------------------------------------------------------
c
        use readarrays
        use fncorpmod
        include "common.h"
c
        dimension   ientmp(npro,nshl),
     &              mattmp(npro),           ien(npro,nshl),
     &              mater(npro)
       integer*8    ienG(npro,nshl)
c
c.... save the element data
c
        do i = 1, nshl
          ien(:,i) = ientmp(:,i)
        enddo
        if(usingpetsc.eq.1) then
          do i = 1, nshl
            if(numpe .ne. 1) then
              ienG(:,i) = fncorp(abs(ientmp(:,i)))
            else
              ienG(:,i) = abs(ientmp(:,i)) 
            endif
          enddo
        endif
c
        mater = mattmp
c
c.... end
c
        return
        end
