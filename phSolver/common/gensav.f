        subroutine gensav (ientmp, mattmp, ien,    mater)
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
        include "common.h"
c
        dimension   ientmp(npro,nshl),
     &              mattmp(npro),           ien(npro,nshl),
     &              mater(npro)
c
c.... save the element data
c
        do i = 1, nshl
          ien(:,i) = ientmp(:,i)
        enddo
c
        mater = mattmp
c
c.... end
c
        return
        end
