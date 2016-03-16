      module solvedata

      integer nsolflow,npermDims, nTmpDims, nPermDimsS, nTmpDimsS

      real*8,  allocatable :: lhsP(:,:), lhsK(:,:), lhsS(:,:)
      real*8,  allocatable :: aperm(:,:), atemp(:,:)
      real*8,  allocatable :: apermS(:,:,:), atempS(:,:)

      end module

c-----------------------------------------------------------------------
c allocate the arrays
c-----------------------------------------------------------------------


      subroutine aSD 

      use solvedata
      include "common.h"
      if(nsolflow.eq.1) then
         allocate (lhsP(4,nnz_tot))
         allocate (lhsK(9,nnz_tot))
         if(leslib.eq.1) then
           allocate (aperm(nshg,nPermDims))
           allocate (atemp(nshg,nTmpDims))
         endif
      endif
 
      return
      end
c-----------------------------------------------------------------------
c delete the arrays
c-----------------------------------------------------------------------

      
      subroutine dSD 

      use solvedata

      deallocate (lhsP)

      return
      end
