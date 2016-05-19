      module solvedata

      integer nsolflow,npermDims, nTmpDims, nPermDimsS, nTmpDimsS

      real*8,  allocatable :: lhsP(:,:), lhsK(:,:), lhsS(:,:)
      real*8,  allocatable :: aperm(:,:), atemp(:,:)
      real*8,  allocatable :: apermS(:,:,:), atempS(:,:)

      end module

c-----------------------------------------------------------------------
c allocate the arrays
c-----------------------------------------------------------------------
      subroutine aSDf 
      use solvedata
      include "common.h"
      allocate (lhsP(4,nnz_tot))
      allocate (lhsK(9,nnz_tot))
      if(leslib.eq.1) then
        allocate (aperm(nshg,nPermDims))
        allocate (atemp(nshg,nTmpDims))
      endif
      return
      end
   
      subroutine aSDs(nsclrsol) 
      use solvedata
      include "common.h"
      allocate (lhsS(nnz_tot,nsclrsol))
      if(leslib.eq.1) then
        allocate (apermS(nshg,nPermDimsS,nsclrsol))
        allocate (atempS(nshg,nTmpDimsS))
      endif
      return
      end
c-----------------------------------------------------------------------
c delete the arrays
c-----------------------------------------------------------------------
      subroutine dSDf 
      use solvedata
      include "common.h"
      deallocate (lhsP)
      deallocate (lhsK)
      if(leslib.eq.1) then
        deallocate (aperm)
        deallocate (atemp)
      endif
      return
      end

      subroutine dSDs 
      use solvedata
      include "common.h"
      deallocate (lhsS)
      if(leslib.eq.1) then
        deallocate (apermS)
        deallocate (atempS)
      endif
      return
      end
