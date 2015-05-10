c-----------------------------------------------------------------------
c
c     module for time averaged statistics (conservative projection).
c
c-----------------------------------------------------------------------
      module stats
      
      integer nResDims, nSolDims, nLhsDims, nTimeStep, stsResFlg
      integer stsCompFreq, stsWriteFreq, stsResetFreq, step1,
     &        stsType
      
      real*8, allocatable :: stsVec(:,:)
      
      real*8, allocatable :: stsReg(:)
      real*8, allocatable :: stsMInv(:,:)
      real*8, allocatable :: stsB(:,:)
      real*8, allocatable :: stsDInv(:,:)
      real*8, allocatable :: stsCInv(:,:)
      
      real*8, allocatable :: stsPres(:), stsPresSqr(:), stsVel(:,:),
     &                       stsVelSqr(:,:), stsVelReg(:,:),
     &                       stsStress(:,:)

      end module
      
c-----------------------------------------------------------------------
c     create the new statistics arrays
c-----------------------------------------------------------------------
      subroutine initStats(x,   iBC,    iper,   ilwork)
      
      use stats
      include "common.h"

      real*8  x(numnp,3)
      integer ilwork(nlwork), iper(nshg), iBC(nshg)
      
      write(*,*) 'Stats are not developed for compressible code'
      return
      end
