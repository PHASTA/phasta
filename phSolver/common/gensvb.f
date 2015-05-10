        subroutine gensvb (ientmp, iBCBtmp, BCBtmp, mattmp,
     &                     ienb,   iBCB,    BCB,    materb)
c
c----------------------------------------------------------------------
c
c  This routine saves the boundary element block.
c
c input:
c  ientmp (npro,nshl)           : boundary nodal connectivity
c  iBCtmp (npro,ndiBCB)         : boundary condition codes
c  BCBtmp (npro,nshlb,ndBCB)    : boundary condition values
c  mattmp (npro)                : material type flag
c
c output:
c  ienb   (npro,nshl)           : boundary nodal connectivity
c  iBCB   (npro,ndiBCB)         : boundary condition codes
c  BCB    (npro,nshlb,ndBCB)    : boundary condition values
c  materb (npro)                : material type flag
c
c
c Zdenek Johan, Winter 1992.
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension   ientmp(npro,nshl),
     &              iBCBtmp(npro,ndiBCB),    BCBtmp(npro,ndBCB)

        dimension   mattmp(npro),           ienb(npro,nshl),
     &              iBCB(npro,ndiBCB),      BCB(npro,nshlb,ndBCB),
     &              materb(npro)
c
c.... generate the boundary element mapping
c
        do i = 1, nshl
          ienb(:,i) = ientmp(:,i)
        enddo
c
c.... save the boundary element data
c
        iBCB   = iBCBtmp
        do i = 1, nenbl ! This is NOT NSHLB as we are just copying the
                        ! piecewise constant data given by NSpre and
                        ! higher order coefficients must be zero
           do j = 1, ndBCB
              BCB(:,i,j)   = BCBtmp(:,j)
           end do
        end do
        do i = nenbl+1, nshlb
           do j = 1, ndBCB
              BCB(:,i,j)   = zero
           end do
        end do

        materb = mattmp
c
c.... return
c
        return
        end
