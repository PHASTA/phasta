c==========================================================
c Initialize the jet mdot when setting the jet mdot
c dynamically based on the contraction mdot, accompanied
c with setBlowing_Duct2.f
c==========================================================

        subroutine setBlowing_Duct(BC,iTurbWall)

        use blowingDuct
        use readarrays ! for qold
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
        integer i, nn
        real*8 BC(nshg,ndofBC)
        real*8 yVel, Temp
        integer isfID
        integer iTurbWall(nshg)

        allocate(jetinletf(nshg))
        isfID=isetBlowing_Duct
        call sfID2np(isfID,njetinlet,jetinletf)
c... above, find jet inlet nodes

        if(njetinlet .gt. 0)then
          do i=1,njetinlet
            nn=jetinletf(i)
c...    from IC (p u v w T nu)
            if (iTurbWall(nn).eq.0)then
               BC(nn,2) = qold(nn,5)  ! set Temp
               BC(nn,3) = 0           ! set xVel=0
               BC(nn,4) = qold(nn,3)  ! set and scale y velocity
               BC(nn,5) = 0           ! set zVel=0
               BC(nn,7) = qold(nn,6)  ! set SA variable
            endif
          enddo
        endif
c... above, take BC from IC, before setting blowing vel dynamically in the time step loop

        return
        end

