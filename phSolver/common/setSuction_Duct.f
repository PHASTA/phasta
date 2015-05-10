c=====================================================================
c Initialize the suction condition when setting the suction dynamically
c based on the main contraction mdot, accompanied by setSuction_Duct2.f
c=====================================================================
        subroutine setSuction_Duct(BC)

        use suctionDuct
        use readarrays
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
        integer i, nn
        real*8 BC(nshg,ndofBC)
        integer isfID

        allocate(suctionf(nshg)) 
        isfID=isetSuction_Duct
        call sfID2np(isfID,nsuctionface,suctionf)
c... above, find suction nodes 

        if(nsuctionface .gt. 0)then !only if current process has nodes on this face
          do i=1,nsuctionface
             nn=suctionf(i) ! map face node index to global node index
             BC(nn,3)=qold(nn,2) 
             BC(nn,4)=qold(nn,3) 
             BC(nn,5)=0 
           enddo
        endif 
c... above, take BC from IC, before set suction vel dynamically in the time step loop

       return
       end  
