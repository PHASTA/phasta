c========================================================================
c Dynamically controll the suction level based on main contraction mdot,
c called by itrdrv.f
c========================================================================
        subroutine setSuction_Duct3(x,BC,y, ilwork)

        use wallData ! wnorm
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        integer i, j, k, nn, nL
        integer ilwork(nlwork)
        integer nSuctionNode
        integer, allocatable :: suctionNodeMap(:)
!        logical inPatch
        real*8 xcoor,ycoor,zcoor
        real*8 BC(nshg,ndofBC)
        real*8 BC3(numnp, 5)                      !used for hack to get suction right on part boundaries
        real*8 x(nshg,nsd), y(nshg,ndof) 
!       real*8 xmin_t, xmax_t, xmin_b, xmax_b     !ramp geometry controls
!       real*8 xmin_su, xmax_su, xmin_sl, xmax_sl !side wall patch limits
!       real*8 t_upper                            !upper thickness of suction area
!       real*8 ymin_su, ymax_sl 
        real*8 nvel_s, nvel_t, nvel_b             !normal velocity on the side walls, top, and bottom
        
!        real*8 delX, delY, X1, Y1     !Variables used to map upper surface or trumpet
!        real*8 Xline(30,2)
!        integer suctionIDs(4)

!        real*8 wallBCout(nshg, 6)

        if(myrank.eq.0)write(*,*)'Setting side suctions'
 
        nvel_b     = suctionVbottom      !normal velocity at the bottom patch
        nvel_sl    = suctionVside_lower  !normal velocity of lower side wall patches
        nvel_su    = suctionVside_upper  !normal velocity of upper side wall patches
        nvel_t     = suctionVtop         !normal velocity at the top patch
       
        allocate(suctionNodeMap(nshg))
        call sfID2np(isetSuctionID_Duct, nSuctionNode, suctionNodeMap)  !to figure out whether it's the top, bottom, or side surface.

        smallNum=1.0e-5
        xsideBotMax=0.0428625+smallNum
        xsideTopMax=0.130175+smallNum 

        do i = 1, nSuctionNode
          nn = suctionNodeMap(i) ! map face node index to global node
          xcoor = x(nn,1)
          ycoor = x(nn,2)
          !zcoor = x(nn,3)
          
          !wnormx = wnorm(nn, 1)
          wnormy = wnorm(nn, 2)
          wnormz = wnorm(nn, 3)

          !Test if the point lies in one of the suction patches 
       
          if((abs(nvel_sl) >  0  .or. abs(nvel_su) > 0) .and.  !if either upper or lower side walls are on...
     &       (abs(wnormz) >=  0.999)) then                     ! and the node normal is in the z direction
            if((ycoor > 0) .and. (xcoor .le. xsideTopMax)) then !hack to prevent edges from turning on. 
              BC(nn, 3) = nvel_su*wnorm(nn, 1) ! set xVel 
              BC(nn, 4) = nvel_su*wnorm(nn, 2) ! set yVel
              BC(nn, 5) = nvel_su*wnorm(nn, 3) ! set zVel
            elseif (xcoor .le. xsideBotMax) then
              BC(nn, 3) = nvel_sl*wnorm(nn, 1) ! set xVel 
              BC(nn, 4) = nvel_sl*wnorm(nn, 2) ! set yVel
              BC(nn, 5) = nvel_sl*wnorm(nn, 3) ! set zVel
            endif
          else if(abs(wnormz) > 0.1) then      ! side walls should have a larger znormal than this - it's ambiguous what to do so do nothing
            cycle  
          else if(abs(nvel_t) > 0      .and. !top suction patch
     &                 wnormy > 0.001) then  !face is pointing up (in the +y direction)
            BC(nn, 3) = nvel_t*wnorm(nn, 1) ! set xVel using wall normal
            BC(nn, 4) = nvel_t*wnorm(nn, 2) ! set yVel
            BC(nn, 5) = nvel_t*wnorm(nn, 3) ! set zVel

          else if(abs(nvel_b) > 0       .and. !bottom suction patches
     &                wnormy  < -0.001) then  !face is pointing down (in the -y direction)
            BC(nn, 3) = nvel_b*wnorm(nn, 1) ! set xVel using wall normal
            BC(nn, 4) = nvel_b*wnorm(nn, 2) ! set yVel
            BC(nn, 5) = nvel_b*wnorm(nn, 3) ! set zVel
          endif

        enddo

        BC3=zero
        if(numpe.gt.1) then
            BC3(:,1:3)=BC(:,3:5)
            do i=1,nshg
                bc3mag=BC3(i,1)**2+BC3(i,2)**2+BC3(i,3)**2
                if(bc3mag.gt.0) BC3(i,4)=one
            enddo
            call commu(BC3,ilwork,5,'in ')
! This accumulated and failed                 call commu(BC(:,3:5),ilwork,3,'in ')
            do i=1,nshg
                bcmag=BC(i,3)**2+BC(i,4)**2+BC(i,5)**2
                if((bcmag.eq.0).and.(BC3(i,4).ne.0)) then !node is slave and has not been correctly set.
                     BC(i,3:5)=BC3(i,1:3)/BC3(i,4)  !division by BC3 is necessary to account for more than two blocks sharing the same node
                endif
            enddo
       endif

        !Debugging loop
!        if(nSuctionNode.gt.0) then
!          icount=0
!          do i=1,nSuctionNode
!             nn=suctionNodeMap(i)
!             dz=dabs(dabs(x(nn,3))-0.0762)
!             if(dz.lt.0.0001) then
!               icount=icount+1
!               if(icount.eq.1)
!     & write(myrank+1000,*) 'writing wnorm and BC setSuctionNG3'
!               write(myrank+1000,1234) nn,i,
!     &           x(nn,1), x(nn,2), x(nn,3),
!     &           wnorm(nn,1),wnorm(nn,2),wnorm(nn,3),
!     &           BC(nn,3),BC(nn,4),BC(nn,5)
!             endif
!          enddo
!         if(icount.gt.0) close(myrank+1000)
!        endif
!1234           format(i5,2x,i5,9(2x,e14.7))
        
      end 






