c========================================================================
c Dynamically controll the suction level based on main contraction mdot,
c called by itrdrv.f
c========================================================================
        subroutine setSuction_Duct2(x,BC,y)

        use suctionDuct ! nsuctionface suctionf wnorm
        use timedata
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
        integer i, nn
        real*8 xcoor,ycoor,zcoor
        real*8 BC(nshg,ndofBC)
        real*8 x(nshg,nsd)
        real*8 y(nshg,ndof)
        integer isfID
        real rx, rz        
        real*8 stripleng, Area
        real*8 xminm,xmaxm,zminm,zmaxm
        real*8 contraPres, contraTemp, contraXVel, contraRho, contraMdot
        integer contraProbeNo
        real*8 nVel, pres, rho, Temp, mdot  

        contraProbeNo=5 ! in xyzts.dat the 5th probe is at contraction inlet
        if(myrank.eq.0)then ! only rank0 process can acess all varts data 
          contraPres=varts(contraProbeNo,1)
          contraTemp=varts(contraProbeNo,5)
          contraXVel=varts(contraProbeNo,2)
          contraArea=1.17475**2 ! m^2
          contraRho=contraPres/287/contraTemp
          contraMdot=contraRho*contraArea*contraXVel !kg/s
          write(*,*)'contraMdot is', contraMdot
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_BCAST(contraMdot,1,MPI_REAL8,master,
     &                 MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if(nsuctionface .gt. 0)then !only if current process has nodes on this face
          mdot=1.0/100*contraMdot   ! 1% of the main mdot     
          stripleng=4*1e-3
          xminm=-0.252*0.0254
          xmaxm=2.486*0.0254
          zminm=-4.5*0.0254/2.0
          zmaxm=4.5*0.0254/2.0
          Area = 0.00965813 ! m^2  
          do i=1,nsuctionface
            nn=suctionf(i) ! map face node index to global node index
            xcoor=x(nn,1)
            zcoor=x(nn,3)
            if(xcoor.ge.xminm.and.xcoor.le.xmaxm)then !only if in the suction patch

             if((xcoor-xminm)/stripleng.lt.1.0)then
                 rx=(xcoor-xminm)/stripleng                          
             elseif((xmaxm-xcoor)/stripleng.lt.1.0)then
                 rx=(xmaxm-xcoor)/stripleng
             else
                 rx=1.0
             endif              
             rx=max(0.0,rx)
             if((zcoor-zminm)/stripleng.lt.1.0)then
                 rz=(zcoor-zminm)/stripleng
             elseif((zmaxm-zcoor)/stripleng.lt.1.0)then
                 rz=(zmaxm-zcoor)/stripleng
             else
                 rz=1.0
             endif
             rz=max(0.0,rz)

             pres=y(nn,4) ! the quantity updating each time step
             Temp=317 ! wall temp 
             rho=pres/287/Temp
             nVel=(mdot/rho/Area)*rx*rz ! normal suction velocity magnitude
             BC(nn,3)=nVel*wnorm(nn,1) ! set xVel 
             BC(nn,4)=nVel*wnorm(nn,2) ! set yVel
            endif ! only if in the suction patch
          enddo

        endif  ! only if current process contains suction face nodes

        end  
