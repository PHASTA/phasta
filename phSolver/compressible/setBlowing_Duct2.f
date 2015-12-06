c================================================================
c Set jet inlet BCs based on contraction mdot dynamically,
c called by itrdrv.f 
c ===============================================================

        subroutine setBlowing_Duct2(x,BC,y,iTurbWall,istp)

        use blowingDuct ! njetinlet, jetinletf 
        use timedataC ! varts data
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
        integer istp
        integer i, nn
        real*8 xcoor,ycoor,zcoor
        real*8 BC(nshg,ndofBC)
        real*8 x(nshg,nsd)
        real*8 y(nshg,ndof)
        real*8 yVel, Temp, pres, mdot, rho, Area
        real*8 xminm,xmaxm,zminm,zmaxm,xlm,zlm
        real rx,rz
        integer isfID
        real*8 contraPres, contraTemp, contraXVel, contraRho, 
     &         contraArea, contraMdot  
        real*8 blowingPres,blowingTemp,blowingRho,blowingArea,
     &         MdotRatio,blowingMdot,blowVel
        integer contraProbeNo,blowingProbeNo        
        integer iTurbWall(nshg)

c--- The following operations must be fulfilled on each proces because even if 
c--- there is no jet mouth on that process, it can be the host process that contains
c--- all the varts data

        if(myrank.eq.0)then ! only rank0 process can acess all varts data 
          contraProbeNo=1 ! in xyzts.dat the FIRST probe is at contraction inlet
          blowingProbeNo=2
          contraArea=1.17475**2 ! meter^2
          blowingArea=0.00237064  !meter^2
 
          contraPres=varts(contraProbeNo,1)
          contraTemp=varts(contraProbeNo,5)
          contraXVel=varts(contraProbeNo,2)
          contraRho =contraPres/(287*contraTemp)
          contraMdot=contraRho*contraArea*contraXVel !kg/s

          blowingPres=varts(blowingProbeNo,1)
          blowingTemp=varts(blowingProbeNo,5)
          blowingRho = blowingPres/(287*blowingTemp)

          MdotRatio =  (BlowingIniMdotDuct) +
     &   real(istp-1)/real(nBlowingStepsDuct-1)
     &   *(BlowingFnlMdotDuct-BlowingIniMdotDuct)

          blowingMdot= (MdotRatio/100.0)*contraMdot   ! The absolute value of mdot
          blowingVel = blowingMdot/(blowingArea*blowingRho)
          write(*,*)'Blowing Velocity:', blowingVel
        endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_BCAST(blowingVel,1,MPI_REAL8,master,
     &                 MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


c--- Only if current process contains jet inlet nodes
        if(njetinlet .gt. 0)then  
          do i=1,njetinlet
            nn=jetinletf(i)
            if(iTurbWall(nn).eq.0)then
c            xcoor=x(nn,1)
c            ycoor=x(nn,2)
c            zcoor=x(nn,3)
c            rx=(xcoor-xminm)*(xmaxm-xcoor)/(xlm/2)**2
c            rz=(zcoor-zminm)*(zmaxm-zcoor)/(zlm/2)**2
c            rx=max(0.0,rx)
c            rz=max(0.0,rz)
c            pres = y(nn,4) ! the quantity updating each time step
c           Temp = y(nn,5) 
c           rho = pres/(287.0*Temp)             
c              yVel = (blowingMdot/rho/Area)*rx*rz 
c............................
              BC(nn,2) = 317          ! Temp
              BC(nn,4) = blowingVel   ! set and scale y velocity
              BC(nn,3) = 0
              BC(nn,5) = 0
              BC(nn,7) = 1.825e-5
            endif
          enddo

        endif   ! only if current process contains jet inlet nodes

        return
        end

