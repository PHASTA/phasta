c..............................................................................
c  set jet inlet velocity, temperature, and SA variable with fixed value on the geometry with truncted jet path
        subroutine setBlowing_Duct3(x,BC,iTurbWall)

        use blowingDuct
        use readarrays ! for qold
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
        integer i, nn
        real*8 BC(nshg,ndofBC)
        real*8 x(nshg,nsd)

        real Temp,veloMag,xcoor,ycoor,zcoor,rz,ratio_s
        real*8 x_begin,y_begin,x_end,y_end,jetnx,jetny
        real*8 edge_length, edge_s
        integer isfID
        integer iTurbWall(nshg)     
   
        allocate(jetinletf(nshg))
        isfID=isetBlowing_Duct
        call sfID2np(isfID,njetinlet,jetinletf)
c... above, find jet inlet nodes
        if(myrank.eq.0)write(*,*)'Blowing Vel:',BlowingVelDuct

c.. jet inlet plane geometric parameters
            x_begin = -0.0100093419234682
            y_begin = -0.0456888187657512
            x_end   = -0.00982113237765022
            y_end   = -0.04627072277068951
            jetnx   = 0.951470142380476
            jetny   = 0.307741073239299
            edge_length =
     &         sqrt(abs((x_end-x_begin)**2+(y_end-y_begin)**2))

        if(njetinlet .gt. 0)then
          do i=1,njetinlet
            nn=jetinletf(i)
            if(iTurbWall(nn).eq.0)then
c---------------------------------------------------
              if(iDuctgeometryType.eq.6)then
                xcoor=x(nn,1)
                ycoor=x(nn,2)
                zcoor=x(nn,3)

c----Geometry6 Condajet                     
c the jet inlet velocity profile is as a plug profile 
c in z direction and a 1-(s-1)^10 profile along the 
c inlet projection in x-y plane, in which s varies 
c from 0 to 2 as the point on the jet inlet moving from one end to the other end
 
                if(zcoor<(-0.05715+1.0e-3))then
                 rz=(zcoor+0.05715)/1.0e-3
                elseif(zcoor>(0.05715-1.0e-3))then
                 rz=(0.05715-zcoor)/1.0e-3
                else
                 rz=1.0
                endif
c... 
                edge_s=2* 
     $    sqrt(abs((xcoor-x_begin)**2+(ycoor-y_begin)**2))/edge_length

                ratio_s=1-(edge_s-1)**20
              
                Temp=317.0
                veloMag=200.0*max(0.0,min(1.0,ratio_s))
     
                BC(nn,2) = Temp   ! Temp
                BC(nn,3) = veloMag*jetnx   ! set x velocity
                BC(nn,4) = veloMag*jetny   ! set y velocity
                BC(nn,5) = 0 ! set zVel=0

                BC(nn,2) = qold(nn,5)   ! Temp
                BC(nn,3) = 0   ! set x velocity
                BC(nn,4) = qold(nn,3)   ! set y velocity
                BC(nn,5) = 0 ! set zVel=0
                BC(nn,7) = qold(nn,6) 
              endif

c-----Geometry8--------------------------------------------
              if (iDuctgeometryType .eq. 8) then
                 BC(nn,2) = 317 ! Temp
                 BC(nn,3) = BlowingVelDuct*0.882947592858928   ! set x velocity
                 BC(nn,4) = BlowingVelDuct*0.46947156278589    ! set y velocity
                 BC(nn,5) = 0   ! set zVel=0
                 BC(nn,7) = 1.825e-5   ! set varSA = 0
              endif
            endif  
          enddo
        endif

        return
        end

