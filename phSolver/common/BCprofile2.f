c..............................................................................

        subroutine SetUniOutPres(BC)
          
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        integer BCfaceNode(nshg)
        integer nBCfaceNode
        integer isfID
        real*8 BC(nshg,ndofBC)

        if(myrank.eq.0)write(*,*)'Outlet pressure:',outPres1
        isfID=isetOutPres
        call sfID2np(isfID,nBCfaceNode,BCfaceNode)  
        if(nBCfaceNode .gt. 0)then
           do i=1,nBCfaceNode
             nn=BCfaceNode(i)
             BC(nn,1)=outPres1
           enddo
        endif
         
        return 
        end

c..............................................................................
        subroutine setInlet_Duct(x,BC,iTurbWall)

        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        real*8 BC(nshg,ndofBC)
        real*8 x(nshg,nsd)
        integer i,nn
        real*8 xcoor,ycoor,zcoor 
        real ry,rz
        real*8 Temp,xVel
        integer BCfaceNode(nshg)
        integer nBCfaceNode
        integer isfID
        integer iTurbWall(nshg)

        if(myrank.eq.0)write(*,*)'Inlet surf:',isetInlet_Duct 
       
        isfID=isetInlet_Duct
        call sfID2np(isfID,nBCfaceNode,BCfaceNode)
        if(nBCfaceNode .gt. 0)then
          do i=1,nBCfaceNode
            nn=BCfaceNode(i)
            if(iTurbWall(nn).eq.0)then

              xcoor=x(nn,1)
              ycoor=x(nn,2)
              zcoor=x(nn,3) 
c...........................    

c Contraction Square Length is 46.25 inch, 0.587375 = 46.25/2*0.0254
c Inlet Vel is computed based on throat Mach 0.43 
c Inlet Temp is 330K, wall temp is 317K
                   
              if(ycoor<(-0.587375+1.0e-4))then
               ry=(ycoor+0.587375)/1.0e-4
              elseif(ycoor>(0.587375-1.0e-4))then
               ry=(0.587375-ycoor)/1.0e-4
              else
               ry=1.0
              endif

              if(zcoor<(-0.587375+1.0e-4))then
               rz=(zcoor+0.587375)/1.0e-4
              elseif(zcoor>(0.587375-1.0e-4))then
               rz=(0.587375-zcoor)/1.0e-4
              else
               rz=1.0
              endif
                 
              ry=max(0.0,ry)
              rz=max(0.0,rz)
          
              xVel = 1.513158*ry*rz
c              xVel = 0.1*ry*rz
              Temp = 317.0+13.0*ry*rz
c..........................................             
              BC(nn,2) = Temp   !Temp 
              BC(nn,3) = xVel   ! set and scale x velocity
              BC(nn,4) = 0
              BC(nn,5) = 0
            endif     
          enddo 
         endif 
       
        return
        end 

