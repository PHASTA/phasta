        subroutine INIprofile(BC,y,x)
        
         include "common.h"
         real*8 BC(nshg,ndofBC)
         real*8 y(nshg,ndof)
         real*8 x(nshg,nsd)

         if(iI2Binlet .gt. 0)then ! useless for geometry with contraction
              call TakeBCfromIC_ScaleInlVelX(BC,y,x)
         endif           

         if(isetInitial .gt.0)then ! obscure
              call setInitial(x,y)
         endif

        return
        end
c..............................................................................
        subroutine TakeBCfromIC_ScaleInlVelX(BC,y,x)

        use BCsfIDmap
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        integer inlf_tmp(nshg)
        real*8 BC(nshg,ndofBC)
        real*8 y(nshg,ndof)
        real*8 x(nshg,nsd)
        real*8 u_max,u_max_all

         call sfID2np(iI2Binlet,ninlet,inlf_tmp)
         if(ninlet.gt.0) then  ! only for process including inlet nodes
            allocate(inlf(ninlet))
            inlf(1:ninlet)=inlf_tmp(1:ninlet) ! generate map mapping inlet nodes to nshg nodes
            do i=1,ninlet
               nn=inlf(i)
               if(y(nn,1).gt.u_max)u_max=y(nn,1) 
            enddo
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE( u_max, u_max_all, 1,
     &        MPI_REAL8, MPI_MAX, MPI_COMM_WORLD,ierr)
         if(myrank.eq.0)then
           write(*,*)''
           write(*,556)'u_max_all=',u_max_all ! find the largest x velocity at inlet
         endif
c apply IC to BC, must be run before applying BC to IC later in genini
         if(ninlet.gt.0) then
            BC(inlf(:),2) = y(inlf(:),5) ! set temp
            BC(inlf(:),3) = y(inlf(:),1)/u_max_all*inletVelX  ! set and scale x velocity
            BC(inlf(:),4) = 0             ! set y velocity
            BC(inlf(:),5) = 0             ! set z velocity
            if(nsclr.eq.1)then 
              BC(inlf(:),7) = y(inlf(:),6)  ! set scalar_1
            endif
         endif
555          format(I3,F10.3)
556          format(A,F10.3)
       
        return
        end 
         
c..............................................................................
c... this subroutine gives the formular of initial conditions, user may want to
c... change this file and recompile the code, just like other commercial code which
c... supports User Defined Function, you need to recompile and relink the code  
        subroutine setInitial(x,y)

        use turbSA   ! this gives us d2wall(1:nshg) nodal distance to the closest wall
        use BCsfIDmap
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        real*8 x(nshg,nsd)
        real*8 y(nshg,ndof)         
        if(isetinitial.eq.1) then

c... user gives the formular of initial conditions from here       
           y(:,1)= xvel_ini 
           y(:,4)= pres_ini
           y(:,2)=yvel_ini
           y(:,3)=zvel_ini
           y(:,5)=temp_ini
           if(nsclr.eq.1)then
              y(:,6)=evis_ini
           endif
!moving previously coded spatially varying pressure to 10 so that 1 is simple
        elseif (isetInitial.eq.10) then   
        if(ninlet.gt.0)then
           xinlet=x(inlf(1),1)
        else
           xinlet=1e18
        endif

        if(noutlet.gt.0)then
           xoutlet=x(outf(1),1)
        else
           xoutlet=-1e18
        endif
             
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(xinlet, xinlet_all, 1,
     &       MPI_REAL8, MPI_MIN, MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(xoutlet, xoutlet_all, 1,
     &       MPI_REAL8, MPI_MAX, MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        do nn=1,nshg
           xdel=(x(nn,1)-xinlet_all)/(xoutlet_all-xinlet_all)
           xdel=1-tanh(7*xdel)
           y(nn,1)= xvel_ini + (inletVelX-xvel_ini)*xdel
           y(nn,4)= outPres1 + (pres_ini-outPres1)*xdel
        enddo
           y(:,2)=yvel_ini
           y(:,3)=zvel_ini
           y(:,5)=temp_ini
           if(nsclr.eq.1)then
              y(:,6)=evis_ini
           endif
        elseif (isetInitial.eq.2) then   
            !----------------------------------
            ! set initial condition for blower
            !----------------------------------
 
!            blthickness=1.0e-4
!            yfloor=-0.0333375 - 1e-7
!            xthroat=-0.3048 + 1e-5  
!            xSlit = 0.025            !x-coordinate of the blower slit
!    
!            p0 = 121800       !total pressure in blower
!            T0 = 305          !total temperature
!                
!            do nn=1,nshg
!                xcoor = x(nn,1)
!
!                if((xcoor .gt. xthroat) .and. (x(nn,2).lt.yfloor)) then
!                    y(nn,3:4)=0  !velocity
!                    y(nn,6)=1e-3 !eddy viscosity
!               
!                    if(xcoor .le. xSlit) then  !test whether the point is infront of or behind the slit. If infront, use M = 0.8
!                        !4th order polynomial best fit for an exit Mach number of 0.8
!                        MEuler=0.0806961 + xcoor*(2.96302 + xcoor*(72.5101 + xcoor*(13914.6 + xcoor*1.09778e6))) 
!                    else
!                        MEuler = 0.8
!                    endif
!    
!                    pEuler = p0/((1 + 0.2*M*M)**3.5)
!                    TEuler = T0/(1 + 0.2*M*M)           
!                    uEuler = MEuler*sqrt(1.4*287.06*TEuler)
!    
!                    bldamper=min(one, d2wall(nn)/blthickness)
!                    y(nn,2)=uEuler*bldamper
!                    y(nn,1)=pEuler
!                    y(nn,5)=TEuler
!                endif
!             enddo
        endif

        return
        end

c=========================================================================
c Set the initial condition based on coordiantes
c=========================================================================

        subroutine setInitial_Duct(x)

        use readarrays
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        real*8 x(nshg,nsd)
        real*8 xcoor,ycoor,zcoor
        integer nn 
        real*8 xVel, pres, Temp
!p u v w T s(calar_1) in array qold 
!x y z in array x 



        a1 =       329.7  
        b1 =    -0.07732  
        c1 =        9.82  
        a2 =       6.431  
        b2 =       1.621  
        c2 =      0.8998  
        a3 =      -4.637  
        b3 =       1.728  
        c3 =      0.2168  


        if (iDuctgeometryType .eq. 6) then
           xContraStart = -74*0.0254;
        elseif (iDuctgeometryType .eq. 8) then
           xContraStart = -84*0.0254;
        endif
        xContraEnd = xContraStart + 72*0.0254
         
        do nn=1,nshg

           xcoor=x(nn,1)          
           ycoor=x(nn,2)
           zcoor=x(nn,3)
c.... [xVel,pres,Temp] are functions of [xcoor,ycoor,zcoor]            

c          if(ycoor<(-0.587375+1.0e-4))then
c             ry=(ycoor+0.587375)/1.0e-4
c          elseif(ycoor>(0.587375+1.0e-4))then
c             ry=(0.587375-ycoor)/1.0e-4
c          else
c             ry=1.0 
c          endif

c          if(zcoor<(-0.587375+1.0e-4))then
c             rz=(zcoor+0.587375)/1.0e-4
c          elseif(zcoor>(0.587375+1.0e-4))then
c             rz=(0.587375-zcoor)/1.0e-4
c          else
c             rz=1.0
c          endif

c          ry=max(0.0,ry)
c          rz=max(0.0,rz)
         
c          if(xcoor.le.0.0.and.xcoor.gt.-12*0.0254)then
c             xVel=120
c          elseif(xcoor.lt.-12*0.0254)then
c             xVel=(1.013158+120)/2+(120-1.013158)/2*
c     &       tanh((xcoor-(-84*0.0254-12*0.0254)/2)*4)
c          else
c             xVel=max(0.0,120*(1-(xcoor/(1.6*4.5*0.0254+0.85))**2))
c          endif           

c          Temp=317+(-tanh(4*(xcoor+0.35+0.254))*
c     &       (13.578+5)*0.5+(13.578-5)*0.5)*ry*rz
c          pres=-tanh(4*(xcoor+0.35+0.254))*
c     &    (116600-97000)*0.5+(116600+97000)*0.5  
c... above, Onkar's method

c         if(xcoor.lt.0.0 .and. 
c     &    xcoor.gt.12*0.0254 .and. 
c     &    ycoor .lt. -0.04446)then
c            xVel=0
c            pres=97800
c            Temp=317
c         else   
c            if(xcoor .ge. -0.0254)then   
c             xVel=158            
c             pres=97800           
c             Temp=318
c            elseif(xcoor .le. -1.0163)then
c             Temp=330
c             pres=111896
c             xVel=1.0132
c            else
c             pres=97800+0.5*(111896-97800)*
c     &       (tanh(-3.4*((xcoor+0.52085)/0.49545))+1)
c             Temp=318+0.5*(330-318)*
c     &       (tanh(-3.4*((xcoor+0.52085)/0.49545))+1)
c             xVel=158+0.5*(1.0132-158)*
c     &       (tanh(-3.4*((xcoor+0.52085)/0.49545))+1)
c            endif
c         endif  
c....................
            

c             pres= 97000
c             xvel= 0.1

c             if (xcoor .gt. xContraEnd)then
c                 Temp = 320
c             elseif (xcoor . le. xContraEnd)then
c                 xSim = xcoor - xContraStart
c                 Temp = a1*exp(-((xSim-b1)/c1)**2) + 
c     &                  a2*exp(-((xSim-b2)/c2)**2) +
c     &                  a3*exp(-((xSim-b3)/c3)**2)
c             endif 
c             varSA = 0.0

c             qold(nn,1)=pres
c             qold(nn,2)=xvel
c             qold(nn,3)=0
c             qold(nn,4)=0
c             qold(nn,5)=Temp
c             qold(nn,6)=varSA

c... above Yi Chen's Method
                 
          if(xcoor .lt. 0.0 .and.
     &       xcoor .gt. (-12*0.0254) .and.
     &       ycoor .lt. -0.044449)then ! this is the jet device
             qold(nn,1)=97000
             qold(nn,2)=0
             qold(nn,3)=0
             qold(nn,4)=0
             qold(nn,5)=320
             qold(nn,6)=1.825e-5
          endif

c... above, used for restart from M2M

        enddo ! end of loop over all nshg points 
c........................... 

        return
        end


