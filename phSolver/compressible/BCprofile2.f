        subroutine INIBCprofile(BC,y,x)
        
         include "common.h"
         real*8 BC(nshg,ndofBC)
         real*8 y(nshg,ndof)
         real*8 x(nshg,nsd)
 
         if(iI2Binlet .gt. 0)then
              call TakeBCfromIC_ScaleInlVelX(BC,y,x)
         endif           

         if(isetOutPres .gt. 0)then
              call SetUniOutPres(BC,y,x)
         endif

         if(isetInitial .gt.0)then
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
        subroutine SetUniOutPres(BC,y,x)
          
        use BCsfIDmap
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        integer outf_tmp(nshg)
        real*8 BC(nshg,ndofBC)
        real*8 y(nshg,ndof)
        real*8 x(nshg,nsd)

        call sfID2np(isetOutPres,noutlet,outf_tmp)  
        if(noutlet .gt. 0)then
           allocate(outf(noutlet))
           outf(1:noutlet)=outf_tmp(1:noutlet)
           BC(outf(:),1)=outPres1
        endif
         
        return 
        end
c..............................................................................
c... this subroutine gives the formular of initial conditions, user may want to
c... change this file and recompile the code, just like other commercial code which
c... supports User Defined Function, you need to recompile and relink the code  
        subroutine setInitial(x,y)

        use BCsfIDmap
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        real*8 x(nshg,nsd)
        real*8 y(nshg,ndof)         

c... user gives the formular of initial conditions from here       
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
c........................... 

        return
        end


