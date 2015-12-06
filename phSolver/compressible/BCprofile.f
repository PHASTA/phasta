c-----------------------------------------------------------------------
c
c  This module conveys ramped BC data. 
c
c-----------------------------------------------------------------------
        module rampBC
        integer kmaxinf, kmaxbot,kmaxtop,nfctop,nfcbot,ninflow
        real*8 ctrl_a(3)
        real*8 ctrl_factor(3)
        end module
        
        subroutine initBCprofileScale(vbc_prof,BC,yold,x)
        use rampBC
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

! assuming three profiles to control (inlet, bottom FC and top FC)
        real*8 vbc_prof(nshg,3), loc_ctrl_pr(3), ctrl_pr(3)
        real*8 BC(nshg,ndofBC),yold(nshg,ndof),x(nshg,3)        
        nstp=nstep(1)
! set constant factors for each profile ((V_peak/V_bar)*R*T/A)
c
c ALL COMPUTABLE WITH READ DATA BUT FOR NOW USER MUST COMPUTE FOR 
c EACH GEOMETRY + PROFILE
C
        ctrl_factor(1) = 1.00*Rgas*330.578/((2*0.587375)**2) ! inlet
c        ctrl_factor(1) = 1.00*Rgas*293/((0.0889*0.1143)) ! inlet
c no chamber       ctrl_factor(2) = 2.25*Rgas*293/((0.02441024589128274-0.008949756886432915)*(2*0.0508)) ! bottom FC
        ctrl_factor(2) = 2.25*Rgas*317/((0.0168910950375201
     &                        -0.01054110785157094)*(2*0.041275)) ! bottom FC
        ctrl_factor(3) = 2.25*Rgas*317/(((-0.126612)
     &                  -(-0.156330))*(2*0.0508)) ! top FC
c        ctrl_factor(2) = 2.25*Rgas*293/(lbot *wbot)           
c        ctrl_factor(1) = CfactInl ! inlet
c        ctrl_factor(2) = CfactBot ! bottom FC
c        ctrl_factor(3) = CfactTop ! top FC
c2D        ctrl_factor(1) = 1.0*Rgas*273/0.046990500 ! inlet
c2D        ctrl_factor(2) = 1.5*Rgas*273/0.000618413 ! bottom FC
c2D        ctrl_factor(3) = 1.5*Rgas*273/0.001188690 ! top FC

!set the first first step
        ninflow=0
        nfctop=0
        nfcbot=0
!assuming max at center line
        kmaxinf=0
        kmaxbot=0
        kmaxtop=0
        rmaxvinflow=-one
        rmaxvbot=-one
        rmaxvtop=-one
c
c find the peak value (and node number) of each inlet
c   SWITCH THIS TO USE SURFID's 
c 
        do kk=1,numnp
           if(vbc_prof(kk,1).ne.zero) then
              if(vbc_prof(kk,2).ne.zero) then ! both of these true means top FC
                 if(abs(vbc_prof(kk,1)).gt.rmaxvtop) then
                    kmaxtop=kk
                    rmaxvtop=abs(vbc_prof(kk,1))
                 endif
		 nfctop=nfctop+1
              else              ! only x true means inflow
                 if(vbc_prof(kk,1).gt.rmaxvinflow) then
                    kmaxinf=kk
                    rmaxvinflow=vbc_prof(kk,1)
                 endif
		 ninflow=ninflow+1
              endif
           else if(vbc_prof(kk,2).ne.zero) then ! y true means low FC
              if(vbc_prof(kk,2).gt.rmaxvbot) then
                 kmaxbot=kk
                 rmaxvbot=vbc_prof(kk,2)
              endif
		 nfcbot=nfcbot+1
           endif
        enddo
c
c  call BCprofileScale. Since called form init, this is going to be used to set
c  the yold value on the first step. Consequently, we will need to set istep to
c  -1 temporarily so that when it looks to the target value at the "end" of the c  step it will actually be looking to time zero. Note we will reset it back to
c  zero after this call
c
        istep=-1 !offset the istep+1 increment so that it gets t0 right
        call BCprofileScale(vbc_prof,BC,yold)
        istep=0  !set back to correct value
        return
        end 

        subroutine BCprofileScale(vbc_prof,BC,yold)
        use rampBC
        use timedataC
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
! assuming three profiles to control (inlet, bottom FC and top FC)
        real*8 vbc_prof(nshg,3), loc_ctrl_pr(3), ctrl_pr(3)
        real*8 BC(nshg,ndofBC),yold(nshg,ndof),mdotNow(3)  
c
c in some way find mdot for this way
c
        call rampedMdot(mdotNow)
c       call flowControl(mdotNow)

c
c compute the pressure from the current solution at the peak of profile
c

        loc_ctrl_pr(1) = yold(kmaxinf,4) ! inlet
        loc_ctrl_pr(2) = yold(kmaxbot,4) ! bottom FC
        loc_ctrl_pr(3) = yold(kmaxtop,4) ! top FC
c
c communicate it so that we are sure all agree
c (in case the profile is split across procs)
c 
       if(numpe > 1) then
           call MPI_ALLREDUCE ( loc_ctrl_pr, ctrl_pr, 3, 
     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr ) 
        else
           ctrl_pr = loc_ctrl_pr
        endif


            ctrl_a(1) = ctrl_factor(1)*mdotNow(1)/ctrl_pr(1) ! inlet
            ctrl_a(2) = ctrl_factor(2)*mdotNow(2)/ctrl_pr(2) ! bottom FC
            ctrl_a(3) = ctrl_factor(3)*mdotNow(3)/ctrl_pr(3) ! top FC
	if(max(ninflow,max(nfctop,nfcbot)).gt. 0) then
            do kk=1,numnp
               if(vbc_prof(kk,1).ne.zero) then
                  if(vbc_prof(kk,2).ne.zero) then ! both of these true means top FC
                     if(ctrl_a(3).ge.zero) then
                        BC(kk,3)=ctrl_a(3)*vbc_prof(kk,1)
                        BC(kk,4)=ctrl_a(3)*vbc_prof(kk,2)
                     endif
                  else ! only x true means inflow
                     if(ctrl_a(1).ge.zero) then
                        BC(kk,3)=ctrl_a(1)*vbc_prof(kk,1)
                     endif
                     if(ctrl_a(1).eq.zero) then
                        BC(kk,7) = 0.0d0 ! scalar_1 set to zero for zero flow
                     endif
                  endif
               else if(vbc_prof(kk,2).ne.zero) then ! y true means low FC
                  if(ctrl_a(2).ge.zero) then
                     BC(kk,4)=ctrl_a(2)*vbc_prof(kk,2)
                  endif
                  if(ctrl_a(2).eq.zero) then
                     BC(kk,7) = 0.0d0 ! scalar_1 set to zero for zero flow
                  endif
               endif
            enddo
c            if(ninflow.gt.0) write(myrank+2000,2000), lstep,kmaxinf,one,
c     &                BC(kmaxinf,3),BC(kmaxinf,4),ctrl_a(1),ctrl_pr(1)
c            if(nfcbot .gt.0) write(myrank+2000,2000), lstep,kmaxbot,two,
c     &                BC(kmaxbot,3),BC(kmaxbot,4),ctrl_a(2),ctrl_pr(2)
c            if(nfctop .gt.0) write(myrank+2000,2000), lstep,kmaxtop,three,
c     &                BC(kmaxtop,3),BC(kmaxtop,4),ctrl_a(3),ctrl_pr(3)
        endif
 2000   format(i6,i6,5(e14.7,2x))
        return
        end

        
        subroutine rampedMdot(mdotNow)
        include "common.h" !gives us rampmndot and nstep

! assuming three profiles to control (inlet, bottom FC and top FC)
        real*8 mdotNow(3),xi,cstep,rnstep

        rnstep=nstep(1)
	cstep=istep+1
	xi=cstep/rnstep

	mdotNow(:)=rampmdot(1,:)+xi*(rampmdot(2,:)-rampmdot(1,:))
	return
        end 
