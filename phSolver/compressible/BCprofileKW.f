c-----------------------------------------------------------------------
        subroutine BCprofileScaleKW(PresBase,VelBase,BC,yold)
       
        use pvsQbi
        include "common.h"
        real*8 PresBase, VelBase
        real*8 BC(nshg,ndofBC), yold(nshg,ndof)
        real*8 PresSin, PresAmp, PresFreq
        real*8 Alpha, AlphaAmp, AlphaFreq
        integer tsBase, BCdtKW

c        PresFreq=1000 ! frequency in Hz of pressure disturbance
c        PresAmp=100 ! Amplitude of pressure disturbance
c        tsBase=100
        PresSin=PresAmp*sin(two*pi*PresFreq*((lstep-tsBase)+1)*Delt(1)) ! pressure disturbance
c        AlphaFreq=1000 ! frequency in Hz of changhe in angle of attack
c        AlphaAmp=1 ! Max angle of attack in degrees
        Alpha=AlphaAmp*sin(two*pi*AlphaFreq*((lstep-tsBase)+1)*Delt(1)) ! pressure disturbance

        if(BCdtKW.eq.1) then
          do kk=1,nshg
            if(ndsurf(kk).eq.801) then ! this means inflow plane of inlet

c              BC(kk,1)=PresBase+PresSin
              BC(kk,1)=11597+PresSin

            endif
          enddo
        else
          do kk=1,nshg
            if(ndsurf(kk).eq.801) then ! this means inflow plane of inlet

              BC(kk,3)=VelBase*cos((pi/180)*Alpha)
              BC(kk,4)=VelBase*sin((pi/180)*Alpha)

            endif
          enddo
        endif

        return
        end 

c--------------------------------------------------------------
        subroutine BCprofileInitKW(PresBase,VelBase,BC)
      
        use pvsQbi
        include "common.h"
        real*8 PresBase, VelBase
        real*8 PresAmp, PresFreq
        real*8 AlphaAmp, AlphaFreq
        real*8 BC(nshg,ndofBC)
        integer iflagKW
        integer tsBase, BCdtKW
 
c        open(unit=789, file='bcprofile.dat',status='unknown')
c        tsBase=lstep
        iflagKW=1
        do kk=1,nshg
          if(iflagKW.eq.1) then
            if(ndsurf(kk).eq.801) then

              PresBase=BC(kk,1)
              VelBase=BC(kk,3)
              iflagKW=0

            endif
          endif
        enddo
c         close(789)
c 987     format(i6,3(2x,e14.7))

        return
        end
