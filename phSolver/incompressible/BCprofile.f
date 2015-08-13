c-----------------------------------------------------------------------
        subroutine BCprofileScale(vbc_prof,BC,yold)
        use pvsQbi
        include "common.h"
        real*8 vbc_prof(nshg,3)
        real*8 BC(nshg,ndofBC),yold(nshg,ndof)
        real*8 offphase
        integer factor
              
        r_amp =rampmdot(1,1)
        r_freq=rampmdot(2,1)
!  Usual sinusoidal in time syn jet is given in the next line.  It assumes that you are NOT changing the time step from previous runs since it computes the 
!  current time in the sin function to be (lstep+1)* Dt where lstep is a running total of all of the runs up to now.  This will be "ok" under the following
!  assumption:  lstep*dt_current = m * T where m is an integer and T is the period of the jet.  That is because this will evaluate to 2*m*Pi and the sin
!  that is zero

        if(abs(rampmdot(1,2)) .lt. 1.0e-5) then
           r_time_factor = r_amp*sin(two*pi*r_freq*(lstep+1)*Delt(1)) ! BC set for next time step - all phases are computed
        else
           r_time_factor = rampmdot(1,2)*r_amp   ! read parameter scales Vmax
        endif

        icount = 0
        do kk=1,nshg
          if(ndsurf(kk).ge.1 .and. ndsurf(kk).le.12) then ! this means all the SJ BC (1->12) for the Beta (Boeing) model

            factor = idnint(rampmdot(2,2))

            if(factor == 0) then
              offphase = 0.d0

            elseif(factor == 1) then
              offphase = 1.d0

            elseif(factor == 2) then
              !Start to count from tip. If factor == 2 -> 1, 3, 5, etc
              if(mod(ndsurf(kk)-1,factor) == 0) then   
                offphase = 1.d0
              else
                offphase = 0.d0
              endif

            elseif(factor == 3) then
              !Start to count from tip. If factor == 3 -> 2, 5, 8, etc
              if(mod(ndsurf(kk)+1,factor) == 0) then   
                offphase = 1.d0
              else
                offphase = 0.d0
              endif

            elseif(factor < 0) then
              !Only one jet blowing. If factor = -5, then jet 5 only
              !blows
              if (ndsurf(kk) == -factor) then 
                offphase = 1.d0
              else 
                offphase = 0.d0
              endif
            endif

            BC(kk,3)=r_time_factor*vbc_prof(kk,1)*offphase
            BC(kk,4)=r_time_factor*vbc_prof(kk,2)*offphase
            BC(kk,5)=r_time_factor*vbc_prof(kk,3)*offphase

            icount = icount + 1

          endif
        enddo

!        if(istep.eq.0 .and. icount.ne.0)
!     &     write(*,*) 'BCprofile count',myrank,icount
    
        return
        end 

c--------------------------------------------------------------
        subroutine BCprofileInit(vbc_prof,x)
      
        use pvsQbi
        include "common.h"
        real*8 vbc_prof(nshg,3), x(numnp,nsd)
        real*8 rcenter(3),rnorml(3)
 
!        open(unit=789, file='bcprofile.dat',status='unknown')
        do kk=1,nshg
          if(ndsurf(kk).eq.1) then     ! SJ1
            rnorml(1)= -0.03805041312991605
            rnorml(2)= 0.99925652732035097
            rnorml(3)= 0.00620956265097128
            rcenter(1)= 2.31206845913690984
            rcenter(2)= 0.39594943486752049
            rcenter(3)= 0.51243427055630897
        
          elseif(ndsurf(kk).eq.2) then     ! SJ2          
            rnorml(1)= -0.03805041312991605
            rnorml(2)= 0.99925652732035097
            rnorml(3)= 0.00620956265097128
            rcenter(1)= 2.28612828683167457
            rcenter(2)= 0.39519558316178066
            rcenter(3)= 0.47479183928994600

          elseif(ndsurf(kk).eq.3) then     ! SJ3
            rnorml(1)= -0.03805041312991605
            rnorml(2)= 0.99925652732035097
            rnorml(3)= 0.00620956265097128
            rcenter(1)= 2.26019111197956502
            rcenter(2)= 0.39444184406245603
            rcenter(3)= 0.43714965468096900

          elseif(ndsurf(kk).eq.4) then     ! SJ4
            rnorml(1)= -0.03806860790482020
            rnorml(2)= 0.99925599805332099
            rnorml(3)= 0.00618315830704482
            rcenter(1)= 2.22406898326471003
            rcenter(2)= 0.39338941696418994
            rcenter(3)= 0.38472970195368150

          elseif(ndsurf(kk).eq.5) then     ! SJ5
            rnorml(1)= -0.03807026719147060
            rnorml(2)= 0.99925594973518905
            rnorml(3)= 0.00618075034236115
            rcenter(1)= 2.19813166796875015
            rcenter(2)= 0.39263411631338135
            rcenter(3)= 0.34709106843832799

          elseif(ndsurf(kk).eq.6) then     ! SJ6
            rnorml(1)= -0.03807076874492700
            rnorml(2)= 0.99925593512835698
            rnorml(3)= 0.00618002248559851
            rcenter(1)= 2.17219205684703009
            rcenter(2)= 0.39187870083238074
            rcenter(3)= 0.30944277157095951

          elseif(ndsurf(kk).eq.7) then     ! SJ7
            rnorml(1)= -0.03786596844535725
            rnorml(2)= 0.99926183452064998
            rnorml(3)= 0.00647722966396847
            rcenter(1)= 2.13621068971991468
            rcenter(2)= 0.39081967107127247
            rcenter(3)= 0.25712280227877149

          elseif(ndsurf(kk).eq.8) then     ! SJ8
            rnorml(1)= -0.03805101631143114
            rnorml(2)= 0.99925650979092495
            rnorml(3)= 0.00620868731108326
            rcenter(1)= 2.11020507860295004
            rcenter(2)= 0.39007468380778049
            rcenter(3)= 0.21948491547550850

          elseif(ndsurf(kk).eq.9) then     ! SJ9
            rnorml(1)= -0.03805101631143114
            rnorml(2)= 0.99925650979092495
            rnorml(3)= 0.00620868731108326
            rcenter(1)= 2.08426231456637501
            rcenter(2)= 0.38931079611337499
            rcenter(3)= 0.18184066095664247

          elseif(ndsurf(kk).eq.10) then    ! SJ10
            rnorml(1)= -0.03427503641589490
            rnorml(2)= 0.99934408637497807
            rnorml(3)= 0.01168840904698280
            rcenter(1)= 2.04819760342475998
            rcenter(2)= 0.38841663756682199
            rcenter(3)= 0.12950419445094899

          elseif(ndsurf(kk).eq.11) then    ! SJ11
            rnorml(1)= -0.03208677268269435
            rnorml(2)= 0.99937455421864496
            rnorml(3)= 0.01486403037855546
            rcenter(1)= 2.02226290188708502
            rcenter(2)= 0.38805554104084300
            rcenter(3)= 0.09186760148390730

          elseif(ndsurf(kk).eq.12) then     ! SJ12
            rnorml(1)= -0.02932216953382100
            rnorml(2)= 0.99939176776025096
            rnorml(3)= 0.01887604055064800
            rcenter(1)= 1.99632633336121512
            rcenter(2)= 0.38789394943524053
            rcenter(3)= 0.05422828360356270

          endif

           rdistfromcsq=(x(kk,1)-rcenter(1))*(x(kk,1)-rcenter(1))+
     &                  (x(kk,2)-rcenter(2))*(x(kk,2)-rcenter(2))+
     &                  (x(kk,3)-rcenter(3))*(x(kk,3)-rcenter(3))

           radius = 0.0183896  ! 0.724 inches in meters

c.............Factors below are negative for desired blowing direction
           if(ndsurf(kk).ge.1 .and. ndsurf(kk).le.12) then  ! All jet blowing profiles
!           if(ndsurf(kk).eq.9) then  ! Only jet 9  blowing for now

              vbc_prof(kk,1)=-rnorml(1)*(1-rdistfromcsq/(radius*radius))
              vbc_prof(kk,2)=-rnorml(2)*(1-rdistfromcsq/(radius*radius))
              vbc_prof(kk,3)=-rnorml(3)*(1-rdistfromcsq/(radius*radius))

!              write(789,987) kk,vbc_prof(kk,1),vbc_prof(kk,2),
!     &                          vbc_prof(kk,3)

           else
              vbc_prof(kk,:)=zero
           endif

        enddo
!        close(789)
!987     format(i6,3(2x,e14.7))

        return
        end
