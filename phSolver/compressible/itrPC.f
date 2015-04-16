c-----------------------------------------------------------------------
c
c    Initialize the predictor multicorrector (set up parameters)
c
c-----------------------------------------------------------------------
      subroutine itrSetup ( y,  acold ) 
      
      include "common.h"
      
      real*8     y(nshg,ndof),  acold(nshg,ndof)
      
c
c  Define the Hulbert parameters
c  second order if between (and including) 0 and 1 otherwise backward Euler
c
      if( rhoinf(itseq).lt.0.or.rhoinf(itseq).gt.1) then ! backward Euler
         almi   = one
         alfi   = one
         gami   = one
         ipred  = 1
      else           !second order family
         almi   = (three-rhoinf(itseq))/(one+rhoinf(itseq))/two
         alfi   = one/(one+rhoinf(itseq))
         gami   = pt5+almi-alfi
      endif
c     
c.... set the jacobian type
c     
      Jactyp=0
c
      if(ipred.eq.4 .and. itseq .eq. 1 ) 
     &           y=y-(one-alfi)*Delt(1)*CFLfl(1)*acold

c
c protect from ipred=4 and rhoi=0
c
      if(ipred.eq.4 .and. rhoinf(itseq) .eq. 0.0 ) ipred=3
      
c
c.... set the global time increment and the CFL data
c
      Dtgl   = one / Delt(itseq)  ! caution: inverse of time step
      CFLfld = CFLfl(itseq)
      CFLsld = CFLsl(itseq)
      
      Dtgl   = Dtgl / 1.0      ! check CFLfld
      return
      end


c-----------------------------------------------------------------------
c
c    Predict solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrPredict ( yold,          acold,
     &                        y,             ac   )
      
      include "common.h"
      
      real*8        yold(nshg,ndof),            acold(nshg,ndof),
     &              y(nshg,ndof),               ac(nshg,ndof)

c
c.... set the solution at t_{n+alfi) and the acceleration at t_{n+almi)
c
       if(ipred.eq.1) then
c
c   yn+1_pred=yn
c
        y  = yold
        ac = acold*(one-almi/gami)
       endif
c
       if(ipred.eq.2) then
c
c   an+1_pred=0
c
        y  = yold+alfi/Dtgl*acold*(one-gami)
        call itrBC (y, ac,  iBC,  BC,  iper,ilwork)
c Elaine-SPEBC
        if((irscale.ge.0).and.(myrank.eq.master)) then
           call genscale(y, x, iBC) 
        endif
        ac = acold*(one-almi)
       endif
c
       if(ipred.eq.3 ) then
c
c   an+1_pred=an
c
        y  = yold+alfi/Dtgl*acold
        call itrBC (y, ac,  iBC,  BC, iper,ilwork)
c Elaine-SPEBC
        if((irscale.ge.0).and.(myrank.eq.master)) then
           call genscale(y, x, iBC)
        endif
        ac = acold
       endif
c
       if(ipred.eq.4 ) then ! protect from DC=4 rho=0
c
c  same dV
c
        fct1=alfi/(one-alfi)
        fct2=one-almi/gami
        fct3=almi/gami/alfi*Dtgl
        y  = yold+fct1*(yold-y)
        call itrBC (y, ac,  iBC,  BC,  iper,ilwork)
c Elaine-SPEBC
        if((irscale.ge.0).and.(myrank.eq.master)) then
           call genscale(y, x, iBC)
        endif
        ac = acold*fct2+(y-yold)*fct3
       endif

c
        if (LCtime .eq. 0) time = time + Delt(itseq)
c
      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrCorrect ( y, ac,yold, acold,Dy)
      
      include "common.h"
      
      real*8   y(nshg,ndof), ac(nshg,ndof), Dy(nshg,nflow)
      real*8   yold(nshg,ndof), acold(nshg,ndof)
c      
c         
c     What we actually solved for here was delta y_{n+alpha_f).
c     This is the variable that we work with in the iteration loop
c     but we have another variable that must be updated
c
      y(:,1:3) = y(:,1:3) - Dy(:,2:4)
      y(:,4)   = y(:,4)   - Dy(:,1)
      y(:,5)   = y(:,5)   - Dy(:,5)
c     
      fct1= (one-almi/gami)
      fct2= almi*Dtgl/gami/alfi ! don't be confused by inverse dt=Dtgl
      ac= acold*fct1+(y-yold)*fct2
      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrCorrectSclr (y, ac, yold, acold, Dyt)
      
      include "common.h"
      
      real*8   y(nshg,ndof), ac(nshg,ndof), Dyt(nshg)
      real*8   yold(nshg,ndof), acold(nshg,ndof)
c      
      is=5+isclr
      y(:,is)  = y(:,is) - Dyt(:)
c    
      fct1= (one-almi/gami)
      fct2= almi*Dtgl/gami/alfi ! don't be confused by inverse dt=Dtgl
      ac(:,is)= acold(:,is)*fct1+(y(:,is)-yold(:,is))*fct2 
c
      return
      end
c-----------------------------------------------------------------------
c
c    Update solution at end of time step
c
c-----------------------------------------------------------------------
      subroutine itrUpdate( yold,          acold,
     &                      y,             ac )

      include "common.h"
      
      real*8        yold(nshg,ndof),            acold(nshg,ndof),
     &              y(nshg,ndof),               ac(nshg,ndof)

      
      if(iLset.eq.2 .and. isclr.eq.2) then  
c ... assigning the redistanced scalar to the first scalar 
c     
         y(:,6)    = y(:,7)
c     
c ... to calculate the right acceleration based on the updated value of the 
c ... scalar, Ref: Eq.25 in generalised alpha method paper by Prof.Jansen
c     
         acold(:,6)= acold(:,6)*(1-1/gami) 
     &             + (y(:,6)-yold(:,6))*Dtgl/gami
         fct2=one/almi
         fct3=one/alfi
         acold(:,1:nflow) = acold(:,1:nflow) 
     &                    + (ac(:,1:nflow)-acold(:,1:nflow))*fct2
         yold(:,1:nflow)  = yold(:,1:nflow)  
     &                    + (y(:,1:nflow)-yold(:,1:nflow))*fct3  
         yold(:,6) = y(:,7)
         acold(:,7)= zero  
      else
         fct2=one/almi
         fct3=one/alfi
         acold = acold + (ac-acold)*fct2
         yold  = yold  + (y-yold)*fct3   
      endif  
c
      return
      end







