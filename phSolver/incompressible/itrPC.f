c-----------------------------------------------------------------------
c
c    Initialize the predictor multicorrector (set up parameters)
c
c    Modified by Alberto Figueroa.  Winter 2004
c-----------------------------------------------------------------------
      subroutine itrSetup ( y,  acold ) 
      
      include "common.h"
      
      real*8     y(nshg,ndof),     acold(nshg,ndof)
      
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
         if(ideformwall.eq.1) then
            betai=1.0
         else
            betai=0.0
         endif
      endif
c     
c.... set the jacobian type
c     
      Jactyp=0
      if (impl(itseq) .eq. 3) then
         Jactyp = 1
         impl(itseq) = 2
      endif
c     
c.... same_Dy predictor special case
c     
      if (ipred.eq.4 .and. itseq .eq. 1 ) then
         y=y-(one-alfi)*Delt(1)*acold
         if ( rhoinf(itseq) .eq. 0.0 ) then
            ipred = 3
         endif
      endif
c
c.... set the global time increment and the CFL data
c
      Dtgl   = one / Delt(itseq)  ! caution: inverse of time step
      CFLfld = CFLfl(itseq)
      CFLsld = CFLsl(itseq)
      
      return
      end


c-----------------------------------------------------------------------
c
c    Predict solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrPredict (yold,  y,  acold,   ac,   uold,   u, iBC)
      
      include "common.h"
      
      real*8        y(nshg,ndof),               ac(nshg,ndof),
     &              u(nshg,nsd),                yold(nshg,ndof),
     &              acold(nshg,ndof),           uold(nshg,nsd)
      dimension     iBC(nshg)

      
      if ( ipred.eq.1) then     ! yn+1_pred=yn
         fct = (gami-one)/gami
         y  = yold
         ac = acold * fct
         if(ideformwall.eq.1) then
            do isd = 1, nsd
               where(btest(iBC,13))
                  u(:,isd) = uold(:,isd) + Delt(itseq)*yold(:,isd) + 
     &                 pt5*((gami-two*betai)/gami)*
     &                 Delt(itseq)*Delt(itseq)*acold(:,isd)
               endwhere
            enddo
         endif
      endif
c     
      if ( ipred.eq.2) then     ! an+1_pred=0
         y  = yold + (one - gami)/Dtgl * acold
         ac = 0.0
      endif
c     
      if(ipred.eq.3 ) then      ! an+1_pred=an
         y  = yold+alfi*Delt(itseq)*acold
         ac = acold
      endif
c     
      if ( ipred.eq.4 ) then    ! protect from DC=4 rho=0, same dV
         fct1 = alfi/(one-alfi)
         fct2 = one-almi/gami
         fct3 = almi/gami/alfi*Dtgl
         y    = yold+fct1*(yold-y)
         ac   = acold*fct2+(y-yold)*fct3
      endif
c     

      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrCorrect ( y,     ac,   u,   solinc ,iBC)
      
      include "common.h"
      
      real*8      y(nshg,ndof),               ac(nshg,ndof),  
     &            u(nshg,nsd),                solinc(nshg,4)
      dimension   iBC(nshg)

      fct1 = gami*Delt(itseq)
      fct2 = gami*alfi*Delt(itseq)
      
      y(:,1:3)  = y(:,1:3)  + fct1 * solinc(:,1:3)
      y(:,4  )  = y(:,4  )  + fct2 * solinc(:,4  )
      ac(:,1:3) = ac(:,1:3) + solinc(:,1:3)
      if(ideformwall.eq.1) then
         do isd = 1,nsd
            where (btest(iBC,13))
               u(:,isd)  = u(:,isd)  + 
     &              Delt(itseq)*Delt(itseq)*betai*solinc(:,isd)
            endwhere
         enddo
      endif
c     
      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1
c
c-----------------------------------------------------------------------
      subroutine itrCorrectSclr ( y,     ac,   solinc )
      
      include "common.h"
      
      real*8      y(nshg,ndof),               ac(nshg,ndof),  
     &            solinc(nshg)
      
      fct1 = gami*Delt(itseq)
      is=5+isclr
      y(:,is)  = y(:,is)  + fct1 * solinc(:)
      ac(:,is) = ac(:,is) + solinc(:)
c     
      return
      end

c-----------------------------------------------------------------------
c
c    Correct solution at time n+1, protecting against negative values
c
c-----------------------------------------------------------------------
      subroutine itrCorrectSclrPos ( y,     ac,   solinc )
      
      include "common.h"
      
      real*8      y(nshg,ndof),               ac(nshg,ndof),  
     &            solinc(nshg),               posinc(nshg)
      
      fct1 = gami*Delt(itseq)
      if(fct1.eq.0) then
         call itrCorrectSclr( y, ac, solinc)
      else
         turbUpdFct = 1.0
         updFct = turbUpdFct
         updFct2 = 1.0
c         if(topHdTscFlag.and.PRM_TCS_PARTIAL_UPDATE)then
c            updFct2 = max(MTH_SQRT_EPS,min(topHdUpdFct,1.0))
c         endif
         fctCr = fct1*updFct*updFct2
         c0 = 0.01
         c1 = -(1.0-c0)/fctCr
         is=5+isclr
c         if(any(updFct*updFct2*solinc(:).lt.c1*y(:,is))) then
c            write(*,*) 'INSTEAD OF GETTING NEGATIVE FIELD'
c            write(*,*) 'BROUGHT FIELD DOWN TO 1 PERCENT'
c            write(*,*) 'FOR SCALAR NUMBER ',isclr
c            write(*,*) '(SEE itrCorrectSclr in itrPC.f)'
c         endif
         posinc = max(updFct*updFct2*solinc,c1*y(:,is))
         y(:,is)  = y(:,is)  + fct1 * posinc(:)
         ac(:,is) = ac(:,is) + posinc(:)
      endif
c     
      return
      end


c-----------------------------------------------------------------------
c
c    Compute solution and acceleration at n+alpha
c
c-----------------------------------------------------------------------
      subroutine itrYAlpha ( uold,        yold,        acold,        
     &                       u,           y,           ac,
     &                       uAlpha,      yAlpha,      acAlpha )

c      use readarrays       !reads in uold and acold     
      include "common.h"
      
      real*8        yold(nshg,ndof),            acold(nshg,ndof),
     &              y(nshg,ndof),               ac(nshg,ndof),
     &              yAlpha(nshg,ndof),          acAlpha(nshg,ndof),
     &              u(nshg,nsd),                uold(nshg,nsd),
     &              uAlpha(nshg,nsd)

      acAlpha(:,4) = zero  !pressure acceleration is never used but....

      yAlpha (:,1:3) = yold(:,1:3) 
     &                  + alfi * (y(:,1:3) - yold(:,1:3))

      acAlpha(:,1:3) = acold(:,1:3)
     &                  + almi * (ac(:,1:3) - acold(:,1:3))

      yAlpha (:,4  ) = y(:,4)

      if(ideformwall.eq.1) uAlpha (:,1:3) = uold(:,1:3) 
     &                  + alfi * (u(:,1:3) - uold(:,1:3))
      
      if(ndof.ge.5) then
c
c  Now take care of temperature, turbulence, what have you
c
      

         yAlpha (:,5:ndof  ) = yold(:,5:ndof) 
     &                       + alfi * (y(:,5:ndof) - yold(:,5:ndof))
         acAlpha(:,5:ndof  ) = acold(:,5:ndof) 
     &                       + almi * (ac(:,5:ndof) - acold(:,5:ndof))

      endif
      return
      end

c-----------------------------------------------------------------------
c
c    Update solution at end of time step
c
c-----------------------------------------------------------------------
      subroutine itrUpdate( yold,          acold,        uold,
     &                      y,             ac,           u )

c      use readarrays            !reads in uold and acold

      include "common.h"
      
      real*8        yold(nshg,ndof),            acold(nshg,ndof),
     &              y(nshg,ndof),               ac(nshg,ndof),
     &              u(nshg,nsd),                uold(nshg,nsd)

            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!	Allow for redistance of the level set function, AD 5/8/00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      if we are re-distancing the levelset function the we need 
!      to update the levelset variable with the corrected value
!      stored in the "next" scalar position.
            
        if(iLSet .eq.2 .and. isclr .eq. 2) then 
     
           y(:,6) = y(:,7) ! will need to also fix the yold ac etc.?????
     	   ac(:,7)=zero	
        endif
        

      
      yold  = y
      acold = ac
      if(ideformwall.eq.1)  uold  = u
      
      return
      end

