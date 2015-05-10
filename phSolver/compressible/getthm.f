        subroutine getthm (pres,    T,     Sclr,    rk,   rho,
     &                     ei,      h,     s,       cv,  cp,
     &                     alfap,   betaT, gamb,    c )
c
c-----------------------------------------------------------------------
c
c  This subroutine calculates the thermodynamic properties.
c
c  The different possibilities are:
c   ipress = 0  : calorifically perfect gas
c   ipress = 1  : thermally perfect gas
c   ipress = 2  : mixture of thermally perfect gases in
c                  thermo-chemical equilibrium
c
c  The options available are:
c
c   ithm = 2    : given rho  and T, compute pres and engBC
c   ithm = 3    : given pres and T
c   ithm = 4    : given pres and T,   engBC
c   ithm = 6    : given pres and T, compute rho,   ei
c                 (also given sclr for levelset)
c   ithm = 7    : given pres and T, compute rho,   ei, h, cv,   cp,
c                                           alfap, betaT, gamb, c
c
c Variables:
c
c  pres   (npro)        : pressure
c  T      (npro)        : temperature
c  rk     (npro)        : specific kinetic energy
c  rho    (npro)        : density
c  ei     (npro)        : internal energy
c  h      (npro)        : enthalpy
c  s      (npro)        : entropy
c  cv     (npro)        : specific heat at constant volume
c  cp     (npro)        : specific heat at constant pressure
c  alfap  (npro)        : expansivity
c  betaT  (npro)        : isothermal compressibility
c  gamb   (npro)        : gamma-bar (defined in paper by Chalot et al.)
c  c      (npro)        : speed of sound
c
c
c Zdenek Johan,    Spring 1990.
c Frederic Chalot, Summer 1990.
c Zdenek Johan,    Winter 1991.  (Fortran 90)
c-----------------------------------------------------------------------
c
        include "common.h"
c
        dimension pres(npro),                Sclr(npro),
     &            T(npro),                   rk(npro),
     &            rho(npro),                 ei(npro),
     &            h(npro),                   s(npro),
     &            cv(npro),                  cp(npro),
     &            alfap(npro),               betaT(npro),
     &            gamb(npro),                c(npro),
     &            rsrhol(npro),              rsrhog(npro),
     &            tmpg(npro),                tmpl(npro)             
c
        dimension Texp1(npro),               Texp2(npro)
        real*8 prop_blend(npro),test_it(npro)

c       ttim(27) = ttim(27) - secs(0.0)
c
c.... get the property type flag
c
        ipress = matflg(1,1)
c
c.... ***********************>  IPRESS = 0  <***************************
c
        if (ipress .eq. 0) then
c
c.... --------------------->  ithm = 1 or 2  <--------------------------
c
c
        if (ithm .eq. 2) then
        pres = Rgas * rho * T
c

c.... compute engBC (internal energy in this case)
c
c         engBC = T * (Rgas / gamma1)
c
          flops = flops + npro
c
        endif
c
c.... --------------------->  ithm = 3 or 4  <--------------------------
c
        if ((ithm .eq. 3) .or. (ithm .eq. 4)) then
c
        endif
c
c        if (ithm .eq. 4) then
c
c.... compute engBC (enthalpy in this case)
c
c          engBC = T * (Rgas * gamma / gamma1)
c
          flops = flops + npro
c
c        endif
c
c.... -------------------->  ithm = 5, 6 or 7  <------------------------
c
c
        if (ithm .ge. 6) then
c
c.... compute density and internal energy
c  
          if (iLSet .eq. 0)then 
             rho  = pres / ( Rgas * T )
c
          else     !  two fluid properties used in this model

!        Smooth the tranistion of properties for a "distance" of epsilon_ls
!        around the interface.  Here "distance" is define as the value of the 
!        levelset function.  If the levelset function is properly defined, 
!        this is the true distance normal from the front.  Of course, the 
!        distance is in a driection normal to the front.

             do i= 1, npro
                if (sclr(i) .lt. - epsilon_ls)then
                   prop_blend(i) = zero
                elseif  (abs(sclr(i)) .le. epsilon_ls)then
                   prop_blend(i) = 0.5*(one + Sclr(i)/epsilon_ls +
     &                  (sin(pi*Sclr(i)/epsilon_ls))/pi )
                elseif (sclr(i) .gt. epsilon_ls) then
                   prop_blend(i) = one
                endif
             enddo
                fact = datmat(1,1,2)/datmat(1,1,1)
c               call eqs(pres,T,rsrhol)
c               rsrhol(:) = pres(:) / ( Rgas * T(:) )
c               rsrhog(:)  = fact* pres(:) / ( Rgas * T(:))
                rsrhog(:)  = pres(:) / ( Rgas * T(:))
                rsrhol(:) = datmat(1,1,1)*
     &                      (1+0.000000000517992*(pres-18.02))
                rho(:) = rsrhol(:)*prop_blend(:)+rsrhog(:)
     &               *(1-prop_blend(:))
c ...        for the VOF case..just in case if we want to run VOF
c$$$         prop_blend(:) = min((max(sclr(:),0.0)),1.0)
c$$$         rho(:)=rsrhol(:) * prop_blend(:) + rsrhog(:) * (1-prop_blend(:))
c 
c
            endif
c           Calculate Internal Energy
            if (ilset .eq. 0) then
               ei   = T * (Rgas / gamma1)
            else
               tmpg = T * (Rgas / gamma1) !for gas phase
               tmpl = 3.264*1000*T !cv*T for liquid phase
c              tmpl = T* 8.314* 18.0/ ((3.598/3.264) -1.0)
               ei(:) = tmpl(:)*prop_blend(:)+tmpg*(1-prop_blend(:))
            endif
            flops = flops + 4*npro
         endif ! end if(ithm==6)
c
        if (ithm .ge. 7) then
c
c.... compute enthalpy, cv, cp, alfap, betaT, gamb and c
c
             if (ilset .eq. 0) then
                h     = T * (Rgas * gamma / gamma1)
c
                cv    = Rgas / gamma1
                cp    = Rgas * gamma / gamma1
c     
                alfap = one / T
                betaT = one / pres
c
                gamb  = gamma1
                c     = sqrt( (gamma * Rgas) * T )
c   
             else
c
                do i= 1, npro
                   if (sclr(i) .lt. - epsilon_ls)then
                      prop_blend(i) = zero
                   elseif  (abs(sclr(i)) .le. epsilon_ls)then
                      prop_blend(i) = 0.5*(one + Sclr(i)/epsilon_ls +
     &                     (sin(pi*Sclr(i)/epsilon_ls))/pi )
                   elseif (sclr(i) .gt. epsilon_ls) then
                      prop_blend(i) = one
                   endif
                enddo
                rsrhol(:) = datmat(1,1,1)*
     &                      (1+0.000000000517992*(pres-18.02))
                tmpg = T * (Rgas * gamma / gamma1) ! enthalpy of gas phase
                tmpl = T*3.264*1000 + pres/rsrhol ! enthalpy of liquid phase
c                tmpl = T*8.314*18.0*3.598/(3.598-3.264)
                h = tmpl(:)*prop_blend(:)+tmpg*(1-prop_blend(:))
c
                tmpg = Rgas / gamma1 ! cv for gas phase
                tmpl    = 3.264*1000.0 ! cv=cp for liquid phase
                   cv = tmpl(:)*prop_blend(:)+tmpg*(1-prop_blend(:))
c
                tmpg = Rgas * gamma / gamma1 ! cp for gas phase
                tmpl = 3.264*1000.0 ! cp=cv for liquid phase
                cp = tmpl(:)*prop_blend(:)+tmpg*(1-prop_blend(:))
c
                tmpg = one / T  ! alfap for gas phase
                tmpl = zero+epsM     ! alfap for nearly incompressible liquid
                alfap = tmpl(:)*prop_blend(:)+tmpg*(1-prop_blend(:))
c
                tmpg = one / pres ! betaT for gas phase
                tmpl = 0.000000000517992 ! betaT for nearly incompressible liquid 
                betaT = tmpl(:)*prop_blend(:)+tmpg*(1-prop_blend(:))
c
                tmpg = gamma1   ! gamb for gas phase 
                tmpl = 0.0 ! gamb for liquid phase 
                gamb = tmpl(:)*prop_blend(:)+tmpg*(1-prop_blend(:))
c
                tmpg = sqrt( (gamma * Rgas) * T ) ! c for gas phase
c                tmpl = sqrt( (3.598/3.264 * 8.314*18) * T ) ! for liquid
                tmpl =  sqrt(pres/rsrhol)
                c = tmpl(:)*prop_blend(:)+tmpg*(1-prop_blend(:))
c
                
             endif
             flops = flops + 12*npro
c
          endif
c
c
c.... end of ipress = 0
c
        endif
c
c.... ***********************>  IPRESS = 1  <***************************
c
        if (ipress .eq. 1) then
c
c.... --------------------->  ithm = 1 or 2  <--------------------------
c
c
        if (ithm .eq. 2) then
c
c.... compute engBC (internal energy in this case)
c
c          engBC = Rgas * T / gamma1
c
        endif
c
c.... --------------------->  ithm = 3 or 4  <--------------------------
c
c
c        if (ithm .eq. 4) then
c
c.... compute engBC (enthalpy in this case)
c
c          engBC = Rgas * T * gamma / gamma1
c
c        endif
c
c.... -------------------->  ithm = 5, 6 or 7  <------------------------
c
c
        if (ithm .ge. 6) then
c
c.... compute density and internal energy
c
          Texp1 = exp ( - Tvib(1)/T )
          Texp2 = exp ( - Tvib(2)/T )
c  
          if (iLSet .eq. 0)then 
             rho  = pres / ( Rgas * T )
c
          else     !  two fluid properties used in this model

!        Smooth the tranistion of properties for a "distance" of epsilon_ls
!        around the interface.  Here "distance" is define as the value of the 
!        levelset function.  If the levelset function is properly defined, 
!        this is the true distance normal from the front.  Of course, the 
!        distance is in a driection normal to the front.
               do i= 1, npro
                  if (sclr(i) .lt. - epsilon_ls)then
                     prop_blend(i) = zero
                  elseif  (abs(sclr(i)) .le. epsilon_ls)then
                     prop_blend(i) = 0.5*(one + Sclr(i)/epsilon_ls +
     &                    (sin(pi*Sclr(i)/epsilon_ls))/pi )
                  elseif (sclr(i) .gt. epsilon_ls) then
                     prop_blend(i) = one
                  endif
               enddo
               fact = datmat(1,1,2)/datmat(1,1,1)
c              rsrhol(:) = pres(:) / ( Rgas * T(:) )
c              call eqs(pres,T,rsrhol)
c              rsrhog(:) = fact*pres(:) / ( Rgas * T(:))
               rsrhog(:) = pres(:) / ( Rgas * T(:))
               rsrhol(:)  = datmat(1,1,1)*
     &              (1+0.000000000517992*(pres-18.02))
               rho(:)=rsrhol(:)*prop_blend(:)
     &              +rsrhog(:)*(1-prop_blend(:))
c ..... for the VOF case .. in case if we want to run VOF
c$$$         prop_blend(:) = min((max(sclr(:),0.0)),1.0)
c$$$         rho(:)=rsrhol(:) * prop_blend(:) + rsrhog(:) * (1-prop_blend(:))
c
        endif
c
          ei    = yN2 * ( cvs(1) * T
     &                   + Rs(1) * Tvib(1) * Texp1 / ( one - Texp1 ) )
     &          + yO2 * ( cvs(2) * T
     &                   + Rs(2) * Tvib(2) * Texp2 / ( one - Texp2 ) )
c
        endif
c
        if (ithm .ge. 7) then
c
c.... compute enthalpy, cp, alfap, betaT, cv, gamb and c
c
          h     = ei + Rgas * T
c
          alfap = one / T
c
          betaT = one / pres
c
          cp    = yN2 * ( cps(1) + Rs(1) * Tvib(1)**2 * Texp1
     &                   / ( ( one - Texp1 ) * T )**2 )
     &          + yO2 * ( cps(2) + Rs(2) * Tvib(2)**2 * Texp2
     &                   / ( ( one - Texp2 ) * T )**2 )
c
          cv    = cp - Rgas
c
          gamb  = Rgas / cv
c
          c     = sqrt( cp * gamb * T )
c
        endif
c
c
c.... end of ipress = 1
c
        endif

c       ttim(27) = ttim(27) + secs(0.0)
c
c.... end
c
        return
        end
