        subroutine getFld (T,      cp,     rmu,    rlm,    rlm2mu,
     &                     con)
c
c----------------------------------------------------------------------
c
c This routine calculates the fluid material properties.
c
c input:
c  T      (npro)        : temperature
c  cp     (npro)        : specific heat at constant pressure
c
c output:
c  rmu    (npro)        : Mu
c  rlm    (npro)        : Lambda
c  rlm2mu (npro)        : Lambda + 2 Mu
c  con    (npro)        : Conductivity
c
c Note: material type flags
c         matflg(2):
c          eq. 0, constant viscosity
c          eq. 1, generalized Sutherland viscosity
c         matflg(3):
c          eq. 0, Stokes approximation
c          eq. 1, shear proportional bulk viscosity
c
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension T(npro),                   cp(npro),
     &            rmu(npro),                 rlm(npro),
     &            rlm2mu(npro),              con(npro)
c
c
c.... constant viscosity
c
        if (matflg(2,1) .eq. 0) then
c
          rmu = datmat(1,2,1)
c
        else
c
c.... generalized Sutherland viscosity
c
          rmu = datmat(1,2,1) * (T/datmat(2,2,1))*sqrt(T/datmat(2,2,1))
     &        * ( datmat(2,2,1) + datmat(3,2,1) ) / (T + datmat(3,2,1))
c
        endif
c
c.... calculate the second viscosity coefficient
c
        if (matflg(3,1) .eq. 0) then
c
          rlm = -pt66 * rmu
c
        else
c
          rlm = (datmat(1,3,1) - pt66) * rmu
c
        endif
c
c.... calculate the remaining quantities
c
        cp     = datmat(1,3,1)
        rlm2mu = rlm + two * rmu
        con    = datmat(1,4,1)
c
c.... return
c
        return
        end
