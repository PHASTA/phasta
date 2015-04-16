        subroutine e3massr (aci,    dui,  ri,     
     &                      rmi,    A0)
c
c----------------------------------------------------------------------
c
c This routine calculates the contribution of the jump condition 
c and the time flux to RHS and LHS.
c
c input:
c  aci     (npro,nflow)           : Y variable acceleration
c  dui     (npro,nflow)           : dU variables at previous step
c  A0      (npro,nflow,nflow)      : weighted Jacobian
c
c output:
c  ri     (npro,nflow*(nsd+1))   : partial residual
c  rmi    (npro,nflow*(nsd+1))   : partial modified residual
c
c
c
c Zdenek Johan, Summer 1990. (Modified from e2jump.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997 Prim Variables
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension aci(npro,nflow),            dui(npro,nflow),
     &            ri(npro,nflow*(nsd+1)),
     &            rmi(npro,nflow*(nsd+1)),    A0(npro,nflow,nflow)
c
c
c.... add contribution of U at previous step
c
       if(ires.eq.1 .or. ires .eq. 3) then

        ri(:,16) = ri(:,16) 
     &           + A0(:,1,1)*aci(:,1)
c    &           + A0(:,1,2)*aci(:,2)
c    &           + A0(:,1,3)*aci(:,3)
c    &           + A0(:,1,4)*aci(:,4)
     &           + A0(:,1,5)*aci(:,5)
c
        ri(:,17) = ri(:,17) 
     &           + A0(:,2,1)*aci(:,1)
     &           + A0(:,2,2)*aci(:,2)
c    &           + A0(:,2,3)*aci(:,3)
c    &           + A0(:,2,4)*aci(:,4)
     &           + A0(:,2,5)*aci(:,5)
c
        ri(:,18) = ri(:,18) 
     &           + A0(:,3,1)*aci(:,1)
c    &           + A0(:,3,2)*aci(:,2)
     &           + A0(:,3,3)*aci(:,3)
c    &           + A0(:,3,4)*aci(:,4)
     &           + A0(:,3,5)*aci(:,5)
c
        ri(:,19) = ri(:,19) 
     &           + A0(:,4,1)*aci(:,1)
c    &           + A0(:,4,2)*aci(:,2)
c    &           + A0(:,4,3)*aci(:,3)
     &           + A0(:,4,4)*aci(:,4)
     &           + A0(:,4,5)*aci(:,5)
c
        ri(:,20) = ri(:,20) 
     &           + A0(:,5,1)*aci(:,1)
     &           + A0(:,5,2)*aci(:,2)
     &           + A0(:,5,3)*aci(:,3)
     &           + A0(:,5,4)*aci(:,4)
     &           + A0(:,5,5)*aci(:,5)
c

         endif
         if(ires.ne.1) then
c
c the modified residual
c
            fct1=almi/gami/alfi*dtgl

            rmi(:,16) = rmi(:,16) + fct1*dui(:,1)
            rmi(:,17) = rmi(:,17) + fct1*dui(:,2)
            rmi(:,18) = rmi(:,18) + fct1*dui(:,3)
            rmi(:,19) = rmi(:,19) + fct1*dui(:,4)
            rmi(:,20) = rmi(:,20) + fct1*dui(:,5)
 
         endif

c
c.... return
c
        return
        end
c
c
c
        subroutine e3massrSclr (acti,   rti,     A0t)
c
c----------------------------------------------------------------------
c
c This routine calculates the contribution of the jump condition 
c and the time flux to RHS and LHS.
c
c input:
c  acti     (npro)           : scalar variable acceleration
c  A0t      (npro)           : weighted Jacobian
c
c output:
c  rti     (npro,nsd+1)   : partial residual
c
c
c
c Zdenek Johan, Summer 1990. (Modified from e2jump.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997 Prim Variables
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension acti(npro),
     &            rti(npro,nsd+1),
     &            rmti(npro,nsd+1),    A0t(npro)
c
c
c.... add contribution of U at previous step
c
       if(ires.eq.1 .or. ires .eq. 3) then

        rti(:,4) = rti(:,4) + A0t(:)*acti(:)

         endif


c         if(ires.ne.1) then
c
c the modified residual
c
c            fct1=almi/gami/alfi*dtgl
c
c            rmi(:,4) = rmi(:,4) + fct1*duti(:)
c 
c         endif

c
c.... return
c
        return
        end

