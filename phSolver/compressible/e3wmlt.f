        subroutine e3wmlt (shp,    shg,    WdetJ,  
     &                     ri,     rmi,    rl,
     &                     rml,    stiff,  EGmass )
c
c----------------------------------------------------------------------
c
c Up to now most of the terms have not been multiplied by the 
c shape function from the weight space.  Now that we have collected
c all the terms that the shape function (and its derivatives for the
c terms that were integrated by parts) multiplies, in this routine
c we carry out that multiplication. After these operations we will
c have this quadrature points contribution to the integral for the
c residual (i.e.  g^e_b(xi^l) D(xi^l) QW(xi^l)  where e is the element,
c b is the local node number on the element, l is the quadrature point
c D is the determinate of the Jacobian of the mapping, and QW is this
c quadrature points weight....WHEW).  When we add up all of the 
c integration points we get G^e_b which we will assemble in the 
c subroutine LOCAL to form G_B.
c
c This routine also forms the LHS contribution from the LS term
c and the diffusion term which has been partially constructed and
c placed in stiff.  Those familiar with elasticity might recognize
c this naming convention since this is like a stiffness matrix that
c if you had a linear problem would be calculated once and saved for
c all time.
c
c    ri:  LS, Diffusion, Convective and mass, 
c   rmi:  LS, Diffusion and mass, 
c stiff:  LS, Diffusion.
c
c input:
c  shp    (nshl)                 : element shape-functions
c  shg    (npro,nshl,nsd)        : element grad-shape-functions
c  WdetJ  (npro)                   : weighted Jacobian
c  ri     (npro,nflow*(nsd+1))      : partial residual
c  rmi    (npro,nflow*(nsd+1))      : partial modified residual
c  stiff  (npro,nsd*nflow,nsd*nflow) : stiffness matrix 
c                                    ( K_ij + A_i tau A_j )
c
c output:
c  rl     (npro,nshl,nflow)       : residual
c  rml    (npro,nshl,nflow)       : modified residual
c  EGmass (npro,nedof,nedof)       : element LHS tangent mass matrix
c
c
c Zdenek Johan, Summer 1990.  (Modified from e2assm.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c Kenneth Jansen, Winter 1997  Prim. variables
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension shp(npro,nshl),  
     &            shg(npro,nshl,nsd),
     &            WdetJ(npro),             ri(npro,nflow*(nsd+1)),
     &            rmi(npro,nflow*(nsd+1)), stiff(npro,3*nflow,3*nflow),
     &            rl(npro,nshl,nflow),     rml(npro,nshl,nflow),
     &            EGmass(npro,nedof,nedof)
c
        dimension shg1(npro),              shg2(npro),
     &            shg3(npro),              stif1(npro,nflow,nflow),
     &            stif2(npro,nflow,nflow), stif3(npro,nflow,nflow)

        ttim(29) = ttim(29) - secs(0.0)
c
c.... ---------------------------->  RHS  <----------------------------
c
c.... add spatial contribution to rl and rml
c
c.... ires = 1 or 3
c
        if ((ires .eq. 1) .or. (ires .eq. 3)) then
c
          do i = 1, nshl
            rl(:,i,1) = rl(:,i,1) + WdetJ * (
     &                  shg(:,i,1) * ri(:, 1) + shg(:,i,2) * ri(:, 6)
     &                                        + shg(:,i,3) * ri(:,11) )
            rl(:,i,2) = rl(:,i,2) + WdetJ * (
     &                  shg(:,i,1) * ri(:, 2) + shg(:,i,2) * ri(:, 7)
     &                                        + shg(:,i,3) * ri(:,12) )
            rl(:,i,3) = rl(:,i,3) + WdetJ * (
     &                  shg(:,i,1) * ri(:, 3) + shg(:,i,2) * ri(:, 8)
     &                                        + shg(:,i,3) * ri(:,13) )
            rl(:,i,4) = rl(:,i,4) + WdetJ * (
     &                  shg(:,i,1) * ri(:, 4) + shg(:,i,2) * ri(:, 9)
     &                                        + shg(:,i,3) * ri(:,14) )
            rl(:,i,5) = rl(:,i,5) + WdetJ * (
     &                  shg(:,i,1) * ri(:, 5) + shg(:,i,2) * ri(:,10)
     &                                        + shg(:,i,3) * ri(:,15) )
          enddo
c
!      flops = flops + 36*nshl*npro
        endif
c
c.... ires = 2 or 3
c
        if ((ires .eq. 2) .or. (ires .eq. 3)) then
          do i = 1, nshl
            rml(:,i,1) = rml(:,i,1) + WdetJ * (
     &                   shg(:,i,1) * rmi(:, 1) + shg(:,i,2) * rmi(:, 6)
     &                 + shg(:,i,3) * rmi(:,11) 
     &                   + shp(:,i) * rmi(:,16)  )
            rml(:,i,2) = rml(:,i,2) + WdetJ * (
     &                   shg(:,i,1) * rmi(:, 2) + shg(:,i,2) * rmi(:, 7)
     &                 + shg(:,i,3) * rmi(:,12) 
     &                   + shp(:,i) * rmi(:,17)    )
            rml(:,i,3) = rml(:,i,3) + WdetJ * (
     &                   shg(:,i,1) * rmi(:, 3) + shg(:,i,2) * rmi(:, 8)
     &                 + shg(:,i,3) * rmi(:,13)
     &                   + shp(:,i) * rmi(:,18)    )
            rml(:,i,4) = rml(:,i,4) + WdetJ * (
     &                   shg(:,i,1) * rmi(:, 4) + shg(:,i,2) * rmi(:, 9)
     &                 + shg(:,i,3) * rmi(:,14) 
     &                   + shp(:,i) * rmi(:,19)    )
            rml(:,i,5) = rml(:,i,5) + WdetJ * (
     &                   shg(:,i,1) * rmi(:, 5) + shg(:,i,2) * rmi(:,10)
     &                 + shg(:,i,3) * rmi(:,15)
     &                   + shp(:,i) * rmi(:,20)    )
          enddo
c
!      flops = flops + (15+45*nshl)*npro
        endif
c
c.... add temporal contribution to rl
c
        if (ngauss .eq. 1 .and. nshl.eq.4) then !already ex. integ mass
         if(matflg(5,1).ge.1) then
          do i = 1, nshl
            rl(:,i,2) = rl(:,i,2) + shp(:,i) * WdetJ * ri(:,17)
            rl(:,i,3) = rl(:,i,3) + shp(:,i) * WdetJ * ri(:,18)
            rl(:,i,4) = rl(:,i,4) + shp(:,i) * WdetJ * ri(:,19)
            rl(:,i,5) = rl(:,i,5) + shp(:,i) * WdetJ * ri(:,20)
          enddo
         endif
       else   !check for a body force
          do i = 1, nshl
            rl(:,i,1) = rl(:,i,1) + shp(:,i) * WdetJ * ri(:,16)  
            rl(:,i,2) = rl(:,i,2) + shp(:,i) * WdetJ * ri(:,17)  
            rl(:,i,3) = rl(:,i,3) + shp(:,i) * WdetJ * ri(:,18)  
            rl(:,i,4) = rl(:,i,4) + shp(:,i) * WdetJ * ri(:,19)  
            rl(:,i,5) = rl(:,i,5) + shp(:,i) * WdetJ * ri(:,20)  
          enddo
 
!      flops = flops + 11*nshl*npro
       endif

c
c.... ---------------------------->  LHS  <----------------------------
c
       if (lhs .eq. 1) then
c
c.... loop through columns (nodes j)
c
          do j = 1, nshl
             j0 = nflow * (j - 1)
c
c.... set up factors
c
             shg1 = WdetJ * shg(:,j,1)
             shg2 = WdetJ * shg(:,j,2)
             shg3 = WdetJ * shg(:,j,3)
c
c.... loop through d.o.f.'s
c
             do jdof = 1, nflow
                do idof = 1, nflow
                   idof2 = idof  + nflow
                   jdof2 = jdof  + nflow

                   idof3 = idof2 + nflow
                   jdof3 = jdof2 + nflow
c
c.... calculate weighted stiffness matrix (first part)
c
                   stif1(:,idof,jdof) = shg1 * stiff(:,idof,jdof)
     &                                + shg2 * stiff(:,idof,jdof2)
     &                                + shg3 * stiff(:,idof,jdof3)
                   stif2(:,idof,jdof) = shg1 * stiff(:,idof2,jdof)
     &                                + shg2 * stiff(:,idof2,jdof2)
     &                                + shg3 * stiff(:,idof2,jdof3)
                   stif3(:,idof,jdof) = shg1 * stiff(:,idof3,jdof)
     &                                + shg2 * stiff(:,idof3,jdof2)
     &                                + shg3 * stiff(:,idof3,jdof3)
                enddo
             enddo
c
c.... loop through rows (nodes i)
c
             do i = 1, nshl
                i0 = nflow * (i - 1)
c
c.... add contribution of stiffness to EGmass
c
                do jdof = 1, nflow
                   EGmass(:,i0+1,j0+jdof) = EGmass(:,i0+1,j0+jdof) 
     &                  + shg(:,i,1) * stif1(:,1,jdof)
     &                  + shg(:,i,2) * stif2(:,1,jdof)
     &                  + shg(:,i,3) * stif3(:,1,jdof)
                   EGmass(:,i0+2,j0+jdof) = EGmass(:,i0+2,j0+jdof) 
     &                  + shg(:,i,1) * stif1(:,2,jdof)
     &                  + shg(:,i,2) * stif2(:,2,jdof)
     &                  + shg(:,i,3) * stif3(:,2,jdof)
                   EGmass(:,i0+3,j0+jdof) = EGmass(:,i0+3,j0+jdof) 
     &                  + shg(:,i,1) * stif1(:,3,jdof)
     &                  + shg(:,i,2) * stif2(:,3,jdof)
     &                  + shg(:,i,3) * stif3(:,3,jdof)
                   EGmass(:,i0+4,j0+jdof) = EGmass(:,i0+4,j0+jdof) 
     &                  + shg(:,i,1) * stif1(:,4,jdof)
     &                  + shg(:,i,2) * stif2(:,4,jdof)
     &                  + shg(:,i,3) * stif3(:,4,jdof)
                   EGmass(:,i0+5,j0+jdof) = EGmass(:,i0+5,j0+jdof) 
     &                  + shg(:,i,1) * stif1(:,5,jdof)
     &                  + shg(:,i,2) * stif2(:,5,jdof)
     &                  + shg(:,i,3) * stif3(:,5,jdof)
                enddo
c
c.... end loop on rows
c
             enddo
c
c.... end loop on columns
c
          enddo
c
c.... end left hand side assembly
c          
       endif
                
       ttim(29) = ttim(29) + secs(0.0)
c
c.... return
c
        return
        end
c
c
c
        subroutine e3wmltSclr(shp,    shg,    WdetJ,  
     &                        rti,    rtl,
     &                        stifft, EGmasst )
c
c----------------------------------------------------------------------
c
c Up to now most of the terms have not been multiplied by the 
c shape function from the weight space.  Now that we have collected
c all the terms that the shape function (and its derivatives for the
c terms that were integrated by parts) multiplies, in this routine
c we carry out that multiplication. After these operations we will
c have this quadrature points contribution to the integral for the
c residual (i.e.  g^e_b(xi^l) D(xi^l) QW(xi^l)  where e is the element,
c b is the local node number on the element, l is the quadrature point
c D is the determinate of the Jacobian of the mapping, and QW is this
c quadrature points weight....WHEW).  When we add up all of the 
c integration points we get G^e_b which we will assemble in the 
c subroutine LOCAL to form G_B.
c
c This routine also forms the LHS contribution from the LS term
c and the diffusion term which has been partially constructed and
c placed in stiff.  Those familiar with elasticity might recognize
c this naming convention since this is like a stiffness matrix that
c if you had a linear problem would be calculated once and saved for
c all time.
c
c    ri:  LS, Diffusion, Convective and mass, 
c stiff:  LS, Diffusion.
c
c input:
c  shp     (npro,nshl)            : element shape-functions
c  shg     (npro,nshl,nsd)        : element grad-shape-functions
c  WdetJ   (npro)                 : weighted Jacobian
c  rti     (npro,nsd+1)           : partial residual
c  stifft  (npro,nsd,nsd)         : stiffness matrix 
c                                    ( K_ij + A_i tau A_j )
c
c output:
c  rtl     (npro,nshl)            : residual  (will end up being G^e_b)
c  EGmasst (npro,nshape,nshape)   : element LHS tangent mass matrix  
c                                  (will end up being dG^e_b/dY_a)
c
c
c Zdenek Johan, Summer 1990.  (Modified from e2assm.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c Kenneth Jansen, Winter 1997  Prim. variables
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension shp(npro,nshl),            shg(npro,nshl,nsd),
     &            WdetJ(npro),               rti(npro,nsd+1),
     &            rmti(npro,nsd+1),          stifft(npro,nsd,nsd),
     &            rtl(npro,nshl),            rmtl(npro,nshl),
     &            EGmasst(npro,nshape,nshape)
c
        dimension shg1(npro),                shg2(npro),
     &            shg3(npro),                stift1(npro),
     &            stift2(npro),              stift3(npro)

        ttim(29) = ttim(29) - tmr(0.0)
c
c.... ---------------------------->  RHS  <----------------------------
c
c.... add spatial contribution to rl and rml
c
c.... ires = 1 or 3
c
        if ((ires .eq. 1) .or. (ires .eq. 3)) then
c
          do i = 1, nshl
            rtl(:,i) = rtl(:,i) + WdetJ * (
     &                  shg(:,i,1) * rti(:,1) + shg(:,i,2) * rti(:,2)
     &                                        + shg(:,i,3) * rti(:,3) )

          enddo
!      flops = flops + 36*nshl*npro
        endif
c
c.... ires = 2 or 3
c
c        if ((ires .eq. 2) .or. (ires .eq. 3)) then
c          do i = 1, nshl
c            rml(:,i,1) = rml(:,i,1) + WdetJ * (
c     &                   shg(:,i,1) * rmi(:, 1) + shg(:,i,2) * rmi(:, 6)
c     &                 + shg(:,i,3) * rmi(:,11) + shp(:,i) * rmi(:,16) )
c            rml(:,i,2) = rml(:,i,2) + WdetJ * (
c     &                   shg(:,i,1) * rmi(:, 2) + shg(:,i,2) * rmi(:, 7)
c     &                 + shg(:,i,3) * rmi(:,12) + shp(:,i) * rmi(:,17) )
c            rml(:,i,3) = rml(:,i,3) + WdetJ * (
c     &                   shg(:,i,1) * rmi(:, 3) + shg(:,i,2) * rmi(:, 8)
c     &                 + shg(:,i,3) * rmi(:,13) + shp(:,i) * rmi(:,18) )
c            rml(:,i,4) = rml(:,i,4) + WdetJ * (
c     &                   shg(:,i,1) * rmi(:, 4) + shg(:,i,2) * rmi(:, 9)
c     &                 + shg(:,i,3) * rmi(:,14) + shp(:,i) * rmi(:,19) )
c            rml(:,i,5) = rml(:,i,5) + WdetJ * (
c     &                   shg(:,i,1) * rmi(:, 5) + shg(:,i,2) * rmi(:,10)
c     &                 + shg(:,i,3) * rmi(:,15) + shp(:,i) * rmi(:,20) )
c          enddo
c
c     !      flops = flops + (15+45*nshl)*npro
c        endif

c
c.... add temporal contribution to rl
c

          do i = 1, nshl
            rtl(:,i) = rtl(:,i) + shp(:,i) * WdetJ * rti(:,4)
          enddo
 
!      flops = flops + 11*nshl*npro

c
c.... ---------------------------->  LHS  <----------------------------
c
       if (lhs .eq. 1) then
c
c.... loop through columns (nodes j)
c
          do j = 1, nshl
c
c.... set up factors
c
             shg1 = WdetJ * shg(:,j,1)
             shg2 = WdetJ * shg(:,j,2)
             shg3 = WdetJ * shg(:,j,3)
c
c.... calculate weighted stiffness matrix (first part)
c
                   stift1(:) = shg1 * stifft(:,1,1)
     &                                + shg2 * stifft(:,1,2)
     &                                + shg3 * stifft(:,1,3)
                   stift2(:) = shg1 * stifft(:,2,1)
     &                                + shg2 * stifft(:,2,2)
     &                                + shg3 * stifft(:,2,3)
                   stift3(:) = shg1 * stifft(:,3,1)
     &                                + shg2 * stifft(:,3,2)
     &                                + shg3 * stifft(:,3,3)
c
c.... loop through rows (nodes i)
c
             do i = 1, nshl
c
c.... add contribution of stiffness to EGmass
c
                   EGmasst(:,i,j) = EGmasst(:,i,j) 
     &                  + shg(:,i,1) * stift1(:)
     &                  + shg(:,i,2) * stift2(:)
     &                  + shg(:,i,3) * stift3(:)

c
c.... end loop on rows
c
             enddo
c
c.... end loop on columns
c
          enddo
c
c.... end left hand side assembly
c          
       endif
                
       ttim(29) = ttim(29) + tmr(0.0)
c
c.... return
c
        return
        end
