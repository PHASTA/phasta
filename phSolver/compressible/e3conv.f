        subroutine e3conv (g1yi,   g2yi,   g3yi,
     &                     A1,     A2,     A3, 
     &                     rho,    pres,   T,
     &                     ei,     rk,     u1,
     &                     u2,     u3,     rLyi,   
     &                     ri,     rmi,    EGmass,
     &                     shg,    shape,  WdetJ )
!
!----------------------------------------------------------------------
!
! This routine calculates the contribution of Galerkin part of the 
! Convective term (Time and Euler fluxes) to both RHS and LHS.
!
! input:
!  g1yi   (npro,nflow)    : grad-y in direction 1
!  g2yi   (npro,nflow)    : grad-y in direction 2
!  g3yi   (npro,nflow)    : grad-y in direction 3
!  A1    (npro,nflow,nflow)  : A-1
!  A2    (npro,nflow,nflow)  : A-2
!  A3    (npro,nflow,nflow)  : A-3 
!  rho    (npro)         : density
!  pres   (npro)         : pressure
!  T      (npro)         : temperature
!  ei     (npro)         : internal energy
!  rk     (npro)         : kinetic energy
!  u1     (npro)         : x1-velocity component
!  u2     (npro)         : x2-velocity component
!  u3     (npro)         : x3-velocity component
!  shg    (npro,nshl,nsd) : global grad's of shape functions
!  shape  (npro,nshl)  : element shape functions      
!  WdetJ  (npro)         : weighted Jacobian determinant
!     
! output:
!  rLyi   (npro,nflow)           : least-squares residual vector
!  ri     (npro,nflow*(nsd+1))   : partial residual
!  rmi    (npro,nflow*(nsd+1))   : partial modified residual
!  EGmass (npro,nedof,nedof)    : partial LHS tangent matrix
!
!
! Zdenek Johan, Summer 1990. (Modified from e2conv.f)
! Zdenek Johan, Winter 1991. (Fortran 90)
! Kenneth Jansen, Winter 1997 Primitive Variables
!----------------------------------------------------------------------
!
        include "common.h"
!
!  passed arrays
!
        dimension g1yi(npro,nflow),           g2yi(npro,nflow),
     &            g3yi(npro,nflow),           
     &            A1(npro,nflow,nflow),
     &            A2(npro,nflow,nflow),       A3(npro,nflow,nflow),
     &            rho(npro),                pres(npro),
     &            T(npro),                  ei(npro),
     &            rk(npro),                 u1(npro),
     &            u2(npro),                 u3(npro),
     &            rLyi(npro,nflow),          ri(npro,nflow*(nsd+1)),
     &            rmi(npro,nflow*(nsd+1)),   EGmass(npro,nedof,nedof),
     &            shg(npro,nshl,nsd),       shape(npro,nshl),
     &            WdetJ(npro)
!
!  local arrays  
!
        dimension AiNbi(npro,nflow,nflow),     fact1(npro),
     &            fact2(npro),               fact3(npro)
        
        ttim(22) = ttim(22) - secs(0.0)

!
!.... ---------------------->  RHS, Euler Flux  <----------------------
!
        if ((ires .eq. 1) .or. (ires .eq. 3)) then
!
!.... calculate integrated by part contribution of Euler flux (Galerkin)
!
          ri(:, 1) = (- u1) * rho
          ri(:, 2) = (- u1) * rho * u1 - pres
          ri(:, 3) = (- u1) * rho * u2
          ri(:, 4) = (- u1) * rho * u3
          ri(:, 5) = (- u1) * rho * (ei + rk) - u1 * pres
!
          ri(:, 6) = (- u2) * rho
          ri(:, 7) = (- u2) * rho * u1
          ri(:, 8) = (- u2) * rho * u2 - pres
          ri(:, 9) = (- u2) * rho * u3
          ri(:,10) = (- u2) * rho * (ei + rk) - u2 * pres
!
          ri(:,11) = (- u3) * rho
          ri(:,12) = (- u3) * rho * u1
          ri(:,13) = (- u3) * rho * u2
          ri(:,14) = (- u3) * rho * u3 - pres
          ri(:,15) = (- u3) * rho * (ei + rk) - u3 * pres
!
!      flops = flops + 28*npro
!
        endif
!
!.... calculate ( A_i Y,i ) --> rLyi   Commented out zeros of A matrices
!
        rLyi(:,1) = 
     &              A1(:,1,1) * g1yi(:,1) 
     &            + A1(:,1,2) * g1yi(:,2)
!    &            + A1(:,1,3) * g1yi(:,3) 
!    &            + A1(:,1,4) * g1yi(:,4)
     &            + A1(:,1,5) * g1yi(:,5)
     &            + A2(:,1,1) * g2yi(:,1) 
!    &            + A2(:,1,2) * g2yi(:,2)
     &            + A2(:,1,3) * g2yi(:,3) 
!    &            + A2(:,1,4) * g2yi(:,4)
     &            + A2(:,1,5) * g2yi(:,5)
     &            + A3(:,1,1) * g3yi(:,1) 
!    &            + A3(:,1,2) * g3yi(:,2)
!    &            + A3(:,1,3) * g3yi(:,3) 
     &            + A3(:,1,4) * g3yi(:,4)
     &            + A3(:,1,5) * g3yi(:,5)
        rLyi(:,2) = 
     &              A1(:,2,1) * g1yi(:,1) 
     &            + A1(:,2,2) * g1yi(:,2)
!    &            + A1(:,2,3) * g1yi(:,3) 
!    &            + A1(:,2,4) * g1yi(:,4)
     &            + A1(:,2,5) * g1yi(:,5)
     &            + A2(:,2,1) * g2yi(:,1) 
     &            + A2(:,2,2) * g2yi(:,2)
     &            + A2(:,2,3) * g2yi(:,3) 
!    &            + A2(:,2,4) * g2yi(:,4)
     &            + A2(:,2,5) * g2yi(:,5)
     &            + A3(:,2,1) * g3yi(:,1) 
     &            + A3(:,2,2) * g3yi(:,2)
!    &            + A3(:,2,3) * g3yi(:,3) 
     &            + A3(:,2,4) * g3yi(:,4)
     &            + A3(:,2,5) * g3yi(:,5)
        rLyi(:,3) = 
     &              A1(:,3,1) * g1yi(:,1) 
     &            + A1(:,3,2) * g1yi(:,2)
     &            + A1(:,3,3) * g1yi(:,3) 
!    &            + A1(:,3,4) * g1yi(:,4)
     &            + A1(:,3,5) * g1yi(:,5)
     &            + A2(:,3,1) * g2yi(:,1) 
!    &            + A2(:,3,2) * g2yi(:,2)
     &            + A2(:,3,3) * g2yi(:,3) 
!    &            + A2(:,3,4) * g2yi(:,4)
     &            + A2(:,3,5) * g2yi(:,5)
     &            + A3(:,3,1) * g3yi(:,1) 
!    &            + A3(:,3,2) * g3yi(:,2)
     &            + A3(:,3,3) * g3yi(:,3) 
     &            + A3(:,3,4) * g3yi(:,4)
     &            + A3(:,3,5) * g3yi(:,5)
        rLyi(:,4) = 
     &              A1(:,4,1) * g1yi(:,1) 
     &            + A1(:,4,2) * g1yi(:,2)
!    &            + A1(:,4,3) * g1yi(:,3) 
     &            + A1(:,4,4) * g1yi(:,4)
     &            + A1(:,4,5) * g1yi(:,5)
     &            + A2(:,4,1) * g2yi(:,1) 
!    &            + A2(:,4,2) * g2yi(:,2)
     &            + A2(:,4,3) * g2yi(:,3) 
     &            + A2(:,4,4) * g2yi(:,4)
     &            + A2(:,4,5) * g2yi(:,5)
     &            + A3(:,4,1) * g3yi(:,1) 
!    &            + A3(:,4,2) * g3yi(:,2)
!    &            + A3(:,4,3) * g3yi(:,3) 
     &            + A3(:,4,4) * g3yi(:,4)
     &            + A3(:,4,5) * g3yi(:,5)
        rLyi(:,5) = 
     &              A1(:,5,1) * g1yi(:,1) 
     &            + A1(:,5,2) * g1yi(:,2)
     &            + A1(:,5,3) * g1yi(:,3) 
     &            + A1(:,5,4) * g1yi(:,4)
     &            + A1(:,5,5) * g1yi(:,5)
     &            + A2(:,5,1) * g2yi(:,1) 
     &            + A2(:,5,2) * g2yi(:,2)
     &            + A2(:,5,3) * g2yi(:,3) 
     &            + A2(:,5,4) * g2yi(:,4)
     &            + A2(:,5,5) * g2yi(:,5)
     &            + A3(:,5,1) * g3yi(:,1) 
     &            + A3(:,5,2) * g3yi(:,2)
     &            + A3(:,5,3) * g3yi(:,3) 
     &            + A3(:,5,4) * g3yi(:,4)
     &            + A3(:,5,5) * g3yi(:,5)
!
!.... add contribution to rmi
!
        if ((ires .eq. 2) .or. (ires .eq. 3))
     &    rmi(:,16:20) = rLyi  ! modified residual uses non i.b.p form of conv.
!
!.... ---------------------->  LHS   <-----------------------
!
        if (lhs .eq. 1) then
!
!.... loop through the columns
!
        do j = 1, nshl
           j0 = nflow * (j - 1)
!
!.... compute some useful factors
!
           fact1 = WdetJ * shg(:,j,1)
           fact2 = WdetJ * shg(:,j,2)
           fact3 = WdetJ * shg(:,j,3)
!
!.... first compute (A_i N_b,i)
!
           AiNbi(:,1,1) = 
     &                    fact1 * A1(:,1,1) 
     &                  + fact2 * A2(:,1,1) 
     &                  + fact3 * A3(:,1,1)
           AiNbi(:,1,2) = 
     &                    fact1 * A1(:,1,2) 
     &                  + fact2 * A2(:,1,2) 
     &                  + fact3 * A3(:,1,2)
           AiNbi(:,1,3) = 
     &                    fact1 * A1(:,1,3) 
     &                  + fact2 * A2(:,1,3) 
     &                  + fact3 * A3(:,1,3)
           AiNbi(:,1,4) = 
     &                    fact1 * A1(:,1,4) 
     &                  + fact2 * A2(:,1,4) 
     &                  + fact3 * A3(:,1,4)
           AiNbi(:,1,5) = 
     &                    fact1 * A1(:,1,5) 
     &                  + fact2 * A2(:,1,5) 
     &                  + fact3 * A3(:,1,5)

           AiNbi(:,2,1) = 
     &                    fact1 * A1(:,2,1) 
     &                  + fact2 * A2(:,2,1) 
     &                  + fact3 * A3(:,2,1)
           AiNbi(:,2,2) = 
     &                    fact1 * A1(:,2,2) 
     &                  + fact2 * A2(:,2,2) 
     &                  + fact3 * A3(:,2,2)
           AiNbi(:,2,3) = 
     &                    fact1 * A1(:,2,3) 
     &                  + fact2 * A2(:,2,3) 
     &                  + fact3 * A3(:,2,3)
           AiNbi(:,2,4) = 
     &                    fact1 * A1(:,2,4) 
     &                  + fact2 * A2(:,2,4) 
     &                  + fact3 * A3(:,2,4)
           AiNbi(:,2,5) = 
     &                    fact1 * A1(:,2,5) 
     &                  + fact2 * A2(:,2,5) 
     &                  + fact3 * A3(:,2,5)

           AiNbi(:,3,1) = 
     &                    fact1 * A1(:,3,1) 
     &                  + fact2 * A2(:,3,1) 
     &                  + fact3 * A3(:,3,1)
           AiNbi(:,3,2) = 
     &                    fact1 * A1(:,3,2) 
     &                  + fact2 * A2(:,3,2) 
     &                  + fact3 * A3(:,3,2)
           AiNbi(:,3,3) = 
     &                    fact1 * A1(:,3,3) 
     &                  + fact2 * A2(:,3,3) 
     &                  + fact3 * A3(:,3,3)
           AiNbi(:,3,4) = 
     &                    fact1 * A1(:,3,4) 
     &                  + fact2 * A2(:,3,4) 
     &                  + fact3 * A3(:,3,4)
           AiNbi(:,3,5) = 
     &                    fact1 * A1(:,3,5) 
     &                  + fact2 * A2(:,3,5) 
     &                  + fact3 * A3(:,3,5)

           AiNbi(:,4,1) = 
     &                    fact1 * A1(:,4,1) 
     &                  + fact2 * A2(:,4,1) 
     &                  + fact3 * A3(:,4,1)
           AiNbi(:,4,2) = 
     &                    fact1 * A1(:,4,2) 
     &                  + fact2 * A2(:,4,2) 
     &                  + fact3 * A3(:,4,2)
           AiNbi(:,4,3) = 
     &                    fact1 * A1(:,4,3) 
     &                  + fact2 * A2(:,4,3) 
     &                  + fact3 * A3(:,4,3)
           AiNbi(:,4,4) = 
     &                    fact1 * A1(:,4,4) 
     &                  + fact2 * A2(:,4,4) 
     &                  + fact3 * A3(:,4,4)
           AiNbi(:,4,5) = 
     &                    fact1 * A1(:,4,5) 
     &                  + fact2 * A2(:,4,5) 
     &                  + fact3 * A3(:,4,5)

           AiNbi(:,5,1) = 
     &                    fact1 * A1(:,5,1) 
     &                  + fact2 * A2(:,5,1) 
     &                  + fact3 * A3(:,5,1)
           AiNbi(:,5,2) = 
     &                    fact1 * A1(:,5,2) 
     &                  + fact2 * A2(:,5,2) 
     &                  + fact3 * A3(:,5,2)
           AiNbi(:,5,3) = 
     &                    fact1 * A1(:,5,3) 
     &                  + fact2 * A2(:,5,3) 
     &                  + fact3 * A3(:,5,3)
           AiNbi(:,5,4) = 
     &                    fact1 * A1(:,5,4) 
     &                  + fact2 * A2(:,5,4) 
     &                  + fact3 * A3(:,5,4)
           AiNbi(:,5,5) = 
     &                    fact1 * A1(:,5,5) 
     &                  + fact2 * A2(:,5,5) 
     &                  + fact3 * A3(:,5,5)
!
!.... now loop through the row nodes and add (N_a A_i N_b,i) to
!     the tangent matrix.
!
           do i = 1, nshl
              i0 = nflow * (i - 1)
!
!.... loop through dof's
!
              do jdof = 1, nflow
                 jl = j0 + jdof
                 
                 EGmass(:,i0+1,jl) = EGmass(:,i0+1,jl) + 
     &                               shape(:,i) * AiNbi(:,1,jdof)

                 EGmass(:,i0+2,jl) = EGmass(:,i0+2,jl) + 
     &                               shape(:,i) * AiNbi(:,2,jdof)

                 EGmass(:,i0+3,jl) = EGmass(:,i0+3,jl) + 
     &                               shape(:,i) * AiNbi(:,3,jdof)

                 EGmass(:,i0+4,jl) = EGmass(:,i0+4,jl) + 
     &                               shape(:,i) * AiNbi(:,4,jdof)

                 EGmass(:,i0+5,jl) = EGmass(:,i0+5,jl) + 
     &                               shape(:,i) * AiNbi(:,5,jdof)
              enddo
!
!.... end loop on rows
!              
           enddo
!
!.... end loop on columns
!
        enddo
!
!.... end of LHS tangent matrix computation
!        
        endif
      
        ttim(22) = ttim(22) + secs(0.0)
!
!.... return
!
        return
        end
!
!
!
        subroutine e3convSclr (g1yti,   g2yti,   g3yti,
     &                         A1t,     A2t,     A3t, 
     &                         rho,     u1,      Sclr,
     &                         u2,      u3,      rLyti,   
     &                         rti,     rmti,    EGmasst,
     &                         shg,     shape,   WdetJ)
!
!----------------------------------------------------------------------
!
! This routine calculates the contribution of Galerkin part of the 
! Convective term (Time and Euler fluxes) to both RHS and LHS.
!
! input:
!  Sclr   (npro)          : Scalar variable
!  g1yti  (npro)          : grad-y in direction 1
!  g2yti  (npro)          : grad-y in direction 2
!  g3yti  (npro)          : grad-y in direction 3
!  A1t    (npro)          : A-1
!  A2t    (npro)          : A-2
!  A3t    (npro)          : A-3 
!  rho    (npro)          : density
!  u1     (npro)          : x1-velocity component
!  u2     (npro)          : x2-velocity component
!  u3     (npro)          : x3-velocity component
!  shg    (npro,nshl,nsd) : global grad's of shape functions
!  shape  (npro,nshl)     : element shape functions      
!  WdetJ  (npro)          : weighted Jacobian determinant
!     
! output:
!  rLyti   (npro)         : least-squares residual vector
!  rti     (npro,nsd+1)   : partial residual
!  rmti    (npro,nsd+1)   : partial modified residual
!  EGmasst (npro,nshape,nshape): partial LHS tangent matrix
!
!
! Zdenek Johan, Summer 1990. (Modified from e2conv.f)
! Zdenek Johan, Winter 1991. (Fortran 90)
! Kenneth Jansen, Winter 1997 Primitive Variables
!----------------------------------------------------------------------
!
        include "common.h"
!
!  passed arrays
!
        dimension g1yti(npro),            g2yti(npro),
     &            g3yti(npro),            Sclr(npro),
     &            A1t(npro),             
     &            A2t(npro),              A3t(npro),
     &            rho(npro),              u1(npro),
     &            u2(npro),               u3(npro),
     &            rLyti(npro),            rti(npro,nsd+1),
     &            rmti(npro,nsd+1),       EGmasst(npro,nshape,nshape),
     &            shg(npro,nshl,nsd),     shape(npro,nshl),
     &            WdetJ(npro)
!
!  local arrays  
!
        dimension AitNbi(npro)
       
        ttim(22) = ttim(22) - tmr() !tmr(0.0)
!
!.... ---------------------->  RHS, Euler Flux  <----------------------
!
        if ((ires .eq. 1) .or. (ires .eq. 3)) then
!
!.... calculate integrated by part contribution of Euler flux (Galerkin)
!     
           if (iconvsclr.eq.2) then ! convective form
!     
              rti(:, 4) = rti(:,4) + ( u1) * g1yti(:)
     &                             + ( u2) * g2yti(:)
     &                             + ( u3) * g3yti(:) 
!     
           else                 ! conservative form     
!     
              rti(:, 1) = rti(:,1) + (- u1) * rho * Sclr
              rti(:, 2) = rti(:,2) + (- u2) * rho * Sclr
              rti(:, 3) = rti(:,3) + (- u3) * rho * Sclr
!     
           endif

!      flops = flops + 28*npro

        endif
!
!.... calculate ( A_i Y,i ) --> rLyi
!
        rLyti(:) = rLyti(:)
     &            + A1t(:) * g1yti(:) 
     &            + A2t(:) * g2yti(:)
     &            + A3t(:) * g3yti(:) 

!
!.... add contribution to rmi
!
!        if ((ires .eq. 2) .or. (ires .eq. 3))
!     &    rmi(:,16:20) = rLyi  ! modified residual uses non i.b.p form of conv
!
!.... ---------------------->  LHS   <-----------------------
!
        if (lhs .eq. 1) then
!
!.... loop through the columns
!
        do j = 1, nshl

!
!.... first compute (A_i N_b,i)
!
           AitNbi(:) = 
     &                    WdetJ * shg(:,j,1) * A1t(:) 
     &                  + WdetJ * shg(:,j,2) * A2t(:) 
     &                  + WdetJ * shg(:,j,3) * A3t(:)

!
!.... now loop through the rows and add (N_a A_i N_b,i) to
!     the tangent matrix.
!
           do i = 1, nshl
                 
             EGmasst(:,i,j) = EGmasst(:,i,j) +  shape(:,i) * AitNbi(:)

            
!
!.... end loop on rows
!              
           enddo
!
!.... end loop on columns
!
        enddo
!
!.... end of LHS tangent matrix computation
!        
        endif
      
        ttim(22) = ttim(22) + tmr()
!
!.... return
!
        return
        end

