        subroutine e3conv (g1yi,   g2yi,   g3yi,
     &                     A1,     A2,     A3, 
     &                     rho,    pres,   T,
     &                     ei,     rk,     u1,
     &                     u2,     u3,     rLyi,   
     &                     ri,     rmi,    EGmass,
     &                     shg,    shape,  WdetJ )
c
c----------------------------------------------------------------------
c
c This routine calculates the contribution of Galerkin part of the 
c Convective term (Time and Euler fluxes) to both RHS and LHS.
c
c input:
c  g1yi   (npro,nflow)    : grad-y in direction 1
c  g2yi   (npro,nflow)    : grad-y in direction 2
c  g3yi   (npro,nflow)    : grad-y in direction 3
c  A1    (npro,nflow,nflow)  : A-1
c  A2    (npro,nflow,nflow)  : A-2
c  A3    (npro,nflow,nflow)  : A-3 
c  rho    (npro)         : density
c  pres   (npro)         : pressure
c  T      (npro)         : temperature
c  ei     (npro)         : internal energy
c  rk     (npro)         : kinetic energy
c  u1     (npro)         : x1-velocity component
c  u2     (npro)         : x2-velocity component
c  u3     (npro)         : x3-velocity component
c  shg    (npro,nshl,nsd) : global grad's of shape functions
c  shape  (npro,nshl)  : element shape functions      
c  WdetJ  (npro)         : weighted Jacobian determinant
c     
c output:
c  rLyi   (npro,nflow)           : least-squares residual vector
c  ri     (npro,nflow*(nsd+1))   : partial residual
c  rmi    (npro,nflow*(nsd+1))   : partial modified residual
c  EGmass (npro,nedof,nedof)    : partial LHS tangent matrix
c
c
c Zdenek Johan, Summer 1990. (Modified from e2conv.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997 Primitive Variables
c----------------------------------------------------------------------
c
        include "common.h"
c
c  passed arrays
c
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
c
c  local arrays  
c
        dimension AiNbi(npro,nflow,nflow),     fact1(npro),
     &            fact2(npro),               fact3(npro)
        
	ttim(22) = ttim(22) - secs(0.0)

c
c.... ---------------------->  RHS, Euler Flux  <----------------------
c
        if ((ires .eq. 1) .or. (ires .eq. 3)) then
c
c.... calculate integrated by part contribution of Euler flux (Galerkin)
c
          ri(:, 1) = (- u1) * rho
          ri(:, 2) = (- u1) * rho * u1 - pres
          ri(:, 3) = (- u1) * rho * u2
          ri(:, 4) = (- u1) * rho * u3
          ri(:, 5) = (- u1) * rho * (ei + rk) - u1 * pres
c
          ri(:, 6) = (- u2) * rho
          ri(:, 7) = (- u2) * rho * u1
          ri(:, 8) = (- u2) * rho * u2 - pres
          ri(:, 9) = (- u2) * rho * u3
          ri(:,10) = (- u2) * rho * (ei + rk) - u2 * pres
c
          ri(:,11) = (- u3) * rho
          ri(:,12) = (- u3) * rho * u1
          ri(:,13) = (- u3) * rho * u2
          ri(:,14) = (- u3) * rho * u3 - pres
          ri(:,15) = (- u3) * rho * (ei + rk) - u3 * pres
c
          flops = flops + 28*npro
c
        endif
c
c.... calculate ( A_i Y,i ) --> rLyi   Commented out zeros of A matrices
c
        rLyi(:,1) = 
     &              A1(:,1,1) * g1yi(:,1) 
     &            + A1(:,1,2) * g1yi(:,2)
c    &            + A1(:,1,3) * g1yi(:,3) 
c    &            + A1(:,1,4) * g1yi(:,4)
     &            + A1(:,1,5) * g1yi(:,5)
     &            + A2(:,1,1) * g2yi(:,1) 
c    &            + A2(:,1,2) * g2yi(:,2)
     &            + A2(:,1,3) * g2yi(:,3) 
c    &            + A2(:,1,4) * g2yi(:,4)
     &            + A2(:,1,5) * g2yi(:,5)
     &            + A3(:,1,1) * g3yi(:,1) 
c    &            + A3(:,1,2) * g3yi(:,2)
c    &            + A3(:,1,3) * g3yi(:,3) 
     &            + A3(:,1,4) * g3yi(:,4)
     &            + A3(:,1,5) * g3yi(:,5)
        rLyi(:,2) = 
     &              A1(:,2,1) * g1yi(:,1) 
     &            + A1(:,2,2) * g1yi(:,2)
c    &            + A1(:,2,3) * g1yi(:,3) 
c    &            + A1(:,2,4) * g1yi(:,4)
     &            + A1(:,2,5) * g1yi(:,5)
     &            + A2(:,2,1) * g2yi(:,1) 
     &            + A2(:,2,2) * g2yi(:,2)
     &            + A2(:,2,3) * g2yi(:,3) 
c    &            + A2(:,2,4) * g2yi(:,4)
     &            + A2(:,2,5) * g2yi(:,5)
     &            + A3(:,2,1) * g3yi(:,1) 
     &            + A3(:,2,2) * g3yi(:,2)
c    &            + A3(:,2,3) * g3yi(:,3) 
     &            + A3(:,2,4) * g3yi(:,4)
     &            + A3(:,2,5) * g3yi(:,5)
        rLyi(:,3) = 
     &              A1(:,3,1) * g1yi(:,1) 
     &            + A1(:,3,2) * g1yi(:,2)
     &            + A1(:,3,3) * g1yi(:,3) 
c    &            + A1(:,3,4) * g1yi(:,4)
     &            + A1(:,3,5) * g1yi(:,5)
     &            + A2(:,3,1) * g2yi(:,1) 
c    &            + A2(:,3,2) * g2yi(:,2)
     &            + A2(:,3,3) * g2yi(:,3) 
c    &            + A2(:,3,4) * g2yi(:,4)
     &            + A2(:,3,5) * g2yi(:,5)
     &            + A3(:,3,1) * g3yi(:,1) 
c    &            + A3(:,3,2) * g3yi(:,2)
     &            + A3(:,3,3) * g3yi(:,3) 
     &            + A3(:,3,4) * g3yi(:,4)
     &            + A3(:,3,5) * g3yi(:,5)
        rLyi(:,4) = 
     &              A1(:,4,1) * g1yi(:,1) 
     &            + A1(:,4,2) * g1yi(:,2)
c    &            + A1(:,4,3) * g1yi(:,3) 
     &            + A1(:,4,4) * g1yi(:,4)
     &            + A1(:,4,5) * g1yi(:,5)
     &            + A2(:,4,1) * g2yi(:,1) 
c    &            + A2(:,4,2) * g2yi(:,2)
     &            + A2(:,4,3) * g2yi(:,3) 
     &            + A2(:,4,4) * g2yi(:,4)
     &            + A2(:,4,5) * g2yi(:,5)
     &            + A3(:,4,1) * g3yi(:,1) 
c    &            + A3(:,4,2) * g3yi(:,2)
c    &            + A3(:,4,3) * g3yi(:,3) 
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
c
c.... add contribution to rmi
c
        if ((ires .eq. 2) .or. (ires .eq. 3))
     &    rmi(:,16:20) = rLyi  ! modified residual uses non i.b.p form of conv.
c
c.... ---------------------->  LHS   <-----------------------
c
        if (lhs .eq. 1) then
c
c.... loop through the columns
c
        do j = 1, nshl
           j0 = nflow * (j - 1)
c
c.... compute some useful factors
c
           fact1 = WdetJ * shg(:,j,1)
           fact2 = WdetJ * shg(:,j,2)
           fact3 = WdetJ * shg(:,j,3)
c
c.... first compute (A_i N_b,i)
c
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
c
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
c
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
c
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
c
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
c
c.... now loop through the row nodes and add (N_a A_i N_b,i) to
c     the tangent matrix.
c
           do i = 1, nshl
              i0 = nflow * (i - 1)
c
c.... loop through dof's
c
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
c
c.... end loop on rows
c              
           enddo
c
c.... end loop on columns
c
        enddo
c
c.... end of LHS tangent matrix computation
c        
        endif
      
        ttim(22) = ttim(22) + secs(0.0)
c
c.... return
c
        return
        end
c
c
c
        subroutine e3convSclr (g1yti,   g2yti,   g3yti,
     &                         A1t,     A2t,     A3t, 
     &                         rho,     u1,      Sclr,
     &                         u2,      u3,      rLyti,   
     &                         rti,     rmti,    EGmasst,
     &                         shg,     shape,   WdetJ)
c
c----------------------------------------------------------------------
c
c This routine calculates the contribution of Galerkin part of the 
c Convective term (Time and Euler fluxes) to both RHS and LHS.
c
c input:
c  Sclr   (npro)          : Scalar variable
c  g1yti  (npro)          : grad-y in direction 1
c  g2yti  (npro)          : grad-y in direction 2
c  g3yti  (npro)          : grad-y in direction 3
c  A1t    (npro)          : A-1
c  A2t    (npro)          : A-2
c  A3t    (npro)          : A-3 
c  rho    (npro)          : density
c  u1     (npro)          : x1-velocity component
c  u2     (npro)          : x2-velocity component
c  u3     (npro)          : x3-velocity component
c  shg    (npro,nshl,nsd) : global grad's of shape functions
c  shape  (npro,nshl)     : element shape functions      
c  WdetJ  (npro)          : weighted Jacobian determinant
c     
c output:
c  rLyti   (npro)         : least-squares residual vector
c  rti     (npro,nsd+1)   : partial residual
c  rmti    (npro,nsd+1)   : partial modified residual
c  EGmasst (npro,nshape,nshape): partial LHS tangent matrix
c
c
c Zdenek Johan, Summer 1990. (Modified from e2conv.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997 Primitive Variables
c----------------------------------------------------------------------
c
        include "common.h"
c
c  passed arrays
c
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
c
c  local arrays  
c
        dimension AitNbi(npro)
       
	ttim(22) = ttim(22) - tmr(0.0)
c
c.... ---------------------->  RHS, Euler Flux  <----------------------
c
        if ((ires .eq. 1) .or. (ires .eq. 3)) then
c
c.... calculate integrated by part contribution of Euler flux (Galerkin)
c     
           if (iconvsclr.eq.2) then ! convective form
c     
              rti(:, 4) = rti(:,4) + ( u1) * g1yti(:)
     &                             + ( u2) * g2yti(:)
     &                             + ( u3) * g3yti(:) 
c     
           else                 ! conservative form     
c     
              rti(:, 1) = rti(:,1) + (- u1) * rho * Sclr
              rti(:, 2) = rti(:,2) + (- u2) * rho * Sclr
              rti(:, 3) = rti(:,3) + (- u3) * rho * Sclr
c     
           endif

           flops = flops + 28*npro

        endif
c
c.... calculate ( A_i Y,i ) --> rLyi
c
        rLyti(:) = rLyti(:)
     &            + A1t(:) * g1yti(:) 
     &            + A2t(:) * g2yti(:)
     &            + A3t(:) * g3yti(:) 

c
c.... add contribution to rmi
c
c        if ((ires .eq. 2) .or. (ires .eq. 3))
c     &    rmi(:,16:20) = rLyi  ! modified residual uses non i.b.p form of conv
c
c.... ---------------------->  LHS   <-----------------------
c
        if (lhs .eq. 1) then
c
c.... loop through the columns
c
        do j = 1, nshl

c
c.... first compute (A_i N_b,i)
c
           AitNbi(:) = 
     &                    WdetJ * shg(:,j,1) * A1t(:) 
     &                  + WdetJ * shg(:,j,2) * A2t(:) 
     &                  + WdetJ * shg(:,j,3) * A3t(:)

c
c.... now loop through the rows and add (N_a A_i N_b,i) to
c     the tangent matrix.
c
           do i = 1, nshl
                 
             EGmasst(:,i,j) = EGmasst(:,i,j) +  shape(:,i) * AitNbi(:)

            
c
c.... end loop on rows
c              
           enddo
c
c.... end loop on columns
c
        enddo
c
c.... end of LHS tangent matrix computation
c        
        endif
      
        ttim(22) = ttim(22) + tmr()
c
c.... return
c
        return
        end

