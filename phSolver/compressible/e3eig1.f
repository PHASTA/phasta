        subroutine e3eig1 (rho,    T,      cp,     gamb,   c,
     &                     u1,     u2,     u3,     a1,     a2,
     &                     a3,     eb1,
     &                     dxidx,  u,      Q)
c
c----------------------------------------------------------------------
c
c This routine performs the first step of the eigenvalue decomposition
c of the Tau matrix.
c
c input:
c  rho    (npro)           : density
c  T      (npro)           : temperature
c  cp     (npro)           : specific heat at constant pressure
c  gamb   (npro)           : gamma_bar (defined in paper by Chalot et al.)
c  c      (npro)           : speed of sound
c  u1     (npro)           : x1-velocity component
c  u2     (npro)           : x2-velocity component
c  u3     (npro)           : x3-velocity component
c  a1     (npro)           : x1-acceleration component
c  a2     (npro)           : x2-acceleration component
c  a3     (npro)           : x3-acceleration component
c  eb1    (npro)           : e1_bar (defined in paper by Chalot et al.)
c  dxidx  (npro,nsd,nsd)   : inverse of deformation gradient
c
c output:
c  u      (npro)           : fluid velocity
c  a1     (npro)           : aspect ratio factor in streamline direction
c  a2     (npro)           : aspect ratio factor in normal_1 direction
c  a3     (npro)           : aspect ratio factor in normal_2 direction
c  Q      (npro,nflow,nflow) : 1st level eigenvectors of Tau matrix
c  
c
c Zdenek Johan, Summer 1990.  (Modified from e2tau.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension rho(npro),                 T(npro),
     &            cp(npro),                  gamb(npro),
     &            c(npro),                   u1(npro),
     &            u2(npro),                  u3(npro),
     &            a1(npro),                  a2(npro),
     &            a3(npro),                  
     &            eb1(npro),                 dxidx(npro,nsd,nsd),
     &            u(npro),                   Q(npro,nflow,nflow)
c
        dimension Rcs(npro,nsd,nsd),         fact(npro),
     &            temp(npro)
c
c.... compute the directional cosines (streamline direction)
c
c

        where (u .ne. zero)
           fact       = one / u
           Rcs(:,1,1) = u1 * fact
           Rcs(:,1,2) = u2 * fact
           Rcs(:,1,3) = u3 * fact
        elsewhere
           Rcs(:,1,1) = one
           Rcs(:,1,2) = zero
           Rcs(:,1,3) = zero
        endwhere
c
c.... compute the directional cosines (normal acceleration direction)
c

        
        fact = a1 * Rcs(:,1,1) + a2 * Rcs(:,1,2) + a3 * Rcs(:,1,3)
        a1   = a1 - fact * Rcs(:,1,1)
        a2   = a2 - fact * Rcs(:,1,2)
        a3   = a3 - fact * Rcs(:,1,3)
        fact = a1**2 + a2**2 + a3**2
c
        where (fact .gt. epsM)
c
          Rcs(:,2,1) = a1
          Rcs(:,2,2) = a2
          Rcs(:,2,3) = a3
c
        elsewhere
c
          Rcs(:,2,1) =   Rcs(:,1,2)
     &                 + sign(one, Rcs(:,1,2)*Rcs(:,1,3)) * Rcs(:,1,3)
          Rcs(:,2,2) = - Rcs(:,1,1)
     &                 - sign(one, Rcs(:,1,2)*Rcs(:,1,3)) * Rcs(:,1,3)
          Rcs(:,2,3) =   sign(one, Rcs(:,1,2)*Rcs(:,1,3)) *
     &                 ( Rcs(:,1,2) - Rcs(:,1,1) )
c
          fact = Rcs(:,2,1)**2 + Rcs(:,2,2)**2 + Rcs(:,2,3)**2
c
        endwhere
c
	fact = one / sqrt(fact)
c
        Rcs(:,2,1) = Rcs(:,2,1) * fact
        Rcs(:,2,2) = Rcs(:,2,2) * fact
        Rcs(:,2,3) = Rcs(:,2,3) * fact
c
c.... compute the directional cosines (last direction)
c
        Rcs(:,3,1) = Rcs(:,1,2) * Rcs(:,2,3) - Rcs(:,1,3) * Rcs(:,2,2)
        Rcs(:,3,2) = Rcs(:,1,3) * Rcs(:,2,1) - Rcs(:,1,1) * Rcs(:,2,3)
        Rcs(:,3,3) = Rcs(:,1,1) * Rcs(:,2,2) - Rcs(:,1,2) * Rcs(:,2,1)

c
c.... calculate the element aspect ratio factors
c
        a1 = Rcs(:,1,1) * dxidx(:,1,1) + Rcs(:,1,2) * dxidx(:,1,2) +
     &       Rcs(:,1,3) * dxidx(:,1,3)
        a2 = Rcs(:,2,1) * dxidx(:,1,1) + Rcs(:,2,2) * dxidx(:,1,2) +
     &       Rcs(:,2,3) * dxidx(:,1,3)
        a3 = Rcs(:,3,1) * dxidx(:,1,1) + Rcs(:,3,2) * dxidx(:,1,2) +
     &       Rcs(:,3,3) * dxidx(:,1,3)
        dxidx(:,1,1) = a1
        dxidx(:,1,2) = a2
        dxidx(:,1,3) = a3
c
        a1 = Rcs(:,1,1) * dxidx(:,2,1) + Rcs(:,1,2) * dxidx(:,2,2) +
     &       Rcs(:,1,3) * dxidx(:,2,3)
        a2 = Rcs(:,2,1) * dxidx(:,2,1) + Rcs(:,2,2) * dxidx(:,2,2) +
     &       Rcs(:,2,3) * dxidx(:,2,3)
        a3 = Rcs(:,3,1) * dxidx(:,2,1) + Rcs(:,3,2) * dxidx(:,2,2) +
     &       Rcs(:,3,3) * dxidx(:,2,3)
        dxidx(:,2,1) = a1
        dxidx(:,2,2) = a2
        dxidx(:,2,3) = a3
c
        a1 = Rcs(:,1,1) * dxidx(:,3,1) + Rcs(:,1,2) * dxidx(:,3,2) +
     &       Rcs(:,1,3) * dxidx(:,3,3)
        a2 = Rcs(:,2,1) * dxidx(:,3,1) + Rcs(:,2,2) * dxidx(:,3,2) +
     &       Rcs(:,2,3) * dxidx(:,3,3)
        a3 = Rcs(:,3,1) * dxidx(:,3,1) + Rcs(:,3,2) * dxidx(:,3,2) +
     &       Rcs(:,3,3) * dxidx(:,3,3)
        dxidx(:,3,1) = a1
        dxidx(:,3,2) = a2
        dxidx(:,3,3) = a3

c
c... original
c

        a1 = dxidx(:,1,1)**2 + dxidx(:,2,1)**2 + dxidx(:,3,1)**2
        a2 = dxidx(:,1,2)**2 + dxidx(:,2,2)**2 + dxidx(:,3,2)**2
        a3 = dxidx(:,1,3)**2 + dxidx(:,2,3)**2 + dxidx(:,3,3)**2

c
c... change from original (analyt., could be error)
c

cc        a1 = dxidx(:,1,1)**2 + dxidx(:,1,2)**2 + dxidx(:,1,3)**2
cc        a2 = dxidx(:,2,1)**2 + dxidx(:,2,2)**2 + dxidx(:,2,3)**2
cc        a3 = dxidx(:,3,1)**2 + dxidx(:,3,2)**2 + dxidx(:,3,3)**2        

c
c.... correct for tetrahedra
c
        if (lcsyst .eq. 1) then
          a1 = ( a1 + (dxidx(:,1,1) + dxidx(:,2,1) +
     &                 dxidx(:,3,1))**2 ) * pt39
          a2 = ( a2 + (dxidx(:,1,2) + dxidx(:,2,2) +
     &                 dxidx(:,3,2))**2 ) * pt39
          a3 = ( a3 + (dxidx(:,1,3) + dxidx(:,2,3) +
     &                 dxidx(:,3,3))**2 ) * pt39
c
!      flops = flops + 15*npro
        endif
c
c.... set up the 1st level Eigenvectors =  R_t*Tau^*
c
        fact     =  (one / sqrt( two * rho * T )) / c
c
        Q(:,1,5) =  fact * ((c + u) * c - eb1 * gamb)
        Q(:,2,5) = -fact * (c + u * gamb) * Rcs(:,1,1)
        Q(:,3,5) = -fact * (c + u * gamb) * Rcs(:,1,2)
        Q(:,4,5) = -fact * (c + u * gamb) * Rcs(:,1,3)
        Q(:,5,5) =  fact * gamb
c
        Q(:,1,1) =  fact * ((c - u) * c - eb1 * gamb)
        Q(:,2,1) =  fact * (c - u * gamb) * Rcs(:,1,1)
        Q(:,3,1) =  fact * (c - u * gamb) * Rcs(:,1,2)
        Q(:,4,1) =  fact * (c - u * gamb) * Rcs(:,1,3)
        Q(:,5,1) =  fact * gamb
c
        Q(:,1,3) =  zero
        Q(:,2,3) = -fact * c * sqt2 * Rcs(:,2,1) ! + in original Jo/Sh code
        Q(:,3,3) = -fact * c * sqt2 * Rcs(:,2,2) ! could be error here,
        Q(:,4,3) = -fact * c * sqt2 * Rcs(:,2,3) ! but unlikely
        Q(:,5,3) =  zero
c
        Q(:,1,2) =  zero
        Q(:,2,2) =  fact * c * sqt2 * Rcs(:,3,1)
        Q(:,3,2) =  fact * c * sqt2 * Rcs(:,3,2)
        Q(:,4,2) =  fact * c * sqt2 * Rcs(:,3,3)
        Q(:,5,2) =  zero
c
        fact     =  (one / sqrt( cp * rho )) / T
c
        Q(:,1,4) =  fact * eb1
        Q(:,2,4) =  fact * u * Rcs(:,1,1)
        Q(:,3,4) =  fact * u * Rcs(:,1,2)
        Q(:,4,4) =  fact * u * Rcs(:,1,3)
        Q(:,5,4) = -fact
c
c.... flop count
c

c        do i = 1, nflow
c           do j = 1, nflow
c              Q(:,i,j) = abs(Q(:,i,j)) !make sure eigenv. are positive
c           enddo
c        enddo

   !      flops = flops + 203*npro
c
c.... return
c
        return
        end


        subroutine e3eig2 (u,      c,      AR1,    AR2,    AR3,
     &                     rlam,   Q,    eigmax)
c
c----------------------------------------------------------------------
c
c  This routine diagonalizes a partially reduced matrix using the 
c Jacobi transformation.  This routine assumes that the original 
c 5x5 matrix is already reduced to 2x2.
c
c input:
c  u      (npro)           : fluid velocity
c  c      (npro)           : speed of sound
c  AR1    (npro)           : aspect ratio factor in streamline direction
c  AR2    (npro)           : aspect ratio factor in normal_1 direction
c  AR3    (npro)           : aspect ratio factor in normal_2 direction
c  Q      (npro,nflow,nflow) : 1st level eigenvectors of Tau matrix
c
c output:
c  rlam   (npro,nflow)      : eigenvalues
c  Q      (npro,nflow,nflow) : eigenvectors
c
c
c                                  T
c  This routine finds S that      S  Atau S = Lam  (Lam is Symm.)
c
c  Then returns       rlam <- Lam     and      Q <- Q S
c
c
c Note: Jacobi transformation is extracted (and modified) from
c       "Numerical Recipes: The Art of Scientific Computing" by
c       Press, Flannery, Teukolsky and Vetterling, pp. 342-349.
c
c
c Farzin Shakib, Spring 1989.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension u(npro),                   c(npro),
     &            AR1(npro),                 AR2(npro),
     &            AR3(npro),                 rlam(npro,nflow),
     &            Q(npro,nflow,nflow),         eigmax(npro,nflow)
c
        dimension offd(npro),                t(npro),
     &            Rcs(npro),                 Rsn(npro)
c
c.... set the reduced eigensystem
c
        offd      = (AR2 + AR3) * pt5 * c**2
c
        rlam(:,1) = AR1 * (u + c)**2 + offd
        rlam(:,2) = AR1 * u**2       + AR3 * c**2
        rlam(:,3) = AR1 * u**2       + AR2 * c**2
        rlam(:,4) = AR1 * u**2
        rlam(:,5) = AR1 * (u - c)**2 + offd
c
c.... modify for time dependent problems
c
! consider time term if iremoveStabTimeTerm is set to zero
        if(iremoveStabTimeTerm.eq.0) then
           tmp  = dtsfct * four * (Dtgl * Dtgl)
           rlam(:,:) = rlam(:,:) + tmp
        endif
c
c.... compute the rotation tangent ( IEEE arithmetic if offd=0 )
c
        t = pt5 * (rlam(:,1) - rlam(:,5)) / offd
        t = sign(one, t) / ( abs(t) + sqrt(one + t**2) )
c
c.... compute cosine and sin, and rotate the eigenvalues
c
        Rcs = one / sqrt( one + t**2 )
        Rsn = t * Rcs
c
        rlam(:,1) = rlam(:,1) - t * offd
        rlam(:,5) = rlam(:,5) + t * offd
c
c.... transform the Eigenvectors (all 5 components)
c
        t        = Rcs * Q(:,1,1) - Rsn * Q(:,1,5)
        Q(:,1,5) = Rsn * Q(:,1,1) + Rcs * Q(:,1,5)
        Q(:,1,1) = t
c
        t        = Rcs * Q(:,2,1) - Rsn * Q(:,2,5)
        Q(:,2,5) = Rsn * Q(:,2,1) + Rcs * Q(:,2,5)
        Q(:,2,1) = t
c
        t        = Rcs * Q(:,3,1) - Rsn * Q(:,3,5)
        Q(:,3,5) = Rsn * Q(:,3,1) + Rcs * Q(:,3,5)
        Q(:,3,1) = t
c
        t        = Rcs * Q(:,4,1) - Rsn * Q(:,4,5)
        Q(:,4,5) = Rsn * Q(:,4,1) + Rcs * Q(:,4,5)
        Q(:,4,1) = t
c
        t        = Rcs * Q(:,5,1) - Rsn * Q(:,5,5)
        Q(:,5,5) = Rsn * Q(:,5,1) + Rcs * Q(:,5,5)
        Q(:,5,1) = t

c
c.... extract maximum eigenvalues
c
        eigmax(:,1) = max(rlam(:,1), rlam(:,2), rlam(:,3),
     &                    rlam(:,4), rlam(:,5))
        eigmax(:,2) = eigmax(:,1)
        eigmax(:,3) = eigmax(:,1)
        eigmax(:,4) = eigmax(:,1)
        eigmax(:,5) = eigmax(:,1)

c        
c.... flop count
c
   !      flops = flops + 85*npro
c
c.... return
c
        return
        end
