        subroutine e3bvar (yl,      ycl,     BCB,     shpb,    shglb,   
     &                     xlb,     lnode,   g1yi,    g2yi,
     &                     g3yi,    WdetJb,  bnorm,   pres,    T,
     &                     u1,      u2,      u3,      rho,     ei,
     &                     cp,      rk,      
     &                     rou,     p,       tau1n,   tau2n,   tau3n,
     &                     heat,    dNadx)
c
c----------------------------------------------------------------------
c
c   This routine computes the variables at integration points for 
c the boundary element routine.
c
c input:
c  yl     (npro,nshl,nflow)   : primitive variables (perturbed, no scalars)
c  ycl    (npro,nshl,ndof)    : primitive variables
c  BCB    (npro,nshlb,ndBCB)  : Boundary Condition values
c  shpb   (npro,nshl)         : boundary element shape-functions
c  shglb  (npro,nsd,nshl)     : boundary element grad-shape-functions
c  xlb    (npro,nenl,nsd)       : nodal coordinates at current step
c  lnode  (nenb)                : local nodes on the boundary
c
c output:
c  g1yi   (npro,nflow)           : grad-v in direction 1
c  g2yi   (npro,nflow)           : grad-v in direction 2
c  g3yi   (npro,nflow)           : grad-v in direction 3
c  WdetJb (npro)                : weighted Jacobian
c  bnorm  (npro,nsd)            : outward normal
c  pres   (npro)                : pressure
c  T      (npro)                : temperature
c  u1     (npro)                : x1-velocity component
c  u2     (npro)                : x2-velocity component
c  u3     (npro)                : x3-velocity component
c  rho    (npro)                : density
c  ei     (npro)                : internal energy
c  cp     (npro)                : specific energy at constant pressure
c  rk     (npro)                : kinetic energy
c  rou    (npro)                : BC mass flux
c  p      (npro)                : BC pressure
c  tau1n  (npro)                : BC viscous flux 1
c  tau2n  (npro)                : BC viscous flux 2
c  tau3n  (npro)                : BC viscous flux 3
c  heat   (npro)                : BC heat flux
c  dNdx   (npro, nsd)           : BC element shape function gradients
c
c
c Zdenek Johan, Summer 1990.  (Modified from e2bvar.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension yl(npro,nshl,nflow),      BCB(npro,nshlb,ndBCB),
     &            ycl(npro,nshl,ndof),
     &            shpb(npro,nshl),
     &            shglb(npro,nsd,nshl), 
     &            xlb(npro,nenl,nsd),        
     &            lnode(27),               g1yi(npro,nflow),
     &            g2yi(npro,nflow),        g3yi(npro,nflow),
     &            WdetJb(npro),            bnorm(npro,nsd),
     &            pres(npro),              T(npro),
     &            u1(npro),                u2(npro),
     &            u3(npro),                rho(npro),
     &            ei(npro),                cp(npro),
     &            rk(npro),                  
     &            rou(npro),               p(npro),
     &            tau1n(npro),             tau2n(npro),
     &            tau3n(npro),             heat(npro)

        dimension gl1yi(npro,nflow),       gl2yi(npro,nflow),
     &            gl3yi(npro,nflow),       dxdxib(npro,nsd,nsd),
     &            dxidxb(npro,nsd,nsd),    temp(npro),
     &            temp1(npro),             temp2(npro),
     &            temp3(npro),
     &            dNadx(npro, nshl, nsd),  dNadxi(npro, nshl, nsd)

        dimension h(npro),                 cv(npro),
     &            alfap(npro),             betaT(npro),
     &            gamb(npro),              c(npro),
     &            tmp(npro),
     &            v1(npro,nsd),            v2(npro,nsd)

        integer   aa
c
c.... ------------------->  integration variables  <--------------------
c
c.... compute the primitive variables at the integration point
c
        pres = zero
        u1   = zero
        u2   = zero
        u3   = zero
        T    = zero
c
        do n = 1, nshlb
          nodlcl = lnode(n)
c
          pres = pres + shpb(:,nodlcl) * yl(:,nodlcl,1)
          u1   = u1   + shpb(:,nodlcl) * yl(:,nodlcl,2)
          u2   = u2   + shpb(:,nodlcl) * yl(:,nodlcl,3)
          u3   = u3   + shpb(:,nodlcl) * yl(:,nodlcl,4)
          T    = T    + shpb(:,nodlcl) * yl(:,nodlcl,5)
        enddo
c
c.... calculate the specific kinetic energy
c
        rk = pt5 * ( u1**2 + u2**2  + u3**2 )
c
c.... get the thermodynamic properties
c
        if (iLSet .ne. 0)then
           temp = zero
           isc=abs(iRANS)+6
           do n = 1, nshlb
              temp = temp + shpb(:,n) * ycl(:,n,isc)
           enddo
        endif

        ithm = 6
        if (Navier .eq. 1) ithm = 7
        call getthm (pres,            T,                  temp,
     &               rk,              rho,                ei,
     &               h,               tmp,                cv,
     &               cp,              alfap,              betaT,
     &               gamb,            c)
c
c.... ---------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
        dxdxib = zero
c
        do n = 1, nenl
           dxdxib(:,1,1) = dxdxib(:,1,1) + xlb(:,n,1) * shglb(:,1,n)
           dxdxib(:,1,2) = dxdxib(:,1,2) + xlb(:,n,1) * shglb(:,2,n)
           dxdxib(:,1,3) = dxdxib(:,1,3) + xlb(:,n,1) * shglb(:,3,n)
           dxdxib(:,2,1) = dxdxib(:,2,1) + xlb(:,n,2) * shglb(:,1,n)
           dxdxib(:,2,2) = dxdxib(:,2,2) + xlb(:,n,2) * shglb(:,2,n)
           dxdxib(:,2,3) = dxdxib(:,2,3) + xlb(:,n,2) * shglb(:,3,n)
           dxdxib(:,3,1) = dxdxib(:,3,1) + xlb(:,n,3) * shglb(:,1,n)
           dxdxib(:,3,2) = dxdxib(:,3,2) + xlb(:,n,3) * shglb(:,2,n)
           dxdxib(:,3,3) = dxdxib(:,3,3) + xlb(:,n,3) * shglb(:,3,n)
        enddo
c
c.... compute the normal to the boundary
c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
        v1 = xlb(:,2,:) - xlb(:,1,:)
        v2 = xlb(:,3,:) - xlb(:,1,:)
c
c.....The following are done in order to correct temp1..3  
c     based on the results from compressible code.  This is done only 
c     for wedges, depending on the bounary face.(tri or quad)  
c        
        if (lcsyst .eq. 4) then
           temp1 = dxdxib(:,2,1) * dxdxib(:,3,3) -
     &             dxdxib(:,2,3) * dxdxib(:,3,1)
           temp2 = dxdxib(:,3,1) * dxdxib(:,1,3) -
     &             dxdxib(:,3,3) * dxdxib(:,1,1)
           temp3 = dxdxib(:,1,1) * dxdxib(:,2,3) -
     &             dxdxib(:,1,3) * dxdxib(:,2,1)
           
        elseif (lcsyst .eq. 1) then
           temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
           temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
           temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
        else 
           temp1 = - v1(:,2) * v2(:,3) + v2(:,2) * v1(:,3)
           temp2 = - v2(:,1) * v1(:,3) + v1(:,1) * v2(:,3)
           temp3 = - v1(:,1) * v2(:,2) + v2(:,1) * v1(:,2)
        endif
c
        temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
        bnorm(:,1) = temp1 * temp
        bnorm(:,2) = temp2 * temp
        bnorm(:,3) = temp3 * temp
c
      
        if (lcsyst .eq. 3) then
           WdetJb     = (1 - Qwtb(lcsyst,intp)) / (four*temp)
        elseif (lcsyst .eq. 4) then
c
c  quad face wedges have a conflict in lnode ordering that makes the
c  normal negative
c
c          bnorm=-bnorm
c
          WdetJb     = Qwtb(lcsyst,intp) / temp
        else
           WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
        endif
c
c.... -------------------------->  Grad-V  <----------------------------
c
c.... compute grad-v for Navier-Stokes terms
c
        if (Navier .eq. 1) then
c
c.... compute the inverse of deformation gradient
c
          dxidxb(:,1,1) =   dxdxib(:,2,2) * dxdxib(:,3,3) 
     &                    - dxdxib(:,3,2) * dxdxib(:,2,3)
          dxidxb(:,1,2) =   dxdxib(:,3,2) * dxdxib(:,1,3) 
     &                    - dxdxib(:,1,2) * dxdxib(:,3,3)
          dxidxb(:,1,3) =   dxdxib(:,1,2) * dxdxib(:,2,3) 
     &                    - dxdxib(:,1,3) * dxdxib(:,2,2)
          temp          = one / ( dxidxb(:,1,1) * dxdxib(:,1,1) 
     &                          + dxidxb(:,1,2) * dxdxib(:,2,1)  
     &                          + dxidxb(:,1,3) * dxdxib(:,3,1) )
          dxidxb(:,1,1) =  dxidxb(:,1,1) * temp
          dxidxb(:,1,2) =  dxidxb(:,1,2) * temp
          dxidxb(:,1,3) =  dxidxb(:,1,3) * temp
          dxidxb(:,2,1) = (dxdxib(:,2,3) * dxdxib(:,3,1) 
     &                   - dxdxib(:,2,1) * dxdxib(:,3,3)) * temp
          dxidxb(:,2,2) = (dxdxib(:,1,1) * dxdxib(:,3,3) 
     &                   - dxdxib(:,3,1) * dxdxib(:,1,3)) * temp
          dxidxb(:,2,3) = (dxdxib(:,2,1) * dxdxib(:,1,3) 
     &                   - dxdxib(:,1,1) * dxdxib(:,2,3)) * temp
          dxidxb(:,3,1) = (dxdxib(:,2,1) * dxdxib(:,3,2) 
     &                   - dxdxib(:,2,2) * dxdxib(:,3,1)) * temp
          dxidxb(:,3,2) = (dxdxib(:,3,1) * dxdxib(:,1,2) 
     &                   - dxdxib(:,1,1) * dxdxib(:,3,2)) * temp
          dxidxb(:,3,3) = (dxdxib(:,1,1) * dxdxib(:,2,2) 
     &                   - dxdxib(:,1,2) * dxdxib(:,2,1)) * temp
c
c.... compute local-grad-Y
c
          gl1yi = zero
          gl2yi = zero
          gl3yi = zero
c
          do n = 1, nshl
            gl1yi(:,1) = gl1yi(:,1) + shglb(:,1,n) * yl(:,n,1)
            gl1yi(:,2) = gl1yi(:,2) + shglb(:,1,n) * yl(:,n,2)
            gl1yi(:,3) = gl1yi(:,3) + shglb(:,1,n) * yl(:,n,3)
            gl1yi(:,4) = gl1yi(:,4) + shglb(:,1,n) * yl(:,n,4)
            gl1yi(:,5) = gl1yi(:,5) + shglb(:,1,n) * yl(:,n,5)
c
            gl2yi(:,1) = gl2yi(:,1) + shglb(:,2,n) * yl(:,n,1)
            gl2yi(:,2) = gl2yi(:,2) + shglb(:,2,n) * yl(:,n,2)
            gl2yi(:,3) = gl2yi(:,3) + shglb(:,2,n) * yl(:,n,3)
            gl2yi(:,4) = gl2yi(:,4) + shglb(:,2,n) * yl(:,n,4)
            gl2yi(:,5) = gl2yi(:,5) + shglb(:,2,n) * yl(:,n,5)
c
            gl3yi(:,1) = gl3yi(:,1) + shglb(:,3,n) * yl(:,n,1)
            gl3yi(:,2) = gl3yi(:,2) + shglb(:,3,n) * yl(:,n,2)
            gl3yi(:,3) = gl3yi(:,3) + shglb(:,3,n) * yl(:,n,3)
            gl3yi(:,4) = gl3yi(:,4) + shglb(:,3,n) * yl(:,n,4)
            gl3yi(:,5) = gl3yi(:,5) + shglb(:,3,n) * yl(:,n,5)
          enddo
c
c.... convert local-grads to global-grads
c
          g1yi(:,2) = dxidxb(:,1,1) * gl1yi(:,2) + 
     &                dxidxb(:,2,1) * gl2yi(:,2) +
     &                dxidxb(:,3,1) * gl3yi(:,2)
          g2yi(:,2) = dxidxb(:,1,2) * gl1yi(:,2) + 
     &                dxidxb(:,2,2) * gl2yi(:,2) +
     &                dxidxb(:,3,2) * gl3yi(:,2)
          g3yi(:,2) = dxidxb(:,1,3) * gl1yi(:,2) + 
     &                dxidxb(:,2,3) * gl2yi(:,2) +
     &                dxidxb(:,3,3) * gl3yi(:,2)
c
          g1yi(:,3) = dxidxb(:,1,1) * gl1yi(:,3) + 
     &                dxidxb(:,2,1) * gl2yi(:,3) +
     &                dxidxb(:,3,1) * gl3yi(:,3)
          g2yi(:,3) = dxidxb(:,1,2) * gl1yi(:,3) + 
     &                dxidxb(:,2,2) * gl2yi(:,3) +
     &                dxidxb(:,3,2) * gl3yi(:,3)
          g3yi(:,3) = dxidxb(:,1,3) * gl1yi(:,3) + 
     &                dxidxb(:,2,3) * gl2yi(:,3) +
     &                dxidxb(:,3,3) * gl3yi(:,3)
c
          g1yi(:,4) = dxidxb(:,1,1) * gl1yi(:,4) + 
     &                dxidxb(:,2,1) * gl2yi(:,4) +
     &                dxidxb(:,3,1) * gl3yi(:,4)
          g2yi(:,4) = dxidxb(:,1,2) * gl1yi(:,4) + 
     &                dxidxb(:,2,2) * gl2yi(:,4) +
     &                dxidxb(:,3,2) * gl3yi(:,4)
          g3yi(:,4) = dxidxb(:,1,3) * gl1yi(:,4) + 
     &                dxidxb(:,2,3) * gl2yi(:,4) +
     &                dxidxb(:,3,3) * gl3yi(:,4)
c
          g1yi(:,5) = dxidxb(:,1,1) * gl1yi(:,5) + 
     &                dxidxb(:,2,1) * gl2yi(:,5) +
     &                dxidxb(:,3,1) * gl3yi(:,5)
          g2yi(:,5) = dxidxb(:,1,2) * gl1yi(:,5) + 
     &                dxidxb(:,2,2) * gl2yi(:,5) +
     &                dxidxb(:,3,2) * gl3yi(:,5)
          g3yi(:,5) = dxidxb(:,1,3) * gl1yi(:,5) + 
     &                dxidxb(:,2,3) * gl2yi(:,5) +
     &                dxidxb(:,3,3) * gl3yi(:,5)
c
c.... end grad-v
c
          !Compute the gradient of the shape function for heat flux's 
          !contribution to lhsk
          if(iLHScond > 0) then
            dNadx = zero
            
            !dNdx(a,i) = dN_a / dx_i

            do aa = 1, nshl  !TODO: get rid of the intermediary dNadxi
                             !shglb(:,nsd,a=1)* N(:,a=1)
              dNadxi(:,aa,1) = shglb(:,1,aa) * 1 !would normally be a sum over
              dNadxi(:,aa,2) = shglb(:,2,aa) * 1 !all nodes, but N = 0 for a /= 1
              dNadxi(:,aa,3) = shglb(:,3,aa) * 1
            enddo 
                        
            do aa = 1, nshl
              dNadx(:,aa,1) = dNadxi(:,aa,1) * dxidxb(:,1,1) + 
     &                        dNadxi(:,aa,2) * dxidxb(:,2,1) +
     &                        dNadxi(:,aa,3) * dxidxb(:,3,1) 
              dNadx(:,aa,2) = dNadxi(:,aa,1) * dxidxb(:,1,2) + 
     &                        dNadxi(:,aa,2) * dxidxb(:,2,2) +
     &                        dNadxi(:,aa,3) * dxidxb(:,3,2) 
              dNadx(:,aa,3) = dNadxi(:,aa,1) * dxidxb(:,1,3) + 
     &                        dNadxi(:,aa,2) * dxidxb(:,2,3) +
     &                        dNadxi(:,aa,3) * dxidxb(:,3,3) 
            enddo
          endif
         
        endif
c
c.... -------------------->  Boundary Conditions  <--------------------
c
c.... compute the Euler boundary conditions
c
        rou = zero
        p   = zero
c
        do n = 1, nshlb
          nodlcl = lnode(n)
c
          rou = rou + shpb(:,nodlcl) * BCB(:,n,1)
          p   = p   + shpb(:,nodlcl) * BCB(:,n,2)
        enddo
c
c.... compute the Navier-Stokes boundary conditions
c
        if (Navier .eq. 1) then
c
          tau1n = zero
          tau2n = zero
          tau3n = zero
          heat  = zero
c
          do n = 1, nshlb
            nodlcl = lnode(n)
c
            tau1n = tau1n + shpb(:,nodlcl) * BCB(:,n,3)
            tau2n = tau2n + shpb(:,nodlcl) * BCB(:,n,4)
            tau3n = tau3n + shpb(:,nodlcl) * BCB(:,n,5)
            heat  = heat  + shpb(:,nodlcl) * BCB(:,n,6)
          enddo
c
c.... flop count
c
!      flops = flops + (184+30*nshl+8*nshlb)*npro
c
        endif
c
c.... flop count
c
!      flops = flops + (27+18*nshl+14*nshlb)*npro
c
c.... return
c
        return
        end
c
c
c
        subroutine e3bvarSclr(ycl,      BCB,     shpb,    shglb, 
     &                        xlb,     lnode,
     &                        u1,      u2,      u3,
     &                        g1yti,   g2yti,   g3yti,   WdetJb,
     &                        bnorm,   T,       rho,     cp,      rou,
     &                        Sclr,    SclrF)
c
c----------------------------------------------------------------------
c
c   This routine computes the variables at integration points for 
c the boundary element routine.
c
c input:
c  ycl     (npro,nshl,ndof)      : Y variables
c  BCB    (npro,nenbl,ndBCB)    : Boundary Condition values
c  shpb   (npro,nen)            : boundary element shape-functions
c  shglb  (nsd,nen)             : boundary element grad-shape-functions
c  xlb    (npro,nshl,nsd)       : nodal coordinates at current step
c  lnode  (nenb)                : local nodes on the boundary
c
c output:
c  g1yti  (npro)
c  g2yti  (npro)
c  g3yti  (npro)
c  WdetJb (npro)                : weighted Jacobian
c  bnorm  (npro,nsd)            : outward normal
c  T      (npro)                : temperature
c  rho    (npro)                : density
c  cp     (npro)                : specific energy at constant pressure
c  rou    (npro)                : BC mass flux
c  SclrF  (npro)                : BC Scalar  flux
c
c Zdenek Johan, Summer 1990.  (Modified from e2bvar.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension ycl(npro,nshl,ndof),        BCB(npro,nshlb,ndBCB),
     &            shpb(npro,nshl),           shglb(npro,nsd,nshl),
     &            xlb(npro,nshl,nsd),        
     &            lnode(27),                 
     &            g1yti(npro),               g2yti(npro),
     &            g3yti(npro),
     &            WdetJb(npro),              bnorm(npro,nsd),
     &            pres(npro),                T(npro),
     &            u1(npro),                  u2(npro),
     &            u3(npro),                  rho(npro),
     &            ei(npro),                  cp(npro),
     &            rk(npro),                  Sclr(npro),
     &            rou(npro),
     &            SclrF(npro)
c
        dimension dxdxib(npro,nsd,nsd),
     &            dxidxb(npro,nsd,nsd),      temp(npro),
     &            temp1(npro),               temp2(npro),
     &            temp3(npro),               gl1yti(npro),
     &            gl2yti(npro),              gl3yti(npro)
c
        dimension h(npro),                   cv(npro),
     &            alfap(npro),               betaT(npro),
     &            gamb(npro),                c(npro),
     &            tmp(npro),                 v1(npro,nsd),   
     &            v2(npro,nsd)
c
c.... ------------------->  integration variables  <--------------------
c
c.... compute the primitive variables at the integration point
c
        pres = zero
        u1   = zero
        u2   = zero
        u3   = zero
        T    = zero
        Sclr = zero

        id  = isclr+5
        ibb = isclr+6
c
        do n = 1, nshlb
          nodlcl = lnode(n)
c
          pres = pres + shpb(:,nodlcl) * ycl(:,nodlcl,1)
          u1   = u1   + shpb(:,nodlcl) * ycl(:,nodlcl,2)
          u2   = u2   + shpb(:,nodlcl) * ycl(:,nodlcl,3)
          u3   = u3   + shpb(:,nodlcl) * ycl(:,nodlcl,4)
          T    = T    + shpb(:,nodlcl) * ycl(:,nodlcl,5)
          Sclr = Sclr + shpb(:,nodlcl) * ycl(:,nodlcl,id)
        enddo
c
c.... calculate the specific kinetic energy
c
        rk = pt5 * ( u1**2 + u2**2  + u3**2 )
c
c.... get the thermodynamic properties
c
        ithm = 6
        if (Navier .eq. 1) ithm = 7
        call getthm (pres,            T,                  Sclr,
     &               rk,              rho,                ei,
     &               h,               tmp,                cv,
     &               cp,              alfap,              betaT,
     &               gamb,            c)
c
       if (iconvsclr.eq.2) rho=one
c
c.... ---------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
        dxdxib = zero
c
        do n = 1, nshl
           dxdxib(:,1,1) = dxdxib(:,1,1) + xlb(:,n,1) * shglb(:,1,n)
           dxdxib(:,1,2) = dxdxib(:,1,2) + xlb(:,n,1) * shglb(:,2,n)
           dxdxib(:,1,3) = dxdxib(:,1,3) + xlb(:,n,1) * shglb(:,3,n)
           dxdxib(:,2,1) = dxdxib(:,2,1) + xlb(:,n,2) * shglb(:,1,n)
           dxdxib(:,2,2) = dxdxib(:,2,2) + xlb(:,n,2) * shglb(:,2,n)
           dxdxib(:,2,3) = dxdxib(:,2,3) + xlb(:,n,2) * shglb(:,3,n)
           dxdxib(:,3,1) = dxdxib(:,3,1) + xlb(:,n,3) * shglb(:,1,n)
           dxdxib(:,3,2) = dxdxib(:,3,2) + xlb(:,n,3) * shglb(:,2,n)
           dxdxib(:,3,3) = dxdxib(:,3,3) + xlb(:,n,3) * shglb(:,3,n)
        enddo
c
c.... compute the normal to the boundary
c
c
        v1 = xlb(:,2,:) - xlb(:,1,:)
        v2 = xlb(:,3,:) - xlb(:,1,:)
c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
c
c.....The following are done in order to correct temp1..3  
c     based on the results from compressible code.  This is done only 
c     for wedges, depending on the bounary face.(tri or quad)  
c        
        if (lcsyst .eq. 4) then
           temp1 = dxdxib(:,2,1) * dxdxib(:,3,3) -
     &             dxdxib(:,2,3) * dxdxib(:,3,1)
           temp2 = dxdxib(:,3,1) * dxdxib(:,1,3) -
     &             dxdxib(:,3,3) * dxdxib(:,1,1)
           temp3 = dxdxib(:,1,1) * dxdxib(:,2,3) -
     &             dxdxib(:,1,3) * dxdxib(:,2,1)
           
        elseif (lcsyst .eq. 1) then
           temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
           temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
           temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
        else 
           temp1 = - v1(:,2) * v2(:,3) + v2(:,2) * v1(:,3)
           temp2 = - v2(:,1) * v1(:,3) + v1(:,1) * v2(:,3)
           temp3 = - v1(:,1) * v2(:,2) + v2(:,1) * v1(:,2)
        endif
c
        temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
        bnorm(:,1) = temp1 * temp
        bnorm(:,2) = temp2 * temp
        bnorm(:,3) = temp3 * temp
c

        if (lcsyst .eq. 3) then
           WdetJb     = (1 - Qwtb(lcsyst,intp)) / (four*temp)
        elseif (lcsyst .eq. 4) then
          WdetJb     = Qwtb(lcsyst,intp)/ temp
        else
           WdetJb     =Qwtb(lcsyst,intp) / (four*temp)
        endif
c
c.... -------------------------->  Grad-V  <----------------------------
c
c.... compute grad-v for Navier-Stokes terms
c
        if (Navier .eq. 1) then
c
c.... compute the inverse of deformation gradient
c
          dxidxb(:,1,1) =   dxdxib(:,2,2) * dxdxib(:,3,3) 
     &                    - dxdxib(:,3,2) * dxdxib(:,2,3)
          dxidxb(:,1,2) =   dxdxib(:,3,2) * dxdxib(:,1,3) 
     &                    - dxdxib(:,1,2) * dxdxib(:,3,3)
          dxidxb(:,1,3) =   dxdxib(:,1,2) * dxdxib(:,2,3) 
     &                    - dxdxib(:,1,3) * dxdxib(:,2,2)
          temp          = one / ( dxidxb(:,1,1) * dxdxib(:,1,1) 
     &                          + dxidxb(:,1,2) * dxdxib(:,2,1)  
     &                          + dxidxb(:,1,3) * dxdxib(:,3,1) )
          dxidxb(:,1,1) =  dxidxb(:,1,1) * temp
          dxidxb(:,1,2) =  dxidxb(:,1,2) * temp
          dxidxb(:,1,3) =  dxidxb(:,1,3) * temp
          dxidxb(:,2,1) = (dxdxib(:,2,3) * dxdxib(:,3,1) 
     &                   - dxdxib(:,2,1) * dxdxib(:,3,3)) * temp
          dxidxb(:,2,2) = (dxdxib(:,1,1) * dxdxib(:,3,3) 
     &                   - dxdxib(:,3,1) * dxdxib(:,1,3)) * temp
          dxidxb(:,2,3) = (dxdxib(:,2,1) * dxdxib(:,1,3) 
     &                   - dxdxib(:,1,1) * dxdxib(:,2,3)) * temp
          dxidxb(:,3,1) = (dxdxib(:,2,1) * dxdxib(:,3,2) 
     &                   - dxdxib(:,2,2) * dxdxib(:,3,1)) * temp
          dxidxb(:,3,2) = (dxdxib(:,3,1) * dxdxib(:,1,2) 
     &                   - dxdxib(:,1,1) * dxdxib(:,3,2)) * temp
          dxidxb(:,3,3) = (dxdxib(:,1,1) * dxdxib(:,2,2) 
     &                   - dxdxib(:,1,2) * dxdxib(:,2,1)) * temp

c
c.... compute local-grad-Y
c
          gl1yti = zero
          gl2yti = zero
          gl3yti = zero
c
          do n = 1, nshl
            gl1yti(:) = gl1yti(:) + shglb(:,1,n) * ycl(:,n,id)
            gl2yti(:) = gl2yti(:) + shglb(:,2,n) * ycl(:,n,id)
            gl3yti(:) = gl3yti(:) + shglb(:,3,n) * ycl(:,n,id)
          enddo
c
c.... convert local-grads to global-grads
c
          g1yti(:) = dxidxb(:,1,1) * gl1yti(:) + 
     &               dxidxb(:,2,1) * gl2yti(:) +
     &               dxidxb(:,3,1) * gl3yti(:)
          g2yti(:) = dxidxb(:,1,2) * gl1yti(:) + 
     &               dxidxb(:,2,2) * gl2yti(:) +
     &               dxidxb(:,3,2) * gl3yti(:)
          g3yti(:) = dxidxb(:,1,3) * gl1yti(:) + 
     &               dxidxb(:,2,3) * gl2yti(:) +
     &               dxidxb(:,3,3) * gl3yti(:)

c
c.... end grad-Sclr
        endif
c
c.... -------------------->  Boundary Conditions  <--------------------
c
c.... compute the Euler boundary conditions
c
        rou = zero
        do n = 1, nshlb
          nodlcl = lnode(n)
          rou = rou + shpb(:,nodlcl) * BCB(:,n,1)
        enddo
c
c.... impose scalar flux boundary conditions
        SclrF = zero
        do n=1,nshlb
          nodlcl = lnode(n)
          SclrF = SclrF + shpb(:,nodlcl) * BCB(:,n,ibb)
        enddo

c
c.... flop count
c
!      flops = flops + (27+18*nshl+14*nenbl)*npro
c
c.... return
c
        return
        end

