       subroutine e3b (yl,      ycl,  iBCB,    BCB,     shpb,    shglb,
     &                 xlb,     rl,      rml,     sgn)
c
c----------------------------------------------------------------------
c
c   This routine calculates the 3D RHS residual of the fluid boundary 
c   elements.
c
c input:
c  yl     (npro,nshl,nflow)     : Y variables  (perturbed, no scalars)
c  ycl    (npro,nshl,ndof)      : Y variables
c  iBCB   (npro,ndiBCB)         : boundary condition code (iBCB(:,1) is
c      a bit tested boundary integral flag i.e.
c                  if set to value of BCB      if set to floating value
c      iBCB(:,1) : convective flux * 1            0  (ditto to all below)
c                  pressure   flux * 2
c                  viscous    flux * 4
c                  heat       flux * 8
c                  turbulence wall * 16
c                  scalarI   flux  * 16*2^I 
c                  (where I is the scalar number)
c
c      iBCB(:,2) is the srfID given by the user in MGI that we will
c                collect integrated fluxes for.
c
c  BCB    (npro,nshlb,ndBCB)    : Boundary Condition values
c                                  BCB (1) : mass flux 
c                                  BCB (2) : pressure 
c                                  BCB (3) : viscous flux in x1-direc.
c                                  BCB (4) : viscous flux in x2-direc.
c                                  BCB (5) : viscous flux in x3-direc.
c                                  BCB (6) : heat flux
c  shpb   (nshl,ngaussb)           : boundary element shape-functions
c  shglb  (nsd,nshl,ngaussb)       : boundary element grad-shape-functions
c  xlb    (npro,nenl,nsd)       : nodal coordinates at current step
c
c output:
c  rl     (npro,nshl,nflow)      : element residual
c  rml    (npro,nshl,nflow)      : element modified residual
c
c
c Note: Always the first side of the element is on the boundary.  
c       However, note that for higher-order elements the nodes on 
c       the boundary side are not the first nshlb nodes, see the 
c       array lnode.
c
c
c Zdenek Johan, Summer 1990.  (Modified from e2b.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c Anilkumar Karanam Spring 2000 (Modified for Hierarchic Hexes)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension yl(npro,nshl,nflow),          iBCB(npro,ndiBCB),
     &            ycl(npro,nshl,ndof),
     &            BCB(npro,nshlb,ndBCB),       shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),         
     &            xlb(npro,nenl,nsd),          
     &            rl(npro,nshl,nflow),          rml(npro,nshl,nflow)
c
        dimension g1yi(npro,nflow),             g2yi(npro,nflow),
     &            g3yi(npro,nflow),             WdetJb(npro),
     &            bnorm(npro,nsd)
c
        dimension un(npro),                    rk(npro),
     &            u1(npro),                    u2(npro),
     &            u3(npro),                    
     &            rho(npro),                   pres(npro),
     &            T(npro),                     ei(npro),
     &            cp(npro)
c
        dimension rou(npro),                   p(npro),
     &            F1(npro),                    F2(npro),
     &            F3(npro),                    F4(npro),
     &            F5(npro),                    Fv2(npro),
     &            Fv3(npro),                   Fv4(npro),
     &            Fv5(npro),                   Fh5(npro)
c
        dimension rmu(npro),                   rlm(npro),
     &            rlm2mu(npro),                con(npro),
     &            tau1n(npro),
     &            tau2n(npro),                 tau3n(npro),
     &            heat(npro)
c
        dimension lnode(27),               sgn(npro,nshl),
     &       shape(npro,nshl),        shdrv(npro,nsd,nshl)
c
        dimension xmudum(npro,ngauss)

        ttim(40) = ttim(40) - secs(0.0)

c
c.... compute the nodes which lie on the boundary
c
        call getbnodes(lnode)

c.... loop through the integration points

        if(lcsyst.eq.3.or.lcsyst.eq.4) then
           ngaussb = nintb(lcsyst)
        else
           ngaussb = nintb(lcsyst)
        endif

        do intp = 1, ngaussb
c
c.... if Det. .eq. 0, do not include this point
c
        if (Qwtb(lcsyst,intp) .eq. zero) cycle         ! precaution
c
c.... create a matrix of shape functions (and derivatives) for each
c     element at this quadrature point. These arrays will contain 
c     the correct signs for the hierarchic basis
c     
c
        call getshpb(shpb,        shglb,        sgn, 
     &              shape,       shdrv)
c     
c.... calculate the integration variables
c
        call e3bvar (yl,   ycl,           BCB,          shape,
     &               shdrv,           xlb,
     &               lnode,           g1yi,
     &               g2yi,            g3yi,         WdetJb,
     &               bnorm,           pres,         T,
     &               u1,              u2,           u3,
     &               rho,             ei,           cp,
     &               rk,              rou,          p,
     &               Fv2,             Fv3,          Fv4,
     &               Fh5)
c
c.... ires = 1 or 3
c
        if ((ires .eq. 1) .or. (ires .eq. 3)) then
c
c.... clear some variables
c
          tau1n = zero
          tau2n = zero
          tau3n = zero
          heat  = zero
c
c.... ------------------------->  convective  <------------------------
c
c
          where (.not.btest(iBCB(:,1),0) )
            un  = bnorm(:,1) * u1 + bnorm(:,2) * u2 + bnorm(:,3) * u3
            rou = rho * ( un )
          elsewhere
            un  = (rou / rho) 
          endwhere
c
c.... ------------------------->  pressure  <--------------------------
c
c.... use one-point quadrature in time
c
          where (.not.btest(iBCB(:,1),1)) p = pres
c
c.... add the Euler contribution
c
          F1 = rou
          F2 = rou * u1        + bnorm(:,1) * p
          F3 = rou * u2        + bnorm(:,2) * p
          F4 = rou * u3        + bnorm(:,3) * p
          F5 = rou * (ei + rk) +         un * p
c
c.... flop count
c
          flops = flops + 23*npro
c
c.... end of ires = 1 or 3
c
        endif
c
c.... ----------------------->  Navier-Stokes  <-----------------------
c
        if (Navier .eq. 1) then

           xmudum = zero

c
c.... get the material properties
c
        call getDiff (T,        cp,    rho,        ycl,
     &                rmu,      rlm,   rlm2mu,     con, shape,
     &                xmudum,   xl)
c
c.... ------------------------>  viscous flux <------------------------
c
c.... floating viscous flux
c
        tau1n = bnorm(:,1) * (rlm2mu* g1yi(:,2) + rlm   *g2yi(:,3) 
     &                                          + rlm   *g3yi(:,4))
     &        + bnorm(:,2) * (rmu   *(g2yi(:,2) + g1yi(:,3)))
     &        + bnorm(:,3) * (rmu   *(g3yi(:,2) + g1yi(:,4)))
        tau2n = bnorm(:,1) * (rmu   *(g2yi(:,2) + g1yi(:,3)))
     &        + bnorm(:,2) * (rlm   * g1yi(:,2) + rlm2mu*g2yi(:,3) 
     &                                          + rlm   *g3yi(:,4))
     &        + bnorm(:,3) * (rmu   *(g3yi(:,3) + g2yi(:,4)))
        tau3n = bnorm(:,1) * (rmu   *(g3yi(:,2) + g1yi(:,4)))
     &        + bnorm(:,2) * (rmu   *(g3yi(:,3) + g2yi(:,4)))
     &        + bnorm(:,3) * (rlm   * g1yi(:,2) + rlm   *g2yi(:,3) 
     &                                          + rlm2mu*g3yi(:,4))
c
        where (.not.btest(iBCB(:,1),2) )
           Fv2 = tau1n          ! wherever traction is not set, use the
           Fv3 = tau2n          ! viscous flux calculated from the field
           Fv4 = tau3n          !
        endwhere
c
        Fv5 = u1 * Fv2 + u2 * Fv3 + u3 * Fv4
c
c.... -------------------------->  heat flux <-------------------------
c
c.... floating heat flux
c
        heat =   con * ( bnorm(:,1) * g1yi(:,5) +
     &                   bnorm(:,2) * g2yi(:,5) +
     &                   bnorm(:,3) * g3yi(:,5) ) 
c
        where (.not.btest(iBCB(:,1),3) ) Fh5 = heat
c
c.... add the Navier-Stokes contribution
c
        F2  = F2 - Fv2
        F3  = F3 - Fv3
        F4  = F4 - Fv4
        F5  = F5 - Fv5 + Fh5
c
c.... flop count
c
        flops = flops + 27*npro
c
c.... end of Navier Stokes part
c
        endif
c
c.... ------------------------->  Residual  <--------------------------
c
c.... add the flux to the residual
c
        if ((ires .eq. 1) .or. (ires .eq. 3)) then
c
c
          do n = 1, nshlb
            nodlcl = lnode(n)
c
            rl(:,nodlcl,1) = rl(:,nodlcl,1)
     &                     + WdetJb * shape(:,nodlcl) * F1
            rl(:,nodlcl,2) = rl(:,nodlcl,2)
     &                     + WdetJb * shape(:,nodlcl) * F2
            rl(:,nodlcl,3) = rl(:,nodlcl,3)
     &                     + WdetJb * shape(:,nodlcl) * F3
            rl(:,nodlcl,4) = rl(:,nodlcl,4)
     &                     + WdetJb * shape(:,nodlcl) * F4
            rl(:,nodlcl,5) = rl(:,nodlcl,5)
     &                     + WdetJb * shape(:,nodlcl) * F5
          enddo
c
          flops = flops + 12*nshlb*npro
c
        endif
c
c.... add the flux to the modified residual
c
        if (((ires .eq. 2) .or. (ires .eq. 3))
     &      .and. (Navier .eq. 1) .and. (Jactyp .eq. 1)) then
c
          do n = 1, nshlb
            nodlcl = lnode(n)
c
            rml(:,nodlcl,2) = rml(:,nodlcl,2) - WdetJb *
     &                        shape(:,nodlcl) *  Fv2
            rml(:,nodlcl,3) = rml(:,nodlcl,3) - WdetJb *
     &                        shape(:,nodlcl) *  Fv3
            rml(:,nodlcl,4) = rml(:,nodlcl,4) - WdetJb *
     &                        shape(:,nodlcl) *  Fv4
            rml(:,nodlcl,5) = rml(:,nodlcl,5) - WdetJb *
     &                        shape(:,nodlcl) * (Fv5 - Fh5)
          enddo
c
          flops = flops + 11*nenbl*npro
c
        endif
c
c  uncomment and run 1 step to get estimate of wall shear vs z
c
c        do i=1,npro
c           tnorm= sqrt(tau1n(i)*tau1n(i)
c     &                +tau2n(i)*tau2n(i)+tau3n(i)*tau3n(i))
c          
c           write(700+myrank,*) xlb(i,1,3),tnorm
c        enddo
        

       do iel = 1, npro
c
c  if we have a nonzero value then
c  calculate the fluxes through this surface 
c
           iface = abs(iBCB(iel,2))
           if (iface .ne. 0 .and. ires.ne.2) then
              flxID(1,iface) =  flxID(1,iface) + WdetJb(iel)! measure area too
c              flxID(2,iface) =  flxID(2,iface) - WdetJb(iel) * un(iel)
              flxID(2,iface) =  flxID(2,iface) - WdetJb(iel) * rou(iel)
              flxID(3,iface) = flxID(3,iface)
     &                   - ( tau1n(iel) - bnorm(iel,1)*pres(iel))
     &                   * WdetJb(iel) 
              flxID(4,iface) = flxID(4,iface)
     &                   - ( tau2n(iel) - bnorm(iel,2)*pres(iel))
     &                   * WdetJb(iel) 
              flxID(5,iface) = flxID(5,iface)
     &                   - ( tau3n(iel) - bnorm(iel,3)*pres(iel))
     &                   * WdetJb(iel) 

           endif
        enddo

c
c.... -------------------->  Aerodynamic Forces  <---------------------
c
        if ((ires .ne. 2) .and. (iter .eq. nitr)) then
c
c.... compute the forces on the body
c
          where (.not.btest(iBCB(:,1),0) )
            tau1n = ( pres * bnorm(:,1) - tau1n ) * WdetJb
            tau2n = ( pres * bnorm(:,2) - tau2n ) * WdetJb
            tau3n = ( pres * bnorm(:,3) - tau3n ) * WdetJb
            heat  = - heat * WdetJb
          elsewhere
            tau1n = zero
            tau2n = zero
            tau3n = zero
            heat  = zero
          endwhere
c
          Force(1) = Force(1) + sum(tau1n)
          Force(2) = Force(2) + sum(tau2n)
          Force(3) = Force(3) + sum(tau3n)
          HFlux    = HFlux    + sum(heat)
c
        endif
c
c.... end of integration loop
c
        enddo

        ttim(40) = ttim(40) + secs(0.0)
c
c.... return
c
        return
        end
c
c
c
        subroutine e3bSclr (ycl,      iBCB,     BCB,    
     &                      shpb,    shglb,    sgn, 
     &                      xlb,     rtl,      rmtl)
c
c----------------------------------------------------------------------
c
c   This routine calculates the 3D RHS residual of the fluid boundary 
c   elements.
c
c input:
c  ycl     (npro,nshl,ndof)      : Y variables
c  iBCB   (npro,ndiBCB)         : boundary condition code & surfID
c  BCB    (npro,nenbl,ndBCB)    : Boundary Condition values
c                                  BCB (1) : mass flux 
c                                  BCB (2) : pressure 
c                                  BCB (3) : viscous flux in x1-direc.
c                                  BCB (4) : viscous flux in x2-direc.
c                                  BCB (5) : viscous flux in x3-direc.
c                                  BCB (6) : heat flux
c                                  BCB (7) : eddy visc flux
c  shpb   (nen,nintg)           : boundary element shape-functions
c  shglb  (nsd,nen,nintg)       : boundary element grad-shape-functions
c  xlb    (npro,nenl,nsd)       : nodal coordinates at current step
c
c output:
c  rtl    (npro,nenl,nflow)      : element residual
c  rmtl   (npro,nenl,nflow)      : element modified residual
c
c
c Note: Always the first side of the element is on the boundary.  
c       However, note that for higher-order elements the nodes on 
c       the boundary side are not the first nenbl nodes, see the 
c       array mnodeb.
c
c
c Zdenek Johan, Summer 1990.  (Modified from e2b.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use turbSA ! for saSigma
        include "common.h"
c
        dimension ycl(npro,nshl,ndof),        iBCB(npro,ndiBCB),
     &            BCB(npro,nshlb,ndBCB),   shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),   sgn(npro,nshl),         
     &            xlb(npro,nenl,nsd),        shape(npro,nshl),  
     &            rtl(npro,nshl),          rmtl(npro,nshl),
     &            shdrv(npro,nsd,nshl)
c
        dimension u1(npro),                    u2(npro),
     &            u3(npro),
     &            g1yti(npro),                 g2yti(npro),
     &            g3yti(npro),                 WdetJb(npro),
     &            bnorm(npro,nsd)
c
        dimension rho(npro),                   Sclr(npro),
     &            T(npro),                     cp(npro)
c
        dimension rou(npro),                   F(npro),
     &            un(npro),                    Sclrn(npro)
c
        dimension rmu(npro),                   rlm(npro),
     &            rlm2mu(npro),                con(npro),
     &            heat(npro),                  srcp(npro)
c
        dimension lnode(27)
        real*8  sign_levelset(npro), sclr_ls(npro), mytmp(npro) 
        ttim(40) = ttim(40) - tmr()
c
c.... get the boundary nodes
c
       id  = isclr + 5
       ib  = isclr + 4
       ibb = isclr + 6
       call getbnodes(lnode)
c
c.... loop through the integration points
c
        ngaussb = nintb(lcsyst)
c        
        do intp = 1, ngaussb
c
c.... if Det. .eq. 0, do not include this point
c
        if (Qwtb(lcsyst,intp) .eq. zero) cycle         ! precaution
c
c.... create a matrix of shape functions (and derivatives) for each
c     element at this quadrature point. These arrays will contain 
c     the correct signs for the hierarchic basis
c 
        call getshpb(shpb,        shglb,        sgn, 
     &       shape,       shdrv)
c
c.... calculate the integraton variables
c
        call e3bvarSclr (ycl,            BCB,
     &                   shape,         shdrv,        
     &                   xlb,           lnode,         u1,
     &                   u2,            u3,            g1yti,
     &                   g2yti,         g3yti,         WdetJb,
     &                   bnorm,         T,             rho,
     &                   cp,            rou,           Sclr,
     &                   F)
c.......********************modification for Ilset=2**********************
          if (ilset.eq.2 .and. isclr.eq.2) then !we are redistancing level-sets

CAD   If Sclr(:,1).gt.zero, result of sign_term function 1
CAD   If Sclr(:,1).eq.zero, result of sign_term function 0
CAD   If Sclr(:,1).lt.zero, result of sign_term function -1

            sclr_ls = zero      !zero out temp variable

            do ii=1,npro

           do jj = 1, nshl  ! first find the value of levelset at point ii
                  
                  sclr_ls(ii) =  sclr_ls(ii) 
     &                        + shape(ii,jj) * ycl(ii,jj,6)

               enddo
               if (sclr_ls(ii) .lt. - epsilon_ls)then
            
                  sign_levelset(ii) = - one
                  
               elseif  (abs(sclr_ls(ii)) .le. epsilon_ls)then
c                 
                 sign_levelset(ii) =sclr_ls(ii)/epsilon_ls 
     &              + sin(pi*sclr_ls(ii)/epsilon_ls)/pi
                  
               elseif (sclr_ls(ii) .gt. epsilon_ls) then
                  
                  sign_levelset(ii) = one
                  
               endif               
c               
               srcp(ii) = sign_levelset(ii)
c               
            enddo  
c
cad   The redistancing equation can be written in the following form
cad
cad   d_{,t} + sign(phi)*( d_{,i}/|d_{,i}| )* d_{,i} = sign(phi)
cad
cad   This is rewritten in the form
cad
cad   d_{,t} + u * d_{,i} = sign(phi)
cad

c$$$CAD   For the redistancing equation the "pseudo velocity" term is
c$$$CAD   calculated as follows



            mytmp = srcp / sqrt   ( g1yti * g1yti 
     &                            + g2yti * g2yti
     &                            + g3yti * g3yti) 

            u1 = mytmp * g1yti 
            u2 = mytmp * g2yti
            u3 = mytmp * g3yti
        endif
         
c
c.... ires = 1 or 3
c
        if ((ires .eq. 1) .or. (ires .eq. 3)) then

c.... ------------------------->  convective  <------------------------
c
           where (.not.btest(iBCB(:,1),0) )
              un  = bnorm(:,1)*u1 + bnorm(:,2)*u2 + bnorm(:,3)*u3
              rou = rho * ( un )
           elsewhere
              un  = (rou / rho) 
           endwhere
c
c.... calculate flux where unconstrained
c
           where (.not.btest(iBCB(:,1),ib) )
              F = Sclr *rou
           endwhere
c
c.... get the material properties
c

        call getDiffSclr (T,          cp,         rmu,
     &                    rlm,        rlm2mu,     con, rho, Sclr)

c
c.... ----------> DiffFlux for Scalar Variable  <--------
c
        if (ilset.ne.2) then

           where (.not.btest(iBCB(:,1),ib) )
              Sclrn  = bnorm(:,1) * g1yti(:) 
     &               + bnorm(:,2) * g2yti(:)
     &               + bnorm(:,3) * g3yti(:)
C
           
c             F = F + rmu*Sclrn  !!!! CHECK THIS 

          F = F + saSigmaInv*rho*((rmu/rho)+Sclr)*Sclrn!confirm the modificationc                                                      in getdiffsclr

c.....this modification of viscosity goes in getdiffsclr

           endwhere
        endif
c
c.... end of ires = 1 or 3
c
        endif
c
c.... ------------------------->  Residual  <--------------------------
c
c.... add the flux to the residual
c
        if ((ires .eq. 1) .or. (ires .eq. 3)) then
           if (iconvsclr.eq.1) then    !conservative boundary integral
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 rtl(:,nodlcl) = rtl(:,nodlcl)
     &                         + WdetJb * shape(:,nodlcl) * F
              enddo
              flops = flops + 12*nshlb*npro
           endif
        endif
c
c.... add the flux to the modified residual
c
c        if (((ires .eq. 2) .or. (ires .eq. 3))
c     &      .and. (Navier .eq. 1) .and. (Jactyp .eq. 1)) then
c
c          do n = 1, nenbl
c            nodlcl = lnode(n)
c
c            rml(:,nodlcl,2) = rml(:,nodlcl,2) - WdetJb *
c     &                        shpb(nodlcl,intp) *  Fv2
c            rml(:,nodlcl,3) = rml(:,nodlcl,3) - WdetJb *
c     &                        shpb(nodlcl,intp) *  Fv3
c            rml(:,nodlcl,4) = rml(:,nodlcl,4) - WdetJb *
c     &                        shpb(nodlcl,intp) *  Fv4
c            rml(:,nodlcl,5) = rml(:,nodlcl,5) - WdetJb *
c     &                        shpb(nodlcl,intp) * (Fv5 - Fh5)
c          enddo
c
c          flops = flops + 11*nenbl*npro
c
c        endif
c

c
c.... end of integration loop
c
        enddo

        ttim(40) = ttim(40) + tmr()
c
c.... return
c
        return
        end

