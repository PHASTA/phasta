        subroutine e3b (ul,      yl,      acl,     iBCB,    BCB,     
     &                  shpb,    shglb,
     &                  xlb,     rl,      sgn,     dwl,     xKebe)
c
c----------------------------------------------------------------------
c
c   This routine calculates the 3D RHS residual of the fluid boundary 
c   elements.
c
c input:
c  yl     (npro,nshl,ndof)      : Y variables
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
c  shpb   (nen,ngaussb)           : boundary element shape-functions
c  shglb  (nsd,nen,ngaussb)       : boundary element grad-shape-functions
c  xlb    (npro,nenl,nsd)       : nodal coordinates at current step
c
c output:
c  rl     (npro,nshl,nflow)      : element residual
c
c Note: Always the first side of the element is on the boundary.  
c       However, note that for higher-order elements the nodes on 
c       the boundary side are not the first nshlb nodes, see the 
c       array mnodeb.
c
c
c Zdenek Johan, Summer 1990.  (Modified from e2b.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c Alberto Figueroa, Winter 2004.  CMM-FSI
c Irene Vignon, Spring 2004
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension yl(npro,nshl,ndof),          iBCB(npro,ndiBCB),
     &            BCB(npro,nshlb,ndBCB),       shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),           
     &            xlb(npro,nenl,nsd),          ul(npro,nshl,nsd),
     &            acl(npro,nshl,ndof),
     &            rl(npro,nshl,nflow)
c
        dimension g1yi(npro,ndof),             g2yi(npro,ndof),
     &            g3yi(npro,ndof),             WdetJb(npro),
     &            bnorm(npro,nsd)
c
        dimension u1(npro),                    u2(npro),
     &            u3(npro),                    rho(npro),
     &            unm(npro),                   pres(npro),
     &            vdot(npro,nsd),              rlKwall(npro,nshlb,nsd)
c
        dimension rou(npro),                   rmu(npro)
c
        dimension tau1n(npro),
     &            tau2n(npro),                 tau3n(npro)
c
        dimension lnode(27),               sgn(npro,nshl),
     &            shape(npro,nshl),        shdrv(npro,nsd,nshl),
     &            rNa(npro,4)

        real*8    xmudmi(npro,ngauss),      dwl(npro,nshl)
c
      	dimension xKebe(npro,9,nshl,nshl),  rKwall_glob(npro,9,nshl,nshl)
      	integer   intp


c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c
c.... loop through the integration points
c
        if(lcsyst.eq.3.or.lcsyst.eq.4) then
           ngaussb = nintb(lcsyst)
        else
           ngaussb = nintb(lcsyst)
        endif
        
        do intp = 1, ngaussb
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(shpb,        shglb,        sgn, 
     &              shape,       shdrv)
c
c     NOTE I DID NOT PASS THE lnode down.  It is not needed
c     since the shape functions are zero on the boundary
c
c     Note that xmudmi is not calculated at these quadrature 
c     points so you give it a zero.  This has implications.
c     the traction calculated by this approach will include 
c     molecular stresses ONLY.  This is why we will use the 
c     consistent flux method to obtain the forces when doing
c     effective viscosity wall modeling.  When doing slip velocity
c     this is not a problem since the traction is given from the
c     log law relation (not the viscosity).
c
        xmudmi=zero
c
c.... get necessary fluid properties (including eddy viscosity)
c
        call getdiff(dwl, yl,     shape,     xmudmi, xlb,   rmu, rho)
c
c.... calculate the integraton variables
c
        call e3bvar (yl,              acl,             ul,              
     &               shape,
     &               shdrv,           xlb,
     &               lnode,           WdetJb,
     &               bnorm,           pres,            
     &               u1,              u2,              u3,
     &               rmu,             unm,
     &               tau1n,           tau2n,           tau3n,
     &               vdot,            rlKwall,         
     &               xKebe,           rKwall_glob)
        
c        
c.... -----------------> boundary conditions <-------------------
c
        do iel = 1, npro
c
c  if we have a nonzero value then
c  calculate the fluxes through this surface 
c
           iface = abs(iBCB(iel,2))
!MR CHANGE
           if(iface>MAXSURF) then
            write(*,*) 'iface>MAXSURF',iface,MAXSURF
            write(*,*) 'Increase MAXSURF or decrease surfID in geom.spj'
            stop !brutal but mpi_finalize will not work without broadcasting the information to other processors.
           endif
!MR CHANGE END

           if (nsrflist(iface) .ne. 0 .and. ires.ne.2) then
              flxID(1,iface) =  flxID(1,iface) + WdetJb(iel)! measure area too
              flxID(2,iface) =  flxID(2,iface) - WdetJb(iel) * unm(iel)
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
c
c
c.... mass flux
c

           if (btest(iBCB(iel,1),0)) then
              unm(iel)  = zero
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 unm(iel) = unm(iel) 
     &                    + shape(iel,nodlcl) * BCB(iel,n,1)
              enddo
           endif
c
c.... pressure
c        
           if (btest(iBCB(iel,1),1)) then
              pres(iel) = zero
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 pres(iel) = pres(iel) 
     &                     + shape(iel,nodlcl) * BCB(iel,n,2)
              enddo
           endif
c
c.... viscous flux
c        
           if (btest(iBCB(iel,1),2)) then
              tau1n(iel) = zero
              tau2n(iel) = zero
              tau3n(iel) = zero
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 tau1n(iel) = tau1n(iel) 
     &                      + shape(iel,nodlcl)*BCB(iel,n,3)
                 tau2n(iel) = tau2n(iel) 
     &                      + shape(iel,nodlcl)*BCB(iel,n,4)
                 tau3n(iel) = tau3n(iel) 
     &                      + shape(iel,nodlcl)*BCB(iel,n,5)
              enddo
           endif
!KJ CHANGE      if(ideformwall.eq.1) then
c
c.... turbulence wall (as a way of checking for deformable wall stiffness)
c
           if (btest(iBCB(iel,1),4)) then
              rlKwall(iel,:,:) = rlKwall(iel,:,:) / ngaussb ! divide by number of gauss points 
              pres(iel) = zero                              ! to avoid the gauss point loop
              tau1n(iel) = zero                             ! and make the traction contribution
              tau2n(iel) = zero                             ! zero
              tau3n(iel) = zero                              
           else
              rlKwall(iel,:,:) = zero                       ! this is not a deformable element
              vdot(iel,:) = zero
              xKebe(iel,:,:,:) = zero
              rKwall_glob(iel,:,:,:) = zero                 ! no stiffness: not a wall element
           endif
!KJ CHANGE      endif

        enddo                                               ! end of bc loop
c
c$$$c.... if we are computing the bdry for the consistent
c$$$c     boundary forces, we must not include the surface elements
c$$$c     in the computataion (could be done MUCH more efficiently!)--->
                                                                  !this
                                                                  !comment should read as for the consistent flux calculation rather than boundary forces
c$$$c
        if (ires .eq. 2) then
           do iel = 1, npro 
              if (nsrflist(iBCB(iel,2)) .ne. 0) then
                 unm(iel) = zero
                 tau1n(iel) = zero
                 tau2n(iel) = zero
                 tau3n(iel) = zero
c                 pres(iel) = zero
c
c whatever is zeroed here will beome part of the post-processed surface
c                 "traction force"
!KJ CHANGE         if(ideformwall.eq.1) then
c
c uncomment the next two lines to get all of the t vector coming from
c                 Alberto's wall motion model.
c                 vdot(iel,:)=zero
c                 rlKwall(iel,:,:)=zero

c
c uncomment the next 8 lines to get only the tangential part
c
                  vn=dot_product(vdot(iel,:),bnorm(iel,:))
                  vdot(iel,:)=vn*bnorm(iel,:)
                  walln1=dot_product(rlkwall(iel,1,:),bnorm(iel,:))
                  walln2=dot_product(rlkwall(iel,2,:),bnorm(iel,:))
                  walln3=dot_product(rlkwall(iel,3,:),bnorm(iel,:))
                  rlkwall(iel,1,:)=walln1*bnorm(iel,:)
                  rlkwall(iel,2,:)=walln2*bnorm(iel,:)
                  rlkwall(iel,3,:)=walln3*bnorm(iel,:)
!KJ CHANGE              endif
             endif
           enddo
        endif
c
c.... assemble the contributions
c
        rNa(:,1) = -WdetJb * ( tau1n - bnorm(:,1) * pres - vdot(:,1))
        rNa(:,2) = -WdetJb * ( tau2n - bnorm(:,2) * pres - vdot(:,2))
        rNa(:,3) = -WdetJb * ( tau3n - bnorm(:,3) * pres - vdot(:,3))
        rNa(:,4) =  WdetJb * unm
c
        if(iconvflow.eq.1) then     ! conservative form was integrated
                                    ! by parts and has a convective 
                                    ! boundary integral
c
c.... assemble the contributions
c
           rou=rho*unm
           rNa(:,1) = rNa(:,1) + WdetJb * rou * u1 
           rNa(:,2) = rNa(:,2) + WdetJb * rou * u2
           rNa(:,3) = rNa(:,3) + WdetJb * rou * u3
        endif
c
c.... ------------------------->  Residual  <--------------------------
c
c.... add the flux to the residual
c
        do n = 1, nshlb
           nodlcl = lnode(n)

           rl(:,nodlcl,1) = rl(:,nodlcl,1) - shape(:,nodlcl) * rNa(:,1)
           rl(:,nodlcl,2) = rl(:,nodlcl,2) - shape(:,nodlcl) * rNa(:,2)
           rl(:,nodlcl,3) = rl(:,nodlcl,3) - shape(:,nodlcl) * rNa(:,3)
           rl(:,nodlcl,4) = rl(:,nodlcl,4) - shape(:,nodlcl) * rNa(:,4)

        enddo
        if(ideformwall.eq.1) then
           rl(:,1,1) = rl(:,1,1) - rlKwall(:,1,1)
           rl(:,1,2) = rl(:,1,2) - rlKwall(:,1,2)
           rl(:,1,3) = rl(:,1,3) - rlKwall(:,1,3)
           
           rl(:,2,1) = rl(:,2,1) - rlKwall(:,2,1)
           rl(:,2,2) = rl(:,2,2) - rlKwall(:,2,2)
           rl(:,2,3) = rl(:,2,3) - rlKwall(:,2,3)
        
           rl(:,3,1) = rl(:,3,1) - rlKwall(:,3,1)
           rl(:,3,2) = rl(:,3,2) - rlKwall(:,3,2)
           rl(:,3,3) = rl(:,3,3) - rlKwall(:,3,3)
        endif 
c
c.... -------------------->  Aerodynamic Forces  <---------------------
c
        if (((ires .ne. 2) .and. (iter .eq. nitr))
     &                     .and. (abs(itwmod).eq.1)) then
c
c.... compute the forces on the body
c
           where (nsrflist(iBCB(:,2)).eq.1)
              tau1n = ( tau1n - bnorm(:,1)*pres)  * WdetJb
              tau2n = ( tau2n - bnorm(:,2)*pres)  * WdetJb
              tau3n = ( tau3n - bnorm(:,3)*pres)  * WdetJb
           elsewhere
              tau1n = zero
              tau2n = zero
              tau3n = zero
           endwhere
c
c Note that the sign has changed from the compressible code to 
c make it consistent with the way the bflux calculates the forces
c Note also that Hflux has moved to e3btemp
c
          Force(1) = Force(1) - sum(tau1n)
          Force(2) = Force(2) - sum(tau2n)
          Force(3) = Force(3) - sum(tau3n)

c
        endif
c
c.... end of integration loop
c
        enddo
        if(ideformwall.eq.1) then
c     
c.... -----> Wall Stiffness and Mass matrices for implicit LHS  <-----------
c     
c.... Now we simply have to add the stiffness contribution in rKwall_glob to 
c.... the mass contribution already contained in xKebe

c.... this line is going to destroy the mass matrix contribution


c      xKebe = zero

         xKebe(:,:,:,:) = ( xKebe(:,:,:,:)*iwallmassfactor
     &           + rKwall_glob(:,:,:,:)*iwallstiffactor )

         endif
c$$$        ttim(40) = ttim(40) + tmr()
c
c.... return
c
        return
        end


c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c*********************************************************************
c*********************************************************************


        subroutine e3bSclr (yl,      iBCB,    BCB,     shpb,    shglb,
     &                      xlb,     rl,      sgn,     dwl)
        include "common.h"
c
        dimension yl(npro,nshl,ndof),          iBCB(npro,ndiBCB),
     &            BCB(npro,nshlb,ndBCB),       shpb(nshl,*),
     &            shglb(nsd,nshl,*),           
     &            xlb(npro,nenl,nsd),          
     &            rl(npro,nshl)
c
        real*8    WdetJb(npro),                bnorm(npro,nsd)
c
        dimension lnode(27),                   sgn(npro,nshl),
     &            shape(npro,nshl),            shdrv(npro,nsd,nshl),
     &            rNa(npro),                   flux(npro)
        real*8    dwl(npro,nshl)

c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c
c.... loop through the integration points
c
        if(lcsyst.eq.3.or.lcsyst.eq.4) then
           ngaussb = nintb(lcsyst)
        else
           ngaussb = nintb(lcsyst)
        endif
        do intp = 1, ngaussb
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(shpb,        shglb,        sgn, 
     &              shape,       shdrv)
c
c.... calculate the integraton variables
c
        call e3bvarSclr (yl,          shdrv,   xlb,
     &                   shape,       WdetJb,  bnorm,
     &                   flux,        dwl )
c        
c.... -----------------> boundary conditions <-------------------
c

c
c.... heat or scalar  flux
c     
        if(isclr.eq.0) then 
           iwalljump=0
        else
           iwalljump=1  !turb wall between heat and scalar flux..jump over
        endif
        ib=4+isclr+iwalljump
        ibb=6+isclr
        do iel=1, npro
c
c  if we have a nonzero value then
c  calculate the fluxes through this surface 
c
           if (nsrflist(iBCB(iel,2)).ne.0 .and. ires.ne.2) then
              iface = abs(iBCB(iel,2))
              flxID(ibb,iface) =  flxID(ibb,iface) 
     &                          - WdetJb(iel) * flux(iel)
           endif

           if (btest(iBCB(iel,1),ib-1)) then
              flux(iel) = zero
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 flux(iel) = flux(iel) 
     &                     + shape(iel,nodlcl) * BCB(iel,n,ibb)
              enddo           
           endif
        enddo
c
c.... assemble the contributions
c
        rNa(:) = -WdetJb * flux
c
c.... ------------------------->  Residual  <--------------------------
c
c.... add the flux to the residual
c
        do n = 1, nshlb
           nodlcl = lnode(n)
 
           rl(:,nodlcl) = rl(:,nodlcl) - shape(:,nodlcl) * rNa(:)
        enddo
c
c.... -------------------->  Aerodynamic Forces  <---------------------
c
        if ((ires .ne. 2) .and. (iter .eq. nitr)) then
c
c.... compute the forces on the body
c
           if(isclr.eq.0)   HFlux    = sum(flux)
c
        endif
c
c.... end of integration loop
c
        enddo

c
c.... return
c
        return
        end
  
