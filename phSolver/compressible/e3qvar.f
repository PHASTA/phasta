        subroutine e3qvar (ycl, shp,     shgl,    rho,
     &                     xl,  g1yi,    g2yi,    g3yi,
     &                     shg, dxidx,   WdetJ,   T, 
     &                     cp,  u1,      u2,      u3)
c
c----------------------------------------------------------------------
c
c  This routine computes the variables at integration point
c  necessary for the computation of the diffusive flux vector.
c
c input:
c  ycl    (npro,nshape,ndof)      : primitive variables
c  shp    (npro,nshape)         : element shape-functions
c  shgl   (npro,nsd,nshape)     : element local-grad-shape-functions
c  xl     (npro,nenl,nsd)       : nodal coordinates at current step
c
c output:
c  g1yi   (npro,nflow)           : grad-y in direction 1
c  g2yi   (npro,nflow)           : grad-y in direction 2
c  g3yi   (npro,nflow)           : grad-y in direction 3
c  shg    (npro,nshape,nsd)       : element global grad-shape-functions
c  dxidx  (npro,nsd,nsd)        : inverse of deformation gradient
c  WdetJ  (npro)                : weighted Jacobian
c  T      (npro)                : temperature
c  cp     (npro)                : specific heat at constant pressure
c  u1     (npro)                : x1-velocity component
c  u2     (npro)                : x2-velocity component
c  u3     (npro)                : x3-velocity component
c
c
c Zdenek Johan, Summer 1990. (Modified from e2ivar.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997. Primitive Variables
c----------------------------------------------------------------------
c
        include "common.h"
c
c  passed arrays
c
        dimension ycl(npro,nshl,ndof),  rho(npro),
     &            shp(npro,nshl),      tmp(npro),
     &            shgl(npro,nsd,nshl), xl(npro,nenl,nsd),
     &            g1yi(npro,nflow),    g2yi(npro,nflow),
     &            g3yi(npro,nflow),    shg(npro,nshl,nsd), 
     &            dxidx(npro,nsd,nsd), WdetJ(npro),
     &            T(npro),             cp(npro),                  
     &            u1(npro),            u2(npro),
     &            u3(npro),            pres(npro)
c
c  local arrays
c
        dimension dxdxi(npro,nsd,nsd), Sclr(npro)

        ttim(34) = ttim(34) - secs(0.0)
c
c.... ------------->  Primitive variables at int. point  <--------------
c
c.... compute primitive variables
c
       pres = zero
       u1   = zero
       u2   = zero
       u3   = zero
       T    = zero


c
       do n = 1, nshl
c
c  y(intp)=SUM_{a=1}^nshape (N_a(intp) Ya) (we don't need pressure here)
c   
          pres = pres + shp(:,n) * ycl(:,n,1)       
          u1   = u1   + shp(:,n) * ycl(:,n,2)
          u2   = u2   + shp(:,n) * ycl(:,n,3)
          u3   = u3   + shp(:,n) * ycl(:,n,4)
          T    = T    + shp(:,n) * ycl(:,n,5)
       enddo
c
c$$$c.... compute cp, only thermodynamic property needed for qdiff
c$$$c
c$$$       cp    = Rgas * gamma / gamma1
c.... get the thermodynamic properties
c

        if (iLSet .ne. 0)then
           Sclr = zero
           isc=abs(iRANS)+6
           do n = 1, nshl
              Sclr = Sclr + shp(:,n) * ycl(:,n,isc)
           enddo
        endif
c
        if (Navier .eq. 1) ithm = 7

        call getthm (pres,            T,                  Sclr,
     &               tmp,             rho,                tmp,
     &               tmp,             tmp,                tmp,
     &               cp,              tmp,                tmp,
     &               tmp,             tmp)
c
c.... --------------------->  Element Metrics  <-----------------------
c
        call e3metric( xl,         shgl,        dxidx,  
     &                 shg,        WdetJ)

c
c.... --------------------->  Global Gradients  <-----------------------
c
        g1yi = zero
        g2yi = zero
        g3yi = zero
c
c
        do n = 1, nshl
c.... compute the global gradient of Y-variables
c
c
c  Y_{,x_i}=SUM_{a=1}^nenl (N_{a,x_i}(intp) Ya)
c
          g1yi(:,1) = g1yi(:,1) + shg(:,n,1) * ycl(:,n,1)
          g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * ycl(:,n,2)
          g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * ycl(:,n,3)
          g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * ycl(:,n,4)
          g1yi(:,5) = g1yi(:,5) + shg(:,n,1) * ycl(:,n,5)
c
          g2yi(:,1) = g2yi(:,1) + shg(:,n,2) * ycl(:,n,1)
          g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * ycl(:,n,2)
          g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * ycl(:,n,3)
          g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * ycl(:,n,4)
          g2yi(:,5) = g2yi(:,5) + shg(:,n,2) * ycl(:,n,5)
c
          g3yi(:,1) = g3yi(:,1) + shg(:,n,3) * ycl(:,n,1)
          g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * ycl(:,n,2)
          g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * ycl(:,n,3)
          g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * ycl(:,n,4)
          g3yi(:,5) = g3yi(:,5) + shg(:,n,3) * ycl(:,n,5)

        enddo

        ttim(34) = ttim(34) + secs(0.0)

c
c.... return
c

       return
       end
