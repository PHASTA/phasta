        subroutine AsIq (y,       x,       shp,
     &                   shgl,    ien,     xmudmi,
     &                   qres,    rmass    )
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c interior elements for the global reconstruction of the diffusive
c flux vector.
c
c input:
c     y     (numnp,ndof)        : Y variables
c     x     (numnp,nsd)         : nodal coordinates
c     shp   (nen,nintg)         : element shape-functions
c     shgl  (nsd,nen,nintg)     : element local shape-function gradients
c     ien   (npro)              : nodal connectivity array
c
c output:
c     qres  (numnp,nsd,nsd)  : residual vector for diffusive flux
c     rmass  (numnp)            : lumped mass matrix
c
c----------------------------------------------------------------------
c
        use turbsa      ! access to d2wall
        include "common.h"
c
        dimension y(nshg,ndof),               x(numnp,nsd),            
     &            shp(nshl,ngauss),         shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),      dwl(npro,nenl),
     &            qres(nshg,idflx),    rmass(nshg)
c
        dimension yl(npro,nshl,ndof),          xl(npro,nenl,nsd),
     &            ql(npro,nshl,idflx),  rmassl(npro,nshl),
     &            xmudmi(npro,ngauss)
c
        dimension sgn(npro,nshl)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        do i=1,nshl
           where ( ien(:,i) < 0 )
              sgn(:,i) = -one
           elsewhere
              sgn(:,i) = one
           endwhere
        enddo

c
c.... gather the variables
c

        call localy(y,      yl,     ien,    ndof,   'gather  ')
        call localx (x,      xl,     ien,    nsd,    'gather  ')
        if (iRANS .eq. -2) then ! kay-epsilon
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
c
c.... get the element residuals 
c
        ql     = zero
        rmassl = zero

        call e3q  (yl,         dwl,      shp,      shgl,    
     &             xl,         ql,       rmassl,
     &             xmudmi,     sgn  )

c
c.... assemble the diffusive flux residual 
c
        call local (qres,   ql,  ien,  idflx,  'scatter ')
        call local (rmass,  rmassl, ien,  1,          'scatter ')
c
c.... end
c
        return
        end


c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c interior elements for the global reconstruction of the diffusive
c flux vector.
c
c----------------------------------------------------------------------
        subroutine AsIqSclr (y,       x,       shp,
     &                       shgl,    ien,     qres,    
     &                       rmass    )
c
        use turbsa      ! access to d2wall
        include "common.h"
c
        dimension y(nshg,ndof),             x(numnp,nsd),            
     &            shp(nshl,ngauss),         shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),      dwl(npro,nenl),
     &            qres(nshg,nsd),           rmass(nshg)
c
        dimension yl(npro,nshl,ndof),       xl(npro,nenl,nsd),         
     &            ql(npro,nshl,nsd),        rmassl(npro,nshl)
c
        dimension sgn(npro,nshl)

        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
c
c.... gather the variables
c
        call localy(y,      yl,     ien,    ndof,   'gather  ')
        call localx (x,      xl,     ien,    nsd,    'gather  ')
        if (iRANS .eq. -2) then ! kay-epsilon
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
c
c.... get the element residuals 
c
        ql     = zero
        rmassl = zero

        call e3qSclr  (yl,      dwl,    shp,    shgl,    
     &                 xl,      ql,     rmassl, 
     &                 sgn             )

c
c.... assemble the temperature diffusive flux residual 
c
        call local (qres,   ql,  ien,  nsd,  'scatter ')
        call local (rmass,  rmassl, ien,  1, 'scatter ')
c
c.... end
c
        return
        end

