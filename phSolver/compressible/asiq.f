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
c     y     (nshg,ndof)        : Y variables
c     x     (numnp,nsd)         : nodal coordinates
c     shp   (nshape,ngauss)     : element shape-functions
c     shgl  (nsd,nshape,ngauss) : element local shape-function gradients
c     ien   (npro)              : nodal connectivity array
c
c output:
c     qres  (nshg,nflow-1,nsd)  : residual vector for diffusive flux
c     rmass  (nshg)            : lumped mass matrix
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,ndof),               x(numnp,nsd),              
     &            shp(nshl,ngauss),  
     &            shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),
     &            qres(nshg,idflx),    rmass(nshg)
c
c.... element level declarations
c
        dimension ycl(npro,nshl,ndof),        xl(npro,nenl,nsd),         
     &            ql(npro,nshl,idflx),       rmassl(npro,nshl),
     &            xmudmi(npro,ngauss)
c
        dimension sgn(npro,nshape)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
c
c.... gather the variables
c

        call localy(y,      ycl,     ien,    ndof,   'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
c
c.... get the element residuals 
c
        ql     = zero
        rmassl = zero

        call e3q  (ycl,      shp,    shgl,    
     &             xl,       ql,     rmassl, 
     &             xmudmi,   sgn )

c
c.... assemble the diffusive flux residual 
c
        call local (qres,   ql,  ien,  idflx,'scatter ')
        call local (rmass,  rmassl, ien,  1,      'scatter ')
c
c.... end
c
        return
        end
