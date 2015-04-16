        subroutine AsIqGradV (y,       x,       shp,
     &                   shgl,    ien,
     &                   qres )
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
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,ndof),               x(numnp,nsd),            
     &            shp(nshl,ngauss),         shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),  
     &            qres(nshg,nsdsq)
c
        dimension yl(npro,nshl,ndof),          xl(npro,nenl,nsd),
     &            ql(npro,nshl,nsdsq)
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
c
c.... get the element residuals 
c
        ql     = zero

        call e3qGradV  (yl,         shp,      shgl,    
     &             xl,         ql,
     &             sgn  )

c
c.... assemble the diffusive flux residual 
c
        call local (qres,   ql,  ien,  nsdsq,  'scatter ')
c
c.... end
c
        return
        end

