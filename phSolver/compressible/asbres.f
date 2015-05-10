        subroutine AsBRes (y, yc,      x,       shpb,    shglb,
     &                     ienb,    materb,  iBCB,    BCB,
     &                     rmes)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  boundary elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,nflow),           x(numnp,nsd),
     &            yc(nshg,ndof),           shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),
     &            ienb(npro,nshl),        materb(npro),
     &            iBCB(npro,ndiBCB),     BCB(npro,nshlb,ndBCB),
     &            rmes(nshg,nflow)
c
        dimension yl(npro,nshl,nflow),    xlb(npro,nenl,nsd),
     &            ycl(npro,nshl,ndof),    rml(npro,nshl,nflow)
c
        dimension sgn(npro,nshl)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c     
c.... gather the variables
c
        call localy(y,      yl,     ienb,   nflow,  'gather  ')
        call localy(yc,     ycl,    ienb,   ndof,  'gather  ')
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
c
c.... get the boundary element residuals
c
        rml = zero
        call e3b  (yl,      ycl,     iBCB,    BCB,     shpb,    shglb,
     &             xlb,     rml,     rml,     sgn)

c
c.... assemble the residual and the modified residual
c
        if (iabres .eq. 1) rml = abs(rml)
c
        call local(rmes,   rml,    ienb,   nflow,  'scatter ')
c
c.... end
c
        return
        end




