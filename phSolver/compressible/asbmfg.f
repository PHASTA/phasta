        subroutine AsBMFG (y,       x,       shpb,    shglb,
     &                     ienb,    materb,  iBCB,    BCB,
     &                     res,     rmes)
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
        dimension y(nshg,ndofl),           x(numnp,nsd),
     &            shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),        
     &            ienb(npro,nshl),          materb(npro),
     &            iBCB(npro,ndiBCB),        BCB(npro,nshlb,ndBCB),
     &            res(nshg,nflow),         rmes(nshg,nflow)
c
        dimension ycl(npro,nshl,ndofl),  xlb(npro,nenl,nsd),
     &            rl(npro,nshl,nflow),
     &            rml(npro,nshl,nflow)
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

        call localy(y,      ycl,     ienb,   ndofl,  'gather  ')
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
c

c.... get the boundary element residuals
c
        rl  = zero
        rml = zero
c
!  pass the memory location of ycl to both yl and ycl in e3b.  This may
 !  seem dangerous since yl in e3b is :,nflow and ycl is :,ndof but they
 !  do not write to yl (out of bounds at least), only use the data there 
 !  so both will access data
 !  properly from this location.
c
        call e3b  (ycl,     ycl,     iBCB,    BCB,     shpb,    shglb,
     &             xlb,     rl,      rml,     sgn)
c
c.... assemble the residual and the modified residual
c
        call local(res,    rl,     ienb,   nflow,  'scatter ')


c
        if (Navier .eq. 1 .and. ires.ne.1 )
     &  call local(rmes,   rml,    ienb,   nflow,  'scatter ')
c
c.... end
c
        return
        end
c
c
c
        subroutine AsBMFGSclr (y,       x,       shpb,    shglb,
     &                         ienb,    materb,  iBCB, 
     &                         BCB,     rest,    rmest)
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
        dimension y(nshg,ndofl),           x(numnp,nsd),
     &            shpb(nshl,maxsh),         
     &            shglb(nsd,nshl,maxsh),         
     &            ienb(npro,nshl),       materb(npro),
     &            iBCB(npro,ndiBCB),   BCB(npro,nshlb,ndBCB),
     &            rest(nshg),         rmest(nshg)
c
        dimension ycl(npro,nshl,ndofl),   xlb(npro,nenl,nsd),
     &            rtl(npro,nshl),      
     &            rmtl(npro,nshl)
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
        call localy (y,      ycl,     ienb,   ndofl,  'gather  ')
        call localx (x,      xlb,    ienb,   nsd,    'gather  ')
c
c.... get the boundary element residuals
c
        rtl  = zero
        rmtl = zero
c
c.... 3D
c
            call e3bSclr (ycl,    iBCB,    BCB,     
     &                    shpb,  shglb,   sgn,
     &                    xlb,   rtl,     rmtl)
c
c.... assemble the residual and the modified residual
c

        call local (rest,    rtl,     ienb,   1,  'scatter ')


c
        if (Navier .eq. 1)
     &  call local (rmest,   rmtl,    ienb,   1,  'scatter ')
c
c.... end
c
        return
        end


