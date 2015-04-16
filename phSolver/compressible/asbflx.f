        subroutine AsBFlx (y,       x,       shpb,    shglb,
     &                     ienb,    materb,  iBCB,    BCB,
     &                     invflx,  flxres,  flxLHS,  flxnrm)
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
        dimension y(nshg,ndof),             x(numnp,nsd),
     &            shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),
     &            ienb(npro,nshl),          materb(npro),
     &            iBCB(npro,ndiBCB),        BCB(npro,nshlb,ndBCB),
     &            invflx(numnp),            flxres(numnp,nflow),
     &            flxLHS(numnp,1),          flxnrm(numnp,nsd)
c
        dimension ycl(npro,nshl,ndof),       xlb(npro,nenl,nsd),
     &            rl(npro,nshl,nflow),      rml(npro,nshl,nflow),
     &            flhsl(npro,nshl,1),       fnrml(npro,nshl,nsd),
     &            lnflx(npro),              lnode(27)
c        
        dimension sgn(npro,nshl)
c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c

c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c
c.... ------------------------>  Residual  <---------------------------
c
c.... gather the variables
c
        call localy(y,      ycl, ienb,   ndofl,  'gather  ')
        call localx(x,      xlb,ienb,   nsd,    'gather  ')
c
c
c.... get the boundary element residual
c
        rl = zero
        call e3b (ycl, ycl,     iBCB,    BCB,     shpb,    shglb,
     &            xlb,     rl,      rml,      sgn)
c
c.... assemble the residual
c 
        call local (flxres, rl, ienb,   nflow,   'scatter ')
c
c.... --------------------->  LHS and Normal  <------------------------
c
c.... compute the boundary LHS and normal
c
        flhsl = zero
        fnrml = zero
c
c.... 2D
c
c       if (nsd .eq. 2) then
c         call f2lhs (shpb,   shglb,  xlb,    flhsl,
c    &                fnrml)
c       endif
c
c.... 3D
c
c       if (nsd .eq. 3) then
          call f3lhs (shpb,   shglb,  xlb,    flhsl,
     &                fnrml, sgn)
c       endif
c
c.... reset the non-contributing element values
c
        lnflx = 0
        do n = 1, nenbl
          lnflx = lnflx + 
     &          min(1, invflx(abs(ienb(:,mnodeb(n,lelCat,nsd)))))
        enddo
c
        do n = 1, nshl
          where (lnflx .ne. nenbl)   flhsl(:,n,1) = zero
c
          do i = 1, nsd
            where (lnflx .ne. nenbl) fnrml(:,n,i) = zero
          enddo
        enddo
c
c.... assemble the boundary LHS and normal
c
        call local (flxLHS, flhsl,  ienb,   1,      'scatter ')
c
        call local (flxnrm, fnrml,  ienb,   nsd,    'scatter ')
c
c.... end
c
        return
        end
