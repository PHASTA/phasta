      subroutine AsBFlx (u,           y,           ac,      
     &                   x,           shpb,    
     &                   shglb,       ienb,        iBCB,    
     &                   BCB,         invflx,      flxres,
     &                   flxLHS,      flxnrm,      xKebe )
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  boundary elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use turbSA ! access to d2wall
        include "common.h"
c
        dimension y(nshg,ndofl),           x(numnp,nsd),
     &            ac(nshg,ndofl),          u(nshg,nsd),
     &            shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),         
     &            ienb(npro,nshl),         
     &            iBCB(npro,ndiBCB),       BCB(npro,nshlb,ndBCB),
     &            invflx(nshg),            flxres(nshg,nflow),
     &            flxLHS(nshg,1),          flxnrm(nshg,nsd)
c
        dimension yl(npro,nshl,ndofl),     xlb(npro,nenl,nsd),
     &            rl(npro,nshl,nflow),     sgn(npro,nshl),
     &            flhsl(npro,nshl,1),      fnrml(npro,nshl,nsd),
     &            lnflx(npro),             lnode(27),
     &            ul(npro,nshl,nsd),       acl(npro,nshl,ndofl)
        real*8 dwl(npro,nshl)
        
        dimension xKebe(npro,9,nshl,nshl) 

c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c     
c.... gather the variables
c
        call localy(y,      yl,     ienb,   ndofl,  'gather  ')
        call localy(ac,     acl,    ienb,   ndofl,  'gather  ')
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
        call localx(u,      ul,     ienb,   nsd,    'gather  ')
        if(iRANS.eq.-2) then
           call local(d2wall, dwl, ienb, 1, 'gather  ')
        endif

        rl    = zero
        flhsl = zero
        fnrml = zero
c
        ires = 2
        call e3b  (ul,      yl,      acl,     iBCB,    BCB,     
     &             shpb,    shglb,
     &             xlb,     rl,      sgn,     dwl,     xKebe)
        ires = 1
c
c.... assemble the residuals
c
        call local (flxres, rl,     ienb,   nflow,  'scatter ')
c
c.... compute the LHS for the flux computation (should only be done
c     once)
c
        call f3lhs (shpb,       shglb,      xlb,
     &              flhsl,      fnrml,      sgn )

c     
c.... reset the non-contributing element values
c
        lnflx = 0
        do n = 1, nshlb
          lnflx = lnflx + min(1, invflx(ienb(:,lnode(n))))
        enddo
c
        do n = 1, nshl
          where (lnflx .ne. nshlb)   flhsl(:,n,1) = zero
          do i = 1, nsd
            where (lnflx .ne. nshlb) fnrml(:,n,i) = zero
          enddo
        enddo
c
c.... assemble the boundary LHS and normal
c
        call local (flxLHS, flhsl,  ienb,   1,      'scatter ')
        call local (flxnrm, fnrml,  ienb,   nsd,    'scatter ')
c     
c.... end
c
        return
        end
