        subroutine AsIFlx (y,       ac,      x,       xmudmi,   shp, 
     &                     shgl,    ien,     mater,   flxres)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use rlssave
        include "common.h"
c
        dimension y(nshg,ndof),            ac(nshg,ndof),
     &            x(numnp,nsd),             
     &            shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            ien(npro,nshl),
     &            mater(npro),            flxres(numnp,nflow)
c
        dimension ycl(npro,nshl,ndof),   acl(npro,nshl,ndof),
     &            xl(npro,nenl,nsd),        
     &            rl(npro,nshl,nflow),   BDiagl(npro,nshl,nflow,nflow)
c
        dimension ql(npro,nshl,(nflow-1)*nsd)
c        
        dimension  xmudmi(npro,ngauss)        
        dimension sgn(npro,nshl)
        dimension EGmass(npro,nedof,nedof), rml(npro,nshl,nflow)

        dimension rlsl(npro,nshl,6) 
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
        call localy(y,      ycl,    ien,    ndof,   'gather  ')
        call localy(ac,    acl,    ien,    ndof,   'gather  ')
        call localx(x,      xl,    ien,    nsd,    'gather  ')

        if( (iLES.gt.10).and.(iLES.lt.20)) then  ! bardina 
           call local (rls, rlsl,     ien,       6, 'gather  ')  
        else
           rlsl = zero
        endif      
c
c
c.... get the element residual
c
        rl = zero

        ql = zero !only used now in elmmfg
        
        EGmassd= one  ! just a dummy real since we don't have a LHS with MFI
        
        if (nsd .eq. 3) then
          call e3 (ycl,      ycl,    acl,   shp(lcsyst,1:nshl,:),
     &             shgl(lcsyst,:,1:nshl,:), xl,      rl,  rml,
     &             xmudmi,   BDiagl,  ql,      sgn, rlsl, EGmassd)
        endif
c
c.... assemble the residual
c
        call local (flxres, rl,   ien,    nflow,   'scatter ')
c
c.... end
c
        return
        end
