        subroutine AsIMFG (y,       ac,      x,     xmudmi,   shp,
     &                     shgl,    ien,     mater,
     &                     res,     rmes,    BDiag,   qres, rerr)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use rlssave   ! Use the resolved Leonard stresses at the nodes.

      include "common.h"
c
        dimension y(nshg,ndofl),            ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(nshl,MAXQPT),  
     &            shgl(nsd,nshl,MAXQPT),
     &            ien(npro,nshl),
     &            mater(npro),               res(nshg,nflow),
     &            rmes(nshg,nflow),         BDiag(nshg,nflow,nflow),
     &            qres(nshg,idflx)

c
        dimension ycl(npro,nshl,ndofl),       acl(npro,nshl,ndofl),
     &            xl(npro,nenl,nsd),         
     &            rl(npro,nshl,nflow),       rml(npro,nshl,nflow),
     &            BDiagl(npro,nshl,nflow,nflow),
     &            ql(npro,nshl,idflx)
c        
        dimension rlsl(npro,nshl,6)          
        dimension  xmudmi(npro,ngauss)
        dimension sgn(npro,nshl)
c
        real*8 rerrl(npro,nshl,6), rerr(nshg,10)
c
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
        call localy(y,      ycl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        
        if (idiff >= 1 .or. isurf .eq. 1)
     &    call local (qres,   ql,  ien, idflx, 'gather  ')

        if( (iLES.gt.10).and.(iLES.lt.20)) then  ! bardina 
           call local (rls, rlsl,     ien,       6, 'gather  ')  
        else
           rlsl = zero
        endif  
c
c.... get the element residuals and preconditioner
c
        rl     = zero
        rml    = zero
        BDiagl = zero
        EGmassd= one  ! just a dummy real since we don't have a LHS with MFI
        if(ierrcalc.eq.1) rerrl = zero
        
        ttim(31) = ttim(31) - secs(0.0)
!  pass the memory location of ycl to both yl and ycl in e3b.  This may
!  seem dangerous since yl in e3b is :,nflow and ycl is :,ndof but they
!  do not write to yl (at least not out of bounds), only use the data
!  there so both will access data
!  properly from this location.

            call e3  (ycl,     ycl,     acl,     shp,
     &                shgl,    xl,      rl,      rml, xmudmi,
     &                BDiagl,  ql,      sgn,     rlsl, EGmassd,
     &                rerrl)

        ttim(31) = ttim(31) + secs(0.0)
c
c.... assemble the residual and the modified residual
c

        call local (res,    rl,     ien,    nflow,  'scatter ')
        call local (rmes,   rml,    ien,    nflow,  'scatter ')
c
c       res is G_A obtained using local  A_{e=1}^n_e G^e_a
c
           if ( ierrcalc .eq. 1 ) then
              call local (rerr, rerrl,  ien, 6, 'scatter ')
           endif
c
c.... assemble the Block-Diagonal
c
        if (iprec .ne. 0)
     &     call local (BDiag,  BDiagl, ien, nflow*nflow, 'scatter ')
        
c
c.... end
c
        return
        end





