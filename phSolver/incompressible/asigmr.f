        subroutine AsIGMR (y,       ac,      x,       xmudmi,
     &                     shp,     shgl,    ien,     
     &                     res,     qres,
     &                     xKebe,   xGoC,    rerr, CFLworst)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use stats
      use rlssave  ! Use the resolved Leonard stresses at the nodes.
      use timedata    ! time series
      use turbsa                ! access to d2wall


      include "common.h"
c
        dimension y(nshg,ndofl),              ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(nshl,ngauss),            shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),
     &            res(nshg,nflow),
     &            qres(nshg,idflx)

c
        dimension yl(npro,nshl,ndofl),         acl(npro,nshl,ndofl),
     &            xl(npro,nenl,nsd),           dwl(npro,nenl),      
     &            rl(npro,nshl,nflow), 
     &            ql(npro,nshl,idflx)
c        
        dimension xKebe(npro,9,nshl,nshl), 
     &            xGoC(npro,4,nshl,nshl)
c
        dimension rlsl(npro,nshl,6) 

c
        real*8    lStsVec(npro,nshl,nResDims)
        
        dimension xmudmi(npro,ngauss)
        dimension sgn(npro,nshl)
        dimension CFLworst(npro)
c
        real*8 rerrl(npro,nshl,6), rerr(nshg,10)
c
c.... gather the variables
c
c
c.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        call local (qres,   ql,     ien,    idflx,  'gather  ')
        if (iRANS .eq. -2) then ! kay-epsilon
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
 
        if( (iLES.gt.10).and.(iLES.lt.20)) then  ! bardina 
           call local (rls, rlsl,     ien,       6, 'gather  ')  
        else
           rlsl = zero
        endif      

c
c.... zero the matrices if they are being recalculated
c
        if (lhs. eq. 1)  then
           xKebe = zero
           xGoC  = zero
        endif   
c
c.... get the element residuals, LHS matrix, and preconditioner
c
        rl     = zero

        if(ierrcalc.eq.1) rerrl = zero

        call e3  (yl,      acl,     dwl,     shp,
     &            shgl,    xl,      rl,      
     &            ql,      xKebe,   xGoC,    xmudmi, 
     &            sgn,     rerrl,  rlsl,     CFLworst)
c
c.... assemble the statistics residual
c
        if ( stsResFlg .eq. 1 ) then
           call e3StsRes ( xl, rl, lStsVec )
           call local( stsVec, lStsVec, ien, nResDims, 'scatter ')
        else
c
c.... assemble the residual
c
           call local (res,    rl,     ien,    nflow,  'scatter ')
           
           if ( ierrcalc .eq. 1 ) then
              call local (rerr, rerrl,  ien, 6, 'scatter ')
           endif
        endif
c
c.... end
c
        if (exts.and.ires.ne.2) then
           if ((iter.eq.1).and.(mod(lstep,freq).eq.0)) then
              call timeseries(yl,xl,ien,sgn)
           endif
        endif
        
        return
        end



C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c-----------------------------------------------------------------------
c=======================================================================


        subroutine AsIGMRSclr(y,       ac,      x,       
     &                     shp,     shgl,    ien,     
     &                     res,     qres,    xSebe, xmudmi )
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use     turbSA  
      include "common.h"
c
        dimension y(nshg,ndofl),              ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(nshl,ngauss),            shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),
     &            res(nshg),                  qres(nshg,nsd)

c
        real*8    yl(npro,nshl,ndofl),        acl(npro,nshl,ndofl),
     &            xl(npro,nenl,nsd),         
     &            rl(npro,nshl),              ql(npro,nshl,nsd),
     &            dwl(npro,nenl)            
c        
        real*8    xSebe(npro,nshl,nshl),      xmudmi(npro,ngauss) 
c
c.... gather the variables
c
        real*8 sgn(npro,nshl)
c
c.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        if(iRANS.lt. 0) 
     &  call localx(d2wall, dwl,    ien,    1,      'gather  ')
        call local (qres,   ql,     ien,    nsd,    'gather  ')
c
c.... zero the matrices if they are being recalculated
c
        if (lhs. eq. 1)  then
           xSebe = zero
        endif   
c
c.... get the element residuals, LHS matrix, and preconditioner
c
      rl = zero
      call e3Sclr  (yl,      acl,     shp,
     &              shgl,    xl,      dwl,
     &              rl,      ql,      xSebe,   
     &              sgn, xmudmi)
c
c.... assemble the residual
c
        call local (res,    rl,     ien,    1,  'scatter ')
c
c.... end
c
        return
        end
