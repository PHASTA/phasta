        subroutine Au2MFG (ypre,      y,    ac,     x,
     &                     rmes,      res,       Dy,
     &                     uBrg,      BDiag,     iBC,       BC,
     &                     iper,      ilwork,
     &                     shp,       shgl,      
     &                     shpb,      shglb)
!
!----------------------------------------------------------------------
!
! This routine performs a matrix-vector product for the Matrix-Free
!  Implicit/Iterative solver using a two-sided scheme and subtracts
!  it from the residual.
!
! input:
!  y      (nshg,ndof)           : Y-variables
!  ypre   (nshg,ndof)           : preconditioned Y-variables
!  x      (numnp,nsd)            : node coordinates
!  rmes   (nshg,nflow)           : modified residual
!  res    (nshg,nflow)           : residual
!  Dy     (nshg,nflow)           : solution step
!  BDiag  (nshg,nflow,nflow)         : block-diagonal preconditioner
!  iBC    (nshg)                : BC codes
!  BC     (nshg,ndofBC)         : BC constraint parameters
!
! output:
!  uBrg   (nshg,nflow)           : Krylov vector ( R - A x )
!
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
        include "common.h"
!
        dimension y(nshg,ndof),               ypre(nshg,nflow),
     &            x(numnp,nsd),               ac(nshg,ndof),
     &            rmes(nshg,nflow),            ytmp(nshg,nflow),
     &            res(nshg,nflow),             Dy(nshg,nflow),
     &            uBrg(nshg,nflow),            BDiag(nshg,nflow,nflow),
     &            iBC(nshg),                  BC(nshg,ndofBC),
     &            iper(nshg),                 Dy2(nshg,nflow)
!
        dimension uBtmp1(nshg,nflow),          uBtmp2(nshg,nflow),
     &            tmpBC(nshg),                ilwork(nlwork)
!
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 

!$$$        dimension shp(nshl,ngauss),    shgl(nsd,nshl,ngauss), 
!$$$     &            shpb(nshl,ngaussb), shglb(nsd,nshl,ngaussb)
!
!.... compute the finite difference interval
!
        Dy2 = Dy**2
        call sumgat (Dy2, nflow, summed, ilwork)
        eps = epsM**pt66 / sqrt(summed)
!
!.... calculate R(y + eps x)
!
        uBtmp1 = zero
!
!        call yshuffle(ypre, 'new2old ')
!
        uBrg = ypre + eps * Dy
!
        call i3LU (BDiag,  uBrg,  'backward')
!
        call yshuffle(uBrg, 'old2new ')
!
        call itrBC (uBrg, uBrg, iBC,  BC,   iper, ilwork)
!
        call itrRes (uBrg,            y,
     &               x,               shp,
     &               shgl,            iBC,
     &               BC,              shpb,
     &               shglb,           uBtmp1,
     &               iper,            ilwork, ac)
!Added ac to the end if itrRes, but not tested - Nicholas
!
        call i3LU (BDiag,  uBtmp1,  'forward ')
!
!.... calculate R(y - eps x)
!
        uBtmp2 = zero
!
!        call yshuffle(ypre, 'new2old ')
!
        uBrg = ypre - eps * Dy
!
        call i3LU (BDiag,  uBrg,  'backward')
!
        call yshuffle(uBrg, 'old2new ')

        call itrBC (uBrg, uBrg, iBC,  BC,   iper, ilwork)
!
        call itrRes (uBrg,            y,
     &               x,               shp,
     &               shgl,            iBC,
     &               BC,              shpb,
     &               shglb,           uBtmp2,
     &               iper,            ilwork, ac)
!
        call i3LU (BDiag,  uBtmp2,  'forward ')
!
!.... compute  R(y) - ( R(y + eps x) - R(y - eps x) ) / 2 eps
!
        uBrg = res - ( uBtmp1 - uBtmp2 ) / (two * eps)
!
!.... flop count
!
   !      flops = flops + 8*nflow*nshg+nshg
!
!.... end
!
        return
        end
