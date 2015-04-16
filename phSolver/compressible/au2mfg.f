        subroutine Au2MFG (ypre,      y,    ac,     x,
     &                     rmes,      res,       Dy,
     &                     uBrg,      BDiag,     iBC,       BC,
     &                     iper,      ilwork,
     &                     shp,       shgl,      
     &                     shpb,      shglb)
c
c----------------------------------------------------------------------
c
c This routine performs a matrix-vector product for the Matrix-Free
c  Implicit/Iterative solver using a two-sided scheme and subtracts
c  it from the residual.
c
c input:
c  y      (nshg,ndof)           : Y-variables
c  ypre   (nshg,ndof)           : preconditioned Y-variables
c  x      (numnp,nsd)            : node coordinates
c  rmes   (nshg,nflow)           : modified residual
c  res    (nshg,nflow)           : residual
c  Dy     (nshg,nflow)           : solution step
c  BDiag  (nshg,nflow,nflow)         : block-diagonal preconditioner
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c
c output:
c  uBrg   (nshg,nflow)           : Krylov vector ( R - A x )
c
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,ndof),               ypre(nshg,nflow),
     &            x(numnp,nsd),               ac(nshg,ndof),
     &            rmes(nshg,nflow),            ytmp(nshg,nflow),
     &            res(nshg,nflow),             Dy(nshg,nflow),
     &            uBrg(nshg,nflow),            BDiag(nshg,nflow,nflow),
     &            iBC(nshg),                  BC(nshg,ndofBC),
     &            iper(nshg),                 Dy2(nshg,nflow)
c
        dimension uBtmp1(nshg,nflow),          uBtmp2(nshg,nflow),
     &            tmpBC(nshg),                ilwork(nlwork)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 

c$$$        dimension shp(nshl,ngauss),    shgl(nsd,nshl,ngauss), 
c$$$     &            shpb(nshl,ngaussb), shglb(nsd,nshl,ngaussb)
c
c.... compute the finite difference interval
c
        Dy2 = Dy**2
	call sumgat (Dy2, nflow, summed, ilwork)
        eps = epsM**pt66 / sqrt(summed)
c
c.... calculate R(y + eps x)
c
        uBtmp1 = zero
c
c        call yshuffle(ypre, 'new2old ')
c
        uBrg = ypre + eps * Dy
c
        call i3LU (BDiag,  uBrg,  'backward')
c
        call yshuffle(uBrg, 'old2new ')
c
        call itrBC (uBrg, uBrg, iBC,  BC,   iper, ilwork)
c
        call itrRes (uBrg,            y,
     &               x,               shp,
     &               shgl,            iBC,
     &               BC,              shpb,
     &               shglb,           uBtmp1,
     &               iper,            ilwork)
c
        call i3LU (BDiag,  uBtmp1,  'forward ')
c
c.... calculate R(y - eps x)
c
        uBtmp2 = zero
c
c        call yshuffle(ypre, 'new2old ')
c
        uBrg = ypre - eps * Dy
c
        call i3LU (BDiag,  uBrg,  'backward')
c
        call yshuffle(uBrg, 'old2new ')

        call itrBC (uBrg, uBrg, iBC,  BC,   iper, ilwork)
c
        call itrRes (uBrg,            y,
     &               x,               shp,
     &               shgl,            iBC,
     &               BC,              shpb,
     &               shglb,           uBtmp2,
     &               iper,            ilwork, ac)
c
        call i3LU (BDiag,  uBtmp2,  'forward ')
c
c.... compute  R(y) - ( R(y + eps x) - R(y - eps x) ) / 2 eps
c
        uBrg = res - ( uBtmp1 - uBtmp2 ) / (two * eps)
c
c.... flop count
c
        flops = flops + 8*nflow*nshg+nshg
c
c.... end
c
        return
        end
