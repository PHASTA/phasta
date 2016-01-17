        subroutine Au1MFG (ypre,      y,    ac,     x,
     &                     rmes,      res,       uBrg,
     &                     BDiag,     iBC,       BC,        
     &                     iper,      ilwork,    shp,
     &                     shgl,      shpb,
     &                     shglb)
c
c----------------------------------------------------------------------
c
c This routine performs a matrix-vector product for the Matrix-Free
c  Implicit/Iterative solver using a one-sided scheme.
c
c input:
c  y      (nshg,ndof)           : Y-variables
c  ypre   (nshg,nflow)           : preconditioned Y-variables 
c                                  (perturbed, no-scalars)
c  x      (numnp,nsd)            : node coordinates
c  rmes   (nshg,nflow)           : modified residual
c  res    (nshg,nflow)           : residual
c  uBrg   (nshg,nflow)           : Krylov space vector
c  BDiag  (nshg,nflow,nflow)         : block-diagonal preconditioner
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  engBC  (nshg)                : energy for BC on density or pressure
c  shp(b) (nshape,ngauss)        : element shape functions (boundary)
c  shgl(b)(nsd,nshape,ngauss)    : local gradients of shape functions
c
c output:
c  uBrg   (nshg,nflow)           : Krylov space vectors
c
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension y(nshg,ndof),             ypre(nshg,nflow),
     &            x(numnp,nsd),             ac(nshg,ndof),
     &            rmes(nshg,nflow),          ytmp(nshg,nflow),
     &            res(nshg,nflow),           uBrg(nshg,nflow),
     &            BDiag(nshg,nflow,nflow),       iBC(nshg),
     &            BC(nshg,ndofBC),          iper(nshg)
c
        dimension uBtmp(nshg,nflow),         tmpBC(nshg),
     &            ilwork(nlwork)
c        

        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
        
c$$$        dimension shp(nshl,ngauss),  shgl(nsd,nshl,ngauss), 
c$$$     &            shpb(nshl,ngaussb),
c$$$     &            shglb(nsd,nshl,ngaussb)
c
c.... calculate Rmod(y + eps u_i)
c
        uBtmp = zero
c
c        call yshuffle(ypre, 'new2old ')
        uBrg = ypre + eGMRES * uBrg
c
        call i3LU (BDiag,  uBrg,  'backward')
c
        call yshuffle(uBrg, 'old2new ')
c
        call itrBC (uBrg, uBrg, iBC,  BC, iper, ilwork)
c
        call itrRes (uBrg,            y,
     &               x,               shp,
     &               shgl,            iBC,
     &               BC,              shpb,
     &               shglb,           uBtmp,
     &               iper,            ilwork, ac)
c
        call i3LU (BDiag,  uBtmp,  'forward ')
c
c.... calculate ( Rmod(y + eps u_i) - Rmod(y) ) / eps
c
        uBrg = ( uBtmp - rmes ) / eGMRES
c ... before returning lets put ypre back in the new format 
c         call yshuffle(ypre, 'old2new ')
c
c.... flop count
c
   !      flops = flops + 4*nflow*nshg
c
c.... end
c
        return
        end




