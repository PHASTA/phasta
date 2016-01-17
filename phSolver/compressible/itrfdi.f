        subroutine itrFDI (ypre,      y,    ac,     x,
     &                     rmes,      uBrg,      BDiag,
     &                     iBC,       BC,        iper,
     &                     ilwork,    shp,       shgl,
     &                     shpb,      shglb)
c
c----------------------------------------------------------------------
c
c  This subroutine computes the "optimum" finite difference
c  interval eps for the forward difference scheme
c
c                 Rmod(y + eps u) - Rmod(y)
c                ---------------------------
c                            eps
c
c  where  u is the step and Rmod is the modified residual.
c
c Note: A good theoretical reference is 'Practical Optimization' by
c        P.E. Gill, W. Murray and M.H. Wright  [1981].
c
c input:
c  y      (nshg,ndof)           : Y-variables
c  ypre   (nshg,ndof)           : preconditioned Y-variables
c  x      (numnp,nsd)            : node coordinates
c  rmes   (nshg,nflow)           : modified residual
c  uBrg   (nshg,nflow)           : step
c  BDiag  (nshg,nflow,nflow)         : block-diagonal preconditioner
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  
c
c
c Zdenek Johan, Winter 1989.
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,ndof),               ypre(nshg,nflow),
     &            x(numnp,nsd),               ac(nshg,ndof),
     &            rmes(nshg,nflow),
     &            uBrg(nshg,nflow),           BDiag(nshg,nflow,nflow),
     &            iBC(nshg),                  BC(nshg,ndofBC),
     &            ilwork(nlwork),
     &            iper(nshg)
c
        dimension ytmp(nshg,nflow),            rtmp(nshg,nflow) 
        dimension tmpy(nshg,ndof)
c        
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)        
c
c.... compute the accuracy (cancellation error)  ->  epsA
c
        rtmp = zero
c
        ytmp = ypre
c
c        call yshuffle(ytmp, 'new2old ')
c
        call i3LU (BDiag, ytmp, 'backward')
c
        call yshuffle(ytmp, 'old2new ')
c
        iabres = 1
c
        call itrRes (ytmp,            y,             
     &               x,               shp,
     &               shgl,            iBC,
     &               BC,              shpb,
     &               shglb,           rtmp,
     &               iper,            ilwork,
     &               ac)
c
        iabres = 0
c
        call i3LU (BDiag, rtmp, 'forward ')
c
        rtmp = rtmp**2
        call sumgat (rtmp, nflow, summed, ilwork)
        epsA = (epsM**2) * sqrt(summed)
c
c.... compute the norm of the second derivative (truncation error)
c
c.... set interval
c
        epsSD = sqrt(epsM)
c
c.... compute the first residual
c
        rtmp = zero
c
c        call yshuffle(ypre, 'new2old ')
c
        ytmp = ypre + epsSD * uBrg
c
        call i3LU (BDiag, ytmp, 'backward')
c
        call yshuffle(ytmp, 'old2new ')
c
        call itrRes (ytmp,            y,
     &               x,               shp,
     &               shgl,            iBC,
     &               BC,              shpb,
     &               shglb,           rtmp,
     &               iper,            ilwork, ac)
!Added ac to the end if itrRes, but not tested - Nicholas
c
c.... compute the second residual and add it to the first one
c
         ytmp = ypre - epsSD * uBrg
c
c         call yshuffle(ypre, 'old2new ')
c
        call i3LU (BDiag, ytmp, 'backward')
c
        call yshuffle(ytmp, 'old2new ')
        call itrRes (ytmp,            y,
     &               x,               shp,
     &               shgl,            iBC,
     &               BC,              shpb,
     &               shglb,           rtmp,
     &               iper,            ilwork, ac)
!Added ac to the end if itrRes, but not tested - Nicholas
c
        call i3LU (BDiag, rtmp, 'forward ')
c
c.... compute the second derivative and its norm
c
        rtmp = (( rtmp - two * rmes ) / epsM)**2
c
        call sumgat (rtmp, nflow, summed, ilwork)
        SDnrm = sqrt(summed)
c
c.... compute the 'optimum' interval
c
        eGMRES = two * sqrt( epsA / SDnrm )
c
c.... flop count
c
   !      flops = flops + 10*nflow*nshg+3*nshg
c
c.... end
c
        return
        end


        subroutine itrFDISclr (y,         ypre,      x,
     &                         rmes,      uBrg,      BDiag,
     &                         iBC,       BC,        engBC,     iper,
     &                         ilwork)
c
c----------------------------------------------------------------------
c
c  This subroutine computes the "optimum" finite difference
c  interval eps for the forward difference scheme
c
c                 Rmod(y + eps u) - Rmod(y)
c                ---------------------------
c                            eps
c
c  where  u is the step and Rmod is the modified residual.
c
c Note: A good theoretical reference is 'Practical Optimization' by
c        P.E. Gill, W. Murray and M.H. Wright  [1981].
c
c input:
c  y      (nshg,ndof)           : Y-variables
c  ypre   (nshg,ndof)           : preconditioned Y-variables
c  x      (numnp,nsd)            : node coordinates
c  rmes   (nshg,nflow)           : modified residual
c  uBrg   (nshg,nflow)           : step
c  BDiag  (nshg,nflow,nflow)         : block-diagonal preconditioner
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  engBC  (nshg)                : energy for BC on density or pressure
c
c
c Zdenek Johan, Winter 1989.
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,ndof),               ypre(nshg,ndof),
     &            x(numnp,nsd),
     &            rmes(nshg,nflow),
     &            uBrg(nshg,nflow),            BDiag(nshg,nflow,nflow),
     &            iBC(nshg),                  BC(nshg,ndofBC),
     &            engBC(nshg),                ilwork(nlwork),
     &            iper(nshg)
c
        dimension ytmp(nshg,ndof),            rtmp(nshg,nflow)
c
c.... compute the accuracy (cancellation error)  ->  epsA
c
        rtmp = zero
c
        ytmp = ypre
c
c       call tnanq(ytmp,ndof,"ytmp       ")
        iabres = 1
c
        call itrRes (ytmp,            y,             
     &               x,               a(mshp),
     &               a(mshgl),        a(mwght),      iBC,
     &               BC,              engBC,         a(mshpb),
     &               a(mshglb),       a(mwghtb),     rtmp,
     &               iper,            ilwork, ac)
!Added ac to the end if itrRes, but not tested - Nicholas
c
        iabres = 0
c
        rtmp = rtmp**2
        call sumgat (rtmp, nflow, summed, ilwork)
        epsA = (epsM**2) * sqrt(summed)
c
c.... compute the norm of the second derivative (truncation error)
c
c.... set interval
c
        epsSD = sqrt(epsM)
c
c.... compute the first residual
c
        rtmp = zero
c
        ytmp = ypre + epsSD * uBrg
c
        call itrRes (ytmp,            y,
     &               x,               a(mshp),
     &               a(mshgl),        a(mwght),      iBC,
     &               BC,              engBC,         a(mshpb),
     &               a(mshglb),       a(mwghtb),     rtmp,
     &               iper,            ilwork, ac)
!Added ac to the end if itrRes, but not tested - Nicholas
c
c.... compute the second residual and add it to the first one
c
        ytmp = ypre - epsSD * uBrg
c
        call itrRes (ytmp,            y,
     &               x,               a(mshp),
     &               a(mshgl),        a(mwght),      iBC,
     &               BC,              engBC,         a(mshpb),
     &               a(mshglb),       a(mwghtb),     rtmp,
     &               iper,            ilwork, ac)
!Added ac to the end if itrRes, but not tested - Nicholas
c
c.... compute the second derivative and its norm
c
        rtmp = (( rtmp - two * rmes ) / epsM)**2
c
        call sumgat (rtmp, nflow, summed, ilwork)
        SDnrm = sqrt(summed)
c
c.... compute the 'optimum' interval
c
        eGMRES = two * sqrt( epsA / SDnrm )
c
c.... flop count
c
   !      flops = flops + 10*nflow*nshg+3*nshg
c
c.... end
c
        return
        end
