        subroutine rstatCheck (res, ilwork,y,ac)
c
c----------------------------------------------------------------------
c
c This subroutine calculates the statistics of the residual.
c
c input:
c  res   (nshg,nflow)   : preconditioned residual
c
c output:
c  The time step, cpu-time and entropy-norm of the residual
c     are printed in the file HISTOR.DAT.
c  
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension res(nshg,nflow)
        dimension rtmp(nshg), nrsmax(1), ilwork(nlwork)
        dimension Forin(4), Forout(4)
!SCATTER        dimension irecvcount(numpe), resvec(numpe)
c        integer TMRC
        real*8 y(nshg,ndof),ac(nshg,ndof)
        save ResLast

        if (numpe == 1) nshgt=nshg   ! global = this processor
c
c.... ----------------------->  Convergence  <-------------------------
c
c.... compute the maximum residual and the corresponding node number
c
        rtmp = zero
        do i = 1, nflow
          rtmp = rtmp + res(:,i)**2
        enddo
 
        call sumgat (rtmp, 1, resnrm, ilwork)
      
        resmaxl = maxval(rtmp)

        irecvcount = 1
        resvec = resmaxl
        if (numpe > 1) then
        call MPI_ALLREDUCE (resvec, resmax, irecvcount,
     &                    MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD,
     &                    ierr)
c        call MPI_REDUCE_SCATTER (resvec, resmax, irecvcount,
c     &                    MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD,
c     &                    ierr)
        else
          resmax=resmaxl
        endif
        nrsmax = maxloc(rtmp)
c
c.... correct the residuals
c
        if (loctim(itseq) .eq. 0) then
          resnrm = resnrm 
          resmax = resmax 
        else
          resnrm = resnrm 
          resmax = resmax 
        endif
c
c.... approximate the number of entries
c
        totres = resnrm / float(nshgt)
        totres = sqrt(totres)
       if((iter.gt.1).and.(totres.gt.10000.0*ResLast)) then !diverging
               call restar('out ',y,res) ! 'res' is used instead of 'ac'
               if(myrank.eq.0) write(*,*) 'ResLast totres', ResLast, totres
               if(myrank.eq.0) write(*,*) 'resmax', resmax
               if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
               call error('rstat    ','Diverge', iter)
       endif
       ResLast=totres
	ttim(68) = ttim(68) + secs(0.0)

c
c.... return
c
        return
c
        end
        subroutine rstatCheckSclr (rest, ilwork,y,ac)
c
c----------------------------------------------------------------------
c
c This subroutine calculates the statistics of the residual.
c
c input:
c  rest   (nshg)   : preconditioned residual
c
c output:
c  The time step, cpu-time and entropy-norm of the residual
c     are printed in the file HISTOR.DAT.
c  
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension rest(nshg)
        dimension rtmp(nshg), nrsmax(1), ilwork(nlwork)
!SCATTER        dimension irecvcount(numpe), resvec(numpe)
c        integer TMRC
        real*8 y(nshg,ndof),ac(nshg,ndof)
        save ResLast
        save lstepLast

	ttim(68) = ttim(68) - secs(0.0)
        if (numpe == 1) nshgt=nshg   ! global = this processor
c
c.... ----------------------->  Convergence  <-------------------------
c
c.... compute the maximum residual and the corresponding node number
c
        rtmp = zero
        rtmp = rtmp + rest**2
 
        call sumgat (rtmp, 1, resnrm, ilwork)
      
        resmaxl = maxval(rtmp)

continue on

        irecvcount = 1
        resvec = resmaxl
        if (numpe > 1) then
        call MPI_ALLREDUCE (resvec, resmax, irecvcount,
     &                    MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD,
     &                    ierr)
c        call MPI_REDUCE_SCATTER (resvec, resmax, irecvcount,
c     &                    MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD,
c     &                    ierr)
        else
          resmax=resmaxl
        endif
        nrsmax = maxloc(rtmp)
c
c.... correct the residuals
c
        if (loctim(itseq) .eq. 0) then
          resnrm = resnrm 
          resmax = resmax 
        else
          resnrm = resnrm 
          resmax = resmax 
        endif
c
c.... approximate the number of entries
c
        totres = resnrm / float(nshgt)
        totres = sqrt(totres)
	if((lstep.gt.0).and.(lstepLast.eq.lstep)) then
           if(totres.gt.10000.0*ResLast) then !diverging
               lstep = lstep+1
               ac(:,5) = rest(:) ! T dot in 'ac' is filled with scl. res
               call restar('out ',y,ac)
               if(myrank.eq.0) write(*,*) 'ResLast totres', ResLast, totres
               if(myrank.eq.0) write(*,*) 'resmax', resmax
               call error('rstatSclr','Diverge', iter)
           endif
	else
		lstepLast=lstep
	endif
       ResLast=totres
c
c.... return
c
        return
c
        end
