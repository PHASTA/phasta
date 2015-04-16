        subroutine rstat (res, ilwork)
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


        real*8  ftots(3,0:MAXSURF),ftot(3),spmasstot(0:MAXSURF),spmasss

	ttim(68) = ttim(68) - secs(0.0)

        if (numpe == 1) nshgt=nshg   ! global = this processor
c
c incompressible style data from flx surface
c
      if (numpe > 1) then
         call MPI_ALLREDUCE (flxID(2,isrfIM), spmasss,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE (flxID(1,isrfIM), Atots,1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE (flxID(3,:), Ftots(1,:),MAXSURF+1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE (flxID(4,:), Ftots(2,:),MAXSURF+1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE (flxID(5,:), Ftots(3,:),MAXSURF+1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
         call MPI_ALLREDUCE (flxID(2,:), spmasstot(:),MAXSURF+1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
      else
         Ftots=flxID(3:5,:)
         Atots=flxID(1,isrfIM)
         spmasss=flxID(2,isrfIM)
         spmasstot(:)=flxID(2,:)
      endif
	if(myrank.eq.0) then	
      write(44,1000)lstep+1,(spmasstot(j),j=1,5)
      call flush(44)
	endif	
      ftot(1)=sum(Ftots(1,0:MAXSURF))
      ftot(2)=sum(Ftots(2,0:MAXSURF))
      ftot(3)=sum(Ftots(3,0:MAXSURF))
c
c end of incompressible style
c
c
c.... -------------------->  Aerodynamic Forces  <----------------------
c
c.... output the forces and the heat flux
c
        if (iter .eq. nitr) then
          Forin = (/ Force(1), Force(2), Force(3), HFlux /)
          if (numpe > 1) then
          call MPI_REDUCE (Forin(1), Forout(1), 4, MPI_DOUBLE_PRECISION,
     &                                   MPI_SUM, master, 
     &                                   MPI_COMM_WORLD,ierr)
          endif
          Force = Forout(1:3)
          HFlux = Forout(4)
          if (myrank .eq. master) then
             write (iforce,1000) lstep+1, (Force(i), i=1,nsd), HFlux, 
     &                           spmasss
             call flush(iforce)
          endif
        endif

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
        resmax = sqrt(resmax)
        if (resfrt .eq. zero) resfrt = totres
        jtotrs = int  ( 10.d0 * log10 ( totres / resfrt ) )
        jresmx = int  ( 10.d0 * log10 ( resmax / totres ) )
c
c.... get the CPU-time
c
        rsec=TMRC()
        cputme = (rsec-ttim(100))
c
c.... output the result
c
        if (myrank .eq. master) then
          print 2000,        lstep+1, cputme, totres, jtotrs, nrsmax,
     &                     jresmx, lGMRES,  iKs, ntotGM
          write (ihist,2000) lstep+1, cputme, totres, jtotrs, nrsmax,
     &                     jresmx, lGMRES,  iKs, ntotGM
          call flush(ihist)
        endif
	ttim(68) = ttim(68) + secs(0.0)

c
c.... return
c
        return
c
1000    format(1p,i6,5e13.5)
2000    format(1p,i6,e10.3,e10.3,3x,'(',i4,')',3x,'<',i6,'|',i4,'>',
     &         ' [',i3,'-',i3,']',i10)
c
        end
        subroutine rstatSclr (rest, ilwork,lgmrest,ikst)
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
        resmax = sqrt(resmax)
        if (resfrt .eq. zero) resfrt = totres
        jtotrs = int  ( 10.d0 * log10 ( totres / resfrt ) )
        jresmx = int  ( 10.d0 * log10 ( resmax / totres ) )
c
c.... get the CPU-time
c
        rsec=TMRC()
        cputme = (rsec-ttim(100))
c
c.... output the result
c
        if (myrank .eq. master) then
          print 2000,        lstep+1, cputme, totres, jtotrs, nrsmax,
     &                     jresmx, lgmrest,  iKst, ntotGM
          write (ihist,2000) lstep+1, cputme, totres, jtotrs, nrsmax,
     &                     jresmx, lgmrest,  iKst, ntotGM
          call flush(ihist)
        endif
        if(totres.gt.1.0e-9) istop=istop-1

	ttim(68) = ttim(68) + secs(0.0)

c
c.... return
c
        return
c
1000    format(1p,i6,4e13.5)
2000    format(1p,i6,e10.3,e10.3,3x,'(',i4,')',3x,'<',i6,'|',i4,'>',
     &         ' [',i3,'-',i3,']',i10)
c
        end
