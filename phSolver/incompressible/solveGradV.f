      subroutine solveGradV( rmass, qres, iBC, iper, ilwork )
c---------------------------------------------------------------------
c
c This routine satisfies the periodic boundary conditions
c on the diffusive flux residual and mass matrix
c
c input:
c     rmass   (nshg)              : mass matrix
c     qres    (nshg,(nflow-1)*nsd) : diffusive flux vector
c 
c output: modified qres and rmass 
c---------------------------------------------------------------------
      include "common.h"
      
      dimension rmass(nshg), qres(nshg,nsdsq),
     &          iBC(nshg), iper(nshg)
c
c
c.... compute qi for node A, i.e., qres <-- qres/rmass
c
       if (numpe > 1) then
          call commu (qres  , ilwork, nsdsq  , 'in ')
       endif
c
c  take care of periodic boundary conditions
c  but not on surface tension terms in qres(:,10-12)
c  that are used to compute normal vector
c
        !write(*,*) 'nflow, nsd, idflx, idflow:',nflow,nsd,idflx,idflow
        !idflow = (nflow-1)*nsd
        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
c            qres(i,:) = qres(i,:) + qres(j,:)
!            qres(i,1:idflow) = qres(i,1:idflow) + qres(j,1:idflow)
            qres(i,1:nsdsq) = qres(i,1:nsdsq) + qres(j,1:nsdsq)
          endif
        enddo

        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
c            qres(j,:) = qres(i,:)
!            qres(j,1:idflow) = qres(i,1:idflow)
            qres(j,1:nsdsq) = qres(i,1:nsdsq)
          endif
        enddo

       ! rmass has already been computed and inversed in qpbc.f
       do i=1, nsdsq
          qres(:,i) = rmass*qres(:,i)
       enddo

       if(numpe > 1) then
          call commu (qres, ilwork, nsdsq, 'out')    
       endif

c
c.... return
c    
        return
        end

