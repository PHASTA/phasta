c     
c--------------
c     drvAllreduce
c--------------
c     
      subroutine drvAllreduce ( eachproc, result, m )
c     
      include "common.h"
      include "mpif.h"
c     
      dimension eachproc(m), result(m)
c     
      if (numpe > 1) then
         if(impistat.eq.1) then
           iAllR = iAllR+1
         elseif(impistat.eq.2) then
           iAllRScal = iAllRScal+1
         endif
         if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         if(impistat.gt.0) rmpitmr = TMRC()
         call MPI_ALLREDUCE ( eachproc, result, m, 
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
         if(impistat.eq.1) then
           rAllR = rAllR+TMRC()-rmpitmr
         elseif(impistat.eq.2) then
           rAllRScal = rAllRScal+TMRC()-rmpitmr
         endif
      else
         result = eachproc
      endif
c     
      return
      end
c     
c------------------
c     drvAllreducesclr
c------------------
c     
      subroutine drvAllreducesclr ( eachproc, result ) 
c     
      include "common.h"
      include "mpif.h"
c     
      if (numpe > 1) then
         if(impistat.eq.1) then
           iAllR = iAllR+1
         elseif(impistat.eq.2) then 
           iAllRScal = iAllRScal+1
         endif
         if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         if(impistat.gt.0) rmpitmr = TMRC()
         call MPI_ALLREDUCE ( eachproc, result, 1,
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
         if(impistat.eq.1) then
           rAllR = rAllR+TMRC()-rmpitmr
         elseif(impistat.eq.2) then
           rAllRScal = rAllRScal+TMRC()-rmpitmr
         endif
      else
         result = eachproc
      endif
c     
      return
      end

c------------------------------------------------------------------------
c
c   sum real*8 array over all processors
c
c------------------------------------------------------------------------
      subroutine sumgat (u, n, summed)

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      dimension u(nshg,n), ilwork(nlwork) 
!SCATTER      dimension sumvec(numpe), irecvcount(numpe)

      summed = sum(u)

      if (numpe > 1) then
         irecvcount = 1
         sumvec = summed
         if(impistat.eq.1) then
           iAllR = iAllR+1
         elseif(impistat.eq.2) then 
            iAllRScal = iAllRScal+1
         endif
         if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         if(impistat.gt.0) rmpitmr = TMRC()
         call MPI_ALLREDUCE (sumvec, summed, irecvcount, 
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         if(impistat.eq.1) then 
           rAllR = rAllR+TMRC()-rmpitmr
         elseif(impistat.eq.2) then
           rAllRScal = rAllRScal+TMRC()-rmpitmr
         endif
c         call MPI_REDUCE_SCATTER (sumvec, summed, irecvcount, 
c     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      endif

      return
      end

c------------------------------------------------------------------------
c
c   sum real*8 array of length nnp over all processors
c
c------------------------------------------------------------------------
      subroutine sumgatN (u, n, summed, nnp)

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      dimension u(nnp,n), ilwork(nlwork) 
!      dimension sumvec(numpe), irecvcount(numpe)

c protect against underflow
c     summed = sum(u)
      summed = sum(u) + 1.e-20

      if (numpe > 1) then
         irecvcount = 1
         sumvec = summed

         if(impistat.eq.1) then
           iAllR = iAllR+1
         elseif(impistat.eq.2) then
            iAllRScal = iAllRScal+1
         endif
         if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         if(impistat.gt.0) rmpitmr = TMRC()
         call MPI_ALLREDUCE (sumvec, summed, irecvcount, 
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
         if(impistat.eq.1) then 
           rAllR = rAllR+TMRC()-rmpitmr
         elseif(impistat.eq.2) then 
           rAllRScal = rAllRScal+TMRC()-rmpitmr
         endif
c         call MPI_REDUCE_SCATTER (sumvec, summed, irecvcount, 
c     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      endif

      return
      end

c------------------------------------------------------------------------
c
c   sum integer array over all processors
c
c------------------------------------------------------------------------
      subroutine sumgatInt (u, n, summed )

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer u(n), summed, sumvec
!SCATTER      integer sumvec(numpe), irecvcount(numpe)

c$$$      ttim(62) = ttim(62) - tmr()

      summed = sum(u)

      if (numpe > 1) then
         irecvcount = 1
         sumvec = summed
         if(impistat.eq.1) then 
           iAllR = iAllR+1
         elseif(impistat.eq.2) then 
           iAllRScal = iAllRScal+1
         endif
         if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         if(impistat.gt.0) rmpitmr = TMRC()
         call MPI_ALLREDUCE (sumvec, summed, irecvcount, 
     &        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         if(impistat.eq.1) then
           rAllR = rAllR+TMRC()-rmpitmr
         elseif(impistat.eq.2) then
           rAllRScal = rAllRScal+TMRC()-rmpitmr
         endif
c         call MPI_REDUCE_SCATTER (sumvec, summed, irecvcount, 
c     &        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      endif
c$$$      ttim(62) = ttim(62) + tmr()

      return
      end

