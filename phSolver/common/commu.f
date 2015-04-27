      subroutine commu (global, ilwork,  n,  code)
c---------------------------------------------------------------------
c 
c This subroutine is responsible for interprocessor communication of
c the residual and solution vectors.
c
c input:   
c     global(nshg,n): global vector to be communicated. Note that
c                      this vector is local to the processor, (i.e.
c                      not distributed across processors)
c     ilwork(nlwork):  this is the local interprocessor work array.
c                      This array is local to the processor, (i.e.
c                      each processor has a unique ilwork array.
c     n:               second dimension of the array to be communicated
c     code:            = 'in' for communicating with the residual
c                      = 'out' for cummunicating the solution 
c
c---------------------------------------------------------------------
c
c The array ilwork describes the details of the communications. 
c Each communication step (call of this routine) consists of a 
c sequence of "tasks", where a task is defined as a communication 
c between two processors where data is exchanged. This would imply 
c that for a given processor, there will be as many tasks as there
c are processors with which it must communicate. Details of the 
c ilwork array appear below.
c
c---------------------------------------------------------------------
c
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      integer status(MPI_STATUS_SIZE), ierr
      integer stat(MPI_STATUS_SIZE, 2*maxtask), req(2*maxtask)
      real*8  rDelISend, rDelIRecv, rDelWaitAll

      dimension global(nshg,n),
     &          rtemp(maxfront*n,maxtask),
     &          ilwork(nlwork)
 
      character*3 code

      if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
      if(impistat.gt.0) rDelIRecv = zero
      if(impistat.gt.0) rDelISend = zero
      if(impistat.gt.0) rDelWaitAll = zero

      if (code .ne. 'in ' .and. code .ne. 'out') 
     &  call error ('commu   ','code    ',0)

      if     (n .eq. 1)      then        ! like a scalar
        kdof = 1
      elseif (n .eq. nsd)    then        ! like the normal vectors
        kdof = 2
      elseif (n .eq. ndof)   then        ! res, y, ac, krylov vectors....
        kdof = 3
      elseif (n .eq. nflow*nflow) then     ! bdiag
        kdof = 4
      elseif (n .eq. (nflow-1)*nsd) then  ! qres
        kdof = 5
      elseif (n .eq. nflow) then
        kdof = 6
      elseif (n .eq. 24 ) then
        kdof = 7
      elseif (n .eq. 9) then
        kdof = 8
      elseif (n .eq. 11 ) then
        kdof = 9
      elseif (n .eq. 7 ) then
        kdof = 10
      elseif (n .eq. 33 ) then
         kdof = 11 
      elseif (n .eq. 22 ) then
         kdof = 12
      elseif (n .eq. 16 ) then
         kdof = 13
      elseif (n .eq. 10 ) then
         kdof = 14
      elseif (n .eq. nflow*nsd ) then   !surface tension + qres
         kdof = 15
      else
        call error ('commu   ','n       ',n)
      endif

c... Note that when adding another kdof to the above set, we must
c... also make changes in ctypes.f and auxmpi.h

c---------------------------------------------------------------------
c  ilwork(1): number of tasks
c
c  The following information is contained in ilwork for each task:
c     itag: tag of the communication
c     iacc: == 0 if task is a send
c           == 1 if task is a recieve
c     iother: rank of processor with which this communication occurs
c     numseg: number of data "segments" to be sent or recieved. A 
c             segment is defined as a continuous section of the global
c             vector to be communicated, (i.e. a group of nodes (or,
c             rather, "shape function coefficients") which occur 
c             sequentially in the array global(nshg,n)).
c     isbeg:  location of the first segment in the array owned by the
c             current processor.
c
c The two types of communication are 'in', where the residual is being
c communicated, and 'out', where the solution is being communicated.
c Note that when the type is 'out', senders recieve and recievers send.
c
c The following comment pertains to a communication of type 'in':
c
c     If the task is a send, then all of the numseg segments are
c     sent with a single call to MPI_SEND. Where these segments live in 
c     the array is built into the array sevsegtype, which is a common 
c     array constructed in the subroutine "ctypes.f". In other words,
c     sevsegtype is a data type that describes the indices of the blocks
c     to be sent, in terms of there beginning index, and the length of 
c     each segment. Using this, we can make a single send to take care of
c     all the segments for this task. 
c      
c     If the task is a recieve, then once the vector is recieved, the
c     recieved segments must be added to the correct locations in the
c     current array. These locations are described in ilwork as the
c     beginning position, then the length of the segment.
c     
c---------------------------------------------------------------------
      numtask = ilwork(1)
      
      itkbeg = 1
      m = 0
      idl=0

      DO itask = 1, numtask
        m      = m + 1
        itag   = ilwork (itkbeg + 1)
        iacc   = ilwork (itkbeg + 2)
        iother = ilwork (itkbeg + 3)
        numseg = ilwork (itkbeg + 4)
        isgbeg = ilwork (itkbeg + 5)
c
c.... if iacc == 0, then this task is a send.
c     slave
c
        if (iacc .EQ. 0) then  
c
c.... residual communication
c
          if (code .eq. 'in ') then
            if(impistat.eq.1) then
              iISend = iISend+1
            elseif(impistat.eq.2) then
               iISendScal = iISendScal+1
            endif
            if(impistat.gt.0) rmpitmr = TMRC()
            call MPI_ISEND(global(isgbeg, 1), 1, sevsegtype(itask,kdof), 
     &                     iother, itag, MPI_COMM_WORLD, req(m), ierr)
            if(impistat.gt.0) rDelISend = TMRC()-rmpitmr
            if(impistat.eq.1) then 
              rISend = rISend+rDelISend
            elseif(impistat.eq.2) then 
              rISendScal = rISendScal+rDelISend
            endif
          endif
c
c.... solution communication
c
          if (code .eq. 'out') then
            if(impistat.eq.1) then
              iIRecv = iIRecv+1
            elseif(impistat.eq.2) then
               iIRecvScal = iIRecvScal+1
            endif
            if(impistat.gt.0) rmpitmr = TMRC()
            call MPI_IRECV(global(isgbeg, 1), 1, sevsegtype(itask,kdof), 
     &                     iother, itag, MPI_COMM_WORLD, req(m), ierr)
            if(impistat.gt.0) rDelIRecv = TMRC()-rmpitmr
            if(impistat.eq.1) then 
              rIRecv = rIRecv+rDelIRecv
            elseif(impistat.eq.2) then 
              rIRecvScal = rIRecvScal+rDelIRecv
            endif
          endif
c
c.... if iacc == 1, then this task is a recieve.
c     master
c
        else
          if (code .eq. 'in ') then
c
c.... determine the number of total number of nodes involved in this
c     communication (lfront), including all segments
c
            lfront = 0
            do is = 1,numseg
              lenseg = ilwork (itkbeg + 4 + 2*is)
              lfront = lfront + lenseg
            enddo
c
c.... recieve all segments for this task in a single step
c
            idl=idl+1 ! stands for i Do Later, the number to fix later
            if(impistat.eq.1) then 
              iIRecv = iIRecv+1
            elseif(impistat.eq.2) then 
              iIRecvScal = iIRecvScal+1
            endif
            if(impistat.gt.0) rmpitmr = TMRC()
            call MPI_IRECV(rtemp(1,idl), lfront*n, MPI_DOUBLE_PRECISION, 
     &                     iother, itag, MPI_COMM_WORLD, req(m), ierr)
            if(impistat.gt.0) rDelIRecv = TMRC()-rmpitmr
            if(impistat.eq.1) then
               rIRecv = rIRecv+rDelIRecv
            elseif(impistat.eq.2) then 
               rIRecvScal = rIRecvScal+rDelIRecv
            endif
          endif
          if (code .eq. 'out') then
            if(impistat.eq.1) then 
              iISend = iISend+1
            elseif(impistat.eq.2) then 
              iISendScal = iISendScal+1
            endif
            if(impistat.gt.0) rmpitmr = TMRC()
            call MPI_ISEND(global(isgbeg, 1), 1, sevsegtype(itask,kdof), 
     &                     iother, itag, MPI_COMM_WORLD, req(m), ierr)
            if(impistat.gt.0) rDelISend = TMRC()-rmpitmr
            if(impistat.eq.1) then
               rISend = rISend+rDelISend
            elseif(impistat.eq.2) then 
               rISendScal = rISendScal+rDelISend
            endif
          endif
        endif

        itkbeg = itkbeg + 4 + 2*numseg

      enddo   !! end tasks loop

      if(impistat.eq.1) then
        iWaitAll = iWaitAll+1
      elseif(impistat.eq.2) then
         iWaitAllScal = iWaitAllScal+1
      endif
      if(impistat.gt.0) rmpitmr = TMRC()
      call MPI_WAITALL(m, req, stat, ierr)
      if(impistat.gt.0) rDelWaitAll = TMRC()-rmpitmr
      if(impistat.eq.1) then
        rWaitAll = rWaitAll+rDelWaitAll
        rCommu = rCommu+rDelIRecv+rDelISend+rDelWaitAll
      elseif(impistat.eq.2) then
        rWaitAllScal = rWaitAllScal+rDelWaitAll
        rCommuScal = rCommuScal+rDelIRecv+rDelISend+rDelWaitAll
      endif

c
c     Stuff added below is a delayed assembly of that which was communicated
c     above but due to the switch to non-blocking receivves could not be
c     assembled until after the waitall.  Only necessary for commu "in"
c

      if(code .eq. 'in ') then
         itkbeg=1
         jdl=0
         do j=1,numtask         ! time to do all the segments that needed to be
                                ! assembled into the global vector

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)
            isgbeg = ilwork (itkbeg + 5)
            if(iacc.eq.1) then
               jdl=jdl+1  ! keep track of order of rtemp's
c
c... add the recieved data to the global array on the current processor.
c    Note that this involves splitting up the chunk of recieved data
c    into its correct segment locations for the current processor.
c
               itemp = 1
               do idof = 1,n
                  do is = 1,numseg
                 isgbeg = ilwork (itkbeg + 3 + 2*is)
                 lenseg = ilwork (itkbeg + 4 + 2*is)
                 isgend = isgbeg + lenseg - 1
                 global(isgbeg:isgend,idof) = global(isgbeg:isgend,idof)
     &                                + rtemp (itemp:itemp+lenseg-1,jdl)
                 itemp = itemp + lenseg
                  enddo
               enddo
            endif ! end of receive (iacc=1)
            itkbeg = itkbeg + 4 + 2*numseg
         enddo
      endif  ! commu "in"
      return
      end



