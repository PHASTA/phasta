      subroutine ctypes (ilwork)

      parameter (maxseg = 60000)

      include "commonM2N.h"
      include "mpif.h"
      include "auxmpiM2N.h"

      integer  sizeofdouble

      dimension ilwork(nlwork)
      dimension isbegin(maxseg), lenseg(maxseg), ioffset(maxseg)

      real*8 timer1,timer2,time_ellapsed
      integer istatRanks(2)
      real*8 rStats(4)

      timer1 = TMRC()

      CALL MPI_TYPE_EXTENT (MPI_DOUBLE_PRECISION,sizeofdouble,ierr)
      lstride = nshg * sizeofdouble
c
c.... maxfront is a common variable being set in this routine
c
      maxfront = 0  
      numtask = ilwork (1)
      itkbeg  = 1

      if (numtask .gt. maxtask) 
     &  call error('ctypes  ','numtask ',numtask)

      nshg0 = nshg

      do itask = 1,numtask
c
c.... iacc = 0 ==> this task is a send
c          = 1 ==> this task is a recieve
c
        iacc   = ilwork (itkbeg + 2)
c
c.... numseg : number of data segments to be communicated
c
        numseg = ilwork (itkbeg + 4)
c
c.... adjust the number of the other processor, since processors
c     are numbered here starting from 0, not 1.
c
        ilwork (itkbeg + 3) = ilwork (itkbeg + 3) - 1 
        if (numseg .gt. maxseg) 
     &    call error('ctypes  ','numseg  ',numseg )
c
c.... lfront = total number of nodes involved in this task
c     
        lfront = 0
        do is = 1,numseg
c
c.... isbegin(is): starting node number for each segment
c
          isbegin (is) = ilwork (itkbeg + 3 + 2*is)
c
c.... lenseg(is): length of each segment (number of nodes)
c
          lenseg  (is) = ilwork (itkbeg + 4 + 2*is)
c
c.... increment the total node counter
c
          lfront       = lfront + lenseg(is)
c
c.... nshg0: number of nodes to be assembled on this processor,
c             i.e. subtract the number of nodes which will be 
c             sent to another processor.
c
	  if (iacc .eq. 0) nshg0 = nshg0 - lenseg(is)
        enddo
c
c.... maxfront: number of nodes which will be communicated, including
c               all segments. Note that after the loop over tasks
c               is complete, maxfront will contain the maximum number
c               of nodes for any of the tasks.
c
        maxfront = MAX(maxfront,lfront)
c
c.... ioffset: array offset from the first node in the first segment
c
        ioffset(1:numseg) = isbegin(1:numseg) - isbegin(1)
c
c.... now set up the MPI data types which will be used in commu.f.
c     These data types represent the indexed sets that will be sent
c     and recieved.
c 
c
c.... the following call to MPI_TYPE_INDEXED will create a new data
c     type which will represent the blocks of data we wish to transfer
c     for this task. A handle to the new type is returned 
c     (sevsegtype(itask,1)). This data type describes the blocks of
c     data to be transferred in terms of segments.
c     Input to this routine:
c          numseg: number of segments in this task
c          lenseg: length of each segment (number of nodes)
c          ioffset: where to begin each block with respect to the
c                   first segment
c          MPI_DOUBLE_PRECISION: type to set for each of the blocks
c
        call MPI_TYPE_INDEXED (numseg, lenseg, ioffset,
     &                  MPI_DOUBLE_PRECISION, sevsegtype(itask,1), ierr)
c
c.... now create a new data type for each of the types of arrays we 
c     may wish to communicate with. For example ndof will be used when
c     communicating the residual vector. Each one of these is derived
c     from the first data type defined above, sevsegtype(itask,1).
c
        call MPI_TYPE_HVECTOR(nsd,    1, lstride, sevsegtype(itask,1),
     &                                        sevsegtype(itask,2), ierr)
c
        call MPI_TYPE_HVECTOR(ndof,   1, lstride, sevsegtype(itask,1),
     &                                        sevsegtype(itask,3), ierr)
c
        call MPI_TYPE_HVECTOR(nflow*nflow,1, lstride,
     &                    sevsegtype(itask,1),sevsegtype(itask,4), ierr)
        call MPI_TYPE_HVECTOR((nflow-1)*nsd,1,lstride,
     &                    sevsegtype(itask,1),sevsegtype(itask,5), ierr)
       call MPI_TYPE_HVECTOR(nflow,1,lstride,sevsegtype(itask,1),
     &                                        sevsegtype(itask,6), ierr)
       call MPI_TYPE_HVECTOR(24,1,lstride,sevsegtype(itask,1),
     &                                        sevsegtype(itask,7), ierr)
       call MPI_TYPE_HVECTOR(9,1,lstride,sevsegtype(itask,1),
     &                                        sevsegtype(itask,8), ierr)
       call MPI_TYPE_HVECTOR(11,1,lstride,sevsegtype(itask,1),
     &                                        sevsegtype(itask,9), ierr)
       call MPI_TYPE_HVECTOR(7,1,lstride,sevsegtype(itask,1),
     &                                   sevsegtype(itask,10), ierr)
       call MPI_TYPE_HVECTOR(33,1,lstride,sevsegtype(itask,1),
     &                                   sevsegtype(itask,11), ierr)
       call MPI_TYPE_HVECTOR(22,1,lstride,sevsegtype(itask,1),
     &                                   sevsegtype(itask,12), ierr)
       call MPI_TYPE_HVECTOR(16,1,lstride,sevsegtype(itask,1),
     &                                   sevsegtype(itask,13), ierr)
       call MPI_TYPE_HVECTOR(10,1,lstride,sevsegtype(itask,1),
     &                                   sevsegtype(itask,14), ierr)
       call MPI_TYPE_HVECTOR(nflow*nsd,1,lstride,sevsegtype(itask,1),
     &                                   sevsegtype(itask,15), ierr)
c
c
c.... now this must be done to make MPI recognize each of the data
c     types that were just defined
c
        do kdof = 1,15
          call MPI_TYPE_COMMIT (sevsegtype(itask,kdof), ierr)
        enddo
c
c.... set the counter to the index in ilwork where the next task
c     begins
c

        itkbeg = itkbeg + 4 + 2*numseg
c
c.... end loop over tasks
c
      enddo

      timer2 = TMRC()
      time_ellapsed = timer2 - timer1
      !write(*,*) 'myrank: ', myrank, '- time in ctype: ',time_ellapsed 
      call rgetMinMaxAvg(time_ellapsed,rStats,istatRanks)
      if(myrank == 0) then
         write(*,811) istatRanks(1),rStats(1),istatRanks(2),rStats(2),
     &             rStats(3),rStats(4)
      endif

811   format('Time in ctype: ',
     &                'min [',I6,',',F10.5,'], ',
     &                'max [',I6,',',F10.5,'], ',
     &                'and avg [.,',F10.5,'], (rms=',F10.5,')')


      return
      end
