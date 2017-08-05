#ifdef USE_CATALYST
c...==============================================================
c... subroutine to do the coprocessing
c... The subroutine is responsible for determining if coprocessing
c... is needed this timestep and if it is needed then the
c... subroutine passes the phasta data structures into
c... the coprocessor. This is meant to be called at the end of
c... every time step.
c... The input is:
c...    itimestep -- the current simulation time step
c...    X -- the coordinates array of the nodes
c...    Y -- the fields array (e.g. velocity, pressure, etc.)
c...    compressibleflow -- flag to indicate whether or not the
c...                         flow is compressible.  if it is then
c...                         temperature will be outputted
c...    computevort -- flag to indicate whether or not vorticity is computed
c...    VORTICITY -- the vorticity array
c... It has no output and should not change any Phasta data.
c...==============================================================

      subroutine phastacoprocessor(itimestep, X, Y, compressibleflow,
     &                      computevort, VORTICITY, dwal )
      use pointer_data
      include "common.h"
      integer iblk, nenl, npro, j, needflag, i
      integer compressibleflow, itimestep, computevort
      dimension x(numnp,nsd), y(nshg,ndof), vorticity(nshg, 5)
      dimension dwal(nshg)
!      dimension ycontainer(nshg,ndof)
      if(docoprocessing .ne. 1) then
        return
      endif
      
      if(nshg .ne. numnp) then
         print *, 'CoProcessing only setup for when nshg equals numnp'
         return
      endif

c  First check if we even need to coprocess this time step/time
      if(myrank.eq.0)  then
            write(6,*) 'Before calling requestdatadescription, '//
     &                  'itimestep:', itimestep, ', time:', time
      endif
      call requestdatadescription(itimestep, time, needflag)
      if(myrank.eq.0)  then
          write(6,*) 'After calling requestdatadescription, needflag:',
     &                  needflag
      endif
      if(needflag .eq. 0) then
c  We don't need to do any coprocessing now so we can return
         return
      endif
c  Check if we need to create the grid for the coprocessor
      call needtocreategrid(needflag)
      if(needflag .ne. 0) then
c     We do need the grid.
         call createpointsandallocatecells(numnp, X, numel)

        do iblk=1,nelblk
            nenl = lcblk(5,iblk) ! no. of vertices per element
            npro = lcblk(1,iblk+1) - lcblk(1,iblk) ! no. of elemens in block
            call insertblockofcells(npro, nenl, mien(iblk)%p(1,1))
        enddo
      endif ! if needflag .ne. 0 --

c  Inside addfields we check to see if we really need the field or not
!      call addfields(nshg, ndof, Y, compressibleflow, computevort,
!     &                             VORTICITY)

!      do i = 1, nshg
!         ycontainer(i,:) = Y(i,:)
!         ycontainer(i,4) = vorticity(i,5) !Replace pressure by Q
!      enddo

      !call addfields(nshg, ndof, Y, compressibleflow)
      !call addfields(nshg, ndof, ycontainer, compressibleflow)
      call addfields(nshg, ndof, Y, compressibleflow,
     &               vorticity,dwal)

      if(myrank.eq.0)  then
         tcorecp5 = TMRC() 
      endif
      call coprocess()
      if(myrank.eq.0)  then
          tcorecp6 = TMRC()
          write(6,*) 'coprocess: ',tcorecp6-tcorecp5 
      endif
      return
      end

      subroutine catalystinit()
!       FIXME - using the old code for now - make configurable        
!       and fix variable allocation
        include "common.h"
        character*200 pyfilename
        integer stringlength

        pyfilename = "../cpscript.py"
        stringlength = 14
        docoprocessing = 1
        call coprocessorinitializewithpython(pyfilename, stringlength)
      end
#endif