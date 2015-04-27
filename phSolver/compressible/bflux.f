        subroutine Bflux (y,         ac,        x,         
     &                    shp,       shgl,      shpb,  
     &                    shglb,     nodflx,    ilwork)
c
c----------------------------------------------------------------------
c
c This routine :
c   1. computes the boundary fluxes
c   2. prints the results in the file [FLUX.lstep]
c   3. prints the primitives variables in the files [OUTPUT.lstep]
c      and [OUTCHM.lstep] (2-D computations only)
c   4. calls the restar routine
c
c
c output:
c  in file FLUX.lstep:
c     run title
c     step number, run time
c     node_number
c     x_1      ... x_nsd                 ! coordinates
c     normal_1 ... normal_nsd            ! outward normal direction
c     tau_1n   ... tau_nsd n             ! boundary viscous flux
c     q_n                                ! boundary head flux
c     ro  u_1  ... u_nsd  t  p  s  Mach  ! primitive variables
c     y_N2 y_O2 y_NO y_N y_O             ! gas composition
c
c  in file OUTPUT.lstep:
c     run title
c     step number, run time, gamma, Rgas
c     density, velocity vector and temperature
c
c  in file OUTCHM.lstep:
c     run title
c     step number, run time, gamma, Rgas
c     pressure, entropy, Mach, y_N2, y_O2, y_NO, y_N and y_O
c
c
c Zdenek Johan, Summer 1991.
c----------------------------------------------------------------------
c
        use pointer_data
c
        include "common.h"
c
        dimension y(nshg,ndof),           ac(nshg,ndof),
     &            x(numnp,nsd),              
     &            shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            shpb(MAXTOP,maxsh,MAXQPT),  
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT),
     &            nodflx(numflx)
c
        dimension invflx(nshg),             flxres(nshg,nflow),
     &            flxLHS(nshg,1),           flxnrm(nshg,nsd),
     &            temp(nshg),               itemp(numflx),
     &            q(nshg,ndof),             qq(nshg,8),
     &            ilwork(nlwork)
c
        character(len=20) fname1,           fname2,
     &                    fmt1,             fmt2
c
c.... ************************>>  Restart  <<***************************
c
c.... set up the timer
c
        call timer ('Output  ')
c
c.... call restart option
c
        q(:,1:ndof)=y(:,1:ndof)

        if (irs .ge. 1) call restar ('out ',  q,ac)
c
c      No Flux output for now KENJ
c
        return
c
c.... scale the primitive variables for output
c
        q(:,1)    = ro     * q(:,1)
        q(:,2)    = vel    * q(:,2)
        q(:,3)    = vel    * q(:,3)
        if (nsd .eq. 3)
     &  q(:,4)    = vel    * q(:,4)
        q(:,nflow) = temper * q(:,nflow)
        qq(:,1)   = press  * qq(:,1)
        qq(:,2)   = entrop * qq(:,2)
c
c.... *************************>>  Output  <<***************************
c
c.... output only for 2-D computations
c
        if ((irs .eq. 2) .and. (nsd .eq. 2)) then
c
c.... get the file names ('output.'lstep) and ('outchm.'lstep)
c
        itmp = 1
        if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
        write (fmt1,"('(''output.'',i',i1,',1x)')") itmp
        write (fmt2,"('(''outchm.'',i',i1,',1x)')") itmp
        write (fname1,fmt1) lstep
        write (fname2,fmt2) lstep
c
c.... open file
c
        if (ioform .eq. 0) then
          open (unit=iout,   file=fname1, status='unknown', err=995)
          if (ichem .eq. 1)
     &    open (unit=ichmou, file=fname2, status='unknown', err=996)
        else
          open (unit=iout,   file=fname1, status='unknown', err=995,
     &          form='unformatted')
          if (ichem .eq. 1)
     &    open (unit=ichmou, file=fname2, status='unknown', err=996,
     &          form='unformatted')
        endif
c
c.... output header
c
        if (ioform .eq. 0) then
          write (iout,  1000) title, lstep, time, gamma, Rgas
          if (ichem .eq. 1)
     &    write (ichmou,1000) title, lstep, time
        else
          write (iout) title
          write (iout) nshg, lstep, time, gamma, Rgas
c
          if (ichem .eq. 1) then
          write (ichmou) title
          write (ichmou) nshg, lstep, time
          endif
        endif
c
        if (ioform .eq. 0) then
          do n = 1, nshg
            write (iout,  2000) n, (q(n,i),  i=1,ndof)
            if (ichem .eq. 1)
     &      write (ichmou,2000) n, (qq(n,i), i=1,8)
          enddo
        else
          write (iout) q
          if (ichem .eq. 1) write (ichmou) qq
        endif
c
c.... close the files
c
        close (iout)
        if (ichem .eq. 1) close (ichmou)
c
        endif
c
c.... stop the timer
c
        call timer ('Back    ')
c
c.... *********************>>  Boundary Fluxes  <<**********************
c
        if ((irs .eq. 2) .and. (numflx .ne. 0)) then
c
c.... set up the timer
c
        call timer ('Bnd_Flux')
c
c.... calculate the inverse of nodflx
c
        itemp  = nodflx
c
        invflx = 0
        invflx(itemp) = (/ (i, i=1,numflx) /)
c
c.... -------------------->   interior elements   <--------------------
c
c.... set up parameters
c
        intrul = intg  (1,itseq)
        intind = intpt (intrul)
c
        jump   = min(iALE + iter-1, 1)
        ires   = 1
        iprec  = 0
c
c.... initialize the arrays
c
        flxres = zero
        flxLHS = zero
        flxnrm = zero
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c
c$$$          iel    = lcblk(1,iblk)
c$$$          nenl   = lcblk(5,iblk)
c$$$          mattyp = lcblk(7,iblk)
c$$$          npro   = lcblk(1,iblk+1) - iel 
c
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          ngauss = nint(lcsyst)
c
c
          if (mattyp .ne. 0) cycle      ! fluid elements only
c
c.... compute and assemble the residuals
c
          call AsIFlx (y,           ac,
     &                 x,           mxmudmi(iblk)%p,        
     &                 shp,         shgl,     mien(iblk)%p,
     &                 mmat(iblk)%p,            flxres)
c
c.... end of interior element loop
c
        enddo
c
c.... -------------------->   boundary elements   <--------------------
c
c.... loop over the elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
c$$$          iel    = lcblkb(1,iblk)
c$$$          nenl   = 4      ! tetrahedral elements (4-vertices)
c$$$          nenbl  = 3      ! triangular boundary element faces
c$$$          mattyp = lcblkb(7,iblk)
c$$$          npro   = lcblkb(1,iblk+1) - iel 
c
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          if(lcsyst.eq.3) lcsyst=nenbl
          ngaussb = nintb(lcsyst)
          
          if (mattyp .ne. 0) cycle      ! fluid elements only
c
c.... compute and assemble the residuals
c
          call AsBFlx (y,                       x,
     &                 shpb(lcsyst,1:nshl,:),
     &                 shglb(lcsyst,:,1:nshl,:),
     &                 mienb(iblk)%p,           mmatb(iblk)%p,
     &                 miBCB(iblk)%p,           mBCB(iblk)%p,
     &                 invflx,                  flxres,
     &                 flxLHS,                  flxnrm)
c
c.... end of boundary element loop
c
        enddo
c
c.... ---------------------------->  Communications <-------------------
c
        if(numpe > 1) then
        call commu (flxres, ilwork, nflow, 'in ')
        call commu (flxLHS, ilwork, 1   , 'in ')
        call commu (flxnrm, ilwork, nsd , 'in ')
        endif
c
c.... ---------------------------->  Solve  <---------------------------
c
c.... compute the viscous and heat fluxes
c
        do i = 2, nflow
          where ( (invflx .ne. 0) .and. (flxLHS(:,1) .ne. zero) )
     &      flxres(:,i) = flxres(:,i) / flxLHS(:,1)
        enddo
c
c.... normalize the outward normal
c
c       if (nsd .eq. 2) then
c
c         temp = sqrt(flxnrm(:,1)**2 + flxnrm(:,2)**2)
c         where ( (invflx .ne. 0) .and. (temp .ne. zero) )
c           flxnrm(:,1) = flxnrm(:,1) / temp
c           flxnrm(:,2) = flxnrm(:,2) / temp
c         endwhere
c
c       else
c
          temp = sqrt(flxnrm(:,1)**2 + flxnrm(:,2)**2 + flxnrm(:,3)**2)
          where ( (invflx .ne. 0) .and. (temp .ne. zero) )
            flxnrm(:,1) = flxnrm(:,1) / temp
            flxnrm(:,2) = flxnrm(:,2) / temp
            flxnrm(:,3) = flxnrm(:,3) / temp
          endwhere
c
c       endif
c
c.... ---------------------------->  Print  <---------------------------
c
c.... get the file name ('flux.'lstep)
c
        itmp = 1
        if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
        write (fmt1,"('(''flux.'',i',i1,',1x)')") itmp
        write (fname1,fmt1) lstep
c
c.... open file
c
        open (unit=iflux, file=fname1, status='unknown', err=997)
c
c.... output the header
c
        write (iflux, 1000) title, lstep, time
c
c.... output the results
c
        do n = 1, numflx
c
          k = nodflx(n)
          m = k                       ! global node number
c
          write (iflux, 2000) k, (x(m,i),       i=1,nsd), 
     &                           (flxnrm(m,i),  i=1,nsd),
     &                           (flxres(m,i),  i=2,nflow),
     &                           (q(k,i),       i=1,ndof),
     &                           (qq(k,i),      i=1,8)
c
        enddo
c
c.... close file
c
        close (iflux)
c
c.... stop the timer
c
        call timer ('Back    ')
c
        endif
c
        return
c
c.... file error handling
c
995     call error ('bflux   ','opening ', iout)
996     call error ('bflux   ','opening ', ichmou)
997     call error ('bflux   ','opening ', iflux)
c
1000    format(' ',a80,/,1x,i10,1p,3e20.7)
2000    format(1p,1x,i6,23e15.7)
c
c.... end
c
        end
