      module rlssave

      real*8, allocatable :: rls(:,:)

      end module

c-----------------------------------------------------------------------------



c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||




      subroutine setrls

      use rlssave

      include "common.h"

      allocate ( rls(nshg,6) )

      return

      end




c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||




      subroutine bardmc (y,      shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x)

      use pointer_data
      use rlssave    ! rls(nshg,6) is retrieved and modified here.

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwtf contains the weights of the 
c                    quad. points.  

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension fres(nshg,33),         fwr(nshg),
     &          strnrm(nshg),         cdelsq(nshg),
     &          xnum(nshg),           xden(nshg),
     &          xmij(nshg,6),         xlij(nshg,6),
     &          xrij(nshg,6),         xhij(nshg,6),
     &          xnude(nfath,2),        xnuder(nfath,2),
     &          nsons(nfath),
     &          strl(numel,maxnint),           
     &          y(nshg,5),            yold(nshg,5),
     &          ifath(nshg),          iper(nshg),
     &          ilwork(nlwork),
     &          x(numnp,3),
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),    
     &          xnutf(nfath),
     &          hfres(nshg,11),       
     &          rho(nshg),           
     &          eps(3),               epst(3),
     &          tracerls(nshg)
     
c
c
c   setup the weights for time averaging of cdelsq (now in quadfilt module)
c
      denom=max((one*lstep),one)
      if(dtavei.lt.0) then
         wcur=one/denom
      else
         wcur=dtavei
      endif  
      whist=1.0-wcur

      fres = zero
      hfres = zero

      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)
      yold(:,5) = y(:,5)

      if(matflg(1,1).eq.0) then ! compressible
      rho(:) = yold(:,1) / (Rgas * yold(:,5))  ! get density at nodes.
      else
      rho(:) = one ! unit density     
      endif

c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      yold(:,5) = 1.0  ! Debugging
c      yold(:,2) = 2.0*x(:,1) - 3*x(:,2) 
c      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
c      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
c                               suitable for the

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss  = nint(lcsyst)
        ngaussf = nintf(lcsyst) 
        
        call hfilter (yold, x, mien(iblk)%p, hfres, 
     &               shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo

      if(numpe>1) call commu (hfres, ilwork, 11, 'in ')
c 
c... account for periodicity in filtered variables
c
      do j = 1,nshg  !    Add on-processor slave contribution to masters
        i = iper(j)
        if (i .ne. j) then
           hfres(i,:) = hfres(i,:) + hfres(j,:)
        endif
      enddo
      do j = 1,nshg ! Set on-processor slaves to be the same as masters
        i = iper(j)
        if (i .ne. j) then
           hfres(j,:) = hfres(i,:)
        endif
      enddo

c... Set off-processor slaves to be the same as their masters

      if(numpe>1)   call commu (hfres, ilwork, 11, 'out')


      hfres(:,5) = one / hfres(:,5) ! one/(volume of hat filter kernel)

      do j = 1, 4
	hfres(:,j) = hfres(:,j) * hfres(:,5)
      enddo	    		
      do j = 6, 11
	hfres(:,j) = hfres(:,j) * hfres(:,5)
      enddo	

c     Uncomment to test dmod

c      do j = 1, 11
c         do i = 1, numnp
c            hfres(i,j) = hfres(14,j)
c         enddo
c      enddo

c     End of dmod test

c... Form the resolved Leonard stress (rls) at each node.

      rls(:,1) = hfres(:,6)  - hfres(:,1)*hfres(:,1)/rho(:) 
      rls(:,2) = hfres(:,7)  - hfres(:,2)*hfres(:,2)/rho(:) 
      rls(:,3) = hfres(:,8)  - hfres(:,3)*hfres(:,3)/rho(:)
      
      tracerls=(rls(:,1)+rls(:,2)+rls(:,3))/three

      rls(:,1) = rls(:,1) - tracerls
      rls(:,2) = rls(:,2) - tracerls
      rls(:,3) = rls(:,3) - tracerls

      rls(:,4) = hfres(:,9)  - hfres(:,1)*hfres(:,2)/rho(:)
      rls(:,5) = hfres(:,10) - hfres(:,1)*hfres(:,3)/rho(:)
      rls(:,6) = hfres(:,11) - hfres(:,2)*hfres(:,3)/rho(:)	

c      rls(:,1) = 4.0d0 ! Debugging
c      rls(:,2) = 5.0d0
c      rls(:,3) = 1.0d0
c      rls(:,4) = 2.1d0
c      rls(:,5) = 1.5d0
c      rls(:,6) = 3.0d0

c... Stored the resolved Leonard stresses in a module.

c...  Compute dissipation of resolved kinetic energy (u_{i} u_{i}) 
c...  due to molecular stresses (epst(3)), to modeled residual stresses 
c...  (epst(1)), and to resolved Leonard stresses (epst(2)).

c      if (lstep .gt. 0) then 

c      eps = zero
c      do iblk = 1,nelblk
c        lcsyst = lcblk(3,iblk) 
c        iel  = lcblk(1,iblk) !Element number where this block begins
c        npro = lcblk(1,iblk+1) - iel
c        lelCat = lcblk(2,iblk)
c        nenl = lcblk(5,iblk)
c        nshl = lcblk(10,iblk)
c        inum  = iel + npro - 1

c        ngauss = nint(lcsyst)
        
c        call disPtion ( yold, x, mien(iblk)%p, 
c     &               shgl(lcsyst,:,1:nshl,:), shp(lcsyst,1:nshl,:),
c     &               eps, mxmudmi(iblk)%p )

c      enddo      

c      if(numpe.gt.1)then
c         call drvAllreduce( eps, epst, 3 )
c      else
c         epst  = eps
c      endif
      
c      if (myrank .eq. master) then
c         write(22,*) one*(lstep+1), epst(1), epst(2)
c         call flush(22)
c      endif

c      endif 


c... Done w/ h-filtering. Begin 2h-filtering.

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss  = nint(lcsyst)
        ngaussf = nintf(lcsyst) 
        
        call twohfilter (yold, x, strl(iel:inum,:), mien(iblk)%p, 
     &               fres, hfres, shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo
c
 
      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss  = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        if (ngauss .ne. ngaussf) then
        call getstrl (yold, x,      mien(iblk)%p,  
     &               strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:))
        endif
      enddo
c
c
C must fix for abc and dynamic model
c      if(any(btest(iBC,12)))   !are there any axisym bc's
c     &      call rotabc(res, iBC, 'in ')
c
      if(numpe>1) call commu (fres, ilwork, 33, 'in ')
c 
c account for periodicity in filtered variables
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:)
           fres(j,:) = fres(i,:)
        endif
      enddo

      if (nfath .eq. 0) then
         fres = fres(iper,:)
         if(numpe>1)   call commu (fres, ilwork, 33, 'out')
c again we will need a  call
c      if(any(btest(iBC,12)))   !are there any axisym bc's
c     &   rotabc(y, iBC, 'out')
      endif

      fres(:,23) = one / fres(:,23)
      do j = 1,22
        fres(:,j) = fres(:,j) * fres(:,23)
      enddo
c     fres(:,24) = fres(:,24) * fres(:,23)
      do j = 25,33
        fres(:,j) = fres(:,j) * fres(:,23)
      enddo
c
c.....at this point fres is really all of our filtered quantities
c     at the nodes
c

      strnrm = sqrt( 
     &  two * (fres(:,10)**2 + fres(:,11)**2 + fres(:,12)**2)
     &  + four * ( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )

      fwr = fwr1 * fres(:,22) * strnrm

      if(matflg(1,1).eq.0) then ! compressible

      xmij(:,1) = -fwr
     &             * pt33 * (two*fres(:,10) - fres(:,11) - fres(:,12))
     &             + pt33 * (two*fres(:,16) - fres(:,17) - fres(:,18))
c
      xmij(:,2) = -fwr
     &             * pt33 * (two*fres(:,11) - fres(:,10) - fres(:,12))
     &             + pt33 * (two*fres(:,17) - fres(:,16) - fres(:,18))
c
      xmij(:,3) = -fwr
     &             * pt33 * (two*fres(:,12) - fres(:,10) - fres(:,11))
     &             + pt33 * (two*fres(:,18) - fres(:,16) - fres(:,17))      

      else

      xmij(:,1) = -fwr
     &             * fres(:,10) + fres(:,16)
      xmij(:,2) = -fwr
     &             * fres(:,11) + fres(:,17) 
      xmij(:,3) = -fwr
     &             * fres(:,12) + fres(:,18) 

      endif

      xmij(:,4) = -fwr * fres(:,13) + fres(:,19)
      xmij(:,5) = -fwr * fres(:,14) + fres(:,20)
      xmij(:,6) = -fwr * fres(:,15) + fres(:,21)

      fres(:,22) = one / fres(:,22)

      xlij(:,1) = fres(:,4) - fres(:,1) * fres(:,1) * fres(:,22)
      xlij(:,2) = fres(:,5) - fres(:,2) * fres(:,2) * fres(:,22)
      xlij(:,3) = fres(:,6) - fres(:,3) * fres(:,3) * fres(:,22)
      xlij(:,4) = fres(:,7) - fres(:,1) * fres(:,2) * fres(:,22)
      xlij(:,5) = fres(:,8) - fres(:,1) * fres(:,3) * fres(:,22)
      xlij(:,6) = fres(:,9) - fres(:,2) * fres(:,3) * fres(:,22)

      xhij(:,1) = fres(:,28) - fres(:,25)*fres(:,25)
      xhij(:,2) = fres(:,29) - fres(:,26)*fres(:,26)
      xhij(:,3) = fres(:,30) - fres(:,27)*fres(:,27)
      xhij(:,4) = fres(:,31) - fres(:,25)*fres(:,26)
      xhij(:,5) = fres(:,32) - fres(:,25)*fres(:,27)
      xhij(:,6) = fres(:,33) - fres(:,26)*fres(:,27)		

      xrij(:,1) = xlij(:,1) - xhij(:,1)
      xrij(:,2) = xlij(:,2) - xhij(:,2)
      xrij(:,3) = xlij(:,3) - xhij(:,3)
      xrij(:,4) = xlij(:,4) - xhij(:,4)
      xrij(:,5) = xlij(:,5) - xhij(:,5)
      xrij(:,6) = xlij(:,6) - xhij(:,6)

      xnum =        xrij(:,1) * xmij(:,1) + xrij(:,2) * xmij(:,2) 
     &                                    + xrij(:,3) * xmij(:,3)
     &     + two * (xrij(:,4) * xmij(:,4) + xrij(:,5) * xmij(:,5)
     &                                    + xrij(:,6) * xmij(:,6))
      xden =        xmij(:,1) * xmij(:,1) + xmij(:,2) * xmij(:,2) 
     &                                    + xmij(:,3) * xmij(:,3)
     &     + two * (xmij(:,4) * xmij(:,4) + xmij(:,5) * xmij(:,5)
     &                                    + xmij(:,6) * xmij(:,6))
      xden = two * xden

c  zero on processor periodic nodes so that they will not be added twice
        do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xnum(j) = zero
            xden(j) = zero
          endif
        enddo

      if (numpe.gt.1 .and. nsons(1).gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xnum(isgbeg:isgend) = zero
                  xden(isgbeg:isgend) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
c
c Description of arrays.   Each processor has an array of length equal
c to the total number of fathers times 2 xnude(nfathers,2). One to collect 
c the numerator and one to collect the denominator.  There is also an array
c of length nshg on each processor which tells the father number of each
c on processor node, ifath(nnshg).  Finally, there is an arry of length
c nfathers to tell the total (on all processors combined) number of sons
c for each father. 
c
c  Now loop over nodes and accumlate the numerator and the denominator
c  to the father nodes.  Only on processor addition at this point.
c  Note that serrogate fathers are collect some for the case where some
c  sons are on another processor
c
      if(nsons(1) .gt. 1) then
         xnude = zero
         do i = 1,nshg
            xnude(ifath(i),1) = xnude(ifath(i),1) + xnum(i)
            xnude(ifath(i),2) = xnude(ifath(i),2) + xden(i)
         enddo
      else                      ! no homogeneity assumed
         xnude(:,1)=xnum(:)
         xnude(:,2)=xden(:)
      endif

c
c Now  the true fathers and serrogates combine results and update
c each other.
c
c ONLY DO THIS IF  1) multiprocessors AND 2) homogenous directions exist 
  
      if(numpe .gt. 1 .and. nsons(1).ne.1 )then
         call drvAllreduce(xnude, xnuder,2*nfath)
c
c  xnude is the sum of the sons for each father on this processor
c
c  xnuder is the sum of the sons for each father on all processor combined
c  (the same as if we had not partitioned the mesh for each processor)
c
c   For each father we have precomputed the number of sons (including
c   the sons off processor). 
c
c   Now divide by number of sons to get the average (not really necessary
c   for dynamic model since ratio will cancel nsons at each father)
c
c         xnuder(:,1) = xnuder(:,1) ! / nsons(:)
c         xnuder(:,2) = xnuder(:,2) ! / nsons(:)
c
c  the next line is c \Delta^2
c
            numNden(:,1) = whist*numNden(:,1)+wcur*xnuder(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnuder(ifath(:),2)
            cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)
      else
c     
c     the next line is c \Delta^2, not nu_T but we want to save the
c     memory
c     
c$$$
c$$$            write(400+myrank,555) (numNden(j*500,1),j=1,5)
c$$$            write(410+myrank,555) (numNden(j*500,2),j=1,5)
c$$$            write(500+myrank,555) (numNden(j+500,1),j=1,5)
c$$$            write(510+myrank,555) (numNden(j+500,2),j=1,5)
c$$$
c$$$            write(430+myrank,555) (xnude(j*500,1),j=1,5)
c$$$            write(440+myrank,555) (xnude(j*500,2),j=1,5)
c$$$            write(530+myrank,555) (xnude(j+500,1),j=1,5)
c$$$            write(540+myrank,555) (xnude(j+500,2),j=1,5)

            numNden(:,1) = whist*numNden(:,1)+wcur*xnude(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnude(ifath(:),2)
            cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)

c 
c  to get started we hold cdelsq fixed
c
c            cdelsq(:) = 2.27e-4

            
            
c$$$            write(450+myrank,555) (cdelsq(j*500),j=1,5)
c$$$            write(460+myrank,555) (y(j*500,1),j=1,5)
c$$$            write(470+myrank,555) (y(j*500,2),j=1,5)
c$$$            write(480+myrank,555) (y(j*500,3),j=1,5)
c$$$            write(490+myrank,555) (strnrm(j*500),j=1,5)
c$$$
c$$$            write(550+myrank,555) (cdelsq(j+500),j=1,5)
c$$$            write(560+myrank,555) (y(j+500,1),j=1,5)
c$$$            write(570+myrank,555) (y(j+500,2),j=1,5)
c$$$            write(580+myrank,555) (y(j+500,3),j=1,5)
c$$$            write(590+myrank,555) (strnrm(j+500),j=1,5)
c$$$
c$$$            call flush(400+myrank)
c$$$            call flush(410+myrank)
c$$$            call flush(430+myrank)
c$$$            call flush(440+myrank)
c$$$            call flush(450+myrank)
c$$$            call flush(460+myrank)
c$$$            call flush(470+myrank)
c$$$            call flush(480+myrank)
c$$$            call flush(490+myrank)
c$$$            call flush(500+myrank)
c$$$            call flush(510+myrank)
c$$$            call flush(530+myrank)
c$$$            call flush(540+myrank)
c$$$            call flush(550+myrank)
c$$$            call flush(560+myrank)
c$$$            call flush(570+myrank)
c$$$            call flush(580+myrank)
c$$$            call flush(590+myrank)
      endif
 555  format(5(2x,e14.7))

c $$$$$$$$$$$$$$$$$$$$$$$$$$$
      tmp1 =  MINVAL(cdelsq)
      tmp2 =  MAXVAL(cdelsq)
      if(numpe>1) then
         call MPI_REDUCE (tmp1, tmp3, 1,MPI_DOUBLE_PRECISION,
     &        MPI_MIN, master, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
     &        MPI_MAX, master, MPI_COMM_WORLD, ierr)
         tmp1=tmp3
         tmp2=tmp4
      endif
      if (myrank .EQ. master) then !print CDelta^2 range
         write(34,*)lstep,tmp1,tmp2
         call flush(34)
      endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$

      if (myrank .eq. master) then
c         write(*,*)'xnut=',sum(cdelsq)/nshg
         write(*,*)'cdelsq=',cdelsq(1), cdelsq(2)
      endif

      cdeslq = zero ! Debugging
      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         
         ngauss = nint(lcsyst)

         call scatnu (mien(iblk)%p, strl(iel:inum,:), 
     &        mxmudmi(iblk)%p,cdelsq,shp(lcsyst,1:nshl,:))
      enddo
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$  tmp1 =  MINVAL(xmudmi)
c$$$  tmp2 =  MAXVAL(xmudmi)
c$$$  if(numpe>1) then
c$$$  call MPI_REDUCE (tmp1, tmp3, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MIN, master, MPI_COMM_WORLD, ierr)
c$$$  call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MAX, master, MPI_COMM_WORLD, ierr)
c$$$      tmp1=tmp3
c$$$  tmp2=tmp4
c$$$  endif
c$$$  if (myrank .EQ. master) then
c$$$  write(35,*) lstep,tmp1,tmp2
c$$$  call flush(35)
c$$$  endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$

c
c  if flag set, write a restart file with info (reuse xmij's memory)
c
      if(irs.eq.11) then
         lstep=999
         xmij(:,1)=xnum(:)
         xmij(:,2)=xden(:)
         xmij(:,3)=cdelsq(:)
         xmij(:,5)=xlij(:,4)    !leave M_{12} in 4 and put L_{12} here
         call restar('out ',xmij,xlij) !also dump all of L_{ij} in ac
         stop
      endif
c
c  local clipping moved to scatnu with the creation of mxmudmi pointers
c
c$$$      rmu=datmat(1,2,1)
c$$$      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
c$$$      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0
c      stop !uncomment to test dmod
c


c  write out the nodal values of xnut (estimate since we don't calc strain
c  there and must use the filtered strain).
c

c      if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
c
c  collect the average strain into xnude(2)
c
c         xnude(:,2) = zero
c         do i = 1,numnp
c            xnude(ifath(i),2) = xnude(ifath(i),2) + strnrm(i)
c         enddo

c
c        get the instantaneous cdelsq of the fathers
c
c         xnuder(:,1)=xnuder(:,1)/(xnuder(:,2)+1.0e-9)

c         if(numpe .gt. 1 .and. nsons(1).ne.1 )then
c            call drvAllreduce(xnude(:,2), xnuder(:,2),nfath)
c         else
c            xnuder=xnude
c         endif
c     
c          nut= cdelsq    * |S|
c 
c         xnutf=xnuder(:,1)*xnuder(:,2)/nsons
c
c  collect the x and y coords into xnude
c
c         xnude = zero
c         do i = 1,numnp
c            xnude(ifath(i),1) = xnude(ifath(i),1) + x(i,1)
c            xnude(ifath(i),2) = xnude(ifath(i),2) + x(i,2)
c         enddo

c         if(numpe .gt. 1 .and. nsons(1).ne.1 )then
c            call drvAllreduce(xnude, xnuder,2*nfath)
c            xnuder(:,1)=xnuder(:,1)/nsons(:)
c            xnuder(:,2)=xnuder(:,2)/nsons(:)
c         else
c            xnuder=xnude
c         endif
c
c  xnude is the sum of the sons for each father on this processor
c
c         if((myrank.eq.master)) then
c            do i=1,nfath      ! cdelsq   * |S|
c               write(444,*) xnuder(i,1),xnuder(i,2),xnutf(i)
c            enddo
c            call flush(444)
c         endif
c      endif

      return
      end







c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



      
      subroutine hfilter (y, x, ien, hfres, shgl, shp, Qwtf)

      include "common.h"

      dimension y(nshg,5),             hfres(nshg,11)
      dimension x(numnp,3),            xl(npro,nenl,3)
      dimension ien(npro,nshl),        yl(npro,nshl,nflow),
     &          fresl(npro,nshl,11),        WdetJ(npro),
     &          u1(npro),              u2(npro),
     &          u3(npro),              rho(npro),
     &          dxdxi(npro,nsd,nsd),   dxidx(npro,nsd,nsd),
     &          shgl(nsd,nshl,maxsh),       shg(npro,nshl,nsd),
     &          shp(nshl,maxsh),           Qwtf(ngaussf),
     &          fresli(npro,nshl,11)

      dimension tmp(npro)

      call local (y,      yl,     ien,    5,  'gather  ')
      call localx (x,      xl,     ien,    3,  'gather  ')
c
c
      if(matflg(1,1).eq.0) then ! compressible
      yl (:,:,1) = yl(:,:,1) / (Rgas * yl(:,:,5))  !get density
      else
      yl (:,:,1) = one
      endif
c
      fresl = zero
c
      do intp = 1, ngaussf

c  calculate the metrics
c
c
c.... --------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
         dxdxi = zero
c
         do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(1,n,intp)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(2,n,intp)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(3,n,intp)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(1,n,intp)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(2,n,intp)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(3,n,intp)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(1,n,intp)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(2,n,intp)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(3,n,intp)
         enddo
c     
c.... compute the inverse of deformation gradient
c
         dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)
     &        - dxdxi(:,3,2) * dxdxi(:,2,3)
         dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)
     &        - dxdxi(:,1,2) * dxdxi(:,3,3)
         dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)
     &        - dxdxi(:,1,3) * dxdxi(:,2,2)
         tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)
     &        + dxidx(:,1,2) * dxdxi(:,2,1)
     &        + dxidx(:,1,3) * dxdxi(:,3,1) )
         dxidx(:,1,1) = dxidx(:,1,1) * tmp
         dxidx(:,1,2) = dxidx(:,1,2) * tmp
         dxidx(:,1,3) = dxidx(:,1,3) * tmp
         dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)
     &        - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
         dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)
     &        - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
         dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)
     &        - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
         dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)
     &        - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
         dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)
     &        - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
         dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)
     &        - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c     
         wght = Qwtf(intp)
         WdetJ(:) = wght / tmp(:)
         
         rho = zero
         do i=1,nshl
            rho = rho+shp(i,intp)*yl(:,i,1) !density at qpt
         enddo
         
         rho = rho * WdetJ      !rho * WdetJ

         u1=zero
         u2=zero
         u3=zero
         do i=1,nshl  ! velocities at qpt
            u1 = u1 + shp(i,intp)*yl(:,i,2)
            u2 = u2 + shp(i,intp)*yl(:,i,3)
            u3 = u3 + shp(i,intp)*yl(:,i,4)
         enddo

         do i = 1, nshl


            fresli(:,i,1) = rho * u1 * shp(i,intp) !rho * u1 * WdetJ
            fresli(:,i,2) = rho * u2 * shp(i,intp) !rho * u2 * WdetJ
            fresli(:,i,3) = rho * u3 * shp(i,intp) !rho * u3 * WdetJ
            fresli(:,i,4) = rho * shp(i,intp) !rho * WdetJ
            fresli(:,i,5) = WdetJ*shp(i,intp) !Integral of filter kernel
c                                                over the element
            fresli(:,i,6) = rho * u1 * u1 * shp(i,intp) !rho * u1 * u1 * WdetJ
            fresli(:,i,7) = rho * u2 * u2 * shp(i,intp) !rho * u2 * u2 * WdetJ
            fresli(:,i,8) = rho * u3 * u3 * shp(i,intp) !rho * u3 * u3 * WdetJ
            fresli(:,i,9) = rho * u1 * u2 * shp(i,intp) !rho * u1 * u2 * WdetJ
            fresli(:,i,10)= rho * u1 * u3 * shp(i,intp) !rho * u1 * u3 * WdetJ
            fresli(:,i,11)= rho * u2 * u3 * shp(i,intp) !rho * u2 * u3 * WdetJ


            do j = 1, 11
               fresl(:,i,j) = fresl(:,i,j) + fresli(:,i,j)
            enddo

         enddo  !end loop over element nodes

      enddo !end of loop over integration points

c      do j = 1, n
c      do i = 1, nenl
c      do nel = 1, npro
c      hfres(ien(nel,i),j) = hfres(ien(nel,i),j) + fresl(nel,i,j) 
c      enddo
c      enddo
c      enddo

      call local (hfres, fresl, ien, 11, 'scatter ')

      return
      end







c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||








      subroutine twohfilter (y, x, strnrm, ien, fres, 
     &     hfres, shgl, shp, Qwtf)

      include "common.h"

      dimension y(nshg,ndof),            fres(nshg,33)
      dimension x(numnp,nsd),            xl(npro,nenl,nsd)
      dimension ien(npro,nshl),        yl(npro,nshl,nflow),
     &          fresl(npro,33),        WdetJ(npro),
     &          u1(npro),              u2(npro),
     &          u3(npro),              dxdxi(npro,nsd,nsd),
     &          strnrm(npro,ngauss),    dxidx(npro,nsd,nsd),
     &          shgl(nsd,nshl,maxsh),       shg(npro,nshl,nsd),
     &          shp(nshl,maxsh),
     &          fresli(npro,33),       Qwtf(ngaussf),
     &          hfres(nshg,11),        hfresl(npro,nshl,11)

      dimension tmp(npro)

      call local (y,      yl,     ien,    5,  'gather  ')
      call localx (x,      xl,     ien,    3,  'gather  ')
      call local (hfres,  hfresl, ien,   11,  'gather  ')
c
      if(matflg(1,1).eq.0) then ! compressible
      yl (:,:,1) = yl(:,:,1) / (Rgas * yl(:,:,5))  !get density
      else
      yl (:,:,1) = one
      endif
c
      fresl = zero

      do intp = 1, ngauss


c  calculate the metrics
c
c
c.... --------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
        dxdxi = zero
c
          do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(1,n,intp)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(2,n,intp)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(3,n,intp)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(1,n,intp)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(2,n,intp)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(3,n,intp)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(1,n,intp)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(2,n,intp)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(3,n,intp)
          enddo
c
c.... compute the inverse of deformation gradient
c
        dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)
     &                 - dxdxi(:,3,2) * dxdxi(:,2,3)
        dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)
     &                 - dxdxi(:,1,2) * dxdxi(:,3,3)
        dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)
     &                 - dxdxi(:,1,3) * dxdxi(:,2,2)
        tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)
     &                       + dxidx(:,1,2) * dxdxi(:,2,1)
     &                       + dxidx(:,1,3) * dxdxi(:,3,1) )
        dxidx(:,1,1) = dxidx(:,1,1) * tmp
        dxidx(:,1,2) = dxidx(:,1,2) * tmp
        dxidx(:,1,3) = dxidx(:,1,3) * tmp
        dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)
     &                - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
        dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)
     &                - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
        dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)
     &                - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
        dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)
     &                - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
        dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)
     &                - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
        dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)
     &                - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
        wght=Qwt(lcsyst,intp)  ! may be different now
c        wght=Qwtf(intp)
        WdetJ = wght / tmp
c
      fresli=zero
c
      do i=1,nshl
        fresli(:,22) = fresli(:,22)+shp(i,intp)*yl(:,i,1)     ! density at qpt
        fresli(:,25) = fresli(:,25)+shp(i,intp)*hfresl(:,i,1) !bar(rho u1)
        fresli(:,26) = fresli(:,26)+shp(i,intp)*hfresl(:,i,2) !bar(rho u2)
        fresli(:,27) = fresli(:,27)+shp(i,intp)*hfresl(:,i,3) !bar(rho u3)
      enddo
c
c...fresli(:,28) = WdetJ * bar(rho u1) * bar(rho u1) / rho
c...fresli(:,29) = WdetJ * bar(rho u2) * bar(rho u2) / rho
c...fresli(:,30) = WdetJ * bar(rho u3) * bar(rho u3) / rho
c...fresli(:,31) = WdetJ * bar(rho u1) * bar(rho u2) / rho
c...fresli(:,32) = WdetJ * bar(rho u1) * bar(rho u3) / rho
c...fresli(:,33) = WdetJ * bar(rho u2) * bar(rho u3) / rho

      fresli(:,28) = WdetJ(:) * fresli(:,25) * 
     &     fresli(:,25) / fresli(:,22)
      fresli(:,29) = WdetJ(:) * fresli(:,26) * 
     &     fresli(:,26) / fresli(:,22) 
      fresli(:,30) = WdetJ(:) * fresli(:,27) * 
     &     fresli(:,27) / fresli(:,22)
      fresli(:,31) = WdetJ(:) * fresli(:,25) * 
     &     fresli(:,26) / fresli(:,22)       
      fresli(:,32) = WdetJ(:) * fresli(:,25) * 
     &     fresli(:,27) / fresli(:,22)
      fresli(:,33) = WdetJ(:) * fresli(:,26) * 
     &     fresli(:,27) / fresli(:,22) 

      fresli(:,25) = fresli(:,25) * WdetJ(:) ! WdetJ*bar(rho u1)
      fresli(:,26) = fresli(:,26) * WdetJ(:) ! WdetJ*bar(rho u2)
      fresli(:,27) = fresli(:,27) * WdetJ(:) ! WdetJ*bar(rho u3)
c
      do n = 1,nshl
        shg(:,n,1) = (shgl(1,n,intp) * dxidx(:,1,1)
     &              + shgl(2,n,intp) * dxidx(:,2,1)
     &              + shgl(3,n,intp) * dxidx(:,3,1))
        shg(:,n,2) = (shgl(1,n,intp) * dxidx(:,1,2)
     &              + shgl(2,n,intp) * dxidx(:,2,2)
     &              + shgl(3,n,intp) * dxidx(:,3,2))
        shg(:,n,3) = (shgl(1,n,intp) * dxidx(:,1,3)
     &              + shgl(2,n,intp) * dxidx(:,2,3)
     &              + shgl(3,n,intp) * dxidx(:,3,3))
      enddo

      do j=10,12  ! normal strainrate u_{i,i} no sum on i
       ig=j-9
       iv=j-8
       do i=1,nshl
        fresli(:,j) = fresli(:,j)+shg(:,i,ig)*yl(:,i,iv)
       enddo
      enddo

c shear stresses  NOTE  there may be faster ways to do this
c                  check agains CM5 code for speed WTP
       
       do i=1,nshl
        fresli(:,13) = fresli(:,13)+shg(:,i,2)*yl(:,i,2)
     &                             +shg(:,i,1)*yl(:,i,3)
        fresli(:,14) = fresli(:,14)+shg(:,i,3)*yl(:,i,2)
     &                             +shg(:,i,1)*yl(:,i,4)
        fresli(:,15) = fresli(:,15)+shg(:,i,3)*yl(:,i,3)
     &                             +shg(:,i,2)*yl(:,i,4)
       enddo

      fresli(:,13) = pt5 * fresli(:,13)
      fresli(:,14) = pt5 * fresli(:,14)
      fresli(:,15) = pt5 * fresli(:,15)

      strnrm(:,intp) = fresli(:,22) * sqrt(
     &   two * (fresli(:,10)**2 + fresli(:,11)**2 + fresli(:,12)**2)
     &  + four * ( fresli(:,13)**2 + fresli(:,14)**2 + 
     &    fresli(:,15)**2 ) )

c
c S_ij
c

      fresli(:,10) = fresli(:,10) * WdetJ ! u_{1,1}*WdetJ
      fresli(:,11) = fresli(:,11) * WdetJ ! u_{2,2}*WdetJ
      fresli(:,12) = fresli(:,12) * WdetJ ! u_{3,3}*WdetJ
      fresli(:,13) = fresli(:,13) * WdetJ ! (1/2)*(u_{1,2}+u_{2,1})*WdetJ
      fresli(:,14) = fresli(:,14) * WdetJ ! (1/2)*(u_{1,3}+u_{3,1})*WdetJ
      fresli(:,15) = fresli(:,15) * WdetJ ! (1/2)*(u_{2,3}+u_{3,2})*WdetJ

      fresli(:,22) = fresli(:,22) * WdetJ   !rho * WdetJ
c     fresli(:,24) = fresli(:,24) * WdetJ
     
      u1=zero
      u2=zero
      u3=zero
      do i=1,nshl
       u1 = u1 + shp(i,intp)*yl(:,i,2)
       u2 = u2 + shp(i,intp)*yl(:,i,3)
       u3 = u3 + shp(i,intp)*yl(:,i,4)
      enddo

      fresli(:,1) = fresli(:,22) * u1   !rho u1 * WdetJ
      fresli(:,2) = fresli(:,22) * u2   !rho u2 * WdetJ
      fresli(:,3) = fresli(:,22) * u3   !rho u3 * WdetJ

      fresli(:,4) = fresli(:,1) * u1    !rho u1 u1 *WdetJ
      fresli(:,5) = fresli(:,2) * u2    !rho u2 u2 *WdetJ
      fresli(:,6) = fresli(:,3) * u3    !rho u3 u3 *WdetJ
      fresli(:,7) = fresli(:,1) * u2    !rho u1 u2 *WdetJ
      fresli(:,8) = fresli(:,1) * u3    !rho u1 u3 *WdetJ
      fresli(:,9) = fresli(:,2) * u3    !rho u2 u3 *WdetJ

      fresli(:,16) = strnrm(:,intp) * fresli(:,10) ! rho *|Eps| *Eps11 *WdetJ
      fresli(:,17) = strnrm(:,intp) * fresli(:,11) ! rho *|Eps| *Eps22 *WdetJ
      fresli(:,18) = strnrm(:,intp) * fresli(:,12) ! rho *|Eps| *Eps33 *WdetJ
      fresli(:,19) = strnrm(:,intp) * fresli(:,13) ! rho *|Eps| *Eps12 *WdetJ
      fresli(:,20) = strnrm(:,intp) * fresli(:,14) ! rho *|Eps| *Eps13 *WdetJ
      fresli(:,21) = strnrm(:,intp) * fresli(:,15) ! rho *|Eps| *Eps23 *WdetJ

      fresli(:,23) = WdetJ   !    Integral of 1 over the element
c
      do i = 1, 33
         fresl(:,i) = fresl(:,i) + fresli(:,i)
      enddo
   
      enddo !end of loop over integration points
c
      do j = 1,nshl
      do nel = 1,npro
        fres(ien(nel,j),:) = fres(ien(nel,j),:) + fresl(nel,:) 
      enddo
      enddo

      return
      end








c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||







      subroutine disPtion (y, x, ien, shgl, shp, eps, xmudmi) 

      use rlssave  ! Use the resolved Leonard stresses at the nodes.

      include "common.h"

      dimension xmudmi(npro,ngauss),         y(nshg,ndof),
     &          x(numnp,nsd),               ien(npro,nshl),
     &          shg(npro,nshl,nsd),
     &          shgl(nsd,nshl,maxsh),       shp(nshl,maxsh),
     &          dxdxi(npro,nsd,nsd),        dxidx(npro,nsd,nsd),
     &          WdetJ(npro),
     &          eps(3),                     fresli(npro,33),
     &          epsli(npro,3),              epsl(npro,3)

      dimension yl(npro,nshl,nflow),         xl(npro,nenl,nsd), 
     &          strnrm(npro,ngauss),         rlsl(npro,nshl,6)          

      dimension tmp(npro)

      call local (y,      yl,     ien,    5,  'gather  ')
      call localx (x,      xl,     ien,    3,  'gather  ')
      call local (rls,    rlsl,   ien,    6,  'gather  ') 

      epsl = zero

      yl (:,:,1) = one ! Unit density

      do intp = 1, ngauss

      fresli=zero

c  calculate the metrics
c
c
c.... --------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
        dxdxi = zero
c
          do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(1,n,intp)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(2,n,intp)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(3,n,intp)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(1,n,intp)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(2,n,intp)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(3,n,intp)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(1,n,intp)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(2,n,intp)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(3,n,intp)
          enddo
c
c.... compute the inverse of deformation gradient
c
        dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)
     &                 - dxdxi(:,3,2) * dxdxi(:,2,3)
        dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)
     &                 - dxdxi(:,1,2) * dxdxi(:,3,3)
        dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)
     &                 - dxdxi(:,1,3) * dxdxi(:,2,2)
        tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)
     &                       + dxidx(:,1,2) * dxdxi(:,2,1)
     &                       + dxidx(:,1,3) * dxdxi(:,3,1) )
        dxidx(:,1,1) = dxidx(:,1,1) * tmp
        dxidx(:,1,2) = dxidx(:,1,2) * tmp
        dxidx(:,1,3) = dxidx(:,1,3) * tmp
        dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)
     &                - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
        dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)
     &                - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
        dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)
     &                - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
        dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)
     &                - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
        dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)
     &                - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
        dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)
     &                - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
        wght=Qwt(lcsyst,intp)  
        WdetJ = wght / tmp      

c
c
      do n = 1,nshl
        shg(:,n,1) = (shgl(1,n,intp) * dxidx(:,1,1)
     &              + shgl(2,n,intp) * dxidx(:,2,1)
     &              + shgl(3,n,intp) * dxidx(:,3,1))
        shg(:,n,2) = (shgl(1,n,intp) * dxidx(:,1,2)
     &              + shgl(2,n,intp) * dxidx(:,2,2)
     &              + shgl(3,n,intp) * dxidx(:,3,2))
        shg(:,n,3) = (shgl(1,n,intp) * dxidx(:,1,3)
     &              + shgl(2,n,intp) * dxidx(:,2,3)
     &              + shgl(3,n,intp) * dxidx(:,3,3))
      enddo


c
c
      do i=1,nshl
        fresli(:,22) = fresli(:,22)+shp(i,intp)    ! unit density at qpt
      enddo

c
c
      do j=10,12  ! normal strainrate u_{i,i} no sum on i
       ig=j-9
       iv=j-8
       do i=1,nshl
        fresli(:,j) = fresli(:,j)+shg(:,i,ig)*yl(:,i,iv)
       enddo
      enddo

c shear stresses  NOTE  there may be faster ways to do this
c                  check agains CM5 code for speed WTP
       
       do i=1,nshl
        fresli(:,13) = fresli(:,13)+shg(:,i,2)*yl(:,i,2)
     &                             +shg(:,i,1)*yl(:,i,3)
        fresli(:,14) = fresli(:,14)+shg(:,i,3)*yl(:,i,2)
     &                             +shg(:,i,1)*yl(:,i,4)
        fresli(:,15) = fresli(:,15)+shg(:,i,3)*yl(:,i,3)
     &                             +shg(:,i,2)*yl(:,i,4)
       enddo

      fresli(:,13) = pt5 * fresli(:,13)
      fresli(:,14) = pt5 * fresli(:,14)
      fresli(:,15) = pt5 * fresli(:,15)

      strnrm(:,intp) = fresli(:,22) * (
     &   two * (fresli(:,10)**2 + fresli(:,11)**2 + fresli(:,12)**2)
     &  + four * ( fresli(:,13)**2 + fresli(:,14)**2 + 
     &    fresli(:,15)**2 ) )

c
c

      fresli(:,10) = fresli(:,10) * WdetJ ! u_{1,1}*WdetJ
      fresli(:,11) = fresli(:,11) * WdetJ ! u_{2,2}*WdetJ
      fresli(:,12) = fresli(:,12) * WdetJ ! u_{3,3}*WdetJ
      fresli(:,13) = fresli(:,13) * WdetJ ! (1/2)*(u_{1,2}+u_{2,1})*WdetJ
      fresli(:,14) = fresli(:,14) * WdetJ ! (1/2)*(u_{1,3}+u_{3,1})*WdetJ
      fresli(:,15) = fresli(:,15) * WdetJ ! (1/2)*(u_{2,3}+u_{3,2})*WdetJ
c
      strnrm(:,intp) =  strnrm(:,intp) * WdetJ !  ( |Eps|^2 )*WdetJ   

      epsli(:,1) = -xmudmi(:,intp)*strnrm(:,intp) 

      epsli(:,2) = rlsl(:,intp,1)*fresli(:,10) +
     &             rlsl(:,intp,2)*fresli(:,11) +     
     &             rlsl(:,intp,3)*fresli(:,12) +
     &             two*( rlsl(:,intp,4)*fresli(:,13)+
     &                   rlsl(:,intp,5)*fresli(:,14)+
     &                   rlsl(:,intp,6)*fresli(:,15) )

      epsl(:,1) = epsl(:,1) + epsli(:,1)
      epsl(:,2) = epsl(:,2) + epsli(:,2)

      enddo  ! end loop over integ. pts

      do i = 1, npro
         eps(1) = eps(1) + epsl(i,1)
         eps(2) = eps(2) + epsl(i,2)
      enddo

      return
      end      



