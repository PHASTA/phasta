      subroutine genscale(y, x, iBC)
c
c----------------------------------------------------------------------
c This subroutine calculate the y^+ and eta at inflow and internal face.
c From these generate the scaling for the inflow data.
c
c input:
c  iBC    (numnp)               : boundary condition code
c  x      (numnp,nsd)           : node coordinates
c
c output:
c  y      (numnp,ndof)          : initial values of Y variables
c
c
c Elaine Bohr december 2001
c----------------------------------------------------------------------
c
       use spebc
c       use pointer_data
       include "common.h"
       include "mpif.h"
       include "auxmpi.h"
c
       dimension y(numnp,ndof),   iBC(numnp),  
     &           x(numnp,nsd), velbarR(nfint,nflow)
       dimension ifath(numnp),  velbarl(nelint,nshl,nflow),
     &           v1(nfint),       ymapped(numnp,ndof),
     &		 shapef(nshl),	shgradl(nshl,nsd),
     &		 xsi(nsd), yintl(nelint,nshl,nflow),
     &		 flucl(nelint,nshl,nflow),
     &		 ubarintl(nelint,nshl,nflow),
     &		 fluc1(npin,nflow), fluc2(npin,nflow),
     &		 ubar1(npin,nflow), ubar2(npin,nflow)
       integer   element, dir

       real*8    ymax, displTi, displTr, correction
       real*8	 freestream(nflow)
       save deltaint

       
c	return
        ymapped(:,2:4)=y(:,1:3)
	ymapped(:,1)=y(:,4)
	ymapped(:,5)=y(:,5)
	
	ubar2 = 0
	fluc2 = 0
	
	ymax = xyn(nfint)

c
c .... Localizing the solution vector on virtual plane
c

        do i = 1, nelint
          do j = 1, 3
            yintl(i,:,j+1) = y(ien2D(i,:),j)
          enddo
          yintl(i,:,1) = y(ien2D(i,:),4)
	  if(nflow.gt.4) then
            do j = 5, nflow
              yintl(i,:,j) = y(ien2D(i,:),j)
            enddo
          endif
	enddo  

c
c .... Finding averaged velocity in spanwise direction
c      for the virtual plane
c

	do i=1,nfint
	  velbarR(i,:)=0
	  do j=1,imax(i)+1
	    call shptet(ipord,xsinfin(i,j,:),shapef(:),shgradl(:,:))
	    do k=1,nshl
	      velbarR(i,:)=velbarR(i,:) 
     &		+ yintl(elcnfin(i,j),k,:)*shapef(k)
	    enddo
	  enddo
	  velbarR(i,:)=velbarR(i,:) / (imax(i)+1)
	enddo
	 
c
c .... Label the nodes that near the BL thickness
c 

       if (thetag.eq.0.0) then
         dir = 2
       else
         dir = 4
       endif

       v1(1)=10.0
       do i=2,nfint+1
          v1(i)=velbarR(i-1,dir)-0.99*vel
          if((v1(i).gt.0).and.(v1(i-1).le.0)) then
             label=i-1
             go to 200
          endif
       enddo
       label=i-1

c     
c.... Find the BL thickness by means of finding the y coord
c     

 200   continue
       dv=velbarR(label,dir)-velbarR(label-1,dir)
       dy=xyn(label)-xyn(label-1)

c     
c .... Current calculation of bl thickness at recycle plane
c

       if(istep.ne.0) then
          dlast=deltaint
          deltaint=xyn(label-1)
     &		+ dy*(0.99*vel-velbarR(label-1,dir))/dv
     
c
c .... Early transients cause jumpy delta, smooth it.
c

          deltaint=min(1.05*dlast,max(deltaint,0.95*dlast))
       else
          deltaint=xyn(label-1)
     &		+ dy*(0.99*vel-velbarR(label-1,dir))/dv
       endif

c
c .... Deltaint is now the ratio of BL thickness at the interior plane
c      to the BL thickness at the inlet plane
c

       deltaint=min(two*rbltin,max(deltaint,pt5*rbltin)) 
       rdelta = deltaint/rbltin
       
c
c .... Finding freestream solutions
c
	
	freestream = 0
	icount = 0
	do i=1, nfint
	  if (xyn(i).ge.deltaint) then
	    freestream(:) = freestream(:) + velbarR(i,:)
	    icount = icount + 1 
	  endif
	enddo
	freestream = freestream / icount
	
c
c .... Putting the freestream values into the average outside the BLT
c

	do i=1, nfint
	  if (xyn(i).ge.deltaint) then
	    velbarR(i,:) = freestream(:)
	  endif
	enddo

c
c .... Localizing the averaged velocity found above
c

	do i=1,nelint
	  do k=1,nshl
	    do j=1,nfint-1
	      if (thetag.eq.0.0) then
	        if ((x(ien2D(i,k),2).ge.xyn(j)) .and.
     &		    (x(ien2D(i,k),2).le.(xyn(j+1)+0.000001))) then
		  tmp = (x(ien2D(i,k),2) - xyn(j)) /
     &			(xyn(j+1) - xyn(j))
     		  do l=1,nflow
		     velbarl(i,k,l) = 
     &		            (velbarR(j+1,l) - velbarR(j,l)) * 
     &                       tmp + velbarR(j,l)
	          enddo
		endif
	      else
	        if ((xcyl(ien2D(i,k),1).ge.xcyl(nrint(j+1),1)) .and.
     &		    (xcyl(ien2D(i,k),1).le.xcyl(nrint(j),1))) then
                  tmp = (xcyl(ien2D(i,k),1) - xcyl(nrint(j+1),1)) / 
     &	              (xcyl(nrint(j),1) - xcyl(nrint(j+1),1))
     		  do l=1,nflow
		     velbarl(i,k,l) = 
     &		            (velbarR(j,l) - velbarR(j+1,l)) * 
     &                       tmp + velbarR(j+1,l)
	          enddo
		endif
     	      endif
	    enddo
	  enddo
	enddo
	
c
c --- For now only Blasius is coded ---
c
       
c
c .... Calculate fluctuations on elements of internal plane
c

c       flucl = yintl - velbarl

c
c .... Calculate mean values on elements of internal plane
c

       ubarintl = velbarl 

c
c .... Calculating the coordinates of the point from where the
c      solution will be projected to the inlet plane
c


	do i=1,npin
	
c
c .... Cartesian coodinate system
c

	  if (thetag.eq.0.0) then
	    xts1 = x(nen1(i),1) + plandist
	    if (xynin(i)*rdelta.gt.ymax) then
	      xts2 = ymax
	    else  
	      xts2 = xynin(i)*rdelta
	    endif
	    xts3 = x(nen1(i),3)

c
c .... Cylindrical coordinate system
c

	  else
	    if (xynin(i).le.0.00001) then
	      xts1 = (radcyl-xynin(i)*rdelta*sang-tolerence)
     & 		*cos(xcyl(nen1(i),2))
	      xts2 = (radcyl-xynin(i)*rdelta*sang-tolerence)
     &		*sin(xcyl(nen1(i),2))
	      xts3 = (aR-(radcyl-xynin(i)*rdelta*sang-tolerence)
     &	       * (xnrml*cos(xcyl(nen1(i),2))
     &	       +  ynrml*sin(xcyl(nen1(i),2))))/znrml
            elseif (xynin(i)*rdelta.gt.ymax) then
	      xts1 = (radcyl-ymax*sang)
     & 		*cos(xcyl(nen1(i),2))
	      xts2 = (radcyl-ymax*sang)
     &		*sin(xcyl(nen1(i),2))
	      xts3 = (aR-(radcyl-ymax*sang)
     &	       * (xnrml*cos(xcyl(nen1(i),2))
     &	       +  ynrml*sin(xcyl(nen1(i),2))))/znrml
	    else
	      xts1 = (radcyl-xynin(i)*rdelta*sang)
     & 		*cos(xcyl(nen1(i),2))
	      xts2 = (radcyl-xynin(i)*rdelta*sang)
     &		*sin(xcyl(nen1(i),2))
	      xts3 = (aR-(radcyl-xynin(i)*rdelta*sang)
     &	       * (xnrml*cos(xcyl(nen1(i),2))
     &	       +  ynrml*sin(xcyl(nen1(i),2))))/znrml
            endif
	  endif
	  
c
c .... Searching for the appropriate element
c

     	  call elem_search(xintl, xts1, xts2, xts3,
     &			   xsi(:), element, 2)
      	  call shptet(ipord,xsi(:),shapef(:),shgradl(:,:))	  

c
c .... Calculating the average velocity and fluctuations
c      for the inlet plane
c

	  do k=1,nshl
	    fluc2(i,:)= 0 !fluc2(i,:) + flucl(element,k,:)*shapef(k)
	    ubar2(i,:)=ubar2(i,:) + ubarintl(element,k,:)*shapef(k)
	  enddo
	enddo  
	        

c$$$c
c$$$c keep freestream values set through averages
c$$$c
c$$$         ubaro=0
c$$$	 tbaro=0
c$$$         icount=0
c$$$         do i=1,nfin
c$$$            if(yin(i).ge.rbltin) then
c$$$               nzl=nsons(i)  !Elaine
c$$$               nzb=ienson1(i,1)
c$$$               nze=nzb+nzl-1
c$$$	       tbaro=tbaro+2.0*ubar2(i,5)+sum(fluc2(nzb:nze,5))
c$$$               ubaro=ubaro               +sum(fluc2(nzb:nze,2))
c$$$               icount=icount+nzl
c$$$            endif
c$$$         enddo
c$$$         
c$$$c     alternative to myway
c$$$c     
c$$$         ubaro=ubaro/icount
c$$$	 tmeaninflow=0.0625097048890964
c$$$	 fact= tmeaninflow/(tbaro/icount)
c$$$	 if (fact.ge. 0.9999999 .and. fact.le.1.0000001) fact = 1.0
c$$$         
c$$$         do i=1,nfin
c$$$            if(yin(i).ge.rbltin) then
c$$$               ubar2(i,2)=1.0-ubaro
c$$$            endif
c$$$         enddo
         fact = 1.0 
	 rvscal = 1.0
	 
c
c .... Putting the freestream value outside the BLT into ubar2
c

	do i = 1, npin
	  if (xynin(i).ge.rbltin) then
	    ubar2(i,:) = freestream(:)
c	    ubar2(i,dir) = 1.0
	    fluc2(i,:) = 0
	  endif
	enddo

c$$$c
c$$$c .... For the cylindrical case the freestream velocity needs
c$$$c      to be corrected for the blockage
c$$$c
c$$$
c$$$	if (thetag.ne.0.0) then
c$$$	  displTi = 0.0
c$$$	  displTr = 0.0
c$$$	  do i = 2, nfint
c$$$
c$$$c
c$$$c .... Displacement thickness for inlet plane
c$$$c
c$$$
c$$$	    displTi = displTi + (1 - y(nrint(i),3))
c$$$     &              * (xyn(i) - xyn(i-1)) * (radcyl - xyn(i))
c$$$
c$$$c
c$$$c .... Displacement thickness for recycle plane
c$$$c
c$$$
c$$$     	    displTr = displTr + (1 - velbarR(i,4))
c$$$     &              * (xyn(i) - xyn(i-1)) * (radcyl - xyn(i))
c$$$	  enddo
c$$$c	  displTi = radcyl - sqrt(radcyl*radcyl - displTi)
c$$$c	  displTr = radcyl - sqrt(radcyl*radcyl - displTr)
c$$$	  correction = (radcyl*radcyl - displTr)
c$$$     &               / (radcyl*radcyl - displTi)
c$$$        else
	  correction = 1.0
c$$$	endif

c     
c .... Scaled plane extraction boundary condition
c     

         ymapped(nen1(1:npin),1)= correction * (ubar2(:,1)+fluc2(:,1))
         ymapped(nen1(1:npin),2)= correction *
     &        (ubar2(:,2)+fluc2(:,2)) !myway *factu
         ymapped(nen1(1:npin),3)= correction *
     &        (ubar2(:,3)+fluc2(:,3))*rvscal
         ymapped(nen1(1:npin),4)= correction * (ubar2(:,4)+fluc2(:,4)) 
c     &			freestream(4)
	 ymapped(nen1(1:npin),5)= correction * fact*(ubar2(:,5)+fluc2(:,5))



c     
c .... Ready to put the solution on the inflow plane 
c

      if(intpres.eq.1) then     !interpolating pressure at inflow
         where (btest(iBC,11))
            y(:,1) = ymapped(:,2)
            y(:,2) = ymapped(:,3)
            y(:,3) = ymapped(:,4)
	    y(:,4) = ymapped(:,1)
	    y(:,5) = ymapped(:,5)
         endwhere
      else                      ! not interpolating pressure at inflow
         where (btest(iBC,11))
            y(:,1) = ymapped(:,2)
            y(:,2) = ymapped(:,3)
            y(:,3) = ymapped(:,4)
	    y(:,5) = ymapped(:,5)
         endwhere
      endif

c     
c     debugging variables
c     
c      if(iter.eq.nitr) then
c         write(555,556)lstep+1,deltaint,label,nfint
c         write(554,556)lstep+1,yplusi(2),ypluso(2),factu,factt,gamt
c         call flush(554)
c         call flush(555)
c      endif
c 556  format(i6,5(2x,e14.7))
c
c.... return
c

      return
      end
