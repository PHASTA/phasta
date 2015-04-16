      subroutine sonfath (ncorp, vnumnpl, ifath, nsons, tfath, 
     &                    numpe, numnp,   maxnode)

c     This program writes to geom.dat.<1..numpe> the number sons and
c     fathers required for averaging in the z-direction (Turbulent 
c     Boundary Layer), (x,y)-space (Turbulent Channel Flow), 
c     or (x,y,z)-space (Isotropic Turbulence)

      integer inform,n1m,n2m,n3m,numpe,numpe1,numnp,maxnode, tfath,
     &     counter,nc8,nc4,node,ifil,itmp
      integer gfath(0:numnp),ifath(numpe,maxnode),ncorp(numpe,maxnode),
     &     nsons(tfath),dummy(numnp),iper(numnp), vnumnpl(numpe)
      integer, allocatable :: imap(:),invmap(:)
      allocate (imap(numnp))
      allocate (invmap(numnp))

      write(*,*)'Enter number of nodes in each direction'
      read(*,*)n1m,n2m,n3m

      write(*,*)'For Z averaging enter 1 For X-Y averaging enter 2'
      write(*,*)'For X-Y-Z averaging enter 3  W=0 SPEBC enter 4'
      write(*,*)'fully local model 5'
      read(*,*)inform

c
c.... read and generate map
c      
      kmap = 15

      open(kmap,file='geom.kmap',status='old',err=183)

      do i=1,numnp
         read(kmap,*,err=183) iorig
         imap(i)=iorig+1  ! imap takes arguement of original and returns new
         invmap(iorig+1)=i ! invmap takes arguemnt of new and returns original
      enddo 
      goto 184
 183  write(*,*) 'geom.kmap not read properly...assuming inform=5'
      inform=5
 184  continue

      nsoncor=1
      if(inform.eq.4) then
         inform=1
         nsoncor=0
      endif

      gfath(:) = 0

      nc4=n3m*(n2m-1)+1
      nc8=n3m*n2m*n1m+1-n3m      
      
      if (numpe.eq.1)then
         numnp   = n1m*n2m*n3m
         maxnode = numnp
      endif

      if (inform.eq.1) then ! For averaging over Z - lines
         
         do i = 1, tfath
            nsons(i) = n3m-nsoncor
         enddo
         
         counter = 1
         do i = 1, nc8, n3m
            do j = i, i+n3m-1
               gfath(j) = counter
            enddo
            counter = counter + 1
         enddo
         
         if(numpe.gt.1)then
            do j = 1, numpe
               do i = 1, maxnode
                  ifath(j,i) = gfath(imap(ncorp(j,i)))
               enddo
            enddo
         else
            do i = 1, maxnode
               ifath(1,i) = gfath(imap(i))
            enddo
         endif

      endif

      if (inform.eq.2) then ! For averaging over X-Z planes
         
         do i = 1, tfath
            nsons(i) = n1m*n3m - (n1m+n3m-1)
         enddo

         counter = 1
         do j = 1, nc4, n3m
            do i = 1, n1m
               do k = j, j+n3m-1
                  node = (i-1)*n2m*n3m + k
                  gfath(node) = counter
               enddo
            enddo
            counter = counter + 1
         enddo

         if(numpe.gt.1)then
            do j = 1, numpe
               do i = 1, maxnode
                  ifath(j,i) = gfath(imap(ncorp(j,i)))
               enddo
            enddo
         else
            do i = 1, maxnode
               ifath(1,i) = gfath(imap(i))
            enddo
         endif
         
      endif

      if( inform.eq.3 )then !For averaging over X-Y-Z   

         open (unit=96,FILE='geom.iper',form='unformatted',
     &        status='old')
         read(96) (iper(l),l=1,numnp)
         close(96)

         dummy(:) = 1

         do i = 1, numnp
            if(iper(i).eq.0)then
               iper(i)=i
            endif
         enddo
         do j = 1, numnp
            i = iper(j)
            if (i .ne. j) then
               dummy(j) = 0
            endif
         enddo
         itmp=sum(dummy)
         nsons(tfath) = itmp
         write(*,*)'nsons(1)=',nsons(1)
         ifath=1
      endif
         
      if(inform.eq.5) then      ! local clipping

         do i = 1, numpe
            ifil = 99
            
            nsons=1
            do j=1,vnumnpl(i)
               ifath(i,j)=j
            enddo

         enddo
      endif

      return
      end
