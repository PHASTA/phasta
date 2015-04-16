      module quadfilt

      real*8, allocatable :: shpf(:,:,:)
      real*8, allocatable :: shglf(:,:,:,:)
      real*8, allocatable :: Qptf(:,:,:)
      real*8, allocatable :: Qwtf(:,:)

      real*8, allocatable :: numNden(:,:)

      real*8, allocatable :: xnd(:,:)
      real*8, allocatable :: xmodcomp(:,:)
      real*8, allocatable :: xmcomp(:,:)
      real*8, allocatable :: xlcomp(:,:)
      real*8, allocatable :: xl1comp(:,:)
      real*8, allocatable :: xl2comp(:,:)
      real*8, allocatable :: ucomp(:,:)
      real*8, allocatable :: scomp(:)

      end module

c------------------------------------------------------------------------------

      subroutine setfilt

      use quadfilt

      include "common.h"

      ifiltrl = mod(iLES,10)

         if (ifiltrl .eq. 1) then
            nintf(1) = 1       ! tets
            nintf(2) = 1       ! hexes
            nintf(3) = 1       ! wedges
            nintf(5) = 1       ! pyramids
         else if (ifiltrl .eq. 2) then
            nintf(1) = 4
            nintf(2) = 8
            nintf(3) = 6
            nintf(5) = 8
         else if (ifiltrl .eq. 3) then
            nintf(1) = 16
            nintf(2) = 27
            nintf(3) = 18
            nintf(5) = 27
         else if (ifiltrl .eq. 4) then
            nintf(1) = 29
            nintf(2) = 64
            nintf(3) = 48
            nintf(5) = 64
         endif


      allocate ( shpf(MAXTOP,maxsh,maxval(nintf)) )
      allocate ( shglf(MAXTOP,nsd,maxsh,maxval(nintf)) )
      allocate ( Qptf(MAXTOP,4,maxval(nintf)) )
      allocate ( Qwtf(MAXTOP,maxval(nintf)) )


      allocate ( numNden(nshg,2) )
C
C In development
C
      allocate ( xnd(70,2) )
      allocate ( xmodcomp(70,5) )
      allocate ( xmcomp(70,6) )
      allocate ( xlcomp(70,6) )
      allocate ( xl1comp(70,6) )
      allocate ( xl2comp(70,6) )
      allocate ( ucomp(70,3) )
      allocate ( scomp(70) )
C
C In development
C

      numNden = zero

c.... return

      return
      end

c------------------------------------------------------------------------------

      subroutine filtprep

      use quadfilt

      include "common.h"
 
      real*8, allocatable :: tmpQptf (:,:), tmpQwtf (:)

      allocate (tmpQptf(4,nint(1)))
      allocate (tmpQwtf(nintf(1)))
      call symtet(nintf(1),tmpQptf,tmpQwtf,nerr)
      Qptf(1,1:4,1:nintf(1)) = tmpQptf(1:4,1:nintf(1))
      Qwtf(1,1:nintf(1))     = tmpQwtf(1:nintf(1))
c
c.... adjust quadrature weights to be consistent with the
c     design of tau. 
c
      do i = 1, nintf(1)
         Qwtf(1,i) = (four/three)*Qwtf(1,i)
      enddo
      
      deallocate (tmpQptf)
      deallocate (tmpQwtf)
         

      allocate (tmpQptf(4,nint(2)))
      allocate (tmpQwtf(nintf(2)))
      call symhex(nintf(2),tmpQptf,tmpQwtf,nerr)
      Qptf(2,1:4,1:nintf(2)) = tmpQptf(1:4,1:nintf(2))
      Qwtf(2,1:nintf(2)) = tmpQwtf(1:nintf(2))
      deallocate (tmpQptf)
      deallocate (tmpQwtf)
      
      allocate (tmpQptf(4,nintf(3)))
      allocate (tmpQwtf(nintf(3)))
      call symwdg(nintf(3),tmpQptf,tmpQwtf,nerr)
      Qptf(3,1:4,1:nintf(3)) = tmpQptf(1:4,1:nintf(3))
      Qwtf(3,1:nintf(3)) = tmpQwtf(1:nintf(3))
      deallocate (tmpQptf)
      deallocate (tmpQwtf)

      allocate (tmpQptf(4,nintf(5)))
      allocate (tmpQwtf(nintf(5)))
      call sympyr(nintf(5),tmpQptf,tmpQwtf,nerr)
      Qptf(5,1:4,1:nintf(5)) = tmpQptf(1:4,1:nintf(5))
      Qwtf(5,1:nintf(5)) = tmpQwtf(1:nintf(5))
      deallocate (tmpQptf)
      deallocate (tmpQwtf)

c
c.... loop through element blocks
c
      do iblk = 1, nelblk
c
c.... get coord. system and element type 
c
         lcsyst = lcblk(3,iblk)
         nshl   = lcblk(10,iblk)
c
c.... generate the shape-functions in local coordinates
c
         if (lcsyst .eq. 1) then ! tets
           
            do i=1,nintf(1)  
               call shpTet(ipord,Qptf(1,1:3,i),shpf(1,:,i),
     &              shglf(1,:,:,i))
            enddo
c
c.... permute to positive volume element
c
            shglf(1,:,1:nshl,1:nintf(1)) = 
     &           shglf(1,:,1:nshl,1:nintf(1))/two

c     
         else if (lcsyst .eq. 2) then ! hexes
c     


            do i=1,nintf(2)
               call shphex  (ipord, Qptf(2,1:3,i),shpf(2,:,i),
     &              shglf(2,:,:,i))
             
            enddo
c
         else if (lcsyst .eq. 3) then ! wedges
c
            do i=1,nintf(3)
               call shp6W  (ipord,Qptf(3,1:3,i),shpf(3,:,i),
     &              shglf(3,:,:,i))
            enddo
c     
         else if (lcsyst .eq. 5) then ! pyramids
c     
            do i=1,nintf(5)
               call shppyr  (ipord,Qptf(5,1:3,i),shpf(5,:,i),
     &              shglf(5,:,:,i))
            enddo

c
c.... nonexistent element
c
         else
c
            call error ('filtprep  ', 'elem Cat', lelCat)
c
         endif
c
c.... end of generation
c
      enddo      

c     
c.... return
c     
      return
      end

