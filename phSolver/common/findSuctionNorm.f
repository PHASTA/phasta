       subroutine findSuctionNorm(x,iBC,ilwork,iper)

       use suctionDuct
       use pointer_data    
       include "common.h"
       include "mpif.h"
       integer iBC(nshg) 
       real*8 x(numnp,nsd)
       integer iper(nshg)
       integer ilwork(nlwork)
       integer, allocatable :: ienb(:)
       real*8 elnrm(nsd),wn(nsd),aedg(nsd),bedg(nsd) 
       real*8 mag
       integer isfID       

      allocate(wnorm(nshg,3))              
      wnorm=0

      do iblk=1,nelblb !loop over all the bounday element blocks at current process
        npro=lcblkb(1,iblk+1)-lcblkb(1,iblk)
        nenbl=lcblkb(6,iblk)
        nshl=lcblkb(9,iblk)
        allocate(ienb(nshl))
        do i=1,npro  ! loop over boundary elements
          isfID=miBCB(iblk)%p(i,2) ! miBCB(2) is the surfID
          if(isfID.ne.isetSuction_Duct) cycle
          ienb(1:nshl)=mienb(iblk)%p(i,1:nshl)
            do j=1,nenbl ! loop over boundary nodes
               nn=ienb(j) !nn is the global index of suction node
               if(j.ne.1.and.j.ne.nenbl)then
                 aedg(:)=x(ienb(j+1),:)-x(nn,:)
                 bedg(:)=x(ienb(j-1),:)-x(nn,:)
               elseif(j.eq.1)then
                 aedg(:)=x(ienb(j+1),:)-x(nn,:)
                 bedg(:)=x(ienb(nenbl),:)-x(nn,:)
               elseif(j.eq.nenbl)then
                 aedg(:)=x(ienb(1),:)-x(nn,:)
                 bedg(:)=x(ienb(j-1),:)-x(nn,:)   
               endif  
               elnrm(1)=aedg(2)*bedg(3)-aedg(3)*bedg(2)
               elnrm(2)=aedg(3)*bedg(1)-aedg(1)*bedg(3)
               elnrm(3)=aedg(1)*bedg(2)-aedg(2)*bedg(1)
               wnorm(nn,:)=wnorm(nn,:)+elnrm(:)
            enddo
        enddo
        deallocate(ienb)
      enddo
      
      if(numpe.gt.1) call commu(wnorm(:,:),ilwork,3,'in ') 
      call bc3per(iBC,wnorm(:,:),iper,ilwork,3)
      if(numpe.gt.1) call commu(wnorm(:,:),ilwork,3,'out')

      do nn=1,nshg ! normalize wnorm
        wn(:)=wnorm(nn,:)
        mag=(wn(1)**2+wn(2)**2+wn(3)**2)**0.5
        if(mag.ne.0) wnorm(nn,:)=wn(:)/mag             
      enddo      
 

c      call write_field(myrank,'a',
c     &      'wnrm',4,wnorm(:,:),'d',nshg,3,-2) 
cc 7 is the number of charaters of the name of array wnrm701
c        call write_field(myrank,'a',
c     &      'wnrm',4,wlnorm(:,:,2),'d',nshg,3,702)
c      stop
 
      return
      end

