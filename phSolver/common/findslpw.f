       subroutine findslpw(x,ilwork,iper,iBC)
       use pointer_data    
       use slpw       
       include "common.h"
       include "mpif.h"
       real*8 x(numnp,nsd)
       integer iper(nshg)
       integer ilwork(nlwork)
       integer ifirstvisit(nshg)
       integer, allocatable :: ienb(:)
       real*8 elnrm(nsd),wn1(nsd),wn2(nsd),wn3(nsd),aedg(nsd),bedg(nsd) 
       real*8 mag
       integer iBC(nshg)
       integer, allocatable :: IDslpw(:)
       real*8 AA(nsd),BB(nsd),tmp
       integer iparall

       allocate(idxslpw(nshg))
       allocate(mpslpw(nshg,2))
       allocate(wlnorm(nshg,nsd,9))
             
              
      nslpwnd=0 
      idxslpw=0
      mpslpw=0
      wlnorm=0
      ifirstvisit(:)=1
      do iblk=1,nelblb 
      npro=lcblkb(1,iblk+1)-lcblkb(1,iblk)
      nenbl=lcblkb(6,iblk)
      nshl=lcblkb(9,iblk)
      allocate(ienb(nshl))
      do i=1,npro 
      isfID=miBCB(iblk)%p(i,2)
      if(isfID.gt.709.or.isfID.lt.701)cycle
      islpw=mod(isfID,700)
      ienb(1:nshl)=mienb(iblk)%p(i,1:nshl)
      do j=1,nenbl
         nn=ienb(j)
         if(ifirstvisit(nn).eq.1)then
            ifirstvisit(nn)=0
            nslpwnd=nslpwnd+1
            idxslpw(nslpwnd)=nn
         endif
         if(.not.btest(mpslpw(nn,2),islpw))then
            mpslpw(nn,1)=mpslpw(nn,1)+1   
            mpslpw(nn,2)=mpslpw(nn,2)+2**islpw
         endif
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
         wlnorm(nn,:,islpw)=wlnorm(nn,:,islpw)+elnrm(:)
      enddo
      enddo
      deallocate(ienb)
      enddo
      
       do k=1,9
         if(numpe.gt.1) call commu(wlnorm(:,:,k),ilwork,3,'in ') 
         call bc3per(iBC,wlnorm(:,:,k),iper,ilwork,3)
         if(numpe.gt.1) call commu(wlnorm(:,:,k),ilwork,3,'out')
       enddo
cccccccccccccccccccccccccccccccccccccccccccccccccc   
      
     
      do id=1,nslpwnd
       nn=idxslpw(id)
       nslpw=mpslpw(nn,1)
       ihis=mpslpw(nn,2)
       nslpwg=nslpw
 
       allocate(IDslpw(nslpw))
       jslpw=0 
       do ni=1,nslpw
          do islpw=1+jslpw,9
             if(btest(ihis,islpw))exit 
          enddo
          IDslpw(ni)=islpw
          jslpw=islpw 
       enddo  

       ired=0   
       do i=1,nslpwg
          if(btest(ired,i))cycle 
       do j=i+1,nslpwg
          AA(:)=wlnorm(nn,:,IDslpw(i))          
          BB(:)=wlnorm(nn,:,IDslpw(j))
          if(iparall(AA,BB).eq.1)then
             nslpw=nslpw-1
             ihis=ihis-2**IDslpw(j)
             wlnorm(nn,:,IDslpw(i))=AA(:)+BB(:) 
             ired=ired+2**j  
          endif
       enddo
       enddo            

       do k=1,9
          wn1(:)=wlnorm(nn,:,k)
          mag=(wn1(1)**2+wn1(2)**2+wn1(3)**2)**0.5
          if(mag.ne.0)wlnorm(nn,:,k)=wn1(:)/mag             
       enddo
            
      mpslpw(nn,1)=nslpw
      mpslpw(nn,2)=ihis
 
      deallocate(IDslpw) 
      enddo
    

c        call write_field(myrank,'a'//char(0),
c     &      'wnrm'//char(0),4,wlnorm(:,:,1),'d'//char(0),nshg,3,701) 
cc 7 is the number of charaters of the name of array wnrm701
c        call write_field(myrank,'a'//char(0),
c     &      'wnrm'//char(0),4,wlnorm(:,:,2),'d'//char(0),nshg,3,702)
c      stop
 
      return
      end

cccccccccccccccccccccc

      integer function iparall(A,B)
      real*8 A(3),B(3)
      real*8 angtol,tmp
 
      tmp=(A(1)**2+A(2)**2+A(3)**2)**0.5
      A(:)=A(:)/tmp
      tmp=(B(1)**2+B(2)**2+B(3)**2)**0.5
      B(:)=B(:)/tmp
      pi=3.1415926535
      angtol=10
      angtol=sin(angtol*(pi/180))
      C1=A(2)*B(3)-A(3)*B(2)
      C2=A(3)*B(1)-A(1)*B(3)   
      C3=A(1)*B(2)-A(2)*B(1)
      Cmag=(C1**2+C2**2+C3**2)**0.5
      dotp=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      if(dotp.gt.0.and.Cmag.lt.angtol)then
         iparall=1
      else
         iparall=0
      endif

      return
      end 
