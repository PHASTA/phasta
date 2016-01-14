      module wallData !was suctionDuct 
        logical wnormCalced
        !integer nsuctionface
        !integer, allocatable :: suctionf(:)

        !Threshold for determining whether a node is on a corner
        real*8 :: cornerCosThresh 
          
        !Wall normal array of size (nshg, 3) 
        real*8, allocatable :: wnorm(:,:)
      contains 
            
        subroutine findWallNorm(x,iBC,ilwork,iper)
   
!         use wallData
          use pointer_data    
          include "common.h"
          include "mpif.h"
   
          integer iBC(nshg) 
          real*8 x(numnp,nsd)
          integer iper(nshg)
          integer ilwork(nlwork)
          integer, allocatable :: ienb(:)
          real*8,  allocatable :: wNormTmp(:,:)
          real*8 eNorm(nsd),en(nsd),wn(nsd),aedg(nsd),bedg(nsd) 
          real*8 mag
!         integer isfID       
          
          real*8 wallBCout(nshg, 6)
          real*8 BC(nshg,ndofBC)
           
          cornerCosThresh = 0.75d00
          !The number 0.75 is somewhat arbitrary. The thinking is that this
          !code is also used for the normal on curved surfaces which will have
          !a non unit dot product from elements surrounding a wall vertex but
          !acos(0.75) is a pretty big angle. On coarse meshes this number may
          !have to be dialed down to get normal on curved walls correct. 
    
          allocate(wNorm(nshg,nsd))              
          wNorm = 0
    
          do iblk=1,nelblb !loop over all the bounday element blocks at current process
            npro=lcblkb(1,iblk+1)-lcblkb(1,iblk)
            nenbl=lcblkb(6,iblk)
            nshl=lcblkb(9,iblk)
            allocate(ienb(nshl))
            do i=1,npro  ! loop over boundary elements
!             isfID=miBCB(iblk)%p(i,2) ! miBCB(2) is the surfID
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
                eNorm(1)=aedg(2)*bedg(3)-aedg(3)*bedg(2)
                eNorm(2)=aedg(3)*bedg(1)-aedg(1)*bedg(3)
                eNorm(3)=aedg(1)*bedg(2)-aedg(2)*bedg(1)
                
                call addElemNormal(wNorm, eNorm, nn)
                 
              enddo  !end do j = 1, nenbl
            enddo  !end do i = 1, npro
            deallocate(ienb)
          enddo  !end do iblk = 1, nelblk
          
          !Now communicate between parts. Note: there still may be a
          !conflict if the boundary edge is on a corner, so a separate array
          !needs to be used to allow us to repeat the z-norm comparison.
          if(numpe.gt.1) then 
            allocate(wNormTmp(nshg,nsd))
            wNormTmp = wNorm
            call commu(wNormTmp(:,:), ilwork, 3, 'in ') 
            wNormTmp = wNormTmp - wNorm   !undo the addition in commu to just
                                          !get contributions from other processors
          
            do i = 1,nshg                 !loop over each node and add the off
              wn(:) = wNormTmp(i,:)       !processor contribution. If this is
              call addElemNormal(wNorm, wn, i) !an interior node, zero gets 
            enddo                              !added to zero. 
          endif
          
          !It is unclear how to treat periodic boundary conditions, i.e. 
          !should periodic walls have zero norm, a correct normal for each
          !side, or a correct normal for the master. With bc3per commented,
          !the wall normal should be calculated as though there is no 
          !periodic BC. 
!         call bc3per(iBC,wNorm(:,:),iper,ilwork,3)
    
          if(numpe.gt.1) call commu(wNorm(:,:),ilwork,3,'out')
          
          do nn = 1,nshg
            wn(:) = wNorm(nn,:)
            mag = (wn(1)**2+wn(2)**2+wn(3)**2)**0.5
           
            if(mag > 0) wNorm(nn,:) = wn(:)/mag
          enddo
    
          zwall=0.0762d00
          do nn=1,nshg ! normalize wNorm
            wn(:)=wNorm(nn,:)
            mag=(wn(1)**2+wn(2)**2+wn(3)**2)**0.5
            if(mag.ne.0) then
              if(((dabs(dabs(x(nn,3))-zwall)).lt.0.000000001d00) 
     &                  .and. (dabs(wn(3)/mag) .lt.0.99)) then 
! logic to trap wall normals in z above failed, possibley in commu? so
! report it and "fix it"
                write(*,*) 
     &             "Warning: still screwing up wall Norm, proc ", 
     &              myrank, ", node ", nn, 
     &              ", x = (", x(nn,1), x(nn,2), x(nn,3), ")", 
     &              ", n = (", wn(1), wn(2), wn(3), ")"
!                write(myrank+1000,1234)
!     &              nn,x(nn,1),x(nn,2),x(nn,3),wn(1),wn(2),wn(3)
                wn(1)=0
                wn(2)=0
                wn(3)=x(nn,3)
                mag=abs(wn(3))
              endif
              wNorm(nn,:)=wn(:)/mag             
            endif
          enddo      
    
1234      format(i5,6(2x,e12.5))
!!DEBUG
!!      dimension   u(nshg,n)
!      do j = 1,3
!        do i = 1,nshg
!         u(i,j)=9.876543e21
!          wallBCout(i, j  ) = BC(   i, j+2)
!          wallBCout(i, j+3) = wNorm(i, j  )
!        enddo
!      enddo
!      call write_restart(myrank,100002,nshg,6,wallBCout,wallBCout)
!!DEBUG
          
          wNormCalced = .true.
    
          return
        end subroutine
    
        subroutine addElemNormal(wNorm, eNorm, nn)
        !Adds contributions from the wall element normal en into the wall
        !normal array wNorm at index nn. 
            
          include "common.h"
    
          real*8 :: wNorm(nshg, nsd)
          real*8 :: eNorm(nsd), en(nsd), wn(nsd)
          real*8 :: mag, oDotn
          integer :: nn
            
          wn(:)=wNorm(nn,:)
          mag=(wn(1)**2+wn(2)**2+wn(3)**2)**0.5  
          
          if(mag < 1.0e-14) then  ! no direction yet so give it this
            wNorm(nn,:) = eNorm(:)
             
          else ! wNorm already has a direction 
            wn(:)=wn(:)/mag
            
            mag=(eNorm(1)**2+eNorm(2)**2+eNorm(3)**2)**0.5  
            en(:)=eNorm(:)/mag
            
            ! when we are in corners we need to establish a precidence since
            ! averaging will produce irregular results.  In this first pass,
            ! let's let the vector with the largest z-normal take precedence. 
            oDotn=en(1)*wn(1)+en(2)*wn(2)+en(3)*wn(3)
            if(abs(oDotn) < cornerCosThresh) then  ! this is a corner-Don't average
              if(abs(en(3)) > abs(wn(3))) then
                 wNorm(nn,:) = eNorm(:)  ! overwrite wNorm
              endif  ! since no else, wNorm unaffected by weaker z
                     ! e.g., for the else, eNorm is ignored
            else !  not a corner so average
              wNorm(nn,:)=wNorm(nn,:) + eNorm(:)
            endif  !end if(abs(oDotn) < cornerCosThresh)
          endif  !end mag < 1e-14
            
        end subroutine
    
      end module
                                                                        
