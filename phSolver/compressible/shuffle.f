        subroutine yshuffle(restmp, code)
c
c ... this subroutine is to shuffle the residual vector from the old format 
c ... of (P,u,v,w,P,T) to new format (u,v,w,T,P)
c
c input:  is the residual in old format
c output: is the residual in new format
c
	include "common.h"
c
        dimension restmp (nshg,nflow), tmp(nshg), tmp1(nshg,nflow)
        character*8 code
c
        if (code .eq. 'old2new ') then
        tmp(:)=restmp(:,1)   ! copying the res of continuity
        restmp(:, 1:3) = restmp(:, 2:4)
        restmp(:, 4)   = tmp(:)
        return
        endif
c
        if( code .eq. 'new2old ') then
           tmp1(:,:) = restmp(:,:)
           do i=1,nsd
              restmp(:,i+1) =  tmp1(:,i)
           enddo
           restmp(:,1) =  tmp1(:,4)
           return
        endif
        return
        end
c
c
      subroutine mshuffle(bdiagtmp)
c
c ... this subroutine is to shuffle the bdiag from new to old format
c
c input:  is the bdiag in old format
c output: is the bdiag in new format
c
      include "common.h"
c
      dimension bdiagtmp (nshg,nflow,nflow), tmp(nshg,nflow)
c
      tmp(:,:)  =  bdiagtmp(:,:,1)         ! reshuffling 1st column with 
      bdiagtmp(:,:,1:3) = bdiagtmp(:,:,2:4)    ! fourth column for all the nodes
      bdiagtmp(:,:,4) = tmp(:,:)
c
     
      tmp(:,:)  =  bdiagtmp(:,1,:)         ! reshuffling 1st row with
      bdiagtmp(:,1:3,:)= bdiagtmp(:,2:4,:)     ! fourth row for all the nodes
      bdiagtmp(:,4,:) = tmp (:,:)
c
      return
      end
