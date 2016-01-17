       subroutine slpwBC(nslpw,ihis,iBCg,BCg,BCtmpg,iBC,BC,wlnorm)
       include "common.h"
       include "mpif.h"
      
       integer nslpw,ihis,iBCg,iBC
       real*8 BCg(ndofBC)
       real*8 BCtmpg(ndof+7)
       real*8 BC(ndofBC)
       real*8 wlnorm(3,9)
       real*8 c4,c5,c6,c7
       real*8 c8,c9,c10,c11
       real*8 c12,c13,c14,c15
       real*8 ck1,ck2,ck3
       real*8 det,det3x3,tmp,u,v,w 
       integer nmax          
       real*8 wn1(3),wn2(3),wn3(3)
       
       icd      = ibits(iBCg,3,3)
       ixset    = ibits(iBCg,3,1)
       iyset    = ibits(iBCg,4,1)
       izset    = ibits(iBCg,5,1)
C
C...   release the previous setting for velocities
C
       iBC      = iBCg-(ixset*8+iyset*16+izset*32)

       
ccc case1: one slipwall and no user-specified velocity BC
       if(nslpw.eq.1.and.icd.eq.0)then 
        do islpw=1,9
           if(btest(ihis,islpw))exit
        enddo
        wn1(:)=wlnorm(:,islpw)
        ii=nmax(wn1(1),wn1(2),wn1(3))
        if(ii.eq.1)then  ! u has the largest constrain, u should be constrained
         iBC  = iBC+8       
         BC(3)= 0
         BC(4)= wn1(2)/wn1(1)
         BC(5)= wn1(3)/wn1(1) 
        elseif(ii.eq.2)then
         iBC=iBC+16       
         BC(3)=0
         BC(4)=wn1(1)/wn1(2)
         BC(5)=wn1(3)/wn1(2) 
        elseif(ii.eq.3)then
         iBC=iBC+32      
         BC(3)=0
         BC(4)=wn1(1)/wn1(3)
         BC(5)=wn1(2)/wn1(3) 
        endif

ccc case2: two slipwall and no user-specified velocity BC
       elseif(nslpw.eq.2.and.icd.eq.0)then 
        do islpw=1,9
           if(btest(ihis,islpw))exit
        enddo         
        do jslpw=islpw+1,9
           if(btest(ihis,jslpw))exit
        enddo      
        wn1(:)=wlnorm(:,islpw)         
        wn2(:)=wlnorm(:,jslpw) 
        c4=wn1(1)
        c5=wn1(2)
        c6=wn1(3)
        c7=0  
        c8=wn2(1)      
        c9=wn2(2) 
        c10=wn2(3)
        c11=0
C
C... ck is the cross product of wn1 and wn2
C
        ck1=c5*c10 - c9*c6
        ck2=c6*c8  - c4*c10
        ck3=c4*c9  - c8*c5

        kk=nmax(ck1,ck2,ck3)

        if(kk.eq.1)then  ! u has the largest freedom, v and w should be constrained
C
          iBC=iBC + 48
          det=c5*c10-c9*c6
          BC(3)=(c10*c7-c6*c11)/det
          BC(4)=(c10*c4-c6*c8)/det          
          BC(5)=(c11*c5-c9*c7)/det
          BC(6)=(c5*c8-c9*c4)/det          
        elseif(kk.eq.2)then
          iBC=iBC+40
          det=c4*c10-c8*c6
          BC(3)=(c10*c7-c6*c11)/det
          BC(4)=(c10*c5-c6*c9)/det          
          BC(5)=(c11*c4-c8*c7)/det
          BC(6)=(c4*c9-c5*c8)/det               
        elseif(kk.eq.3)then  
          iBC=iBC+24
          det=c4*c9-c8*c5
          BC(3)=(c9*c7-c5*c11)/det
          BC(4)=(c9*c6-c5*c10)/det          
          BC(5)=(c4*c11-c8*c7)/det
          BC(6)=(c4*c10-c8*c6)/det        
        endif

ccc case3: one slipwall and one user-specified velocity BC
       elseif(nslpw.eq.1.and.(icd.eq.1.or.icd.eq.2.or.icd.eq.4))then    
        do islpw=1,9
           if(btest(ihis,islpw))exit
        enddo                  
        wn1(:)=wlnorm(:,islpw)                
        c4=BCtmpg(4)
        c5=BCtmpg(5)
        c6=BCtmpg(6)
        c7=BCtmpg(7)
        c8=wn1(1)
        c9=wn1(2)
        c10=wn1(3)
        c11=0

        ck1=c5*c10 - c9*c6
        ck2=c6*c8  - c4*c10
        ck3=c4*c9  - c8*c5

        kk=nmax(ck1,ck2,ck3)
        if(kk.eq.1)then  ! u has the largest freedom
          iBC=iBC+48
          det=c5*c10-c9*c6
          BC(3)=(c10*c7-c6*c11)/det
          BC(4)=(c10*c4-c6*c8)/det          
          BC(5)=(c11*c5-c9*c7)/det
          BC(6)=(c5*c8-c9*c4)/det          
        elseif(kk.eq.2)then
          iBC=iBC+40
          det=c4*c10-c8*c6
          BC(3)=(c10*c7-c6*c11)/det
          BC(4)=(c10*c5-c6*c9)/det          
          BC(5)=(c11*c4-c8*c7)/det
          BC(6)=(c4*c9-c5*c8)/det               
        elseif(kk.eq.3)then  
          iBC=iBC+24
          det=c4*c9-c8*c5
          BC(3)=(c9*c7-c5*c11)/det
          BC(4)=(c9*c6-c5*c10)/det          
          BC(5)=(c4*c11-c8*c7)/det
          BC(6)=(c4*c10-c8*c6)/det        
        endif         

ccc case4: two slipwalls and one user-specified velocity BC
       elseif(nslpw.eq.2.and.(icd.eq.1.or.icd.eq.2.or.icd.eq.4))then    
        iBC=iBC+56
        c4=BCtmpg(4)
        c5=BCtmpg(5)
        c6=BCtmpg(6)
        c7=BCtmpg(7)                
        do islpw=1,9
           if(btest(ihis,islpw))exit
        enddo         
        do jslpw=islpw+1,9
           if(btest(ihis,jslpw))exit
        enddo      
        wn1(:)=wlnorm(:,islpw)         
        wn2(:)=wlnorm(:,jslpw)   
        c8=wn1(1)
        c9=wn1(2)
        c10=wn1(3)
        c11=0
        c12=wn2(1)
        c13=wn2(2)
        c14=wn2(3)
        c15=0
        tmp=det3x3(c4,c5,c6,c8,c9,c10,c12,c13,c14)
        BC(4)=det3x3(c4,c7,c6,c8,c11,c10,c12,c15,c14)/tmp
        BC(3)=det3x3(c7,c5,c6,c11,c9,c10,c15,c13,c14)/tmp
        BC(5)=det3x3(c4,c5,c7,c8,c9,c11,c12,c13,c15)/tmp

ccc case5: one slipwall and two user-specified velocity BC
       elseif(nslpw.eq.1.and.(icd.eq.3.or.icd.eq.5.or.icd.eq.6))then
        iBC=iBC+56
        c4=BCtmpg(4)
        c5=BCtmpg(5)
        c6=BCtmpg(6)
        c7=BCtmpg(7)        
        c8=BCtmpg(8)
        c9=BCtmpg(9)
        c10=BCtmpg(10)
        c11=BCtmpg(11)     
        do islpw=1,9
           if(btest(ihis,islpw))exit
        enddo
        wn1(:)=wlnorm(:,islpw)        
        c12=wn1(1)
        c13=wn1(2)
        c14=wn1(3)
        c15=0
        tmp=det3x3(c4,c5,c6,c8,c9,c10,c12,c13,c14)
        BC(4)=det3x3(c4,c7,c6,c8,c11,c10,c12,c15,c14)/tmp
        BC(3)=det3x3(c7,c5,c6,c11,c9,c10,c15,c13,c14)/tmp
        BC(5)=det3x3(c4,c5,c7,c8,c9,c11,c12,c13,c15)/tmp

ccc case6: two slipwalls and two user-specified velocity BC

ccc case7: one slipwall and three user-specified velocity BC
       elseif(nslpw.eq.1.and.icd.eq.7)then
        iBC=iBC+56
        u=BC(3)
        v=BC(4)
        w=BC(5)
        do islpw=1,9
           if(btest(ihis,islpw))exit
        enddo
        wn1(:)=wlnorm(:,islpw)
        ii=nmax(wn1(1),wn1(2),wn1(3))         
        if(ii.eq.1)then  ! u has the largest constrain       
         BC(3)=-wn1(2)/wn1(1)*v-wn1(3)/wn1(1)*w  ! u is determined by v and w
        elseif(ii.eq.2)then ! v has the largest constrain
         BC(4)=-wn1(1)/wn1(2)*u-wn1(3)/wn1(2)*w  ! v is determined by u and w
        elseif(ii.eq.3)then ! w has the largest constrain 
         BC(5)=-wn1(1)/wn1(3)*u-wn1(2)/wn1(3)*v  ! w is determined by u and v
        endif
  
ccc case8: two slipwall and three user-specified velocity BC
       elseif(nslpw.eq.2.and.icd.eq.7)then
        iBC=iBC+56
        u=BC(3)
        v=BC(4)
        w=BC(5)
        do islpw=1,9
           if(btest(ihis,islpw))exit
        enddo         
        do jslpw=islpw+1,9
           if(btest(ihis,jslpw))exit
        enddo      
        wn1(:)=wlnorm(:,islpw)         
        wn2(:)=wlnorm(:,jslpw) 
        c4=wn1(1)
        c5=wn1(2)
        c6=wn1(3)
        c7=0  
        c8=wn2(1)      
        c9=wn2(2) 
        c10=wn2(3)
        c11=0

        ck1=c5*c10 - c9*c6     ! ck is the cross product of wn1 and wn2
        ck2=c6*c8  - c4*c10
        ck3=c4*c9  - c8*c5

        kk=nmax(ck1,ck2,ck3)
        if(kk.eq.1)then  ! u has the largest freedom
          det=c5*c10-c9*c6
          BC(4)=(c10*c7-c6*c11)/det-(c10*c4-c6*c8)/det*u ! v is determined by u          
          BC(5)=(c11*c5-c9*c7)/det-(c5*c8-c9*c4)/det*u   ! w is determined by u          
        elseif(kk.eq.2)then ! v has the largest freedom
          det=c4*c10-c8*c6
          BC(3)=(c10*c7-c6*c11)/det-(c10*c5-c6*c9)/det*v ! u is determined by v          
          BC(5)=(c11*c4-c8*c7)/det-(c4*c9-c5*c8)/det*v   ! w is determined by v            
        elseif(kk.eq.3)then  ! w has the largest freedom
          det=c4*c9-c8*c5
          BC(3)=(c9*c7-c5*c11)/det-(c9*c6-c5*c10)/det*w  ! u is determined by w          
          BC(4)=(c4*c11-c8*c7)/det-(c4*c10-c8*c6)/det*w  ! v is determined by w       
        endif

ccc case9: more than three slipwalls is acutally non-slip wall, regardless of user-specified velocity BC 
       elseif(nslpw.ge.3)then 
          iBC=iBC+56
          BC(3:5)=0
ccc

       endif
        
       return
       end 

C-------------------------------------------------
C function nmax
C     returns the index at which the input number
C     has the largest absolute value. 
C     For example nmax( 0.1, -3.4, 1.2  ) = 2
C-------------------------------------------------

       integer function nmax(nx,ny,nz)
       real*8 nx,ny,nz
       real*8 lnx,lny,lnz
          lnx=abs(nx)
          lny=abs(ny)
          lnz=abs(nz)          
          if(lnx.gt.lny.and.lnx.gt.lnz)nmax=1
          if(lny.gt.lnx.and.lny.gt.lnz)nmax=2
          if(lnz.gt.lnx.and.lnz.gt.lny)nmax=3  
          return
       end

C-----------------------------------------------------------
C function det3x3 computes the determinant of a 3x3 matrix
C-----------------------------------------------------------
      
       real*8 function det3x3(a1,a2,a3,a4,a5,a6,a7,a8,a9)
       real*8 a1,a2,a3,a4,a5,a6,a7,a8,a9        
         det3x3=a1*a5*a9+a4*a8*a3+a7*a6*a2-a7*a5*a3-a8*a6*a1-a9*a4*a2
       return
       end 








