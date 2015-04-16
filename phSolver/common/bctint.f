c-----------------------------------------------------------------------
c
c  This module conveys temporal BC data.  Below functions read in the data
c  and interpolate it to the current time level. 
c
c-----------------------------------------------------------------------
      module specialBC
      use pointer_data

      real*8, allocatable ::  BCt(:,:), acs(:,:), spamp(:)
      real*8, allocatable ::  ytarget(:,:)
      real*8, allocatable ::  BCperiod(:)
      integer, allocatable :: nBCt(:), numBCt(:)
     
      type(r2d),dimension(:),allocatable :: BCtptr 
        ! BCtptr%p(:,:) is to replace BCt(ic,:,:), in which 
        ! the second index is dynamically decided in
        ! the setup stage.

      integer ntv,nptsmax
c$$$      integer itvn
      end module

c-----------------------------------------------------------------------
c
c  This module conveys flow rate history for the different impedance outlets
c  over one period. Below functions read in the data and store it for the
c  current time level. 
c
c-----------------------------------------------------------------------
      module convolImpFlow

      real*8, allocatable ::  QHistImp(:,:), ValueImpt(:,:,:)
      real*8, allocatable ::  ValueListImp(:,:), ConvCoef(:,:)
      real*8, allocatable ::  ImpConvCoef(:,:), poldImp(:)
      integer ntimeptpT,numDataImp
      integer, allocatable :: nImpt(:), numImpt(:)
      integer nptsImpmax
      real*8, allocatable ::  QHistTry(:,:), QHistTryF(:,:) !filter
      integer cutfreq !filter
      end module
c-----------------------------------------------------------------------
c
c  This module conveys the parameters for the different RCR outlets.
c  Below functions read in the inputs (proximal resistance, capacitance, 
c  distal resistance and distal pressure) and store it for the
c  current time level. 
c
c-----------------------------------------------------------------------
      module convolRCRFlow

      real*8, allocatable ::  ValueListRCR(:,:), ValuePdist(:,:,:) !inputs
      real*8, allocatable ::  QHistRCR(:,:), HopRCR(:) !calc
      real*8, allocatable ::  RCRConvCoef(:,:), poldRCR(:) !calc
      real*8, allocatable ::  dtRCR(:) !scaled timestep: deltat/RdC
      real*8, allocatable ::  RCRic(:) !(P(0)-RQ(0)-Pd(0))
      integer nptsRCRmax,numDataRCR !to read inputs
      integer, allocatable :: numRCRt(:) !to read inputs
      end module

c-----------------------------------------------------------------------
c
c     Initialize:
c
c-----------------------------------------------------------------------
      subroutine initSponge( y,x)
      
      use     specialBC
      include "common.h"
      
      real*8   y(nshg,nflow), x(numnp,3)
      allocate (ytarget(nshg,nflow))  
      
      if(matflg(5,1).eq.5) then
         write(*,*) 'calculating IC sponge'
         ytarget = y
      else
         write(*,*) 'calculating Analytic sponge'

c
c OLD style sponge pushed onto target.  You need to be sure that your
c solver.inp entries for start and stop of sponge match as well as the
c growth rates
c
      vcl=datmat(1,5,1)         ! velocity on centerline
      rslc=datmat(2,5,1)        ! shear layer center radius
      bfz=datmat(3,5,1)
      we=3.0*29./682.
      rsteep=3.0
      zstart=30.0
      radst=10.0
      radsts=radst*radst
      do id=1,numnp
         radsqr=x(id,2)**2+x(id,1)**2
c         if((x(id,3).gt. zstart) .or. (radsqr.gt.radsts))  then
            rad=sqrt(radsqr)
            radc=max(rad,radst)
            zval=max(x(id,3),zstart)
            utarget=(tanh(rsteep*(rslc-rad))+one)/two*
     &                    (vcl-we) + we
            Ttarget  = press/(ro*Rgas)
            ptarget= press
            ytarget(id,1) = zero
            ytarget(id,2) = zero
            ytarget(id,3) = utarget
            ytarget(id,4) = ptarget
            ytarget(id,5) = Ttarget            
c         endif
      enddo
      endif
      return
      end


c-----------------------------------------------------------------------
c
c     Initialize:time varying boundary condition
c
c-----------------------------------------------------------------------
      subroutine initBCt( x, iBC, BC )
      
      use     specialBC
      include "common.h"
 
      real*8   x(numnp,nsd), BC(nshg,ndofBC), rj1,rj2,rj3,rj4,distd,epsd
      integer  iBC(nshg)
      character*80 card
      real*8 distds
      real*8 dd

      integer irstart
      real(kind=8),allocatable,dimension(:) ::  bcttest

      real*8 t0,tlen,t1,tstart,tend
      integer i0,ilen,i1,nper,istart,iend
c
c  This one should be used for boundary layer meshes where bct.dat must
c  be given to greater precision than is currently being generated.
c
c      epsd=1.0d-12    ! this is distance SQUARED to save square root

      epsd=1.0d-8              ! this is distance SQUARED to save square root

      ic=0                      !count the number on this processor

! ************ Chun Sun: limit the memory use if required time steps
! ************ is smaller than total bct.dat length.

      ! readin the start timestep
      
      open(unit=72,file='numstart.dat',status='old')
      read(72,*) irstart
      close(72)

      ! use irstart-1 and nstep+1 to avoid machine tolerance issues
      t0 = max(zero,(irstart-1)*Delt(1))
      tlen = (nstep(1)+1)*Delt(1) 
      ! Assumption: constant time step in one run. If you use variable
      ! time step in one run, this should be modified.
      t1 = t0 + tlen

      if(myrank.eq.master)
     &   write(*,*) 'necessary bct timing: from ',t0,' to ',t1

! **************************************************************      
      
    
      if(any(ibits(iBC,3,3).eq.7)) then
         if(myrank.eq.master) write(*,*) 'opening bct.dat'
c         open(unit=567, file='bct.dat',status='old')
         open(unit=567, file='bct.dat',ACTION='READ',STATUS='old')
         read (567,'(a80)') card
           read (card,*) ntv, nptsmax
c        read(567,*) ntv,nptsmax
         allocate (nBCt(numnp))  
         allocate (numBCt(ntv))  
         allocate (BCt(nptsmax,4))
         allocate (BCperiod(ntv))
         allocate (BCtptr(ntv))  ! dynamic bct allocation
         do k=1,ntv
            read(567,*) x1,x2,x3,ntpts
c
c Find the point on the boundary (if it is on this processor)
c that matches this point
c
            do i=1,numnp
               if(ibits(ibc(i),3,3) .eq.7) then
                  dd= distds(x1,x2,x3,x(i,1),x(i,2),x(i,3))
                  if(dd.lt.epsd) then
                     ic=ic+1
                     nBCt(ic)=i ! the pointer to this point
!                     numBCt(ic)=ntpts ! the number of time series
                     do j=1,ntpts
c                        read(567,*) BCt(ic,j,4),(BCt(ic,j,n),n=1,3)
                        read(567,*) (BCt(j,n),n=1,4)
                     enddo
              ! at this point we have all bc data of spacial point
              ! ic/i in all time domain. now we figure out which subset
              ! is necessary.
                     if (tlen.ge.BCt(ntpts,4)) then
                         ! whole run is larger than whole period
                         ! should take all data
                         ilen = ntpts
                         allocate(BCtptr(ic)%p(ilen,4))
                         BCtptr(ic)%p = BCt(1:ilen,:)
                     else
                         nper = t0/BCt(ntpts,4)
                         tstart = t0-nper*BCt(ntpts,4)
                         nper = t1/BCt(ntpts,4)
                         tend = t1-nper*BCt(ntpts,4)
                         istart = iBfind(BCt(:,4),ntpts,tstart)
                         iend = iBfind(BCt(:,4),ntpts,tend)
                         if (istart>1) istart=istart-1
                         if (iend<ntpts) iend=iend+1
                         !write(*,*)'bcstart:',BCt(istart,4),tstart
                         !write(*,*)'bcend:',BCt(iend,4),tend
                         if (tstart.lt.tend) then ! does not loop
                             ilen = iend-istart+1
                             allocate(BCtptr(ic)%p(ilen,4))
                             BCtptr(ic)%p = BCt(istart:iend,:)
                         else ! loop within two consequetive time period
                             i0 = (ntpts-istart+1) ! first segment
                             ilen = iend + i0
                             allocate(BCtptr(ic)%p(ilen,4))
                      BCtptr(ic)%p(1:i0,:) = BCt(istart:ntpts,:)
                      BCtptr(ic)%p(i0+1:ilen,:) = BCt(1:iend,:)
                      BCtptr(ic)%p(i0+1:ilen,4) =
     &                 BCtptr(ic)%p(i0+1:ilen,4) + BCt(ntpts,4) 
                         endif
                     endif
                     numBCt(ic) = ilen
              BCtptr(ic)%p(:,4) = BCtptr(ic)%p(:,4)*bcttimescale
              BCperiod(ic) = BCt(ntpts,4)*bcttimescale
                     exit
                  endif
               endif
            enddo
            if(i.eq.numnp+1) then
c
c  if we get here the point was not found.  It must be on another
c  processor so we read past this record and move on
c
               do j=1,ntpts
                  read(567,*) rj1,rj2,rj3,rj4
               enddo
            endif
         enddo                  ! end of the loop over ntv
         !BCt(:,:,4)=BCt(:,:,4)*bcttimescale
         deallocate(BCt)
      endif                     ! any 3 component nodes
      itvn=ic
      close(567)
      if (ic.gt.0) then
         write(*,*)'myrank=',myrank,' and I found ',ic,' nodes.'
c      else
c         deallocate(nBCt)
c         deallocate(numBCt)
c         deallocate(BCt)
      endif

      return
      end

!***************************************************************
!      A Binary search return an index of a real array.
!      This array should be an ascending sorted array.
!***************************************************************      
      function iBfind(bArray,bLen,bVal)
      integer,intent(in) :: bLen
      real(kind=8),intent(in),dimension(bLen) :: bArray
      real(kind=8),intent(in) :: bVal
      integer mlen,indx,bstart,bend
      bstart = 1
      bend = bLen
      indx = (bLen+1)/2
      do while ((bstart+1).lt.bend)
         if (bVal.gt.bArray(indx)) then 
             bstart = indx
         else
             bend = indx
         endif
         indx = (bstart+bend)/2
      enddo
      iBfind = indx
      return
      end function


      subroutine BCint(timel,shp,shgl,shpb,shglb,x,BC,iBC)

      use     specialBC ! brings in itvn,nbct, bct, numbct, nptsmax

      include "common.h"

      real*8   BC(nshg,ndofBC), timel,t
      real*8   x(numnp,nsd),   
     &         shp(MAXTOP,maxsh,MAXQPT),
     &         shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &         shpb(MAXTOP,maxsh,MAXQPT),
     &         shglb(MAXTOP,nsd,maxsh,MAXQPT)

      integer  iBC(numnp),nlast,i,j,nper 

      do i =1,itvn ! itvn is the number of varying nodes on this proc 

         nlast=numBCt(i)     ! number of time series to interpolate from
         nper=timel/BCperiod(i) 
        ! number of periods completed to shift off


         t=timel-nper*BCperiod(i)  ! now time in periodic domain

         do j=2,nlast   !loop to find the interval that we are in

            if(BCtptr(i)%p(j,4).gt.t) then  
            ! this is upper bound, j-1 is lower

      wr=(t-BCtptr(i)%p(j-1,4))/(BCtptr(i)%p(j,4)-BCtptr(i)%p(j-1,4))
               BC(nbct(i),3:5)= BCtptr(i)%p(j-1,1:3)*(one-wr) 
     &                        + BCtptr(i)%p(j,1:3)*wr
               exit

            endif
         enddo
      enddo
      return
      end

      function distds(x1,y1,z1,x2,y2,z2)
      real*8 distds 
      real*8 x1,y1,z1,x2,y2,z2,x,y,z
      x=x1-x2
      y=y1-y2
      z=z1-z2
      distds=x*x+y*y+z*z
      return
      end
c-----------------------------------------------------------------------
c   initialize the impedance boundary condition:
c   read the data in initImpt
c   interpolate the data to match the process time step in Impint
c-----------------------------------------------------------------------
      subroutine initImpt()
      
      use convolImpFlow
      include "common.h"

      open(unit=817, file='impt.dat',status='old')
         read (817,*) nptsImpmax
         allocate (numImpt(numImpSrfs))  
         allocate (ValueImpt(nptsImpmax,2,numImpSrfs))
         ValueImpt=0
         do k=1,numImpSrfs
            read (817,*) numDataImp
            numImpt(k) = numDataImp
            do j=1,numDataImp
               read(817,*) (ValueImpt(j,n,k),n=1,2) ! n=1 time, 2 value
            enddo
         enddo
      close(817)
      
      allocate (ValueListImp(ntimeptpT+1,numImpSrfs))
      ValueListImp = zero 
      ValueListImp(ntimeptpT+1,:) = ValueImpt(1,2,:) !Z(time=0), last entry
      ValueListImp(1,:) = ValueImpt(1,2,:) !Z(time=0)=Z(time=T)
      return
      end
      
      
      
      subroutine Impint(ctime,jstep)
      
      use convolImpFlow
      include "common.h"
      
      real*8 ctime, ptime
      integer nlast, nper, k, j , jstep
      
         
      do k =1,numImpSrfs
         nlast=numImpt(k)     ! number of time series to interpolate from
         nper=ctime/ValueImpt(nlast,1,k)!number of periods completed to shift off
         ptime = ctime-nper*ValueImpt(nlast,1,k)  ! now time in periodic domain
            
         do j=2,nlast   !loop to find the interval that we are in

            if(ValueImpt(j,1,k).gt.ptime) then  ! this is upper bound, j-1 is lower
               wr=(ptime-ValueImpt(j-1,1,k))
     &             / ( ValueImpt(j,1,k)-ValueImpt(j-1,1,k) )
               ValueListImp(jstep,k)= ValueImpt(j-1,2,k)*(one-wr) 
     &                        + ValueImpt(j,2,k)*wr
               exit
            endif

         enddo
      enddo
      return
      end

c-----------------------------------------------------------------------------
c     time filter for a periodic function (sin cardinal + window function)     
c     is used for the impedance and the flow rate history
c-----------------------------------------------------------------------------
      subroutine Filter(Filtered,DataHist,nptf,timestep,cutfreq)
      
      include "common.h"

      integer nptf, cutfreq, j, k, m, s, Filtime(nptf)
      real*8  DataHist(nptf,numImpSrfs), Window(nptf)
      real*8  Sinc(nptf), FilterSW(nptf), Filtered(nptf,numImpSrfs)
      real*8  windK, timestep

      windK = cutfreq*2 + 1
      do j=1,nptf
         Filtime(j) = j-1
         Window(j) = 0.42+0.5*cos(2*pi*Filtime(j)/nptf)
     &              +0.08*cos(4*pi*Filtime(j)/nptf)
         Sinc(j) = sin(pi*Filtime(j)*windK/nptf)/sin(pi*Filtime(j)/nptf) 
      enddo          
      Sinc(1) = windK
            
      do j=1,nptf     
         FilterSW(j) = Window(nptf+1-j)*Sinc(nptf+1-j) !filter for convolution
      enddo     
      
      Filtered = zero
      do m=1,nptf
         do j=1,nptf
            s=modulo(m-nptf+j,nptf)
            if(s.eq.zero) then
               s=nptf
            endif
            Filtered(m,:) = Filtered(m,:)
     &              +FilterSW(j)*DataHist(s,:)/nptf !filter convolution
         enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
c   initialize the RCR boundary condition:
c   read the data in initRCRt
c   interpolate the data to match the process time step in RCRint
c-----------------------------------------------------------------------
      subroutine initRCRt()
      
      use convolRCRFlow
      include "common.h"

      open(unit=818, file='rcrt.dat',status='old')
         read (818,*) nptsRCRmax
         allocate (numRCRt(numRCRSrfs))  
         allocate (ValuePdist(nptsRCRmax,2,numRCRSrfs))
         allocate (ValueListRCR(3,numRCRSrfs))
         ValuePdist=0
         ValueListRCR=0
         do k=1,numRCRSrfs
            read (818,*) numDataRCR
            numRCRt(k) = numDataRCR
            do j=1,3
               read(818,*) ValueListRCR(j,k) ! reads Rp,C,Rd
            enddo
            do j=1,numDataRCR
               read(818,*) (ValuePdist(j,n,k),n=1,2) ! n=1 time, 2 value
            enddo
         enddo
      close(818)

      allocate (dtRCR(numRCRSrfs))

      return
      end
           
      
      subroutine RCRint(ctime,Pdist)
      
      use convolRCRFlow ! brings numRCRSrfs, ValuePdist
      include "common.h"
      
      real*8  ctime, ptime
      integer nlast, nper, k, j
      real*8  Pdist(0:MAXSURF)      
         
      do k =1,numRCRSrfs
         nlast=numRCRt(k)     ! number of time series to interpolate from
         nper=ctime/ValuePdist(nlast,1,k)!number of periods completed to shift off
         ptime = ctime-nper*ValuePdist(nlast,1,k)  ! now time in periodic domain
            
         do j=2,nlast   !loop to find the interval that we are in

            if(ValuePdist(j,1,k).gt.ptime) then  ! this is upper bound, j-1 is lower
               wr=(ptime-ValuePdist(j-1,1,k))
     &             / ( ValuePdist(j,1,k)-ValuePdist(j-1,1,k) )
               Pdist(k)= ValuePdist(j-1,2,k)*(one-wr) 
     &                        + ValuePdist(j,2,k)*wr
               exit
            endif

         enddo
      enddo
      return
      end

c----------------------------------------------------------------------- 
c     returns in pold the history dependent part of the pressure in the
c     impedance/flow rate convolution for the impedance and RCR BC      
c-----------------------------------------------------------------------      
      subroutine pHist(pressHist,QHist,betas,nTimePoint,nSrfs)

      include "common.h"
      
      integer  nTimePoint,nSrfs
      real*8   pressHist(0:MAXSURF)
      real*8   QHist(nTimePoint+1,nSrfs),betas(nTimePoint+2,nSrfs)
      !don't need here betas(ntimePoint+2)
      !but pb of array passing if cut at nTimePoint+1
      pressHist=zero
      do k=1,nSrfs
        do j=1,nTimePoint+1
            pressHist(k) = pressHist(k) + QHist(j,k)*betas(j,k)
        enddo
      enddo
      return
      end
