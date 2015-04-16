c-----------------------------------------------------------------------
c
c     module for time averaged statistics (conservative projection).
c
c-----------------------------------------------------------------------
      module stats
      
      integer nResDims, nSolDims, nLhsDims, nTimeStep, stsResFlg
      integer stsCompFreq, stsWriteFreq, stsResetFreq, step1,
     &        stsType
      
      real*8, allocatable :: stsVec(:,:)
      
      real*8, allocatable :: stsReg(:)
      real*8, allocatable :: stsMInv(:,:)
      real*8, allocatable :: stsB(:,:)
      real*8, allocatable :: stsDInv(:,:)
      real*8, allocatable :: stsCInv(:,:)
      
      real*8, allocatable :: stsPres(:), stsPresSqr(:), stsVel(:,:),
     &                       stsVelSqr(:,:), stsVelReg(:,:),
     &                       stsStress(:,:)

      end module
      
c-----------------------------------------------------------------------
c     create the new statistics arrays
c-----------------------------------------------------------------------
      subroutine initStats(x,   iBC,    iper,   ilwork)
      
      use stats
      include "common.h"

      real*8  x(numnp,3)
      integer ilwork(nlwork), iper(nshg), iBC(nshg)
      
      if (ipord .eq. 1) then
         stsType      = 1
      else
         stsType      = 2
      endif
      
      stsWriteFreq = 200
      
      nResDims = 11
      nSolDims = 10
      nLhsDims = 19
      
      allocate ( stsVec(nshg,nResDims) )

      if (stsType .eq. 1) then
         allocate ( stsReg(nshg)          )
         allocate ( stsMInv(nshg,6)       )
         allocate ( stsB(nshg,3)          )
         allocate ( stsDInv(nshg,3)       )
         allocate ( stsCInv(nshg,6)       )
      endif

      allocate ( stsPres(nshg)         )
      allocate ( stsPresSqr(nshg)      )
      allocate ( stsVel(nshg,3)        )
      allocate ( stsVelSqr(nshg,6)     )
      allocate ( stsVelReg(nshg,3)     )
      allocate ( stsStress(nshg,6)     )
      
      stsPres    = 0.0
      stsPresSqr = 0.0
      stsVel     = 0.0
      stsVelSqr  = 0.0
      stsVelReg  = 0.0
      stsStress  = 0.0

      step1      = lstep+1
      nTimeStep  = 0
      stsResFlg  = 0

      if (stsType .eq. 1) then
         call elmStatsLhs( x,   iBC,   iper,   ilwork )
         call stsInitLhs(  nshg )
      endif

      return
      end
      
c-----------------------------------------------------------------------
c     compute the Lhs matrices needed for the conservative projection
c     of the statistics
c-----------------------------------------------------------------------
      subroutine stsInitLhs(nshg)

      use     stats
      integer nshg

      real*8  res(nResDims), reg, det, r0, r1, r2, r3, r4, r5,
     &        d0, d1, d2, c0, c1, c2, c3, c4, c5
      integer i
      
c
c.... build the regularization
c      
      do i = 1, nshg
         res = stsVec(i,:)
         reg = res(1) * res(2) * res(3)
         det = res(1) * (res(2) * res(3) - res(5) * res(5))
     &       + res(4) * (res(5) * res(6) - res(4) * res(3))
     &       + res(6) * (res(4) * res(5) - res(2) * res(6)) 
	
         if ( det .gt. 1.d-10*reg .and. reg .ne. 0 ) then
	    stsReg(i) = 0 
         else 
	    stsReg(i) = (res(1) + res(2) + res(3)) / 1000. 
         endif
      enddo
      
c
c.... form M and factorize
c
      do i = 1, nshg
         res   = stsVec(i,:)
         reg   = stsReg(i)
         r0    = res(1) + reg
         r1    = res(2) + reg         
         r2    = res(3) + reg         
         r3    = res(4)
         r4    = res(5)
         r5    = res(6)
         
         det   = r0 * (r1 * r2 - r4 * r4)
     &         + r3 * (r4 * r5 - r3 * r2)
     &         + r5 * (r3 * r4 - r1 * r5)
         det   = 1.0/det
         
         stsMInv(i,1) = det * (r1 * r2 - r4 * r4)
         stsMInv(i,2) = det * (r0 * r2 - r5 * r5)
         stsMInv(i,3) = det * (r0 * r1 - r3 * r3)
         stsMInv(i,4) = det * (r4 * r5 - r2 * r3)
         stsMInv(i,5) = det * (r3 * r5 - r0 * r4)
         stsMInv(i,6) = det * (r3 * r4 - r1 * r5)
      enddo

c
c.... form B, DInv and CInv      
c
      do i = 1, nshg
	res          = stsVec(i,:)
	reg          = stsReg(i) 
	r0           = res(1) 
	r1           = res(2) 
	r2           = res(3) 
	r3           = res(4) 
	r4           = res(5) 
	r5           = res(6) 
	d0           = 1. / ( reg + r0 ) 
	d1           = 1. / ( reg + r1 ) 
	d2           = 1. / ( reg + r2 ) 
	stsDInv(i,1) = d0 
	stsDInv(i,2) = d1 
	stsDInv(i,3) = d2 
	stsB(i,1)    = r3 
	stsB(i,2)    = r4 
	stsB(i,3)    = r5 
	c0           = r0 + r1 - r3 * r3 * (d0 + d1) + reg 
	c1           = r1 + r2 - r4 * r4 * (d1 + d2) + reg 
	c2           = r2 + r0 - r5 * r5 * (d2 + d0) + reg 
	c3           = r5      - r3 * r4 * d1 
	c4           = r3      - r4 * r5 * d2 
	c5           = r4      - r5 * r3 * d0 
	det          = c0 * (c1 * c2 - c4 * c4)
     &               + c3 * (c4 * c5 - c3 * c2)
     &               + c5 * (c3 * c4 - c1 * c5) 
	det          = 1. / det 
	stsCInv(i,1) = det * (c1 * c2 - c4 * c4) 
	stsCInv(i,2) = det * (c0 * c2 - c5 * c5) 
	stsCInv(i,3) = det * (c0 * c1 - c3 * c3) 
	stsCInv(i,4) = det * (c4 * c5 - c2 * c3) 
	stsCInv(i,5) = det * (c3 * c5 - c0 * c4) 
	stsCInv(i,6) = det * (c3 * c4 - c1 * c5) 
      enddo
      
      return
      end
      
         
c-----------------------------------------------------------------------
c     collect the desired statistics 
c-----------------------------------------------------------------------
      subroutine stsGetStats( y,      yold,   ac,     acold, 
     &                        u,      uold,   x,
     &                        shp,    shgl,   shpb,   shglb,
     &                        iBC,    BC,     iper,   ilwork,
     &                        rowp,   colm,   lhsK,   lhsP )
      
      use     stats

      include "common.h"
      
      real*8  y(nshg,ndof),             yold(nshg,ndof),
     &        ac(nshg,ndof),            acold(nshg,ndof),
     &        u(nshg,nsd),              uold(nshg,nsd),
     &        shp(MAXTOP,maxsh,MAXQPT),  shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &        shpb(MAXTOP,maxsh,MAXQPT),
     &        shglb(MAXTOP,nsd,maxsh,MAXQPT),
     &        BC(nshg,ndofBC),          lhsK(9,nnz_tot),
     &        lhsP(4,nnz_tot),         x(numnp,nsd)

      integer iBC(nshg),                iper(nshg),
     &        ilwork(nlwork),           rowp(nshg*nnz),
     &        colm(nshg+1)
      
      
      real*8  yAlpha(nshg,ndof),      acAlpha(nshg,ndof),
     &        uAlpha(nshg,nsd),
     &        res(nResDims),          MInv(6),
     &        DInv(3),                B(3),
     &        CInv(6)
      
      real*8 u1, u2, u3, r0, r1, r2, r3, r4, r5, t3, t4, t5

      integer i
      
      nTimeStep = nTimeStep + 1
c
c.... compute solution at n+alpha
c
      call itrYAlpha( uold,     yold,     acold,
     &                u,        y,        ac,  
     &                uAlpha,   yAlpha,   acAlpha)
      
c
c.... assemble the residual
c
      if (stsType .eq. 1) then
         call elmStatsRes( yAlpha,   acAlpha,     x,       shp,   shgl, 
     &                     shpb,     shglb,       iBC,     BC, 
     &                     iper,     ilwork,      rowp,    colm,
     &                     lhsK,     lhsP  )

c
c.... compute the statistics
c
         do i = 1, nshg
            res   = stsVec(i,:)
            reg   = stsReg(i)
      
            MInv  = stsMInv(i,:)
            DInv  = stsDInv(i,:)
            B     = stsB(i,:)
            CInv  = stsCInv(i,:)
            
            u1    = yAlpha(i,1)
            u2    = yAlpha(i,2)
            u3    = yAlpha(i,3)
            
            stsPres(i)    = stsPres(i)    +  y(i,4) 
            stsPresSqr(i) = stsPresSqr(i) +  y(i,4)*y(i,4)  

            r0            = res(1) + reg * u1 
            r1            = res(2) + reg * u2 
            r2            = res(3) + reg * u3 
         
            stsVel(i,1)   = stsVel(i,1) 
     &                    + MInv(1) * r0 + MInv(4) * r1 + MInv(6) * r2 
            stsVel(i,2)   = stsVel(i,2)
     &                    + MInv(4) * r0 + MInv(2) * r1 + MInv(5) * r2 
            stsVel(i,3)   = stsVel(i,3)
     &                    + MInv(6) * r0 + MInv(5) * r1 + MInv(3) * r2 
            
            stsVelReg(i,1) = stsVelReg(i,1) + u1 
            stsVelReg(i,2) = stsVelReg(i,2) + u2 
            stsVelReg(i,3) = stsVelReg(i,3) + u3 
            
            r0	        = res(1) * u1               + reg * u1 * u1 
            r1		= res(2) * u2               + reg * u2 * u2 
            r2		= res(3) * u3               + reg * u3 * u3 
            r3	        = res(1) * u2 + res(2) * u1 + reg * u2 * u1 
            r4		= res(2) * u3 + res(3) * u2 + reg * u3 * u2 
            r5		= res(3) * u1 + res(1) * u3 + reg * u1 * u3 
            r0          = DInv(1) * r0 
            r1          = DInv(2) * r1 
            r2          = DInv(3) * r2 
            r3          = r3 - B(1) * (r0 + r1) 
            r4          = r4 - B(2) * (r1 + r2) 
            r5          = r5 - B(3) * (r2 + r0) 
            t3          = CInv(1) * r3 + CInv(4) * r4 + CInv(6) * r5 
            t4          = CInv(4) * r3 + CInv(2) * r4 + CInv(5) * r5 
            t5          = CInv(6) * r3 + CInv(5) * r4 + CInv(3) * r5 
            
            stsVelSqr(i,1) = stsVelSqr(i,1)  
     &                  + r0 - DInv(1) * (B(1) * t3 + B(3) * t5) 
            stsVelSqr(i,2) = stsVelSqr(i,2)  
     &                  + r1 - DInv(2) * (B(2) * t4 + B(1) * t3) 
            stsVelSqr(i,3) = stsVelSqr(i,3)  
     &                  + r2 - DInv(3) * (B(3) * t5 + B(2) * t4) 

            stsVelSqr(i,4) = stsVelSqr(i,4) + t3 
            stsVelSqr(i,5) = stsVelSqr(i,5) + t4 
            stsVelSqr(i,6) = stsVelSqr(i,6) + t5 

            r0	        = res(4) 
            r1		= res(5) 
            r2		= res(6) 
            r3		= res(7) 
            r4		= res(8) 
            r5		= res(9) 
            
            r0          = DInv(1) * r0 
            r1          = DInv(2) * r1 
            r2          = DInv(3) * r2 

            r3          = r3 - B(1) * (r0 + r1) 
            r4          = r4 - B(2) * (r1 + r2) 
            r5          = r5 - B(3) * (r2 + r0) 

            t3          = CInv(1) * r3 + CInv(4) * r4 + CInv(6) * r5 
            t4          = CInv(4) * r3 + CInv(2) * r4 + CInv(5) * r5 
            t5          = CInv(6) * r3 + CInv(5) * r4 + CInv(3) * r5 

            stsStress(i,1) = stsStress(i,1)
     &                  + r0 - DInv(1) * (B(1) * t3 + B(3) * t5) 
            stsStress(i,2) = stsStress(i,2)
     &                  + r1 - DInv(2) * (B(2) * t4 + B(1) * t3) 
            stsStress(i,3) = stsStress(i,3)
     &                  + r2 - DInv(3) * (B(3) * t5 + B(2) * t4) 
            stsStress(i,4) = stsStress(i,4) + t3 
            stsStress(i,5) = stsStress(i,5) + t4 
            stsStress(i,6) = stsStress(i,6) + t5 
         enddo
      else if (stsType .eq. 2) then
         
         call evalAtInterp( yAlpha,     stsVec,         x, 
     &                      nResDims,   nshape)
         
         do i = 1, nshg
            
            u1    = stsVec(i,1)
            u2    = stsVec(i,2)
            u3    = stsVec(i,3)

            stsPres(i)    = stsPres(i)    +  stsVec(i,4) 
            stsPresSqr(i) = stsPresSqr(i) +  stsVec(i,4)*stsVec(i,4)  
            
            stsVel(i,1) = stsVel(i,1) + u1 
            stsVel(i,2) = stsVel(i,2) + u2 
            stsVel(i,3) = stsVel(i,3) + u3 

            stsVelSqr(i,1) = stsVelSqr(i,1) + u1*u1
            stsVelSqr(i,2) = stsVelSqr(i,2) + u2*u2
            stsVelSqr(i,3) = stsVelSqr(i,3) + u3*u3
            stsVelSqr(i,4) = stsVelSqr(i,4) + u1*u2
            stsVelSqr(i,5) = stsVelSqr(i,5) + u2*u3
            stsVelSqr(i,6) = stsVelSqr(i,6) + u3*u1
            
            stsStress(i,1) = stsStress(i,1) + stsVec(i,6)
            stsStress(i,2) = stsStress(i,2) + stsVec(i,7)
            stsStress(i,3) = stsStress(i,3) + stsVec(i,8)
            stsStress(i,4) = stsStress(i,4) + stsVec(i,9)
            stsStress(i,5) = stsStress(i,5) + stsVec(i,10)
            stsStress(i,6) = stsStress(i,6) + stsVec(i,11)

         enddo
      endif
      
      if ( mod(nTimeStep,stsWriteFreq) .eq. 0 .or. 
     &     nTimeStep .eq. nstep(itseq) ) then
         call stsWriteStats()
      endif

c$$$      if ( mod( nTimeStep, stsResetFreq) .eq. 0 ) then
c$$$         stsPres    = 0.0
c$$$         stsPresSqr = 0.0
c$$$         stsVel     = 0.0
c$$$         stsVelSqr  = 0.0
c$$$         stsVelReg  = 0.0
c$$$         stsStress  = 0.0
c$$$      
c$$$         nTimeStep  = 0
c$$$      endif
      
      return
      end
         
c-----------------------------------------------------------------------
c     collect the desired statistics 
c-----------------------------------------------------------------------
      subroutine stsWriteStats()
      
      use     stats
      include "common.h"

      integer      iofile, step2, itmp1, itmp,iofile2
      character*30 fname,  fmt1
      character*5  cname
      character*1  dash
      real*8       outvec(nshg,19)
c
c.... open the output file
c
      iofile = 39
      step2 = lstep+1  ! current time step
      itmp  = 1
      itmp1 = 1
      if (step1 .gt. 0) itmp  = int(log10(float(step1)))+1
      if (step2 .gt. 0) itmp1 = int(log10(float(step2)))+1
      dash = '-'
      write (fmt1,
     &     "('(''stats.'',i',i1,',a1,i',i1,',1x)')")
     &     itmp,itmp1
      write (fname,fmt1) step1,dash,step2
      fname = trim(fname) // cname(myrank+1)
      open ( unit = iofile, file = fname, status = 'unknown',
     &       form = 'unformatted')
c
c.... write the statistics
c
      outvec(:,1)     = stsPres(:)
      outvec(:,2:4)   = stsVel(:,:)
c      outvec(:,2:4)   = stsVelReg(:,:)
      outvec(:,5)     = zero   ! later will be temperature
      outvec(:,6)     = stsPresSqr(:)
      outvec(:,7:12)  = stsVelSqr(:,:)
      outvec(:,13)    = zero   ! later wil be tempSqr
      outvec(:,14:19) = stsStress(:,:)
      
      write (iofile) numnp, nshg, nTimeStep
      write (iofile) outvec(1:nshg,:)
      close (iofile)

c$$$      write (iofile) numnp, numnp, nTimeStep
c$$$      write (iofile) outvec(1:numnp,:)
c$$$      close (iofile)
      
      iofile2 = 40
c$$$      open (unit=iofile2,file='stats.asc',status='unknown')
c$$$c$$$c
c$$$c$$$c.... pressure, velocity, temperature
c$$$c$$$c
c$$$c$$$      write (iofile2) outvec(:,1:5)
c$$$c$$$c
c$$$c$$$c.... pressSqr, u_i u_j, tempSqr
c$$$c$$$c
c$$$c$$$      write (iofile2) outvec(:,6:13)
c$$$c$$$c
c$$$c$$$c.... viscous stress
c$$$c$$$c
c$$$c$$$      write (iofile2) outvec(:,14:19)
c$$$c$$$
c$$$
c$$$c
c$$$c.... write the statistics
c$$$c
c$$$      write (iofile2,*) 'nNodes ',numnp
c$$$      write (iofile2,*) 'power ',2
c$$$      write (iofile2,*) 'nTimeSteps = ', nTimeStep
c$$$c
c$$$c.... velocity
c$$$c      
c$$$      write (iofile2,*) 'vel 3 ',numnp
c$$$      do i = 1, numnp
c$$$         write (iofile2,111) stsVel(i,1), stsVel(i,2), stsVel(i,3)
c$$$      enddo
c$$$c
c$$$c.... pressure
c$$$c
c$$$      write (iofile2,*) 'pres 1 ',numnp
c$$$      do i = 1, numnp
c$$$         write (iofile2,112) stsPres(i)
c$$$      enddo
c$$$c
c$$$c.... velSqr
c$$$c
c$$$      write (iofile2,*) 'velSqr 6 ',numnp
c$$$      do i = 1, numnp
c$$$         write (iofile2,113) (stsVelSqr(i,j),j=1,6)
c$$$      enddo
c$$$c
c$$$c.... presSqr
c$$$c
c$$$      write (iofile2,*) 'presSqr 1 ',numnp
c$$$      do i = 1, numnp
c$$$         write (iofile2,112) stsPresSqr(i)
c$$$      enddo
c$$$c
c$$$c.... stress      
c$$$c
c$$$      write (iofile2,*) 'stress 6 ',numnp
c$$$      do i = 1, numnp
c$$$         write (iofile2,113) (stsStress(i,j),j=1,6)
c$$$      enddo
c$$$c
c$$$c.... velocity
c$$$c
c$$$      write (iofile2,*) 'vel 3 ',numnp
c$$$      do i = 1, numnp
c$$$         write (iofile2,111) (stsVelReg(i,j),j=1,3)
c$$$      enddo
c$$$      
c$$$      close (iofile2)

 111  format(1p,3e24.16)
 112  format(1p, e24.16)
 113  format(1p,6e24.16)
      
      return
      end
