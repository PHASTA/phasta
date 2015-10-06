c-----------------------------------------------------------------------
c
c  Natural pressure boundary condition can be calculated with p, the pressure,
c  related (in some prescribed manner) to Q, the flow rate, through the same 
c  boundary.  To do this efficiently requires us to precompute the integral 
c  of N_A over the boundary for each node A and store it in a vector of length
c  nshg (a bit wasteful since only the nodes on the boundary will be non zero
c  in this vector but it is probably slower to index it than to multiply and 
c  add the extra zeros....check later).
c
c-----------------------------------------------------------------------
      module pvsQbi

      real*8, allocatable ::  NABI(:,:)
      real*8, allocatable ::  NASC(:)
      integer, allocatable :: ndsurf(:)
      
      end module

      subroutine finalizeNABI
        use pvsQbi
        if( allocated(NABI) ) then
          deallocate(NABI)
        endif
        if( allocated(NASC) ) then
          deallocate(NASC)
        endif
        if( allocated(ndsurf) ) then
          deallocate(ndsurf)
        endif
      end

c-----------------------------------------------------------------------
c
c     Initialize:
c
c-----------------------------------------------------------------------
      subroutine initNABI( x, shpb )
      
      use     pointer_data
      use     pvsQbi
      include "common.h"
      
      real*8   x(numnp,nsd)
c
c use is like
c 
c      NABI=pvsQbi -> NABI
c
        dimension   shpb(MAXTOP,maxsh,MAXQPT)
        real*8, allocatable :: tmpshpb(:,:)
        allocate ( NABI(nshg,3) ) 
        allocate ( NASC(nshg)   )
        allocate ( ndsurf(nshg) ) 

c
c....  calculate NABI
c
      NABI=zero
      NASC=zero
      ndsurf=0
c
c.... -------------------->   boundary elements   <--------------------
c
c.... set up parameters
c
c        intrul = intg   (2,itseq)
c        intind = intptb (intrul)
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel 


          if(lcsyst.eq.3) lcsyst=nenbl
c
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          allocate (tmpshpb(nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)

          call AsBNABI (                       x,
     &                 tmpshpb,
     &                 mienb(iblk)%p,
     &                 miBCB(iblk)%p)

          call AsBNASC(                       x,
     &                 tmpshpb,
     &                 mienb(iblk)%p,
     &                 miBCB(iblk)%p)

          deallocate (tmpshpb)

      enddo 

c
c     note that NABI has NOT been communicated.  It
C     is the on processor contribution to this vector.  It will used to
C     build the on processor contribution to res that will then be made
C     complete via a call to commu.  Similarly the LHS usage will create
C     the on-processor contribution to the lhsK. Same for NASC
c
      return
      end


