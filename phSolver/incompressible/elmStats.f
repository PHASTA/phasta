      subroutine elmStatsLhs( x,  iBC,   iper,  ilwork )
c-----------------------------------------------------------------------
c     compute the necessary terms for the statistics projection
c     matrices.
c-----------------------------------------------------------------------
      use     stats
      use     pointer_data
      
      include "common.h"

      real*8  x(numnp,3)
      integer iBC(nshg), iper(nshg), ilwork(nlwork)
      
      real*8, allocatable :: xl(:,:,:)
      real*8, allocatable :: lStsVec(:,:,:)

c
c.... loop over element blocks
c
      stsVec = zero
      
      do iblk = 1, nelblk
         iel    = lcblk(1,iblk)
         lcsyst = lcblk(3,iblk)
         nenl   = lcblk(5,iblk) ! no. of vertices per element
         nshl   = lcblk(10,iblk)
         ndofl  = lcblk(8,iblk)
         npro   = lcblk(1,iblk+1) - iel 

         allocate ( xl(npro,nenl,3)             )
         allocate ( lStsVec(npro,nshl,nResDims) )
c
c.... localize needed data
c
         call localx ( x,    xl,  mien(iblk)%p, nsd,   'gather  ' )
c
c.... form the Lhs
c
         call e3StsLhs( xl, lStsVec )
c
c.... assemble
c
         call local (stsVec, lStsVec, mien(iblk)%p,
     &               nResDims, 'scatter ' ) 

         deallocate ( xl       )
         deallocate ( lStsVec  )
c
c.... end loop over element blocks
c
      enddo

      if (numpe > 1) then
        call commu (stsVec, ilwork, nResDims  , 'in ')
      endif
c
c.... local periodic boundary conditions (no communications)
c
      do j = 1,nshg
         if (btest(iBC(j),10)) then
            i = iper(j)
            stsVec(i,:) = stsVec(i,:) + stsVec(j,:)
         endif
      enddo
c
      do i = 1,nshg
         stsVec(i,:) = stsVec(iper(i),:)
      enddo
      if (numpe > 1) then
        call commu (stsVec, ilwork, nResDims  , 'out')
      endif

      return
      end
      
c-----------------------------------------------------------------------
c  Assemble the residual for the statistics
c-----------------------------------------------------------------------
      subroutine elmStatsRes( y,        ac,    x,      shp,     shgl, 
     &                        shpb,     shglb,       iBC,     BC, 
     &                        iper,     ilwork,      rowp,    colm,
     &                        lhsK,     lhsP )
      
      use     stats
      include "common.h"
      
      
      real*8  y(nshg,ndof),             ac(nshg,ndof), x(numnp,nsd),
     &        shp(MAXTOP,maxsh,MAXQPT),  shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &        shpb(MAXTOP,maxsh,MAXQPT),
     &        shglb(MAXTOP,nsd,maxsh,MAXQPT),
     &        BC(nshg,ndofBC),          lhsK(9,nnz_tot),
     &        lhsP(4,nnz_tot),         res(nshg,ndof)

      integer iBC(nshg),                iper(nshg),
     &        ilwork(nlwork),           rowp(nshg,nnz),
     &        colm(nshg+1)
      

      lhs    = 0
      stsVec = zero
      
      stsResFlg = 1
      ierrcalctmp=ierrcalc ! save the current value of ierrcalc
      ierrcalc=0           ! set it to zero so that we don't calculate error
                           ! here (which would overflow memory around rjunk)
      call ElmGMR (y,         ac,         x,
     &             shp,       shgl,       iBC,       
     &             BC,        shpb,       shglb,
     &             res,       iper,       ilwork,   
     &             rowp,      colm,       lhsK,      
     &             lhsP,      rjunk   )
      stsResFlg = 0
      ierrcalc=ierrcalctmp  ! reset it back to original value
      return 
      end

