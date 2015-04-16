!***********************************************************
!      ramg_cheby_apply: u = Cheby(v,(lhs,r))
!      v : approx. solution before
!      u : approx. solution (after)
!      r : residual vector, r=Ax-b
!***********************************************************


      subroutine ramg_cheby_apply(u,v,r,level,colm,rowp,lhsK,lhsP,
     &                          ilwork,BC,iBC,iper,fwdbck)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in),dimension(nlwork)            :: ilwork
      integer,intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)  :: BC
      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp

      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP


      real(kind=8),intent(in),dimension(amg_nshg(level)) :: v,r
      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: u
      integer,intent(in) :: fwdbck

      ! Local 
      real(kind=8),dimension(amg_nshg(level)) :: pAux
      real(kind=8),dimension(amg_nshg(level)) :: epdk,epx

      real(kind=8),dimension(amg_nshg(level)) :: cy1,cz1

      real(kind=8) :: cf,beta,alpha,delta,theta,s1,rhok,invtheta
      real(kind=8) :: rhokp1,dtemp1,dtemp2,twoodelta
      real(kind=8) :: chebyniu,chebyomega,chebyck1,chebyck2,chebygamma
      integer :: i,j,k

      beta = 1.1*mlsCf(level,1) ! 1.1 * Ev_max
      alpha = mlsCf(level,1)*ramg_chebyratio

      IF (.FALSE.) THEN
      ! mine ( matrix computing )
      chebyniu = 1.0+2.0*(1.0-beta)/(beta-alpha)
      chebygamma = 2.0/(2.0-alpha-beta)
      chebyomega = chebyniu/(2.0*chebyniu*chebyniu-1.0)
      chebyomega = chebyomega*2.0*(2.0-beta-alpha)/(beta-alpha)
      if (fwdbck.eq.2) then ! initial guess x0 = 0
          cy1 = r
      else
          call ramg_calcAv_g(level,pAux,v,colm,rowp,lhsK,lhsP,
     &         ilwork,BC,iBC,iper,1)
          cy1 = v-1.1*3.0/beta*(pAux-r)
      endif

      call ramg_calcAv_g(level,pAux,cy1,colm,rowp,lhsK,lhsP,
     &                   ilwork,BC,iBC,iper,1)
      cz1 = r-pAux

      if (fwdbck.eq.2) then
          u = chebyomega*(cy1+chebygamma*cz1)
      else
          u = chebyomega*(cy1-v+chebygamma*cz1)+v
      endif
      
      ELSE
      ! trilinos
      delta = (beta-alpha)/2.0
      twoodelta = 4.0/(beta-alpha)
      theta = (beta+alpha)/2.0
      invtheta = 2.0/(beta+alpha)
      s1 = 2.0*theta/delta
      rhok = delta/theta

      !write(*,*)delta,twoodelta,2.0/delta
      if (fwdbck.eq.2) then ! pre-smoothing, initial guess v=0
          epdk = r*invtheta!/theta
          epx  = epdk
      else
          call ramg_calcAv_g(level,pAux,v,colm,rowp,lhsK,lhsP,
     &         ilwork,BC,iBC,iper,1)
          epdk = r-pAux
          epdk = epdk*invtheta!/theta
          epx = v+epdk
      endif

      do k=1,mlsDeg-1
          call ramg_calcAv_g(level,pAux,epx,colm,rowp,lhsK,lhsP,
     &         ilwork,BC,iBC,iper,1)
          rhokp1 = 1.0/(s1-rhok)
          dtemp1 = rhokp1*rhok
          dtemp2 = rhokp1*twoodelta
          rhok = rhokp1
          epdk = dtemp1*epdk+dtemp2*(r-pAux)
          epx = epx + epdk
      enddo
      u = epx

      ENDIF

      end subroutine ! ramg_cheby_apply


      subroutine ramg_cheby_setup(colm,rowp,lhsK,lhsP,
     &                          ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in),dimension(nshg+1)            :: colm
      integer,intent(in),dimension(nnz_tot)           :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP
      integer,intent(in),dimension(nlwork)            :: ilwork
      integer,intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)  :: BC

      integer :: i,j,k,lvl,ideg

      if (ramg_setup_flag .ne. 0 ) return

      if (myrank.eq.master) then
         write(*,*)'Start Chebyshev setup'
      endif
      do lvl=1,ramg_levelx
      call ramg_mls_eigen(mlsCf(lvl,1),lvl,1,colm,rowp,lhsK,lhsP,
     &                    ilwork,BC,iBC,iper)
!      call ramg_mls_min_eigen(mlsCf(lvl,2),lvl,1,colm,rowp,lhsK,lhsP,
!     &                    ilwork,BC,iBC,iper)
!      mlsCf(lvl,2) = mlsCf(lvl,1)
      if (myrank.eq.master) then
          write(*,7000)lvl,mlsCf(lvl,1)
7000      format('Level: ',I2,T15,'Ev_max: ',E10.4)
      endif
      enddo
      if (myrank.eq.master) then
          write(*,*)'End Chebyshev setup'
      endif

      end subroutine !ramg_cheby_setup
