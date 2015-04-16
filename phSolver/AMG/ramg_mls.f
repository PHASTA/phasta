!**************************************************************
!      mls_eigen: get maximum eigenvalue of A
!      mls_setup: set coefficients of A^n
!      
!      These routine uses ARPACK and auxilary routines in 
!      ramg_ggb.f
!**************************************************************
      subroutine ramg_mls_eigen(evmax,level,flagfb,
     &                          colm,rowp,lhsK,lhsP,
     &                          ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      include 'debug.h'

      real(kind=8),intent(inout)                      :: evmax
      integer,intent(in)                              :: level
      integer,intent(in)                              :: flagfb
      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp

      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP

      integer,intent(in),dimension(nlwork)            :: ilwork
      integer,intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)  :: BC

!      ! parameter for PPE Ap-product
!      real(kind=8),dimension(nshg,4) :: diag
      ! parameter for eigenvalue 

      integer   iparam(11),ipntr(14)
      integer gcomm

!      logical   selectncv(maxncv)
!      real(kind=8) :: d(maxncv,3),dr(maxncv),di(maxncv),
!     &                workev(3*maxncv),
!     &                workl(3*maxncv*maxncv+6*maxncv)
      logical,dimension(:),allocatable :: selectncv
      real(kind=8),dimension(:,:),allocatable :: d
      real(kind=8),dimension(:),allocatable :: dr,di,workev,workl

      real(kind=8),allocatable,dimension(:) ::  resid,workd

      character  bmat*1,which*2
      integer   :: ido,n,nx,nev,ncv,lworkl,info,ierr,
     &             i,j,ishifts,maxitr,model,nconv
      real(kind=8)  :: tol,sigmar,sigmai
      logical     ::  first,rvec

      real(kind=8) :: res_i,res_o,rtemp,tnorm,tnorm1
      real(kind=8),allocatable,dimension(:,:):: tramg_ev
      real(kind=8),dimension(amg_nshg(level)) :: twork1,twork2

      integer,allocatable,dimension(:) :: gmap,grevmap

      integer :: p,q,nsize,step,asize,gsize

!      ! if ( freeze parameter ) return
!      call drvLesPrepDiag(diag,ilwork,iBC,BC,iper,rowp,colm,lhsK,lhsP)

      asize = amg_nshg(level)
      allocate(gmap(asize))
      allocate(grevmap(asize))

      call ramg_generate_gmap(ilwork,asize,nsize,gmap,grevmap,level)

      call MPI_AllReduce(nsize,gsize,1,MPI_INTEGER,MPI_MAX,
     &                   MPI_COMM_WORLD,ierr)

      ncv = min(maxncv,gsize)
      !write(*,*)'mcheck ncv=',ncv

      if (ncv.eq.1) then
         evmax = 1
         return
      endif

      allocate(selectncv(ncv))
      allocate(d(ncv,3))
      allocate(dr(ncv))
      allocate(di(ncv))
      allocate(workev(3*ncv))
      allocate(workl(3*ncv*ncv+6*ncv))

      allocate(resid(gsize))
      allocate(workd(3*gsize))
      !***********   Start of Parameter setting *******
      
      ndigit = -3
      logfil = 6
      mnaitr = 0
      mnapps = 0
      mnaupd = 0 ! controls the output (verbose) mode
      mnaup2 = 0
      mneigh = 0
      mneupd = 0

      nev = 1
      bmat = 'I'
      which = 'LM'

      lworkl = 3*ncv**2+6*ncv
      tol = 1.0E-3
      ido = 0
      info = 0

      iparam(1) = 1
      iparam(3) = 500
      iparam(7) = 1

      !********** End of parameter setting ********

      allocate(tramg_ev(nsize,maxncv))

      step = 1
      tramg_ev = 0
      gcomm = MPI_COMM_WORLD

      call pdnaupd(gcomm,
     &               ido,bmat,nsize,which,nev,tol,resid,ncv,
     &               tramg_ev,
     &               nsize,iparam,ipntr,workd,workl,lworkl,
     &               info)
      do while ((ido.eq.1).or.(ido.eq.-1))
        call ramg_ggb_G2A(twork1,workd(ipntr(1)),asize,nsize,1,
     &                    gmap,grevmap)
!        call ramg_PPEAp(twork2,twork1,diag,
!     &                colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
!        call ramg_calcAv(amg_A_colm(level)%p,amg_A_rowp(level)%p,
!     &                     amg_A_lhs(level)%p,amg_nshg(level),
!     &                     amg_nnz(level),twork2,twork1,1)
        if (flagfb.eq.1) then
        call ramg_calcAv_g(level,twork2,twork1,colm,rowp,lhsK,lhsP,
     &                 ilwork,BC,iBC,iper,1)
        else 
        call ramg_mls_calcPAv(level,twork2,twork1,
     &           colm,rowp,lhsK,lhsP,
     &           ilwork,BC,iBC,iper)
        endif
        call ramg_ggb_A2G(twork2,workd(ipntr(2)),asize,nsize,1,
     &                    gmap,grevmap)
        call pdnaupd(gcomm,
     &               ido,bmat,nsize,which,nev,tol,resid,ncv,
     &               tramg_ev,
     &               nsize,iparam,ipntr,workd,workl,lworkl,
     &               info)
         step = step + 1
      enddo
      if (info.lt.0) then
          if (myrank.eq.master) then
          write(*,*)'mcheck: mls: info:',info,iparam(5)
          endif
          ramg_levelx = level-1
      else
          rvec = .false.
          call pdneupd (gcomm,rvec,'A',selectncv,dr,di,
     &         tramg_ev,nsize,
     &         sigmar,sigmai,workev,bmat,nsize,which,nev,tol,
     &         resid,ncv,tramg_ev,
     &         nsize,iparam,ipntr,workd,workl,
     &         lworkl,ierr) 
      end if
      evmax = dr(1)
!      if (myrank.eq.master) then
!         write(*,*)'MLS: largest eigenvalue of A at level ',level,
!     &    ' is ',evmax
!      endif
      !write(*,*)'mcheck: ggb over pdnaupd'

      deallocate(tramg_ev)
      deallocate(gmap)
      deallocate(grevmap)
      deallocate(resid)
      deallocate(workd)

      deallocate(selectncv)
      deallocate(d)
      deallocate(dr)
      deallocate(di)
      deallocate(workev)
      deallocate(workl)
   
      end subroutine ! ramg_mls_eigen

      subroutine ramg_mls_min_eigen(evmax,level,flagfb,
     &                          colm,rowp,lhsK,lhsP,
     &                          ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      include 'debug.h'

      real(kind=8),intent(inout)                      :: evmax
      integer,intent(in)                              :: level
      integer,intent(in)                              :: flagfb
      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp

      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP

      integer,intent(in),dimension(nlwork)            :: ilwork
      integer,intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)  :: BC

!      ! parameter for PPE Ap-product
!      real(kind=8),dimension(nshg,4) :: diag
      ! parameter for eigenvalue 

      integer   iparam(11),ipntr(14)
      integer gcomm

!      logical   selectncv(maxncv)
!      real(kind=8) :: d(maxncv,3),dr(maxncv),di(maxncv),
!     &                workev(3*maxncv),
!     &                workl(3*maxncv*maxncv+6*maxncv)
      logical,dimension(:),allocatable :: selectncv
      real(kind=8),dimension(:,:),allocatable :: d
      real(kind=8),dimension(:),allocatable :: dr,di,workev,workl

      real(kind=8),allocatable,dimension(:) ::  resid,workd

      character  bmat*1,which*2
      integer   :: ido,n,nx,nev,ncv,lworkl,info,ierr,
     &             i,j,ishifts,maxitr,model,nconv
      real(kind=8)  :: tol,sigmar,sigmai
      logical     ::  first,rvec

      real(kind=8) :: res_i,res_o,rtemp,tnorm,tnorm1
      real(kind=8),allocatable,dimension(:,:):: tramg_ev
      real(kind=8),dimension(amg_nshg(level)) :: twork1,twork2

      integer,allocatable,dimension(:) :: gmap,grevmap

      integer :: p,q,nsize,step,asize,gsize

!      ! if ( freeze parameter ) return
!      call drvLesPrepDiag(diag,ilwork,iBC,BC,iper,rowp,colm,lhsK,lhsP)

      asize = amg_nshg(level)
      allocate(gmap(asize))
      allocate(grevmap(asize))

      call ramg_generate_gmap(ilwork,asize,nsize,gmap,grevmap,level)

      call MPI_AllReduce(nsize,gsize,1,MPI_INTEGER,MPI_MAX,
     &                   MPI_COMM_WORLD,ierr)

      ncv = min(maxncv,gsize)
      !write(*,*)'mcheck ncv=',ncv

      if (ncv.eq.1) then
         evmax = 1
         return
      endif

      allocate(selectncv(ncv))
      allocate(d(ncv,3))
      allocate(dr(ncv))
      allocate(di(ncv))
      allocate(workev(3*ncv))
      allocate(workl(3*ncv*ncv+6*ncv))

      allocate(resid(gsize))
      allocate(workd(3*gsize))
      !***********   Start of Parameter setting *******
      
      ndigit = -3
      logfil = 6
      mnaitr = 0
      mnapps = 0
      mnaupd = 0 ! controls the output (verbose) mode
      mnaup2 = 0
      mneigh = 0
      mneupd = 0

      nev = 1
      bmat = 'I'
      which = 'SM'

      lworkl = 3*ncv**2+6*ncv
      tol = 1.0E-3
      ido = 0
      info = 0

      iparam(1) = 1
      iparam(3) = 500
      iparam(7) = 1

      !********** End of parameter setting ********

      allocate(tramg_ev(nsize,maxncv))

      step = 1
      tramg_ev = 0
      gcomm = MPI_COMM_WORLD

      call pdnaupd(gcomm,
     &               ido,bmat,nsize,which,nev,tol,resid,ncv,
     &               tramg_ev,
     &               nsize,iparam,ipntr,workd,workl,lworkl,
     &               info)
      do while ((ido.eq.1).or.(ido.eq.-1))
        call ramg_ggb_G2A(twork1,workd(ipntr(1)),asize,nsize,1,
     &                    gmap,grevmap)
!        call ramg_PPEAp(twork2,twork1,diag,
!     &                colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
!        call ramg_calcAv(amg_A_colm(level)%p,amg_A_rowp(level)%p,
!     &                     amg_A_lhs(level)%p,amg_nshg(level),
!     &                     amg_nnz(level),twork2,twork1,1)
        if (flagfb.eq.1) then
        call ramg_calcAv_g(level,twork2,twork1,colm,rowp,lhsK,lhsP,
     &                 ilwork,BC,iBC,iper,1)
        else 
        call ramg_mls_calcPAv(level,twork2,twork1,
     &           colm,rowp,lhsK,lhsP,
     &           ilwork,BC,iBC,iper)
        endif
        call ramg_ggb_A2G(twork2,workd(ipntr(2)),asize,nsize,1,
     &                    gmap,grevmap)
        call pdnaupd(gcomm,
     &               ido,bmat,nsize,which,nev,tol,resid,ncv,
     &               tramg_ev,
     &               nsize,iparam,ipntr,workd,workl,lworkl,
     &               info)
         step = step + 1
      enddo
      if (info.lt.0) then
          if (myrank.eq.master) then
          write(*,*)'mcheck: mls: info:',info,iparam(5)
          endif
          ramg_levelx = level-1
      else
          rvec = .false.
          call pdneupd (gcomm,rvec,'A',selectncv,dr,di,
     &         tramg_ev,nsize,
     &         sigmar,sigmai,workev,bmat,nsize,which,nev,tol,
     &         resid,ncv,tramg_ev,
     &         nsize,iparam,ipntr,workd,workl,
     &         lworkl,ierr) 
      end if
      evmax = dr(1)
!      if (myrank.eq.master) then
!         write(*,*)'MLS: smallest eigenvalue of A is ',evmax
!      endif
      !write(*,*)'mcheck: ggb over pdnaupd'

      deallocate(tramg_ev)
      deallocate(gmap)
      deallocate(grevmap)
      deallocate(resid)
      deallocate(workd)

      deallocate(selectncv)
      deallocate(d)
      deallocate(dr)
      deallocate(di)
      deallocate(workev)
      deallocate(workl)
   
      end subroutine ! ramg_mls_min_eigen

      subroutine ramg_mls_coeff(level,ilwork,BC,iBC,iper)
      end subroutine !ramg_mls_coeff


      subroutine ramg_mls_setup(colm,rowp,lhsK,lhsP,
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

      
      real(kind=8)    :: evmax,ddeg,aux1,aux0,auxom,rho,rho2
      real(kind=8),dimension(10) :: omloc
      real(kind=8)    :: mlsOver  ! magic number! 
      integer         :: i,j,k,lvl,ideg
      integer  :: nSample,nGrid
      real(kind=8) :: gridStep,coord,samplej

      if (ramg_setup_flag .ne. 0 ) return
      if ((iamg_smoother.eq.1).and.(numpe.eq.1)) return
      ! do not setup if serial Gauss-Seidel
      ! magic numbers....
      mlsOver = 1.1

      do lvl=1,ramg_levelx
      call ramg_mls_eigen(evmax,lvl,1,colm,rowp,lhsK,lhsP,
     &                    ilwork,BC,iBC,iper)

      mlsCf(lvl,:) = 0
      omloc = 0

      if (iamg_smoother.eq.1) then !Gauss-Seidel in parallel
          ideg = 1
      else
        ideg = iabs(mlsDeg)
      endif
      ddeg = ideg
      aux1 = 1.0/(2.0*ddeg+1.0)
      rho = evmax!*mlsOver

      do i=1,ideg
         aux0 = 2.0*i*PI
         auxom = rho/2.0*(1.0-cos(aux0*aux1))
         omloc(i) = 1.0/auxom
         !write(*,*)'mcheck: ',omloc(i)
      enddo

      ! 4th order maximum
      mlsCf(lvl,1) = omloc(1)+omloc(2)+omloc(3)+omloc(4)+omloc(5)
      
      mlsCf(lvl,2) = -( omloc(1)*omloc(2)+omloc(1)*omloc(3)
     &             +omloc(1)*omloc(4)+omloc(1)*omloc(5)
     &             +omloc(2)*omloc(3)+omloc(2)*omloc(4)
     &             +omloc(2)*omloc(5)+omloc(3)*omloc(4)
     &             +omloc(3)*omloc(5)+omloc(4)*omloc(5) )

      mlsCf(lvl,3) = ( omloc(1)*omloc(2)*omloc(3)
     &            +omloc(1)*omloc(2)*omloc(4)
     &            +omloc(1)*omloc(2)*omloc(5)
     &            +omloc(1)*omloc(3)*omloc(4)
     &            +omloc(1)*omloc(3)*omloc(5)
     &            +omloc(1)*omloc(4)*omloc(5)
     &            +omloc(2)*omloc(3)*omloc(4)
     &            +omloc(2)*omloc(3)*omloc(5)
     &            +omloc(2)*omloc(4)*omloc(5)
     &            +omloc(3)*omloc(4)*omloc(5) )

      mlsCf(lvl,4) = -( omloc(1)*omloc(2)*omloc(3)*omloc(4)
     &             +omloc(1)*omloc(2)*omloc(3)*omloc(5)
     &             +omloc(1)*omloc(2)*omloc(4)*omloc(5)
     &             +omloc(1)*omloc(3)*omloc(4)*omloc(5)
     &             +omloc(2)*omloc(3)*omloc(4)*omloc(5) )

      mlsCf(lvl,5) = omloc(1)*omloc(2)*omloc(3)*omloc(4)*omloc(5)

      do i=1,ideg
         mlsOm(lvl,i) = omloc(i)
      enddo

      ! setup coefficients for post smoothing
      ! ref: trilinos ml

      nSample=20000
      gridStep = rho/nSample
      rho2=zero
      do j=1,nSample
         coord=j*gridStep
         samplej=1.0-omloc(1)*coord
         do i=2,ideg
            samplej=samplej*(1.0-omloc(i)*coord)
         enddo
         samplej=samplej*(samplej*coord)
         if (samplej.gt.rho2) rho2=samplej
      enddo
      smlsOm(lvl,1) = 2.0/(1.025*rho2)!*mlsCf(lvl,1)
      !write(*,*)'mcheckpost:',rho2,smlsOm(lvl,1)
    
      !do i=1,mlsDeg
      !   write(*,*)'mcheck: ',lvl,i,mlsCf(lvl,i)
      !enddo
      if (myrank.eq.master) then
          write(*,*)'MLS: Setup finished at level:',lvl,mlsCf(lvl,1)
      endif
      enddo !lvl

      end subroutine ! ramg_mls_setup

!***********************************************************
!     Construct Sfc = Ifc*(I-\lambda_1 A-\lambda_2 A^2)
!      
!***********************************************************      
!      Additional note: This is useless which results in
!      a super dense S instead of I, using this interpolator
!      is not effective won't save time ignore!
      
!***********************************************************
!     u = (S^2*A) v
!     S = polynomial smoothing      
!      
!**********************************************************
      subroutine ramg_mls_calcPAv(level,u,v,colm,rowp,lhsK,lhsP,
     &                            ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      
      integer,intent(in),dimension(nlwork)  :: ilwork
      integer,intent(in),dimension(nshg)    :: iBC,iper
      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp

      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP

      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: v
      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: u

      real(kind=8),dimension(amg_nshg(level)) :: aux,aux2,r0

!      r0 = 0
      call ramg_calcAv_g(level,aux,v,colm,rowp,lhsK,lhsP,
     &                   ilwork,BC,iBC,iper,1)
      call ramg_mls_sandw_pre(aux2,aux,level,colm,rowp,lhsK,lhsP,
     &                    ilwork,BC,iBC,iper)!,
      call ramg_mls_sandw_post(u,aux2,level,colm,rowp,lhsK,lhsP,
     &                    ilwork,BC,iBC,iper)!,
     
      end subroutine ! ramg_mls_calcPav

      subroutine ramg_mls_setupPost(colm,rowp,lhsK,lhsP,
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

      real(kind=8)    :: evmax,ddeg,aux1,aux0,auxom,rho
      real(kind=8),dimension(10) :: omloc
      integer :: i,j,k,lvl

      do lvl = 1,ramg_levelx
      call ramg_mls_eigen(evmax,lvl,2,colm,rowp,lhsK,lhsP,
     &                    ilwork,BC,iBC,iper)
      smlsOm(lvl,2) = 1.0/(1.025*evmax)!*mlsCf(lvl,1)
      !mlsOm(lvl,2) = -1.0*evmax*mlsCf(lvl,1)*mlsCf(lvl,1)
      write(*,*)'mcheck: post,',evmax,smlsOm(lvl,2)
      enddo
      
      end subroutine ! ramg_mls_setupPost

!***********************************************************
!      ramg_mls_apply: u = MLS(v,(lhs,r))
!      v : approx. solution before
!      u : approx. solution (after)
!      r : residual vector, r=Ax-b
!***********************************************************

      subroutine ramg_mls_apply(u,v,r,level,colm,rowp,lhsK,lhsP,
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
      
      real(kind=8),dimension(amg_nshg(level)) :: aux

!     calculate residual first if necessary
!     do not calculate residual at pre-smoothing, 
!     only calculate it at post-smoothing.
      if (fwdbck.eq.2) then
          aux = -r
      else
          call ramg_calcAv_g(level,aux,v,colm,rowp,lhsK,lhsP,
     &         ilwork,BC,iBC,iper,1)
          aux = aux-r
      endif
!      if (.false.) then ! always post smoother
      if (.true.) then ! always pre smoother
          call ramg_mls_apply_fwd(u,v,aux,level,colm,rowp,lhsK,lhsP,
     &                  ilwork,BC,iBC,iper)
      else
          call ramg_mls_apply_post(u,v,aux,level,colm,rowp,lhsK,lhsP,
     &                  ilwork,BC,iBC,iper)
      endif
      
      u = v-1.1*u
 
      end subroutine ! ramg_mls_apply


      subroutine ramg_mls_apply_fwd(u,v,r,level,colm,rowp,lhsK,lhsP,
     &                          ilwork,BC,iBC,iper)
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

      ! Local 
      real(kind=8),dimension(amg_nshg(level)) :: pAux,res
      real(kind=8) :: cf
      integer :: i,j,k

      pAux = r
      ! c_1*A*r
      
      cf = mlsCf(level,1)
      u = cf*pAux
      do i=2,mlsDeg
         ! c_i*A*(A*...*r)
         call ramg_calcAv_g(level,res,pAux,colm,rowp,lhsK,lhsP,
     &        ilwork,BC,iBC,iper,1)
         pAux = res
         cf = mlsCf(level,i)
         u = u+cf*res
      enddo

      end subroutine ! ramg_mls_apply_fwd

      subroutine ramg_mls_apply_post(u,v,r,level,colm,rowp,lhsK,lhsP,
     &                          ilwork,BC,iBC,iper)
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

      ! Local 
      real(kind=8),dimension(amg_nshg(level)) :: pAux,y,res,y1
      real(kind=8) :: cf
      integer :: i,j,k

      call ramg_mls_apply_fwd(y1,v,r,level,colm,rowp,lhsK,lhsP,
     &     ilwork,BC,iBC,iper)
!      y1=1.1*y1
      call ramg_calcAv_g(level,pAux,y1,colm,rowp,lhsK,lhsP,
     &     ilwork,BC,iBC,iper,1)
      res = b-pAux!b-pAux ! rhs to feed in post smoothing, keep y

      call ramg_mls_sandw_post(pAux,res,level,colm,rowp,lhsK,lhsP,
     &     ilwork,BC,iBC,iper)
      call ramg_mls_sandw_pre(y,pAux,level,colm,rowp,lhsK,lhsP,
     &    ilwork,BC,iBC,iper)
      u = smlsOm(level,1)*y+y1

      end subroutine ! ramg_mls_apply_post

      subroutine ramg_mls_sandw_pre(u,v,level,colm,rowp,lhsK,lhsP,
     &                          ilwork,BC,iBC,iper)
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


      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: v
      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: u

      ! Local 
      real(kind=8),dimension(amg_nshg(level)) :: res
      real(kind=8) :: cf
      integer :: i,j,k

      u = v
      ! c_1*A*r

      do i=mlsDeg,1,-1
         call ramg_calcAv_g(level,res,u,colm,rowp,lhsK,lhsP,
     &        ilwork,BC,iBC,iper,1)
         u = u-mlsOm(level,i)*res
      enddo
      
      end subroutine ! ramg_mls_sandw_pre

      subroutine ramg_mls_sandw_post(u,v,level,colm,rowp,lhsK,lhsP,
     &                          ilwork,BC,iBC,iper)
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


      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: v
      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: u

      ! Local 
      real(kind=8),dimension(amg_nshg(level)) :: res
      real(kind=8) :: cf
      integer :: i,j,k

      u = v
      ! c_1*A*r

      do i=1,mlsDeg,1
         call ramg_calcAv_g(level,res,u,colm,rowp,lhsK,lhsP,
     &        ilwork,BC,iBC,iper,1)
         u = u-mlsOm(level,i)*res
      enddo
      
      end subroutine ! ramg_mls_sandw_post

