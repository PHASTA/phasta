!************************************************************
!      ramg_calcPPEv
!      calculate u = PPE*v ( parallelly using lesLib calls)
!************************************************************
      subroutine ramg_calcPPEv(colm,rowp,lhsK,lhsP,
     &                   ilwork,BC,iBC,iper,
     &                         u,v)
      use ramg_data
      include "common.h"
      include "mpif.h"
      
      integer,intent(in),dimension(nshg+1) :: colm
      integer,intent(in),dimension(nnz_tot) :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot) :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot) :: lhsP
      integer,intent(in),dimension(nlwork)  :: ilwork
      integer,intent(in),dimension(nshg)    :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      
      real(kind=8),intent(in),dimension(nshg) :: v
      real(kind=8),intent(inout),dimension(nshg) :: u

      real(kind=8),dimension(nshg,3) :: utmp
      real(kind=8),dimension(nshg,4) :: utmp4
      integer                         :: i,j

      call commOut(v,ilwork,1,iper,iBC,BC)
      call fLesSparseApG(colm,rowp,lhsP,v,utmp,nshg,nnz_tot)
      call commIn(utmp,ilwork,3,iper,iBC,BC)

      do i=1,3
         utmp(:,i) = utmp(:,i)*ramg_flowDiag%p(:,i)**2
      enddo

      do i=1,3
         utmp4(:,i) = -utmp(:,i)
      enddo

      utmp4(:,4) = v

      call commOut(utmp4,ilwork,4,iper,iBC,BC)
      call fLesSparseApNGtC(colm,rowp,lhsP,utmp4,u,nshg,nnz_tot)
      call commIn(u,ilwork,1,iper,iBC,BC)
!     There is a slight modify at lesSparse.f:
!     row(20*nNodes) ==> row(nnz_tot)
      
      end subroutine ! ramg_calcPPEv


!************************************************************
!      ramg_jacobi
!      calculate u = D^-1 * [ -( L+U ) v + r ]
!
!      and also the residule in and out resout = | r-Au |
!                                       resin = | r-Av | 
!      
!************************************************************
      subroutine ramg_jacobi(Acolm,Arowp,Alhs,Anshg,Annz_tot,!A
     &                      r,u,v)
     
      include "common.h"
      integer,intent(in)                   :: Anshg,Annz_tot
      integer,intent(in),dimension(Anshg+1) :: Acolm
      integer,intent(in),dimension(Annz_tot) :: Arowp
      real(kind=8),intent(in),dimension(Annz_tot) :: Alhs
      real(kind=8),intent(in),dimension(Anshg) :: v,r
      real(kind=8),intent(inout),dimension(Anshg) :: u
      real(kind=8)                        :: tmp,tmpin,tmpout
      
      integer                              :: i,j,k,p

      real(kind=8)                     :: damp_jacobi

      ! Useless function, will be removed

      u = 0
      damp_jacobi = 1.0/ramg_relax
      do i=1,Anshg
         do k = Acolm(i)+1,Acolm(i+1)-1
            j = Arowp(k)
            tmp = Alhs(k)*v(j)
            u(i) = u(i) - tmp
         enddo
         u(i) = damp_jacobi*u(i) + r(i)
         u(i) = u(i)/Alhs(Acolm(i))
      enddo
       
      end subroutine ! ramg_jacobi

!***************************************************************
!     ramg_boundary_mls1: 
!     Do polynomial smoothing on boundary nodes only.
!     remember this can always be simplified to boundary jacobi
!**************************************************************      
      subroutine ramg_boundary_mls1(acolm,arowp,alhs,
     &                      anshg,annz_tot,!A
     &                      r,u,v,clevel,pflag,
     &                      ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in)                   :: anshg,annz_tot
      integer,intent(in),dimension(anshg+1) :: acolm
      integer,intent(in),dimension(annz_tot) :: arowp
      real(kind=8),intent(in),dimension(annz_tot) :: alhs
      real(kind=8),intent(in),dimension(anshg) :: v,r
      real(kind=8),intent(inout),dimension(anshg) :: u
      integer,intent(in)                      :: clevel
      integer,intent(in),dimension(nlwork)        :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      logical,dimension(anshg),intent(inout) :: pflag

      integer i,j,k,p,q
      real(kind=8),dimension(anshg)    :: tmp,v2,aux
      real(kind=8) :: cf1
      real(kind=8) :: cpusec(2)
      
      if (numpe.eq.1) return

      call cpu_time(cpusec(1))
      
      tmp = 0
      v2 = v

      ! tmp = A*x
      call ramg_calcAv_b(clevel,tmp,v2,ilwork,BC,iBC,iper,1,0)
!      call ramg_calcAv(acolm,arowp,alhs,anshg,annz_tot,
!     &                 tmp,v2,1)
      ! tmp = A*x-b=r
      tmp = tmp-r
      ! aux = A*r ! not necessary
      ! call ramg_calcAv_b(clevel,aux,tmp,ilwork,BC,iBC,iper,1,0)

      cf1 = 1.1*mlsCf(clevel,1)
      
      !v2 = mlsCf(clevel,1)*tmp!+mlsCf(clevel,2)*aux
      !v2 = 1.1*v2
      !v2 = tmp
      
      do i=1,anshg
      p = amg_paraext(clevel)%p(i)
      q = amg_paramap(clevel)%p(i)
      if (p.ne.(myrank+1)) then
          pflag(i) = .true.
      if ((p.gt.0).or.(p.ne.q)) then ! master boundary or extended
         u(i) = v(i)-cf1*tmp(i)
      endif
      endif
      enddo

      call cpu_time(cpusec(2))
      ramg_time(20) = ramg_time(20)+cpusec(2)-cpusec(1)

      end subroutine !ramg_boundary_mls1
 

!************************************************************
!      ramg_gauss
!      
!      forward and backward, defined in fwdbck
!      Gauss-Seidel smoothing, u = gaussseidel(M,r)
!
!************************************************************
      subroutine ramg_gauss(acolm,arowp,alhs,anshg,annz_tot,!A
     &                      r,u,v,fwdbck,clevel,
     &                      ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in)                   :: anshg,annz_tot
      integer,intent(in),dimension(anshg+1) :: acolm
      integer,intent(in),dimension(annz_tot) :: arowp
      real(kind=8),intent(in),dimension(annz_tot) :: alhs
      real(kind=8),intent(in),dimension(anshg) :: v,r
      real(kind=8),intent(inout),dimension(anshg) :: u
      integer,intent(in)                      :: fwdbck!1=fwd,2=bck
      integer,intent(in)                      :: clevel
      integer,intent(in),dimension(nlwork)        :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

     
      real(kind=8)                        :: tmp
      real(kind=8),dimension(anshg) :: u2
      integer                             :: istr,iend,iint
      integer                              :: i,j,k,p,ki,kj,kp
      logical                             :: postsmooth,p2
      logical,dimension(anshg) :: pflag
!      integer            :: itag,iacc,iother,numseg,isgbeg,itkbeg
!      integer            :: numtask

      u = v

      postsmooth = (fwdbck.eq.1)
      if (postsmooth) then ! post smoothing : 1->nshg
          istr = 1
          iend = anshg
          iint = 1
      else ! pre smoothing : nshg->1
          istr = anshg
          iend = 1
          iint = -1
      end if

      pflag = .false.
      ! boundary 
!      if (postsmooth) then
          call ramg_boundary_mls1(acolm,arowp,alhs,anshg,annz_tot,
     &         r,u,v,clevel,pflag,ilwork,BC,iBC,iper)
!      endif
!      if (.false.) then

      do i=istr,iend,iint
         ki = CF_map(clevel)%p(i)
         p2 = .not.(pflag(ki)) ! nodes that have not been treated
!         if (numpe.gt.1) then
!            p2=p2.and.(amg_paraext(clevel)%p(ki).eq.(myrank+1))
!         endif
         if (p2) then
         tmp = 0
         do k = acolm(ki)+1,acolm(ki+1)-1
            j = arowp(k)
            if ((postsmooth).or.(pflag(j))) then 
            ! save some if u(j)=0
            tmp = tmp + alhs(k)*u(j)
            endif
         enddo
         u(ki) = r(ki) - tmp
         pflag(ki) = .true.
         !u(ki) = u(ki)/alhs(acolm(ki))
         ! The above line is commented out because we know that
         ! A(i,i)=1.0, this is the result of prescaling of the matrix.
         endif !p2
      enddo
!      if (.not.postsmooth) then
      u2 = u
          call ramg_boundary_mls1(acolm,arowp,alhs,anshg,annz_tot,
     &         r,u,u2,clevel,pflag,ilwork,BC,iBC,iper)
!      endif

      end subroutine ! ramg_gauss

!**********************************************
!      ramg dense direct solver routines.
!      
!      ramg_sparse2dense
!      Sparse to dense form
!      
!      ramg_direct_LU ( setup only )
!      Outside package for direct solve sparse
!        matrix by LU
!**********************************************

      subroutine ramg_direct_LU(Acolm,Arowp,Alhs,Anshg,Annz_tot)
      use ramg_data

      integer,intent(in)                   :: Anshg,Annz_tot
      integer,intent(in),dimension(Anshg+1) :: Acolm
      integer,intent(in),dimension(Annz_tot) :: Arowp
      real(kind=8),intent(in),dimension(Annz_tot) :: Alhs

!      real(kind=8),dimension(Anshg,Anshg)      :: mtxA
!      integer,dimension(Anshg)             :: indx
      real(kind=8)                     :: d

      integer :: i,j

      !write(*,*) "ramg_setup_flag",ramg_setup_flag
      if (ramg_setup_flag .ne. 0) return

      !write(*,*)'again'
      if (allocated(cmtxA)) deallocate(cmtxA)
      if (allocated(cindx)) deallocate(cindx)

      allocate(cmtxA(Anshg,Anshg))
      allocate(cindx(Anshg))

      if (myrank.eq.master) then
      write(*,*)"Start to setup LU decomposition at coarsest level"
      endif
      call ramg_sparse2dense(Acolm,Arowp,Alhs,Anshg,Annz_tot,cmtxA)

!      Asol = Arhs
      
      call ludcmp(cmtxA,Anshg,Anshg,cindx,d)
!      call lubksb(mtxA,Anshg,Anshg,indx,Asol)
      if (myrank.eq.master) then
          write(*,*)"End of setup LU decomposition at coarsest level"
      endif
     
      end subroutine ! ramg_direct_LU

      subroutine ramg_sparse2dense(Acolm,Arowp,Alhs,Anshg,Annz_tot,
     &                             mtxA)

      integer,intent(in)                   :: Anshg,Annz_tot
      integer,intent(in),dimension(Anshg+1) :: Acolm
      integer,intent(in),dimension(Annz_tot) :: Arowp
      real(kind=8),intent(in),dimension(Annz_tot) :: Alhs
      real(kind=8),intent(inout),dimension(Anshg,Anshg) :: mtxA

      integer                       :: i,j,k,ip,jp,kp

      mtxA = 0
      do i=1,Anshg
         do j= Acolm(i),Acolm(i+1)-1
            k = Arowp(j)
            mtxA(i,k) = Alhs(j)
         enddo
      enddo

      end subroutine !ramg_sparse2dense
      
 
!**********************************************
!      ramg_direct_solve
!**********************************************
      subroutine ramg_direct_solve(
     &             Arhs,Asol,
     &             colm,rowp,lhsK,lhsP,
     &             ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      real(kind=8),intent(in),dimension(amg_nshg(ramg_levelx))
     &                                                     ::Arhs
      real(kind=8),intent(inout),dimension(amg_nshg(ramg_levelx)) 
     &                                                     :: Asol
      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)  :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)  :: lhsP

      integer,intent(in),dimension(nlwork)        :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      real(kind=8)              :: mres_i,mres_n,mres_o
      real(kind=8),dimension(amg_nshg(ramg_levelx)) :: myV
      integer                      :: titer,cl

!      Asol=Arhs
!      return
      cl = ramg_levelx

      Asol = 0
      call ramg_L2_norm_p(cl,Arhs,1,mres_i)
      !mres_i=sqrt(flesDot1(Arhs,amg_nshg(cl),1))
      mres_o = 1e-7*mres_i
      titer = 0

      !write(*,*)'mcheck direct solve begin,',myrank
      if (iamg_c_solver.eq.0) then
          ! Do two smoothing steps
          call ramg_smoother(cl,myV,Asol,Arhs,
     &              colm,rowp,lhsK,lhsP,
     &              ilwork,BC,iBC,iper,2)
          call ramg_smoother(cl,Asol,myV,Arhs,
     &              colm,rowp,lhsK,lhsP,
     &              ilwork,BC,iBC,iper,1) 
      else if ( iamg_c_solver .eq. 1) then ! direct solve using G-S
          IF (.FALSE.) THEN ! smoother solve to exact
          mres_n = mres_i
          do while ( (mres_n.gt.mres_o).and.(titer.lt.1000) )
          call ramg_smoother(cl,myV,Asol,Arhs,
     &              colm,rowp,lhsK,lhsP,
     &              ilwork,BC,iBC,iper,2)
          call ramg_smoother(cl,Asol,myV,Arhs,
     &              colm,rowp,lhsK,lhsP,
     &              ilwork,BC,iBC,iper,1) 
          call ramg_calcAv_g(cl,myV,Asol,colm,rowp,lhsK,lhsP,
     &         ilwork,BC,iBC,iper,1)
          myV = myV-Arhs
          mres_n=sqrt(flesDot1(myV,amg_nshg(cl),1))
          titer = titer + 1
          enddo
          ELSE ! cg solve to exact
              call ramg_CG(Asol,Arhs,cl,
     &        colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper,titer)
          ENDIF
      else ! direct solve using dense-matrix LU-decomposition
          Asol = Arhs
          call lubksb(cmtxA,amg_nshg(cl),amg_nshg(cl),cindx,Asol)
      endif

      end subroutine ! ramg_direct_solve


      subroutine ramg_smoother(level,xnew,xold,xrhs,
     &                         colm,rowp,lhsK,lhsP,
     &                         ilwork,BC,iBC,iper,fwdbck)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      
      integer level
      real(kind=8),dimension(amg_nshg(level)),intent(in) ::xrhs,xold
      real(kind=8),dimension(amg_nshg(level)),intent(inout) :: xnew
      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp

      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP

      integer,intent(in),dimension(nlwork)   :: ilwork
      integer,intent(in),dimension(nshg)     :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      integer fwdbck

      integer i,k
      logical :: jacobi

      ! NOTICE: fwdbck=2 : pre-smoothing
      ! NOTICE: fwdbck=1 : post-smoothing
      
      jacobi = .false. !(ramg_relax.gt.0)
      ! Disable Jacobi, if needed, using 1-st order MLS instead

      k = 1
      do i=1,k
          if ( iamg_smoother .eq. 3 ) then !.or. (level.gt.1) ) then
              call ramg_mls_apply(xnew,xold,xrhs,level,
     &                  colm,rowp,lhsK,lhsP,
     &                  ilwork,BC,iBC,iper,fwdbck)
          else if ( iamg_smoother .eq. 1 ) then 
          call ramg_gauss(amg_A_colm(level)%p,amg_A_rowp(level)%p,
     &                    amg_A_lhs(level)%p,amg_nshg(level),
     &                    amg_nnz(level),xrhs,xnew,xold,
     &                    fwdbck,level,
     &                    ilwork,BC,iBC,iper)
          else if ( ( iamg_smoother .eq. 2) ) then !.and. (level.eq.1) ) then
              call ramg_cheby_apply(xnew,xold,xrhs,level,
     &             colm,rowp,lhsK,lhsP,
     &             ilwork,BC,iBC,iper,fwdbck)
          else if (jacobi) then
             call ramg_jacobi(amg_A_colm(level)%p,
     &                 amg_A_rowp(level)%p,amg_A_lhs(level)%p,
     &                 amg_nshg(level),amg_nnz(level),
     &                 xrhs,xnew,xold)
          endif
      enddo

      end subroutine !ramg_smoother
