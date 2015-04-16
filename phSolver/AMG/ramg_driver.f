!!*****************************************************
!
!         Turbo AMG (ramg) interface functions
!         1. V_cycle
!         2. plain
!         3. driver      
!
!*****************************************************

      !*******************************************
      !    ramg_V_cycle:
      !     One V_cycle call. 
      !     Now use SAMG, later will use own code
      !*******************************************
      recursive subroutine ramg_V_cycle(
     &               ramg_sol,ramg_rhs,clevel,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper
     &           )!ramg_V_cycle
      use ramg_data
      include "common.h"
      include "mpif.h"

      !Variable Declaration
      ! arguments
      integer, intent(in)         :: clevel
      real(kind=8),intent(inout),dimension(amg_nshg(clevel))
     &               :: ramg_sol,ramg_rhs

      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)  :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)  :: lhsP
      integer,intent(in),dimension(nlwork)          :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC


      ! local
      real(kind=8),dimension(amg_nshg(clevel))  :: myvF,myvE
      real(kind=8),dimension(:),allocatable :: myvC,myvCS
      real(kind=8)          :: cpusec(10)

      ! In scale
      
      if (clevel.ne.1) then
      ramg_rhs = ramg_rhs*amg_ppeDiag(clevel)%p
      endif

      if (ramg_levelx-clevel .eq. 0) then 
      ! finest level solver
          call ramg_direct_solve(ramg_rhs,ramg_sol,
     &    colm,rowp,lhsK,lhsP,
     &    ilwork,BC,iBC,iper)
      !**********************
      ! higher level smoother
      !**********************
      else if (ramg_levelx-clevel .gt. 0) then
          allocate(myvC(amg_nshg(clevel+1)))
          allocate(myvCS(amg_nshg(clevel+1)))
          
          ! smoothing (pre)
          myvF = 0 ! initial guess is zero
          call ramg_smoother(clevel,ramg_sol,myvF,ramg_rhs,
     &              colm,rowp,lhsK,lhsP,
     &              ilwork,BC,iBC,iper,2) !dx1
          
          ! restriction
          call ramg_calcAv_g(clevel,myvF,ramg_sol,colm,rowp,lhsK,lhsP,
     &               ilwork,BC,iBC,iper,1)
          myvF =  myvF - ramg_rhs! r2
          call ramg_calcIvFC(clevel,clevel+1,myvF,myvC)
          ! up one level, solve
          myvCS = 0
          call ramg_V_cycle(myvCS,myvC,clevel+1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
          ! prolongation
          call ramg_calcIvCF(clevel,clevel+1,myvCS,myvF,
     &         ilwork,BC,iBC,iper) ! v2hat
          ramg_sol = ramg_sol - myvF
          
          ! smoothing (post)
          call ramg_smoother(clevel,myvF,ramg_sol,ramg_rhs,
     &              colm,rowp,lhsK,lhsP,
     &              ilwork,BC,iBC,iper,1)
          ramg_sol =  myvF

          deallocate(myVC)
          deallocate(myVCS)
      end if
      
      ! Out scale
      if (clevel.ne.1) then
      ramg_sol = ramg_sol * amg_ppeDiag(clevel)%p
      endif

      return

      end subroutine ramg_V_cycle
      

!*****************************************************
!      ramg Cg solver
!      using ramg_V_cycle as preconditioner
!*****************************************************
      subroutine ramg_CG(
     &              sol,rhs,level,
     &              colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper,iterNum)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in) :: level
      real(kind=8),intent(in),dimension(amg_nshg(level)) :: rhs
      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: sol

      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp

      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP

      integer,intent(in),dimension(nlwork)   :: ilwork
      integer,intent(in),dimension(nshg)     :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      integer,intent(inout)                    :: iterNum

      ! local
      real(kind=8),dimension(amg_nshg(level))   :: cgP,cgQ,cgZ,cgR
      real(kind=8)                       :: rz,rz_0,pq,alpha,beta
      real(kind=8)               :: tmp,norm_0,norm_1,norm_e,norm_c
      real(kind=8)               :: mres_0,mres_n

      cgR = rhs
      call ramg_L2_norm_p(level,cgR,1,norm_0)
      norm_c = norm_0
      norm_e = norm_0*1e-7

      cgZ = cgR

      call ramg_dot_p(level,cgZ,cgR,1,rz)
      rz_0 = rz
      cgP = cgZ
      sol = 0

      iterNum = 1

      do while ( (norm_c .gt. norm_e).and.(iterNum.lt.500) )
         call ramg_calcAv_g(level,cgQ,cgP,colm,rowp,lhsK,lhsP,
     &                      ilwork,BC,iBC,iper,1)
         call ramg_dot_p(level,cgP,cgQ,1,pq)
         alpha = rz/pq
         sol = sol + alpha*cgP
         cgR = cgR - alpha*cgQ
         call ramg_L2_norm_p(level,cgR,1,tmp)
         norm_c = tmp
         cgZ = cgR
         tmp = rz
         call ramg_dot_p(level,cgZ,cgR,1,rz)
         beta = rz/tmp
         cgP = beta*cgP+cgZ
         iterNum = iterNum + 1
      enddo
     
      end subroutine ! ramg_CG
      
      subroutine ramg_PCG(
     &              sol,rhs,run_mode,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper
     &           )!ramg_CG
      use ramg_data
      include "common.h"
      include "mpif.h"
      
      real(kind=8),intent(inout),dimension(nshg) :: sol
      real(kind=8),intent(in),dimension(nshg)    :: rhs
      integer,intent(in)                                :: run_mode

      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)  :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)  :: lhsP
      integer,intent(in),dimension(nlwork)          :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC


      ! local
      real(kind=8),dimension(amg_nshg(1))   :: cgP,cgQ,cgZ,cgR
      real(kind=8)                       :: rz,rz_0,pq,alpha,beta
      real(kind=8)               :: tmp,norm_0,norm_p,norm_e,norm_c
      real(kind=8)               :: mres_0,mres_n
      real(kind=8),dimension(nshg,4)               :: diag
      integer                    :: iterNum

      ! controls

      real(kind=8),dimension(maxiters) :: rconv
      real(kind=8) :: avgconv
      logical :: vcycle,gcycle,restart
      character cflagv,cflagg

      cgR = rhs 
      call ramg_L2_norm_p(1,cgR,1,norm_0)
      norm_c = norm_0
      norm_e = norm_0 * epstol(2)
      write(*,*)'norm_0 ', norm_0

      vcycle = .true.
      gcycle = .true.
      restart = .false.

      ! amg preconditioning control         

      if (vcycle) then
         call ramg_V_cycle(cgZ,cgR,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         if (gcycle) then
         call ramg_G_cycle(cgZ,cgR,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         endif
      else
          cgZ = cgR
      endif

      call ramg_dot_p(1,cgZ,cgR,1,rz)
      rz_0 = rz

      cgP = cgZ
      sol = 0

      iterNum = 1
      rconv(1) = 1.0
      norm_p=norm_c
      avgconv = 0

      do while ( (norm_c.gt.norm_e).and.(iterNum.lt.maxiters))
      !do while ((iterNum.lt.maxiters))
         call ramg_calcAv_g(1,cgQ,cgP,colm,rowp,lhsK,lhsP,
     &                      ilwork,BC,iBC,iper,1)
         call ramg_dot_p(1,cgP,cgQ,1,pq)
         alpha = rz/pq
         sol = sol + alpha*cgP
         cgR = cgR - alpha*cgQ
         call ramg_L2_norm_p(1,cgR,1,tmp)
         norm_c = tmp
         if (iterNum.eq.1) rconv(1) = norm_p/norm_c

      ! amg preconditioning control         
      if (vcycle) then
         call ramg_V_cycle(cgZ,cgR,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         if (gcycle) then
         call ramg_G_cycle(cgZ,cgR,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         endif
      else
          cgZ = cgR
      endif
      
         tmp = rz
         call ramg_dot_p(1,cgZ,cgR,1,rz)
         beta = rz/tmp
         cgP = beta*cgP + cgZ
         iterNum = iterNum + 1
         rconv(iterNum) = norm_p/norm_c
         norm_p = norm_c
         cflagv='n'
         cflagg='n'
         if (vcycle) cflagv='v'
         if (gcycle) cflagg='g'
         if (myrank.eq.master) then 
         write(*,"((A7)(I4)(E14.4)(T30)(F8.4)(T40)(A1)(A1))")
     &    'AMGCG: ',
     &    iterNum,norm_c/norm_0,rconv(iterNum),cflagv,cflagg
         end if
      enddo

      end subroutine !ramg_PCG

!*******************************************
!      ramg Interface
!      Interface for libLES
!******************************************* 
      subroutine ramg_interface(
     &           colm,rowp,lhsK,lhsP,flowDiag,
     &           mcgR,mcgZ,
     &           ilwork,BC,iBC,iper
     &           )
      use ramg_data
      include "common.h"
      include "mpif.h"
      integer,intent(in),dimension(nshg+1)             :: colm
      integer,intent(in),dimension(nnz_tot)            :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)     :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)     :: lhsP
      real(kind=8),intent(in),dimension(nshg)          :: mcgR
      real(kind=8),intent(inout),dimension(nshg)       :: mcgZ
      integer, intent(in), dimension(nlwork)           :: ilwork
      integer, intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)   :: BC
      real(kind=8),dimension(nshg)                     :: myR
      real(kind=8),intent(in),dimension(nshg,4)        :: flowDiag
      
      integer::i,j,k,amgmode
      logical :: precflag
      !character*10                                     :: fname

      mcgZ = 0
      amgmode = 11

      ramg_winct = mod(ramg_winct+1,2)

      if (irun_amg_prec.eq.0) then
          precflag = .false.
      else if ((irun_amg_prec.ge.2).and.(ramg_flag.eq.1)) then
      ! first solve attempt with CG if using smart solver
          precflag=.false.
      else 
          precflag = .true.
      endif
          
      if (precflag) then !.and.(ramg_redo.ne.0)) then
          myR = mcgR*amg_ppeDiag(1)%p
          call ramg_driver(colm,rowp,lhsK,lhsP,
     &           nshg,nnz_tot,nflow,
     &           myR,
     &           nlwork,ilwork,ndofBC,BC,iBC,iper,mcgZ,1,
     &           ramg_eps,amgmode)
          mcgZ = mcgZ*amg_ppeDiag(1)%p
      else
          mcgZ = mcgR
      endif
      
      end subroutine ramg_interface

!*******************************************      
!        ramg Driver
!        1. extract PPE / system
!        2. call coarsening
!        3. solve     
!********************************************

      !*************************************
      !    ramg_driver
      !     Input:  global matrix,# of systems
      !             control params
      !     Output: solution
      !*************************************
      subroutine ramg_driver(
     &           colm,rowp,lhsK,lhsP,
     &           nshg,nnz_tot,nflow,
     &           pperhs,
     &           nlwork,ilwork,ndofBC,BC,iBC,iper,
     &           ramg_sol,n_sol,
     &           amg_eps,amg_mode
     &           )
      
      use ramg_data
      implicit none
      
      !***********parameters**************
      !the matrix
      integer,intent(in)                               :: nshg
      integer,intent(in)                        :: nnz_tot
      integer,intent(in)                               :: nflow
      !the matrix
      integer,intent(in),dimension(nshg+1)             :: colm
      integer,intent(in),dimension(nnz_tot)            :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)     :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)     :: lhsP
      !the forcing term, rhs
      real(kind=8),intent(in),dimension(nshg)       :: pperhs
      ! the boundary info
      integer, intent(in)                              :: nlwork
      integer, intent(in), dimension(nlwork)           :: ilwork
      integer, intent(in)                              :: ndofBC
      integer, intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)   :: BC
      !the output solution
      integer, intent(in)                              :: n_sol
      real(kind=8),intent(inout),dimension(nshg,n_sol) :: ramg_sol

      !the tolerance
      real(kind=8),intent(in)                          :: amg_eps
      !control run mode
      integer,intent(inout)                            :: amg_mode
      ! my AMG parameter
      !*********parameters end**************

      !*****local variables*****************
      real(kind=8)                      :: cpusec(10)

      call cpu_time(cpusec(1))

      if (amg_mode .eq. 11 ) then     ! solve PPE les Precondition ramg
         call cpu_time(cpusec(2))
         call ramg_V_cycle(ramg_sol,pperhs,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         call ramg_G_cycle(ramg_sol,pperhs,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         call cpu_time(cpusec(3))
         ramg_time(4)=ramg_time(4)+cpusec(3)-cpusec(2)
      else if (amg_mode .eq. 3) then ! solve PPE CG
          call cpu_time(cpusec(2))
          call ramg_PCG(ramg_sol,pperhs,1,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
          call cpu_time(cpusec(3))
          ramg_time(4)=ramg_time(4)+cpusec(3)-cpusec(2)
      end if

      return

      end subroutine ramg_driver 
!*******************************************
!      <EOF> ramg Driver
!*******************************************      
