!*********************************************************
!      ramg control 
!      Control AMG preparation/solve process
!      will have dynamic control (based on solver vecters)
!*********************************************************

      subroutine ramg_control(colm,rowp,lhsK,lhsP,
     &                         ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"

      integer,intent(in),dimension(nshg+1)             :: colm
      integer,intent(in),dimension(nnz_tot)            :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)     :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)     :: lhsP
      ! the boundary info
      integer, intent(in), dimension(nlwork)           :: ilwork
      integer, intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)   :: BC

      integer i

      if (iamg_init .eq. 0 ) then
          ! The overall initialization
          ! do init
          ramg_winct = 0
          ramg_setup_flag = 0
          iamg_init = 1
          ramg_time = 0
      else
         ramg_setup_flag = mod(ramg_setup_flag+1 ,iamg_setup_frez)
      end if 
 
      ! Extract PPE
          call ramg_extract_ppe(colm,rowp,lhsK,lhsP,
     &               ilwork,BC,iBC,iper)
      ! Coarsening
          call ramg_prep(ilwork,BC,iBC,iper)
          call ramg_init_ilwork(ilwork,BC,iBC,iper)
      ! Prepare for MLS smoothing
      if (iamg_smoother.eq.2) then
          call ramg_cheby_setup(colm,rowp,lhsK,lhsP,
     &                        ilwork,BC,iBC,iper)
      else
          call ramg_mls_setup(colm,rowp,lhsK,lhsP,
     &                        ilwork,BC,iBC,iper)
      endif
      if (iamg_c_solver.eq.2) then
      ! Setup Coarsest Direct solver
          call ramg_direct_LU(amg_A_colm(ramg_levelx)%p,
     &          amg_A_rowp(ramg_levelx)%p,
     &          amg_A_lhs(ramg_levelx)%p,amg_nshg(ramg_levelx),
     &          amg_nnz(ramg_levelx))
      endif
      if (maxnev.gt.0) then
      ! Setup GGB
         call ramg_ggb_setup(colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
      endif

     
      end subroutine !ramg_control

!*********************************************************
!     Ramg preparation
!*********************************************************

      subroutine ramg_prep(ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in),dimension(nlwork)      :: ilwork
      integer,intent(in),dimension(nshg)        :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      logical                  :: maxstopsign,mxs2
      integer                  :: i,p,p2

      if (ramg_setup_flag.ne.0) return
      maxstopsign = .false.
      i = 1
      do while ((i+1.le.iamg_nlevel) .and. (.not.maxstopsign))
         if (amg_nshg(i+1).ne.0) then
             call ramg_deallocate(i+1)
         end if
         call ramg_coarse_setup(i,i+1,strong_eps,iamg_interp,
     &        ramg_trunc,
     &        ilwork,BC,iBC,iper,nshg,nlwork,ndofBC)
         call ramg_calcITAI(i,i+1,maxstopsign)
         if ((iamg_verb.gt.2).and.(myrank.eq.master)) then
         write(*,*)'COARSEN: level:',i+1,' nshg:',amg_nshg(i+1),
     &          ' nnz:',amg_nnz(i+1)
         endif
         i = i+1
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         maxstopsign = (amg_nshg(i).eq.amg_nshg(i-1)) 
         IF (.true.) THEN
         !call ramg_checkcoarse(i,ilwork,BC,iBC,iper,maxstopsign)
         p = 1
         if (maxstopsign) p = 0
         call MPI_AllReduce(p,p2,1,MPI_INTEGER,MPI_SUM,
     &                MPI_COMM_WORLD,ierr)
         if (p2.eq.0) then
           maxstopsign = .true.
           !write(*,*)'mcheck stopped'
         else
           maxstopsign = .false.
         endif
         ENDIF
      enddo
      ramg_levelx=i
      if (maxstopsign) then 
         ramg_levelx = ramg_levelx-1
      endif

!      if ( (irun_amg_prec.eq.1).and.(mlsDeg.gt.0) ) then
      if (.false.) then
      deallocate(amg_A_colm(1)%p)
      deallocate(amg_A_rowp(1)%p)
      deallocate(amg_A_lhs(1)%p)
      endif

      allocate(CF_map(ramg_levelx)%p(amg_nshg(ramg_levelx)))
      allocate(CF_revmap(ramg_levelx)%p(amg_nshg(ramg_levelx)))
      do i=1,amg_nshg(ramg_levelx)
         CF_map(ramg_levelx)%p(i) = i
         CF_revmap(ramg_levelx)%p(i) = i
      enddo

      !call ramg_output_coarsening

      end subroutine !ramg_prep

      subroutine ramg_checkcoarse(level1,ilwork,BC,iBC,iper,
     &           cfstop)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in) :: level1
      integer,intent(in),dimension(nlwork)      :: ilwork
      integer,intent(in),dimension(nshg)        :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      logical,intent(inout)                  :: cfstop
      integer                  :: i,numtask,itkbeg,numseg,iacc
      integer                  :: j,k,p
      integer,allocatable,dimension(:) 
     &            :: subcfstop,subnnz,subcf,subcfrev,subnei
      real(kind=8) :: rhoratio
      logical      :: subneireduce

      IF (( iamg_reduce.le.1).or.(numpe.gt.1)) THEN
!      IF (.TRUE.) THEN
      if (numpe.ge.2) then 
          numtask = ilwork(1)+1
      else 
          numtask = 1
      endif
      
      allocate(subcfrev(numpe))
      subcfrev = 0
      allocate(subnei(numtask))
      subnei = 0
 
      subnei(1) = myrank+1
      subcfrev(subnei(1))=1
      itkbeg = 1
      do i=2,numtask
         subnei(i)=ilwork(itkbeg+3)+1
         subcfrev(subnei(i)) = i
         itkbeg = itkbeg + 4 + 2*ilwork(itkbeg+4)
      enddo
   
      ELSE !reduced
          numtask = rmapmax-1
          allocate(subcfrev(numtask))
          allocate(subnei(numtask))
          do i=1,numtask
             subcfrev(i) = i
          enddo
      ENDIF

      allocate(subcfstop(numtask))
      allocate(subcf(numtask))
      allocate(subnnz(numtask))
      subcfstop = 0
      subcf = 0
      subnnz = 0


      do i = 1,amg_nnz(level1)
         k = amg_A_rowp(level1)%p(i)
         p = iabs(amg_paramap(level1)%p(k))
!         if (iamg_reduce.gt.1) then
!            p = 1 ! for reduced only
!         endif
         p = subcfrev(p)
         subnnz(p) = subnnz(p) + 1
      enddo
      do i = 1,amg_nshg(level1)
         p = iabs(amg_paramap(level1)%p(i))
!         if (iamg_reduce.gt.1) then
!            p = 1 ! for reduced only
!         endif
         p = subcfrev(p)
         subcf(p) = subcf(p) + 1
      enddo

      do i=1,numtask
         IF ((iamg_reduce.le.1).or.(numpe.gt.1)) THEN
!         IF (.TRUE.) THEN
         p = subnei(i)
         subneireduce = (p.ne.(myrank+1))
         ELSE !reduced
             subneireduce = (i.gt.iamg_reduce)
         END IF
         if ((subcf(i).lt.30).and.subneireduce) then
             subcfstop(i) = 1
         endif
         if (.not.subneireduce) then
             if (subcf(i).lt.200) then
                 subcfstop(i) = 1
             endif
             if ((subnnz(i)/(subcf(i)**2)).gt.0.6) then
                 subcfstop(i) = 1
             endif
         end if
      enddo
  
      cfstop = .true.
      do i=1,numtask
         if (subcfstop(i).eq.0) then
             cfstop = .false.
             exit
         endif
      enddo
    

      IF ((iamg_reduce.le.1).or.(numpe.gt.1)) THEN
      ! if interior is coarsened, boundary should be fixed whatever
      if (subcfstop(1).eq.1) then
          cfstop = .true.
      endif
      ENDIF

      if (.not.cfstop) then
      do i=1,amg_nshg(level1)
         p = iabs(amg_paramap(level1)%p(i))
!         if (iamg_reduce.gt.1) then
!             p = 1 ! reduced only
!         endif
         p = subcfrev(p)
         if (subcfstop(p).eq.1) then
         amg_paramap(level1)%p(i) = -iabs(amg_paramap(level1)%p(i))
         endif
      enddo
      endif

      deallocate(subcfrev)
      deallocate(subnei)
      deallocate(subcfstop)
      deallocate(subcf)
      deallocate(subnnz)
     
      end subroutine ! ramg_checkcoarse

!******************************************************************
!     A bunch of code that hook to leslib's cg solve
!     Force it to restart with AMG if meets plateau
!******************************************************************

!********
!      sc writes:
!      irun_amg_prec is a variable controling how to use AMG wisely.
!      0 : no run
!      1 : always use AMG prec'd CG
!      2 : restart CG a ( if plain CG hits plateau, do again with AMG )
!      3 : restart CG b ( if plain exceeds maxiter, do again with AMG )
!      also refer to input.config for a detailed description.
!********      

      
      subroutine ramg_normcheck(tmpnorm)
      use ramg_data

      include "common.h"
      include "mpif.h"

      real(kind=8),intent(inout) ::  tmpnorm
      real(kind=8) :: sqnorm

      !write(*,*)myrank,' winct:',ramg_winct,sqrt(tmpnorm)
      !return
      if (irun_amg_prec.ne.2) return ! No control at all
      if (ramg_winct.eq.0) then
          !ramg_window = sqrt(tmpnorm)
          return ! Not relevant
      endif
      ramg_winct=0
      if (ramg_flag.gt.1) return ! Second run does not need control
     
      ! Following are for the first run
      !write(*,*)'normcheck: ',sqnorm
      !return;

      if (ramg_redo.ge.100) then
          ! The rest of GMRES should be terminated
          if (ramg_redo.lt.104) then
             ramg_winct = 1
             ramg_redo=ramg_redo+1
!             write(*,*)'***'
          else
             ramg_redo = 0
             tmpnorm = 0
          endif
          return
      endif
      sqnorm = sqrt(tmpnorm)

      if (ramg_redo.le.7) then
          ramg_redo = ramg_redo + 1
          call ramg_winpushin(sqnorm)
          return
      endif
 
      call ramg_winpushin(sqnorm)

      if (ramg_window(7).gt.ramg_window(1)) then ! Stop point
!          write(*,*)'myrank:',myrank,ramg_window(7),ramg_window(1)
!      if (ramg_window(7).lt.1e-3) then
          ! stop point in cg
          if (myrank.eq.master) then
          write(*,*)"Prepare to restart CG"
          endif
          tmpnorm = 0
          ramg_redo = 100
          ramg_winct = 1
          ramg_flag = ramg_flag -1
      endif

      end subroutine ! ramg_normcheck

      subroutine ramg_winpushin(sqnorm)
      use ramg_data
      include "common.h"
      include "mpif.h"

      real(kind=8),intent(in) :: sqnorm
      integer i
      integer totlen

      totlen = 7
      
      do i=1,totlen-1
        ramg_window(i)=ramg_window(i+1)
      enddo
      ramg_window(totlen) = sqnorm
      
      end subroutine ! ramg_winpushin


!*********************************************************
!      <EOF> RAMG Control
!*********************************************************
