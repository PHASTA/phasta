!**********************************************************
!    RAMG tools, functions needed in ramg_driver and
!    ramg_interface. 
!    Module defined in ramg_data.f      
!
!    ramg_calcIAI
!    ramg_calcIvFC/CF
!    ramg_direct
!    ramg_allocate
!    ramg_deallocate
!    ramg_dump
!      
!***********************************************************
      
      
!!***********************************************************
!      ramg_calcIvCF
!       calculate V coarse to fine, 
!       v = I * VC
!***********************************************************
      subroutine ramg_calcIvCF(level1,level2,VC,vf,
     &           ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      
      !implicit none
      integer,intent(in)                :: level1,level2
      real(kind=8),intent(inout),dimension(amg_nshg(level1)) :: vf
      real(kind=8),intent(in),dimension(amg_nshg(level2)) :: VC

      integer,intent(in),dimension(nlwork)          :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      integer                           :: i,j,k,p

      real(kind=8) :: cpusec(2)

      call cpu_time(cpusec(1))

      vf = 0
      do i=1,amg_nshg(level1)
         do k=I_cf_colm(level1)%p(i),I_cf_colm(level1)%p(i+1)-1
            j = I_cf_rowp(level1)%p(k)
            vf(i) = vf(i) + VC(j)*I_cf(level1)%p(k)
         enddo
      enddo

      if ( level1 .eq. -1 ) then
          call commIn(vf,ilwork,1,iper,iBC,BC)
          call commOut(vf,ilwork,1,iper,iBC,BC)
      end if

      call cpu_time(cpusec(2))
      ramg_time(13)=ramg_time(13)+cpusec(2)-cpusec(1)

      end subroutine ! ramg_calcIvCF
      
!***********************************************************
!      ramg_calcIvFC
!       calculate v fine to coarse, 
!       VC = IT * v
!***********************************************************
      subroutine ramg_calcIvFC(level1,level2,vf,VC)
      use ramg_data
      implicit none
      integer,intent(in)                :: level1,level2
      real(kind=8),intent(in),dimension(amg_nshg(level1)) :: vf
      real(kind=8),intent(inout),dimension(amg_nshg(level2)) :: VC

      integer                           :: i,j,k,p
      real(kind=8) :: cpusec(2)

      call cpu_time(cpusec(1))


      VC = 0
      do i=1,amg_nshg(level2)
         do k=I_fc_colm(level1)%p(i),I_fc_colm(level1)%p(i+1)-1
            j = I_fc_rowp(level1)%p(k)
            VC(i) = VC(i) + vf(j)*I_fc(level1)%p(k)
         enddo
      enddo
      call cpu_time(cpusec(2))
      ramg_time(13)=ramg_time(13)+cpusec(2)-cpusec(1)

      end subroutine ! ramg_calcIvFC
     
!!***********************************************************
!      ramg_calcSvCF
!       calculate V coarse to fine, 
!       v = S * VC
!***********************************************************
      
!***********************************************************
!      ramg_calcIvFC
!       calculate v fine to coarse, 
!       VC = ST * v
!***********************************************************

!************************************************************
!      ramg_calcAv
!      calculate u = A*v 
!************************************************************
      subroutine ramg_calcAv(gcolm,growp,glhs,gnshg,gnnz_tot,
     &                         u,v,gcomm)
      use ramg_data
      integer,intent(in)                 :: gnshg,gnnz_tot
      integer,intent(in),dimension(gnshg+1) :: gcolm
      integer,intent(in),dimension(gnnz_tot) :: growp
      real(kind=8),intent(in),dimension(gnnz_tot) :: glhs
      real(kind=8),intent(in),dimension(gnshg) :: v
      real(kind=8),intent(inout),dimension(gnshg) :: u
      integer,intent(in) :: gcomm

      integer                             :: i,j,k,p

      u = 0
      do i=1,gnshg
         do k=gcolm(i),gcolm(i+1)-1
            j = growp(k)
            u(i) = u(i)+glhs(k)*v(j)
         enddo
      enddo

!      if (gcomm.eq.1) then
!      end if
      
      end subroutine ! ramg_calcAv

!************************************************************
!      ramg_calcAv_g
!      calculate u = A*v 
!      Globally: commout, do product, commin, (zeroout)
!************************************************************
      subroutine ramg_calcAv_g(level,u,v,colm,rowp,lhsK,lhsP,
     &           ilwork,BC,iBC,iper,gcomm)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in),dimension(nlwork) :: ilwork
      integer,intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)  :: BC
      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp

      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP
      integer,intent(in) :: level
      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: v
      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: u
      integer,intent(in) :: gcomm

      integer                             :: i,j,k,p,mi,mj
      real(kind=8)  :: cpusec(2)

      call cpu_time(cpusec(1))
 
      IF (level.eq.1) THEN
      !IF (.FALSE.) THEN
      call ramg_PPEAp(u,v,colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)    
      ELSE
      if (gcomm.eq.1) then
      call ramg_commOut(v,level,ilwork,1,iper,iBC,BC)
      endif
      u = 0
      do i=1,amg_nshg(level)
         mi = amg_paramap(level)%p(i)
         do k=amg_A_colm(level)%p(i),amg_A_colm(level)%p(i+1)-1
            j = amg_A_rowp(level)%p(k)
            mj = amg_paramap(level)%p(j)
            if (.not.( (mi.eq.mj).and.(mi.lt.0) )) then
                u(i) = u(i)+amg_A_lhs(level)%p(k,1)*v(j)
            endif
         enddo
      enddo
      if (gcomm.eq.1) then
      call ramg_commIn(u,level,ilwork,1,iper,iBC,BC)
      endif
      ENDIF
      call cpu_time(cpusec(2))
      if (level.eq.1) then
      ramg_time(11) = ramg_time(11) + cpusec(2)-cpusec(1)
      else
      ramg_time(12) = ramg_time(12) + cpusec(2)-cpusec(1)
      endif
     
      end subroutine ! ramg_calcAv_g

!************************************************************
!      ramg_calcAv_b: same as calcAv_g, but only apply
!                     on boundary nodes.
!      calculate u = A*v 
!      Globally: commout, do product, commin, (zeroout)
!************************************************************
      subroutine ramg_calcAv_b(level,u,v,
     &           ilwork,BC,iBC,iper,gcomm,diag)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in),dimension(nlwork) :: ilwork
      integer,intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)  :: BC

      integer,intent(in) :: level
      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: v
      real(kind=8),intent(inout),dimension(amg_nshg(level)) :: u
      integer,intent(in) :: gcomm
      integer,intent(in) :: diag !=1: NOT include diagonal A_(i,i)
                                 !=0: include diagonal 

      integer                             :: i,j,k,p,mi,mj,mk
      real(kind=8)  :: cpusec(2)

      call cpu_time(cpusec(1))
 
      if (gcomm.eq.1) then
      call ramg_commOut(v,level,ilwork,1,iper,iBC,BC)
      endif
      u = 0
      do i=1,amg_nshg(level)
         mi = amg_paramap(level)%p(i)
         mk = amg_paraext(level)%p(i)
         IF (mk.ne.(myrank+1)) THEN ! only treat boundary nodes
         do k=amg_A_colm(level)%p(i)+diag,amg_A_colm(level)%p(i+1)-1
            j = amg_A_rowp(level)%p(k)
            mj = amg_paramap(level)%p(j)
            if (.not.( (mi.eq.mj).and.(mi.lt.0) )) then
                u(i) = u(i)+amg_A_lhs(level)%p(k,1)*v(j)
            endif
         enddo
         ELSE
             u(i) = v(i)
         ENDIF
      enddo
      if (gcomm.eq.1) then
      call ramg_commIn(u,level,ilwork,1,iper,iBC,BC)
      endif
     
      end subroutine ! ramg_calcAv_b


!*************************************************
!      check_paracomm
!*************************************************
      subroutine ramg_check_paracomm(ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in),dimension(nlwork)        :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      integer i,j
      type(r1d),dimension(ramg_levelx) :: tarray
      character :: fname*80

      do i=1,ramg_levelx
         allocate(tarray(i)%p(amg_nshg(i)))
         do j=1,amg_nshg(i)
            call random_number(tarray(i)%p)
         enddo
      enddo

      do i=1,ramg_levelx
         write(fname,'((A8)(I1))')'tarray_a',i
         call ramg_dump(tarray(i)%p,amg_nshg(i),fname)
      enddo

      do i=1,ramg_levelx
         call ramg_commIn(tarray(i)%p,i,ilwork,1,iper,iBC,BC)
         write(fname,'((A8)(I1))')'tarray_i',i
         call ramg_dump(tarray(i)%p,amg_nshg(i),fname)
      enddo

      do i=1,ramg_levelx
         call ramg_commOut(tarray(i)%p,i,ilwork,1,iper,iBC,BC)
         write(fname,'((A8)(I1))')'tarray_b',i
         call ramg_dump(tarray(i)%p,amg_nshg(i),fname)
      enddo

!      call commOut(tarray(1)%p,ilwork,1,iper,iBC,BC)

!      do i=1,ramg_levelx
!         write(fname,'((A8)(I1))')'tarray_i',i
!         call ramg_dump_rn_map(tarray(i)%p,amg_nshg(i),1,fname)
!      enddo


      do i=1,ramg_levelx
         deallocate(tarray(i)%p)
      enddo

      end subroutine ! ramg_check_paracomm
     
!**************************************************
!      ramg_zeroOut: zero out slave values
!**************************************************
      subroutine ramg_zeroOut(u,ilwork,BC,iBC,iper)
      include "common.h"
      include "auxmpi.h"
      include "mpif.h"
      
      real(kind=8),dimension(nshg),intent(inout)  :: u
      integer,intent(in),dimension(nlwork)        :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
     
      integer i,j,k,p
      integer            :: itag,iacc,iother,numseg,isgbeg,itkbeg
      integer            :: numtask,lenseg

      if (numpe.lt.2) return
      !call commOut(u,ilwork,1,iper,iBC,BC)
      numtask = ilwork(1)
      itkbeg = 1
      do i=1,numtask
         iacc = ilwork(itkbeg+2)
         numseg = ilwork(itkbeg+4)
         if (iacc.eq.0) then
         do j=1,numseg
            isgbeg = ilwork(itkbeg+3+j*2)
            lenseg = ilwork(itkbeg+4+j*2)
            isgend = isgbeg + lenseg -1
            u(isgbeg:isgend) = 0
         enddo
         endif
         itkbeg = itkbeg+4+2*numseg
      enddo

      end subroutine ! ramg_zeroOut

!************************************************
!      ramg_freeBC: set "fine" on boundary nodes
!************************************************
      subroutine ramg_freeBC(amg_F,ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "auxmpi.h"
      include "mpif.h"
      integer,intent(inout),dimension(nshg)  :: amg_F
      integer,intent(in),dimension(nlwork)     :: ilwork
      integer,intent(in),dimension(nshg)       :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      integer      :: itag,iacc,iother,numseg,isgbeg,itkbeg,numtask
      integer      :: i,p
      integer,dimension(nshg,2)  :: tmpout

      character     :: filename*80
      integer       :: nentry,nshg,iglobal,ilocal,icf,iproc

      return 
      
      if (numpe.ge.2) then
      
          numtask = ilwork(1)
          itkbeg = 1
          do i=1,numtask
             itag = ilwork(itkbeg+1)
             iacc = ilwork(itkbeg+2)
             iother = ilwork(itkbeg+3)
             numseg = ilwork(itkbeg+4)
             isgbeg = ilwork(itkbeg+5)
             do j=1,numseg
                p = ilwork(itkbeg+3+j*2)
                amg_F(p) = 2
             enddo
          enddo

      tmpout(:,1) = ncorp_map
      tmpout(:,2) = amg_F

      call ramg_dump_i(tmpout,nshg,2,'CFsplit   ')
      else  ! DEBUGGING PART, READ 2-PROC CASE INTO 1-proc
          filename = "cfs.dat"
          filename = trim(filename)
          open(unit=999,file=filename,status='unknown')
          read(999,*)nentry
          do i=1,nentry
             read(999,*)iglobal,ilocal,icf,iproc
             !amg_F(iglobal) = icf
          enddo
          close(999)
      end if

      end subroutine ! ramg_freeBC


!**********************************************
!      ramg_read_map:
!       read in ncorp array from local to global
!       mapping
!**********************************************
      subroutine ramg_prep_reduce_map
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      character*5 cname
      character*255 fname1
      integer igeom,map_nshg,i

      allocate(ncorp_map(nshg))
      if (numpe.lt.2) then! construct dummy-array
          do i=1,nshg
             ncorp_map(i) = i
          enddo
      else
      fname1='geombc.dat'
      fname1=trim(fname1)//cname(myrank+1)

      call openfile(fname1,"read",igeom)
      !write(*,*)'mcheck:',myrank,fname1,igeom

      fname1="mode number map from partition to global"
      ione=1
      call readheader(igeom,
     & 'mode number map from partition to global' // char(0),
     &  map_nshg,ione,"integer",iotype)
      
      call readdatablock(igeom,
     & 'mode number map from partition to global' // char(0),
     & ncorp_map,map_nshg, "integer",iotype)

      call closefile(igeom,"read")

      endif
      
      end subroutine ! ramg_read_map

!***********************************************
!      ramg_check_diag
!***********************************************
      subroutine ramg_check_diag(colm,rowp,lhs,nshg,nnz_tot)
      implicit none
      integer,intent(in)                   :: nshg,nnz_tot
      integer,intent(in),dimension(nshg+1) :: colm
      integer,intent(in),dimension(nnz_tot) :: rowp
      real(kind=8),intent(in),dimension(nnz_tot) :: lhs

      integer                              :: i,j,p
      real(kind=8)                         :: diagrow
      logical                              :: diagokay
      
      diagokay = .true.
      
      do i=1,nshg
         p = colm(i)
         if (rowp(p).ne.i) then
             write(*,*)'matrix first row entry is not diagonal',i
             diagokay = .false.
         end if
         diagrow = lhs(p)
         if (diagrow.lt.0) then
             write(*,*)'matrix diagonal < 0 at',i,diagrow
             diagokay = .false.
         end if
         do j = colm(i)+1,colm(i+1)-1
            p = rowp(j)
            if (lhs(p).gt.diagrow) then
                write(*,*)'matrix not diagonal dominant at row',i
                diagokay = .false.
                write(*,*)'diag:',diagrow,p,lhs(p)
                exit
            end if
         end do
      enddo

      if (diagokay) then
          write(*,*)'matrix check diagonal dominance okay'
      end if
      
      end subroutine

      subroutine ramg_dot_p(level,v,u,redun,norm)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      
      integer :: redun,level
      real(kind=8),intent(in),dimension(amg_nshg(level),redun) ::v,u
      real(kind=8),intent(inout)    :: norm
      integer                       :: i,k
      real(kind=8) :: normt
      normt = 0
      do i=1,amg_nshg(level)
         do k=1,redun
         normt = normt + v(i,k)*u(i,k)
         enddo
      enddo
      if (numpe.gt.1) then
      call MPI_AllReduce(normt,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     & MPI_COMM_WORLD,ierr)
      else
          norm=normt
      endif
      end subroutine !ramg_dot_p
   
!************************************************************
!      ramg_L2_norm
!      calculate norm = | r | 
!************************************************************      
      subroutine ramg_L2_norm(nshg,v,norm)
      implicit none
      integer,intent(in)           :: nshg
      real(kind=8),intent(in),dimension(nshg) :: v
      real(kind=8),intent(inout)    :: norm
      integer                       :: i
      norm = 0
      do i=1,nshg
         norm = norm + v(i)*v(i)
      enddo
      norm = sqrt(norm)
      end subroutine !ramg_L2_norm
  
      subroutine ramg_L2_norm_p(level,v,vflow,norm)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      
      integer,intent(in) :: vflow,level
      real(kind=8),intent(in),dimension(amg_nshg(level),vflow) :: v
      real(kind=8),intent(inout)    :: norm
      integer                       :: i,k
      real(kind=8) :: normt
      normt = 0
      do i=1,amg_nshg(level)
         do k=1,vflow
         normt = normt + v(i,k)*v(i,k)
         enddo
      enddo
      if (numpe.gt.1) then
      call MPI_AllReduce(normt,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     & MPI_COMM_WORLD,ierr)
      else
          norm = normt
      endif
      norm = sqrt(norm)
      end subroutine !ramg_L2_norm_p
      
!************************************************************
!   ramg_allocate & deallocate
!    (de)allocate memory for level l, lhs matrix, rhs vector
!************************************************************
      subroutine ramg_allocate(
     &                level,lnshg,lnnz_tot,n_sol)
      use ramg_data
      implicit none

      integer,intent(in)             :: level,lnshg,lnnz_tot
      integer,intent(in)             :: n_sol
      integer                        :: mem_err,mem_err_s
      mem_err_s = 0
      if (lnnz_tot.ne.0) then ! zero means manual alloc later
      amg_nnz(level) = lnnz_tot
      allocate(amg_A_lhs(level)%p(amg_nnz(level),n_sol),
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      allocate(amg_A_rowp(level)%p(amg_nnz(level)),
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      endif
      if (lnshg.ne.0) then
      amg_nshg(level) = lnshg
      allocate(amg_A_colm(level)%p(amg_nshg(level)+1),
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      allocate(amg_A_rhs(level)%p(amg_nshg(level)),
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      allocate(amg_ppeDiag(level)%p(amg_nshg(level)),
     &        stat=mem_err)
      endif
      mem_err_s = mem_err_s + mem_err
      if (mem_err_s .ne. 0 ) then
          write(6,7001)level
      end if
7001  format(/' **** ramg: Allocation error at level',i5)
      ! zero out
      if (lnnz_tot.ne.0) then
      amg_A_lhs(level)%p(:,:)  = 0
      amg_A_rowp(level)%p(:) = 0
      endif
      if (lnshg.ne.0) then
      amg_A_colm(level)%p(:) = 0
      amg_A_rhs(level)%p(:)  = 0
      endif
      return

      end subroutine ! ramg_allocate

      subroutine ramg_deallocate(level)
      
      use ramg_data
      include "common.h"
      
      integer,intent(inout)            :: level
      integer                        :: mem_err,mem_err_s,i
      

      if (level.eq.1) then

      if (maxnev.gt.0) then
        deallocate(ramg_ggb_ev)
        deallocate(ramg_ggb_eA)
        deallocate(cmtxA)
        deallocate(cindx)
      endif


        if (iamg_reduce.gt.1) then
            deallocate(reducemap)
            deallocate(rmap1d)
        endif
 
      end if


      level = abs(level)
      
      mem_err_s = 0
      if (associated(amg_A_lhs(level)%p)) then
      deallocate(amg_A_lhs(level)%p,
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      endif
      if (associated(amg_A_rowp(level)%p)) then
      deallocate(amg_A_rowp(level)%p,
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      endif
      if (associated(amg_A_colm(level)%p)) then
      deallocate(amg_A_colm(level)%p,
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      endif
      if (associated(amg_ppeDiag(level)%p)) then
      deallocate(amg_ppeDiag(level)%p,
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      endif
      if (associated(amg_A_rhs(level)%p)) then
      deallocate(amg_A_rhs(level)%p,
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      endif
      if (associated(amg_ilwork(level)%p)) then
      deallocate(amg_ilwork(level)%p,
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      write(*,*)'mcheck deallocate ilwork,',level,myrank
      endif

      if (associated(amg_paramap(level)%p)) then
      deallocate(amg_paramap(level)%p,
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      endif
      if (associated(amg_paraext(level)%p)) then
      deallocate(amg_paraext(level)%p,
     &         stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      endif

      if (mem_err_s .ne. 0 ) then
          write(6,7002)level
      end if

      if (associated(I_fc_colm(level)%p)) then
      deallocate(CF_map(level)%p,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      deallocate(CF_revmap(level)%p,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      deallocate(I_cf_colm(level)%p,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      deallocate(I_cf_rowp(level)%p,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      deallocate(I_cf(level)%p,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      deallocate(I_fc_colm(level)%p,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      deallocate(I_fc_rowp(level)%p,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      deallocate(I_fc(level)%p,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      end if
      amg_nnz(level) = 0
      amg_nshg(level) = 0
7002  format(/' **** ramg: Deallocation error at level',i5)

      return
     
      end subroutine ! ramg_deallocate

      subroutine ramg_readin_i(iarray,rn1,rfname)
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      integer rn1
      integer,dimension(rn1) ::  iarray
      character*10     rfname
      character*5      cname
      
      character*20     mfname

      integer i,t
      mfname = trim(rfname)//'.dat'//cname(myrank+1)

      write(*,*)'RAMG READ: ',mfname
      open(264,file=trim(mfname),status='unknown')
      do i=1,rn1
         read(264,*)t,iarray(i)
      enddo
      close(264)
      end subroutine ! ramg_readin_i
      

      subroutine ramg_dump(rarray,rn1,rfname)
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      integer rn1
      real(kind=8),dimension(rn1) ::  rarray
      character*10     rfname
      character*5      cname
      
      character*20     mfname

      integer i
      mfname = trim(rfname)//'.dat'//cname(myrank+1)

      write(*,*)'RAMG DUMP: ',mfname
      open(264,file=trim(mfname),status='unknown')
      do i=1,rn1
         !write(264,'((I10)(A)(D10.3))')i,'  ',rarray(i)
         write(264,*)i,rarray(i)
      enddo
      close(264)
      
      end subroutine ! ramg_dump
 
      subroutine ramg_dump_ic(rarray,rn1,rn2,rfname)
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      integer rn1,rn2
      integer,dimension(rn1,rn2) ::  rarray
      character*10     rfname
      character*5      cname
      
      character*20     mfname

      integer i
      mfname = trim(rfname)//'.dat'//cname(myrank+1)

      write(*,*)'RAMG DUMP: ',mfname
      open(264,file=trim(mfname),status='unknown')
      write(264,*)rn2
      do j=1,rn2
         write(264,*)j,(rarray(i,j),i=1,rn1)
      enddo
      close(264)
      end subroutine
     
      subroutine ramg_dump_i(rarray,rn1,rn2,rfname)
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      integer rn1,rn2
      integer,dimension(rn1,rn2) ::  rarray
      character*10     rfname
      character*5      cname
      
      character*20     mfname

      integer i
      mfname = trim(rfname)//'.dat'//cname(myrank+1)

      write(*,*)'RAMG DUMP: ',mfname
      open(264,file=trim(mfname),status='unknown')
      write(264,*)rn1
      do i=1,rn1
         write(264,*)i,(rarray(i,j),j=1,rn2)
      enddo
      close(264)
      
      end subroutine ! ramg_dump_i
 
      subroutine ramg_dump_ir(iarray,rarray,rn1,rn2,rfname)
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      integer rn1,rn2
      integer,dimension(rn1) ::  iarray
      real(kind=8),dimension(rn1,rn2) ::  rarray
      character*10     rfname
      character*5      cname
      
      character*20     mfname
      
      character*20 mformat

      integer i
      mfname = trim(rfname)//'.dat'//cname(myrank+1)

      write(mformat,'((A6)(I1)(A7))')'((2I)(',rn2,'E14.3))'
      mformat = trim(mformat)
      write(*,*)'RAMG DUMP: ',mfname
      open(264,file=trim(mfname),status='unknown')
      do i=1,rn1
         write(264,mformat)
     &   i,iarray(i),(rarray(i,j),j=1,rn2)
      enddo
      close(264)
      
      end subroutine ! ramg_dump_ir
 
      subroutine ramg_dump_rn_map(rarray,rn1,rn2,rfname)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      integer rn1,rn2
      real(kind=8),dimension(rn1,rn2) ::  rarray
      character*10     rfname
      character*5      cname
      
      character*20     mfname
      character*20 mformat

      integer i,j,ii
      mfname = trim(rfname)//'.dat'//cname(myrank+1)

      write(*,*)'RAMG DUMP: ',mfname
      open(264,file=trim(mfname),status='unknown')
      write(mformat,'((A6)(I1)(A7))')'((1I)(',rn2,'E14.4))'
      mformat=trim(mformat)
      do i=1,rn1
         if (numpe.gt.1) then
             ii = ncorp_map(i)
         else 
             ii = i
         endif
         write(264,mformat)ii,(rarray(i,j),j=1,rn2)
      enddo
      close(264)
      end subroutine ! ramg_dump_rn_map

      subroutine ramg_dump_rn(rarray,rn1,rn2,rfname)
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      integer rn1,rn2
      real(kind=8),dimension(rn1,rn2) ::  rarray
      character*10     rfname
      character*5      cname
      
      character*20     mfname
      character*20 mformat

      integer i,j
      mfname = trim(rfname)//'.dat'//cname(myrank+1)

      write(*,*)'RAMG DUMP: ',mfname
      open(264,file=trim(mfname),status='unknown')
      write(mformat,'((A6)(I1)(A7))')'((1I)(',rn2,'E14.4))'
      mformat=trim(mformat)
      do i=1,rn1
         write(264,mformat)i,(rarray(i,j),j=1,rn2)
      enddo
      close(264)
      
      end subroutine ! ramg_dump_rn

      subroutine ramg_dump_matlab_A(acolm,arowp,alhs,an,annz,redun,
     &           fname)

      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer :: an,annz,redun
      integer,dimension(an+1) :: acolm
      integer,dimension(annz) :: arowp
      real(kind=8),dimension(redun,annz) :: alhs
      character fname*10,mtxname*5
      character cname*5

      character mfname*15,mAname*5
      character mformat*20

      integer i,j,p,k
      
      mfname = trim(fname)//'.dat'//cname(myrank+1)

      write(*,*)'RAMG Dump to matlab  ',mfname,myrank
      open(264,file=trim(mfname),status='unknown')
      !write(264,*)an
      !write(264,*)annz
      !write(264,*)'1'
      write(mformat,'((A6)(I1)(A7))')'((2I)(',redun,'E14.4))'
      do i=1,an
         do p=acolm(i),acolm(i+1)-1
            j = arowp(p)
            !write(264,"((I6)(I7)(A)(E8.2))")i,j,' ',alhs(p)
!            if ( (amg_paramap(1)%p(i).eq.amg_paramap(2)%p(j)).and.
!           if   (iabs(amg_paramap(1)%p(i)).ne.(myrank+1)) then
!     &          (amg_paramap(1)%p(i).ne.1)) then
            write(264,mformat)i,j,(alhs(k,p),k=1,redun)
!            endif
         enddo
      enddo
      close(264)
      end subroutine

      subroutine ramg_dump_matlab_map(acolm,arowp,alhs,an,annz,redun,
     &           fname)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer :: an,annz,redun
      integer,dimension(an+1) :: acolm
      integer,dimension(annz) :: arowp
      real(kind=8),dimension(redun,annz) :: alhs
      !integer,dimension(an) :: ipmap
      character fname*10

      character mfname*20
      character cname*5
      character mformat*20

      integer i,j,p,k
      integer ii,jj

      if (numpe.eq.1) then
          call ramg_dump_matlab_A(acolm,arowp,alhs,an,annz,redun,
     &         fname)
         return;
      endif
      
      mfname = trim(fname)//'.dat'//cname(myrank+1)

      write(*,*)'RAMG Dump to matlab  ',mfname,nshg,myrank
      open(264,file=trim(mfname),status='unknown')
      !write(264,*)an
      !write(264,*)annz
      !write(264,*)'1'
      write(mformat,'((A6)(I1)(A7))')'((2I)(',redun,'E14.4))'
      !call ramg_dump_i(ncorp_map,nshg,1,'pam_corp  ')
      do i=1,an
         ii = ncorp_map(i)!ipmap(i))
         do p=acolm(i),acolm(i+1)-1
            j = arowp(p)
            jj = ncorp_map(j)!ipmap(j))
            !write(264,"((I6)(I7)(A)(E8.2))")i,j,' ',alhs(p)
            write(264,mformat)ii,jj,(alhs(k,p),k=1,redun)
         enddo
      enddo
      close(264)
      end subroutine

      
      subroutine ramg_dump_mupad_A(acolm,arowp,alhs,an,annz,
     &           fname,mtxname)
      integer :: an,annz
      integer,dimension(an+1) :: acolm
      integer,dimension(annz) :: arowp
      real(kind=8),dimension(annz) :: alhs
      character fname*10,mtxname*5

      character mfname*15,mAname*5

      integer i,j,p,k
      
      mfname = trim(fname)//'.mws'
      mAname = trim(mtxname)

      write(*,*)'RAMG Dump to Mupad  ',mfname
      open(264,file=trim(mfname),status='unknown')
      write(264,*)mAname,':= matrix(',an,',',an,'):'
      do i=1,an
         do p=acolm(i),acolm(i+1)-1
            j = arowp(p)
            write(264,*)mAname,'[',i,',',j,']:=',alhs(p),':'
         enddo
      enddo
      close(264)
      
      end subroutine

      subroutine ramg_print_time(str,v)

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      character*(*) str
      real(kind=8) :: v,vave,vmax

      if (numpe.gt.1) then
          call MPI_AllReduce(v,vmax,1,
     &    MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
          call MPI_AllReduce(v,vave,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
          vave = vave/numpe
      else
          vave = v
          vmax = v
      endif

      if ((iamg_verb.gt.1).and.(myrank.eq.master)) then
      write(*,7101)trim(str),vave,vmax
7101  format(T1,A,T40,f9.2,' s (ave), ',f9.2,' s (max)')
      endif

      end subroutine

      subroutine ramg_output_coarsening
      use readarrays
      use ramg_data
      include "common.h"
      include "mpif.h"

      character*20 mfname
      character*20 mformat

      integer i

      write(*,*)'ramg dump coarsening'
      mfname = 'amgcoarsen.dat'
      open(265,file=trim(mfname),status='unknown')
      write(265,*)nshg
      do i=1,nshg
      write(265,'((2I)(3E14.5))')i,amg_cfmap(i),
     & point2x(i,1),point2x(i,2),point2x(i,3)
      enddo
      close(265)

      end subroutine
     
!***********************************************************
!      <EOF> ramg TOOLS
!**********************************************************      
