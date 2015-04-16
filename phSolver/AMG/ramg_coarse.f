!*************************************************************
!     RAMG Coarse: 
!      do C/F spliting and setup coarse matrix
!**************************************************************
!     ramg: ramg_coarse_setup
!      Input: amg_A matrix level1
!      Output: I_h_H matrix level2 to level1
!              I_H_h matrix level1 to level2
!**************************************************************
      subroutine ramg_coarse_setup(level1,level2,eps_str,
     &           interp_flag,trunc_tol,
     &           ilwork,BC,iBC,iper)!,nshg,nlwork,ndofBC)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      
!      implicit none
      
      !******* parameters *******
      integer,intent(in)                        :: level1,level2
      real(kind=8),intent(in)                   :: eps_str
      integer,intent(in)                        :: interp_flag
      real(kind=8),intent(in)                   :: trunc_tol
!      integer,intent(in)                       :: nshg,nlwork,ndofBC
      integer,intent(in),dimension(nlwork)     :: ilwork
      integer,intent(in),dimension(nshg)       :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      !******* temp variables *******
      integer,dimension(amg_nnz(level1))        :: amg_S,amg_Ip
      integer,dimension(amg_nshg(level1)) :: amg_F,aLoc
      real(kind=8),dimension(amg_nshg(level1))  :: amg_la
      integer,dimension(amg_nshg(level1))       :: amg_Fn,amg_Fp
      real(kind=8)              :: alpha,beta,alphac,betac,diag,rtp
      real(kind=8)              :: tmp1,tmp2,tmp3,tmp4
      ! AMG_F : 0: UNDECIDED, 1: COARSE, 2: FINE
      integer                                   :: I_nnz,ki,kj,jj
      integer                                   :: i,j,k,m,n,p,q
      character                                 :: fname*80

      integer            :: numtask,itkbeg,numseg,iacc
      integer,allocatable,dimension(:) :: subcf,subcfrev,subnei
    
      integer                          :: mnnz
      integer                          :: cfilter
         type(i2dd)     :: amg_I_rowp
         type(r2dd)     :: amg_I
         type(i1d)      :: amg_I_colm

    
      real,dimension(10)                        :: cpusec
      logical                                   :: ti_flag
      integer                                   :: mem_err,mem_err_s
      !******* end of variables ********

      call cpu_time(cpusec(1))

      if (abs(trunc_tol).lt.1e-4) then
          ti_flag = .false.
      else
          ti_flag = .true.
      end if

      amg_S = 0
      amg_F = 0

      if (numpe.ge.2) then 
          numtask = ilwork(1)+1
      else
          numtask = 1
      endif

      IF ((iamg_reduce.le.1).or.(numpe.gt.1)) THEN
!      IF (.TRUE.) THEN
      allocate(subcfrev(numpe))
      subcfrev = 0
      allocate(subcf(numtask))
      subcf = 0
      allocate(subnei(numtask))
      subnei = 0

      subnei(1) = myrank+1
      subcfrev(subnei(1))=1
      itkbeg = 1
      do i=2,numtask
         subnei(i)=ilwork(itkbeg+3)+1 ! iother+1
         subcfrev(subnei(i)) = i
         itkbeg = itkbeg + 4 + 2*ilwork(itkbeg+4)
      enddo

      ELSE !reduced
          numtask = rmapmax-1
          allocate(subcfrev(numtask))
          allocate(subcf(numtask))
          allocate(subnei(numtask))
          subcfrev = 0
          subcf = 0
          subnei = 0
          do i=1,numtask
             subcfrev(i) = i
          enddo
      END IF

      call ramg_CFsplit(amg_F,amg_S,amg_nshg(level1),amg_nnz(level1),
     &          amg_A_colm(level1)%p,amg_A_rowp(level1)%p,
     &          amg_A_lhs(level1)%p,amg_paramap(level1)%p,
     &          ilwork,BC,iBC,iper,
     &          eps_str,-1,0)
      ! read FC spliting to file
      if (level1.eq.1) then
      !    write(fname,'((A7)(I1))')'GTG_FC_',level1
      !    call ramg_readin_i(amg_F,amg_nshg(level1),fname)
      ! write FC spliting to file
      !write(fname,'((A7)(I1))')'amg_FC_',level1
      !call ramg_dump_i(amg_F,amg_nshg(level1),1,fname)
      endif

      call ramg_update_cfmap(amg_F,level1)
      ! communication goes here
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call commOut_i(amg_cfmap,ilwork,1,iper,iBC,BC)
      call ramg_readin_cfmap(amg_F,level1)

      IF (.true.) THEN
      subcf = 0
      subnei = 0
      do i = 1,amg_nshg(level1)
         p = iabs(amg_paramap(level1)%p(i))
         p = subcfrev(p)
         subnei(p) = subnei(p) + 1
         if (amg_F(i).eq.1) then
         subcf(p) = subcf(p) + 1
         endif
      enddo
      do i = 1,numtask
      enddo
      ENDIF

      call cpu_time(cpusec(3))

      deallocate(subcf)
      deallocate(subcfrev)
      deallocate(subnei)

      ! SETUP Smoothing ordering
      allocate(CF_map(level1)%p(amg_nshg(level1)))
      allocate(CF_revmap(level1)%p(amg_nshg(level1)))

      k = 1
      do i=1,amg_nshg(level1)
         if (amg_F(i).eq.1) then ! fine/coarse first
             CF_map(level1)%p(k) = i
             CF_revmap(level1)%p(i) = k
             k = k+1
         end if
      enddo
      do i=1,amg_nshg(level1)
         if (amg_F(i).eq.2) then
             CF_map(level1)%p(k) = i
             CF_revmap(level1)%p(i) = k
             k = k+1
         end if
      enddo

      if (level1.eq.-1) then
      do i=1,amg_nshg(level1)
         CF_map(level1)%p(i) = i
         CF_revmap(level1)%p(i) = i
      enddo
      endif

      ! generate I := I_H_h, from coarse to fine
      ! then TRANSPOSE and get I^T, from fine to coarse

      ! determine the size of second level
      amg_nshg(level2) = 0
      amg_Ip = 0
      I_nnz = 0
      do i = 1, amg_nshg(level1)
         if (amg_F(i).eq.1) then
             amg_nshg(level2) = amg_nshg(level2) + 1
             I_nnz = I_nnz + 1
             amg_Ip(i) = I_nnz
         end if
      enddo

      allocate(amg_paramap(level2)%p(amg_nshg(level2)))
      allocate(amg_paraext(level2)%p(amg_nshg(level2)))
      j = 1
      do i=1,amg_nshg(level1)
         if (amg_F(i).eq.1) then
             amg_paramap(level2)%p(j) = amg_paramap(level1)%p(i)
             amg_paraext(level2)%p(j) = amg_paraext(level1)%p(i)
             j = j + 1
         end if
      enddo

      ! ALLOCATE I_H_h from coarse to fine
      mem_err_s = 0
      allocate(amg_I_colm%p(amg_nshg(level1)),stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      allocate(amg_I_rowp%pp(amg_nshg(level1)),stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      allocate(amg_I%pp(amg_nshg(level1)),stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      if (mem_err_s.ne.0) then
          write(*,*)'Memory allocation error while geting I_H_h'
      end if

      cpusec(4) = 0
      cpusec(5) = 0
      cpusec(8) = 0
      mnnz = 0
      aLoc = 0
      loop_i: do i = 1, amg_nshg(level1)
         call cpu_time(cpusec(6))
         cfilter = amg_paramap(level1)%p(i)
         if (amg_F(i).eq.1) then ! i as a coarse node do trivial work
             allocate(amg_I_rowp%pp(i)%p(1))
             allocate(amg_I%pp(i)%p(1))
             amg_I_colm%p(i) = 1
             amg_I_rowp%pp(i)%p(1) = amg_Ip(i)
             amg_I%pp(i)%p(1) = 1
             mnnz = mnnz + 1
             cycle loop_i
         end if
         n = 0
         ! a-hat, P-hat, N-hat
         ! amg_la  ! a-hat
         ! amg_Fn  ! N-hat
         ! amg_Fp  ! P-hat
         ! diagonal part
         diag = -1.0d0/amg_A_lhs(level1)%p(amg_A_colm(level1)%p(i),1)
         ! non-diagonal 
         do k = amg_A_colm(level1)%p(i)+1,amg_A_colm(level1)%p(i+1)-1
            j = amg_A_rowp(level1)%p(k)
         if (cfilter.eq.amg_paramap(level1)%p(j)) then
            n = n + 1
            if ( (amg_F(j).eq.1) .and. (mod(amg_S(k),2).eq.1) ) then
                amg_Fp(n) = 1
            else
                amg_Fp(n) = 0
            endif
            amg_Fn(n) = j
            amg_la(n) = amg_A_lhs(level1)%p(k,1)
            aLoc(j) = n
         end if
         enddo
         
         ! standard interpolation here ! for parallel dont use
         if ( interp_flag .eq. 1) then
         do k=amg_A_colm(level1)%p(i)+1,amg_A_colm(level1)%p(i+1)-1
          j = amg_A_rowp(level1)%p(k)
          if ( (amg_F(j).eq.2) .and. (mod(amg_S(k),2).eq.1)) then
             rtp = amg_A_lhs(level1)%p(amg_A_colm(level1)%p(j),1)
             if (rtp.ne.0) then
             rtp = amg_la(aLoc(j))/rtp
             endif
           do kj=amg_A_colm(level1)%p(j),amg_A_colm(level1)%p(j+1)-1
             jj = amg_A_rowp(level1)%p(kj)
             m = aLoc(jj)
             if (m.eq.0) then
               m = n+1
               aLoc(jj) = m
               amg_Fn(m) = jj
               amg_Fp(m) = 0
               amg_la(m) = 0
               n = n+1
             end if
             amg_la(m) = amg_la(m) - rtp*amg_A_lhs(level1)%p(kj,1)
             if ( (amg_F(jj).eq.1) .and. (mod(amg_S(kj),2).eq.1) ) then
                amg_Fp(m) = 1
             else
                amg_Fp(m) = 0
             endif
           enddo
          end if
         enddo
         endif
         
         call cpu_time(cpusec(7))
         cpusec(4) = cpusec(4) + cpusec(7)-cpusec(6)
         ! calculate alpha, beta, nnz-in-row
         alpha = 0
         beta = 0
         alphac = 0
         betac = 0
         I_nnz = 0
         do k = 1,n ! all non diagonal
            aLoc(amg_Fn(k))=0
            ! in N-hat
                if ( amg_la(k).lt.0 ) then ! a-hat < 0
                     alpha = alpha + amg_la(k)
                else                       ! a-hat > 0
                     beta = beta + amg_la(k)
                end if
                ! in P-hat
                if (amg_Fp(k).eq.1) then
                    I_nnz = I_nnz + 1     ! nonzero in I_h_H row
                    if (amg_la(k).lt.0) then ! a-hat < 0 in coarse
                        alphac = alphac + amg_la(k)
                    else                     ! a-hat > 0 in coarse
                        betac = betac + amg_la(k)
                    end if
                    amg_Fn(I_nnz) = amg_Fn(k)
                    amg_la(I_nnz) = amg_la(k)
                end if
         enddo
         alpha = diag*alpha/alphac
         beta = diag*beta/betac
         rtp = 0
         do k = 1,I_nnz
            if ( amg_la(k).lt.0) then 
                amg_la(k) = amg_la(k)*alpha
            else
                amg_la(k) = amg_la(k)*beta
            end if
            rtp = max(rtp,abs(amg_la(k)))
         enddo
         if (ti_flag) then
             n = I_nnz
             I_nnz = 0
             tmp1 = 0
             tmp2 = 0
             tmp3 = 0
             tmp4 = 0
             rtp = rtp*trunc_tol
             do k=1,n
                tmp3 = tmp3 + min(amg_la(k),zero)
                tmp4 = tmp4 + max(amg_la(k),zero)
                if (abs(amg_la(k)).gt.rtp) then
                    tmp1 = tmp1 + min(amg_la(k),zero)
                    tmp2 = tmp2 + max(amg_la(k),zero)
                    I_nnz = I_nnz + 1
                    amg_Fn(I_nnz) = amg_Fn(k)
                    amg_la(I_nnz) = amg_la(k)
                end if
             enddo
             if ( abs(tmp1).gt.1e-5) tmp1 = tmp3/tmp1
             if ( abs(tmp2).gt.1e-5) tmp2 = tmp4/tmp2
             do k=1,I_nnz
                if (amg_la(k).lt.0) then
                    amg_la(k) = amg_la(k)*tmp1
                else
                    amg_la(k) = amg_la(k)*tmp2
                end if
             enddo
         end if
         ! put up I matrix row
         call cpu_time(cpusec(6))
         cpusec(5) = cpusec(5) + cpusec(6)-cpusec(7)
         mnnz = mnnz + I_nnz
         amg_I_colm%p(i) = I_nnz
         allocate(amg_I_rowp%pp(i)%p(I_nnz))
         allocate(amg_I%pp(i)%p(I_nnz))
         do k = 1,I_nnz ! scan over P-hat
                amg_I_rowp%pp(i)%p(k) = amg_Ip(amg_Fn(k))
                ! put in rowp
                if (amg_la(k).lt.0) then
                   amg_I%pp(i)%p(k) = amg_la(k)
                else
                    amg_I%pp(i)%p(k) = amg_la(k)
                endif
         enddo
         call cpu_time(cpusec(7))
         cpusec(8) = cpusec(8) + cpusec(7) - cpusec(6)
      enddo loop_i ! Now we have I_H_h

      ! copy the data
      mem_err_s = 0
      allocate(I_cf_colm(level1)%p(amg_nshg(level1)+1),stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      allocate(I_cf_rowp(level1)%p(mnnz),stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      allocate(I_cf(level1)%p(mnnz),stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      if (mem_err_s .ne. 0 ) then
          write(*,*)'allocation error in I_cf'
      end if
      mnnz = 0
      do i=1,amg_nshg(level1)
         I_cf_colm(level1)%p(i)=mnnz+1
         do j=1,amg_I_colm%p(i)
            I_cf_rowp(level1)%p(mnnz+j)=amg_I_rowp%pp(i)%p(j)
            I_cf(level1)%p(mnnz+j)=amg_I%pp(i)%p(j)
         enddo
         mnnz = mnnz + amg_I_colm%p(i)
      enddo
      I_cf_colm(level1)%p(amg_nshg(level1)+1) = mnnz+1

      ! dump interpolator
!      if (level1.eq.1) then
!      write(fname,'((A4)(I1))')'Icf_',level1
!      write(*,*)amg_nshg(level1),mnnz
!      call ramg_dump_matlab_A(I_cf_colm(level1)%p,I_cf_rowp(level1)%p,
!     &     I_cf(level1)%p,amg_nshg(level1),mnnz,1,fname)
!      endif
      
      mem_err_s = 0
      do i = 1,amg_nshg(level1)
         if (amg_I_colm%p(i).ne.0) then
             deallocate(amg_I_rowp%pp(i)%p,stat=mem_err)
             mem_err_s = mem_err_s + mem_err
             deallocate(amg_I%pp(i)%p,stat=mem_err)
             mem_err_s = mem_err_s + mem_err
         end if
      enddo
      deallocate(amg_I_colm%p,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      deallocate(amg_I_rowp%pp,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      deallocate(amg_I%pp,stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      if (mem_err_s .ne. 0) then
          write(*,*)'deallocation error in I_cf'
      end if

      ! now build I_fc, inverse of I_cf
      call cpu_time(cpusec(6))

      mem_err_s = 0
      allocate(I_fc_colm(level1)%p(amg_nshg(level2)+1),stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      allocate(I_fc_rowp(level1)%p(mnnz),stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      allocate(I_fc(level1)%p(mnnz),stat=mem_err)
      mem_err_s = mem_err_s + mem_err
      if (mem_err_s .ne. 0 ) then
          write(*,*)'allocation error in I_fc'
      end if

      I_fc_colm(level1)%p(1:amg_nshg(level2)+1) = 0
      do i=1,mnnz
         I_fc_colm(level1)%p(I_cf_rowp(level1)%p(i)) =
     &   I_fc_colm(level1)%p(I_cf_rowp(level1)%p(i)) + 1
      enddo
      mnnz = 1
      do i=1,amg_nshg(level2)
         j = I_fc_colm(level1)%p(i)
         I_fc_colm(level1)%p(i) = mnnz
         mnnz = mnnz + j
      enddo
      I_fc_colm(level1)%p(amg_nshg(level2)+1) = mnnz
      
      do i=1,amg_nshg(level1)
         do k=I_cf_colm(level1)%p(i),I_cf_colm(level1)%p(i+1)-1
            j = I_cf_rowp(level1)%p(k)
            kj = I_fc_colm(level1)%p(j)
            I_fc_colm(level1)%p(j) = I_fc_colm(level1)%p(j) + 1
            I_fc_rowp(level1)%p(kj) = i
            I_fc(level1)%p(kj) = I_cf(level1)%p(k)
         enddo
      enddo

      do i=amg_nshg(level2),2,-1
         I_fc_colm(level1)%p(i) = I_fc_colm(level1)%p(i-1)
      enddo
      I_fc_colm(level1)%p(1) = 1

      call cpu_time(cpusec(10))

      end subroutine!ramg_coarse_setup

!**************************************************************
! Heap routines
!**************************************************************
      subroutine ramg_initheap(heap,invMap,wght,nheaps,ilen)
      implicit none
      integer :: nheaps,ilen
      integer,dimension(0:ilen-1),intent(inout) :: heap,invMap
      real(kind=8),dimension(0:ilen-1),intent(in)    :: wght

      integer :: i,j,k,t

      do i=0,nheaps-1
        heap(i) = i
      enddo

      do i=1,nheaps-1
         k = i
         j = ishft(k-1,-1)
         do while ( ( k.gt.0 ).and.(wght(heap(k)).gt.wght(heap(j))) )
            t = heap(j)
            heap(j) = heap(k)
            heap(k) = t
            k = j
            j = ishft(j-1,-1)
         enddo
      enddo

      do i=0,nheaps-1
        invmap(heap(i)) = i
      enddo
      
      end subroutine ! ramg_initheap

      subroutine ramg_popheap(heap,invmap,wght,nheaps,popid,ilen)
      implicit none
      integer,intent(inout) :: nheaps
      integer,intent(in)    :: popid,ilen
      real(kind=8),dimension(0:ilen-1),intent(in) :: wght
      integer,dimension(0:ilen-1),intent(inout) :: heap,invmap

      integer i,j,k
      nheaps = nheaps-1
      heap(popid) = heap(nheaps)
      invmap(heap(popid)) = popid
      call ramg_adjheap(heap,invmap,wght,nheaps,popid,ilen)

      end subroutine !ramg_popheap

      subroutine ramg_adjheap(heap,invmap,wght,nheaps,popid,ilen)
      implicit none
      integer,intent(in) :: nheaps,ilen
      integer,dimension(0:ilen-1),intent(inout) :: heap,invmap
      real(kind=8),dimension(0:ilen-1),intent(in)    :: wght
      integer,intent(in) :: popid

      integer i,j,k,t
      
      i = popid
      j = ishft(i-1,-1);
      do while ( (i.gt.0).and.(wght(heap(j)).lt.wght(heap(i))))
        t = heap(i)
        heap(i) = heap(j)
        heap(j) = t
        invmap(heap(i)) = i
        invmap(heap(j)) = j
        i = j
        j = ishft(i-1,-1)
      enddo

      i = popid
      do while (i.lt.(nheaps/2))
         j = 2*i+1
      if ((j.lt.(nheaps-1)).and.(wght(heap(j)).lt.wght(heap(j+1)))) 
     &   then
         j = j + 1
      end if
      if (wght(heap(i)).gt.wght(heap(j))) then 
          exit
      end if
      t = heap(i)
      heap(i) = heap(j)
      heap(j) = t
      invmap(heap(i)) = i
      invmap(heap(j)) = j
      i = j
      enddo

      end subroutine !ramg_adjheap

!****************************************************
!     ramg_CFsplit: split coarse/fine nodes
!     Input: matrix in sparse form
!            parameters: (filter array)(aggressive/standard)
!     Output: C/F splitting result, Strong corelation set amg_S
!            within given filter subset
!****************************************************

      subroutine ramg_CFsplit(amg_CF,amg_S,anshg,annz,acolm,arowp,
     &                        alhs,afmap,ilwork,BC,iBC,iper,
     &           eps_str,afilter,spflag)
      use ramg_data

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      
      ! parameters 
      integer,intent(in)                         :: anshg,annz
      integer,intent(inout),dimension(anshg)     :: amg_CF
      integer,intent(in),dimension(anshg+1)      :: acolm
      integer,intent(in),dimension(annz)         :: arowp
      integer,intent(inout),dimension(annz)      :: amg_S
      real(kind=8),intent(in),dimension(annz)    :: alhs
      integer,intent(in),dimension(anshg)        :: afmap
      integer,intent(in),dimension(nlwork)     :: ilwork
      integer,intent(in),dimension(nshg)       :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
 
      real(kind=8)                               :: eps_str
      integer,intent(in)                         :: afilter,spflag

      ! afilter: indicating which subset of afmap for coarsening
         !          -1 for self
         ! spflag: aggressive or not

      integer       :: i,j,k,n,m
      integer,dimension(anshg)  :: aheap,ainvheap
      real(kind=8)  :: Lmax
      real(kind=8),dimension(anshg) :: amg_L
      real(kind=8),dimension(anshg) :: rowmax,prowmax
      real(kind=8) :: MSCALE,LSCALE,deltal
      integer       :: flagS,nheaps,cfilter
      ! setup S, S^T matrix in amg_S
      ! after this, S(i,j) = 1,3 : S
      !             S(i,j) = 2,3 : S^T

      integer,allocatable,dimension(:) :: oneck
    
      !amg_S = 0

      allocate(oneck(numpe+1))
      oneck = 0

      LSCALE = 0.005
      MSCALE = 0.1*LSCALE/anshg

      do i=1,anshg
         cfilter = afmap(i) ! each row only deal with its own group
         amg_CF(i) = 0
         amg_L(i) = 0
         rowmax(i) = 1E+10
         prowmax(i) = -1e+10
         do j=acolm(i)+1,acolm(i+1)-1
            k = arowp(j)
          if ((cfilter.eq.afmap(k)).and.(alhs(j).lt.rowmax(i))) then
                rowmax(i) = alhs(j)
          end if
          if ((cfilter.eq.afmap(k)).and.(alhs(j).gt.prowmax(i))) then
                prowmax(i) = alhs(j)
          end if
         enddo
         rowmax(i) = eps_str*rowmax(i)
         prowmax(i) = eps_str*prowmax(i)
       !  if (anshg.lt.1000) then
       !write(*,"((I5)(3E14.5))")i,rowmax(i),alhs(acolm(i)),prowmax(i)
       !  endif
      enddo
      do i=1,anshg
         cfilter = afmap(i)
         p = iabs(cfilter)
         oneck(p) = oneck(p)+1
         do j=acolm(i)+1,acolm(i+1)-1
            k = arowp(j)
            if (cfilter.eq.afmap(k)) then ! same subset
                if (alhs(j).lt.0) then ! negative entries
            if (alhs(j) .lt. rowmax(i)) then
                amg_S(j) = amg_S(j) + 1
            end if
            if (alhs(j) .lt. rowmax(k)) then
                amg_S(j) = amg_S(j) + 2
            end if 
                else ! positive entries
            IF (.FALSE.) THEN
            if (alhs(j) .gt. prowmax(i)) then
                amg_S(j) = amg_S(j) + 1
            end if
            if (alhs(j) .gt. prowmax(k)) then
                amg_S(j) = amg_S(j) + 2
            end if 
            ENDIF
                endif ! if negative or positive
            endif
         enddo
      enddo
      ! mark isolated unknowns "coarse", S^T=0,
      ! keep territory which should not be coarsened "coarse"
      ! keep slave coarse temporarily
      do i=1,anshg
         cfilter = afmap(i)
         flagS = 0
         do j=acolm(i),acolm(i+1)-1
            k = arowp(j)
            if ((amg_S(j).ne.0).and.(afmap(k).eq.cfilter)) then
                flagS = 1
                exit
            end if
         enddo
         if ( flagS.eq.0) then
             p = iabs(cfilter)
             oneck(p) = oneck(p)-1
!             if (oneck(p).gt.0) then
             if (.false.) then
             amg_CF(i) = 2 ! mark fine ! NO, should mark isolated coarse
             else
                 amg_CF(i) = 1 ! last coarse ! everybody coarse
             endif
         end if
      enddo

      k = 0
      do i=1,anshg
         if (amg_CF(i).eq.2) then 
             k = k+1
         endif
      enddo
      if (k.eq.anshg) then
          write(*,*)'mcheck coarsening error here'
          stop
      endif
      ! setup initial lamda value
      Lmax = -1.0e+6 ! Lmax stores max value of lambda, m points to it
      do i=1,anshg
         cfilter = afmap(i)
         if (amg_CF(i).eq.0) then 
         ! unassigned with in the filter
         amg_L(i) = i*MSCALE+iabs(cfilter)
         do j=acolm(i),acolm(i+1)-1
            k = arowp(j)
            if ((afmap(k).eq.cfilter).and.(amg_S(j).ge.2)) then ! in S^T
                deltal = 0
                if (amg_CF(k).eq.0) then ! undecided
                    deltal = LSCALE
                else if (amg_CF(k).eq.2) then ! free
                    deltal = 2*LSCALE
                end if!filter
                if (alhs(j).gt.0) deltal=-deltal
                amg_L(i) = amg_L(i) + deltal

                !if (alhs(j).gt.0) then
                !    write(*,*)'hooh: ',alhs(j)
                !endif
                ! this is to seperate nodes in each section
                ! in the heap. e.g. 0-1000 for neighbour with proc 0, 
                ! 2000-3000 for neighbour with proc 1... 
            end if
         enddo
         if (amg_L(i).gt.Lmax) then
             Lmax = amg_L(i)
             m = i
         end if
         end if
      enddo
      ! aLoc: heap, amg_Fn: invheap
      nheaps = anshg
      call ramg_initheap(aheap,ainvheap,amg_L,nheaps,anshg)
      !write(*,*)oneck
      ! setup coarse/fine grid
      do while (nheaps.gt.0)
         m = aheap(1)+1
         if (amg_CF(m).ne.0) then
             amg_L(m) = 0
             call ramg_popheap(aheap,ainvheap,amg_L,nheaps,
     &            ainvheap(m),anshg)
           cycle
         endif
         Lmax = amg_L(m)-iabs(afmap(m))-m*MSCALE ! reuse Lmax
         if (Lmax.lt.LSCALE) then
             k = iabs(afmap(m))
             if (oneck(k).eq.1) then
                 amg_CF(m) = 1
             else
                 amg_CF(m) = 2
             endif
             oneck(iabs(afmap(m)))=oneck(iabs(afmap(m)))-1
 
             amg_L(m) = 0!afmap(m)+m*MSCALE
             call ramg_popheap(aheap,ainvheap,amg_L,nheaps,
     &            ainvheap(m),anshg)
           cycle
         end if
        ! set coarse grid
         amg_CF(m) = 1 ! mark coarse
         amg_L(m) = 0!afmap(m)+m*MSCALE!0 !floor(amg_L(m))
         call ramg_popheap(aheap,ainvheap,amg_L,nheaps,ainvheap(m),
     &        anshg)
         ! set fine grid
         do j=acolm(m),acolm(m+1)-1
            k = arowp(j)
            if (afmap(k).eq.afmap(m)) then
             if (amg_S(j).ge.2) then ! in m's S^T if !B
                if (amg_CF(k).eq.0) then !not defined yet !if A
                   ! update flag
                   amg_CF(k) = 2  ! mark fine
                  amg_L(k) = 0!afmap(k)+k*MSCALE!0 !floor(amg_L(k))
         call ramg_popheap(aheap,ainvheap,amg_L,nheaps,ainvheap(k),
     &        anshg)
                  ! update lambda in k's S
                   do i = acolm(k),acolm(k+1)-1
                      n = arowp(i)
                    if (((amg_S(i).eq.1).or.(amg_S(i).eq.3)).and.
     &             (amg_CF(n).eq.0)) then
                        ! in k's S , U==>F
                        if (alhs(i).le.0) then ! negative
                        amg_L(n) = amg_L(n) + LSCALE 
                        else ! positive
                        amg_L(n) = amg_L(n) - LSCALE 
                        endif
         call ramg_adjheap(aheap,ainvheap,amg_L,nheaps,ainvheap(n),
     &         anshg)
                      end if
                   enddo
                end if ! if A
             end if ! if B
             if (mod(amg_S(j),2).eq.1) then !if B2
                if (amg_CF(k).eq.0) then
                    if (alhs(j).le.0) then
                    amg_L(k) = amg_L(k) - LSCALE
                    else
                    amg_L(k) = amg_L(k) + LSCALE
                    endif
         call ramg_adjheap(aheap,ainvheap,amg_L,nheaps,ainvheap(k),
     &         anshg)
                end if
             end if ! if B2
            end if ! subset
         enddo
      enddo
      k = 0
      do i=1,anshg
         if (amg_CF(i).eq.2) then 
             k = k+1
         endif
      enddo
      if (k.eq.anshg) then
          write(*,*)'mcheck coarsening error here,2'
          stop
      endif
      !write(*,*)oneck

      ! any U left
      j = 0
      do i = 1,anshg
         if (amg_CF(i).eq.0) then
             amg_CF(i) = 1 ! Mark every node coarse if isolated
!             amg_CF(i) = 2 ! Mark every node fine if isolated
         end if
         if (amg_CF(i).eq.1) then 
             j = j+1
         endif
      enddo
      !write(*,*)'mcheck rank2b, j=',myrank,j
      if (.false.) then
      k = 0
      do i=1,anshg
         if (amg_CF(i).eq.1) then 
             k = k+1
             do j=acolm(i)+1,acolm(i+1)-1
                m = arowp(j)
                if (amg_CF(m).eq.1) then ! two coarse connected
                if (amg_S(j).ge.1) then
                !    write(*,*)'mcheck err',alhs(j)
                endif
                endif
             enddo
         endif
      enddo
      if (k.eq.anshg) then
          write(*,*)'mcheck coarsening error here,3,'
          stop
      endif
      endif  
      !write(*,*)'mcheck end of cfsplit'
      !call ramg_dump_i(amg_CF,anshg,1,'splitcf   ')
      deallocate(oneck)

      end subroutine !ramg_CFsplit

      subroutine ramg_update_cfmap(amg_CF,slevel)
      use ramg_data
      integer,intent(in)                       :: slevel
      integer,intent(in),dimension(amg_nshg(slevel)) :: amg_CF
      integer :: i,j,p
      integer,dimension(amg_nshg(slevel))      :: revmap

      revmap = 0
      j = 1
      do i=1,amg_nshg(1)
         if (amg_cfmap(i).eq.slevel) then
             revmap(j) = i
             j = j+1
         end if
      enddo
      do i=1,j-1
         if (amg_CF(i).eq.1) then
             amg_cfmap(revmap(i)) = amg_cfmap(revmap(i)) + 1
         end if
      enddo
      
      end subroutine ! ramg_update_cfmap

      subroutine ramg_readin_cfmap(amg_CF,slevel)
      use ramg_data
      integer,intent(in)                     :: slevel
      integer,intent(inout),dimension(amg_nshg(slevel)) :: amg_CF
      integer :: i,j,p

      j = 1
      do i=1,amg_nshg(1)
         if (amg_cfmap(i).eq.slevel) then
             amg_CF(j) = 2 ! fine
             j = j + 1
         else if (amg_cfmap(i).eq.(slevel+1)) then
             amg_CF(j) = 1 ! coarse
             j = j + 1
         endif
      enddo

      end subroutine ! ramg_readin_cfmap
