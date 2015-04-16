
!***********************************************************
!      mtx_colm,mtx_rowp,mtx_mtx
!**********************************************************
      subroutine mtx_colm(a_colm,a_rowp,
     &                    b_colm,b_rowp,
     &                    c_colm,m,l,n)
      integer,intent(in)                :: m,l,n
      integer,intent(in) :: a_colm(m+1),b_colm(l+1)
      integer,intent(inout) :: c_colm(m+1)
      integer,intent(in) :: a_rowp(a_colm(m+1)-1)
      integer,intent(in) :: b_rowp(b_colm(l+1)-1)

      integer      :: i,j,k,p,q,mnnz
      integer      :: itemp(max(m,l))
      
      mnnz = 1
      itemp = 0
      do i=1,m
         c_colm(i) = mnnz
         do k=a_colm(i),a_colm(i+1)-1
            j=a_rowp(k)
            do p=b_colm(j),b_colm(j+1)-1
               q=b_rowp(p)
               if (itemp(q).ne.i) then
                   itemp(q)=i
                   mnnz = mnnz + 1
               end if
            enddo
         enddo 
      enddo
      c_colm(m+1) = mnnz
      
      end subroutine!mtx_colm

      subroutine mtx_rowp(a_colm,a_rowp,b_colm,b_rowp,
     &                    c_colm,c_rowp,m,l,n,diagflag)
      integer,intent(in)                :: m,l,n
      integer,intent(in) :: a_colm(m+1),b_colm(l+1)
      integer,intent(in) :: c_colm(m+1)
      integer,intent(inout) :: c_rowp(c_colm(m+1)-1)
      integer,intent(in) :: a_rowp(a_colm(m+1)-1)
      integer,intent(in) :: b_rowp(b_colm(l+1)-1)
      logical,intent(in) :: diagflag

      integer      :: i,j,k,p,q,mnnz
      integer      :: itemp(max(m,l))
      
      mnnz = 0
      itemp = 0
      do i=1,m
         if (diagflag) then
         mnnz = mnnz + 1
         c_rowp(mnnz) = i
         itemp(i)=i
         endif
         do k=a_colm(i),a_colm(i+1)-1
            j=a_rowp(k)
            do p=b_colm(j),b_colm(j+1)-1
               q=b_rowp(p)
               if (itemp(q).ne.i) then
                   itemp(q)=i
                   mnnz = mnnz + 1
                   c_rowp(mnnz) = q
               end if
            enddo
         enddo 
      enddo
      
      end subroutine !mtx_rowp 
      
      subroutine mtx_mtx(a_colm,a_rowp,a_mtx,
     &                   b_colm,b_rowp,b_mtx,
     &                   c_colm,c_rowp,c_mtx,
     &                   m,l,n)
      integer,intent(in)                :: m,l,n
      integer,intent(in) :: a_colm(m+1),b_colm(l+1)
      integer,intent(in) :: c_colm(m+1)
      integer,intent(in) :: c_rowp(c_colm(m+1)-1)
      integer,intent(in) :: a_rowp(a_colm(m+1)-1)
      integer,intent(in) :: b_rowp(b_colm(l+1)-1)

      real(kind=8),intent(in) :: a_mtx(a_colm(m+1)-1)
      real(kind=8),intent(in) :: b_mtx(b_colm(l+1)-1)
      real(kind=8),intent(inout) :: c_mtx(c_colm(m+1)-1)

      integer      :: i,j,k,p,q
      integer      :: itemp(max(m,l))
      real(kind=8) :: rtemp

      c_mtx = 0
      do i=1,m
         do k=c_colm(i),c_colm(i+1)-1
            itemp(c_rowp(k)) = k
         enddo
         do k=a_colm(i),a_colm(i+1)-1
            j = a_rowp(k)
            rtemp = a_mtx(k)
            do p = b_colm(j),b_colm(j+1)-1
               q = b_rowp(p)
               c_mtx(itemp(q)) = c_mtx(itemp(q)) + rtemp*b_mtx(p)
            enddo
         enddo
      enddo
      
      end subroutine !mtx_mtx
 
!************************************************************
!     ramg_calcIAI
!        do matrix multiplication for I^T A I
!        actual operation is: A2_ij = I^T_ip A_pq I_qj
!        where p,q are repeated indicies
!      
!************************************************************
      subroutine ramg_calcITAI(level1,level2,maxstopsign)
      use ramg_data
      implicit none
      
      integer,intent(in)             :: level1,level2
      logical,intent(inout)          :: maxstopsign
      
      integer                        :: i,j,k,m,p,q,r,n,s,p0
      integer,dimension(amg_nshg(level1))  :: itemp
      integer,dimension(amg_nshg(level2))  :: itemp2
      integer                              :: mnnz

      integer,dimension(amg_nshg(level2)+1)   :: ITA_colm
      integer,dimension(:),allocatable        :: ITA_rowp
      real(kind=8),dimension(:),allocatable   :: ITA

      integer                              :: mem_err, mem_err_s
      
      real(kind=8)                             :: rtp
      character                            :: fname*80
      real,dimension(10)                   :: cpusec

      call cpu_time(cpusec(1))

      IF (level1.eq.0) THEN
          ! reduced serial test here, DO NOT USE IT unless RECOVER
          ! from extract.
      call mtx_colm(I_fc_colm(level1)%p,I_fc_rowp(level1)%p,
!     &              amg_A_colm(level1)%p,amg_A_rowp(level1)%p,
     &              reducecolm,reducerowp,
     &              ITA_colm,amg_nshg(level2),amg_nshg(level1),
     &              amg_nshg(level1))
      mnnz = ITA_colm(amg_nshg(level2)+1)-1
      allocate(ITA_rowp(mnnz),stat=mem_err)
      allocate(ITA(mnnz),stat=mem_err)
      ! IT * A
      ! rowp
      call mtx_rowp(I_fc_colm(level1)%p,I_fc_rowp(level1)%p,
!     &              amg_A_colm(level1)%p,amg_A_rowp(level1)%p,
     &              reducecolm,reducerowp,
     &              ITA_colm,ITA_rowp,
     &              amg_nshg(level2),amg_nshg(level1),
     &              amg_nshg(level1),.false.)
      ! IT * A
      ! value
      call mtx_mtx(I_fc_colm(level1)%p,I_fc_rowp(level1)%p,
     &              I_fc(level1)%p,
!     &              amg_A_colm(level1)%p,amg_A_rowp(level1)%p,
!     &              amg_A_lhs(level1)%p,
     &              reducecolm,reducerowp,reducelhs,
     &              ITA_colm,ITA_rowp,ITA,
     &              amg_nshg(level2),amg_nshg(level1),
     &              amg_nshg(level1))
      call cpu_time(cpusec(2))
      
      ELSE
          
      ! IT * A
      ! colm and nnz
      call mtx_colm(I_fc_colm(level1)%p,I_fc_rowp(level1)%p,
     &              amg_A_colm(level1)%p,amg_A_rowp(level1)%p,
     &              ITA_colm,amg_nshg(level2),amg_nshg(level1),
     &              amg_nshg(level1))
      mnnz = ITA_colm(amg_nshg(level2)+1)-1
      allocate(ITA_rowp(mnnz),stat=mem_err)
      allocate(ITA(mnnz),stat=mem_err)
      ! IT * A
      ! rowp
      call mtx_rowp(I_fc_colm(level1)%p,I_fc_rowp(level1)%p,
     &              amg_A_colm(level1)%p,amg_A_rowp(level1)%p,
     &              ITA_colm,ITA_rowp,
     &              amg_nshg(level2),amg_nshg(level1),
     &              amg_nshg(level1),.false.)
      ! IT * A
      ! value
      call mtx_mtx(I_fc_colm(level1)%p,I_fc_rowp(level1)%p,
     &              I_fc(level1)%p,
     &              amg_A_colm(level1)%p,amg_A_rowp(level1)%p,
     &              amg_A_lhs(level1)%p,
     &              ITA_colm,ITA_rowp,ITA,
     &              amg_nshg(level2),amg_nshg(level1),
     &              amg_nshg(level1))
      call cpu_time(cpusec(2))
      
      ENDIF

      ! ( IT*A ) * I
      ! colm and nnz
      call ramg_allocate(level2,amg_nshg(level2),0,1)
      call mtx_colm(ITA_colm,ITA_rowp,I_cf_colm(level1)%p,
     &              I_cf_rowp(level1)%p,amg_A_colm(level2)%p,
     &              amg_nshg(level2),amg_nshg(level1),
     &              amg_nshg(level2))
      
      call cpu_time(cpusec(3))
      mnnz = amg_A_colm(level2)%p(amg_nshg(level2)+1)-1
      call ramg_allocate(level2,0,mnnz,1)
      ! ( IT*A ) * I
      ! rowp
      call mtx_rowp(ITA_colm,ITA_rowp,I_cf_colm(level1)%p,
     &              I_cf_rowp(level1)%p,amg_A_colm(level2)%p,
     &              amg_A_rowp(level2)%p,
     &              amg_nshg(level2),amg_nshg(level1),
     &              amg_nshg(level2),.true.)
      call cpu_time(cpusec(4))
      ! ( IT*A ) * I
      ! value
      call mtx_mtx(ITA_colm,ITA_rowp,ITA,I_cf_colm(level1)%p,
     &              I_cf_rowp(level1)%p,I_cf(level1)%p,
     &              amg_A_colm(level2)%p,
     &              amg_A_rowp(level2)%p,amg_A_lhs(level2)%p,
     &              amg_nshg(level2),amg_nshg(level1),
     &              amg_nshg(level2))
      
      deallocate(ITA_rowp,stat=mem_err)
      deallocate(ITA,stat=mem_err)

      call cpu_time(cpusec(5))
     
      amg_nshg(level2+1) = 0

      ! put on right hand side
      call ramg_calcIvFC(level1,level2,amg_A_rhs(level1)%p,
     &               amg_A_rhs(level2)%p)

      call cpu_time(cpusec(10))

      do i=1,amg_nshg(level2)
         amg_ppeDiag(level2)%p(i) = 
     &   1.0/sqrt(amg_A_lhs(level2)%p(amg_A_colm(level2)%p(i),1))
      enddo
      do i=1,amg_nshg(level2)
         do m=amg_A_colm(level2)%p(i),amg_A_colm(level2)%p(i+1)-1
            j = amg_A_rowp(level2)%p(m)
            amg_A_lhs(level2)%p(m,1)=amg_A_lhs(level2)%p(m,1)*
     &        amg_ppeDiag(level2)%p(i)*amg_ppeDiag(level2)%p(j)
         enddo
      enddo

!      do i=1,amg_nshg(level1)
!         do j=I_cf_colm(level1)%p(i),I_cf_colm(level1)%p(i+1)-1
!            k = I_cf_rowp(level1)%p(j)
!            I_cf(level1)%p(j) = 
!     &      I_cf(level1)%p(j)*amg_ppeDiag(level2)%p(k)
!         enddo
!      enddo

!      do i=1,amg_nshg(level2)
!         do j=I_fc_colm(level1)%p(i),I_fc_colm(level1)%p(i+1)-1
!         I_cf(level1)%p(j)=
!     &   I_cf(level1)%p(j)*amg_ppeDiag(level2)%p(i)
!         enddo
!      enddo

      end subroutine ! ramg_calcIAI


