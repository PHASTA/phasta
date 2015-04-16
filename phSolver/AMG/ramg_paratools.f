
!**************************************************************
!      ramg_global_PPE
!       Make lhsP global complete 
!**************************************************************

      subroutine ramg_global_lhs(
     &               acolm,arowp,alhsP,annz_tot,
     &               lhsGPcolm,lhsGProwp,lhsGP,
     &               redun,
     &               ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in)            :: annz_tot
      integer,intent(in),dimension(nshg+1)            :: acolm
      integer,intent(in),dimension(annz_tot)           :: arowp
      integer,intent(in)                              :: redun
      real(kind=8),intent(in),dimension(redun,annz_tot)    :: alhsP
!      real(kind=8),dimension(:,:),allocatable          :: lhsGP
!      integer,dimension(:),allocatable                 :: lhsGProwp
!      integer,dimension(:),allocatable                 :: lhsGPcolm
      type(r2d),intent(inout) :: lhsGP
      type(i1d),intent(inout) :: lhsGProwp,lhsGPcolm
      integer,intent(in),dimension(nlwork)            :: ilwork
      integer,intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)  :: BC

      integer           :: numtask,i,j,k,m,p,ki,kj,krowp,ierr
      integer           :: gi,gj,gk
      integer           :: itag,iacc,iother,numseg,isgbeg,itkbeg
      integer    :: stat(MPI_STATUS_SIZE,redun*numpe),req(redun*numpe)
      integer           :: mcheck

      real(kind=8) :: swaptemp(redun)

      character fname*10
      
      ! temp variables, used for communication

      ! each task has a set of storage rooms for the submatrix
      ! to communicate

      integer,dimension(nshg) :: tmp_rowmap,tmp_revrmap
      real(kind=8),dimension(redun,nshg) :: tmp_rmtx

      integer,dimension(nshg) :: Pflag
      integer,dimension(nshg+1) :: Pcolm
      type(i2dd) :: Prowp
      type(r12dd) :: Pmtx
      integer :: rownnz
      

      if (numpe .le. 1) then
          allocate(lhsGP%p(redun,annz_tot),stat=ierr)
          allocate(lhsGProwp%p(annz_tot))
          allocate(lhsGPcolm%p(nshg+1))
          lhsGP%p(:,:) = alhsP(:,:)
          lhsGProwp%p(:) = arowp(:)
          lhsGPcolm%p(:) = acolm(:)
          return
      end if

      numtask = ilwork(1)
      m = 0
      itkbeg = 1

      revmap = 0

      allocate(sub_map%pp(numtask))
      allocate(sub_revmap%pp(numtask))
      allocate(sub_colm%pp(numtask))
      allocate(sub_colm2%pp(numtask))
      allocate(sub_rowp%pp(numtask))
      allocate(sub_rowp2%pp(numtask))
      allocate(sub_rowpmap%pp(numtask))
      allocate(sub_mtx%pp(numtask))
      allocate(sub_mtx2%pp(numtask))
      allocate(sub_nnz%p(numtask))
      allocate(sub_nnz2%p(numtask))
      allocate(sub_nshg%p(numtask))

      ! submatrix : colm
      do i=1,numtask
      
         ! Beginning of each task
         m = m+1
         ! Task head information, ref commu.f
         itag = ilwork(itkbeg+1)
         iacc = ilwork(itkbeg+2)
         iother = ilwork(itkbeg+3)
         numseg = ilwork(itkbeg+4)
         isgbeg = ilwork(itkbeg+5)

         sub_nshg%p(i) = numseg

         allocate(sub_map%pp(i)%p(numseg))
         allocate(sub_revmap%pp(i)%p(nshg))

         sub_map%pp(i)%p = 0
         sub_revmap%pp(i)%p = 0
         
         ! prepare vector mapping
         do j=1,numseg
            k = ilwork(itkbeg+3+2*j) ! row idx
            sub_map%pp(i)%p(j) = k
            sub_revmap%pp(i)%p(k) = j
!            if ((myrank.eq.1).and.(iother.eq.3)) then
!            write(*,*)'mcheck:',myrank,iother,j,k,ncorp_map(k)
!            endif
         enddo

         allocate(sub_colm%pp(i)%p(numseg+1))
         allocate(sub_colm2%pp(i)%p(numseg+1))

         ! prepare matrix mapping, sub matrix colm
         sub_nnz%p(i) = 0
         do j=1,numseg
            sub_colm%pp(i)%p(j) = sub_nnz%p(i) + 1
            ki = sub_map%pp(i)%p(j)
            do kj = acolm(ki),acolm(ki+1)-1
               krowp = arowp(kj)
               if (sub_revmap%pp(i)%p(krowp).ne.0) then
!      if ((ncorp_map(ki).eq.964)) then
!      if ((ncorp_map(krowp).eq.884)) then
!          write(*,*)'964/884 prepared in',myrank,i,iother
!      else
!          write(*,*)'line 964:',myrank,i,iother,ncorp_map(krowp)
!      endif
!      endif
                   sub_nnz%p(i) = sub_nnz%p(i) + 1
               end if
            enddo
            !write(*,*)'mcheck: ',myrank,j,sub_nnz
         enddo
         sub_colm%pp(i)%p(numseg+1) = sub_nnz%p(i) + 1

         if (iacc.eq.0) then ! this task is a send, master
            call MPI_ISEND(sub_colm%pp(i)%p(1),numseg+1,MPI_INTEGER,
     &            iother,itag,MPI_COMM_WORLD,req(m),ierr)
            call MPI_IRECV(sub_colm2%pp(i)%p(1),numseg+1,MPI_INTEGER,
     &            iother,itag,MPI_COMM_WORLD,req(m+1),ierr)
 
         else if (iacc.eq.1) then ! this task is a receive, slave
            call MPI_IRECV(sub_colm2%pp(i)%p(1),numseg+1,MPI_INTEGER,
     &            iother,itag,MPI_COMM_WORLD,req(m),ierr)
            call MPI_ISEND(sub_colm%pp(i)%p(1),numseg+1,MPI_INTEGER,
     &            iother,itag,MPI_COMM_WORLD,req(m+1),ierr)

         end if

         m = m + 1

         ! Task end, next task
         itkbeg = itkbeg + 4 + 2*numseg
         
      enddo

      call MPI_WAITALL(m,req,stat,ierr)

      do i=1,numtask
         sub_nnz2%p(i) = sub_colm2%pp(i)%p(sub_nshg%p(i)+1)-1
      enddo

      ! sub matrix : rowp,mtx
      m = 0
      itkbeg = 1
      do i=1,numtask
      
         ! Beginning of each task
         m = m+1
         ! Task head information, ref commu.f
         itag = ilwork(itkbeg+1)
         iacc = ilwork(itkbeg+2)
         iother = ilwork(itkbeg+3)
         numseg = ilwork(itkbeg+4)
         isgbeg = ilwork(itkbeg+5)

         allocate(sub_rowp%pp(i)%p(sub_nnz%p(i)))
         allocate(sub_rowp2%pp(i)%p(sub_nnz2%p(i)))
         allocate(sub_rowpmap%pp(i)%p(sub_nnz%p(i)))
         allocate(sub_mtx%pp(i)%p(redun,sub_nnz%p(i)))
         allocate(sub_mtx2%pp(i)%p(redun,sub_nnz2%p(i)))

         ! prepare matrix mapping, sub matrix rowp
         k = 0
         do j=1,numseg
            ki = sub_map%pp(i)%p(j)
            do kj = acolm(ki),acolm(ki+1)-1
               krowp = arowp(kj)
               if (sub_revmap%pp(i)%p(krowp).ne.0) then
                   k = k + 1
                   sub_rowp%pp(i)%p(k) = sub_revmap%pp(i)%p(krowp)
                   sub_rowpmap%pp(i)%p(k) = kj
                   do p=1,redun
                      sub_mtx%pp(i)%p(p,k) = alhsP(p,kj)
                   enddo
               end if
            enddo
         enddo

         if (iacc.eq.0) then ! this task is a send, master
        call MPI_ISEND(sub_rowp%pp(i)%p(1),sub_nnz%p(i),MPI_INTEGER,
     &            iother,itag,MPI_COMM_WORLD,req(m),ierr)
       call MPI_IRECV(sub_rowp2%pp(i)%p(1),sub_nnz2%p(i),MPI_INTEGER,
     &            iother,itag,MPI_COMM_WORLD,req(m+1),ierr)
 
         else if (iacc.eq.1) then ! this task is a receive, slave
       call MPI_IRECV(sub_rowp2%pp(i)%p(1),sub_nnz2%p(i),MPI_INTEGER,
     &            iother,itag,MPI_COMM_WORLD,req(m),ierr)
        call MPI_ISEND(sub_rowp%pp(i)%p(1),sub_nnz%p(i),MPI_INTEGER,
     &            iother,itag,MPI_COMM_WORLD,req(m+1),ierr)
         end if
         m = m + 2
         if (iacc.eq.0) then ! this task is a send, master
        call MPI_ISEND(sub_mtx%pp(i)%p(1,1),redun*sub_nnz%p(i),
     &            MPI_DOUBLE_PRECISION,
     &            iother,itag,MPI_COMM_WORLD,req(m),ierr)
       call MPI_IRECV(sub_mtx2%pp(i)%p(1,1),redun*sub_nnz2%p(i),
     &            MPI_DOUBLE_PRECISION,
     &            iother,itag,MPI_COMM_WORLD,req(m+1),ierr)
 
         else if (iacc.eq.1) then ! this task is a receive, slave
        call MPI_IRECV(sub_mtx2%pp(i)%p(1,1),redun*sub_nnz2%p(i),
     &            MPI_DOUBLE_PRECISION,
     &            iother,itag,MPI_COMM_WORLD,req(m),ierr)
        call MPI_ISEND(sub_mtx%pp(i)%p(1,1),redun*sub_nnz%p(i),
     &            MPI_DOUBLE_PRECISION,
     &            iother,itag,MPI_COMM_WORLD,req(m+1),ierr)
         end if

         m = m + 1

         ! Task end, next task
         itkbeg = itkbeg + 4 + 2*numseg
         
      enddo

      call MPI_WAITALL(m,req,stat,ierr)

      if (.false.) then
      ! dump some matrices for matlab
      do i=1,numtask
       write(fname,'((A4)(I1)(A5))')'subP',i,'     '
       !call ramg_dump_i(ncorp_map,nshg,1,'pam_corpmd')
!       call ramg_dump_matlab_map(sub_colm%pp(i)%p,sub_rowp%pp(i)%p,
!     &      sub_mtx%pp(i)%p,sub_nshg%p(i),sub_nnz%p(i),4,
!     &      sub_map%pp(i)%p,fname)
!      enddo
!      do i=1,numtask
!       write(fname,'((A5)(I1)(A4))')'subPs',i,'    '
!       !call ramg_dump_i(ncorp_map,nshg,1,'pam_corpmd')
!       call ramg_dump_matlab_map(sub_colm2%pp(i)%p,sub_rowp2%pp(i)%p,
!     &      sub_mtx2%pp(i)%p,sub_nshg%p(i),sub_nnz2%p(i),4,
!     &      sub_map%pp(i)%p,fname)
      enddo
      endif ! dump? no dump?

!      call ramg_dump_matlab_map(acolm,arowp,alhsP,nshg,annz_tot,redun,
!     &     'locallhsPb')

      ! merge with extra entries

      Pflag = 0
      allocate(Prowp%pp(nshg))
      allocate(Pmtx%pp(nshg))
      Pcolm = 0

      do i=1,numtask ! each task
         do j=1,sub_nshg%p(i) ! each line
            
            tmp_rowmap = 0
            tmp_revrmap = 0
            tmp_rmtx = 0
            rownnz = 0
            
            ! Prepare
            ki = sub_map%pp(i)%p(j)
            if (Pflag(ki).eq.0) then ! a new line
                Pflag(ki) = 1
                do k=acolm(ki),acolm(ki+1)-1
                   kj = arowp(k)
                   tmp_rowmap(kj) = k-acolm(ki)+1
                   tmp_revrmap(k-acolm(ki)+1)=kj
                   tmp_rmtx(:,kj) = alhsP(:,k)
                enddo
                rownnz = rownnz+acolm(ki+1)-acolm(ki)
            else ! existing line
                ! expand sparse line to full line
                do k=1,Pcolm(ki)
                   kj = Prowp%pp(ki)%p(k)
                   tmp_rowmap(kj) = k
                   tmp_revrmap(k) = kj
                   tmp_rmtx(:,kj) = Pmtx%pp(ki)%p(:,k)
                enddo
                rownnz = rownnz+Pcolm(ki)
                deallocate(Prowp%pp(ki)%p)
                deallocate(Pmtx%pp(ki)%p)
            endif

            ! Merge
            do ki = sub_colm2%pp(i)%p(j),sub_colm2%pp(i)%p(j+1)-1
            ! each entry in second mtx /slave
               krowp = sub_rowp2%pp(i)%p(ki)
               kj = sub_map%pp(i)%p(krowp)
               if (tmp_rowmap(kj).eq.0) then ! not yet occupied
               rownnz = rownnz+1
               tmp_rowmap(kj) = rownnz
               tmp_revrmap(rownnz) = kj
               endif
!         if ((ncorp_map(sub_map%pp(i)%p(j)).eq.964)) then
!         if ((ncorp_map(kj).eq.883)) then
!             write(*,*)'paramcheck',myrank,i,sub_mtx2%pp(i)%p(1,ki)
!         else
!             write(*,*)'paramcheck2',myrank,ncorp_map(kj)
!         endif
!         endif
              do p=1,redun
               tmp_rmtx(p,kj) = tmp_rmtx(p,kj)+sub_mtx2%pp(i)%p(p,ki)
               enddo
            enddo

            ! Store
            ki = sub_map%pp(i)%p(j)
            Pcolm(ki) = rownnz
            allocate(Prowp%pp(ki)%p(rownnz))
            allocate(Pmtx%pp(ki)%p(redun,rownnz))
            do k=1,rownnz
               kj = tmp_revrmap(k)
               Prowp%pp(ki)%p(k) = kj
               Pmtx%pp(ki)%p(:,k) = tmp_rmtx(:,kj)
!               if ((ncorp_map(ki).eq.964)) then
!                   write(*,*)'paramcheck',myrank,i,k,kj,ncorp_map(kj)
!               endif
!         if ((ncorp_map(ki).eq.964).and.(ncorp_map(kj).eq.883)) then
!             write(*,*)'paramcheck',myrank,i,tmp_rmtx(1,kj)
!         endif
            enddo

            !sort
            do gi=1,rownnz
            do gj=gi+1,rownnz
               if (Prowp%pp(ki)%p(gi).gt.Prowp%pp(ki)%p(gj)) then
                   gk = Prowp%pp(ki)%p(gj)
                   Prowp%pp(ki)%p(gj) = Prowp%pp(ki)%p(gi)
                   Prowp%pp(ki)%p(gi) = gk
                   swaptemp(:) = Pmtx%pp(ki)%p(:,gj)
                   Pmtx%pp(ki)%p(:,gj) = Pmtx%pp(ki)%p(:,gi)
                   Pmtx%pp(ki)%p(:,gi) = swaptemp(:)
               endif
            enddo
            enddo

         enddo
      enddo

      allocate(lhsGPcolm%p(nshg+1))
      
      rownnz = 0
      k = 0
      lhsGPcolm%p(1) = 1
      do i=1,nshg
      ! colm 
         if (Pflag(i).eq.0) then ! original lhsP
             rownnz = rownnz+acolm(i+1)-acolm(i)
         else ! new lhsP
             rownnz = rownnz+Pcolm(i)
             k = k+1
         endif
         lhsGPcolm%p(i+1)=rownnz+1
      enddo

!      call ramg_dump_i(lhsGPcolm%p,nshg,1,'lhsGPcolm ')
      
      allocate(lhsGProwp%p(rownnz))
      allocate(lhsGP%p(redun,rownnz))

      rownnz = 1
      do i=1,nshg
      ! rowp & mtx.
         if (Pflag(i).eq.0) then ! original lhsP
             do j=acolm(i),acolm(i+1)-1
                lhsGProwp%p(rownnz) = arowp(j)
                lhsGP%p(:,rownnz) = alhsP(:,j)
                rownnz = rownnz + 1
             enddo
         else ! new lhsP
             do j=1,Pcolm(i)
                lhsGProwp%p(rownnz) = Prowp%pp(i)%p(j)
                lhsGP%p(:,rownnz) = Pmtx%pp(i)%p(:,j)
                rownnz = rownnz + 1
             enddo
         endif
      enddo

      ! First entry be diagonal for PPE
      if (redun.eq.1) then
      loop_i: do i=1,nshg
         gi = lhsGPcolm%p(i)
         if (lhsGProwp%p(gi).ne.i) then
         do j=gi+1,lhsGPcolm%p(i+1)-1
            k = lhsGProwp%p(j)
            if (k.eq.i) then !swap first and k(diag)
!                 gj = lhsGProwp%p(gi)
                lhsGProwp%p(j) = lhsGProwp%p(gi)
                lhsGProwp%p(gi) = i
                
                swaptemp(:) = lhsGP%p(:,gi)
                lhsGP%p(:,gi) = lhsGP%p(:,j)
                lhsGP%p(:,j) = swaptemp(:)
                cycle loop_i
            endif
         enddo
         endif
      enddo loop_i
      endif

      if (redun.eq.1) then ! check PPE
      do i=1,nshg
         do j=lhsGPcolm%p(i),lhsGPcolm%p(i+1)-1
            k = lhsGProwp%p(j)
            gi = amg_paramap(1)%p(i)
            gk = amg_paramap(1)%p(k)
            if ((gi.eq.gk).and.(gi.lt.0) ) then
            !lhsGP%p(:,j) = 0
            endif
         enddo
      enddo
      endif

!      call ramg_dump_matlab_map(lhsGPcolm%p,lhsGProwp%p,lhsGP%p,
!     &      nshg,rownnz-1,redun,'locallhsP ')

!      write(*,*)'mcheck',myrank,nnz_tot,rownnz,k

!      write(*,*)'mcheck,',myrank,'okay here'

      ! Deallocate temporary arrays
      do i=1,nshg
         if (Pflag(i) .eq. 1) then
             deallocate(Prowp%pp(i)%p)
             deallocate(Pmtx%pp(i)%p)
         endif
      enddo
!          write(*,*)'mcheck deallocate,',myrank

      deallocate(Prowp%pp)
      deallocate(Pmtx%pp)

        do i=1,numtask

           deallocate(sub_map%pp(i)%p)
           deallocate(sub_revmap%pp(i)%p)
           deallocate(sub_colm%pp(i)%p)
           deallocate(sub_colm2%pp(i)%p)
           deallocate(sub_rowp%pp(i)%p)
           deallocate(sub_rowp2%pp(i)%p)
           deallocate(sub_rowpmap%pp(i)%p)
           deallocate(sub_mtx%pp(i)%p)
           deallocate(sub_mtx2%pp(i)%p)

         enddo

        deallocate(sub_map%pp)
        deallocate(sub_revmap%pp)
        deallocate(sub_rowpmap%pp)
        deallocate(sub_nnz%p)
        deallocate(sub_nnz2%p)
        deallocate(sub_nshg%p)
        deallocate(sub_colm%pp)
        deallocate(sub_colm2%pp)
        deallocate(sub_rowp%pp)
        deallocate(sub_rowp2%pp)
        deallocate(sub_mtx%pp)
        deallocate(sub_mtx2%pp)

      end subroutine ! ramg_global_lhs
      
!*******************************************************************
!      ramg_PPEAp
!      produce parallel A-p product for PPE correctly
!      q = PPE * p without scaling
!*******************************************************************
      subroutine ramg_PPEAp(q,p,
     &                      colm,rowp,lhsK,lhsP,
     &                      ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"


      real(kind=8),intent(in),dimension(nshg)         :: p
      real(kind=8),intent(inout),dimension(nshg)      :: q
      integer,intent(in),dimension(nshg+1)            :: colm
      integer,intent(in),dimension(nnz_tot)           :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP
      integer,intent(in),dimension(nlwork)            :: ilwork
      integer,intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)  :: BC

      real(kind=8),dimension(nshg,1) :: t1,t1a
      real(kind=8),dimension(nshg,3) :: t2,t2a
      real(kind=8),dimension(nshg,4) :: t2b
      real(kind=8),dimension(nshg,4) :: diag

      diag = ramg_flowDiag%p

      ! scaling
      call fMtxVdimVecMult(p,amg_ppeDiag(1)%p,t1a,1,1,1,1,nshg)
      call fMtxVdimVecMult(t1a,diag(:,4),t1,1,1,1,1,nshg)
      !t1(:,1) = p
      ! G
      call commOut(t1,ilwork,1,iper,iBC,BC)
      call fLesSparseApG(colm,rowp,lhsP,t1,t2,nshg,nnz_tot)
      call commIn(t2,ilwork,3,iper,iBC,BC)
      ! K^{-1}
      call fMtxVdimVecMult(t2,diag,t2a,3,4,3,3,nshg)
      call fMtxVdimVecMult(t2a,diag,t2,3,4,3,3,nshg)
         ! note: different from lestools.c (3,4,4,3)
         ! t1 should be good to use by now
      ! G^T ... + C
      t2b(:,1:3) = -t2(:,1:3)
      t2b(:,4) = t1(:,1)
      call commOut(t2b,ilwork,4,iper,iBC,BC)
      call fLesSparseApNGtC(colm,rowp,lhsP,t2b,t1,nshg,nnz_tot)
      call commIn(t1,ilwork,1,iper,iBC,BC)
      !q = t1(:,1)
      ! scaling again
      call fMtxVdimVecMult(t1,diag(:,4),t1a,1,1,1,1,nshg)
      call fMtxVdimVecMult(t1a,amg_ppeDiag(1)%p,q,1,1,1,1,nshg)

      end subroutine


!*******************************************************************
!      ramg_PPEAps 
!      produce parallel A-p product for PPE correctly
!      q = PPE * p, WITH SCALING!
!*******************************************************************
      subroutine ramg_PPEAps(q,p,diag,
     &                      colm,rowp,lhsK,lhsP,
     &                      ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"


      real(kind=8),intent(in),dimension(nshg)         :: p
      real(kind=8),intent(inout),dimension(nshg)      :: q
      real(kind=8),intent(in),dimension(nshg,4)            :: diag
      integer,intent(in),dimension(nshg+1)            :: colm
      integer,intent(in),dimension(nnz_tot)           :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP
      integer,intent(in),dimension(nlwork)            :: ilwork
      integer,intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)  :: BC

      real(kind=8),dimension(nshg,1) :: t1,t1a
      real(kind=8),dimension(nshg,3) :: t2,t2a
      real(kind=8),dimension(nshg,4) :: t2b
      
      ! scaling
      call fMtxVdimVecMult(p,amg_ppeDiag(1)%p,t1a,1,1,1,1,nshg)
      call fMtxVdimVecMult(t1a,diag(:,4),t1,1,1,1,1,nshg)
      !call ramg_dump(t1,nshg,'ppe_t1_a  ')
      ! G
      call commOut(t1,ilwork,1,iper,iBC,BC)
      call fLesSparseApG(colm,rowp,lhsP,t1,t2,nshg,nnz_tot)
      call commIn(t2,ilwork,3,iper,iBC,BC)
      ! K^{-1}
      call fMtxVdimVecMult(t2,diag,t2a,3,4,3,3,nshg)
      call fMtxVdimVecMult(t2a,diag,t2,3,4,3,3,nshg)
         ! note: different from lestools.c (3,4,4,3)
         ! t1 should be good to use by now
      ! G^T ... + C
      t2b(:,1:3) = -t2(:,1:3)
      t2b(:,4) = t1(:,1)
      call commOut(t2b,ilwork,4,iper,iBC,BC)
      call fLesSparseApNGtC(colm,rowp,lhsP,t2b,t1,nshg,nnz_tot)
      call commIn(t1,ilwork,1,iper,iBC,BC)
      ! scaling again
      call fMtxVdimVecMult(t1,diag(:,4),t1a,1,1,1,1,nshg)
      call fMtxVdimVecMult(t1a,amg_ppeDiag(1)%p,q,1,1,1,1,nshg)
 
      end subroutine

!*******************************************************************
!      ramg_PPErhs
!      produce a globally correct rhs
!*******************************************************************
      subroutine ramg_PPErhs(rhsp,rhsg,diag,
     &                      colm,rowp,lhsK,lhsP,
     &                      ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"


      real(kind=8),intent(inout),dimension(nshg)         :: rhsp
      real(kind=8),intent(in),dimension(nshg,4)      :: rhsg
      real(kind=8),intent(in),dimension(nshg,4)            :: diag
      integer,intent(in),dimension(nshg+1)            :: colm
      integer,intent(in),dimension(nnz_tot)           :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)    :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)    :: lhsP
      integer,intent(in),dimension(nlwork)            :: ilwork
      integer,intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)  :: BC

      real(kind=8),dimension(nshg,1) :: t1
      real(kind=8),dimension(nshg,3) :: t2,t2a
      
      ! scaling
      ! call fMtxVdimVecMult(rhsg,diag(:,4),t1,1,4,1,1,nshg)
      ! call ramg_dump(t1,nshg,'ppe_t1_a  ')
      ! G
      t2 = rhsg(:,1:3)
      t1 = rhsg(:,4:4)
      call commIn(t1,ilwork,1,iper,iBC,BC)
      call commOut(t1,ilwork,1,iper,iBC,BC)
      call commIn(t2,ilwork,3,iper,iBC,BC)
      call commOut(t2,ilwork,3,iper,iBC,BC)
      ! K^{-1}
      call fMtxVdimVecMult(t2,diag,t2a,3,4,3,3,nshg)
      call fMtxVdimVecMult(t2a,diag,t2,3,4,3,3,nshg)
      call commOut(t2,ilwork,3,iper,iBC,BC)
      call fLesSparseApNGt(colm,rowp,lhsP,t2,rhsp,nshg,nnz_tot)
      call commIn(rhsp,ilwork,1,iper,iBC,BC)
      call commOut(rhsp,ilwork,1,iper,iBC,BC)
      call ramg_zeroOut(rhsp,ilwork,BC,iBC,iper)
      ! -G^T ... + Rc
      rhsp = rhsp - t1(:,1)
      ! scaling again
      call fMtxVdimVecMult(rhsp,diag(:,4),t1,1,1,1,1,nshg)
      rhsp = t1(:,1)
 
      end subroutine ! ramg_PPErhs

!********************************************************
!     ramg_init_ilwork(ilwork,BC,iBC,iper)
!     Initialize amg_ilwork(level)%p(:) array      
!     For each level, same structure as ilwork
!           
!********************************************************      
      subroutine ramg_init_ilwork(ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      integer, intent(in), dimension(nlwork)           :: ilwork
      integer, intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)   :: BC

      integer :: lvl,i,j,iacc,iother,numseg,numtask,itkbeg,sigi,k
      integer :: tasknnztot,itask,amgseg
      integer,dimension(numpe) :: revp
      integer,dimension(ilwork(1)) :: tasknnz,iotherl
      type(i1d),dimension(ilwork(1)) :: taskfill

      integer,dimension(nshg) :: lvlmap

      character :: fname*80

!      call ramg_dump_i(amg_cfmap,nshg,1,'amgcfmap  ')
      !write(*,*)'mcheck ilwork,',myrank,amg_nshg

      if (numpe.le.1) return

      if (ramg_setup_flag.ne.0) return

      numtask = ilwork(1)

      do lvl = 1,ramg_levelx
!         lvl = 1 ! test for 1 level only
      
         ! Create mapping from base to coarse lvl
         ! lvlmap(baselevelindex) = coarselevelindex
         lvlmap=0
         k = 0
         do i=1,nshg
            if (amg_cfmap(i).ge.lvl) then
                k = k+1
                lvlmap(i) = k
            endif
         enddo
         
         ! Count nnz for each task
         itkbeg=1
         tasknnz = 0
         do i=1,numtask
            iacc=ilwork(itkbeg+2)
            iother=ilwork(itkbeg+3) ! starts from 0
            numseg=ilwork(itkbeg+4)
            do j=1,numseg
               k=ilwork(itkbeg+3+2*j) ! row index
               if (amg_cfmap(k).ge.lvl) then
                   tasknnz(i) = tasknnz(i) + 1
               endif
            enddo
            itkbeg=itkbeg+4+2*numseg
         enddo
         tasknnztot = sum(tasknnz)
         !write(*,*)'mcheck ilwork',myrank,lvl,tasknnz

         ! Fill in ilwork array at each lvl
         
         amg_nlwork(lvl)=tasknnztot+1+2*numtask
         allocate(amg_ilwork(lvl)%p(tasknnztot+1+2*numtask))
         amg_ilwork(lvl)%p = 0

         amg_ilwork(lvl)%p(1) = numtask
         itkbeg = 1
         kk = 1 ! pointer to amg_ilwork array
         do i=1,numtask
            iacc = ilwork(itkbeg+2)
            iother = ilwork(itkbeg+3)+1
            if (iacc.eq.0) iother=-iother
            kk = kk+1
            amg_ilwork(lvl)%p(kk) = iother  ! first put iother
            kk = kk+1
            amg_ilwork(lvl)%p(kk) = tasknnz(i) ! then numseg
            numseg = ilwork(itkbeg+4)
            do j=1,numseg
               k = ilwork(itkbeg+3+2*j)
               if (amg_cfmap(k).ge.lvl) then
                   kk = kk+1
                   amg_ilwork(lvl)%p(kk)=lvlmap(k)
               endif
            enddo
            itkbeg=itkbeg+4+2*numseg
         enddo
         !write(fname,'((A9)(I1))')'amgilwork',lvl
         !call ramg_dump_i(amg_ilwork(lvl)%p,amg_nlwork(lvl),1,fname)
      enddo
     
      end subroutine ! ramg_init_ilwork


!*******************************************************
!      ramg_initBCflag: Setup amg_paramap on PPE level(level0)
!      -1: self, n: neighbour proc number
!*******************************************************
      subroutine ramg_initBCflag(flag,ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(inout),dimension(nshg)            :: flag

      integer, intent(in), dimension(nlwork)           :: ilwork
      integer, intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)   :: BC

      integer   :: numtask,itag,iacc,iother,numseg,isgbeg,itkbeg
      integer   :: i,j,k,lenseg,nr,nj,mi,mj
      integer,allocatable,dimension(:)    :: neimap
      character*5 cname
      character*255 fname,fname1

      flag = myrank+1!note*

      if (numpe.lt.2) then
          IF (iamg_reduce.gt.1) then 
          ! this block of code create reduced case
              nr = iamg_reduce
              allocate(reducemap(nr,nr))
              allocate(rmap1d(2,nr*nr)) ! rmap1d(master,slave,rankid)
              reducemap = 0
              rmap1d = 0
              rmapmax = 1
              do i=1,nr
                 reducemap(i,i) = i
                 rmap1d(1,i) = i
                 rmap1d(2,i) = i
              enddo
              rmapmax = nr+1
              write(fname,"((I1)(A))")nr,'-procs_case/neimap.dat'
              !fname=trim(cname(nr))//'-proc_case/map_nproc.dat'
              fname=trim(fname)
              do i=1,nr
                 fname1 = trim(fname)//cname(i)
                 write(*,*)'mcheck reduced case:',fname1
                 open(264,file=trim(fname1),status='unknown')
                 read(264,*)nj
                 do j=1,nj
                    read(264,*)mi,mj
                    if (mj.gt.0) then
                        reducemap(i,mj) = rmapmax
                        reducemap(mj,i) = -rmapmax
                        write(*,*)'mcheck reducemap,',i,mj,rmapmax
                        rmap1d(1,rmapmax) = i
                        rmap1d(2,rmapmax) = mj
                        rmapmax = rmapmax + 1
                    end if
                 enddo
                 close(264)
              enddo
!              call ramg_dump_i(reducemap,nr,nr,'reducemap ')
              call ramg_dump_ic(rmap1d,2,rmapmax,'rmap1d    ')
              flag = 0
              write(fname,"((I1)(A))")nr,'-procs_case/map_nproc.dat'
              fname=trim(fname)
              do i=nr,1,-1
                 fname1 = trim(fname)//cname(i)
                 open(264,file=trim(fname1),status='unknown')
                 read(264,*)nj
                 do j=1,nj
                    read(264,*)mi,mj
                    if (flag(mj).eq.0) then
                        flag(mj) = i
                    else if (flag(mj).le.nr) then
                        !if ((flag(mj).eq.2).and.(i.eq.4)) then
                      !write(*,*)'mcheck!',mj,reducemap(flag(mj),i)
                        !endif
                        flag(mj) = iabs(reducemap(flag(mj),i))
                    end if
                 enddo
                 close(264)
              enddo
              call ramg_dump_i(flag,nshg,1,'initflagbc')
!              flag = myrank + 1
          ENDIF
          return
      end if

      IF (numpe.gt.1) THEN
      numtask = ilwork(1)
      allocate(neimap(numtask))
      itkbeg = 1
      do i=1,numtask
         iacc = ilwork(itkbeg+2)
         iother = ilwork(itkbeg+3)
         numseg = ilwork(itkbeg+4)
         if (iacc.eq.0) then !slave
             neimap(i) = -(iother+1)
         else !master
             neimap(i) = (iother+1)
         endif
         do j=1,numseg
            isgbeg = ilwork(itkbeg+3+j*2)
            lenseg = ilwork(itkbeg+4+j*2)
            isgend = isgbeg + lenseg - 1
            flag(isgbeg:isgend) = neimap(i)!iother+1!note*
         enddo
         itkbeg = itkbeg + 4 + 2*numseg
      enddo

      !call ramg_dump_i(flag,nshg,1,'amgparamap')

      !if (iamg_reduce.gt.0) then
          !call ramg_dump_i(neimap,numtask,1,'neimap    ')
      !endif
      deallocate(neimap)
      ENDIF

      ! note*: +1 to make paramap array range from 1 to numprocs
      ! this will avoid 0.
      ! n=proc indicates neighbour
      ! n=-proc indicates no coarsening necessary or other info


      end subroutine ! ramg_initBCflag

!****************************************************************
!      commOut_i : commOut for integer arrays, pack commOut
!****************************************************************
      subroutine commOut_i(global,ilwork,n,iper,iBC,BC)
      include "common.h"
      integer global(nshg,n)
      real*8  BC(nshg,ndofBC)
      integer ilwork(nlwork),iper(nshg),iBC(nshg)
      ! temp array
      real(kind=8)  aglobal(nshg,n)
      integer i,j
      do i=1,nshg
      do j=1,n
         aglobal(i,j) = global(i,j)
      enddo
      enddo
      call commOut(aglobal,ilwork,n,iper,iBC,BC)
      do i=1,nshg
      do j=1,n
         global(i,j)=aglobal(i,j)
      enddo
      enddo
      end subroutine ! commOut_i

!****************************************************************
!      ramg_commOut: commOut for amg array in higher level
!****************************************************************
      subroutine ramg_commOut(global,level,ilwork,redun,iper,iBC,BC)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in)                            :: level
      integer,intent(in) :: redun
      real(kind=8),intent(inout),dimension(amg_nshg(level),redun) 
     &           :: global
      integer,intent(in),dimension(nlwork)        :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      integer :: numtask,iacc,itag(ilwork(1)),i,j,m,k
      integer :: iother(ilwork(1)),ierr,iotheri
      integer :: numseg(ilwork(1))
      integer :: mstat(MPI_STATUS_SIZE,ilwork(1))
      integer :: req(ilwork(1))
      type(r2d),dimension(ilwork(1)) :: tmparray

      if (numpe.le.1) return

      if (level.eq.1) then
         call commOut(global,ilwork,1,iper,iBC,BC)
         return
      endif
      
      numtask = amg_ilwork(level)%p(1)

      lother = -1
      itkbeg=1
      do i=1,numtask
         itag(i)=ilwork(itkbeg+1)
         itkbeg=itkbeg+4+2*ilwork(itkbeg+4)
      enddo
      itkbeg=1
      do i=1,numtask
         iother(i)=amg_ilwork(level)%p(itkbeg+1)
         numseg(i)=amg_ilwork(level)%p(itkbeg+2)
         allocate(tmparray(i)%p(numseg(i),redun))
         itkbeg=itkbeg+2
         do j=1,numseg(i)
            itkbeg=itkbeg+1
         tmparray(i)%p(j,:)=global(amg_ilwork(level)%p(itkbeg),:)
         enddo
      enddo

      if ((myrank.eq.1).and.(level.eq.2)) then ! debug
!          write(*,*)'mcheck debug ilwork',amg_nlwork(2)
!          call ramg_dump(global,amg_nshg(level),'debuglobal')
!      call ramg_dump_i(amg_ilwork(2)%p,amg_nlwork(2),1,'debugilwok')
!          call ramg_dump(tmparray(1)%p,numseg(1),'debuglvl2 ')
      endif

      IF (.TRUE.) THEN
      m =0
      do i=1,numtask
         m = m+1
         iotheri = iabs(iother(i))-1
         if (iother(i)>0) then ! master send
!      write(*,*)'mcheck commou send',myrank,iotheri,numseg(i),itag(i)
            call MPI_ISEND(tmparray(i)%p(1,1),redun*numseg(i),
     &      MPI_DOUBLE_PRECISION,iotheri,itag(i),MPI_COMM_WORLD,
     &      req(m),ierr)
         else ! slave receive
!      write(*,*)'mcheck commou recv',myrank,iotheri,numseg(i),itag(i)
            call MPI_IRECV(tmparray(i)%p(1,1),redun*numseg(i),
     &      MPI_DOUBLE_PRECISION,iotheri,itag(i),MPI_COMM_WORLD,
     &      req(m),ierr)
         endif
      enddo

!      write(*,*)'mcheck commout m,',m

      call MPI_WAITALL(m,req,mstat,ierr)
      ENDIF

      ! put back

      itkbeg=1
      do i=1,numtask
         numseg(i)=amg_ilwork(level)%p(itkbeg+2)
         itkbeg=itkbeg+2
         do j=1,numseg(i)
            itkbeg=itkbeg+1
            global(amg_ilwork(level)%p(itkbeg),:)=tmparray(i)%p(j,:)
         enddo
      enddo

      do i=1,numtask
         deallocate(tmparray(i)%p)
      enddo
      
      end subroutine ! ramg_commOut

      subroutine ramg_commIn(global,level,ilwork,redun,iper,iBC,BC)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer,intent(in)                            :: level
      integer,intent(in) :: redun
      real(kind=8),intent(inout),dimension(amg_nshg(level),redun) 
     &           :: global
      integer,intent(in),dimension(nlwork)        :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      integer :: numtask,iacc,itag(ilwork(1)),i,j,m,k
      integer :: iother(ilwork(1)),ierr,iotheri
      integer :: numseg(ilwork(1))
      integer :: mstat(MPI_STATUS_SIZE,ilwork(1))
      integer :: req(ilwork(1))
      type(r2d),dimension(ilwork(1)) :: tmparray

      if (numpe.le.1) return

      if (level.eq.1) then
         call commIn(global,ilwork,1,iper,iBC,BC)
         return
      endif
      
      numtask = amg_ilwork(level)%p(1)

      lother = -1
      itkbeg=1
      do i=1,numtask
         itag(i)=ilwork(itkbeg+1)
         itkbeg=itkbeg+4+2*ilwork(itkbeg+4)
      enddo
      itkbeg=1
      do i=1,numtask
         iother(i)=amg_ilwork(level)%p(itkbeg+1)
         numseg(i)=amg_ilwork(level)%p(itkbeg+2)
         allocate(tmparray(i)%p(numseg(i),redun))
         itkbeg=itkbeg+2
         do j=1,numseg(i)
            itkbeg=itkbeg+1
         tmparray(i)%p(j,:)=global(amg_ilwork(level)%p(itkbeg),:)
         enddo
      enddo

      IF (.TRUE.) THEN
      m =0
      do i=1,numtask
         m = m+1
         iotheri = iabs(iother(i))-1
         if (iother(i)>0) then ! master receive
            call MPI_IRECV(tmparray(i)%p(1,1),redun*numseg(i),
     &      MPI_DOUBLE_PRECISION,iotheri,itag(i),MPI_COMM_WORLD,
     &      req(m),ierr)
         else ! slave send
            call MPI_ISEND(tmparray(i)%p(1,1),redun*numseg(i),
     &      MPI_DOUBLE_PRECISION,iotheri,itag(i),MPI_COMM_WORLD,
     &      req(m),ierr)
         endif
      enddo

!      write(*,*)'mcheck commout m,',m

      call MPI_WAITALL(m,req,mstat,ierr)
      ENDIF

      ! put back

      itkbeg=1
      do i=1,numtask
         numseg(i)=amg_ilwork(level)%p(itkbeg+2)
         itkbeg=itkbeg+2
         do j=1,numseg(i)
            itkbeg=itkbeg+1
            if (iother(i)>0) then ! master
            global(amg_ilwork(level)%p(itkbeg),:)=
     &      global(amg_ilwork(level)%p(itkbeg),:)+tmparray(i)%p(j,:)
            else
            global(amg_ilwork(level)%p(itkbeg),:)=0
            endif
         enddo
      enddo

      do i=1,numtask
         deallocate(tmparray(i)%p)
      enddo

      end subroutine ! ramg_commIn

      subroutine ramg_mapv2g(level,carray,garray)
      ! map to level 0 array
      use ramg_data
      include "common.h"
      integer,intent(in)            :: level
      real(kind=8),intent(inout),dimension(nshg) :: garray
      real(kind=8),intent(in),dimension(amg_nshg(level)) :: carray

      integer :: i,j,k
      integer,dimension(amg_nshg(level)) :: revmap
      revmap = 0
      garray(:) = 0
      j = 1
      do i=1,nshg
         if (amg_cfmap(i).ge.level) then
             revmap(j) = i
             j = j+1
         end if
      enddo
      do i=1,j-1
         garray(revmap(i)) = carray(i)
      enddo
      
      end subroutine !ramg_mapv2g

      subroutine ramg_mapg2v(level,carray,garray)
      use ramg_data
      include "common.h"
      integer,intent(in)            :: level
      real(kind=8),intent(in),dimension(nshg) :: garray
      real(kind=8),intent(inout),dimension(amg_nshg(level)) ::carray
      integer :: i,j,k
      carray(:) = 0
      j = 1
      do i=1,nshg
         if (amg_cfmap(i).ge.level) then
             carray(j) = garray(i)
             j = j+1
         end if
      enddo
      end subroutine !ramg_mapg2v
