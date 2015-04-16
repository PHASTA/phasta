!************************************
!     RAMG Extract:
!      extract ppe and scale ppe
!**********************************************************
     
!***********************************************************
!    ramg_extract_ppe
!    prepare PPE matrix
!     input: boundary, lhsK,lhsP
!     output: PPE & RHS at level 1
!**********************************************************
      subroutine ramg_extract_ppe(
     &              colm,rowp,lhsK,lhsP,
     &              ilwork,BC,iBC,iper)
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
!      implicit none
      
     
      !***********parameters**************
      !the matrix
!      integer,intent(in)                               :: nshg
!      integer,intent(in)                        :: nnz_tot
!      integer,intent(in)                               :: nflow
      !the matrix
      integer,intent(in),dimension(nshg+1)             :: colm
      integer,intent(in),dimension(nnz_tot)            :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)     :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)     :: lhsP
!      real(kind=8),dimension(:,:),allocatable          :: lhsGP
!      integer,dimension(:),allocatable                 :: lhsGProwp
!      integer,dimension(:),allocatable                 :: lhsGPcolm
      type(r2d) :: lhsGP
      type(i1d) :: lhsGProwp,lhsGPcolm
      ! the boundary info
      integer, intent(in), dimension(nlwork)           :: ilwork
      integer, intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)   :: BC

      !*********parameters end**************

      !****** local variables **********
      real(kind=8),dimension(nshg,4)     :: mflowDiag,pflowDiag
      real(kind=8)                             :: rtemp,rt,rtp,rtn
      real(kind=8),dimension(nshg)                    :: rtest
      integer                     :: i,j,k,m,n,p,q,ki,kj,ni,nj,qq
      integer :: cki,ckj,ckk
      integer,dimension(nshg)                  :: temprow
      integer                                         :: mnnz

      logical  :: extentri

      character                                       :: fname*80
      !****** end local variables ******
    
      if (ramg_setup_flag .ne. 0) return

      !*** calculate memory for PPE ***
      !*** using the fact that matrix is symmetric ***
      !*** lhs = GT_ik * KI_kk * G_kj + C_ij ***
      !*** rhs = -GT_ik * KI_kk * Rm_k - Rc_k ***
      !*** where G and GT has same sparsity pattern if 
      !*** only consider the main matrix (M,3)

      if ( amg_nshg(1) .gt.0 ) then
          call ramg_deallocate(1)
          deallocate(ramg_flowDiag%p)
          deallocate(amg_cfmap)
      end if

      ! colm
      call ramg_allocate(1,nshg,0,1)
      allocate(amg_paramap(1)%p(nshg))
      allocate(amg_paraext(1)%p(nshg))

      mflowDiag(:,:)=0
 
!      if ( (numpe.gt.1) ) then !.and.(iamg_reduce.gt.0) ) then
!          call ramg_prep_reduce_map
          !call ramg_dump_i(ncorp_map,nshg,1,'map_nproc ')
          !deallocate(ncorp_map)
!      endif

      call drvLesPrepDiag(mflowDiag,ilwork,iBC,BC,iper,
     &                    rowp,colm,lhsK,lhsP)

      !*** Complete diagonal values in mflowdiag,pflowdiag ***
      !*** lhs K, G, C should have "global" values on diagonal ***
      !call ramg_dump_rn_map(mflowDiag,nshg,4,'mflowDiagb')
      call commIn(mflowdiag,ilwork,4,iper,iBC,BC)
      call commOut(mflowdiag,ilwork,4,iper,iBC,BC)
      !call ramg_dump_rn_map(mflowDiag,nshg,4,'mflowDiaga')
      ! call ramg_dump(mflowDiag,nshg,'flowdiag  ')

      allocate(ramg_flowDiag%p(nshg,4)) 
      allocate(amg_cfmap(nshg))

      amg_cfmap = 1
      ! amg_cfmap is such a variable that it keeps
      ! coarsening information through the coarsest level
      ! in a finest level array to enable communication.
      ! The structure is : 111223322113231122...
      ! obviously the nodes that kept into next level got
      ! an extra 1

      do i=1,4
      ramg_flowDiag%p(:,i) = mflowDiag(:,i)
      do j=1,nshg
         if (mflowdiag(j,4).eq.0) then
             write(*,*)'mflowdiag zero at ',j,i
             stop
         end if
      enddo
      enddo

      call ramg_initBCflag(amg_paramap(1)%p,ilwork,BC,iBC,iper)

      call ramg_global_lhs(colm,rowp,lhsP,nnz_tot,
     &                    lhsGPcolm,lhsGProwp,lhsGP,4,
     &                    ilwork,BC,iBC,iper)

      ! prepare extended boundary in amg_paraext
      ni = 0
      nj = 0
      amg_paraext(1)%p = amg_paramap(1)%p
      do i=1,nshg
         if (amg_paramap(1)%p(i).ne.(myrank+1)) then
             ni = ni+1
             nj = nj+1
             do m=lhsGPcolm%p(i),lhsGPcolm%p(i+1)-1
                k = lhsGProwp%p(m)
                if (amg_paraext(1)%p(k).eq.(myrank+1)) then
                    nj = nj+1
                    amg_paraext(1)%p(k)=amg_paramap(1)%p(i)
                endif
             enddo
         endif
      enddo
      !write(*,*)'proc',myrank,ni,nj,nshg
      
      ! end of preparing extended boundary 
    
      extentri = .false. ! HAVE extra entries
      !extentri = .true. ! NO extra entries

      ! output K^{hat}^-1 and lhsP (G,C) to external file
      !write(fname,*)'Khatinv'
      !call ramg_dump_rn(mflowdiag,nshg,4,fname)
      !write(fname,*)'lhsP'
      !call ramg_dump_matlab_A(colm,rowp,lhsGP,nshg,nnz_tot,4,fname)
      
      mnnz = 1
      temprow = 0
      do i = 1,nshg                ! each row
        amg_A_colm(1)%p(i) = mnnz
        !q = iabs(amg_paramap(1)%p(i)) ! used only in reduced serial
        do m = lhsGPcolm%p(i),lhsGPcolm%p(i+1)-1
           k = lhsGProwp%p(m)
           do j = lhsGPcolm%p(k),lhsGPcolm%p(k+1)-1
              p = lhsGProwp%p(j)
              !qq = iabs(amg_paramap(1)%p(p)) ! used only in reduced
              !serial
              if (temprow(p).ne.i) then
!        if ((iamg_reduce.gt.1).and.(numpe.eq.1).and.(extentri))then
!                     if  ( (rmap1d(1,qq).eq.rmap1d(1,q)).or.
!     &                     (rmap1d(2,qq).eq.rmap1d(2,q)).or.
!     &                     (rmap1d(1,qq).eq.rmap1d(2,q)).or.
!     &                     (rmap1d(2,qq).eq.rmap1d(1,q)) ) then
!                       temprow(p) = i
!                       mnnz = mnnz + 1
!                     endif
!                   else
                   temprow(p) = i
                   mnnz = mnnz + 1
!                   endif
              end if
           end do
        end do
      end do
      amg_A_colm(1)%p(nshg+1) = mnnz
      mnnz = mnnz - 1
      !*** now we have amg_nnz(1) as nnz_tot for PPE
      call ramg_allocate(1,0,mnnz,1)
      !allocate(levelbaseflag(mnnz))
      !levelbaseflag = 0
     
      if (iamg_verb.gt.5) then
      write(6,7003) nshg,amg_nnz(1)
7003  format(/'nshg='i10',  nnz='i10)
      endif

      !*** end of PPE memory alloc ***
      
      !**** begin putting rowp and values into PPE **
      mnnz = 0
      temprow = 0


      ! rowp
      do i = 1,nshg
         mnnz = mnnz + 1
         amg_A_rowp(1)%p(mnnz) = i
         temprow(i) = i
         !levelbaseflag(mnnz) = 1
         !q = iabs(amg_paramap(1)%p(i)) !used in reduced serial
         do m = lhsGPcolm%p(i),lhsGPcolm%p(i+1)-1
            j = lhsGProwp%p(m)
            do n = lhsGPcolm%p(j),lhsGPcolm%p(j+1)-1
               k = lhsGProwp%p(n)
               !qq = iabs(amg_paramap(1)%p(k)) !reduced serial
               if ( temprow(k) .ne. i) then
!         if ((iamg_reduce.gt.1).and.(numpe.eq.1).and.(.true.))then
!                     if  ( (rmap1d(1,qq).eq.rmap1d(1,q)).or.
!     &                     (rmap1d(2,qq).eq.rmap1d(2,q)).or.
!     &                     (rmap1d(1,qq).eq.rmap1d(2,q)).or.
!     &                     (rmap1d(2,qq).eq.rmap1d(1,q)) ) then
!                       mnnz = mnnz + 1
!                       amg_A_rowp(1)%p(mnnz) = k
!                       temprow(k) = i
!                        levelbaseflag(mnnz+1) = 1
!                     endif
!         endif
!                   else
                   mnnz = mnnz + 1
                   amg_A_rowp(1)%p(mnnz) = k
                   temprow(k) = i
!                   endif
               end if
            enddo
         enddo
         do m=amg_A_colm(1)%p(i)+1,amg_A_colm(1)%p(i+1)-1
            ki = amg_A_rowp(1)%p(m)
            do n=m+1,amg_A_colm(1)%p(i+1)-1
               kj = amg_A_rowp(1)%p(n)
               if (kj.lt.ki) then
                   amg_A_rowp(1)%p(m) = kj
                   amg_A_rowp(1)%p(n) = ki
                   ki = kj
!                   p = levelbaseflag(n) 
!                   levelbaseflag(n) = levelbaseflag(m)
!                   levelbaseflag(m) = p
               end if
            enddo
         enddo
      enddo

!      j = 0
!      do i=1,mnnz
!         j = j + levelbaseflag(i)
!      enddo
!      write(*,*)'mcheck incomplete: ',j

      rt = 0

      ! matrix value
      ! cki,ckj,ckk, this is to avoid double summation on PPE,
      ! since both master and slave have complete lhsP, parts of it
      ! should be removed when constructing PPE. if (i,j) are both from
      ! slave
      do i=1,nshg
         do m = amg_A_colm(1)%p(i),amg_A_colm(1)%p(i+1)-1
            j = amg_A_rowp(1)%p(m)
            !if (levelbaseflag(m).eq.0) then
            !    amg_A_lhs(1)%p(m,1) = 0
            !    cycle
            !endif
            ki = lhsGPcolm%p(i)
            ni = lhsGPcolm%p(i+1)
            kj = lhsGPcolm%p(j)
            nj = lhsGPcolm%p(j+1)
            rtemp = 0
            ! question: the original lhsK,lhsP, are they symmetric?
            ! symmetric in nnz structure but not in value?
            kloop: do while ( ( ki.lt.ni ).and.(kj.lt.nj) )
              do while ((ki.lt.ni).and.
     &           (lhsGProwp%p(ki).lt.lhsGProwp%p(kj)) )
                 ki = ki + 1
              enddo
              if (ki.eq.ni) exit kloop
              do while ( (kj.lt.nj) .and. 
     &                   (lhsGProwp%p(kj).lt.lhsGProwp%p(ki)) )
                 kj = kj + 1
              enddo
              if (kj.eq.nj) exit kloop
              p = lhsGProwp%p(ki)
              q = lhsGProwp%p(kj)
              if (p.eq.q) then
c              if (amg_paramap(1)%p(p).eq.amg_paramap(1)%p(q)) then
                  k = q
               cki = amg_paramap(1)%p(i)
               ckj = amg_paramap(1)%p(j)
               ckk = amg_paramap(1)%p(k)
               if (.not.((iabs(cki).eq.iabs(ckj)).and.
     &             (iabs(ckj).eq.iabs(ckk)).and.
     &             (cki.lt.0).and.(ckj.lt.0).and.(ckk.lt.0))) then
                  rtemp = rtemp  
     &            + lhsGP%p(1,ki)*lhsGP%p(1,kj)*mflowDiag(k,1)**2
     &            + lhsGP%p(2,ki)*lhsGP%p(2,kj)*mflowDiag(k,2)**2
     &            + lhsGP%p(3,ki)*lhsGP%p(3,kj)*mflowDiag(k,3)**2
               endif
                  ki = ki+1
                  kj = kj+1
c              endif
              end if
            enddo kloop
            amg_A_lhs(1)%p(m,1)=rtemp
         enddo
      enddo
      do i=1,nshg
         ki = amg_A_colm(1)%p(i)+1
         mloop: do m=lhsGPcolm%p(i),lhsGPcolm%p(i+1)-1
            j = lhsGProwp%p(m)
            if (j.eq.i) then
            cki = amg_paramap(1)%p(i)
            ckj = amg_paramap(1)%p(j)
            if (.not.((iabs(cki).eq.iabs(ckj)).and.
     &                (cki.lt.0).and.(ckj.lt.0))) then
             amg_A_lhs(1)%p(amg_A_colm(1)%p(i),1) = 
     &          + amg_A_lhs(1)%p(amg_A_colm(1)%p(i),1) + lhsGP%p(4,m)
            endif
                cycle mloop
            end if
            do while(amg_A_rowp(1)%p(ki).lt.j)
               ki = ki+1
            enddo
            cki = amg_paramap(1)%p(i)
            ckj = amg_paramap(1)%p(j)
            ckk = amg_paramap(1)%p(amg_A_rowp(1)%p(ki))
            if (.not.((iabs(cki).eq.iabs(ckj)).and.
     &                (iabs(ckj).eq.iabs(ckk)).and.
     &                (cki.lt.0).and.(ckj.lt.0).and.(ckk.lt.0))) then
            amg_A_lhs(1)%p(ki,1) = amg_A_lhs(1)%p(ki,1) + lhsGP%p(4,m)
            endif
            ki = ki+1
         enddo mloop
      enddo

      !gtg: lhsgp(4,m), mflowdiag ignored

      deallocate(lhsGP%p)
      deallocate(lhsGProwp%p)
      deallocate(lhsGPcolm%p)

      if (.true.) then
          
      call ramg_global_lhs(amg_A_colm(1)%p,amg_A_rowp(1)%p,
     &                    amg_A_lhs(1)%p,amg_nnz(1),
     &                    lhsGPcolm,lhsGProwp,lhsGP,1,
     &                    ilwork,BC,iBC,iper)

      amg_A_colm(1)%p = lhsGPcolm%p
      amg_nnz(1) = lhsGPcolm%p(amg_nshg(1)+1)-1
      deallocate(amg_A_rowp(1)%p)
      deallocate(amg_A_lhs(1)%p)
      allocate(amg_A_rowp(1)%p(amg_nnz(1)))
      allocate(amg_A_lhs(1)%p(amg_nnz(1),1))
      amg_A_rowp(1)%p = lhsGProwp%p
      amg_A_lhs(1)%p = lhsGP%p
      deallocate(lhsGP%p)
      deallocate(lhsGProwp%p)
      deallocate(lhsGPcolm%p)

      endif
      
      if (.false.) then
          ! Dump for Matlab
         call ramg_dump_matlab_map(amg_A_colm(1)%p,amg_A_rowp(1)%p,
     &             amg_A_lhs(1)%p,
     &             amg_nshg(1),amg_nnz(1),1,
     &             'A0        ')
!         call ramg_dump_matlab_map(lhsGPcolm%p,lhsGProwp%p,
!     &             lhsGP%p,
!     &             amg_nshg(1),lhsGPcolm%p(amg_nshg(1)+1)-1,1,
!     &             'A0        ')
          !call ramg_dump(amg_ppeDiag(1)%p,nshg,'ppediag   ')
      end if


      if (.false.) then ! Check Diagonal Scaling, M-matrix
          open(264,file='mtcheck.dat',status='unknown')
          do i=1,amg_nshg(1)
             m = 0
             rtemp=amg_A_lhs(1)%p(amg_A_colm(1)%p(i),1)
             rtp = 0
             rtn = 0
             do j=amg_A_colm(1)%p(i)+1,amg_A_colm(1)%p(i+1)-1
                rtemp = rtemp + amg_A_lhs(1)%p(j,1)
                rt = amg_A_lhs(1)%p(j,1)
                if (rt.gt.0) then
                   if ( rt.gt.rtp ) rtp = rt
                end if
                if (rt.lt.0) then
                    if (rt.le.rtn) rtn = rt
                end if
             enddo
             rtp = rtp/ (amg_A_lhs(1)%p(amg_A_colm(1)%p(i),1))
             rtn = rtn/amg_A_lhs(1)%p(amg_A_colm(1)%p(i),1)
             write(264,"((I5)(A)(E9.3)(A)(E9.3)(A)(E9.3))")
     &         i,'  ',rtp,'  ',rtn,'  ',rtemp
          enddo
          close(264)
      end if

    ! scaled by diag(PPE) 1...1
       do i=1,nshg
         amg_ppeDiag(1)%p(i) = 
     &      1.0/sqrt(amg_A_lhs(1)%p(amg_A_colm(1)%p(i),1))
       enddo
     
      do i=1,nshg
         do m=amg_A_colm(1)%p(i),amg_A_colm(1)%p(i+1)-1
            j = amg_A_rowp(1)%p(m)
            amg_A_lhs(1)%p(m,1) = amg_A_lhs(1)%p(m,1) *
     &      amg_ppeDiag(1)%p(i)*amg_ppeDiag(1)%p(j)
         end do
         amg_A_rhs(1)%p(i) = amg_A_rhs(1)%p(i)*amg_ppeDiag(1)%p(i)
      end do

      ! If solve heat/scalar, should disable following 3 lines
      ! because we only have 1 scaling.
      do i=1,nshg
         amg_ppeDiag(1)%p(i) = amg_ppeDiag(1)%p(i)/mflowDiag(i,4)
      enddo

      if (.false.) then
          ! Dump for Matlab
         call ramg_dump_matlab_map(amg_A_colm(1)%p,amg_A_rowp(1)%p,
     &             amg_A_lhs(1)%p,
     &             amg_nshg(1),amg_nnz(1),1,
     &             'A0scaled  ')
!         call ramg_dump_matlab_map(lhsGPcolm%p,lhsGProwp%p,
!     &             lhsGP%p,
!     &             amg_nshg(1),lhsGPcolm%p(amg_nshg(1)+1)-1,1,
!     &             'A0        ')
          !call ramg_dump(amg_ppeDiag(1)%p,nshg,'ppediag   ')
      end if


      return

      end subroutine !ramg_extract_ppe
     
