      subroutine ramg_ggb_setup(colm,rowp,lhsK,lhsP,
     &                          ilwork,BC,iBC,iper)
      
      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      include 'debug.h'
 
      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp

      real(kind=8),intent(in),dimension(9,nnz_tot)  :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)  :: lhsP

      integer,intent(in),dimension(nlwork)      :: ilwork
      integer,intent(in),dimension(nshg)        :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC
      integer   iparam(11),ipntr(14)

      integer gcomm

      logical   selectncv(maxncv)
      real(kind=8) :: d(maxncv,3),dr(maxncv),di(maxncv),
     &                workev(3*maxncv),
     &                workl(3*maxncv*maxncv+6*maxncv)

      real(kind=8),allocatable,dimension(:) ::  resid,workd

      character  bmat*1,which*2
      integer   :: ido,n,nx,nev,ncv,lworkl,info,ierr,
     &             i,j,ishifts,maxitr,model,nconv
      real(kind=8)  :: tol,sigmar,sigmai
      logical     ::  first,rvec

      real(kind=8) :: rtemp,tnorm,tnorm1
      real(kind=8),allocatable,dimension(:) :: tv,tw
      real(kind=8),allocatable,dimension(:,:):: teA,tramg_ev

      integer,allocatable,dimension(:) :: gmap,grevmap

      integer :: p,q,nsize,step,asize,gsize
      character fname*20

      if (maxnev.le.0) return
      if (ramg_setup_flag.ne.0) return

      if (numpe.gt.1) then
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end if

      asize = amg_nshg(1)
      allocate(gmap(asize))
      allocate(grevmap(asize))

      call ramg_generate_gmap(ilwork,asize,nsize,gmap,grevmap,1)

      call MPI_AllReduce(nsize,gsize,1,MPI_INTEGER,MPI_MAX,
     &                   MPI_COMM_WORLD,ierr)

      !write(*,*)'mcheck ggb:',myrank,asize,nsize,gsize

      allocate(resid(gsize))
      allocate(workd(3*gsize))
      !***********   Start of Parameter setting *******
      
      ndigit = -3
      logfil = 6
      mnaitr = 0
      mnapps = 0
      mnaupd = 1 ! controls the output (verbose) mode
      mnaup2 = 0
      mneigh = 0
      mneupd = 1

      nev = maxnev
      ncv = maxncv
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
      allocate(tv(asize))
      allocate(tw(asize))

      step = 1
      if (myrank.eq.master) then
      write(*,*)'AMG GGB: calculating eigenvalue/vectors'
      endif
      tramg_ev = 0
      gcomm = MPI_COMM_WORLD

      call pdnaupd(gcomm,
     &               ido,bmat,nsize,which,nev,tol,resid,ncv,
     &               tramg_ev,
     &               nsize,iparam,ipntr,workd,workl,lworkl,
     &               info)
      do while ((ido.eq.1).or.(ido.eq.-1))
         call ramg_ggb_av(workd(ipntr(2)),workd(ipntr(1)),
     &                     asize,nsize,gmap,grevmap,
     &                     1,
     &                     colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
         call pdnaupd(gcomm,
     &               ido,bmat,nsize,which,nev,tol,resid,ncv,
     &               tramg_ev,
     &               nsize,iparam,ipntr,workd,workl,lworkl,
     &               info)
         step = step + 1
      enddo
      !write(*,*)'mcheck: ggb over pdnaupd'

      if (info.lt.0) then
          write(*,*)'mcheck: ggb: info:',info
      else
          rvec = .false.
          call pdneupd (gcomm,rvec,'A',selectncv,dr,di,
     &         tramg_ev,nsize,
     &         sigmar,sigmai,workev,bmat,nsize,which,nev,tol,
     &         resid,ncv,tramg_ev,
     &         nsize,iparam,ipntr,workd,workl,
     &         lworkl,ierr) 
          !write(*,*)'mcheck:ggb:',ierr,'iconv:',iparam(5)
      end if

      allocate(ramg_ggb_ev(asize,maxnev))
      do i=1,maxnev
        call ramg_ggb_G2A(ramg_ggb_ev(:,i),tramg_ev(:,i),
     &                    asize,nsize,1,gmap,grevmap)
      enddo


      ! Set matrix eA for GGB V-cycle

      allocate(ramg_ggb_eA(maxnev,maxnev))
      allocate(teA(asize,maxnev))

      ramg_ggb_eA = 0
      teA = 0

      do i=1,maxnev
         call ramg_calcAv_g(1,teA(:,i),ramg_ggb_ev(:,i),
     &                      colm,rowp,lhsK,lhsP,
     &                      ilwork,BC,iBC,iper,1)
      enddo
      do i=1,maxnev
            do j=1,maxnev
               do p=1,asize
                  ramg_ggb_eA(j,i) = ramg_ggb_eA(j,i) + 
     &            ramg_ggb_ev(p,j)*teA(p,i)
               enddo
            enddo
      enddo

      ! communicate to complete GGB matrix
      if (numpe.gt.1) then
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         do i=1,maxnev
            do j=1,maxnev
           call MPI_Allreduce(ramg_ggb_eA(j,i),rtemp,1,
     &          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
           ramg_ggb_eA(j,i) = rtemp      
            enddo
          enddo
      end if

      !call MPI_Barrier(MPI_COMM_WORLD,ierr)

      deallocate(teA)

      deallocate(tramg_ev)
      deallocate(tv)
      deallocate(tw)
      deallocate(gmap)
      deallocate(grevmap)
      deallocate(resid)
      deallocate(workd)

      ramg_time(11:30) = 0
      !write(*,*)'mcheck: over at g-cycle'
     
      end subroutine !ramg_ggb_setup

      subroutine ramg_generate_gmap(ilwork,asize,nsize,gmap,grevmap,
     &           level)

      use ramg_data
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      
      integer,intent(inout) :: nsize
      integer,intent(in),dimension(nlwork) :: ilwork
      integer,intent(in) :: asize
      integer,intent(inout),dimension(asize) :: gmap,grevmap
      integer,intent(in) :: level

      integer  :: numtask,i,j,k,m,p,newn
      integer  :: itag,iacc,iother,numseg,isgbeg,itkbeg
      integer,allocatable,dimension(:) :: neimap

      do i=1,asize
         gmap(i) = i
         grevmap(i) = i
      enddo

      if (numpe.le.1) then
          nsize = asize
          return 
      end if

      numtask = ilwork(1)
      m = 0
      itkbeg = 1
      allocate(neimap(numpe))
      neimap = 0

      newn = 0

      do i=1,numtask
         m = m+1
         itag = ilwork(itkbeg+1)
         iacc = ilwork(itkbeg+2)
         iother = ilwork(itkbeg+3)
         numseg = ilwork(itkbeg+4)
         isgbeg = ilwork(itkbeg+5)
         if (iacc.eq.0) then !slave
             neimap(iother+1) = 1
!             do j=1,numseg
!                k = ilwork(itkbeg+3+2*j)
!                gmap(k) = -iother
!             enddo
         end if
         m = m + 1
         itkbeg = itkbeg + 4 + 2*numseg
      enddo

      k = 1
      do i=1,asize
         m = iabs(amg_paramap(level)%p(i))
         if (neimap(m).eq.0) then ! self or master
            gmap(i) = k
            grevmap(k) = i
            k = k+1
         end if
      enddo
      nsize = k-1
      deallocate(neimap)

      !write(*,*)'mcheck gmap,',myrank,asize,nsize

      !call ramg_dump_i(gmap,asize,1,'gmapdump  ')
      !call ramg_dump_i(grevmap,asize,1,'grevmap   ')

      end subroutine !ramg_generate_gmap

      subroutine ramg_ggb_av(
     &               gramg_sol,gramg_rhs,
     &               asize,nsize,gmap,grevmap,
     &               clevel,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper
     &           )
      use ramg_data
      include "common.h"
      include "mpif.h"
      
      !arguments
      integer, intent(in)         :: asize,nsize
      integer, intent(in)         :: clevel
      real(kind=8),intent(inout),dimension(nsize)
     &               :: gramg_sol
      real(kind=8),intent(in),dimension(nsize) :: gramg_rhs
      integer, intent(in),dimension(asize) :: gmap,grevmap

      integer,intent(in),dimension(nshg+1)       :: colm
      integer,intent(in),dimension(nnz_tot)      :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)  :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)  :: lhsP
      integer,intent(in),dimension(nlwork)          :: ilwork
      integer,intent(in),dimension(nshg)            :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC) :: BC

      !local
      real(kind=8),dimension(amg_nshg(clevel)) :: vF
      real(kind=8),dimension(asize) :: ramg_sol,ramg_rhs

      call ramg_ggb_G2A(ramg_rhs,gramg_rhs,asize,nsize,1,gmap,grevmap)

      vF = 0
      
      call ramg_calcAv_g(clevel,vF,ramg_rhs,colm,rowp,lhsK,lhsP,
     &                   ilwork,BC,iBC,iper,1)

      call ramg_V_cycle(ramg_sol,vF,1,
     &                  colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper)
      ramg_sol = ramg_rhs-ramg_sol

      call ramg_ggb_A2G(ramg_sol,gramg_sol,asize,nsize,1,
     &                  gmap,grevmap)

      end subroutine !ramg_ggb_av

      
      subroutine ramg_G_cycle(
     &               ramg_sol,ramg_rhs,
     &               clevel,
     &      colm,rowp,lhsK,lhsP,ilwork,BC,iBC,iper
     &           )
      use ramg_data
      include "common.h"
      include "mpif.h"
      
      !arguments
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

      !local
      real(kind=8),dimension(maxnev) :: vR,rtemp
      real(kind=8),dimension(amg_nshg(clevel)) :: ggberr
      real(kind=8),dimension(maxnev,maxnev) :: mtxggbA
      integer,dimension(maxnev) :: indx
      real(kind=8) :: d
      integer :: i,j

      if (maxnev.le.0) return

!      ggberr = 0
      call ramg_calcAv_g(clevel,ggberr,ramg_sol,colm,rowp,lhsK,lhsP,
     &                   ilwork,BC,iBC,iper,1)
      ggberr = ggberr - ramg_rhs

      ! restriction
      vR = 0
      do i=1,maxnev
         do j=1,amg_nshg(clevel)
            vR(i) = vR(i)+ramg_ggb_eV(j,i)*ggberr(j)
         enddo
      enddo
      rtemp = 0
      if (numpe.gt.1) then
          call MPI_Allreduce(vR,rtemp,maxnev,
     &            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

      vR = rtemp
      end if

      ! solve maxnev by maxnev matrix
      mtxggbA = ramg_ggb_eA
      call ludcmp(mtxggbA,maxnev,maxnev,indx,d)
      call lubksb(mtxggbA,maxnev,maxnev,indx,vR)

      ! prolongation
      ggberr = 0
      do i=1,maxnev
         do j=1,amg_nshg(clevel)
            ggberr(j) = ggberr(j)+ramg_ggb_eV(j,i)*vR(i)
         enddo
      enddo

      ramg_sol = ramg_sol - ggberr

      end subroutine !ramg_G_cycle

      subroutine ramg_ggb_G2A(avec,gvec,asize,gsize,redun,
     &                           gmap,grevmap)
      integer,intent(in) :: asize, gsize
      integer,intent(in) :: redun
      integer,intent(in),dimension(asize) :: gmap,grevmap
      real(kind=8),intent(inout),dimension(asize,redun) :: avec
      real(kind=8),intent(in),dimension(gsize,redun) :: gvec
      integer :: i

      avec = 0
      do i=1,gsize
         avec(grevmap(i),:) = gvec(i,:)
      enddo

      end subroutine ! ramg_ggb_G2A

      
      subroutine ramg_ggb_A2G(avec,gvec,asize,gsize,redun,
     &           gmap,grevmap)
      integer,intent(in) :: asize, gsize
      integer,intent(in) :: redun
      integer,intent(in),dimension(asize) :: gmap,grevmap
      real(kind=8),intent(in),dimension(asize,redun) :: avec
      real(kind=8),intent(inout),dimension(gsize,redun) :: gvec
      integer :: i

      gvec = 0
      do i=1,asize
         if (gmap(i).gt.0) then
             gvec(gmap(i),:) = avec(i,:)
         end if
      enddo

      end subroutine ! ramg_ggb_A2G
 
