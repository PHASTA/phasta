      module ramg_data

         parameter (MAXAMGLVL=15)
         
         type r1d
           real*8, pointer :: p(:)
         end type
c
         type r2d
           real*8, pointer :: p(:,:)
         end type
c
         type r3d
           real*8, pointer :: p(:,:,:)
         end type
c
         type i1d
           integer, pointer :: p(:)
         end type
c
         type i2d
           integer, pointer :: p(:,:)
         end type
c
         type i3d
           integer, pointer :: p(:,:,:)
         end type
c
         type i2dd
           type(i1d),pointer :: pp(:)
         end type
c
         type r2dd
           type(r1d),pointer :: pp(:)
         end type         
c
         type r12dd
            type(r2d),pointer :: pp(:)
         end type
c
         ! For communication pattern
         type(i2dd)            :: sub_map,sub_revmap,sub_rowpmap
         type(i2dd)       :: sub_colm,sub_colm2,sub_rowp,sub_rowp2
         type(r12dd)           :: sub_mtx,sub_mtx2
         type(i1d)             :: sub_nnz,sub_nnz2,sub_nshg


         ! For AMG

         integer,allocatable,dimension(:)          :: levelbaseflag
         integer,allocatable,dimension(:)     :: reducecolm
         integer,allocatable,dimension(:)     :: reducerowp
         real(kind=8),allocatable,dimension(:) :: reducelhs
         
         type(r2d),dimension(MAXAMGLVL)        :: amg_A_lhs
         type(i1d),dimension(MAXAMGLVL)        :: amg_A_rowp
         type(i1d),dimension(MAXAMGLVL)        :: amg_A_colm
         type(r1d),dimension(MAXAMGLVL)        :: amg_A_rhs

         integer,dimension(MAXAMGLVL)          :: amg_nshg
         integer,dimension(MAXAMGLVL)          :: amg_nnz
 
         type(i1d),dimension(MAXAMGLVL)        :: amg_paramap
         type(i1d),dimension(MAXAMGLVL)        :: amg_paraext
         type(i1d),dimension(MAXAMGLVL)        :: amg_ilwork
         integer,dimension(MAXAMGLVL)          :: amg_nlwork
         integer,dimension(:),allocatable      :: amg_cfmap

         integer                               :: ramg_flag
         real(kind=8),dimension(10)            :: ramg_window
         integer                               :: ramg_winct
         integer                               :: ramg_redo
         real(kind=8),dimension(30)            :: ramg_time

         integer                               :: ramg_verb
         integer                               :: ramg_setup_flag
         integer                               :: ramg_levelx

         integer,allocatable,dimension(:)      :: ncorp_map

         type(r2d)        :: ramg_flowDiag
         type(r1d),dimension(MAXAMGLVL)        :: amg_ppeDiag

         type(i1d),dimension(MAXAMGLVL)        :: CF_map,CF_revmap
         type(r1d),dimension(MAXAMGLVL)        :: I_fc
         type(i1d),dimension(MAXAMGLVL)        :: I_fc_rowp
         type(i1d),dimension(MAXAMGLVL)        :: I_fc_colm
         type(r1d),dimension(MAXAMGLVL)        :: I_cf
         type(i1d),dimension(MAXAMGLVL)        :: I_cf_rowp
         type(i1d),dimension(MAXAMGLVL)        :: I_cf_colm

         type(r2d),dimension(MAXAMGLVL)        :: sI_fc
         type(r2d),dimension(MAXAMGLVL)        :: sI_cf

         ! GGB parameters

         real(kind=8),dimension(:,:),allocatable :: ramg_ggb_ev
         real(kind=8),dimension(:,:),allocatable :: ramg_ggb_eA
         ! dense coarsest
         real(kind=8),dimension(:,:),allocatable :: cmtxA
         real(kind=8),dimension(:),allocatable :: cindx

         ! MLS paramters
         real(kind=8),dimension(MAXAMGLVL,10)    :: mlsCf,mlsOm
            ! for system:
         real(kind=8),dimension(MAXAMGLVL,10)    :: smlsCf,smlsOm

         ! reduced case
         integer,allocatable,dimension(:,:) :: reducemap,rmap1d
         integer :: rmapmax

      end module


