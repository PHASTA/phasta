       module pointer_data
c
c.... maximum number of blocks
c
         parameter ( MAXBLK2 = 50000 ) ! Note compiler was complaining 
c                                       because MAXBLK in common.h be careful
c    					to chang both places
c
c.... data type definitions
c
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
         type i2d64
           integer*8, pointer :: p(:,:)
         end type
c
         type i3d
           integer, pointer :: p(:,:,:)
         end type
c
c.... pointer declarations
c
         type (i1d), dimension(MAXBLK2) ::  mmat,  mmatb
         type (i2d), dimension(MAXBLK2) ::  mien
         type (i2d64), dimension(MAXBLK2) ::  mienG
         type (i2d), dimension(MAXBLK2) ::  mienb,  miBCB
         type (r2d), dimension(MAXBLK2) ::  mxmudmi
         type (r3d), dimension(MAXBLK2) ::  mBCB
c
         real*8, allocatable :: gmass(:)
       end module
c
c
c
