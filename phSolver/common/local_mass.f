      module local_mass
      
      parameter (MAXBLK2 = 5000)
      
      integer :: iblock = 0
      integer :: have_local_mass = 0
      
      type r2da
      real*8, pointer :: p(:,:,:)
      end type
      
      type (r2da), dimension(MAXBLK2) :: lmassinv
      
      end module
 
