      module fncorpmod

      integer*8, allocatable :: fncorp(:)
      integer, allocatable :: ltg(:)

      contains

        subroutine destroyfncorp()
          if ( allocated(fncorp) ) then
            deallocate(fncorp)
          endif
          if ( allocated(ltg) ) then
            deallocate(ltg)
          endif
        end subroutine destroyfncorp

      end module
