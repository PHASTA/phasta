          module BCsfIDmap 
              integer ninlet
              integer, allocatable, dimension(:) :: inlf
              real*8, allocatable, dimension(:,:) :: prof_inlet ! profile for inlet velocity
              integer noutlet
              integer, allocatable, dimension(:) :: outf 
          end module

