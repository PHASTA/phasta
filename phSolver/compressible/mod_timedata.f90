      module timedataC
            
        integer ntspts, freq, iterat, varcod, nbuff
        integer iblkts
        real*8  tolpt
        logical exts
              
        integer,  allocatable :: statptts(:,:)
        real*8,  allocatable :: ptts(:,:)
        real*8,  allocatable :: parptts(:,:)
        real*8,  allocatable :: varts(:,:)
        
        !Extra variables for buffering and writting data     
        character(len=60) vartsIOFrmtStr
        integer :: ivartsBuff 
        integer, allocatable, dimension(:) :: ivarts
        integer, allocatable, dimension(:) :: ivartsg
        real*8, allocatable, dimension(:) :: vartssoln      
        real*8, allocatable, dimension(:) :: vartssolng
        real*8, allocatable, dimension(:,:,:) :: vartsbuff
        integer, allocatable, dimension(:) :: vartsbuffStep

        logical :: vartsResetBuffer     !Used so that the buffer does
                                        !not need to be reset 
                                        !immediately after writing. 

      end module





      module pp_data

      integer numppnodes, ppfreq

      integer,  allocatable :: ppnodes(:,:)

      end module






