      module streamio
      use :: iso_c_binding
      type(c_ptr) :: grstream
      interface 
        subroutine streamio_setup(stream, handle) 
     &   bind(C, NAME='streamio_setup')
        use :: iso_c_binding
          type(c_ptr), value :: stream
          type(c_ptr) :: handle
        end subroutine
      end interface
      end module
