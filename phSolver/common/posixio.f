      module posixio
      interface 
        subroutine posixio_setup(handle, mode) 
     &   bind(C, NAME='posixio_setup')
        use :: iso_c_binding
          type(c_ptr) :: handle
          character(c_char), value :: mode
        end subroutine
      end interface
      end module
