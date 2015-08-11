      module syncio
      use :: iso_c_binding
      interface 
        subroutine syncio_setup_read(nfiles, handle) 
     &   bind(C, NAME='syncio_setup_read')
        use :: iso_c_binding
          integer(c_int), value :: nfiles
          type(c_ptr) :: handle
        end subroutine
        subroutine syncio_setup_write(nfiles, nfields, nppf, handle) 
     &   bind(C, NAME='syncio_setup_write')
        use :: iso_c_binding
          integer(c_int), value :: nfiles
          integer(c_int), value :: nfields
          integer(c_int), value :: nppf
          type(c_ptr) :: handle
        end subroutine
      end interface
      end module
