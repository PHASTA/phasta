      module phio
      interface 
        subroutine phio_openfile_read(fname, nfiles, handle) 
     &   bind(C, NAME='phio_openfile_read')
        use :: iso_c_binding
          character(c_char), intent(in) :: fname(*)
          integer(c_int), intent(in) :: nfiles
          type(c_ptr) :: handle
        end subroutine
      end interface
      interface 
        subroutine phio_closefile_read(handle) 
     &   bind(C, NAME='phio_closefile_read')
        use :: iso_c_binding
          type(c_ptr), value :: handle
        end subroutine
      end interface
      interface 
        subroutine phio_readheader(handle, phrase, vals, nvals, 
     &                             datatype, iotype) 
     &   bind(C, NAME='phio_readheader')
        use :: iso_c_binding
          type(c_ptr), value :: handle
          character(c_char), intent(in) :: phrase(*)
          type(c_ptr), value :: vals
          integer(c_int), intent(in) :: nvals
          character(c_char), intent(in) :: datatype(*)
          character(c_char), intent(in) :: iotype(*)
        end subroutine
      end interface
      interface 
        subroutine phio_readdatablock(handle, phrase, vals, nvals, 
     &                                datatype, iotype) 
     &   bind(C, NAME='phio_readdatablock')
        use :: iso_c_binding
          type(c_ptr), value :: handle
          character(c_char), intent(in) :: phrase(*)
          type(c_ptr), value :: vals
          integer(c_int), intent(in) :: nvals
          character(c_char), intent(in) :: datatype(*)
          character(c_char), intent(in) :: iotype(*)
        end subroutine
      end interface
      end module

