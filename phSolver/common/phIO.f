      module phio
      use :: iso_c_binding
      type(c_ptr) :: fhandle
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
        subroutine phio_openfile_write(fname, nfiles, nfields,
     &   nppf, handle) 
     &   bind(C, NAME='phio_openfile_write')
        use :: iso_c_binding
          character(c_char), intent(in) :: fname(*)
          integer(c_int), intent(in) :: nfiles
          integer(c_int), intent(in) :: nfields
          integer(c_int), intent(in) :: nppf
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
        subroutine phio_closefile_write(handle) 
     &   bind(C, NAME='phio_closefile_write')
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
        subroutine phio_writeheader(handle, phrase, vals, nitems, ndata,
     &                             datatype, iotype) 
     &   bind(C, NAME='phio_writeheader')
        use :: iso_c_binding
          type(c_ptr), value :: handle
          character(c_char), intent(in) :: phrase(*)
          type(c_ptr), value :: vals
          integer(c_int), intent(in) :: nitems
          integer(c_int), intent(in) :: ndata
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
      interface 
        subroutine phio_writedatablock(handle, phrase, vals, nvals,
     &                                datatype, iotype)
     &   bind(C, NAME='phio_writedatablock')
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
        subroutine phio_appendStep(str, val)
     &   bind(C, NAME='phio_appendStep')
        use :: iso_c_binding
          character(c_char) :: str(*)
          integer(c_int), value, intent(in) :: val
        end subroutine
      end interface
      end module
