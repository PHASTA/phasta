      module phio
      use :: iso_c_binding
      type(c_ptr) :: fhandle
      interface 
        subroutine phio_openfile(fname, handle) 
     &   bind(C, NAME='phio_openfile')
        use :: iso_c_binding
          character(c_char), intent(in) :: fname(*)
          type(c_ptr), value :: handle
        end subroutine
      end interface
      interface 
        subroutine phio_closefile(handle) 
     &   bind(C, NAME='phio_closefile')
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
        subroutine phio_appendInt(str, val)
     &   bind(C, NAME='phio_appendInt')
        use :: iso_c_binding
          character(c_char) :: str(*)
          integer(c_int), value, intent(in) :: val
        end subroutine
      end interface
      interface
        subroutine phio_constructName(handle, inName, outName)
     &   bind(C, NAME='phio_constructName')
        use :: iso_c_binding
          type(c_ptr), value :: handle
          character(c_char) :: inName(*)
          character(c_char) :: outName(*)
        end subroutine
      end interface

      end module
