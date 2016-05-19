      module mkdir_mod
        implicit none
        interface 
          integer function c_mkdir(path) bind(C, name="ph_mkdir")
            use iso_c_binding
            character(c_char) :: path(*)
          end function
        end interface
        contains 
        subroutine mkdir(path, err)
          use iso_c_binding
          character(len=*) :: path
          integer, optional, intent(out) :: err
          integer :: loc_err
          loc_err = c_mkdir(path//c_null_char)
          if( present(err) ) err = loc_err
        end subroutine
      end module mkdir_mod
