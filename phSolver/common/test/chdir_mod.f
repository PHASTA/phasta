      module chdir_mod
        implicit none
        interface 
          integer function c_chdir(path) bind(C, name="chdir")
            use iso_c_binding
            character(c_char) :: path(*)
          end function
        end interface
        contains 
        subroutine chdir(path, err)
          use iso_c_binding
          character(len=*) :: path
          integer, optional, intent(out) :: err
          integer :: loc_err
          loc_err = c_chdir(path//c_null_char)
          if( present(err) ) err = loc_err
        end subroutine
      end module chdir_mod
