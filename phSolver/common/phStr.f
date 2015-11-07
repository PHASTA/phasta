      module phstr
      use :: iso_c_binding
      interface
        subroutine phstr_appendInt(str, val)
     &   bind(C, NAME='phstr_appendInt')
          use :: iso_c_binding
          character(c_char) :: str(*)
          integer(c_int), value, intent(in) :: val
        end subroutine
      end interface

      interface
        subroutine phstr_appendDbl(str, val)
     &   bind(C, NAME='phstr_appendDbl')
          use :: iso_c_binding
          character(c_char) :: str(*)
          real(c_double), value, intent(in) :: val
        end subroutine
      end interface

      interface
        subroutine phstr_appendStr(dest, src)
     &   bind(C, NAME='phstr_appendStr')
          use :: iso_c_binding
          character(c_char) :: dest(*)
          character(c_char) :: src(*)
        end subroutine
      end interface
      
      end module
