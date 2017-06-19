      module phiotimer
      enum, bind(C)
        enumerator :: GEOMBC_READ, RESTART_READ, RESTART_WRITE
      end enum
      interface 
        subroutine phastaio_setfile(fileidx)
     &   bind(C, NAME='phastaio_setfile')
        use :: iso_c_binding
          integer(c_int), intent(in), value:: fileidx
        end subroutine
      end interface
      end module
