      module streamio
      use :: iso_c_binding
      type(c_ptr) :: geomRestartStream 
      bind(C, name='geomRestartStream') :: geomRestartStream
      type(c_ptr) :: restartStream
      bind(C, name='restartStream') :: restartStream
      interface 
        subroutine streamio_setup_read(handle, stream) 
     &   bind(C, NAME='streamio_setup_read')
        use :: iso_c_binding
          type(c_ptr) :: handle
          type(c_ptr), value :: stream
        end subroutine
      end interface
      interface 
        subroutine streamio_setup_write(handle, stream) 
     &   bind(C, NAME='streamio_setup_write')
        use :: iso_c_binding
        type(c_ptr) :: handle
        type(c_ptr), value :: stream
        end subroutine
      end interface
      end module
