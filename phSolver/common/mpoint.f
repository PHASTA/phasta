        function mpoint (name,  ndim1,  ndim2,  ndim3)
c
c----------------------------------------------------------------------
c
c This function dynamically allocates memory for the arrays.
c
c input:
c  name                 : name of the array
c  ndim1, ndim2, ndim3  : dimensions of the array
c
c output:
c  mpoint               : memory location
c
c Farzin Shakib, Summer 1985.
c Taken from dlearn (modified summer 1985)
c----------------------------------------------------------------------
c
        include "common.h"
c
        character*8 name
c
c.... calculate the array size
c
        idim1  = ndim1
        idim2  = ndim2 * min(idim1, 1)
        idim3  = ndim3 * min(idim2, 1)
c
c.... store the array information
c
        mpoint = mbeg
c
c.... set the memory pointer 
c
        mbeg   = mpoint + max(1,idim1) * max(1,idim2) * max(1,idim3)
c
c.... if past the end of total memory allocated, set the error message
c
        if (mbeg .gt. mend) call error ('mpoint  ', name, mbeg-mend)
c
c.... return
c
        return
        end
