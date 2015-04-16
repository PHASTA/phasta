        subroutine input()
c
c----------------------------------------------------------------------
c This routine inputs all the necessary data, allocates required array 
c storage, and sets up the appropriate parameters for the processing.
c 
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "commonM2NFixBnd.h"
        include "mpif.h"

        external endata

        integer, allocatable :: nsons(:)
c
        character*8  date
        character*80 card

c assigned in phasta.cc
c        numpe=npe
c        myrank=mrank

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        epsM = sqrt(epsilon(one))
c
c.... read in and block all data
c
        call readnblk()
c
c....return
c
        return
        end
