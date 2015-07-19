program readwrite
    use iso_c_binding
      include "common.h"
      include "mpif.h"

    integer :: rank, ierror, ione, nfiles
    type(c_ptr), TARGET :: igeom

    call MPI_Init(ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    ione = 1
    nfiles = 2
    write (*,*) rank, 'rank numfiles', nfiles, 'c_loc(igeom)', c_loc(igeom), 'igeom', igeom
    call phio_openfile_read('geombc-dat.' // char(0), nfiles, c_loc(igeom))
    write (*,*) rank, 'rank numfiles', nfiles, 'c_loc(igeom)', c_loc(igeom), 'igeom', igeom
    !      call phio_readheader(igeom,'number of nodes' // char(0),numnp,ione,
    !     & 'integer' // char(0), iotype)
    write (*,*) rank, ' calling closefile_read'

    call phio_closefile_read(igeom)


    call MPI_Finalize(ierror)
end
