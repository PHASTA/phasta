 
        subroutine D2wall (dwall,       x)
c
c----------------------------------------------------------------------
c
c This routine calculates the distance to the nearest wall for all 
c field points.
c
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use pointer_data
        include "common.h"
        include "mpif.h"
c
        dimension x(numnp,nsd),        dwall(numnp),
     &            nwall(numpe),        idisp(numpe)

c
c  Allocatable arrays
c
        allocatable xwi(:,:,:),  xw(:,:,:)

c
c.... Count the welts (wall-elements)
c
        nwalli=0
        do iblk = 1, nelblb ! loop over boundary elt blocks
           npro = lcblkb(1,iblk+1) - lcblkb(1,iblk)
           do j = 1, npro
              if(miBCB(iblk)%p(j,5).eq.2) nwalli=nwalli+1
           enddo
        enddo
c
c.... Create wallnode-coord list for welts on processor
c
        if (nwalli.eq.0) nwalli=1  !  patch for mpi's lack of imagination
        allocate (xwi(nwalli,nenb+1,nsd))
        xwi = 1.0d32
        xwi(:,nenb+1,:)=zero
        nwalli = 0
        do iblk = 1, nelblb ! loop over boundary elt blocks
c
           iel    = lcblkb(1,iblk)
           nenbl  = lcblkb(6,iblk) ! no. of vertices per bdry. face
           npro   = lcblkb(1,iblk+1) - iel 
c
           do j = 1, npro ! loop over belts in this blk
              if(miBCB(iblk)%p(j,5).eq.2) then
                 nwalli=nwalli+1
c assemble local coord list
                 do node = 1, nenbl
                    xwi(nwalli,node,1:3)=x(mienb(iblk)%p(j,node),:)
                 enddo
                 do node = 1, nenbl
                    xwi(nwalli,nenb+1,:)=xwi(nwalli,nenb+1,:)
     &                                   +xwi(nwalli,node,:)
                 enddo
                 xwi(nwalli,nenb+1,:)=xwi(nwalli,nenb+1,:)/nenbl
c
              endif
           enddo          ! loop over belts in this blk
c
        enddo               ! loop over boundary elt blocks
c
        if (nwalli.eq.0) xwi=1.0e32 ! fix for mpi's lack of imagination
        if (nwalli.eq.0) nwalli=1  !  patch for mpi's lack of imagination
c
c  Pool "number of welts" info from all processors
c
        if (numpe.gt.1) then
           call MPI_ALLGATHER(nwalli,1,MPI_INTEGER,nwall,1,
     &          MPI_INTEGER,MPI_COMM_WORLD,ierr)
           nwallt=0
           do j=1,numpe
              nwallt=nwallt+nwall(j)
           enddo
        else
           nwall=nwalli
           nwallt=nwalli
        endif
c
c  Make all-processor wallnode-coord collage
c
        allocate (xw(nwallt,nenb+1,nsd))
        if (numpe.gt.1) then
           idisp(:)=0
           do j=2,numpe
              idisp(j)=idisp(j-1)+nwall(j-1)
           enddo
           do j=1,nenb+1
              do k=1,nsd
                 call MPI_ALLGATHERV(xwi(:,j,k),nwalli,
     &                MPI_DOUBLE_PRECISION,xw(:,j,k),nwall,idisp,
     &                MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
              enddo
           enddo
        else
           xw=xwi
        endif

c
c  For each point, calculate distance to centroid; 
c  save the distance in this node's position of dwall if it's 
c  shorter than currently stored distance
c
        dwall=1.0e32
        do i=1,numnp
           do j=1, nwallt
              do k=1,nenb+1
                 distance =  ( x(i,1) - xw(j,k,1) )**2
     &                      +( x(i,2) - xw(j,k,2) )**2
     &                      +( x(i,3) - xw(j,k,3) )**2
                 if ( dwall(i).gt.distance ) dwall(i) = distance
              enddo
           enddo
        enddo
        dwall=sqrt(dwall)
c
        deallocate(xwi)
        deallocate(xw)
c
        return
        end



