        subroutine genshp (shp,    shgl, nshp, nblk)  
c
c----------------------------------------------------------------------
c
c This subroutine generates shape functions for triangular,
c quadrilateral, tetrahedron, wedge and brick elements and pyramids.
c
c Christian Whiting, Winter 1999
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension shp(MAXTOP,maxsh,MAXQPT), 
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT)
c
c.... loop through element blocks
c
        maxnint=1
          do iblk = 1, nblk
c
c.... get coord. system and element type 
c
            lcsyst = lcblk(3,iblk)
            nshl   = lcblk(10,iblk)
c
c.... generate the shape-functions in local coordinates
c
            select case ( lcsyst )
            case ( 1 )          ! tets
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)  
               call shpTet(ipord,Qpt(1,1:3,i),shp(1,:,i),shgl(1,:,:,i))
            enddo
            shgl(1,:,1:nshl,1:nint(lcsyst)) = 
     &      shgl(1,:,1:nshl,1:nint(lcsyst))/two
c     
            case ( 2 )          ! hexes
c     
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)
               call shphex  (ipord, Qpt(2,1:3,i),shp(2,:,i),
     &                       shgl(2,:,:,i))
            enddo
c
            case ( 3 )          ! wedges
c
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)
               call shp6w  (ipord,Qpt(3,1:3,i),shp(3,:,i),
     &                       shgl(3,:,:,i))
            enddo

         case ( 5)              ! pyramids
            
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)
               call shppyr (ipord,Qpt(5,1:3,i),shp(5,:,i),shgl(5,:,:,i))
               
            enddo
c
c.... nonexistent element
c
            case default
c
            call error ('genshp  ', 'elem Cat', lelCat)
c
            end select
c
c.... end of generation
c
          enddo
c
c.... return
c
        return
        end
