      subroutine genshpb (shpb,    shglb, nshpb, nblk)  
c
c----------------------------------------------------------------------
c
c This subroutine generates shape functions for triangular,
c quadrilateral, tetrahedron, wedge and brick elements.
c
c Christian Whiting, Winter 1999
c----------------------------------------------------------------------
c
      include "common.h"
c
      dimension shpb(MAXTOP,maxsh,MAXQPT), 
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT)
c
c.... loop through element blocks
c
      do iblk = 1, nblk
c
c.... get coord. system and element type 
c

         lcsyst = lcblkb(3,iblk)

c.... generate the shape-functions in local coordinates
c
         select case ( lcsyst )
         case ( 1 )             ! tets
            nshl=lcblkb(9,iblk)
            do i=1,nintb(lcsyst)  
               call shpTet(ipord,Qptb(1,1:3,i),shpb(1,:,i),
     &              shglb(1,:,:,i))
            enddo
            shglb(1,:,1:nshl,1:nintb(lcsyst)) = 
     &           shglb(1,:,1:nshl,1:nintb(lcsyst))/two
c     
         case ( 2 )             ! hexes
c     
            do i=1,nintb(lcsyst)
               call shpHex  (ipord, Qptb(2,1:3,i),shpb(2,:,i),
     &              shglb(2,:,:,i))
            enddo
c     
         case ( 3 )             ! wedges with tri bd face
            
            do i=1,nintb(lcsyst)
               call shp6w (ipord,Qptb(3,1:3,i),
     &              shpb(3,:,i),shglb(3,:,:,i))
            enddo
c     
         case ( 4 )             ! wedges with quad bd face
c     
            do i=1,nintb(lcsyst)
               call shp6w (ipord,Qptb(4,1:3,i),
     &              shpb(4,:,i),shglb(4,:,:,i))
            enddo
         case ( 5 )             ! pyramids with quad bd face
c     
            do i=1,nintb(lcsyst)
               call shppyr (ipord,Qptb(5,1:3,i),
     &              shpb(5,:,i),shglb(5,:,:,i))
            enddo
c
         case ( 6 )             ! pyramids with quad bd face
c     
            do i=1,nintb(lcsyst)
               call shppyr (ipord,Qptb(6,1:3,i),
     &              shpb(6,:,i),shglb(6,:,:,i))
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
