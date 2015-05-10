      subroutine sfID2np(isfID,nsfnp,nmap)

ccc
c     input: isfID, the surface ID of the face which user want to play with, it is specified in Simaker
c     output: nsfnp, number of the nodal points which lie on the surface with this ID, local for one process
c             nmap, a map from nsfnp to nshg, map the node on the surface to the global nodal numbering for one process
c
ccc
      use pointer_data
      include "common.h"
      include "mpif.h"
      integer isfID, nsfnp
      integer nmap(nshg)
      integer ifirstvisit(nshg)
      integer, allocatable :: ienb(:)

c--- Start to find matching surfID --- c
      ifirstvisit(:)=1
      nsfnp=0
      do iblk=1,nelblb 
         npro=lcblkb(1,iblk+1)-lcblkb(1,iblk)
         nenbl=lcblkb(6,iblk)
         nshl=lcblkb(9,iblk)
         allocate(ienb(nshl))
         do i=1,npro 

c If surfID does not match, do not consider
            iBCB2=miBCB(iblk)%p(i,2)
            if(isfID.ne.iBCB2)cycle

            ienb(1:nshl)=mienb(iblk)%p(i,1:nshl)
            do j=1,nenbl
               nn=ienb(j)
               if(ifirstvisit(nn).eq.1)then
                  ifirstvisit(nn)=0
                  nsfnp=nsfnp+1
                  nmap(nsfnp)=nn
               endif
            enddo
         enddo
         deallocate(ienb)
      enddo 
      
      return
      end 
