      subroutine npSucPatch(nsfnp,nmap)

ccc
c     output: nsfnp, number of the nodal points which lie on the surface, local for one process
c             nmap, a map from nsfnp to nshg, map the node on the surface to the global nodal numbering for one process
c
ccc
      use pointer_data
      use wallData
      include "common.h"
      include "mpif.h"
      integer nsfnp
      integer nmap(nshg)
      

c--- Start to find surface nodes by looking for nonzero wall normals --- c
      nsfnp=0
      do ind = 1, nshg
         wnrmsq = wnorm(ind,1)**2 + wnorm(ind,2)**2 + wnorm(ind,3)**2
         if(wnrmsq.gt.pt5) then
                  nsfnp=nsfnp+1
                  nmap(nsfnp)=ind
         endif
      enddo
      
      return
      end 
