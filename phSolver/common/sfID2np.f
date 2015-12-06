      subroutine sfID2np(isurfID,nsfnp,nmap)
      !Finds the nodes associated with the surface ID and returns the 
      !count in nsfnp and the nodal mapping in nmap. 
      !
      !input: 
      !  isurfID	(1)		the surface ID of the face which user want 
      !                     to play with, it is specified in Simaker
      !output: 
      ! nsfnp		(1)  	number of the nodal points which lie on the 
      !                     surface with this ID, local for one process
      ! nmap		(nshg)	a map from nsfnp to nshg; map the node on 
      !                     the surface to the global nodal numbering 
      !                     for one process

      use pointer_data
      include "common.h"
      include "mpif.h"
      integer isurfID, nsfnp
      integer nmap(nshg)
      integer ifirstvisit(nshg)
      integer, allocatable :: ienb(:)

      !Start to find matching surfID
      ifirstvisit(:)=1
      nsfnp=0
      do iblk=1,nelblb    !loop over boundary element blocks 
         npro=lcblkb(1,iblk+1)-lcblkb(1,iblk)
         nenbl=lcblkb(6,iblk)    !nenbl = number of element nodes on the boundary local?
         nshl=lcblkb(9,iblk)     !nshl = number of shape [functions] local
         allocate(ienb(nshl))
         do i=1,npro             !loop over boundary elements

         !  If surfID does not match, do not consider
            iBCB2=miBCB(iblk)%p(i,2)   
            if(isurfID.ne.iBCB2)cycle   

            ienb(1:nshl)=mienb(iblk)%p(i,1:nshl) !ienb = index of element nodes local
            do j=1,nenbl     !loop over nodes in an elmeent 
               nn=ienb(j)    !nn is a global node index. The array ienb is probably the global indeces of the nodes on the boundary for a given element
               if(ifirstvisit(nn).eq.1)then  !need this check because nodes are shared by other elements
                  ifirstvisit(nn)=0 
                  nsfnp=nsfnp+1    !nsfnp = Number of SurFace Nodal Points
                  nmap(nsfnp)=nn   !nmap is the mapping from nodal surface index to the global nodal numbering for one process. The ith index of nmap gives you the global node numbering. Thus, nn is a global node number (on a processor)
               endif
            enddo
         enddo
         deallocate(ienb)
      enddo 
      
      return
      end

      subroutine asfID2np(isurfID, nSurfNP, nodeMap)
      !Same as sfID2np, except that nmap is allocated in the subroutine
      !call to be of size nsfnp. 
      !
      !input: 
      !  isurfID	(1)		the surface ID of the face which user want to play with, it is specified in Simaker
      !output: 
      ! nSurfNP		(1)  	number of the nodal points which lie on the surface with this ID, local for one process
      ! nmap		(nshg)	a map from nsfnp to nshg, map the node on the surface to the global nodal numbering for one process
       
        include "common.h"
 
        integer :: isurfID, nSurfNP
        integer, allocatable :: nodeMap(:), tmp(:)
        
        allocate(tmp(nshg)) 
        call sfID2np(isurfID, nSurfNP, tmp)

        allocate(nodeMap(max(nSurfNP, 1)))
        nodeMap = tmp(1:nSurfNP)

!        if(nSurfNP > 0) then
!          allocate(nmap(nSurfNP))
!          nmap = tmp(1:nSurfNP)
!        else
!          allocate(nmap(1))
!        endif

        deallocate(tmp)
      end subroutine 
