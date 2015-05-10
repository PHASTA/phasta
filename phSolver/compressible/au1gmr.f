        subroutine Au1GMR (EGmass,   uBrg,   ilwork,iBC,iper )
c
c----------------------------------------------------------------------
c
c This routine performs a matrix-vector product for the EBE - 
c preconditioned GMRES solver.
c
c input:
c     EGmass   (numel, nedof, nedof)  : element mass matrices
c     ilwork   (nlwork)               : local MPI communication array
c     
c output:
c     uBrg   (nshg,nflow)              : next Krylov vector
c
c----------------------------------------------------------------------
c
      use pointer_data
      
      include "common.h"
      include "mpif.h"
c     
      dimension EGmass(numel,nedof,nedof),  uBrg(nshg,nflow),
     &          uBtmp(nshg,nflow),          ilwork(nlwork),
     &          iBC(nshg),
     &          iper(nshg)
c     
c.... communicate:: copy the master's portion of uBrg to each slave
c
      if (numpe > 1) then
         call commu (uBrg, ilwork, nflow  , 'out')
      endif
c
c.... local periodic boundary conditions (no communications)
c
        do j=1,nflow
           uBrg(:,j)=uBrg(iper(:),j)
        enddo
c
c       slave has masters value, for abc we need to rotate it
c        (if this is a vector only no SCALARS)
        if((iabc==1)) !are there any axisym bc's
     &     call rotabc(uBrg(1,2), iBC,  'out')

c
c.... initialize
c
      uBtmp = zero
c     
c.... loop over element blocks
c
      do iblk = 1, nelblk
         iel   = lcblk(1,iblk)
         nenl  = lcblk(5,iblk)
         npro  = lcblk(1,iblk+1) - iel
         inum = iel + npro - 1
         nshl = lcblk(10,iblk)
c
c.... compute and assemble the Au product
c
         call asAuGMR (mien(iblk)%p,  EGmass(iel:inum,:,:), uBrg,
     &                 uBtmp )
c
      enddo
      
      uBrg = uBtmp
c
c.... -------------------->   communications <-------------------------
c
c
        if((iabc==1)) !are there any axisym bc's
     &       call rotabc(uBrg(1,2), iBC,  'in ')
c
      if (numpe > 1) then
c
c.... send slave's copy of uBrg to the master
c
        call commu (uBrg  , ilwork, nflow  , 'in ')
c     
c.... nodes treated on another processor are eliminated
c     
         numtask = ilwork(1)
         itkbeg = 1

         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  uBrg(isgbeg:isgend,:) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg

         enddo
      endif
c
c.... end
c
      return
      end
c
c
c
        subroutine Au1GMRSclr (EGmasst,   uBrg,   ilwork, iper )
c
c----------------------------------------------------------------------
c
c This routine performs a matrix-vector product for the EBE - 
c preconditioned GMRES solver.
c
c input:
c     EGmasst  (numel, nshape, nshape)  : element mass matrices
c     ilwork   (nlwork)                 : local MPI communication array
c     
c output:
c     uBrg   (nshg)              : next Krylov vector
c
c----------------------------------------------------------------------
c
      use pointer_data
      
      include "common.h"
      include "mpif.h"
c     
      dimension EGmasst(numel,nshape,nshape),uBrg(nshg),
     &          uBtmp(nshg),  ilwork(nlwork), iper(nshg)
c     
c.... communicate:: copy the master's portion of uBrg to each slave
c
      if (numpe > 1) then
         call commu (uBrg, ilwork, 1, 'out')
      endif
c ... changed
c.... local periodic boundary conditions (no communications)
c
           uBrg(:)=uBrg(iper(:))
c
c
c.... initialize
c
      uBtmp = zero
c     
c.... loop over element blocks
c
      do iblk = 1, nelblk
         iel   = lcblk(1,iblk)
         nenl  = lcblk(5,iblk)
         npro  = lcblk(1,iblk+1) - iel
         inum = iel + npro - 1
         nshl = lcblk(10,iblk)
c
c.... compute and assemble the Au product
c
         call asAuGMRSclr (mien(iblk)%p,  EGmassT(iel:inum,:,:), uBrg,
     &                 uBtmp )
c
      enddo
      
      uBrg = uBtmp
c
c.... -------------------->   communications <-------------------------
c
      if (numpe > 1) then
c
c.... send slave's copy of uBrg to the master
c
        call commu (uBrg  , ilwork, 1, 'in ')
c     
c.... nodes treated on another processor are eliminated
c     
         numtask = ilwork(1)
         itkbeg = 1

         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  uBrg(isgbeg:isgend) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg

         enddo
      endif
c
c.... end
c
      return
      end


