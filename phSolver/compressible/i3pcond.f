      subroutine i3Pcond ( Binv,  uBrg,  ilwork,  code )
c
c---------------------------------------------------------------------
c This routine is the preconditioner driver which calls
c local routines to perform the right or left EBE preconditioning
c of a vector.
c
c input: 
c     Binv   (numel,nedof,nedof)     : element preconditioners
c     uBrg  (nshg, nflow)            : vector to be preconditioned
c     code                           : preconditioning code
c                                        .eq. 'R_Pcond ', Right precond.
c                                        .eq. 'L_Pcond ', Left precond.
c
c output:
c     uBrg   (nshg, nflow)            : preconditioned vector
c
c---------------------------------------------------------------------
c     
      use pointer_data
      
      include "common.h"
c
      dimension Binv(numel,nedof,nedof),   uBrg(nshg,nflow)
c
      dimension uBtmp(nshg,nflow), ilwork(nlwork)
c
      character*8 code
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
         inum  = iel + npro - 1
c     
c.... right precondition the vector
c     
         if (code .eq. 'R_Pcond ') then
c            
            if (iPcond .eq. 1) then
               call itrPr1 (mien(iblk)%p, Binv(iel:inum,:,:),  uBrg,  
     &                      uBtmp,        'R_Pcond ')
            endif
c     
c            if (iPcond .eq. 2) then
c               call itrPr2 (mien(iblk)%p, Binv(iel:inum,:,:),  uBrg,  
c     &              'R_Pcond ')
c            endif
         endif
c     
c.... left precondition the vector
c     
         if (code .eq. 'L_Pcond ') then
c            
            if (iPcond .eq. 1) then
               call itrPr1 (mien(iblk)%p, Binv(iel:inum,:,:),  uBrg,  
     &                      uBtmp,        'L_Pcond ')
            endif
c     
c            if (iPcond .eq. 2) then
c               call itrPr2 (mien(iblk)%p, Binv(iel:inum,:,:),  uBrg,  
c     &              'L_Pcond ')
c            endif
         endif
c
      enddo
c
c.... update the vector
c
c      if (iPcond .ne. 0) uBrg = uBtmp
c
c.... return
c
      return
      end
