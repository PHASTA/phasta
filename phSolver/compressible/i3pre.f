      subroutine i3pre (BDtmp,  Binv,  EGmass,  ilwork )
c
c---------------------------------------------------------------------
c This is the initialization routine for the EBE-GMRES solver.
c It pre-preconditions the LHS mass matrix and sets up the 
c EBE preconditioners. (pre-preconditioning is block diagonal
c scaling). 
c
c input:
c     BDtmp  (nshg,nflow,nflow)  : block diagonal scaling matrix 
c                                  which is already LU factored.
c     EGmass (numel,nedof,nedof) : element mass matrix
c
c output:
c     EGmass (numel,nedof,nedof) : pre-preconditioned (scaled) mass matrix
c     Binv   (numel,nedof,nedof) : EBE preconditioner for each element
c
c---------------------------------------------------------------------
c
      use pointer_data

      include "common.h"
      
      dimension  BDtmp(nshg,nflow,nflow),
     &           EGmass(numel,nedof,nedof),
ctoomuchmemory     &           Binv(numel,nedof,nedof),  
     &           ilwork(nlwork)
c     
      dimension  BDiagl(numel,nshape,nflow,nflow), 
     &           BDiag(nshg,nflow,nflow)

      BDiag = BDtmp
      
      if (numpe > 1) then
        call commu (BDiag  , ilwork, nflow*nflow  , 'out')
      endif
c
c.... block-diagonal pre-precondition LHS 
c
c
c.... loop over element blocks
c
      do iblk = 1, nelblk
         iel    = lcblk(1,iblk)
         nenl   = lcblk(5,iblk)
         npro   = lcblk(1,iblk+1) - iel 
         n      = iel - 1
         inum   = iel + npro - 1
         nshl   = lcblk(10,iblk)
c     
c.... localize block diagonal scaling matrices
c     
         call local (BDiag,  BDiagl(iel:inum,:,:,:), abs(mien(iblk)%p),  
     &               nflow*nflow,  'gather  ' ) 
c     
c.... loop through local nodes and reduce all columns and rows
c     
         do inode = 1, nshl
            i = (inode - 1) * nflow ! EGmass dof offset
c     
c.... reduce by columns, (left block diagonal preconditioning)
c
c     EGmass <-- inverse(L^tilde) EGmass   
c     
            do j = 1, nedof
               do iv = 1, npro
                  EGmass(n+iv,i+1,j) = EGmass(n+iv,i+1,j)
c     
                  EGmass(n+iv,i+2,j) = (EGmass(n+iv,i+2,j) 
     &                 - BDiagl(n+iv,inode,2,1) * EGmass(n+iv,i+1,j))
c     
                  EGmass(n+iv,i+3,j) = (EGmass(n+iv,i+3,j)
     &                 - BDiagl(n+iv,inode,3,1) * EGmass(n+iv,i+1,j) 
     &                 - BDiagl(n+iv,inode,3,2) * EGmass(n+iv,i+2,j))
c     
                  EGmass(n+iv,i+4,j) = (EGmass(n+iv,i+4,j)
     &                 - BDiagl(n+iv,inode,4,1) * EGmass(n+iv,i+1,j) 
     &                 - BDiagl(n+iv,inode,4,2) * EGmass(n+iv,i+2,j) 
     &                 - BDiagl(n+iv,inode,4,3) * EGmass(n+iv,i+3,j))
c     
                  EGmass(n+iv,i+5,j) = (EGmass(n+iv,i+5,j)
     &                 - BDiagl(n+iv,inode,5,1) * EGmass(n+iv,i+1,j) 
     &                 - BDiagl(n+iv,inode,5,2) * EGmass(n+iv,i+2,j) 
     &                 - BDiagl(n+iv,inode,5,3) * EGmass(n+iv,i+3,j)
     &                 - BDiagl(n+iv,inode,5,4) * EGmass(n+iv,i+4,j))
               enddo
            enddo
         enddo
            
         do inode = 1, nshl
            i = (inode - 1) * nflow ! EGmass dof offset
            
c     
c.... reduce by rows, (right block diagonal preconditioning)
c
c     EGmass <-- EGmass inverse(U^tilde)
c     
            do j = 1, nedof
               do iv = 1, npro
                  EGmass(n+iv,j,i+1) = BDiagl(n+iv,inode,1,1) * 
     &                 (EGmass(n+iv,j,i+1))
c
                  EGmass(n+iv,j,i+2) = BDiagl(n+iv,inode,2,2) * (
     &                 EGmass(n+iv,j,i+2) 
     &                 - BDiagl(n+iv,inode,1,2) * EGmass(n+iv,j,i+1))
c
                  EGmass(n+iv,j,i+3) = BDiagl(n+iv,inode,3,3) * (
     &                 EGmass(n+iv,j,i+3)
     &                 - BDiagl(n+iv,inode,1,3) * EGmass(n+iv,j,i+1) 
     &                 - BDiagl(n+iv,inode,2,3) * EGmass(n+iv,j,i+2))
c
                  EGmass(n+iv,j,i+4) = BDiagl(n+iv,inode,4,4) * (
     &                 EGmass(n+iv,j,i+4)
     &                 - BDiagl(n+iv,inode,1,4) * EGmass(n+iv,j,i+1) 
     &                 - BDiagl(n+iv,inode,2,4) * EGmass(n+iv,j,i+2) 
     &                 - BDiagl(n+iv,inode,3,4) * EGmass(n+iv,j,i+3))
c
                  EGmass(n+iv,j,i+5) = BDiagl(n+iv,inode,5,5) * (
     &                 EGmass(n+iv,j,i+5)
     &                 - BDiagl(n+iv,inode,1,5) * EGmass(n+iv,j,i+1) 
     &                 - BDiagl(n+iv,inode,2,5) * EGmass(n+iv,j,i+2) 
     &                 - BDiagl(n+iv,inode,3,5) * EGmass(n+iv,j,i+3)
     &                 - BDiagl(n+iv,inode,4,5) * EGmass(n+iv,j,i+4))
               enddo
            enddo
c     
c.... end loops over row and column nodes
c     
         enddo
c     
c.... end of element blocks loop
c     
      enddo
c     
c.... calculate non-symmetric Cholesky EBE preconditioners
c     
ctoomuchmemory      Binv = EGmass
c$$$      if (iPcond .eq. 2) then
c$$$         call itrPr2 (ieneg, lcblk, Binv, ubBgl, ubBgl, 'LDU_Fact')
c$$$      endif
c
c.... end of (pre process), return
c     
      return

      end
c
c
c
      subroutine i3preSclr (Diag,  Dinv,  EGmassT,  ilwork )
c
c---------------------------------------------------------------------
c This is the initialization routine for the EBE-GMRES solver.
c It pre-preconditions the LHS mass matrix and sets up the 
c EBE preconditioners. (pre-preconditioning is block diagonal
c scaling). 
c
c input:
c     Diag  (numnp,nflow,nflow)    : diagonal scaling matrix 
c     EGmass (numel,nedof,nedof) : element mass matrix
c
c output:
c     EGmass (numel,nedof,nedof) : pre-preconditioned (scaled) mass matrix
c     Dinv   (numel,nedof,nedof) : EBE preconditioner for each element
c
c---------------------------------------------------------------------
c
      use pointer_data

      include "common.h"
      
      dimension  EGmassT(numel,nshape,nshape),
     &           Dinv(nshg), ilwork(nlwork)
c     
      dimension  Dinvl(numel,nshape), Diag(nshg)
c
      Dinv = one/Diag
c      
      if (numpe > 1) then
        call commu (Dinv  , ilwork, 1  , 'out')
      endif
c
c.... loop over element blocks
c
      do iblk = 1, nelblk
         iel    = lcblk(1,iblk)
         nenl   = lcblk(5,iblk)
         npro   = lcblk(1,iblk+1) - iel 
         n      = iel - 1
         inum   = iel + npro - 1
         nshl   = lcblk(10,iblk)
c
c.... localize diagonal scaling matrices
c     
         call local (Dinv,  Dinvl(iel:inum,:), mien(iblk)%p,  
     &               1,  'gather  ' ) 
c     
c.... loop through and reduce all columns
c     
         do icol = 1, nshl
            do irow = 1, nshl
               do iv = 1, npro
                  EGmassT(n+iv,irow,icol) = EGmassT(n+iv,irow,icol)
     &                                      *Dinvl(n+iv,icol)
               enddo
            enddo
         enddo
c     
c.... end of element blocks loop
c     
      enddo
c
c.... end of (pre process), return
c     
      return

      end

      

      
