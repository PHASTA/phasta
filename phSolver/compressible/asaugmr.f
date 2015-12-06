        subroutine AsAuGMR (ien,  EGmass,  uBrg, uBtmp )
c
c----------------------------------------------------------------------
c This routine computes and assembles the Au product for the 
c GMRES solver.
c
c input:
c     ien    (npro,nshl)       : nodal connectivity
c     EGmass (npro,nedof,nedof)  : element mass matrix
c     uBrg   (nshg,nflow)         : u_i before product
c
c output:
c     uBrg   (nshg,nflow)         : result of product ( u_{i+1} )
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension ien(npro,nshl),     EGmass(npro,nedof,nedof),
     &            uBrg(nshg,nflow),    uBtmp(nshg,nflow)
c     
        dimension uBrgl(npro,nedof),  ubBgl(npro,nedof)
c
c.... localize this K-vector for the EBE product
c        
        call localt (uBrg,  uBrgl,  abs(ien),  nflow,  'gather  ')
        
        ubBgl = zero
c
c.... ----------------------->  Au product  <---------------------------
c       
        do i = 1, nflow*nshl, nflow
           do j = 1, nflow*nshl, nflow
              ubBgl(:,i  ) = ubBgl(:,i  )  
     &                     + EGmass(:,i  ,j  ) * uBrgl(:,j  )
     &                     + EGmass(:,i  ,j+1) * uBrgl(:,j+1)
     &                     + EGmass(:,i  ,j+2) * uBrgl(:,j+2)
     &                     + EGmass(:,i  ,j+3) * uBrgl(:,j+3)
     &                     + EGmass(:,i  ,j+4) * uBrgl(:,j+4)
c     
              ubBgl(:,i+1) = ubBgl(:,i+1)  
     &                     + EGmass(:,i+1,j  ) * uBrgl(:,j  )
     &                     + EGmass(:,i+1,j+1) * uBrgl(:,j+1)
     &                     + EGmass(:,i+1,j+2) * uBrgl(:,j+2)
     &                     + EGmass(:,i+1,j+3) * uBrgl(:,j+3)
     &                     + EGmass(:,i+1,j+4) * uBrgl(:,j+4)
c     
              ubBgl(:,i+2) = ubBgl(:,i+2)  
     &                     + EGmass(:,i+2,j  ) * uBrgl(:,j  )
     &                     + EGmass(:,i+2,j+1) * uBrgl(:,j+1)
     &                     + EGmass(:,i+2,j+2) * uBrgl(:,j+2)
     &                     + EGmass(:,i+2,j+3) * uBrgl(:,j+3)
     &                     + EGmass(:,i+2,j+4) * uBrgl(:,j+4)
c     
              ubBgl(:,i+3) = ubBgl(:,i+3)  
     &                     + EGmass(:,i+3,j  ) * uBrgl(:,j  )
     &                     + EGmass(:,i+3,j+1) * uBrgl(:,j+1)
     &                     + EGmass(:,i+3,j+2) * uBrgl(:,j+2)
     &                     + EGmass(:,i+3,j+3) * uBrgl(:,j+3)
     &                     + EGmass(:,i+3,j+4) * uBrgl(:,j+4)
c
              ubBgl(:,i+4) = ubBgl(:,i+4)  
     &                     + EGmass(:,i+4,j  ) * uBrgl(:,j  )
     &                     + EGmass(:,i+4,j+1) * uBrgl(:,j+1)
     &                     + EGmass(:,i+4,j+2) * uBrgl(:,j+2)
     &                     + EGmass(:,i+4,j+3) * uBrgl(:,j+3)
     &                     + EGmass(:,i+4,j+4) * uBrgl(:,j+4)
c
           enddo
        enddo
c
c.... assemble the result of the product
c    
        call localt (uBtmp,  ubBgl,  abs(ien),  nflow,  'scatter ')
c     
c.... end
c
        return
        end
c
c
c
        subroutine AsAuGMRSclr (ien,  EGmass,  uBrg, uBtmp )
c
c----------------------------------------------------------------------
c This routine computes and assembles the Au product for the 
c GMRES solver.
c
c input:
c     ien    (npro,nshl)       : nodal connectivity
c     EGmass (npro,nen,nen)    : element mass matrix
c     uBrg   (nshg)           : u_i before product
c
c output:
c     uBtmp   (nshg)           : result of product ( u_{i+1} )
c
c----------------------------------------------------------------------
c
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
       dimension ien(npro,nshl),   EGmass(npro,nshape,nshape),
     &           uBrg(nshg),       uBtmp(nshg)
c     
       dimension uBrgl(npro,nshl), ubBgl(npro,nshl)
c
c.... localize this K-vector for the EBE product
c
        uBrgl = zero
        call localtSclr(uBrg,  uBrgl,  ien,  'gather  ')
        
        ubBgl = zero
c
c.... ----------------------->  Au product  <---------------------------
c       
        do i = 1, nshl
           do j = 1, nshl
              ubBgl(:,i  ) = ubBgl(:,i  )  
     &                     + EGmass(:,i  ,j  ) * uBrgl(:,j  )
c
           enddo
         enddo
c
c.... assemble the result of the product
c    
        call localtSclr(uBtmp,  ubBgl,  ien,  'scatter ')
c     
c.... end
c
        return
        end






