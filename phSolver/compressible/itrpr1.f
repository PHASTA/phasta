        subroutine itrPr1 (ien,  Binv,   uBrg,    uBtmp,  code)
c
c----------------------------------------------------------------------
c
c This routine preconditions a given vector, element-by-element.
c The preconditioner used is Gauss-Siedel.
c
c input:
c  ien    (npro,nshl)         : element nodal connectivity
c  Binv   (npro,nedof,nedof)	: LHS element preconditioner matrices
c  code				: preconditioning code
c				    .eq. 'R_Pcond ', Right precond.
c				    .eq. 'L_Pcond ', Left  precond.
c
c output:
c  uBrg    (nshg,nflow)	        : preconditioned vector (uBrg) 
c
c Farzin Shakib, Winter 1987.
c----------------------------------------------------------------------
c
        include "common.h"
c
	dimension Binv(npro,nedof,nedof),  uBrg(nshg,nflow),
     &		  uBrgl(npro,nshl*nflow), ien(npro,nshl),
     &            uBtmp(nshg,nflow)
c
	character*8 code
c
c.... -------------------->  Right Pre-condition  <--------------------
c
	if (code .eq. 'R_Pcond ') then
c
c.... perform the upper triangular solve
c
	   call localt (uBrg,   uBrgl,   abs(ien),  nflow,   'gather  ' )
c
           do i = nedof-1, 1, -1
              do j = i+1, nedof
                 uBrgl(:,i) = uBrgl(:,i) - Binv(:,i,j) * uBrgl(:,j)
              enddo
           enddo
c     
	   call localt (uBrg,   uBrgl,   abs(ien),  nflow,   'globaliz')
c
	return
c
	endif
c
c.... -------------------->  Left Pre-condition  <---------------------
c
	if (code .eq. 'L_Pcond ') then
c
c.... perform the lower triangular solve (in reverse order)
c
           call localt (uBrg,   uBrgl,   abs(ien),  nflow, 'gather  ')
c
           do  i = 2, nedof
              do  j = 1, i-1
                 uBrgl(:,i) = uBrgl(:,i) - Binv(:,i,j) * uBrgl(:,j)
              enddo
           enddo
              
           call localt (uBrg,   uBrgl,   abs(ien),  nflow, 'globaliz')
c
	return
c
      endif
c
c.... error handling
c
      call error ('itrPr1  ', code, iGMRES)
c     
c.... end
c
      end
