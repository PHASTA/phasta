      subroutine rotabc (global, iBC, code)
c---------------------------------------------------------------------
c 
c This subroutine is responsible for rotating 
c the residual and solution vectors for axisymmetric BC's.
c
c input:   
c     global(nshg,n): global vector to be rotated.
c     code:            = 'in' for rotating with the residual
c                      = 'out' for rotating the solution 
c
c  note that the cos and sin of the rotation angles are preprocessed and
c  stored in acs(1 and 2) respectively.
c
c---------------------------------------------------------------------
c
      use specialBC  ! gives us acs, contains (:,1)=cos(theta) (:,2)=sin(theta)
      include "common.h"
 
      dimension global(nshg,2),             iBC(nshg),
     &          tmp(nshg)
 
      character*3 code

      if (code .eq. 'in ') then
         where( btest(iBC,10))
            tmp         =  global(:,1)*acs(:,1) - global(:,2)*acs(:,2)
            global(:,2) =  global(:,1)*acs(:,2) + global(:,2)*acs(:,1)
            global(:,1) = tmp
         endwhere
      else  if (code .eq. 'out') then
         where( btest(iBC,10))
            tmp         =  global(:,1)*acs(:,1) + global(:,2)*acs(:,2)
            global(:,2) = -global(:,1)*acs(:,2) + global(:,2)*acs(:,1)
            global(:,1) = tmp
         endwhere
      else 
         call error ('rotabc  ','code    ',0)
      endif

      return
      end
