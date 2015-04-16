        subroutine gettab  (mut,  rhot,  xst)
c
c-----------------------------------------------------------------------
c
c  This subroutine reads the three tables for equilibrium chemistry.
c
c
c output:
c
c    mut  (71,451)      : specific chemical potential function of (p,T)
c    rhot (71,451)      : density function of (p,T)
c    xst  (5,71,451)    : mole fractions functions of (p,T)
c
c Note: These three arrays are always in double precision.
c 
c Frederic Chalot and Zdenek Johan, Fall 1990.
c-----------------------------------------------------------------------
c
        include "common.h"
c
        real*8  mut(71,451),  rhot(71,451),  xst(5,71,451)
c
c.... open table file
c
        open (unit=itable, file=ftable, form='unformatted',
     &                                  status='unknown')
c
c.... read tables
c
        read (itable) mut
c
        read (itable) rhot
c
        read (itable) xst
c
c.... close table file
c
        close(unit=itable)
c
c.... end
c
        return
        end
