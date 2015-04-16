        subroutine Asadj (   row_fill_list,
     &                    iens,        adjcnt   )
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        integer row_fill_list(nshg,15*nnz),
     &          ien(npro,nshl),
     &          adjcnt(nshg), ndlist(nshl)

	integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
	ien=abs(iens)

        do i=1,npro
           do j=1,nshl
              ndlist(j)=ien(i,j)
           enddo
           do j=1,nshl
              jnd=ndlist(j)
              jlngth=adjcnt(jnd) ! current length of j's list
              do k=1,nshl 
                 knd=ndlist(k)
                 ibroke=zero
                 do l= 1,jlngth
                    if(row_fill_list(jnd,l).eq. knd) then
                       ibroke=1
                       exit
                    endif
                 enddo
                 
c
c  to get here k was not in  j's list so add it
c
                 if(ibroke.eq.0) then
                    jlngth=jlngth+1 ! lenthen list
                    if(jlngth.gt.15*nnz) then
                       write(*,*) 'increase overflow factor in genadj'
                       stop
                    endif
                    row_fill_list(jnd,jlngth)=knd ! add unique entry to list
                 endif
              enddo ! finished checking all the k's for this j
              adjcnt(jnd)=jlngth  ! update the counter
           enddo                  ! done with j's
        enddo                   ! done with elements in this block
c
c
c.... end
c
        return
        end
