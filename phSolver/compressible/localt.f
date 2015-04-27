        subroutine localt (global, rlocal, ien, n, code)
c
c----------------------------------------------------------------------
c
c This subroutine performs a vector gather/scatter operation. This 
c is the transpose of local.f, i.e., recieves a global and returns
c a transposed local, or the oposite.
c
c input:
c  global (nshg,n)             : global array
c  rlocal (npro,n,nenl)         : local array
c  ien    (npro,nshl)      : nodal connectivity
c  n                            : number of d.o.f.'s to be copied
c  code                         : the transfer code
c                                  .eq. 'gather  ', from global to local
c                                  .eq. 'scatter ', add  local to global 
c                                  .eq. 'globaliz', from local to global
c
c
c Zdenek Johan, Winter 1992.
c----------------------------------------------------------------------
c
        include "common.h"

        dimension global(nshg,n),           rlocal(npro,n,nshl),
     &            ien(npro,nshl)
c
        character*8 code
c
c.... ------------------------>  'localization  '  <--------------------
c
        if (code .eq. 'gather  ') then
c
c.... set timer
c
          call timer ('Gather  ')
c
c.... gather the data
c
          ttim(3) = ttim(3) - secs(0.0)

          do j = 1, nshl
            do i = 1, n
              rlocal(:,i,j) = global(ien(:,j),i)
            enddo
          enddo

	  ttim(3) = ttim(3) + secs(0.0)

c
c.... transfer count
c
c          gbytes = gbytes + n*nenl*npro
c
c.... return
c
          call timer ('Back    ')
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
c
c.... set timer
c
          call timer ('Scatter ')
c
c.... scatter the data (possible collisions)
c
          ttim(4) = ttim(4) - secs(0.0)

          do j = 1, nshl
            do i = 1, n
              do nel = 1,npro
                global(ien(nel,j),i) = global(ien(nel,j),i) 
     &                               + rlocal(nel,i,j)
              enddo
            enddo
          enddo

	  ttim(4) = ttim(4) + secs(0.0)

c
c.... transfer and flop counts
c
c          sbytes = sbytes + n*nenl*npro
c          flops  = flops  + n*nenl*npro
c
c.... return
c
          call timer ('Back    ')
          return
        endif
c
c.... ------------------------->  'globalizing '  <----------------------
c
        if (code .eq. 'globaliz') then
c
c.... scatter the data (possible collisions)
c
          do j = 1, nshl
            do i = 1, n
              do nel = 1,npro
                global(ien(nel,j),i) = rlocal(nel,i,j)
              enddo
            enddo
          enddo
c
c.... return
c
          call timer ('Back    ')
          return
        endif
c
c.... --------------------------->  error  <---------------------------
c
        call error ('local   ', code, 0)
c
c.... end
c
        end
c
c
c
        subroutine localtSclr (global, rlocal, ien, code)
c
c----------------------------------------------------------------------
c
c This subroutine performs a vector gather/scatter operation. This 
c is the transpose of local.f, i.e., recieves a global and returns
c a transposed local, or the oposite.
c
c input:
c  global (nshg)              : global array
c  rlocal (npro,nshl)         : local array
c  ien    (npro,nshl)         : nodal connectivity
c  n                            : number of d.o.f.'s to be copied
c  code                         : the transfer code
c                                  .eq. 'gather  ', from global to local
c                                  .eq. 'scatter ', add  local to global 
c                                  .eq. 'globaliz', from local to global
c
c
c Zdenek Johan, Winter 1992.
c----------------------------------------------------------------------
c
        include "common.h"

        dimension global(nshg),           rlocal(npro,nshl),
     &            ien(npro,nshl)
c
        character*8 code
c
c.... ------------------------>  'localization  '  <--------------------
c
        if (code .eq. 'gather  ') then
c
c.... set timer
c
          call timer ('Gather  ')
c
c.... gather the data
c
          ttim(3) = ttim(3) - tmr()

          do j = 1, nshl
              rlocal(:,j) = global(ien(:,j))
          enddo

	  ttim(3) = ttim(3) + tmr()

c
c.... transfer count
c
c          gbytes = gbytes + n*nshl*npro
c
c.... return
c
          call timer ('Back    ')
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
c
c.... set timer
c
          call timer ('Scatter ')
c
c.... scatter the data (possible collisions)
c
          ttim(4) = ttim(4) - tmr()

          do j = 1, nshl
             do nel = 1,npro
                global(ien(nel,j)) = global(ien(nel,j)) 
     &               + rlocal(nel,j)
             enddo
          enddo

	  ttim(4) = ttim(4) + tmr()

c
c.... transfer and flop counts
c
c          sbytes = sbytes + n*nshl*npro
c          flops  = flops  + n*nshl*npro
c
c.... return
c
          call timer ('Back    ')
          return
        endif
c
c.... ------------------------->  'globalizing '  <----------------------
c
        if (code .eq. 'globaliz') then
c
c.... scatter the data (possible collisions)
c
          do j = 1, nshl
              do nel = 1,npro
                global(ien(nel,j)) = rlocal(nel,j)
              enddo
          enddo
c
c.... return
c
          call timer ('Back    ')
          return
        endif
c
c.... --------------------------->  error  <---------------------------
c
        call error ('local   ', code, 0)
c
c.... end
c
        end
c



