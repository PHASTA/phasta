        subroutine local (global, rlocal, ientmp, n, code)
c
c----------------------------------------------------------------------
c
c This subroutine performs a vector gather/scatter operation.
c
c input:
c  global (nshg,n)             : global array
c  rlocal (npro,nshl,n)         : local array
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

        dimension global(nshg,n),           rlocal(npro,nshl,n),
     &            ien(npro,nshl),           ientmp(npro,nshl)
c
        character*8 code
        
c
c.... cubic basis has negatives in ien
c
        if (ipord > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
c
c.... ------------------------>  'localization  '  <--------------------
c
        if (code .eq. 'gather  ') then
c
c.... gather the data
c

          do j = 1, n
            do i = 1, nshl
              rlocal(:,i,j) = global(ien(:,i),j)
            enddo
          enddo


c
c.... transfer count
c
          gbytes = gbytes + n*nshl*npro
c
c.... return
c
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
c
c.... scatter the data (possible collisions)
c
          do j = 1, n
            do i = 1, nshl
              do nel = 1,npro
                global(ien(nel,i),j) = global(ien(nel,i),j) 
     &                               + rlocal(nel,i,j)
              enddo
            enddo
          enddo

c
c.... transfer and flop counts
c
          sbytes = sbytes + n*nshl*npro
          flops  = flops  + n*nshl*npro
c
c.... return
c
          return
        endif
c
c.... ------------------------->  'globalizing '  <----------------------
c
        if (code .eq. 'globaliz') then
c
c.... scatter the data (possible collisions)
c
          do j = 1, n
            do i = 1, nshl
              do nel = 1,npro
                global(ien(nel,i),j) = rlocal(nel,i,j)
              enddo
            enddo
          enddo
c
c.... return
c
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
        subroutine localx (global, rlocal, ien, n, code)
c
c----------------------------------------------------------------------
c
c This subroutine performs a vector gather/scatter operation for the
c nodal coordinates array.
c
c input:
c  global (numnp,n)             : global array
c  rlocal (npro,nenl,n)         : local array
c  ien    (npro,nshl)      : nodal connectivity
c  n                            : number of d.o.f.'s to be copied
c  code                         : the transfer code
c                                  .eq. 'gather  ', from global to local
c                                  .eq. 'scatter ', add  local to global 
c
c
c Zdenek Johan, Winter 1992.
c----------------------------------------------------------------------
c
        include "common.h"

        dimension global(numnp,n),           rlocal(npro,nenl,n),
     &            ien(npro,nshl)
c
        character*8 code
c
c.... ------------------------>  'localization  '  <--------------------
c
        if (code .eq. 'gather  ') then
c
c.... gather the data
c
          do j = 1, n
            do i = 1, nenl
              rlocal(:,i,j) = global(ien(:,i),j)
            enddo
          enddo


c
c.... transfer count
c
          gbytes = gbytes + n*nenl*npro
c
c.... return
c
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
c
c.... scatter the data (possible collisions)
c

          do j = 1, n
            do i = 1, nenl
              do nel = 1,npro
                global(ien(nel,i),j) = global(ien(nel,i),j) 
     &                               + rlocal(nel,i,j)
              enddo
            enddo
          enddo


c
c.... transfer and flop counts
c
          sbytes = sbytes + n*nenl*npro
          flops  = flops  + n*nenl*npro
c
c.... return
c
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

c        subroutine localM (global, xKebe, xGoC, ien)
cc
cc----------------------------------------------------------------------
cc This routine assembles a global tangent matrix from the element
cc matrices.
cc
cc
cc 
cc
cc
cc                         |  C      G^T |
cc           globalK   =   |             |
cc                         |  G      K   |   
cc
cc
cc
cc
cc Chris Whiting,  Winter '98
cc----------------------------------------------------------------------
cc
c        include "common.h"
c
c        dimension global(nshg*4,nshg*4),xKebe(npro,3*nshl,3*nshl),
c     &            xGoC(npro,4*nshl,nshl),
c     &            ien(npro,nshape)
cc
c        character*8 code
cc
cc.... ------------------------->  'assembling '  <----------------------
cc
c
cc     
cc.... scatter the data (possible collisions)
cc
c
cc
cc.... k
cc          
c          do iel = 1, numel
c
c             do i = 1, nshl
c                i0 = (i-1)*3
cc                
c                do j = 1, nshl
c                   j0 = (j-1)*3 
cc
c                   ia = (ien(iel,i)-1)*4 + 1
c                   ib = (ien(iel,j)-1)*4 + 1 
cc                      
c                   global(ia+1,ib+1) = global(ia+1,ib+1) 
c     &                                       + xKebe(iel,i0+1,j0+1)
c                   global(ia+1,ib+2) = global(ia+1,ib+2) 
c     &                                       + xKebe(iel,i0+1,j0+2)
c                   global(ia+1,ib+3) = global(ia+1,ib+3) 
c     &                                       + xKebe(iel,i0+1,j0+3)
c                   global(ia+2,ib+1) = global(ia+2,ib+1) 
c     &                                       + xKebe(iel,i0+2,j0+1)
c                   global(ia+2,ib+2) = global(ia+2,ib+2) 
c     &                                       + xKebe(iel,i0+2,j0+2)
c                   global(ia+2,ib+3) = global(ia+2,ib+3) 
c     &                                       + xKebe(iel,i0+2,j0+3)
c                   global(ia+3,ib+1) = global(ia+3,ib+1) 
c     &                                       + xKebe(iel,i0+3,j0+1)
c                   global(ia+3,ib+2) = global(ia+3,ib+2) 
c     &                                       + xKebe(iel,i0+3,j0+2)
c                   global(ia+3,ib+3) = global(ia+3,ib+3) 
c     &                                       + xKebe(iel,i0+3,j0+3)
cc                   
c                enddo
cc
c             enddo
cc
c          enddo
c
cc
cc.... G and G^T
cc          
c          do iel = 1, numel
c
c             do i = 1, nshl
c                i0 = (i-1)*3
c                do j = 1, nshl 
c                
c                   ia = (ien(iel,i)-1)*4 + 1
c                   ib = (ien(iel,j)-1)*4 + 1 
cc                      
c                global(ia+1,ib  ) = global(ia+1,ib  )+ xGoC(iel,i0+1,j)
c                global(ia+2,ib  ) = global(ia+2,ib  )+ xGoC(iel,i0+2,j)
c                global(ia+3,ib  ) = global(ia+3,ib  )+ xGoC(iel,i0+3,j)
c                global(ia  ,ib+1) = global(ia  ,ib+1)+ xGoC(iel,i0+1,j)
c                global(ia  ,ib+2) = global(ia  ,ib+2)+ xGoC(iel,i0+2,j)
c                global(ia  ,ib+3) = global(ia  ,ib+3)+ xGoC(iel,i0+3,j)
c
cc
c             enddo
cc
c          enddo
c       enddo
c       
cc
cc.... C
cc
c          do iel = 1, numel
c             do i = 1, nshl
c                i0 = 3*nshl + i
c                do j = 1, nshl
c                   ia = (ien(iel,i)-1)*4 + 1
c                   ib = (ien(iel,j)-1)*4 + 1
cc                      
c                   global(ia,ib) = global(ia,ib) + xGoC(iel,i0,j)
cc
c                enddo
c             enddo
c             
cc
c          enddo
c       
c          
c       
ccad	  ttim(4) = ttim(4) + secs(0.0)
c
cc
cc.... transfer and flop counts
cc
c          sbytes = sbytes + nshl*nenl*npro
c          flops  = flops  + nshl*nenl*npro
cc
cc.... return
cc
ccad          call timer ('Back    ')
c          return
cc
cc.... --------------------------->  error  <---------------------------
cc
c        call error ('local   ', code, 0)
cc
cc.... end
cc
c        end
cc
c


        subroutine localSum (global, rlocal, ientmp, nHits, n)
c
c----------------------------------------------------------------------
c
c  sum the data from the local array to the global degrees of
c  freedom and keep track of the number of locals contributing
c  to each global dof. This may be used to find the average.
c
c----------------------------------------------------------------------
c
        include "common.h"

        dimension global(nshg,n),           rlocal(npro,nshl,n),
     &            ien(npro,nshl),           ientmp(npro,nshl),
     &            nHits(nshg)
c
c.... cubic basis has negatives in ien
c
        if (ipord > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        do j = 1, n
           do i = 1, nshl
              do nel = 1,npro
                 idg = ien(nel,i)
                 global(idg,j) = global(idg,j) + rlocal(nel,i,j)
              enddo
           enddo
        enddo
        do i = 1, nshl
           do nel = 1,npro
              idg = ien(nel,i)
              nHits(idg) = nHits(idg) + 1
           enddo
        enddo
c
c.... end
c
        end
 
      subroutine localb (global, rlocal, ientmp, n, code)
c
c----------------------------------------------------------------------
c
c This subroutine performs a vector gather/scatter operation on boundary only.
c
c input:
c  global (nshg,n)             : global array
c  rlocal (npro,nshl,n)         : local array
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

        dimension global(nshg,n),           rlocal(npro,nshlb,n),
     &            ien(npro,nshl),           ientmp(npro,nshl)
c
        character*8 code
        
c
c.... cubic basis has negatives in ien
c
        if (ipord > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
c
c.... ------------------------>  'localization  '  <--------------------
c
        if (code .eq. 'gather  ') then
c
c.... set timer
c
cad          call timer ('Gather  ')
c
c.... gather the data
c

          do j = 1, n
            do i = 1, nshlb
              rlocal(:,i,j) = global(ien(:,i),j)
            enddo
          enddo


c
c.... transfer count
c
          gbytes = gbytes + n*nshl*npro
c
c.... return
c
cad          call timer ('Back    ')
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
c
c.... set timer
c
cad          call timer ('Scatter ')
c
c.... scatter the data (possible collisions)
c
          do j = 1, n
            do i = 1, nshlb
              do nel = 1,npro
                global(ien(nel,i),j) = global(ien(nel,i),j) 
     &                               + rlocal(nel,i,j)
              enddo
            enddo
          enddo

c
c.... transfer and flop counts
c
          sbytes = sbytes + n*nshlb*npro
          flops  = flops  + n*nshlb*npro
c
c.... return
c
CAD          call timer ('Back    ')
          return
        endif
c
c.... ------------------------->  'globalizing '  <----------------------
c
        if (code .eq. 'globaliz') then
c
c.... scatter the data (possible collisions)
c
          do j = 1, n
            do i = 1, nshlb
              do nel = 1,npro
                global(ien(nel,i),j) = rlocal(nel,i,j)
              enddo
            enddo
          enddo
c
c.... return
c
cad          call timer ('Back    ')
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




