        subroutine genblk (IBKSZ)
c
c----------------------------------------------------------------------
c
c  This routine reads the interior elements and generates the
c  appropriate blocks.
c
c Zdenek Johan, Fall 1991.
c----------------------------------------------------------------------
c
        use pointer_data
c
        include "common.h"
c
        integer, allocatable :: ientp(:,:)
        integer mater(ibksz)
        integer intfromfile(50) ! integers read from headers
        character*255 fname1
c
        iel=1
        itpblk=nelblk

        nelblk=0
        mattyp = 0
        ndofl = ndof
        nsymdl = nsymdf
        do iblk = 1, itpblk
c
c           read(igeom) neltp,nenl,ipordl,nshl, ijunk, ijunk, lcsyst
           iseven=7
c           call creadlist(igeom,iseven,
c     &          neltp,nenl,ipordl,nshl, ijunk, ijunk, lcsyst)
           iseven=7
           fname1='connectivity interior?'
           call readheader(igeom,fname1,intfromfile,iseven,
     &                     "integer", iotype)
           neltp  =intfromfile(1)
           nenl   =intfromfile(2)
           ipordl =intfromfile(3)
           nshl   =intfromfile(4)
           ijunk  =intfromfile(5)
           ijunk  =intfromfile(6)
           lcsyst =intfromfile(7)
           allocate (ientp(neltp,nshl))
c           read(igeom) ientp
           iientpsiz=neltp*nshl
           call readdatablock(igeom,fname1,ientp,iientpsiz,
     &                     "integer", iotype)

           do n=1,neltp,ibksz 
             
              nelblk=nelblk+1
              npro= min(IBKSZ, neltp - n + 1)
c
              lcblk(1,nelblk)  = iel
c              lcblk(2,nelblk)  = iopen ! available for later use
              lcblk(3,nelblk)  = lcsyst
              lcblk(4,nelblk)  = ipordl
              lcblk(5,nelblk)  = nenl
              lcblk(6,nelblk)  = nfacel
              lcblk(7,nelblk)  = mattyp
              lcblk(8,nelblk)  = ndofl
              lcblk(9,nelblk)  = nsymdl 
              lcblk(10,nelblk) = nshl ! # of shape functions per elt
c
c.... allocate memory for stack arrays
c
              allocate (mmat(nelblk)%p(npro))
c
              allocate (mien(nelblk)%p(npro,nshl))
              allocate (mxmudmi(nelblk)%p(npro,maxsh))
c
c.... save the element block
c
              n1=n
              n2=n+npro-1
              mater=1   ! all one material for now
              call gensav (ientp(n1:n2,1:nshl),
     &                     mater,           mien(nelblk)%p,
     &                     mmat(nelblk)%p)
              iel=iel+npro
c
           enddo
           deallocate(ientp)
        enddo
        lcblk(1,nelblk+1) = iel
c
c.... return
c
CAD        call timer ('Back    ')
c
        return
c
1000    format(a80,//,
     &  ' N o d a l   C o n n e c t i v i t y',//,
     &  '   Elem  ',/,
     &  '  Number  ',7x,27('Node',i2,:,2x))
1100    format(2x,i5,6x,27i8)
        end
