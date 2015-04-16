        subroutine SolDir (y,          ac,        yold,      acold,
     &			   x, 
     &                     iBC,       BC,        
     &                     res,             
     &                     iper,       ilwork,
     &                     shp, shgl, shpb, shglb)
c
c----------------------------------------------------------------------
c     direct solver
c----------------------------------------------------------------------
c
        use pointer_data
        
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension y(nshg,ndof),             ac(nshg,ndof),
     &            yold(nshg,ndof),          acold(nshg,ndof),
     &            x(numnp,nsd),
     &            iBC(nshg),                BC(nshg,ndofBC),
     &            res(nshg,nflow),
     &            ilwork(nlwork),            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
        
        real*8    yAlpha(nshg,5),           acAlpha(nshg,5)

c
        dimension dyf(4*nshg),  indx(4*nshg), solinc(nshg,4)

        dimension globMas(4*nshg,4*nshg)

        write (*,*) 'Warning: using direct solver...'
c
c.... set the element matrix flag
c        
        lhs = 1   ! always
c
c.... compute solution at n+alpha
c
      call itrYAlpha( yold,  acold,  y,  ac,  yAlpha,  acAlpha)
c     
c.... *******************>> Element Data Formation <<******************
c
c.... form the LHS matrices, the residual vector 
c
        call ElmGMR (yAlpha,    acAlpha,            x,
     &               shp,       shgl,             iBC,
     &               BC,        shpb,           shglb,
     &               res,       iper,          ilwork,   
     &              rowp,       colm,            lhsK,      
     &              lhsP,       rerr   )

        globMas = zero
        npro = numel
c
cccc need to assemble here!
c
        call bc3Global(globMas,   iBC)
c
cDEBUG: write the global matrix (nonzero blocks)
c
        do i=1,nshg
           i0 = (i-1)*4
           
           do j=1,nshg
              j0 = (j-1)*4

              if (globMas(i0+1,j0+1) .ne. 0) then
                 write (544,21) i,j
                 do ii=1,4
                    write(543,20) (globMas(i0+ii,j0+kk), kk=1,4)
                 enddo
              endif
           enddo
        enddo

 20     format (4(2x,e14.7))
 21     format (2(2x,i8))
        
c$$$        stop
        
c
c.... LU factor the mass matrix
c
        indx = 0
        call ludcmp(globMas,   4*nshg,    4*nshg,  indx,  d)
        write(543,*) 'rhs'
        do i=1, nshg
           i0 = 4*(i-1)
           dyf(i0+1) = res(i,1)
           dyf(i0+2) = res(i,2)
           dyf(i0+3) = res(i,3)
           dyf(i0+4) = res(i,4)
           write(543,20) (dyf(i0+j),j=1,4)
        enddo
c
c.... back-substitute to find dY
c
        call lubksb(globMas,   4*nshg,    4*nshg,  indx,  dyf)
c     
        write(543,*) 'soln'

        do i=1,nshg
           i0 = 4*(i-1)
           solinc(i,1) = dyf(i0+1)
           solinc(i,2) = dyf(i0+2)
           solinc(i,3) = dyf(i0+3)
           solinc(i,4) = dyf(i0+4)
        enddo
c
c
c.... Now, you satisfy the boundary conditions to newly
c     obtained p,u,v,w 
c
c     
c     You have to set boundary conditions first so Dy distributes
c     
        call itrCorrect ( y,    ac,    solinc, iBC)
  	call itrBC (y,  ac,  iBC,  BC,  iper, ilwork)
c
c.... output the statistics
c
      call rstatic (res, y, solinc)
c     
c.... end
c     
      	return
        end
