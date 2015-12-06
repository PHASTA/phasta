      subroutine errsmooth(rerr,   x,     iper,   ilwork, 
     &                     shp,    shgl,  iBC)
c
        use pointer_data
c
        include "common.h"
        include "mpif.h"
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension rerrsm(nshg, 10), rerr(nshg,10), rmass(nshg),
     &            x(nshg,3)
c
        dimension ilwork(nlwork), iBC(nshg), iper(nshg)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

c
c loop over element blocks for the global reconstruction
c of the smoothed error and lumped mass matrix, rmass
c
        rerrsm = zero
        rmass = zero
        
        do iblk = 1, nelblk
c
c.... set up the parameters
c
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          ngauss = nint(lcsyst)
c
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call smooth (rerr,                x,                       
     &               tmpshp,              
     &               tmpshgl,
     &               mien(iblk)%p,
     &               rerrsm,                   
     &               rmass)

          deallocate ( tmpshp )
          deallocate ( tmpshgl ) 
       enddo
c
       if (numpe > 1) then
          call commu (rerrsm , ilwork,  10   , 'in ')
          call commu (rmass  , ilwork,  1    , 'in ')
       endif       
c
c.... take care of periodic boundary conditions
c
        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(i) = rmass(i) + rmass(j)
            rerrsm(i,:) = rerrsm(i,:) + rerrsm(j,:)
          endif
        enddo

        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(j) = rmass(i)
            rerrsm(j,:) = rerrsm(i,:)
          endif
        enddo
c
c.... invert the diagonal mass matrix and find q
c
        rmass = one/rmass
       
       do i=1, 10
          rerrsm(:,i) = rmass*rerrsm(:,i)
       enddo
       if(numpe > 1) then
          call commu (rerrsm, ilwork, 10, 'out')    
       endif
c
c      copy the smoothed error overwriting the original error.
c

       rerr = rerrsm 

       return
       end

        subroutine smooth (rerr,       x,       shp,
     &                     shgl,       ien,          
     &                     rerrsm,     rmass    )
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c interior elements for the global reconstruction of the diffusive
c flux vector.
c
c input:
c     y     (nshg,ndof)        : Y variables
c     x     (numnp,nsd)         : nodal coordinates
c     shp   (nshape,ngauss)     : element shape-functions
c     shgl  (nsd,nshape,ngauss) : element local shape-function gradients
c     ien   (npro)              : nodal connectivity array
c
c output:
c     qres  (nshg,nflow-1,nsd)  : residual vector for diffusive flux
c     rmass  (nshg)            : lumped mass matrix
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension rerr(nshg,10),               x(numnp,nsd),     
     &            shp(nshl,maxsh),  
     &            shgl(nsd,nshl,maxsh),
     &            ien(npro,nshl),
     &            rerrsm(nshg,10),    rmass(nshg)
c
c.... element level declarations
c
        dimension rerrl(npro,nshl,10),        xl(npro,nenl,nsd),         
     &            rerrsml(npro,nshl,10),       rmassl(npro,nshl)
c
        dimension sgn(npro,nshl),          shape(npro,nshl),
     &            shdrv(npro,nsd,nshl),    WdetJ(npro),
     &            dxidx(npro,nsd,nsd),     shg(npro,nshl,nsd)
c
        dimension error(npro,10)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
c
c.... gather the variables
c

        call local(rerr,   rerrl,  ien,    10,   'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
c
c.... get the element residuals 
c
        rerrsml     = zero
        rmassl      = zero

c
c.... loop through the integration points
c
        
                
        do intp = 1, ngauss
        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
c
c.... create a matrix of shape functions (and derivatives) for each
c     element at this quadrature point. These arrays will contain 
c     the correct signs for the hierarchic basis
c
        call getshp(shp,          shgl,      sgn, 
     &              shape,        shdrv)
c
        call e3metric( xl,         shdrv,        dxidx,  
     &                 shg,        WdetJ)
        error=zero
        do n = 1, nshl
           do i=1,10
              error(:,i)=error(:,i) + shape(:,n) * rerrl(:,n,i)
           enddo
        enddo
        do i=1,nshl
           do j=1,10
              rerrsml(:,i,j)  = rerrsml(:,i,j)  
     &                       + shape(:,i)*WdetJ*error(:,j)
           enddo

           rmassl(:,i) = rmassl(:,i) + shape(:,i)*WdetJ
        enddo
 
c.... end of the loop over integration points
c
      enddo
c
c.... assemble the diffusive flux residual 
c
        call local (rerrsm,   rerrsml,  ien,  10,'scatter ')
        call local (rmass,   rmassl,  ien,  1,  'scatter ')
c

      return
      end
