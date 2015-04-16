        subroutine e3qGradV (yl,      shp,     shgl,
     &                  xl,      ql,      
     &                  sgn )
c                                                                      
c----------------------------------------------------------------------
c
c This routine computes the element contribution to the 
c diffusive flux vector and the lumped mass matrix.
c
c input: 
c  yl     (npro,nshl,ndof)       : Y variables
c  shp    (nen,ngauss)            : element shape-functions
c  shgl   (nsd,nen,ngauss)        : element local-grad-shape-functions
c  xl     (npro,nshl,nsd)        : nodal coordinates at current step
c  
c output:
c  ql     (npro,nshl,nsd*nsd) : element RHS diffusion residual 
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension yl(npro,nshl,ndof),
     &            shp(nshl,ngauss),      shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),
     &            ql(npro,nshl,nsdsq)  
c
c local arrays
c
        dimension g1yi(npro,nflow),           g2yi(npro,nflow),
     &            g3yi(npro,nflow),           shg(npro,nshl,nsd),
     &            dxidx(npro,nsd,nsd),       WdetJ(npro)
c
c
        dimension sgn(npro,nshl),          shape(npro,nshl),
     &            shdrv(npro,nsd,nshl)

c
c.... loop through the integration points
c
        
        do intp = 1, ngauss
        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
c     
        call getshp(shp,          shgl,      sgn, 
     &              shape,        shdrv)
        
c
c
c.... calculate the integration variables necessary for the
c     formation of q
c

        call e3qvar   (yl,        shdrv,   
     &                 xl,           g1yi,
     &                 g2yi,      g3yi,         shg,
     &                 dxidx,     WdetJ )      
c  
c
c     each element node
c     
        do i=1,nshl
           ql(:,i,1 ) = ql(:,i,1 )+ shape(:,i)*WdetJ*g1yi(:,2 ) ! du/dx
           ql(:,i,2 ) = ql(:,i,2 )+ shape(:,i)*WdetJ*g2yi(:,2 ) ! du/dy
           ql(:,i,3 ) = ql(:,i,3 )+ shape(:,i)*WdetJ*g3yi(:,2 ) ! du/dz

           ql(:,i,4 ) = ql(:,i,4 )+ shape(:,i)*WdetJ*g1yi(:,3 ) ! dv/dx
           ql(:,i,5 ) = ql(:,i,5 )+ shape(:,i)*WdetJ*g2yi(:,3 ) ! dv/dy
           ql(:,i,6 ) = ql(:,i,6 )+ shape(:,i)*WdetJ*g3yi(:,3 ) ! dv/dz

           ql(:,i,7 ) = ql(:,i,7 )+ shape(:,i)*WdetJ*g1yi(:,4 ) ! dw/dx
           ql(:,i,8 ) = ql(:,i,8 )+ shape(:,i)*WdetJ*g2yi(:,4 ) ! dw/dy
           ql(:,i,9 ) = ql(:,i,9 )+ shape(:,i)*WdetJ*g3yi(:,4 ) ! dw/dz

        enddo
c
c
c.... end of the loop over integration points
c
      enddo
      

c
c.... return
c
       return
       end
