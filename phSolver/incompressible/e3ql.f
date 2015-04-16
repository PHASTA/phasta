      subroutine e3ql (yl,      dwl,     shp,     shgl,
     &                 xl,      ql,      xmudmi,
     &                 sgn )
c                                                                      
c----------------------------------------------------------------------
c
c This routine computes the local diffusive flux vector using a 
c local projection algorithm
c
c input: 
c  yl     (npro,nshl,ndof)       : Y variables
c  shp    (nen,ngauss)           : element shape-functions
c  shgl   (nsd,nen,ngauss)       : element local-grad-shape-functions
c  xl     (npro,nshape,nsd)      : nodal coordinates at current step
c  sgn    (npro,nshl)            : signs for reversed shape functions
c  
c output:
c  ql     (npro,nshl,nsd*nsd) : element RHS diffusion residual 
c
c----------------------------------------------------------------------
c
      use local_mass
      include "common.h"
c
      dimension yl(npro,nshl,ndof),        dwl(npro,nshl),
     &          shp(nshl,ngauss),          shgl(nsd,nshl,ngauss),
     &          xl(npro,nenl,nsd),         sgn(npro,nshl),
     &          ql(npro,nshl,idflx), xmudmi(npro,ngauss)
c
c local arrays
c
      dimension g1yi(npro,ndof),           g2yi(npro,ndof),
     &          g3yi(npro,ndof),           shg(npro,nshl,nsd),
     &          dxidx(npro,nsd,nsd),       WdetJ(npro),
     &          rmu(npro),
     &          rminv(npro,nshl,nshl),
     &          qrl(npro,nshl,nsd*nsd)
c
      dimension qdi(npro,nsd*nsd),    shape(npro,nshl),
     &          shdrv(npro,nsd,nshl),      indx(nshl),
     &          rmass(npro,nshl,nshl)


        real*8 tmp(npro)
c
c.... loop through the integration points
c
      rminv = zero
      rmass = zero
      qrl   = zero
        
      do intp = 1, ngauss

         call getshp(shp, shgl, sgn, shape, shdrv)

         qdi = zero
c
c.... calculate the integration variables 
c    
c
         call e3qvar   (yl,           shdrv,   
     &                  xl,           g1yi,
     &                  g2yi,         g3yi,         shg,
     &                  dxidx,        WdetJ )

         call getdiff(dwl,  yl, shape, xmudmi, xl,rmu, tmp)
c
c.... diffusive flux in x1-direction
c
         qdi(:,1) =  two * rmu *  g1yi(:,2)
         qdi(:,4) =        rmu * (g1yi(:,3) + g2yi(:,2))
         qdi(:,7) =        rmu * (g1yi(:,4) + g3yi(:,2))
c
c.... diffusive flux in x2-direction
c
         qdi(:,2) =        rmu * (g1yi(:,3) + g2yi(:,2))
         qdi(:,5) =  two * rmu *  g2yi(:,3)
         qdi(:,8) =        rmu * (g2yi(:,4) + g3yi(:,3))
c     
c.... diffusive flux in x3-direction
c
         qdi(:,3) =        rmu * (g1yi(:,4) + g3yi(:,2))
         qdi(:,6)=        rmu * (g2yi(:,4) + g3yi(:,3))
         qdi(:,9)=  two * rmu *  g3yi(:,4)
c
c
c.... assemble contribution of qdi to qrl,i.e., contribution to 
c     each element shape function
c
         tmp = Qwt(lcsyst,intp)
         if (lcsyst .eq. 1) then 
            tmp = tmp*(three/four)
         endif
c
c reconsider this when hierarchic wedges come into code WDGCHECK
c
        
         do i=1,nshl
            qrl(:,i,1 ) = qrl(:,i,1 )+ shape(:,i)*tmp*qdi(:,1 )
            qrl(:,i,2 ) = qrl(:,i,2 )+ shape(:,i)*tmp*qdi(:,2 )
            qrl(:,i,3 ) = qrl(:,i,3 )+ shape(:,i)*tmp*qdi(:,3 )
            
            qrl(:,i,4 ) = qrl(:,i,4 )+ shape(:,i)*tmp*qdi(:,4 )
            qrl(:,i,5 ) = qrl(:,i,5 )+ shape(:,i)*tmp*qdi(:,5 )
            qrl(:,i,6 ) = qrl(:,i,6 )+ shape(:,i)*tmp*qdi(:,6 )

            qrl(:,i,7 ) = qrl(:,i,7 )+ shape(:,i)*tmp*qdi(:,7 )
            qrl(:,i,8 ) = qrl(:,i,8 )+ shape(:,i)*tmp*qdi(:,8 )
            qrl(:,i,9 ) = qrl(:,i,9 )+ shape(:,i)*tmp*qdi(:,9 )
         enddo
c
c.... add contribution to local mass matrix
c

         if (have_local_mass .eq. 0) then
            do i=1,nshl
               do j=1,nshl
                  rmass(:,i,j) = rmass(:,i,j)+shape(:,i)*shape(:,j)*tmp
              enddo
           enddo
        endif
c
c.... end of the loop over integration points
c
      enddo

c
c.... find the inverse of the local mass matrix for each element


         if (have_local_mass .eq. 0) then
            allocate (lmassinv(iblock)%p(npro,nshl,nshl))

            do iel=1,npro
               do i=1,nshl      ! form the identy matrix
                  do j=1,nshl
                     lmassinv(iblock)%p(iel,i,j) = 0.0
                  enddo
                  lmassinv(iblock)%p(iel,i,i)=1.0
               enddo
c     
c.... LU factor the mass matrix
c
               call ludcmp(rmass(iel,:,:),nshl,nshl,indx,d)
c     
c.... back substitute with the identy matrix to find the
c     matrix inverse
c          
               do j=1,nshl
                  call lubksb(rmass(iel,:,:),nshl,nshl,indx,
     &                        lmassinv(iblock)%p(iel,:,j))
               enddo
            enddo
            rminv(:,:,:) = lmassinv(iblock)%p(:,:,:)
         else
            rminv(:,:,:) = lmassinv(iblock)%p(:,:,:)
         endif
c
c.... find the modal coefficients of ql by multiplying by the inverse of
c     the local mass matrix
c
      do iel=1,npro
        do j=1,9
c         do j=1, 3*nsd
            ql(iel,:,j) = matmul( rminv(iel,:,:),qrl(iel,:,j) )
         enddo
      enddo
c
c.... return
c
      return
      end




      subroutine e3qlSclr (yl,      dwl,     shp,     shgl,
     &                     xl,      ql,      sgn )
c                                                                      
c----------------------------------------------------------------------
c
c This routine computes the local diffusive flux vector using a 
c local projection algorithm: 
c     diffus * phi,i
c
c----------------------------------------------------------------------
c
      use local_mass
      include "common.h"
c
      dimension yl(npro,nshl,ndof),        dwl(npro,nshl),
     &          shp(nshl,ngauss),          shgl(nsd,nshl,ngauss),
     &          xl(npro,nenl,nsd),         sgn(npro,nshl),
     &          ql(npro,nshl,nsd)
c
c local arrays
c
      dimension dxidx(npro,nsd,nsd),       WdetJ(npro),
     &          diffus(npro),
     &          rminv(npro,nshl,nshl),
     &          qrl(npro,nshl,nsd)
c
      dimension qdi(npro,nsd),    shape(npro,nshl),
     &          shdrv(npro,nsd,nshl),      indx(nshl),
     &          rmass(npro,nshl,nshl),     gradT(npro,nsd),
     &          eviscv(npro)

c
c.... loop through the integration points
c
      rminv = zero
      rmass = zero
      qrl   = zero
        
      do intp = 1, ngauss

         call getshp(shp, shgl, sgn, shape, shdrv)

         qdi = zero
c
c.... calculate the integration variables 
c    
c
         call e3qvarSclr  (yl,           shdrv,        xl,           
     &                     gradT,        dxidx,        WdetJ )
c
c....  call function to sort out diffusivity (at end of this file)
c
         call getdiffsclr(dwl,shape,yl, diffus)
c
c.... diffusive flux in x1-direction
c
         qdi(:,1) =  diffus * gradT(:,1)
         qdi(:,2) =  diffus * gradT(:,2)
         qdi(:,3) =  diffus * gradT(:,3)

c
c.... assemble contribution of qdi to qrl,i.e., contribution to 
c     each element shape function
c
         tmp = Qwt(lcsyst,intp)
         if (lcsyst .eq. 1) then 
            tmp = tmp*(three/four)
         endif
        
         do i=1,nshl
            qrl(:,i,1 ) = qrl(:,i,1 )+ shape(:,i)*tmp*qdi(:,1 )
            qrl(:,i,2 ) = qrl(:,i,2 )+ shape(:,i)*tmp*qdi(:,2 )
            qrl(:,i,3 ) = qrl(:,i,3 )+ shape(:,i)*tmp*qdi(:,3 )
         enddo
c
c.... add contribution to local mass matrix
c
         if (have_local_mass .eq. 0) then
            do i=1,nshl
               do j=1,nshl
                  rmass(:,i,j)=rmass(:,i,j)+shape(:,i)*shape(:,j)*tmp
               enddo
            enddo
         endif

c.... end of the loop over integration points
c
      enddo

c
c.... find the inverse of the local mass matrix for each element
c
       qrl   = qrl/6.d0
c
c.... Assuming that lmassinv was already computed for flow equations
c     
       rmass = rmass/6.0
c
c.... for cubics, it cannot be precomputed, so compute and
c     save it the first time it is needed
c
         if (have_local_mass .eq. 0) then
            allocate (lmassinv(iblock)%p(npro,nshl,nshl))

            do iel=1,npro
               do i=1,nshl      ! form the identy matrix
                  do j=1,nshl
                     lmassinv(iblock)%p(iel,i,j) = 0.0
                  enddo
                  lmassinv(iblock)%p(iel,i,i)=1.0
               enddo
c     
c.... LU factor the mass matrix
c
               call ludcmp(rmass(iel,:,:),nshl,nshl,indx,d)
c     
c.... back substitute with the identy matrix to find the
c     matrix inverse
c          
               do j=1,nshl
                  call lubksb(rmass(iel,:,:),nshl,nshl,indx,
     &                        lmassinv(iblock)%p(iel,:,j))
               enddo
            enddo
            rminv(:,:,:) = lmassinv(iblock)%p(:,:,:)
         else
            rminv(:,:,:) = lmassinv(iblock)%p(:,:,:)
         endif
c
c.... find the modal coefficients of ql by multiplying by the inverse of
c     the local mass matrix
c
      do iel=1,npro
         do j=1,nsd
            ql(iel,:,j) = matmul( rminv(iel,:,:),qrl(iel,:,j) )
         enddo
      enddo
c
c.... return
c
      return
      end
