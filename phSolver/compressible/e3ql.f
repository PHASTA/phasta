        subroutine e3ql (ycl,      shp,     shgl,
     &                   xl,      ql,      xmudmi,
     &                   sgn )
c                                                                      
c----------------------------------------------------------------------
c
c This routine computes the local diffusive flux vector using a 
c local projection algorithm
c
c input: 
c  ycl     (npro,nshl,ndof)      : Y variables
c  shp    (nen,ngauss)           : element shape-functions
c  shgl   (nsd,nen,ngauss)       : element local-grad-shape-functions
c  xl     (npro,nshape,nsd)      : nodal coordinates at current step
c  sgn    (npro,nshl)            : signs for reversed shape functions
c  
c output:
c  ql     (npro,nshape,(nflow-1)*nsd) : element RHS diffusion residual 
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension ycl(npro,nshl,ndof),  
     &            shp(nshl,ngauss),  
     &            shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),       sgn(npro,nshl),
     &            ql(npro,nshl,idflx), xmudmi(npro,ngauss)
c
c local arrays
c
        dimension g1yi(npro,nflow),          g2yi(npro,nflow),
     &            g3yi(npro,nflow),          shg(npro,nshl,nsd),
     &            dxidx(npro,nsd,nsd),       WdetJ(npro),
     &            T(npro),                   cp(npro),
     &            u1(npro),                  u2(npro),
     &            u3(npro),                  rmu(npro),
     &            rlm(npro),                 rlm2mu(npro),
     &            con(npro),                 rminv(npro,nshl,nshl),
     &            qrl(npro,nshl,(nflow-1)*nsd)
c
        dimension qdi(npro,nsd*(nflow-1)),   shape(npro,nshl),
     &            shdrv(npro,nsd,nshl),      indx(nshl),
     &            rmass(npro,nshl,nshl),     rho(npro)

        real*8 tmp(npro)
c
c.... loop through the integration points
c
        rminv = zero
        rmass = zero
        qrl = zero
        
        do intp = 1, ngauss
        if (Qwt(lcsyst,intp) .eq. zero) cycle ! precaution
c
c.... create a matrix of shape functions (and derivatives) for each
c     element at this quadrature point. These arrays will contain 
c     the correct signs for the hierarchic basis
c
        call getshp(shp,          shgl,      sgn, 
     &              shape,        shdrv)
c
c
c.... initialize
c
        qdi = zero
c
c
c.... calculate the integration variables necessary for the
c     formation of q
c

        call e3qvar   (ycl,       shape,        shdrv,   
     &                 rho,       xl,           g1yi,
     &                 g2yi,      g3yi,         shg,
     &                 dxidx,     WdetJ,        T,
     &                 cp,        u1,           u2,
     &                 u3                                 )              
c
c
c.... compute diffusive flux vector at this integration point
c
c
c.... get material properties
c
        call getDiff (T,        cp,    rho,      ycl,
     &               rmu,      rlm,    rlm2mu,   con, shape,
     &               xmudmi,   xl)
c
c.... compute diffusive fluxes 
c
c.... diffusive flux in x1-direction
c
        qdi(:,1) =  rlm2mu      * g1yi(:,2) 
     &               +      rlm * g2yi(:,3) 
     &               +      rlm * g3yi(:,4)
        qdi(:,2) =  rmu         * g1yi(:,3) 
     &               +      rmu * g2yi(:,2) 
        qdi(:,3) =  rmu         * g1yi(:,4)
     &               +      rmu * g3yi(:,2)
        qdi(:,4) =  rlm2mu * u1 * g1yi(:,2) + rmu * u2 * g1yi(:,3)
     &                                   +    rmu * u3 * g1yi(:,4)
     &               +    rmu * u2 * g2yi(:,2) + rlm * u1 * g2yi(:,3)
     &               +    rmu * u3 * g3yi(:,2) + rlm * u1 * g3yi(:,4)
     &               +    con      * g1yi(:,5)

c
c.... diffusive flux in x2-direction
c
        qdi(:,5) =       rmu * g1yi(:,3) 
     &            +      rmu * g2yi(:,2)
        qdi(:,6) =       rlm * g1yi(:,2)
     &            +   rlm2mu * g2yi(:,3)
     &            +      rlm * g3yi(:,4)
        qdi(:,7) =       rmu * g2yi(:,4)
     &            +      rmu * g3yi(:,3)
        qdi(:,8) =  rlm * u2 * g1yi(:,2) +    rmu * u1 * g1yi(:,3)
     &            + rmu * u1 * g2yi(:,2) + rlm2mu * u2 * g2yi(:,3)
     &            + rmu * u3 * g2yi(:,4)
     &            + rmu * u3 * g3yi(:,3) +    rlm * u2 * g3yi(:,4)
     &            +    con      * g2yi(:,5)
c
c.... diffusive flux in x3-direction
c
        qdi(:,9 ) =       rmu * g1yi(:,4)
     &            +      rmu * g3yi(:,2)
        qdi(:,10) =       rmu * g2yi(:,4)
     &            +      rmu * g3yi(:,3)
        qdi(:,11) =       rlm * g1yi(:,2)
     &            +      rlm * g2yi(:,3)
     &            +   rlm2mu * g3yi(:,4)
        qdi(:,12) =     rlm * u3 * g1yi(:,2) + rmu * u1 * g1yi(:,4)
     &            +    rlm * u3 * g2yi(:,3) + rmu * u2 * g2yi(:,4)
     &            +    rmu * u1 * g3yi(:,2) + rmu * u2 * g3yi(:,3)
     &            + rlm2mu * u3 * g3yi(:,4)
     &            +    con      * g3yi(:,5) 

c
c
c.... assemble contribution of qdi to qrl,i.e., contribution to 
c     each element shape function
c
        tmp = Qwt(lcsyst,intp)
        if ( lcsyst .eq. 1) then
     	   tmp = tmp*(three/four)
        end if
        
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
           qrl(:,i,10) = qrl(:,i,10)+ shape(:,i)*tmp*qdi(:,10)
           qrl(:,i,11) = qrl(:,i,11)+ shape(:,i)*tmp*qdi(:,11)
           qrl(:,i,12) = qrl(:,i,12)+ shape(:,i)*tmp*qdi(:,12)
        enddo
c
c.... add contribution to local mass matrix
c
        do i=1,nshl
           do j=1,nshl
              rmass(:,i,j) = rmass(:,i,j)+shape(:,i)*shape(:,j)*tmp
           enddo
        enddo
c
c.... end of the loop over integration points
c
       enddo
       if ( lcsyst .eq. 1) then
          qrl = qrl/6.d0
          rmass = rmass/6.0     	  
       end if
c
c.... find the inverse of the local mass matrix for each element
c
       do iel=1,npro

          do i=1,nshl   ! form the identy matrix
             do j=1,nshl
                rminv(iel,i,j) = 0.0
             enddo
             rminv(iel,i,i)=1.0
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
             call lubksb(rmass(iel,:,:),nshl,nshl,
     &                                  indx,rminv(iel,:,j))
          enddo
       enddo

c
c.... find the modal coefficients of ql by multiplying by the inverse of
c     the local mass matrix
c
       do iel=1,npro
          do j=1,12
             ql(iel,:,j) = matmul( rminv(iel,:,:),qrl(iel,:,j) )
          enddo
       enddo
c
c.... return
c
       return
       end
