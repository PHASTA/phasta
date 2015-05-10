        subroutine e3q (ycl,      shp,     shgl,
     &                  xl,       ql,      rmassl, 
     &                  xmudmi,   sgn)
c                                                                      
c----------------------------------------------------------------------
c
c This routine computes the element contribution to the 
c diffusive flux vector and the lumped mass matrix.
c
c input: 
c  ycl     (npro,nshl,ndof)       : Y variables
c  shp    (nen,ngauss)            : element shape-functions
c  shgl   (nsd,nen,ngauss)        : element local-grad-shape-functions
c  xl     (npro,nshl,nsd)        : nodal coordinates at current step
c  
c output:
c  ql     (npro,nshl,idflx) : element RHS diffusion residual 
c  rmassl     (npro,nshl)        : element lumped mass matrix
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension ycl(npro,nshl,ndof),  
     &            shp(nshl,ngauss),  
     &            shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),
     &            ql(npro,nshl,idflx),  rmassl(npro,nshl),
     &            xmudmi(npro,ngauss)
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
     &            con(npro),                 rho(npro)
c
        dimension qdi(npro,idflx),alph1(npro),alph2(npro)
c
        dimension sgn(npro,nshl),          shape(npro,nshl),
     &            shdrv(npro,nsd,nshl),    shpsum(npro)
c
c.... for surface tension
c     
        dimension g1yti(npro),          g2yti(npro),
     &            g3yti(npro)
        integer idflow
c
c.... loop through the integration points
c
        
        ttim(33) = ttim(33) - secs(0.0)
        
        alph1 = 0.d0
        alph2 = 0.d0
        
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
c.... initialize
c
        qdi = zero
c
c
c.... calculate the integration variables necessary for the
c     formation of q
c

        call e3qvar   (ycl,        shape,        shdrv,   
     &                 rho,       xl,           g1yi,
     &                 g2yi,      g3yi,         shg,
     &                 dxidx,     WdetJ,        T,
     &                 cp,        u1,           u2,
     &                 u3                                 )              
c
c.... compute diffusive flux vector at this integration point
c
c
c.... get material properties
c
        call getDiff (T,        cp,       rho,        ycl,
     &                rmu,      rlm,      rlm2mu,     con,  shape,
     &                xmudmi,   xl)
          
        idflow = 0
        if(idiff >= 1) then   !so taking care of all the idiff=1,2,3
        idflow = idflow+12
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
c.... assemble contribution of qdi to ql,i.e., contribution to 
c     each element node
c
        do i=1,nshl
           ql(:,i,1 ) = ql(:,i,1 )+ shape(:,i)*WdetJ*qdi(:,1 )
           ql(:,i,2 ) = ql(:,i,2 )+ shape(:,i)*WdetJ*qdi(:,2 )
           ql(:,i,3 ) = ql(:,i,3 )+ shape(:,i)*WdetJ*qdi(:,3 )
           ql(:,i,4 ) = ql(:,i,4 )+ shape(:,i)*WdetJ*qdi(:,4 )
           ql(:,i,5 ) = ql(:,i,5 )+ shape(:,i)*WdetJ*qdi(:,5 )
           ql(:,i,6 ) = ql(:,i,6 )+ shape(:,i)*WdetJ*qdi(:,6 )
           ql(:,i,7 ) = ql(:,i,7 )+ shape(:,i)*WdetJ*qdi(:,7 )
           ql(:,i,8 ) = ql(:,i,8 )+ shape(:,i)*WdetJ*qdi(:,8 )
           ql(:,i,9 ) = ql(:,i,9 )+ shape(:,i)*WdetJ*qdi(:,9 )
           ql(:,i,10) = ql(:,i,10)+ shape(:,i)*WdetJ*qdi(:,10)
           ql(:,i,11) = ql(:,i,11)+ shape(:,i)*WdetJ*qdi(:,11)
           ql(:,i,12) = ql(:,i,12)+ shape(:,i)*WdetJ*qdi(:,12)
        enddo
c
c.... compute and assemble the element contribution to the lumped
c     mass matrix
c
c
c.... row sum technique
c
        if ( idiff == 1 ) then
           do i=1,nshl
              rmassl(:,i) = rmassl(:,i) + shape(:,i)*WdetJ
           enddo
        endif
c
c.... "special lumping technique" (Hughes p. 445)
c
        if ( idiff == 3 ) then
           shpsum = zero
           do i=1,nshl
              shpsum = shpsum + shape(:,i)*shape(:,i)
              rmassl(:,i)=rmassl(:,i)+shape(:,i)*shape(:,i)*WdetJ
           enddo
           alph1 = alph1+WdetJ
           alph2 = alph2+shpsum*WdetJ
        endif
      endif  ! end of idiff=1 .or. 3 
c
      if(isurf .eq. 1) then
c
c.... initialize
c
        g1yti   = zero
        g2yti   = zero
        g3yti   = zero
c
c.... calculate the integration variables necessary for the
c     formation of q
c
c.... compute the global gradient of Yt-variables, assuming 6th entry as 
c.... the phase indicator function 
c
c  Yt_{,x_i}=SUM_{a=1}^nshl (N_{a,x_i}(int) Yta)
c
        do n = 1, nshl
          g1yti(:)  = g1yti(:)  + shg(:,n,1) * ycl(:,n,6)
          g2yti(:)  = g2yti(:)  + shg(:,n,2) * ycl(:,n,6)
          g3yti(:)  = g3yti(:)  + shg(:,n,3) * ycl(:,n,6)
        enddo
c
c    computing N_{b}*N_{a,x_i)*yta*WdetJ
c
        do i=1,nshl
           ql(:,i,idflow+1)  = ql(:,i,idflow+1)  
     &                       + shape(:,i)*WdetJ*g1yti
           ql(:,i,idflow+2)  = ql(:,i,idflow+2)  
     &                       + shape(:,i)*WdetJ*g2yti
           ql(:,i,idflow+3)  = ql(:,i,idflow+3)  
     &                       + shape(:,i)*WdetJ*g3yti
           rmassl(:,i) = rmassl(:,i) + shape(:,i)*WdetJ
        enddo
      endif  !end of the isurf
c
c.... end of the loop over integration points
c
      enddo
c
c.... normalize the mass matrix for idiff == 3
c
      if ( idiff == 3 ) then
         do i=1,nshl
            rmassl(:,i) = rmassl(:,i)*alph1/alph2
         enddo
      endif
      
      ttim(33) = ttim(33) + secs(0.0)

c
c.... return
c
       return
       end
