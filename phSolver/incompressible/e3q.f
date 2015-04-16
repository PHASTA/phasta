        subroutine e3q (yl,      dwl,     shp,     shgl,
     &                  xl,      ql,      rmassl, 
     &                  xmudmi,  sgn )
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
c  ql     (npro,nshl,idflx) : element RHS diffusion residual 
c  rmassl     (npro,nshl)        : element lumped mass matrix
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension yl(npro,nshl,ndof),     dwl(npro,nenl),
     &            shp(nshl,ngauss),      shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),
     &            ql(npro,nshl,idflx),  rmassl(npro,nshl),
     &            xmudmi(npro,ngauss)
c
c local arrays
c
        dimension g1yi(npro,nflow),           g2yi(npro,nflow),
     &            g3yi(npro,nflow),           shg(npro,nshl,nsd),
     &            dxidx(npro,nsd,nsd),       WdetJ(npro),
     &            rmu(npro) 
c
        dimension qdi(npro,idflx),alph1(npro),alph2(npro)
c
        dimension sgn(npro,nshl),          shape(npro,nshl),
     &            shdrv(npro,nsd,nshl),    shpsum(npro)

        real*8 tmp(npro)
c
c.... for surface tension
c     
        dimension g1yti(npro),          g2yti(npro),
     &            g3yti(npro)
        integer idflow
c
c.... loop through the integration points
c
        
        
        alph1 = 0.d0
        alph2 = 0.d0
        
        do intp = 1, ngauss
        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
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

        call e3qvar   (yl,        shdrv,   
     &                 xl,           g1yi,
     &                 g2yi,      g3yi,         shg,
     &                 dxidx,     WdetJ )      
c  
        idflow = 9   ! we ALWAYS save space for tau_{ij} in q_i 
                     ! even if idiff is not greater than 1

        if(idiff >= 1) then   !so taking care of all the idiff=1,3
c
c.... compute diffusive fluxes 
c
c.... compute the viscosity
c
        call getdiff(dwl, yl, shape, xmudmi, xl, rmu, tmp)
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
        qdi(:,6)=         rmu * (g2yi(:,4) + g3yi(:,3))
        qdi(:,9)=  two * rmu *  g3yi(:,4)
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
      endif                     ! end of idiff=1 .or. 3 
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
          g1yti(:)  = g1yti(:)  + shg(:,n,1) * yl(:,n,6)
          g2yti(:)  = g2yti(:)  + shg(:,n,2) * yl(:,n,6)
          g3yti(:)  = g3yti(:)  + shg(:,n,3) * yl(:,n,6)
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
      

c
c.... return
c
       return
       end


        subroutine e3qSclr (yl,      dwl,     shp,     shgl,
     &                      xl,      ql,      rmassl, 
     &                      sgn )
c                                                                      
c----------------------------------------------------------------------
c
c This routine computes the element contribution to the 
c diffusive flux vector and the lumped mass matrix.
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension yl(npro,nshl,ndof),    dwl(npro,nshl),
     &            shp(nshl,ngauss),      shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),
     &            ql(npro,nshl,nsd),     rmassl(npro,nshl)
c
c local arrays
c
        dimension gradT(npro,nsd),     
     &            dxidx(npro,nsd,nsd),       WdetJ(npro)
c
        dimension qdi(npro,nsd),alph1(npro),alph2(npro)
c
        dimension sgn(npro,nshl),          shape(npro,nshl),
     &            shdrv(npro,nsd,nshl),    shpsum(npro)

        real*8 diffus(npro)
c
c.... loop through the integration points
c
        
        
        alph1 = 0.d0
        alph2 = 0.d0
        
        do intp = 1, ngauss
        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
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
        call e3qvarSclr   (yl,        shdrv,   
     &                     xl,        gradT,
     &                     dxidx,     WdetJ )        

c
c.... compute diffusive flux vector at this integration point
c
        call getdiffsclr(shape, dwl, yl, diffus)

c
c.... diffusive flux 
c
        qdi(:,1) =  diffus * gradT(:,1)
        qdi(:,2) =  diffus * gradT(:,2)
        qdi(:,3) =  diffus * gradT(:,3)
c
c
c.... assemble contribution of qdi to ql,i.e., contribution to 
c     each element node
c     
        do i=1,nshl
           ql(:,i,1 ) = ql(:,i,1 )+ shape(:,i)*WdetJ*qdi(:,1 )
           ql(:,i,2 ) = ql(:,i,2 )+ shape(:,i)*WdetJ*qdi(:,2 )
           ql(:,i,3 ) = ql(:,i,3 )+ shape(:,i)*WdetJ*qdi(:,3 )

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
c
c.... return
c
       return
       end

