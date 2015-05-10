        subroutine e3massl (bcool,    shape,  WdetJ,   A0, 
     &                      EGmass )
c
c----------------------------------------------------------------------
c This routine calculates the time contribution to the LHS tangent
c mass matrix.
c
c input:
c  shape   (npro,nshl)        : element shape functions
c  WdetJ   (npro)               : weighted Jacobian determinant
c  A0      (npro)               : weighted Jacobian matrix
c  EGmass  (npro,nedof,nedof)   : partial LHS matrix
c
c output:
c  EGmass  (npro,nedof,nedof)   : partial LHS matrix
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension  A0(npro,nflow,nflow),   shape(npro,nshl),
     &             WdetJ(npro),          EGmass(npro,nedof,nedof)
c     
        dimension  fact(npro,5),           temp(npro),
     &             shpij(npro),          bcool(npro)
c
c.... ---------------------------->  LHS  <----------------------------
c
        temp = WdetJ * almi/gami/alfi*Dtgl
        bcool = WdetJ * bcool  ! note that bcool is not used after this
                               ! (if we got in here at least).
c
c.... loop through columns (nodes j) and rows (nodes i)
c
        do j = 1, nshl
           j0 = nflow * (j - 1)
           
           do i = 1, nshl
              i0 = nflow * (i - 1)
c
c.... get the factor
c
              shpij = shape(:,i) * shape(:,j)
              fact(:,1)  = shpij * (temp + bcool*spongeContinuity)
              fact(:,2)  = shpij * (temp + bcool*spongeMomentum1)
              fact(:,3)  = shpij * (temp + bcool*spongeMomentum2)
              fact(:,4)  = shpij * (temp + bcool*spongeMomentum3)
              fact(:,5)  = shpij * (temp + bcool*spongeEnergy)
c
c.... loop through d.o.f.'s
c
              do jdof = 1, nflow
                 EGmass(:,i0+1,j0+jdof) = EGmass(:,i0+1,j0+jdof) 
     &                                  + fact(:,1) * A0(:,1,jdof)
                 EGmass(:,i0+2,j0+jdof) = EGmass(:,i0+2,j0+jdof) 
     &                                  + fact(:,2) * A0(:,2,jdof)
                 EGmass(:,i0+3,j0+jdof) = EGmass(:,i0+3,j0+jdof) 
     &                                  + fact(:,3) * A0(:,3,jdof)
                 EGmass(:,i0+4,j0+jdof) = EGmass(:,i0+4,j0+jdof) 
     &                                  + fact(:,4) * A0(:,4,jdof)
                 EGmass(:,i0+5,j0+jdof) = EGmass(:,i0+5,j0+jdof) 
     &                                  + fact(:,5) * A0(:,5,jdof)
              enddo
c
c.... end loop on row nodes
c
           enddo
c
c.... end loop on column nodes
c
        enddo
c
c.... return
c
        return
        end
c
c
c
        subroutine e3masslSclr (shape,  WdetJ,   A0t, 
     &                      EGmasst,srcp)
c
c----------------------------------------------------------------------
c This routine calculates the time contribution to the LHS tangent
c mass matrix.
c
c input:
c  shape   (npro,nshl)          : element shape functions
c  WdetJ   (npro)               : weighted Jacobian determinant
c  A0t     (npro)               : weighted Jacobian matrix
c  EGmasst (npro,nshape, nshape): partial LHS matrix
c
c output:
c  EGmasst (npro,nshape,nshape) : partial LHS matrix
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension  A0t(npro),         shape(npro,nshl),
     &             WdetJ(npro),       EGmasst(npro,nshape,nshape)
c     
        dimension  fact(npro),        temp(npro),
     &             shpij(npro),       srcp(npro)
c
c.... ---------------------------->  LHS  <----------------------------
c
        temp = WdetJ * almi/gami/alfi*Dtgl
c
c.... loop through columns (nodes j) and rows (nodes i)
c
        do j = 1, nshl
           do i = 1, nshl
c
c.... get the factor
c
              shpij = shape(:,i) * shape(:,j)
              fact  = shpij * temp
              EGmasst(:,i,j) = EGmasst(:,i,j) + fact * A0t(:)
              EGmasst(:,i,j) = EGmasst(:,i,j) 
     &                       - shpij*WdetJ(:) * srcp(:)
c
c.... end loop on row nodes
c
           enddo
c
c.... end loop on column nodes
c
        enddo
c
c.... return
c
        return
        end

