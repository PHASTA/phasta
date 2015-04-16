        subroutine e3bdg (shp, 	   shg,    WdetJ,
     &		          A1, 	   A2, 	   A3, 	  
     &                    A0,      bcool,  tau,    
     &		          u1,      u2,     u3,  
     &			  BDiagl,
     &                    rmu,     rlm2mu, con)
        include "common.h"

c
c  passed arrays
c
        dimension BDiagl(npro,nshl,nflow,nflow), 
     &            shp(npro,nshl),          bcool(npro),
     &            shg(npro,nshl,nsd),      WdetJ(npro),
     &            A1tauA0(npro,nflow,nflow), A2tauA0(npro,nflow,nflow),
     &            A3tauA0(npro,nflow,nflow), Atau(npro,nflow,nflow),
     &            A1(npro,nflow,nflow),      A2(npro,nflow,nflow),
     &            A3(npro,nflow,nflow),      A0(npro,nflow,nflow),
     &            tau(npro,5),             u1(npro),
     &            u2(npro),                u3(npro),                
     &            rlm2mu(npro),
     &            rmu(npro),               con(npro)
c
c
c  passed work arrays for local variables
c
	dimension tmp1(npro),             tmp2(npro)          
	dimension tmp(npro)

	ttim(30) = ttim(30) - secs(0.0)

c $$ Ex-E3conv
c
c.... calculate the contribution in non-integrated by part form
c.... add (N_b A_tilde_i N_b,i) to BDiagl
c
       do j = 1, nshl   ! May be worth eliminating zeros in A(prim) matrices
        tmp=shp(:,j)*WdetJ
            BDiagl(:,j,1,1) = BDiagl(:,j,1,1)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,1,1) +
     &                                   shg(:,j,2) * A2(:,1,1) +
     &                                   shg(:,j,3) * A3(:,1,1) 
     &                                                            )
            BDiagl(:,j,1,2) = BDiagl(:,j,1,2)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,1,2) 
c    &                                  +shg(:,j,2) * A2(:,1,2)
c    &                                  +shg(:,j,3) * A3(:,1,2) 
     &                                                            )
            BDiagl(:,j,1,3) = BDiagl(:,j,1,3)
     &                          + tmp * (
c    &                                   shg(:,j,1) * A1(:,1,3) 
     &                                  +shg(:,j,2) * A2(:,1,3) 
c    &                                  +shg(:,j,3) * A3(:,1,3) 
     &                                                            )
            BDiagl(:,j,1,4) = BDiagl(:,j,1,4)
     &                          + tmp * (
c    &                                   shg(:,j,1) * A1(:,1,4) +
c    &                                   shg(:,j,2) * A2(:,1,4) +
     &                                   shg(:,j,3) * A3(:,1,4) 
     &                                                            )
            BDiagl(:,j,1,5) = BDiagl(:,j,1,5)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,1,5) +
     &                                   shg(:,j,2) * A2(:,1,5) +
     &                                   shg(:,j,3) * A3(:,1,5) 
     &                                                            )
            BDiagl(:,j,2,1) = BDiagl(:,j,2,1)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,2,1) +
     &                                   shg(:,j,2) * A2(:,2,1) +
     &                                   shg(:,j,3) * A3(:,2,1) 
     &                                                            )
            BDiagl(:,j,2,2) = BDiagl(:,j,2,2)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,2,2) +
     &                                   shg(:,j,2) * A2(:,2,2) +
     &                                   shg(:,j,3) * A3(:,2,2) 
     &                                                            )
            BDiagl(:,j,2,3) = BDiagl(:,j,2,3)
     &                          + tmp * (
c    &                                   shg(:,j,1) * A1(:,2,3) 
     &                                  +shg(:,j,2) * A2(:,2,3) 
c    &                                  +shg(:,j,3) * A3(:,2,3) 
     &                                                            )
            BDiagl(:,j,2,4) = BDiagl(:,j,2,4)
     &                          + tmp * (
c    &                                   shg(:,j,1) * A1(:,2,4) +
c    &                                   shg(:,j,2) * A2(:,2,4) +
     &                                   shg(:,j,3) * A3(:,2,4) 
     &                                                            )
            BDiagl(:,j,2,5) = BDiagl(:,j,2,5)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,2,5) +
     &                                   shg(:,j,2) * A2(:,2,5) +
     &                                   shg(:,j,3) * A3(:,2,5) 
     &                                                            )
            BDiagl(:,j,3,1) = BDiagl(:,j,3,1)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,3,1) +
     &                                   shg(:,j,2) * A2(:,3,1) +
     &                                   shg(:,j,3) * A3(:,3,1) 
     &                                                            )
            BDiagl(:,j,3,2) = BDiagl(:,j,3,2)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,3,2) 
c    &                                  +shg(:,j,2) * A2(:,3,2) 
c    &                                  +shg(:,j,3) * A3(:,3,2) 
     &                                                            )
            BDiagl(:,j,3,3) = BDiagl(:,j,3,3)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,3,3) +
     &                                   shg(:,j,2) * A2(:,3,3) +
     &                                   shg(:,j,3) * A3(:,3,3) 
     &                                                            )
            BDiagl(:,j,3,4) = BDiagl(:,j,3,4)
     &                          + tmp * (
c    &                                   shg(:,j,1) * A1(:,3,4) +
c    &                                   shg(:,j,2) * A2(:,3,4) +
     &                                   shg(:,j,3) * A3(:,3,4) 
     &                                                            )
            BDiagl(:,j,3,5) = BDiagl(:,j,3,5)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,3,5) +
     &                                   shg(:,j,2) * A2(:,3,5) +
     &                                   shg(:,j,3) * A3(:,3,5) 
     &                                                            )
            BDiagl(:,j,4,1) = BDiagl(:,j,4,1)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,4,1) +
     &                                   shg(:,j,2) * A2(:,4,1) +
     &                                   shg(:,j,3) * A3(:,4,1) 
     &                                                            )
            BDiagl(:,j,4,2) = BDiagl(:,j,4,2)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,4,2) 
c    &                                  +shg(:,j,2) * A2(:,4,2) 
c    &                                  +shg(:,j,3) * A3(:,4,2) 
     &                                                            )
            BDiagl(:,j,4,3) = BDiagl(:,j,4,3)
     &                          + tmp * (
c    &                                   shg(:,j,1) * A1(:,4,3) 
     &                                  +shg(:,j,2) * A2(:,4,3) 
c    &                                  +shg(:,j,3) * A3(:,4,3) 
     &                                                            )
            BDiagl(:,j,4,4) = BDiagl(:,j,4,4)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,4,4) +
     &                                   shg(:,j,2) * A2(:,4,4) +
     &                                   shg(:,j,3) * A3(:,4,4) 
     &                                                            )
            BDiagl(:,j,4,5) = BDiagl(:,j,4,5)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,4,5) +
     &                                   shg(:,j,2) * A2(:,4,5) +
     &                                   shg(:,j,3) * A3(:,4,5) 
     &                                                            )
            BDiagl(:,j,5,1) = BDiagl(:,j,5,1)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,5,1) +
     &                                   shg(:,j,2) * A2(:,5,1) +
     &                                   shg(:,j,3) * A3(:,5,1) 
     &                                                            )
            BDiagl(:,j,5,2) = BDiagl(:,j,5,2)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,5,2) +
     &                                   shg(:,j,2) * A2(:,5,2) +
     &                                   shg(:,j,3) * A3(:,5,2) 
     &                                                            )
            BDiagl(:,j,5,3) = BDiagl(:,j,5,3)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,5,3) +
     &                                   shg(:,j,2) * A2(:,5,3) +
     &                                   shg(:,j,3) * A3(:,5,3) 
     &                                                            )
            BDiagl(:,j,5,4) = BDiagl(:,j,5,4)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,5,4) +
     &                                   shg(:,j,2) * A2(:,5,4) +
     &                                   shg(:,j,3) * A3(:,5,4) 
     &                                                            )
            BDiagl(:,j,5,5) = BDiagl(:,j,5,5)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1(:,5,5) +
     &                                   shg(:,j,2) * A2(:,5,5) +
     &                                   shg(:,j,3) * A3(:,5,5) 
     &                                                            )
       enddo
        if (ngauss .eq. 1 .and. nshl .eq. 4) then  ! Exact integ of mass term for tets

c $$ Ex-e3mael
c
c.... calculate factor   V             * sum(column of M)
c                     =WdetJ/(3/4)Qwt(lcsyst,intp) * (2+1+1+1)/20
c
	tmp =  0.4d0 * WdetJ * Dtgl/Qwt(lcsyst,intp)/three*almi/gami/alfi
c
        do j=1,nshl  ! take advantage of zeros in A0(Prim)
	  BDiagl(:,j,2,2) = BDiagl(:,j,2,2) + tmp * A0(:,2,2)
	  BDiagl(:,j,3,3) = BDiagl(:,j,3,3) + tmp * A0(:,3,3)
	  BDiagl(:,j,4,4) = BDiagl(:,j,4,4) + tmp * A0(:,4,4)
	  BDiagl(:,j,1,5) = BDiagl(:,j,1,5) + tmp * A0(:,1,5)
	  BDiagl(:,j,2,5) = BDiagl(:,j,2,5) + tmp * A0(:,2,5)
	  BDiagl(:,j,3,5) = BDiagl(:,j,3,5) + tmp * A0(:,3,5)
	  BDiagl(:,j,4,5) = BDiagl(:,j,4,5) + tmp * A0(:,4,5)
	  BDiagl(:,j,1,1) = BDiagl(:,j,1,1) + tmp * A0(:,1,1)
	  BDiagl(:,j,2,1) = BDiagl(:,j,2,1) + tmp * A0(:,2,1)
	  BDiagl(:,j,3,1) = BDiagl(:,j,3,1) + tmp * A0(:,3,1)
	  BDiagl(:,j,4,1) = BDiagl(:,j,4,1) + tmp * A0(:,4,1)
	  BDiagl(:,j,5,1) = BDiagl(:,j,5,1) + tmp * A0(:,5,1)
	  BDiagl(:,j,5,2) = BDiagl(:,j,5,2) + tmp * A0(:,5,2)
	  BDiagl(:,j,5,3) = BDiagl(:,j,5,3) + tmp * A0(:,5,3)
	  BDiagl(:,j,5,4) = BDiagl(:,j,5,4) + tmp * A0(:,5,4)
	  BDiagl(:,j,5,5) = BDiagl(:,j,5,5) + tmp * A0(:,5,5)
	enddo

        else

c $$ Ex-e3mass
c
        ff = almi / gami / alfi
        tmp = WdetJ * (Dtgl * ff + bcool)
c
        do j = 1, nshl
c
          tmp2 = (shp(:,j)*shp(:,j)) * tmp   
c
	  BDiagl(:,j,2,2) = BDiagl(:,j,2,2) + tmp2 * A0(:,2,2)
	  BDiagl(:,j,3,3) = BDiagl(:,j,3,3) + tmp2 * A0(:,3,3)
	  BDiagl(:,j,4,4) = BDiagl(:,j,4,4) + tmp2 * A0(:,4,4)
	  BDiagl(:,j,1,5) = BDiagl(:,j,1,5) + tmp2 * A0(:,1,5)
	  BDiagl(:,j,2,5) = BDiagl(:,j,2,5) + tmp2 * A0(:,2,5)
	  BDiagl(:,j,3,5) = BDiagl(:,j,3,5) + tmp2 * A0(:,3,5)
	  BDiagl(:,j,4,5) = BDiagl(:,j,4,5) + tmp2 * A0(:,4,5)
	  BDiagl(:,j,1,1) = BDiagl(:,j,1,1) + tmp2 * A0(:,1,1)
	  BDiagl(:,j,2,1) = BDiagl(:,j,2,1) + tmp2 * A0(:,2,1)
	  BDiagl(:,j,3,1) = BDiagl(:,j,3,1) + tmp2 * A0(:,3,1)
	  BDiagl(:,j,4,1) = BDiagl(:,j,4,1) + tmp2 * A0(:,4,1)
	  BDiagl(:,j,5,1) = BDiagl(:,j,5,1) + tmp2 * A0(:,5,1)
	  BDiagl(:,j,5,2) = BDiagl(:,j,5,2) + tmp2 * A0(:,5,2)
	  BDiagl(:,j,5,3) = BDiagl(:,j,5,3) + tmp2 * A0(:,5,3)
	  BDiagl(:,j,5,4) = BDiagl(:,j,5,4) + tmp2 * A0(:,5,4)
	  BDiagl(:,j,5,5) = BDiagl(:,j,5,5) + tmp2 * A0(:,5,5)
c
        enddo

        endif

c
c.... calculate (Atau) <-- (A_1 tau) (Recall that we are using a 
c                                     diagonal tau here)          
c     
       if (ivart .ge. 2) then
       do i = 1, nflow
          Atau(:,i,1) = A1(:,i,1)*tau(:,1)
          Atau(:,i,2) = A1(:,i,2)*tau(:,2)
          Atau(:,i,3) = A1(:,i,3)*tau(:,2)
          Atau(:,i,4) = A1(:,i,4)*tau(:,2)
          Atau(:,i,5) = A1(:,i,5)*tau(:,3)
       enddo
c
c.... calculate (A_1 tau A_0) (for L.S. time term of EGmass)
c
       do j = 1, nflow
          do i = 1, nflow
             A1tauA0(:,i,j) = 
     &            Atau(:,i,1)*A0(:,1,j) +
     &            Atau(:,i,2)*A0(:,2,j) +
     &            Atau(:,i,3)*A0(:,3,j) +
     &            Atau(:,i,4)*A0(:,4,j) +
     &            Atau(:,i,5)*A0(:,5,j)
          enddo
       enddo
c
c.... calculate (Atau) <-- (A_2 tau) (Recall that we are using a 
c                                     diagonal tau here)          
c     
       do i = 1, nflow
          Atau(:,i,1) = A2(:,i,1)*tau(:,1)
          Atau(:,i,2) = A2(:,i,2)*tau(:,2)
          Atau(:,i,3) = A2(:,i,3)*tau(:,2)
          Atau(:,i,4) = A2(:,i,4)*tau(:,2)
          Atau(:,i,5) = A2(:,i,5)*tau(:,3)
       enddo
c
c.... calculate (A_2 tau A_0) (for L.S. time term of EGmass)
c
       do j = 1, nflow
          do i = 1, nflow
             A2tauA0(:,i,j) = 
     &            Atau(:,i,1)*A0(:,1,j) +
     &            Atau(:,i,2)*A0(:,2,j) +
     &            Atau(:,i,3)*A0(:,3,j) +
     &            Atau(:,i,4)*A0(:,4,j) +
     &            Atau(:,i,5)*A0(:,5,j)
          enddo
       enddo
c     
       do i = 1, nflow
          Atau(:,i,1) = A3(:,i,1)*tau(:,1)
          Atau(:,i,2) = A3(:,i,2)*tau(:,2)
          Atau(:,i,3) = A3(:,i,3)*tau(:,2)
          Atau(:,i,4) = A3(:,i,4)*tau(:,2)
          Atau(:,i,5) = A3(:,i,5)*tau(:,3)
       enddo
c
c.... calculate (A_3 tau A_0) (for L.S. time term of EGmass)
c
       do j = 1, nflow
          do i = 1, nflow
             A3tauA0(:,i,j) =
     &            Atau(:,i,1)*A0(:,1,j) +
     &            Atau(:,i,2)*A0(:,2,j) +
     &            Atau(:,i,3)*A0(:,3,j) +
     &            Atau(:,i,4)*A0(:,4,j) +
     &            Atau(:,i,5)*A0(:,5,j)
          enddo
       enddo
c
c
c.... add least squares time term to BDiagl
c
c
c.... loop through rows (nodes i)
c
       do i = 1, nshl
c
c.... loop through column nodes, add (N_a,i A_i tau N_b) to BDiagl
c
          tmp1 = shp(:,i) * WdetJ * almi/gami/alfi*dtgl
          do idof = 1, nflow
             do jdof = 1, nflow
                BDiagl(:,i,idof,jdof) = BDiagl(:,i,idof,jdof)
     &                                + tmp1*(
     &               shg(:,i,1) * A1tauA0(:,idof,jdof) +
     &               shg(:,i,2) * A2tauA0(:,idof,jdof) +
     &               shg(:,i,3) * A3tauA0(:,idof,jdof))
             enddo
          enddo
       enddo
       endif
c
c.... add contribution of stiffness to BDiagl
c
c....  we no longer build stiff so we have to get it in here.
c      recall that this is the (A_i tau A_j + K_{ij}) that
c      multiplies N_{a,i} N_{b,j} (Actually for BDiag a=b).
c
c...   we are through with A0 so use it now for the sub-blocks
c      of what used to be stiff
c
c
     
c
c  start with 1 1 contribution
c
c
c  Note that I comment out lines where zeros appear (for Pvariables)
c 
          if(ivart .ge. 2) then 
c
c  row one
c
          A0(:,1,1)= 
     &             tau(:,1) *    A1(:,1,1) * A1(:,1,1)
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A1(:,2,1)
c    &                         + A1(:,1,3) * A1(:,3,1)
c    &                         + A1(:,1,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A1(:,1,5) * A1(:,5,1)
c
          A0(:,1,2)= 
     &             tau(:,1) *    A1(:,1,1) * A1(:,1,2)
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A1(:,2,2)
c    &                         + A1(:,1,3) * A1(:,3,2)
c    &                         + A1(:,1,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) *    A1(:,1,5) * A1(:,5,2)
c
          A0(:,1,3)= 
c    &             tau(:,1) *    A1(:,1,1) * A1(:,1,3)
c    &           + tau(:,2) * ( 
c    &                         + A1(:,1,2) * A1(:,2,3)
c    &                         + A1(:,1,3) * A1(:,3,3)
c    &                         + A1(:,1,4) * A1(:,4,3)
c    &                                                  )
     &           + tau(:,3) *    A1(:,1,5) * A1(:,5,3)
c
          A0(:,1,4)= 
c    &             tau(:,1) *    A1(:,1,1) * A1(:,1,4)
c    &           + tau(:,2) * ( 
c    &                         + A1(:,1,2) * A1(:,2,4)
c    &                         + A1(:,1,3) * A1(:,3,4)
c    &                         + A1(:,1,4) * A1(:,4,4)
c    &                                                  )
     &           + tau(:,3) *    A1(:,1,5) * A1(:,5,4)
c
          A0(:,1,5)= 
     &             tau(:,1) *    A1(:,1,1) * A1(:,1,5)
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A1(:,2,5)
c    &                         + A1(:,1,3) * A1(:,3,5)
c    &                         + A1(:,1,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A1(:,1,5) * A1(:,5,5)
c
c  row two
c
          A0(:,2,1)= 
     &             tau(:,1) *    A1(:,2,1) * A1(:,1,1)
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A1(:,2,1)
c    &                         + A1(:,2,3) * A1(:,3,1)
c    &                         + A1(:,2,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A1(:,2,5) * A1(:,5,1)
c
          A0(:,2,2)= 
     &             tau(:,1) *    A1(:,2,1) * A1(:,1,2)
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A1(:,2,2)
c    &                         + A1(:,2,3) * A1(:,3,2)
c    &                         + A1(:,2,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) *    A1(:,2,5) * A1(:,5,2)
c
          A0(:,2,3)= 
c    &             tau(:,1) *    A1(:,2,1) * A1(:,1,3)
c    &           + tau(:,2) * ( 
c    &                         + A1(:,2,2) * A1(:,2,3)
c    &                         + A1(:,2,3) * A1(:,3,3)
c    &                         + A1(:,2,4) * A1(:,4,3)
c    &                                                  )
     &           + tau(:,3) *    A1(:,2,5) * A1(:,5,3)
c
          A0(:,2,4)= 
c    &             tau(:,1) *    A1(:,2,1) * A1(:,1,4)
c    &           + tau(:,2) * ( 
c    &                         + A1(:,2,2) * A1(:,2,4)
c    &                         + A1(:,2,3) * A1(:,3,4)
c    &                         + A1(:,2,4) * A1(:,4,4)
c    &                                                  )
     &           + tau(:,3) *    A1(:,2,5) * A1(:,5,4)
c
          A0(:,2,5)= 
     &             tau(:,1) *    A1(:,2,1) * A1(:,1,5)
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A1(:,2,5)
c    &                         + A1(:,2,3) * A1(:,3,5)
c    &                         + A1(:,2,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A1(:,2,5) * A1(:,5,5)
c
c  row three
c
          A0(:,3,1)= 
     &             tau(:,1) *    A1(:,3,1) * A1(:,1,1)
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A1(:,2,1)
     &                         + A1(:,3,3) * A1(:,3,1)
c    &                         + A1(:,3,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A1(:,3,5) * A1(:,5,1)
c
          A0(:,3,2)= 
     &             tau(:,1) *    A1(:,3,1) * A1(:,1,2)
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A1(:,2,2)
     &                         + A1(:,3,3) * A1(:,3,2)
c    &                         + A1(:,3,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) *    A1(:,3,5) * A1(:,5,2)
c
          A0(:,3,3)= 
c    &             tau(:,1) *    A1(:,3,1) * A1(:,1,3)
     &           + tau(:,2) * ( 
c    &                         + A1(:,3,2) * A1(:,2,3)
     &                         + A1(:,3,3) * A1(:,3,3)
c    &                         + A1(:,3,4) * A1(:,4,3)
     &                                                  )
     &           + tau(:,3) *    A1(:,3,5) * A1(:,5,3)
c
          A0(:,3,4)= 
c    &             tau(:,1) *    A1(:,3,1) * A1(:,1,4)
c    &           + tau(:,2) * ( 
c    &                         + A1(:,3,2) * A1(:,2,4)
c    &                         + A1(:,3,3) * A1(:,3,4)
c    &                         + A1(:,3,4) * A1(:,4,4)
c    &                                                  )
     &           + tau(:,3) *    A1(:,3,5) * A1(:,5,4)
c
          A0(:,3,5)= 
     &             tau(:,1) *    A1(:,3,1) * A1(:,1,5)
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A1(:,2,5)
     &                         + A1(:,3,3) * A1(:,3,5)
c    &                         + A1(:,3,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A1(:,3,5) * A1(:,5,5)
c
c  row four
c
          A0(:,4,1)= 
     &             tau(:,1) *    A1(:,4,1) * A1(:,1,1)
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A1(:,2,1)
c    &                         + A1(:,4,3) * A1(:,3,1)
     &                         + A1(:,4,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A1(:,4,5) * A1(:,5,1)
c
          A0(:,4,2)= 
     &             tau(:,1) *    A1(:,4,1) * A1(:,1,2)
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A1(:,2,2)
c    &                         + A1(:,4,3) * A1(:,3,2)
     &                         + A1(:,4,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) *    A1(:,4,5) * A1(:,5,2)
c
          A0(:,4,3)= 
c    &             tau(:,1) *    A1(:,4,1) * A1(:,1,3)
c    &           + tau(:,2) * ( 
c    &                         + A1(:,4,2) * A1(:,2,3)
c    &                         + A1(:,4,3) * A1(:,3,3)
c    &                         + A1(:,4,4) * A1(:,4,3)
c    &                                                  )
     &           + tau(:,3) *    A1(:,4,5) * A1(:,5,3)
c
          A0(:,4,4)= 
c    &             tau(:,1) *    A1(:,4,1) * A1(:,1,4)
     &           + tau(:,2) * ( 
c    &                         + A1(:,4,2) * A1(:,2,4)
c    &                         + A1(:,4,3) * A1(:,3,4)
     &                         + A1(:,4,4) * A1(:,4,4)
     &                                                  )
     &           + tau(:,3) *    A1(:,4,5) * A1(:,5,4)
c
          A0(:,4,5)= 
     &             tau(:,1) *    A1(:,4,1) * A1(:,1,5)
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A1(:,2,5)
c    &                         + A1(:,4,3) * A1(:,3,5)
     &                         + A1(:,4,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A1(:,4,5) * A1(:,5,5)
c
c  row two
c
          A0(:,5,1)= 
     &             tau(:,1) *    A1(:,5,1) * A1(:,1,1)
     &           + tau(:,2) * ( 
     &                         + A1(:,5,2) * A1(:,2,1)
     &                         + A1(:,5,3) * A1(:,3,1)
     &                         + A1(:,5,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A1(:,5,5) * A1(:,5,1)
c
          A0(:,5,2)= 
     &             tau(:,1) *    A1(:,5,1) * A1(:,1,2)
     &           + tau(:,2) * ( 
     &                         + A1(:,5,2) * A1(:,2,2)
     &                         + A1(:,5,3) * A1(:,3,2)
     &                         + A1(:,5,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) *    A1(:,5,5) * A1(:,5,2)
c
          A0(:,5,3)= 
c    &             tau(:,1) *    A1(:,5,1) * A1(:,1,3)
     &           + tau(:,2) * ( 
c    &                         + A1(:,5,2) * A1(:,2,3)
     &                         + A1(:,5,3) * A1(:,3,3)
c    &                         + A1(:,5,4) * A1(:,4,3)
     &                                                  )
     &           + tau(:,3) *    A1(:,5,5) * A1(:,5,3)
c
          A0(:,5,4)= 
c    &             tau(:,1) *    A1(:,5,1) * A1(:,1,4)
     &           + tau(:,2) * ( 
c    &                         + A1(:,5,2) * A1(:,2,4)
c    &                         + A1(:,5,3) * A1(:,3,4)
     &                         + A1(:,5,4) * A1(:,4,4)
     &                                                  )
     &           + tau(:,3) *    A1(:,5,5) * A1(:,5,4)
c
          A0(:,5,5)= 
     &             tau(:,1) *    A1(:,5,1) * A1(:,1,5)
     &           + tau(:,2) * ( 
     &                         + A1(:,5,2) * A1(:,2,5)
     &                         + A1(:,5,3) * A1(:,3,5)
     &                         + A1(:,5,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A1(:,5,5) * A1(:,5,5)
c
         else
          A0=zero
         endif
c
         if (Navier .eq. 1) then
           A0(:,2,2) = A0(:,2,2) + rlm2mu !rlm2mu
           A0(:,3,3) = A0(:,3,3) + rmu !rmu
           A0(:,4,4) = A0(:,4,4) + rmu !rmu
           A0(:,5,2) = A0(:,5,2) + rlm2mu * u1
           A0(:,5,3) = A0(:,5,3) + rmu * u2
           A0(:,5,4) = A0(:,5,4) + rmu * u3
           A0(:,5,5) = A0(:,5,5) + con !con
        endif
        do j = 1, nshl
          tmp = WdetJ * shg(:,j,1) * shg(:,j,1)
           do i=1,nflow
           do k=1,nflow
             BDiagl(:,j,i,k) = BDiagl(:,j,i,k) + tmp * A0(:,i,k)
           enddo
           enddo
        enddo
c
c  start with 2 2 contribution
c
          if(ivart.ge.2) then
c
c  row one
c
          A0(:,1,1)= 
     &             tau(:,1) *    A2(:,1,1) * A2(:,1,1)
     &           + tau(:,2) * ( 
c    &                         + A2(:,1,2) * A2(:,2,1)
     &                         + A2(:,1,3) * A2(:,3,1)
c    &                         + A2(:,1,4) * A2(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A2(:,1,5) * A2(:,5,1)
c
          A0(:,1,2)= 
c    &             tau(:,1) *    A2(:,1,1) * A2(:,1,2)
c    &           + tau(:,2) * ( 
c    &                         + A2(:,1,2) * A2(:,2,2)
c    &                         + A2(:,1,3) * A2(:,3,2)
c    &                         + A2(:,1,4) * A2(:,4,2)
c    &                                                  )
     &           + tau(:,3) *    A2(:,1,5) * A2(:,5,2)
c
          A0(:,1,3)= 
     &             tau(:,1) *    A2(:,1,1) * A2(:,1,3)
     &           + tau(:,2) * ( 
c    &                         + A2(:,1,2) * A2(:,2,3)
     &                         + A2(:,1,3) * A2(:,3,3)
c    &                         + A2(:,1,4) * A2(:,4,3)
     &                                                  )
     &           + tau(:,3) *    A2(:,1,5) * A2(:,5,3)
c
          A0(:,1,4)= 
c    &             tau(:,1) *    A2(:,1,1) * A2(:,1,4)
c    &           + tau(:,2) * ( 
c    &                         + A2(:,1,2) * A2(:,2,4)
c    &                         + A2(:,1,3) * A2(:,3,4)
c    &                         + A2(:,1,4) * A2(:,4,4)
c    &                                                  )
     &           + tau(:,3) *    A2(:,1,5) * A2(:,5,4)
c
          A0(:,1,5)= 
     &             tau(:,1) *    A2(:,1,1) * A2(:,1,5)
     &           + tau(:,2) * ( 
c    &                         + A2(:,1,2) * A2(:,2,5)
     &                         + A2(:,1,3) * A2(:,3,5)
c    &                         + A2(:,1,4) * A2(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A2(:,1,5) * A2(:,5,5)
c
c  row two
c
          A0(:,2,1)= 
     &             tau(:,1) *    A2(:,2,1) * A2(:,1,1)
     &           + tau(:,2) * ( 
     &                         + A2(:,2,2) * A2(:,2,1)
     &                         + A2(:,2,3) * A2(:,3,1)
c    &                         + A2(:,2,4) * A2(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A2(:,2,5) * A2(:,5,1)
c
          A0(:,2,2)= 
c    &             tau(:,1) *    A2(:,2,1) * A2(:,1,2)
     &           + tau(:,2) * ( 
     &                         + A2(:,2,2) * A2(:,2,2)
c    &                         + A2(:,2,3) * A2(:,3,2)
c    &                         + A2(:,2,4) * A2(:,4,2)
     &                                                  )
     &           + tau(:,3) *    A2(:,2,5) * A2(:,5,2)
c
          A0(:,2,3)= 
     &             tau(:,1) *    A2(:,2,1) * A2(:,1,3)
     &           + tau(:,2) * ( 
     &                         + A2(:,2,2) * A2(:,2,3)
     &                         + A2(:,2,3) * A2(:,3,3)
c    &                         + A2(:,2,4) * A2(:,4,3)
     &                                                  )
     &           + tau(:,3) *    A2(:,2,5) * A2(:,5,3)
c
          A0(:,2,4)= 
c    &             tau(:,1) *    A2(:,2,1) * A2(:,1,4)
c    &           + tau(:,2) * ( 
c    &                         + A2(:,2,2) * A2(:,2,4)
c    &                         + A2(:,2,3) * A2(:,3,4)
c    &                         + A2(:,2,4) * A2(:,4,4)
c    &                                                  )
     &           + tau(:,3) *    A2(:,2,5) * A2(:,5,4)
c
          A0(:,2,5)= 
     &             tau(:,1) *    A2(:,2,1) * A2(:,1,5)
     &           + tau(:,2) * ( 
     &                         + A2(:,2,2) * A2(:,2,5)
     &                         + A2(:,2,3) * A2(:,3,5)
c    &                         + A2(:,2,4) * A2(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A2(:,2,5) * A2(:,5,5)
c
c  row three
c
          A0(:,3,1)= 
     &             tau(:,1) *    A2(:,3,1) * A2(:,1,1)
     &           + tau(:,2) * ( 
c    &                         + A2(:,3,2) * A2(:,2,1)
     &                         + A2(:,3,3) * A2(:,3,1)
c    &                         + A2(:,3,4) * A2(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A2(:,3,5) * A2(:,5,1)
c
          A0(:,3,2)= 
c    &             tau(:,1) *    A2(:,3,1) * A2(:,1,2)
c    &           + tau(:,2) * ( 
c    &                         + A2(:,3,2) * A2(:,2,2)
c    &                         + A2(:,3,3) * A2(:,3,2)
c    &                         + A2(:,3,4) * A2(:,4,2)
c    &                                                  )
     &           + tau(:,3) *    A2(:,3,5) * A2(:,5,2)
c
          A0(:,3,3)= 
     &             tau(:,1) *    A2(:,3,1) * A2(:,1,3)
     &           + tau(:,2) * ( 
c    &                         + A2(:,3,2) * A2(:,2,3)
     &                         + A2(:,3,3) * A2(:,3,3)
c    &                         + A2(:,3,4) * A2(:,4,3)
     &                                                  )
     &           + tau(:,3) *    A2(:,3,5) * A2(:,5,3)
c
          A0(:,3,4)= 
c    &             tau(:,1) *    A2(:,3,1) * A2(:,1,4)
c    &           + tau(:,2) * ( 
c    &                         + A2(:,3,2) * A2(:,2,4)
c    &                         + A2(:,3,3) * A2(:,3,4)
c    &                         + A2(:,3,4) * A2(:,4,4)
c    &                                                  )
     &           + tau(:,3) *    A2(:,3,5) * A2(:,5,4)
c
          A0(:,3,5)= 
     &             tau(:,1) *    A2(:,3,1) * A2(:,1,5)
     &           + tau(:,2) * ( 
c    &                         + A2(:,3,2) * A2(:,2,5)
     &                         + A2(:,3,3) * A2(:,3,5)
c    &                         + A2(:,3,4) * A2(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A2(:,3,5) * A2(:,5,5)
c
c  row four
c
          A0(:,4,1)= 
     &             tau(:,1) *    A2(:,4,1) * A2(:,1,1)
     &           + tau(:,2) * ( 
c    &                         + A2(:,4,2) * A2(:,2,1)
     &                         + A2(:,4,3) * A2(:,3,1)
     &                         + A2(:,4,4) * A2(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A2(:,4,5) * A2(:,5,1)
c
          A0(:,4,2)= 
c    &             tau(:,1) *    A2(:,4,1) * A2(:,1,2)
c    &           + tau(:,2) * ( 
c    &                         + A2(:,4,2) * A2(:,2,2)
c    &                         + A2(:,4,3) * A2(:,3,2)
c    &                         + A2(:,4,4) * A2(:,4,2)
c    &                                                  )
     &           + tau(:,3) *    A2(:,4,5) * A2(:,5,2)
c
          A0(:,4,3)= 
     &             tau(:,1) *    A2(:,4,1) * A2(:,1,3)
     &           + tau(:,2) * ( 
c    &                         + A2(:,4,2) * A2(:,2,3)
     &                         + A2(:,4,3) * A2(:,3,3)
     &                         + A2(:,4,4) * A2(:,4,3)
     &                                                  )
     &           + tau(:,3) *    A2(:,4,5) * A2(:,5,3)
c
          A0(:,4,4)= 
c    &             tau(:,1) *    A2(:,4,1) * A2(:,1,4)
     &           + tau(:,2) * ( 
c    &                         + A2(:,4,2) * A2(:,2,4)
c    &                         + A2(:,4,3) * A2(:,3,4)
     &                         + A2(:,4,4) * A2(:,4,4)
     &                                                  )
     &           + tau(:,3) *    A2(:,4,5) * A2(:,5,4)
c
          A0(:,4,5)= 
     &             tau(:,1) *    A2(:,4,1) * A2(:,1,5)
     &           + tau(:,2) * ( 
c    &                         + A2(:,4,2) * A2(:,2,5)
     &                         + A2(:,4,3) * A2(:,3,5)
     &                         + A2(:,4,4) * A2(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A2(:,4,5) * A2(:,5,5)
c
c  row two
c
          A0(:,5,1)= 
     &             tau(:,1) *    A2(:,5,1) * A2(:,1,1)
     &           + tau(:,2) * ( 
     &                         + A2(:,5,2) * A2(:,2,1)
     &                         + A2(:,5,3) * A2(:,3,1)
     &                         + A2(:,5,4) * A2(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A2(:,5,5) * A2(:,5,1)
c
          A0(:,5,2)= 
c    &             tau(:,1) *    A2(:,5,1) * A2(:,1,2)
     &           + tau(:,2) * ( 
     &                         + A2(:,5,2) * A2(:,2,2)
c    &                         + A2(:,5,3) * A2(:,3,2)
c    &                         + A2(:,5,4) * A2(:,4,2)
     &                                                  )
     &           + tau(:,3) *    A2(:,5,5) * A2(:,5,2)
c
          A0(:,5,3)= 
     &             tau(:,1) *    A2(:,5,1) * A2(:,1,3)
     &           + tau(:,2) * ( 
     &                         + A2(:,5,2) * A2(:,2,3)
     &                         + A2(:,5,3) * A2(:,3,3)
     &                         + A2(:,5,4) * A2(:,4,3)
     &                                                  )
     &           + tau(:,3) *    A2(:,5,5) * A2(:,5,3)
c
          A0(:,5,4)= 
c    &             tau(:,1) *    A2(:,5,1) * A2(:,1,4)
     &           + tau(:,2) * ( 
c    &                         + A2(:,5,2) * A2(:,2,4)
c    &                         + A2(:,5,3) * A2(:,3,4)
     &                         + A2(:,5,4) * A2(:,4,4)
     &                                                  )
     &           + tau(:,3) *    A2(:,5,5) * A2(:,5,4)
c
          A0(:,5,5)= 
     &             tau(:,1) *    A2(:,5,1) * A2(:,1,5)
     &           + tau(:,2) * ( 
     &                         + A2(:,5,2) * A2(:,2,5)
     &                         + A2(:,5,3) * A2(:,3,5)
     &                         + A2(:,5,4) * A2(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A2(:,5,5) * A2(:,5,5)
          else
             A0 = zero
          endif
c
          if (Navier .eq. 1) then
           A0(:,2,2) = A0(:,2,2) + rmu !rmu
           A0(:,3,3) = A0(:,3,3) + rlm2mu !rlm2mu
           A0(:,4,4) = A0(:,4,4) + rmu !rmu
           A0(:,5,2) = A0(:,5,2) + rmu * u1
           A0(:,5,3) = A0(:,5,3) + rlm2mu * u2
           A0(:,5,4) = A0(:,5,4) + rmu * u3
           A0(:,5,5) = A0(:,5,5) + con !con
        endif
c
        do j = 1, nshl
          tmp = WdetJ * shg(:,j,2) * shg(:,j,2)
           do i=1,nflow
           do k=1,nflow
             BDiagl(:,j,i,k) = BDiagl(:,j,i,k) + tmp * A0(:,i,k)
           enddo
           enddo
        enddo
c
c  start with 3 3 contribution
c
       if(ivart.ge.2) then
c
c  row one
c
          A0(:,1,1)= 
     &             tau(:,1) *    A3(:,1,1) * A3(:,1,1)
     &           + tau(:,2) * ( 
c    &                         + A3(:,1,2) * A3(:,2,1)
c    &                         + A3(:,1,3) * A3(:,3,1)
     &                         + A3(:,1,4) * A3(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A3(:,1,5) * A3(:,5,1)
c
          A0(:,1,2)= 
c    &             tau(:,1) *    A3(:,1,1) * A3(:,1,2)
c    &           + tau(:,2) * ( 
c    &                         + A3(:,1,2) * A3(:,2,2)
c    &                         + A3(:,1,3) * A3(:,3,2)
c    &                         + A3(:,1,4) * A3(:,4,2)
c    &                                                  )
     &           + tau(:,3) *    A3(:,1,5) * A3(:,5,2)
c
          A0(:,1,3)= 
c    &             tau(:,1) *    A3(:,1,1) * A3(:,1,3)
c    &           + tau(:,2) * ( 
c    &                         + A3(:,1,2) * A3(:,2,3)
c    &                         + A3(:,1,3) * A3(:,3,3)
c    &                         + A3(:,1,4) * A3(:,4,3)
c    &                                                  )
     &           + tau(:,3) *    A3(:,1,5) * A3(:,5,3)
c
          A0(:,1,4)= 
     &             tau(:,1) *    A3(:,1,1) * A3(:,1,4)
     &           + tau(:,2) * ( 
c    &                         + A3(:,1,2) * A3(:,2,4)
c    &                         + A3(:,1,3) * A3(:,3,4)
     &                         + A3(:,1,4) * A3(:,4,4)
     &                                                  )
     &           + tau(:,3) *    A3(:,1,5) * A3(:,5,4)
c
          A0(:,1,5)= 
     &             tau(:,1) *    A3(:,1,1) * A3(:,1,5)
     &           + tau(:,2) * ( 
c    &                         + A3(:,1,2) * A3(:,2,5)
c    &                         + A3(:,1,3) * A3(:,3,5)
     &                         + A3(:,1,4) * A3(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A3(:,1,5) * A3(:,5,5)
c
c  row two
c
          A0(:,2,1)= 
     &             tau(:,1) *    A3(:,2,1) * A3(:,1,1)
     &           + tau(:,2) * ( 
     &                         + A3(:,2,2) * A3(:,2,1)
c    &                         + A3(:,2,3) * A3(:,3,1)
     &                         + A3(:,2,4) * A3(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A3(:,2,5) * A3(:,5,1)
c
          A0(:,2,2)= 
c    &             tau(:,1) *    A3(:,2,1) * A3(:,1,2)
     &           + tau(:,2) * ( 
     &                         + A3(:,2,2) * A3(:,2,2)
c    &                         + A3(:,2,3) * A3(:,3,2)
c    &                         + A3(:,2,4) * A3(:,4,2)
     &                                                  )
     &           + tau(:,3) *    A3(:,2,5) * A3(:,5,2)
c
          A0(:,2,3)= 
c    &             tau(:,1) *    A3(:,2,1) * A3(:,1,3)
c    &           + tau(:,2) * ( 
c    &                         + A3(:,2,2) * A3(:,2,3)
c    &                         + A3(:,2,3) * A3(:,3,3)
c    &                         + A3(:,2,4) * A3(:,4,3)
c    &                                                  )
     &           + tau(:,3) *    A3(:,2,5) * A3(:,5,3)
c
          A0(:,2,4)= 
     &             tau(:,1) *    A3(:,2,1) * A3(:,1,4)
     &           + tau(:,2) * ( 
     &                         + A3(:,2,2) * A3(:,2,4)
c    &                         + A3(:,2,3) * A3(:,3,4)
     &                         + A3(:,2,4) * A3(:,4,4)
     &                                                  )
     &           + tau(:,3) *    A3(:,2,5) * A3(:,5,4)
c
          A0(:,2,5)= 
     &             tau(:,1) *    A3(:,2,1) * A3(:,1,5)
     &           + tau(:,2) * ( 
     &                         + A3(:,2,2) * A3(:,2,5)
c    &                         + A3(:,2,3) * A3(:,3,5)
     &                         + A3(:,2,4) * A3(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A3(:,2,5) * A3(:,5,5)
c
c  row three
c
          A0(:,3,1)= 
     &             tau(:,1) *    A3(:,3,1) * A3(:,1,1)
     &           + tau(:,2) * ( 
c    &                         + A3(:,3,2) * A3(:,2,1)
     &                         + A3(:,3,3) * A3(:,3,1)
     &                         + A3(:,3,4) * A3(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A3(:,3,5) * A3(:,5,1)
c
          A0(:,3,2)= 
c    &             tau(:,1) *    A3(:,3,1) * A3(:,1,2)
c    &           + tau(:,2) * ( 
c    &                         + A3(:,3,2) * A3(:,2,2)
c    &                         + A3(:,3,3) * A3(:,3,2)
c    &                         + A3(:,3,4) * A3(:,4,2)
c    &                                                  )
     &           + tau(:,3) *    A3(:,3,5) * A3(:,5,2)
c
          A0(:,3,3)= 
c    &             tau(:,1) *    A3(:,3,1) * A3(:,1,3)
     &           + tau(:,2) * ( 
c    &                         + A3(:,3,2) * A3(:,2,3)
     &                         + A3(:,3,3) * A3(:,3,3)
c    &                         + A3(:,3,4) * A3(:,4,3)
     &                                                  )
     &           + tau(:,3) *    A3(:,3,5) * A3(:,5,3)
c
          A0(:,3,4)= 
     &             tau(:,1) *    A3(:,3,1) * A3(:,1,4)
     &           + tau(:,2) * ( 
c    &                         + A3(:,3,2) * A3(:,2,4)
     &                         + A3(:,3,3) * A3(:,3,4)
     &                         + A3(:,3,4) * A3(:,4,4)
     &                                                  )
     &           + tau(:,3) *    A3(:,3,5) * A3(:,5,4)
c
          A0(:,3,5)= 
     &             tau(:,1) *    A3(:,3,1) * A3(:,1,5)
     &           + tau(:,2) * ( 
c    &                         + A3(:,3,2) * A3(:,2,5)
     &                         + A3(:,3,3) * A3(:,3,5)
     &                         + A3(:,3,4) * A3(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A3(:,3,5) * A3(:,5,5)
c
c  row four
c
          A0(:,4,1)= 
     &             tau(:,1) *    A3(:,4,1) * A3(:,1,1)
     &           + tau(:,2) * ( 
c    &                         + A3(:,4,2) * A3(:,2,1)
c    &                         + A3(:,4,3) * A3(:,3,1)
     &                         + A3(:,4,4) * A3(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A3(:,4,5) * A3(:,5,1)
c
          A0(:,4,2)= 
c    &             tau(:,1) *    A3(:,4,1) * A3(:,1,2)
c    &           + tau(:,2) * ( 
c    &                         + A3(:,4,2) * A3(:,2,2)
c    &                         + A3(:,4,3) * A3(:,3,2)
c    &                         + A3(:,4,4) * A3(:,4,2)
c    &                                                  )
     &           + tau(:,3) *    A3(:,4,5) * A3(:,5,2)
c
          A0(:,4,3)= 
c    &             tau(:,1) *    A3(:,4,1) * A3(:,1,3)
c    &           + tau(:,2) * ( 
c    &                         + A3(:,4,2) * A3(:,2,3)
c    &                         + A3(:,4,3) * A3(:,3,3)
c    &                         + A3(:,4,4) * A3(:,4,3)
c    &                                                  )
     &           + tau(:,3) *    A3(:,4,5) * A3(:,5,3)
c
          A0(:,4,4)= 
     &             tau(:,1) *    A3(:,4,1) * A3(:,1,4)
     &           + tau(:,2) * ( 
c    &                         + A3(:,4,2) * A3(:,2,4)
c    &                         + A3(:,4,3) * A3(:,3,4)
     &                         + A3(:,4,4) * A3(:,4,4)
     &                                                  )
     &           + tau(:,3) *    A3(:,4,5) * A3(:,5,4)
c
          A0(:,4,5)= 
     &             tau(:,1) *    A3(:,4,1) * A3(:,1,5)
     &           + tau(:,2) * ( 
c    &                         + A3(:,4,2) * A3(:,2,5)
c    &                         + A3(:,4,3) * A3(:,3,5)
     &                         + A3(:,4,4) * A3(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A3(:,4,5) * A3(:,5,5)
c
c  row two
c
          A0(:,5,1)= 
     &             tau(:,1) *    A3(:,5,1) * A3(:,1,1)
     &           + tau(:,2) * ( 
     &                         + A3(:,5,2) * A3(:,2,1)
     &                         + A3(:,5,3) * A3(:,3,1)
     &                         + A3(:,5,4) * A3(:,4,1)
     &                                                  )
     &           + tau(:,3) *    A3(:,5,5) * A3(:,5,1)
c
          A0(:,5,2)= 
c    &             tau(:,1) *    A3(:,5,1) * A3(:,1,2)
     &           + tau(:,2) * ( 
     &                         + A3(:,5,2) * A3(:,2,2)
c    &                         + A3(:,5,3) * A3(:,3,2)
c    &                         + A3(:,5,4) * A3(:,4,2)
     &                                                  )
     &           + tau(:,3) *    A3(:,5,5) * A3(:,5,2)
c
          A0(:,5,3)= 
c    &             tau(:,1) *    A3(:,5,1) * A3(:,1,3)
     &           + tau(:,2) * ( 
c    &                         + A3(:,5,2) * A3(:,2,3)
     &                         + A3(:,5,3) * A3(:,3,3)
c    &                         + A3(:,5,4) * A3(:,4,3)
     &                                                  )
     &           + tau(:,3) *    A3(:,5,5) * A3(:,5,3)
c
          A0(:,5,4)= 
     &             tau(:,1) *    A3(:,5,1) * A3(:,1,4)
     &           + tau(:,2) * ( 
     &                         + A3(:,5,2) * A3(:,2,4)
     &                         + A3(:,5,3) * A3(:,3,4)
     &                         + A3(:,5,4) * A3(:,4,4)
     &                                                  )
     &           + tau(:,3) *    A3(:,5,5) * A3(:,5,4)
c
          A0(:,5,5)= 
     &             tau(:,1) *    A3(:,5,1) * A3(:,1,5)
     &           + tau(:,2) * ( 
     &                         + A3(:,5,2) * A3(:,2,5)
     &                         + A3(:,5,3) * A3(:,3,5)
     &                         + A3(:,5,4) * A3(:,4,5)
     &                                                  )
     &           + tau(:,3) *    A3(:,5,5) * A3(:,5,5)
          else
             A0 = zero
          endif
c
          if (Navier .eq. 1) then
           A0(:,2,2) = A0(:,2,2) + rmu !rmu
           A0(:,3,3) = A0(:,3,3) + rmu !rmu
           A0(:,4,4) = A0(:,4,4) + rlm2mu !rlm2mu
           A0(:,5,2) = A0(:,5,2) + rmu * u1
           A0(:,5,3) = A0(:,5,3) + rmu * u2
           A0(:,5,4) = A0(:,5,4) + rlm2mu * u3
           A0(:,5,5) = A0(:,5,5) + con !con
          endif
c
        do j = 1, nshl
          tmp = WdetJ * shg(:,j,3) * shg(:,j,3)
           do i=1,nflow
           do k=1,nflow
             BDiagl(:,j,i,k) = BDiagl(:,j,i,k) + tmp * A0(:,i,k)
           enddo
           enddo
        enddo
c
c now for the 1 2  plus 2 1
c
          if(ivart.ge.2) then
c
c  row one
c
          A0(:,1,1)= 
     &             tau(:,1) * (
     &                         + A1(:,1,1) * A2(:,1,1)
     &                         + A2(:,1,1) * A1(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A2(:,2,1)
c    &                         + A1(:,1,3) * A2(:,3,1)
c    &                         + A1(:,1,4) * A2(:,4,1)
c    &                         + A2(:,1,2) * A1(:,2,1)
     &                         + A2(:,1,3) * A1(:,3,1)
c    &                         + A2(:,1,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,1,5) * A2(:,5,1)
     &		               + A2(:,1,5) * A1(:,5,1)
     &		                                        )
c
          A0(:,1,2)= 
     &             tau(:,1) * (
c    &                         + A1(:,1,1) * A2(:,1,2)
     &                         + A2(:,1,1) * A1(:,1,2)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A2(:,2,2)
c    &                         + A1(:,1,3) * A2(:,3,2)
c    &                         + A1(:,1,4) * A2(:,4,2)
c    &                         + A2(:,1,2) * A1(:,2,2)
     &                         + A2(:,1,3) * A1(:,3,2)
c    &                         + A2(:,1,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,1,5) * A2(:,5,2)
     &		               + A2(:,1,5) * A1(:,5,2)
     &		                                        )
c
          A0(:,1,3)= 
     &             tau(:,1) * (
     &                         + A1(:,1,1) * A2(:,1,3)
c    &                         + A2(:,1,1) * A1(:,1,3)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A2(:,2,3)
c    &                         + A1(:,1,3) * A2(:,3,3)
c    &                         + A1(:,1,4) * A2(:,4,3)
c    &                         + A2(:,1,2) * A1(:,2,3)
     &                         + A2(:,1,3) * A1(:,3,3)
c    &                         + A2(:,1,4) * A1(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,1,5) * A2(:,5,3)
     &		               + A2(:,1,5) * A1(:,5,3)
     &		                                        )
c
          A0(:,1,4)= 
c    &             tau(:,1) * (
c    &                         + A1(:,1,1) * A2(:,1,4)
c    &                         + A2(:,1,1) * A1(:,1,4)
c    &		                                        )
c    &           + tau(:,2) * ( 
c    &                         + A1(:,1,2) * A2(:,2,4)
c    &                         + A1(:,1,3) * A2(:,3,4)
c    &                         + A1(:,1,4) * A2(:,4,4)
c    &                         + A2(:,1,2) * A1(:,2,4)
c    &                         + A2(:,1,3) * A1(:,3,4)
c    &                         + A2(:,1,4) * A1(:,4,4)
c    &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,1,5) * A2(:,5,4)
     &		               + A2(:,1,5) * A1(:,5,4)
     &		                                        )
c
          A0(:,1,5)= 
     &             tau(:,1) * (
     &                         + A1(:,1,1) * A2(:,1,5)
     &                         + A2(:,1,1) * A1(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A2(:,2,5)
c    &                         + A1(:,1,3) * A2(:,3,5)
c    &                         + A1(:,1,4) * A2(:,4,5)
c    &                         + A2(:,1,2) * A1(:,2,5)
     &                         + A2(:,1,3) * A1(:,3,5)
c    &                         + A2(:,1,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,1,5) * A2(:,5,5)
     &		               + A2(:,1,5) * A1(:,5,5)
     &		                                        )
c
c  row two
c
          A0(:,2,1)= 
     &             tau(:,1) * (
     &                         + A1(:,2,1) * A2(:,1,1)
     &                         + A2(:,2,1) * A1(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A2(:,2,1)
c    &                         + A1(:,2,3) * A2(:,3,1)
c    &                         + A1(:,2,4) * A2(:,4,1)
     &                         + A2(:,2,2) * A1(:,2,1)
     &                         + A2(:,2,3) * A1(:,3,1)
c    &                         + A2(:,2,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,2,5) * A2(:,5,1)
     &		               + A2(:,2,5) * A1(:,5,1)
     &		                                        )
c
          A0(:,2,2)= 
     &             tau(:,1) * (
c    &                         + A1(:,2,1) * A2(:,1,2)
     &                         + A2(:,2,1) * A1(:,1,2)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A2(:,2,2)
c    &                         + A1(:,2,3) * A2(:,3,2)
c    &                         + A1(:,2,4) * A2(:,4,2)
     &                         + A2(:,2,2) * A1(:,2,2)
     &                         + A2(:,2,3) * A1(:,3,2)
c    &                         + A2(:,2,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,2,5) * A2(:,5,2)
     &		               + A2(:,2,5) * A1(:,5,2)
     &		                                        )
c
          A0(:,2,3)= 
     &             tau(:,1) * (
     &                         + A1(:,2,1) * A2(:,1,3)
c    &                         + A2(:,2,1) * A1(:,1,3)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A2(:,2,3)
c    &                         + A1(:,2,3) * A2(:,3,3)
c    &                         + A1(:,2,4) * A2(:,4,3)
c    &                         + A2(:,2,2) * A1(:,2,3)
     &                         + A2(:,2,3) * A1(:,3,3)
c    &                         + A2(:,2,4) * A1(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,2,5) * A2(:,5,3)
     &		               + A2(:,2,5) * A1(:,5,3)
     &		                                        )
c
          A0(:,2,4)= 
c    &             tau(:,1) * (
c    &                         + A1(:,2,1) * A2(:,1,4)
c    &                         + A2(:,2,1) * A1(:,1,4)
c    &		                                        )
c    &           + tau(:,2) * ( 
c    &                         + A1(:,2,2) * A2(:,2,4)
c    &                         + A1(:,2,3) * A2(:,3,4)
c    &                         + A1(:,2,4) * A2(:,4,4)
c    &                         + A2(:,2,2) * A1(:,2,4)
c    &                         + A2(:,2,3) * A1(:,3,4)
c    &                         + A2(:,2,4) * A1(:,4,4)
c    &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,2,5) * A2(:,5,4)
     &		               + A2(:,2,5) * A1(:,5,4)
     &		                                        )
c
          A0(:,2,5)= 
     &             tau(:,1) * (
     &                         + A1(:,2,1) * A2(:,1,5)
     &                         + A2(:,2,1) * A1(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A2(:,2,5)
c    &                         + A1(:,2,3) * A2(:,3,5)
c    &                         + A1(:,2,4) * A2(:,4,5)
     &                         + A2(:,2,2) * A1(:,2,5)
     &                         + A2(:,2,3) * A1(:,3,5)
c    &                         + A2(:,2,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,2,5) * A2(:,5,5)
     &		               + A2(:,2,5) * A1(:,5,5)
     &		                                        )
c
c  row three
c
          A0(:,3,1)= 
     &             tau(:,1) * (
     &                         + A1(:,3,1) * A2(:,1,1)
     &                         + A2(:,3,1) * A1(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A2(:,2,1)
     &                         + A1(:,3,3) * A2(:,3,1)
c    &                         + A1(:,3,4) * A2(:,4,1)
c    &                         + A2(:,3,2) * A1(:,2,1)
     &                         + A2(:,3,3) * A1(:,3,1)
c    &                         + A2(:,3,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,3,5) * A2(:,5,1)
     &		               + A2(:,3,5) * A1(:,5,1)
     &		                                        )
c
          A0(:,3,2)= 
     &             tau(:,1) * (
c    &                         + A1(:,3,1) * A2(:,1,2)
     &                         + A2(:,3,1) * A1(:,1,2)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A2(:,2,2)
c    &                         + A1(:,3,3) * A2(:,3,2)
c    &                         + A1(:,3,4) * A2(:,4,2)
c    &                         + A2(:,3,2) * A1(:,2,2)
     &                         + A2(:,3,3) * A1(:,3,2)
c    &                         + A2(:,3,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,3,5) * A2(:,5,2)
     &		               + A2(:,3,5) * A1(:,5,2)
     &		                                        )
c
          A0(:,3,3)= 
     &             tau(:,1) * (
     &                         + A1(:,3,1) * A2(:,1,3)
c    &                         + A2(:,3,1) * A1(:,1,3)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A2(:,2,3)
     &                         + A1(:,3,3) * A2(:,3,3)
c    &                         + A1(:,3,4) * A2(:,4,3)
c    &                         + A2(:,3,2) * A1(:,2,3)
     &                         + A2(:,3,3) * A1(:,3,3)
c    &                         + A2(:,3,4) * A1(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,3,5) * A2(:,5,3)
     &		               + A2(:,3,5) * A1(:,5,3)
     &		                                        )
c
          A0(:,3,4)= 
c    &             tau(:,1) * (
c    &                         + A1(:,3,1) * A2(:,1,4)
c    &                         + A2(:,3,1) * A1(:,1,4)
c    &		                                        )
c    &           + tau(:,2) * ( 
c    &                         + A1(:,3,2) * A2(:,2,4)
c    &                         + A1(:,3,3) * A2(:,3,4)
c    &                         + A1(:,3,4) * A2(:,4,4)
c    &                         + A2(:,3,2) * A1(:,2,4)
c    &                         + A2(:,3,3) * A1(:,3,4)
c    &                         + A2(:,3,4) * A1(:,4,4)
c    &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,3,5) * A2(:,5,4)
     &		               + A2(:,3,5) * A1(:,5,4)
     &		                                        )
c
          A0(:,3,5)= 
     &             tau(:,1) * (
     &                         + A1(:,3,1) * A2(:,1,5)
     &                         + A2(:,3,1) * A1(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A2(:,2,5)
     &                         + A1(:,3,3) * A2(:,3,5)
c    &                         + A1(:,3,4) * A2(:,4,5)
c    &                         + A2(:,3,2) * A1(:,2,5)
     &                         + A2(:,3,3) * A1(:,3,5)
c    &                         + A2(:,3,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,3,5) * A2(:,5,5)
     &		               + A2(:,3,5) * A1(:,5,5)
     &		                                        )
c
c  row four
c
          A0(:,4,1)= 
     &             tau(:,1) * (
     &                         + A1(:,4,1) * A2(:,1,1)
     &                         + A2(:,4,1) * A1(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A2(:,2,1)
c    &                         + A1(:,4,3) * A2(:,3,1)
     &                         + A1(:,4,4) * A2(:,4,1)
c    &                         + A2(:,4,2) * A1(:,2,1)
     &                         + A2(:,4,3) * A1(:,3,1)
     &                         + A2(:,4,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,4,5) * A2(:,5,1)
     &		               + A2(:,4,5) * A1(:,5,1)
     &		                                        )
c
          A0(:,4,2)= 
     &             tau(:,1) * (
c    &                         + A1(:,4,1) * A2(:,1,2)
     &                         + A2(:,4,1) * A1(:,1,2)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A2(:,2,2)
c    &                         + A1(:,4,3) * A2(:,3,2)
c    &                         + A1(:,4,4) * A2(:,4,2)
c    &                         + A2(:,4,2) * A1(:,2,2)
     &                         + A2(:,4,3) * A1(:,3,2)
     &                         + A2(:,4,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,4,5) * A2(:,5,2)
     &		               + A2(:,4,5) * A1(:,5,2)
     &		                                        )
c
          A0(:,4,3)= 
     &             tau(:,1) * (
     &                         + A1(:,4,1) * A2(:,1,3)
c    &                         + A2(:,4,1) * A1(:,1,3)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A2(:,2,3)
c    &                         + A1(:,4,3) * A2(:,3,3)
     &                         + A1(:,4,4) * A2(:,4,3)
c    &                         + A2(:,4,2) * A1(:,2,3)
     &                         + A2(:,4,3) * A1(:,3,3)
c    &                         + A2(:,4,4) * A1(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,4,5) * A2(:,5,3)
     &		               + A2(:,4,5) * A1(:,5,3)
     &		                                        )
c
          A0(:,4,4)= 
c    &             tau(:,1) * (
c    &                         + A1(:,4,1) * A2(:,1,4)
c    &                         + A2(:,4,1) * A1(:,1,4)
c    &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A1(:,4,2) * A2(:,2,4)
c    &                         + A1(:,4,3) * A2(:,3,4)
     &                         + A1(:,4,4) * A2(:,4,4)
c    &                         + A2(:,4,2) * A1(:,2,4)
c    &                         + A2(:,4,3) * A1(:,3,4)
     &                         + A2(:,4,4) * A1(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,4,5) * A2(:,5,4)
     &		               + A2(:,4,5) * A1(:,5,4)
     &		                                        )
c
          A0(:,4,5)= 
     &             tau(:,1) * (
     &                         + A1(:,4,1) * A2(:,1,5)
     &                         + A2(:,4,1) * A1(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A2(:,2,5)
c    &                         + A1(:,4,3) * A2(:,3,5)
     &                         + A1(:,4,4) * A2(:,4,5)
c    &                         + A2(:,4,2) * A1(:,2,5)
     &                         + A2(:,4,3) * A1(:,3,5)
     &                         + A2(:,4,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,4,5) * A2(:,5,5)
     &		               + A2(:,4,5) * A1(:,5,5)
     &		                                        )
c
c  row five
c
          A0(:,5,1)= 
     &             tau(:,1) * (
     &                         + A1(:,5,1) * A2(:,1,1)
     &                         + A2(:,5,1) * A1(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,5,2) * A2(:,2,1)
     &                         + A1(:,5,3) * A2(:,3,1)
     &                         + A1(:,5,4) * A2(:,4,1)
     &                         + A2(:,5,2) * A1(:,2,1)
     &                         + A2(:,5,3) * A1(:,3,1)
     &                         + A2(:,5,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,5,5) * A2(:,5,1)
     &		               + A2(:,5,5) * A1(:,5,1)
     &		                                        )
c
          A0(:,5,2)= 
     &             tau(:,1) * (
c    &                         + A1(:,5,1) * A2(:,1,2)
     &                         + A2(:,5,1) * A1(:,1,2)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,5,2) * A2(:,2,2)
c    &                         + A1(:,5,3) * A2(:,3,2)
c    &                         + A1(:,5,4) * A2(:,4,2)
     &                         + A2(:,5,2) * A1(:,2,2)
     &                         + A2(:,5,3) * A1(:,3,2)
     &                         + A2(:,5,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,5,5) * A2(:,5,2)
     &		               + A2(:,5,5) * A1(:,5,2)
     &		                                        )
c
          A0(:,5,3)= 
     &             tau(:,1) * (
     &                         + A1(:,5,1) * A2(:,1,3)
c    &                         + A2(:,5,1) * A1(:,1,3)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,5,2) * A2(:,2,3)
     &                         + A1(:,5,3) * A2(:,3,3)
     &                         + A1(:,5,4) * A2(:,4,3)
c    &                         + A2(:,5,2) * A1(:,2,3)
     &                         + A2(:,5,3) * A1(:,3,3)
c    &                         + A2(:,5,4) * A1(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,5,5) * A2(:,5,3)
     &		               + A2(:,5,5) * A1(:,5,3)
     &		                                        )
c
          A0(:,5,4)= 
c    &             tau(:,1) * (
c    &                         + A1(:,5,1) * A2(:,1,4)
c    &                         + A2(:,5,1) * A1(:,1,4)
c    &		                                        )
     &           + tau(:,2) * (
c    &                         + A1(:,5,2) * A2(:,2,4)
c    &                         + A1(:,5,3) * A2(:,3,4)
     &                         + A1(:,5,4) * A2(:,4,4)
c    &                         + A2(:,5,2) * A1(:,2,4)
c    &                         + A2(:,5,3) * A1(:,3,4)
     &                         + A2(:,5,4) * A1(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,5,5) * A2(:,5,4)
     &		               + A2(:,5,5) * A1(:,5,4)
     &		                                        )
c
          A0(:,5,5)= 
     &             tau(:,1) * (
     &                         + A1(:,5,1) * A2(:,1,5)
     &                         + A2(:,5,1) * A1(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,5,2) * A2(:,2,5)
     &                         + A1(:,5,3) * A2(:,3,5)
     &                         + A1(:,5,4) * A2(:,4,5)
     &                         + A2(:,5,2) * A1(:,2,5)
     &                         + A2(:,5,3) * A1(:,3,5)
     &                         + A2(:,5,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,5,5) * A2(:,5,5)
     &		               + A2(:,5,5) * A1(:,5,5)
     &		                                        )
        else
           A0 = zero
        endif
c
        if (Navier .eq. 1) then
           A0(:,2,3) = A0(:,2,3) + rlm2mu - rmu
           A0(:,3,2) = A0(:,3,2) + rlm2mu - rmu
           A0(:,5,2) = A0(:,5,2) + (rlm2mu - rmu) * u2
           A0(:,5,3) = A0(:,5,3) + (rlm2mu- rmu) * u1
        endif
c
        do j = 1, nshl
          tmp = WdetJ * shg(:,j,1) * shg(:,j,2)
           do i=1,nflow
           do k=1,nflow
             BDiagl(:,j,i,k) = BDiagl(:,j,i,k) + tmp * A0(:,i,k)
           enddo
           enddo
         enddo
c
c now for the 1 3  plus 3 1
c
          if(ivart.ge.2) then
c
c  row one
c
          A0(:,1,1)= 
     &             tau(:,1) * (
     &                         + A1(:,1,1) * A3(:,1,1)
     &                         + A3(:,1,1) * A1(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A3(:,2,1)
c    &                         + A1(:,1,3) * A3(:,3,1)
c    &                         + A1(:,1,4) * A3(:,4,1)
c    &                         + A3(:,1,2) * A1(:,2,1)
c    &                         + A3(:,1,3) * A1(:,3,1)
     &                         + A3(:,1,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,1,5) * A3(:,5,1)
     &		               + A3(:,1,5) * A1(:,5,1)
     &		                                        )
c
          A0(:,1,2)= 
     &             tau(:,1) * (
c    &                         + A1(:,1,1) * A3(:,1,2)
     &                         + A3(:,1,1) * A1(:,1,2)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A3(:,2,2)
c    &                         + A1(:,1,3) * A3(:,3,2)
c    &                         + A1(:,1,4) * A3(:,4,2)
c    &                         + A3(:,1,2) * A1(:,2,2)
c    &                         + A3(:,1,3) * A1(:,3,2)
     &                         + A3(:,1,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,1,5) * A3(:,5,2)
     &		               + A3(:,1,5) * A1(:,5,2)
     &		                                        )
c
          A0(:,1,3)= 
c    &             tau(:,1) * (
c    &                         + A1(:,1,1) * A3(:,1,3)
c    &                         + A3(:,1,1) * A1(:,1,3)
c    &		                                        )
c    &           + tau(:,2) * ( 
c    &                         + A1(:,1,2) * A3(:,2,3)
c    &                         + A1(:,1,3) * A3(:,3,3)
c    &                         + A1(:,1,4) * A3(:,4,3)
c    &                         + A3(:,1,2) * A1(:,2,3)
c    &                         + A3(:,1,3) * A1(:,3,3)
c    &                         + A3(:,1,4) * A1(:,4,3)
c    &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,1,5) * A3(:,5,3)
     &		               + A3(:,1,5) * A1(:,5,3)
     &		                                        )
c
          A0(:,1,4)= 
     &             tau(:,1) * (
     &                         + A1(:,1,1) * A3(:,1,4)
c    &                         + A3(:,1,1) * A1(:,1,4)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A3(:,2,4)
c    &                         + A1(:,1,3) * A3(:,3,4)
c    &                         + A1(:,1,4) * A3(:,4,4)
c    &                         + A3(:,1,2) * A1(:,2,4)
c    &                         + A3(:,1,3) * A1(:,3,4)
     &                         + A3(:,1,4) * A1(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,1,5) * A3(:,5,4)
     &		               + A3(:,1,5) * A1(:,5,4)
     &		                                        )
c
          A0(:,1,5)= 
     &             tau(:,1) * (
     &                         + A1(:,1,1) * A3(:,1,5)
     &                         + A3(:,1,1) * A1(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,1,2) * A3(:,2,5)
c    &                         + A1(:,1,3) * A3(:,3,5)
c    &                         + A1(:,1,4) * A3(:,4,5)
c    &                         + A3(:,1,2) * A1(:,2,5)
c    &                         + A3(:,1,3) * A1(:,3,5)
     &                         + A3(:,1,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,1,5) * A3(:,5,5)
     &		               + A3(:,1,5) * A1(:,5,5)
     &		                                        )
c
c  row two
c
          A0(:,2,1)= 
     &             tau(:,1) * (
     &                         + A1(:,2,1) * A3(:,1,1)
     &                         + A3(:,2,1) * A1(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A3(:,2,1)
c    &                         + A1(:,2,3) * A3(:,3,1)
c    &                         + A1(:,2,4) * A3(:,4,1)
     &                         + A3(:,2,2) * A1(:,2,1)
c    &                         + A3(:,2,3) * A1(:,3,1)
     &                         + A3(:,2,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,2,5) * A3(:,5,1)
     &		               + A3(:,2,5) * A1(:,5,1)
     &		                                        )
c
          A0(:,2,2)= 
     &             tau(:,1) * (
c    &                         + A1(:,2,1) * A3(:,1,2)
     &                         + A3(:,2,1) * A1(:,1,2)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A3(:,2,2)
c    &                         + A1(:,2,3) * A3(:,3,2)
c    &                         + A1(:,2,4) * A3(:,4,2)
     &                         + A3(:,2,2) * A1(:,2,2)
c    &                         + A3(:,2,3) * A1(:,3,2)
     &                         + A3(:,2,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,2,5) * A3(:,5,2)
     &		               + A3(:,2,5) * A1(:,5,2)
     &		                                        )
c
          A0(:,2,3)= 
c    &             tau(:,1) * (
c    &                         + A1(:,2,1) * A3(:,1,3)
c    &                         + A3(:,2,1) * A1(:,1,3)
c    &		                                        )
c    &           + tau(:,2) * ( 
c    &                         + A1(:,2,2) * A3(:,2,3)
c    &                         + A1(:,2,3) * A3(:,3,3)
c    &                         + A1(:,2,4) * A3(:,4,3)
c    &                         + A3(:,2,2) * A1(:,2,3)
c    &                         + A3(:,2,3) * A1(:,3,3)
c    &                         + A3(:,2,4) * A1(:,4,3)
c    &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,2,5) * A3(:,5,3)
     &		               + A3(:,2,5) * A1(:,5,3)
     &		                                        )
c
          A0(:,2,4)= 
     &             tau(:,1) * (
     &                         + A1(:,2,1) * A3(:,1,4)
c    &                         + A3(:,2,1) * A1(:,1,4)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A3(:,2,4)
c    &                         + A1(:,2,3) * A3(:,3,4)
c    &                         + A1(:,2,4) * A3(:,4,4)
c    &                         + A3(:,2,2) * A1(:,2,4)
c    &                         + A3(:,2,3) * A1(:,3,4)
     &                         + A3(:,2,4) * A1(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,2,5) * A3(:,5,4)
     &		               + A3(:,2,5) * A1(:,5,4)
     &		                                        )
c
          A0(:,2,5)= 
     &             tau(:,1) * (
     &                         + A1(:,2,1) * A3(:,1,5)
     &                         + A3(:,2,1) * A1(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,2,2) * A3(:,2,5)
c    &                         + A1(:,2,3) * A3(:,3,5)
c    &                         + A1(:,2,4) * A3(:,4,5)
     &                         + A3(:,2,2) * A1(:,2,5)
c    &                         + A3(:,2,3) * A1(:,3,5)
     &                         + A3(:,2,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,2,5) * A3(:,5,5)
     &		               + A3(:,2,5) * A1(:,5,5)
     &		                                        )
c
c  row three
c
          A0(:,3,1)= 
     &             tau(:,1) * (
     &                         + A1(:,3,1) * A3(:,1,1)
     &                         + A3(:,3,1) * A1(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A3(:,2,1)
     &                         + A1(:,3,3) * A3(:,3,1)
c    &                         + A1(:,3,4) * A3(:,4,1)
c    &                         + A3(:,3,2) * A1(:,2,1)
     &                         + A3(:,3,3) * A1(:,3,1)
     &                         + A3(:,3,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,3,5) * A3(:,5,1)
     &		               + A3(:,3,5) * A1(:,5,1)
     &		                                        )
c
          A0(:,3,2)= 
     &             tau(:,1) * (
c    &                         + A1(:,3,1) * A3(:,1,2)
     &                         + A3(:,3,1) * A1(:,1,2)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A3(:,2,2)
c    &                         + A1(:,3,3) * A3(:,3,2)
c    &                         + A1(:,3,4) * A3(:,4,2)
c    &                         + A3(:,3,2) * A1(:,2,2)
     &                         + A3(:,3,3) * A1(:,3,2)
     &                         + A3(:,3,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,3,5) * A3(:,5,2)
     &		               + A3(:,3,5) * A1(:,5,2)
     &		                                        )
c
          A0(:,3,3)= 
c    &             tau(:,1) * (
c    &                         + A1(:,3,1) * A3(:,1,3)
c    &                         + A3(:,3,1) * A1(:,1,3)
c    &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A1(:,3,2) * A3(:,2,3)
     &                         + A1(:,3,3) * A3(:,3,3)
c    &                         + A1(:,3,4) * A3(:,4,3)
c    &                         + A3(:,3,2) * A1(:,2,3)
     &                         + A3(:,3,3) * A1(:,3,3)
c    &                         + A3(:,3,4) * A1(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,3,5) * A3(:,5,3)
     &		               + A3(:,3,5) * A1(:,5,3)
     &		                                        )
c
          A0(:,3,4)= 
     &             tau(:,1) * (
     &                         + A1(:,3,1) * A3(:,1,4)
c    &                         + A3(:,3,1) * A1(:,1,4)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A3(:,2,4)
     &                         + A1(:,3,3) * A3(:,3,4)
c    &                         + A1(:,3,4) * A3(:,4,4)
c    &                         + A3(:,3,2) * A1(:,2,4)
c    &                         + A3(:,3,3) * A1(:,3,4)
     &                         + A3(:,3,4) * A1(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,3,5) * A3(:,5,4)
     &		               + A3(:,3,5) * A1(:,5,4)
     &		                                        )
c
          A0(:,3,5)= 
     &             tau(:,1) * (
     &                         + A1(:,3,1) * A3(:,1,5)
     &                         + A3(:,3,1) * A1(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,3,2) * A3(:,2,5)
     &                         + A1(:,3,3) * A3(:,3,5)
c    &                         + A1(:,3,4) * A3(:,4,5)
c    &                         + A3(:,3,2) * A1(:,2,5)
     &                         + A3(:,3,3) * A1(:,3,5)
     &                         + A3(:,3,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,3,5) * A3(:,5,5)
     &		               + A3(:,3,5) * A1(:,5,5)
     &		                                        )
c
c  row four
c
          A0(:,4,1)= 
     &             tau(:,1) * (
     &                         + A1(:,4,1) * A3(:,1,1)
     &                         + A3(:,4,1) * A1(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A3(:,2,1)
c    &                         + A1(:,4,3) * A3(:,3,1)
     &                         + A1(:,4,4) * A3(:,4,1)
c    &                         + A3(:,4,2) * A1(:,2,1)
c    &                         + A3(:,4,3) * A1(:,3,1)
     &                         + A3(:,4,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,4,5) * A3(:,5,1)
     &		               + A3(:,4,5) * A1(:,5,1)
     &		                                        )
c
          A0(:,4,2)= 
     &             tau(:,1) * (
c    &                         + A1(:,4,1) * A3(:,1,2)
     &                         + A3(:,4,1) * A1(:,1,2)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A3(:,2,2)
c    &                         + A1(:,4,3) * A3(:,3,2)
c    &                         + A1(:,4,4) * A3(:,4,2)
c    &                         + A3(:,4,2) * A1(:,2,2)
c    &                         + A3(:,4,3) * A1(:,3,2)
     &                         + A3(:,4,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,4,5) * A3(:,5,2)
     &		               + A3(:,4,5) * A1(:,5,2)
     &		                                        )
c
          A0(:,4,3)= 
c    &             tau(:,1) * (
c    &                         + A1(:,4,1) * A3(:,1,3)
c    &                         + A3(:,4,1) * A1(:,1,3)
c    &		                                        )
c    &           + tau(:,2) * ( 
c    &                         + A1(:,4,2) * A3(:,2,3)
c    &                         + A1(:,4,3) * A3(:,3,3)
c    &                         + A1(:,4,4) * A3(:,4,3)
c    &                         + A3(:,4,2) * A1(:,2,3)
c    &                         + A3(:,4,3) * A1(:,3,3)
c    &                         + A3(:,4,4) * A1(:,4,3)
c    &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,4,5) * A3(:,5,3)
     &		               + A3(:,4,5) * A1(:,5,3)
     &		                                        )
c
          A0(:,4,4)= 
     &             tau(:,1) * (
     &                         + A1(:,4,1) * A3(:,1,4)
c    &                         + A3(:,4,1) * A1(:,1,4)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A3(:,2,4)
c    &                         + A1(:,4,3) * A3(:,3,4)
     &                         + A1(:,4,4) * A3(:,4,4)
c    &                         + A3(:,4,2) * A1(:,2,4)
c    &                         + A3(:,4,3) * A1(:,3,4)
     &                         + A3(:,4,4) * A1(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,4,5) * A3(:,5,4)
     &		               + A3(:,4,5) * A1(:,5,4)
     &		                                        )
c
          A0(:,4,5)= 
     &             tau(:,1) * (
     &                         + A1(:,4,1) * A3(:,1,5)
     &                         + A3(:,4,1) * A1(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,4,2) * A3(:,2,5)
c    &                         + A1(:,4,3) * A3(:,3,5)
     &                         + A1(:,4,4) * A3(:,4,5)
c    &                         + A3(:,4,2) * A1(:,2,5)
c    &                         + A3(:,4,3) * A1(:,3,5)
     &                         + A3(:,4,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,4,5) * A3(:,5,5)
     &		               + A3(:,4,5) * A1(:,5,5)
     &		                                        )
c
c  row five
c
          A0(:,5,1)= 
     &             tau(:,1) * (
     &                         + A1(:,5,1) * A3(:,1,1)
     &                         + A3(:,5,1) * A1(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,5,2) * A3(:,2,1)
     &                         + A1(:,5,3) * A3(:,3,1)
     &                         + A1(:,5,4) * A3(:,4,1)
     &                         + A3(:,5,2) * A1(:,2,1)
     &                         + A3(:,5,3) * A1(:,3,1)
     &                         + A3(:,5,4) * A1(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,5,5) * A3(:,5,1)
     &		               + A3(:,5,5) * A1(:,5,1)
     &		                                        )
c
          A0(:,5,2)= 
     &             tau(:,1) * (
c    &                         + A1(:,5,1) * A3(:,1,2)
     &                         + A3(:,5,1) * A1(:,1,2)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,5,2) * A3(:,2,2)
c    &                         + A1(:,5,3) * A3(:,3,2)
c    &                         + A1(:,5,4) * A3(:,4,2)
     &                         + A3(:,5,2) * A1(:,2,2)
     &                         + A3(:,5,3) * A1(:,3,2)
     &                         + A3(:,5,4) * A1(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,5,5) * A3(:,5,2)
     &		               + A3(:,5,5) * A1(:,5,2)
     &		                                        )
c
          A0(:,5,3)= 
c    &             tau(:,1) * (
c    &                         + A1(:,5,1) * A3(:,1,3)
c    &                         + A3(:,5,1) * A1(:,1,3)
c    &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A1(:,5,2) * A3(:,2,3)
     &                         + A1(:,5,3) * A3(:,3,3)
c    &                         + A1(:,5,4) * A3(:,4,3)
c    &                         + A3(:,5,2) * A1(:,2,3)
     &                         + A3(:,5,3) * A1(:,3,3)
c    &                         + A3(:,5,4) * A1(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,5,5) * A3(:,5,3)
     &		               + A3(:,5,5) * A1(:,5,3)
     &		                                        )
c
          A0(:,5,4)= 
     &             tau(:,1) * (
     &                         + A1(:,5,1) * A3(:,1,4)
c    &                         + A3(:,5,1) * A1(:,1,4)
     &		                                        )
     &           + tau(:,2) * (
     &                         + A1(:,5,2) * A3(:,2,4)
     &                         + A1(:,5,3) * A3(:,3,4)
     &                         + A1(:,5,4) * A3(:,4,4)
c    &                         + A3(:,5,2) * A1(:,2,4)
c    &                         + A3(:,5,3) * A1(:,3,4)
     &                         + A3(:,5,4) * A1(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,5,5) * A3(:,5,4)
     &		               + A3(:,5,5) * A1(:,5,4)
     &		                                        )
c
          A0(:,5,5)= 
     &             tau(:,1) * (
     &                         + A1(:,5,1) * A3(:,1,5)
     &                         + A3(:,5,1) * A1(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A1(:,5,2) * A3(:,2,5)
     &                         + A1(:,5,3) * A3(:,3,5)
     &                         + A1(:,5,4) * A3(:,4,5)
     &                         + A3(:,5,2) * A1(:,2,5)
     &                         + A3(:,5,3) * A1(:,3,5)
     &                         + A3(:,5,4) * A1(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A1(:,5,5) * A3(:,5,5)
     &		               + A3(:,5,5) * A1(:,5,5)
     &		                                        )
        else
           A0 = zero
        endif
c
        if (Navier .eq. 1) then
           A0(:,2,4) = A0(:,2,4) + rlm2mu - rmu
           A0(:,4,2) = A0(:,4,2) + rlm2mu - rmu
           A0(:,5,2) = A0(:,5,2) + (rlm2mu-rmu) * u3
           A0(:,5,4) = A0(:,5,4) + (rlm2mu-rmu) * u1
        endif
c
        do j = 1, nshl
          tmp = WdetJ * shg(:,j,1) * shg(:,j,3)
           do i=1,nflow
           do k=1,nflow
             BDiagl(:,j,i,k) = BDiagl(:,j,i,k) + tmp * A0(:,i,k)
           enddo
           enddo
         enddo
c
c now for the 2 3  plus 3 2
c
          if(ivart.ge.2) then
c
c  row one
c
          A0(:,1,1)= 
     &             tau(:,1) * (
     &                         + A2(:,1,1) * A3(:,1,1)
     &                         + A3(:,1,1) * A2(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,1,2) * A3(:,2,1)
     &                         + A2(:,1,3) * A3(:,3,1)
c    &                         + A2(:,1,4) * A3(:,4,1)
c    &                         + A3(:,1,2) * A2(:,2,1)
c    &                         + A3(:,1,3) * A2(:,3,1)
     &                         + A3(:,1,4) * A2(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,1,5) * A3(:,5,1)
     &		               + A3(:,1,5) * A2(:,5,1)
     &		                                        )
c
          A0(:,1,2)= 
c    &             tau(:,1) * (
c    &                         + A2(:,1,1) * A3(:,1,2)
c    &                         + A3(:,1,1) * A2(:,1,2)
c    &		                                        )
c    &           + tau(:,2) * ( 
c    &                         + A2(:,1,2) * A3(:,2,2)
c    &                         + A2(:,1,3) * A3(:,3,2)
c    &                         + A2(:,1,4) * A3(:,4,2)
c    &                         + A3(:,1,2) * A2(:,2,2)
c    &                         + A3(:,1,3) * A2(:,3,2)
c    &                         + A3(:,1,4) * A2(:,4,2)
c    &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,1,5) * A3(:,5,2)
     &		               + A3(:,1,5) * A2(:,5,2)
     &		                                        )
c
          A0(:,1,3)= 
     &             tau(:,1) * (
c    &                         + A2(:,1,1) * A3(:,1,3)
     &                         + A3(:,1,1) * A2(:,1,3)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,1,2) * A3(:,2,3)
     &                         + A2(:,1,3) * A3(:,3,3)
c    &                         + A2(:,1,4) * A3(:,4,3)
c    &                         + A3(:,1,2) * A2(:,2,3)
c    &                         + A3(:,1,3) * A2(:,3,3)
     &                         + A3(:,1,4) * A2(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,1,5) * A3(:,5,3)
     &		               + A3(:,1,5) * A2(:,5,3)
     &		                                        )
c
          A0(:,1,4)= 
     &             tau(:,1) * (
     &                         + A2(:,1,1) * A3(:,1,4)
c    &                         + A3(:,1,1) * A2(:,1,4)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,1,2) * A3(:,2,4)
     &                         + A2(:,1,3) * A3(:,3,4)
c    &                         + A2(:,1,4) * A3(:,4,4)
c    &                         + A3(:,1,2) * A2(:,2,4)
c    &                         + A3(:,1,3) * A2(:,3,4)
     &                         + A3(:,1,4) * A2(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,1,5) * A3(:,5,4)
     &		               + A3(:,1,5) * A2(:,5,4)
     &		                                        )
c
          A0(:,1,5)= 
     &             tau(:,1) * (
     &                         + A2(:,1,1) * A3(:,1,5)
     &                         + A3(:,1,1) * A2(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,1,2) * A3(:,2,5)
     &                         + A2(:,1,3) * A3(:,3,5)
c    &                         + A2(:,1,4) * A3(:,4,5)
c    &                         + A3(:,1,2) * A2(:,2,5)
c    &                         + A3(:,1,3) * A2(:,3,5)
     &                         + A3(:,1,4) * A2(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,1,5) * A3(:,5,5)
     &		               + A3(:,1,5) * A2(:,5,5)
     &		                                        )
c
c  row two
c
          A0(:,2,1)= 
     &             tau(:,1) * (
     &                         + A2(:,2,1) * A3(:,1,1)
     &                         + A3(:,2,1) * A2(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A2(:,2,2) * A3(:,2,1)
     &                         + A2(:,2,3) * A3(:,3,1)
c    &                         + A2(:,2,4) * A3(:,4,1)
     &                         + A3(:,2,2) * A2(:,2,1)
c    &                         + A3(:,2,3) * A2(:,3,1)
     &                         + A3(:,2,4) * A2(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,2,5) * A3(:,5,1)
     &		               + A3(:,2,5) * A2(:,5,1)
     &		                                        )
c
          A0(:,2,2)= 
c    &             tau(:,1) * (
c    &                         + A2(:,2,1) * A3(:,1,2)
c    &                         + A3(:,2,1) * A2(:,1,2)
c    &		                                        )
     &           + tau(:,2) * ( 
     &                         + A2(:,2,2) * A3(:,2,2)
c    &                         + A2(:,2,3) * A3(:,3,2)
c    &                         + A2(:,2,4) * A3(:,4,2)
     &                         + A3(:,2,2) * A2(:,2,2)
c    &                         + A3(:,2,3) * A2(:,3,2)
c    &                         + A3(:,2,4) * A2(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,2,5) * A3(:,5,2)
     &		               + A3(:,2,5) * A2(:,5,2)
     &		                                        )
c
          A0(:,2,3)= 
     &             tau(:,1) * (
c    &                         + A2(:,2,1) * A3(:,1,3)
     &                         + A3(:,2,1) * A2(:,1,3)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,2,2) * A3(:,2,3)
     &                         + A2(:,2,3) * A3(:,3,3)
c    &                         + A2(:,2,4) * A3(:,4,3)
     &                         + A3(:,2,2) * A2(:,2,3)
c    &                         + A3(:,2,3) * A2(:,3,3)
     &                         + A3(:,2,4) * A2(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,2,5) * A3(:,5,3)
     &		               + A3(:,2,5) * A2(:,5,3)
     &		                                        )
c
          A0(:,2,4)= 
     &             tau(:,1) * (
     &                         + A2(:,2,1) * A3(:,1,4)
c    &                         + A3(:,2,1) * A2(:,1,4)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A2(:,2,2) * A3(:,2,4)
     &                         + A2(:,2,3) * A3(:,3,4)
c    &                         + A2(:,2,4) * A3(:,4,4)
c    &                         + A3(:,2,2) * A2(:,2,4)
c    &                         + A3(:,2,3) * A2(:,3,4)
     &                         + A3(:,2,4) * A2(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,2,5) * A3(:,5,4)
     &		               + A3(:,2,5) * A2(:,5,4)
     &		                                        )
c
          A0(:,2,5)= 
     &             tau(:,1) * (
     &                         + A2(:,2,1) * A3(:,1,5)
     &                         + A3(:,2,1) * A2(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A2(:,2,2) * A3(:,2,5)
     &                         + A2(:,2,3) * A3(:,3,5)
c    &                         + A2(:,2,4) * A3(:,4,5)
     &                         + A3(:,2,2) * A2(:,2,5)
c    &                         + A3(:,2,3) * A2(:,3,5)
     &                         + A3(:,2,4) * A2(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,2,5) * A3(:,5,5)
     &		               + A3(:,2,5) * A2(:,5,5)
     &		                                        )
c
c  row three
c
          A0(:,3,1)= 
     &             tau(:,1) * (
     &                         + A2(:,3,1) * A3(:,1,1)
     &                         + A3(:,3,1) * A2(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,3,2) * A3(:,2,1)
     &                         + A2(:,3,3) * A3(:,3,1)
c    &                         + A2(:,3,4) * A3(:,4,1)
c    &                         + A3(:,3,2) * A2(:,2,1)
     &                         + A3(:,3,3) * A2(:,3,1)
     &                         + A3(:,3,4) * A2(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,3,5) * A3(:,5,1)
     &		               + A3(:,3,5) * A2(:,5,1)
     &		                                        )
c
          A0(:,3,2)= 
c    &             tau(:,1) * (
c    &                         + A2(:,3,1) * A3(:,1,2)
c    &                         + A3(:,3,1) * A2(:,1,2)
c    &		                                        )
c    &           + tau(:,2) * ( 
c    &                         + A2(:,3,2) * A3(:,2,2)
c    &                         + A2(:,3,3) * A3(:,3,2)
c    &                         + A2(:,3,4) * A3(:,4,2)
c    &                         + A3(:,3,2) * A2(:,2,2)
c    &                         + A3(:,3,3) * A2(:,3,2)
c    &                         + A3(:,3,4) * A2(:,4,2)
c    &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,3,5) * A3(:,5,2)
     &		               + A3(:,3,5) * A2(:,5,2)
     &		                                        )
c
          A0(:,3,3)= 
     &             tau(:,1) * (
c    &                         + A2(:,3,1) * A3(:,1,3)
     &                         + A3(:,3,1) * A2(:,1,3)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,3,2) * A3(:,2,3)
     &                         + A2(:,3,3) * A3(:,3,3)
c    &                         + A2(:,3,4) * A3(:,4,3)
c    &                         + A3(:,3,2) * A2(:,2,3)
     &                         + A3(:,3,3) * A2(:,3,3)
     &                         + A3(:,3,4) * A2(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,3,5) * A3(:,5,3)
     &		               + A3(:,3,5) * A2(:,5,3)
     &		                                        )
c
          A0(:,3,4)= 
     &             tau(:,1) * (
     &                         + A2(:,3,1) * A3(:,1,4)
c    &                         + A3(:,3,1) * A2(:,1,4)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,3,2) * A3(:,2,4)
     &                         + A2(:,3,3) * A3(:,3,4)
c    &                         + A2(:,3,4) * A3(:,4,4)
c    &                         + A3(:,3,2) * A2(:,2,4)
c    &                         + A3(:,3,3) * A2(:,3,4)
     &                         + A3(:,3,4) * A2(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,3,5) * A3(:,5,4)
     &		               + A3(:,3,5) * A2(:,5,4)
     &		                                        )
c
          A0(:,3,5)= 
     &             tau(:,1) * (
     &                         + A2(:,3,1) * A3(:,1,5)
     &                         + A3(:,3,1) * A2(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,3,2) * A3(:,2,5)
     &                         + A2(:,3,3) * A3(:,3,5)
c    &                         + A2(:,3,4) * A3(:,4,5)
c    &                         + A3(:,3,2) * A2(:,2,5)
     &                         + A3(:,3,3) * A2(:,3,5)
     &                         + A3(:,3,4) * A2(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,3,5) * A3(:,5,5)
     &		               + A3(:,3,5) * A2(:,5,5)
     &		                                        )
c
c  row four
c
          A0(:,4,1)= 
     &             tau(:,1) * (
     &                         + A2(:,4,1) * A3(:,1,1)
     &                         + A3(:,4,1) * A2(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,4,2) * A3(:,2,1)
     &                         + A2(:,4,3) * A3(:,3,1)
     &                         + A2(:,4,4) * A3(:,4,1)
c    &                         + A3(:,4,2) * A2(:,2,1)
c    &                         + A3(:,4,3) * A2(:,3,1)
     &                         + A3(:,4,4) * A2(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,4,5) * A3(:,5,1)
     &		               + A3(:,4,5) * A2(:,5,1)
     &		                                        )
c
          A0(:,4,2)= 
c    &             tau(:,1) * (
c    &                         + A2(:,4,1) * A3(:,1,2)
c    &                         + A3(:,4,1) * A2(:,1,2)
c    &		                                        )
c    &           + tau(:,2) * ( 
c    &                         + A2(:,4,2) * A3(:,2,2)
c    &                         + A2(:,4,3) * A3(:,3,2)
c    &                         + A2(:,4,4) * A3(:,4,2)
c    &                         + A3(:,4,2) * A2(:,2,2)
c    &                         + A3(:,4,3) * A2(:,3,2)
c    &                         + A3(:,4,4) * A2(:,4,2)
c    &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,4,5) * A3(:,5,2)
     &		               + A3(:,4,5) * A2(:,5,2)
     &		                                        )
c
          A0(:,4,3)= 
     &             tau(:,1) * (
c    &                         + A2(:,4,1) * A3(:,1,3)
     &                         + A3(:,4,1) * A2(:,1,3)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,4,2) * A3(:,2,3)
     &                         + A2(:,4,3) * A3(:,3,3)
c    &                         + A2(:,4,4) * A3(:,4,3)
c    &                         + A3(:,4,2) * A2(:,2,3)
c    &                         + A3(:,4,3) * A2(:,3,3)
     &                         + A3(:,4,4) * A2(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,4,5) * A3(:,5,3)
     &		               + A3(:,4,5) * A2(:,5,3)
     &		                                        )
c
          A0(:,4,4)= 
     &             tau(:,1) * (
     &                         + A2(:,4,1) * A3(:,1,4)
c    &                         + A3(:,4,1) * A2(:,1,4)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,4,2) * A3(:,2,4)
     &                         + A2(:,4,3) * A3(:,3,4)
     &                         + A2(:,4,4) * A3(:,4,4)
c    &                         + A3(:,4,2) * A2(:,2,4)
c    &                         + A3(:,4,3) * A2(:,3,4)
     &                         + A3(:,4,4) * A2(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,4,5) * A3(:,5,4)
     &		               + A3(:,4,5) * A2(:,5,4)
     &		                                        )
c
          A0(:,4,5)= 
     &             tau(:,1) * (
     &                         + A2(:,4,1) * A3(:,1,5)
     &                         + A3(:,4,1) * A2(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,4,2) * A3(:,2,5)
     &                         + A2(:,4,3) * A3(:,3,5)
     &                         + A2(:,4,4) * A3(:,4,5)
c    &                         + A3(:,4,2) * A2(:,2,5)
c    &                         + A3(:,4,3) * A2(:,3,5)
     &                         + A3(:,4,4) * A2(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,4,5) * A3(:,5,5)
     &		               + A3(:,4,5) * A2(:,5,5)
     &		                                        )
c
c  row five
c
          A0(:,5,1)= 
     &             tau(:,1) * (
     &                         + A2(:,5,1) * A3(:,1,1)
     &                         + A3(:,5,1) * A2(:,1,1)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A2(:,5,2) * A3(:,2,1)
     &                         + A2(:,5,3) * A3(:,3,1)
     &                         + A2(:,5,4) * A3(:,4,1)
     &                         + A3(:,5,2) * A2(:,2,1)
     &                         + A3(:,5,3) * A2(:,3,1)
     &                         + A3(:,5,4) * A2(:,4,1)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,5,5) * A3(:,5,1)
     &		               + A3(:,5,5) * A2(:,5,1)
     &		                                        )
c
          A0(:,5,2)= 
c    &             tau(:,1) * (
c    &                         + A2(:,5,1) * A3(:,1,2)
c    &                         + A3(:,5,1) * A2(:,1,2)
c    &		                                        )
     &           + tau(:,2) * ( 
     &                         + A2(:,5,2) * A3(:,2,2)
c    &                         + A2(:,5,3) * A3(:,3,2)
c    &                         + A2(:,5,4) * A3(:,4,2)
     &                         + A3(:,5,2) * A2(:,2,2)
c    &                         + A3(:,5,3) * A2(:,3,2)
c    &                         + A3(:,5,4) * A2(:,4,2)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,5,5) * A3(:,5,2)
     &		               + A3(:,5,5) * A2(:,5,2)
     &		                                        )
c
          A0(:,5,3)= 
     &             tau(:,1) * (
c    &                         + A2(:,5,1) * A3(:,1,3)
     &                         + A3(:,5,1) * A2(:,1,3)
     &		                                        )
     &           + tau(:,2) * ( 
c    &                         + A2(:,5,2) * A3(:,2,3)
     &                         + A2(:,5,3) * A3(:,3,3)
c    &                         + A2(:,5,4) * A3(:,4,3)
     &                         + A3(:,5,2) * A2(:,2,3)
     &                         + A3(:,5,3) * A2(:,3,3)
     &                         + A3(:,5,4) * A2(:,4,3)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,5,5) * A3(:,5,3)
     &		               + A3(:,5,5) * A2(:,5,3)
     &		                                        )
c
          A0(:,5,4)= 
     &             tau(:,1) * (
     &                         + A2(:,5,1) * A3(:,1,4)
c    &                         + A3(:,5,1) * A2(:,1,4)
     &		                                        )
     &           + tau(:,2) * (
     &                         + A2(:,5,2) * A3(:,2,4)
     &                         + A2(:,5,3) * A3(:,3,4)
     &                         + A2(:,5,4) * A3(:,4,4)
c    &                         + A3(:,5,2) * A2(:,2,4)
c    &                         + A3(:,5,3) * A2(:,3,4)
     &                         + A3(:,5,4) * A2(:,4,4)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,5,5) * A3(:,5,4)
     &		               + A3(:,5,5) * A2(:,5,4)
     &		                                        )
c
          A0(:,5,5)= 
     &             tau(:,1) * (
     &                         + A2(:,5,1) * A3(:,1,5)
     &                         + A3(:,5,1) * A2(:,1,5)
     &		                                        )
     &           + tau(:,2) * ( 
     &                         + A2(:,5,2) * A3(:,2,5)
     &                         + A2(:,5,3) * A3(:,3,5)
     &                         + A2(:,5,4) * A3(:,4,5)
     &                         + A3(:,5,2) * A2(:,2,5)
     &                         + A3(:,5,3) * A2(:,3,5)
     &                         + A3(:,5,4) * A2(:,4,5)
     &                                                  )
     &           + tau(:,3) * (
     &		               + A2(:,5,5) * A3(:,5,5)
     &		               + A3(:,5,5) * A2(:,5,5)
     &		                                        )
        else
           A0 = zero
        endif
c
        if (Navier .eq. 1) then
           A0(:,3,4) = A0(:,3,4) + rlm2mu - rmu
           A0(:,4,3) = A0(:,4,3) + rlm2mu - rmu
           A0(:,5,3) = A0(:,5,3) + (rlm2mu-rmu) * u3
           A0(:,5,4) = A0(:,5,4) + (rlm2mu-rmu) * u2
        endif

        do j = 1, nshl
          tmp = WdetJ * shg(:,j,2) * shg(:,j,3)
           do i=1,nflow
           do k=1,nflow
             BDiagl(:,j,i,k) = BDiagl(:,j,i,k) + tmp * A0(:,i,k)
           enddo
           enddo
        enddo

	ttim(30) = ttim(30) + tmr()
 
        return
        end


        subroutine e3bdgSclr  (	shp, 	shg, 	WdetJ,
     &				A1t, 	A2t, 	A3t, 	A0t, 
     &                        	taut,
     &   			Sclr, 
     &				u1,  	u2,  	u3,  
     &				rmu,    BDiagtl)

        include "common.h"

c
c  passed arrays
c
        dimension BDiagtl(npro,nshl),  
     &            shp(npro,nshl),      Sclr(npro),
     &            shg(npro,nshl,nsd),  WdetJ(npro),
     &            A1tautA0(npro),      A2tautA0(npro),
     &            A3tautA0(npro),      Ataut(npro),
     &            A1t(npro),           A2t(npro),
     &            A3t(npro),           A0t(npro),
     &            taut(npro),          u1(npro),
     &            u2(npro),            u3(npro),
     &            rmu(npro)
c
c
c  passed work arrays for local variables
c
	dimension tmp1(npro),             tmp2(npro)              
	dimension tmp(npro)

	ttim(30) = ttim(30) - tmr()

c $$ Ex-E3conv
c
c.... calculate the contribution in non-integrated by part form
c.... add (N_b A_tilde_i N_b,i) to BDiagl
c
       do j = 1, nshl   ! May be worth eliminating zeros in A(prim) matrices
        tmp=shp(:,j)*WdetJ
            BDiagtl(:,j) = BDiagtl(:,j)
     &                          + tmp * (
     &                                   shg(:,j,1) * A1t(:) +
     &                                   shg(:,j,2) * A2t(:) +
     &                                   shg(:,j,3) * A3t(:) 
     &                                                            )

       enddo
        if (ngauss .eq. 1 .and. nshl .eq. 4) then  ! Exact integ of mass term for tets

c $$ Ex-e3mael
c
c.... calculate factor   V             * sum(column of M)
c                     =WdetJ/(3/4)Qwt(lcsyst,intp)* (2+1+1+1)/20
c
	tmp =  0.4d0 * WdetJ * Dtgl/Qwt(lcsyst,intp)/three*almi/gami/alfi
c
        do j=1,nshl  ! take advantage of zeros in A0t(Prim)
	  BDiagtl(:,j) = BDiagtl(:,j) + tmp * A0t(:)
	enddo

        else

c $$ Ex-e3mass
c
        ff = almi / gami / alfi
        tmp = WdetJ * Dtgl * ff
c
        do j = 1, nshl
c
          tmp2 = (shp(:,j)*shp(:,j)) * tmp   
c
	  BDiagtl(:,j) = BDiagtl(:,j) + tmp2 * A0t(:)

        enddo

        endif

c
c.... calculate (Atau) <-- (A_1 tau) 
c                                   
c     
       if (ivart .ge. 2) then
       do i = 1, nflow
          Ataut(:) = A1t(:)*taut(:)
       enddo
c
c.... calculate (A_1 tau A_0) (for L.S. time term of EGmass)
c
             A1tautA0(:) = Ataut(:)*A0t(:) 
c
c.... calculate (Atau) <-- (A_2 tau) 
c                                           
c     
          Ataut(:) = A2t(:)*taut(:)

c
c.... calculate (A_2 tau A_0) (for L.S. time term of EGmass)
c

             A2tautA0(:) = Ataut(:)*A0t(:)
c
c.... calculate (Atau) <-- (A_3 tau)  
c                                           
c     

          Ataut(:) = A3t(:)*taut(:)
c
c.... calculate (A_3 tau A_0) (for L.S. time term of EGmass)
c

             A3tautA0(:) = Ataut(:)*A0t(:) 


c.... add least squares time term to BDiagl
c
c
c.... loop through rows (nodes i)
c
       do i = 1, nshl
c
c.... loop through column nodes, add (N_a,i A_i tau N_b) to BDiagl
c
          tmp1 = shp(:,i) * WdetJ * almi/gami/alfi*dtgl

                BDiagtl(:,i) = BDiagtl(:,i)  + tmp1*(
     &               shg(:,i,1) * A1tautA0(:) +
     &               shg(:,i,2) * A2tautA0(:) +
     &               shg(:,i,3) * A3tautA0(:))

       enddo
       endif
       
c
c $$ Ex-e3wmlt
c
c.... add contribution of stiffness to BDiagl
c
c....  we no longer build stiff so we have to get it in here.
c      recall that this is the (A_i tau A_j + K_{ij}) that
c      multiplies N_{a,i} N_{b,j} (Actually for BDiag a=b).
c
c...   we are through with A0 so use it now for the sub-blocks
c      of what used to be stiff
c
c
     
c
c  start with 1 1 contribution
c
        if(ivart .ge. 2) then 


          A0t(:)= A1t(:)*taut(:) *A1t(:)
        else
          A0t=zero
        endif
c
        if (Navier .eq. 1) then
c.... K11
c
           A0t(:) = A0t(:) +1.5*(rmu+Sclr)

        endif
        do j = 1, nshl
          tmp = WdetJ * shg(:,j,1) * shg(:,j,1)
  
             BDiagtl(:,j) = BDiagtl(:,j) + tmp * A0t(:)
        enddo
c
c  start with 2 2 contribution
c
        if(ivart .ge. 2) then 


          A0t(:)= A2t(:)*taut(:) *A2t(:)
        else
          A0t=zero
        endif
c
        if (Navier .eq. 1) then
c.... K22
c
           A0t(:) = A0t(:) +1.5*(rmu+Sclr)

        endif
        do j = 1, nshl
          tmp = WdetJ * shg(:,j,2) * shg(:,j,2)
  
             BDiagtl(:,j) = BDiagtl(:,j) + tmp * A0t(:)
        enddo

c
c  start with 3 3 contribution
c
        if(ivart .ge. 2) then 


          A0t(:)= A3t(:)*taut(:) *A3t(:)
        else
          A0t=zero
        endif
c
        if (Navier .eq. 1) then
c.... K33
c
           A0t(:) = A0t(:) +1.5*(rmu+Sclr)

        endif
        do j = 1, nshl
          tmp = WdetJ * shg(:,j,3) * shg(:,j,3)
  
             BDiagtl(:,j) = BDiagtl(:,j) + tmp * A0t(:)
        enddo



c
c now for the 1 2  plus 2 1
c
        if(ivart.ge.2) then

          A0t(:)= taut(:) * two * A1t(:) * A2t(:)

        else
           A0t = zero
        endif

        do j = 1, nshl
          tmp = WdetJ * shg(:,j,1) * shg(:,j,2)

             BDiagtl(:,j) = BDiagtl(:,j) + tmp * A0t(:)
        enddo
c
c now for the 1 3  plus 3 1
c
        if(ivart.ge.2) then
          A0t(:)= taut(:) * two * A1t(:) * A3t(:)

        else
           A0t = zero
        endif

        do j = 1, nshl
          tmp = WdetJ * shg(:,j,1) * shg(:,j,3)

             BDiagtl(:,j) = BDiagtl(:,j) + tmp * A0t(:)
        enddo
c
c now for the 2 3  plus 3 2
c
        if(ivart.ge.2) then

          A0t(:)= taut(:) * two * A2t(:) * A3t(:)

        else
           A0t = zero
        endif

        do j = 1, nshl
          tmp = WdetJ * shg(:,j,2) * shg(:,j,3)

             BDiagtl(:,j) = BDiagtl(:,j) + tmp * A0t(:)
        enddo


	ttim(30) = ttim(30) + tmr()
 
        return
        end
