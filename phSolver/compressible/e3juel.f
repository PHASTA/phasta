        subroutine e3juel (yl,     acl,    sls,    A0,    
     &			   WdetJ,  rl,     rml)
c
c----------------------------------------------------------------------
c
c This routine calculates Exactly integrated Linear Tetrahedra
c Mass term (assuming U(Y) is linear it is not).
c
c input:
c  WdetJ  (npro)                : weighted Jacobian
c
c output:
c  rl     (npro,nshl,nflow)      : residual
c  rml    (npro,nshl,nflow)      : modified residual
c
c
c  note that this routine wipes out yl by putting ul into it
c  and then (in ires=1 case ) it is used again
c
c Kenneth Jansen, Winter 1997, Primitive Variables
c----------------------------------------------------------------------
c
	include "common.h"
c
        dimension yl(npro,nshl,nflow),        acl(npro,nshl,ndof),
     &            WdetJ(npro),               A0(npro,nflow,nflow),
     &            rl(npro,nshl,nflow),        rml(npro,nshl,nflow)
c
        dimension rk(npro),                  rho(npro),
     &            ei(npro),                  tmp(npro),
     &            ub(npro,nflow),             fact(npro),
     &		  fddt(npro)

	ttim(28) = ttim(28) - secs(0.0)

c
c.... --------------------->  Time term   <--------------------
c
c.... add contribution of U to rml
c
c.... compute conservative variables
c
c
c   multiply by exact mass matrix   integral N_aN_b=(I+1)*V/20
c   where 1 is a matrix with every element=1
c
c   note that the wght has 4/3 multiplier so 3/4*20=15
c

        fact=WdetJ/(Qwt(lcsyst,intp)*15.0d0)
         fct1=almi/(gami*alfi)*Dtgl  ! factor for predictor (scalar)
        if(ires.ne.1) then
         fddt=fact*fct1
         do inod=1,nshl
	  rk = pt5 * (yl(:,inod,2)**2 + yl(:,inod,3)**2 + yl(:,inod,4)**2)
c
          ithm = 6
          call getthm (yl(:,inod,1), yl(:,inod,5), sls,
     &                 rk,           rho,          ei,
     &                 tmp,          tmp,          tmp,
     &                 tmp,          tmp,          tmp,
     &                 tmp,          tmp)
c
         yl(:,inod,1) = rho 
         yl(:,inod,2) = rho * yl(:,inod,2)
         yl(:,inod,3) = rho * yl(:,inod,3)
         yl(:,inod,4) = rho * yl(:,inod,4)
         yl(:,inod,5) = rho * (ei + rk)
        enddo
	    ub(:,:)=yl(:,1,:)+yl(:,2,:)
     &             +yl(:,3,:)+yl(:,4,:)

c
c  what we have now in yl is the U_b^e
c  we want to get Resm(mass)=M^e_{ab} (U^e_b(Y+epsilon P) 
c                                              - U^e_b(Y))*fact/dt,
c  since only the difference between resm's is important we
c  do not have to subtract off the unperturbed U(Y) vector
c  This term is meant to carry through the effect of a perturbation
c  in Y upon dY/dt (through the predictor into the term fact)
c
c

	  do i = 1, nshl
	    rml(:,i,1) = rml(:,i,1) + fddt * (yl(:,i,1)+ub(:,1))
	    rml(:,i,2) = rml(:,i,2) + fddt * (yl(:,i,2)+ub(:,2))
	    rml(:,i,3) = rml(:,i,3) + fddt * (yl(:,i,3)+ub(:,3))
	    rml(:,i,4) = rml(:,i,4) + fddt * (yl(:,i,4)+ub(:,4))
	    rml(:,i,5) = rml(:,i,5) + fddt * (yl(:,i,5)+ub(:,5))
	  enddo
c
!      flops = flops + 35*nshl*npro
	endif

c
c.... ires = 2 or 3
c
	if ((ires .eq. 1) .or. (ires .eq. 3)) then

	    ub(:,:)=acl(:,1,:)+acl(:,2,:)
     &             +acl(:,3,:)+acl(:,4,:)

	  do i = 1, nshl
	    yl(:,i,1) = fact*(acl(:,i,1)+ub(:,1))
	    yl(:,i,2) = fact*(acl(:,i,2)+ub(:,2))
	    yl(:,i,3) = fact*(acl(:,i,3)+ub(:,3))
	    yl(:,i,4) = fact*(acl(:,i,4)+ub(:,4))
	    yl(:,i,5) = fact*(acl(:,i,5)+ub(:,5))
	  enddo
c
c  what we have now in yl is the dY_a^e/dt=M^e_{ab} Y_{b,t},  must multiply by
c  A0 to get dU^e_a/dt= dU/dY(centroid) dY_a^e/dt,  take advantage of zeros 
c  in A0(Prim) with comments
c
          do i = 1, nshl
            rl(:,i,1) = rl(:,i,1) 
     &     + A0(:,1,1)*yl(:,i,1)
c    &     + A0(:,1,2)*yl(:,i,2)
c    &     + A0(:,1,3)*yl(:,i,3)
c    &     + A0(:,1,4)*yl(:,i,4)
     &     + A0(:,1,5)*yl(:,i,5)
c
            rl(:,i,2) = rl(:,i,2) 
     &     + A0(:,2,1)*yl(:,i,1)
     &     + A0(:,2,2)*yl(:,i,2)
c    &     + A0(:,2,3)*yl(:,i,3)
c    &     + A0(:,2,4)*yl(:,i,4)
     &     + A0(:,2,5)*yl(:,i,5)
c
            rl(:,i,3) = rl(:,i,3) 
     &     + A0(:,3,1)*yl(:,i,1)
c    &     + A0(:,3,2)*yl(:,i,2)
     &     + A0(:,3,3)*yl(:,i,3)
c    &     + A0(:,3,4)*yl(:,i,4)
     &     + A0(:,3,5)*yl(:,i,5)
c
            rl(:,i,4) = rl(:,i,4) 
     &     + A0(:,4,1)*yl(:,i,1)
c    &     + A0(:,4,2)*yl(:,i,2)
c    &     + A0(:,4,3)*yl(:,i,3)
     &     + A0(:,4,4)*yl(:,i,4)
     &     + A0(:,4,5)*yl(:,i,5)
c
            rl(:,i,5) = rl(:,i,5) 
     &     + A0(:,5,1)*yl(:,i,1)
     &     + A0(:,5,2)*yl(:,i,2)
     &     + A0(:,5,3)*yl(:,i,3)
     &     + A0(:,5,4)*yl(:,i,4)
     &     + A0(:,5,5)*yl(:,i,5)
c
          enddo
c
!      flops = flops + 45*nshl*npro
	endif

	ttim(28) = ttim(28) + tmr()
c
c.... return
c
	return
	end

c$$$
c$$$        subroutine e3juelSclr (ycl,     acl,    A0t,    
c$$$     &			       WdetJ,  rtl,     rmtl)
c$$$c
c$$$c----------------------------------------------------------------------
c$$$c
c$$$c This routine calculates Exactly integrated Linear Tetrahedra
c$$$c Mass term (assuming U(Y) is linear it is not).
c$$$c
c$$$c input:
c$$$c  WdetJ  (npro)                : weighted Jacobian
c$$$c
c$$$c output:
c$$$c  rtl     (npro,nshl,nflow)      : residual
c$$$c  rmtl    (npro,nshl,nflow)      : modified residual
c$$$c
c$$$c
c$$$c  note that this routine wipes out ycl by putting ul into it
c$$$c  and then (in ires=1 case ) it is used again
c$$$c
c$$$c Kenneth Jansen, Winter 1997, Primitive Variables
c$$$c----------------------------------------------------------------------
c$$$c
c$$$	include "common.h"
c$$$c
c$$$        dimension ycl(npro,nshl,ndof),    acl(npro,nshl,ndof),
c$$$     &            WdetJ(npro),           A0t(npro),
c$$$     &            rtl(npro,nshl),        rmtl(npro,nshl)
c$$$
c$$$c
c$$$        dimension rk(npro),              rho(npro),
c$$$     &            ei(npro),              tmp(npro),
c$$$     &            ubt(npro),             fact(npro),
c$$$     &		  fddt(npro)
c$$$
c$$$	ttim(28) = ttim(28) - tmr()
c$$$
c$$$c
c$$$c.... --------------------->  Time term   <--------------------
c$$$c
c$$$c.... add contribution of U to rml
c$$$c
c$$$c.... compute conservative variables
c$$$c
c$$$c
c$$$c   multiply by exact mass matrix   integral N_aN_b=(I+1)*V/20
c$$$c   where 1 is a matrix with every element=1
c$$$c
c$$$c   note that the wght has 4/3 multiplier so 3/4*20=15
c$$$c
c$$$
c$$$        fact=WdetJ/(Qwt(lcsyst,intp)*15.0d0)
c$$$         fct1=almi/(gami*alfi)*Dtgl  ! factor for predictor (scaler)
c$$$c
c$$$c
c$$$c.... ires = 2 or 3
c$$$c
c$$$
c$$$	if ((ires .eq. 1) .or. (ires .eq. 3)) then
c$$$	    ubt(:)=acl(:,1,id)+acl(:,2,id)
c$$$     &            +acl(:,3,id)+acl(:,4,id)
c$$$	  do i = 1, nshl
c$$$	    ycl(:,i,id) = fact*(acl(:,i,id)+ubt(:))
c$$$
c$$$	  enddo
c$$$c
c$$$c  what we have now in ycl is the dY_a^e/dt=M^e_{ab} Y_{b,t},  must multiply by
c$$$c  A0 to get dU^e_a/dt= dU/dY(centroid) dY_a^e/dt,  take advantage of zeros 
c$$$c  in A0(Prim) with comments
c$$$c
c$$$          do i = 1, nshl
c$$$            rtl(:,i) = rtl(:,i) 
c$$$     &     + A0t(:)*ycl(:,i,id)
c$$$
c$$$c
c$$$          enddo
c$$$c
c$$$     !      flops = flops + 45*nenl*npro
c$$$	endif
c$$$
c$$$	ttim(28) = ttim(28) + tmr()
c$$$c
c$$$c.... return
c$$$c
c$$$	return
c$$$	end
