c$$$	subroutine e3DC (g1yi,   g2yi,   g3yi,  
c$$$     &                   giju,   tauM,    A0,     raLS,
c$$$     &			 rtLS,   giju,   DC,     ri,
c$$$     &                   rmi,    stiff, A0DC)
c$$$c
c$$$c----------------------------------------------------------------------
c$$$c
c$$$c This routine calculates the contribution of the Discontinuity-
c$$$c Capturing operator to RHS and preconditioner.
c$$$c
c$$$c  g1yi   (nflow,npro)           : grad-y in direction 1
c$$$c  g2yi   (nflow,npro)           : grad-y in direction 2
c$$$c  g3yi   (nflow,npro)           : grad-y in direction 3
c$$$c  A0     (nsymdf,npro)          : A0 matrix (Symm. storage)
c$$$c  raLS   (npro)                 : square of LS residual (A0inv norm)
c$$$c  rtLS   (npro)                 : square of LS residual (Tau norm)
c$$$c  giju    (6,npro)              : metric matrix
c$$$c  DC     (ngauss,npro)          : discontinuity-capturing factor
c$$$c  intp				 : integration point number
c$$$c
c$$$c output:
c$$$c  ri     (nflow*(nsd+1),npro)   : partial residual
c$$$c  rmi    (nflow*(nsd+1),npro)   : partial modified residual
c$$$c  stiff  (nsymdf,6,npro)       : diffusivity matrix
c$$$c  DC     (npro)                : discontinuity-capturing factor
c$$$c
c$$$c
c$$$c Zdenek Johan, Summer 1990. (Modified from e2dc.f)
c$$$c Zdenek Johan, Winter 1991. (Recoded)
c$$$c Zdenek Johan, Winter 1991. (Fortran 90)
c$$$c----------------------------------------------------------------------
c$$$c
c$$$	include "common.h"
c$$$c
c$$$        dimension g1yi(npro,nflow),          g2yi(npro,nflow),
c$$$     &            g3yi(npro,nflow),          A0(npro,5,5),
c$$$     &            raLS(npro),                rtLS(npro),
c$$$     &            giju(npro,6),              DC(npro,ngauss),
c$$$     &            ri(npro,nflow*(nsd+1)),    rmi(npro,nflow*(nsd+1)),
c$$$     &            stiff(npro,3*nflow,3*nflow),itmp(npro)
c$$$c
c$$$
c$$$        dimension ggyi(npro,nflow),         gAgyi(npro,15),
c$$$     &            gnorm(npro),              A0gyi(npro,15),
c$$$     &            yiA0DCyj(npro,6),         A0DC(npro,4)
c$$$c
c$$$c ... -----------------------> initialize <----------------------------
c$$$c
c$$$        A0gyi    = zero
c$$$        gAgyi    = zero
c$$$        yiA0DCyj = zero
c$$$        DC       = zero
c$$$c.... ----------------------->  global gradient  <----------------------
c$$$c
c$$$c.... calculate (A0 y_,j) --> A0gyi
c$$$c
c$$$c  A0 Y_{,1}
c$$$c
c$$$        A0gyi( :,1) = A0(:,1,1)*g1yi(:,1)
c$$$     &              + A0(:,1,2)*g1yi(:,2)
c$$$     &              + A0(:,1,3)*g1yi(:,3)
c$$$     &              + A0(:,1,4)*g1yi(:,4)
c$$$     &              + A0(:,1,5)*g1yi(:,5)      
c$$$        A0gyi( :,2) = A0(:,2,1)*g1yi(:,1)
c$$$     &              + A0(:,2,2)*g1yi(:,2)
c$$$     &              + A0(:,2,3)*g1yi(:,3)
c$$$     &              + A0(:,2,4)*g1yi(:,4)
c$$$     &              + A0(:,2,5)*g1yi(:,5)
c$$$        A0gyi( :,3) = A0(:,3,1)*g1yi(:,1)
c$$$     &              + A0(:,3,2)*g1yi(:,2)
c$$$     &              + A0(:,3,3)*g1yi(:,3)
c$$$     &              + A0(:,3,4)*g1yi(:,4)
c$$$     &              + A0(:,3,5)*g1yi(:,5)
c$$$        A0gyi( :,4) = A0(:,4,1)*g1yi(:,1)
c$$$     &              + A0(:,4,2)*g1yi(:,2)
c$$$     &              + A0(:,4,3)*g1yi(:,3)
c$$$     &              + A0(:,4,4)*g1yi(:,4)
c$$$     &              + A0(:,4,5)*g1yi(:,5)
c$$$        A0gyi( :,5) = A0(:,5,1)*g1yi(:,1)
c$$$     &              + A0(:,5,2)*g1yi(:,2)
c$$$     &              + A0(:,5,3)*g1yi(:,3)
c$$$     &              + A0(:,5,4)*g1yi(:,4)
c$$$     &              + A0(:,5,5)*g1yi(:,5)
c$$$c
c$$$c  A0 Y_{,2}
c$$$c
c$$$        A0gyi( :,6) = A0(:,1,1)*g2yi(:,1)
c$$$     &              + A0(:,1,2)*g2yi(:,2)
c$$$     &              + A0(:,1,3)*g2yi(:,3)
c$$$     &              + A0(:,1,4)*g2yi(:,4)
c$$$     &              + A0(:,1,5)*g2yi(:,5)
c$$$        A0gyi( :,7) = A0(:,2,1)*g2yi(:,1)
c$$$     &              + A0(:,2,2)*g2yi(:,2)
c$$$     &              + A0(:,2,3)*g2yi(:,3)
c$$$     &              + A0(:,2,4)*g2yi(:,4)
c$$$     &              + A0(:,2,5)*g2yi(:,5)
c$$$        A0gyi( :,8) = A0(:,3,1)*g2yi(:,1)
c$$$     &              + A0(:,3,2)*g2yi(:,2)
c$$$     &              + A0(:,3,3)*g2yi(:,3)
c$$$     &              + A0(:,3,4)*g2yi(:,4)
c$$$     &              + A0(:,3,5)*g2yi(:,5)
c$$$        A0gyi( :,9) = A0(:,4,1)*g2yi(:,1)
c$$$     &              + A0(:,4,2)*g2yi(:,2)
c$$$     &              + A0(:,4,3)*g2yi(:,3)
c$$$     &              + A0(:,4,4)*g2yi(:,4)
c$$$     &              + A0(:,4,5)*g2yi(:,5)
c$$$        A0gyi(:,10) = A0(:,5,1)*g2yi(:,1)
c$$$     &              + A0(:,5,2)*g2yi(:,2)
c$$$     &              + A0(:,5,3)*g2yi(:,3)
c$$$     &              + A0(:,5,4)*g2yi(:,4)
c$$$     &              + A0(:,5,5)*g2yi(:,5)
c$$$c
c$$$c  A0 Y_{,3}
c$$$c
c$$$        A0gyi(:,11) = A0(:,1,1)*g3yi(:,1)
c$$$     &              + A0(:,1,2)*g3yi(:,2)
c$$$     &              + A0(:,1,3)*g3yi(:,3)
c$$$     &              + A0(:,1,4)*g3yi(:,4)
c$$$     &              + A0(:,1,5)*g3yi(:,5)
c$$$        A0gyi(:,12) = A0(:,2,1)*g3yi(:,1)
c$$$     &              + A0(:,2,2)*g3yi(:,2)
c$$$     &              + A0(:,2,3)*g3yi(:,3)
c$$$     &              + A0(:,2,4)*g3yi(:,4)
c$$$     &              + A0(:,2,5)*g3yi(:,5)
c$$$        A0gyi(:,13) = A0(:,3,1)*g3yi(:,1)
c$$$     &              + A0(:,3,2)*g3yi(:,2)
c$$$     &              + A0(:,3,3)*g3yi(:,3)
c$$$     &              + A0(:,3,4)*g3yi(:,4)
c$$$     &              + A0(:,3,5)*g3yi(:,5)
c$$$        A0gyi(:,14) = A0(:,4,1)*g3yi(:,1)
c$$$     &              + A0(:,4,2)*g3yi(:,2)
c$$$     &              + A0(:,4,3)*g3yi(:,3)
c$$$     &              + A0(:,4,4)*g3yi(:,4)
c$$$     &              + A0(:,4,5)*g3yi(:,5)
c$$$        A0gyi(:,15) = A0(:,5,1)*g3yi(:,1)
c$$$     &              + A0(:,5,2)*g3yi(:,2)
c$$$     &              + A0(:,5,3)*g3yi(:,3)
c$$$     &              + A0(:,5,4)*g3yi(:,4)
c$$$     &              + A0(:,5,5)*g3yi(:,5)
c$$$c
c$$$c.... calculate (giju A0 y_,j) --> gAgyi
c$$$c
c$$$
c$$$        gAgyi( :,1) = giju(:,1)*A0gyi( :,1)
c$$$     &              + giju(:,4)*A0gyi( :,6)
c$$$     &              + giju(:,5)*A0gyi(:,11)
c$$$
c$$$        gAgyi( :,2) = giju(:,1)*A0gyi( :,2)
c$$$     &              + giju(:,4)*A0gyi( :,7)
c$$$     &              + giju(:,5)*A0gyi(:,12)
c$$$
c$$$	gAgyi( :,3) = giju(:,1)*A0gyi( :,3)
c$$$     &              + giju(:,4)*A0gyi( :,8)
c$$$     &              + giju(:,5)*A0gyi(:,13)
c$$$
c$$$	gAgyi( :,4) = giju(:,1)*A0gyi( :,4)
c$$$     &              + giju(:,4)*A0gyi( :,9)
c$$$     &              + giju(:,5)*A0gyi(:,14)
c$$$
c$$$	gAgyi( :,5) = giju(:,1)*A0gyi( :,5)
c$$$     &              + giju(:,4)*A0gyi(:,10)
c$$$     &              + giju(:,5)*A0gyi(:,15)
c$$$
c$$$	gAgyi( :,6) = giju(:,4)*A0gyi( :,1)
c$$$     &              + giju(:,2)*A0gyi( :,6)
c$$$     &              + giju(:,6)*A0gyi(:,11)
c$$$
c$$$	gAgyi( :,7) = giju(:,4)*A0gyi( :,2)
c$$$     &              + giju(:,2)*A0gyi( :,7)
c$$$     &              + giju(:,6)*A0gyi(:,12)
c$$$
c$$$	gAgyi( :,8) = giju(:,4)*A0gyi( :,3)
c$$$     &              + giju(:,2)*A0gyi( :,8)
c$$$     &              + giju(:,6)*A0gyi(:,13)
c$$$
c$$$	gAgyi( :,9) = giju(:,4)*A0gyi( :,4)
c$$$     &              + giju(:,2)*A0gyi( :,9)
c$$$     &              + giju(:,6)*A0gyi(:,14)
c$$$
c$$$	gAgyi(:,10) = giju(:,4)*A0gyi( :,5)
c$$$     &              + giju(:,2)*A0gyi(:,10)
c$$$     &              + giju(:,6)*A0gyi(:,15)
c$$$
c$$$	gAgyi(:,11) = giju(:,5)*A0gyi( :,1)
c$$$     &              + giju(:,6)*A0gyi( :,6)
c$$$     &              + giju(:,3)*A0gyi(:,11)
c$$$
c$$$	gAgyi(:,12) = giju(:,5)*A0gyi( :,2)
c$$$     &              + giju(:,6)*A0gyi( :,7)
c$$$     &              + giju(:,3)*A0gyi(:,12)
c$$$
c$$$	gAgyi(:,13) = giju(:,5)*A0gyi( :,3)
c$$$     &              + giju(:,6)*A0gyi( :,8)
c$$$     &              + giju(:,3)*A0gyi(:,13)
c$$$
c$$$	gAgyi(:,14) = giju(:,5)*A0gyi( :,4)
c$$$     &              + giju(:,6)*A0gyi( :,9)
c$$$     &              + giju(:,3)*A0gyi(:,14)
c$$$
c$$$	gAgyi(:,15) = giju(:,5)*A0gyi( :,5)
c$$$     &              + giju(:,6)*A0gyi(:,10)
c$$$     &              + giju(:,3)*A0gyi(:,15)
c$$$c	
c$$$c... the denominator term of the DC factor
c$$$c... evaluation of the term  Y,i.A0DC Y,j 
c$$$c
c$$$        yiA0DCyj(:,1) = A0DC(:,1)*g1yi(:,1)**2
c$$$     &                + two*g1yi(:,1)*A0DC(:,2)*g1yi(:,5)
c$$$     &                + A0DC(:,3)*g1yi(:,2)**2
c$$$     &                + A0DC(:,3)*g1yi(:,3)**2
c$$$     &                + A0DC(:,3)*g1yi(:,4)**2
c$$$     &                + A0DC(:,4)*g1yi(:,5)**2
c$$$
c$$$        yiA0DCyj(:,2) = A0DC(:,1)*g2yi(:,1)**2
c$$$     &                + two*g2yi(:,1)*A0DC(:,2)*g2yi(:,5)
c$$$     &                + A0DC(:,3)*g2yi(:,2)**2
c$$$     &                + A0DC(:,3)*g2yi(:,3)**2
c$$$     &                + A0DC(:,3)*g2yi(:,4)**2
c$$$     &                + A0DC(:,4)*g2yi(:,5)**2
c$$$
c$$$        yiA0DCyj(:,3) = A0DC(:,1)*g3yi(:,1)**2
c$$$     &                + two*g3yi(:,1)*A0DC(:,2)*g3yi(:,5)
c$$$     &                + A0DC(:,3)*g3yi(:,2)**2
c$$$     &                + A0DC(:,3)*g3yi(:,3)**2
c$$$     &                + A0DC(:,3)*g3yi(:,4)**2
c$$$     &                + A0DC(:,4)*g3yi(:,5)**2
c$$$
c$$$        yiA0DCyj(:,4) = g1yi(:,1)*A0DC(:,1)*g2yi(:,1)
c$$$     &                + g1yi(:,1)*A0DC(:,2)*g2yi(:,5)
c$$$     &                + g1yi(:,2)*A0DC(:,3)*g2yi(:,2)
c$$$     &                + g1yi(:,3)*A0DC(:,3)*g2yi(:,3)
c$$$     &                + g1yi(:,4)*A0DC(:,3)*g2yi(:,4)
c$$$     &                + g1yi(:,5)*A0DC(:,2)*g2yi(:,1)
c$$$     &                + g1yi(:,5)*A0DC(:,4)*g2yi(:,5)
c$$$
c$$$        yiA0DCyj(:,5) = g1yi(:,1)*A0DC(:,1)*g3yi(:,1)
c$$$     &                + g1yi(:,1)*A0DC(:,2)*g3yi(:,5)
c$$$     &                + g1yi(:,2)*A0DC(:,3)*g3yi(:,2)
c$$$     &                + g1yi(:,3)*A0DC(:,3)*g3yi(:,3)
c$$$     &                + g1yi(:,4)*A0DC(:,3)*g3yi(:,4)
c$$$     &                + g1yi(:,5)*A0DC(:,2)*g3yi(:,1)
c$$$     &                + g1yi(:,5)*A0DC(:,4)*g3yi(:,5)
c$$$
c$$$        yiA0DCyj(:,6) = g2yi(:,1)*A0DC(:,1)*g3yi(:,1)
c$$$     &                + g2yi(:,1)*A0DC(:,2)*g3yi(:,5)
c$$$     &                + g2yi(:,2)*A0DC(:,3)*g3yi(:,2)
c$$$     &                + g2yi(:,3)*A0DC(:,3)*g3yi(:,3)
c$$$     &                + g2yi(:,4)*A0DC(:,3)*g3yi(:,4)
c$$$     &                + g2yi(:,5)*A0DC(:,2)*g3yi(:,1)
c$$$     &                + g2yi(:,5)*A0DC(:,4)*g3yi(:,5)
c$$$c
c$$$c.... ------------------------->  DC factor  <--------------------------
c$$$c
c$$$	if ((ires .ne. 2) .or. (Jactyp .eq. 1)) then
c$$$c
c$$$c.... calculate 2-norm of Grad-local-V with respect to A0
c$$$c
c$$$c.... DC-mallet
c$$$c
c$$$	  if (iDC .eq. 1) then
c$$$c
c$$$	    fact = one
c$$$	    if (ipord .eq. 2)  fact = 0.9
c$$$	    if (ipord .eq. 3) fact = 0.75
c$$$	
c$$$c
c$$$            gnorm = one / (
c$$$     &              giju(:,1)*yiA0DCyj(:,1)
c$$$     &            + two*giju(:,4)*yiA0DCyj(:,4)
c$$$     &            + two*giju(:,5)*yiA0DCyj(:,5)
c$$$     &            + giju(:,2)*yiA0DCyj(:,2) 
c$$$     &            + two*giju(:,6)*yiA0DCyj(:,6)
c$$$     &            + giju(:,3)*yiA0DCyj(:,3) 
c$$$     &            + epsM  )
c$$$c
c$$$	    DC(:,intp)=dim((fact*sqrt(raLS*gnorm)),(rtLS*gnorm))
c$$$c
c$$$c.... flop count
c$$$c
c$$$	    flops = flops + 46*npro
c$$$c
c$$$	  endif
c$$$c
c$$$c.... DC-quadratic
c$$$c
c$$$	  if (iDC .eq. 2) then
c$$$c
c$$$            gnorm = one / (
c$$$     &              giju(:,1)*yiA0DCyj(:,1)
c$$$     &            + two*giju(:,4)*yiA0DCyj(:,4)
c$$$     &            + two*giju(:,5)*yiA0DCyj(:,5)
c$$$     &            + giju(:,2)*yiA0DCyj(:,2) 
c$$$     &            + two*giju(:,6)*yiA0DCyj(:,6)
c$$$     &            + giju(:,3)*yiA0DCyj(:,3) 
c$$$     &            + epsM  )
c$$$         
c$$$c
c$$$	    DC(:,intp) = two * rtLS * gnorm
c$$$c
c$$$c.... flop count
c$$$c
c$$$	    flops = flops + 36*npro
c$$$c
c$$$	  endif
c$$$c
c$$$c.... DC-min
c$$$c
c$$$	  if (iDC .eq. 3) then
c$$$c
c$$$	    fact = one
c$$$	    if (ipord .eq. 2)  fact = pt5
c$$$c
c$$$            gnorm = one / (
c$$$     &              giju(:,1)*yiA0DCyj(:,1)
c$$$     &            + two*giju(:,4)*yiA0DCyj(:,4)
c$$$     &            + two*giju(:,5)*yiA0DCyj(:,5)
c$$$     &            + giju(:,2)*yiA0DCyj(:,2) 
c$$$     &            + two*giju(:,6)*yiA0DCyj(:,6)
c$$$     &            + giju(:,3)*yiA0DCyj(:,3) 
c$$$     &            + epsM  )
c$$$
c$$$c
c$$$	    DC(:,intp) = min( dim(fact * sqrt(raLS * gnorm),
c$$$     &                       rtLS * gnorm), two * rtLS * gnorm )
c$$$c
c$$$c.... flop count
c$$$c
c$$$	    flops = flops + 48*npro
c$$$c
c$$$	  endif
c$$$c
c$$$	endif
c$$$c
c$$$c.... ---------------------------->  RHS  <----------------------------
c$$$c
c$$$c.... add the contribution of DC to ri and/or rmi
c$$$c
c$$$c.... ires = 1 or 3
c$$$c
c$$$	if ((ires .eq. 1) .or. (ires .eq. 3)) then
c$$$c
c$$$	  ri ( :,1) = ri ( :,1) + DC(:,intp) * gAgyi( :,1)
c$$$	  rmi( :,1) = rmi( :,1) + DC(:,intp) * gAgyi( :,1)
c$$$	  ri ( :,2) = ri ( :,2) + DC(:,intp) * gAgyi( :,2)
c$$$	  rmi( :,2) = rmi( :,2) + DC(:,intp) * gAgyi( :,2)
c$$$	  ri ( :,3) = ri ( :,3) + DC(:,intp) * gAgyi( :,3)
c$$$	  rmi( :,3) = rmi( :,3) + DC(:,intp) * gAgyi( :,3)
c$$$	  ri ( :,4) = ri ( :,4) + DC(:,intp) * gAgyi( :,4)
c$$$	  rmi( :,4) = rmi( :,4) + DC(:,intp) * gAgyi( :,4)
c$$$	  ri ( :,5) = ri ( :,5) + DC(:,intp) * gAgyi( :,5)
c$$$	  rmi( :,5) = rmi( :,5) + DC(:,intp) * gAgyi( :,5)
c$$$c
c$$$	  ri ( :,6) = ri ( :,6) + DC(:,intp) * gAgyi( :,6)
c$$$	  rmi( :,6) = rmi( :,6) + DC(:,intp) * gAgyi( :,6)
c$$$	  ri ( :,7) = ri ( :,7) + DC(:,intp) * gAgyi( :,7)
c$$$	  rmi( :,7) = rmi( :,7) + DC(:,intp) * gAgyi( :,7)
c$$$	  ri ( :,8) = ri ( :,8) + DC(:,intp) * gAgyi( :,8)
c$$$	  rmi( :,8) = rmi( :,8) + DC(:,intp) * gAgyi( :,8)
c$$$	  ri ( :,9) = ri ( :,9) + DC(:,intp) * gAgyi( :,9)
c$$$	  rmi( :,9) = rmi( :,9) + DC(:,intp) * gAgyi( :,9)
c$$$	  ri (:,10) = ri (:,10) + DC(:,intp) * gAgyi(:,10)
c$$$	  rmi(:,10) = rmi(:,10) + DC(:,intp) * gAgyi(:,10)
c$$$c
c$$$	  ri (:,11) = ri (:,11) + DC(:,intp) * gAgyi(:,11)
c$$$	  rmi(:,11) = rmi(:,12) + DC(:,intp) * gAgyi(:,12)
c$$$	  ri (:,12) = ri (:,12) + DC(:,intp) * gAgyi(:,12)
c$$$	  rmi(:,12) = rmi(:,12) + DC(:,intp) * gAgyi(:,12)
c$$$	  ri (:,13) = ri (:,13) + DC(:,intp) * gAgyi(:,13)
c$$$	  rmi(:,13) = rmi(:,13) + DC(:,intp) * gAgyi(:,13)
c$$$	  ri (:,14) = ri (:,14) + DC(:,intp) * gAgyi(:,14)
c$$$	  rmi(:,14) = rmi(:,14) + DC(:,intp) * gAgyi(:,14)
c$$$	  ri (:,15) = ri (:,15) + DC(:,intp) * gAgyi(:,15)
c$$$	  rmi(:,15) = rmi(:,15) + DC(:,intp) * gAgyi(:,15)
c$$$c
c$$$	  flops = flops + 45*npro
c$$$c
c$$$	endif
c$$$c
c$$$c.... ires = 2
c$$$c
c$$$	if (ires .eq. 2) then
c$$$c
c$$$	  rmi( :,1) = rmi( :,1) + DC(:,intp) * gAgyi( :,1)
c$$$	  rmi( :,2) = rmi( :,2) + DC(:,intp) * gAgyi( :,2)
c$$$	  rmi( :,3) = rmi( :,3) + DC(:,intp) * gAgyi( :,3)
c$$$	  rmi( :,4) = rmi( :,4) + DC(:,intp) * gAgyi( :,4)
c$$$	  rmi( :,5) = rmi( :,5) + DC(:,intp) * gAgyi( :,5)
c$$$c
c$$$	  rmi( :,6) = rmi( :,6) + DC(:,intp) * gAgyi( :,6)
c$$$	  rmi( :,7) = rmi( :,7) + DC(:,intp) * gAgyi( :,7)
c$$$	  rmi( :,8) = rmi( :,8) + DC(:,intp) * gAgyi( :,8)
c$$$	  rmi( :,9) = rmi( :,9) + DC(:,intp) * gAgyi( :,9)
c$$$	  rmi(:,10) = rmi(:,10) + DC(:,intp) * gAgyi(:,10)
c$$$c
c$$$	  rmi(:,11) = rmi(:,11) + DC(:,intp) * gAgyi(:,11)
c$$$	  rmi(:,12) = rmi(:,12) + DC(:,intp) * gAgyi(:,12)
c$$$	  rmi(:,13) = rmi(:,13) + DC(:,intp) * gAgyi(:,13)
c$$$	  rmi(:,14) = rmi(:,14) + DC(:,intp) * gAgyi(:,14)
c$$$	  rmi(:,15) = rmi(:,15) + DC(:,intp) * gAgyi(:,15)
c$$$c
c$$$	  flops = flops + 30*npro
c$$$c
c$$$	endif
c$$$c
c$$$c.... ------------------------->  Stiffness  <--------------------------
c$$$c
c$$$c.... add the contribution of DC to stiff
c$$$c
c$$$	if (iprec .eq. 1) then
c$$$	     nflow2=two*nflow
c$$$       do j = 1, nflow
c$$$          do i = 1, nflow
c$$$             itmp(:)=A0(:,i,j)*DC(:,intp)
c$$$c
c$$$c.... add (DC g^1 A0) to stiff [1,1]
c$$$c
c$$$             stiff(:,i,j) = stiff(:,i,j) 
c$$$     &                    + itmp(:)*giju(:,1)
c$$$c
c$$$c.... add (DC g^1 A0) to stiff [1,2]
c$$$c
c$$$
c$$$             stiff(:,i,j+nflow) = stiff(:,i,j+nflow) 
c$$$     &                    + itmp(:)*giju(:,4)
c$$$c
c$$$c.... add (DC g^1 A0) to stiff [1,3]
c$$$c
c$$$
c$$$             stiff(:,i,j+nflow2) = stiff(:,i,j+nflow2) 
c$$$     &                    + itmp(:)*giju(:,5)
c$$$
c$$$c.... add (DC g^1 A0) to stiff [2,1] (similarly below)
c$$$c
c$$$
c$$$             stiff(:,i+nflow,j) = stiff(:,i+nflow,j) 
c$$$     &                    + itmp(:)*giju(:,4)
c$$$
c$$$             stiff(:,i+nflow,j+nflow) = stiff(:,i+nflow,j+nflow) 
c$$$     &                    + itmp(:)*giju(:,2)
c$$$
c$$$             stiff(:,i+nflow,j+nflow2) = stiff(:,i+nflow,j+nflow2) 
c$$$     &                    + itmp(:)*giju(:,6)
c$$$
c$$$             stiff(:,i+nflow2,j) = stiff(:,i+nflow2,j) 
c$$$     &                    + itmp(:)*giju(:,5)
c$$$
c$$$             stiff(:,i+nflow2,j+nflow) = stiff(:,i+nflow2,j+nflow) 
c$$$     &                    + itmp(:)*giju(:,6)
c$$$
c$$$             stiff(:,i+nflow2,j+nflow2) = stiff(:,i+nflow2,j+nflow2) 
c$$$     &                    + itmp(:)*giju(:,3)
c$$$          enddo
c$$$       enddo
c$$$c
c$$$c.... flop count
c$$$c
c$$$	  flops = flops + 210*npro
c$$$c
c$$$c.... end of stiffness
c$$$c
c$$$	endif
c$$$c
c$$$c.... return
c$$$c
c$$$	return
c$$$	end
c$$$c
c
        subroutine e3dcSclr ( gradS,    giju,     gGradS,
     &                        rLS,      tauS,     srcR,
     &                        dcFct)
c
c
c----------------------------------------------------------------------
c
c This routine calculates the contribution of the Discontinuity-
c Capturing operator to RHS and preconditioner for the scalar solve.
c
c  g1yti   (nflow,npro)           : grad-y in direction 1
c  g2yti   (nflow,npro)           : grad-y in direction 2
c  g3yti   (nflow,npro)           : grad-y in direction 3
c  A0     (nsymdf,npro)          : A0 matrix (Symm. storage)
c  raLS   (npro)                 : square of LS residual (A0inv norm)
c  rtLS   (npro)                 : square of LS residual (Tau norm)
c  giju    (6,npro)              : metric matrix
c  DC     (ngauss,npro)          : discontinuity-capturing factor
c  intp				 : integration point number
c
c output:
c  ri     (nflow*(nsd+1),npro)   : partial residual
c  rmi    (nflow*(nsd+1),npro)   : partial modified residual
c  stiff  (nsymdf,6,npro)       : diffusivity matrix
c  DC     (npro)                : discontinuity-capturing factor
c
c
c Zdenek Johan, Summer 1990. (Modified from e2dc.f)
c Zdenek Johan, Winter 1991. (Recoded)
c Zdenek Johan, Winter 1991. (Fortran 90)
c----------------------------------------------------------------------
c
	include "common.h"
c
        dimension gradS(npro,nsd),            gGradS(npro,nsd),
     &            rLS(npro),                  tauS(npro),
     &            giju(npro,6),               dcFct(npro),
     &            srcR(npro)
c
c.... Form GijUp gradS and  gradS . GijUp gradS (store in dcFct)
c
	
	    gGradS(:,1) = GijU(:,1) * gradS(:,1)
     1			+ GijU(:,4) * gradS(:,2)
     2			+ GijU(:,6) * gradS(:,3)
	    gGradS(:,2) = GijU(:,4) * gradS(:,1)
     1			+ GijU(:,2) * gradS(:,2)
     2			+ GijU(:,5) * gradS(:,3)
	    gGradS(:,3) = GijU(:,6) * gradS(:,1)
     1			+ GijU(:,5) * gradS(:,2)
     2			+ GijU(:,3) * gradS(:,3)
c
	    dcFct(:)    = gradS(:,1) * gGradS(:,1)
     1		        + gradS(:,2) * gGradS(:,2)
     2		        + gradS(:,3) * gGradS(:,3)
     3		        + epsM
	
	    dcFct(:) = 1.0/ dcFct(:)
c
c.... Form pdeRes 2-norm / gradT 2-norm
c

	    dcFct  = dcFct * (rLS - srcR) ** 2 
c
c.... ------------------------->  DC factor  <------------------------
c
c.... DC-mallet
c
	    if (idcsclr(1) .eq. 1) then
c       
	       fact = one
	       if (ipord .eq. 2)  fact = 0.9
	       if (ipord .eq. 3) fact = 0.75
c       
c$$$  dcFct(:)=dim((fact*sqrt(dcFct(:))),(tauS(:)*dcFct(:))) !not work
                                                          !with all compilers
	       dcFct(:)=max(zero,(fact*sqrt(dcFct(:)))-(tauS(:)*dcFct(:)))
c
	    endif
c       
c       
c....   DC-quadratic
c       
	    if (idcsclr(1) .eq. 2) then
c       
	       dcFct(:) = two * tauS(:) * dcFct(:)
c       
	    endif
c       
c....   DC-min
c       
	    if (idcsclr(1) .eq. 3) then
c       
	       fact = one
	       if (ipord .eq. 2)  fact = 0.9
c       
          dcFct(:) = min( max(zero, (fact * sqrt(dcFct(:)) -
     &	             tauS(:)*dcFct(:)) ), two * tauS(:) * dcFct(:))
c       
	    endif
c
c.... Scale the gGradT for residual formation
c	
	    gGradS(:,1) = dcFct(:) * gGradS(:,1)
	    gGradS(:,2) = dcFct(:) * gGradS(:,2)
	    gGradS(:,3) = dcFct(:) * gGradS(:,3)
	


	return
	end
c
