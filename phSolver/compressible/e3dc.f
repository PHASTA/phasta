	subroutine e3DC (g1yi,   g2yi,   g3yi,   A0,     raLS,
     &			 rtLS,   giju,   DC,     ri,
     &                   rmi,    stiff, A0DC)
c
c----------------------------------------------------------------------
c
c This routine calculates the contribution of the Discontinuity-
c Capturing operator to RHS and preconditioner.
c
c  g1yi   (nflow,npro)           : grad-y in direction 1
c  g2yi   (nflow,npro)           : grad-y in direction 2
c  g3yi   (nflow,npro)           : grad-y in direction 3
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
        dimension g1yi(npro,nflow),          g2yi(npro,nflow),
     &            g3yi(npro,nflow),          A0(npro,5,5),
     &            raLS(npro),                rtLS(npro),
     &            giju(npro,6),              DC(npro,ngauss),
     &            ri(npro,nflow*(nsd+1)),    rmi(npro,nflow*(nsd+1)),
     &            stiff(npro,3*nflow,3*nflow),dtmp(npro)
c

        dimension ggyi(npro,nflow),         gAgyi(npro,15),
     &            gnorm(npro),              A0gyi(npro,15),
     &            yiA0DCyj(npro,6),         A0DC(npro,4)
c
c ... -----------------------> initialize <----------------------------
c
        A0gyi    = zero
        gAgyi    = zero
        yiA0DCyj = zero
        DC       = zero
c.... ----------------------->  global gradient  <----------------------
c
c.... calculate (A0 y_,j) --> A0gyi
c
c  A0 Y_{,1}
c
        A0gyi( :,1) = A0(:,1,1)*g1yi(:,1)
     &              + A0(:,1,2)*g1yi(:,2)
     &              + A0(:,1,3)*g1yi(:,3)
     &              + A0(:,1,4)*g1yi(:,4)
     &              + A0(:,1,5)*g1yi(:,5)      
        A0gyi( :,2) = A0(:,2,1)*g1yi(:,1)
     &              + A0(:,2,2)*g1yi(:,2)
     &              + A0(:,2,3)*g1yi(:,3)
     &              + A0(:,2,4)*g1yi(:,4)
     &              + A0(:,2,5)*g1yi(:,5)
        A0gyi( :,3) = A0(:,3,1)*g1yi(:,1)
     &              + A0(:,3,2)*g1yi(:,2)
     &              + A0(:,3,3)*g1yi(:,3)
     &              + A0(:,3,4)*g1yi(:,4)
     &              + A0(:,3,5)*g1yi(:,5)
        A0gyi( :,4) = A0(:,4,1)*g1yi(:,1)
     &              + A0(:,4,2)*g1yi(:,2)
     &              + A0(:,4,3)*g1yi(:,3)
     &              + A0(:,4,4)*g1yi(:,4)
     &              + A0(:,4,5)*g1yi(:,5)
        A0gyi( :,5) = A0(:,5,1)*g1yi(:,1)
     &              + A0(:,5,2)*g1yi(:,2)
     &              + A0(:,5,3)*g1yi(:,3)
     &              + A0(:,5,4)*g1yi(:,4)
     &              + A0(:,5,5)*g1yi(:,5)
c
c  A0 Y_{,2}
c
        A0gyi( :,6) = A0(:,1,1)*g2yi(:,1)
     &              + A0(:,1,2)*g2yi(:,2)
     &              + A0(:,1,3)*g2yi(:,3)
     &              + A0(:,1,4)*g2yi(:,4)
     &              + A0(:,1,5)*g2yi(:,5)
        A0gyi( :,7) = A0(:,2,1)*g2yi(:,1)
     &              + A0(:,2,2)*g2yi(:,2)
     &              + A0(:,2,3)*g2yi(:,3)
     &              + A0(:,2,4)*g2yi(:,4)
     &              + A0(:,2,5)*g2yi(:,5)
        A0gyi( :,8) = A0(:,3,1)*g2yi(:,1)
     &              + A0(:,3,2)*g2yi(:,2)
     &              + A0(:,3,3)*g2yi(:,3)
     &              + A0(:,3,4)*g2yi(:,4)
     &              + A0(:,3,5)*g2yi(:,5)
        A0gyi( :,9) = A0(:,4,1)*g2yi(:,1)
     &              + A0(:,4,2)*g2yi(:,2)
     &              + A0(:,4,3)*g2yi(:,3)
     &              + A0(:,4,4)*g2yi(:,4)
     &              + A0(:,4,5)*g2yi(:,5)
        A0gyi(:,10) = A0(:,5,1)*g2yi(:,1)
     &              + A0(:,5,2)*g2yi(:,2)
     &              + A0(:,5,3)*g2yi(:,3)
     &              + A0(:,5,4)*g2yi(:,4)
     &              + A0(:,5,5)*g2yi(:,5)
c
c  A0 Y_{,3}
c
        A0gyi(:,11) = A0(:,1,1)*g3yi(:,1)
     &              + A0(:,1,2)*g3yi(:,2)
     &              + A0(:,1,3)*g3yi(:,3)
     &              + A0(:,1,4)*g3yi(:,4)
     &              + A0(:,1,5)*g3yi(:,5)
        A0gyi(:,12) = A0(:,2,1)*g3yi(:,1)
     &              + A0(:,2,2)*g3yi(:,2)
     &              + A0(:,2,3)*g3yi(:,3)
     &              + A0(:,2,4)*g3yi(:,4)
     &              + A0(:,2,5)*g3yi(:,5)
        A0gyi(:,13) = A0(:,3,1)*g3yi(:,1)
     &              + A0(:,3,2)*g3yi(:,2)
     &              + A0(:,3,3)*g3yi(:,3)
     &              + A0(:,3,4)*g3yi(:,4)
     &              + A0(:,3,5)*g3yi(:,5)
        A0gyi(:,14) = A0(:,4,1)*g3yi(:,1)
     &              + A0(:,4,2)*g3yi(:,2)
     &              + A0(:,4,3)*g3yi(:,3)
     &              + A0(:,4,4)*g3yi(:,4)
     &              + A0(:,4,5)*g3yi(:,5)
        A0gyi(:,15) = A0(:,5,1)*g3yi(:,1)
     &              + A0(:,5,2)*g3yi(:,2)
     &              + A0(:,5,3)*g3yi(:,3)
     &              + A0(:,5,4)*g3yi(:,4)
     &              + A0(:,5,5)*g3yi(:,5)
c
c.... calculate (giju A0 y_,j) --> gAgyi
c

        gAgyi( :,1) = giju(:,1)*A0gyi( :,1)
     &              + giju(:,4)*A0gyi( :,6)
     &              + giju(:,5)*A0gyi(:,11)

        gAgyi( :,2) = giju(:,1)*A0gyi( :,2)
     &              + giju(:,4)*A0gyi( :,7)
     &              + giju(:,5)*A0gyi(:,12)

	gAgyi( :,3) = giju(:,1)*A0gyi( :,3)
     &              + giju(:,4)*A0gyi( :,8)
     &              + giju(:,5)*A0gyi(:,13)

	gAgyi( :,4) = giju(:,1)*A0gyi( :,4)
     &              + giju(:,4)*A0gyi( :,9)
     &              + giju(:,5)*A0gyi(:,14)

	gAgyi( :,5) = giju(:,1)*A0gyi( :,5)
     &              + giju(:,4)*A0gyi(:,10)
     &              + giju(:,5)*A0gyi(:,15)

	gAgyi( :,6) = giju(:,4)*A0gyi( :,1)
     &              + giju(:,2)*A0gyi( :,6)
     &              + giju(:,6)*A0gyi(:,11)

	gAgyi( :,7) = giju(:,4)*A0gyi( :,2)
     &              + giju(:,2)*A0gyi( :,7)
     &              + giju(:,6)*A0gyi(:,12)

	gAgyi( :,8) = giju(:,4)*A0gyi( :,3)
     &              + giju(:,2)*A0gyi( :,8)
     &              + giju(:,6)*A0gyi(:,13)

	gAgyi( :,9) = giju(:,4)*A0gyi( :,4)
     &              + giju(:,2)*A0gyi( :,9)
     &              + giju(:,6)*A0gyi(:,14)

	gAgyi(:,10) = giju(:,4)*A0gyi( :,5)
     &              + giju(:,2)*A0gyi(:,10)
     &              + giju(:,6)*A0gyi(:,15)

	gAgyi(:,11) = giju(:,5)*A0gyi( :,1)
     &              + giju(:,6)*A0gyi( :,6)
     &              + giju(:,3)*A0gyi(:,11)

	gAgyi(:,12) = giju(:,5)*A0gyi( :,2)
     &              + giju(:,6)*A0gyi( :,7)
     &              + giju(:,3)*A0gyi(:,12)

	gAgyi(:,13) = giju(:,5)*A0gyi( :,3)
     &              + giju(:,6)*A0gyi( :,8)
     &              + giju(:,3)*A0gyi(:,13)

	gAgyi(:,14) = giju(:,5)*A0gyi( :,4)
     &              + giju(:,6)*A0gyi( :,9)
     &              + giju(:,3)*A0gyi(:,14)

	gAgyi(:,15) = giju(:,5)*A0gyi( :,5)
     &              + giju(:,6)*A0gyi(:,10)
     &              + giju(:,3)*A0gyi(:,15)
c	
c... the denominator term of the DC factor
c... evaluation of the term  Y,i.A0DC Y,j 
c
        yiA0DCyj(:,1) = A0DC(:,1)*g1yi(:,1)**2
     &                + two*g1yi(:,1)*A0DC(:,2)*g1yi(:,5)
     &                + A0DC(:,3)*g1yi(:,2)**2
     &                + A0DC(:,3)*g1yi(:,3)**2
     &                + A0DC(:,3)*g1yi(:,4)**2
     &                + A0DC(:,4)*g1yi(:,5)**2

        yiA0DCyj(:,2) = A0DC(:,1)*g2yi(:,1)**2
     &                + two*g2yi(:,1)*A0DC(:,2)*g2yi(:,5)
     &                + A0DC(:,3)*g2yi(:,2)**2
     &                + A0DC(:,3)*g2yi(:,3)**2
     &                + A0DC(:,3)*g2yi(:,4)**2
     &                + A0DC(:,4)*g2yi(:,5)**2

        yiA0DCyj(:,3) = A0DC(:,1)*g3yi(:,1)**2
     &                + two*g3yi(:,1)*A0DC(:,2)*g3yi(:,5)
     &                + A0DC(:,3)*g3yi(:,2)**2
     &                + A0DC(:,3)*g3yi(:,3)**2
     &                + A0DC(:,3)*g3yi(:,4)**2
     &                + A0DC(:,4)*g3yi(:,5)**2

        yiA0DCyj(:,4) = g1yi(:,1)*A0DC(:,1)*g2yi(:,1)
     &                + g1yi(:,1)*A0DC(:,2)*g2yi(:,5)
     &                + g1yi(:,2)*A0DC(:,3)*g2yi(:,2)
     &                + g1yi(:,3)*A0DC(:,3)*g2yi(:,3)
     &                + g1yi(:,4)*A0DC(:,3)*g2yi(:,4)
     &                + g1yi(:,5)*A0DC(:,2)*g2yi(:,1)
     &                + g1yi(:,5)*A0DC(:,4)*g2yi(:,5)

        yiA0DCyj(:,5) = g1yi(:,1)*A0DC(:,1)*g3yi(:,1)
     &                + g1yi(:,1)*A0DC(:,2)*g3yi(:,5)
     &                + g1yi(:,2)*A0DC(:,3)*g3yi(:,2)
     &                + g1yi(:,3)*A0DC(:,3)*g3yi(:,3)
     &                + g1yi(:,4)*A0DC(:,3)*g3yi(:,4)
     &                + g1yi(:,5)*A0DC(:,2)*g3yi(:,1)
     &                + g1yi(:,5)*A0DC(:,4)*g3yi(:,5)

        yiA0DCyj(:,6) = g2yi(:,1)*A0DC(:,1)*g3yi(:,1)
     &                + g2yi(:,1)*A0DC(:,2)*g3yi(:,5)
     &                + g2yi(:,2)*A0DC(:,3)*g3yi(:,2)
     &                + g2yi(:,3)*A0DC(:,3)*g3yi(:,3)
     &                + g2yi(:,4)*A0DC(:,3)*g3yi(:,4)
     &                + g2yi(:,5)*A0DC(:,2)*g3yi(:,1)
     &                + g2yi(:,5)*A0DC(:,4)*g3yi(:,5)
c
c.... ------------------------->  DC factor  <--------------------------
c
c	if ((ires .ne. 2) .or. (Jactyp .eq. 1)) then
c
c.... calculate 2-norm of Grad-local-V with respect to A0
c
c.... DC-mallet
c
	  if (iDC .eq. 1) then
c
	    fact = one
	    if (ipord .eq. 2)  fact = 0.9
	    if (ipord .eq. 3) fact = 0.75
	
c
            gnorm = one / (
     &              giju(:,1)*yiA0DCyj(:,1)
     &            + two*giju(:,4)*yiA0DCyj(:,4)
     &            + two*giju(:,5)*yiA0DCyj(:,5)
     &            + giju(:,2)*yiA0DCyj(:,2) 
     &            + two*giju(:,6)*yiA0DCyj(:,6)
     &            + giju(:,3)*yiA0DCyj(:,3) 
     &            + epsM  )
c
c	    DC(:,intp)=dim((fact*sqrt(raLS*gnorm)),(rtLS*gnorm))
	    DC(:,intp)=max(zero,(fact*sqrt(raLS*gnorm))-(rtLS*gnorm))
c
c.... flop count
c
!	    flops = flops + 46*npro
c
	  endif
c
c.... DC-quadratic
c
	  if (iDC .eq. 2) then
c
            gnorm = one / (
     &              giju(:,1)*yiA0DCyj(:,1)
     &            + two*giju(:,4)*yiA0DCyj(:,4)
     &            + two*giju(:,5)*yiA0DCyj(:,5)
     &            + giju(:,2)*yiA0DCyj(:,2) 
     &            + two*giju(:,6)*yiA0DCyj(:,6)
     &            + giju(:,3)*yiA0DCyj(:,3) 
     &            + epsM  )
         
c
	    DC(:,intp) = two * rtLS * gnorm
c
c.... flop count
c
!	    flops = flops + 36*npro
c
	  endif
c
c.... DC-min
c
	  if (iDC .eq. 3) then
c
	    fact = one
	    if (ipord .eq. 2)  fact = pt5
c
            gnorm = one / (
     &              giju(:,1)*yiA0DCyj(:,1)
     &            + two*giju(:,4)*yiA0DCyj(:,4)
     &            + two*giju(:,5)*yiA0DCyj(:,5)
     &            + giju(:,2)*yiA0DCyj(:,2) 
     &            + two*giju(:,6)*yiA0DCyj(:,6)
     &            + giju(:,3)*yiA0DCyj(:,3) 
     &            + epsM  )

c
c	    DC(:,intp) = min( dim(fact * sqrt(raLS * gnorm),
	    DC(:,intp) = min( max(zero,fact * sqrt(raLS * gnorm)-
     &                       rtLS * gnorm), two * rtLS * gnorm )
c
c.... flop count
c
!	    flops = flops + 48*npro
c
	  endif
c
c	endif
c
c.... ---------------------------->  RHS  <----------------------------
c
c.... add the contribution of DC to ri and/or rmi
c
c.... ires = 1 or 3
c
	if ((ires .eq. 1) .or. (ires .eq. 3)) then
c
	  ri ( :,1) = ri ( :,1) + DC(:,intp) * gAgyi( :,1)
	  rmi( :,1) = rmi( :,1) + DC(:,intp) * gAgyi( :,1)
	  ri ( :,2) = ri ( :,2) + DC(:,intp) * gAgyi( :,2)
	  rmi( :,2) = rmi( :,2) + DC(:,intp) * gAgyi( :,2)
	  ri ( :,3) = ri ( :,3) + DC(:,intp) * gAgyi( :,3)
	  rmi( :,3) = rmi( :,3) + DC(:,intp) * gAgyi( :,3)
	  ri ( :,4) = ri ( :,4) + DC(:,intp) * gAgyi( :,4)
	  rmi( :,4) = rmi( :,4) + DC(:,intp) * gAgyi( :,4)
	  ri ( :,5) = ri ( :,5) + DC(:,intp) * gAgyi( :,5)
	  rmi( :,5) = rmi( :,5) + DC(:,intp) * gAgyi( :,5)
c
	  ri ( :,6) = ri ( :,6) + DC(:,intp) * gAgyi( :,6)
	  rmi( :,6) = rmi( :,6) + DC(:,intp) * gAgyi( :,6)
	  ri ( :,7) = ri ( :,7) + DC(:,intp) * gAgyi( :,7)
	  rmi( :,7) = rmi( :,7) + DC(:,intp) * gAgyi( :,7)
	  ri ( :,8) = ri ( :,8) + DC(:,intp) * gAgyi( :,8)
	  rmi( :,8) = rmi( :,8) + DC(:,intp) * gAgyi( :,8)
	  ri ( :,9) = ri ( :,9) + DC(:,intp) * gAgyi( :,9)
	  rmi( :,9) = rmi( :,9) + DC(:,intp) * gAgyi( :,9)
	  ri (:,10) = ri (:,10) + DC(:,intp) * gAgyi(:,10)
	  rmi(:,10) = rmi(:,10) + DC(:,intp) * gAgyi(:,10)
c
	  ri (:,11) = ri (:,11) + DC(:,intp) * gAgyi(:,11)
	  rmi(:,11) = rmi(:,12) + DC(:,intp) * gAgyi(:,12)
	  ri (:,12) = ri (:,12) + DC(:,intp) * gAgyi(:,12)
	  rmi(:,12) = rmi(:,12) + DC(:,intp) * gAgyi(:,12)
	  ri (:,13) = ri (:,13) + DC(:,intp) * gAgyi(:,13)
	  rmi(:,13) = rmi(:,13) + DC(:,intp) * gAgyi(:,13)
	  ri (:,14) = ri (:,14) + DC(:,intp) * gAgyi(:,14)
	  rmi(:,14) = rmi(:,14) + DC(:,intp) * gAgyi(:,14)
	  ri (:,15) = ri (:,15) + DC(:,intp) * gAgyi(:,15)
	  rmi(:,15) = rmi(:,15) + DC(:,intp) * gAgyi(:,15)
c
!	  flops = flops + 45*npro
c
	endif
c
c.... ires = 2
c
	if (ires .eq. 2) then
c
	  rmi( :,1) = rmi( :,1) + DC(:,intp) * gAgyi( :,1)
	  rmi( :,2) = rmi( :,2) + DC(:,intp) * gAgyi( :,2)
	  rmi( :,3) = rmi( :,3) + DC(:,intp) * gAgyi( :,3)
	  rmi( :,4) = rmi( :,4) + DC(:,intp) * gAgyi( :,4)
	  rmi( :,5) = rmi( :,5) + DC(:,intp) * gAgyi( :,5)
c
	  rmi( :,6) = rmi( :,6) + DC(:,intp) * gAgyi( :,6)
	  rmi( :,7) = rmi( :,7) + DC(:,intp) * gAgyi( :,7)
	  rmi( :,8) = rmi( :,8) + DC(:,intp) * gAgyi( :,8)
	  rmi( :,9) = rmi( :,9) + DC(:,intp) * gAgyi( :,9)
	  rmi(:,10) = rmi(:,10) + DC(:,intp) * gAgyi(:,10)
c
	  rmi(:,11) = rmi(:,11) + DC(:,intp) * gAgyi(:,11)
	  rmi(:,12) = rmi(:,12) + DC(:,intp) * gAgyi(:,12)
	  rmi(:,13) = rmi(:,13) + DC(:,intp) * gAgyi(:,13)
	  rmi(:,14) = rmi(:,14) + DC(:,intp) * gAgyi(:,14)
	  rmi(:,15) = rmi(:,15) + DC(:,intp) * gAgyi(:,15)
c
!	  flops = flops + 30*npro
c
	endif
c
c.... ------------------------->  Stiffness  <--------------------------
c
c.... add the contribution of DC to stiff
c
	if (iprec .eq. 1) then ! leave out of LHS, when called from itrres
	     nflow2=two*nflow
       do j = 1, nflow
          do i = 1, nflow
             dtmp(:)=A0(:,i,j)*DC(:,intp)
c
c.... add (DC g^1 A0) to stiff [1,1]
c
             stiff(:,i,j) = stiff(:,i,j) 
     &                    + dtmp(:)*giju(:,1)
c
c.... add (DC g^1 A0) to stiff [1,2]
c

             stiff(:,i,j+nflow) = stiff(:,i,j+nflow) 
     &                    + dtmp(:)*giju(:,4)
c
c.... add (DC g^1 A0) to stiff [1,3]
c

             stiff(:,i,j+nflow2) = stiff(:,i,j+nflow2) 
     &                    + dtmp(:)*giju(:,5)

c.... add (DC g^1 A0) to stiff [2,1] (similarly below)
c

             stiff(:,i+nflow,j) = stiff(:,i+nflow,j) 
     &                    + dtmp(:)*giju(:,4)

             stiff(:,i+nflow,j+nflow) = stiff(:,i+nflow,j+nflow) 
     &                    + dtmp(:)*giju(:,2)

             stiff(:,i+nflow,j+nflow2) = stiff(:,i+nflow,j+nflow2) 
     &                    + dtmp(:)*giju(:,6)

             stiff(:,i+nflow2,j) = stiff(:,i+nflow2,j) 
     &                    + dtmp(:)*giju(:,5)

             stiff(:,i+nflow2,j+nflow) = stiff(:,i+nflow2,j+nflow) 
     &                    + dtmp(:)*giju(:,6)

             stiff(:,i+nflow2,j+nflow2) = stiff(:,i+nflow2,j+nflow2) 
     &                    + dtmp(:)*giju(:,3)
          enddo
       enddo
c
c.... flop count
c
!	  flops = flops + 210*npro
c
c.... end of stiffness
c
	endif
c
c.... return
c
	return
	end
c
        subroutine e3dcSclr ( g1yti,    g2yti,       g3yti,
     &                        A0t,      raLSt,       rTLSt,
     &                        DCt,      giju,       
     &                        rti,      rmti,        stifft)
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
        dimension g1yti(npro),                g2yti(npro),
     &            g3yti(npro),                A0t(npro),
     &            raLSt(npro),                rtLSt(npro),
     &            giju(npro,6),               DCt(npro,ngauss),
     &            rti(npro,nsd+1),            rmti(npro,nsd+1),
     &            stifft(npro,nsd,nsd),       dtmp(npro)
c

        dimension ggyit(npro,nflow),        gAgyit(npro,3),
     &            gnormt(npro),             A0gyit(npro,3),
     &            yiA0DCyjt(npro,6),        A0DCt(npro)
c
c ... -----------------------> initialize <----------------------------
c
        A0gyit    = zero
        gAgyit    = zero
        yiA0DCyjt = zero
        DCt       = zero
        A0DCt    = A0t
c.... ----------------------->  global gradient  <----------------------
c
c.... calculate (A0 y_,j) --> A0gyit
c
c  A0 Y_{,1}
c
        A0gyit( :,1) = A0t(:)*g1yti(:)
c  A0 Y_{,2}      
        A0gyit( :,2) = A0t(:)*g2yti(:)
c  A0 Y_{,3} 
        A0gyit( :,3) = A0t(:)*g3yti(:)
c
c.... calculate (giju A0 y_,j) --> gAgyit
c

        gAgyit( :,1) = giju(:,1)*A0gyit( :,1)
     &               + giju(:,4)*A0gyit( :,2)
     &               + giju(:,5)*A0gyit( :,3)

        gAgyit( :,2) = giju(:,4)*A0gyit( :,1)
     &               + giju(:,2)*A0gyit( :,2)
     &               + giju(:,6)*A0gyit( :,3)

	gAgyit( :,3) = giju(:,5)*A0gyit( :,1)
     &               + giju(:,6)*A0gyit( :,2)
     &               + giju(:,3)*A0gyit( :,3)
c	
c... the denominator term of the DC factor
c... evaluation of the term  Y,i.A0DC Y,j 
c
        yiA0DCyjt(:,1) = A0DCt(:)*g1yti(:)**2
c    
        yiA0DCyjt(:,2) = A0DCt(:)*g2yti(:)**2
c     
        yiA0DCyjt(:,3) = A0DCt(:)*g3yti(:)**2
c
        yiA0DCyjt(:,4) = A0DCt(:)*g1yti(:)*g2yti(:)
c
        yiA0DCyjt(:,5) = A0DCt(:)*g1yti(:)*g3yti(:)
c
        yiA0DCyjt(:,6) = A0DCt(:)*g2yti(:)*g3yti(:)
c
c
c.... ------------------------->  DC factor  <--------------------------
c
c	if ((ires .ne. 2) .or. (Jactyp .eq. 1)) then
c
c.... calculate 2-norm of Grad-local-V with respect to A0
c
c.... DC-mallet
c
	  if (iDCsclr(1) .eq. 1) then
c
	    fact = one
	    if (ipord .eq. 2)  fact = 0.9
	    if (ipord .eq. 3) fact = 0.75
	
c
            gnormt = one / (
     &              giju(:,1)*yiA0DCyjt(:,1)
     &            + two*giju(:,4)*yiA0DCyjt(:,4)
     &            + two*giju(:,5)*yiA0DCyjt(:,5)
     &            + giju(:,2)*yiA0DCyjt(:,2) 
     &            + two*giju(:,6)*yiA0DCyjt(:,6)
     &            + giju(:,3)*yiA0DCyjt(:,3) 
     &            + epsM  )
c
c	    DCt(:,intp)=dim((fact*sqrt(raLSt*gnormt)),(rtLSt*gnormt))
	    DCt(:,intp)=max(zero,(fact*sqrt(raLSt*gnormt))-(rtLSt*gnormt))
c
c.... flop count
c
!	    flops = flops + 46*npro
c
	  endif
c
c.... DC-quadratic
c
	  if (iDCSclr(1) .eq. 2) then  
c
            gnormt = one / (
     &              giju(:,1)*yiA0DCyjt(:,1)
     &            + two*giju(:,4)*yiA0DCyjt(:,4)
     &            + two*giju(:,5)*yiA0DCyjt(:,5)
     &            + giju(:,2)*yiA0DCyjt(:,2) 
     &            + two*giju(:,6)*yiA0DCyjt(:,6)
     &            + giju(:,3)*yiA0DCyjt(:,3) 
     &            + epsM  )
         
c
	    DCt(:,intp) = two * rtLSt * gnormt
c
c.... flop count
c
!	    flops = flops + 36*npro
c
	  endif
c
c.... DC-min
c
	  if (iDCSclr(1) .eq. 3) then  
c
	    fact = one
	    if (ipord .eq. 2)  fact = pt5
c
            gnormt = one / (
     &              giju(:,1)*yiA0DCyjt(:,1)
     &            + two*giju(:,4)*yiA0DCyjt(:,4)
     &            + two*giju(:,5)*yiA0DCyjt(:,5)
     &            + giju(:,2)*yiA0DCyjt(:,2) 
     &            + two*giju(:,6)*yiA0DCyjt(:,6)
     &            + giju(:,3)*yiA0DCyjt(:,3) 
     &            + epsM  )

c
c	    DCt(:,intp) = min( dim(fact * sqrt(raLSt * gnormt),
	    DCt(:,intp) = min( max(zero,fact * sqrt(raLSt * gnormt)-
     &                       rtLSt * gnormt), two * rtLSt * gnormt )
c
c.... flop count
c
!	    flops = flops + 48*npro
c
	  endif
c
c	endif
c	DCt=DCt*two
c
c.... ---------------------------->  RHS  <----------------------------
c
c.... add the contribution of DC to ri and/or rmi
c
c.... ires = 1 or 3
c
	if ((ires .eq. 1) .or. (ires .eq. 3)) then
c
	  rti ( :,1) = rti ( :,1) + DCt(:,intp) * gAgyit( :,1)
	  rmti( :,1) = rmti( :,1) + DCt(:,intp) * gAgyit( :,1)
	  rti ( :,2) = rti ( :,2) + DCt(:,intp) * gAgyit( :,2)
	  rmti( :,2) = rmti( :,2) + DCt(:,intp) * gAgyit( :,2)
	  rti ( :,3) = rti ( :,3) + DCt(:,intp) * gAgyit( :,3)
	  rmti( :,3) = rmti( :,3) + DCt(:,intp) * gAgyit( :,3)
	 
c
!	  flops = flops + 45*npro
c
	endif
c
c.... ires = 2
c
	if (ires .eq. 2) then
c
	  rmti( :,1) = rmti( :,1) + DCt(:,intp) * gAgyit( :,1)
	  rmti( :,2) = rmti( :,2) + DCt(:,intp) * gAgyit( :,2)
	  rmti( :,3) = rmti( :,3) + DCt(:,intp) * gAgyit( :,3)
	  
c
!	  flops = flops + 30*npro
c
	endif
c
c.... ------------------------->  Stiffness  <--------------------------
c
c.... add the contribution of DC to stiff
c$$$c
c	if (iprec .eq. 1) then !leave out of LHS, if called from itrres
	                       !anyway matrix free not implemented for scalar
             dtmp(:)=A0t(:)*DCt(:,intp)	    
c
c.... add (DC g^1 A0) to stifft [1,1]
c
             stifft(:,1,1) = stifft(:,1,1) 
     &                    + dtmp(:)*giju(:,1)
c
c.... add (DC g^1 A0) to stifft [1,2]
c
             stifft(:,1,2) = stifft(:,1,2) 
     &                    + dtmp(:)*giju(:,4)
c
c.... add (DC g^1 A0) to stifft [1,3]
c
             stifft(:,1,3) = stifft(:,1,3) 
     &                    + dtmp(:)*giju(:,5)

c.... add (DC g^1 A0) to stifft [2,1] (similarly below)
c
             stifft(:,2,1) = stifft(:,2,1) 
     &                    + dtmp(:)*giju(:,4)

             stifft(:,2,2) = stifft(:,2,2) 
     &                    + dtmp(:)*giju(:,2)

             stifft(:,2,3) = stifft(:,2,3) 
     &                    + dtmp(:)*giju(:,6)

             stifft(:,3,1) = stifft(:,3,1) 
     &                    + dtmp(:)*giju(:,5)

             stifft(:,3,2) = stifft(:,3,2) 
     &                    + dtmp(:)*giju(:,6)

             stifft(:,3,3) = stifft(:,3,3) 
     &                    + dtmp(:)*giju(:,3)
c
c.... flop count
c
!	  flops = flops + 210*npro
c
c.... end of stiffness
c
c	endif
c
c.... return
c
	return
	end
c
