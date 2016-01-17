      subroutine bc3Res (y,  iBC,  BC,  res, iper, ilwork)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the residual vector for 3D elements.
c
c input:
c  y     (nshg,ndof)   : Y variables 
c  iBC   (nshg)        : Boundary Condition Code
c  BC    (nshg,ndofBC) : the boundary condition constraint parameters
c  res   (nshg,nflow)   : residual before BC is applied
c
c output:
c  res   (nshg,nflow)   : residual after satisfaction of BC
c  
c
c Thuc Bui,      Winter 1989.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      include "common.h"
c     
      dimension y(nshg,ndof),             iBC(nshg),
     &            BC(nshg,ndofBC),         
     &            res(nshg,nflow),           ilwork(nlwork),
     &            iper(nshg)
c
c.... density 
c
      where (btest(iBC,0))
         res(:,5) = res(:,5) + BC(:,1)*Rgas* res(:,1) !IDEAL GAS ASSUMED
         res(:,1) = zero
      endwhere
c 
c.... pressure 
c
      if(EntropyPressure.eq.1) then
         where (btest(iBC,2))

c Thought this would be correct if W was in the tangent space of V
c as described in Shakib's thesis it does not work for primitive
c variables
c
c Instead we use the entropy variable W
c
            res(:,2) = res(:,2) -y(:,1)*res(:,1)
            res(:,3) = res(:,3) -y(:,2)*res(:,1)
            res(:,4) = res(:,4) -y(:,3)*res(:,1)
            res(:,5) = res(:,5) -
     &           (gamma*Rgas/gamma1*y(:,5) 
     &           + pt5 * ( y(:,1)**2 + y(:,2)**2 + y(:,3)**2))*res(:,1)
            res(:,1) = zero
         endwhere
      else
         where (btest(iBC,2))
            res(:,1) = zero
         endwhere
      endif
c
c.... velocities
c
c ibits(n1,n2,n3) extracts bits n2+1 through n2+n3 (extending to the left
c as is traditional in binary) of the integer n1
c and returns the base 10 integer. In examples below x y z a b can 
c be 1 or zero without any effect.
c
c.... x1-velocity
c
c if iBC=4   bits of ibc =00000100 => ibits(4,3,3)=0
c if iBC=40  bits of ibc =00101000 => ibits(4,3,3)=5
c if iBC=40  bits of ibc =00101000 => ibits(4,3,2)=1
c
        where (ibits(iBC,3,3) .eq. 1)   ! bits of iBC= xy001zab 
c
c     notice that the extracted 3 bits form the number 1.  below
c     you will see the combinations which make up 2-7, all of the
c     possible velocity combinations
c
          res(:,3) = res(:,3) - BC(:,4) * res(:,2)
          res(:,4) = res(:,4) - BC(:,5) * res(:,2)
          res(:,2) = zero
        endwhere
c
c.... x2-velocity
c
        where (ibits(iBC,3,3) .eq. 2)   ! bits of iBC= xy010zab 
          res(:,2) = res(:,2) - BC(:,4) * res(:,3)
          res(:,4) = res(:,4) - BC(:,5) * res(:,3)
          res(:,3) = zero
        endwhere
c
c.... x1-velocity and x2-velocity
c
        where (ibits(iBC,3,3) .eq. 3)  ! bits of iBC= xy011zab 
          res(:,4) = res(:,4) - BC(:,4) * res(:,2) - BC(:,6) * res(:,3)
          res(:,2) = zero
          res(:,3) = zero
        endwhere
c
c.... x3-velocity
c
        where (ibits(iBC,3,3) .eq. 4)  ! bits of iBC= xy100zab 
          res(:,2) = res(:,2) - BC(:,4) * res(:,4)
          res(:,3) = res(:,3) - BC(:,5) * res(:,4)
          res(:,4) = zero
        endwhere
c
c.... x1-velocity and x3-velocity
c
        where (ibits(iBC,3,3) .eq. 5)  ! bits of iBC= xy101zab 
          res(:,3) = res(:,3) - BC(:,4) * res(:,2) - BC(:,6) * res(:,4)
          res(:,2) = zero
          res(:,4) = zero
        endwhere
c
c.... x2-velocity and x3-velocity
c
        where (ibits(iBC,3,3) .eq. 6)  ! bits of iBC= xy110zab 
          res(:,2) = res(:,2) - BC(:,4) * res(:,3) - BC(:,6) * res(:,4)
          res(:,3) = zero
          res(:,4) = zero
        endwhere
c
c.... x1-velocity, x2-velocity and x3-velocity
c
        where (ibits(iBC,3,3) .eq. 7)  ! bits of iBC= xy111zab 
          res(:,2) = zero
          res(:,3) = zero
          res(:,4) = zero
        endwhere
c
c.... scaled plane extraction boundary condition
c
        if(intpres.eq.1) then  ! interpolating pressure so zero continuity res 
           where (btest(iBC,11))
              res(:,1) = zero
              res(:,2) = zero
              res(:,3) = zero
              res(:,4) = zero
	      res(:,5) = zero ! added to correspond to genscale (Elaine)
           endwhere
        else  ! leave residual in continuity equation
           where (btest(iBC,11))
              res(:,2) = zero
              res(:,3) = zero
              res(:,4) = zero
	      res(:,5) = zero ! added to correspond to genscale (Elaine)
           endwhere
        endif
c
c.... temperature
c
        where (btest(iBC,1)) res(:,5) = zero
c
c.... local periodic boundary conditions (no communications)
c
        do j = 1,nshg
          if (btest(iBC(j),10)) then
            i = iper(j)
            res(i,:) = res(i,:) + res(j,:)
            res(j,:) = zero          
          endif
        enddo
c
c.... periodic slaves get the residual values of the masters
c
c        do i = 1,nshg
c           if (btest(iBC(i),10)) then
c              res(i,:) = res(iper(i),:)
c           endif
c        enddo       
       if(numpe.gt.1) then
       if(usingPETSc.eq.0) then !kill this code for petsc
c
c.... nodes treated on another processor are eliminated
c
        numtask = ilwork(1)
        itkbeg = 1

        do itask = 1, numtask

          iacc   = ilwork (itkbeg + 2)
          numseg = ilwork (itkbeg + 4)

          if (iacc .eq. 0) then
            do is = 1,numseg
              isgbeg = ilwork (itkbeg + 3 + 2*is)
              lenseg = ilwork (itkbeg + 4 + 2*is)
              isgend = isgbeg + lenseg - 1
              res(isgbeg:isgend,:) = zero
            enddo
          endif
          
          itkbeg = itkbeg + 4 + 2*numseg

        enddo
        endif
        endif
c
c.... return
c
        return
        end
c
c
c
        subroutine bc3ResSclr (y,  iBC,  BC,   rest, 
     &                         iper, ilwork)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the residual vector for 3D elements.
c
c input:
c  Y     (nshg,ndof)   : Y Variables
c  iBC   (nshg)        : Boundary Condition Code
c  BC    (nshg,ndofBC) : the boundary condition constraint parameters
c  rest  (nshg)        : residual before BC is applied
c
c output:
c  rest  (nshg)        : residual after satisfaction of BC
c  
c
c Thuc Bui,      Winter 1989.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,ndof),      iBC(nshg),
     &            BC(nshg,ndofBC),     
     &            rest(nshg),        ilwork(nlwork),
     &            iper(nshg) 
c
c    
         id = isclr+5
c.... Scalar Variable
c
        where (btest(iBC,id))
          rest(:) = zero
        endwhere
c
c.... local periodic boundary conditions (no communications)
c
        do j = 1,nshg
          if (btest(iBC(j),10)) then
            i = iper(j)
            rest(i) = rest(i) + rest(j)
            rest(j) = zero   !changed
          endif
        enddo
c
c.... periodic slaves get the residual values of the masters
c
c$$$         do i = 1,nshg
c$$$            if (btest(iBC(i),10)) then
c$$$               rest(i) = rest(iper(i))
c$$$            endif
c$$$         enddo 
c
c removed for impl=4 as we have set the loops over ntopsh
c
        if(numpe.gt.1 ) then
         if(usingPETSc.eq.0) then !kill this code for petsc
c     
c.... nodes treated on another processor are eliminated
c     
           numtask = ilwork(1)
           itkbeg = 1
           
           do itask = 1, numtask
              
              iacc   = ilwork (itkbeg + 2)
              numseg = ilwork (itkbeg + 4)
              
              if (iacc .eq. 0) then
                 do is = 1,numseg
                    isgbeg = ilwork (itkbeg + 3 + 2*is)
                    lenseg = ilwork (itkbeg + 4 + 2*is)
                    isgend = isgbeg + lenseg - 1
                    rest(isgbeg:isgend) = zero
                 enddo
              endif
              
              itkbeg = itkbeg + 4 + 2*numseg
              
           enddo
         endif
        endif
c
c
c.... return
c
        return
        end

