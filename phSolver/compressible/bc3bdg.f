        subroutine bc3BDg (y,  iBC,  BC, BDiag, iper, ilwork)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the block-diagonal preconditioning
c   matrix for 3D elements.
c
c input:
c  y      (nshg,ndof)   : Y variables
c  iBC    (nshg)        : boundary condition code
c  BC     (nshg,ndofBC) : Dirichlet BC constraint parameters
c  BDiag   (nshg,nflow,nflow) : preconditionning matrix before BC
c                          (only upper part)
c
c output:
c  BDiag   (nshg,nflow,nflow) : preconditionning matrix after BC
c                          is satisfied
c
c
c Zdenek Johan, Summer 1990. (Modified from g3bce.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,ndof),             iBC(nshg),
     &            BC(nshg,ndofBC),          
     &            BDiag(nshg,nflow,nflow),        ilwork(nlwork),
     &            iper(nshg)
c
        real*8 a5(nshg)
c
c.... density 
c
          do i = 1, nshg 
             a5(i) = - y(i,5) * (Rgas * gamma / gamma1) !IDEAL GAS ASSUMED 
          end do

        where (btest(iBC,0))
c
c.... engbc was replaced for a5 by following

          BDiag(:,5,5) = BDiag(:,5,5) + a5 * a5 * BDiag(:,1,1) +
     &                                       a5 * BDiag(:,1,5) +
     &                                       a5 * BDiag(:,5,1)
          BDiag(:,4,5) = BDiag(:,4,5) +  a5 * BDiag(:,4,1) 
          BDiag(:,3,5) = BDiag(:,3,5) +  a5 * BDiag(:,3,1) 
          BDiag(:,2,5) = BDiag(:,2,5) +  a5 * BDiag(:,2,1) 
          BDiag(:,5,4) = BDiag(:,5,4) +  a5 * BDiag(:,1,4) 
          BDiag(:,5,3) = BDiag(:,5,3) +  a5 * BDiag(:,1,3) 
          BDiag(:,5,2) = BDiag(:,5,2) +  a5 * BDiag(:,1,2) 
          BDiag(:,1,2) = zero
          BDiag(:,1,3) = zero
          BDiag(:,1,4) = zero
          BDiag(:,1,5) = zero
          BDiag(:,2,1) = zero
          BDiag(:,3,1) = zero
          BDiag(:,4,1) = zero
          BDiag(:,5,1) = zero
          BDiag(:,1,1) = one
        endwhere

c       where (btest(iBC,11)) ! pressure  to deactivate
        where (btest(iBC,2)) ! pressure

          BDiag(:,1,2) = zero
          BDiag(:,1,3) = zero
          BDiag(:,1,4) = zero
          BDiag(:,1,5) = zero
          BDiag(:,2,1) = zero
          BDiag(:,3,1) = zero
          BDiag(:,4,1) = zero
          BDiag(:,5,1) = zero
          BDiag(:,1,1) = one
        endwhere

c
c.... velocities
c
c.... x1-velocity   
c
        where (ibits(iBC,3,3) .eq. 1)
          BDiag(:,5,4) = BDiag(:,5,4) - BC(:,5) * BDiag(:,5,2)
          BDiag(:,5,3) = BDiag(:,5,3) - BC(:,4) * BDiag(:,5,2)

          BDiag(:,4,5) = BDiag(:,4,5) - BC(:,5) * BDiag(:,2,5)
          BDiag(:,3,5) = BDiag(:,3,5) - BC(:,4) * BDiag(:,2,5)

          BDiag(:,4,1) = BDiag(:,4,1) - BC(:,5) * BDiag(:,2,1)
          BDiag(:,3,1) = BDiag(:,3,1) - BC(:,4) * BDiag(:,2,1)

          BDiag(:,1,4) = BDiag(:,1,4) - BC(:,5) * BDiag(:,1,2)
          BDiag(:,1,3) = BDiag(:,1,3) - BC(:,4) * BDiag(:,1,2)

          BDiag(:,4,4) = BDiag(:,4,4) + BC(:,5) * BC(:,5) * BDiag(:,2,2)
     &                                -           BC(:,5) * BDiag(:,2,4)
     &                                -           BC(:,5) * BDiag(:,4,2)
          BDiag(:,3,4) = BDiag(:,3,4) + BC(:,4) * BC(:,5) * BDiag(:,2,2)
     &                                -           BC(:,5) * BDiag(:,3,2)
     &                                -           BC(:,4) * BDiag(:,2,4)
          BDiag(:,4,3) = BDiag(:,4,3) + BC(:,4) * BC(:,5) * BDiag(:,2,2)
     &                                -           BC(:,5) * BDiag(:,2,3)
     &                                -           BC(:,4) * BDiag(:,4,2)
          BDiag(:,3,3) = BDiag(:,3,3) + BC(:,4) * BC(:,4) * BDiag(:,2,2)
     &                                -           BC(:,4) * BDiag(:,2,3)
     &                                -           BC(:,4) * BDiag(:,3,2)
          BDiag(:,2,1) = zero
          BDiag(:,1,2) = zero
          BDiag(:,2,3) = zero
          BDiag(:,2,4) = zero
          BDiag(:,2,5) = zero
          BDiag(:,3,2) = zero
          BDiag(:,4,2) = zero
          BDiag(:,5,2) = zero
          BDiag(:,2,2) = one
        endwhere
c
c.... x2-velocity
c
        where (ibits(iBC,3,3) .eq. 2)
          BDiag(:,5,4) = BDiag(:,5,4) - BC(:,5) * BDiag(:,5,3)
          BDiag(:,5,2) = BDiag(:,5,2) - BC(:,4) * BDiag(:,5,3)

          BDiag(:,4,5) = BDiag(:,4,5) - BC(:,5) * BDiag(:,3,5)
          BDiag(:,2,5) = BDiag(:,2,5) - BC(:,4) * BDiag(:,3,5)

          BDiag(:,4,1) = BDiag(:,4,1) - BC(:,5) * BDiag(:,3,1)
          BDiag(:,2,1) = BDiag(:,2,1) - BC(:,4) * BDiag(:,3,1)

          BDiag(:,1,4) = BDiag(:,1,4) - BC(:,5) * BDiag(:,1,3)
          BDiag(:,1,2) = BDiag(:,1,2) - BC(:,4) * BDiag(:,1,3)

          BDiag(:,4,4) = BDiag(:,4,4) + BC(:,5) * BC(:,5) * BDiag(:,3,3)
     &                                -           BC(:,5) * BDiag(:,3,4)
     &                                -           BC(:,5) * BDiag(:,4,3)
          BDiag(:,2,4) = BDiag(:,2,4) + BC(:,4) * BC(:,5) * BDiag(:,3,3)
     &                                - 	  BC(:,5) * BDiag(:,2,3)
     &                                - 	  BC(:,4) * BDiag(:,3,4) 
          BDiag(:,4,2) = BDiag(:,4,2) + BC(:,4) * BC(:,5) * BDiag(:,3,3)
     &                                - 	  BC(:,5) * BDiag(:,3,2)
     &                                - 	  BC(:,4) * BDiag(:,4,3)
          BDiag(:,2,2) = BDiag(:,2,2) + BC(:,4) * BC(:,4) * BDiag(:,3,3)
     &                                - 	  BC(:,4) * BDiag(:,2,3)
     &                                - 	  BC(:,4) * BDiag(:,3,2)
          BDiag(:,3,1) = zero
          BDiag(:,3,2) = zero
          BDiag(:,3,4) = zero
          BDiag(:,3,5) = zero
          BDiag(:,1,3) = zero
          BDiag(:,2,3) = zero
          BDiag(:,4,3) = zero
          BDiag(:,5,3) = zero
          BDiag(:,3,3) = one
        endwhere
c
c.... x1-velocity and x2-velocity 
c
      where (ibits(iBC,3,3) .eq. 3) 
      BDiag(:,4,4) = BDiag(:,4,4) + BC(:,4) * BC(:,4) * BDiag(:,2,2) 
     &                            + BC(:,6) * BC(:,6) * BDiag(:,3,3)
     &           + BC(:,4) * BC(:,6) * ( BDiag(:,2,3) * BDiag(:,3,2))
     &           -           BC(:,6) * ( BDiag(:,4,3) * BDiag(:,3,4))
     &           -           BC(:,4) * ( BDiag(:,4,2) * BDiag(:,2,4))
      BDiag(:,1,4) = BDiag(:,1,4) -           BC(:,4) * BDiag(:,1,2) 
     &                            -           BC(:,6) * BDiag(:,1,3)
      BDiag(:,4,1) = BDiag(:,4,1) -           BC(:,4) * BDiag(:,2,1) 
     &                            -           BC(:,6) * BDiag(:,3,1)
      BDiag(:,5,4) = BDiag(:,5,4) -           BC(:,4) * BDiag(:,5,2) 
     &                            -           BC(:,6) * BDiag(:,5,3)
      BDiag(:,4,5) = BDiag(:,4,5) -           BC(:,4) * BDiag(:,2,5) 
     &                            -           BC(:,6) * BDiag(:,3,5)
          BDiag(:,2,1) = zero
          BDiag(:,2,3) = zero
          BDiag(:,2,4) = zero
          BDiag(:,2,5) = zero
          BDiag(:,3,1) = zero
          BDiag(:,3,2) = zero
          BDiag(:,3,4) = zero
          BDiag(:,3,5) = zero
          BDiag(:,1,2) = zero
          BDiag(:,4,2) = zero
          BDiag(:,5,2) = zero
          BDiag(:,1,3) = zero
          BDiag(:,4,3) = zero
          BDiag(:,5,3) = zero
          BDiag(:,3,3) = one
          BDiag(:,2,2) = one
        endwhere
c
c.... x3-velocity
c
        where (ibits(iBC,3,3) .eq. 4)
          BDiag(:,5,3) = BDiag(:,5,3) - BC(:,5) * BDiag(:,5,4)
          BDiag(:,5,2) = BDiag(:,5,2) - BC(:,4) * BDiag(:,5,4)

          BDiag(:,3,5) = BDiag(:,3,5) - BC(:,5) * BDiag(:,4,5)
          BDiag(:,2,5) = BDiag(:,2,5) - BC(:,4) * BDiag(:,4,5)

          BDiag(:,3,1) = BDiag(:,3,1) - BC(:,5) * BDiag(:,4,1)
          BDiag(:,2,1) = BDiag(:,2,1) - BC(:,4) * BDiag(:,4,1)

          BDiag(:,1,3) = BDiag(:,1,3) - BC(:,5) * BDiag(:,1,4)
          BDiag(:,1,2) = BDiag(:,1,2) - BC(:,4) * BDiag(:,1,4)

          BDiag(:,3,3) = BDiag(:,3,3) + BC(:,5) * BC(:,5) * BDiag(:,4,4)
     &                                - 	  BC(:,5) * BDiag(:,3,4)
     &                                - 	  BC(:,5) * BDiag(:,4,3)
          BDiag(:,2,3) = BDiag(:,2,3) + BC(:,4) * BC(:,5) * BDiag(:,4,4)
     &                                - 	  BC(:,5) * BDiag(:,2,4)
     &                                -	          BC(:,4) * BDiag(:,4,3)
          BDiag(:,3,2) = BDiag(:,3,2) + BC(:,4) * BC(:,5) * BDiag(:,4,4)
     &                                -   	  BC(:,5) * BDiag(:,4,2)
     &                                - 	  BC(:,4) * BDiag(:,3,4)
          BDiag(:,2,2) = BDiag(:,2,2) + BC(:,4) * BC(:,4) * BDiag(:,4,4)
     &                                - 	  BC(:,4) * BDiag(:,2,4)
     &                                - 	  BC(:,4) * BDiag(:,4,2)
          BDiag(:,4,1) = zero
          BDiag(:,4,2) = zero
          BDiag(:,4,3) = zero
          BDiag(:,4,5) = zero
          BDiag(:,1,4) = zero
          BDiag(:,2,4) = zero
          BDiag(:,3,4) = zero
          BDiag(:,5,4) = zero
          BDiag(:,4,4) = one
        endwhere
c
c.... x1-velocity and x3-velocity
c
        where (ibits(iBC,3,3) .eq. 5)
          BDiag(:,3,3) = BDiag(:,3,3) + BC(:,4) * BC(:,4) * BDiag(:,2,2)
     &                                + BC(:,6) * BC(:,6) * BDiag(:,4,4)
     &                + BC(:,4) * BC(:,6) *(BDiag(:,2,4) + BDiag(:,4,2))
     &                -           BC(:,4) *(BDiag(:,2,3) + BDiag(:,3,2))
     &                -           BC(:,6) *(BDiag(:,4,3) + BDiag(:,3,4))
          BDiag(:,1,3) = BDiag(:,1,3) -           BC(:,4) * BDiag(:,1,2) 
     &                                -           BC(:,6) * BDiag(:,1,4)
          BDiag(:,3,1) = BDiag(:,3,1) -           BC(:,4) * BDiag(:,2,1) 
     &                                -           BC(:,6) * BDiag(:,4,1)
          BDiag(:,5,3) = BDiag(:,5,3) -           BC(:,4) * BDiag(:,5,2) 
     &                                -           BC(:,6) * BDiag(:,5,4)
          BDiag(:,3,5) = BDiag(:,3,5) -           BC(:,4) * BDiag(:,2,5) 
     &                                -           BC(:,6) * BDiag(:,4,5)
          BDiag(:,2,1) = zero
          BDiag(:,2,3) = zero
          BDiag(:,2,4) = zero
          BDiag(:,2,5) = zero
          BDiag(:,4,1) = zero
          BDiag(:,4,2) = zero
          BDiag(:,4,3) = zero
          BDiag(:,4,5) = zero
          BDiag(:,1,2) = zero
          BDiag(:,4,2) = zero
          BDiag(:,5,2) = zero
          BDiag(:,1,4) = zero
          BDiag(:,3,4) = zero
          BDiag(:,5,4) = zero
          BDiag(:,4,4) = one
          BDiag(:,2,2) = one
        endwhere
c
c.... x2-velocity and x3-velocity
c
        where (ibits(iBC,3,3) .eq. 6)
          BDiag(:,2,2) = BDiag(:,2,2) + BC(:,4) * BC(:,4) * BDiag(:,3,3)
     &                                + BC(:,6) * BC(:,6) * BDiag(:,4,4)
     &               + BC(:,4) * BC(:,6) * (BDiag(:,3,4) + BDiag(:,4,3))
     &               -            BC(:,4) *(BDiag(:,2,3) + BDiag(:,3,2))
     &               -            BC(:,6) *(BDiag(:,4,2) + BDiag(:,2,4))
          BDiag(:,1,2) = BDiag(:,1,2) -           BC(:,4) * BDiag(:,1,3) 
     &                                -           BC(:,6) * BDiag(:,1,4)
          BDiag(:,2,1) = BDiag(:,2,1) -           BC(:,4) * BDiag(:,3,1) 
     &                                -           BC(:,6) * BDiag(:,4,1)
          BDiag(:,5,2) = BDiag(:,5,2) -           BC(:,4) * BDiag(:,5,3) 
     &                                -           BC(:,6) * BDiag(:,5,4)
          BDiag(:,2,5) = BDiag(:,2,5) -           BC(:,4) * BDiag(:,3,5) 
     &                                -           BC(:,6) * BDiag(:,4,5)
          BDiag(:,3,1) = zero
          BDiag(:,3,2) = zero
          BDiag(:,3,4) = zero
          BDiag(:,3,5) = zero
          BDiag(:,4,1) = zero
          BDiag(:,4,2) = zero
          BDiag(:,4,3) = zero
          BDiag(:,4,5) = zero
          BDiag(:,1,3) = zero
          BDiag(:,2,3) = zero
          BDiag(:,5,3) = zero
          BDiag(:,1,4) = zero
          BDiag(:,2,4) = zero
          BDiag(:,5,4) = zero
          BDiag(:,4,4) = one
          BDiag(:,3,3) = one
        endwhere
c
c.... x1-velocity and x2-velocity and x3-velocity
c
        where (ibits(iBC,3,3) .eq. 7)
          BDiag(:,2,1) = zero
          BDiag(:,2,3) = zero
          BDiag(:,2,4) = zero
          BDiag(:,2,5) = zero
          BDiag(:,3,1) = zero
          BDiag(:,3,2) = zero
          BDiag(:,3,4) = zero
          BDiag(:,3,5) = zero
          BDiag(:,4,1) = zero
          BDiag(:,4,2) = zero
          BDiag(:,4,3) = zero
          BDiag(:,4,5) = zero
          BDiag(:,1,2) = zero
          BDiag(:,5,2) = zero
          BDiag(:,1,3) = zero
          BDiag(:,5,3) = zero
          BDiag(:,1,4) = zero
          BDiag(:,5,4) = zero
          BDiag(:,2,2) = one
          BDiag(:,3,3) = one
          BDiag(:,4,4) = one
        endwhere
c
c.... temperature
c
        where (btest(iBC,1))
          BDiag(:,5,5) = one
          BDiag(:,1,5) = zero
          BDiag(:,2,5) = zero
          BDiag(:,3,5) = zero
          BDiag(:,4,5) = zero
          BDiag(:,5,1) = zero
          BDiag(:,5,2) = zero
          BDiag(:,5,3) = zero
          BDiag(:,5,4) = zero
        endwhere
c
c.... local periodic boundary conditions (no communications)
c

           do j = 1,nshg
              if (btest(iBC(j),10)) then
                 i = iper(j)
                 BDiag(i, :,:) = BDiag(i,:,:) + BDiag(j,:,:)
              endif
           enddo
c     
c.... periodic slaves get the residual values of the masters
c
           do j = 1,nshg
              if (btest(iBC(j),10)) then
                 i=iper(j)
                 BDiag(j,:,:) = BDiag(i,:,:)
              endif
           enddo  
c$$$        endif

        if(numpe.gt.1) then
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
              BDiag(isgbeg:isgend,:,:) = zero
              BDiag(isgbeg:isgend,1,1) = one
              BDiag(isgbeg:isgend,2,2) = one
              BDiag(isgbeg:isgend,3,3) = one
              BDiag(isgbeg:isgend,4,4) = one
              BDiag(isgbeg:isgend,5,5) = one
            enddo
          endif
          
          itkbeg = itkbeg + 4 + 2*numseg

        enddo
        endif
c
c.... return
c
        return
        end
c
c
c
        subroutine bc3BDgSclr (iBC, Diag, iper, ilwork)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the block-diagonal preconditioning
c   matrix for 3D elements.
c
c input:
c  iBC    (numnp)        : boundary condition code
c  BC     (numnp,ndofBC) : Dirichlet BC constraint parameters
c  Diag   (numnp) : preconditionning diagonal matrix before BC
c
c output:
c  Diag   (numnp) : preconditionning matrix after BC
c                          is satisfied
c
c
c Zdenek Johan, Summer 1990. (Modified from g3bce.f)
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension iBC(nshg),
     &            Diag(nshg),             ilwork(nlwork),
     &            iper(nshg)
c
c
      id = isclr+5
c
c.... scalar variable 
c
        where (btest(iBC,id)) 
          Diag(:) = one
        endwhere
c
c.... local periodic boundary conditions (no communications)
c
        do j = 1,nshg
          if (btest(iBC(j),10)) then
            i = iper(j)
            Diag(i) = Diag(i) + Diag(j)
          endif
        enddo
c     
c.... periodic slaves get the residual values of the masters
c
      do i = 1,nshg
         if (btest(iBC(i),10)) then
            Diag(i) = Diag(iper(i))
         endif
      enddo       
c
c.... nodes treated on another processor are eliminated
c
      if(numpe.gt.1) then
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
                  Diag(isgbeg:isgend) = one
               enddo
            endif
          
            itkbeg = itkbeg + 4 + 2*numseg

         enddo
      endif
c
c.... return
c
        return
        end



