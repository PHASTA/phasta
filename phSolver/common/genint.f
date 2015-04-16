      subroutine genint
c
c----------------------------------------------------------------------
c
c This subroutine inputs the integration information.
c
c----------------------------------------------------------------------
c
      include "common.h"

      real*8, allocatable :: tmpQpt (:,:), tmpQwt (:) 
      real*8, allocatable :: tmpQptb(:,:), tmpQwtb(:) 

c
c.... compute the shape function parameters
c
        

c 
c.... get quadrature data for interior and boundary elements
c
c
c  Tets
c             
      if (nen.eq.4) then 
         nshape  = (ipord+1)*(ipord+2)*(ipord+3)/6
         nshapeb = (ipord+1)*(ipord+2)/2
      endif

      select case (intg(1,1))
      case (1)
         nint(1) = 1
      case (2)
         nint(1) = 4
      case (3)
         nint(1) = 16
      case (4)
         nint(1) = 29
      end select
      select case (intg(2,1))
      case (1)
         nintb(1) = 1
      case (2)
         nintb(1) = 3
      case (3)
         nintb(1) = 6
      case (4)
         nintb(1) = 12
      end select
      
      allocate (tmpQpt (4,nint(1)))           
      allocate (tmpQwt (nint(1)))           
      allocate (tmpQptb(4,nintb(1)))           
      allocate (tmpQwtb(nintb(1)))           
      
      call symtet(nint(1),tmpQpt,tmpQwt,nerr) ! interior elements
      Qpt(1,1:4,1:nint(1)) = tmpQpt(1:4,1:nint(1))
      Qwt(1,1:nint(1))     = tmpQwt(1:nint(1))
      
      
      call symtri(nintb(1),tmpQptb,tmpQwtb,nerr) ! boundary elements
      Qptb(1,1:4,1:nintb(1)) = tmpQptb(1:4,1:nintb(1))
      Qwtb(1,1:nintb(1))     = tmpQwtb(1:nintb(1))
      
      deallocate (tmpQpt)
      deallocate (tmpQwt)
      deallocate (tmpQptb)
      deallocate (tmpQwtb)
        
c
c.... adjust quadrature weights to be consistent with the
c     design of tau. 
c
      Qwt(1,:) = (four/three)*Qwt(1,:)
      Qwtb(1,:) = two*Qwtb(1,:)
c     
c     Hexes now
c     
      if (nen.eq.8) then
         nshape = nen
         if ( ipord .gt. 1 ) then 
            nshape = nshape + 12*(ipord - 1)
         endif
         
         if ( ipord .gt. 3 ) then
            nshape = nshape + 3*(ipord -2)*(ipord - 3)
         endif
         
         if ( ipord .gt. 5 ) then
            nshape = nshape + (ipord - 3)*(ipord - 4)*(ipord - 5)/6
         endif
         
         nshapeb = nenb 
         
         if ( ipord .gt. 1 ) then
            nshapeb = nshapeb + 4*(ipord - 1)
         endif
         
         if ( ipord .gt. 3 ) then 
            nshapeb = nshapeb +(ipord -2)*(ipord - 3)/2
         endif
      endif
      
      select case (intg(1,1))
      case (1)
         nint(2) = 1
      case (2)
         nint(2) = 8
      case (3)
         nint(2) =27
      case (4)
         nint(2) = 64
      case (5)
         nint(2) = 125
      end select
      
      select case (intg(2,1))
      case (1)
         nintb(2) = 1 
      case (2) 
         nintb(2) = 4
      case(3)
         nintb(2) = 9
      case (4)
         nintb(2) = 16
      case (5) 
         nintb(2) = 25
      end select
      
      
      allocate (tmpQpt (4,nint(2)))           
      allocate (tmpQwt (nint(2)))           
      allocate (tmpQptb(4,nintb(2)))           
      allocate (tmpQwtb(nintb(2)))      
      
      call symhex(nint(2),tmpQpt,tmpQwt,nerr)
      Qpt(2,1:4,1:nint(2)) = tmpQpt(1:4,1:nint(2))
      Qwt(2,1:nint(2)) = tmpQwt(1:nint(2))
      
      call symquad(nintb(2),tmpQptb,tmpQwtb,nerr)
      Qptb(2,1:4,1:nintb(2)) = tmpQptb(1:4,1:nintb(2))
      Qwtb(2,1:nintb(2)) = tmpQwtb(1:nintb(2))
      
      deallocate (tmpQpt)
      deallocate (tmpQwt)
      deallocate (tmpQptb)
      deallocate (tmpQwtb)
      
      if (nen.eq.5) then        ! for pyramid
         
         nshape = nen           ! counting the nodal functions
         
         if ( ipord .gt. 1 ) then 
            nshape = nshape + 8*(ipord - 1) ! counting the edge functions
         endif
           
         if ( ipord .gt. 2 ) then            
                                ! counting the triangular face functions
            nshape = nshape + 2*(ipord -1)*(ipord - 2)
         endif
         
         if ( ipord .gt. 3 ) then
                                ! counting the quadrilateral face functions
            nshape = nshape + (ipord - 2)*(ipord - 3)/2
         endif
         
         if ( ipord .gt. 5 ) then
            nshape = nshape + (ipord - 3)*(ipord - 4)*(ipord - 5)/6
         endif
         
         nshapeb = nenb
                                ! assume only the quadrilateral face on boundary
         if ( ipord .gt. 1 ) then
            nshapeb = nshapeb + 4*(ipord - 1)
         endif
         if ( ipord .gt. 3 ) then 
            nshape = nshape +(ipord - 2)*(ipord - 3)/2
         endif  
         
      endif
      
      select case (intg(1,1))
      case (1)
         nint(5) = 1
      case (2)
         nint(5) = 8
      case (3)
         nint(5) = 27
      case (4)
         nint(5) = 64
      case (5)
         nint(5) = 125
      end select
      select case (intg(2,1))
      case (1)
         nintb(5) = 1
      case (2)
         nintb(5) = 4
      case (3)
         nintb(5) = 9
      case (4)
         nintb(5) = 16
      case (5)
         nintb(5) = 25
      end select
      select case (intg(2,1))
      case (1)
         nintb(6) = 1
      case (2)
         nintb(6) = 4
      case (3)
         nintb(6) = 9
      case (4)
         nintb(6) = 12
      case (5)
         nintb(6) = 12 ! tri face boundary integration rules not coded
      end select
      
      allocate (tmpQpt (4,nint(5)))           
      allocate (tmpQwt (nint(5)))
      
      call sympyr(nint(5),tmpQpt,tmpQwt,nerr) ! interior elements
      Qpt(5,1:4,1:nint(5)) = tmpQpt(1:4,1:nint(5))
      Qwt(5,1:nint(5)) = tmpQwt(1:nint(5))
      
      deallocate (tmpQpt)
      deallocate (tmpQwt)
      
      allocate (tmpQptb(4,nintb(5)))
      allocate (tmpQwtb(nintb(5))) 
      
      call symquad (nintb(5),tmpQptb,tmpQwtb,nerr) ! quad boundary elements
      Qptb(5,1:4,1:nintb(5)) = tmpQptb(1:4,1:nintb(5))
      Qwtb(5,1:nintb(5)) = tmpQwtb(1:nintb(5))    
      
      deallocate (tmpQptb)
      deallocate (tmpQwtb)
      
      allocate (tmpQptb(4,nintb(6)))
      allocate (tmpQwtb(nintb(6))) 
      
      call symtripyr (nintb(6),tmpQptb,tmpQwtb,nerr) ! tri boundary elements
      Qptb(6,1:4,1:nintb(6)) = tmpQptb(1:4,1:nintb(6))
      Qwtb(6,1:nintb(6)) = tmpQwtb(1:nintb(6))    
      
      deallocate (tmpQptb)
      deallocate (tmpQwtb)
      
      if (nen.eq.6) then
c     
c     later for hierarchic basis nshape formulas will be inserted
c     
         nshape = nen
         
         if ( ipord .gt. 1 ) then 
            nshape = nshape + 9*(ipord - 1)
         endif
         
         if ( ipord .gt. 2 ) then
            nshape = nshape + (ipord -1)*(ipord - 2)
         endif
         
         if ( ipord .gt. 3 ) then
            nshape = nshape + 3*(ipord - 2)*(ipord - 3)/2
         endif
         
         if ( ipord .gt. 4 ) then
            nshape = nshape + (ipord - 2)*(ipord - 3)*(ipord - 4)/6
         endif
         
         nshapeb = nenb
c     
         if (nenb .eq. 3) then   
            if ( ipord .gt. 1 ) then
               nshapeb = nshapeb + 3*(ipord - 1)
            endif
            if ( ipord .gt. 2 ) then 
               nshape = nshape +(ipord - 1)*(ipord - 2)/2
            endif 
         endif
c     
         if (nenb .eq. 4) then
            if ( ipord .gt. 1 ) then
               nshapeb = nshapeb + 4*(ipord - 1)
            endif
            if ( ipord .gt. 2 ) then 
               nshape = nshape +(ipord - 2)*(ipord - 3)/2
            endif  
         endif
      endif
      
      select case (intg(1,1))
      case (2)
         nint(3) = 6
      case (3)
         nint(3) = 18
      case (4)
         nint(3) = 48
      end select
      select case (intg(2,1))
      case (2)
         nintb(3) = 3
         nintb(4) = 4
      case (3)
         nintb(3) = 6
         nintb(4) = 9
      case (4)
         nintb(3) = 12
         nintb(4) = 16
      end select
      
      allocate (tmpQpt (4,nint(3)))           
      allocate (tmpQwt (nint(3)))           
      
      call symwdg(nint(3),tmpQpt,tmpQwt,nerr) ! interior elements
      Qpt(3,1:4,1:nint(3)) = tmpQpt(1:4,1:nint(3))
      Qwt(3,1:nint(3)) = tmpQwt(1:nint(3))
      
      deallocate (tmpQpt)
      deallocate (tmpQwt)
      
      allocate (tmpQptb(4,nintb(3)))           
      allocate (tmpQwtb(nintb(3))) 
      
      call symtri(nintb(3),tmpQptb,tmpQwtb,nerr) ! boundary elements
      Qptb(3,1:2,2:nintb(3)) = tmpQptb(1:2,1:nintb(3)-1)
      Qptb(3,1:2,1) = tmpQptb(1:2,nintb(3))
c     
c     wedges want the third entry to be zeta=-1 (not t=1-r-s)
c     4th entry not used
c     
      Qptb(3,3:4,1:nintb(3)) = -1
c$$$  Qptb(3,1:4,1:nintb(3)) = tmpQptb(1:4,1:nintb(3))
      Qwtb(3,1:nintb(3)) = tmpQwtb(1:nintb(3))
      
      deallocate (tmpQptb)
      deallocate (tmpQwtb)
      
      allocate (tmpQptb(4,nintb(4)))           
      allocate (tmpQwtb(nintb(4)))
      
      call symquadw(nintb(4),tmpQptb,tmpQwtb,nerr) ! boundary elements
      Qptb(4,1:4,1:nintb(4)) = tmpQptb(1:4,1:nintb(4))
      Qwtb(4,1:nintb(4)) = tmpQwtb(1:nintb(4))
      
      deallocate (tmpQptb)
      deallocate (tmpQwtb)
      
c     
c.... return
c     
      return
      end
