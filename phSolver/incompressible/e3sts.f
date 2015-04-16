      subroutine e3StsLhs( xl,  lStsVec )
c-----------------------------------------------------------------------
c 
c  compute the terms needed for the left hand side matrices 
c  needed for the conservative projection       
c
c-----------------------------------------------------------------------
      use     stats
      include "common.h"
      
      integer i
      real*8  lDir(npro,nshl,3), lStsVec(npro,nshl,nResDims),
     &        xl(npro,nenl,3)

      call e3StsDir( xl,  lDir )
      
      do i = 1, nshl
         lStsVec(:,i,1) = lDir(:,i,1) * lDir(:,i,1)
         lStsVec(:,i,2) = lDir(:,i,2) * lDir(:,i,2)
         lStsVec(:,i,3) = lDir(:,i,3) * lDir(:,i,3)

         lStsVec(:,i,4) = lDir(:,i,1) * lDir(:,i,2)
         lStsVec(:,i,5) = lDir(:,i,2) * lDir(:,i,3)
         lStsVec(:,i,6) = lDir(:,i,3) * lDir(:,i,1)

         lStsVec(:,i,7) = 0.0
         lStsVec(:,i,8) = 0.0
         lStsVec(:,i,9) = 0.0
         lStsVec(:,i,10) = 0.0
         lStsVec(:,i,11) = 0.0
      enddo
      
      return
      end


      subroutine e3StsRes( xl, rl, lStsVec )
c-----------------------------------------------------------------------
c  
c  compute the residual terms for the consistent projection
c
c-----------------------------------------------------------------------      
      use     stats
      include "common.h"
      
      real*8  xl(npro,nenl,3),  rl(npro,nshl,ndof)
      real*8  lDir(npro,nshl,3), lStsVec(npro,nshl,nResDims)
      
      call e3StsDir( xl,  lDir )
      
      do i = 1, nshl
         lStsVec(:,i,1) = lDir(:,i,1) * rl(:,i,4)
         lStsVec(:,i,2) = lDir(:,i,2) * rl(:,i,4)
         lStsVec(:,i,3) = lDir(:,i,3) * rl(:,i,4)
      
         lStsVec(:,i,4) = lDir(:,i,1) * rl(:,i,1)
         lStsVec(:,i,5) = lDir(:,i,2) * rl(:,i,2)
         lStsVec(:,i,6) = lDir(:,i,3) * rl(:,i,3)
      
         lStsVec(:,i,7) = lDir(:,i,1) * rl(:,i,2)
     &                  + lDir(:,i,2) * rl(:,i,1)
         lStsVec(:,i,8) = lDir(:,i,2) * rl(:,i,3)
     &                  + lDir(:,i,3) * rl(:,i,2)
         lStsVec(:,i,9) = lDir(:,i,3) * rl(:,i,1)
     &        + lDir(:,i,1) * rl(:,i,3)
         lStsVec(:,i,10) = 0
         lStsVec(:,i,11) = 0
      enddo
      

      return
      end

      subroutine e3StsDir( xl,  lDir )
c-----------------------------------------------------------------------
c
c  compute the normal to each of the nodes
c
c-----------------------------------------------------------------------
      include "common.h"
      
      real*8  xl(npro,nenl,3), lDir(npro,nshl,3)
      integer e

c
c.... linear tets
c
      if (nshl .eq. 4 ) then
         fct = 1.d0 / 6.d0
         do e = 1, npro

            x12         = xl(e,2,1) - xl(e,1,1)
            x13         = xl(e,3,1) - xl(e,1,1)
            x14         = xl(e,4,1) - xl(e,1,1)
            x23         = xl(e,3,1) - xl(e,2,1)
            x24         = xl(e,4,1) - xl(e,2,1)
            x34         = xl(e,4,1) - xl(e,3,1)

            y12         = xl(e,2,2) - xl(e,1,2)
            y13         = xl(e,3,2) - xl(e,1,2)
            y14         = xl(e,4,2) - xl(e,1,2)
            y23         = xl(e,3,2) - xl(e,2,2)
            y24         = xl(e,4,2) - xl(e,2,2)
            y34         = xl(e,4,2) - xl(e,3,2)

            z12         = xl(e,2,3) - xl(e,1,3)
            z13         = xl(e,3,3) - xl(e,1,3)
            z14         = xl(e,4,3) - xl(e,1,3)
            z23         = xl(e,3,3) - xl(e,2,3)
            z24         = xl(e,4,3) - xl(e,2,3)
            z34         = xl(e,4,3) - xl(e,3,3)
            
c
c.. The calculation of the direction of a vertex is based on the average of
c.. the normals of the neighbor faces(3); And the calculation of the direction
c.. of the edge is based on the neighbor faces(2);
c
	    lDir(e,1,1) = fct * (y14*(z12 - z13) + y12*(z13 - z14) 
     &                         + y13*(-z12 + z14))
	    lDir(e,1,2) = fct * ( x14*(-z12 + z13) + x13*(z12 - z14) 
     &                          + x12*(-z13 + z14))
	    lDir(e,1,3) = fct * ( x14*(y12 - y13)
     &                          + x12*(y13 - y14) + x13*(-y12 + y14))
c
	    lDir(e,2,1) = fct * (-(y13*z12) + y14*z12 + y12*z13
     &                          - y12*z14 + y24*z23 - y23*z24)
            lDir(e,2,2) = fct * (x13*z12 - x14*z12-x12*z13 + x12*z14
     &                          - x24*z23 + x23*z24 )
            lDir(e,2,3) = fct * (-(x13*y12) + x14*y12 + x12*y13
     &                         - x12*y14 + x24*y23 - x23*y24)
c
            lDir(e,3,1) = fct * (y12*z13 - y14*z13 + y13*(-z12 + z14)
     &                         + y24*z23 - y23*z24)
            lDir(e,3,2) = fct * (-(x12*z13) + x14*z13 + x13*(z12 - z14)
     &                         - x24*z23 + x23*z24)
            lDir(e,3,3) = fct * (x12*y13 - x14*y13 + x13*(-y12 + y14)
     &                         + x24*y23 - x23*y24)
c
            lDir(e,4,1) = fct * (y14*(z12 - z13) - y12*z14 + y13*z14
     &                         + y24*z23 - y23*z24)
            lDir(e,4,2) = fct * (x14*(-z12 + z13) + x12*z14 - x13*z14
     &                         - x24*z23 + x23*z24)
            lDir(e,4,3) = fct * (x14*(y12 - y13) - x12*y14 + x13*y14
     &                         + x24*y23 - x23*y24)
c
         enddo
c
c.... quadratic tets
c
      else if (nshl .eq. 10 ) then
         fct = 1.d0 / 6.d0
         do e = 1, npro

            x12         = xl(e,2,1) - xl(e,1,1)
            x13         = xl(e,3,1) - xl(e,1,1)
            x14         = xl(e,4,1) - xl(e,1,1)
            x23         = xl(e,3,1) - xl(e,2,1)
            x24         = xl(e,4,1) - xl(e,2,1)
            x34         = xl(e,4,1) - xl(e,3,1)

            y12         = xl(e,2,2) - xl(e,1,2)
            y13         = xl(e,3,2) - xl(e,1,2)
            y14         = xl(e,4,2) - xl(e,1,2)
            y23         = xl(e,3,2) - xl(e,2,2)
            y24         = xl(e,4,2) - xl(e,2,2)
            y34         = xl(e,4,2) - xl(e,3,2)

            z12         = xl(e,2,3) - xl(e,1,3)
            z13         = xl(e,3,3) - xl(e,1,3)
            z14         = xl(e,4,3) - xl(e,1,3)
            z23         = xl(e,3,3) - xl(e,2,3)
            z24         = xl(e,4,3) - xl(e,2,3)
            z34         = xl(e,4,3) - xl(e,3,3)
            
c
c.... vertex modes
c            
	    lDir(e,1,1) = fct * (y14*(z12 - z13) + y12*(z13 - z14) 
     &                         + y13*(-z12 + z14))
	    lDir(e,1,2) = fct * ( x14*(-z12 + z13) + x13*(z12 - z14) 
     &                          + x12*(-z13 + z14))
	    lDir(e,1,3) = fct * ( x14*(y12 - y13)
     &                          + x12*(y13 - y14) + x13*(-y12 + y14))
c
	    lDir(e,2,1) = fct * (-(y13*z12) + y14*z12 + y12*z13
     &                          - y12*z14 + y24*z23 - y23*z24)
            lDir(e,2,2) = fct * (x13*z12 - x14*z12-x12*z13 + x12*z14
     &                          - x24*z23 + x23*z24 )
            lDir(e,2,3) = fct * (-(x13*y12) + x14*y12 + x12*y13
     &                         - x12*y14 + x24*y23 - x23*y24)
c
            lDir(e,3,1) = fct * (y12*z13 - y14*z13 + y13*(-z12 + z14)
     &                         + y24*z23 - y23*z24)
            lDir(e,3,2) = fct * (-(x12*z13) + x14*z13 + x13*(z12 - z14)
     &                         - x24*z23 + x23*z24)
            lDir(e,3,3) = fct * (x12*y13 - x14*y13 + x13*(-y12 + y14)
     &                         + x24*y23 - x23*y24)
c
            lDir(e,4,1) = fct * (y14*(z12 - z13) - y12*z14 + y13*z14
     &                         + y24*z23 - y23*z24)
            lDir(e,4,2) = fct * (x14*(-z12 + z13) + x12*z14 - x13*z14
     &                         - x24*z23 + x23*z24)
            lDir(e,4,3) = fct * (x14*(y12 - y13) - x12*y14 + x13*y14
     &                         + x24*y23 - x23*y24)
c
c.... edge modes (quadratic)
c
            lDir(e,5,1) = pt25*(-(y13*z12) + y14*z12 + y12*(z13-z14))
            lDir(e,5,2) = pt25*(x13*z12 - x14*z12 + x12*(-z13 + z14))
            lDir(e,5,3) = pt25*(-(x13*y12) + x14*y12 + x12*(y13 - y14))

            lDir(e,6,1) = pt25*(-(y13*z12) + y12*z13 + y24*z23-y23*z24)
            lDir(e,6,2) = pt25*(x13*z12 - x12*z13 - x24*z23 + x23*z24)
            lDir(e,6,3) = pt25*(-(x13*y12) + x12*y13 + x24*y23-x23*y24)

            lDir(e,7,1) = pt25*((y12 - y14)*z13 + y13*(-z12 + z14))
            lDir(e,7,2) = pt25*((-x12 + x14)*z13 + x13*(z12 - z14))
            lDir(e,7,3) = pt25*((x12 - x14)*y13 + x13*(-y12 + y14))

            lDir(e,8,1) = pt25*(y14*(z12 - z13) + (-y12 + y13)*z14)
            lDir(e,8,2) = pt25*(x14*(-z12 + z13) + (x12 - x13)*z14)
            lDir(e,8,3) = pt25*(x14*(y12 - y13) + (-x12 + x13)*y14)

            lDir(e,9,1) = pt25*(y14*z12 - y12*z14 + y24*z23 - y23*z24)
            lDir(e,9,2) = pt25*(-(x14*z12) + x12*z14 - x24*z23+x23*z24)
            lDir(e,9,3) = pt25*(x14*y12 - x12*y14 + x24*y23 - x23*y24)

            lDir(e,10,1) = pt25*(-(y14*z13) + y13*z14+y24*z23-y23*z24)
            lDir(e,10,2) = pt25*(x14*z13 - x13*z14-x24*z23 + x23*z24)
            lDir(e,10,3) = pt25*(-(x14*y13) + x13*y14+x24*y23-x23*y24)


         enddo
c
c.... cubic tets
c
      else if (nshl .eq. 20 ) then
         fct = 1.d0 / 6.d0
         do e = 1, npro

            x12         = xl(e,2,1) - xl(e,1,1)
            x13         = xl(e,3,1) - xl(e,1,1)
            x14         = xl(e,4,1) - xl(e,1,1)
            x23         = xl(e,3,1) - xl(e,2,1)
            x24         = xl(e,4,1) - xl(e,2,1)
c$$$            x34         = xl(e,4,1) - xl(e,3,1)

            y12         = xl(e,2,2) - xl(e,1,2)
            y13         = xl(e,3,2) - xl(e,1,2)
            y14         = xl(e,4,2) - xl(e,1,2)
            y23         = xl(e,3,2) - xl(e,2,2)
            y24         = xl(e,4,2) - xl(e,2,2)
c$$$            y34         = xl(e,4,2) - xl(e,3,2)

            z12         = xl(e,2,3) - xl(e,1,3)
            z13         = xl(e,3,3) - xl(e,1,3)
            z14         = xl(e,4,3) - xl(e,1,3)
            z23         = xl(e,3,3) - xl(e,2,3)
            z24         = xl(e,4,3) - xl(e,2,3)
c$$$            z34         = xl(e,4,3) - xl(e,3,3)
            
c
c.... vertex modes
c            
	    lDir(e,1,1) = fct * (y14*(z12 - z13) + y12*(z13 - z14) 
     &                         + y13*(-z12 + z14))
	    lDir(e,1,2) = fct * ( x14*(-z12 + z13) + x13*(z12 - z14) 
     &                          + x12*(-z13 + z14))
	    lDir(e,1,3) = fct * ( x14*(y12 - y13)
     &                          + x12*(y13 - y14) + x13*(-y12 + y14))
c
	    lDir(e,2,1) = fct * (-(y13*z12) + y14*z12 + y12*z13
     &                          - y12*z14 + y24*z23 - y23*z24)
            lDir(e,2,2) = fct * (x13*z12 - x14*z12-x12*z13 + x12*z14
     &                          - x24*z23 + x23*z24 )
            lDir(e,2,3) = fct * (-(x13*y12) + x14*y12 + x12*y13
     &                         - x12*y14 + x24*y23 - x23*y24)
c
            lDir(e,3,1) = fct * (y12*z13 - y14*z13 + y13*(-z12 + z14)
     &                         + y24*z23 - y23*z24)
            lDir(e,3,2) = fct * (-(x12*z13) + x14*z13 + x13*(z12 - z14)
     &                         - x24*z23 + x23*z24)
            lDir(e,3,3) = fct * (x12*y13 - x14*y13 + x13*(-y12 + y14)
     &                         + x24*y23 - x23*y24)
c
            lDir(e,4,1) = fct * (y14*(z12 - z13) - y12*z14 + y13*z14
     &                         + y24*z23 - y23*z24)
            lDir(e,4,2) = fct * (x14*(-z12 + z13) + x12*z14 - x13*z14
     &                         - x24*z23 + x23*z24)
            lDir(e,4,3) = fct * (x14*(y12 - y13) - x12*y14 + x13*y14
     &                         + x24*y23 - x23*y24)
c
c.... edge modes (quadratic and cubic)
c
            lDir(e,5,1) = pt25*(-(y13*z12) + y14*z12 + y12*(z13-z14))
            lDir(e,5,2) = pt25*(x13*z12 - x14*z12 + x12*(-z13 + z14))
            lDir(e,5,3) = pt25*(-(x13*y12) + x14*y12 + x12*(y13 - y14))
            lDir(e,6,1) = lDir(e,5,1)
            lDir(e,6,2) = lDir(e,5,2)
            lDir(e,6,3) = lDir(e,5,3)

            lDir(e,7,1) = pt25*(-(y13*z12) + y12*z13 + y24*z23-y23*z24)
            lDir(e,7,2) = pt25*(x13*z12 - x12*z13 - x24*z23 + x23*z24)
            lDir(e,7,3) = pt25*(-(x13*y12) + x12*y13 + x24*y23-x23*y24)
            lDir(e,8,1) = lDir(e,7,1)
            lDir(e,8,2) = lDir(e,7,2)
            lDir(e,8,3) = lDir(e,7,3)

            lDir(e,9,1) = pt25*((y12 - y14)*z13 + y13*(-z12 + z14))
            lDir(e,9,2) = pt25*((-x12 + x14)*z13 + x13*(z12 - z14))
            lDir(e,9,3) = pt25*((x12 - x14)*y13 + x13*(-y12 + y14))
            lDir(e,10,1) = lDir(e,9,1)
            lDir(e,10,2) = lDir(e,9,2)
            lDir(e,10,3) = lDir(e,9,3)

            lDir(e,11,1) = pt25*(y14*(z12 - z13) + (-y12 + y13)*z14)
            lDir(e,11,2) = pt25*(x14*(-z12 + z13) + (x12 - x13)*z14)
            lDir(e,11,3) = pt25*(x14*(y12 - y13) + (-x12 + x13)*y14)
            lDir(e,12,1) = lDir(e,11,1)
            lDir(e,12,2) = lDir(e,11,2)
            lDir(e,12,3) = lDir(e,11,3)


            lDir(e,13,1) = pt25*(y14*z12 - y12*z14 + y24*z23 - y23*z24)
            lDir(e,13,2) = pt25*(-(x14*z12) + x12*z14-x24*z23+x23*z24)
            lDir(e,13,3) = pt25*(x14*y12 - x12*y14 + x24*y23 - x23*y24)
            lDir(e,14,1) = lDir(e,13,1)
            lDir(e,14,2) = lDir(e,13,2)
            lDir(e,14,3) = lDir(e,13,3)

            lDir(e,15,1) = pt25*(-(y14*z13) + y13*z14+y24*z23-y23*z24)
            lDir(e,15,2) = pt25*(x14*z13 - x13*z14-x24*z23 + x23*z24)
            lDir(e,15,3) = pt25*(-(x14*y13) + x13*y14+x24*y23-x23*y24)
            lDir(e,16,1) = lDir(e,15,1)
            lDir(e,16,2) = lDir(e,15,2)
            lDir(e,16,3) = lDir(e,15,3)
c
c.... face modes (cubic)
c
            lDir(e,17,1) = pt5*(-(y13*z12) + y12*z13)
            lDir(e,17,2) = pt5*(x13*z12 - x12*z13)
            lDir(e,17,3) = pt5*(-(x13*y12) + x12*y13)

            lDir(e,18,1) = pt5*(y14*z12 - y12*z14)
            lDir(e,18,2) = pt5*(-(x14*z12) + x12*z14)
            lDir(e,18,3) = pt5*(x14*y12 - x12*y14)

            lDir(e,19,1) = pt5*(y24*z23 - y23*z24)
            lDir(e,19,2) = pt5*(-(x24*z23) + x23*z24)
            lDir(e,19,3) = pt5*(x24*y23 - x23*y24)

            lDir(e,20,1) = pt5*(-(y14*z13) + y13*z14)
            lDir(e,20,2) = pt5*(x14*z13 - x13*z14)
            lDir(e,20,3) = pt5*(-(x14*y13) + x13*y14)

         enddo
c
c.... hexes
c     
      else if (nenl .eq. 8) then
         fct = 1.d0 / 12.d0
         do e = 1, npro
	    x13		= xl(e,1,1) - xl(e,3,1)
	    x16		= xl(e,1,1) - xl(e,6,1)
	    x18		= xl(e,1,1) - xl(e,8,1)
	    x24		= xl(e,2,1) - xl(e,4,1)
	    x25		= xl(e,2,1) - xl(e,5,1)
	    x27		= xl(e,2,1) - xl(e,7,1)
	    x36		= xl(e,3,1) - xl(e,6,1)
	    x38		= xl(e,3,1) - xl(e,8,1)
	    x45		= xl(e,4,1) - xl(e,5,1)
	    x47		= xl(e,4,1) - xl(e,7,1)
	    x57		= xl(e,5,1) - xl(e,7,1)
	    x68		= xl(e,6,1) - xl(e,8,1)
c
	    y13		= xl(e,1,2) - xl(e,3,2)
	    y16		= xl(e,1,2) - xl(e,6,2)
	    y18		= xl(e,1,2) - xl(e,8,2)
	    y24		= xl(e,2,2) - xl(e,4,2)
	    y25		= xl(e,2,2) - xl(e,5,2)
	    y27		= xl(e,2,2) - xl(e,7,2)
	    y36		= xl(e,3,2) - xl(e,6,2)
	    y38		= xl(e,3,2) - xl(e,8,2)
	    y45		= xl(e,4,2) - xl(e,5,2)
	    y47		= xl(e,4,2) - xl(e,7,2)
	    y57		= xl(e,5,2) - xl(e,7,2)
	    y68		= xl(e,6,2) - xl(e,8,2)
c
	    z13		= xl(e,1,3) - xl(e,3,3)
	    z16		= xl(e,1,3) - xl(e,6,3)
	    z18		= xl(e,1,3) - xl(e,8,3)
	    z24		= xl(e,2,3) - xl(e,4,3)
	    z25		= xl(e,2,3) - xl(e,5,3)
	    z27		= xl(e,2,3) - xl(e,7,3)
	    z36		= xl(e,3,3) - xl(e,6,3)
	    z38		= xl(e,3,3) - xl(e,8,3)
	    z45		= xl(e,4,3) - xl(e,5,3)
	    z47		= xl(e,4,3) - xl(e,7,3)
	    z57		= xl(e,5,3) - xl(e,7,3)
	    z68		= xl(e,6,3) - xl(e,8,3)
c
	x31= -x13
	x61= -x16
	x81= -x18
	x42= -x24
	x52= -x25
	x72= -x27
	x63= -x36
	x83= -x38
	x54= -x45
	x74= -x47
	x75= -x57
	x86= -x68
	y31= -y13
	y61= -y16
	y81= -y18
	y42= -y24
	y52= -y25
	y72= -y27
	y63= -y36
	y83= -y38
	y54= -y45
	y74= -y47
	y75= -y57
	y86= -y68
	z31= -z13
	z61= -z16
	z81= -z18
	z42= -z24
	z52= -z25
	z72= -z27
	z63= -z36
	z83= -z38
	z54= -z45
	z74= -z47
	z75= -z57
	z86= -z68

	    lDir(e,1,1) = fct * (-y24 * z45 + y36 * z24 - y68 * z45 
     1				+ z24 * y45 - z36 * y24 + z68 * y45 )
	    lDir(e,2,1) = fct * (-y16 * z63 + y54 * z16 - y47 * z63 
     1				+ z16 * y63 - z54 * y16 + z47 * y63 )
	    lDir(e,3,1) = fct * (-y42 * z27 + y18 * z42 - y86 * z27 
     1				+ z42 * y27 - z18 * y42 + z86 * y27 )
	    lDir(e,4,1) = fct * (-y38 * z81 + y72 * z38 - y25 * z81 
     1				+ z38 * y81 - z72 * y38 + z25 * y81 )
	    lDir(e,5,1) = fct * (-y61 * z18 + y27 * z61 - y74 * z18 
     1				+ z61 * y18 - z27 * y61 + z74 * y18 )
	    lDir(e,6,1) = fct * (-y57 * z72 + y81 * z57 - y13 * z72 
     1				+ z57 * y72 - z81 * y57 + z13 * y72 )
	    lDir(e,7,1) = fct * (-y83 * z36 + y45 * z83 - y52 * z36 
     1				+ z83 * y36 - z45 * y83 + z52 * y36 )
	    lDir(e,8,1) = fct * (-y75 * z54 + y63 * z75 - y31 * z54 
     1				+ z75 * y54 - z63 * y75 + z31 * y54 )
c
	    lDir(e,1,2) = fct * (-z24 * x45 + z36 * x24 - z68 * x45 
     1				+ x24 * z45 - x36 * z24 + x68 * z45 )
	    lDir(e,2,2) = fct * (-z16 * x63 + z54 * x16 - z47 * x63 
     1				+ x16 * z63 - x54 * z16 + x47 * z63 )
	    lDir(e,3,2) = fct * (-z42 * x27 + z18 * x42 - z86 * x27 
     1				+ x42 * z27 - x18 * z42 + x86 * z27 )
	    lDir(e,4,2) = fct * (-z38 * x81 + z72 * x38 - z25 * x81 
     1				+ x38 * z81 - x72 * z38 + x25 * z81 )
	    lDir(e,5,2) = fct * (-z61 * x18 + z27 * x61 - z74 * x18 
     1				+ x61 * z18 - x27 * z61 + x74 * z18 )
	    lDir(e,6,2) = fct * (-z57 * x72 + z81 * x57 - z13 * x72 
     1				+ x57 * z72 - x81 * z57 + x13 * z72 )
	    lDir(e,7,2) = fct * (-z83 * x36 + z45 * x83 - z52 * x36 
     1				+ x83 * z36 - x45 * z83 + x52 * z36 )
	    lDir(e,8,2) = fct * (-z75 * x54 + z63 * x75 - z31 * x54 
     1				+ x75 * z54 - x63 * z75 + x31 * z54 )
c
	    lDir(e,1,3) = fct * (-x24 * y45 + x36 * y24 - x68 * y45 
     1				+ y24 * x45 - y36 * x24 + y68 * x45 )
	    lDir(e,2,3) = fct * (-x16 * y63 + x54 * y16 - x47 * y63 
     1				+ y16 * x63 - y54 * x16 + y47 * x63 )
	    lDir(e,3,3) = fct * (-x42 * y27 + x18 * y42 - x86 * y27 
     1				+ y42 * x27 - y18 * x42 + y86 * x27 )
	    lDir(e,4,3) = fct * (-x38 * y81 + x72 * y38 - x25 * y81 
     1				+ y38 * x81 - y72 * x38 + y25 * x81 )
	    lDir(e,5,3) = fct * (-x61 * y18 + x27 * y61 - x74 * y18 
     1				+ y61 * x18 - y27 * x61 + y74 * x18 )
	    lDir(e,6,3) = fct * (-x57 * y72 + x81 * y57 - x13 * y72 
     1				+ y57 * x72 - y81 * x57 + y13 * x72 )
	    lDir(e,7,3) = fct * (-x83 * y36 + x45 * y83 - x52 * y36 
     1				+ y83 * x36 - y45 * x83 + y52 * x36 )
	    lDir(e,8,3) = fct * (-x75 * y54 + x63 * y75 - x31 * y54 
     1				+ y75 * x54 - y63 * x75 + y31 * x54 )
c
         enddo
      else
         write(*,*) 'Error in e3sts: elt type not impl.'
         stop
      endif
      
      return
      end
         
      
         
