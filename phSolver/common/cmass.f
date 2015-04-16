      module lhsGkeep

      real*8, allocatable :: lhsG(:)

      end module

c----------------------------------------------------------------------------

      subroutine keeplhsG

      use lhsGkeep

      include "common.h"

      allocate ( lhsG(nnz*nshg) )

      return

      end
      subroutine cmass (shp, shgl, xl, em)
c       
c----------------------------------------------------------------------
c     
c     This subroutine computes the consistent mass matrices
c     
c     Ken Jansen, Spring 2000
c----------------------------------------------------------------------
c     
c
      include "common.h"
c     
      integer ne, na, nb, nodlcla, nodlclb, iel
      dimension 
     &     shp(nshl,MAXQPT),   shgl(nsd,nshl,MAXQPT),
     &     em(npro,nshl,nshl),
     &     xl(npro,nenl,nsd)
c
      dimension shape(npro,nshl),   shdrv(npro,nsd,nshl),
     &          sgn(npro,nshl),     dxidx(npro,nsd,nsd),  
     &          shg(npro,nshl,nsd), 
     &          WdetJ(npro)
c
      em = zero
c     
c.... loop through the integration points
c     
      do intp = 1, ngauss      ! (these are in common.h)
c
c.... get the hierarchic shape functions at this int point
c
         call getshp(shp,         shgl,         sgn, 
     &               shape,       shdrv,        intp)
c     
c.... calculate the determinant of the jacobian and weight it
c     
         call e3metric( xl, shdrv,dxidx,shg,WdetJ)
c     
         do iel = 1, npro
            do  na  = 1, nshl
               do  nb = 1, nshl
                  shp2 = shape(iel,na) * shape(iel,nb)
                  em(iel,na,nb) = em(iel,na,nb) + shp2*WdetJ(iel)
               enddo
            enddo
         enddo
      enddo
c     
c.... return
c     
      return
      end

      subroutine cmassl (shp, shgl, shpf, shglf, xl, em, Qwtf)
c       
c----------------------------------------------------------------------
c     
c     This subroutine computes the consistent mass matrices
c     
c     Ken Jansen, Spring 2000
c----------------------------------------------------------------------
c     
c
      include "common.h"
c     
      integer ne, na, nb, nodlcla, nodlclb, iel
      dimension 
     &     shp(nshl,MAXQPT),   shgl(nsd,nshl,MAXQPT),
     &     shpf(nshl,MAXQPT),  shglf(nsd,nshl,MAXQPT),
     &     em(npro,nshl,nshl), eml(npro,nshl),
     &     xl(npro,nenl,nsd)
c
      dimension shape(npro,nshl),   shdrv(npro,nsd,nshl),
     &          sgn(npro,nshl),     dxidx(npro,nsd,nsd),  
     &          shg(npro,nshl,nsd), Qwtf(ngaussf),
     &          WdetJ(npro)
c
      em = zero
      eml= zero

      if (ifproj.eq.1)then
         nods = nshl
      else
         nods = nenl
      endif

c----------------> Get the lumped mass matrix <-----------------------

c     
c.... loop through the integration points
c     
      do intp = 1, ngaussf      ! (these are in common.h)
c
c.... get the hierarchic shape functions at this int point
c
         call getshp(shpf,         shglf,         sgn, 
     &               shape,       shdrv,        intp)
c     
c.... calculate the determinant of the jacobian and weight it
c     
         call e3metricf( xl, shdrv,dxidx,shg,WdetJ,Qwtf)
c     
         do i=1,nods !nenl !nshl
            eml(:,i) = eml(:,i) + shape(:,i)*WdetJ(:)
         enddo         

      enddo ! End loop over quad points.

          
c--------------> Get the consistent mass matrix <------------------------


      shape = zero
      shdrv = zero
      dxidx = zero
      WdetJ = zero
      shg   = zero

c     
c.... loop through the integration points
c     
      do intp = 1, ngauss       ! (these are in common.h)

c.... get the hierarchic shape functions at this int point
c
c
c.... for the mass matrix to be consistent shp and shgl must be
c.... evaluated with at least higher quadrature than one-pt. quad. 

         call getshp(shp,         shgl,         sgn,  
     &               shape,       shdrv,        intp)

c     
c.... calculate the determinant of the jacobian and weight it
c     
         call e3metric( xl, shdrv,dxidx,shg,WdetJ)
c     

         do iel = 1, npro
            do  na  = 1, nods !nenl !nshl
               do  nb = 1, nods !nenl !nshl
                  shp2 = shape(iel,na) * shape(iel,nb)
                  em(iel,na,nb) = em(iel,na,nb) + shp2*WdetJ(iel)
               enddo
            enddo
         enddo

      enddo    ! End loop over quadrature points



c----------> Obtain a mixed (lumped/consistent) mass matrix <------------

c... Different combinations of the lump and mass matrices yield
c... filters of varying widths. In the limiting case were
c... the entire matrix is lumped, we obtain the same filter as 
c... in getdmc.f. Note that in these equivalent ways of
c... filtering one-point quadrature is used for shpf and shgl..


      em = (one-flump)*em

      do iel = 1, npro
         do  na  = 1, nods !nenl !nshl
            em(iel,na,na) = em(iel,na,na)+flump*eml(iel,na)
         enddo
      enddo

c     
c.... return
c     
      return
      end



      subroutine cmasstl (shp, shgl, shpf, shglf, xl, em, Qwtf)
c       
c----------------------------------------------------------------------
c     
c     This subroutine computes the consistent mass matrices
c     
c     Ken Jansen, Spring 2000
c----------------------------------------------------------------------
c     
c
      include "common.h"
c     
      integer ne, na, nb, nodlcla, nodlclb, iel
      dimension 
     &     shp(nshl,MAXQPT),   shgl(nsd,nshl,MAXQPT),
     &     shpf(nshl,MAXQPT),  shglf(nsd,nshl,MAXQPT),
     &     em(npro,nshl,nshl), eml(npro,nshl),
     &     xl(npro,nenl,nsd)
c
      dimension shape(npro,nshl),   shdrv(npro,nsd,nshl),
     &          sgn(npro,nshl),     dxidx(npro,nsd,nsd),  
     &          shg(npro,nshl,nsd), Qwtf(ngaussf),
     &          WdetJ(npro)
c
      em = zero
      eml= zero

c----------------> Get the lumped mass matrix <-----------------------

c     
c.... loop through the integration points
c     
      do intp = 1, ngaussf      ! (these are in common.h)
c
c.... get the hierarchic shape functions at this int point
c
         call getshp(shpf,         shglf,         sgn, 
     &               shape,       shdrv,        intp)
c     
c.... calculate the determinant of the jacobian and weight it
c     
         call e3metricf( xl, shdrv,dxidx,shg,WdetJ,Qwtf)
c     
         do i=1,nshl
            eml(:,i) = eml(:,i) + shape(:,i)*WdetJ(:)
         enddo         

      enddo ! End loop over quad points.

          
c--------------> Get the consistent mass matrix <------------------------


      shape = zero
      shdrv = zero
      dxidx = zero
      WdetJ = zero
      shg   = zero

c     
c.... loop through the integration points
c     
      do intp = 1, ngauss       ! (these are in common.h)

c.... get the hierarchic shape functions at this int point
c
c
c.... for the mass matrix to be consistent shp and shgl must be
c.... evaluated with at least higher quadrature than one-pt. quad. 

         call getshp(shp,         shgl,         sgn,  
     &               shape,       shdrv,        intp)

c     
c.... calculate the determinant of the jacobian and weight it
c     
         call e3metric( xl, shdrv,dxidx,shg,WdetJ)
c     

         do iel = 1, npro
            do  na  = 1, nshl
               do  nb = 1, nshl
                  shp2 = shape(iel,na) * shape(iel,nb)
                  em(iel,na,nb) = em(iel,na,nb) + shp2*WdetJ(iel)
               enddo
            enddo
         enddo

      enddo    ! End loop over quadrature points



c----------> Obtain a mixed (lumped/consistent) mass matrix <------------

c... Different combinations of the lump and mass matrices yield
c... filters of varying widths. In the limiting case were
c... the entire matrix is lumped, we obtain the same filter as 
c... in getdmc.f. Note that in these equivalent ways of
c... filtering one-point quadrature is used for shpf and shgl..


      do iel = 1, npro
         do  na  = 1, nshl
            em(iel,na,na) = eml(iel,na)
         enddo
      enddo

c     
c.... return
c     
      return
      end

c-----------------------------------------------------------------------
c
c  compute the metrics of the mapping from global to local 
c  coordinates and the jacobian of the mapping (weighted by 
c  the quadrature weight
c
c-----------------------------------------------------------------------
      subroutine e3metricf(  xl,      shgl,     dxidx,
     &                      shg,     WdetJ, Qwtf)

      include "common.h"
      
      real*8     xl(npro,nenl,nsd),    shgl(npro,nsd,nshl),
     &           dxidx(npro,nsd,nsd),  shg(npro,nshl,nsd), 
     &           WdetJ(npro),          Qwtf(ngaussf)

      real*8     dxdxi(npro,nsd,nsd),  tmp(npro)

c
c.... compute the deformation gradient
c
      dxdxi = zero
c
       do n = 1, nenl
          dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(:,1,n)
          dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(:,2,n)
          dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(:,3,n)
          dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(:,1,n)
          dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(:,2,n)
          dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(:,3,n)
          dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(:,1,n)
          dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(:,2,n)
          dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(:,3,n)
       enddo
c
c.... compute the inverse of deformation gradient
c
       dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3) 
     &                - dxdxi(:,3,2) * dxdxi(:,2,3)
       dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3) 
     &                - dxdxi(:,1,2) * dxdxi(:,3,3)
       dxidx(:,1,3) =  dxdxi(:,1,2) * dxdxi(:,2,3) 
     &                - dxdxi(:,1,3) * dxdxi(:,2,2)
       tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1) 
     &                       + dxidx(:,1,2) * dxdxi(:,2,1)  
     &                       + dxidx(:,1,3) * dxdxi(:,3,1) )
       dxidx(:,1,1) = dxidx(:,1,1) * tmp
       dxidx(:,1,2) = dxidx(:,1,2) * tmp
       dxidx(:,1,3) = dxidx(:,1,3) * tmp
       dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1) 
     &                - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
       dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3) 
     &                - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
       dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3) 
     &                - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
       dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2) 
     &                - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
       dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2) 
     &                - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
       dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2) 
     &                - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
c       WdetJ = Qwt(lcsyst,intp) / tmp

       WdetJ = Qwtf(intp) / tmp
c
c.... compute the global gradient of shape-functions
c
       do n = 1, nshl
          shg(:,n,1) = shgl(:,1,n) * dxidx(:,1,1) + 
     &                 shgl(:,2,n) * dxidx(:,2,1) +
     &                 shgl(:,3,n) * dxidx(:,3,1)
          shg(:,n,2) = shgl(:,1,n) * dxidx(:,1,2) + 
     &                 shgl(:,2,n) * dxidx(:,2,2) +
     &                 shgl(:,3,n) * dxidx(:,3,2) 
          shg(:,n,3) = shgl(:,1,n) * dxidx(:,1,3) + 
     &                 shgl(:,2,n) * dxidx(:,2,3) +
     &                 shgl(:,3,n) * dxidx(:,3,3) 
       enddo

       return
       end



