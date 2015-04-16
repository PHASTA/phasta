
c---------------------------------------------------------------------------

      subroutine SUPGstress (y, ac, x, qres, ien, xmudmi,
     &     cdelsq, shgl, shp, Qwtf, shglo, shpo, stress, diss, vol)

      use stats
      use rlssave   ! Use the resolved Leonard stresses at the nodes.

      include "common.h"

      dimension y(nshg,5),                  ac(nshg,5),
     &          x(numnp,nsd),               ien(npro,nshl),
     &          shp(nshl,ngauss),            shpfun(npro,nshl),
     &          shgl(nsd,nshl,ngauss),       shg(npro,nshl,nsd),
     &          shglo(nsd,nshl,ngauss),      shpo(nshl,ngauss),
     &          Qwtf(ngaussf),              acl(npro,nshl,ndof), 
     &          yl(npro,nshl,ndof),         xl(npro,nenl,nsd)
      dimension stress(nshg,9),             stressl(npro,9),
     &          stressli(npro,9),
     &          dxdxi(npro,nsd,nsd),        dxidx(npro,nsd,nsd),
     &          WdetJ(npro),                rho(npro),
     &          tmp(npro),                  aci(npro,nsd),
     &          pres(npro),                 u1(npro),
     &          u2(npro),                   u3(npro)
      dimension qres(nshg,nsd*nsd),         ql(npro,nshl,nsd*nsd),
     &          g1yi(npro,ndof),            g2yi(npro,ndof),
     &          g3yi(npro,ndof),            divqi(npro,3),
     &          src(npro,nsd),             Temp(npro),
     &          xx(npro,nsd),
     &          rlsl(npro,nshl,6),         rlsli(npro,6),
     &          rLui(npro,3),
     &          tauC(npro),                tauM(npro),
     &          tauBar(npro),              uBar(npro,nsd)
      dimension Sij(npro,6),               Snorm(npro),
     &          Snorm2(npro),              cdelsq(nshg),              
     &          xmudmi(npro,ngauss),        xmudmif(npro,ngauss),
     &          dissi(npro,3),             dissl(npro,3),
     &          voli(npro),                voll(npro),
     &          vol(nshg),                 diss(nshg,3),
     &          rmu(npro)

      real*8    omega(3), divu(npro)

c.... Note that the xmudmi passed in here is 
c.... evaluated at quadrature points of the flow. xmudmif will
c...  be evaluated at the ngaussf quad. pts.

      xmudmif = zero

c.... Debuggin

c      xmudmi = zero

c.... Localization
      
      call localy(y,      yl,     ien,    ndofl,  'gather  ')
      call localy(ac,    acl,     ien,    ndofl,  'gather  ')
      call localx(x,      xl,     ien,    nsd,    'gather  ')

      if (idiff==1 .or. idiff==3) then ! global reconstruction of qdiff
         call local (qres,   ql,     ien, nsd*nsd, 'gather  ')
      endif

      if ( idiff==2 .and. ires .eq. 1 ) then
         call e3ql (yl,        shpo,       shglo, 
     &              xl,        ql,        xmudmi, 
     &              sgn)
      endif

      if( (iLES.gt.10).and.(iLES.lt.20)) then ! bardina 
         call local (rls, rlsl,     ien,       6, 'gather  ')  
      else
         rlsl = zero
      endif      

c... Now that everything is localized, begin loop over ngaussf quad. pts.


      stressl = zero
      dissl   = zero
      voll    = zero

      do intp = 1, ngaussf

c
c.... ------------->  Primitive variables at int. point  <--------------
c
       pres = zero
       u1   = zero
       u2   = zero
       u3   = zero
c
       do n = 1, nshl 
          pres(:) = pres(:) + shp(n,intp) * yl(:,n,1)
          u1(:)   = u1(:)   + shp(n,intp) * yl(:,n,2)
          u2(:)   = u2(:)   + shp(n,intp) * yl(:,n,3)
          u3(:)   = u3(:)   + shp(n,intp) * yl(:,n,4)
       enddo

c
c.... ----------------------->  accel. at int. point  <----------------------
c
       aci = zero
       do n = 1, nshl
          aci(:,1) = aci(:,1) + shp(n,intp) * acl(:,n,2)
          aci(:,2) = aci(:,2) + shp(n,intp) * acl(:,n,3)
          aci(:,3) = aci(:,3) + shp(n,intp) * acl(:,n,4)
       enddo


c
c.... --------------->  Element Metrics at int. point <-------------
c
c.... compute the deformation gradient
c
        dxdxi = zero
c
          do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(1,n,intp)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(2,n,intp)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(3,n,intp)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(1,n,intp)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(2,n,intp)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(3,n,intp)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(1,n,intp)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(2,n,intp)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(3,n,intp)
          enddo
c
c.... compute the inverse of deformation gradient
c
        dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)
     &                 - dxdxi(:,3,2) * dxdxi(:,2,3)
        dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)
     &                 - dxdxi(:,1,2) * dxdxi(:,3,3)
        dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)
     &                 - dxdxi(:,1,3) * dxdxi(:,2,2)
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

        wght=Qwtf(intp)
        WdetJ = wght / tmp

c    Obtain the global gradient of the shape functions at current qpt.

      do n = 1,nshl
        shg(:,n,1) = (shgl(1,n,intp) * dxidx(:,1,1)
     &              + shgl(2,n,intp) * dxidx(:,2,1)
     &              + shgl(3,n,intp) * dxidx(:,3,1))
        shg(:,n,2) = (shgl(1,n,intp) * dxidx(:,1,2)
     &              + shgl(2,n,intp) * dxidx(:,2,2)
     &              + shgl(3,n,intp) * dxidx(:,3,2))
        shg(:,n,3) = (shgl(1,n,intp) * dxidx(:,1,3)
     &              + shgl(2,n,intp) * dxidx(:,2,3)
     &              + shgl(3,n,intp) * dxidx(:,3,3))
      enddo

c
c.... compute the global gradient of u and P
c
c
       g1yi = zero
       g2yi = zero
       g3yi = zero
       do n = 1, nshl
          g1yi(:,1) = g1yi(:,1) + shg(:,n,1) * yl(:,n,1)
          g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * yl(:,n,2)
          g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * yl(:,n,3)
          g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * yl(:,n,4)
c
          g2yi(:,1) = g2yi(:,1) + shg(:,n,2) * yl(:,n,1)
          g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * yl(:,n,2)
          g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * yl(:,n,3)
          g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * yl(:,n,4)
c
          g3yi(:,1) = g3yi(:,1) + shg(:,n,3) * yl(:,n,1)
          g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * yl(:,n,2)
          g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * yl(:,n,3)
          g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * yl(:,n,4)        
       enddo

c.... Let us build the Sij tensor and its norms

       Sij(:,1) = g1yi(:,2)
       Sij(:,2) = g2yi(:,3)
       Sij(:,3) = g3yi(:,4)
       Sij(:,4) = (g2yi(:,2)+g1yi(:,3))*pt5
       Sij(:,5) = (g3yi(:,2)+g1yi(:,4))*pt5
       Sij(:,6) = (g3yi(:,3)+g2yi(:,4))*pt5

       Snorm(:) = Sij(:,1)**2 + Sij(:,2)**2 + Sij(:,3)**2 
     &      + two*(Sij(:,4)**2 + Sij(:,5)**2 + Sij(:,6)**2)

       Snorm2(:) = sqrt( two*(Sij(:,1)**2 + Sij(:,2)**2 + Sij(:,3)**2) 
     &      + four*(Sij(:,4)**2 + Sij(:,5)**2 + Sij(:,6)**2) )

c... Let us build xmudmif at current quad pt. a la scatnu.f

       do n = 1,nshl
          xmudmif(:,intp) = xmudmif(:,intp) + 
     &         cdelsq(ien(:,n)) * Snorm2(:)*shp(n,intp)
       enddo          

      rmu=datmat(1,2,1)
      xmudmif(:,intp)=min(xmudmif(:,intp),1000.0*rmu(:)) !
c                                don't let it get larger than 1000 mu
      xmudmif(:,intp)=max(xmudmif(:,intp), zero) ! don't let (xmudmi) < 0

c.... Debugging

c      xmudmif(:,intp) = rmu(:)


c
c.... get necessary fluid properties (including the updated viscosity)
c
       do i = 1, npro
          do n = 1, nshl
             shpfun(i,n) = shp(n,intp)
          enddo
       enddo

        call getdiff(yl, shpfun, xmudmif,xl, rmu, rho)


       divqi = zero
       if ( idiff >= 1 ) then
c
c.... compute divergence of diffusive flux vector, qi,i
c
          do n=1, nshl
             divqi(:,1) = divqi(:,1) + shg(:,n,1)*ql(:,n,1 ) 
     &                               + shg(:,n,2)*ql(:,n,4 )
     &                               + shg(:,n,3)*ql(:,n,7 )

             divqi(:,2) = divqi(:,2) + shg(:,n,1)*ql(:,n,2 ) 
     &                               + shg(:,n,2)*ql(:,n,5 )
     &                               + shg(:,n,3)*ql(:,n,8)

             divqi(:,3) = divqi(:,3) + shg(:,n,1)*ql(:,n,3 ) 
     &                               + shg(:,n,2)*ql(:,n,6 )
     &                               + shg(:,n,3)*ql(:,n,9 )

          enddo

       endif                    ! diffusive flux computation

c
c.... take care of the body force term here
c
       src = zero
       if(matflg(5,1) .ge. 1) then
c
         bfx      = datmat(1,5,1) ! Boussinesq, g*alfap
         bfy      = datmat(2,5,1)
         bfz      = datmat(3,5,1)
         
         select case ( matflg(5,1) )
            case ( 1 )               ! standard linear body force
               src(:,1) = bfx
               src(:,2) = bfy
               src(:,3) = bfz
            case ( 2 )               ! boussinesq body force
               Temp = zero
               do n = 1, nshl 
                  Temp = Temp + shp(n,intp) * yl(:,n,5)
               enddo
               Tref = datmat(2,2,1)
               src(:,1) = bfx * (Temp(:)-Tref)
               src(:,2) = bfy * (Temp(:)-Tref)
               src(:,3) = bfz * (Temp(:)-Tref)
            case ( 3 )               ! user specified f(x,y,z)
               xx = zero
               do n  = 1,nenl
                  xx(:,1) = xx(:,1)  + shp(n,intp) * xl(:,n,1)
                  xx(:,2) = xx(:,2)  + shp(n,intp) * xl(:,n,2)
                  xx(:,3) = xx(:,3)  + shp(n,intp) * xl(:,n,3)
               enddo

               call e3source(xx, src)
          end select
            
       endif
c
c.... -------------------> Coriolis force  <-----------------
c
      omag=datmat(3,5,1)  ! frame rotation rate
       if(omag.ne.0) then
c
c.... unit vector of axis of rotation currently selecting the i,j,k
c
          e1 = one/sqrt(3.0d0)
          e2 = e1
          e3 = e1

          omega(1)=omag*e1
          omega(2)=omag*e2
          omega(3)=omag*e3

          if(matflg(5,1) .ne. 3) then ! we need to calculate the int pt. coords
             xx = zero
             do n  = 1,nenl
                xx(:,1) = xx(:,1)  + shp(n,intp) * xl(:,n,1)
                xx(:,2) = xx(:,2)  + shp(n,intp) * xl(:,n,2)
                xx(:,3) = xx(:,3)  + shp(n,intp) * xl(:,n,3)
             enddo

          endif
c
c  note that we calculate f as if it contains the usual source
c  plus the Coriolis and the centrifugal forces taken to the rhs (sign change)
c  as long as we are doing SUPG with no accounting for these terms in the
c  LHS this is the only change (which will find its way to the RHS momentum
c  equation (both Galerkin and SUPG parts)).
c
c  uncomment later if you want rotation always about z axis
c                 orig_src - om x om x r       - two om x u
c
c$$$          src(:,1)=src(:,1)+omega(3)*omega(3)*xx(:,1)+two*omega(3)*u2
c$$$          src(:,2)=src(:,2)+omega(3)*omega(3)*xx(:,2)-two*omega(3)*u1
c
c more general for testing
c
          src(:,1)=src(:,1)
     &            -(omega(2)*(omega(1)*xx(:,2)-omega(2)*xx(:,1))
     &             -omega(3)*(omega(3)*xx(:,1)-omega(1)*xx(:,3)))
     &            -two*(omega(2)*u3-omega(3)*u2)
          src(:,2)=src(:,2)
     &            -(omega(3)*(omega(2)*xx(:,3)-omega(3)*xx(:,2))
     &             -omega(1)*(omega(1)*xx(:,2)-omega(2)*xx(:,1)))
     &            -two*(omega(3)*u1-omega(1)*u3)
          src(:,3)=src(:,3)
     &            -(omega(1)*(omega(3)*xx(:,1)-omega(1)*xx(:,3))
     &             -omega(2)*(omega(2)*xx(:,3)-omega(3)*xx(:,2)))
     &            -two*(omega(1)*u2-omega(2)*u1)
       endif
c
c.... -------------------> momentum residual  <-----------------
c
       rLui(:,1) =(aci(:,1) + u1 * g1yi(:,2)
     &                      + u2 * g2yi(:,2)
     &                      + u3 * g3yi(:,2) - src(:,1) ) * rho
     &           + g1yi(:,1)
     &           - divqi(:,1)
       rLui(:,2) =(aci(:,2) + u1 * g1yi(:,3)
     &                      + u2 * g2yi(:,3)
     &                      + u3 * g3yi(:,3) - src(:,2) ) * rho
     &           + g2yi(:,1)
     &           - divqi(:,2)
       rLui(:,3) =(aci(:,3) + u1 * g1yi(:,4)
     &                      + u2 * g2yi(:,4)
     &                      + u3 * g3yi(:,4) - src(:,3) ) * rho
     &           + g3yi(:,1)
     &           - divqi(:,3)
       if(iconvflow.eq.1) then
          divu(:)  = (g1yi(:,2) + g2yi(:,3) + g3yi(:,4))*rho
          rLui(:,1)=rlui(:,1)+u1*divu
          rLui(:,2)=rlui(:,2)+u2*divu
          rLui(:,3)=rlui(:,3)+u3*divu
       endif

c
c.... compute the stabilization terms
c
        call e3stab (rho,          u1,       u2,
     &               u3,           dxidx,    rLui,   
     &               rmu,          tauC,     tauM,   
     &               tauBar,       uBar )         
c
c... Compute the SUPG stress at the current quad point multiplied
c... by the quadrature point weight.
c
        stressli(:,1) = u1(:)*rLui(:,1) 
        stressli(:,2) = u1(:)*rLui(:,2) 
        stressli(:,3) = u1(:)*rLui(:,3)
        stressli(:,4) = u2(:)*rLui(:,1)        
        stressli(:,5) = u2(:)*rLui(:,2)
        stressli(:,6) = u2(:)*rLui(:,3)        
        stressli(:,7) = u3(:)*rLui(:,1)
        stressli(:,8) = u3(:)*rLui(:,2)
        stressli(:,9) = u3(:)*rLui(:,3)

        if (iconvflow .eq. 1) then
           stressli(:,1) = stressli(:,1) + u1(:)*rLui(:,1) 
           stressli(:,2) = stressli(:,2) + u2(:)*rLui(:,1) 
           stressli(:,3) = stressli(:,3) + u3(:)*rLui(:,1)
           stressli(:,4) = stressli(:,4) + u1(:)*rLui(:,2)        
           stressli(:,5) = stressli(:,5) + u2(:)*rLui(:,2)
           stressli(:,6) = stressli(:,6) + u3(:)*rLui(:,2)        
           stressli(:,7) = stressli(:,7) + u1(:)*rLui(:,3)
           stressli(:,8) = stressli(:,8) + u2(:)*rLui(:,3)
           stressli(:,9) = stressli(:,9) + u3(:)*rLui(:,3)           
        endif

c.... Debugging

c        stressli = two
c        tauM     = one

c.... Multiply  ui*Luj times tauM and times WdetJ       

        do l = 1, 9
           do k = 1, npro
              stressli(k,l) = stressli(k,l)*WdetJ(k)*tauM(k)
           enddo
        enddo

c.... Obtain the SUPG energy dissipation (tau_{ij} S_{ij}) at the
c.... current qpt.

        dissi(:,1) = stressli(:,1)*Sij(:,1) + stressli(:,5)*Sij(:,2)
     &       + stressli(:,9)*Sij(:,3) + stressli(:,4)*Sij(:,4)
     &       + stressli(:,7)*Sij(:,5) + stressli(:,8)*Sij(:,6)
     &       + stressli(:,2)*Sij(:,4) + stressli(:,3)*Sij(:,5)
     &       + stressli(:,6)*Sij(:,6)

c.... Obtain the eddy viscosity dissipation multiplied by WdetJ

        dissi(:,2) = xmudmif(:,intp)*Snorm(:)*rho(:)*WdetJ(:)
        dissi(:,3) = rmu(:)*Snorm(:)*rho(:)*WdetJ(:) ! Total dissipation 
c                                             from molec. and eddy

c.... Debugging

c        dissi(:,1) = two*WdetJ(:)
c        dissi(:,2) = two*WdetJ(:)
c        dissi(:,3) = two*WdetJ(:)        

c..... Volume of element

        voli = WdetJ  ! Volume of element patch
c
c.... For debugging purposes let us keep track of rLui

c        rLui(:,1) = rLui(:,1)*WdetJ(:)
c        rLui(:,1) = rLui(:,2)*WdetJ(:)
c        rLui(:,1) = rLui(:,3)*WdetJ(:)

c.... Acumulate integration point contributions for each each element

        do l = 1, 9
           stressl(:,l) = stressl(:,l) + stressli(:,l)
        enddo
        
        do l = 1, 3
           dissl(:,l) = dissl(:,l) + dissi(:,l)
        enddo

        voll = voll + voli

      enddo    ! End loop over quadrature points.

      do j = 1,nshl
      do nel = 1,npro
        stress(ien(nel,j),:) = stress(ien(nel,j),:) + stressl(nel,:)
      enddo
      enddo
      
      do j = 1,nshl
      do nel = 1,npro
        diss(ien(nel,j),:) = diss(ien(nel,j),:) + dissl(nel,:)
      enddo
      enddo

      do j = 1,nshl
      do nel = 1,npro
         vol(ien(nel,j)) = vol(ien(nel,j)) + voll(nel)      
      enddo
      enddo

      return
      end
      subroutine cpjdmcnoi (y,      shgl,      shp, 
     &                   iper,   ilwork,       x,
     &                   rowp,   colm, 
     &                   iBC,    BC)

      use pointer_data

      use lhsGkeep ! This module stores the mass (Gram) matrix.

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwt contains the weights of the 
c                    quad. points.  

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension fres(nshg,24),         fwr(nshg),
     &          strnrm(nshg),         cdelsq(nshg),
     &          xnum(nshg),           xden(nshg),
     &          xmij(nshg,6),         xlij(nshg,6),
     &          xnude(nfath,2),        xnuder(nfath,2),
     &          strl(numel,ngauss),           
     &          y(nshg,5),            yold(nshg,5),
     &          iper(nshg),
     &          ilwork(nlwork),
     &          x(numnp,3),
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),    
     &          pfres(nshg,22),                ifath(nshg),
     &          nsons(nshg),                   iBC(nshg),
     &          BC(nshg,ndofBC),               xnutf(nfath),
     &          xnut(nshg)

        integer   rowp(nshg*nnz),         colm(nshg+1)

        real*8, allocatable, dimension(:,:,:) :: em 


      denom=max(1.0d0*(lstep),one)
      if(dtavei.lt.0) then
         wcur=one/denom
      else
         wcur=dtavei
      endif  
      whist=1.0-wcur
     
      if (istep .eq. 0) then
         lhsG = zero
      endif

      fres = zero
      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)
c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      do i = 1, nshg  ! No periodicity for testing
c      iper(i) = i
c      enddo
c      yold(:,5) = 1.0
c      yold(:,2) = 3.0d0
c      yold(:,2) = 2.0*x(:,1) - 3*x(:,2) 
c      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
c      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
c                               suitable for the


      intrul=intg(1,itseq)
      intind=intpt(intrul)

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        call hfilterBB (yold, x, mien(iblk)%p, fres, 
     &               shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))


        if ( istep.eq.0 ) then

           allocate ( em(npro,nshl,nshl) )

           call getgram2 (x, mien(iblk)%p, 
     &          shgl(lcsyst,:,1:nshl,:),  shp(lcsyst,1:nshl,:),     
     &          shglf(lcsyst,:,1:nshl,:), shpf(lcsyst,1:nshl,:), em, 
     &          Qwtf(lcsyst,1:ngaussf))

           call fillsparseSclr (mien(iblk)%p, 
     &                          em,            lhsG,
     &                          rowp,          colm)


           deallocate ( em )

        endif

      enddo   ! End loop over element blocks
c

c      write(*,*)'Im here'
 
      if(numpe>1) call commu (fres, ilwork, 24, 'in ')
c 
c account for periodicity in filtered variables
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:)
        endif
      enddo

      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = zero
        endif
      enddo

c     Need to zero off-processor slaves as well.

      if (numpe.gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  

         do itask = 1, numtask
            
            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  fres(isgbeg:isgend,:) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif

c... At this point fres has the right hand side vector (b) and lhsG has
c... the Gram matrix (M_{AB}) (in sparse storage). Now we need to solve
c... Ax = b using the conjugate gradient method to finish off the 
c... L2-projection.


      do i = 1, 21
         call sparseCG (fres(:,i), pfres(:,i), lhsG, 
     &        rowp, colm, iper, ilwork,
     &        iBC,  BC)
      enddo


      write(*,*)'Done with least-squares projection'
      
      do i = 1, 21
         fres(:,i) = pfres(:,i)
      enddo

      fres(:,22) = one

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)
 
        call getstrl (yold, x,      mien(iblk)%p,  
     &               strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:))

      enddo


      strnrm = sqrt( 
     &  two * (fres(:,10)**2 + fres(:,11)**2 + fres(:,12)**2)
     &  + four * ( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )

      fwr = fwr1 * fres(:,22) * strnrm

      xmij(:,1) = -fwr
     &             * fres(:,10) + fres(:,16)
      xmij(:,2) = -fwr
     &             * fres(:,11) + fres(:,17) 
      xmij(:,3) = -fwr
     &             * fres(:,12) + fres(:,18) 

      xmij(:,4) = -fwr * fres(:,13) + fres(:,19)
      xmij(:,5) = -fwr * fres(:,14) + fres(:,20)
      xmij(:,6) = -fwr * fres(:,15) + fres(:,21)

      fres(:,22) = one / fres(:,22)

      xlij(:,1) = fres(:,4) - fres(:,1) * fres(:,1) * fres(:,22)
      xlij(:,2) = fres(:,5) - fres(:,2) * fres(:,2) * fres(:,22)
      xlij(:,3) = fres(:,6) - fres(:,3) * fres(:,3) * fres(:,22)
      xlij(:,4) = fres(:,7) - fres(:,1) * fres(:,2) * fres(:,22)
      xlij(:,5) = fres(:,8) - fres(:,1) * fres(:,3) * fres(:,22)
      xlij(:,6) = fres(:,9) - fres(:,2) * fres(:,3) * fres(:,22)

      xnum =        xlij(:,1) * xmij(:,1) + xlij(:,2) * xmij(:,2) 
     &                                    + xlij(:,3) * xmij(:,3)
     &     + two * (xlij(:,4) * xmij(:,4) + xlij(:,5) * xmij(:,5)
     &                                    + xlij(:,6) * xmij(:,6))
      xden =        xmij(:,1) * xmij(:,1) + xmij(:,2) * xmij(:,2) 
     &                                    + xmij(:,3) * xmij(:,3)
     &     + two * (xmij(:,4) * xmij(:,4) + xmij(:,5) * xmij(:,5)
     &                                    + xmij(:,6) * xmij(:,6))
      xden = two * xden

c 
c don't account for periodic nodes twice
c
      do j = 1,numnp
        i = iper(j)
        if (i .ne. j) then
           xden(j) = zero
           xnum(j) = zero
        endif
      enddo

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
                     xnum(isgbeg:isgend) = zero
                     xden(isgbeg:isgend) = zero
                  enddo
               endif
               
               itkbeg = itkbeg + 4 + 2*numseg
               
            enddo

c            if (myrank.eq.0)then
c               do i = 1, numnp
c                  write(253,*)xnum(i),xden(i),myrank
c               enddo
c            endif
c            if (myrank.eq.1)then
c               do i = 1, numnp
c                  write(254,*)xnum(i),xden(i),myrank
c               enddo
c            endif

c            xnuml = sum(xnum)
c            xdenl = sum(xden)

            xnuml = zero
            xdenl = zero
            do i = 1, numnp
               xnuml = xnuml + xnum(i)
               xdenl = xdenl + xden(i)
            enddo

c            write(*,*)xnuml,xdenl,myrank
            
            call drvAllreducesclr ( xnuml, xnumt )
            call drvAllreducesclr ( xdenl, xdent )
cd 
         else

c            xnumt = sum(xnum)
c            xdent = sum(xden)
            xnumt = zero
            xdent = zero
            do i = 1, numnp
               xnumt = xnumt + xnum(i)
               xdent = xdent + xden(i)
            enddo


         endif

         scalar = xnumt / (xdent + 1.d-09)      
         xnut = scalar


      if (myrank .eq. 0)then
         write(*,*) 'xnut=', xnut(100)
      endif
c      do i = 1, numnp
c         write(*,*)xnumt/xdent,myrank
c      enddo
c
      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         
         ngauss = nint(lcsyst)

         call scatnu (mien(iblk)%p, strl(iel:inum,:), 
     &        mxmudmi(iblk)%p,cdelsq,shp(lcsyst,1:nshl,:))
      enddo
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$  tmp1 =  MINVAL(xmudmi)
c$$$  tmp2 =  MAXVAL(xmudmi)
c$$$  if(numpe>1) then
c$$$  call MPI_REDUCE (tmp1, tmp3, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MIN, master, MPI_COMM_WORLD, ierr)
c$$$  call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MAX, master, MPI_COMM_WORLD, ierr)
c$$$      tmp1=tmp3
c$$$  tmp2=tmp4
c$$$  endif
c$$$  if (myrank .EQ. master) then
c$$$  write(35,*) lstep,tmp1,tmp2
c$$$  call flush(35)
c$$$  endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$

c
c  if flag set, write a restart file with info (reuse xmij's memory)
c
      if(irs.eq.11) then
         lstep=999
         xmij(:,1)=xnum(:)
         xmij(:,2)=xden(:)
         xmij(:,3)=cdelsq(:)
         xmij(:,5)=xlij(:,4)    !leave M_{12} in 4 and put L_{12} here
         call restar('out ',xmij,xlij) !also dump all of L_{ij} in ac
         stop
      endif
c
c  local clipping moved to scatnu with the creation of mxmudmi pointers
c
c$$$      rmu=datmat(1,2,1)
c$$$      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
c$$$      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0
c      stop !uncomment to test dmod
c


c  write out the nodal values of xnut (estimate since we don't calc strain
c  there and must use the filtered strain).
c



      return
      end
      
c-----------------------------------------------------

      subroutine getgram (x, ien, shgl, shp, em, Qwtf)

      include "common.h"

      dimension x(numnp,nsd),            xl(npro,nenl,nsd)      
      dimension ien(npro,nshl),
     &          shgl(nsd,nshl,ngauss),    shp(nshl,ngauss),
     &          em(npro,nshl,nshl),      Qwtf(ngaussf)
      
      call localx(x,      xl,     ien,    nsd,    'gather  ')

      call cmass(shp,shgl,xl,em)
         

      return

      end

c----------------------------------------------------------------------


      subroutine getgram2 (x, ien, shgl, shp, shglf, shpf, em, Qwtf)

      include "common.h"

      dimension x(numnp,nsd),            xl(npro,nenl,nsd)      
      dimension ien(npro,nshl),
     &          shgl(nsd,nshl,ngauss),    shp(nshl,ngauss),
     &          shglf(nsd,nshl,ngauss),   shpf(nshl,ngauss),
     &          em(npro,nshl,nshl),      Qwtf(ngaussf) 

      
      call localx(x,      xl,     ien,    nsd,    'gather  ')

      call cmassl(shp,shgl,shpf,shglf,xl,em,Qwtf)
         

      return

      end

c-----------------------------------------------------------------------

      subroutine getgram3 (x, ien, shgl, shp, shglf, shpf, em, Qwtf)

      include "common.h"

      dimension x(numnp,nsd),            xl(npro,nenl,nsd)      
      dimension ien(npro,nshl),
     &          shgl(nsd,nshl,ngauss),    shp(nshl,ngauss),
     &          shglf(nsd,nshl,ngauss),   shpf(nshl,ngauss),
     &          em(npro,nshl,nshl),      Qwtf(ngaussf) 

      
      call localx(x,      xl,     ien,    nsd,    'gather  ')

      call cmasstl(shp,shgl,shpf,shglf,xl,em,Qwtf)
         

      return

      end
      subroutine cdelBHsq (y,      shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x, cdelsq1)

      use pointer_data

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwtf contains the weights of the 
c                    quad. points.  

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension fres(nshg,33),         fwr(nshg),
     &          strnrm(nshg),         cdelsq1(nfath),
     &          xnum(nshg),           xden(nshg),
     &          xmij(nshg,6),         xlij(nshg,6),
     &          xnude(nfath,2),        xnuder(nfath,2),
     &          nsons(nshg),
     &          strl(numel,ngauss),           
     &          y(nshg,5),            yold(nshg,5),
     &          ifath(nshg),          iper(nshg),
     &          ilwork(nlwork),
     &          x(numnp,3),
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),    
     &          xnutf(nfath),
     &          hfres(nshg,16)      

c
     
      fres = zero
      hfres = zero

      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)

c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      yold(:,5) = 1.0  ! Debugging
c      yold(:,2) = 2.0*x(:,1) - 3.0*x(:,2) 
c      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
c      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
c                               suitable for the

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
c        call hfilterB (yold, x, mien(iblk)%p, hfres, 
c     &               shglf(lcsyst,:,1:nshl,:),
c     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

        call hfilterC (yold, x, mien(iblk)%p, hfres, 
     &               shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo

      if(numpe>1) call commu (hfres, ilwork, 16, 'in ')
c 
c... account for periodicity in filtered variables
c
      do j = 1,nshg  !    Add on-processor slave contribution to masters
        i = iper(j)
        if (i .ne. j) then
           hfres(i,:) = hfres(i,:) + hfres(j,:)
        endif
      enddo
      do j = 1,nshg ! Set on-processor slaves to be the same as masters
        i = iper(j)
        if (i .ne. j) then
           hfres(j,:) = hfres(i,:)
        endif
      enddo

c... Set off-processor slaves to be the same as their masters

      if(numpe>1)   call commu (hfres, ilwork, 16, 'out')


      hfres(:,16) = one / hfres(:,16) ! one/(volume of hat filter kernel)

      do j = 1, 15
	hfres(:,j) = hfres(:,j) * hfres(:,16)
      enddo	    		

c... For debugging

c      hfres(:,1) = 2.0*x(:,1) - 3.0*x(:,2) 
c      hfres(:,2) = 3.0*x(:,1) + 4.0*x(:,2)
c      hfres(:,3) = 4.0*x(:,1) + x(:,2) + x(:,3)

c... Done w/ h-filtering. Begin 2h-filtering.

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        call twohfilterB (yold, x, strl(iel:inum,:), mien(iblk)%p, 
     &               fres, hfres, shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo
c
 

      if(numpe>1) call commu (fres, ilwork, 33, 'in ')
c 
c account for periodicity in filtered variables
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:)
        endif
      enddo

      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = fres(i,:)
        endif
      enddo

      if(numpe>1)then
         call commu (fres, ilwork, 33, 'out')
      endif

      fres(:,22) = one / fres(:,22)
      do j = 1,21
        fres(:,j) = fres(:,j) * fres(:,22)
      enddo
      do j = 23,33
        fres(:,j) = fres(:,j) * fres(:,22)
      enddo

      
      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)
 
        call getstrl (yold, x,      mien(iblk)%p,  
     &               strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:))

      enddo

c
c... Obtain the hat-tilde strain rate norm at the nodes 
c

      strnrm = sqrt( 
     &  two * (fres(:,10)**2 + fres(:,11)**2 + fres(:,12)**2)
     &  + four * ( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )

      fwr = fwr1 * strnrm

      xmij(:,1) = -fwr
     &             * fres(:,10) + fres(:,16)
      xmij(:,2) = -fwr
     &             * fres(:,11) + fres(:,17) 
      xmij(:,3) = -fwr
     &             * fres(:,12) + fres(:,18) 

      xmij(:,4) = -fwr * fres(:,13) + fres(:,19)
      xmij(:,5) = -fwr * fres(:,14) + fres(:,20)
      xmij(:,6) = -fwr * fres(:,15) + fres(:,21)


      xlij(:,1) = fres(:,4) - fres(:,1) * fres(:,1) 
      xlij(:,2) = fres(:,5) - fres(:,2) * fres(:,2) 
      xlij(:,3) = fres(:,6) - fres(:,3) * fres(:,3) 
      xlij(:,4) = fres(:,7) - fres(:,1) * fres(:,2) 
      xlij(:,5) = fres(:,8) - fres(:,1) * fres(:,3) 
      xlij(:,6) = fres(:,9) - fres(:,2) * fres(:,3) 

      xnum =        xlij(:,1) * xmij(:,1) + xlij(:,2) * xmij(:,2) 
     &                                    + xlij(:,3) * xmij(:,3)
     &     + two * (xlij(:,4) * xmij(:,4) + xlij(:,5) * xmij(:,5)
     &                                    + xlij(:,6) * xmij(:,6))
      xden =        xmij(:,1) * xmij(:,1) + xmij(:,2) * xmij(:,2) 
     &                                    + xmij(:,3) * xmij(:,3)
     &     + two * (xmij(:,4) * xmij(:,4) + xmij(:,5) * xmij(:,5)
     &                                    + xmij(:,6) * xmij(:,6))
      xden = two * xden

c  zero on processor periodic nodes so that they will not be added twice
        do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xnum(j) = zero
            xden(j) = zero
          endif
        enddo

      if (numpe.gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xnum(isgbeg:isgend) = zero
                  xden(isgbeg:isgend) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
c
c Description of arrays.   Each processor has an array of length equal
c to the total number of fathers times 2 xnude(nfathers,2). One to collect 
c the numerator and one to collect the denominator.  There is also an array
c of length nshg on each processor which tells the father number of each
c on processor node, ifath(nnshg).  Finally, there is an arry of length
c nfathers to tell the total (on all processors combined) number of sons
c for each father. 
c
c  Now loop over nodes and accumlate the numerator and the denominator
c  to the father nodes.  Only on processor addition at this point.
c  Note that serrogate fathers are collect some for the case where some
c  sons are on another processor
c
      xnude = zero
      do i = 1,nshg
         xnude(ifath(i),1) = xnude(ifath(i),1) + xnum(i)
         xnude(ifath(i),2) = xnude(ifath(i),2) + xden(i)
      enddo

c
c Now  the true fathers and serrogates combine results and update
c each other.
c       
      if(numpe .gt. 1)then
         call drvAllreduce(xnude, xnuder,2*nfath)
c
c  xnude is the sum of the sons for each father on this processor
c
c  xnuder is the sum of the sons for each father on all processor combined
c  (the same as if we had not partitioned the mesh for each processor)
c
c   For each father we have precomputed the number of sons (including
c   the sons off processor). 
c
c   Now divide by number of sons to get the average (not really necessary
c   for dynamic model since ratio will cancel nsons at each father)
c
c         xnuder(:,1) = xnuder(:,1) ! / nsons(:)
c         xnuder(:,2) = xnuder(:,2) ! / nsons(:)
c
c  the next line is c \Delta^2
c
         xnuder(:,1) = xnuder(:,1) / (xnuder(:,2) + 1.d-09)
         do i = 1,nfath
            cdelsq1(i) = xnuder(i,1)
         enddo
      else
c     
c     the next line is c \Delta^2, not nu_T but we want to save the
c     memory
c     
         xnude(:,1) = xnude(:,1) / (xnude(:,2) + 1.d-09)
         do i = 1,nfath
            cdelsq1(i) = xnude(i,1)
         enddo
      endif

      if (myrank .eq. master) then
         if (numpe .gt. 1) then
            do i = 1, nfath
               write(22,*)i, xnuder(i,1)
            enddo
         else
            do i = 1, nfath
               write(22,*)i, xnude(i,1)           
            enddo             
         endif
      endif
      call flush(22)

      do i = 1, nfath
         if (cdelsq1(i) .lt. zero) then
            cdelsq1(i) = zero
         endif
      enddo

      return
      end
      subroutine SUPGdis (y,           ac,         shgl,      
     &                  shp,         iper,       ilwork,    
     &                  nsons,       ifath,      x,
     &                  iBC,    BC,  stabdis,    xavegt)


      use stats            !  
      use pointer_data     ! brings in the pointers for the blocked arrays
      use local_mass
      use rlssave  ! Use the resolved Leonard stresses at the nodes.      
      use quadfilt ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwt contains the weights of the 
c                    quad. points.  



      include "common.h"
      include "mpif.h"
      include "auxmpi.h"


      dimension y(nshg,ndof),                  ac(nshg,ndof), 
     &          yold(nshg,ndof),
     &          ifath(nshg),                   nsons(nshg),
     &          iper(nshg),                    ilwork(nlwork),        
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),
     &          x(numnp,3),           
     &          qres(nshg,nsd*nsd),             rmass(nshg),  
     &          iBC(nshg),                      BC(nshg,ndofBC),
     &          cdelsq(nshg),                   vol(nshg),
     &          stress(nshg,9),                 diss(nshg,3),
     &          xave(nshg,12),                  xaveg(nfath,12),
     &          xavegr(nfath,12),               stabdis(nfath),
     &          dmodc(nfath),                   strl(numel,ngauss),
     &          xavegt(nfath,12)
 
      character*5  cname
      character*30 fname 
      
      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)     

c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      yold(:,5) = 1.0  ! Debugging
c      yold(:,2) = 2.0*x(:,1) - 3.0*x(:,2) 
c      yold(:,2) = 2.0
c      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      yold(:,3) = 3.0
c      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
c      yold(:,4) = 4.0
c      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
c                               suitable for the

c.... First let us obtain cdelsq at each node in the domain.
c.... We use numNden which lives in the quadfilt module.

      if ( (istep .eq. 0) ) then
         fname =  'dmodc.dat' // cname (myrank+1)
         open (99,file=fname,form='unformatted',status='unknown')
         read(99) dmodc
         close(99)
         cdelsq(:) = dmodc(ifath(:))
      else
         cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)       
      endif

c      if (myrank .eq. master) then
c         do i = 1, nfath
c            write(*,*)'dmod=', dmodc(i)
c         enddo
c      endif

      if ( istep .eq. (nstep(1)-1) ) then
         dmodc(ifath(:)) = cdelsq(:)
         fname =  'dmodc.dat' // cname (myrank+1)
         open (99,file=fname,form='unformatted', status='replace')
         write(99) dmodc         
         close(99)
c         if (myrank .eq. master) then
c            do i = 1, nfath
c               write(*,*)'dmod=', dmodc(i)
c            enddo
c         endif

      endif

c      if (istep .eq. 0)
      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)
 
        call getstrl (yold, x,      mien(iblk)%p,  
     &               strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:))

      enddo

      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         
         ngauss = nint(lcsyst)

         call scatnu (mien(iblk)%p, strl(iel:inum,:), 
     &        mxmudmi(iblk)%p,cdelsq,shp(lcsyst,1:nshl,:))
      enddo
c      endif



        if (idiff==1 .or. idiff==3) then ! global reconstruction of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
           qres = zero
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)
c     
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

              call AsIq (y,                x,                       
     &                   shp(lcsyst,1:nshl,:), 
     &                   shgl(lcsyst,:,1:nshl,:),
     &                   mien(iblk)%p,     mxmudmi(iblk)%p,  
     &                   qres,             rmass )
           enddo
       
c
c.... form the diffusive flux approximation
c
           call qpbc( rmass, qres, iBC, BC, iper, ilwork )       
c
        endif 


c.... form the SUPG stresses well as dissipation due to eddy viscosity,
c...  and SUPG stabilization.


        stress = zero
        vol    = zero
        diss   = zero

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        call SUPGstress (y, ac, x, qres, mien(iblk)%p, mxmudmi(iblk)%p, 
     &                   cdelsq, shglf(lcsyst,:,1:nshl,:),
     &                   shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf),
     &                   shgl(lcsyst,:,1:nshl,:), shp(lcsyst,1:nshl,:),
     &                   stress, diss, vol)

      enddo

      if(numpe>1) call commu (stress, ilwork, 9, 'in ')      
      if(numpe>1) call commu (diss, ilwork, 3, 'in ')       
      if(numpe>1) call commu (vol, ilwork, 1, 'in ')      

c 
c account for periodicity 
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           stress(i,:) = stress(i,:) + stress(j,:)
           diss(i,:)   = diss(i,:)   + diss(j,:)
           vol(i)      = vol(i)      + vol(j)
        endif
      enddo

      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           stress(j,:) = stress(i,:) 
           diss(j,:)   = diss(i,:)   
           vol(j)      = vol(i)      
        endif
      enddo      

      if(numpe>1) call commu (stress, ilwork, 9, 'out ')      
      if(numpe>1) call commu (diss, ilwork, 3, 'out ')       
      if(numpe>1) call commu (vol, ilwork, 1, 'out ')      

      vol = one / vol
      do i = 1, 9
         stress(:,i) = stress(:,i)*vol(:)
      enddo
      do i = 1, 3
         diss(:,i) = diss(:,i)*vol(:)
      enddo

c---------- > Begin averaging dissipations and SUPG stress <--------------

      do i = 1, 9
         xave(:,i) = stress(:,i)
      enddo
      xave(:,10) = diss(:,1)
      xave(:,11) = diss(:,2)
      xave(:,12) = diss(:,3)

c  zero on processor periodic nodes so that they will not be added twice
        do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xave(j,:) = zero
          endif
        enddo

      if (numpe.gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xave(isgbeg:isgend,:) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
c

      xaveg = zero
      do i = 1,nshg      
         xaveg(ifath(i),:) = xaveg(ifath(i),:) + xave(i,:)
      enddo

      if(numpe .gt. 1)then
         call drvAllreduce(xaveg, xavegr,12*nfath)         

         do m = 1, 12
            xavegr(:,m) = xavegr(:,m)/nsons(:)
         enddo

c         if (myrank .eq. master) then
c            write(*,*)'diss=', xavegt(14,11), xavegr(14,11)
c         endif

         do m = 1, 12
            xavegt(:,m) = xavegt(:,m) + xavegr(:,m)
         enddo

         stabdis(:) = xavegr(:,10)

      else

         do m = 1, 12
            xaveg(:,m) = xaveg(:,m)/nsons(:)
         enddo         
         
         do m = 1, 12
            xavegt(:,m) = xavegt(:,m) + xaveg(:,m)
         enddo         
         
         stabdis(:) = xaveg(:,10)

      endif

c      if (myrank .eq. master) then
c         write(*,*)'diss=', xavegt(14,11), xavegr(14,11)
c      endif

       if ( istep .eq. (nstep(1)-1) ) then
          if ( myrank .eq. master) then

             do i = 1, nfath
c               write(376,*)xavegt(i,1),xavegt(i,2),xavegt(i,3)
c               write(377,*)xavegt(i,4),xavegt(i,5),xavegt(i,6)
c               write(378,*)xavegt(i,7),xavegt(i,8),xavegt(i,9)
                write(380,*)xavegt(i,10),xavegt(i,11),xavegt(i,12)
            enddo

c            call flush(376)
c            call flush(377)
c            call flush(378)
            call flush(380)
            
         endif
      endif      


      return

      end
      subroutine dmcSUPG(y,           ac,         shgl,      
     &                  shp,         iper,       ilwork,    
     &                  nsons,       ifath,      x,
     &                  iBC,    BC,  rowp,       colm,
     &                  xavegt, stabdis) 

      use lhsGkeep ! This module stores the mass (Gram) matrix.
      use stats            !  
      use pointer_data     ! brings in the pointers for the blocked arrays
      use local_mass
      use rlssave  ! Use the resolved Leonard stresses at the nodes.      
      use quadfilt ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwt contains the weights of the 
c                    quad. points.  



      include "common.h"
      include "mpif.h"
      include "auxmpi.h"


      dimension y(nshg,ndof),                  ac(nshg,ndof), 
     &          ifath(nshg),                   nsons(nshg),
     &          iper(nshg),                    ilwork(nlwork),        
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),
     &          x(numnp,3),           
     &          qres(nshg,nsd*nsd),             rmass(nshg),  
     &          iBC(nshg),                      BC(nshg,ndofBC),
     &          cdelsq(nshg),                   vol(nshg),
     &          stress(nshg,9),                 diss(nshg,3),
     &          xave(nshg,12),                  xaveg(nfath,12),
     &          xavegr(nfath,12),               stabdis(nfath),
     &          yold(nshg,ndof),                xavegt(nfath,12),
     &          fres(nshg,24),                  pfres(nshg,24),
     &          cdel(nfath),                    xnume(nfath),
     &          xdeno(nfath),                    strl(numel,ngauss),
     &          rden(nshg),                     rnum(nshg)


      integer   rowp(nshg*nnz),         colm(nshg+1)
      
      real*8, allocatable, dimension(:,:,:) :: em 

      real*8, allocatable, dimension(:,:) :: fakexmu      


      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)
      fres = zero

c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      yold(:,5) = 1.0  ! Debugging
c      yold(:,2) = 2.0*x(:,1) - 3.0*x(:,2) 
c      yold(:,2) = 2.0
c      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      yold(:,3) = 3.0
c      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
c      yold(:,4) = 4.0
c      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
c                               suitable for the


      intrul=intg(1,itseq)
      intind=intpt(intrul)

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        call resSij (yold, x, mien(iblk)%p, fres, 
     &               shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

        if ( istep.eq.0 ) then

           allocate ( em(npro,nshl,nshl) )

           call getgram2 (x, mien(iblk)%p, 
     &          shgl(lcsyst,:,1:nshl,:),  shp(lcsyst,1:nshl,:),     
     &          shglf(lcsyst,:,1:nshl,:), shpf(lcsyst,1:nshl,:), em, 
     &          Qwtf(lcsyst,1:ngaussf))

c           call getgram (x, mien(iblk)%p, 
c     &          shgl(lcsyst,:,1:nshl,:),  shp(lcsyst,1:nshl,:),     
c     &          em, Qwtf(lcsyst,1:ngaussf))

           call fillsparseSclr (mien(iblk)%p, 
     &                          em,            lhsG,
     &                          rowp,          colm)


           deallocate ( em )

        endif

      enddo   ! End loop over element blocks
c

      if(numpe>1) call commu (fres, ilwork, 24, 'in ')
c 
c account for periodicity in filtered variables
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:)
        endif
      enddo
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = fres(i,:)
        endif
      enddo

      if(numpe>1)   call commu (fres, ilwork, 24, 'out')

      fres(:,22) = one / fres(:,22)
      do j = 1,21
        fres(:,j) = fres(:,j) * fres(:,22)
      enddo
      pfres = fres

c---- Needed for consistent projection
 
c      if(numpe>1) call commu (fres, ilwork, 24, 'in ')
c 
c account for periodicity in filtered variables
c
c      do j = 1,nshg
c        i = iper(j)
c        if (i .ne. j) then
c           fres(i,:) = fres(i,:) + fres(j,:)
c        endif
c      enddo

c      do j = 1,nshg
c        i = iper(j)
c        if (i .ne. j) then
c           fres(j,:) = zero
c        endif
c      enddo

c     Need to zero off-processor slaves as well.

c      if (numpe.gt.1 .and. nsons(1).gt.1) then

c         numtask = ilwork(1)
c         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  

c         do itask = 1, numtask
            
c            iacc   = ilwork (itkbeg + 2)
c            numseg = ilwork (itkbeg + 4)

c            if (iacc .eq. 0) then
c               do is = 1,numseg
c                  isgbeg = ilwork (itkbeg + 3 + 2*is)
c                  lenseg = ilwork (itkbeg + 4 + 2*is)
c                  isgend = isgbeg + lenseg - 1
c                  fres(isgbeg:isgend,:) = zero
c               enddo
c            endif
            
c            itkbeg = itkbeg + 4 + 2*numseg
            
c         enddo
         
c      endif

c... At this point fres has the right hand side vector (b) and lhsG has
c... the Gram matrix (M_{AB}) (in sparse storage). Now we need to solve
c... Ax = b using the conjugate gradient method to finish off the 
c... L2-projection.


c      do i = 16, 16
c         call sparseCG (fres(:,i), pfres(:,i), lhsG, 
c     &        rowp, colm, iper, ilwork, 
c     &        iBC,  BC)
c         write(*,*) 'i=', i
c      enddo


c      write(*,*)'Done with least-squares projection'
      




c.... First let us obtain cdelsq at each node in the domain.
c.... We use numNden which lives in the quadfilt module.

      cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)       
c      cdelsq(:) = zero ! Debugging

      if (istep .eq. 0) then
         xavegt = zero  ! For averaging dissipations and SUPG stresses
      endif

        if (idiff==1 .or. idiff==3) then ! global reconstruction of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
           qres = zero
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)

              allocate ( fakexmu(npro,ngauss) )              
              fakexmu = zero

c     
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

              call AsIq (y,                x,                       
     &                   shp(lcsyst,1:nshl,:), 
     &                   shgl(lcsyst,:,1:nshl,:),
     &                   mien(iblk)%p,     mxmudmi(iblk)%p,
     &                   qres,             rmass )

              deallocate ( fakexmu )
           enddo
       
c
c.... form the diffusive flux approximation
c
           call qpbc( rmass, qres, iBC, BC, iper, ilwork )       
c
        endif 


c.... form the SUPG stresses well as dissipation due to eddy viscosity,
c...  and SUPG stabilization.


        stress = zero
        vol    = zero
        diss   = zero

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        allocate ( fakexmu(npro,ngauss) )              
        fakexmu = zero

        call SUPGstress (y, ac, x, qres, mien(iblk)%p, fakexmu, 
     &                   cdelsq, shglf(lcsyst,:,1:nshl,:),
     &                   shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf),
     &                   shgl(lcsyst,:,1:nshl,:), shp(lcsyst,1:nshl,:),
     &                   stress, diss, vol)
        
        deallocate ( fakexmu )
      enddo

      if(numpe>1) call commu (stress, ilwork, 9, 'in ')      
      if(numpe>1) call commu (diss, ilwork, 3, 'in ')       
      if(numpe>1) call commu (vol, ilwork, 1, 'in ')      

c 
c account for periodicity 
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           stress(i,:) = stress(i,:) + stress(j,:)
           diss(i,:)   = diss(i,:)   + diss(j,:)
           vol(i)      = vol(i)      + vol(j)
        endif
      enddo

      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           stress(j,:) = stress(i,:) 
           diss(j,:)   = diss(i,:)   
           vol(j)      = vol(i)      
        endif
      enddo      

      if(numpe>1) call commu (stress, ilwork, 9, 'out ')      
      if(numpe>1) call commu (diss, ilwork, 3, 'out ')       
      if(numpe>1) call commu (vol, ilwork, 1, 'out ')      

      vol = one / vol
      do i = 1, 9
         stress(:,i) = stress(:,i)*vol(:)
      enddo
      do i = 1, 3
         diss(:,i) = diss(:,i)*vol(:)
      enddo

c---------- > Begin averaging dissipations and SUPG stress <--------------

      do i = 1, 9
         xave(:,i) = stress(:,i)
      enddo
      xave(:,10) = diss(:,1)
      xave(:,11) = diss(:,2)
      xave(:,12) = pfres(:,16)

c  zero on processor periodic nodes so that they will not be added twice
        do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xave(j,:) = zero
          endif
        enddo

      if (numpe.gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xave(isgbeg:isgend,:) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
c

      xaveg = zero
      do i = 1,nshg      
         xaveg(ifath(i),:) = xaveg(ifath(i),:) + xave(i,:)
      enddo

      if(numpe .gt. 1)then
         call drvAllreduce(xaveg, xavegr,12*nfath)         

         do m = 1, 12
            xavegr(:,m) = xavegr(:,m)/nsons(:)
         enddo

c         if (myrank .eq. master) then
c            write(*,*)'diss=', xavegt(14,11), xavegr(14,11)
c         endif

         do m = 1, 12
            xavegt(:,m) = xavegt(:,m) + xavegr(:,m)
         enddo

      else

         do m = 1, 12
            xaveg(:,m) = xaveg(:,m)/nsons(:)
         enddo         
         
         do m = 1, 12
            xavegt(:,m) = xavegt(:,m) + xaveg(:,m)
         enddo         
         
      endif

      if (myrank .eq. master) then
         write(*,*)'diss0=', xavegt(14,11), xavegr(14,11)
      endif

      if ( istep .eq. (nstep(1)-1) ) then
         if ( myrank .eq. master) then

            do i = 1, nfath
c               write(376,*)xavegt(i,1),xavegt(i,2),xavegt(i,3)
c               write(377,*)xavegt(i,4),xavegt(i,5),xavegt(i,6)
c               write(378,*)xavegt(i,7),xavegt(i,8),xavegt(i,9)
               write(381,*)xavegt(i,10),xavegt(i,11),xavegt(i,12)
            enddo

c            call flush(376)
c            call flush(377)
c            call flush(378)
c            call flush(379)
            call flush(381)
         endif
      endif      

      rnum(ifath(:)) = numNden(:,1)
      rden(ifath(:)) = numNden(:,2)

      if (numpe .gt. 1) then
      do i = 1, nfath
         if (stabdis(i) .gt. zero) then
            cdel(i) = (two*xavegr(i,11)-stabdis(i))/xavegr(i,12)
            xnume(i) = two*xavegr(i,11)-stabdis(i)
            xdeno(i) = xavegr(i,12)
         else
            xnume(i) = rnum(i)
            xdeno(i) = rden(i)
         endif
      enddo
      else
      do i = 1, nfath
         if (stabdis(i) .gt. zero) then
            cdel(i) = (two*xaveg(i,11)-stabdis(i))/xaveg(i,12)
            xnume(i) = two*xaveg(i,11)-stabdis(i)
            xdeno(i) = xaveg(i,12)
         else
            xnume(i) = rnum(i)
            xdeno(i) = rden(i)
         endif
      enddo
      endif

      do i = 1, nfath
         if (xnume(i) .lt. zero) then
            xnume(i) = rnum(i)
            xdeno(i) = rden(i)
         endif
      enddo

      do i = 1, nshg
            numNden(i,1) = xnume(ifath(i))
            numNden(i,2) = xdeno(ifath(i))
      enddo

      cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)
 
        call getstrl (yold, x,      mien(iblk)%p,  
     &               strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:))

      enddo


      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         
         ngauss = nint(lcsyst)

         call scatnu (mien(iblk)%p, strl(iel:inum,:), 
     &        mxmudmi(iblk)%p,cdelsq,shp(lcsyst,1:nshl,:))
      enddo

      return

      end
      subroutine FiltRat (y,      shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x,   cdelsq1, fwr4, 
     &                   fwr3)

      use pointer_data

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwt contains the weights of the 
c                    quad. points.  

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension fres(nshg,24),         fwr(nshg),
     &          strnrm(nshg),         cdelsq1(nfath),
     &          xnum(nshg),           xden(nshg),
     &          xmij(nshg,6),         xlij(nshg,6),
     &          xnude(nfath,5),        xnuder(nfath,5),
     &          nsons(nshg),           xfac(nshg,5),
     &          strl(numel,ngauss),     xa(nfath,3),      
     &          y(nshg,5),            yold(nshg,5),
     &          ifath(nshg),          iper(nshg),
     &          ilwork(nlwork),!        xmudmi(numel,ngauss),
     &          x(numnp,3),
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),    
     &          xnutf(nfath),          xkap(nfath),
     &          fwr2(nshg),            fwr3(nshg),
     &          xlamb1(nfath),         xlamb2(nfath),
     &          fwr4(nshg)
c
     
      fres = zero
      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)
c
c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      yold(:,5) = 1.0  ! Debugging
c      yold(:,2) = 2.0*x(:,1) - 3.0*x(:,2) 
c      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
c      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
c                               suitable for the


      intrul=intg(1,itseq)
      intind=intpt(intrul)

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)

        call asithf (yold, x, strl(iel:inum,:), mien(iblk)%p, fres, 
     &               shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo
c
 
      if(numpe>1) call commu (fres, ilwork, 24, 'in ')
c 
c account for periodicity in filtered variables
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:)
        endif
      enddo
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = fres(i,:)
        endif
      enddo
      
      if(numpe>1)   call commu (fres, ilwork, 24, 'out')

      fres(:,23) = one / fres(:,23)
      do j = 1,22
        fres(:,j) = fres(:,j) * fres(:,23)
      enddo
c
c.....at this point fres is really all of our filtered quantities
c     at the nodes
c

      xlij(:,1) = fres(:,4) - fres(:,1)*fres(:,1)
      xlij(:,2) = fres(:,5) - fres(:,2)*fres(:,2)
      xlij(:,3) = fres(:,6) - fres(:,3)*fres(:,3)
      xlij(:,4) = fres(:,7) - fres(:,1)*fres(:,2)
      xlij(:,5) = fres(:,8) - fres(:,1)*fres(:,3)
      xlij(:,6) = fres(:,9) - fres(:,2)*fres(:,3)

      strnrm = sqrt( 
     &  two * (fres(:,10)**2 + fres(:,11)**2 + fres(:,12)**2)
     &  + four * ( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )

      xfac(:,1) = strnrm*strnrm*( fres(:,10)**2 + fres(:,11)**2 + 
     &     fres(:,12)**2
     &     + two*( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )
      
      xfac(:,2) = strnrm*( xlij(:,1)*fres(:,10) + xlij(:,2)*fres(:,11) 
     &     + xlij(:,3)*fres(:,12) + 
     &     two*(xlij(:,4)*fres(:,13) + xlij(:,5)*fres(:,14) +
     &     xlij(:,6)*fres(:,15)) )

      xfac(:,3) = strnrm*( fres(:,10)*fres(:,16) + fres(:,11)*fres(:,17)
     &     + fres(:,12)*fres(:,18) + 
     &     two*(fres(:,13)*fres(:,19) + fres(:,14)*fres(:,20) +
     &     fres(:,15)*fres(:,21)) )

      xfac(:,4) = xlij(:,1)*fres(:,16) + xlij(:,2)*fres(:,17)
     &     + xlij(:,3)*fres(:,18) + 
     &     two*(xlij(:,4)*fres(:,19) + xlij(:,5)*fres(:,20) +
     &     xlij(:,6)*fres(:,21))

      xfac(:,5) = fres(:,16)*fres(:,16) + fres(:,17)*fres(:,17)
     &     + fres(:,18)*fres(:,18) + 
     &     two*(fres(:,19)*fres(:,19) + fres(:,20)*fres(:,20) +
     &     fres(:,21)*fres(:,21))


c      xfac(:,1) = one ! Debugging
c      xfac(:,2) = one 
c      xfac(:,3) = two
c      xfac(:,4) = one 
c      xfac(:,5) = one       

c  zero on processor periodic nodes so that they will not be added twice
      
      do j = 1, nshg
          i = iper(j)
          if (i .ne. j) then
            xfac(j,1) = zero
            xfac(j,2) = zero
            xfac(j,3) = zero
            xfac(j,4) = zero
            xfac(j,5) = zero
          endif
       enddo

      if (numpe.gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xfac(isgbeg:isgend,:) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif

c... Debugging

      xatm1 = sum(xfac(:,1))
      xatm2 = sum(xfac(:,2))
      xatm3 = sum(xfac(:,3))
      xatm4 = sum(xfac(:,4))      
      xatm5 = sum(xfac(:,5))


c
c Description of arrays.   Each processor has an array of length equal
c to the total number of fathers times 2 xnude(nfathers,2). One to collect 
c the numerator and one to collect the denominator.  There is also an array
c of length nshg on each processor which tells the father number of each
c on processor node, ifath(nnshg).  Finally, there is an arry of length
c nfathers to tell the total (on all processors combined) number of sons
c for each father. 
c
c  Now loop over nodes and accumlate the numerator and the denominator
c  to the father nodes.  Only on processor addition at this point.
c  Note that serrogate fathers are collect some for the case where some
c  sons are on another processor
c
      xnude = zero
      do i = 1,nshg
         xnude(ifath(i),1) = xnude(ifath(i),1) + xfac(i,1)
         xnude(ifath(i),2) = xnude(ifath(i),2) + xfac(i,2)
         xnude(ifath(i),3) = xnude(ifath(i),3) + xfac(i,3)
         xnude(ifath(i),4) = xnude(ifath(i),4) + xfac(i,4)
         xnude(ifath(i),5) = xnude(ifath(i),5) + xfac(i,5)
      enddo

c
c Now  the true fathers and serrogates combine results and update
c each other.
c       
      if(numpe .gt. 1)then
         call drvAllreduce(xnude, xnuder,5*nfath)
c
c  xnude is the sum of the sons for each father on this processor
c
c  xnuder is the sum of the sons for each father on all processor combined
c  (the same as if we had not partitioned the mesh for each processor)
c
c   For each father we have precomputed the number of sons (including
c   the sons off processor). 
c
c   Now divide by number of sons to get the average (not really necessary
c   for dynamic model since ratio will cancel nsons at each father)
c
c         xnuder(:,1) = xnuder(:,1)  / nsons(:)
c         xnuder(:,2) = xnuder(:,2)  / nsons(:)
c         xnuder(:,3) = xnuder(:,3)  / nsons(:)
c         xnuder(:,4) = xnuder(:,4)  / nsons(:)
c         xnuder(:,5) = xnuder(:,5)  / nsons(:)
c
c  the next line are the  a, b, c coefficients in the quadratic eq.
c

         do i = 1,nfath
            xa(i,1) = two*cdelsq1(i)*xnuder(i,1) + 
     &           xnuder(i,2)  
            xa(i,2) = four*cdelsq1(i)*xnuder(i,3) +
     &           xnuder(i,4)
            xa(i,3) = two*cdelsq1(i)*xnuder(i,5)

c            xa(i,1) = xnuder(ifath(i),1) + ! Debugging
c     &           xnuder(ifath(i),2)  
c            xa(i,2) = xnuder(ifath(i),3) +
c     &           xnuder(ifath(i),4)
c            xa(i,3) = xnuder(ifath(i),5)


         enddo
      else

c         xnude(:,1) = xnude(:,1)  / nsons(:)
c         xnude(:,2) = xnude(:,2)  / nsons(:)
c         xnude(:,3) = xnude(:,3)  / nsons(:)
c         xnude(:,4) = xnude(:,4)  / nsons(:)
c         xnude(:,5) = xnude(:,5)  / nsons(:)
         
         do i = 1,nfath
            xa(i,1) = two*cdelsq1(i)*xnude(i,1) + 
     &           xnude(i,2)  
            xa(i,2) = four*cdelsq1(i)*xnude(i,3) + 
     &           xnude(i,4)
            xa(i,3) = two*cdelsq1(i)*xnude(i,5)       

c            xa(i,1) = xnude(ifath(i),1) + ! Debugging
c     &           xnude(ifath(i),2)  
c            xa(i,2) = xnude(ifath(i),3) + 
c     &           xnude(ifath(i),4)
c            xa(i,3) = xnude(ifath(i),5)       

         enddo
      endif

c... Solve a*x*x - b*x + c


      do i = 1, nfath

      xdisc = xa(i,2)**2 - four*xa(i,1)*xa(i,3)

      if (xdisc .lt. zero) then
         write(*,*) '*********Warning on filter width ratio********'
      xlamb1(i) = fwr1
      xlamb2(i) = fwr1
      if (xdisc .lt. -0.5d0) then
         write(*,*) '*********Warning on filter width ratio********'
      endif
      endif

      if (xdisc .eq. zero) then
      xlamb1(i) = xa(i,2) / (two*xa(i,1))
      xlamb2(i) = xa(i,2) / (two*xa(i,1))      
      endif

      if (xdisc .gt. zero) then
      xlamb1(i)= ( xa(i,2) + sqrt( xa(i,2)**2 - four*xa(i,1)*xa(i,3) ) )
     &     / (two*xa(i,1))
      xlamb2(i)= ( xa(i,2) - sqrt( xa(i,2)**2 - four*xa(i,1)*xa(i,3) ) )
     &     / (two*xa(i,1))
      endif

      enddo

      do i = 1, nshg
         fwr2(i) = xlamb1(ifath(i))
         fwr3(i) = xlamb2(ifath(i))
      enddo

      if (myrank .eq. master) then
         do i = 1, nfath
            write(23,*)i,xlamb1(i), xlamb2(i)
         enddo
      endif
      call flush(23)
      


      do i = 1, nfath
         xkap(i) = cdelsq1(i) / xlamb2(i)
         xa(i,1) = two*xkap(i)*xnuder(i,1)
         xa(i,2) = four*xkap(i)*xnuder(i,3) - xnuder(i,2)
         xa(i,3) = two*xkap(i)*xnuder(i,5) - xnuder(i,4)

         xlamb1(i)= ( xa(i,2) + sqrt( xa(i,2)**2 - four*xa(i,1)*xa(i,3)
     &        ) )/ (two*xa(i,1))
         xlamb2(i)= ( xa(i,2) - sqrt( xa(i,2)**2 - four*xa(i,1)*xa(i,3)
     &        ) )/ (two*xa(i,1))         

      enddo


      if (myrank .eq. master) then
         do i = 1, nfath
            write(255,*)i, xlamb1(i), xlamb2(i)
         enddo 
      endif
      call flush(255)

      fwr4(:) = xlamb1(ifath(:))

      return
      end
      subroutine DFWRsfdmc (y,      shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x, fwr2, fwr3)

      use pointer_data

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwt contains the weights of the 
c                    quad. points.  



      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension fres(nshg,24),         fwr(nshg),
     &          strnrm(nshg),         cdelsq(nshg),
     &          cdelsq2(nshg),
     &          xnum(nshg),           xden(nshg),
     &          xmij(nshg,6),         xlij(nshg,6),
     &          xnude(nfath,2),        xnuder(nfath,2),
     &          ynude(nfath,6),        ynuder(nfath,6),
     &          ui(nfath,3),           snorm(nfath),
     &          uir(nfath,3),          snormr(nfath),
     &          xm(nfath,6),           xl(nfath,6)
      dimension xl1(nfath,6),          xl2(nfath,6),
     &          xl1r(nfath,6),         xl2r(nfath,6),
     &          xmr(nfath,6),          xlr(nfath,6),
     &          nsons(nshg),
     &          strl(numel,ngauss),           
     &          y(nshg,5),            yold(nshg,5),
     &          ifath(nshg),          iper(nshg),
     &          ilwork(nlwork),!        xmudmi(numel,ngauss),
     &          x(numnp,3)
      dimension shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),    
     &          xnutf(nfath),         xfac(nshg,5),
     &          fwr2(nshg),           fwr3(nshg)

      character*10 cname
      character*30 fname1, fname2, fname3, fname4, fname5, fname6,
     &             fname0
c
c
c   setup the weights for time averaging of cdelsq (now in quadfilt module)
c
      denom=max(1.0d0*(lstep),one)
      if(dtavei.lt.0) then
         wcur=one/denom
      else
         wcur=dtavei
      endif  
      whist=1.0-wcur

      if (istep .eq. 0) then
         xnd      = zero
         xmodcomp = zero
         xmcomp  = zero
         xlcomp  = zero
         xl1comp  = zero
         xl2comp  = zero
         ucomp    = zero
         scomp    = zero
      endif

     
      fres = zero
      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)
c

c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      yold(:,5) = 1.0  ! Debugging
c      yold(:,2) = 2.0*x(:,1) - 3.0*x(:,2) 
c      yold(:,2) = 2.0
c      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      yold(:,3) = 3.0
c      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
c      yold(:,4) = 4.0
c      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
c                               suitable for the



      intrul=intg(1,itseq)
      intind=intpt(intrul)

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        
        call asithf (yold, x, strl(iel:inum,:), mien(iblk)%p, fres, 
     &               shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo
c
 
      if (ngaussf .ne. ngauss) then
      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)
 
        call getstrl (yold, x,      mien(iblk)%p,  
     &               strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:))

      enddo
      endif
c
c
C must fix for abc and dynamic model
c      if(iabc==1)   !are there any axisym bc's
c     &      call rotabc(res, iBC, BC,nflow, 'in ')
c
      if(numpe>1) call commu (fres, ilwork, 24, 'in ')
c 
c account for periodicity in filtered variables
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:)
        endif
      enddo
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = fres(i,:)
        endif
      enddo

      if(numpe>1)   call commu (fres, ilwork, 24, 'out')

      fres(:,23) = one / fres(:,23)
      do j = 1,22
        fres(:,j) = fres(:,j) * fres(:,23)
      enddo
c     fres(:,24) = fres(:,24) * fres(:,23)
c
c.....at this point fres is really all of our filtered quantities
c     at the nodes
c

      strnrm = sqrt( 
     &  two * (fres(:,10)**2 + fres(:,11)**2 + fres(:,12)**2)
     &  + four * ( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )

c      fwr = fwr1 * fres(:,22) * strnrm
      fwr = fwr3 * fres(:,22) * strnrm

      xmij(:,1) = -fwr
     &             * fres(:,10) + fres(:,16)
      xmij(:,2) = -fwr
     &             * fres(:,11) + fres(:,17) 
      xmij(:,3) = -fwr
     &             * fres(:,12) + fres(:,18) 

      xmij(:,4) = -fwr * fres(:,13) + fres(:,19)
      xmij(:,5) = -fwr * fres(:,14) + fres(:,20)
      xmij(:,6) = -fwr * fres(:,15) + fres(:,21)

      fres(:,22) = one / fres(:,22)

      xlij(:,1) = fres(:,4) - fres(:,1) * fres(:,1) * fres(:,22)
      xlij(:,2) = fres(:,5) - fres(:,2) * fres(:,2) * fres(:,22)
      xlij(:,3) = fres(:,6) - fres(:,3) * fres(:,3) * fres(:,22)
      xlij(:,4) = fres(:,7) - fres(:,1) * fres(:,2) * fres(:,22)
      xlij(:,5) = fres(:,8) - fres(:,1) * fres(:,3) * fres(:,22)
      xlij(:,6) = fres(:,9) - fres(:,2) * fres(:,3) * fres(:,22)

      xnum =        xlij(:,1) * xmij(:,1) + xlij(:,2) * xmij(:,2) 
     &                                    + xlij(:,3) * xmij(:,3)
     &     + two * (xlij(:,4) * xmij(:,4) + xlij(:,5) * xmij(:,5)
     &                                    + xlij(:,6) * xmij(:,6))
      xden =        xmij(:,1) * xmij(:,1) + xmij(:,2) * xmij(:,2) 
     &                                    + xmij(:,3) * xmij(:,3)
     &     + two * (xmij(:,4) * xmij(:,4) + xmij(:,5) * xmij(:,5)
     &                                    + xmij(:,6) * xmij(:,6))
      xden = two * xden

c... For collectection of statistics on dyn. model components 

      xfac(:,1) = strnrm*strnrm*( fres(:,10)**2 + fres(:,11)**2 + 
     &     fres(:,12)**2
     &     + two*( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )
      
      xfac(:,2) = strnrm*( xlij(:,1)*fres(:,10) + xlij(:,2)*fres(:,11) 
     &     + xlij(:,3)*fres(:,12) + 
     &     two*(xlij(:,4)*fres(:,13) + xlij(:,5)*fres(:,14) +
     &     xlij(:,6)*fres(:,15)) )

      xfac(:,3) = strnrm*( fres(:,10)*fres(:,16) + fres(:,11)*fres(:,17)
     &     + fres(:,12)*fres(:,18) + 
     &     two*(fres(:,13)*fres(:,19) + fres(:,14)*fres(:,20) +
     &     fres(:,15)*fres(:,21)) )

      xfac(:,4) = xlij(:,1)*fres(:,16) + xlij(:,2)*fres(:,17)
     &     + xlij(:,3)*fres(:,18) + 
     &     two*(xlij(:,4)*fres(:,19) + xlij(:,5)*fres(:,20) +
     &     xlij(:,6)*fres(:,21))

      xfac(:,5) = fres(:,16)*fres(:,16) + fres(:,17)*fres(:,17)
     &     + fres(:,18)*fres(:,18) + 
     &     two*(fres(:,19)*fres(:,19) + fres(:,20)*fres(:,20) +
     &     fres(:,21)*fres(:,21))

c  zero on processor periodic nodes so that they will not be added twice
        do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xnum(j) = zero
            xden(j) = zero
            xfac(j,:) = zero
            xmij(j,:) = zero
            xlij(j,:) = zero
            fres(j,:) = zero
            strnrm(j) = zero
          endif
        enddo

      if (numpe.gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xnum(isgbeg:isgend) = zero
                  xden(isgbeg:isgend) = zero
                  strnrm(isgbeg:isgend) = zero
                  xfac(isgbeg:isgend,:) = zero
                  xmij(isgbeg:isgend,:) = zero
                  xlij(isgbeg:isgend,:) = zero
                  fres(isgbeg:isgend,:) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
c
c Description of arrays.   Each processor has an array of length equal
c to the total number of fathers times 2 xnude(nfathers,2). One to collect 
c the numerator and one to collect the denominator.  There is also an array
c of length nshg on each processor which tells the father number of each
c on processor node, ifath(nnshg).  Finally, there is an arry of length
c nfathers to tell the total (on all processors combined) number of sons
c for each father. 
c
c  Now loop over nodes and accumlate the numerator and the denominator
c  to the father nodes.  Only on processor addition at this point.
c  Note that serrogate fathers are collect some for the case where some
c  sons are on another processor
c
      xnude = zero
      ynude = zero
      xm    = zero
      xl    = zero
      xl1   = zero
      xl2   = zero
      ui    = zero
      snorm = zero

      do i = 1,nshg
         xnude(ifath(i),1) = xnude(ifath(i),1) + xnum(i)
         xnude(ifath(i),2) = xnude(ifath(i),2) + xden(i)

         ynude(ifath(i),1) = ynude(ifath(i),1) + xfac(i,1)
         ynude(ifath(i),2) = ynude(ifath(i),2) + xfac(i,2)
         ynude(ifath(i),3) = ynude(ifath(i),3) + xfac(i,3)
         ynude(ifath(i),4) = ynude(ifath(i),4) + xfac(i,4)
         ynude(ifath(i),5) = ynude(ifath(i),5) + xfac(i,5)

         xm(ifath(i),1) = xm(ifath(i),1) + xmij(i,1)
         xm(ifath(i),2) = xm(ifath(i),2) + xmij(i,2)
         xm(ifath(i),3) = xm(ifath(i),3) + xmij(i,3)
         xm(ifath(i),4) = xm(ifath(i),4) + xmij(i,4)
         xm(ifath(i),5) = xm(ifath(i),5) + xmij(i,5)
         xm(ifath(i),6) = xm(ifath(i),6) + xmij(i,6)

         xl(ifath(i),1) = xl(ifath(i),1) + xlij(i,1)
         xl(ifath(i),2) = xl(ifath(i),2) + xlij(i,2)
         xl(ifath(i),3) = xl(ifath(i),3) + xlij(i,3)
         xl(ifath(i),4) = xl(ifath(i),4) + xlij(i,4)
         xl(ifath(i),5) = xl(ifath(i),5) + xlij(i,5)
         xl(ifath(i),6) = xl(ifath(i),6) + xlij(i,6)         

         xl1(ifath(i),1) = xl1(ifath(i),1) + fres(i,4)
         xl1(ifath(i),2) = xl1(ifath(i),2) + fres(i,5)
         xl1(ifath(i),3) = xl1(ifath(i),3) + fres(i,6)
         xl1(ifath(i),4) = xl1(ifath(i),4) + fres(i,7)
         xl1(ifath(i),5) = xl1(ifath(i),5) + fres(i,8)
         xl1(ifath(i),6) = xl1(ifath(i),6) + fres(i,9)         

         xl2(ifath(i),1) = xl2(ifath(i),1) + fres(i,1)*fres(i,1) 
         xl2(ifath(i),2) = xl2(ifath(i),2) + fres(i,2)*fres(i,2)
         xl2(ifath(i),3) = xl2(ifath(i),3) + fres(i,3)*fres(i,3) 
         xl2(ifath(i),4) = xl2(ifath(i),4) + fres(i,1)*fres(i,2)
         xl2(ifath(i),5) = xl2(ifath(i),5) + fres(i,1)*fres(i,3)
         xl2(ifath(i),6) = xl2(ifath(i),6) + fres(i,2)*fres(i,3)

         ui(ifath(i),1) = ui(ifath(i),1) + fres(i,1)
         ui(ifath(i),2) = ui(ifath(i),2) + fres(i,2)
         ui(ifath(i),3) = ui(ifath(i),3) + fres(i,3)

         snorm(ifath(i)) = snorm(ifath(i)) + strnrm(i)

      enddo

c
c Now  the true fathers and serrogates combine results and update
c each other.
c       
      if(numpe .gt. 1)then
         call drvAllreduce(xnude, xnuder,2*nfath)
         call drvAllreduce(ynude, ynuder,6*nfath)
         call drvAllreduce(xm, xmr,6*nfath)
         call drvAllreduce(xl, xlr,6*nfath)
         call drvAllreduce(xl1, xl1r,6*nfath)
         call drvAllreduce(xl2, xl2r,6*nfath)
         call drvAllreduce(ui, uir,3*nfath)
         call drvAllreduce(snorm, snormr,nfath)

         do i = 1, nfath
            ynuder(i,6) = ( ynuder(i,4) - fwr1*ynuder(i,2) ) /
     &           ( two*ynuder(i,5) - four*fwr1*ynuder(i,3)
     &           + two*fwr1*fwr1*ynuder(i,1) )
         enddo

         cdelsq2(:) = ynuder(ifath(:),6)  ! For comparison w/ cdelsq
c
c  xnude is the sum of the sons for each father on this processor
c
c  xnuder is the sum of the sons for each father on all processor combined
c  (the same as if we had not partitioned the mesh for each processor)
c
c   For each father we have precomputed the number of sons (including
c   the sons off processor). 
c
c   Now divide by number of sons to get the average (not really necessary
c   for dynamic model since ratio will cancel nsons at each father)
c
         xnuder(:,1) = xnuder(:,1) / nsons(:)
         xnuder(:,2) = xnuder(:,2) / nsons(:)

         do m = 1, 5
         ynuder(:,m) = ynuder(:,m)/nsons(:)
         enddo
         do m = 1,6
         xmr(:,m) = xmr(:,m)/nsons(:)
         xlr(:,m) = xlr(:,m)/nsons(:)
         xl1r(:,m) = xl1r(:,m)/nsons(:)
         xl2r(:,m) = xl2r(:,m)/nsons(:)
         enddo

         uir(:,1) = uir(:,1)/nsons(:)
         uir(:,2) = uir(:,2)/nsons(:)
         uir(:,3) = uir(:,3)/nsons(:)

         snormr(:) = snormr(:)/nsons(:)
c
cc  the next line is c \Delta^2
cc
cc         xnuder(:,1) = xnuder(:,1) / (xnuder(:,2) + 1.d-09)
cc         do i = 1,nshg
cc            cdelsq(i) = xnuder(ifath(i),1)
cc         enddo

            numNden(:,1) = whist*numNden(:,1)+wcur*xnuder(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnuder(ifath(:),2)
            cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)
            
c            cdelsq(:) = xnuder(ifath(:),1)/(xnuder(ifath(:),2)+1.d-09)

            xnd(:,1) = xnd(:,1) + xnuder(:,1)
            xnd(:,2) = xnd(:,2) + xnuder(:,2)

            xmodcomp(:,1) = xmodcomp(:,1)+ynuder(:,1)
            xmodcomp(:,2) = xmodcomp(:,2)+ynuder(:,2)            
            xmodcomp(:,3) = xmodcomp(:,3)+ynuder(:,3)
            xmodcomp(:,4) = xmodcomp(:,4)+ynuder(:,4)
            xmodcomp(:,5) = xmodcomp(:,5)+ynuder(:,5)

            xmcomp(:,:) = xmcomp(:,:)+xmr(:,:)
            xlcomp(:,:) = xlcomp(:,:)+xlr(:,:)

            xl1comp(:,:) = xl1comp(:,:)+xl1r(:,:)
            xl2comp(:,:) = xl2comp(:,:)+xl2r(:,:)

            ucomp(:,:) = ucomp(:,:)+uir(:,:)
            u1 = uir(32,1)
            scomp(:)   = scomp(:)+snormr(:)

      else

         xnude(:,1) = xnude(:,1)/nsons(:)
         xnude(:,2) = xnude(:,2)/nsons(:)

         do m = 1, 5
         ynude(:,m) = ynude(:,m)/nsons(:)
         enddo
         do m = 1,6
         xm(:,m) = xm(:,m)/nsons(:)
         xl(:,m) = xl(:,m)/nsons(:)
         xl1(:,m) = xl1(:,m)/nsons(:)
         xl2(:,m) = xl2(:,m)/nsons(:)
         enddo

         ui(:,1) = ui(:,1)/nsons(:)
         ui(:,2) = ui(:,2)/nsons(:)
         ui(:,3) = ui(:,3)/nsons(:)

         snorm(:) = snorm(:)/nsons(:)

c     
c     the next line is c \Delta^2, not nu_T but we want to save the
c     memory
c     

cc         xnude(:,1) = xnude(:,1) / (xnude(:,2) + 1.d-09)
cc        do i = 1,nshg
cc            cdelsq(i) = xnude(ifath(i),1)
cc         enddo
cc      endif

         do i = 1, nfath
            ynude(i,6) = ( ynude(i,4) - fwr1*ynude(i,2) ) /
     &           ( two*ynude(i,5) - four*fwr1*ynude(i,3)
     &           + fwr1*fwr1*ynude(i,1) )
         enddo

            numNden(:,1) = whist*numNden(:,1)+wcur*xnude(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnude(ifath(:),2)

            xnd(:,1) = xnd(:,1)+xnude(:,1)
            xnd(:,2) = xnd(:,2)+xnude(:,2)

            cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)

c            cdelsq(:) = xnude(ifath(:),1)/(xnude(ifath(:),2))!+1.d-09)
            

          cdelsq2(:) = ynude(ifath(:),6)  ! For comparison w/ cdelsq

            xmodcomp(:,1) = xmodcomp(:,1)+ynude(:,1)
            xmodcomp(:,2) = xmodcomp(:,2)+ynude(:,2)            
            xmodcomp(:,3) = xmodcomp(:,3)+ynude(:,3)
            xmodcomp(:,4) = xmodcomp(:,4)+ynude(:,4)
            xmodcomp(:,5) = xmodcomp(:,5)+ynude(:,5)

            xmcomp(:,:) = xmcomp(:,:)+xm(:,:)
            xlcomp(:,:) = xlcomp(:,:)+xl(:,:)

            xl1comp(:,:) = xl1comp(:,:)+xl1(:,:)
            xl2comp(:,:) = xl2comp(:,:)+xl2(:,:)

            ucomp(:,:) = ucomp(:,:)+ui(:,:)
            scomp(:)   = scomp(:)+snorm(:)

         endif

c         do i = 1, nfath
c            xmodcomp(i,:) = xmodcomp(i,:)/nsons(i)
c            xmcomp(i,:) = xmcomp(i,:)/nsons(i)         
c            xlcomp(i,:) = xlcomp(i,:)/nsons(i)
c            xl2comp(i,:) = xl2comp(i,:)/nsons(i)         
c            xl1comp(i,:) = xl1comp(i,:)/nsons(i)
c            xnd(i,:) = xnd(i,:)/nsons(i)
c            scomp(i) = scomp(i)/nsons(i)
c            ucomp(i,:) = ucomp(i,:)/nsons(i)
c         enddo

         if (myrank .eq. master) then
            write(*,*)'istep, nstep=', istep, nstep(1)
         endif

         if ( istep .eq. (nstep(1)-1) ) then
         if ( myrank .eq. master) then

            do i = 1, nfath
            write(365,*)xmodcomp(i,1),xmodcomp(i,2),xmodcomp(i,3),
     &              xmodcomp(i,4),xmodcomp(i,5)

            write(366,*)xmcomp(i,1),xmcomp(i,2),xmcomp(i,3)
            write(367,*)xmcomp(i,4),xmcomp(i,5),xmcomp(i,6)            

            write(368,*)xlcomp(i,1),xlcomp(i,2),xlcomp(i,3)
            write(369,*)xlcomp(i,4),xlcomp(i,5),xlcomp(i,6)

            write(370,*)xl1comp(i,1),xl1comp(i,2),xl1comp(i,3)
            write(371,*)xl1comp(i,4),xl1comp(i,5),xl1comp(i,6) 

            write(372,*)xl2comp(i,1),xl2comp(i,2),xl2comp(i,3)
            write(373,*)xl2comp(i,4),xl2comp(i,5),xl2comp(i,6)

            write(374,*)xnd(i,1),xnd(i,2),scomp(i)
            write(375,*)ucomp(i,1),ucomp(i,2),ucomp(i,3) 
            enddo

            call flush(365)
            call flush(366)
            call flush(367)
            call flush(368)
            call flush(369)
            call flush(370)
            call flush(371)
            call flush(372)
            call flush(373)
            call flush(374)
            call flush(375)

c            close(852)
c            close(853)
c            close(854)

         endif
         endif

            if (myrank .eq. master) then
               write(*,*)'uit uic=', ucomp(32,1),u1
            endif

 555     format(e14.7,4(2x,e14.7))
 556     format(e14.7,5(2x,e14.7))

         


c $$$$$$$$$$$$$$$$$$$$$$$$$$$
      tmp1 =  MINVAL(cdelsq)
      tmp2 =  MAXVAL(cdelsq)
      if(numpe>1) then
         call MPI_REDUCE (tmp1, tmp3, 1,MPI_DOUBLE_PRECISION,
     &        MPI_MIN, master, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
     &        MPI_MAX, master, MPI_COMM_WORLD, ierr)
         tmp1=tmp3
         tmp2=tmp4
      endif
      if (myrank .EQ. master) then !print CDelta^2 range
         write(34,*)lstep,tmp1,tmp2
         call flush(34)
      endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$
      
      if (myrank .eq. master) then
         write(*,*) 'cdelsq=', cdelsq(1),cdelsq(2)
         write(*,*) 'cdelsq=', cdelsq2(1),cdelsq2(2)
         write(22,*) lstep, cdelsq(1)
         call flush(22)
      endif

      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         
         ngauss = nint(lcsyst)

         call scatnu (mien(iblk)%p, strl(iel:inum,:), 
     &        mxmudmi(iblk)%p,cdelsq,shp(lcsyst,1:nshl,:))
      enddo
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$  tmp1 =  MINVAL(xmudmi)
c$$$  tmp2 =  MAXVAL(xmudmi)
c$$$  if(numpe>1) then
c$$$  call MPI_REDUCE (tmp1, tmp3, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MIN, master, MPI_COMM_WORLD, ierr)
c$$$  call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MAX, master, MPI_COMM_WORLD, ierr)
c$$$      tmp1=tmp3
c$$$  tmp2=tmp4
c$$$  endif
c$$$  if (myrank .EQ. master) then
c$$$  write(35,*) lstep,tmp1,tmp2
c$$$  call flush(35)
c$$$  endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$

c
c  if flag set, write a restart file with info (reuse xmij's memory)
c
      if(irs.eq.11) then
         lstep=999
         xmij(:,1)=xnum(:)
         xmij(:,2)=xden(:)
         xmij(:,3)=cdelsq(:)
         xmij(:,5)=xlij(:,4)    !leave M_{12} in 4 and put L_{12} here
         call restar('out ',xmij,xlij) !also dump all of L_{ij} in ac
         stop
      endif
c
c  local clipping moved to scatnu with the creation of mxmudmi pointers
c
c$$$      rmu=datmat(1,2,1)
c$$$      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
c$$$      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0
c      stop !uncomment to test dmod
c


c  write out the nodal values of xnut (estimate since we don't calc strain
c  there and must use the filtered strain).
c

      if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
c
c  collect the average strain into xnude(2)
c
         xnude(:,2) = zero
         do i = 1,numnp
            xnude(ifath(i),2) = xnude(ifath(i),2) + strnrm(i)
         enddo

         if(numpe .gt. 1) then
             call drvAllreduce(xnude(:,2), xnuder(:,2),nfath)
          else
             xnuder=xnude
          endif
c     
c          nut= cdelsq    * |S|
c 
         xnutf=xnuder(:,1)*xnuder(:,2)/nsons(:)
c
c  collect the x and y coords into xnude
c
         xnude = zero
         do i = 1,numnp
            xnude(ifath(i),1) = xnude(ifath(i),1) + x(i,1)
            xnude(ifath(i),2) = xnude(ifath(i),2) + x(i,2)
         enddo

         if(numpe .gt. 1) 
     &        call drvAllreduce(xnude, xnuder,2*nfath)
         xnuder(:,1)=xnuder(:,1)/nsons(:)
         xnuder(:,2)=xnuder(:,2)/nsons(:)
c
c  xnude is the sum of the sons for each father on this processor
c
         if((myrank.eq.master)) then
            do i=1,nfath      ! cdelsq   * |S|
               write(444,*) xnuder(i,1),xnuder(i,2),xnutf(i)
            enddo
            call flush(444)
         endif
      endif

      return
      end
      subroutine DFWRwfdmc (y,      shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x,    fwr2, fwr3)

      use pointer_data

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwtf contains the weights of the 
c                    quad. points.  

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension fres(nshg,33),         fwr(nshg),
     &          strnrm(nshg),         cdelsq(nshg),
     &          cdelsq2(nshg),
     &          xnum(nshg),           xden(nshg),
     &          xmij(nshg,6),         xlij(nshg,6),
     &          xnude(nfath,2),        xnuder(nfath,2),
     &          ynude(nfath,6),        ynuder(nfath,6),
     &          ui(nfath,3),           snorm(nfath),
     &          uir(nfath,3),          snormr(nfath)
      dimension xm(nfath,6),           xl(nfath,6),
     &          xl1(nfath,6),          xl2(nfath,6),
     &          xl1r(nfath,6),         xl2r(nfath,6),
     &          xmr(nfath,6),          xlr(nfath,6),
     &          nsons(nshg),
     &          strl(numel,ngauss),           
     &          y(nshg,5),            yold(nshg,5),
     &          ifath(nshg),          iper(nshg),
     &          ilwork(nlwork),
     &          x(numnp,3),
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT),    
     &          xnutf(nfath),
     &          hfres(nshg,22),
     &          xfac(nshg,5),         fwr2(nshg),
     &          fwr3(nshg)

      real*8 u1

      character*10 cname
      character*30 fname1, fname2, fname3, fname4, fname5, fname6
c
     
c
c
c   setup the weights for time averaging of cdelsq (now in quadfilt module)
c

      denom=max(1.0d0*(lstep),one)
      if(dtavei.lt.0) then
         wcur=one/denom
      else
         wcur=dtavei
      endif  
      whist=1.0-wcur

      if (myrank .eq. master) then
         write(*,*)'istep=', istep
      endif

      if (istep .eq. 0) then
         xnd      = zero
         xmodcomp = zero
         xmcomp  = zero
         xlcomp  = zero
         xl1comp  = zero
         xl2comp  = zero
         ucomp    = zero
         scomp    = zero
      endif


      fres = zero
      hfres = zero

      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)

c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      yold(:,5) = 1.0  ! Debugging
c      yold(:,2) = 2.0*x(:,1) - 3.0*x(:,2) 
c      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
c      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
c                               suitable for the

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
c        call hfilterBB (yold, x, mien(iblk)%p, hfres, 
c     &               shglf(lcsyst,:,1:nshl,:),
c     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

        call hfilterCC (yold, x, mien(iblk)%p, hfres, 
     &               shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo

      if(numpe>1) call commu (hfres, ilwork, 22, 'in ')
c 
c... account for periodicity in filtered variables
c
      do j = 1,nshg  !    Add on-processor slave contribution to masters
        i = iper(j)
        if (i .ne. j) then
           hfres(i,:) = hfres(i,:) + hfres(j,:)
        endif
      enddo
      do j = 1,nshg ! Set on-processor slaves to be the same as masters
        i = iper(j)
        if (i .ne. j) then
           hfres(j,:) = hfres(i,:)
        endif
      enddo

c... Set off-processor slaves to be the same as their masters

      if(numpe>1)   call commu (hfres, ilwork, 22, 'out')


      hfres(:,16) = one / hfres(:,16) ! one/(volume filter kernel)

      do j = 1, 15
	hfres(:,j) = hfres(:,j) * hfres(:,16)
      enddo	    		
      do j = 17, 22
	hfres(:,j) = hfres(:,j) * hfres(:,16)
      enddo	

c... For debugging

c      hfres(:,1) = 2.0*x(:,1) - 3.0*x(:,2) 
c      hfres(:,2) = 3.0*x(:,1) + 4.0*x(:,2)
c      hfres(:,3) = 4.0*x(:,1) + x(:,2) + x(:,3)

c... Done w/ h-filtering. Begin 2h-filtering.

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        call twohfilterBB (yold, x, strl(iel:inum,:), mien(iblk)%p, 
     &               fres, hfres, shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo
c
 

      if(numpe>1) call commu (fres, ilwork, 33, 'in ')
c 
c account for periodicity in filtered variables
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:)
        endif
      enddo

      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = fres(i,:)
        endif
      enddo

      if(numpe>1)then
         call commu (fres, ilwork, 33, 'out')
      endif

      fres(:,22) = one / fres(:,22)
      do j = 1,21
        fres(:,j) = fres(:,j) * fres(:,22)
      enddo
      do j = 23,33
        fres(:,j) = fres(:,j) * fres(:,22)
      enddo

      
      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)
 
        call getstrl (yold, x,      mien(iblk)%p,  
     &               strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:))

      enddo

c
c... Obtain the hat-tilde strain rate norm at the nodes 
c

      strnrm = sqrt( 
     &  two * (fres(:,10)**2 + fres(:,11)**2 + fres(:,12)**2)
     &  + four * ( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )

c      fwr = fwr1 * strnrm

      fwr = fwr1 * fwr3 * strnrm

      xmij(:,1) = -fwr
     &             * fres(:,10) + fres(:,16)
      xmij(:,2) = -fwr
     &             * fres(:,11) + fres(:,17) 
      xmij(:,3) = -fwr
     &             * fres(:,12) + fres(:,18) 

      xmij(:,4) = -fwr * fres(:,13) + fres(:,19)
      xmij(:,5) = -fwr * fres(:,14) + fres(:,20)
      xmij(:,6) = -fwr * fres(:,15) + fres(:,21)


      xlij(:,1) = fres(:,4) - fres(:,1) * fres(:,1) 
      xlij(:,2) = fres(:,5) - fres(:,2) * fres(:,2) 
      xlij(:,3) = fres(:,6) - fres(:,3) * fres(:,3) 
      xlij(:,4) = fres(:,7) - fres(:,1) * fres(:,2) 
      xlij(:,5) = fres(:,8) - fres(:,1) * fres(:,3) 
      xlij(:,6) = fres(:,9) - fres(:,2) * fres(:,3) 

      xnum =        xlij(:,1) * xmij(:,1) + xlij(:,2) * xmij(:,2) 
     &                                    + xlij(:,3) * xmij(:,3)
     &     + two * (xlij(:,4) * xmij(:,4) + xlij(:,5) * xmij(:,5)
     &                                    + xlij(:,6) * xmij(:,6))
      xden =        xmij(:,1) * xmij(:,1) + xmij(:,2) * xmij(:,2) 
     &                                    + xmij(:,3) * xmij(:,3)
     &     + two * (xmij(:,4) * xmij(:,4) + xmij(:,5) * xmij(:,5)
     &                                    + xmij(:,6) * xmij(:,6))
      xden = two * xden

c... For collectection of statistics on dyn. model components 

      xfac(:,1) = strnrm*strnrm*( fres(:,10)**2 + fres(:,11)**2 + 
     &     fres(:,12)**2
     &     + two*( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )
      
      xfac(:,2) = strnrm*( xlij(:,1)*fres(:,10) + xlij(:,2)*fres(:,11) 
     &     + xlij(:,3)*fres(:,12) + 
     &     two*(xlij(:,4)*fres(:,13) + xlij(:,5)*fres(:,14) +
     &     xlij(:,6)*fres(:,15)) )

      xfac(:,3) = strnrm*( fres(:,10)*fres(:,16) + fres(:,11)*fres(:,17)
     &     + fres(:,12)*fres(:,18) + 
     &     two*(fres(:,13)*fres(:,19) + fres(:,14)*fres(:,20) +
     &     fres(:,15)*fres(:,21)) )

      xfac(:,4) = xlij(:,1)*fres(:,16) + xlij(:,2)*fres(:,17)
     &     + xlij(:,3)*fres(:,18) + 
     &     two*(xlij(:,4)*fres(:,19) + xlij(:,5)*fres(:,20) +
     &     xlij(:,6)*fres(:,21))

      xfac(:,5) = fres(:,16)*fres(:,16) + fres(:,17)*fres(:,17)
     &     + fres(:,18)*fres(:,18) + 
     &     two*(fres(:,19)*fres(:,19) + fres(:,20)*fres(:,20) +
     &     fres(:,21)*fres(:,21))

c  zero on processor periodic nodes so that they will not be added twice
        do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xnum(j) = zero
            xden(j) = zero
            xfac(j,:) = zero
            xmij(j,:) = zero
            xlij(j,:) = zero
            fres(j,:) = zero
            strnrm(j) = zero
          endif
        enddo

      if (numpe.gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xnum(isgbeg:isgend) = zero
                  xden(isgbeg:isgend) = zero
                  strnrm(isgbeg:isgend) = zero
                  xfac(isgbeg:isgend,:) = zero
                  xmij(isgbeg:isgend,:) = zero
                  xlij(isgbeg:isgend,:) = zero
                  fres(isgbeg:isgend,:) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
c
c Description of arrays.   Each processor has an array of length equal
c to the total number of fathers times 2 xnude(nfathers,2). One to collect 
c the numerator and one to collect the denominator.  There is also an array
c of length nshg on each processor which tells the father number of each
c on processor node, ifath(nnshg).  Finally, there is an arry of length
c nfathers to tell the total (on all processors combined) number of sons
c for each father. 
c
c  Now loop over nodes and accumlate the numerator and the denominator
c  to the father nodes.  Only on processor addition at this point.
c  Note that serrogate fathers are collect some for the case where some
c  sons are on another processor
c
      xnude = zero
      ynude = zero
      xm    = zero
      xl    = zero
      xl1   = zero
      xl2   = zero
      ui    = zero
      snorm = zero

      do i = 1,nshg
         xnude(ifath(i),1) = xnude(ifath(i),1) + xnum(i)
         xnude(ifath(i),2) = xnude(ifath(i),2) + xden(i)

         ynude(ifath(i),1) = ynude(ifath(i),1) + xfac(i,1)
         ynude(ifath(i),2) = ynude(ifath(i),2) + xfac(i,2)
         ynude(ifath(i),3) = ynude(ifath(i),3) + xfac(i,3)
         ynude(ifath(i),4) = ynude(ifath(i),4) + xfac(i,4)
         ynude(ifath(i),5) = ynude(ifath(i),5) + xfac(i,5)

         xm(ifath(i),1) = xm(ifath(i),1) + xmij(i,1)
         xm(ifath(i),2) = xm(ifath(i),2) + xmij(i,2)
         xm(ifath(i),3) = xm(ifath(i),3) + xmij(i,3)
         xm(ifath(i),4) = xm(ifath(i),4) + xmij(i,4)
         xm(ifath(i),5) = xm(ifath(i),5) + xmij(i,5)
         xm(ifath(i),6) = xm(ifath(i),6) + xmij(i,6)

         xl(ifath(i),1) = xl(ifath(i),1) + xlij(i,1)
         xl(ifath(i),2) = xl(ifath(i),2) + xlij(i,2)
         xl(ifath(i),3) = xl(ifath(i),3) + xlij(i,3)
         xl(ifath(i),4) = xl(ifath(i),4) + xlij(i,4)
         xl(ifath(i),5) = xl(ifath(i),5) + xlij(i,5)
         xl(ifath(i),6) = xl(ifath(i),6) + xlij(i,6)         

         xl1(ifath(i),1) = xl1(ifath(i),1) + fres(i,4)
         xl1(ifath(i),2) = xl1(ifath(i),2) + fres(i,5)
         xl1(ifath(i),3) = xl1(ifath(i),3) + fres(i,6)
         xl1(ifath(i),4) = xl1(ifath(i),4) + fres(i,7)
         xl1(ifath(i),5) = xl1(ifath(i),5) + fres(i,8)
         xl1(ifath(i),6) = xl1(ifath(i),6) + fres(i,9)         

         xl2(ifath(i),1) = xl2(ifath(i),1) + fres(i,1)*fres(i,1) 
         xl2(ifath(i),2) = xl2(ifath(i),2) + fres(i,2)*fres(i,2)
         xl2(ifath(i),3) = xl2(ifath(i),3) + fres(i,3)*fres(i,3) 
         xl2(ifath(i),4) = xl2(ifath(i),4) + fres(i,1)*fres(i,2)
         xl2(ifath(i),5) = xl2(ifath(i),5) + fres(i,1)*fres(i,3)
         xl2(ifath(i),6) = xl2(ifath(i),6) + fres(i,2)*fres(i,3)

         ui(ifath(i),1) = ui(ifath(i),1) + fres(i,1)
         ui(ifath(i),2) = ui(ifath(i),2) + fres(i,2)
         ui(ifath(i),3) = ui(ifath(i),3) + fres(i,3)

         snorm(ifath(i)) = snorm(ifath(i)) + strnrm(i)

      enddo

c
c Now  the true fathers and serrogates combine results and update
c each other.
c       
      if(numpe .gt. 1)then
         call drvAllreduce(xnude, xnuder,2*nfath)
         call drvAllreduce(ynude, ynuder,6*nfath)
         call drvAllreduce(xm, xmr,6*nfath)
         call drvAllreduce(xl, xlr,6*nfath)
         call drvAllreduce(xl1, xl1r,6*nfath)
         call drvAllreduce(xl2, xl2r,6*nfath)
         call drvAllreduce(ui, uir,3*nfath)
         call drvAllreduce(snorm, snormr,nfath)

         do i = 1, nfath
            ynuder(i,6) = ( ynuder(i,4) - fwr1*ynuder(i,2) ) /
     &           ( two*ynuder(i,5) - four*fwr1*ynuder(i,3)
     &           + two*fwr1*fwr1*ynuder(i,1) )
         enddo

         cdelsq2(:) = ynuder(ifath(:),6)  ! For comparison w/ cdelsq
c
c  xnude is the sum of the sons for each father on this processor
c
c  xnuder is the sum of the sons for each father on all processor combined
c  (the same as if we had not partitioned the mesh for each processor)
c
c   For each father we have precomputed the number of sons (including
c   the sons off processor). 
c
c   Now divide by number of sons to get the average (not really necessary
c   for dynamic model since ratio will cancel nsons at each father)
c
         xnuder(:,1) = xnuder(:,1) / nsons(:)
         xnuder(:,2) = xnuder(:,2) / nsons(:)

         do m = 1, 5
         ynuder(:,m) = ynuder(:,m)/nsons(:)
         enddo
         do m = 1,6
         xmr(:,m) = xmr(:,m)/nsons(:)
         xlr(:,m) = xlr(:,m)/nsons(:)
         xl1r(:,m) = xl1r(:,m)/nsons(:)
         xl2r(:,m) = xl2r(:,m)/nsons(:)
         enddo

         uir(:,1) = uir(:,1)/nsons(:)
         uir(:,2) = uir(:,2)/nsons(:)
         uir(:,3) = uir(:,3)/nsons(:)

         snormr(:) = snormr(:)/nsons(:)

c
cc  the next line is c \Delta^2
cc
cc         xnuder(:,1) = xnuder(:,1) / (xnuder(:,2) + 1.d-09)
cc         do i = 1,nshg
cc            cdelsq(i) = xnuder(ifath(i),1)
cc         enddo

            numNden(:,1) = whist*numNden(:,1)+wcur*xnuder(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnuder(ifath(:),2)
            cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)
            
c            cdelsq(:) = xnuder(ifath(:),1)/(xnuder(ifath(:),2)+1.d-09)

            xnd(:,1) = xnd(:,1) + xnuder(:,1)
            xnd(:,2) = xnd(:,2) + xnuder(:,2)

            xmodcomp(:,1) = xmodcomp(:,1)+ynuder(:,1)
            xmodcomp(:,2) = xmodcomp(:,2)+ynuder(:,2)            
            xmodcomp(:,3) = xmodcomp(:,3)+ynuder(:,3)
            xmodcomp(:,4) = xmodcomp(:,4)+ynuder(:,4)
            xmodcomp(:,5) = xmodcomp(:,5)+ynuder(:,5)

            xmcomp(:,:) = xmcomp(:,:)+xmr(:,:)
            xlcomp(:,:) = xlcomp(:,:)+xlr(:,:)

            xl1comp(:,:) = xl1comp(:,:)+xl1r(:,:)
            xl2comp(:,:) = xl2comp(:,:)+xl2r(:,:)

            ucomp(:,:) = ucomp(:,:)+uir(:,:)
            u1 = uir(32,1)
            scomp(:)   = scomp(:)+snormr(:)

      else

         xnude(:,1) = xnude(:,1)/nsons(:)
         xnude(:,2) = xnude(:,2)/nsons(:)

         do m = 1, 5
         ynude(:,m) = ynude(:,m)/nsons(:)
         enddo
         do m = 1,6
         xm(:,m) = xm(:,m)/nsons(:)
         xl(:,m) = xl(:,m)/nsons(:)
         xl1(:,m) = xl1(:,m)/nsons(:)
         xl2(:,m) = xl2(:,m)/nsons(:)
         enddo

         ui(:,1) = ui(:,1)/nsons(:)
         ui(:,2) = ui(:,2)/nsons(:)
         ui(:,3) = ui(:,3)/nsons(:)

         snorm(:) = snorm(:)/nsons(:)
c     
c     the next line is c \Delta^2, not nu_T but we want to save the
c     memory
c     

cc         xnude(:,1) = xnude(:,1) / (xnude(:,2) + 1.d-09)
cc        do i = 1,nshg
cc            cdelsq(i) = xnude(ifath(i),1)
cc         enddo
cc      endif

         do i = 1, nfath
            ynude(i,6) = ( ynude(i,4) - fwr1*ynude(i,2) ) /
     &           ( two*ynude(i,5) - four*fwr1*ynude(i,3)
     &           + fwr1*fwr1*ynude(i,1) )
         enddo

            numNden(:,1) = whist*numNden(:,1)+wcur*xnude(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnude(ifath(:),2)

            xnd(:,1) = xnd(:,1)+xnude(:,1)
            xnd(:,2) = xnd(:,2)+xnude(:,2)

            cdelsq(:) = numNden(:,1) / (numNden(:,2)) ! + 1.d-09)

c            cdelsq(:) = xnude(ifath(:),1)/(xnude(ifath(:),2))!+1.d-09)
            

          cdelsq2(:) = ynude(ifath(:),6)  ! For comparison w/ cdelsq

            xmodcomp(:,1) = xmodcomp(:,1)+ynude(:,1)
            xmodcomp(:,2) = xmodcomp(:,2)+ynude(:,2)            
            xmodcomp(:,3) = xmodcomp(:,3)+ynude(:,3)
            xmodcomp(:,4) = xmodcomp(:,4)+ynude(:,4)
            xmodcomp(:,5) = xmodcomp(:,5)+ynude(:,5)

            xmcomp(:,:) = xmcomp(:,:)+xm(:,:)
            xlcomp(:,:) = xlcomp(:,:)+xl(:,:)

            xl1comp(:,:) = xl1comp(:,:)+xl1(:,:)
            xl2comp(:,:) = xl2comp(:,:)+xl2(:,:)

            ucomp(:,:) = ucomp(:,:)+ui(:,:)
            u1 = ui(32,1)
            scomp(:)   = scomp(:)+snorm(:)

         endif


c         do i = 1, nfath
c            xmodcomp(i,:) = xmodcomp(i,:)/nsons(i)
c            xmcomp(i,:) = xmcomp(i,:)/nsons(i)         
c            xlcomp(i,:) = xlcomp(i,:)/nsons(i)
c            xl2comp(i,:) = xl2comp(i,:)/nsons(i)         
c            xl1comp(i,:) = xl1comp(i,:)/nsons(i)
c            xnd(i,:) = xnd(i,:)/nsons(i)
c            scomp(i) = scomp(i)/nsons(i)
c            ucomp(i,:) = ucomp(i,:)/nsons(i)
c         enddo

         if ( istep .eq. (nstep(1)-1) ) then
         if ( myrank .eq. master) then

            do i = 1, nfath
            write(365,*)xmodcomp(i,1),xmodcomp(i,2),xmodcomp(i,3),
     &              xmodcomp(i,4),xmodcomp(i,5)

            write(366,*)xmcomp(i,1),xmcomp(i,2),xmcomp(i,3)
            write(367,*)xmcomp(i,4),xmcomp(i,5),xmcomp(i,6)            

            write(368,*)xlcomp(i,1),xlcomp(i,2),xlcomp(i,3)
            write(369,*)xlcomp(i,4),xlcomp(i,5),xlcomp(i,6)

            write(370,*)xl1comp(i,1),xl1comp(i,2),xl1comp(i,3)
            write(371,*)xl1comp(i,4),xl1comp(i,5),xl1comp(i,6) 

            write(372,*)xl2comp(i,1),xl2comp(i,2),xl2comp(i,3)
            write(373,*)xl2comp(i,4),xl2comp(i,5),xl2comp(i,6)

            write(374,*)xnd(i,1),xnd(i,2),scomp(i)
            write(375,*)ucomp(i,1),ucomp(i,2),ucomp(i,3) 

c            write(*,*)'uit uic=', ucomp(32,1),u1
            enddo


            call flush(365)
            call flush(366)
            call flush(367)
            call flush(368)
            call flush(369)
            call flush(370)
            call flush(371)
            call flush(372)
            call flush(373)
            call flush(374)
            call flush(375)

c            if (myrank .eq. master) then
c               write(*,*)'uit uic=', ucomp(32,1),u1
c            endif


c            close(852)
c            close(853)
c            close(854)

         endif
         endif

            if (myrank .eq. master) then
               write(*,*)'uit uic=', ucomp(32,1),u1
            endif


 555     format(e14.7,4(2x,e14.7))
 556     format(e14.7,5(2x,e14.7))

c         close(849)
c         close(850)         
c         close(851)
c         close(852)  
c         close(853)
c         close(854)          

c $$$$$$$$$$$$$$$$$$$$$$$$$$$
      tmp1 =  MINVAL(cdelsq)
      tmp2 =  MAXVAL(cdelsq)
      if(numpe>1) then
         call MPI_REDUCE (tmp1, tmp3, 1,MPI_DOUBLE_PRECISION,
     &        MPI_MIN, master, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
     &        MPI_MAX, master, MPI_COMM_WORLD, ierr)
         tmp1=tmp3
         tmp2=tmp4
      endif
      if (myrank .EQ. master) then !print CDelta^2 range
         write(34,*)lstep,tmp1,tmp2
         call flush(34)
      endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$
      
      if (myrank .eq. master) then
         write(*,*) 'cdelsq=', cdelsq(1),cdelsq(2)
         write(*,*) 'cdelsq=', cdelsq2(1),cdelsq2(2)
         write(22,*) lstep, cdelsq(1)
         call flush(22)
      endif

      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         
         ngauss = nint(lcsyst)

         call scatnu (mien(iblk)%p, strl(iel:inum,:), 
     &        mxmudmi(iblk)%p,cdelsq,shp(lcsyst,1:nshl,:))
      enddo
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$  tmp1 =  MINVAL(xmudmi)
c$$$  tmp2 =  MAXVAL(xmudmi)
c$$$  if(numpe>1) then
c$$$  call MPI_REDUCE (tmp1, tmp3, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MIN, master, MPI_COMM_WORLD, ierr)
c$$$  call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
c$$$  &                 MPI_MAX, master, MPI_COMM_WORLD, ierr)
c$$$      tmp1=tmp3
c$$$  tmp2=tmp4
c$$$  endif
c$$$  if (myrank .EQ. master) then
c$$$  write(35,*) lstep,tmp1,tmp2
c$$$  call flush(35)
c$$$  endif
c $$$$$$$$$$$$$$$$$$$$$$$$$$$

c
c  if flag set, write a restart file with info (reuse xmij's memory)
c
      if(irs.eq.11) then
         lstep=999
         xmij(:,1)=xnum(:)
         xmij(:,2)=xden(:)
         xmij(:,3)=cdelsq(:)
         xmij(:,5)=xlij(:,4)    !leave M_{12} in 4 and put L_{12} here
         call restar('out ',xmij,xlij) !also dump all of L_{ij} in ac
         stop
      endif
c
c  local clipping moved to scatnu with the creation of mxmudmi pointers
c
c$$$      rmu=datmat(1,2,1)
c$$$      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
c$$$      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0
c      stop !uncomment to test dmod
c


c  write out the nodal values of xnut (estimate since we don't calc strain
c  there and must use the filtered strain).
c

      if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
c
c  collect the average strain into xnude(2)
c
         xnude(:,2) = zero
         do i = 1,numnp
            xnude(ifath(i),2) = xnude(ifath(i),2) + strnrm(i)
         enddo

         if(numpe .gt. 1) then
             call drvAllreduce(xnude(:,2), xnuder(:,2),nfath)
          else
             xnuder=xnude
          endif
c     
c          nut= cdelsq    * |S|
c 
         xnutf=xnuder(:,1)*xnuder(:,2)/nsons(:)
c
c  collect the x and y coords into xnude
c
         xnude = zero
         do i = 1,numnp
            xnude(ifath(i),1) = xnude(ifath(i),1) + x(i,1)
            xnude(ifath(i),2) = xnude(ifath(i),2) + x(i,2)
         enddo

         if(numpe .gt. 1) 
     &        call drvAllreduce(xnude, xnuder,2*nfath)
         xnuder(:,1)=xnuder(:,1)/nsons(:)
         xnuder(:,2)=xnuder(:,2)/nsons(:)
c
c  xnude is the sum of the sons for each father on this processor
c
         if((myrank.eq.master)) then
            do i=1,nfath      ! cdelsq   * |S|
               write(444,*) xnuder(i,1),xnuder(i,2),xnutf(i)
            enddo
            call flush(444)
         endif
      endif

      return
      end
