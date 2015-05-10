      subroutine solvecon(y,       x,      iBC,  BC, 
     &                     iper,    ilwork, shp,  shgl)

c---------------------------------------------------------------------
c This subroutine is to calculate the constarint for the redistancing 
c scalar of the level set method. This is to prevent interface from 
c moving by applying the condition that the volume must stay constant 
c in each element when the redisatnce step is applied.
c--------------------------------------------------------------------
c
c
      use pointer_data
c     
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"      
c      
      dimension y(nshg,ndof),                    
     &          x(numnp,nsd),            iBC(nshg),
     &          BC(nshg,ndofBC),         ilwork(nlwork),
     &          iper(nshg)
c
c     
      dimension shp(MAXTOP,maxsh,MAXQPT),  
     &           shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &           v_lambda(nshg),   hprime(nshg),
     &           v_lambda1(nshg), v_lambda2(nshg),
     &           rmass(nshg) 
c
        real*8, allocatable :: tmpshp(:,:),  tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

c
c ... intialize
c       
      rmass  = zero
      v_lambda = zero
      v_lambda1 = zero
      v_lambda2 = zero
      hprime = zero
c
c ... loop over element blocks
c
      do iblk = 1, nelblk
c
c.... set up the parameters
c
         nenl   = lcblk(5,iblk) ! no. of vertices per element
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
c.... compute and assemble the constarint factor, and mass matrix
c     

         allocate (tmpshp(nshl,MAXQPT))
         allocate (tmpshgl(nsd,nshl,MAXQPT))

         tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
         tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
c   
         call volcon (y,          x,             tmpshp,              
     &                tmpshgl,    mien(iblk)%p,  rmass,     
     &                v_lambda1,  hprime,        v_lambda2)

         deallocate ( tmpshp )
         deallocate ( tmpshgl ) 
      enddo

c
c ... multiple processor communication
c
      if (numpe > 1) then
         call commu (v_lambda1  , ilwork, 1  , 'in ')
         call commu (v_lambda2  , ilwork, 1  , 'in ')
         call commu (hprime     , ilwork, 1  , 'in ')
         call commu (rmass      , ilwork, 1  , 'in ')
      endif
c
c.... take care of periodic boundary conditions
c
      do j= 1,nshg
         if (btest(iBC(j),10)) then
            i = iper(j)
            rmass(i) = rmass(i) + rmass(j)
            v_lambda1(i) =  v_lambda1(i) +  v_lambda1(j)
            v_lambda2(i) =  v_lambda2(i) +  v_lambda2(j)
            hprime(i)   =  hprime(i) + hprime(j)
         endif
      enddo
c
      do j= 1,nshg
         if (btest(iBC(j),10)) then
            i = iper(j)
            rmass(j) = rmass(i)
            v_lambda1(j) =  v_lambda1(i)
            v_lambda2(j) =  v_lambda2(i)
            hprime(j) = hprime(i)
         endif
      enddo
c
c ... calculation of constraint factor
c
      rmass  = one/rmass
      v_lambda1 = v_lambda1*rmass ! numerator of lambda
      v_lambda2 = v_lambda2*rmass ! denominator of lambda
      v_lambda  = v_lambda1/(v_lambda2+epsM**2)
      hprime    = hprime*rmass
      v_lambda  = v_lambda*hprime  
c
c ... commu out for the multiple processor
c
      if(numpe > 1) then
         call commu (v_lambda, ilwork, 1, 'out')    
      endif
c          
c ... the following commented lines are for the different way of getting 
c     the denominator of constraint (lambda) calculation 
c$$$                hprime=zero
c$$$                do kk=1, nshg
c$$$                   if (abs (y(kk,6)) .le. epsilon_ls) then
c$$$                      hprime(kk) = (0.5/epsilon_ls) * (1 
c$$$     &                   + cos(pi*y(kk,6)/epsilon_ls))
c$$$                   endif
c$$$                enddo
c$$$                y(:,7) = y(:,7)+v_lambda*hprime/dtgl
c
c ... the vlome constraint applied on the second scalar
c
      y(:,7) = y(:,7)+v_lambda/dtgl   
c  
      return
      end
c
c
c

      subroutine volcon (y,         x,      shp,      
     &                   shgl,      ien,    rmass, 
     &                   v_lambda1, hprime, v_lambda2)

c---------------------------------------------------------------------
c
c This subroutine is to calculate the element contribution to the 
c constraint factor and mass matrix.
c
c---------------------------------------------------------------------
      include "common.h"
c     
      dimension y(nshg,ndof),               x(numnp,nsd),              
     &            shp(nshl,maxsh),  
     &            shgl(nsd,nshl,maxsh),
     &            ien(npro,nshl),
     &            qres(nshg,idflx),         rmass(nshg)
c
c.... element level declarations
c
      dimension ycl(npro,nshl,ndof),      xl(npro,nenl,nsd),         
     &          rmassl(npro,nshl)     
      dimension sgn(npro,nshape),         v_lambdal1(npro,nshl),
     &          v_lambda1(nshg),          hprimel(npro,nshl),
     &          hprime(nshg),             v_lambdal2(npro,nshl),
     &          v_lambda2(nshg)
c
c local arrays
c
      dimension shg(npro,nshl,nsd),
     &          dxidx(npro,nsd,nsd),      WdetJ(npro)
c
      dimension shape(npro,nshl),
     &          shdrv(npro,nsd,nshl)
c
c
c.... for volume constraint calculation of redistancing step
c
      dimension Sclr(npro),              Sclrtmp(npro),
     &          h_prime(npro),           tmp1(npro), 
     &          tmp2(npro),
     &          v_lambdatmp(npro),       v_lambdal(npro,nshl)
c$$$     &          ,hprimel(npro,nshl),      v_lambdal1(npro,nshl),
c$$$     &          v_lambdal2(npro,nshl)
c above arrays must be uncommented for alternate method included below (commented)
      real epsilon_tmp

c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
      if (ipord .gt. 1) then
         call getsgn(ien,sgn)
      endif
c
c.... gather the variables
c

      call localy(y,      ycl,     ien,    ndof,   'gather  ')
      call localx(x,      xl,      ien,    nsd,    'gather  ')
c
c.... get the element contributions of the numerator and denominator 
c     of lambda 
c
      rmassl = zero
      v_lambdal1= zero
      v_lambdal2= zero
      hprimel = zero
c
c.... loop through the integration points
c
      do intp = 1, ngauss
         if (Qwt(lcsyst,intp) .eq. zero) cycle ! precaution
c
c.... create a matrix of shape functions (and derivatives) for each
c     element at this quadrature point. These arrays will contain 
c     the correct signs for the hierarchic basis
c
         call getshp(shp,          shgl,      sgn, 
     &               shape,        shdrv)
c
c.... initialize
c     
         h_prime   = zero        
	 sclr      = zero
	 sclrtmp   = zero
         v_lambdatmp=zero
         tmp1       =zero
         tmp2       =zero

c
c.... --------------------->  Element Metrics  <-----------------------
c
         call e3metric( xl,         shdrv,        dxidx,  
     &                  shg,        WdetJ)

c
         do i = 1, nshl 
c
c  y(intp)=SUM_{a=1}^nshl (N_a(intp) Ya)
c     
            Sclr    = Sclr    + shape(:,i) * ycl(:,i,7) !d^kbar
            sclrtmp = sclrtmp + shape(:,i) * ycl(:,i,6) !d^0
         enddo

         if (isclr .eq. 2) then
            epsilon_tmp = epsilon_lsd
         else
            epsilon_tmp = epsilon_ls
         endif

         do i=1,npro
            if (abs (Sclrtmp(i)) .le. epsilon_tmp) then
               h_prime(i) = (0.5/epsilon_tmp) * (1 
     &                    + cos(pi*Sclrtmp(i)/epsilon_tmp))
               tmp1(i)=-h_prime(i)*(sclr(i)-sclrtmp(i))*dtgl
               tmp2(i)=h_prime(i)**2
c              v_lambdatmp(i)=tmp1(i)/(tmp2(i)+ epsM)
            endif
         enddo
c
         do i=1,nshl
            v_lambdal1(:,i)=v_lambdal1(:,i)+ shape(:,i)*WdetJ*tmp1
            v_lambdal2(:,i)=v_lambdal2(:,i)+ shape(:,i)*WdetJ*tmp2
c           v_lambdal(:,i)=v_lambdal(:,i)+ shape(:,i)*WdetJ*v_lambdatmp
            hprimel(:,i) = hprimel(:,i) + shape(:,i)*WdetJ*h_prime
            rmassl(:,i)  =rmassl(:,i) + shape(:,i)*WdetJ
         enddo
c
c.... end of the loop over integration points
c
      enddo
c
      call local (v_lambda1,  v_lambdal1, ien,  1,  'scatter ')
      call local (v_lambda2,  v_lambdal2, ien,  1,  'scatter ')
      call local (hprime,     hprimel,    ien,  1,  'scatter ')
      call local (rmass,      rmassl,     ien,  1,  'scatter ') 
c
      return
      end

