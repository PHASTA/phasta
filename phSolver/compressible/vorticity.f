        subroutine vortGLB  (y, x, shp, shgl, ilwork, vortG)

        use pointer_data
        include "common.h"
        include "mpif.h"
c
        real*8  y(nshg,ndof),  x(numnp,nsd)              
c
        real*8  shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT) 

        integer ilwork(nlwork)
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8  vortG(nshg,4) ! The first three components are vorticity vector, the 4th one is Q criterion 
        real*8  vortIntG(nshg,4)
        real*8  lpmassG(nshg)
        real*8  tmp(nshg)     

        vortG=0
        lpmassG=0
        vortIntG=0

        do iblk = 1, nelblk
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
c
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          if(myrank.eq.0)then
            if(mod(iblk,100).eq.0)then
               write(*,*)'Compute vortIntG:', real(iblk)/nelblk
            endif
          endif
        call vortELE(y,x,tmpshp,tmpshgl,mien(iblk)%p,vortIntG,lpmassG)
          
          deallocate(tmpshp)
          deallocate(tmpshgl)
        enddo

         if(myrank.eq.0)write(*,*)'Communicating...'
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)         
         call commu(lpmassG,ilwork,1,'in ')
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         do i=1,4
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
           tmp(:)=vortIntG(:,i)
           call commu(tmp,ilwork,1,'in ')
           vortIntG(:,i)=tmp(:)
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)           
         enddo
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         do i=1,4
            vortG(:,i)=vortIntG(:,i)/lpmassG(:)  
         enddo

!         call write_field(myrank,'a'//char(0),'vortG'//char(0),5,
!     &                    vortG, 'd'//char(0), nshg,4,lstep)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine vortELE(y,x,shp1,shgl1,ien,vortIntG,lpmassG)

  
       include "common.h"

       real*8 y(nshg,ndof), x(numnp,nsd)
       real*8 shp1(nshl,MAXQPT), shgl1(nsd,nshl,MAXQPT)
       real*8 shp2(nshl,ngauss), shgl2(nsd,nshl,ngauss) 
       real*8 shp(npro,nshl), shgl(npro,nsd,nshl)
       real*8 sgn(npro,nshl)
       integer ien(npro,nshl)
       real*8 yl(npro,nshl,ndof)
       real*8 xl(npro,nenl,nsd)
       real*8 dxidx(npro,nsd,nsd) 
       real*8 WdetJ(npro)  
       real*8 vortIntL(npro,nshl,4)
       real*8 lpmass(npro,nshl)  
       real*8 volum(npro)
       real*8 alph2(npro)
       real*8 vortIntG(nshg,4)
       real*8 lpmassG(nshg)
       real*8 shg(npro,nshl,nsd)
       real*8 g1yi(npro,nflow)
       real*8 g2yi(npro,nflow)
       real*8 g3yi(npro,nflow)     

       call getsgn(ien,sgn)
       shp2(:,1:ngauss)=shp1(:,1:ngauss)
       shgl2(:,:,1:ngauss)=shgl1(:,:,1:ngauss)
  
       call localy (y, yl, ien, ndof,  'gather  ')
       call localx (x, xl, ien, nsd,   'gather  ')

      
       vortIntL=0
       lpmass=0
       volum=0
       alph2=0


       do intp=1,ngauss
                
         call getshp (shp2, shgl2, sgn, shp, shgl )
    
         call e3metric( xl,         shgl,        dxidx,  
     &                shg,        WdetJ) 
      
         g1yi = zero
         g2yi = zero
         g3yi = zero
       
       
         do n=1,nshl
           g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * yl(:,n,2) ! du/dx
           g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * yl(:,n,2) ! du/dy
           g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * yl(:,n,2) ! du/dz

           g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * yl(:,n,3) ! dv/dx
           g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * yl(:,n,3) ! dv/dy
           g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * yl(:,n,3) ! dv/dz

           g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * yl(:,n,4) ! dw/dx
           g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * yl(:,n,4) ! dw/dy       
           g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * yl(:,n,4) ! dw/dz

           alph2(:) = alph2(:) + shp(:,n)*shp(:,n)*WdetJ(:)
         
         enddo
         volum(:)  = volum(:) + WdetJ(:)
 
c--------- compute vorticity and Q criterion     
         do n=1,nshl

            lpmass(:,n) = lpmass(:,n) + WdetJ(:)*shp(:,n)*shp(:,n)

            vortIntL(:,n,1) = vortIntL(:,n,1) +
     &              (g2yi(:,4)-g3yi(:,3)) !dw/dy-dv/dz
     &               *WdetJ(:)*shp(:,n)

            vortIntL(:,n,2) = vortIntL(:,n,2) +
     &              (g3yi(:,2)-g1yi(:,4)) !du/dz-dw/dx
     &               *WdetJ(:)*shp(:,n)

            vortIntL(:,n,3) = vortIntL(:,n,3) +
     &              (g1yi(:,3)-g2yi(:,2)) !dv/dx-du/dy
     &               *WdetJ(:)*shp(:,n)


            vortIntL(:,n,4) = vortIntL(:,n,4) + 
     &              ( 0.5*(g2yi(:,4)-g3yi(:,3))**2 !dw/dy-dv/dz
     &               +0.5*(g3yi(:,2)-g1yi(:,4))**2 !du/dz-dw/dx
     &               +0.5*(g1yi(:,3)-g2yi(:,2))**2 !dv/dx-du/dy
     &               -g1yi(:,2)**2-g2yi(:,3)**2-g3yi(:,4)**2
     &               -0.5*(g2yi(:,4)+g3yi(:,3))**2
     &               -0.5*(g3yi(:,2)+g1yi(:,4))**2
     &               -0.5*(g1yi(:,3)+g2yi(:,2))**2 )*0.5   
     &               *WdetJ(:)*shp(:,n) 
 
         enddo 
  
       enddo  ! end of loop over Gaussian points

       do n=1,nshl
          lpmass(:,n)=lpmass(:,n)*volum(:)/alph2(:)
       enddo  
        
       call local(vortIntG,vortIntL,ien,4,'scatter ')
       call local(lpmassG,lpmass,ien,1,'scatter ')
        
       return  
       end  

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
       subroutine PRcalcDuct  (y)

        include "common.h"
        include "mpif.h"

       real*8 y(nshg,ndof)  
       real*8 Mach(nshg),PresRec(nshg),PresCoeff(nshg),velocity(nshg,3)
       real*8 Ptotal_contra, Pstatic_throat
       open(911,file='PresConst.inp')
       read(911,*) Ptotal_contra
       read(911,*) Pstatic_throat
       close(911)
       write(*,*)'Ptotal_contra=',Ptotal_contra
       write(*,*)'Pstatic_throat=',Pstatic_throat
c in y the variables are u v w p T
       Mach(:)=sqrt(y(:,1)**2+y(:,2)**2+y(:,3)**2)/sqrt(1.4*287*y(:,5)) 
       PresRec(:)=y(:,4)*(1+0.2*Mach(:)**2)**3.5/Ptotal_contra 
       PresCoeff(:)=(y(:,4)-Pstatic_throat)/(Ptotal_contra-Pstatic_throat) 
       velocity(:,1:3)=y(:,1:3)
       write(*,*)'writing Mach, PR and Cp'
       call write_field(myrank,'a'//char(0),'M'//char(0),1,Mach,
     &                         'd'//char(0),nshg,1,lstep)
       call write_field(myrank,'a'//char(0),'PR'//char(0),2,PresRec,
     &                         'd'//char(0),nshg,1,lstep)
       call write_field(myrank,'a'//char(0),'Cp'//char(0),2,PresCoeff,
     &                         'd'//char(0),nshg,1,lstep)
       call write_field(myrank,'a'//char(0),'V'//char(0),1,velocity,
     &                         'd'//char(0),nshg,3,lstep) 
      end
