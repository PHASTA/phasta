      subroutine Bflux ( y,          ac,        u,      x,
     &                   shp,       shgl,       shpb,   
     &                   shglb,     ilwork,     iBC,
     &                   BC,        iper  ,     wallssVec)
c
c----------------------------------------------------------------------
c
c This routine :
c   1. computes the boundary fluxes
c   2. prints the results in the file [FLUX.lstep]
c
c output:
c  in file flux.<lstep>.n (similar to restart file):
c     machin  nshg  lstep 
c     normal_1 ... normal_nsd            ! outward normal direction
c     tau_1n   ... tau_nsd n             ! boundary viscous flux
c
c Zdenek Johan, Summer 1991.
c----------------------------------------------------------------------
c
      
      use pointer_data
      
      include "common.h"
      include "mpif.h"

      character*10  cname2
      
      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          u(nshg,nsd),              x(numnp,nsd)
      dimension iBC(nshg),           
     &            BC(nshg,ndofBC),  
     &            iper(nshg)
     
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c    
c
      real*8    flxres(nshg,nflow),
     &          flxLHS(nshg,1),           flxnrm(nshg,nsd),
     &          temp(nshg),               rtmp(nshg,ndof),
     &          flxTot(nflow),            wallssVec(nshg,3)

      real*8    qres(nshg,nsd*nsd)

c
      integer   ilwork(nlwork),
     &          invflx(nshg),             nodflx(nshg)             
c
      character*60 fname1,  fmt1, fmt2, fnamer, fname2
      character*13 fieldbflux 
      character*19 fieldwss 
      integer irstin, isize, nitems, ndofwss
      integer iarray(50)  ! integers read from headers

      real*8, allocatable, dimension(:,:,:,:) :: xKebe, xGoC
      integer, allocatable, dimension(:,:)    :: ien2
      integer, allocatable, dimension(:)      :: map
      real*8, allocatable, dimension(:,:)     :: xmu2
c
c....  calculate the flux nodes
c
      numflx  = 0
      invflx  = 0
      nodflx  = 0
      do iblk = 1, nelblb 
         iel    = lcblkb(1,iblk)
         lcsyst = lcblkb(3,iblk)
         nenl   = lcblkb(5,iblk) 
         nenbl  = lcblkb(6,iblk) 
         ndofl  = lcblkb(8,iblk)
         nshl   = lcblkb(9,iblk)
         nshlb  = lcblkb(10,iblk)
         npro   = lcblkb(1,iblk+1) - iel 
         call flxNode (mienb(iblk)%p,   miBCB(iblk)%p,  invflx) 
      enddo 
c      do i = 1, nshg
      do i = 1, numnp
         if (invflx(i) .ne. 0) then
            numflx = numflx + 1
            nodflx(numflx) = i
         endif
      enddo
c     
c.... -------------------->   interior elements   <--------------------
c     
c.... initialize the arrays
c     
      flxres = zero
      flxLHS = zero
      flxnrm = zero
      if (numflx .ne. 0)  then !we have flux nodes
      qres   = zero
c     
c.... loop over the element-blocks
c
      lhs    = 0

      ires=2  ! shield e3ql from an unmapped lmassinv
      ierrcalcsave=ierrcalc
      ierrcalc=0
      do iblk = 1, nelblk
c     
c.... set up the parameters
c     
         iel    = lcblk(1,iblk)
         nenl   = lcblk(5,iblk) ! no. of vertices per element
         nshl   = lcblk(10,iblk)
         ndofl  = lcblk(8,iblk)
         lcsyst = lcblk(3,iblk)
         npro   = lcblk(1,iblk+1) - iel 
         ngauss = nint(lcsyst)       
         allocate ( ien2(npro,nshl) )
         allocate ( xmu2(npro,maxsh))
         allocate ( map(npro) )
c
c.... get the elements touching the boundary
c         
         call mapConn( mien(iblk)%p,    ien2,    invflx,
     &                 map,             nshl,    npro,    
     &                 npro2,           nshg )

         nprold = npro
         npro = npro2
         
         if (npro .ne. 0) then

            call mapArray( mxmudmi(iblk)%p, xmu2,    map,
     &                     maxsh,           nprold)
c
c.... allocate the element matrices (though they're not needed)
c
            allocate ( xKebe(npro,9,nshl,nshl) )
            allocate ( xGoC (npro,4,nshl,nshl) )
c     
c.... compute and assemble the residuals
c     
            call AsIGMR (y,                    ac,
     &                   x,                    xmu2(1:npro,:),
     &                   shp(lcsyst,1:nshl,:),
     &                   shgl(lcsyst,:,1:nshl,:),
     &                   ien2(1:npro,:),       
     &                   flxres,               qres,
     &                   xKebe,                xGoC,
     &                   rtmp)
c     
            deallocate ( xKebe )
            deallocate ( xGoC  )
         endif
         deallocate ( ien2  )
         deallocate ( xmu2  )
         deallocate ( map   )
c     
      enddo
      ierrcalc=ierrcalcsave
c     
c.... -------------------->   boundary elements   <--------------------
c     
      do iblk = 1, nelblb
c     
c.... set up the parameters
c
         iel    = lcblkb(1,iblk)
         lcsyst = lcblkb(3,iblk)
         nenl   = lcblkb(5,iblk)
         nshl   = lcblkb(9,iblk)
         nenbl  = lcblkb(6,iblk)
         nshlb  = lcblkb(10,iblk)
         npro   = lcblkb(1,iblk+1) - iel 
 
         if(lcsyst.eq.3) lcsyst=nenbl
c     
         if(lcsyst.eq.3 .or. lcsyst.eq.4) then
            ngaussb = nintb(lcsyst)
         else
            ngaussb = nintb(lcsyst)
         endif
c
c.... allocate the element matrices (though they're not needed)
c
            allocate ( xKebe(npro,9,nshl,nshl) )   
c.... compute and assemble the residuals
c
         call AsBFlx (u,                      y,
     &                ac,                     x,
     &                shpb(lcsyst,1:nshl,:),
     &                shglb(lcsyst,:,1:nshl,:),
     &                mienb(iblk)%p,
     &                miBCB(iblk)%p,           mBCB(iblk)%p,
     &                invflx,                  flxres,
     &                flxLHS,                  flxnrm,
     &                xKebe  )
c     
            deallocate ( xKebe )
c     
c.... end of boundary element loop
c
      enddo
c.... Communication needed before we take care of periodicity and
c     division of RHS by LHS ???
c
      if(numpe > 1) then
         call commu (flxres, ilwork, nflow, 'in ')
         call commu (flxLHS, ilwork, 1   , 'in ')
         call commu (flxnrm, ilwork, nsd , 'in ')
      endif
c
c  take care of periodic boundary conditions
c
        do j= 1,nshg
          if ((btest(iBC(j),10))) then
             i = iper(j)
             flxLHS(i,1) = flxLHS(i,1) + flxLHS(j,1)
             flxres(i,:) =  flxres(i,:) + flxres(j,:)
          endif
        enddo

        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            flxLHS(j,1) = flxLHS(i,1)
            flxres(j,:) = flxres(i,:)
          endif
        enddo
c
c        call bc3per(iBC,  flxres, iper, ilwork, nflow)

c
c.... integrated fluxes (aerodynamic forces update)
c
        flxTot = zero
        do n = 1, numflx
           flxTot = flxTot + flxres(nodflx(n),:)
        enddo
        Force(1) = flxTot(1)
        Force(2) = flxTot(2)
        Force(3) = flxTot(3)
      else
c       need to call commu for procs. with no flux nodes
        if(numpe > 1) then
           call commu (flxres, ilwork, nflow, 'in ')
           call commu (flxLHS, ilwork, 1   , 'in ')
           call commu (flxnrm, ilwork, nsd , 'in ')
        endif
      endif  ! make sure the zero numflux processors still commu
c
c.... only need to commu if we are going to print surface flux (or compute avg) since
c     the force calculation just sums flxres (and each on processor node
c     has his "piece" of the sum already).
c
      ntoutv=ntout
      if ((irs .ge. 1) .and. (mod(lstep, ntoutv) .eq. 0) .or. 
     &     istep .eq. nstep(itseq) .or. ioybar == 1) then


c
c  need to zero the slaves to prevent counting twice
c  (actually unnecessary since flxres of boundary nodes will be counted n
c  times while flxlhs will be counted n times-> the ratio is still
c  correct
c      

      ndofwss = 3;
      wallssVec=rtmp(:,1:ndofwss)

      if (numflx .eq. 0) then   !no flux nodes
         rtmp=zero
         wallssVec = zero

c        need to call commu for procs. with no flux nodes 
         if(numpe > 1) then
            call commu (flxres, ilwork, nflow, 'out')
            call commu (flxLHS, ilwork, 1   , 'out')
            call commu (flxnrm, ilwork, nsd , 'out')
         endif
      else
c     
c.... ---------------------------->  Solve  <---------------------------
c
c.... compute the viscous and heat fluxes
c     
c
c.... ---------------------------->  Print  <---------------------------
c
c.... nodal fluxes
c
      do i = 1, 3
         where ( (invflx .ne. 0) .and. (flxLHS(:,1) .ne. zero) )
     &        flxres(:,i) = flxres(:,i) / flxLHS(:,1)
      enddo
c     
c.... normalize the outward normal
c     
      temp = sqrt(flxnrm(:,1)**2 + flxnrm(:,2)**2 + flxnrm(:,3)**2)
      where ( (invflx .ne. 0) .and. (temp .ne. zero) )
         flxnrm(:,1) = flxnrm(:,1) / temp
         flxnrm(:,2) = flxnrm(:,2) / temp
         flxnrm(:,3) = flxnrm(:,3) / temp
      endwhere
         
c     
c.... ---------------------------->  Communications <-------------------
c
      if(numpe > 1) then
         call commu (flxres, ilwork, nflow, 'out')
         call commu (flxLHS, ilwork, 1   , 'out')
         call commu (flxnrm, ilwork, nsd , 'out')
      endif
c

         rtmp = zero
         wallssVec  = zero

         do i=1, numnp
            if (invflx(i) .ne. 0) then
               rtmp(i,2:4) = flxres(i,1:3) !viscous flux
c     calculate the WSS
               tn=    flxres(i,1) * flxnrm(i,1)
     &              + flxres(i,2) * flxnrm(i,2)
     &              + flxres(i,3) * flxnrm(i,3)

               wallssVec(i,1) = flxres(i,1) - tn * flxnrm(i,1)
               wallssVec(i,2) = flxres(i,2) - tn * flxnrm(i,2)
               wallssVec(i,3) = flxres(i,3) - tn * flxnrm(i,3)

            endif
         enddo
      endif
cc      itmp = 1
cc      if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
cc      write (fmt1,"('(''flux.'',i',i1,',1x)')") itmp
cc      write (fname1,fmt1) lstep
      
cc         fname1 = trim(fname1) // cname2(myrank+1)

cc      open (unit=iflux, file=fname1, status='unknown', 
cc     &         form='formatted',err=997)

c      write (iflux) machin, nshg, lstep
c      write (iflux) rtmp(:,1:6)
c
c.... output the results
c     
cc         do n = 1, numflx
cc            k = nodflx(n)
cc            write (iflux,2000) k,invflx(k),flxLHS(k,1),(x(k,i),i=1,3), 
cc     &           (flxnrm(k,i),  i=1,3),
cc     &           (flxres(k,i),  i=1,3)
cc         enddo
cc         close (iflux)
c... output the results in the new format in restart.step#.proc# file

cc         itmp = 1
cc         if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
cc         write (fmt2,"('(''restart.'',i',i1,',1x)')") itmp
cc         write (fname2,fmt2) lstep

cc         fname2 = trim(fname2) // cname2(myrank+1)
c
c.... open input files
c
cc         call openfile(  fname2,  'append?', irstin )
         
cc         fieldbflux = 'boundary flux' 
cc         isize = nshg*ndof
cc         nitems = 3;
cc         iarray(1) = nshg
cc         iarray(2) = ndof
cc         iarray(3) = lstep
cc         call writeheader(irstin, fieldbflux,iarray, nitems, isize, 
cc     &        'double', iotype )
    
c         fieldbflux = 'boundary flux'        
cc         nitems = nshg*ndof
cc         call writedatablock(irstin, fieldbflux,rtmp, nitems, 
cc     &        'double', iotype)
        
cc         call closefile( irstin, "append" )
c         call Write_boundaryflux(myrank,lstep,nshg,ndof,rtmp(:,1:ndof))

c     wallss vectors into the restart file(s)

        ntoutv=ntout
        if ((irs .ge. 1) .and. (mod(lstep, ntoutv) .eq. 0) .or. 
     &       istep .eq. nstep(itseq) ) then

            call write_field(myrank,'a','wss',3,
     &                       wallssVec,'d',nshg,ndofwss,lstep)
         endif

      endif
c     
      return
c
c.... file error handling
c
997     call error ('bflux   ','opening ', iflux)
c
c$$$1000    format(' ',a80,/,1x,i10,1p,3e20.7)
 2000   format(i6,i6,10(2x,E12.5e2))
c$$$2001    format(1p,1x,i6,3e15.7)
c
c.... end
c
        end


      subroutine flxNode(ienb, iBCB, flg)
c---------------------------------------------------------------------
c
c     This routine flags the flux nodes
c
c----------------------------------------------------------------------
      include "common.h"
c
      integer   flg(nshg),        iBCB(npro,ndiBCB),     
     &          ienb(npro, nshl), lnode(27)

c
c.... compute the nodes which lie on the boundary (hierarchic)
c
      call getbnodes(lnode)

      do i=1, npro 
         if (nsrflist(iBCB(i,2)).eq.1) then
            do j=1, nshlb
               nodelcl = lnode(j)
               flg(abs(ienb(i,nodelcl)))=flg(abs(ienb(i,nodelcl)))+1  
            enddo
         endif
      enddo
c
      return
      end

      
      subroutine mapConn( ien,      ien2,    mask,
     &                    map,      nshl,    npro,    
     &                    npro2,    nshg )
c-----------------------------------------------------------------------
c
c  Create a condensed connectivity array based on the nodes in
c  mask.
c
c-----------------------------------------------------------------------
      
      integer ien(npro,nshl),  ien2(npro,nshl), mask(nshg),
     &        map(npro)

      integer nshl, nshg, npro, npro2, i, iel

c
c.... first build the map
c      
      map = 0
      do i = 1, nshl
         do iel = 1, npro
            map(iel) = map(iel) + mask( abs(ien(iel,i)) )
         enddo
      enddo
      
      npro2 = 0
      do iel = 1, npro
         if ( map(iel) .gt. 0 ) then
            npro2 = npro2 + 1
            map(iel) = npro2
         else
            map(iel) = npro
         endif
      enddo
c
c.... create the condensed connectivity array
c
      if ( npro2 .gt. 0 ) then
         do i = 1, nshl
            do iel = 1, npro
               ien2(map(iel),i) = ien(iel,i)
            enddo
         enddo
      endif
      
      return 
      end
         
      
      subroutine mapArray( x,      x2,      map,
     &                     nshl,   nprold)
c-----------------------------------------------------------------------
c
c  Maps array x into array x2 based on the given map
c
c-----------------------------------------------------------------------
      real*8   x(nprold,nshl),    x2(nprold,nshl)
      
      integer  map(nprold)

      integer   nprold, nshl,  i
      
c
c.... map the array
c
      do i = 1, nshl
         x2(map(:),i) = x(:,i)
      enddo

      return 
      end
