      subroutine AsBNABI ( x,       shpb,
     &                   ienb,  iBCB)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  boundary elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use pvsQbi
        include "common.h"
c
        dimension xlb(npro,nenl,nsd),    bnorm(npro,nsd),
     &            rl(npro,nshl,nsd),     WdetJb(npro),
     &            Wfactor(npro)

        dimension x(numnp,nsd),
     &            shpb(nshl,ngaussb),      shglb(nsd,nshl,ngaussb),
     &            ienb(npro,nshl),
     &            iBCB(npro,ndiBCB)

        dimension lnode(27),               sgn(npro,nshl),
     &            shpfun(npro,nshl),        shdrv(npro,nsd,nshl)

c
        dimension dxdxib(npro,nsd,nsd),      temp(npro),
     &            temp1(npro),               temp2(npro),
     &            temp3(npro),
     &            v1(npro,nsd),              v2(npro,nsd)

c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c
c.... gather the variables
c
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
c
c.... get the boundary element residuals
c
        rl  = zero
c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
c
  	if(lcsyst.eq.1) then	  ! set to curl into element all others out
  	   ipt2=2
  	   ipt3=3
  	elseif(lcsyst.eq.2) then
  	   ipt2=4
  	   ipt3=2
  	elseif(lcsyst.eq.3) then
  	   ipt2=3
  	   ipt3=2
  	elseif(lcsyst.eq.4) then
  	   ipt2=2
  	   ipt3=4
  	elseif(lcsyst.eq.5) then
  	   ipt2=4
  	   ipt3=2
  	elseif(lcsyst.eq.6) then
  	   ipt2=2
  	   ipt3=5
  	endif
  	v1 = xlb(:,ipt2,:) - xlb(:,1,:)
  	v2 = xlb(:,ipt3,:) - xlb(:,1,:)
c
c compute cross product
c
        temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
        temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
        temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
c     
c mag is area for quads, twice area for tris
c 
  	temp	   = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
  	bnorm(:,1) = temp1 * temp
  	bnorm(:,2) = temp2 * temp
  	bnorm(:,3) = temp3 * temp
c
  	if (lcsyst .eq. 1) then
  	   Wfactor(:) = one / (four*temp(:))
  	elseif (lcsyst .eq. 2) then
  	   Wfactor(:) = one / (four*temp(:))
  	elseif (lcsyst .eq. 3) then
  	   Wfactor(:) = one / (two*temp(:))
  	elseif (lcsyst .eq. 4) then
  	   Wfactor(:) = one / (four*temp(:))
  	elseif (lcsyst .eq. 5) then
  	   Wfactor(:) = one / (four*temp(:))
  	elseif (lcsyst .eq. 6) then
  	   Wfactor(:) = one / (two*temp(:))
  	endif
c
        do iel=1,npro
           if (.not.btest(iBCB(iel,1),1)) then
              Wfactor(iel) = zero  ! we want zeros where we are not integrating
           endif
        enddo
c
        ngaussb = nintb(lcsyst)
c
c.... loop through the integration points
c
        do intp = 1, ngaussb
c
c.... get the hierarchic shape functions at this int point
c
           shglb=zero  ! protect debugger 
           call getshpb(shpb,        shglb,        sgn, 
     &              shpfun,       shdrv)
c
           WdetJb(:) = Qwtb(lcsyst,intp) * Wfactor(:)
c
c  Now lets calculate Integral N_(a:e)^i n_i d Gamma
c
           do n = 1, nshlb
              nodlcl = lnode(n)
              rl(:,nodlcl,1) = rl(:,nodlcl,1) + 
     &             shpfun(:,nodlcl) * bnorm(:,1) * WdetJb(:)
              rl(:,nodlcl,2) = rl(:,nodlcl,2) + 
     &             shpfun(:,nodlcl) * bnorm(:,2) * WdetJb(:)
              rl(:,nodlcl,3) = rl(:,nodlcl,3) + 
     &             shpfun(:,nodlcl) * bnorm(:,3) * WdetJb(:)
           enddo

        enddo  ! quadrature point loop
c
c.... assemble the NABI vector
c
        call local (NABI,    rl,     ienb,   3,  'scatter ')
c
c     push the surf number which we have associated with boundary
C     elements up to the global level in the array ndsurf
c
        do iel=1,npro
           if (iBCB(iel,2) .ne. 0) then
              iface = iBCB(iel,2)
              ndsurf(ienb(iel,1:nshlb))=iface   
           endif
        enddo
c     
c.... end
c
        return
        end


      subroutine AsBNASC ( x,       shpb,
     &                   ienb,  iBCB)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  boundary elements.
c
c Irene Vignon - copied from AsBNABI, Fall 2005.  (Fortran 90)
c----------------------------------------------------------------------
c
      use pvsQbi
        include "common.h"
c
        dimension xlb(npro,nenl,nsd),
     &            rl(npro,nshl),         WdetJb(npro),
     &            Wfactor(npro)

        dimension x(numnp,nsd),
     &            shpb(nshl,ngaussb),    shglb(nsd,nshl,ngaussb),
     &            ienb(npro,nshl),
     &            iBCB(npro,ndiBCB)

        dimension lnode(27),             sgn(npro,nshl),
     &            shpfun(npro,nshl),     shdrv(npro,nsd,nshl)

c
        dimension dxdxib(npro,nsd,nsd),  temp(npro),
     &            temp1(npro),           temp2(npro),
     &            temp3(npro),
     &            v1(npro,nsd),          v2(npro,nsd)

c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c
c.... gather the variables
c
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
c
c.... get the boundary element residuals
c
        rl  = zero
c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
c
  	if(lcsyst.eq.1) then	  ! set to curl into element all others out
  	   ipt2=2
  	   ipt3=3
  	elseif(lcsyst.eq.2) then
  	   ipt2=4
  	   ipt3=2
  	elseif(lcsyst.eq.3) then
  	   ipt2=3
  	   ipt3=2
  	elseif(lcsyst.eq.4) then
  	   ipt2=2
  	   ipt3=4
  	elseif(lcsyst.eq.5) then
  	   ipt2=4
  	   ipt3=2
  	elseif(lcsyst.eq.6) then
  	   ipt2=2
  	   ipt3=5
  	endif
  	v1 = xlb(:,ipt2,:) - xlb(:,1,:)
  	v2 = xlb(:,ipt3,:) - xlb(:,1,:)
c
c compute cross product
c
        temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
        temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
        temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
c     
c mag is area for quads, twice area for tris
c 
  	temp	   = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
c
c......here I only want the d Gamma, not the n_i
   	if (lcsyst .eq. 1) then
   	   Wfactor(:) = one / (four*temp(:))
   	elseif (lcsyst .eq. 2) then
   	   Wfactor(:) = one / (four*temp(:))
   	elseif (lcsyst .eq. 3) then
   	   Wfactor(:) = one / (two*temp(:))
   	elseif (lcsyst .eq. 4) then
   	   Wfactor(:) = one / (four*temp(:))
   	elseif (lcsyst .eq. 5) then
   	   Wfactor(:) = one / (four*temp(:))
   	elseif (lcsyst .eq. 6) then
   	   Wfactor(:) = one / (two*temp(:))
   	endif
c
        do iel=1,npro
           if (.not.btest(iBCB(iel,1),1)) then
              Wfactor(iel) = zero  ! we want zeros where we are not integrating
           endif
        enddo
c
        ngaussb = nintb(lcsyst)
c
c.... loop through the integration points
c
        do intp = 1, ngaussb

           WdetJb(:) = Qwtb(lcsyst,intp) * Wfactor(:)
c
c.... get the hierarchic shape functions at this int point
c
           shglb=zero  ! protect debugger 
           call getshpb(shpb,        shglb,        sgn, 
     &              shpfun,       shdrv)

c
c  Now lets calculate Integral N_(a:e) d Gamma
c
           do n = 1, nshlb
              nodlcl = lnode(n)
              rl(:,nodlcl) = rl(:,nodlcl) + shpfun(:,nodlcl) * WdetJb(:)         
           enddo

        enddo  ! quadrature point loop
c
c.... assemble the NASC vector
c
        call local (NASC,    rl,     ienb,   1,  'scatter ')
c
c     push the surf number which we have associated with boundary
C     elements up to the global level in the array ndsurf
c
        do iel=1,npro
           if (iBCB(iel,2) .ne. 0) then
              iface = iBCB(iel,2)
              ndsurf(ienb(iel,1:nshlb))=iface   
           endif
        enddo
c     
c.... end
c
        return
        end
