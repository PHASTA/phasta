      subroutine SolFlow(y,          ac,         u,
     &                   yold,       acold,      uold,
     &                   x,          iBC,
     &                   BC,         res,        iper,       
     &                   ilwork,     shp,        shgl, 
     &                   shpb,       shglb,      rowp,     
     &                   colm,       
     &                   solinc,     rerr,       tcorecp,
     &                   GradV,       sumtime
#ifdef HAVE_SVLS     
     &                   ,svLS_lhs,  svLS_ls,   svLS_nFaces)
#else
     &                   )
#endif
c
c----------------------------------------------------------------------
c
c This is the 2nd interface routine to the  linear equation 
c solver library that uses the CGP and GMRES methods.
c
c input:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  yold   (nshg,ndof)           : Y-variables at beginning of step
c  acold   (nshg,ndof)          : Primvar. accel. at beginning of step
c  x      (numnp,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodic nodal information
c
c output:
c  res    (nshg,nflow)           : preconditioned residual
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c
c
c The followings are preliminary steps required to use Farzin's
c solver library.  New way of writing has to be used such as
c
c          |  K     G | | du |    | Rmom  |
c          |          | |    | =  |       |
c          | G^t    C | | dp |    | Rcon  |
c
c          |     E    | | dT | =  | Rtemp |
c
c     where
c
c      xKebe : K_ab = dRmom_a/du_b    xTe : E_ab = dRtemp_a/dT_b 
c
c              G_ab = dRmom_a/dp_b
c      xGoC  :
c              C_ab = dRcon_a/dp_b       
c
c              resf = Rmon Rcon       rest = Rtemp
c
c  
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c Juin Kim, Summer 1998. (Incompressible flow solver interface)
c Alberto Figueroa.  CMM-FSI
c----------------------------------------------------------------------
c
      use pointer_data
      use solvedata
#ifdef AMG      
      use ramg_data
#endif     
        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
#ifdef HAVE_SVLS      
        include "svLS.h"
#endif        
C
      REAL*8                rdtmp
C    
#ifdef HAVE_SVLS
      TYPE(svLS_lhsType), INTENT(INOUT) :: svLS_lhs
      TYPE(svLS_lsType), INTENT(INOUT) ::  svLS_ls
#endif      
       
      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          yold(nshg,ndof),          acold(nshg,ndof),
     &          u(nshg,nsd),              uold(nshg,nsd),
     &          x(numnp,nsd),             BC(nshg,ndofBC),
     &          res(nshg,nflow),          tmpres(nshg,nflow),
     &          flowDiag(nshg,4),         
     &          sclrDiag(nshg,1),         
     &          GradV(nshg,nsdsq)
c
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
      integer   usr(100),                 eqnType,temp,
     &          rowp(nshg*nnz),           colm(nshg+1),
     &          iBC(nshg),                ilwork(nlwork),
     &          iper(nshg)
c
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof),
     &          uAlpha(nshg,nsd),
     &          lesP(nshg,4),             lesQ(nshg,4),
     &          solinc(nshg,ndof),        CGsol(nshg)

      real*8     tcorecp(2)
      
      real*8    rerr(nshg,10),            rtmp(nshg,4),rtemp
      
      real*8    msum(4),mval(4),cpusec(10)
      REAL*8 sumtime
#ifdef HAVE_SVLS      
      INTEGER svLS_nFaces
#endif      
      INTEGER dof, i, j, k, l
      INTEGER, ALLOCATABLE :: incL(:)
      REAL*8, ALLOCATABLE :: faceRes(:), Res4(:,:), Val4(:,:)

c     
c.... *******************>> Element Data Formation <<******************
c
c
c.... set the parameters for flux and surface tension calculations
c
c
        temp = npro
        

        idflx = 0 
        if(idiff >= 1 )  idflx= (nflow-1) * nsd
        if (isurf == 1) idflx=nflow*nsd
c.... compute solution at n+alpha
c
      call itrYAlpha( uold,    yold,    acold,
     &                u,       y,       ac,  
     &                uAlpha,  yAlpha,  acAlpha)

c
c.... form the LHS matrices, the residual vector (at alpha)
c
c      call summary_start()
      impistat=1
      impistat2=1
      telmcp1 = TMRC()
      call ElmGMR (uAlpha,    yAlpha,     acAlpha,    x,
     &             shp,       shgl,       iBC,       
     &             BC,        shpb,       shglb,
     &             res,       iper,       ilwork,   
     &             rowp,      colm,       lhsK,      
     &             lhsP,      rerr,       GradV   )
      telmcp2 = TMRC()
      impistat=0
      impistat2=0
c      call summary_stop()


            tmpres(:,:) = res(:,:)
            iblk = 1
#ifdef HAVE_SVLS            
      IF (svLSFlag .EQ. 1) THEN

c####################################################################
!     Here calling svLS

      ALLOCATE(faceRes(svLS_nFaces), incL(svLS_nFaces))
      faceRes=zero  ! function to compute this left out at this time but would be needed to enable adnvanced p vs. Q BC's
      incL = 1
      dof = 4
      IF (.NOT.ALLOCATED(Res4)) THEN
         ALLOCATE (Res4(dof,nshg), Val4(dof*dof,nnz_tot))
      END IF

      DO i=1, nshg
         Res4(1:dof,i) = res(i,1:dof)
      END DO

      DO i=1, nnz_tot
         Val4(1:3,i)   = lhsK(1:3,i)
         Val4(5:7,i)   = lhsK(4:6,i)
         Val4(9:11,i)  = lhsK(7:9,i)
         Val4(13:15,i) = lhsP(1:3,i)
         Val4(16,i)    = lhsP(4,i)
      END DO

      !Val4(4:12:4,:) = -lhsP(1:3,:)^t
      DO i=1, nshg
         Do j=colm(i), colm(i+1) - 1
            k = rowp(j)
            DO l=colm(k), colm(k+1) - 1
               IF (rowp(l) .EQ. i) THEN
                  Val4(4:12:4,l) = -lhsP(1:3,j)
                  EXIT
               END IF
            END DO
         END DO
      END DO
      CALL svLS_SOLVE(svLS_lhs, svLS_ls, dof, Res4, Val4, incL, 
     2   faceRes)
      
      if(myrank.eq.master) write(*,*) 'svLS outer iterations', svLS_ls%RI%itr
      statsflow(1)=1.0*svLS_ls%GM%itr
      statsflow(4)=1.0*svLS_ls%CG%itr
      DO i=1, nshg
         solinc(i,1:dof) = Res4(1:dof,i)
      END DO
      ENDIF 
#endif

#ifdef HAVE_LESLIB  
      if(leslib.eq.1) then    
c
c.... lesSolve : main matrix solver
c
      lesId   = numeqns(1)
      eqnType = 1

c      call summary_start()
      impistat=1
      impistat2=1
      tlescp1 = TMRC()
#ifdef AMG      
      ! Initial Time Profiling
      call cpu_time(cpusec(1))
      if (irun_amg_prec.gt.0) then
          call ramg_control(colm,rowp,lhsK,lhsP,
     &         ilwork,BC,iBC,iper)
      end if

      call cpu_time(cpusec(6))
      if (irun_amg_prec.gt.0) then
      ramg_flag = 1
      if (irun_amg_prec.eq.2) then ! Some setup work (mode a)
        ramg_window = 1.0
        ramg_redo = 0
      endif
      do while (ramg_flag.le.irun_amg_prec) 
      ! if smart solve, possible run solve twice
      ! restart only if meets plateau
#endif      
      
c
c.... setup the linear algebra solver
c
      rtmp = res(:,1:4)
      call usrNew ( usr,        eqnType,          aperm,
     &              atemp,      rtmp,             solinc,          
     &              flowDiag,   sclrDiag,         lesP,   
     &              lesQ,       iBC,              BC,
     &              iper,       ilwork,           numpe,
     &              nshg,       nshl,             nPermDims,  
     &              nTmpDims,   rowp,             colm,     
     &              lhsK,       lhsP,             rdtmp,      
     &              nnz_tot,    CGsol )
c
c.... solve linear system
c
      call myfLesSolve ( lesId, usr )
#ifdef AMG
      ramg_flag = ramg_flag + 2 ! Default no second run in mode a
      if (irun_amg_prec.eq.3) then
          if (maxIters.gt.int(statsflow(4))) then
          ramg_flag = ramg_flag + 1 ! Default no second run in mode b
          endif
      endif
      enddo
      else
c
c.... setup the linear algebra solver
c
      rtmp = res(:,1:4)
      call usrNew ( usr,        eqnType,          aperm,
     &              atemp,      rtmp,             solinc,          
     &              flowDiag,   sclrDiag,         lesP,   
     &              lesQ,       iBC,              BC,
     &              iper,       ilwork,           numpe,
     &              nshg,       nshl,             nPermDims,  
     &              nTmpDims,   rowp,             colm,     
     &              lhsK,       lhsP,             rdtmp,      
     &              nnz_tot,    CGsol )

          call myfLesSolve( lesId, usr )
      endif

      call cpu_time(cpusec(3))

      ramg_time(1) = ramg_time(1) + cpusec(3)-cpusec(1)
      ramg_time(7) = ramg_time(7) + cpusec(6)-cpusec(1)

      ! ramg_time: 1 : local total
      !            4 : local VG-cycle
      !            7 : local setup time
      !           11 : Ap-product level 1
      !           12 : Ap-product level >1
      !           13 : Prolongation/Restriction
      !           20 : local boundary MLS time
     
      if (myrank.eq.master) then
      write(*,*)
      endif
      call ramg_print_time(" == AMG == Total ACUSIM Solver:",
     &                    ramg_time(1))
      call ramg_print_time(" == AMG == Setup: ",ramg_time(7))
      call ramg_print_time(" == AMG == Prec'd cycle: ",ramg_time(4))
      call ramg_print_time(" == AMG == Ap product(level=1): ",
     &                      ramg_time(11))
      call ramg_print_time(" == AMG == Ap product(level>=2): ",
     &                      ramg_time(12))
      call ramg_print_time(" == AMG == Pro/Restr ",
     &                      ramg_time(13))
      call ramg_print_time(" == AMG == Boundary Ap (GS only)",
     &                      ramg_time(20))
      if (myrank.eq.master) then
      write(*,*)
      endif

#endif     
      
      ! End Time profiling output
      
      call getSol ( usr, solinc )

      if (numpe > 1) then
         call commu ( solinc, ilwork, nflow, 'out')
      endif
      ENDIF ! end of leslib flow solve
#endif   
      tlescp2 = TMRC()
      impistat=0
      impistat2=0
c      call summary_stop()

      tcorecp(1) = tcorecp(1) + telmcp2-telmcp1 ! elem. formation
      tcorecp(2) = tcorecp(2) + tlescp2-tlescp1 ! linear alg. solution
      call rstatic (res, y, solinc) ! output flow stats
c     
c.... end
c     
      return
      end

      subroutine SolSclr(y,          ac,         u,
     &                   yold,       acold,      uold,
     &                   x,          iBC,
     &                   BC,         iper,       
     &                   ilwork,     shp,        shgl, 
     &                   shpb,       shglb,      rowp,     
     &                   colm,       solinc,
     &                   tcorecpscal
#ifdef HAVE_SVLS     
     &                   ,svLS_lhs,  svLS_ls,   svLS_nFaces)
#else
     &                   )      
#endif      
c
c----------------------------------------------------------------------
c
c This is the 2nd interface routine to the Farzin's linear equation 
c solver library.
c
c input:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  yold   (nshg,ndof)           : Y-variables at beginning of step
c  x      (numnp,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodic nodal information
c
c output:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c
c
c The followings are preliminary steps required to use Farzin's
c solver library.  New way of writing has to be used such as
c
c          |     E    | | dS | =  | RScal |
c
c----------------------------------------------------------------------
c
      use pointer_data
      use solvedata
        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
#ifdef HAVE_SVLS      
        include "svLS.h"
#endif        
c     
      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          yold(nshg,ndof),          acold(nshg,ndof),
     &          u(nshg,nsd),              uold(nshg,nsd),
     &          x(numnp,nsd),             BC(nshg,ndofBC),
     &          res(nshg,1),
     &          flowDiag(nshg,4),
     &          sclrDiag(nshg,1)           
c
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
      integer   usr(100),                 eqnType,
     &          rowp(nshg*nnz),           colm(nshg+1),
     &          iBC(nshg),                ilwork(nlwork),
     &          iper(nshg)
c
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof),
     &          uAlpha(nshg,nsd),
     &          lesP(nshg,1),               lesQ(nshg,1),
     &          solinc(nshg,1),           CGsol(nshg),
     &          tcorecpscal(2)
#ifdef HAVE_SVLS     
      TYPE(svLS_lhsType), INTENT(INOUT) :: svLS_lhs
      TYPE(svLS_lsType), INTENT(INOUT) ::  svLS_ls
      INTEGER svLS_nFaces
#endif      
      REAL*8 sumtime
      INTEGER dof, i, j, k, l
      INTEGER, ALLOCATABLE :: incL(:)
      REAL*8, ALLOCATABLE :: faceRes(:), Res1(:,:), Val1(:,:)
      
c     
c.... *******************>> Element Data Formation <<******************
c
c.... compute solution at n+alpha
c
      call itrYAlpha( uold,    yold,    acold, 
     &                u,       y,       ac,  
     &                uAlpha,  yAlpha,  acAlpha)
c
c.... form the LHS matrices, the residual vector (at alpha)
c
      impistat=2
      impistat2=1
      telmcp1 = TMRC()
      jsol=nsolt+isclr
      call ElmGMRSclr(yAlpha,acAlpha,    x,
     &             shp,       shgl,       iBC,       
     &             BC,        shpb,       shglb,
     &             res,       iper,       ilwork,   
     &             rowp,      colm,       lhsS   )
      telmcp2 = TMRC()
      impistat=0
      impistat2=0
      statssclr(1)=0
#ifdef HAVE_SVLS      
      IF (svLSFlag .EQ. 1) THEN

c####################################################################
!     Here calling svLS

      ALLOCATE(faceRes(svLS_nFaces), incL(svLS_nFaces))
      faceRes=zero  
      incL = 1
      dof = 1
      IF (.NOT.ALLOCATED(Res1)) THEN
         ALLOCATE (Res1(dof,nshg), Val1(dof*dof,nnz_tot))
      END IF

      DO i=1, nshg
         Res1(1,i) = res(i,1)
      END DO

      DO i=1, nnz_tot
         Val1(1,i)    = lhsS(i,jsol) ! see above jsol indexs for scalars
      END DO

      CALL svLS_SOLVE(svLS_lhs, svLS_ls, dof, Res1, Val1, incL, 
     2   faceRes)
      statssclr(1)=1.0*svLS_ls%RI%itr
      DO i=1, nshg
         solinc(i,1) = Res1(1,i)
      END DO
      ENDIF
#endif          
#ifdef HAVE_LESLIB
      if(leslib.eq.1) then
c
c.... lesSolve : main matrix solver
c
      lesId   = numeqns(1+nsolt+isclr)
      eqnType = 2
c
c.... setup the linear algebra solver
c

      impistat=2
      impistat2=1
      tlescp1 = TMRC()
      call usrNew ( usr,        eqnType,          apermS(1,1,jsol),
     &              atempS,     res,              solinc,          
     &              flowDiag,   sclrDiag,         lesP,   
     &              lesQ,       iBC,              BC,
     &              iper,       ilwork,           numpe,
     &              nshg,       nshl,             nPermDimsS,  
     &              nTmpDimsS,  rowp,             colm,     
     &              rlhsK,      rlhsP,            lhsS,      
     &              nnz_tot,    CGsol )
c
c.... solve linear system
c
      call myfLesSolve ( lesId, usr )
      call getSol ( usr, solinc )

      if (numpe > 1) then
         call commu ( solinc, ilwork, 1, 'out')
      endif      
      ENDIF ! leslib conditional
#endif      
      tlescp2 = TMRC()
      impistat=0
      impistat2=0

      tcorecpscal(1) = tcorecpscal(1) + telmcp2-telmcp1 ! elem. formation
      tcorecpscal(2) = tcorecpscal(2) + tlescp2-tlescp1 ! linear alg. solution
 
      nsolsc=5+isclr
      call rstaticSclr (res, y, solinc, nsolsc) ! output scalar stats
c     
c.... end
c     
      return
      end
