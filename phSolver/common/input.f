        subroutine input()
c
c----------------------------------------------------------------------
c This routine inputs all the necessary data, allocates required array 
c storage, and sets up the appropriate parameters for the processing.
c 
c
c Farzin Shakib, Winter 1987.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
        include "mpif.h"

        external endata

        integer, allocatable :: nsons(:)
c
        character*8  date
        character*80 card

c assigned in phasta.cc
c        numpe=npe
c        myrank=mrank

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        rsec=TMRC()
        ttim(100) = rsec

        epsM = sqrt(epsilon(one))
c
c.... read in and block all data
c
        call readnblk()
c
c.... open the echo file (echo closed at exit)
c
        if (myrank == master) 
     &  open (unit=iecho, file=fecho, status='unknown',   err=996)

c
c.... -------------------->  Control Parameters  <---------------------
c
c.... echo the global information
c

        title = 'Default Ensa Case'     
        call date_and_time (date)
        title  = title(1:69) // ' ' // date(7:8) // '/' // date(5:6)
     &                                           // '/' // date(3:4)
        ititle = char(12) // title(1:78)

        if (myrank == master) then
          write (iecho,1100) ititle, numpe,  numnp,  numel,  numelb,
     &                               nen,    nfaces, nsd,    numflx
          write (iecho,1200)         iALE,   icoord, navier, irs,
     &                               iexec,  necho
c
c.... check the input parameters
c
          if (iALE .lt. 0 .or. iALE .gt. 1)
     &                     call error ('input   ','iALE    ',iALE)
c
          if (icoord .lt. 0 .or. icoord .gt. 1)
     &                     call error ('input   ','icoord  ',icoord)
c
          if (navier .lt. 0 .or. navier .gt. 1)
     &                     call error ('input   ','navier  ',navier)

          if (irs    .lt. 0 .or. irs    .gt. 3)
     &                     call error ('input   ','irs     ',irs)
c
          if (iexec  .lt. 0 .or. iexec  .gt. 1)
     &                     call error ('input   ','iexec   ',iexec)
c
          if (necho  .lt. 0 .or. necho  .gt. 3)
     &                     call error ('input   ','necho   ',necho)
        endif

        if (myrank == master) then
          write (iecho,1300) ititle, ntseq,  imap,   ivart,  iDC,
     &                             Kspace, nGMRES
c
c.... check the input parameters
c
          if (ntseq  .gt. 100) call error ('input   ','ntseq   ',ntseq)
c
          if (imap   .lt. 0 .or. imap  .gt. 1) 
     &                       call error ('input   ','imap    ',imap)
c
          if (ivart  .lt. 1 .or. ivart .gt. 3) 
     &                       call error ('input   ','ivart   ',ivart)
c
          if (iDC    .lt. 0 .or. iDC   .gt. 4) 
     &                       call error ('input   ','iDC     ',iDC)
c
          if (Kspace .lt. 1)   call error ('input   ','Kspace  ',Kspace)
c
          if (nGMRES .lt. 1)   call error ('input   ','nGMRES  ',nGMRES)
        endif
c
c.... ----------------->  Time Sequence Parameters  <-----------------
c
c.... echo the solver information
c
        iprev = 0
        do i = 1, ntseq
          if (mod(i,50).eq.1 .and. myrank .eq. master) 
     &      write(iecho,1400) ititle

          if (myrank .eq. master) 
     &      write (iecho,1500)      i, nstep(i),  niter(i),  impl(i),
     &                                 LHSupd(i), epstol(i)
c
          if ((iALE .eq. 1) .or. (niter(i) .gt. 1)) iprev = 1
        enddo
c
c.... echo the spatial and time integration information
c
        do i = 1, ntseq
          if (mod(i,50).eq.1 .and. myrank .eq. master) 
     &      write(iecho,1600) ititle
          if (myrank .eq. master) 
     &      write (iecho,1700)      i, intg(1,i), intg(2,i), rhoinf(i),
     &                                 loctim(i), Delt(i),   CFLfl(i),
     &                                 CFLsl(i)
c
        enddo
c
        if (myrank .eq. master) 
     &    write (iecho,1800) ititle, ntout,  ioform, ro,     vel,    
     &                               temper, press,  entrop
        
        if (myrank .eq. master) then
           write (*,*) 'Element block size = ',ibksiz
        endif

        if (iLSet .gt. 0 .and. myrank .eq. master)then 
           write(iecho,1900)iLSet, epsilon_ls
        endif
c
c.... generate the spatial integration rules
c
        call genint

        ichem = 0
c
c.... estimate number of nonzero global entries:
c....       nnonzero ~ nnz * nshg
c
        if (ipord .eq. 1) then
           nnz = 35
        else if (ipord .eq. 2) then
           nnz = 85
        else  !assumed cubic
           nnz = 300
        endif


c
c.... compute fluid thermodynamic properties
c
        Boltzm = Rh / Nh
c
        do i = 1, 5
          Rs(i)   = Rh / Msh(i)
          h0s(i)  = h0sh(i) / Msh(i)
          cpsh(i) = ( pt5 * dofs(i) + one ) * Rh
          cps(i)  = ( pt5 * dofs(i) + one ) * Rs(i)
          cvs(i)  = pt5 * dofs(i) * Rs(i)
        enddo
c
        do i = 1, 5
          s0sh(i) = Rh * ( pt5*( log( (two*pi*Msh(i)/(Nh*Planck**2))**3
     &                  * Boltzm**5 ) + five ) + log(g0s(i)) )
        enddo
c
        do i = 1, 3
          s0sh(i) = s0sh(i) + Rh * ( one - log(sigs(i)*Trot(i)) )
        enddo
c
        Rgas  = one / ( xN2 / Rs(1) + xO2 / Rs(2) ) 
        if(myrank.eq.0) write(*,*) 'input.f computes Rgas to be', Rgas
c        Rgas  = 0.4*716.5
c        Rgas = 8314/28.95
        yN2   = xN2 * Rgas / Rs(1)
        yO2   = xO2 * Rgas / Rs(2)
c
        s0    =     yN2 * s0sh(1) / Msh(1) + yO2 * s0sh(2) / Msh(2)
        const = - ( yN2 * Rs(1) * log(xN2) + yO2 * Rs(2) * log(xO2) )
c
c.... stop CPU-timer
c
c        call timer ('Back    ')
cc
c..dumping common (useful for checking differences with
c        old format input
c
        if(myrank.eq.-1) then
        mxats=1
        open (unit=23,   file="dumpnew.dat",   status='unknown')
        write (23,*)" master, numpe, myrank"
        write (23,*) master, numpe, myrank
        write (23,*)" maxfront, nlwork"
        write (23,*) maxfront, nlwork
        write (23,*)"  numper, nshgt, nshg0"
        write (23,*)  numper, nshgt, nshg0
        write (23,*) " birth, death, comtim"
        write (23,*)  birth, death, comtim
        write (23,*)"  pzero, wtavei,dtavei, dke, ierrcalc,"
        write (23,*)  pzero, wtavei,dtavei, dke, ierrcalc,
     &                   itwmod,taucfct
        write (23,*)"irscale, intpres,rxinlt, rxrecy,"
        write (23,*)irscale, intpres,rxinlt, rxrecy,
     &            rbltin,rvscal,  xlngth, ylngth, zlngth 
 
        write (23,*)"  scdiff(5),nsclr,isclr,nsolt"
        write (23,*)  scdiff(5),nsclr,isclr,nsolt
        write (23,*) " flxID(10,20), Force(3),HFlux, nsrflist(0:20)"
        write (23,*)  flxID(10,20), Force(3),HFlux, nsrflist(0:20)
        write (23,*) " numnp,  numel,  numelb, numpbc, nen,    nfaces,"
        write (23,*)  numnp,  numel,  numelb, numpbc, nen,    nfaces,
     &                  numflx, ndof,   iALE,   icoord, navier,
     &                  irs,    iexec,  necho,  ichem,  iRK,    nedof,
     &                  nshg,   nnz,    istop,  nflow,  nnz_tot, idtn,
     &                  iLSet
        write (23,*)"  epsilon_ls, epsilon_lsd, dtlset"
        write (23,*)  epsilon_ls, epsilon_lsd, dtlset
        write (23,*)" nshape, nshapeb, maxshb,"
        write (23,*) nshape, nshapeb, maxshb,
     &                  nshl, nshlb,nfath,  ntopsh,  nsonmax
        write (23,*)" mcsyst, melCat, nenCat(8,3),    nfaCat(8,3)"
        write (23,*) mcsyst, melCat, nenCat(8,3),    nfaCat(8,3)
 
        write (23,*)" lelCat, lcsyst, iorder, nenb, "  
        write (23,*) lelCat, lcsyst, iorder, nenb,   
     &                  nelblk, nelblb, ndofl,  nsymdl, nenl,   nfacel,
     &                  nenbl,  intind, mattyp 
        write (23,*)" E3nsd,  I3nsd,  nsymdf, ndofBC, ndiBCB, ndBCB,"
        write (23,*) E3nsd,  I3nsd,  nsymdf, ndofBC, ndiBCB, ndBCB,
     &                  Jactyp, jump,   ires,   iprec,  ibound,
     &                  idiff,  lhs,    itau,   ipord,  ipred,  lstres,
     &                  iepstm, dtsfct, ibksiz, iabc
        write (23,*)"  epstol(mxats),  Delt(mxats),"
        write (23,*)  epstol(mxats),  Delt(mxats),     nstep(mxats), 
     &                  impl(mxats),    rhoinf(mxats),
     &                  LHSupd(mxats),  loctim(mxats),  deltol(mxats,2)

        write (23,*)" intg(2,mxats),  intpt(3),       intptb(3)"
        write (23,*) intg(2,mxats),  intpt(3),       intptb(3)
        write (23,*) " indQpt(3,3,4),  numQpt(3,3,4),"
        write (23,*)  indQpt(3,3,4),  numQpt(3,3,4),
     &                  intmax
        write (23,*)" iin,    igeom,  ipar,   ibndc,  imat,   iecho,"
        write (23,*) iin,    igeom,  ipar,   ibndc,  imat,   iecho,
     &                  iout,   ichmou, irstin, irstou, ihist,  iflux,
     &                  ierror, itable, iforce, igraph, itime 
        write (23,*)"fwr1,ngaussf,idim,nlist "
        write (23,*)fwr1,ngaussf,idim,nlist 
        write (23,*) " fin,    fgeom,  fpar,   fbndc,  fmat,   fecho,"
        write (23,*)  fin,    fgeom,  fpar,   fbndc,  fmat,   fecho,
     &                  frstin, frstou, fhist,  ferror, ftable, fforce,
     &                  fgraph, ftime 
        write (23,*)" fin,    fgeom,  fpar,   fbndc,  fmat,   fecho,"
        write (23,*) fin,    fgeom,  fpar,   fbndc,  fmat,   fecho,
     &                  frstin, frstou, fhist,  ferror, ftable, fforce,
     &                  fgraph, ftime 
        write (23,*)" eGMRES, lGMRES, iKs,    ntotGM "
        write (23,*) eGMRES, lGMRES, iKs,    ntotGM 
        write (23,*) " mHBrg,  meBrg,  myBrg,  mRcos,  mRsin"
        write (23,*)  mHBrg,  meBrg,  myBrg,  mRcos,  mRsin
c
        write (23,*)" pr,     Planck, Stefan, Nh,     Rh,     Rgas,"
        write (23,*) pr,     Planck, Stefan, Nh,     Rh,     Rgas,
     &                  gamma,  gamma1, s0,     const,  xN2,    xO2,
     &                  yN2,    yO2,    Msh(5), cpsh(5),s0sh(5),h0sh(5),
     &                  Rs(5),  cps(5), cvs(5), h0s(5), Trot(5),sigs(5),
     &                  Tvib(5),g0s(5), dofs(5),ithm 
        write (23,*) " mexist"
        write (23,*)  mexist
        write (23,*) " datmat(3,5,mxats),      matflg(5,mxats),"
        write (23,*)  datmat(3,5,mxats),      matflg(5,mxats),
     &                  nummat,                 mexist
        write (23,*)"ro,     vel,    temper, press,  entrop, ntout,"
        write (23,*)ro,     vel,    temper, press,  entrop, ntout,
     &                  ioform 
        write (23,*)"mbeg,   mend,   mprec "
        write (23,*)mbeg,   mend,   mprec 
        write (23,*)"epsM,   iabres, npro,resfrt"
        write (23,*)epsM,   iabres, npro,resfrt
        write (23,*)"  imap,   ivart,  iDC,    iPcond, Kspace, nGMRES"
        write (23,*)  imap,   ivart,  iDC,    iPcond, Kspace, nGMRES

        write (23,*)" indsym(5,5) "
        write (23,*) indsym(5,5) 
        write (23,*) " time,   CFLfld, CFLsld, Dtgl,   Dtmax,  alpha,"
        write (23,*)  time,   CFLfld, CFLsld, Dtgl,   Dtmax,  alpha,
     &                  etol,   lstep,  ifunc,  itseq,  istep,  iter,
     &                  nitr,   almi,   alfi,   gami,   flmpl,  flmpr,
     &                  dtol(2) 
        write (23,*) "LCtime, ntseq"
        write (23,*) LCtime, ntseq
        write (23,*) " numeqns(100), minIters, maxIters," 
        write (23,*)  numeqns(100), minIters, maxIters, 
     &                  iprjFlag,     nPrjs,    ipresPrjFlag, nPresPrjs,
     &                  prestol,      statsflow(6), statssclr(6),
     &                  iverbose
        write (23,*) " ccode" 
        write (23,*)  ccode 
        write (23,*) " flops,  gbytes, sbytes, iclock, icd,    icode,"
        write (23,*)  flops,  gbytes, sbytes, iclock, icd,    icode,
     &                  icode2, icode3
        write (23,*) " cpu(11),        cpu0(11),       nacess(11)"
        write (23,*)  cpu(11),        cpu0(11),       nacess(11)
        write (23,*) " title,  ititle"
        write (23,*)  title,  ititle
        close (23)
        endif
c
c....return
c
        return
c
c.... end of file error handling
c
992     call error ('input   ','opening ', imat)
993     call error ('input   ','opening ', iin)
996     call error ('input   ','opening ', iecho)
999     call error ('input   ','end file', iin)
c
1000    format(a69)
1100    format(a80,//,
     &  ' M a i n   C o n t r o l   P a r a m e t e r s        '   //,
     &  ' number of processing elements . . . . . . . (numpe )=',i10//,
     &  ' number of mesh nodes  . . . . . . . . . . . (numnp )=',i10//,
     &  ' number of elements  . . . . . . . . . . . . (numel )=',i10//,
     &  ' number of boundary elements . . . . . . . . (numelb)=',i10//,
     &  ' number of element nodes . . . . . . . . . . (nen   )=',i10//,
     &  ' number of element faces . . . . . . . . . . (nfaces)=',i10//,
     &  ' number of space dimensions  . . . . . . . . (nsd   )=',i10//,
     &  ' number of boundary flux nodes . . . . . . . (numflx)=',i10/)
1200    format(
     &  ' frame of reference  . . . . . . . . . . . . (iALE  )=',i10//,
     &  '    eq. 0, Eulerian                                   ',  / ,
     &  '    eq. 1, arbitrary Lagrangian-Eulerian              ',  //,
     &  ' coordinate system . . . . . . . . . . . . . (icoord)=',i10//,
     &  '    eq. 0, cartesian                                  ',  / ,
     &  '    eq. 1, axisymmetric                               ',  //,
     &  ' equation type . . . . . . . . . . . . . . . (navier)=',i10//,
     &  '    eq. 0, Euler (inviscid)                           ',  / ,
     &  '    eq. 1, Navier-Stokes (viscous)                    ',  //,
     &  ' restart option  . . . . . . . . . . . . . . (irs   )=',i10//,
     &  '    eq. 0, no restart nor solution written            ',  / ,
     &  '    eq. 1, restart written                            ',  / ,
     &  '    eq. 2, restart and solution written               ',  //,
     &  ' execution code  . . . . . . . . . . . . . . (iexec )=',i10//,
     &  '    eq. 0, data check only                            ',  / ,
     &  '    eq. 1, execution                                  ',  //,
     &  ' input echo parameter  . . . . . . . . . . . (necho )=',i10)
1300    format(a80,//,
     &  ' S o l u t i o n   P a r a m e t e r s                '   //,
     &  ' number of time sequences  . . . . . . . . . (ntseq )=',i10//,
     &  ' blocking algorithm  . . . . . . . . . . . . (imap  )=',i10//,
     &  '    eq. 0, ordered blocking                           ',  / ,
     &  '    eq. 1, disjoint element blocking                  ',  //,
     &  ' variational formulation . . . . . . . . . . (ivart )=',i10//,
     &  '    eq. 1, Galerkin                                   ',  / ,
     &  '    eq. 2, Galerkin/least-squares                     ',  / ,
     &  '    eq. 3, plus discontinuity-capturing operator      ',  //,
     &  ' discontinuity-capturing type  . . . . . . . (iDC   )=',i10//,
     &  '    eq. 1, DC-mallet                                  ',  / ,
     &  '    eq. 2, quadratic DC                               ',  / ,
     &  '    eq. 3, smallest of the previous two DCs           ',  //,
     &  ' dimension of Krylov space . . . . . . . . . (kspace)=',i10//,
     &  ' maximum number of GMRES cycles  . . . . . . (ngmres)=',i10)
1400    format(a80,//,
     &  ' S o l v e r   I n f o r m a t i o n                    ',//,
     &  ' Seq num    Nstep    Niter    Impl      Nupdate',
     &  '     Eps_Tol')
1500    format(i6,i10,i9,i8,i11,2x,e15.5)
1600    format(a80,//,
     &  ' S p a t i a l   a n d   T i m e   I n t e g r a t i o n',//,
     &  ' Seq num  Elem Int.  Bound Int.  Level  LCtime',
     &  '    Delt       CFLfld    CFLsld')
1700    format(i6,i8,i12,e13.4,1p,i8,1p,e13.4,0p,2f10.4)
1800    format(a80,//,
     &  ' O u t p u t   I n f o r m a t i o n                  ',1p,//,
     &  ' number of time steps per output . . . . . . (ntout )=',i10//,
     &  ' I/O format  . . . . . . . . . . . . . . . . (ioform)=',i10//,
     &  '    eq. 0, ASCII                                      ',  / ,
     &  '    eq. 1, binary                                     ',  //,
     &' scaling factor for density  . . . . . . . . (ro    )=',e15.5//,
     &' scaling factor for velocity . . . . . . . . (vel   )=',e15.5//,
     &' scaling factor for temperature. . . . . . . (temper)=',e15.5//,
     &' scaling factor for pressure . . . . . . . . (press )=',e15.5//,
     &' scaling factor for entropy  . . . . . . . . (entrop)=',e15.5)
c

1900    format(//,
     &  ' L e v e l   S e t   P a r a m e t e r s               '   //,
     &  ' Level Set Switch        . . . . . . . . . . (iLSet )=',i10//,
     &  '    eq. 0, No Level Set Solution Calculated            ',  / ,
     &  '    eq. 1, Level Set Calculated, 2 Fluid Props Read    ',  / ,
     &  '    eq. 2, Level Set and Redistancing Calcuations      ',  //,
     &  ' Property Smearing Band Width  . . . . . .(epsilon_ls)=',e15.5)



        end
