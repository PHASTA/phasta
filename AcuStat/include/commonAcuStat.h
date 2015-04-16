c----------------------------------------------------------------------
c
c This file contains the common blocks and the data declaration needed
c for the routines.
c
c Input variables that have been previously declared in common_c.h have to be
c re-declared here, in a consistant block. 
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------

	IMPLICIT REAL*8 (a-h,o-z)
c
c.... parameters  IF YOU CHANGE THES YOU HAVE TO CHANGE THEM IN
c                  common_c.h ALSO
c
        parameter     ( MAXBLK = 50000, MAXTS = 100)
        parameter     ( MAXSH = 32, NSD = 3 )
c
c  The five types of region topology are  1= Tet, 2=Hex, 3= Wedge (tri-start),
c                                         4= Wedge (quad-first) 5=pyramid
c
c  The two types of face topology are  1= tri, 2=quad
c
        parameter     ( MAXTOP = 6, MAXSURF=30 )
 
c the common block nomodule holds all the things which have been removed
c from different modules
     
        integer seqsize, stepseq
        integer consrv_sclr_conv_vel
        integer spongecontinuity, spongemomentum1, spongemomentum2
        integer spongeenergy, spongemomentum3
	common /amgvarr/strong_eps,ramg_eps,ramg_relax,ramg_trunc,
     &              ramg_chebyratio
	common /amgvari/irun_amg,irun_amg_prec,
     &                  iamg_verb,
     &                  iamg_neg_sten,iamg_nlevel,
     &                  iamg_c_solver,
     &                  iamg_init,
     &        iamg_setup_frez,
     &        iamg_interp,maxnev,maxncv,iamg_smoother,mlsdeg,
     &        iamg_reduce

        common /nomodule/ bcttimescale,ValueListResist(0:MAXSURF),
     &  rhovw,thicknessvw, evw, rnuvw, rshearconstantvw, betai,
     &  icardio, itvn, ipvsq, numResistSrfs, nsrflistResist(0:MAXSURF),
     &  numImpSrfs, nsrflistImp(0:MAXSURF),impfile,
     &  numRCRSrfs, nsrflistRCR(0:MAXSURF),ircrfile,
     &  ideformwall, iwallmassfactor, iwallstiffactor, iviscflux 
        common /sequence/ seqsize, stepseq(100)
	common /workfc/ master, numpe, myrank
	common /fronts/ maxfront, nlwork
	common /newdim/ numper, nshgt, nshg0
	common /timer4/ birth, death, comtim
        common /extrat/ ttim(100)
        common /spongevar/ zoutSponge, radSponge, zinSponge,
     &            grthOSponge,grthISponge,betamax,
     &            spongecontinuity, spongemomentum1, spongemomentum2,
     &            spongeenergy, spongemomentum3
        common /turbvar/ eles,ylimit(3,9), rampmdot(2,3),
     &                   rmutarget, pzero,  wtavei, 
     &                   dtavei, dke,  fwr1, flump, DES_SA_hmin,
     &                   ierrcalc, ihessian, itwmod, ngaussf,idim,
     &                   nlist, nintf(MAXTOP)
        common /turbvari/iRANS, iLES, idistcalc, isubmod, ifproj,
     &                   i2filt, modlstats, idis, nohomog,
     &                   ierrsmooth, iramp
        common /spebcvr/ irscale, intpres, plandist,
     &            thetag, ds, tolerence, radcyl, rbltin, rvscal

        common /sclrs/ scdiff(5),tdecay,nsclr,isclr,nsolt,nosource,
     &            consrv_sclr_conv_vel
c 
c.... common blocks
c
      parameter (MAXQPT = 125)
c
c.... common blocks for hierarchic basis functions
c
      common /intpt/  Qpt (MAXTOP ,4,MAXQPT), Qwt (MAXTOP ,MAXQPT), 
     &                Qptb(MAXTOP,4,MAXQPT),  Qwtb(MAXTOP,MAXQPT), 
     &                nint(MAXTOP),           nintb(MAXTOP),
     &                ngauss,                 ngaussb,   intp,
     &                   maxnint
 
c nsrflist is a binary switch that tells us if a given srfID should be
c included in the consistent flux calculation.  It starts from zero
c since we need to be able to handle/ignore surfaces with no srfID attached
c
c flxID(numfluxes,nIDs+1)
c numfluxes = area, mass, fx, fy, fz, heat, scalar_flux_{1,2,3,4}
c nIDs currently set to MAXSURF, each surface has its own
c
        common /aerfrc/ flxID(10,0:MAXSURF), Force(3),HFlux, 
     &                  nsrflist(0:MAXSURF), isrfIM, 
     &                  flxIDsclr(4,MAXSURF)
c
        common /astore/ a(100000)
c
        common /blkdat/ lcblk(10,MAXBLK+1),      lcblkb(10,MAXBLK+1)
c
        common /mbndnod/ mnodeb(9,8,3)
c
        common /conpar/ numnp,  numel,  numelb, numpbc, nen,    nfaces,
     &                  numflx, ndof,   iALE,   icoord, navier,
     &                  irs,    iexec,  necho,  ichem,  iRK,    nedof,
     &                  nshg,   nnz,    istop,  nflow,  nnz_tot, idtn

c...........................................................................
        common /ctrlvari/ iI2Binlet, isetOutPres, isetInitial

        real*8 inletVelX
        common /ctrlvar/  inletVelX,   outPres1, 
     &                    xvel_ini,    yvel_ini,    zvel_ini,
     &                    temp_ini,    pres_ini,    evis_ini
c...........................................................................

c
        common /levlset/ epsilon_ls, epsilon_lsd, dtlset, iLSet, 
     &                   ivconstraint, iExpLSSclr1, iExpLSSclr2

c 
        common /shpdat/ nshape, nshapeb, maxshb,
     &                  nshl, nshlb,nfath,  ntopsh,  nsonmax
c
        common /datpnt/ mshp,   mshgl,  mwght,  mshpb,  mshglb, mwghtb,
     &                  mmut,   mrhot,  mxst
c
        common /melmcat/ mcsyst, melCat, nenCat(8,3),    nfaCat(8,3)
c
        common /elmpar/ lelCat, lcsyst, iorder, nenb,   
     &                  nelblk, nelblb, ndofl,  nsymdl, nenl,   nfacel,
     &                  nenbl,  intind, mattyp
c

        integer EntropyPressure

        common /genpar/ E3nsd,  I3nsd,  nsymdf, ndofBC, ndiBCB, ndBCB,
     &                  Jactyp, jump,   ires,   iprec,  iprev,  ibound,
     &                  idiff,  lhs,    itau,   ipord,  ipred,  lstres,
     &                  iepstm, dtsfct, taucfct, ibksiz, iabc, isurf,
     &                  idflx,  Bo,     EntropyPressure, irampViscOutlet,
     &                  istretchOutlet, iremoveStabTimeTerm

c
        common /inpdat/ epstol(6),  Delt(MAXTS),    CFLfl(MAXTS),
     &                  CFLsl(MAXTS),   nstep(MAXTS),   niter(MAXTS),
     &                  impl(MAXTS),    rhoinf(MAXTS),
     &                  LHSupd(6),  loctim(MAXTS),  deltol(MAXTS,2)
c
        common /intdat/ intg(2,MAXTS),  intpt(3),       intptb(3)
c
        common /mintpar/ indQpt(3,3,4),  numQpt(3,3,4),
     &                  intmax
c
        common /mio    / iin,    igeom,  ipar,   ibndc,  imat,   iecho,
     &                  iout,   ichmou, irstin, irstou, ihist,  iflux,
     &                  ierror, itable, iforce, igraph, itime
c
c /*         common /andres/ fwr1,ngaussf,idim,nlist */

        character*80    fin,    fgeom,  fpar,   fbndc,  fmat,   fecho,
     &                  frstin, frstou, fhist,  ferror, ftable, fforce,
     &                  fgraph, ftime,  iotype     
        common /mioname/ fin,    fgeom,  fpar,   fbndc,  fmat,   fecho,
     &                  frstin, frstou, fhist,  ferror, ftable, fforce,
     &                  fgraph, ftime
c
        common /itrpar/ eGMRES, lGMRES, iKs,    ntotGM
c
        common /itrpnt/ mHBrg,  meBrg,  myBrg,  mRcos,  mRsin
c
        REAL*8          Nh, Msh
        common /mmatpar/ pr,     Planck, Stefan, Nh,     Rh,     Rgas,
     &                  gamma,  gamma1, s0,     const,  xN2,    xO2,
     &                  yN2,    yO2,    Msh(5), cpsh(5),s0sh(5),h0sh(5),
     &                  Rs(5),  cps(5), cvs(5), h0s(5), Trot(5),sigs(5),
     &                  Tvib(5),g0s(5), dofs(5),ithm
c
        logical         mexist
        common /matdat/ datmat(3,7,MAXTS),      matflg(6,MAXTS),
     &                  nummat,                 mexist
c
cMR CHANGE
c         common /outpar/ ro,     vel,    temper, press,  entrop, ntout,
c      &                  ioform, iowflux, iofieldv, iotype, ioybar,
c      &                  nstepsincycle, nphasesincycle, ncycles_startphaseavg
        common /outpar/ ro,     vel,    temper, press,  entrop, ntout,
     &                  ioform, iowflux, iofieldv, iotype, ioybar,
     &                  nstepsincycle, nphasesincycle, 
     &                  ncycles_startphaseavg,
     &                  nsynciofiles, nsynciofieldswriterestart
cMR CHANGE END

c
        common /point / mbeg,   mend,   mprec
c
        common /precis/ epsM,   iabres
c
        common /propar/ npro
c
        common /resdat/ resfrt
c
        common /solpar/ imap,   ivart,  iDC,    iPcond, Kspace, nGMRES,
     &                  iconvflow, iconvsclr, idcsclr(2)
c
        common /msympar/ indsym(5,5)
c
        common /timdat/ time,   CFLfld, CFLsld, Dtgl,   Dtmax,  alpha,
     &                  etol,   lstep,  ifunc,  itseq,  istep,  iter,
     &                  nitr,   almi,   alfi,   gami,   flmpl,  flmpr,
     &                  dtol(2), iCFLworst, lskeep
c
        common /timpar/ LCtime, ntseq
c
        common /incomp/ numeqns(100), minIters, maxIters, 
     &                  iprjFlag,     nPrjs,    ipresPrjFlag, nPresPrjs,
     &                  prestol,      statsflow(6), statssclr(6),
     &                  iverbose
c
        character*8     ccode
        common /mtimer1/ ccode(13)
c
        integer       flops,  gbytes, sbytes
        common /mtimer2/ flops,  gbytes, sbytes, iclock, icd,    icode,
     &                  icode2, icode3
c
        common /timer3/ cpu(11),        cpu0(11),       nacess(11)
c
        character*80    title,  ititle
        common /title / title,  ititle
c
        character*8     machin
        parameter     ( machin = 'RS/6000 ' )
        parameter     ( machfl = 4 )
 
        parameter 
     &           ( zero   = 0.0000000000000000000000000000000d0,
     &             pt125  = 0.1250000000000000000000000000000d0,
     &             pt25   = 0.2500000000000000000000000000000d0,
     &             pt33   = 0.3333333333333333333333333333333d0,
     &             pt39   = 0.3968502629920498686879264098181d0,
     &             pt5    = 0.5000000000000000000000000000000d0,
     &             pt57   = 0.5773502691896257645091487805020d0,
     &             pt66   = 0.6666666666666666666666666666667d0,
     &             pt75   = 0.7500000000000000000000000000000d0,
     &             one    = 1.0000000000000000000000000000000d0,
     &             sqt2   = 1.4142135623730950488016887242097d0,
     &             onept5 = 1.5000000000000000000000000000000d0,
     &             two    = 2.0000000000000000000000000000000d0,
     &             three  = 3.0000000000000000000000000000000d0,
     &             four   = 4.0000000000000000000000000000000d0,
     &             five   = 5.0000000000000000000000000000000d0,
     &             pi     = 3.1415926535897932384626433832795d0)

c
c----------------------------------------------------------------------
c
c.... element pointers
c
c mmat   (MAXBLK)  : pointer to interior element material number
c mmatb  (MAXBLK)  : pointer to boundary element material number
c mien   (MAXBLK)  : pointer to ien array
c mienb  (MAXBLK)  : pointer to ienb array
c miBCB  (MAXBLK)  : pointer to iBCB array
c mDt    (MAXBLK)  : pointer to Dt array
c mDC    (MAXBLK)  : pointer to DC array
c mBCB   (MAXBLK)  : pointer to BCB array
c mstiff (MAXBLK)  : pointer to stiff array
c
c----------------------------------------------------------------------
c
c.... common /aerfrc/   : aerodynamic forces
c
c Force(3)      : components of the aerodynamic forces
c HFlux         : total heat flux
c
c----------------------------------------------------------------------
c
c.... common /astore/   : the dynamic memory allocation area
c
c a(...)        : the blank array used for front-end data storage
c
c----------------------------------------------------------------------
c
c.... common /blkdat/   : blocking data
c
c lcblk  (10,MAXBLK+1) : blocking data for the interior elements
c lcblkb (10,MAXBLK+1) : blocking data for the boundary elements
c
c----------------------------------------------------------------------
c
c.... common /bndnod/   : boundary nodes of boundary elements
c
c mnodeb (9,8,3) : boundary nodes of each element category and dimension
c
c----------------------------------------------------------------------
c
c.... common /conpar/   : input constants
c
c numnp         : number of nodal points
c numel         : number of elements
c numelb        : number of boundary elements
c numpbc        : number of nodes having a boundary condition
c nen           : maximum number of element nodes
c nfaces        : maximum number of element faces
c nsd           : number of space dimensions
c numflx        : number of flux boundary nodes
c ndof          : number of degrees of freedom per node
c iALE          : ALE formulation flag
c icoord        : coordinate system flag
c navier        : Navier-Stokes calculation flag
c irs           : restart option 
c iexec         : execute flag
c necho         : input echo parameter
c ichem         : equilibrium chemistry flag (for outchem.step dump)
c iRK           : Runge-Kutta flag
c nshg          : global number of shape functions (degrees of freedom,
c                 or equations). Computed from the specified p-order,
c                 the number of edges, and the number of faces (in the
c                 entire mesh)
c
c----------------------------------------------------------------------
c
c.... common /datpnt/   : front-end data pointers
c
c mshp          : pointer to shape-functions 
c mshgl         : pointer to local-grad-shape-functions
c mwght         : pointer to quadrature weights
c mshpb         : pointer to shape-functions of boundary elements
c mshglb        : pointer to local-grad-shape-functions of bound. elem.
c mwghtb        : pointer to quadrature weights of bound. elements
c mmut          : pointer to table mu  = mu  (p,T)
c mrhot         : pointer to table rho = rho (p,T)
c mxst          : pointer to table xs  = xs  (p,T)
c
c----------------------------------------------------------------------
c
c.... common /elmcat/   : element category information
c
c mcsyst        : maximum number of element coordinate system
c melCat        : maximum number of element categories
c nenCat (8,3)  : number of nodes for each category and dimension
c nfaCat (8,3)  : number of faces for each category and dimension
c
c----------------------------------------------------------------------
c
c.... common /elmpar/   : element parameters
c
c lelCat        : element category (P1, Q1, P2, Q2, etc.)
c lcsyst        : element coordinate system
c iorder        : element order (=k for Pk and Qk)
c nenb          : number of element nodes per boundary sides
c maxsh         : total number integration points
c maxshb        : total number integration points of boundary elements
c nelblk        : number of element blocks
c nelblb        : number of boundary element blocks
c ndofl         : number of degrees of freedom (for current block)
c nsymdl        : number of d.o.f for symm. storage (for current block)
c nenl          : number of element nodes (for current block)
c nfacel        : number of element faces (for current block)
c nenbl         : number of boundary element nodes
c intind        : integration data index
c nintg         : number of integration points
c mattyp        : material type ( = 0 for fluid; = 1 for solid )
c
c----------------------------------------------------------------------
c
c.... common /genpar/   : control parameters
c
c E3nsd         : NSD .eq. 3 flag; 0. for 2D, 1. for 3D
c I3nsd         : NSD .eq. 3 flag; 0  for 2D, 1  for 3D
c nsymdf        : number of d.o.f.'s in symm. storage (= ndof*(ndof+1)/2)
c ndofBC        : dimension size of the boundary condition array BC
c ndiBCB        : dimension size of the boundary condition array iBCB
c ndBCB         : dimension size of the boundary condition array BCB
c Jactyp        : Jacobian type flag
c jump          : jump term computation flag
c ires          : residual type computation flag
c iprec         : block-diagonal preconditioner flag
c iprev         : ypl array allocation flag
c ibound        : boundary element flag
c idiff         : diffusive flux vector flag
c                 ( = 0 not used; = 1 global reconstruction )
c itau          : type of tau to be used
c
c----------------------------------------------------------------------
c
c.... common /inpdat/   : time sequence input data
c
c epstol (MAXTS)  : tolerance for GMRES solvers
c Delt   (MAXTS)  : global time step
c CFLfl  (MAXTS)  : CFL number for fluid flow
c CFLsl  (MAXTS)  : CFL number for structural heating
c nstep  (MAXTS)  : number of time steps
c niter  (MAXTS)  : number of iterations per time step
c impl   (MAXTS)  : solver flag
c iturb  (MAXTS)  : turbulence model flag
c rhoinf (MAXTS)  : time integration spectral radius paramter
c                             (0=Gears       1= trapezoidal rule)
c LHSupd (MAXTS)  : LHS/preconditioner update
c loctim (MAXTS)  : local time stepping flag
c
c----------------------------------------------------------------------
c
c.... common /intdat/   : integration data
c
c intg  (2,MAXTS) : integration parameters
c intpt (3)       : integration pointers
c intptb(3)       : integration pointers of boundary elements
c
c----------------------------------------------------------------------
c
c.... common /shpdat/   : hierarchic shape function quadrature data
c
c Qpt  (3,MAXQPT)  : interior element quadrature points (xi,eta,zeta)
c Qwt  (MAXQPT)    : interior element quad. weights
c Qptb (2,MAXQPT)  : boundary element quad. pnts.
c Qwtb (MAXQPT)    : boundary element quad. weights
c nshape           : number of interior element shape functions
c nshapeb          :   "    "  boundary  "        "       "
c ngauss           : number of interior element integration points
c ngaussb          :   "    "  boundary  "        "           "
c----------------------------------------------------------------------
c
c.... common /intpar/   : integration parameters
c
c Qpt   (4,*)   : xi, eta, zeta, weight of quadrature points
c indQpt(3,3,4) : index to quadrature points for a given rule
c numQpt(3,3,4) : number of quadrature points for a given rule
c intmax        : number of allowable spatial integ. points per nsd
c
c----------------------------------------------------------------------
c
c.... common /io    /   : io channels
c
c iin           : input  (main parameters)          [INPUT.DAT]
c igeom         : input  (problem geometry)         [GEOM.DAT]
c ipar          : in/out (spectral mapping)         [PARTITION.DAT]
c ibndc         : input  (problem boundary cond.)   [BC.DAT]
c imat          : input  (element material types)   [MATERIAL.DAT]
c iecho         : output (echo of input)            [ECHO.DAT]
c iout          : output (result output)            [OUTPUT.lstep]
c ichmou        : output (chemistry output)         [OUTCHM.lstep]
c irstin        : input  (input restart)            [RESTAR.INP]
c irstou        : output (output restart)           [RESTAR.OUT]
c ihist         : output (history output)           [HISTOR.DAT]
c iflux         : output (boundary flux)            [FLUX.lstep]
c ierror        : output (error messages)           [ERROR.DAT]
c itable        : input  (equilibrium chemistry)    [TABLE.DAT]
c iforce        : output (aerodynamic forces)       [FORCES.DAT]
c
c----------------------------------------------------------------------
c
c.... common /ioname/   : io file names
c
c fin           : input.dat
c fgeom         : geom.dat
c fpar          : partition.dat
c fbndc         : bc.dat
c fmat          : material.dat
c fecho         : echo.dat
c frstin        : restar.inp
c frstou        : restar.out
c fhist         : histor.dat
c ferror        : error.dat
c ftable        : table.dat
c fforce        : forces.dat
c
c----------------------------------------------------------------------
c
c.... common /itrpar/   : Preconditioned GMRES parameters
c
c eGMRES        : finite difference interval
c lGMRES        : number of GMRES cycles
c iKs           : current Krylov vector
c ntotGM        : total number of GMRES iterations
c
c----------------------------------------------------------------------
c
c.... common /itrpnt/   : Preconditioned GMRES array pointers
c
c mHBrg         : pointer to Hessenberg matrix
c meBrg         : pointer to Hessenberg's RHS matrix
c myBrg         : pointer to minimize solution matrix
c mRcos         : pointer to Rotation Cosine of QR algorithm
c mRsin         : pointer to Rotation Sine   of QR algorithm
c
c----------------------------------------------------------------------
c
c.... common /matpar/   : material constants
c
c pr            : Prandtl number
c Planck        : Planck's constant
c Stefan        : Stefan's constant (for radiation)
c Nh            : Avogadro's number
c Rh            : universal gas constant
c Rgas          : specific gas constant
c gamma         : specific heat ratio
c gamma1        : gamma - 1
c s0            : reference specific entropy
c const         : special constant
c xN2           : mole fraction of diatomic nitrogen
c xO2           : mole fraction of diatomic oxygen
c yN2           : mole fraction of diatomic nitrogen
c yO2           : mole fraction of diatomic oxygen
c Msh  (5)      : molar mass of species
c cpsh (5)      : molar heat at constant pressure of species
c s0sh (5)      : molar reference entropy of species
c h0sh (5)      : molar heat of formation of species
c Rs   (5)      : specific gas constant of species
c cps  (5)      : specific heat at constant pressure of species
c cvs  (5)      : specific heat at constant volume of species
c h0s  (5)      : specific heat of formation of species
c Trot (5)      : characteristic rotational temperature of species
c sigs (5)      : symmetry factor of species
c Tvib (5)      : characteristic vibrational temperature of species
c g0s  (5)      : ground degeneracy of electronic energy
c dofs (5)      : degrees of freedom of species
c ithm          : thermodynamic property flag
c
c----------------------------------------------------------------------
c
c.... common /matdat/   : material data
c
c datmat (3,5,2) : material data
c matflg (5,100)   : material type flag
c nummat           : number of materials
c mexist           : flag indicating the presence of MATERIAL.DAT
c
c----------------------------------------------------------------------
c
c.... common /outpar/   : output parameters
c
c ro            : density     rescaling factor for output
c vel           : velocity    rescaling factor for output
c temper        : temperature rescaling factor for output
c press         : pressure    rescaling factor for output
c entrop        : entropy     rescaling factor for output
c ntout         : number of steps between consecutive printouts
c ioform        : output I/O format
c
c----------------------------------------------------------------------
c
c.... common /point /   : dynamic storage pointer management data
c
c mbeg          : pointer to the beginning of the free storage
c mend          : pointer to the end of the storage
c mprec         : precision of the floating point data
c
c----------------------------------------------------------------------
c
c.... common /precis/   : finite difference interval data
c
c epsM          : square root of machine precision
c iabres        : absolute value residual flag
c
c----------------------------------------------------------------------
c
c....common /propar/    : processor related information
c
c npro          : number of virtual processors for the current block
c
c----------------------------------------------------------------------
c
c....common /resdat/    : residual statistics data
c
c resfrt        : first residual of convergence
c
c----------------------------------------------------------------------
c
c.... common /solpar/   : solution parameters
c
c imap          : permutation mapping flag
c ivart         : variational formulation type
c iDC           : DC type
c iPcond        : type of preconditioner
c Kspace        : dimension of Krylov space
c nGMRES        : maximum number of GMRES iterations
c
c----------------------------------------------------------------------
c
c.... common /sympar/   : symmetric storage parameters
c
c indsym (5,5)  : mapping from 2D storage to symmetric one
c
c----------------------------------------------------------------------
c
c.... common /timdat/   : time data
c
c time          : current run time
c CFLfld        : CFL number for fluid flow
c CFLsld        : CFL number for structural heating
c Dtgl          : inverse of global time step
c Dtmax         : maximum delta-time
c alpha         : trapezoidal rule parameter
c etol          : epsilon tolerance for GMRES
c lstep         : current time step
c ifunc         : func. eval. counter (=niter*(lstep-lstep0) + iter)
c itseq         : sequence number
c istep         : step number (reseted at the beginning of the run)
c iter          : iteration number
c nitr          : number of multi-corrector iterations for this sequence
c
c----------------------------------------------------------------------
c
c.... common /timpar/   : time integration parameters
c
c LCtime        : local time stepping flag
c ntseq         : number of time sequences
c
c----------------------------------------------------------------------
c
c.... common /timer1/   : timer parameters
c.... common /timer2/   : timer parameters
c.... common /timer3/   : timer parameters
c
c ccode(13)     : timing entities codes
c flops         : flop counter
c gbytes        : byte counter for gather operation
c sbytes        : byte counter for scatter operation
c iclock        : wall-clock time (in milliseconds)
c icd           : number of timing entities
c icode         : current timer code
c icode2        : last timer code
c icode3        : next-to-last timer code
c cpu(11)       : cpu time of each entity
c cpu0(11)      : initial cpu time of each entity
c nacess(11)    : number of times each entity is accessed
c
c----------------------------------------------------------------------
c
c.... common /title /   : problem title
c
c title         : problem title
c ititle        : problem title (with form feed)
c
c----------------------------------------------------------------------
c
c.... common /avging / : nfath
c 
c nfath         : total number of global fathers over which certain
c                 quantities will be averaged
c 
c----------------------------------------------------------------------
c
c.... parameters        : machine data
c
c machin        : machine type
c                  (set parameter)
c machfl        : single precision floating point lenght in bytes
c                  (set parameter)
c
c----------------------------------------------------------------------
c
c.... parameters        : useful constants
c
c zero          : 0.0
c pt125         : 0.125
c pt25          : 0.25
c pt33          : 0.33 (1/3)
c pt39          : 2^(-4/3)
c pt5           : 0.5
c pt57          : 1/sqrt(3)
c pt66          : 0.66 (2/3)
c pt75          : 0.75
c one           : 1.0
c sqt2          : sqrt(2)
c onept5        : 1.5
c two           : 2.0
c three         : 3.0
c four          : 4.0
c five          : 5.0
c pi            : the magical number :-)
c 
c---------------------------------------------------------------------- 
c 
c Zdenek Johan, Winter 1991.
c 
c----------------------------------------------------------------------
