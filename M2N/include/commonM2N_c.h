// Routine contains the structures for reading the user input through
// input_fform.cpp. The default values for all these variables are defined in
// input.config.
//
// Input variables that have been previously declared in common.h have to be
// re-declared here, in a consistant structure.   

#include <FCMangle.h>

#define workfc FortranCInterface_GLOBAL_(workfc,WORKFC)
#define fronts FortranCInterface_GLOBAL_(fronts,FRONTS)
#define newdim FortranCInterface_GLOBAL_(newdim,NEWDIM)
#define timer4 FortranCInterface_GLOBAL_(timer4,TIMER4)
#define extrat FortranCInterface_GLOBAL_(extrat,EXTRAT)
#define spongevar FortranCInterface_GLOBAL_(spongevar,SPONGEVAR)
#define turbvar FortranCInterface_GLOBAL_(turbvar,TURBVAR)
#define turbvari FortranCInterface_GLOBAL_(turbvari,TURBVARI)
#define spebcvr FortranCInterface_GLOBAL_(spebcvr,SPEBCVR)
#define aerfrc FortranCInterface_GLOBAL_(aerfrc,AERFRC)
#define astore FortranCInterface_GLOBAL_(astore,ASTORE)
#define conpar FortranCInterface_GLOBAL_(conpar,CONPAR)
#define ctrlvari FortranCInterface_GLOBAL_(ctrlvari,CTRLVARI)
#define ctrlvar  FortranCInterface_GLOBAL_(ctrlvar,CTRLVAR)
#define shpdat FortranCInterface_GLOBAL_(shpdat,SHPDAT)
#define datpnt FortranCInterface_GLOBAL_(datpnt,DATPNT)
#define elmpar FortranCInterface_GLOBAL_(elmpar,ELMPAR)
#define genpar FortranCInterface_GLOBAL_(genpar,GENPAR)
#define inpdat FortranCInterface_GLOBAL_(inpdat,INPDAT)
#define intdat FortranCInterface_GLOBAL_(intdat,INTDAT)
#define mio FortranCInterface_GLOBAL_(mio,MIO)
#define mioname FortranCInterface_GLOBAL_(mioname,MIONAME)
#define itrpar FortranCInterface_GLOBAL_(itrpar,ITRPAR)
#define itrpnt FortranCInterface_GLOBAL_(itrpnt,ITRPNT)
#define matdat FortranCInterface_GLOBAL_(matdat,MATDAT)
#define mmatpar FortranCInterface_GLOBAL_(mmatpar,MMATPAR)
#define outpar FortranCInterface_GLOBAL_(outpar,OUTPAR)
#define point FortranCInterface_GLOBAL_(point,POINT)
#define precis FortranCInterface_GLOBAL_(precis,PRECIS)
#define propar FortranCInterface_GLOBAL_(propar,PROPAR)
#define resdat FortranCInterface_GLOBAL_(resdat,RESDAT)
#define solpar FortranCInterface_GLOBAL_(solpar,SOLPAR)
#define timdat FortranCInterface_GLOBAL_(timdat,TIMDAT)
#define timpar FortranCInterface_GLOBAL_(timpar,TIMPAR)
#define incomp FortranCInterface_GLOBAL_(incomp,INCOMP)
#define mtimer1 FortranCInterface_GLOBAL_(mtimer1,MTIMER1)
#define mtimer2 FortranCInterface_GLOBAL_(mtimer2,MTIMER2)
#define timer3 FortranCInterface_GLOBAL_(timer3,TIMER3)
#define title FortranCInterface_GLOBAL_(title,TITLE)
#define sclrs FortranCInterface_GLOBAL_(sclrs,SCLRS)
#define levlset FortranCInterface_GLOBAL_(levlset,LEVLSET)
#define nomodule FortranCInterface_GLOBAL_(nomodule,NOMODULE)
#define sequence FortranCInterface_GLOBAL_(sequence,SEQUENCE)
#define amgvarr FortranCInterface_GLOBAL_(amgvarr,AMGVARR)
#define amgvari FortranCInterface_GLOBAL_(amgvari,AMGVARI)

#define MAXBLK   50000
#define MAXSURF  30  
#define MAXTS   100
#define MAXTOP   5
#define MAXQPT   125
#define MAXSH    125
#define NSD      3
#define machin   'RS/6000'
#define machfl   4
#define zero   0.0000000000000000000000000000000d0
#define pt125   0.1250000000000000000000000000000d0
#define pt25   0.2500000000000000000000000000000d0
#define pt33   0.3333333333333333333333333333333d0
#define pt39   0.3968502629920498686879264098181d0
#define pt5   0.5000000000000000000000000000000d0
#define pt57   0.5773502691896257645091487805020d0
#define pt66   0.6666666666666666666666666666667d0
#define pt75   0.7500000000000000000000000000000d0
#define one   1.0000000000000000000000000000000d0
#define sqt2   1.4142135623730950488016887242097d0
#define onept5   1.5000000000000000000000000000000d0
#define two   2.0000000000000000000000000000000d0
#define three   3.0000000000000000000000000000000d0
#define four   4.0000000000000000000000000000000d0
#define five   5.0000000000000000000000000000000d0
#define pi   3.1415926535897932384626433832795d0

#ifdef __cplusplus
extern "C" {
#endif
  extern struct { 
    int master;
    int numpe;
    int myrank;
  } workfc ;

  extern struct { 
    int maxfront;
    int nlwork;
  } fronts ;

  extern struct { 
    int numper;
    int nshgt;
    int nshg0;
  } newdim ;

  extern struct { 
    double birth;
    double death;
    double comtim;
  } timer4 ;

  extern struct { 
    double ttim[100];
  } extrat ;

  extern struct {
    double zoutsponge, radsponge, zinsponge, grthosponge, grthisponge;
    double betamax;
    int spongecontinuity, spongemomentum1, spongemomentum2;
    int spongeenergy, spongemomentum3;
  } spongevar ;

  extern struct {
    double eles;
    double ylimit[9][3]; /* 9 = 5 + 4 = puvwT + 4Scalars */
    double rampmdot[3][2];
    double rmutarget;
    double pzero;
    double wtavei;
    double dtavei;
    double dke;
    double fwr1;
    double flump;
    double DES_SA_hmin;
    int ierrcalc;
    int ihessian;
    int itwmod;
    int ngaussf;
    int idim;
    int nlist;
    int nintf[MAXTOP];
  } turbvar ;

  extern struct {
    int irans, iles, idistcalc, isubmod;
    int ifproj;
    int i2filt;
    int modlstats;
    int idis;
    int nohomog;
    int ierrsmooth;
    int iramp;

/*      int itwmod; */
/*      double rtavei; */
/*      int ierrcalc; */
  } turbvari ;

  extern struct { 
    int irscale;
    int intpres;
    double plandist;
    double thetag;
    double ds;
    double tolerence;
    double radcyl;
    double rbltin;
    double rvscal;
  } spebcvr ;

  extern struct {
    double scdiff[5];
    double tdecay;
    int nsclr, isclr,nsolt, nosource;
    int consrv_sclr_conv_vel;
  } sclrs;

  extern struct { 
    double flxID[MAXSURF+1][10] ;
    double Force[3];
    double HFlux;
    int nsrflist[MAXSURF+1];
    int isrfIM;
    double flxIDsclr[MAXSURF][4];
  } aerfrc ;

  extern struct { 
    double a[100000];
  } astore ;

  extern struct { 
    int numnp;
    int numel;
    int numelb;
    int numpbc;
    int nen;
    int nfaces;
    int numflx;
    int ndof;
    int iALE;
    int icoord;
    int navier;
    int irs;
    int iexec;
    int necho;
    int ichem;
    int iRK;
    int nedof;
    int nshg;
    int nnz;
    int istop;
    int nflow;
    int nnz_tot;
    int idtn;
  } conpar ;
 
/*chen Sep 25 2009  Flow Control Parameters*/
  extern struct{
    int iI2Binlet;
    int isetOutPres;
    int isetInitial;
  } ctrlvari;

  extern struct{
    double inletVelX;
    double outPres1; 
    double xvel_ini;
    double yvel_ini;
    double zvel_ini;
    double temp_ini;
    double pres_ini;
    double evis_ini;
  } ctrlvar;
//////////////////////////////////////////

 
  extern struct { 
    double epsilon_ls;
    double epsilon_lsd;
    double dtlset;
    int iLSet;
    int ivconstraint;
    int iExpLSSclr1;
    int iExpLSSclr2;
  } levlset;

  extern struct { 
    int nshape;
    int nshapeb;
    int maxshb;
    int nshl;
    int nshlb;
    int nfath;
    int ntopsh;
    int nsonmax;
  } shpdat ;

  extern struct { 
    int mshp;
    int mshgl;
    int mwght;
    int mshpb;
    int mshglb;
    int mwghtb;
    int mmut;
    int mrhot;
    int mxst;
  } datpnt ;

  extern struct { 
    int lelCat;
    int lcsyst;
    int iorder;
    int nenb;
    int nelblk;
    int nelblb;
    int ndofl;
    int nsymdl;
    int nenl;
    int nfacel;
    int nenbl;
    int intind;
    int mattyp;
  } elmpar ;

  extern struct { 
    double E3nsd;
    int I3nsd;
    int nsymdf;
    int ndofBC;
    int ndiBCB;
    int ndBCB;
    int Jactyp;
    int jump;
    int ires;
    int iprec;
    int iprev;
    int ibound;
    int idiff;
    int lhs;
    int itau;
    int ipord;
    int ipred;
    int lstres;
    int iepstm;
    double dtsfct;
    double taucfct;
    int ibksiz;
    int iabc;
    int isurf;
    int idflx;
    double Bo;
    int EntropyPressure;
    int irampViscOutlet;
    int istretchOutlet;
    int iremoveStabTimeTerm;
  } genpar ;

  extern struct { 
    double epstol[6];  /* 1+ max number of scalars  (beginning of the
                          end of time sequences) */
    double Delt[MAXTS];
    double CFLfl[MAXTS];
    double CFLsl[MAXTS];
    int nstep[MAXTS];
    int niter[MAXTS];
    int impl[MAXTS];
    double rhoinf[MAXTS];
    int LHSupd[6];
    int loctim[MAXTS];
    double deltol[2][MAXTS];
  } inpdat ;

  extern struct { 
    int iin;
    int igeom;
    int ipar;
    int ibndc;
    int imat;
    int iecho;
    int iout;
    int ichmou;
    int irstin;
    int irstou;
    int ihist;
    int iflux;
    int ierror;
    int itable;
    int iforce;
    int igraph;
    int itime;
  } mio ;

  extern struct { 
    double fin;
    double fgeom;
    double fpar;
    double fbndc;
    double fmat;
    double fecho;
    double frstin;
    double frstou;
    double fhist;
    double ferror;
    double ftable;
    double fforce;
    double fgraph;
    double ftime;
  } mioname ;

  extern struct { 
    double eGMRES;
    int lGMRES;
    int iKs;
    int ntotGM;
  } itrpar ;

  extern struct { 
    int mHBrg;
    int meBrg;
    int myBrg;
    int mRcos;
    int mRsin;
  } itrpnt ;

  extern struct { 
    double datmat[MAXTS][7][3];
    int matflg[MAXTS][6];
    int nummat;
    int mexist;
  } matdat ;

  extern struct { 
    double pr, Planck, Stephan, Nh, Rh, Rgas;
    double gamma, gamma1, s0;
    //, const, xN2, xO2;
    //double yN2,    yO2,    Msh[5], cpsh[5],s0sh[5],h0sh[5];
    //double Rs[5],  cps[5], cvs[5], h0s[5], Trot[5],sigs[5];
    //double Tvib[5],g0s[5], dofs[5],ithm;
  } mmatpar ;

  extern struct { 
    double ro;
    double vel;
    double temper;
    double press;
    double entrop;
    int ntout;
    int ioform;
    int iowflux;
    int iofieldv;
    char iotype[80];
    int ioybar;
    int nstepsincycle;
    int nphasesincycle;
    int ncycles_startphaseavg;
//MR CHANGE
    int nsynciofiles;
    int nsynciofilesred;
//    int nsynciofieldsreadgeombc;
//    int nsynciofieldsreadrestart;
    int nsynciofieldswriterestart;
//MR CHANGE END
    /*  int iostats; */
/*      int ipresref; */
  } outpar ;

  extern struct { 
    int mbeg;
    int mend;
    int mprec;
  } point ;

  extern struct { 
    double epsM;
    int iabres;
  } precis ;

  extern struct { 
    int npro;
  } propar ;

  extern struct { 
    double resfrt;
  } resdat ;

  extern struct { 
    int imap;
    int ivart;
    int iDC;
    int iPcond;
    int Kspace;
    int nGMRES;
    int iconvflow;
    int iconvsclr;
    int idcsclr[2];
  } solpar ;

  extern struct { 
    double time;
    double CFLfld;
    double CFLsld;
    double Dtgl;
    double Dtmax;
    double alpha;
    double etol;
    int lstep;
    int ifunc;
    int itseq;
    int istep;
    int iter;
    int nitr;
    double almi;
    double alfi;
    double gami;
    double flmpl;
    double flmpr;
    double dtol[2];
    int iCFLworst;
    int lskeep;
  } timdat ;

  extern struct { 
    int LCtime;
    int ntseq;
  } timpar ;

  extern struct { 
    int numeqns[100];
    int minIters;
    int maxIters;
    int iprjFlag;
    int nPrjs;
    int ipresPrjFlag;
    int nPresPrjs;
    double prestol;
    double statsflow[6];
    double statssclr[6];
    int iverbose;
  } incomp ;

  extern struct { 
    double ccode[13];
  } mtimer1 ;

  extern struct { 
    double flops;
    double gbytes;
    double sbytes;
    int iclock;
    int icd;
    int icode;
    int icode2;
    int icode3;
  } mtimer2 ;

  extern struct { 
    double cpu[11];
    double cpu0[11];
    int nacess[11];
  } timer3 ;

  extern struct { 
    double title;
    int ititle;
  } title ;

  extern struct {
    int intg[MAXTS][2];
  }intdat;

  extern struct {
    double bcttimescale;    
    double ValueListResist[MAXSURF+1];
    double rhovw;
    double thicknessvw;
    double evw;
    double rnuvw;
    double rshearconstantvw;
    double betai;
    int icardio;
    int itvn;
    int ipvsq;
    int numResistSrfs;
    int nsrflistResist[MAXSURF+1];
    int numImpSrfs;
    int nsrflistImp[MAXSURF+1];
    int impfile;
    int numRCRSrfs;
    int nsrflistRCR[MAXSURF+1];
    int ircrfile;
    int ideformwall;  
    int iwallmassfactor;
    int iwallstiffactor;
    int iviscflux;   
 } nomodule;

  extern struct {
    int seqsize;
    int stepseq[100];
  } sequence;

  extern struct {
    double strong_eps;      /* strong criterion Stuben factor    */
    double ramg_eps;        /* AMG convergence eps               */
    double ramg_relax;       /* relaxation factor Gauss-Seidel/Jac*/
    double ramg_trunc;      /* truncation select */
    double ramg_chebyratio; /* Eigen ratio for chebyshev smoothing */
 } amgvarr ;
  
  extern struct {
    int irun_amg;           /* Employ AMG feature solfar.f      */
    int irun_amg_prec;      /* Run AMG as preconditioner to CG */
    int iamg_verb;          /* amg verbosity flag                */
    int iamg_neg_sten;      /* neg only stencil or neg and pos   */
    int iamg_nlevel;        /* number of levels 2-V etc.         */
    int iamg_c_solver;     /* solve fine level iter. method     */
    int iamg_init;           /* setup flag */
    int iamg_setup_frez;    /* how many solfars to re setup amg */
    int iamg_interp;        /* interpolation select */
    int maxnev;             /* total eigenvectors used for ggb*/
    int maxncv;             /* total iterative vectors for ggb*/
    int iamg_smoother;      /* Smoother type */
    int mlsdeg;             /* Polynomial Smoothing (MLS) degree */
    int iamg_reduce;        /* Run a reduced case */
 } amgvari ;

#ifdef __cplusplus
}
#endif
