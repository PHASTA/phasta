#ifdef HAVE_PETSC
#include <stdio.h>

#include "petscsys.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscksp.h"

#include "assert.h"
#include "common_c.h"

#include <FCMangle.h>
#define SolGMRp FortranCInterface_GLOBAL_(solgmrp,SOLGMRP)
#define SolGMRpSclr FortranCInterface_GLOBAL_(solgmrpsclr,SOLGMRPSCLR)
#define ElmGMRPETSc FortranCInterface_GLOBAL_(elmgmrpetsc, ELMGMRPETSC)
#define ElmGMRPETScSclr FortranCInterface_GLOBAL_(elmgmrpetscsclr, ELMGMRPETSCSCLR)
#define rstatp FortranCInterface_GLOBAL_(rstatp, RSTATP)
#define rstatpSclr FortranCInterface_GLOBAL_(rstatpsclr, RSTATPSCLR)
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define get_time FortranCInterface_GLOBAL_(get_time,GET_TIME) 
#define get_max_time_diff FortranCInterface_GLOBAL_(get_max_time_diff,GET_MAX_TIME_DIFF) 

typedef long long int gcorp_t;

void get_time(uint64_t* rv, uint64_t* cycle);
void get_max_time_diff(uint64_t* first, uint64_t* last, uint64_t* c_first, uint64_t* c_last, char* lbl);

//#include "auxmpi.h"
//      static PetscOffset poff;

      static Mat lhsP;
      static PC pc;
      static KSP ksp;
      static Vec DyP, resP, DyPLocal;
      static PetscErrorCode ierr;
      static PetscInt PetscOne, PetscRow, PetscCol, LocalRow, LocalCol;
      static IS LocalIndexSet;
      static ISLocalToGlobalMapping VectorMapping;
      static  VecScatter scatter7;
      static int firstpetsccall = 1;
   
      static Mat lhsPs;
      static PC pcs;
      static KSP ksps;
      static Vec DyPs, resPs, DyPLocals;
      static IS LocalIndexSets;
      static ISLocalToGlobalMapping VectorMappings;
      static VecScatter scatter7s;
      static int firstpetsccalls = 1;

      static int rankdump=-1;  // 8121 was the problem rank with 3.5.3
      PetscReal resNrm;

   
void     SolGMRp(double* y,         double* ac,        double* yold,      
     	double* x,         int* iBC,       double* BC,  
     	int* col,       int* row,       double* lhsk,         
     	double* res,       double* BDiag,     
     	int* iper,
     	int* ilwork,    double* shp,       double* shgl,      double* shpb,
     	double* shglb,     double* Dy,        double* rerr, long long int* fncorp)
{
// 
// ----------------------------------------------------------------------
// 
//  This is the preconditioned GMRES driver routine.
// 
// input:
//  y      (nshg,ndof)           : Y-variables at n+alpha_v
//  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
//  yold   (nshg,ndof)           : Y-variables at beginning of step
//  acold  (nshg,ndof)           : Primvar. accel. variable at begng step
//  x      (numnp,nsd)            : node coordinates
//  iBC    (nshg)                : BC codes
//  BC     (nshg,ndofBC)         : BC constraint parameters
//  shp(b) (nen,maxsh,melCat)     : element shape functions (boundary)
//  shgl(b)(nsd,nen,maxsh,melCat) : local gradients of shape functions
// 
// output:
//  res    (nshg,nflow)           : preconditioned residual
//  BDiag  (nshg,nflow,nflow)      : block-diagonal preconditioner
// 
//  
// Zdenek Johan,  Winter 1991.  (Fortran 90)
// ----------------------------------------------------------------------
// 


// Get variables from common_c.h
      int nshg, nflow, nsd, iownnodes;
      nshg  = conpar.nshg; 
      nflow = conpar.nflow; 
      nsd = NSD; 
      iownnodes = conpar.iownnodes;

      gcorp_t nshgt;
      gcorp_t mbeg;
      gcorp_t mend;
      nshgt = (gcorp_t) newdim.nshgt; //Fix nshgt in common_c.h
      mbeg = (gcorp_t) newdim.minowned; 
      mend = (gcorp_t) newdim.maxowned; 


      int node, element, var, eqn;
      double valtoinsert;
      int nenl, iel, lelCat, lcsyst, iorder, nshl;
      int mattyp, ndofl, nsymdl, npro, ngauss, nppro;
      double DyFlat[nshg*nflow];
      double DyFlatPhasta[nshg*nflow];
      double rmes[nshg*nflow];
// DEBUG
      int i,j,k,l,m;

      // FIXME: PetscScalar
      double  real_rtol, real_abstol, real_dtol;
// /DEBUG
      //double parray[1]; // Should be a PetscScalar 
      double *parray; // Should be a PetscScalar 
      PetscInt petsc_bs, petsc_m, petsc_M,petsc_PA;
      PetscInt petsc_n;
      PetscOne = 1;
      uint64_t duration[4];

      gcorp_t glbNZ;

      if(firstpetsccall == 1) {
//Everthing in this conditional block should be moved to a function to estimate the size of PETSc's matrix which improves time on the first matrix fill
//
        PetscInt* idiagnz= (PetscInt*) malloc(sizeof(PetscInt)*iownnodes);
        PetscInt* iodiagnz= (PetscInt*) malloc(sizeof(PetscInt)*iownnodes);
        for(i=0;i<iownnodes;i++) {
             idiagnz[i]=0;
             iodiagnz[i]=0;
        }
        i=0;
        for(k=0;k<nshg;k++) {
          if((fncorp[k] < mbeg) || (fncorp[k] >mend)){
// this node is not owned by this rank so we skip 
          } else { 
           for(j=col[i]-1;j<col[i+1]-1;j++) {
//          assert(row[j]<=nshg);
//          assert(fncorp[row[j]-1]<=nshgt);
             glbNZ=fncorp[row[j]-1];
             if((glbNZ < mbeg) || (glbNZ > mend)) {
                iodiagnz[i]++;
             } else {
                idiagnz[i]++;
             }
           }
           i++; 
          }
        }
        gcorp_t mind=1000;
        gcorp_t mino=1000;
        gcorp_t maxd=0;
        gcorp_t maxo=0;
        for(i=0;i<iownnodes;i++) {
           mind=min(mind,idiagnz[i]);
           mino=min(mino,iodiagnz[i]);
           maxd=max(maxd,idiagnz[i]);
           maxo=max(maxo,iodiagnz[i]);
//           iodiagnz[i]=max(iodiagnz[i],10);
//           idiagnz[i]=max(idiagnz[i],10);
//           iodiagnz[i]=2*iodiagnz[i];   //  estimate a bit higher for off-part interactions
//           idiagnz[i]=2*idiagnz[i];   //  estimate a bit higher for off-part interactions
        }
// the above was pretty good but below is faster and not too much more memory...of course once you do this 
// could just use the constant fill parameters in create but keep it alive for potential later optimization

        for(i=0;i<iownnodes;i++) {
           iodiagnz[i]=1.3*maxd;
           idiagnz[i]=1.3*maxd;
        }


        
        if(workfc.numpe < 200){
          printf("myrank,i,iownnodes,nshg %d %d %d %d \n",workfc.myrank,i,iownnodes,nshg);
          printf("myrank,mind,maxd,mino,maxo %d %d %d %d %d \n",workfc.myrank,mind,maxd,mino,maxo);
        }
        // Print debug info
        if(nshgt < 200){
          int irank;
          for(irank=0;irank<workfc.numpe;irank++) {
            if(irank == workfc.myrank){
              printf("mbeg,mend,myrank,idiagnz, iodiagnz %d %d %d \n",mbeg,mend,workfc.myrank);
              for(i=0;i<iownnodes;i++) {
                printf("%d %ld %ld \n",workfc.myrank,idiagnz[i],iodiagnz[i]);
              }
            }
           MPI_Barrier(MPI_COMM_WORLD); 
          }
        } 
        petsc_bs = (PetscInt) nflow;
        petsc_m  = (PetscInt) nflow* (PetscInt) iownnodes;
        petsc_M  = (PetscInt) nshgt * (PetscInt) nflow;
        petsc_PA  = (PetscInt) 40;
        ierr = MatCreateBAIJ(PETSC_COMM_WORLD, petsc_bs, petsc_m, petsc_m, petsc_M, petsc_M,
                             0, idiagnz, 0, iodiagnz, &lhsP);

        ierr = MatSetOption(lhsP, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

// Next is Jed Brown's improvement to imprint Assembly to make that stage scalable after the first call
#ifdef JEDBROWN        
        ierr = MatSetOption(lhsP, MAT_SUBSET_OFF_PROC_ENTRIES, PETSC_TRUE);
#endif        
        ierr = MatSetUp(lhsP);
      
      PetscInt myMatStart, myMatEnd;
      ierr = MatGetOwnershipRange(lhsP, &myMatStart, &myMatEnd);
//debug
      if(workfc.myrank == rankdump) printf("Flow myrank,myMatStart,myMatEnd %d,%d,%d, \n", workfc.myrank,myMatStart,myMatEnd);
      }
      if(genpar.lhs ==1) ierr = MatZeroEntries(lhsP);

//     
// .... *******************>> Element Data Formation <<******************
// 
// 
// .... set the parameters for flux and surface tension calculations
// 
// 
      genpar.idflx = 0 ;
      if(genpar.idiff >= 1)  genpar.idflx= genpar.idflx + (nflow-1) * nsd;
      if (genpar.isurf == 1) genpar.idflx=genpar.idflx + nsd;

      get_time(duration, (duration+1)); 
         
      ElmGMRPETSc(y,             ac,            x,
                  shp,           shgl,          iBC,
                  BC,            shpb,
                  shglb,         res,
                  rmes,                   iper,
                  ilwork,        rerr ,         &lhsP);
      
      get_time((duration+2), (duration+3));
      get_max_time_diff((duration), (duration+2), 
                        (duration+1), (duration+3),
                        "ElmGMRPETSc \0");  // char(0))
       if(firstpetsccall == 1) {
      // Setup IndexSet. For now, we mimic vector insertion procedure
      // Since we always reference by global indexes this doesn't matter
      // except for cache performance)
      // TODO: Better arrangment?
      PetscInt* indexsetary = malloc(sizeof(PetscInt)*nflow*nshg);
      PetscInt nodetoinsert;
        nodetoinsert = 0;
        k=0;
//debug
         if(workfc.myrank == rankdump) {
             printf("myrank,i,iownnodes,nshg %d %d %d %d \n",workfc.myrank,i,iownnodes,nshg);
             printf("myrank,mbeg,mend %d %d %d \n",workfc.myrank,mbeg,mend);
           }
        if(workfc.numpe > 1) {
          for (i=0; i<nshg ; i++) {
            nodetoinsert = fncorp[i]-1;
//debug
         if(workfc.myrank == rankdump) {
             printf("myrank,i,nodetoinsert %d %d %d \n",workfc.myrank,i,nodetoinsert);
         }
            
//            assert(fncorp[i]>=0);
            for (j=1; j<=nflow; j++) {
              indexsetary[k] = nodetoinsert*nflow+(j-1);
              assert(indexsetary[k]>=0);
//              assert(fncorp[i]>=0);
              k = k+1;
            }
          }
        }
        else {
          for (i=0; i<nshg ; i++) {
            nodetoinsert = i;
            for (j=1; j<=nflow; j++) {
              indexsetary[k] = nodetoinsert*nflow+(j-1);
              k = k+1;
            }
          }
        } 
        
//  Create Vector Index Maps
        petsc_n  = (PetscInt) nshg * (PetscInt) nflow;
        ierr = ISCreateGeneral(PETSC_COMM_SELF, petsc_n, indexsetary,
     	       PETSC_COPY_VALUES, &LocalIndexSet);
      free(indexsetary);
      }
      if(genpar.lhs == 1) {
      get_time((duration), (duration+1));
      ierr = MatAssemblyBegin(lhsP, MAT_FINAL_ASSEMBLY);
      ierr = MatAssemblyEnd(lhsP, MAT_FINAL_ASSEMBLY);
      get_time((duration+2),(duration+3));
      get_max_time_diff((duration), (duration+2), 
                        (duration+1), (duration+3),
                        "MatAssembly \0"); // char(0))
      get_time(duration, (duration+1));
      } 
      if(firstpetsccall == 1) {
        ierr = MatGetLocalSize(lhsP, &LocalRow, &LocalCol);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, LocalRow, petsc_M, &resP);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, LocalRow, petsc_M, &DyP);
      }
      ierr = VecZeroEntries(resP);
      if(firstpetsccall == 1) {
        ierr = VecCreateSeq(PETSC_COMM_SELF, petsc_n, &DyPLocal);
      }

      PetscRow=0;
      k = 0;
      int index;
      for (i=0; i<nshg; i++ ){
        for (j = 1; j<=nflow; j++){
          index = i + (j-1)*nshg;
          valtoinsert = res[index];
          if(workfc.numpe > 1) {
            PetscRow = (fncorp[i]-1)*nflow+(j-1);
          }
          else {
            PetscRow = i*nflow+(j-1);
          }
          assert(fncorp[i]<=nshgt);
          assert(fncorp[i]>0);
          assert(PetscRow>=0);
          assert(PetscRow<=nshgt*nflow);
          ierr =  VecSetValue(resP, PetscRow, valtoinsert, ADD_VALUES);
        }
      }
      ierr = VecAssemblyBegin(resP);
      ierr = VecAssemblyEnd(resP);
      ierr = VecNorm(resP,NORM_2,&resNrm);
      get_time((duration+2), (duration+3));
      get_max_time_diff((duration), (duration+2), 
                        (duration+1), (duration+3),
                        "VectorWorkPre \0"); // char(0))
     
      get_time((duration),(duration+1));
      
      if(firstpetsccall == 1) {
        ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
        ierr = KSPSetOperators(ksp, lhsP, lhsP);
        ierr = KSPGetPC(ksp, &pc);
        ierr = PCSetType(pc, PCPBJACOBI);
        PetscInt maxits;
        maxits = (PetscInt)  solpar.nGMRES * (PetscInt) solpar.Kspace;
        ierr = KSPSetTolerances(ksp, timdat.etol, PETSC_DEFAULT, PETSC_DEFAULT, maxits);
        ierr = KSPSetFromOptions(ksp); 
      }
      ierr = KSPSolve(ksp, resP, DyP);
      get_time((duration+2),(duration+3));
      get_max_time_diff((duration), (duration+2), 
                        (duration+1), (duration+3),
                        "KSPSolve \0"); // char(0))
      get_time((duration),(duration+1));
      if(firstpetsccall == 1) {
      ierr = VecScatterCreate(DyP, LocalIndexSet, DyPLocal, NULL, &scatter7);
      }
      ierr = VecScatterBegin(scatter7, DyP, DyPLocal, INSERT_VALUES, SCATTER_FORWARD);
      ierr = VecScatterEnd(scatter7, DyP, DyPLocal, INSERT_VALUES, SCATTER_FORWARD);
      ierr = VecGetArray(DyPLocal, &parray);
      PetscRow = 0;
      for ( node = 0; node< nshg; node++) {
        for (var = 1; var<= nflow; var++) {
          index = node + (var-1)*nshg;
          Dy[index] = parray[PetscRow];
          PetscRow = PetscRow+1;
        }
      }
      ierr = VecRestoreArray(DyPLocal, &parray);
      
      firstpetsccall = 0;

// .... output the statistics
// 
      itrpar.iKs=0; // see rstat()
      PetscInt its;
      ierr = KSPGetIterationNumber(ksp, &its);
      itrpar.iKs = (int) its;
      get_time((duration+2),(duration+3));
      get_max_time_diff((duration), (duration+2), 
                        (duration+1), (duration+3),
                        "solWork \0"); // char(0))
      itrpar.ntotGM += (int) its;
      rstatp (&resNrm);
//    
// .... end
//     
}

void     SolGMRpSclr(double* y,         double* ac, 
     	double* x,      double* elDw,   int* iBC,       double* BC,  
     	int* col,       int* row,      
     	int* iper,
     	int* ilwork,    double* shp,       double* shgl,      double* shpb,
     	double* shglb,  double* res,   double* Dy,     long long int* fncorp)
{
// 
// ----------------------------------------------------------------------
// 
//  This is the preconditioned GMRES driver routine.
// 
// input:
//  y      (nshg,ndof)           : Y-variables at n+alpha_v
//  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
//  yold   (nshg,ndof)           : Y-variables at beginning of step
//  acold  (nshg,ndof)           : Primvar. accel. variable at begng step
//  x      (numnp,nsd)            : node coordinates
//  iBC    (nshg)                : BC codes
//  BC     (nshg,ndofBC)         : BC constraint parameters
//  shp(b) (nen,maxsh,melCat)     : element shape functions (boundary)
//  shgl(b)(nsd,nen,maxsh,melCat) : local gradients of shape functions
// 
// ----------------------------------------------------------------------
// 


// Get variables from common_c.h
      int nshg, nflow, nsd, iownnodes;
      nshg  = conpar.nshg; 
      nsd = NSD; 
      iownnodes = conpar.iownnodes;

      gcorp_t nshgt;
      gcorp_t mbeg;
      gcorp_t mend;
      nshgt = (gcorp_t) newdim.nshgt; //Fix nshgt in common_c.h
      mbeg = (gcorp_t) newdim.minowned; 
      mend = (gcorp_t) newdim.maxowned; 


      int node, element, var, eqn;
      double valtoinsert;
      int nenl, iel, lelCat, lcsyst, iorder, nshl;
      int mattyp, ndofl, nsymdl, npro, ngauss, nppro;
      double DyFlats[nshg];
      double DyFlatPhastas[nshg];
      double rmes[nshg];
// DEBUG
      int i,j,k,l,m;

      double  real_rtol, real_abstol, real_dtol;
      double *parray; 
      PetscInt petsc_bs, petsc_m, petsc_M,petsc_PA;
      PetscInt petsc_n;
      PetscOne = 1;
      uint64_t duration[4];


//      
//     
// .... *******************>> Element Data Formation <<******************
// 
// 
// .... set the parameters for flux and surface tension calculations
// 
// 
      genpar.idflx = 0 ;
      if(genpar.idiff >= 1)  genpar.idflx= genpar.idflx + (nflow-1) * nsd;
      if (genpar.isurf == 1) genpar.idflx=genpar.idflx + nsd;



      gcorp_t glbNZ;

      if(firstpetsccalls == 1) {
        PetscInt* idiagnz= (PetscInt*) malloc(sizeof(PetscInt)*iownnodes);
        PetscInt* iodiagnz= (PetscInt*) malloc(sizeof(PetscInt)*iownnodes);
        for(i=0;i<iownnodes;i++) {
             idiagnz[i]=0;
             iodiagnz[i]=0;
        }
        i=0;
        for(k=0;k<nshg;k++) {
          if((fncorp[k] < mbeg) || (fncorp[k] >mend)){
// this node is not owned by this rank so we skip 
          } else { 
           for(j=col[i]-1;j<col[i+1]-1;j++) {
             glbNZ=fncorp[row[j]-1];
             if((glbNZ < mbeg) || (glbNZ > mend)) {
                iodiagnz[i]++;
             } else {
                idiagnz[i]++;
             }
           }
           i++; 
          }
        }
        gcorp_t mind=1000;
        gcorp_t mino=1000;
        gcorp_t maxd=0;
        gcorp_t maxo=0;
        for(i=0;i<iownnodes;i++) {
           mind=min(mind,idiagnz[i]);
           mino=min(mino,iodiagnz[i]);
           maxd=max(maxd,idiagnz[i]);
           maxo=max(maxo,iodiagnz[i]);
        }

        for(i=0;i<iownnodes;i++) {
           iodiagnz[i]=1.3*maxd;
           idiagnz[i]=1.3*maxd;
        }


        
//       }
//       if(firstpetsccalls == 1) {

        petsc_m  = (PetscInt) iownnodes;
        petsc_M  = (PetscInt) nshgt;
        petsc_PA  = (PetscInt) 40;
        
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, petsc_m, petsc_m, petsc_M, petsc_M,
                            0, idiagnz, 0, iodiagnz, &lhsPs);

        ierr = MatSetOption(lhsPs, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

// Next is Jed Brown's improvement to imprint Assembly to make that stage scalable after the first call
#ifdef HIDEJEDBROWN        
        ierr = MatSetOption(lhsPs, MAT_SUBSET_OFF_PROC_ENTRIES, PETSC_TRUE);
#endif        
        ierr = MatSetUp(lhsPs);
      
      PetscInt myMatStart, myMatEnd;
      ierr = MatGetOwnershipRange(lhsPs, &myMatStart, &myMatEnd);
      if(workfc.myrank ==rankdump) printf("Sclr myrank,myMatStart,myMatEnd %d,%d,%d, \n", workfc.myrank,myMatStart,myMatEnd);
      }
//      MPI_Barrier(MPI_COMM_WORLD); 
 //     if(workfc.myrank ==0) printf("Before MatZeroEntries  \n");
      if(genpar.lhs == 1) ierr = MatZeroEntries(lhsPs);

      get_time(duration, (duration+1)); 
         
//      MPI_Barrier(MPI_COMM_WORLD); 
 //     if(workfc.myrank ==0) printf("Before elmgmr  \n");
//      for (i=0; i<nshg ; i++) {
//            assert(fncorp[i]>0);
//      }

      ElmGMRPETScSclr(y,             ac,            x, elDw,
                  shp,           shgl,          iBC,
                  BC,            shpb,
                  shglb,         res,
                  rmes,                   iper,
                  ilwork,        &lhsPs);
       if(firstpetsccalls == 1) {
      // Setup IndexSet. For now, we mimic vector insertion procedure
      // Since we always reference by global indexes this doesn't matter
      // except for cache performance)
      // TODO: Better arrangment?

      PetscInt* indexsetarys = malloc(sizeof(PetscInt)*nshg);
      PetscInt nodetoinsert;

        nodetoinsert = 0;
        k=0;

//debug
         if(workfc.myrank == rankdump) {
             printf("myrank,i,iownnodes,nshg %d %d %d %d \n",workfc.myrank,i,iownnodes,nshg);
             printf("myrank,mbeg,mend %d %d %d \n",workfc.myrank,mbeg,mend);
           }

        if(workfc.numpe > 1) {
          for (i=0; i<nshg ; i++) {
            nodetoinsert = fncorp[i]-1;
//debug
         if(workfc.myrank == rankdump) {
             printf("myrank,i,nodetoinsert %d %d %d \n",workfc.myrank,i,nodetoinsert);
         }
            
            indexsetarys[k] = nodetoinsert;
            k = k+1;
          }
        }
        else {
          for (i=0; i<nshg ; i++) {
            nodetoinsert = i;
            indexsetarys[k] = nodetoinsert;
            k = k+1;
          }
        } 
        
//  Create Vector Index Maps
        petsc_n  = (PetscInt) nshg;
        ierr = ISCreateGeneral(PETSC_COMM_SELF, petsc_n, indexsetarys,
     	       PETSC_COPY_VALUES, &LocalIndexSets);
      free(indexsetarys);
      }
      if(genpar.lhs ==1) {
      ierr = MatAssemblyBegin(lhsPs, MAT_FINAL_ASSEMBLY);
      ierr = MatAssemblyEnd(lhsPs, MAT_FINAL_ASSEMBLY);
      }
      if(firstpetsccalls == 1) {
        ierr = MatGetLocalSize(lhsPs, &LocalRow, &LocalCol);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, LocalRow, petsc_M, &resPs);
        ierr = VecCreateMPI(PETSC_COMM_WORLD, LocalRow, petsc_M, &DyPs);
      }
      ierr = VecZeroEntries(resPs);
      if(firstpetsccalls == 1) {
        ierr = VecCreateSeq(PETSC_COMM_SELF, petsc_n, &DyPLocals);
      }

      PetscRow=0;
      k = 0;
      int index;
      for (i=0; i<nshg; i++ ){
          valtoinsert = res[i];
          if(workfc.numpe > 1) {
            PetscRow = (fncorp[i]-1);
          }
          else {
            PetscRow = i-1;
          }
          assert(fncorp[i]<=nshgt);
          assert(fncorp[i]>0);
          assert(PetscRow>=0);
          assert(PetscRow<=nshgt);
          ierr =  VecSetValue(resPs, PetscRow, valtoinsert, ADD_VALUES);
      }
      ierr = VecAssemblyBegin(resPs);
      ierr = VecAssemblyEnd(resPs);
      ierr = VecNorm(resPs,NORM_2,&resNrm);
      
      if(firstpetsccalls == 1) {
        ierr = KSPCreate(PETSC_COMM_WORLD, &ksps);
        ierr = KSPSetOperators(ksps, lhsPs, lhsPs);
        ierr = KSPGetPC(ksps, &pcs);
        ierr = PCSetType(pcs, PCPBJACOBI);
        PetscInt maxits;
        maxits = (PetscInt)  solpar.nGMRES * (PetscInt) solpar.Kspace;
        ierr = KSPSetTolerances(ksps, timdat.etol, PETSC_DEFAULT, PETSC_DEFAULT, maxits);
        ierr = KSPSetFromOptions(ksps); 
      }
      ierr = KSPSolve(ksps, resPs, DyPs);
      if(firstpetsccalls == 1) {
      ierr = VecScatterCreate(DyPs, LocalIndexSets, DyPLocals, NULL, &scatter7s);
      }
      ierr = VecScatterBegin(scatter7s, DyPs, DyPLocals, INSERT_VALUES, SCATTER_FORWARD);
      ierr = VecScatterEnd(scatter7s, DyPs, DyPLocals, INSERT_VALUES, SCATTER_FORWARD);
      ierr = VecGetArray(DyPLocals, &parray);
      PetscRow = 0;
      for ( node = 0; node< nshg; node++) {
          index = node;
          Dy[index] = parray[PetscRow];
          PetscRow = PetscRow+1;
      }
      ierr = VecRestoreArray(DyPLocals, &parray);
      
      firstpetsccalls = 0;

// .... output the statistics
// 
//      itrpar.iKss=0; // see rstat()
      PetscInt its;
      ierr = KSPGetIterationNumber(ksps, &its);
      itrpar.iKss = (int) its;
      itrpar.ntotGMs += (int) its;
      rstatpSclr (&resNrm);
//     
// .... end
//     
}
#endif
