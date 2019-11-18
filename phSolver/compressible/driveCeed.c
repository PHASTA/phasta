#include <ceed.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "assert.h"
#include "common_c.h"
#include "lccommon.h"
#include "advection.h"
#include <FCMangle.h>
#define driveceed  FortranCInterface_GLOBAL_(driveceed ,DRIVECEED)
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define get_time FortranCInterface_GLOBAL_(get_time,GET_TIME) 
#define get_max_time_diff FortranCInterface_GLOBAL_(get_max_time_diff,GET_MAX_TIME_DIFF) 

void get_time(uint64_t* rv, uint64_t* cycle);
void get_max_time_diff(uint64_t* first, uint64_t* last, uint64_t* c_first, uint64_t* c_last, char* lbl);

   
int    driveceed(double* y,   double* ac,  
     	double* x,         int* ien,   
     	double* rest)     
{
//
//  Can test by replacing passed data above by data written from PHASTA into the file data4libCEED.dat 
// written as follows (fortran)
//      open(unit=777, file='data4libCEED.dat',status='unknown')
//      write(777,*) npro,nen 
//      do i =1,npro
//        write(777,*) (ien(i,j), j=1,nen)
//      enddo
//      write(777,*) numnp,nsd
//      do i=1,numnp
//        write(777,*) (point2x(i,j),j=1,nsd)
//      enddo
//      write(777,*) numnp,ndof
//      do i=1,numnp
//         write(777,*) (y(i,j),j=1,ndof)
//      enddo
//      write(777,*) numnp,ndof
//      do i=1,numnp
//         write(777,*) (ac(i,j),j=1,ndof)
//      enddo
//      close(777)
//   note, rest is output data thus not needed
// numbers comming from common.h (PHASTA's not libCEEDs)   
// nshg=numnp
// nshl=nen (for this case)


// 
// ----------------------------------------------------------------------
// 
//  This is the libCEED driver routine.
// 
// input:
//  y      (nshg,ndof)           : Y-variables at n+alpha_v
//  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
//  acold  (nshg,ndof)           : Primvar. accel. variable at begng step
//  x      (numnp,nsd)           : node coordinates
//  ien    (npro,nshl)           : connectivity
// 
// output:
//  rest    (nshg)           : residual
// 
//  
// Jansen 2019
// ----------------------------------------------------------------------
// 


// Get variables from common_c.h
      int tmp32,nen,numnp,nshg, nflow, nsd, iownnodes;
      nshg  = conpar.nshg; 
      nflow = conpar.nflow; 
      numnp = conpar.numnp; 
      nen = conpar.nen; 
      nsd = NSD; 
      int node, element, var, eqn;
      double valtoinsert;
      int nenl, iel, lelCat, lcsyst, iorder;
      int mattyp, ndofl, nsymdl, npro, ngauss, nppro;
      npro = propar.npro; 
      double qfp[5*numnp]; // libCeed will get this array with scalar thrown into the temperature slot since advection works that way later we will need to make a real scalar equation. 
      double qdotfp[5*numnp];
      double q0[5*numnp]; // just for debugging using ics
      double x0[3*numnp]; // just for debugging using ics
// DEBUG
      int i,j,k,l,m;

      // FIXME: PetscScalar
      double  real_rtol, real_abstol, real_dtol;
// /DEBUG
//

  Ceed ceed;
  CeedElemRestriction restrictx, restrictq, restrictxi, restrictqdi,restrictxFake, restrictxcoord;
  CeedBasis basisxc, bx, bq; 
  CeedQFunction qf_setup, qf_ifunction, qf_ics;
  CeedOperator op_setup, op_ifunction,  op_ics;
  CeedVector qdata, X, U, Udot, V, Xfake;
  CeedVector xceed, q0ceed, qceed, qdotceed, gceed;
  const CeedScalar *hv;
  const CeedScalar *hvt;
  CeedInt nelem = npro, P = 2, Q = 2, qpownsd=Q*Q*Q;
  CeedInt nshl=P*P*P, qdatasize=10;
  CeedInt indx[npro*nen], indq[npro*nshl];
  CeedScalar xref[24];
  CeedScalar theta0     = 300.;     // K
  CeedScalar thetaC     = -15.;     // K
  CeedScalar P0         = 1.e5;     // Pa
  CeedScalar N          = 0.01;     // 1/s
  CeedScalar cv         = 717.;     // J/(kg K)
  CeedScalar cp         = 1004.;    // J/(kg K)
  CeedScalar g          = 9.81;     // m/s^2
  CeedScalar lx        = 8000.;    // m
  CeedScalar ly        = 8000.;    // m
  CeedScalar lz        = 4000.;    // m
  CeedScalar rc         = 1000.;    // m (Radius of bubble)
  CeedScalar Rd;
  CeedInt periodicity[3];
  Rd=288.29438; //PHASTA VALUE cp-cv;

  CeedScalar ctxSetup[] = {theta0, thetaC, P0, N, cv, cp, Rd, g, rc,
                           lx, ly, lz,
                           periodicity[0], periodicity[1], periodicity[2],
                          };


//  ! [Ceed Init]
//  const char* intStr="/cpu/self/ref/memcheck";
  const char* intStr="/cpu/self/ref/serial";
  CeedInit(intStr, &ceed);
//! [Ceed Init]
  for (CeedInt i=0; i<nelem; i++) {
    for (CeedInt j=0; j<nen; j++) 
      indx[j+i*nen] = ien[i+j*nelem]-1; // transpose and shift to C-based node numbering
  }
  for (CeedInt i=0; i<nelem; i++) {
      tmp32=indx[3+i*nen];
      indx[3+i*nen]=indx[2+i*nen];
      indx[2+i*nen]=tmp32;
      tmp32=indx[7+i*nen];
      indx[7+i*nen]=indx[6+i*nen];
      indx[6+i*nen]=tmp32;
  }
// jiggle the coordinates with a shift scaling off their node number to "encode" node number 
  for (CeedInt j=0; j<0*numnp; j++) {
      x[j]+=0.000001*j;
      x[j+numnp]+=0.000001*j;
      x[j+numnp*2]+=0.000001*j;
   }
  CeedVectorCreate(ceed, numnp*3, &X); //ns757
  CeedVectorSetArray(X, CEED_MEM_HOST, CEED_USE_POINTER, x); // the jiggled  coordinates will go INTO IC and solution

//! [Basis Create]
  CeedBasisCreateTensorH1Lagrange(ceed, 3, 3, 2, Q, CEED_GAUSS, &bx);
  CeedBasisCreateTensorH1Lagrange(ceed, 3, 5, P, Q, CEED_GAUSS, &bq);

//! [ElemRestr Create]
  CeedElemRestrictionCreate(ceed, npro, nen, numnp, 3, CEED_MEM_HOST, CEED_USE_POINTER, indx, &restrictx); // coordinates
  CeedElemRestrictionCreate(ceed, npro, nshl, nshg, 5, CEED_MEM_HOST, CEED_USE_POINTER, indx, &restrictq); // solution: change to indq if HO
  CeedElemRestrictionCreateIdentity(ceed, npro, qpownsd, qpownsd*npro, qdatasize, &restrictqdi); //metrics shared from setup to residual
  CeedElemRestrictionCreateIdentity(ceed, npro, qpownsd, qpownsd*npro, 1, &restrictxi); // weight
//! [QFunction Create]
  CeedQFunctionCreateInterior(ceed, 1, Setup, Setup_loc, &qf_setup);
  CeedQFunctionAddInput(qf_setup, "dx", nsd*nsd, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(qf_setup, "weight", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddOutput(qf_setup, "qdata", qdatasize, CEED_EVAL_NONE);
  // Create the operator that builds the quadrature data for the NS operator
  CeedOperatorCreate(ceed, qf_setup, NULL, NULL, &op_setup);
  CeedOperatorSetField(op_setup, "dx", restrictx, CEED_NOTRANSPOSE, //K AddInput says gradient using basisx of operatorApply defined from xcorners
                       bx, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(op_setup, "weight", restrictxi, CEED_NOTRANSPOSE, //K get the weight to setup Q function 
                       bx, CEED_VECTOR_NONE);
  CeedOperatorSetField(op_setup, "qdata", restrictqdi, CEED_NOTRANSPOSE,//K output of setup Q function metric data at quadrature points to share with op
                       CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);
  CeedVectorCreate(ceed, qdatasize*npro*qpownsd, &qdata);

  CeedOperatorApply(op_setup, X, qdata, CEED_REQUEST_IMMEDIATE); //K Apply setup for the selected case.  Creates qdata: WdetJ and dxidx at qpts


  CeedQFunctionCreateInterior(ceed, 1, IFunction_Advection, IFunction_Advection_loc, &qf_ifunction);
  CeedQFunctionAddInput(qf_ifunction, "q", 5, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(qf_ifunction, "dq", 5*nsd, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(qf_ifunction, "qdot", 5, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(qf_ifunction, "qdata", qdatasize, CEED_EVAL_NONE);
  CeedQFunctionAddOutput(qf_ifunction, "v", 5, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(qf_ifunction, "dv", 5*nsd, CEED_EVAL_GRAD);


  CeedVectorCreate(ceed, 5*nshg, &U);
  for (CeedInt i=0; i< nshg; i++) {
    qfp[i]=y[i+3*nshg]/Rd/y[i+4*nshg]; // density
    qfp[i+  nshg]=qfp[i]*y[i+0*nshg]; // density*u1
    qfp[i+2*nshg]=qfp[i]*y[i+1*nshg]; // density*u2
    qfp[i+3*nshg]=qfp[i]*y[i+2*nshg]; // density*u3
    qfp[i+4*nshg]=qfp[i]*y[i+5*nshg]; // PHASTA scalar
    qdotfp[i+4*nshg]= ac[i+5*nshg]; // PHASTA scalar
  }
  CeedVectorSetArray(U, CEED_MEM_HOST, CEED_USE_POINTER, qfp);

  CeedVectorCreate(ceed, 5*nshg, &Udot);
  CeedVectorSetArray(Udot, CEED_MEM_HOST, CEED_USE_POINTER, qdotfp);
  CeedVectorCreate(ceed, 5*nshg, &V);


  { // Create the IFunction operator  
    CeedOperatorCreate(ceed, qf_ifunction, NULL, NULL, &op_ifunction);
    CeedOperatorSetField(op_ifunction, "q", restrictq, CEED_NOTRANSPOSE, //K Active input is current solution vector Q set on OperatorApply line q=B_q_i G_q Q  
                         bq, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_ifunction, "dq", restrictq, CEED_NOTRANSPOSE, //K Active input is current solution vector Q set on OperatorApply line q=B_q_{gi} G_q Q  
                         bq, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_ifunction, "qdot", restrictq, CEED_NOTRANSPOSE, //K not an active vector but the is qdot (like ac in PHASTA)
                         bq, Udot);
    CeedOperatorSetField(op_ifunction, "qdata", restrictqdi, CEED_NOTRANSPOSE, //K shared data from setup is "set"  
                         CEED_BASIS_COLLOCATED, qdata);
    CeedOperatorSetField(op_ifunction, "v", restrictq, CEED_NOTRANSPOSE, //K Output 
                         bq, CEED_VECTOR_ACTIVE);
    CeedOperatorSetField(op_ifunction, "dv", restrictq, CEED_NOTRANSPOSE, //K Output
                         bq, CEED_VECTOR_ACTIVE);
  }
  double CtauS=1.0;
  int strong_form=0;
  int stab=2;
  struct Advection2dContext_ ctxAdvection2d = { //K struct that passes data needed at quadrature points for both advection
    .CtauS = CtauS,
    .strong_form = strong_form,
    .stabilization = stab,
  };
  CeedQFunctionSetContext(qf_ifunction, &ctxAdvection2d, sizeof ctxAdvection2d); //K This function associates the struct with qf_rhs and in next line qf_ifunction
// Calculate qdata 
  // Apply Setup Ceed Operators

  CeedOperatorApply(op_ifunction, U, V, CEED_REQUEST_IMMEDIATE); //K Apply setup for the selected case.  Creates qdata: WdetJ and dxidx at qpts

  CeedVectorGetArrayRead(V, CEED_MEM_HOST, &hv);
  const CeedInt shft=4*nshg;
  for (CeedInt i=0; i<nshg; i++)
    rest[i]=hv[i+shft];
  CeedVectorRestoreArrayRead(V, &hv);

// pausing here to line 897 of copying from nsplex.c

//  I know I  have not finished the cleanup below
  CeedQFunctionDestroy(&qf_setup);
  CeedQFunctionDestroy(&qf_ifunction);
  CeedOperatorDestroy(&op_setup);
  CeedOperatorDestroy(&op_ifunction);
  CeedElemRestrictionDestroy(&restrictq);
  CeedElemRestrictionDestroy(&restrictx);
  CeedElemRestrictionDestroy(&restrictxi);
  CeedElemRestrictionDestroy(&restrictqdi);
  CeedBasisDestroy(&bq);
  CeedBasisDestroy(&bx);
  CeedVectorDestroy(&U);
  CeedVectorDestroy(&Udot);
  CeedVectorDestroy(&V);
  CeedVectorDestroy(&qdata);
  CeedVectorDestroy(&X);
  CeedDestroy(&ceed);
  return 0;
// .... end
}
