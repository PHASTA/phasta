#include <ceed.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "assert.h"
#include "common_c.h"
#include "lccommon.h"
#include "advection.h"
#include "densitycurrent_primitive.h"
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
     	double* res,      int lcmode)     
{

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
//  res    (lcmode*nshg)           : residual
// 
//  
// Jansen 2019
// ----------------------------------------------------------------------
// 


// Get variables from common_c.h
      int npro,tmp32,nen,numnp,nshg, nflow, nsd ;
      nshg  = conpar.nshg; 
      nflow = conpar.nflow; 
      numnp = conpar.numnp; 
      nen = conpar.nen; 
      nsd = NSD; 
      npro = propar.npro; 

  Ceed ceed;
  CeedElemRestriction restrictx, restrictq, restrictxi, restrictqdi;
  CeedBasis basisxc, bx, bq; 
  CeedQFunction qf_setup, qf_ifunction;
  CeedOperator op_setup, op_ifunction;
  CeedVector qdata, X, U, Udot, V;
  const CeedScalar *hv;
  CeedInt P = 2, Q = 2, qpownsd=Q*Q*Q;
  CeedInt nshl=P*P*P, qdatasize=10;
  CeedInt indx[npro*nen], indq[npro*nshl];
  CeedScalar theta0     = 300.;     // K
  CeedScalar thetaC     = -15.;     // K
  CeedScalar P0         = 1.e5;     // Pa
  CeedScalar N          = 0.01;     // 1/s
  CeedScalar Rd=mmatpar.Rgas; //  Rd=288.29438; //PHASTA VALUE cp-cv;
  CeedScalar gamma=mmatpar.gamma;
  CeedScalar cv         = Rd/(gamma-1.0);
  CeedScalar cp         = cv*gamma;    // J/(kg K)
  CeedScalar g          = 9.81;     // m/s^2
  CeedScalar lx        = 8000.;    // m
  CeedScalar ly        = 8000.;    // m
  CeedScalar lz        = 4000.;    // m
  CeedScalar rc         = 1000.;    // m (Radius of bubble)
  CeedInt periodicity[3];


  CeedScalar ctxSetup[] = {theta0, thetaC, P0, N, cv, cp, Rd, g, rc,
                           lx, ly, lz,
                           periodicity[0], periodicity[1], periodicity[2],
                          };


//  ! [Ceed Init]
//  const char* intStr="/cpu/self/ref/memcheck";
  const char* intStr="/cpu/self/ref/serial";
  CeedInit(intStr, &ceed);
//! [Ceed Init]
  for (CeedInt i=0; i<npro; i++) {
    for (CeedInt j=0; j<nen; j++) 
      indx[j+i*nen] = ien[i+j*npro]-1; // transpose and shift to C-based node numbering
  }
  for (CeedInt i=0; i<npro; i++) {
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


  if(lcmode==1) {
    CeedQFunctionCreateInterior(ceed, 1, IFunction_Advection, IFunction_Advection_loc, &qf_ifunction);
  } else {
    CeedQFunctionCreateInterior(ceed, 1, IFunction_DCPrim, IFunction_DCPrim_loc, &qf_ifunction);
  }
  CeedQFunctionAddInput(qf_ifunction, "q", 5, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(qf_ifunction, "dq", 5*nsd, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(qf_ifunction, "qdot", 5, CEED_EVAL_INTERP);
  CeedQFunctionAddInput(qf_ifunction, "qdata", qdatasize, CEED_EVAL_NONE);
  CeedQFunctionAddOutput(qf_ifunction, "v", 5, CEED_EVAL_INTERP);
  CeedQFunctionAddOutput(qf_ifunction, "dv", 5*nsd, CEED_EVAL_GRAD);


  CeedVectorCreate(ceed, 5*nshg, &U);
  CeedVectorCreate(ceed, 5*nshg, &Udot);
    double qfp[5*numnp]; // libCeed will get this array with scalar thrown into the temperature slot since advection works that way later we will need to make a real scalar equation. 
  if(lcmode==1) {
    double qdotfp[5*numnp];
    for (CeedInt i=0; i< nshg; i++) {
      qfp[i]=y[i+3*nshg]/Rd/y[i+4*nshg]; // density
      qfp[i+  nshg]=qfp[i]*y[i+0*nshg]; // density*u1
      qfp[i+2*nshg]=qfp[i]*y[i+1*nshg]; // density*u2
      qfp[i+3*nshg]=qfp[i]*y[i+2*nshg]; // density*u3
      qfp[i+4*nshg]=qfp[i]*y[i+5*nshg];  // shift in scalar if lcmode=1 
      qdotfp[i+4*nshg]= ac[i+5*nshg]; // PHASTA scalar
    }
    CeedVectorSetArray(Udot, CEED_MEM_HOST, CEED_USE_POINTER, qdotfp);
  } else {
    for (CeedInt i=0; i< nshg; i++) {
      qfp[i]=y[i+3*nshg]; // pressure
      qfp[i+1*nshg]=y[i+0*nshg]; // u1
      qfp[i+2*nshg]=y[i+1*nshg]; // u2
      qfp[i+3*nshg]=y[i+2*nshg]; // u3
      qfp[i+4*nshg]=y[i+4*nshg]; // T 
    }
    CeedVectorSetArray(Udot, CEED_MEM_HOST, CEED_USE_POINTER, ac);
  } 
  CeedVectorSetArray(U, CEED_MEM_HOST, CEED_USE_POINTER, qfp);
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
  double mu= matdat.datmat[0][1][0];
  double lambda= -2.0/3.0;
  double dt=1.0/timdat.Dtgl;
  double k=mu*cp/0.72;
  struct Advection2dContext_ ctxAdvection2d = { //K struct that passes data needed at quadrature points for both advection
     .CtauS = CtauS,
     .strong_form = strong_form,
     .stabilization = stab
  };
/*
  struct DCPrimContext_ ctxDCPrim = { 
     .lambda = lambda,
     .mu     = mu,
     .k      = k,
     .cv     = cv,
     .cp     = cp,
     .g      = g,
     .Rd     = Rd,
     .dt     = dt,
     .CtauS = CtauS,
     .strong_form = strong_form,
     .stabilization = stab
  };
*/
  CeedScalar ctxDCPrim[9]={lambda, mu, k, cv, cp, g, Rd, dt, stab};

  if(lcmode==1) {
    CeedQFunctionSetContext(qf_ifunction, &ctxAdvection2d, sizeof ctxAdvection2d); //K This function associates the struct with qf_rhs and in next line qf_ifunction
  } else {
    CeedQFunctionSetContext(qf_ifunction, &ctxDCPrim, sizeof ctxDCPrim); //K This function associates the struct with qf_rhs and in next line qf_ifunction
  }
// Calculate qdata 
  // Apply Setup Ceed Operators

  CeedOperatorApply(op_ifunction, U, V, CEED_REQUEST_IMMEDIATE); //K Apply setup for the selected case.  Creates qdata: WdetJ and dxidx at qpts

  CeedVectorGetArrayRead(V, CEED_MEM_HOST, &hv);
  CeedInt shft=4*nshg;
  if (lcmode==5) shft=0; // we only needed to shift if we solved a scalar hiding in the 5th equation like advection.h does 
  for (CeedInt i=0; i<lcmode*nshg; i++)
    res[i]=hv[i+shft];
  CeedVectorRestoreArrayRead(V, &hv);

// pausing here to line 897 of copying from nsplex.c

//  I am not sure if I finished the cleanup below
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
