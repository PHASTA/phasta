#include "shapeFuncInternals.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Calculate shape functions and derivative for tet elements */
int TetShapeAndDrv(int p,double par[3],double N[],double dN[][3]) {
  static int TetEMAP[6][2]={{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
  static int TetFMAP[4][3]={{0,1,2},{0,3,1},{1,3,2},{0,2,3}};
  int NS,is,i,j,ip,nshp=0,(*fpdrv)[2];
  double L[4],Le[2],Lf[3],mode,blend,bdrv[3],mdrv[3],epdrv[3][2],mfdrv[2];
  double tmp,rst,rsw,stw,rstw,rtw;
  if(p<1)
    return nshp;
  L[0]=par[0];
  L[1]=par[1];
  L[2]=par[2];
  L[3]=1.0e0-par[0]-par[1]-par[2];
  /* collect all vertex modes */
  for(i=0; i <4; i++) {
    N[i] = L[i];
    if(i==3) 
      dN[i][0]=dN[i][1]=dN[i][2]=-1.0e0;
    else {
      for(j=0; j <3; j++) {
	if(i==j)
	  dN[i][j] = 1.0e0;
	else
	  dN[i][j] = 0.0e0;
      }
    }
  }
  nshp=4;
  if( p > 1 ) {
    /* collect all edge modes, for all p */
    for(i=0; i <6; i++) {
      Le[0]=L[TetEMAP[i][0]];
      Le[1]=L[TetEMAP[i][1]];
      blend = R_edgeBlendTet(TetEMAP[i],L);
      R_edgeBlendTetDrv(TetEMAP[i],L,bdrv);
      E_parDrv(TetEMAP[i][0],TetEMAP[i][1],Stet,epdrv);
      for(ip=2; ip <= p; ip++) {
	mode = E_modeShape(ip,Le);
	E_modeShapeDrv(ip, Le, mdrv);
	N[nshp] = blend*mode;
	dN[nshp][0] = bdrv[0]*mode + 
	        blend*(mdrv[0]*epdrv[0][0]+mdrv[1]*epdrv[0][1]);
        dN[nshp][1] = bdrv[1]*mode + 
	        blend*(mdrv[0]*epdrv[1][0]+mdrv[1]*epdrv[1][1]);
        dN[nshp][2] = bdrv[2]*mode + 
	        blend*(mdrv[0]*epdrv[2][0]+mdrv[1]*epdrv[2][1]);
	nshp++;
      }
    }
  }
  if( p > 2 ) {
    /* collect all face modes for all p */
    for(i=0; i <4; i++) {
      Lf[0]=L[TetFMAP[i][0]];
      Lf[1]=L[TetFMAP[i][1]];
      Lf[2]=L[TetFMAP[i][2]];
      blend=Lf[0]*Lf[1]*Lf[2];
      R_faceBlendTetDrv(TetFMAP[i],L,bdrv);
      F_parDrv(TetFMAP[i][0],TetFMAP[i][1],TetFMAP[i][2],Stet,&fpdrv);
      for(ip=3; ip <= p; ip++) {
	NS = ip-2;
	for(is=0; is < NS; is++) {
	  mode = F_modeShapeTri(ip,is,Lf);
	  N[nshp] = blend*mode;
	  F_modeShapeTriDrv(ip,is,Lf,mfdrv);
	  mdrv[0]=mfdrv[0]*(double)fpdrv[0][0]+ mfdrv[1]*(double)fpdrv[0][1];
	  mdrv[1]=mfdrv[0]*(double)fpdrv[1][0]+ mfdrv[1]*(double)fpdrv[1][1];
	  mdrv[2]=mfdrv[0]*(double)fpdrv[2][0]+ mfdrv[1]*(double)fpdrv[2][1];
	  dN[nshp][0] = bdrv[0]*mode + blend*mdrv[0] ;      
	  dN[nshp][1] = bdrv[1]*mode + blend*mdrv[1] ;      
	  dN[nshp][2] = bdrv[2]*mode + blend*mdrv[2] ;
	  nshp++;
	}
      }
    }
  }
  if( p > 3 ) {
    /* collect all region modes */
    for(ip=4; ip <= p; ip++) {
      NS = (ip-3)*(ip-2)/2 ;
      blend = L[0]*L[1]*L[2]*L[3];
      tmp = 1.0e0-L[0]-L[1]-L[2] ;
      rsw = L[0]*L[1]*tmp;
      stw = L[1]*L[2]*tmp ;
      rtw = L[0]*L[2]*tmp ;
      rst = L[0]*L[1]*L[2] ;
      rstw = rst*L[3] ;
      for(is=0; is < NS; is++) {
	mode = R_modeShapeTet(ip,is,L);
	N[nshp] = blend*mode;
	R_modeShapeTetDrv(ip, is, L, mdrv);
	dN[nshp][0] = mode*(stw-rst) + rstw*mdrv[0] ; 
	dN[nshp][1] = mode*(rtw-rst) + rstw*mdrv[1] ;
	dN[nshp][2] = mode*(rsw-rst) + rstw*mdrv[2] ;
	nshp++;
      }
    }
  }
  return nshp;
}

/* calculate the shape functions and their derivatives for
   triangular faces
   */

int TriShapeAndDrv(int p,double par[2],double N[],double dN[][2]){
  int i,j,nshp=0;
  double L[3];
  
  if(p > 2) /* not supported */
    return nshp;
  L[0]=par[0];
  L[1]=par[1];
  L[2]=1.0e0-par[0]-par[1];
    
  /* define shape functions for a quadratic triangle */

  /* collect all vertex modes */
  for(i=0; i<3; i++) {
    N[i] = L[i];
    if(i==2)
      dN[i][0]=dN[i][1]=-1.0e0;
    else {
      for(j=0; j<2; j++) {
	if(i==j)
	  dN[i][j] = 1.0e0;
	else
	  dN[i][j] = 0.0e0;
      }
    }
  }
  nshp=3;
  if( p > 1 ){
    /* collect edge modes (only quadratic for now) */
    N[3] = -2.0e0*L[0]*L[1];
    N[4] = -2.0e0*L[1]*L[2];
    N[5] = -2.0e0*L[0]*L[2];
    dN[3][0] = -2.0e0*L[1];
    dN[3][1] = -2.0e0*L[0];
    dN[4][0] = 2.0e0*L[1];
    dN[4][1] = -2.0e0+2.0e0*L[0]+4.0e0*L[1];
    dN[5][0] = -2.0e0+4.0e0*L[0]+2.0e0*L[1];
    dN[5][1] = 2.0e0*L[0];
    nshp=6;
  }
  return nshp;
}
    
    
  
  


#ifdef __cplusplus
}
#endif

