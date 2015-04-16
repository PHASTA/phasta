/* calculate the shape functions and their derivatives for
   a triangular face
   */

int TriShapeAndDrv(int p,double par[2],double N[],double dN[][2]){
  int i,j,nshp=0;
  double L[3];
  
  if(p != 2)
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
    
    
  
  
