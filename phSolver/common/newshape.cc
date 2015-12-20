/* This file contains the routines which generate hierarchic shape
   functions for a any (right now hexahedral) element, using the
   concept of Blend functions and entity level shapefunctions.

   The shapefunction for a mesh entity Mj^dj is defined as

              N = psi(Mj^di,Me^de)*phi(Mj^dj)

   where psi(...) is the blending function and phi(..) is the entity
   level function which takes care of the polynomial order .
*/



#include <math.h>
  
#include "topo_shapedefs.h"

int nshp =0;
double LP(int j, double x)
{
    /* Bonnet's recursion formula for Legendre Polynomials */
    if( j==0) return 1;
    if( j==1) return x;

    double ply;
    ply = ((2*j-1)*x*LP(j-1,x)-(j-1)*LP(j-2,x))/j;
    return ply;

}
/***************************************************************************/
double LPdrv(int j, double x)
{
    if(j == 0) return 0;
    if(j == 1) return 1;
  
    double plydrv;
    plydrv = (2*j -1)*LP(j-1,x)+LPdrv(j-2,x);
    return plydrv;
}
  
/***************************************************************************/
double phi(int p, double x)
{
    double PHI;

    PHI = LP(p,x)-LP(p-2,x);
    PHI = PHI / (2*p-1);
    /*  PHI = PHI/sqrt(2*(2*p-1)); */
    return PHI;
}  

double phiDrv(int p, double x)
{
    double Phidrv;
    Phidrv = LP(p-1,x);
    /*  Phidrv = sqrt(0.5*(2*p-1))*LP(p-1,x); */
    return Phidrv;
}


/*******************************************************************/
/* Construction of Blending functions */

/* Line element edge blend and its derivatives */

double Line_eB(double xi1)
{
    /* edge parameterization I defined in Saikat's thesis */
    return 0.5*(xi1*xi1 - 1);
}

double dLEBdxi1(double xi1){  return xi1; }
double dLEBdxi2(double xi1){return 0;}
double dLEBdxi3(double xi1){return 0;}

/* Quadrilateral edge blend and its derivatives */


double Quad_eB(double xi1, double xi2, int sign)
{
    /* Quad element edge blend, for edge along xi1 */
    return 0.25*(xi1*xi1 - 1)*(1+sign*xi2);
}

double dQEBdxi1(double xi1, double xi2, int sign)
{ return 0.5*xi1*(1+sign*xi2); }

double dQEBdxi2(double xi1, double xi2, int sign)
{ return 0.25*sign*(xi1*xi1 - 1); }

double dQEBdxi3(double xi1, double xi2, int sign)
{ return 0; }


/* Quadrilateral face blend and its derivatives */

double Quad_fB(double xi1, double xi2)
{
    return 0.25*(xi1*xi1 - 1)*(xi2*xi2 - 1);
}

double dQFBdxi1(double xi1, double xi2)
{ return 0.5*xi1*(xi2*xi2 - 1); }

double dQFBdxi2(double xi1, double xi2)
{ return 0.5*xi2*(xi1*xi1 - 1); }

double dQFBdxi3(double xi1, double xi2)
{ return 0.0 ; }


/* The hexahedral element edge blend and its derivatives */

double Hex_eB(double xi[3], int sign2, int sign3)
{
    /* For edge along xi1 */
  
    //  return 0.125*(xi[0]*xi[0] - 1)*(1+sign2*xi[1])*(1+sign3*xi[2]);
    return 0.25*(1+sign2*xi[1])*(1+sign3*xi[2]);
}

double dHEBdxi1(double xi[3], int sign2, int sign3)
{
    /* For edge along xi1 */

    //  return 0.25*xi[0]*(1+sign2*xi[1])*(1+sign3*xi[2]);
    return 0;
}

double dHEBdxi2(double xi[3], int sign2, int sign3)
{

    /* For edge along xi1 */

    //  return 0.125*sign2*(xi[0]*xi[0] - 1)*(1+sign3*xi[2]);

    return 0.25*sign2*(1+sign3*xi[2]);
}

double dHEBdxi3(double xi[3], int sign2, int sign3)
{
    /* For edge along xi1 */ 

    return 0.25*sign3*(1+sign2*xi[1]);
}


/* The hexahedral face blend and its derivatives */

double Hex_fB(double xi[3], int  sign3)
{
    /* for face perpendicular to xi3 */
  
    return 0.5*(1+sign3*xi[2]);
    //  return 0.125*(xi[0]*xi[0] - 1)*(xi[1]*xi[1] - 1)*(1+sign3*xi[2]);
}

double dHFBdxi1(double xi[3], int sign3)
    //{ return 0.25*xi[0]*(xi[1]*xi[1] - 1)*(1+sign3*xi[2]);}
{ return 0;}

double dHFBdxi2(double xi[3], int sign3)
    //{ return 0.25*xi[1]*(xi[0]*xi[0] - 1)*(1+sign3*xi[2]);}
{ return 0; }

double dHFBdxi3(double xi[3], int sign3)
    //{ return 0.125*(xi[0]*xi[0] - 1)*(xi[1]*xi[1] - 1)*sign3; }
{ return 0.5*sign3; }


// The Construction of higher-order blending functions for edge and face of pyramid */

/* The pyramid element edge blend and its derivatives */
double Pyr_eB(double xi[3], int sign[3], int k, int m, int along)
{
    double psi;
    double tmp0 =1/(1-xi[2]); // xi3->xi[2]
  
    k--; m--; along--;
  
    if (along==2) { // if actually along=3
        // For edge along xi3, which bounds the triangular face
        psi =0.125*(xi[2]*xi[2]-1)*(1+sign[k]*(2*xi[k]*tmp0))*(1+sign[m]*(2*xi[m]*tmp0));
    } else {        // if actually along=k=1 or 2
        // For edge along xik, k=1, 2, m=2, 1 which bounds the quadrilater face
        double tmp1 =2*xi[k]*tmp0;
        psi =0.125*(tmp1*tmp1-1)*(1+sign[m]*(2*xi[m]*tmp0))*(1-xi[2]);
    }
    return psi;
}

double dPeBdxi(double xi[3], int sign[3], int k, int m, int along, int byWhich)
{
    double dPsidxi = 0;
    double tmp0 =1/(1-xi[2]); // xi3->xi[2]
  
    // byWhich was decreased by 1.  When it comes in as j it is 0 initially
    // then it becomes -1, wrong.  Do not decrease byWhich.
    k--; m--; along--; 
    if (along==2) {  // for edges along xi3
        if (byWhich==k)
            dPsidxi =-0.25*(1+xi[2])*sign[k]*(1+sign[m]*(2*xi[m]*tmp0));
        else if (byWhich==m)
            dPsidxi =-0.25*(1+xi[2])*sign[m]*(1+sign[k]*(2*xi[k]*tmp0));
        else if (byWhich==2)        // actually byWhich =2
            dPsidxi =0.25*(xi[2]*(1+sign[k]*2*xi[k]*tmp0)*(1+sign[m]*2*xi[m]*tmp0)-
                           (1+xi[2])*tmp0*(sign[k]*xi[k]*(1+sign[m]*2*xi[m]*tmp0)+
                                           sign[m]*xi[m]*(1+sign[k]*2*xi[k]*tmp0)));
    } else {         // for edges along xik with k=1 or 2
        double tmp1 =2*xi[k]*tmp0;
        if (byWhich==k)
            dPsidxi =0.5*tmp1*(1+sign[m]*(2*xi[m]*tmp0));
        else if (byWhich==m)
            dPsidxi =0.25*(tmp1*tmp1-1)*sign[m];
        else if (byWhich==2)       // actually byWhich =2
            dPsidxi =0.25*(tmp1*tmp1)*(1+sign[m]*(2*xi[m]*tmp0))-0.125*(tmp1*tmp1-1);
    }
    return dPsidxi;
}

/* The pyramid element face blend and its derivatives */
/* Needs to be corrected as we did for the edge blend */
double Pyr_fB (double xi[3], int sign[3], int k, int m, int faceType)
{
    double tmp0 =1/(1-xi[1]); // xi2->xi[1]
    double psi;

    if (faceType==4) {
        // for the quadrilater face
        double tmp1 =2*xi[0]*tmp0;
        double tmp2 =2*xi[2]*tmp0;
    
        psi =0.125*(1-tmp1*tmp1)*(1-tmp2*tmp2)*(1-xi[1]);
    } else {
        // for the triangular faces with k=1, 3 and m=3, 1
        k--; m--;
        double tmp1 =2*xi[m]*tmp0;
    
        psi =0.125*(1+sign[k]*2*xi[k]*tmp0)*(1-tmp1*tmp1)*(1-xi[1]*xi[1]);
    }
    return psi;
}

double dPfBdxi(double xi[3], int sign[3], int k, int m, int faceType, int byWhich)
{
    double dPsidxi = 0;
    double tmp0 =1/(1-xi[2]); // xi3->xi[2]
  
    if (faceType==4) {  
        // for the quadrilater face
        double tmp1 =2*xi[0]*tmp0;
        double tmp2 =2*xi[2]*tmp0;
    
        switch (byWhich) {
        case 1:
            dPsidxi =-0.5*tmp1*(1-tmp2*tmp2);
            break;
        case 2:
            dPsidxi = -0.125* (1+tmp1*tmp1+tmp2*tmp2-3*tmp1*tmp1*tmp2*tmp2);
            break;
        case 3:
            dPsidxi =-0.5*(1-tmp1*tmp1)*tmp2;
            break;
        }
    } else {         // for the triangular face with k=1,3 and m=3,1
        k--; m--; byWhich--;
        if (byWhich==k) {
            double tmp1 =2*xi[m]*tmp0;
            dPsidxi =0.25*sign[k]*(1-tmp1*tmp1)*(1+xi[1]);
        }
        else if (byWhich==m)
            dPsidxi =-(1+sign[k]*2*xi[k]*tmp0)*(xi[m]*tmp0)*(1+xi[1]);
        else if (byWhich==1) {      // actually byWhich =2
            double tmp1 =xi[k]*tmp0;
            double tmp2 =xi[m]*tmp0;
            dPsidxi =(1+xi[1])*(0.25*sign[k]*tmp1*(1-4*tmp2*tmp2)-
                                (1+2*sign[k]*tmp1)*tmp2*tmp2)-
                0.25*(1+2*sign[k]*tmp1)*(1-4*tmp2*tmp2)*xi[1];
        }
    }
    return dPsidxi;
}


/* Add further Blending functions and their derivatives before this line */
/*******************************************************************/

/* Entity level functions  phi(....) */


int mesh_edge(double xi1, int gOrd[3], int p,double* entfn,double** edrv)
{
    int nem = p-1;
//    double leb = Line_eB(xi1);
//    double dlebdxi1 = dLEBdxi1(xi1);

    if(nem > 0){

        //      entfn = new double [nem];
        //      edrv = new double* [nem];
    
        for(int i =0; i< nem; i++) {

            // edrv[i] = new double [3];

//        entfn[i] = phi(i+2, xi1)/leb;
//        edrv[i][gOrd[0]] = (leb*phiDrv(i+2,xi1)-phi(i+2,xi1)*dlebdxi1)/leb*leb;
//        edrv[i][gOrd[1]] = 0.0;
//        edrv[i][gOrd[2]] = 0.0;

            entfn[i] = phi(i+2, xi1);
            edrv[i][gOrd[0]] = phiDrv(i+2,xi1);
            edrv[i][gOrd[1]] = 0.0;      
            edrv[i][gOrd[2]] = 0.0;

        }
    }
  
    return nem;
}

int quad_face(double xi[3], int gOrd[3], int p, double* entfn, double** edrv)
{
    int nfm;
    int a1, a2;
    int mc=0; /* mode counter*/
//    double temp1, temp2;
  
    double xi1 = xi[0];
    double xi2 = xi[1];

//    double qfb = Quad_fB(xi1,xi2);
//    double dqfbdxi1 = dQFBdxi1(xi1,xi2);
//    double dqfbdxi2 = dQFBdxi2(xi1,xi2);
  
    if(p > 3){

        nfm = (p-2)*(p-3)/2;
        //      entfn = new double [nfm];
        //      edrv = new double* [nfm];
        //      for(int i=0; i< nfm; i++){
        //        edrv[i]=new double [3];
        //      }
    
        for(int ip =3; ip <p+1; ip++){     /* for each p */

            for(a1 = 2; a1 < ip-1; a1++){    /* a1,a2 = 2....p-2 */
                for(a2 = 2; a2 < ip-1; a2++){  /* a1+a2 = p        */
                    if( a1+a1 == ip ){

//  	    entfn[mc] = phi(a1,xi1)*phi(a2,xi2)/qfb;

//  	    temp1 = (phi(a2,xi2)/qfb*qfb);
//  	    temp2 = (qfb*phiDrv(a1,xi1)- dqfbdxi1*phi(a1,xi1));
//  	    edrv[mc][gOrd[0]]= temp1*temp2;
	    
//  	    temp1 = (phi(a1,xi1)/qfb*qfb);
//  	    temp2 = (qfb*phiDrv(a2,xi2)- dqfbdxi2*phi(a2,xi2));
//  	    edrv[mc][gOrd[1]]= temp1*temp2;
	    
//  	    edrv[mc++][gOrd[2]]=0.0;

                        entfn[mc] = phi(a1,xi1)*phi(a2,xi2);
                        edrv[mc][gOrd[0]] = phiDrv(a1, xi1)*phi(a2, xi2);
                        edrv[mc][gOrd[1]] = phi(a1, xi1)*phiDrv(a2, xi2);
                        edrv[mc++][gOrd[2]]=0.0;
	    
                    }
                }
            }
        }
    } else nfm =0;
    return nfm;
} 

int hex_regn(double xi[3],int p, double* entfn, double** edrv)
{
    int a1, a2, a3;
    int nrm, mc=0;

    double xi1 = xi[0];
    double xi2 = xi[1];
    double xi3 = xi[2];

    if( p > 5){
    
        nrm = (p-3)*(p-4)*(p-5)/6;
        //      entfn = new double [nrm];
        //      edrv = new double* [nrm];
        //      for(int i=0; i< nrm; i++){
        //        edrv[i]=new double [3];
        //      }

        for(int ip =6; ip< p+1; ip++){
      
            for(a1=2; a1 < ip-3; a1++){
                for(a2=2; a2 < ip-3; a2++){
                    for(a3=2; a3 < ip-3; a3++){
                        if(a1+a2+a3 == ip){
                            entfn[mc]=phi(a1,xi1)*phi(a2,xi2)*phi(a3,xi3);  
	      
                            edrv[mc][0]=phiDrv(a1,xi1)*phi(a2,xi2)*phi(a3,xi3);

                            edrv[mc][1]=phi(a1,xi1)*phiDrv(a2,xi2)*phi(a3,xi3);
	      
                            edrv[mc++][2]=phi(a1,xi1)*phi(a2,xi2)*phiDrv(a3,xi3);
                        }
                    }
                }
            }
        }
    }else nrm =0;
    return nrm;
}



/* hex hierarchic shape function */
int HexShapeAndDrv(int p,double par[3],double N[],double dN[][3])
{
    int nshp = 0;
    int tmp1[4];
    int a,b;
    double EdgeBlend,dEBdxi,dEBdeta,dEBdzeta;
    double arg[3];
    int arg2[3];
    double* entfn;
    double** endrv;
    int num_e_modes, num_f_modes, num_r_modes;

    int** edge[12];
    int n[8][3]={{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
                 {-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};

    int face[6][4] = {{0,3,2,1},{0,1,5,4},{1,2,6,5},
                      {0,4,7,3},{2,3,7,6},{4,5,6,7}};
    int vrt[4], z;
    int normal = 0,sign;
    double FaceBlend, dFBdxi, dFBdeta, dFBdzeta;
  
    if(p<1) return nshp;
  
    double xi = par[0];
    double eta = par[1];
    double zeta = par[2];
  
    double xim=1-xi;
    double etam=1-eta;
    double zetam=1-zeta;
  
    double xip=1+xi;
    double etap=1+eta;
    double zetap=1+zeta;
  
    /* Shape functions for the Nodes. 
     *  There are eight nodal shapefunctions. These are same as the
     *  standard shape functions used in the eight-noded hexahedral
     *  elements   
     */
  
    N[0]= 0.125* xim * etam * zetam ;
    N[1]= 0.125* xip * etam * zetam ;
    N[2]= 0.125* xip * etap * zetam ;
    N[3]= 0.125* xim * etap * zetam ;
    N[4]= 0.125* xim * etam * zetap ;
    N[5]= 0.125* xip * etam * zetap ;
    N[6]= 0.125* xip * etap * zetap ;
    N[7]= 0.125* xim * etap * zetap ;
  
    /* Derivative of the above Shape Functions */
  
    dN[0][0]=-0.125*etam*zetam;
    dN[0][1]=-0.125*xim*zetam;
    dN[0][2]=-0.125*xim*etam;
  
    dN[1][0]=0.125 * etam * zetam;
    dN[1][1]=-0.125 * xip * zetam;
    dN[1][2]=-0.125 * xip * etam;
  
    dN[2][0]= 0.125 * etap * zetam;
    dN[2][1] = 0.125 * xip * zetam;
    dN[2][2] = -0.125 * xip * etap;
  
    dN[3][0] = -0.125 * etap * zetam;
    dN[3][1] = 0.125 * xim * zetam;
    dN[3][2] =-0.125 * xim * etap;
  
    dN[4][0] = -0.125 * etam * zetap;
    dN[4][1] = -0.125 * xim * zetap;
    dN[4][2] = 0.125 * xim * etam;
  
  
    dN[5][0] = 0.125 * etam * zetap;
    dN[5][1] =-0.125 * xip * zetap;
    dN[5][2] = 0.125 * xip * etam;
  
    dN[6][0] = 0.125 * etap * zetap;
    dN[6][1] = 0.125 * xip * zetap;
    dN[6][2] = 0.125 * xip * etap;
  
    dN[7][0] = -0.125 * etap * zetap ;
    dN[7][1] = 0.125 * xim * zetap;
    dN[7][2] = 0.125 * xim * etap;
  
    nshp = 8;

    if( p > 1) {
    
        /* Generate Shape Functions for Edge Modes.
         * For a polynomial Order of p there will be 12*(p-1)
         * edge modes for the entire element.
         */ 
    
        /*
         * edge order description;
         */
    
        for(int y=0;y<12;y++)
        {
            edge[y]=new int * [2];
        }
    
        edge[0][0]=n[0];edge[0][1]=n[1];
        edge[1][0]=n[1];edge[1][1]=n[2];
        edge[2][0]=n[2];edge[2][1]=n[3];
        edge[3][0]=n[3];edge[3][1]=n[0];
        edge[4][0]=n[4];edge[4][1]=n[5];
        edge[5][0]=n[5];edge[5][1]=n[6];
        edge[6][0]=n[6];edge[6][1]=n[7];
        edge[7][0]=n[7];edge[7][1]=n[4];
        edge[8][0]=n[0];edge[8][1]=n[4];
        edge[9][0]=n[1];edge[9][1]=n[5];
        edge[10][0]=n[2];edge[10][1]=n[6];
        edge[11][0]=n[3];edge[11][1]=n[7];


        int nem = p-1;

        for(int e=0; e < 12; e++){

            for(z=0;z<3;z++){
                if(!(edge[e][0][z]==edge[e][1][z]))
                    tmp1[3]=z;
                tmp1[z]=edge[e][1][z];
            }
      
            if((arg2[0] = tmp1[3]) == 0) { arg2[1]=1;arg2[2]=2 ;}
            else if(tmp1[3]==1) { arg2[1]=2;arg2[2]=0;}
            else { arg2[1]=0;arg2[2]=1;}
      
            /* arg[0]=par[arg2[0]]*tmp1[arg2[0]]; */
            arg[0]=par[arg2[0]];
            arg[1]=par[arg2[1]];  
            arg[2]=par[arg2[2]];

            a = arg2[1];
            b = arg2[2];

            EdgeBlend = Hex_eB(arg,tmp1[a],tmp1[b]);

            if( tmp1[3] == 0){
                dEBdxi    = dHEBdxi1(arg,tmp1[a],tmp1[b]); 
                dEBdeta   = dHEBdxi2(arg,tmp1[a],tmp1[b]);
                dEBdzeta  = dHEBdxi3(arg,tmp1[a],tmp1[b]);
            } else if( tmp1[3] == 1 ){
                dEBdeta   = dHEBdxi1(arg,tmp1[a],tmp1[b]);
                dEBdzeta  = dHEBdxi2(arg,tmp1[a],tmp1[b]);
                dEBdxi    = dHEBdxi3(arg,tmp1[a],tmp1[b]);
            } else {
                dEBdzeta  = dHEBdxi1(arg,tmp1[a],tmp1[b]);
                dEBdxi    = dHEBdxi2(arg,tmp1[a],tmp1[b]);
                dEBdeta   = dHEBdxi3(arg,tmp1[a],tmp1[b]);
            }

            entfn = new double [nem];
            endrv = new double* [nem];
            for(z =0; z < nem ; z++){
                endrv[z]=new double [3];
            }

            num_e_modes = mesh_edge(arg[0],arg2,p,entfn,endrv);

            for(int nm=0; nm < num_e_modes; nm++){
	
                N[nshp]= EdgeBlend * entfn[nm];
                dN[nshp][0]= EdgeBlend * endrv[nm][0] + dEBdxi * entfn[nm];
                dN[nshp][1]= EdgeBlend * endrv[nm][1] + dEBdeta * entfn[nm];
                dN[nshp][2]= EdgeBlend * endrv[nm][2] + dEBdzeta * entfn[nm];

                nshp++;
            }

            delete [] entfn;
            for(z=0; z< num_e_modes; z++){
                delete [] endrv[z]; }
            delete [] endrv;

        }
    }
  
    if( p > 3 ) {

        /* Generating Shape functions for the face modes .
         * For a Hexahedral Element there are 3*(p-2)*(p-3) shape
         * functions associated with the face modes [ p>=4]
         */
    
        int nfm = (p-2)*(p-3)/2;

        for(int f=0; f< 6; f++) {
            /* Identitfying the normal to the face */
      
            for(int u=0;u<4;u++){
                vrt[u]=face[f][u];
            }
            for(int v=0;v<3;v++)
            {
                if((n[vrt[0]][v]==n[vrt[1]][v])&&(n[vrt[1]][v]==n[vrt[2]][v]))
                    normal=v;
                sign=n[vrt[0]][v];
            }

      
            arg2[2] = normal;

            if( normal==0) { arg2[0]=1;arg2[1]=2 ;}
            else if(normal==1) { arg2[0]=2;arg2[1]=0;}
            else { arg2[0]=0;arg2[1]=1;}

            arg[0] = par[arg2[0]];
            arg[1] = par[arg2[1]];
            arg[2] = par[arg2[2]];
      
            FaceBlend = Hex_fB(arg,sign);

            if( normal == 0){
                dFBdxi    = dHFBdxi3(arg,sign); 
                dFBdeta   = dHFBdxi1(arg,sign);
                dFBdzeta  = dHFBdxi2(arg,sign);
            } else if( normal == 1 ){	 
                dFBdeta   = dHFBdxi3(arg,sign);
                dFBdzeta  = dHFBdxi1(arg,sign);
                dFBdxi    = dHFBdxi2(arg,sign);
            } else {			 
                dFBdzeta  = dHFBdxi3(arg,sign);
                dFBdxi    = dHFBdxi1(arg,sign);
                dFBdeta   = dHFBdxi2(arg,sign);
            }

            entfn = new double [nfm];
            endrv = new double* [nfm];
            for(z=0; z < nfm ; z++){
                endrv[z]=new double [3];
            }

            num_f_modes = quad_face(arg,arg2,p,entfn,endrv);

            for(int nm =0; nm < num_f_modes ; nm++){
	
                N[nshp] =  FaceBlend * entfn[nm]; 
                dN[nshp][0] = FaceBlend * endrv[nm][0] + dFBdxi * entfn[nm];
                dN[nshp][1] = FaceBlend * endrv[nm][1] + dFBdeta * entfn[nm];
                dN[nshp][2] = FaceBlend * endrv[nm][2] + dFBdzeta * entfn[nm];

                nshp++;
            }

            delete [] entfn;
            for( z=0; z< num_f_modes; z++){
                delete [] endrv[z]; }
            delete [] endrv;
        }
    }

    if ( p > 5 ) {

        int nrm = (p-3)*(p-4)*(p-5)/6;

        entfn = new double [nrm];
        endrv = new double* [nrm];
        for(z =0; z < nrm ; z++){
            endrv[z]=new double [3];
        }

        num_r_modes = hex_regn(par,p,entfn,endrv);

        for(int nm =0; nm < num_r_modes ; nm++){
            N[nshp] = entfn[nm];
            for(int d=0; d < 3; d++){
                dN[nshp][d] = endrv[nm][d];
            }
            nshp++;
        }
        delete [] entfn;
        for(z=0; z< num_r_modes; z++){
            delete [] endrv[z]; }
        delete [] endrv;
    
    }
    return nshp;
}

/* hierarchic wedge element shape function */

int WedgeShapeAndDrv( int p, double Inputpar[3], double N[], double dN[][3] )
{
    int i,j;
    double FaceBlend, FaceBlendDrv[4];
    double FaceEnt; //, FaceEntDrv[2][3];
    double par[4];
    //  int sign;
    //  int temp[4]={0,0,0,0};
    double EdgeBlend[9], EdgeBlendDrv[9][3];
    //  double arg[3];
    //  double arg2[2];
    //  int ** edge[9];
    double entfn[9];
    double endrv[9][3];
    //  int num_e_modes;
    // int n[6][4]={{1,0,0,-1},{0,1,0,-1},{0,0,1,-1},{1,0,0,1},{0,1,0,1},{0,0,1,1}};
  
    int nshp = 0;
    par[0] = 1.0 - Inputpar[0]- Inputpar[1];
    par[1] = Inputpar[0];
    par[2] = Inputpar[1];
    par[3] = Inputpar[2];


    if(p<1) return nshp;

   /* there are six nodal shape functions same as the standard
      shape functions used in the six-noded wedge element */
  
   N[0]=0.5*par[0]*(1.0-par[3]);
   N[1]=0.5*par[1]*(1.0-par[3]);
   N[2]=0.5*par[2]*(1.0-par[3]);
   N[3]=0.5*par[0]*(1.0+par[3]);
   N[4]=0.5*par[1]*(1.0+par[3]);
   N[5]=0.5*par[2]*(1.0+par[3]);

   /* Derivative of the above Shape functions */
   dN[0][0]=-0.25*(1.0-par[3]);
   dN[0][1]=-0.25*(1.0-par[3]);
   dN[0][2]=-0.5*par[0];

   dN[1][0]=0.25*(1.0-par[3]);
   dN[1][1]=0.0;
   dN[1][2]=-0.5*par[1];

   dN[2][0]=0.0;
   dN[2][1]=0.25*(1.0-par[3]);
   dN[2][2]=-0.5*par[2];

   dN[3][0]=-0.25*(1.0+par[3]);
   dN[3][1]=-0.25*(1.0+par[3]);
   dN[3][2]=0.5*par[0];

   dN[4][0]=0.25*(1.0+par[3]);
   dN[4][1]=0.0;
   dN[4][2]=0.5*par[1];
  
   dN[5][0]=0.0;
   dN[5][1]=0.25*(1.0+par[3]);
   dN[5][2]=0.5*par[2];

  
   nshp = 6;

//   if( p > 1 ){

//     /* calculate the blending shape functions and their drivitative */
//      EdgeBlend[0] = -par[0]*par[1]*(1.0-par[3]);
//      EdgeBlendDrv[0][0] = -0.5*(par[0]-par[1])*(1.0-par[3]);
//      EdgeBlendDrv[0][1] = 0.5*par[1]*(1.0-par[3]);
//      EdgeBlendDrv[0][2] = par[0]*par[1];

//      EdgeBlend[1] = -par[1]*par[2]*(1.0-par[3]);
//      EdgeBlendDrv[1][0] = -0.5*par[2]*(1.0-par[3]);
//      EdgeBlendDrv[1][1] = -0.5*par[1]*(1.0-par[3]);
//      EdgeBlendDrv[1][2] = par[1]*par[2];
     
//      EdgeBlend[2] = -par[0]*par[2]*(1.0-par[3]);
//      EdgeBlendDrv[2][0] = 0.5*par[2]*(1.0-par[3]);
//      EdgeBlendDrv[2][1] = -0.5*(par[0]-par[2])*(1.0-par[3]);
//      EdgeBlendDrv[2][2] = par[0]*par[2];

//      EdgeBlend[3] = -par[0]*par[1]*(1.0+par[3]);
//      EdgeBlendDrv[3][0] = -0.5*(par[0]-par[1])*(1.0+par[3]);
//      EdgeBlendDrv[3][1] = 0.5*par[1]*(1.0+par[3]);
//      EdgeBlendDrv[3][2] = -par[0]*par[1];

//      EdgeBlend[4] = -par[1]*par[2]*(1.0+par[3]);
//      EdgeBlendDrv[4][0] = -0.5*par[2]*(1.0+par[3]);
//      EdgeBlendDrv[4][1] = -0.5*par[1]*(1.0+par[3]);
//      EdgeBlendDrv[4][2] = -par[1]*par[2];
     
//      EdgeBlend[5] = -par[0]*par[2]*(1.0+par[3]);  
//      EdgeBlendDrv[5][0] = 0.5*par[2]*(1.0+par[3]);
//      EdgeBlendDrv[5][1] = -0.5*(par[0]-par[2])*(1.0+par[3]);
//      EdgeBlendDrv[5][2] = -par[0]*par[2];

//      EdgeBlend[6] = 0.5*par[0]*(par[3]*par[3]-1.0);
//      EdgeBlendDrv[6][0] = -0.25*(par[3]*par[3] - 1.0);
//      EdgeBlendDrv[6][1] = -0.25*(par[3]*par[3] - 1.0);
//      EdgeBlendDrv[6][2] = par[0]*par[3];

//      EdgeBlend[7] = 0.5*par[1]*(par[3]*par[3]-1.0);
//      EdgeBlendDrv[7][0] = 0.25*(par[3]*par[3] - 1.0);
//      EdgeBlendDrv[7][1] = 0.0;
//      EdgeBlendDrv[7][2] = par[1]*par[3];
     
//      EdgeBlend[8] = 0.5*par[2]*(par[3]*par[3]-1.0);
//      EdgeBlendDrv[8][0] = 0.0;
//      EdgeBlendDrv[8][1] = 0.25*(par[3]*par[3] - 1.0);
//      EdgeBlendDrv[8][2] = par[2]*par[3];
      
//       /* calculate the mesh entity shape function */
     
//      /* for the edge in the two triangle faces */
//      for(j=0; j<9; j++){
//        /* this is where I think the change in phi is for wedges */
//        entfn[j] = 1.0;
// 	 /*    entfn[j] = sqrt(1.5); */  
//        for(i=0; i<3; i++)
// 	 endrv[j][i] = 0.0;
//      }
     
// //       /* for the edge in the perpendicular to triangle faces */
// //       for(j=6; j<9; j++){
// //         entfn[j] = phi(2, par[3]);
// //         for(i=0; i<3; i++)
// //  	 endrv[j][i] = 0.0;	      
// //       }
     
//      for(j=0; j<9; j++){
//        N[nshp] = EdgeBlend[j]*entfn[j];
//        dN[nshp][0] = EdgeBlend[j]*endrv[j][0] + EdgeBlendDrv[j][0]*entfn[j];
//        dN[nshp][1] = EdgeBlend[j]*endrv[j][1] + EdgeBlendDrv[j][1]*entfn[j];
//        dN[nshp][2] = EdgeBlend[j]*endrv[j][2] + EdgeBlendDrv[j][2]*entfn[j];
	
//        ++nshp;
//      }
      
//   }
  
//   /* generate the tri face shape fuction */
//   if( p > 2 ){
//     for(i=0;i<11;i++){
//       /* get the face blending functions */
//       if(par[3]>0){
// 	FaceBlend = 0.5*par[0]*par[1]*par[2]*(1.0-par[3]);
	
// 	/* derivate the face blending functions */
// 	FaceBlendDrv[0]= 0.5*par[1]*par[2]*(1.0-par[3]);
// 	FaceBlendDrv[1]= 0.5*par[0]*par[2]*(1.0-par[3]);
// 	FaceBlendDrv[2]= 0.5*par[0]*par[1]*(1.0-par[3]);
// 	FaceBlendDrv[3]=-0.5*par[0]*par[1]*par[2];
//       } 
//       else{
// 	FaceBlend = 0.5*par[0]*par[1]*par[2]*(1.0+par[3]);
	
// 	/* derivate the face blending functions */
// 	FaceBlendDrv[0]= 0.5*par[1]*par[2]*(1.0+par[3]);
// 	FaceBlendDrv[1]= 0.5*par[0]*par[2]*(1.0+par[3]);
// 	FaceBlendDrv[2]= 0.5*par[0]*par[1]*(1.0+par[3]);
// 	FaceBlendDrv[3]= 0.5*par[0]*par[1]*par[2];
//       }
      
//       /*calculate the mesh entity function for cubic triangular face */
//       FaceEnt = 1.0;
      
//       /*calculate the shape functions*/
//       N[nshp] = FaceBlend*FaceEnt;
//       dN[nshp][0]=FaceBlendDrv[0];
//       dN[nshp][1]=FaceBlendDrv[1];
//       dN[nshp][2]=FaceBlendDrv[2];
//       dN[nshp][3]=FaceBlendDrv[3];
      
//       nshp++;
//     }
//   }

    return nshp;
}

/* Pyramid hierarchic shape function up to order 3*/
int PyrShapeAndDrv(int p,double par[3],double N[],double dN[][3])
{
    int i;
    double EdgeBlend, EdgeBlendDrv[3];
    double EdgeEntity[2], EdgeEntityDrv[2][3];
    double FaceBlend, FaceBlendDrv[3];
    double FaceEntity[2], FaceEntityDrv[2][3];

    //dEBdxi,dEBdeta,dEBdzeta;
    //  double arg[3];
    // int arg2[3];
    // double* entfn;
    //  double** endrv;
    // int num_e_modes, num_f_modes, num_r_modes;

    int** edge[8];  // total 8 edges
    int n[5][3]={{-1,-1, -1},{1,-1,-1},{1,1,-1},{-1,1, -1}, {0, 0, 1}}; // vertex coordinates

//    int n[5][3]={{-1,-1, 1},{-1,-1,-1},{1,-1,-1},{1,-1, 1}, {0, 1, 0}}; // vertex coordinates

    //  int face[5][4] = {{0,3,2,1},{0,1,4, -1},{1,2,4, -1},{2,3,4,-1},{3,0,4,-1}}; // face vertices
    //  int face[5][4] = {{0,1,2,3},{1,0,4, -1},{0,3,4, -1},{3,2,4,-1},{2,1,4,-1}}; // face vertices
    //  int vrt[4];
    //  int normal;
    int sign[3];
  
    if(p<1) return nshp;
  
 
    // when p=1, there are only nodal shapefunctions
    double zeta = par[2];
    double xi = par[0];
    double eta = par[1];

    //  double xim=1-xi;
    //  double etam=1-eta;
    double zetam=1-zeta;

    double zetamap=2.0/zetam;
  
    //  double xip=1+xi;
    // double etap=1+eta;
    double zetap=1+zeta;

    double xipp=1+xi*zetamap;
    double etapp=1+eta*zetamap;

    double ximp=1-xi*zetamap;
    double etamp=1-eta*zetamap;
  
    // Shape functions for the Nodes. 
    // There are five nodal shapefunctions.

    N[0]= 0.125* ximp * etamp * zetam ;
    N[1]= 0.125* xipp * etamp * zetam ;
    N[2]= 0.125* xipp * etapp * zetam ;
    N[3]= 0.125* ximp * etapp * zetam ;
    N[4]= 0.5* zetap; 
  
    // Derivative of the above Shape Functions
    dN[0][0] =-0.25 *  etamp;
    dN[0][1] =-0.25 *  ximp;
    dN[0][2] = 0.125 * (xi*eta*zetamap*zetamap-1);
                         
    dN[1][0] = 0.25 * etamp;
    dN[1][1] =-0.25 * xipp;
    dN[1][2] =-0.125 * (xi*eta*zetamap*zetamap+1);
                         
    dN[2][0] = 0.25 * etapp;
    dN[2][1] = 0.25 * xipp;
    dN[2][2] = dN[0][2];

    dN[3][0] =-0.25 * etapp;
    dN[3][1] = 0.25 * ximp;
    dN[3][2] = dN[1][2];
  
    dN[4][0] = 0;
    dN[4][1] = 0;
    dN[4][2] = 0.5;
  
    nshp = 5;

    if( p>1) {
    
        // Generate Shape Functions for Edges.
        // For a polynomial Order of p there are 8*(p-1)
        // edge modes for the entire element.
    
        // allocate spaces for edge order description
        for(i=0; i<8; i++)
            edge[i]=new int * [2];
    
        // for the edges on the quadrilateral face
        edge[0][0]=n[0]; edge[0][1]=n[1];
        edge[1][0]=n[1]; edge[1][1]=n[2];
        edge[2][0]=n[2]; edge[2][1]=n[3];
        edge[3][0]=n[3]; edge[3][1]=n[0];
    
        // for other edges on triangular faces
        edge[4][0]=n[0]; edge[4][1]=n[4];
        edge[5][0]=n[1]; edge[5][1]=n[4];
        edge[6][0]=n[2]; edge[6][1]=n[4];
        edge[7][0]=n[3]; edge[7][1]=n[4];
    
    
        // calculate the edge blending functions and their derivatives
        for(i=0; i < 8; i++) {
            int k, m, along, j;
            switch (i) {
                // first four are on quadrilateral face

                // P.N anil fix this 
            case 0:
                k =1;	m =2;   sign[m-1] =-1;  along =k;
                break;
            case 1:
                k =2;   m =1;   sign[m-1] =1;   along =k;
                break;
            case 2:
                k =1;   m =2;   sign[m-1] =1;   along =k;
                break;
            case 3:
                k =2;   m =1;   sign[m-1] =-1;  along =k;
                break;
                // next four are on triangular faces
            case 4:
                k =1;   m =2;   sign[k-1] =-1;    sign[m-1] =-1;  along =3;
                break;
            case 5:
                k =2;   m =1;   sign[k-1] =-1;    sign[m-1] =1; along =3;
                break;
            case 6:
                k =1;   m =2;   sign[k-1] =1;     sign[m-1] =1; along =3;
                break;
            case 7:
                k =2;   m =1;   sign[k-1] =1;     sign[m-1] =-1; along =3;
            }
            EdgeBlend =Pyr_eB (par, sign, k, m, along);
            for (j=0; j<3; j++)
                EdgeBlendDrv[j] =dPeBdxi(par, sign, k, m, along, j);
 
            // calculate the mesh entity shape function
    
            // calculate the edge entity function for p=2
            EdgeEntity[0] =1;
            EdgeEntityDrv[0][0] =EdgeEntityDrv[0][1] =EdgeEntityDrv[0][2] =0;
      
            if (p>2) {
                // calculate the edge entity function for p=3
                if (i<4) {
                    // for the edges of the quadrilateral face with parameterization I
                    // equation (33)
                    EdgeEntity[1] =par[k-1];
                    for (j=0; j<3; j++)
                        EdgeEntityDrv[1][j] =(int)(k-1==j);
                } else {
                    // for the edges of the triangular faces with parameterization II
                    // First of all, we need to represent these edges use xi's
                    // In our definition, u[0] points to the quadrilateral base
                    //                    u[1] points to the top
                    // then, we get	u[0] =(1-xi[1])/2;	  u[1] =(1+xi[1])/2;
                    //                    EdgeEntity[1] =u[1]-u[0] =xi[1];           
                    EdgeEntity[1] =par[1];
                    EdgeEntityDrv[1][1] =1; EdgeEntityDrv[1][0]=EdgeEntityDrv[1][2]=0;
                }
            }
      
            for(j=0; j < p-1; j++){
	
                N[nshp]= EdgeBlend * EdgeEntity[j];
                dN[nshp][0]= EdgeBlend * EdgeEntityDrv[j][0] + EdgeBlendDrv[0] * EdgeEntity[j];
                dN[nshp][1]= EdgeBlend * EdgeEntityDrv[j][1] + EdgeBlendDrv[1] * EdgeEntity[j];
                dN[nshp][2]= EdgeBlend * EdgeEntityDrv[j][2] + EdgeBlendDrv[2] * EdgeEntity[j];
                nshp++;
            }
        }

        // calculate the shape function for triangular faces 
        if (p==3) {
            for (i=1; i<5; i++) {
                // calculate the face blending functions
                int k, m, j;
                switch (i) {
                case 1:
                    k =1;   m =3;   sign[k-1] =-1;
                    break;
                case 2:
                    k =3;   m =1;   sign[k-1] =-1;
                    break;
                case 3:
                    k =1;   m =3;   sign[k-1] =1;
                    break;
                case 4:
                    k =3;   m =1;   sign[k-1] =1;
                    break;
                }
                FaceBlend =Pyr_fB (par, sign, k, m, 3);
                for (j=0; j<3; j++)
                    FaceBlendDrv[j] =dPfBdxi(par, sign, k, m, 3, j);

                // for p=3
                FaceEntity[0] =1;
                FaceEntityDrv[0][0] =FaceEntityDrv[0][1] =FaceEntityDrv[0][2] =0;
	
                for(j=0; j < p-2; j++){
	  
                    N[nshp]= FaceBlend * EdgeEntity[j];
                    dN[nshp][0]= FaceBlend * FaceEntityDrv[j][0] + FaceBlendDrv[0] * FaceEntity[j];
                    dN[nshp][1]= FaceBlend * FaceEntityDrv[j][1] + FaceBlendDrv[1] * FaceEntity[j];
                    dN[nshp][2]= FaceBlend * FaceEntityDrv[j][2] + FaceBlendDrv[2] * FaceEntity[j];
                    nshp++;
                }
            }
        }
    }
  
    return nshp;
}




