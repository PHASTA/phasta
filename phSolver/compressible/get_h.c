/* fortran wrapper for C program length which returns the length
   of an element in a given direction.
   */
#include "vec_func.h"


double length(double [][3],double []);
int intersect(double x[4][3],double u[],double points[4][3] );
int isPntInTri(double [][3],double [], double [], double );

#ifdef sun4_5
get_h_(double xtmp[][4], double u[], double *h)
#endif
#ifdef LINUX
get_h_(double xtmp[][4], double u[], double *h)
#endif
#ifdef ibm
get_h(double xtmp[][4], double u[], double *h)
#endif
#ifdef sgi
void get_h(double xtmp[][4], double u[], double *h)
#endif
#ifdef decalp
void get_h(double xtmp[][4], double u[], double *h)
#endif
#ifdef intel
void GET_H(double xtmp[][4], double u[], double *h)
#endif
{
  double x[4][3];
  int i,j;
  /*recall that arrays in fortran and C are transposed */
  for(i=0; i < 4; i++)
    for(j=0; j < 3; j++)
      x[i][j] = xtmp[j][i];

  *h = length(x,u);
}
  
double length(double x[4][3],double u[] ){
  int nint;
  double tmp[3],he;
  double points[4][3];

  /* compute the intersection points. nint is returned
     as the number of faces that this line intersected
     within the tet */
  nint = intersect(x,u,points);

  /* if nint is larger than one, then the vector intersected
     an edge (nint=3) or corner (nint=4) of the tet. In this
     case, we need to make sure that we are computing the
     length between two different points. */

  if(nint > 1){			/* vecEqual returns 1 if
				   these two vectors are the same.
				   */
    if( vecEqual(points[0],points[1]) ){ 
      if( vecEqual(points[0],points[2]) ){
	diffVt(points[0],points[3],tmp);
      }
      else{
	diffVt(points[0],points[2],tmp);
      }
    }
    else{
      diffVt(points[0],points[1],tmp);
    }
  }
  else{
    diffVt(points[0],points[1],tmp);
  }

  /* compute the length of the element */
  he = vecMag(tmp);
  return he;
}

/* 
   This function returns the intersection points of a line
   through the centroid of a tet in a given direction.
*/


/* u_line is a function that defines the parametric equation
   of a line given a point on the lineg, a direction, and a
   parameter value */

void u_line(double s, double xc[3], double udir[3], double pnt[3]);

/*
  All points which intersect the tet are returned, however
  only two of them need be used. The value of the function on
  exit is the number of points of intersection, ie. the number
  of faces not parallel to udir */

int intersect(double x[4][3],double u[],double points[4][3] ){
  int pos = 0;
  double tol = 0.000001;
  double eps = 0.000001;
  int face,i,j,k;
  double xc[3];
  double xbar[3][3];
  double xref[3],v1[3],v2[3],normal[3];
  /* define local vertex numbers of four faces of the tet */
  int f_index[4][3] = {{0,1,2},{1,3,2},{0,2,3},{0,3,1}};
  double n_dot_u;
  double xr_xc[3], s, pt[3];

  /* compute centroid coordinates */
  for(i=0; i<3; i++)
    xc[i] = 0.25*(x[0][i]+x[1][i]+x[2][i]+x[3][i]);

  /* loop over each face */
  for( face= 0 ; face< 4; face++ ){
    /* compute coordinates of this face */
    for( j= 0 ; j<3 ; j++)
      for( k= 0; k< 3; k++)
	xbar[j][k] = x[f_index[face][j]][k];
    
    /* set reference coordinate */
    if( face != 2 ){
      for( i= 0; i<3; i++ )
	xref[i] = xbar[0][i];
    }
    else{			/* face 2 is the only face that
				 can't use xbar[0] as a reference */
      for( i= 0; i<3; i++ )
	xref[i] = xbar[1][i];
    }

    /* compute two vectors which describe the plane of this
       face */
    for( i=0; i< 3; i++){
      v1[i] = xbar[1][i] - xbar[0][i];
      v2[i] = xbar[2][i] - xbar[0][i];
    }
    /* get the normal to this plane */
    crossProd(v1,v2,normal);

    n_dot_u = dotProd(normal,u);
    
    /* if u is not parallel to this face, then find
       the point of intersection */
    if( fabs(n_dot_u) > eps ){
      diffVt( xref, xc, xr_xc );

      /* s is the parameter value where the line parallel to
	 to the vector u intersects this face*/
      s = dotProd(normal,xr_xc)/n_dot_u;

      u_line(s,xc,u,pt); /* find intersection point */
      
      /* if this point lies on the tet, then save it */
      if( isPntInTri(xbar,normal,pt,tol) ){
	for(i=0; i<3; i++)
	  points[pos][i] = pt[i];
	pos += 1;
      }
    }
    else{			/* do nothing */

    }
  }

  return pos;
}

/* parametric equation of the line passing through xc
   in the direction of u.
*/

void u_line(double s, double xc[3], double udir[3], double pnt[3]){
  int i;
  for(i=0; i<3; i++)
    pnt[i] = xc[i]+s*udir[i];
}


void diffVt(double vec1[3],double vec2[3], double ans[3]){
  ans[0] = vec1[0]-vec2[0];
  ans[1] = vec1[1]-vec2[1];
  ans[2] = vec1[2]-vec2[2];
}

void crossProd(double a[3], double b[3], double ans[3]){
  ans[0] = a[1]*b[2]-b[1]*a[2];
  ans[1] = a[2]*b[0]-b[2]*a[0];
  ans[2] = a[0]*b[1]-b[0]*a[1];
}
  
double dotProd(double a[3], double b[3]){
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double vecMag(double a[3]){
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

int vecEqual(double a[3], double b[3]){
  return (a[0]==b[0] && a[1]==b[1] && a[2]==b[2]);
}
    
/*
Determines if a point is inside or out of a triangle
Pt will be declared in triangle if distance to side < tol
if u want the tolerancing to be relative to 't_xyz',
pass tol (absolute tolerance) * size of 't_xyz'
*/

int isPntInTri(
 double t_xyz[3][3],
 double t_normal[3], /* Normal to triangle (not normed) */
 double p_xyz[3],
 double tol
)

{

 double cross_p[3];
 double d;
 int side;
 int status;
 double cross_p_sqnorm;
 double dnbr;
 double dtol;
 double t_vec[3];
 double t_vec_norm;
 double tp_vec[3];
 double tp_vec_sqnorm;

 status= 1;

 for ( side= 0 ; side< 3 ; side++ ) {

    diffVt(t_xyz[(side+1)%3],t_xyz[side],t_vec);
    crossProd(t_normal,t_vec,cross_p);
    cross_p_sqnorm= dotProd(cross_p,cross_p);

    diffVt(p_xyz,t_xyz[side],tp_vec);

    d= dotProd(cross_p,tp_vec);
    dnbr= d * d;
    dtol= tol*tol*cross_p_sqnorm;

    if ( dnbr < dtol ) {

       /*
       Pt is considered to be on the side (extended to infinity)
       */

       /*
       Make sure distance to point 0 or 1 less than ...
       */

       t_vec_norm= vecMag(t_vec);
       dtol= t_vec_norm + tol;
       dtol *= dtol;

       tp_vec_sqnorm= dotProd(tp_vec,tp_vec);
       if ( tp_vec_sqnorm > dtol ) {
          status= 0;
          break;
       }

       diffVt(p_xyz,t_xyz[(side+1)%3],tp_vec); 
       tp_vec_sqnorm= dotProd(tp_vec,tp_vec); 
       if ( tp_vec_sqnorm > dtol ) { 
          status= 0; 
          break; 
       }

       continue;
    }

    if ( d > 0.0 ) {

       /*
       Pt is strictly in the interior of the side
       */
    }
    else if ( d < 0.0 ) {

       /*
       Pt is strictly in the outisde of the side
       */

       status= 0;
       break;
    }

 }

 return status;

}





	
  
  
  

	 
