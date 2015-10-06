/* fortran wrapper for C function PyrShapeAndDrv which returns
   the shape functions and their derivatives for pyramid
   */
#include <FCMangle.h>
#define shppyr FortranCInterface_GLOBAL_(shppyr,SHPPYR)

int PyrShapeAndDrv(int p,double par[3],double N[],double dN[][3]);
/* p:      the order of the element shape function
   par[3]: xi[3]
   N[]:    shape function
   dN[][3]:derivative of shape function
   */

void shppyr(int *p, double par[3], double N[], double dN[][3])
{

  PyrShapeAndDrv(*p,par,N,dN);
  
}








