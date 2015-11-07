/* fortran wrapper for C function HexShapeAndDrv which returns
   the shape functions and their derivatives for tets
   */
#include <FCMangle.h>

int HexShapeAndDrv(int p,double par[3],double N[],double dN[][3]);

#define shphex FortranCInterface_GLOBAL_(shphex, SHPHEX)
void shphex(int *p, double par[], double N[], double dN[][3])
{
  
  HexShapeAndDrv(*p,par,N,dN);

}
