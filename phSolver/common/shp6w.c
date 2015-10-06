/* fortran wrapper for C function WedgeShapeAndDrv which returns
   the shape functions and their derivatives for wedge
   */
#include <FCMangle.h>
#define shp6w FortranCInterface_GLOBAL_(shp6w, SHP6W)

int WedgeShapeAndDrv(int p,double par[3],double N[],double dN[][3]);

void shp6w(int *p, double par[], double N[], double dN[][3])
{
  WedgeShapeAndDrv(*p,par,N,dN);

}








