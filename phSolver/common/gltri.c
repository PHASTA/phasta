#include <FCMangle.h>
int GaussLegendreTri(int,int,double [][3],double[]);

#define gltri FortranCInterface_GLOBAL_(gltri, GLTRI)
void gltri(int *n1, int *n2, double pt[][3], double wt[], int *npt)
{
  *npt = GaussLegendreTri(*n1,*n2,pt,wt);
}
