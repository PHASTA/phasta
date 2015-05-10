/* fortran wrapper for C function TetShapeAndDrv which returns
   the shape functions and their derivatives for tets
   */

int TriShapeAndDrv(int p,double par[3],double N[],double dN[][3]);

#ifdef sun4_5
shptri_(int *p, double par[], double N[], double dN[][3])
#elif LINUX
shptri_(int *p, double par[], double N[], double dN[][3])
#elif ibm
shptri(int *p, double par[], double N[], double dN[][3])
#elif sgi
void shptri_(int *p, double par[], double N[], double dN[][3])
#elif decalp
void shptri_(int *p, double par[], double N[], double dN[][3])
#elif intel
void  SHPTRI(int *p, double par[], double N[], double dN[][3])
#endif
{
 TriShapeAndDrv(*p,par,N,dN);

}
