int GaussLegendreTri(int,int,double [][3],double[]);

#ifdef sun4_5
gltri_(int *n1, int *n2, double pt[][3], double wt[], int *npt)
#endif
#ifdef LINUX 
gltri_(int *n1, int *n2, double pt[][3], double wt[], int *npt)
#endif
#ifdef ibm
gltri(int *n1, int *n2, double pt[][3], double wt[], int *npt)
#endif
#ifdef sgi
void gltri_(int *n1, int *n2, double pt[][3], double wt[], int *npt)
#endif
#ifdef decalp
void gltri_(int *n1, int *n2, double pt[][3], double wt[], int *npt)
#endif
#ifdef intel
void GLTRI(int *n1, int *n2, double pt[][3], double wt[], int *npt)
#endif
{
  *npt = GaussLegendreTri(*n1,*n2,pt,wt);
}
