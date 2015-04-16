/* Functions to generate the Legandre Polynomials and their derivatives */

extern "C" double LP(int j, double x);
extern "C" double LPdrv(int j, double x);
extern "C" double phi(int p, double x);
extern "C" double phiDrv(int p,double x);
extern "C" int HexShapeAndDrv(int p, double par[3], double N[], double
			      dN[][3]);
extern "C" int WedgeShapeAndDrv(int p, double Inputpar[3], double N[], double 
				dN[][3]);
extern "C" int PyrShapeAndDrv (int p, double Inputpar[3], double N[], double 
				dN[][3]);

/* Blending functions */

extern "C" double Line_eB(double xi1);
extern "C" double dLEBdxi1(double xi1);
extern "C" double dLEBdxi2(double xi1);
extern "C" double dLEBdxi3(double xi1);

extern "C" double Quad_eB(double xi1, double xi2, int sign);
extern "C" double dQEBdxi1(double xi1, double xi2, int sign);
extern "C" double dQEBdxi2(double xi1, double xi2, int sign);
extern "C" double dQEBdxi3(double xi1, double xi2, int sign);

extern "C" double Quad_fB(double xi1, double xi2);
extern "C" double dQFBdxi1(double xi1, double xi2);
extern "C" double dQFBdxi2(double xi1, double xi2);
extern "C" double dQFBdxi3(double xi1, double xi2);

extern "C" double Hex_eB(double xi[3], int sign2, int sign3);
extern "C" double dHEBdxi1(double xi[3], int sign2, int sign3);
extern "C" double dHEBdxi2(double xi[3], int sign2, int sign3);
extern "C" double dHEBdxi3(double xi[3], int sign2, int sign3);

extern "C" double Hex_fB(double xi[3], int sign3);
extern "C" double dHFBdxi1(double xi[3], int sign3);
extern "C" double dHFBdxi2(double xi[3], int sign3);
extern "C" double dHFBdxi3(double xi[3], int sign3);

/* Entity Level functions */

extern "C" int mesh_edge(double xi1,int gOrd[3], int p, double* entfn,
			 double** edrv);
extern "C" int quad_face(double xi[3], int gOrd[3], int p, double*
			 entfn, double** edrv); 
extern "C" int hex_regn(double xi[3], int p, double*
			entfn, double** edrv);

