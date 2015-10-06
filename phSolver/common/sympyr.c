#include <stdlib.h>

#include <FCMangle.h>
#define sympyr FortranCInterface_GLOBAL_(sympyr, SYMPYR)

typedef double DARR4[4];

int pyrIntPnt(int,DARR4**,double **);

void sympyr(int *n1, double pt[][4], double wt[], int *err)
{
  double *lwt;
  DARR4 *lpt;
  int i,j;
  *err = pyrIntPnt(*n1, &lpt, &lwt);
  for(i=0; i < *n1; i++) {
    wt[i] = lwt[i];
    for(j=0; j <4; j++) 
      pt[i][j]=lpt[i][j];
  }
}

/*$Id$*/
#include <stdio.h>

/* these are the rule 2 int points and weights */

#define mpt455 -0.455341801261479548
#define ppt455  0.455341801261479548
#define mpt577 -0.577350269189625764
#define ppt577  0.577350269189625764
#define mpt122 -0.122008467928146216
#define ppt122  0.122008467928146216
#define ppt622  0.622008467928146212
#define ppt044  0.0446581987385204510

/* these are the rule 3 int points and weights */

#define zero   0.0000000000000000
#define mpt687 -0.6872983346207417
#define ppt687  0.6872983346207417
#define mpt774 -0.7745966692414834
#define ppt774  0.7745966692414834
#define mpt387 -0.3872983346207416
#define ppt387  0.3872983346207416
#define mpt087 -0.08729833462074170
#define ppt087  0.08729833462074170
#define ppt134 0.1349962850858610
#define ppt215 0.2159940561373775
#define ppt345 0.3455904898198040
#define ppt068 0.06858710562414265
#define ppt109 0.1097393689986282
#define ppt175 0.1755829903978052
#define ppt002 0.002177926162424265
#define ppt003 0.003484681859878818
#define ppt005 0.005575490975806120

/* these are the rule 4 integration points */
/*NOT FIXED YET*/
#define Qp41  -0.801346029378026
#define Qp42  -0.316375532749811
#define Qp43   0.316375532749811
#define Qp44   0.801346029378026
#define Qp45  -0.8611363116

#define Qp46  -0.576953166749811
#define Qp47  -0.227784076803672
#define Qp48   0.227784076803672
#define Qp49   0.576953166749811
#define Qp410 -0.3399810436

#define Qp411 -0.284183144850188
#define Qp412 -0.112196966796327
#define Qp413  0.112196966796327
#define Qp414  0.284183144850188
#define Qp415  0.3399810436

#define Qp416 -0.0597902822219738
#define Qp417 -0.0236055108501886
#define Qp418  0.0236055108501886
#define Qp419  0.0597902822219738
#define Qp420  0.8611363116

/* Rule 1*/

static double  rstw1[][4] = {
  { 0.0, 0.0, 0.0, 0.0 }
};

static double twt1[] = { 8.0000000 };


/* Rule 2*/


static double rstw8[][4] = {
  {mpt455,mpt455,mpt577,0.0},
  {ppt455,mpt455,mpt577,0.0},
  {mpt455,ppt455,mpt577,0.0},
  {ppt455,ppt455,mpt577,0.0},
  {mpt122,mpt122,ppt577,0.0},
  {ppt122,mpt122,ppt577,0.0},
  {mpt122,ppt122,ppt577,0.0},
  {ppt122,ppt122,ppt577,0.0}
};

static double twt8[] = {ppt622,ppt622,ppt622,ppt622,
                        ppt044,ppt044,ppt044,ppt044};

/* Rule 3*/

static double rstw27[][4] = {
  { mpt687,  mpt687,  mpt774,  0.0 },
  {   zero,  mpt687,  mpt774,  0.0 },
  { ppt687,  mpt687,  mpt774,  0.0 },
  { mpt687,    zero,  mpt774,  0.0 },
  {   zero,    zero,  mpt774,  0.0 },
  { ppt687,    zero,  mpt774,  0.0 },
  { mpt687,  ppt687,  mpt774,  0.0 },
  {   zero,  ppt687,  mpt774,  0.0 },
  { ppt687,  ppt687,  mpt774,  0.0 },
  { mpt387,  mpt387,    zero,  0.0 },
  {   zero,  mpt387,    zero,  0.0 },
  { ppt387,  mpt387,    zero,  0.0 },
  { mpt387,    zero,    zero,  0.0 },
  {   zero,    zero,    zero,  0.0 },
  { ppt387,    zero,    zero,  0.0 },
  { mpt387,  ppt387,    zero,  0.0 },
  {   zero,  ppt387,    zero,  0.0 },
  { ppt387,  ppt387,    zero,  0.0 },
  { mpt087,  mpt087,  ppt774,  0.0 },
  {   zero,  mpt087,  ppt774,  0.0 },
  { ppt087,  mpt087,  ppt774,  0.0 },
  { mpt087,    zero,  ppt774,  0.0 },
  {   zero,    zero,  ppt774,  0.0 },
  { ppt087,    zero,  ppt774,  0.0 },
  { mpt087,  ppt087,  ppt774,  0.0 },
  {   zero,  ppt087,  ppt774,  0.0 },
  { ppt087,  ppt087,  ppt774,  0.0 }
};


static double twt27[] = 
{0.1349962850858610,   0.2159940561373775,   0.1349962850858610,
 0.2159940561373775,   0.3455904898198040,   0.2159940561373775,
 0.1349962850858610,   0.2159940561373775,   0.1349962850858610,
 0.06858710562414265,  0.1097393689986282,   0.06858710562414265,
 0.1097393689986282,   0.1755829903978052,   0.1097393689986282,
 0.06858710562414265,  0.1097393689986282,   0.06858710562414265,
 0.002177926162424265, 0.003484681859878818, 0.002177926162424265,
 0.003484681859878822, 0.005575490975806120, 0.003484681859878825,
 0.002177926162424265, 0.003484681859878822, 0.002177926162424265};

/* Rule 4 */

static double rstw64[][4] = {
  { Qp41,  Qp41,  Qp45,  0.0 },
  { Qp42,  Qp41,  Qp45,  0.0 },
  { Qp43,  Qp41,  Qp45,  0.0 },
  { Qp44,  Qp41,  Qp45,  0.0 },
  { Qp41,  Qp42,  Qp45,  0.0 },
  { Qp42,  Qp42,  Qp45,  0.0 },
  { Qp43,  Qp42,  Qp45,  0.0 },
  { Qp44,  Qp42,  Qp45,  0.0 },
  { Qp41,  Qp43,  Qp45,  0.0 },
  { Qp42,  Qp43,  Qp45,  0.0 },
  { Qp43,  Qp43,  Qp45,  0.0 },
  { Qp44,  Qp43,  Qp45,  0.0 },
  { Qp41,  Qp44,  Qp45,  0.0 },
  { Qp42,  Qp44,  Qp45,  0.0 },
  { Qp43,  Qp44,  Qp45,  0.0 },
  { Qp44,  Qp44,  Qp45,  0.0 },             	
  { Qp46,  Qp46,  Qp410,  0.0 },
  { Qp47,  Qp46,  Qp410,  0.0 },
  { Qp48,  Qp46,  Qp410,  0.0 },
  { Qp49,  Qp46,  Qp410,  0.0 },
  { Qp46,  Qp47,  Qp410,  0.0 },
  { Qp47,  Qp47,  Qp410,  0.0 },
  { Qp48,  Qp47,  Qp410,  0.0 },
  { Qp49,  Qp47,  Qp410,  0.0 },
  { Qp46,  Qp48,  Qp410,  0.0 },
  { Qp47,  Qp48,  Qp410,  0.0 },
  { Qp48,  Qp48,  Qp410,  0.0 },
  { Qp49,  Qp48,  Qp410,  0.0 },
  { Qp46,  Qp49,  Qp410,  0.0 },
  { Qp47,  Qp49,  Qp410,  0.0 },
  { Qp48,  Qp49,  Qp410,  0.0 },
  { Qp49,  Qp49,  Qp410,  0.0 },
  { Qp411,  Qp411,  Qp415,  0.0 },
  { Qp412,  Qp411,  Qp415,  0.0 },
  { Qp413,  Qp411,  Qp415,  0.0 },
  { Qp414,  Qp411,  Qp415,  0.0 },
  { Qp411,  Qp412,  Qp415,  0.0 },
  { Qp412,  Qp412,  Qp415,  0.0 },
  { Qp413,  Qp412,  Qp415,  0.0 },
  { Qp414,  Qp412,  Qp415,  0.0 },
  { Qp411,  Qp413,  Qp415,  0.0 },
  { Qp412,  Qp413,  Qp415,  0.0 },
  { Qp413,  Qp413,  Qp415,  0.0 },
  { Qp414,  Qp413,  Qp415,  0.0 },
  { Qp411,  Qp414,  Qp415,  0.0 },
  { Qp412,  Qp414,  Qp415,  0.0 },
  { Qp413,  Qp414,  Qp415,  0.0 },
  { Qp414,  Qp414,  Qp415,  0.0 },
  { Qp416,  Qp416,  Qp420,  0.0 },
  { Qp417,  Qp416,  Qp420,  0.0 },
  { Qp418,  Qp416,  Qp420,  0.0 },
  { Qp419,  Qp416,  Qp420,  0.0 },
  { Qp416,  Qp417,  Qp420,  0.0 },
  { Qp417,  Qp417,  Qp420,  0.0 },
  { Qp418,  Qp417,  Qp420,  0.0 },
  { Qp419,  Qp417,  Qp420,  0.0 },
  { Qp416,  Qp418,  Qp420,  0.0 },
  { Qp417,  Qp418,  Qp420,  0.0 },
  { Qp418,  Qp418,  Qp420,  0.0 },
  { Qp419,  Qp418,  Qp420,  0.0 },
  { Qp416,  Qp419,  Qp420,  0.0 },
  { Qp417,  Qp419,  Qp420,  0.0 },
  { Qp418,  Qp419,  Qp420,  0.0 },
  { Qp419,  Qp419,  Qp420,  0.0 }
};
  
#ifdef __cplusplus
extern "C" {
#endif

  int pyrIntPnt(int nint, DARR4 **bcord, double **wt)
{
  int retval = 1;
  int i,j,k,l;
  DARR4 *rstw;
  double *twt;

  /* Rule 4 & 5*/
  
  double bp5[]={-0.9061798459,-0.5384693101,0,0.5384693101,0.9061798459};
  double bw4[]={0.3478548446,0.6521451548,0.6521451548,0.3478548446};
  double bw5[]={0.2369268850,0.4786286708,0.5019607843,0.4786286708,0.2369268850};

  if( nint == 1){ *bcord = rstw1; *wt = twt1;}

  else if( nint == 8 ){*bcord = rstw8 ; *wt = twt8; }

  else if( nint == 27 ) {*bcord = rstw27 ; *wt = twt27; }

  else if( nint == 64 ) {

 /*     rstw = (DARR4 *)malloc(64*sizeof(DARR4)); */
    twt   = (double *)malloc(64*sizeof(double));

    i=j=k=0;

    for(l=0;l<64;l++){
      twt[l]=bw4[i]*bw4[j]*bw4[k];
      i++;
      if(i == 4){
	i=0 ; j++ ;
	if(j == 4) { j=0; k++;}
      }
    }
    *bcord = rstw64; *wt = twt; 

  } else if( nint == 125) {
    
    rstw = (DARR4 *)malloc(125*sizeof(DARR4));
    twt   = (double *)malloc(125*sizeof(double));
    
    i=j=k=0;
    
    for(l=0;l<125;l++){
      rstw[l][0]=bp5[i];
      rstw[l][1]=bp5[j];
      rstw[l][2]=bp5[k];
      rstw[l][3]=0.0;
      twt[l]=bw5[i]*bw5[j]*bw5[k];
      i++;
      if(i == 5){
	i=0 ; j++ ;
	if(j == 5) { j=0; k++;}
      }
    }
    *bcord = rstw; *wt = twt; 

  } else {

    fprintf(stderr,"\n%d integration points unsupported; give {8,27,64,125,.}\n",nint);
    retval = 0;
  }
  return retval ;
}

#ifdef __cplusplus
}
#endif
