#include <FCMangle.h>

typedef double DARR2[2];

int lineIntPnt(int, DARR2**,double**);

#define symline FortranCInterface_GLOBAL_(symline, SYMLINE)
void symline(int *n1, double pt[][4], double wt[], int *err)
{
  double *lwt = 0;
  DARR2 *lpt = 0;
  int i,j;
  *err = lineIntPnt(*n1, &lpt, &lwt);
  for(i=0; i < *n1; i++){
    wt[i] = lwt[i];
    for(j=0; j < 2; j++)
      pt[i][j] = lpt[i][j];
  }
}

#include <stdio.h>

#define QpNon -1.000000000000000
#define Qp11   0.000000000000000
#define Qw1    2.000000000000000

#define Qp21  -0.577350269189626
#define Qp22   0.577350269189626
#define Qw2    1.000000000000000

#define Qp31  -0.774596669241483
#define Qp32   0.000000000000000
#define Qp33   0.774596669241483
#define Qw31   0.555555555555556
#define Qw32   0.888888888888889

/* typedef double DARR2[2] ; */

static double rst1[][2] = {{Qp11,QpNon}};
static double wt1[] = {Qw1};

static double rst2[][2] = {{Qp21,QpNon},{Qp22,QpNon}};
static double wt2[] = {Qw2,Qw2};

static double rst3[][2] = {{Qp31,QpNon},{Qp32,QpNon},{Qp33,QpNon}};
static double wt3[] = {Qw31,Qw32,Qw31};

                             
int lineIntPnt(int nint, DARR2 **bcord, double **wt)
{
  int retval = 1 ;
  if( nint == 1 ) {*bcord = rst1 ; *wt = wt1; }
  else if( nint == 2 ){*bcord = rst2 ; *wt = wt2; }
  else if( nint == 3 ){*bcord = rst3 ; *wt = wt3; }
/*  else if( nint == 4 ){*bcord = rst4 ; *wt = wt4; }
  else if( nint == 5 ){*bcord = rst5 ; *wt = wt5; }*/
  else
  {
    fprintf(stderr,"\n%d integration points unsupported in symline.c; give {1,2,3}\n",nint);
    retval = 0;
  }
  return retval ;
}



