/*-------------------------------------------------------------------------
  Scientific Computation Research Center, RPI, Troy NY
  (C) Copyright 1995, RPI-SCOREC
 
  Project  : shapeFuntions
  Author(s): Saikat Dey
  Creation : Oct., 95
  Modifi.  : 
  Function :
             evaluate the blend funtion for a given mesh entity over another
             mesh entity
-------------------------------------------------------------------------*/

#include <math.h>
#include "shapeFuncInternals.h"

#ifdef __cplusplus
extern "C" {
#endif

double V_blendOnEntity(int vid, int etype, double *L) {
  /* blend a vertex mode on a mesh entity */
  /* vid = local vertex index 
     etype = entity type
     L = par. coords
  */
  if( etype == Sedge )
    return V_blendIndexedOnEdge(vid,L);
  else if(etype == Stri || etype == Stet)
    return V_blendIndexed(vid,L);
  else
    return 0.0e0;
}

double V_blendIndexed(int i, double *L) {
  return L[i] ;
}

double V_blendIndexedOnEdge(int i, double *L) {
  if( i == 0 )
    return 0.5*(1.0-L[0]);
  else
    return 0.5*(1.0+L[0]);
}

double E_blendOnFace(int eindex[], int etype, double *L) {
  /* blend a vertex mode on a mesh face */
  /* vid = local vertex index 
     etype = entity type
     L = par. coords
  */
  if( etype == Stri)
    return F_edgeBlendTri(eindex, L);
  else if(etype == Squad)
    return F_edgeBlendQuad(eindex, L);
  else
    return 0.0e0;
}

double F_edgeBlendTri(int index[2], double *L) {
  return -2.0*L[index[0]]*L[index[1]];
}

double F_edgeBlendQuad(int *index, double *L) {
  return 0.0;
}

double E_blendOnRegion(int eindex[], int etype, double *L) {
  /* blend a mesh edge mode on a tetra. region */
 if( etype == Stet ) 
    return R_edgeBlendTet(eindex, L);   
 else
  return 0.0;
}

double R_edgeBlendTet(int index[2], double *L) {
  return -2.0*L[index[0]]*L[index[1]];
}

double F_blendOnRegion(int index[], int etype, double *L) {
  /* blend a face mode on a tet. region */
  if( etype == Stet ) {
    return L[index[0]]*L[index[1]]*L[index[2]] ;
  } else
    return 0.0;
}

#ifdef __cplusplus
}
#endif
