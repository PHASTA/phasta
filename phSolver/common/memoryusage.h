/*  Primary interface for the Phasta Binary read and write routines these*/
/*  functions are 'C' callable.( All arguments have been kept as pointers to*/
/*  facilitate calling from Fortran )*/
/*  Michel Rasquin  Spring 2013, inspired from phastaIO.h*/
#ifndef __MEMORYUSAGE_H__
#define __MEMORYUSAGE_H__

#include <FCMangle.h>

#define initmemstat  FortranCInterface_GLOBAL_(initmemstat,INITMEMSTAT)
#define printmeminfo  FortranCInterface_GLOBAL_(printmeminfo,PRINTMEMINFO)

#if defined (__cplusplus)
extern "C" {
#endif

  void initmemstat( void);
  void printmeminfo( const char *msg);

#ifdef __cplusplus
} // end of extern "C".

#endif // __cplusplus

#endif // __MEMORYUSAGE_H__
